/* dd/internal_fma.c
 * 
 * Copyright (C) 2026 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* this module implements double-double arithmetic using C99 fma() */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>

#include "internal.h"

static gsl_dd two_prod_fma(const double a, const double b);        /* * */
static gsl_dd dd_mul_fma(const gsl_dd a, const gsl_dd b);          /* * */
static gsl_dd dd_mul_d_fma(const gsl_dd a, const double b);        /* * */
static gsl_dd dd_div_fma(const gsl_dd a, const gsl_dd b);          /* * */
static gsl_dd dd_div_d_fma(const gsl_dd a, const double b);        /* * */

/*
 * Error-free transformation of a*b into a leading double and its
 * rounding residual using fused multiply-add.
 *
 * This is the Fast2Mult algorithm (Algorithm 2) of:
 *
 *   M. Joldes, J.-M. Muller, and V. Popescu,
 *   "Tight and Rigorous Error Bounds for Basic Building Blocks of
 *   Double-Word Arithmetic," ACM Transactions on Mathematical
 *   Software, 44(2), 2017.
 *
 * In the absence of underflow,
 *
 *   r.hi + r.lo = a*b
 *
 * exactly.
 */
static gsl_dd
two_prod_fma(const double a, const double b)
{
  gsl_dd r;

  r.hi = a * b;
  r.lo = fma(a, b, -r.hi);

  return r;
}

/*
 * Multiply two double-double numbers using FMA
 *
 * This is Algorithm 12 (DWTimesDW3) of:
 *
 *   M. Joldes, J.-M. Muller, and V. Popescu,
 *   "Tight and Rigorous Error Bounds for Basic Building Blocks of
 *   Double-Word Arithmetic," ACM Transactions on Mathematical
 *   Software, 44(2), 2017.
 */
static gsl_dd
dd_mul_fma(const gsl_dd a, const gsl_dd b)
{
  gsl_dd p = two_prod_fma(a.hi, b.hi);
  double e;

  e = fma(a.hi, b.lo, a.lo * b.lo);
  e += a.lo * b.hi;
  e += p.lo;

  return fast_two_sum(p.hi, e);
}

/*
 * Multiply a double-double number by a double using FMA.
 *
 * This is Algorithm 9 (DWTimesFP3) of:
 *
 *   M. Joldes, J.-M. Muller, and V. Popescu,
 *   "Tight and Rigorous Error Bounds for Basic Building Blocks of
 *   Double-Word Arithmetic," ACM Transactions on Mathematical
 *   Software, 44(2), 2017.
 */
static gsl_dd
dd_mul_d_fma(const gsl_dd a, const double b)
{
  gsl_dd c = two_prod_fma(a.hi, b);
  return fast_two_sum(c.hi, fma(a.lo, b, c.lo));
}

#if 0
/*
 * Divide two double-double numbers using FMA.
 *
 * Algorithm 17 (DWDivDW2) of:
 *
 *   M. Joldes, J.-M. Muller, and V. Popescu,
 *   "Tight and Rigorous Error Bounds for Basic Building Blocks of
 *   Double-Word Arithmetic," ACM Transactions on Mathematical
 *   Software, 44(2), 2017.
 *
 * The DWTimesFP operation uses Algorithm 9 (DWTimesFP3).
 */
static gsl_dd
dd_div_fma(const gsl_dd a, const gsl_dd b)
{
  double q1;
  gsl_dd p;
  double pi_h, delta_l, delta, q2;

  q1 = a.hi / b.hi;

  p = dd_mul_d_fma(b, q1);

  pi_h = a.hi - p.hi;
  delta_l = a.lo - p.lo;
  delta = pi_h + delta_l;

  q2 = delta / b.hi;

  return fast_two_sum(q1, q2);
}
#else
/*
 * Divide two double-double numbers using FMA
 *
 * This is Algorithm 18 (DWDivDW3) of:
 *
 *   M. Joldes, J.-M. Muller, and V. Popescu,
 *   "Tight and Rigorous Error Bounds for Basic Building Blocks of
 *   Double-Word Arithmetic," ACM Transactions on Mathematical
 *   Software, 44(2), 2017.
 *
 * The algorithm uses:
 *
 *   Algorithm 9  (DWTimesFP3)
 *   Algorithm 4  (DWPlusFP)
 *   Algorithm 12 (DWTimesDW3)
 *
 * Sep 5 2026: Algorithm 17 above is faster (lower operation count),
 * but we select this Algorithm 18 as the default since it is a bit
 * higher accuracy.
 */
static gsl_dd
dd_div_fma(const gsl_dd x, const gsl_dd y)
{
  double t_h;
  double r_h;
  double r_l;
  gsl_dd e;
  gsl_dd delta;
  gsl_dd m;

  t_h = 1.0 / y.hi;

  r_h = fma(-y.hi, t_h, 1.0);
  r_l = -(y.lo * t_h);

  e = fast_two_sum(r_h, r_l);

  delta = dd_mul_d_fma(e, t_h);

  m = dd_add_d(delta, t_h);

  return dd_mul_fma(x, m);
}
#endif

/*
 * Divide a double-double number by a double using FMA.
 *
 * This is Algorithm 15 (DWDivFP3) of:
 *
 *   M. Joldes, J.-M. Muller, and V. Popescu,
 *   "Tight and Rigorous Error Bounds for Basic Building Blocks of
 *   Double-Word Arithmetic," ACM Transactions on Mathematical
 *   Software, 44(2), 2017.
 */
static gsl_dd
dd_div_d_fma(const gsl_dd a, const double b)
{
  double q1;
  gsl_dd p;
  double delta_h, delta_t, delta, q2;

  q1 = a.hi / b;

  p = two_prod_fma(q1, b);

  delta_h = a.hi - p.hi;
  delta_t = delta_h - p.lo;
  delta = delta_t + a.lo;

  q2 = delta / b;

  return fast_two_sum(q1, q2);
}

#undef DD_SUFFIX
#undef DD_FUNC
#define DD_SUFFIX     _fma
#define DD_FUNC(a)    DD_CAT(a, DD_SUFFIX)
#include "internal_exp_source.c"
