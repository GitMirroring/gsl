/* dd/internal_p.c
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

/* this module implements double-double arithmetic portable for C89 systems */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>

#include "internal.h"

static void split_p(const double a, double * ahi, double * alo); /* * */
static gsl_dd two_prod_p(const double a, const double b);        /* * */

static gsl_dd dd_mul_p(const gsl_dd a, const gsl_dd b);          /* * */
static gsl_dd dd_mul_d_p(const gsl_dd a, const double b);        /* * */
static gsl_dd dd_div_p(const gsl_dd a, const gsl_dd b);          /* * */
static gsl_dd dd_div_d_p(const gsl_dd a, const double b);        /* * */

#if 0
static gsl_dd dd_sqrt_p(const gsl_dd a);
static gsl_dd dd_expm1_p(const gsl_dd z);
#endif

/*
 * Split a double into high and low parts for Dekker multiplication.
 *
 * This is Dekker/Veltkamp splitting with scaling for large operands
 * to prevent overflow in DD_SPLIT*a.
 */
static void
split_p(const double a, double * ahi, double * alo)
{
  const double split_thresh =
    6.696928794914170755245419009244910293412e299;

  double c;
  double m;
  double n;
  int e;

  if (fabs(a) > split_thresh)
    {
      /*
       * For large operands, a direct DD_SPLIT*a can be represented,
       * but scaling the resulting high part back by 2^28 can overflow
       * at the top of the binary64 range.  Instead, split the binary64
       * significand directly, retaining at most 26 significant bits in
       * the high part.
       */
      m = frexp(a, &e);
      if (a > 0.0)
        n = floor(ldexp(m, 26));
      else
        n = ceil(ldexp(m, 26));

      *ahi = ldexp(n, e - 26);
      *alo = a - *ahi;
    }
  else
    {
      c = DD_SPLIT * a;

      *ahi = c - (c - a);
      *alo = a - *ahi;
    }
}

/*
 * Error-free transformation of a*b into a leading double and its
 * rounding residual using Dekker splitting.
 *
 * This is the 2Prod algorithm (Algorithm 3) of:
 *
 *   M. Joldes, J.-M. Muller, and V. Popescu,
 *   "Tight and Rigorous Error Bounds for Basic Building Blocks of
 *   Double-Word Arithmetic," ACM Transactions on Mathematical
 *   Software, 44(2), 2017.
 *
 * The product splitting is based on Dekker's multiplication algorithm:
 *
 *   T. J. Dekker, "A Floating-Point Technique for Extending the
 *   Available Precision," Numerische Mathematik, 18, 1971.
 *
 * In the absence of overflow and underflow in the required
 * intermediate quantities,
 *
 *   r.hi + r.lo = a*b
 *
 * exactly.
 */
static gsl_dd
two_prod_p(const double a, const double b)
{
  gsl_dd r;
  double ah, al, bh, bl;

  split_p(a, &ah, &al);
  split_p(b, &bh, &bl);

  r.hi = a * b;
  r.lo = ((ah * bh - r.hi) + ah * bl + al * bh) + al * bl;

  return r;
}

/*
 * Multiply two double-double numbers using non-FMA double-word
 * multiplication.
 *
 * This is Algorithm 10 (DWTimesDW1) of:
 *
 *   M. Joldes, J.-M. Muller, and V. Popescu,
 *   "Tight and Rigorous Error Bounds for Basic Building Blocks of
 *   Double-Word Arithmetic," ACM Transactions on Mathematical
 *   Software, 44(2), 2017.
 *
 * The algorithm uses an error-free 2Prod of the leading components,
 * adds the two cross products, and renormalizes the result with
 * Fast2Sum.  The product of the two low components is omitted as a
 * higher-order contribution.
 */
static gsl_dd
dd_mul_p(const gsl_dd a, const gsl_dd b)
{
  gsl_dd p = two_prod_p(a.hi, b.hi);
  double e = a.hi * b.lo + a.lo * b.hi;
  return fast_two_sum(p.hi, p.lo + e);
}

/*
 * Multiply a double-double number by a double
 *
 * This is Algorithm 7 (DWTimesFP1) of:
 *
 *   M. Joldes, J.-M. Muller, and V. Popescu,
 *   "Tight and Rigorous Error Bounds for Basic Building Blocks of
 *   Double-Word Arithmetic," ACM Transactions on Mathematical
 *   Software, 44(2), 2017.
 */
static gsl_dd
dd_mul_d_p(const gsl_dd a, const double b)
{
  gsl_dd p = two_prod_p(a.hi, b);
  gsl_dd s = fast_two_sum(p.hi, a.lo * b);
  return fast_two_sum(s.hi, s.lo + p.lo);
}

/*
 * Divide two double-double numbers.
 *
 * This is Algorithm 17 (DWDivDW2) of:
 *
 *   M. Joldes, J.-M. Muller, and V. Popescu,
 *   "Tight and Rigorous Error Bounds for Basic Building Blocks of
 *   Double-Word Arithmetic," ACM Transactions on Mathematical
 *   Software, 44(2), 2017.
 */
static gsl_dd
dd_div_p(const gsl_dd a, const gsl_dd b)
{
  double q1;
  gsl_dd p;
  double pi_h, delta_l, delta, q2;

  q1 = a.hi / b.hi;

  p = dd_mul_d_p(b, q1);

  pi_h = a.hi - p.hi;
  delta_l = a.lo - p.lo;
  delta = pi_h + delta_l;

  q2 = delta / b.hi;

  return fast_two_sum(q1, q2);
}

/*
 * Divide a double-double number by a double.
 *
 * This is Algorithm 15 (DWDivFP3) of:
 *
 *   M. Joldes, J.-M. Muller, and V. Popescu,
 *   "Tight and Rigorous Error Bounds for Basic Building Blocks of
 *   Double-Word Arithmetic," ACM Transactions on Mathematical
 *   Software, 44(2), 2017.
 */
static gsl_dd
dd_div_d_p(const gsl_dd a, const double b)
{
  double q1;
  gsl_dd p;
  double delta_h, delta_t, delta, q2;

  q1 = a.hi / b;

  p = two_prod_p(q1, b);

  delta_h = a.hi - p.hi;
  delta_t = delta_h - p.lo;
  delta = delta_t + a.lo;

  q2 = delta / b;

  return fast_two_sum(q1, q2);
}

#if 0
/*
 * Compute the square root of a nonnegative double-double number using
 * a double-precision initial approximation followed by a compensated
 * residual correction.
 */
static gsl_dd
dd_sqrt_p(const gsl_dd a)
{
  gsl_dd yy, neg, rem, r;
  double av = a.hi + a.lo;
  double y;

  if (av <= 0.0)
    {
      gsl_dd r;
      r.hi = 0.0;
      r.lo = 0.0;
      return r;
    }

  y = sqrt(av);
  yy = two_prod_p(y, y);
  neg.hi = -yy.hi;
  neg.lo = -yy.lo;
  rem = dd_add(a, neg);
  r.hi = y;
  r.lo = (rem.hi + rem.lo) / (2.0*y);

  return r;
}

/*
 * Compute exp(z)-1 for a double-double argument z.
 *
 * For |z| < 1/2, evaluate the expm1 Taylor series directly in
 * double-double arithmetic to avoid cancellation.  For larger
 * arguments, evaluate exp(z) in double-double arithmetic and subtract
 * one.
 */
static gsl_dd
dd_expm1_p(const gsl_dd z)
{
  const double av = fabs(z.hi + z.lo);
  gsl_dd term;

  if (av < 0.5)
    {
      gsl_dd sum = z;
      int k;

      term = z;
      for (k = 2; k < 200; k++)
        {
          term = dd_mul_p(term, z);
          term = dd_div_d_p(term, (double)k);
          sum = dd_add(sum, term);

          if (fabs(term.hi) < 1.0e-32 * fabs(sum.hi))
            break;
        }

      return sum;
    }

  term.hi = -1.0; term.lo = 0.0;
  return dd_add(dd_exp_p(z), term);
}
#endif

#undef DD_SUFFIX
#undef DD_WRAP
#define DD_SUFFIX     _p
#define DD_WRAP(a)    DD_CAT(a, DD_SUFFIX)
#include "internal_exp_source.c"
#include "internal_log_source.c"
#include "internal_sqrt_source.c"
