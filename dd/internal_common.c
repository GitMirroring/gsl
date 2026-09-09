/* dd/internal_common.c
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

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>

#include "internal.h"

/*
 * Error-free transformation of a+b into a leading double and its
 * rounding residual.
 */
static gsl_dd
two_sum(const double a, const double b)
{
  gsl_dd r;
  double s = a + b, bb = s - a;
  r.hi = s;
  r.lo = (a - (s - bb)) + (b - bb);
  return r;
}

/*
 * Error-free transformation of a+b into a leading double and its
 * rounding residual, assuming |a| >= |b|.
 *
 * This is the Fast2Sum algorithm.
 */
static gsl_dd
fast_two_sum(const double a, const double b)
{
  gsl_dd r;
  r.hi = a + b;
  r.lo = b - (r.hi - a);
  return r;
}

/*
 * Add two double-double numbers using accurate double-word addition
 * and renormalize the result.
 *
 * This is the accurate DWPlusDW algorithm.
 */
static gsl_dd
dd_add(const gsl_dd a, const gsl_dd b)
{
  gsl_dd s, t, v;
  double c, w;

  s = two_sum(a.hi, b.hi);
  t = two_sum(a.lo, b.lo);

  c = s.lo + t.hi;
  v = fast_two_sum(s.hi, c);

  w = t.lo + v.lo;

  return fast_two_sum(v.hi, w);
}

/* add a double to a double-double number and renormalize the result */
static gsl_dd
dd_add_d(const gsl_dd a, const double b)
{
  gsl_dd s;
  s = two_sum(a.hi, b);
  s.lo += a.lo;
  return fast_two_sum(s.hi, s.lo);
}

/*
 * Subtract two double-double numbers using accurate double-word
 * subtraction and renormalize the result.
 */
static gsl_dd
dd_sub(const gsl_dd a, const gsl_dd b)
{
  gsl_dd s, t, v;
  double c, w;

  s = two_sum(a.hi, -b.hi);
  t = two_sum(a.lo, -b.lo);

  c = s.lo + t.hi;
  v = fast_two_sum(s.hi, c);

  w = t.lo + v.lo;

  return fast_two_sum(v.hi, w);
}

/*
 * Subtract a double from a double-double number and renormalize
 * the result.
 */
static gsl_dd
dd_sub_d(const gsl_dd a, const double b)
{
  gsl_dd s;

  s = two_sum(a.hi, -b);
  s.lo += a.lo;

  return fast_two_sum(s.hi, s.lo);
}

/*
 * Scale a double-double number by an integral power of two.
 *
 * Apply ldexp() independently to the high and low components.  This
 * preserves the represented value exactly provided the component
 * scalings do not overflow or underflow.
 */
static gsl_dd
dd_ldexp(const gsl_dd a, const int e)
{
  gsl_dd r;

  r.hi = ldexp(a.hi, e);
  r.lo = ldexp(a.lo, e);

  return r;
}

/*
 * Accumulate a binary64 value into a three-term expansion.
 *
 * Add q to the leading component with TwoSum, then propagate the
 * resulting residual through the second component.  Any residual from
 * that addition is accumulated in the third component.
 */
static void
ex3_add_d(dd_ex3 * x, const double q)
{
  gsl_dd s, t;

  s = two_sum(x->a, q);
  t = two_sum(x->b, s.lo);

  x->a = s.hi;
  x->b = t.hi;
  x->c += t.lo;
}

/*
 * Renormalize a three-term floating-point expansion.
 *
 * Fold the third component into the second with TwoSum, then propagate
 * the leading part of that result into the first component.  The
 * remaining residual is retained in the second component and the third
 * component is cleared.
 *
 * The resulting first two components are suitable for subsequent
 * compression to a normalized double-double value.
 */
static void
ex3_renorm(dd_ex3 * x)
{
  gsl_dd t, u;

  t = two_sum(x->b, x->c);
  u = two_sum(x->a, t.hi);

  x->a = u.hi;
  x->b = u.lo + t.lo;
  x->c = 0.0;
}

/*
 * Compress a three-term floating-point expansion to double-double.
 *
 * Combine the two trailing components with TwoSum, add their leading
 * contribution to the first component, and renormalize the resulting
 * two-word value with FastTwoSum.
 */
static gsl_dd
ex3_to_dd(const dd_ex3 x)
{
  gsl_dd t, u;

  t = two_sum(x.b, x.c);
  u = two_sum(x.a, t.hi);
  u.lo += t.lo;

  return fast_two_sum(u.hi, u.lo);
}

/*
 * Accumulate a binary64 value into a four-term expansion.
 *
 * Add q to the leading component with TwoSum, then propagate each
 * resulting residual through the second and third components.  Any
 * residual from the final addition is accumulated in the fourth
 * component.
 */
static void
ex4_add_d(dd_ex4 * x, const double q)
{
  gsl_dd s;

  s = two_sum(x->a, q);
  x->a = s.hi;
  s = two_sum(x->b, s.lo);
  x->b = s.hi;
  s = two_sum(x->c, s.lo);
  x->c = s.hi;
  x->d += s.lo;
}

/*
 * Compress a four-term floating-point expansion to double-double.
 *
 * Combine the components from smallest to largest with TwoSum,
 * retaining the residuals, then renormalize the resulting two-word
 * value with FastTwoSum.
 */
static gsl_dd
ex4_to_dd(const dd_ex4 x)
{
  gsl_dd s, t, u;

  s = two_sum(x.c, x.d);
  t = two_sum(x.b, s.hi);
  u = two_sum(x.a, t.hi);
  u.lo += t.lo + s.lo;

  return fast_two_sum(u.hi, u.lo);
}

/* initialize a four-component expansion to zero */
static int
ex4_zero(dd_ex4 * x)
{
  x->a = 0.0;
  x->b = 0.0;
  x->c = 0.0;
  x->d = 0.0;
  return GSL_SUCCESS;
}

#include "internal_exp.c"
#include "internal_log.c"
#include "internal_sqrt.c"
