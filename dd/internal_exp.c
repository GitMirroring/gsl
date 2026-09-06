/* dd/internal_exp.c
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

/* common routines for double-double dd_exp() */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>

#include "internal.h"

/* three-term expansion used during range reduction */
typedef struct
{
  double a;
  double b;
  double c;
} dd_exp3;

static gsl_dd mul_pwr2_extred(const gsl_dd a, const double c);
static gsl_dd ldexp_extred(const gsl_dd a, const int e);
static void exp3_add_double(dd_exp3 * x, const double q);
static void exp3_renorm(dd_exp3 * x);
static gsl_dd exp3_to_dd(const dd_exp3 x);

/*
 * Multiply a double-double number by a power of two.
 *
 * The multiplier c is assumed to be an exactly representable power of
 * two.  Scaling the two components separately is therefore exact,
 * provided neither component overflows or underflows.
 */
static gsl_dd
mul_pwr2_extred(const gsl_dd a, const double c)
{
  gsl_dd r;
  r.hi = a.hi * c;
  r.lo = a.lo * c;
  return r;
}

/*
 * Scale a double-double number by an integral power of two.
 *
 * Apply ldexp() independently to the high and low components.  This
 * preserves the represented value exactly provided the component
 * scalings do not overflow or underflow.
 */
static gsl_dd
ldexp_extred(const gsl_dd a, const int e)
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
 *
 * This expansion is used only while forming the high-accuracy exp()
 * range reduction; it is not a general triple-double arithmetic type.
 */
static void
exp3_add_double(dd_exp3 * x, const double q)
{
  gsl_dd s, t;

  s = two_sum(x->a, q);
  t = two_sum(x->b, s.lo);

  x->a = s.hi;
  x->b = t.hi;
  x->c += t.lo;
}

/*
 * Renormalize a three-term range-reduction expansion.
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
exp3_renorm(dd_exp3 * x)
{
  gsl_dd t, u;

  t = two_sum(x->b, x->c);
  u = two_sum(x->a, t.hi);

  x->a = u.hi;
  x->b = u.lo + t.lo;
  x->c = 0.0;
}

/*
 * Compress a three-term range-reduction expansion to double-double.
 *
 * Combine the two trailing components with TwoSum, add their leading
 * contribution to the first component, and renormalize the resulting
 * two-word value with FastTwoSum.
 */
static gsl_dd
exp3_to_dd(const dd_exp3 x)
{
  gsl_dd t, u;

  t = two_sum(x.b, x.c);
  u = two_sum(x.a, t.hi);
  u.lo += t.lo;

  return fast_two_sum(u.hi, u.lo);
}
