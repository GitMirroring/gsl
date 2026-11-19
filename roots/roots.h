/* roots/roots.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gough
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

/* roots.h -- declarations for internal root finding and RF support stuff. */

#ifndef __ROOTS_H__
#define __ROOTS_H__

#include <math.h>
#include <gsl/gsl_math.h>

/* Call the pointed-to function with argument x, put its result in y, and 
   return an error if the function value is Inf/Nan. */

#define SAFE_FUNC_CALL(f, x, yp) \
do { \
  *yp = GSL_FN_EVAL(f,x); \
  if (!gsl_finite(*yp)) \
    GSL_ERROR("function value is not finite", GSL_EBADFUNC); \
} while (0)

/* (a + b)/2 and (b - a)/2 without overflow. For finite results these
   are identical to the usual expressions, since scaling by 1/2 is
   exact; the alternative forms are used only when a + b or b - a
   overflows, e.g. for the bracket [-GSL_DBL_MAX, GSL_DBL_MAX]. */

static inline double
midpoint (const double a, const double b)
{
  const double m = 0.5 * (a + b);
  return gsl_finite (m) ? m : 0.5 * a + 0.5 * b;
}

static inline double
half_diff (const double a, const double b)
{
  const double h = 0.5 * (b - a);
  return gsl_finite (h) ? h : 0.5 * b - 0.5 * a;
}

/* Return the trial point c moved, if necessary, to lie at least
   delta = 2 GSL_DBL_EPSILON max(|a|,|b|) inside the bracket (a,b), as
   in Brent's method, or the midpoint if the bracket is narrower than
   2 delta or c is not a number. This keeps every new point distinct
   from the endpoints, so that the bracket collapses once one endpoint
   has converged to the root. */

static inline double
interior_point (const double c, const double a, const double b)
{
  const double delta = 2.0 * GSL_DBL_EPSILON * GSL_MAX (fabs (a), fabs (b));

  if (half_diff (a, b) <= delta || gsl_isnan (c))
    {
      return midpoint (a, b);
    }
  else if (c < a + delta)
    {
      return a + delta;
    }
  else if (c > b - delta)
    {
      return b - delta;
    }

  return c;
}

#endif /* __ROOTS_H__ */


