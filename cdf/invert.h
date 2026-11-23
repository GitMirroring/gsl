/* cdf/invert.h
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

/* Safeguarded Halley iteration for inverting a continuous distribution
   function.

   The caller supplies fn(x, params, &g, &u, &d) which returns, for the
   current x,

     g  a residual which is an increasing function of x and is zero at
        the solution (used only for its sign, to maintain a bracket),
     u  the Newton correction f/f' of the equation f(x) = 0 being solved,
     d  the ratio f''/f'.

   Each iteration takes the Halley step x - u / (1 - u d / 2) (or the
   Newton step x - u if the Halley denominator is not close to 1).  The
   iteration keeps a bracket [lo, hi] around the solution, starting from
   the support of the distribution, and replaces any step which leaves
   the bracket by a bisection: geometric when the bracket spans more than
   a factor of 4, otherwise arithmetic; when lo = 0 the upper end is
   divided by 1024, and when hi is infinite the current point is
   multiplied by 4.  The iteration stops when the step is below a few
   ulps of x.  The support [lo, hi] must satisfy 0 <= lo < hi. */

#ifndef __GSL_CDF_INVERT_H__
#define __GSL_CDF_INVERT_H__

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#define CDF_INVERT_MAXITER 300

typedef int (*cdf_invert_fn) (const double x, void *params,
                              double *g, double *u, double *d);

static int
cdf_invert (cdf_invert_fn fn, void *params, const double x0,
            double lo, double hi, double *result)
{
  double x = x0;
  int iter;

  if (!(x > lo && x < hi))
    {
      /* starting value outside the support: start inside the bracket */
      if (lo > 0.0 && gsl_finite (hi))
        x = 0.5 * lo + 0.5 * hi;
      else if (gsl_finite (hi))
        x = hi / 1024.0;
      else
        x = (lo > 0.0) ? 4.0 * lo : 1.0;
    }

  for (iter = 0; iter < CDF_INVERT_MAXITER; iter++)
    {
      double g, u, d, h, xn;
      int status = fn (x, params, &g, &u, &d);

      if (status)
        return status;

      if (g == 0.0)
        {
          *result = x;
          return GSL_SUCCESS;
        }

      /* tighten the bracket */
      if (g < 0.0)
        lo = x;
      else
        hi = x;

      /* Halley step, or Newton step when the correction is large */
      h = u;
      if (fabs (0.5 * u * d) < 0.5)
        h = u / (1.0 - 0.5 * u * d);

      /* converged: the step is below a few ulps of x (x - h may round
         to x itself, which is now an end of the bracket) */
      if (fabs (h) <= 4.0 * GSL_DBL_EPSILON * fabs (x))
        {
          *result = x - h;
          return GSL_SUCCESS;
        }

      xn = x - h;

      if (!(xn > lo && xn < hi) || !gsl_finite (xn))
        {
          /* bisection */
          if (lo > 0.0 && gsl_finite (hi))
            xn = (hi > 4.0 * lo) ? sqrt (lo) * sqrt (hi) : 0.5 * lo + 0.5 * hi;
          else if (gsl_finite (hi))
            xn = (lo == 0.0) ? hi / 1024.0 : 0.5 * lo + 0.5 * hi;
          else
            xn = 4.0 * GSL_MAX_DBL (x, lo);

          if (!(xn > lo && xn < hi))
            {
              /* bracket cannot be split further */
              *result = x;
              return GSL_SUCCESS;
            }
        }

      if (fabs (xn - x) <= 4.0 * GSL_DBL_EPSILON * fabs (xn))
        {
          *result = xn;
          return GSL_SUCCESS;
        }

      x = xn;
    }

  *result = x;
  return GSL_EMAXITER;
}

#endif /* __GSL_CDF_INVERT_H__ */
