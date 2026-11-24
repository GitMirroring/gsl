/* cdf/tdist.c
 *
 * Copyright (C) 2002 Jason H. Stover.
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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * Computes the Student's t cumulative distribution function from the
 * incomplete beta function (Abramowitz & Stegun 26.7.1):
 *
 *   P(T < -|x|) = 1/2 I_e(nu/2, 1/2),       e = nu / (nu + x^2),
 *   P(|T| < |x|) = I_f(1/2, nu/2),          f = x^2 / (nu + x^2),
 *
 * using the first form for |x| >= sqrt(nu) and the second for
 * |x| < sqrt(nu).
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>


/* below this ln(sqrt(nu)/|x|), e = nu/(nu + x^2) may underflow and the
   tail is given by the leading term of the series for I_e(nu/2, 1/2),
   whose relative error is O(e) */
#define TDIST_LOG_R_SMALL (-333.0)

double
gsl_cdf_tdist_P (const double x, const double nu)
{
  double P;

  if (fabs (x) < sqrt (nu))
    {
      double u = x * x / nu;
      double f = u / (1 + u);

      if (x >= 0)
        {
          P = 0.5 + 0.5 * gsl_sf_beta_inc (0.5, nu / 2.0, f);
        }
      else
        {
          P = 0.5 * gsl_sf_beta_inc_Q (0.5, nu / 2.0, f);
        }
    }
  else
    {
      /* r = sqrt(nu)/|x| <= 1, e = r^2 / (1 + r^2), computed without
         forming x^2, which can overflow */
      const double logr = 0.5 * log (nu) - log (fabs (x));

      if (logr < TDIST_LOG_R_SMALL)
        {
          /* I_e(a,b) ~ e^a / (a B(a,b)) with a = nu/2, b = 1/2, e = r^2;
             r^nu is formed with pow when r and a are normal, since the
             error of exp(nu ln r) grows with |nu ln r| */
          const double a = nu / 2.0;
          const double r = sqrt (nu) / fabs (x);
          const double tail = (r >= GSL_DBL_MIN && a >= GSL_DBL_MIN)
            ? 0.5 * pow (r, nu) / (a * gsl_sf_beta (a, 0.5))
            : 0.5 * exp (nu * logr - log (a) - gsl_sf_lnbeta (a, 0.5));

          P = (x < 0) ? tail : 1.0 - tail;
        }
      else
        {
          const double r = sqrt (nu) / fabs (x);
          const double v = r * r;
          const double e = v / (1 + v);

          if (x >= 0)
            {
              P = 1.0 - 0.5 * gsl_sf_beta_inc (nu / 2.0, 0.5, e);
            }
          else
            {
              P = 0.5 * gsl_sf_beta_inc (nu / 2.0, 0.5, e);
            }
        }
    }

  return P;
}

double
gsl_cdf_tdist_Q (const double x, const double nu)
{
  /* the distribution is symmetric */
  return gsl_cdf_tdist_P (-x, nu);
}
