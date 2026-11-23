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
 * |x| < sqrt(nu).  For large nu the Cornish-Fisher expansion of
 *
 * G.W. Hill and A.W. Davis. "Generalized asymptotic expansions
 * of Cornish-Fisher type." Annals of Mathematical Statistics,
 * vol. 39, 1264-1273. 1968.
 *
 * is used instead (see TDIST_NU_CF below).
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include "beta_inc.c"

static double
poly_eval (const double c[], unsigned int n, double x)
{
  unsigned int i;
  double y = c[0] * x;

  for (i = 1; i < n; i++)
    {
      y = x * (y + c[i]);
    }

  y += c[n];

  return y;
}

/* 
 * Use the Cornish-Fisher asymptotic expansion to find a point u such
 * that gsl_cdf_gauss(y) = tcdf(t).
 * 
 */

static double
cornish_fisher (double t, double n)
{
  const double coeffs6[10] = {
    0.265974025974025974026,
    5.449696969696969696970,
    122.20294372294372294372,
    2354.7298701298701298701,
    37625.00902597402597403,
    486996.1392857142857143,
    4960870.65,
    37978595.55,
    201505390.875,
    622437908.625
  };
  const double coeffs5[8] = {
    0.2742857142857142857142,
    4.499047619047619047619,
    78.45142857142857142857,
    1118.710714285714285714,
    12387.6,
    101024.55,
    559494.0,
    1764959.625
  };
  const double coeffs4[6] = {
    0.3047619047619047619048,
    3.752380952380952380952,
    46.67142857142857142857,
    427.5,
    2587.5,
    8518.5
  };
  const double coeffs3[4] = {
    0.4,
    3.3,
    24.0,
    85.5
  };

  double a = n - 0.5;
  double b = 48.0 * a * a;

  double z2 = a * gsl_sf_log1p (t * t / n);
  double z = sqrt (z2);

  double p5 = z * poly_eval (coeffs6, 9, z2);
  double p4 = -z * poly_eval (coeffs5, 7, z2);
  double p3 = z * poly_eval (coeffs4, 5, z2);
  double p2 = -z * poly_eval (coeffs3, 3, z2);
  double p1 = z * (z2 + 3.0);
  double p0 = z;

  double y = p5;
  y = (y / b) + p4;
  y = (y / b) + p3;
  y = (y / b) + p2;
  y = (y / b) + p1;
  y = (y / b) + p0;

  if (t < 0)
    y *= -1;

  return y;
}

/* For nu >= TDIST_NU_CF the Cornish-Fisher expansion is used for
   |x| < sqrt(10 nu).  Its terms are powers of z^2/nu with |z| < 38 for
   P >= GSL_DBL_MIN, and compared with high precision values it is
   accurate to about 1e3 ulps there for nu >= 5000, better than the
   incomplete beta function in beta_inc.c, which loses accuracy for large
   parameters (up to 2e4 ulps near nu = 1e5).  Below 5000 the expansion
   loses accuracy rapidly in the tails (1.8e6 ulps at nu = 1000) and the
   incomplete beta function is used.  tdistinv.c uses the same value. */
#define TDIST_NU_CF 5e3

/* below this ln(sqrt(nu)/|x|), e = nu/(nu + x^2) may underflow and the
   tail is given by the leading term of the series for I_e(nu/2, 1/2),
   whose relative error is O(e) */
#define TDIST_LOG_R_SMALL (-333.0)

double
gsl_cdf_tdist_P (const double x, const double nu)
{
  double P;

  if (nu >= TDIST_NU_CF && fabs (x) < sqrt (10.0 * nu))
    {
      return gsl_cdf_ugaussian_P (cornish_fisher (x, nu));
    }

  if (fabs (x) < sqrt (nu))
    {
      double u = x * x / nu;
      double f = u / (1 + u);

      if (x >= 0)
        {
          P = beta_inc_AXPY (0.5, 0.5, 0.5, nu / 2.0, f);
        }
      else
        {
          P = beta_inc_AXPY (-0.5, 0.5, 0.5, nu / 2.0, f);
        }
    }
  else
    {
      /* r = sqrt(nu)/|x| <= 1, e = r^2 / (1 + r^2), computed without
         forming x^2, which can overflow */
      const double logr = 0.5 * log (nu) - log (fabs (x));

      if (logr < TDIST_LOG_R_SMALL)
        {
          /* I_e(a,b) ~ e^a / (a B(a,b)) with a = nu/2, b = 1/2, e = r^2 */
          const double a = nu / 2.0;
          const double tail = 0.5 * exp (2.0 * a * logr - log (a) - gsl_sf_lnbeta (a, 0.5));

          P = (x < 0) ? tail : 1.0 - tail;
        }
      else
        {
          const double r = sqrt (nu) / fabs (x);
          const double v = r * r;
          const double e = v / (1 + v);

          if (x >= 0)
            {
              P = beta_inc_AXPY (-0.5, 1.0, nu / 2.0, 0.5, e);
            }
          else
            {
              P = beta_inc_AXPY (0.5, 0.0, nu / 2.0, 0.5, e);
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
