/* specfunc/gamma/test_beta_inc.c
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

#include "../test_sf.h"

/* reference values computed at 50 digits with mpmath (series and continued
   fraction), checked against mpmath.betainc, or for large parameters
   against binomial sums and quadrature */

int
test_beta_inc(void)
{
  int s = 0;

  /* Domain errors. */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (1.0, 1.0, -0.1, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (1.0, 1.0, 1.1, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (-2.0, 1.0, 0.5, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  /* Endpoints. */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (2.0, 3.0, 0.0, &r),
          1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (2.0, 3.0, 1.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  /* I_{1/2}(a,a) = 1/2 exactly. */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (0.3, 0.3, 0.5, &r),
          0.5, TEST_TOL0, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (3000.0, 3000.0, 0.5, &r),
          0.5, TEST_TOL0, GSL_SUCCESS);

  /* Q = 1 - P for nonpositive parameters. */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (-0.1, -0.2, 0.5, &r),
          0.324747998041610028008665,
          TEST_TOL2, GSL_SUCCESS);

  /* Subnormal result: no underflow error, as for gsl_sf_gamma_inc_P. */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (128.16604656167866, 120.87446678101553, 0.001, &r),
          5.7822410806596158761e-312,
          TEST_TOL6, GSL_SUCCESS);

  /* Result below the smallest subnormal. */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (500.0, 500.0, 0.01, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  /* power series */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (0.5, 0.5, 0.1, &r),
          0.20483276469913345754,
          TEST_TOL1, GSL_SUCCESS);

  /* power series */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (2, 3, 0.2, &r),
          0.18080000000000001705,
          TEST_TOL1, GSL_SUCCESS);

  /* power series, small a */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (0.001, 2, 0.3, &r),
          0.99949590940632457545,
          TEST_TOL1, GSL_SUCCESS);

  /* power series, small b */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (5, 0.5, 0.01, &r),
          2.4712578086395446223e-11,
          TEST_TOL1, GSL_SUCCESS);

  /* complement of the power series */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (0.5, 0.5, 0.1, &r),
          0.79516723530086654246,
          TEST_TOL1, GSL_SUCCESS);

  /* continued fraction */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (10, 20, 0.3, &r),
          0.36400408107194422775,
          TEST_TOL1, GSL_SUCCESS);

  /* continued fraction */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (50, 60, 0.5, &r),
          0.83090729390166941434,
          TEST_TOL1, GSL_SUCCESS);

  /* continued fraction near the transition */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (500, 500, 0.49, &r),
          0.26357560835658528519,
          TEST_TOL2, GSL_SUCCESS);

  /* continued fraction, small result */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (5000, 5000, 0.4, &r),
          4.5185437272610854309e-91,
          TEST_TOL1, GSL_SUCCESS);

  /* continued fraction, small result */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (5000, 2000, 0.6, &r),
          8.4453887739027059415e-89,
          TEST_TOL1, GSL_SUCCESS);

  /* continued fraction, small complement */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (2000, 5000, 0.4, &r),
          8.4453887739027059415e-89,
          TEST_TOL1, GSL_SUCCESS);

  /* continued fraction, large parameters */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (15966.3, 52454.3, 0.232506, &r),
          0.30014718924392706511,
          TEST_TOL1, GSL_SUCCESS);

  /* continued fraction, large parameters */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (6784.84, 82879.7, 0.0756661, &r),
          0.50004502047848345574,
          TEST_TOL1, GSL_SUCCESS);

  /* continued fraction, large parameters */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (73807.9, 69740.1, 0.513477, &r),
          0.29998439465801664386,
          TEST_TOL1, GSL_SUCCESS);

  /* large a, small b expansion */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (100, 0.5, 0.99, &r),
          0.15677586542444078337,
          TEST_TOL1, GSL_SUCCESS);

  /* large a, small b expansion */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (10000, 0.1, 0.9999, &r),
          0.024127150309053508845,
          TEST_TOL1, GSL_SUCCESS);

  /* large a, small b expansion */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (513.495, 0.0488652, 0.997911, &r),
          0.0099961925494060508472,
          TEST_TOL2, GSL_SUCCESS);

  /* large a, small b expansion */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (603926, 0.122301, 0.999997, &r),
          0.0092277118072517134303,
          TEST_TOL2, GSL_SUCCESS);

  /* small b expansion with a shifted */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (2, 0.5, 0.95, &r),
          0.67017997331878087828,
          TEST_TOL1, GSL_SUCCESS);

  /* shifted expansion with near-zero odd coefficients */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (1.65291e-06, 0.00293736, 8.54994e-116, &r),
          0.00099999804439519395252,
          TEST_TOL1, GSL_SUCCESS);

  /* small b: the smaller tail is not the one a/(a+b) suggests */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (32, 0.002, 0.99999, &r),
          0.014863631718544292863,
          TEST_TOL1, GSL_SUCCESS);

  /* small a and b */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (0.00402207, 0.0460761, 0.00423218, &r),
          0.10000000200607816645,
          TEST_TOL2, GSL_SUCCESS);

  /* small a */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (0.0076923, 0.836189, 1.52817e-06, &r),
          0.099999985587033940037,
          TEST_TOL1, GSL_SUCCESS);

  /* subnormal x */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (0.000961682, 46946.7, 1.13172e-318, &r),
          0.50000016107364777214,
          TEST_TOL1, GSL_SUCCESS);

  /* subnormal x */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (0.000961682, 46946.7, 1.13172e-318, &r),
          0.49999983892635222786,
          TEST_TOL1, GSL_SUCCESS);

  /* tiny x */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (0.000599239, 27322.1, 8.98781e-82, &r),
          0.10000004846964818677,
          TEST_TOL2, GSL_SUCCESS);

  /* small b, x near 1 */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (0.451732, 0.0130053, 0.999, &r),
          0.10496090777253726024,
          TEST_TOL1, GSL_SUCCESS);

  /* tiny parameters */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (1e-300, 1e-300, 0.5, &r),
          0.50000000000000000000,
          TEST_TOL1, GSL_SUCCESS);

  /* tiny parameters and x */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (1e-06, 1e-06, 1e-300, &r),
          0.49965473150211685875,
          TEST_TOL1, GSL_SUCCESS);

  /* symmetric */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (90, 90, 0.5, &r),
          0.50000000000000000000,
          TEST_TOL1, GSL_SUCCESS);

  /* symmetric, large */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (1000000, 1000000, 0.5, &r),
          0.50000000000000000000,
          TEST_TOL1, GSL_SUCCESS);

  /* large parameters, deep tail */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_e, (1000000, 1000000, 0.4999, &r),
          0.38864871786232204188,
          TEST_TOL2, GSL_SUCCESS);

  /* Q = 2^-50 exactly */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (1, 50, 0.5, &r),
          8.8817841970012523234e-16,
          TEST_TOL1, GSL_SUCCESS);

  /* Q = 1 - 2^-10 */
  TEST_SF_RESULT(s, gsl_sf_beta_inc_Q_e, (10, 1, 0.5, &r),
          0.99902343750000000000,
          TEST_TOL1, GSL_SUCCESS);

  return s;
}
