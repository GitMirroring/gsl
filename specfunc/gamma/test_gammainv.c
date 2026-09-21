/* specfunc/gamma/test_gammainv.c
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

#include "../test_sf_new.h"

int
test_gammainv(void)
{
  int s = 0;

  /*
   * Negative reflection region: x < -10.
   *
   * Exact negative integers are zeros of reciprocal Gamma.  The first
   * nonintegral point also exercises the negative-overflow branch.
   */
  TEST_SF(s, gsl_new_sf_gammainv_e, (-200.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-176.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-175.75, &r),
          GSL_POSINF, GSL_POSINF, GSL_EOVRFLW);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-100.25, &r),
          -6.6529717048286683418e+157, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-20.5, &r),
          -3.5277642061974392980e+18, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-10.5, &r),
          -3.7877040075087946942e+06, TEST_TOL4, GSL_SUCCESS);

  /*
   * Moderate negative recurrence region: -10 <= x < -1/2.
   */
  TEST_SF(s, gsl_new_sf_gammainv_e, (-10.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-9.25, &r),
          1.4387788462519049915e+05, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-5.5, &r),
          9.1636730015295728169e+01, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-3.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-2.25, &r),
          -5.7378441030346615857e-01, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-1.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-0.75, &r),
          -2.0686174712265698577e-01, TEST_TOL4, GSL_SUCCESS);

  /*
   * Central H representation: -1/2 <= x <= 1.
   */
  TEST_SF(s, gsl_new_sf_gammainv_e, (-0.5, &r),
          -2.8209479177387814347e-01, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (-0.25, &r),
          -2.0401223477456574527e-01, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (0.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (0.25, &r),
          2.7581566283020931436e-01, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (0.5, &r),
          5.6418958354775628695e-01, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (1.0, &r),
          1.0, TEST_TOL0, GSL_SUCCESS);

  /*
   * Shifted H representation: 1 < x <= 2.
   */
  TEST_SF(s, gsl_new_sf_gammainv_e, (1.25, &r),
          1.1032626513208372574e+00, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (1.5, &r),
          1.1283791670955125739e+00, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (2.0, &r),
          1.0, TEST_TOL0, GSL_SUCCESS);

  /*
   * Positive low-Gamma kernel: 2 < x <= 10.
   */
  TEST_SF(s, gsl_new_sf_gammainv_e, (2.25, &r),
          8.8261012105666980595e-01, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (5.5, &r),
          1.9104832458760001251e-02, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (10.0, &r),
          2.7557319223985890653e-06, TEST_TOL4, GSL_SUCCESS);

  /*
   * Positive high-Gamma kernel.  10 < x < 12 exercises the shifted
   * Stirling branch; x >= 12 uses the direct scaled reconstruction.
   */
  TEST_SF(s, gsl_new_sf_gammainv_e, (10.5, &r),
          8.8239572002038009055e-07, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (11.5, &r),
          8.4037687620988580052e-08, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (12.0, &r),
          2.5052108385441718775e-08, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (20.0, &r),
          8.2206352466243297170e-18, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (100.0, &r),
          1.0715102881254669232e-156, TEST_TOL4, GSL_SUCCESS);

  /*
   * Normal/subnormal reciprocal-Gamma transition and subnormal range.
   */
  TEST_SF(s, gsl_new_sf_gammainv_e, (171.0, &r),
          1.3779009677917705867e-307, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (171.5, &r),
          1.0544777400574992603e-308, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (175.0, &r),
          1.5563171257343448644e-316, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammainv_e, (178.0, &r),
          2.8547896502574379345e-323, TEST_TOL4, GSL_SUCCESS);

  /*
   * Positive underflow branch: x > GAMMA_SCALED_X_MAX.
   * The exact reciprocal Gamma is below the binary64 range.
   */
  TEST_SF(s, gsl_new_sf_gammainv_e, (201.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  return s;
}
