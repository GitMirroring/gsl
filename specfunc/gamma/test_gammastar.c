/* specfunc/gamma/test_gammastar.c
 *
 * Copyright (C) 2026 Patrick Alken
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
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
test_gammastar(void)
{
  int s = 0;

  /*
   * Small-x regularized kernel: 0 < x < 0.4.
   *
   * The smallest positive binary64 value is included to exercise the
   * extreme subnormal/singular scaling of Gamma*(x).
   */

  TEST_SF(s, gsl_new_sf_gammastar_e, (5e-324, &r),
          1.7948069285245254e+161, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1e-300, &r),
          3.9894228040143267e+149, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1e-200, &r),
          3.9894228040143267e+99, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1e-100, &r),
          3.989422804014327e+49, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1e-20, &r),
          3989422804.014327, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1e-10, &r),
          39894.22813368979, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (0.0001, &r),
          39.93267749222754, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (0.01, &r),
          4.195480832378425, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (0.1, &r),
          1.669860485516386, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (0.25, &r),
          1.3132571183090294, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (0.39, &r),
          1.2094254740208996, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (0.39999999999999997, &r),
          1.2045550086474184, TEST_TOL4, GSL_SUCCESS);

  /*
   * Finite Chebyshev kernel: 0.4 <= x < 2.
   *
   * Values immediately below the two dispatcher boundaries exercise the
   * adjacent branch without relying on a non-binary64 decimal boundary.
   */

  TEST_SF(s, gsl_new_sf_gammastar_e, (0.4, &r),
          1.2045550086474184, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (0.5, &r),
          1.165821990798562, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (0.75, &r),
          1.1121142533175508, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1.0, &r),
          1.0844375514192275, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1.25, &r),
          1.06762010132462, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1.5, &r),
          1.0563442442685598, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1.9999999999999998, &r),
          1.042207120816673, TEST_TOL4, GSL_SUCCESS);

  /*
   * Finite Chebyshev kernel: 2 <= x < 16.
   */

  TEST_SF(s, gsl_new_sf_gammastar_e, (2.0, &r),
          1.042207120816673, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (2.25, &r),
          1.0374910051074155, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (3.0, &r),
          1.0280645179187893, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (5.0, &r),
          1.0167839858278083, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (8.0, &r),
          1.010465651061946, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (10.0, &r),
          1.0083653591324002, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (12.0, &r),
          1.0069669974961886, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (14.0, &r),
          1.0059691146999385, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (15.0, &r),
          1.005570189689369, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (15.999999999999998, &r),
          1.005221239328702, TEST_TOL4, GSL_SUCCESS);

  /*
   * Direct asymptotic kernel: x >= 16.
   *
   * Include the crossover, moderate values, and very large arguments.
   */

  TEST_SF(s, gsl_new_sf_gammastar_e, (16.0, &r),
          1.005221239328702, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (16.000000000000004, &r),
          1.005221239328702, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (20.0, &r),
          1.0041750108677654, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (32.0, &r),
          1.0026074754856418, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (50.0, &r),
          1.0016680340707353, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (100.0, &r),
          1.0008336778720122, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1000.0, &r),
          1.000083336802874, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1000000.0, &r),
          1.000000083333337, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1e+20, &r),
          1.0, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1e+100, &r),
          1.0, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1e+300, &r),
          1.0, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gammastar_e, (1.7976931348623155e+308, &r),
          1.0, TEST_TOL4, GSL_SUCCESS);

  /*
   * Domain x <= 0.
   */
  TEST_SF(s, gsl_new_sf_gammastar_e, (0.0, &r),
          GSL_NAN, TEST_TOL0, GSL_EDOM);

  TEST_SF(s, gsl_new_sf_gammastar_e, (-0.25, &r),
          GSL_NAN, TEST_TOL0, GSL_EDOM);

  TEST_SF(s, gsl_new_sf_gammastar_e, (-1.0, &r),
          GSL_NAN, TEST_TOL0, GSL_EDOM);

  return s;
}
