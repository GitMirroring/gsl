/* specfunc/gamma/test_taylorcoeff.c
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
test_taylorcoeff(void)
{
  int s = 0;

  /*
   * Domain errors.
   */
  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (-1, 1.0, &r),
          GSL_NAN, TEST_TOL0, GSL_EDOM);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (2, -1.0, &r),
          GSL_NAN, TEST_TOL0, GSL_EDOM);

  /*
   * Exact special cases: n = 0, n = 1, and x = 0.
   */
  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (0, 123.0, &r),
          1.00000000000000000000e+00, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (1, 0.125, &r),
          1.25000000000000000000e-01, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (7, 0.0, &r),
          0.00000000000000000000e+00, TEST_TOL0, GSL_SUCCESS);

  /*
   * Scaled-power kernel: 2 <= n <= 100.
   *
   * These points exercise low, moderate, and upper-end orders in the
   * precomputed reciprocal-factorial branch.
   */
  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (2, 0.1, &r),
          5.00000000000000097145e-03, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (5, 1.25, &r),
          2.54313151041666678231e-02, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (10, 20.0, &r),
          2.82186948853615531698e+06, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (20, 5.5, &r),
          2.63711563134619796610e-04, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (30, 30.0, &r),
          7.76207020879728149414e+11, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (50, 25.0, &r),
          2.59373475312613067217e+05, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (100, 55.5, &r),
          2.87935131741831720000e+16, TEST_TOL2, GSL_SUCCESS);

  /*
   * Generic large-order branch: n > 100.
   */
  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (101, 25.0, &r),
          1.65050121123519168813e-19, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (150, 50.0, &r),
          1.22632968045484534447e-08, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (200, 100.0, &r),
          1.26797695348096246921e+25, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (500, 200.0, &r),
          2.68280617236282960000e+16, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (1000, 300.0, &r),
          3.28556828371631311966e-91, TEST_TOL4, GSL_SUCCESS);

  /*
   * Underflow and overflow in the large-order branch.
   */
  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (200, 1e-10, &r),
          0.0, TEST_TOL0, GSL_EUNDRFLW);

  TEST_SF(s, gsl_new_sf_taylorcoeff_e, (200, 10000000000.0, &r),
          GSL_POSINF, GSL_POSINF, GSL_EOVRFLW);

  return s;
}
