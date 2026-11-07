/* specfunc/gamma/test_beta.c
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
test_beta(void)
{
  int s = 0;

  /*
   * 0 < a <= b < 10: both arguments reduced to (0,1] by the exact
   * recurrence.  Reference values from 300-bit arithmetic.
   */
  TEST_SF(s, gsl_sf_beta_e, (1.0, 1.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (0.25, 0.75, &r), 4.44288293815836624702, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (2.5, 3.5, &r), 3.68155389092553895132e-02, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (4.0, 5.0, &r), 3.57142857142857142857e-03, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (9.75, 9.5, &r), 1.8597635976851754746e-06, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (1e-200, 3e-150, &r), 1.0000000000000000179e+200, TEST_TOL1, GSL_SUCCESS);

  /* a < 10 <= b: Stirling ratio Gamma(b)/Gamma(a+b) */
  TEST_SF(s, gsl_sf_beta_e, (0.5, 10.0, &r), 5.67546385503041849791e-01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (2.5, 300.0, &r), 8.47470784778353801499e-07, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (9.5, 1e6, &r), 1.19287645673135061725e-52, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (1e-5, 1e300, &r), 9.93110316097755390005e+04, TEST_TOL1, GSL_SUCCESS);

  /*
   * 10 <= a <= b: Stirling form with exactly formed powers; formerly
   * exp(log B), with relative error ~ |log B| eps (here ~700 eps).
   */
  TEST_SF(s, gsl_sf_beta_e, (10.0, 10.0, &r), 1.08250882244690294226e-06, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (200.0, 2000.0, &r), 1.60526424228496071453e-292, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (8035.125, 90.5, &r), 3.75229240073801339037e-217, TEST_TOL1, GSL_SUCCESS);

  /* negative arguments, by reflection to positive arguments */
  TEST_SF(s, gsl_sf_beta_e, (-0.5, 0.75, &r), -1.19814023473559220744, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (-2.5, -1.25, &r), -1.38385197111960899959e+01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (-100.1, 99.3, &r), 4.4742581995800501513e-04, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (-100.1, -99.1, &r), -1.00518191779769781601e+60, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (-1e-5, 2.0, &r), -1.00001000010000091821e+05, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (-10.999999, -6.25, &r), 1.58027236304663200885e+10, TEST_TOL1, GSL_SUCCESS);

  /* zero (a+b a non-positive integer), poles, overflow and underflow */
  TEST_SF(s, gsl_sf_beta_e, (-1.5, 0.5, &r), 0.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_beta_e, (0.0, 1.0, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);
  TEST_SF(s, gsl_sf_beta_e, (-3.0, 0.5, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);
  TEST_SF(s, gsl_sf_beta_e, (1e-310, 1e-310, &r), GSL_POSINF, GSL_POSINF, GSL_EOVRFLW);
  TEST_SF(s, gsl_sf_beta_e, (1e6, 1e7, &r), 0.0, TEST_TOL0, GSL_EUNDRFLW);

  return s;
}
