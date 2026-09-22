/* specfunc/gamma/test_fact.c
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
test_fact(void)
{
  int s = 0;

  /*
   * Exact factorials.
   */
  TEST_SF(s, gsl_new_sf_fact_e, (0U, &r),
          1.00000000000000000000e+00, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_fact_e, (1U, &r),
          1.00000000000000000000e+00, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_fact_e, (5U, &r),
          1.20000000000000000000e+02, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_fact_e, (22U, &r),
          1.12400072777760768000e+21, TEST_TOL0, GSL_SUCCESS);

  /*
   * Rounded finite factorials.
   */
  TEST_SF(s, gsl_new_sf_fact_e, (23U, &r),
          2.58520167388849782129e+22, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_fact_e, (50U, &r),
          3.04140932017133755764e+64, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_fact_e, (100U, &r),
          9.33262154439441509656e+157, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_fact_e, (170U, &r),
          7.25741561530799903500e+306, TEST_TOL1, GSL_SUCCESS);

  /*
   * Overflow.
   */
  TEST_SF(s, gsl_new_sf_fact_e, (171U, &r),
          GSL_POSINF, GSL_POSINF, GSL_EOVRFLW);

  return s;
}
