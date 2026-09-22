/* specfunc/gamma/test_doublefact.c
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
test_doublefact(void)
{
  int s = 0;

  /*
   * Exact special cases and small exact values.
   */
  TEST_SF(s, gsl_new_sf_doublefact_e, (0U, &r),
          1.00000000000000000000e+00, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_doublefact_e, (1U, &r),
          1.00000000000000000000e+00, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_doublefact_e, (20U, &r),
          3.71589120000000000000e+09, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_doublefact_e, (21U, &r),
          1.37493105750000000000e+10, TEST_TOL0, GSL_SUCCESS);

  /*
   * Finite even branch.
   */
  TEST_SF(s, gsl_new_sf_doublefact_e, (50U, &r),
          5.20469842636666625086e+32, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_doublefact_e, (150U, &r),
          9.37262840292216559752e+131, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_doublefact_e, (298U, &r),
          2.71813802312686482634e+305, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_doublefact_e, (300U, &r),
          8.15441406938059447903e+307, TEST_TOL2, GSL_SUCCESS);

  /*
   * Finite odd branch.
   */
  TEST_SF(s, gsl_new_sf_doublefact_e, (73U, &r),
          1.74865344543353982608e+53, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_doublefact_e, (183U, &r),
          3.61797606666319711405e+168, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_doublefact_e, (299U, &r),
          3.75327411157192600059e+306, TEST_TOL2, GSL_SUCCESS);

  /*
   * Overflow.
   */
  TEST_SF(s, gsl_new_sf_doublefact_e, (301U, &r),
          GSL_POSINF, GSL_POSINF, GSL_EOVRFLW);

  return s;
}
