/* specfunc/gamma/test_lnfact.c
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
test_lnfact(void)
{
  int s = 0;

  /*
   * Shared factorial-table branch: n <= 300.
   */
  TEST_SF(s, gsl_new_sf_lnfact_e, (0U, &r),
          0.00000000000000000000e+00, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnfact_e, (1U, &r),
          0.00000000000000000000e+00, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnfact_e, (2U, &r),
          6.93147180559945286227e-01, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnfact_e, (22U, &r),
          4.84711813518352272467e+01, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnfact_e, (100U, &r),
          3.63739375555563469788e+02, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnfact_e, (300U, &r),
          1.41490584994506798466e+03, TEST_TOL2, GSL_SUCCESS);

  /*
   * Large-n lngamma fallback: n > 300.
   */
  TEST_SF(s, gsl_new_sf_lnfact_e, (301U, &r),
          1.42061296020981694710e+03, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnfact_e, (333U, &r),
          1.60493470919185779167e+03, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnfact_e, (1000U, &r),
          5.91212817848816303012e+03, TEST_TOL2, GSL_SUCCESS);

  return s;
}
