/* specfunc/gamma/test_lndoublefact.c
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
test_lndoublefact(void)
{
  int s = 0;

  /*
   * Exact special cases.
   */
  TEST_SF(s, gsl_new_sf_lndoublefact_e, (0U, &r),
          0.00000000000000000000e+00, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lndoublefact_e, (1U, &r),
          0.00000000000000000000e+00, TEST_TOL0, GSL_SUCCESS);

  /*
   * Shared-table even branch: n <= 300.
   */
  TEST_SF(s, gsl_new_sf_lndoublefact_e, (2U, &r),
          6.93147180559945286227e-01, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lndoublefact_e, (100U, &r),
          1.83135125979770293725e+02, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lndoublefact_e, (300U, &r),
          7.08992182933415506341e+02, TEST_TOL2, GSL_SUCCESS);

  /*
   * Shared-table odd branch: n <= 300.
   */
  TEST_SF(s, gsl_new_sf_lndoublefact_e, (3U, &r),
          1.09861228866810978211e+00, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lndoublefact_e, (299U, &r),
          7.05913667011652478323e+02, TEST_TOL2, GSL_SUCCESS);

  /*
   * Large-n odd lngamma branch.
   */
  TEST_SF(s, gsl_new_sf_lndoublefact_e, (301U, &r),
          7.11620777276401440758e+02, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lndoublefact_e, (401U, &r),
          1.00463323610615339021e+03, TEST_TOL2, GSL_SUCCESS);

  /*
   * Large-n even lngamma branch.
   */
  TEST_SF(s, gsl_new_sf_lndoublefact_e, (302U, &r),
          7.14702609950790360926e+02, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lndoublefact_e, (1000U, &r),
          2.95790404874012892833e+03, TEST_TOL2, GSL_SUCCESS);

  return s;
}
