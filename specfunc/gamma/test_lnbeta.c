/* specfunc/gamma/test_lnbeta.c
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
test_lnbeta(void)
{
  int s = 0;

  /* B1: 0 < a < 1, b < 10. */
  TEST_SF(s, gsl_new_sf_lnbeta_e, (0.25, 0.75, &r),
          1.4913034761293729, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (0.75, 1.5, &r),
          -0.042373001096346445, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (0.99999999999999989, 1.0000000000000002, &r),
          -1.110223024625156e-16, TEST_TOL2, GSL_SUCCESS);

  /* B2: 0 < a < 1, b >= 10. */
  TEST_SF(s, gsl_new_sf_lnbeta_e, (0.25, 10.0, &r),
          0.72182633330350077, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (0.75, 100.0, &r),
          -3.2496599722271737, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (0.007085510039214777, 4.3737812239094046e+295, &r),
          0.12227713696720754, TEST_TOL2, GSL_SUCCESS);

  /* B3: 1 <= a <= 2, b <= 2. */
  TEST_SF(s, gsl_new_sf_lnbeta_e, (1.0, 1.0, &r),
          0.0, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (1.25, 1.75, &r),
          -0.87582013800224401, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (2.0, 2.0, &r),
          -1.791759469228055, TEST_TOL2, GSL_SUCCESS);

  /* B4: 1 <= a <= 2, 2 < b < 10. */
  TEST_SF(s, gsl_new_sf_lnbeta_e, (1.0, 3.0, &r),
          -1.0986122886681098, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (1.5, 5.0, &r),
          -2.605290467144441, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (2.0, 9.5, &r),
          -4.6026670557699729, TEST_TOL2, GSL_SUCCESS);

  /* B5: 1 <= a <= 2, b >= 10. */
  TEST_SF(s, gsl_new_sf_lnbeta_e, (1.0, 10.0, &r),
          -2.3025850929940459, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (1.5, 100.0, &r),
          -7.0322750633365985, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (2.0, 10000000000.0, &r),
          -46.051701859980916, TEST_TOL2, GSL_SUCCESS);

  /* B6: 2 < a < 10, b < 10. */
  TEST_SF(s, gsl_new_sf_lnbeta_e, (2.5, 3.5, &r),
          -3.3018352699620528, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (5.0, 7.0, &r),
          -7.7450028035158391, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (9.0, 9.5, &r),
          -12.6493794533343, TEST_TOL2, GSL_SUCCESS);

  /* B7: 2 < a < 10, 10 <= b <= 1000. */
  TEST_SF(s, gsl_new_sf_lnbeta_e, (2.5, 10.0, &r),
          -5.647837161382057, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (5.0, 100.0, &r),
          -19.946329573136683, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (9.0, 1000.0, &r),
          -51.601093037913266, TEST_TOL2, GSL_SUCCESS);

  /* B8: 2 < a < 10, b > 1000. */
  TEST_SF(s, gsl_new_sf_lnbeta_e, (2.5, 1001.0, &r),
          -16.989075958351606, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (5.0, 10000000000.0, &r),
          -111.95120082035434, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (9.0, 1e+100, &r),
          -2061.7219807918959, TEST_TOL2, GSL_SUCCESS);

  /* B9: a >= 10. */
  TEST_SF(s, gsl_new_sf_lnbeta_e, (10.0, 10.0, &r),
          -13.736229227036555, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (20.0, 100.0, &r),
          -54.550806681721241, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lnbeta_e, (10000000000.0, 1e+20, &r),
          -240258509310.49857, TEST_TOL2, GSL_SUCCESS);

  return s;
}
