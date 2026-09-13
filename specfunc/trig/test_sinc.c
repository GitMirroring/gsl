/* specfunc/trig/test_sinc.c
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
#include <gsl/gsl_sf_trig.h>

#include "../test_sf_new.h"

int
test_sinc(void)
{
  int s = 0;

  /*
   * Small-argument polynomial region: |x| <= 0.5.
   * Include the origin, representative interior points, and the
   * transition point itself.
   */
  TEST_SF(s, gsl_sf_sinc_e, (0.0, &r),
          1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (1.0/1024.0, &r),
          0.9999984312693665404344614320988369838505, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (0.1, &r),
          0.9836316430834659655392050846235283603732, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (0.25, &r),
          0.9003163161571060695551991910067405826646, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (0.49, &r),
          0.6492914683060573829375459187899089702951, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (0.5, &r),
          0.6366197723675813430755350534900574481378, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (-0.5, &r),
          0.6366197723675813430755350534900574481378, TEST_TOL0, GSL_SUCCESS);

  /*
   * Reduced-argument region: |x| > 0.5.
   * These points exercise positive and negative reduced arguments,
   * alternating signs, and points near the extrema |r| = 0.5.
   */
  TEST_SF(s, gsl_sf_sinc_e, (0.75, &r),
          0.3001054387190353565183997303355801942215, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (1.25, &r),
          -0.1800632632314212139110398382013481165329, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (1.5, &r),
          -0.2122065907891937810251783511633524827126, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (2.25, &r),
          0.1000351462396784521727999101118600647405, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (5.4712, &r),
          -0.05794120000086965577374709280237705351242, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (5.5288, &r),
          -0.05733755850180111719818141974533404703245, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (-5.5288, &r),
          -0.05733755850180111719818141974533404703245, TEST_TOL0, GSL_SUCCESS);

  /*
   * Exact integer zeros.  The new reduction is designed to preserve
   * these zeros exactly.
   */
  TEST_SF(s, gsl_sf_sinc_e, (1.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (-1.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (10.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (-10.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  /*
   * Large arguments.  Include cases from the original GSL sinc tests
   * and additional points exercising the reduced-argument formulation.
   */
  TEST_SF(s, gsl_sf_sinc_e, (80.5, &r),
          0.003954160076817275422829410270124580423216, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (100.5, &r),
          0.003167262549092444492913109718856007204666, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (1000.25, &r),
          0.0002250228233334431565996498852803650544025, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (1000000.5, &r),
          0.000000318309727028927157074188989650533898802, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (22909.473649078136, &r),
          -0.00001384666440873361480204870222190494207555, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (1000000000000.25, &r),
          2.25079079039220247619039992689780385668e-13, TEST_TOL0, GSL_SUCCESS);

  /*
   * For IEEE binary64 all finite representable values with
   * |x| >= 2^52 are integral, so sinc(x) is exactly zero.
   */
  TEST_SF(s, gsl_sf_sinc_e, (4503599627370496.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_sinc_e, (-4503599627370496.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  return s;
}
