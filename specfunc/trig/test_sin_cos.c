/* specfunc/trig/test_sin_cos.c
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
#include <gsl/gsl_nan.h>
#include <gsl/gsl_sf_trig.h>

#include "../test_sf.h"

/* reference values: mpmath with 4000-bit arithmetic (enough to reduce
   any double argument exactly) */

int
test_sin_cos(void)
{
  gsl_sf_result r;
  int s = 0;

  /* moderate arguments (from specfunc/test_sf.c) */
  TEST_SF(s, gsl_sf_sin_e, (-10.0, &r), 0.5440211108893698134047, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (-10.0, &r), -0.8390715290764524522589, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (1.0, &r), 0.8414709848078965066525, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1.0, &r), 0.5403023058681397174009, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (1000.0, &r), 0.8268795405320025602559, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1000.0, &r), 0.5623790762907029910782, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (1048576.75, &r), 0.8851545351115651914008, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1048576.75, &r), 0.4652971620066351799007, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (62831853.75, &r), 0.6273955953485000827111, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (62831853.75, &r), 0.7787006914966116436461, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (1073741822.5, &r), -0.8284043541754465987557, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1073741822.5, &r), -0.5601305436977716102466, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (1073741824.0, &r), -0.6173264150460421707542, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1073741824.0, &r), 0.7867071229411881196183, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (1073741825.5, &r), 0.7410684679436226925804, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1073741825.5, &r), 0.671429464515589951389, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (1099511627776.0, &r), -0.4057050115328287198206, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1099511627776.0, &r), -0.9140040719915570022997, TEST_TOL0, GSL_SUCCESS);

  /* large arguments: the previous argument reduction lost all accuracy
     beyond about 1e8 and gave values far outside [-1,1] beyond about
     1e16 (bugs #45746, #52570) */
  TEST_SF(s, gsl_sf_sin_e, (1.0e16, &r), 0.7796880066069787502355, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1.0e16, &r), -0.6261681981330861717628, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (1.0e17, &r), -0.4645301048353726961545, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1.0e17, &r), -0.8855573282976306850496, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (1.0e22, &r), -0.8522008497671888017727, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1.0e22, &r), 0.5232147853951389454976, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (8.1e50, &r), 0.07771142392161757197073, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (8.1e50, &r), -0.9969758946895730458433, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (1.0e300, &r), -0.8178819121159085970459, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1.0e300, &r), -0.5753861119575490466882, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (GSL_DBL_MAX, &r), 0.004961954789184061790503, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (GSL_DBL_MAX, &r), -0.9999876894265599374649, TEST_TOL0, GSL_SUCCESS);

  /* near zeros of sin and cos, and tiny arguments */
  TEST_SF(s, gsl_sf_sin_e, (M_PI, &r), 1.224646799147353177226e-16, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (M_PI, &r), -1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (M_PI_2, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (M_PI_2, &r), 6.12323399573676588613e-17, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (0.0, &r), 0.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (0.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (1.0e-300, &r), 1.0e-300, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (-1.0e-300, &r), -1.0e-300, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_e, (1.0e-300, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_e, (4.9406564584124654e-324, &r), 4.9406564584124654e-324, TEST_TOL0, GSL_SUCCESS);

  /* arguments with an uncertainty dx: the error includes |cos x| dx
     (resp. |sin x| dx) */
  TEST_SF(s, gsl_sf_sin_err_e, (1.0, 1.0e-10, &r), 0.8414709848078965066525, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_cos_err_e, (1.0, 1.0e-10, &r), 0.5403023058681397174009, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_sin_err_e, (1.0e22, 0.0, &r), -0.8522008497671888017727, TEST_TOL0, GSL_SUCCESS);

  /* non-finite arguments: NaN with GSL_EDOM */
  TEST_SF(s, gsl_sf_sin_e, (GSL_POSINF, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);
  TEST_SF(s, gsl_sf_sin_e, (GSL_NEGINF, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);
  TEST_SF(s, gsl_sf_sin_e, (GSL_NAN, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);
  TEST_SF(s, gsl_sf_cos_e, (GSL_POSINF, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);
  TEST_SF(s, gsl_sf_cos_e, (GSL_NEGINF, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);
  TEST_SF(s, gsl_sf_cos_e, (GSL_NAN, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);

  return s;
}
