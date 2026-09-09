/* specfunc/log/test_log1pmx.c
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
#include <gsl/gsl_sf_log.h>

#include "../test_sf_new.h"

int
test_log1pmx(void)
{
  int s = 0;

  /*
   * Domain endpoint and negative outer region, -1 < x < -0.5.
   * This region uses the logmxp1(1+x) transformation.  Include a
   * point adjacent to -1, representative interior points, and a
   * point immediately below the transition at -0.5.
   */
  TEST_SF(s, gsl_sf_log1pmx_e, (-1.0, &r),
          GSL_NEGINF, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_log1pmx_e, (-0.9999999999999999, &r),
          -35.7368005696771015101356048997990121503646739, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (-0.99, &r),
          -3.61517018598809048873934740624514283013853235, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (-0.9, &r),
          -1.4025850929940458838581358872125661357580899, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (-0.75, &r),
          -0.636294361119890618834464242916353136151000269, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (-0.5000000000000001, &r),
          -0.1931471805599454204395345839738552623419551, TEST_TOL0, GSL_SUCCESS);

  /*
   * Legacy GSL Chebyshev region, -0.5 <= x <= -0.425.
   * Exercise both endpoints and an interior point.
   */
  TEST_SF(s, gsl_sf_log1pmx_e, (-0.5, &r),
          -0.193147180559945309417232121458176568075500134, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (-0.49, &r),
          -0.183344553263765587857724120199324094336447239, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (-0.425, &r),
          -0.128385238184786603839706596821973069451901651, TEST_TOL0, GSL_SUCCESS);

  /*
   * Julia log1pmx polynomial region, -0.425 < x < 0.7.
   * Include both sides of zero, tiny arguments, representative
   * interior points, the observed portable-Horner stress point, and
   * values immediately inside both transitions.
   */
  TEST_SF(s, gsl_sf_log1pmx_e, (-0.42499999999999993, &r),
          -0.12838523818478656280972525197923788227697014, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (-0.25, &r),
          -0.0376820724517809274392190059938274315035097109, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (-0.1, &r),
          -0.00536051565782630184429155007551089534077939241, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (-9.5367431640625e-07, &r),
          -0.000000000000454747640007250243736279163253167028009406259, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (-9.094947017729282e-13, &r),
          -4.13590306276764609563879938426956203269161793e-25, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (0.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (9.094947017729282e-13, &r),
          -4.13590306276263065307529095756889542097693024e-25, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (9.5367431640625e-07, &r),
          -0.000000000000454747061766091584485036663856396851929020697, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (0.1, &r),
          -0.00468982019567514046069470609430607525083014919, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (0.25, &r),
          -0.0268564486857902442337049096901654966253989145, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (0.5, &r),
          -0.0945348918918356180219868845356508634280095765, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (0.50841875, &r),
          -0.0973568332853746442724701581103883781385526427, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (0.69, &r),
          -0.165271471065017874171264490388870034808353747, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (0.6999999999999998, &r),
          -0.169371748937829539767364829008100109679187505, TEST_TOL0, GSL_SUCCESS);

  /*
   * Dedicated degree-16 polynomial region, 0.7 <= x <= 1.
   * Exercise the lower transition from both sides, an interior point,
   * and the upper endpoint.
   */
  TEST_SF(s, gsl_sf_log1pmx_e, (0.7, &r),
          -0.169371748937829585482430548867483229767661149, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (0.7000000000000001, &r),
          -0.169371748937829631197496268726870614891305755, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (0.85, &r),
          -0.234814360909766538865078489636682231263979674, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (0.9999999999999999, &r),
          -0.306852819440054635071616647283997951486871971, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (1.0, &r),
          -0.306852819440054690582767878541823431924499866, TEST_TOL0, GSL_SUCCESS);

  /*
   * Direct log1p(x) - x region, x > 1.  Include the value immediately
   * above the transition and arguments spanning many orders of magnitude.
   */
  TEST_SF(s, gsl_sf_log1pmx_e, (1.0000000000000002, &r),
          -0.306852819440054801605070341057483637263488714, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (2.0, &r),
          -0.901387711331890308604754763077474295352509442, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (10.0, &r),
          -7.60210472720162945593805642203487070017829315, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (50.0, &r),
          -46.0681743672756742283552201452043477597643064, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (1000000.0, &r),
          -999986.184488442036225895558718188560281421227, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (1e+100, &r),
          -1.00000000000000001590289110975991804683608086e+100, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1pmx_e, (1e+300, &r),
          -1.00000000000000005250476025520442024870446858e+300, TEST_TOL0, GSL_SUCCESS);

  /*
   * Domain error below x = -1.
   */
  TEST_SF(s, gsl_sf_log1pmx_e, (-1.0000000000000002, &r),
          GSL_NAN, TEST_TOL0, GSL_EDOM);

  return s;
}
