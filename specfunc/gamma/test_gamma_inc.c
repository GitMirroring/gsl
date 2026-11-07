/* specfunc/gamma/test_gamma_inc.c
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

#include "../test_sf.h"

int
test_gamma_inc(void)
{
  int s = 0;

  /* Domain errors: x < 0, and x = 0 for a <= 0 (divergent integral). */
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (1.0, -1.0, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.5, 0.0, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (0.0, 0.0, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  /* Gamma(a,0) = Gamma(a) for a > 0. */
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (4.0, 0.0, &r),
          6.0, TEST_TOL0, GSL_SUCCESS);

  /* Overflow and underflow. */
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-300.0, 0.001, &r),
          GSL_POSINF, GSL_POSINF, GSL_EOVRFLW);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (1.0, 800.0, &r),
          0.0, TEST_TOL0, GSL_EUNDRFLW);

  /* Recurrence: a <= 0, x < 1/2 (including a just below an integer). */
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-1e-10, 0.25, &r),
          1.044282634497776086685133390637808803097,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.5, 0.001, &r),
          5.976388051594219632350583920886792577545e+1,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.9999999999, 0.3, &r),
          1.563717417158352381128731626468308552213,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-1.0000000001, 0.4999, &r),
          6.535303975952015344344723000425528179735e-1,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-1.0, 1e-06, &r),
          9.999857617046069825939327402560036889169e+5,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-3.0, 0.1, &r),
          2.877360907483771821177057126857290035585e+2,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-7.5, 0.01, &r),
          1.318039257856408542069631771776075330701e+14,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-25.25, 0.2, &r),
          1.433198873610303755921782730540177228213e+16,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-100.0, 0.45, &r),
          3.029330123420625434810003334923301909011e+32,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-300.5, 0.49, &r),
          2.539330514868830708153846358066026180366e+90,
          TEST_TOL1, GSL_SUCCESS);

  /* a = 0: Gamma(0,x) = E1(x). */
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (0.0, 1e-08, &r),
          1.784346508905083256561485898204520147986e+1,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (0.0, 0.3, &r),
          9.056766516758467398461090442311552443662e-1,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (0.0, 0.75, &r),
          3.403408129112300078065357686796446251898e-1,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (0.0, 5.0, &r),
          1.148295591275325797330561969819722076266e-3,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (0.0, 700.0, &r),
          1.406518766234032922774410685110643005117e-307,
          TEST_TOL1, GSL_SUCCESS);

  /* Continued fraction and asymptotic expansion, x >= 1/2, x > a. */
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.25, 0.5, &r),
          5.712673030999203793826448389339991874966e-1,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-2.5, 1.0, &r),
          9.655664863127516026409178753300635318331e-2,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-50.0, 30.0, &r),
          1.62171926313707205899995654825337129981e-89,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (0.5, 2.0, &r),
          8.064711796031769078862607302130517570136e-2,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (1.5, 4.0, &r),
          4.077681246780469426903874830300103859855e-2,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (5.0, 60.0, &r),
          1.214411055802217745436665029992304454012e-19,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (30.0, 90.0, &r),
          5.652235162325394271836978455806863493497e+17,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (150.0, 1000.0, &r),
          5.963475480925828848198187522868988420224e+12,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (1e-05, 700.0, &r),
          1.406610931463870547047658946589444404054e-307,
          TEST_TOL1, GSL_SUCCESS);

  /* Gamma(a) Q(a,x), x <= a or near the transition (Temme). */
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (1e-08, 0.1, &r),
          1.822923938566501501184481387484165529097,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (0.4, 0.3, &r),
          7.952320599257721747604732700538973709038e-1,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (3.0, 2.0, &r),
          1.353352832366126918939994949724844034076,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (20.0, 21.0, &r),
          4.674368351547204627239846001725071793468e+16,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (100.0, 95.0, &r),
          6.370850531210829575333620302396088603108e+155,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (170.5, 100.0, &r),
          5.562092414042237649427081943589917548372e+305,
          TEST_TOL1, GSL_SUCCESS);

  /* Gamma(a) overflows but Gamma(a,x) does not. */
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (171.8, 200.0, &r),
          8.549337615668314927260698432212549075895e+306,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (175.0, 500.0, &r),
          4.556048177524615848267617858151662494557e+252,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (500.0, 3500.0, &r),
          3.357727708057382509115685790071158244146e+248,
          TEST_TOL1, GSL_SUCCESS);

  /* Near the underflow threshold. */
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (53.323558148631186, 1056.5870930796477, &r),
          2.363609961735573306272581288402176223408e-301,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (0.5, 740.0, &r),
          1.538773092722771099724092856681373050309e-323,
          TEST_TOL1, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-3.0, 714.0, &r),
          3.137069944989783126755117861303906007085e-322,
          TEST_TOL1, GSL_SUCCESS);


  /* cases carried over from the original specfunc test suite */
  /* non-normalized "Q" function */
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-1.0/1048576.0, 1.0/1048576.0, &r), 13.285819596290624271, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.001, 1.0/1048576.0, &r), 13.381275128625328858, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-1.0,   1.0/1048576.0, &r), 1.0485617142715768655e+06, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.00001,0.001, &r), 6.3317681434563592142, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.0001,0.001, &r), 6.3338276439767189385, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.001, 0.001, &r), 6.3544709102510843793, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-1.0,   0.001, &r), 992.66896046923884234, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-3.5,   0.001, &r), 9.0224404490639003706e+09, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-10.5,  0.001, &r), 3.0083661558184815656e+30, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.001, 0.1,   &r), 1.8249109609418620068, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.5,   0.1,   &r), 3.4017693366916154163, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-10.0,  0.1,   &r), 8.9490757483586989181e+08, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-10.5,  0.1,   &r), 2.6967403834226421766e+09, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.001, 1.0,   &r), 0.21928612679072766340, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.5,   1.0,   &r), 0.17814771178156069019, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-1.0,   1.0,   &r), 0.14849550677592204792, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-1.0,   10.0,  &r), 3.8302404656316087616e-07, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.001, 10.0,  &r), 4.1470562324807320961e-06, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-0.5,   10.0,  &r), 1.2609042613241570681e-06, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-1.0,   10.0,  &r), 3.8302404656316087616e-07, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-10.5,  10.0,  &r), 6.8404927328441566785e-17, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-100.0, 10.0,  &r), 4.1238327669858313997e-107, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (-200.0, 10.0,  &r), 2.1614091830529343423e-207, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (  0.0,     0.001, &r), 6.3315393641361493320, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (  0.001,   0.001, &r), 6.3087159394864007261, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (  1.0,     0.001, &r), 0.99900049983337499167, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, ( 10.0,     0.001, &r), 362880.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (  0.0,     1.0,   &r), 0.21938393439552027368, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (  0.001,   1.0,   &r), 0.21948181320730279613, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (  1.0,     1.0,   &r), 0.36787944117144232160, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, ( 10.0,     1.0,   &r), 362879.95956592242045, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (100.0,     1.0,   &r), 9.3326215443944152682e+155, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (  0.0,   100.0, &r), 3.6835977616820321802e-46, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (  0.001, 100.0, &r), 3.7006367674063550631e-46, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (  1.0,   100.0, &r), 3.7200759760208359630e-44, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, ( 10.0,   100.0, &r), 4.0836606309106112723e-26, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_gamma_inc_e, (100.0,   100.0, &r), 4.5421981208626694294e+155, TEST_TOL1, GSL_SUCCESS);
  return s;
}
