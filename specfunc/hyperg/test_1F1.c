/* specfunc/hyperg/test_1F1.c
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
#include <gsl/gsl_test.h>
#include <gsl/gsl_sf_hyperg.h>

#include "../test_sf.h"

/* expected status only (domain errors) */
#define TEST_1F1_STATUS(stat, func, args, expected_status)                 \
  do {                                                                    \
    gsl_sf_result r;                                                      \
    stat += test_sf_return(func args, (expected_status), #func #args);    \
  } while (0)

/*
 * Reference values: the Maclaurin series summed in 200-bit arithmetic
 * with the precision raised to cover its cancellation (mpmath).  The
 * tests carried over from specfunc/test_hyperg.c keep their arguments;
 * some of their reference values were inaccurate, and some tests had
 * been disabled because the previous implementation failed them.
 */

int
test_1F1(void)
{
  int s = 0;

  /* integer parameters */
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (1, 1, 0.5, &r), 1.6487212707001281468, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (1, 2, 500.0, &r), 2.8071844357056748215e+214, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (1, 2, -500.0, &r), 0.002, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (8, 1, 0.5, &r), 13.108875178030540372, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 1, 1.0, &r), 131.63017574352619931, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 1, 10.0, &r), 8.5146254765462807963e+9, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 1, 100.0, &r), 1.567136364680035332e+56, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 20, 1.0, &r), 1.6585618002669675465, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 20, 10.0, &r), 265.26686430340188871, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 20, 100.0, &r), 3.6404773550632271294e+34, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 1.0, &r), 1.1056660194025527099, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 10.0, &r), 2.8491063634727594206, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 40.0, &r), 133.85880835831230986, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 80.0, &r), 3.1036116228011433406e+5, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 100.0, &r), 8.0321713367541682821e+7, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 500.0, &r), 7.6339612025287314256e+123, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 1, 1.0, &r), 6.8928427290464699646e+7, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 1, 10.0, &r), 2.4175917112200409098e+28, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 1, 100.0, &r), 1.9303216896309102993e+110, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 1.0, &r), 1.6497469106162459226, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 10.0, &r), 157.93286197349321981, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 100.0, &r), 2.181957750125507524e+24, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 400.0, &r), 3.7289755299265733002e+119, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 10.0, &r), 12.473087623658878813, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 100.0, &r), 9.0712303768185502414e+11, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 150.0, &r), 7.1609495157421707754e+18, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 200.0, &r), 2.7406690412731576823e+26, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 300.0, &r), 6.1751106134732761934e+43, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 400.0, &r), 1.180741766271137144e+64, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 600.0, &r), 2.407607635488888603e+112, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 1, -1.0, &r), 0.1139485482464454281, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 1, -10.0, &r), 6.7155063653961278629e-4, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 1, -100.0, &r), -4.2081385374802698683e-32, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 50, -1.0, &r), 0.82000619607937954587, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -10.0, &r), 0.38378859043466242587, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -100.0, &r), 8.4601434014641890614e-4, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -500.0, &r), 1.1090822141973655929e-8, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -10000.0, &r), 5.1737835080882722922e-21, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (50, 1, -90.0, &r), -1.6624258547648311554e-21, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (50, 1, -100.0, &r), 4.0696617751220482035e-24, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (50, 1, -110.0, &r), 1.0072444993946236025e-25, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 10, -100.0, &r), -2.7819353611733941962e-37, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 1, -90.0, &r), 7.501705041159802854e-22, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 1, -100.0, &r), 6.3051288931522911871e-25, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (100, 1, -110.0, &r), -7.0071221154224397547e-26, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (200, 50, -1.0, &r), 0.016087060191732290813, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (200, 50, -300.0, &r), -4.2949759797064214707e-121, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -1.0, &r), 0.13397521083325179687, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -10.0, &r), 5.835134393749807387e-10, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -100.0, &r), 4.888460453078914804e-74, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -500.0, &r), -1.4478509059582015053e-195, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-1, 1, 2.0, &r), -1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-1, -2, 2.0, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-2, -3, 2.0, &r), 3.0, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, 1, 1.0, &r), 0.41894593253968253968, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, 1, 10.0, &r), 27.984126984126984127, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, 1, 100.0, &r), 9.0512837954295714286e+12, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-100, 20, 1.0, &r), 0.0020203016320697069566, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, 1.0, &r), 1.6379141878548080173, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, 10.0, &r), 78.652024045212899702, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, 100.0, &r), 4.4161697132626243153e+8, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, 1.0, &r), 1.1046713999681950919, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, 10.0, &r), 2.6035952191039006838, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, 100.0, &r), 1151.6852040836932392, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, 1.0, &r), 1.6476859702535324743, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, 10.0, &r), 139.3802682954068727, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, 100.0, &r), 1.1669433576237933752e+19, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, -1.0, &r), 0.60255495611480357352, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, -10.0, &r), 0.0035720118279457804742, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, -100.0, &r), 4.9303272262405369733e+7, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, -1.0, &r), 0.90442397250313891333, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, -10.0, &r), 0.35061515251367212539, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, -100.0, &r), 4.4148250205403363517e-7, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, -1.0, &r), 0.60614979396289526293, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, -10.0, &r), 0.0063278543908877662249, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, -100.0, &r), 4.3411179500733645697e-25, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-1, -1, 0.1, &r), 1.1000000000000000056, TEST_TOL0, GSL_SUCCESS);

  /* real parameters */
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.0, 1.5, 1.0, &r), 2.0300784692787049755, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.0, 1.5, 10.0, &r), 6172.8595610784068551, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.0, 1.5, 100.0, &r), 2.3822817898485692114e+42, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.0, 1.5, 500.0, &r), 5.5628953517235135812e+215, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.5, 2.5, 1.0, &r), 1.8834451238277954398, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.5, 2.5, 10.0, &r), 3128.7352996840916381, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 1.1, 1.0, &r), 110.17623733873887904, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 1.1, 10.0, &r), 6.1466579752683836998e+9, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 1.1, 100.0, &r), 9.3318338972303080967e+55, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 1.1, 500.0, &r), 4.5194033687957131718e+235, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 50.1, 2.0, &r), 1.5001295507968071614, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 50.1, 10.0, &r), 8.7133858492650443415, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 50.1, 100.0, &r), 5.9094239322733733367e+18, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 50.1, 500.0, &r), 9.7400606184571666497e+165, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 1.1, 1.0, &r), 5.1835310671168077544e+7, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 1.1, 10.0, &r), 1.603264911009697372e+28, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 1.1, 100.0, &r), 1.1045151213192274661e+110, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 50.1, 1.0, &r), 7.2229531332166033628, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 50.1, 10.0, &r), 1.0998696410887166862e+8, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 50.1, 100.0, &r), 7.2353048623222701485e+63, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.0, 1.5, -1.0, &r), 0.53807950691276841914, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.0, 1.5, -10.0, &r), 0.053037580992901644848, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.0, 1.5, -100.0, &r), 0.0050253847187598528033, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.0, 1.5, -500.0, &r), 0.0010010030151059555322, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.0, 1.1, -500.0, &r), 2.0036137599690226057e-4, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 1.1, -1.0, &r), 0.072276456489359348919, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 1.1, -10.0, &r), 3.1924154096955855995e-4, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 1.1, -100.0, &r), -8.2934253161231640204e-16, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (10.0, 1.1, -500.0, &r), -3.4003792167077035166e-23, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (50.0, 1.1, -90.0, &r), -7.8431294118029157482e-22, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (50.0, 1.1, -100.0, &r), 4.632883869540640306e-24, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (50.0, 1.1, -110.0, &r), 5.6426846513053069792e-26, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 1.1, -1.0, &r), 0.081163734409604191647, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 1.1, -10.0, &r), 2.5945610092231559649e-4, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 1.1, -50.0, &r), 2.4284830988994064948e-13, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 1.1, -90.0, &r), 2.4468224638378396569e-22, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 1.1, -99.0, &r), 1.05070962726176047e-23, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 1.1, -100.0, &r), 1.8315497474210143069e-24, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 1.1, -101.0, &r), -2.3916306291344435243e-24, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 1.1, -110.0, &r), -4.5175819860377304495e-26, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, 10.1, -220.0, &r), -4.2961303000216930026e-64, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-10.0, -10.1, 10.0, &r), 1.095960320463306408e+4, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-10.0, -10.1, 1000.0, &r), 2.0942691895502263477e+23, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-10.0, -100.1, 10.0, &r), 2.6012036337980079419, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-1000.0, -1000.1, 10.0, &r), 2.2004341698908626609e+4, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-1000.0, -1000.1, 200.0, &r), 7.066514294896209194e+86, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-8.1, -10.1, -10.0, &r), 1.846968527634719994e-4, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-8.1, -1000.1, -10.0, &r), 0.921828018508003607, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-10.0, -5.1, 1.0, &r), 16.936141866089636273, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-10.0, -5.1, 10.0, &r), 7.7153403495438500776e+5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-10.0, -5.1, 100.0, &r), 2.2733956505085043617e+17, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, -50.1, -1.0, &r), 0.13854540373629276338, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, -50.1, -10.0, &r), -9.1422603143532253404e+19, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, -50.1, -100.0, &r), -1.7437371339223674694e+87, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, -50.1, 1.0, &r), 7.5168317481703507286, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, -50.1, 10.0, &r), 1.0551632286359513663e+11, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, -50.1, 50.0, &r), -7.5647556009402111986e+36, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, -50.1, 100.0, &r), 4.2187769626759775901e+55, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-10.5, -8.1, 0.1, &r), 1.1387201443786421872, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-10.5, -11.1, 1.0, &r), 2.5682766147138453174, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.5, -80.1, 10.0, &r), 3.5514545173052242003e+5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.5, -102.1, 10.0, &r), 1.8678558725244375766e+4, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.5, -500.1, 10.0, &r), 7.3422090111014536622, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.5, -500.1, 100.0, &r), 1.207744307536716803e+8, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-500.5, -80.1, 2.0, &r), 7.7405785413253508381e+5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, -10.1, 1.0, &r), -2.1213846338338643907e+12, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, -10.1, 10.0, &r), -6.6248493461451325087e+39, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, -10.1, 100.0, &r), -1.2413466759089202619e+129, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, -10.1, -1.0, &r), 3.4456294053055615764e+4, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, -10.1, -10.0, &r), -7.8092242514677371776e+7, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (100.0, -10.1, -100.0, &r), -5.2140654527540013557e-7, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 1.1, 1.0, &r), 0.21519810496314433931, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 1.1, 10.0, &r), 8.1961237155978663877, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 1.1, 100.0, &r), -1.4612966715976522654e+20, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 20.1, 1.0, &r), 0.0021267655527278471746, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 20.1, 10.0, &r), 2.0908665169032182314e-11, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 20.1, 100.0, &r), -0.04159447537001313131, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 1.1, -1.0, &r), 2.1214770215694680221e+7, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 1.1, -10.0, &r), 1.0258848879387569246e+24, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 1.1, -100.0, &r), 1.1811367147091755133e+67, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 50.1, -1.0, &r), 6.9652593172714270236, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 50.1, -10.0, &r), 1.0690052487716994723e+7, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-100.0, 50.1, -100.0, &r), 6.8896444357770894225e+36, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-2.05, 1.0, 5.05, &r), 3.7939338951678502897, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-26.0, 2.0, 100.0, &r), 1.4447867811074369871e+19, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-26.1, 2.0, 100.0, &r), 1.341557199575984828e+19, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.2, 1.1e-15, 1.5, &r), 8.2545031596724290466e+15, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.0, 1000000.5, 800000.5, &r), 4.9999225050994475524, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.0, 1000000.5, 1001000.5, &r), 3480.3699557431856166, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.1, 1000000.5, 1001000.5, &r), 7304.6126942641397987, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (0.9, 1000000.5, 1001000.5, &r), 1645.4879293475413734, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-1.1, 1000000.5, 1001000.5, &r), -5.3006648869744553677e-4, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1.5, 1000000.5, 800000.5, &r), 11.180012889778946505, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-1.5, 1.5, -100.0, &r), 456.44010011787485545, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-1.5, 1.5, 99.0, &r), 4.1336043601464330976e+36, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-1.5, 1.5, 100.0, &r), 1.0893724312430935129e+37, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-1.5, 1.5, 709.0, &r), 8.73968041602649e+298, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-1.5, 1.5, 710.0, &r), 2.3656318721741789817e+299, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-37.8, 2.01, 103.58, &r), -6.2192721100846497496e+17, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-1.0, -1.0, 0.1, &r), 1.1000000000000000056, TEST_TOL0, GSL_SUCCESS);

  /* bug #66826 (cases as reported; values for the terminating polynomial
     where a and b are non-positive integers) */
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (5, 5, -8.0, &r), 3.3546262790251183882e-4, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (5, 5, -2.0, &r), 0.13533528323661269189, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (5, 5, -1.0, &r), 0.3678794411714423216, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (5, 5, 1.0, &r), 2.7182818284590452354, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (5, 5, 2.0, &r), 7.3890560989306502272, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (5, 6, -3.0, &r), 0.091228027395668115414, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (5, 6, 7.0, &r), 480.54877204888424963, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (5, 6, 8.0, &r), 1204.4605636118311266, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (1, 1, -1.0, &r), 0.3678794411714423216, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-1, -1, -1.0, &r), 0.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-1, -1, 1.0, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (2, 2, 1.0, &r), 2.7182818284590452354, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (2, 3, -1.0, &r), 0.52848223531423071362, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (2, 3, 1.0, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-1, 1, 8.0, &r), -7.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-1, 1, -1.0, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-2, 2, 1.0, &r), 0.16666666666666666667, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-2, -2, 2.0, &r), 5.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-2, -2, 1.0, &r), 2.5, TEST_TOL0, GSL_SUCCESS);

  /*
   * a, b non-positive integers with b <= a: the terminating polynomial
   * sum_{k=0}^{-a} (a)_k/(b)_k x^k/k! (as DLMF, Mathematica and mpmath; the
   * previous implementation returned e^x M(b-a,b,-x) for some of these)
   */
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (-3, -5, 2.5, &r), 3.6979166666666666667, TEST_TOL0, GSL_SUCCESS);

  /* bug #28267: a << 0, x >> 0 (oscillatory region) */
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-37.8, 2.0, 103.5, &r), 3.5275730880693920985e+17, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-37.8, 2.0, 103.58, &r), -9.2455271368295069688e+17, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-37.8, 2.0, 103.57946018789922, &r), -9.1559757070422175517e+17, TEST_TOL3, GSL_SUCCESS);

  /* recurrences in a and b, Kummer transformation with an inexact b - a,
     large parameters */
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-300.3, 1.5, 12.992632226094093, &r), -3.8334000866796407679, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-139.07991113016945, 1.5, 487.36386299015754, &r), 1.3350489124648239991e+102, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-40.3, 1.7, 37.460317460317455, &r), -2.4805098441821372308e+5, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-310.9485358798302, 10.54582202411016, 6.734419054687314, &r), 1.7840837738248241915e-11, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (28.376190476190477, -29.33809523809524, -25.0, &r), -9.7826210759951095824e+7, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (1946.4513995365796, 89.43692354089157, -1.7014409523195615, &r), -1.2761832381096384841e-22, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-20.345459483508215, -1657.2760884319734, 90.75, &r), 2.9593650166183391667, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (517.168123737633, 10.0, 200.0, &r), 1.1335459801804294726e+307, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-4.519851867770328, -1.6719001628409487, 10.0, &r), -915.50685988800754978, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_e, (-11.0, 20.0, 7.5, &r), 1.0257392929852525468e-4, TEST_TOL0, GSL_SUCCESS);

  /* M(a,b,x)/Gamma(b), including b a non-positive integer (DLMF 13.2.5) */
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_renorm_e, (0.37, 1.71, 3.5, &r), 3.8819963479484633477, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_renorm_e, (2.5, -3.5, 1.2, &r), 106.6354763935302966, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_renorm_e, (1.5, -3.0, 2.0, &r), 344.88294749634114525, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_renorm_e, (-2.0, -3.0, 2.0, &r), 0.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_renorm_e, (1.0, 0.0, 0.5, &r), 0.82436063535006407342, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_renorm_e, (-0.5, -10.0, -3.0, &r), 107.21471843869280632, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_renorm_e, (4.0, 1e-08, 2.0, &r), 374.37883778066845525, TEST_TOL1, GSL_SUCCESS);

  /* bugs #55687, #58032, #58062: NaN arguments and b a non-positive integer */
  TEST_1F1_STATUS(s, gsl_sf_hyperg_1F1_e, (1.0, GSL_NAN, -1.0, &r), GSL_EDOM);
  TEST_1F1_STATUS(s, gsl_sf_hyperg_1F1_e, (GSL_NAN, 1.0, -1.0, &r), GSL_EDOM);
  TEST_1F1_STATUS(s, gsl_sf_hyperg_1F1_e, (1.0, 2.0, GSL_NAN, &r), GSL_EDOM);
  TEST_1F1_STATUS(s, gsl_sf_hyperg_1F1_e, (2.0, 0.0, 0.0, &r), GSL_EDOM);
  TEST_1F1_STATUS(s, gsl_sf_hyperg_1F1_int_e, (2, 0, 0.0, &r), GSL_EDOM);
  TEST_1F1_STATUS(s, gsl_sf_hyperg_1F1_e, (2.5, -3.0, 1.0, &r), GSL_EDOM);
  TEST_1F1_STATUS(s, gsl_sf_hyperg_1F1_renorm_e, (1.0, GSL_NAN, 1.0, &r), GSL_EDOM);
  /* M(0,b,x) = 1 for all x (bug #58062, not a bug) */
  TEST_SF_RESULT(s, gsl_sf_hyperg_1F1_int_e, (0, 1, GSL_NAN, &r), 1.0, TEST_TOL0, GSL_SUCCESS);

  return s;
}
