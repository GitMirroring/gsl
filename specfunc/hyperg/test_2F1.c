/* specfunc/hyperg/test_2F1.c
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
#include <gsl/gsl_test.h>
#include <gsl/gsl_sf_hyperg.h>

#include "../test_sf.h"

/* expected status only (domain errors, loss of significance) */
#define TEST_2F1_STATUS(stat, func, args, expected_status)                 do {                                                                      gsl_sf_result r;                                                        stat += test_sf_return(func args, (expected_status), #func #args);    } while (0)

int
test_2F1(void)
{
  int s = 0;

  /*
   * Regions of x: series in x (0 < x <= 1/2), Pfaff transformation
   * (-1 <= x < 0), connection formula in 1-x (1/2 < x < 1), Pfaff plus
   * connection formula in 1/(1-x) (x < -1).  Reference values from
   * 300-bit arithmetic (mpmath).
   */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, -50.0, &r), 0.0062954360187886893100, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, -3.0, &r), 0.19971166974921137706, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, -1.0, &r), 0.46076590921890255476, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, -0.7, &r), 0.55572066443896334228, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, -0.3, &r), 0.75143732100577077333, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, 0.1, &r), 1.1186584708889406040, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, 0.4, &r), 1.6980813749917111478, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, 0.6, &r), 2.5206336687353790532, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, 0.9, &r), 8.3054121721322190635, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, 0.999, &r), 146.93712897700675067, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, 0.9999999999, &r), 499988.19854245463519, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-2.75, 1.25, 0.5, -20.0, &r), 19255.772800976439566, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (3.25, -4.5, -2.75, -20.0, &r), -961438867.09624393344, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-2.75, 1.25, 0.5, -0.6, &r), 8.9820942643819004111, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (3.25, -4.5, -2.75, -0.6, &r), -230.85795237084704187, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-2.75, 1.25, 0.5, 0.3, &r), -0.33155578550006047877, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (3.25, -4.5, -2.75, 0.3, &r), 0.11729554465872049465, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-2.75, 1.25, 0.5, 0.8, &r), -0.35224886509831802961, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (3.25, -4.5, -2.75, 0.8, &r), -37.110662319464911899, TEST_TOL1, GSL_SUCCESS);

  /* x = 1: Gauss summation, c-a-b > 0 (bug #39056) */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (0.5, 0.25, 1.5, 1.0, &r), 1.1981402347355922074, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.0, 2.0, 4.5, 1.0, &r), 2.3333333333333333333, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-2.5, 0.75, 1.25, 1.0, &r), 0.23123860548545750616, TEST_TOL1, GSL_SUCCESS);

  /* c-a-b an integer, x in (1/2,1): evaluated uniformly in e = c-a-b-m */
  /* (2/pi) E(sqrt x) (bug #21835) */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-0.5, 0.5, 1.0, 0.5, &r), 0.85984660010223779136, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-0.5, 0.5, 1.0, 0.7, &r), 0.79047203432119513987, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-0.5, 0.5, 1.0, 0.9, &r), 0.70332143885152267280, TEST_TOL1, GSL_SUCCESS);
  /* sign of Gamma at a negative argument (bug #54998) */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-0.25, 0.25, 1.0, 0.9, &r), 0.92058925138209276481, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-0.25, 0.25, 1.0, 0.25, &r), 0.98334265077516523252, TEST_TOL1, GSL_SUCCESS);
  /* bug #39056 */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (3.5, -0.5, 5.0, 0.9, &r), 0.59239812843706533177, TEST_TOL1, GSL_SUCCESS);
  /* wrong sign for x >= 1/2 (bug #32306) */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-0.5, 1.5, 1.0, 0.5, &r), 0.53935260118837935667, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-0.5, 1.5, 1.0, 0.8, &r), 0.063518433378569439130, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-0.5, 1.5, 1.0, 0.9, &r), -0.23462153663932545420, TEST_TOL1, GSL_SUCCESS);
  /* GSL_EMAXITER before (bug #21835) */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.0, 13.0, 14.0, 0.999227196008978, &r), 53.464514418791908495, TEST_TOL1, GSL_SUCCESS);
  /* c-a-b near an integer */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (0.5, 0.5, 1.000000000001, 0.95, &r), 1.8515049970711781023, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (0.5, 0.5, 0.999999999, 0.95, &r), 1.8515049988230225640, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.25, 2.5, 6.75, 0.97, &r), 2.0000460688231835146, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.25, 2.5, 0.75, 0.97, &r), 74340.052438823675217, TEST_TOL1, GSL_SUCCESS);

  /* x < -1 (bug #30324), b-a an integer */
  /* log(4)/3 */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.0, 1.0, 2.0, -3.0, &r), 0.46209812037329687294, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (0.5, 0.25, 1.5, -3.0, &r), 0.86687088900112860890, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (2.0, 4.0, 3.5, -10.0, &r), 0.0055756193231225999731, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (2.0, 4.0000000001, 3.5, -10.0, &r), 0.0055756193227159658046, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (0.75, -1.25, 2.0, -1000000.0, &r), 10122873.532503759610, TEST_TOL1, GSL_SUCCESS);

  /* terminating series: polynomials, any real x (bugs #53905, #50711, #21835, #39056) */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-10.0, 2.0, 0.5, 0.5, &r), 0.077832384333932321548, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-20.0, 2.0, 0.5, 0.5, &r), 0.012059082243344950701, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-1.0, -13.0, 1.0, 0.651439, &r), 9.4687069999999998737, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-3.0, -14.0, 2.0, 0.651439, &r), 78.455361713711291748, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-1.0, -1.0, -0.5, 1.5, &r), -2.0, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-1.0, -10.0, 1.0, 0.5, &r), 6.0, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-20.0, 3.5, 2.25, -4.0, &r), 1238657831404021.2301, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-30.0, 10.5, 7.25, 0.75, &r), 1.6735813002789132498e-11, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-12.0, -15.625230111567706, 10.20713844289294, -2.5789846574473607, &r), -1.2018274053900162206, TEST_TOL1, GSL_SUCCESS);
  /* c a negative integer below the degree */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-5.0, 2.5, -7.0, 3.0, &r), 315.88448660714285714, TEST_TOL1, GSL_SUCCESS);
  /* c-a a non-positive integer: (1-x)^(c-a-b) times a polynomial */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (4.5, 1.5, -3.5, 0.75, &r), 425422262.85714285714, TEST_TOL1, GSL_SUCCESS);

  /* transformed parameters (c-a, c-b, b-c-n+1, ...) near a non-positive
     integer, where their rounding matters (value and error estimate) */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-2.0130078883480818, -5.0130078883480813, 3.6645626135354181, -36.829744673182056, &r), 266.2784542861888966, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (4.136800937994893, 2.1951381036270874, -1.8037892210681306, -334.7965332536889, &r), 6.8672298890991802915e-08, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (0.26683566950992965, -3.5917461680233851, -3.732177628893163, 0.81626975745088148, &r), 1.5554215285047883637, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-38.0, 15.820632106994474, 9.893293791793436, 1.2906655465822192, &r), 2.6897017314605569276e-14, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-34.0, -13.969783974451877, 1.4625069569488538, 3.5607011744759376, &r), 166446993222578073.1, TEST_TOL3, GSL_SUCCESS);

  /* 1-x connection formula with a+m+e = c-b or b+m+e = c-a next to a pole
     of Gamma, which x+e formed in double does not resolve */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (0.2624320771614812, 2.3499999999999974, 2.35, 0.8, &r), 1.5255700361653603772, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (0.2624320771614812, 2.34999999999, 2.35, 0.8, &r), 1.5255700361609519098, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-2.0000000000000044, 6.007339179098662, 2.35, 0.8, &r), 0.33207474662482539979, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-3.1227573348747697, 5.349999999999572, 2.35, 0.8, &r), 0.075914064485542676392, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (9.094947017729282e-13, 2.4999999999990905, 2.5, 0.8, &r), 1.0000000000014637753, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (9.094947017729282e-13, 2.4999999999990905, 2.5, 0.3, &r), 1.000000000000324394, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-1.9999999999999, 3.7, 1.2, 0.9, &r), 0.78556818181825490296, TEST_TOL1, GSL_SUCCESS);

  /* 2F1(a+1,b;a;x) (bug #66850) and c = a */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (3.5, 0.7, 2.5, 0.6, &r), 2.6967851649099269163, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (2.2, 1.3, 2.2, 0.7, &r), 4.7834624472214907140, TEST_TOL1, GSL_SUCCESS);

  /* larger parameters */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (20.0, -15.5, 35.0, 0.95, &r), 0.000064517310611484181275, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-50.5, 30.0, 20.0, 0.3, &r), 1.0387590071656132052e-11, TEST_TOL2, GSL_SUCCESS);
  /* strong cancellation: summed in double-double */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (15.396738555932103, 15.226563731581493, -14.582407970441443, -0.8865763203540584, &r), 93.256691115332563490, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.652281748192408, -19.379030960881483, -25.899690545546754, -0.9112222183977783, &r), 0.42160999572895354941, TEST_TOL1, GSL_SUCCESS);

  /* invalid arguments */
  TEST_2F1_STATUS(s, gsl_sf_hyperg_2F1_e, (1.0, GSL_NAN, 3.0, 0.5, &r), GSL_EDOM);
  TEST_2F1_STATUS(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, 1.5, &r), GSL_EDOM);         /* branch cut */
  TEST_2F1_STATUS(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, 3.5, 1.0, &r), GSL_EDOM);         /* divergent: c-a-b < 0 */
  TEST_2F1_STATUS(s, gsl_sf_hyperg_2F1_e, (1.5, 2.5, -2.0, 0.5, &r), GSL_EDOM);        /* c a pole */
  TEST_2F1_STATUS(s, gsl_sf_hyperg_2F1_e, (1.0, 1.0, 0.0, 0.0, &r), GSL_EDOM);         /* undefined also at x = 0 (bug #66850) */
  TEST_2F1_STATUS(s, gsl_sf_hyperg_2F1_e, (-5.0, 2.5, -3.0, 0.5, &r), GSL_EDOM);       /* pole before termination */

  /* no significant digit can be obtained: GSL_ELOSS, not a wrong value */
  TEST_2F1_STATUS(s, gsl_sf_hyperg_2F1_e, (78.6623289153141, 181.9965807473751, -25.13452167592419, -0.3112423999189009, &r), GSL_ELOSS);

  /*
   * F/Gamma(c); for c a non-positive integer, the limit
   * (a)_{n+1} (b)_{n+1} x^{n+1} F(a+n+1,b+n+1;n+2;x)/(n+1)! (bug #53876)
   */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (1.0, 2.0, 3.5, 0.4, &r), 0.39578371013649785288, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (1.0, 2.0, -2.5, 0.4, &r), 34.335231118802926340, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (5.0, 5.0, -1.0, 0.5, &r), 1237440.0, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (5.0, 5.0, -10.0, 0.5, &r), 68070553334784000.000, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (2.5, -1.5, -3.0, 0.3, &r), 0.13445002339815395684, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (2.5, -1.5, -3.0, -0.7, &r), 0.22983448482571312928, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (-1.0, 2.5, -3.0, 0.3, &r), 0.0, TEST_TOL0, GSL_SUCCESS);

  /* conjugate parameters a = aR + i aI, b = aR - i aI */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (1.0, 2.0, 3.0, 0.4, &r), 2.0719501020923322856, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (1.0, 2.0, 3.0, 0.9, &r), 8.3098337118941490538, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (1.0, 2.0, 3.0, 0.9999, &r), 16.983550913803348388, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (1.0, 2.0, 3.0, -0.8, &r), 0.29044828453612014974, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (1.0, 2.0, 3.0, -5.0, &r), -0.035659696106304395926, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (-2.5, 1.5, 0.75, 0.6, &r), 14.107684387440519283, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (0.5, 1.5, 1.0, 0.9, &r), 16.739572938994562151, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (0.25, 0.5, 0.50000000001, 0.95, &r), 3.2902435537654156332, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (1.0, 2.0, 3.5, 1.0, &r), 9.4498856994264078440, TEST_TOL2, GSL_SUCCESS);

  /* conjugate parameters, renormalized */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_renorm_e, (1.0, 2.0, 3.0, 0.4, &r), 1.0359750510461661428, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_renorm_e, (5.0, 5.0, -1.0, 0.5, &r), 37212765.838394834063, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_renorm_e, (5.0, 5.0, -10.0, 0.5, &r), 3890859063196462448.5, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_renorm_e, (1.5, 0.5, -2.0, -0.6, &r), -0.64652501997084340459, TEST_TOL2, GSL_SUCCESS);

  /* cases from the original specfunc test suite */
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1, 1, 1, 0.5, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (8, 8, 1, 0.5, &r), 12451584.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (8, -8, 1, 0.5, &r), 0.13671875, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (8, -8.1, 1, 0.5, &r), 0.14147385378899930422, TEST_TOL4, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (8, -8, 1, -0.5, &r), 4945.136718750000000, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (8, -8, -5.5, 0.5, &r),  -906.6363636363636364, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (8, -8, -5.5, -0.5, &r), 24565.363636363636364, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (8, 8, 1, -0.5, &r), -0.006476312098196747669, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (8, 8, 5, 0.5, &r), 4205.714285714285714, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (8, 8, 5, -0.5, &r), 0.0028489656290296436616, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (9, 9, 1, 0.99, &r), 1.2363536673577259280e+38 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (9, 9, -1.5, 0.99, &r), 3.796186436458346579e+46, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (9, 9, -1.5, -0.99, &r), 0.14733409946001025146, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (9, 9, -8.5, 0.99, &r), -1.1301780432998743440e+65, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (9, 9, -8.5, -0.99, &r), -8.856462606575344483, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (9, 9, -21.5, 0.99, &r), 2.0712920991876073253e+95, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (9, 9, -21.5, -0.99, &r), -74.30517015382249216, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (9, 9, -100.5, 0.99, &r),  -3.1867780614279306568e+262, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (9, 9, -100.5, -0.99, &r),  2.4454358338375677520, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (25, 25, 1, -0.5, &r), 5.1696944095663206268e-06, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/64.0, &r), 3.17175539044729373926, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/128.0, &r), 3.59937243502024563424, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/256.0, &r), 4.03259299524392504369, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/1024.0, &r), 4.90784159359675398250, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/65536.0, &r), 7.552266033399683914, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/16777216.0, &r), 11.08235454026043830363, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, -1.0+1.0/1024.0, &r), 0.762910940909954974527, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, -1.0+1.0/65536.0, &r), 0.762762124908845424449, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, -1.0+1.0/1048576.0, &r), 0.762759911089064738044, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 3.0, 1.0, &r), 1.6976527263135502482014268 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (1.5, -4.2, 3.0, 1.0, &r), .15583601560025710649555254 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-7.4, 0.7, -1.5, 1.0, &r), -.34478866959246584996859 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (0.1, -2.7, -1.5, 1.0, &r), 1.059766766063610122925 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (0, -2, -4, 0.5, &r), 1.0 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-10.34, 2.05, 3.05, 0.1725, &r), 0.310473552213130010351006093079548, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-9.99999999999, 2.05, 3.05, 0.1725, &r),0.32141934630197487540298837643890, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (11, -1, 11.0/2.0, 0.125 , &r), 0.75, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-0.2, 8.8, 10.0, 0.8, &r), 0.77998971427681563, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-0.2, 9.8, 11.0, 0.8, &r), 0.77574573497387233396, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (3.5, -0.5, 5.0, 0.9, &r), 0.5923981284370653465208973272, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (-1.0, -10.0, 1.0, 0.5, &r), 6.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_e, (3.23191, -4.0229, 8.02291, 0.5, &r), 0.4300243900348170646, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (1, 1, 1, 0.5, &r), 3.352857095662929028, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 1, 0.5, &r), 1.7078067538891293983e+09, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 5, 0.5, &r), 285767.15696901140627, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 1, -0.5, &r), 0.007248196261471276276, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 5, -0.5, &r), 0.00023301916814505902809, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_e, (25, 25, 1, -0.5, &r), 5.1696944095663206268e-06, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (1, 1, 1, 0.5, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (8, 8, 1, 0.5, &r), 12451584.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, 1, 0.5, &r), 0.13671875, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, 1, -0.5, &r), 4945.13671875, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, -5.5, 0.5, &r), -83081.19167659493609, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, -5.5, -0.5, &r), 2.2510895952730178518e+06, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (8, 8, 5, 0.5, &r), 175.2380952380952381, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (9, 9, -1.5, 0.99, &r), 1.6063266334913066551e+46, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (9, 9, -1.5, -0.99, &r), 0.06234327316254516616, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (5, 5, -1, 0.5, &r), 1237440.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (5, 5, -10, 0.5, &r), 68070553334784000.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_renorm_e, (5, 5, -100, 0.5, &r), 1.191184208770032773e+176, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_renorm_e, (9, 9, -1.5, 0.99, &r), 5.912269095984229412e+49, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_renorm_e, (9, 9, -1.5, -0.99, &r), 0.10834020229476124874, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_renorm_e, (5, 5, -1, 0.5, &r), 37212765.838394834063, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_renorm_e, (5, 5, -10, 0.5, &r), 3890859063196462448.5, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_hyperg_2F1_conj_renorm_e, (5, 5, -100, 0.5, &r), 1.2271985759108486647e+178, TEST_TOL2, GSL_SUCCESS);

  return s;
}
