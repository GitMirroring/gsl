/* specfunc/gamma/test_gamma_inc_P.c
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
test_gamma_inc_P(void)
{
  int s = 0;

  /* Domain errors. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (-1.0, 0.5, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.0, 0.5, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.5, -1.0, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  /* Exact x = 0 limit: P(a,0) = 0. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.1, 0.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.5, 0.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.9999999999999999, 0.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  /*
   * Taylor-x region: 0 < a < 1, 0 < x < 1.1.
   *
   * Small a.
   */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (1.0e-8, 1.0e-8, &r),
          0.99999982156536494670635569207705201,
          TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (1.0e-4, 1.0e-6, &r),
          0.99867703807315909596139559293657603,
          TEST_TOL3, GSL_SUCCESS);

  /*
   * Very small x.  These points exercise the pow(x,a) direct-P
   * evaluation used to retain accuracy when x is far below a.
   */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.25, 1.0e-14, &r),
          0.00034888228355700128012707464422634103,
          TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.75, 1.0e-14, &r),
          3.4407644396193763776316191024870257e-11,
          TEST_TOL3, GSL_SUCCESS);

  /* Interior of the Taylor-x region. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.1, 0.1, &r),
          0.82755175958585053664032521721696758,
          TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.25, 0.25, &r),
          0.7436779447314610416708665005619282,
          TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.5, 0.5, &r),
          0.68268949213708589717046509126407584,
          TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.75, 0.5, &r),
          0.52793710983467178605328824058882304,
          TEST_TOL3, GSL_SUCCESS);

  /* x > a, including P close to one. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.1, 0.5, &r),
          0.94140244589013351865712777699107035,
          TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.01, 1.0, &r),
          0.9977837653767720096568733155691142,
          TEST_TOL3, GSL_SUCCESS);

  /*
   * Upper part of the Taylor-x region, where dense validation found
   * the largest ULP errors.
   */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.75, 1.05, &r),
          0.75453241518899625667827878814416893,
          TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.9, 1.05, &r),
          0.69214101206478296539768616005281304,
          TEST_TOL4, GSL_SUCCESS);

  /* Points close to the current dispatcher boundaries. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e,
          (0.9999999999999999, 0.5, &r),
          0.39346934028736663073690752161685254,
          TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e,
          (0.5, 1.0999999999999999, &r),
          0.86198926243134041650517411437966101,
          TEST_TOL4, GSL_SUCCESS);


  /*
   * Positive P-series kernel.
   *
   * These points lie outside the Taylor and Temme regions and exercise
   * the positive series together with the rgammax scaling primitive.
   */

  /* series: small/moderate a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (1, 0.25, &r),
          0.2211992169285951317548297330216793527,
          TEST_TOL5, GSL_SUCCESS);

  /* series: lower tail. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (2.5, 0.5, &r),
          0.037434226752703631042917910102580765413,
          TEST_TOL5, GSL_SUCCESS);

  /* series: moderate a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (10, 5, &r),
          0.031828057306204811737186574180837104177,
          TEST_TOL5, GSL_SUCCESS);

  /* series: rgammax direct branch. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (20, 10, &r),
          0.0034543419758568076821662580527984202044,
          TEST_TOL5, GSL_SUCCESS);

  /* series: rgammax large-a branch. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (21, 10, &r),
          0.0015882606618580481603012858923114564642,
          TEST_TOL5, GSL_SUCCESS);

  /* series: large a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (100, 50, &r),
          0.00000000032000653245851252937786610660316508243,
          TEST_TOL5, GSL_SUCCESS);

  /* series: deep lower tail. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (1000, 300, &r),
          2.4149201482967856031764299161202925847e-221,
          TEST_TOL5, GSL_SUCCESS);

  /*
   * Temme uniform asymptotic kernel.
   *
   * Exercise both sides of x=a, the exact transition x=a, the lower
   * a boundary, and large-a cases.  All points satisfy |eta| <= 0.5.
   */

  /* Temme: lower a boundary, transition. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (12, 12, &r),
          0.53840266693638179935675347015248864341,
          TEST_TOL6, GSL_SUCCESS);

  /* Temme: x<a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (20, 18, &r),
          0.34908387201929826333694991682198185345,
          TEST_TOL6, GSL_SUCCESS);

  /* Temme: transition. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (50, 50, &r),
          0.51880831547204328189090570091126561057,
          TEST_TOL6, GSL_SUCCESS);

  /* Temme: x>a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (100, 110, &r),
          0.84172132993991290619829962098296885319,
          TEST_TOL6, GSL_SUCCESS);

  /* Temme: large a, x<a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (1000, 900, &r),
          0.00054990226571178292301303702674925853711,
          TEST_TOL6, GSL_SUCCESS);

  /* Temme: large a, x>a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (1000, 1100, &r),
          0.99894067674607002265112506662004963676,
          TEST_TOL6, GSL_SUCCESS);

  /* Temme: very large a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (100000, 100000, &r),
          0.50042052211036517669331257904382628615,
          TEST_TOL6, GSL_SUCCESS);



  /*
   * Additional tests for the updated x > a dispatcher.
   *
   * These tests exercise the Q dispatcher used by P_e() on the upper
   * side of the distribution: the small-a kernel, P complementation,
   * the continued fraction, and the large-x asymptotic kernel.  They
   * also exercise the principal Q dispatcher boundaries x = 40,
   * x = 200, and (x-a)/a = 1.5.
   *
   * All reference values below were generated with mpmath at 110-digit
   * precision from the exact binary64 input values.
   */

  /* P_e -> Q_small_a: a <= 2 and x >= 1.1. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (1.0e-12, 1.1, &r),
          0.9999999999998140090954637562090841601373578699546159218812200786646815,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (1.0e-15, 1.1, &r),
          0.9999999999999998140090954639596406675594076255388165152182615578658721,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (0.5, 1.5, &r),
          0.9167354833364495981450806795009063640671589448016793451359957626672847,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (2.0, 2.1, &r),
          0.6203850724157560783221931341748585093895731550429172939479206556759825,
          TEST_TOL6, GSL_SUCCESS);

  /* P-series / Q-complement dispatcher boundary at x = a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (10.0, 10.0, &r),
          0.5420702855281477916858351429406654182473646400324224861641016997854368,
          TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (10.0, 10.000000000000002, &r),
          0.5420702855281480419059065852072383657701783088101746311111932273128815,
          TEST_TOL5, GSL_SUCCESS);

  /* P_e -> Q P-complement branch. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (2.01, 2.01201, &r),
          0.5943041989809751455757828178481641969111555068148362379302354434279444,
          TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (11.99, 12.1099, &r),
          0.5520661932824699535896397861254902756242562441514677404580039400522753,
          TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (11.99, 12.5895, &r),
          0.6050298251544352551037869580463921020899400579303404068674605912262133,
          TEST_TOL5, GSL_SUCCESS);

  /* P_e -> Q continued-fraction branch. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (2.01, 2.814, &r),
          0.7691267436085769547172170244631922867426169110869124149989138776504187,
          TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (2.1, 2.94, &r),
          0.7709831136580132310089200829010705467277363098322397787228153099234216,
          TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (12.0, 30.0, &r),
          0.999936122974600726635009341709838434432235186267933952001872295737161,
          TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (50.0, 87.5, &r),
          0.9999947845865511028181962912811272303375267428700956838809860922736591,
          TEST_TOL5, GSL_SUCCESS);

  /* Q continued-fraction / asymptotic boundary at x = 40 for a < 12. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (10.0, 39.99999999999999, &r),
          0.9999999960740677737137811152360000523016833828511069733400796654023459,
          TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (10.0, 40.0, &r),
          0.9999999960740677737138118052871987726120092158064136085421823581330804,
          TEST_TOL5, GSL_SUCCESS);

  /* Q asymptotic boundary at x = 200 for a >= 20. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (100.0, 199.99999999999997, &r),
          0.9999999999999981561063502883976524246564988513988377874173460490457667,
          TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (100.0, 200.0, &r),
          0.9999999999999981561063502884258486381388505861514731554918610543123582,
          TEST_TOL5, GSL_SUCCESS);

  /* Q asymptotic boundary at (x-a)/a = 1.5 for a >= 20. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (20.0, 49.99999999999999, &r),
          0.9999995208642699661868613066727036674129539081701401667582405193175603,
          TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (20.0, 50.0, &r),
          0.9999995208642699661898855126400524430595791595033287733047207889828299,
          TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (20.0, 50.00000000000001, &r),
          0.9999995208642699661929097186074011999561274132741664960133352362278058,
          TEST_TOL5, GSL_SUCCESS);

  /* Large-a asymptotic upper-tail case from the Q qualification. */
  TEST_SF(s, gsl_new_sf_gamma_inc_P_e, (124.6205714888511, 215.54338189725357, &r),
          0.9999999999930090683591903734857500877601992852096545827354403449412772,
          TEST_TOL5, GSL_SUCCESS);

  return s;
}
