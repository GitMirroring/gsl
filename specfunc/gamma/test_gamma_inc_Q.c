/* specfunc/gamma/test_gamma_inc_Q.c
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
test_gamma_inc_Q(void)
{
  int s = 0;

  /* Domain errors. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (-1.0, 0.5, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (0.0, 0.5, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (0.5, -1.0, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  /* Exact x = 0 limit: Q(a,0) = 1. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (0.10000000000000001, 0.0, &r),
          1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (2, 0.0, &r),
          1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (12, 0.0, &r),
          1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (1000, 0.0, &r),
          1.0, TEST_TOL0, GSL_SUCCESS);

  /*
   * small-a direct Taylor/GRAT1 Q.
   */

  /* small a, tiny x. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (1e-08, 1e-08, &r),
          1.7843463505329364e-07,
          TEST_TOL6, GSL_SUCCESS);

  /* small a, x=1 worst-area. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (9.9999999999999998e-13, 1, &r),
          2.1938393439574474e-13,
          TEST_TOL6, GSL_SUCCESS);

  /* small a, near x=1.1. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (9.9999999999999998e-13, 1.099999999999, &r),
          1.8599090453654639e-13,
          TEST_TOL6, GSL_SUCCESS);

  /* Taylor interior. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (0.10000000000000001, 0.10000000000000001, &r),
          0.17244824041414947,
          TEST_TOL6, GSL_SUCCESS);

  /* Taylor interior. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (0.5, 0.5, &r),
          0.31731050786291409,
          TEST_TOL6, GSL_SUCCESS);

  /* Taylor upper. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (0.75, 1.05, &r),
          0.24546758481100375,
          TEST_TOL6, GSL_SUCCESS);

  /* GRAT1 just above 1.1. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (0.5, 1.1000000000000001, &r),
          0.13801073756865953,
          TEST_TOL6, GSL_SUCCESS);

  /* GRAT1. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (1, 2, &r),
          0.1353352832366127,
          TEST_TOL6, GSL_SUCCESS);

  /* GRAT1 upper a boundary. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (2, 2, &r),
          0.40600584970983805,
          TEST_TOL6, GSL_SUCCESS);

  /* GRAT1 upper a boundary upper tail. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (2, 10, &r),
          0.00049939922738733333,
          TEST_TOL6, GSL_SUCCESS);

  /*
   * P-series complement: x <= a.
   */

  /* just above small-a boundary. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (2.0000000000000009, 2, &r),
          0.40600584970983833,
          TEST_TOL6, GSL_SUCCESS);

  /* transition lower side. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (2.5, 2.5, &r),
          0.41588018699550794,
          TEST_TOL6, GSL_SUCCESS);

  /* moderate. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (4, 4, &r),
          0.43347012036670896,
          TEST_TOL6, GSL_SUCCESS);

  /* moderate. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (10, 5, &r),
          0.96817194269379514,
          TEST_TOL6, GSL_SUCCESS);

  /* near a=12. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (11.99, 11.869999999999999, &r),
          0.47537749071619217,
          TEST_TOL6, GSL_SUCCESS);

  /*
   * moderate-a P complement, x>a and Q>=0.34.
   */

  /* near a=2 boundary. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (2.0099999999999998, 2.0120100000000001, &r),
          0.40569580101902475,
          TEST_TOL6, GSL_SUCCESS);

  /* moderate. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (4, 4.0999999999999996, &r),
          0.41418154152825021,
          TEST_TOL6, GSL_SUCCESS);

  /* moderate. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (8, 8.1999999999999993, &r),
          0.42540929947339051,
          TEST_TOL6, GSL_SUCCESS);

  /* near a=12. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (11.99, 12.1099, &r),
          0.44793380671753014,
          TEST_TOL6, GSL_SUCCESS);

  /* historical Pcomp high error. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (11.99, 12.589499999999999, &r),
          0.39497017484556485,
          TEST_TOL6, GSL_SUCCESS);

  /*
   * Julia DD continued fraction.
   */

  /* near a=2 boundary CF. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (2.0099999999999998, 2.8140000000000001, &r),
          0.23087325639142298,
          TEST_TOL6, GSL_SUCCESS);

  /* moderate upper. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (2.1000000000000001, 2.9399999999999999, &r),
          0.2290168863419868,
          TEST_TOL6, GSL_SUCCESS);

  /* moderate upper. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (5, 8, &r),
          0.099632400487046011,
          TEST_TOL6, GSL_SUCCESS);

  /* large-a upper. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (15, 30, &r),
          0.00092068239614866631,
          TEST_TOL6, GSL_SUCCESS);

  /* large-a upper. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (50, 87.5, &r),
          5.2154134488971821e-06,
          TEST_TOL6, GSL_SUCCESS);

  /* below x=200 boundary. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (100, 190, &r),
          2.6632899618434939e-13,
          TEST_TOL6, GSL_SUCCESS);

  /*
   * Temme.
   */

  /* lower a boundary transition. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (12, 12, &r),
          0.46159733306361822,
          TEST_TOL6, GSL_SUCCESS);

  /* x<a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (20, 18, &r),
          0.65091612798070175,
          TEST_TOL6, GSL_SUCCESS);

  /* transition. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (50, 50, &r),
          0.48119168452795674,
          TEST_TOL6, GSL_SUCCESS);

  /* x>a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (100, 110, &r),
          0.15827867006008708,
          TEST_TOL6, GSL_SUCCESS);

  /* large a x<a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (1000, 900, &r),
          0.99945009773428817,
          TEST_TOL6, GSL_SUCCESS);

  /* large a x>a. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (1000, 1100, &r),
          0.0010593232539299773,
          TEST_TOL6, GSL_SUCCESS);

  /* very large a transition. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (100000, 100000, &r),
          0.4995794778896348,
          TEST_TOL6, GSL_SUCCESS);

  /* known Temme difficult. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (381.79190914000145, 597.3106540677893, &r),
          1.4603195565153375e-21,
          TEST_TOL6, GSL_SUCCESS);

  /*
   * large-x asymptotic.
   */

  /* moderate a x=40 boundary. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (3.6299999999999999, 40, &r),
          1.9290667987644042e-14,
          TEST_TOL6, GSL_SUCCESS);

  /* moderate a far upper. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (10, 100, &r),
          1.1253473960842733e-31,
          TEST_TOL6, GSL_SUCCESS);

  /* a=15 upper. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (15, 150, &r),
          2.6480484856308778e-46,
          TEST_TOL6, GSL_SUCCESS);

  /* x=200 boundary. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (100, 200, &r),
          1.8438936497115741e-15,
          TEST_TOL6, GSL_SUCCESS);

  /* s=1.5 boundary. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (20, 50, &r),
          4.7913573003381007e-07,
          TEST_TOL6, GSL_SUCCESS);

  /* known asym difficult. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (124.6205714888511, 215.54338189725357, &r),
          6.9909316408096273e-12,
          TEST_TOL6, GSL_SUCCESS);

  /* deep upper. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (228.10573094666549, 1092.3258664353218, &r),
          4.7610010498191436e-223,
          TEST_TOL6, GSL_SUCCESS);

  /* far-tail DD prefactor. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (15, 701, &r),
          2.9373826275874304e-276,
          TEST_TOL6, GSL_SUCCESS);

  /* Known small-a worst cases from exhaustive qualification. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (1e-15, 1.0, &r),
          0.0000000000000002193839343955205151995209016738667675472749624906021601604058850730977,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (1e-15, 1.0999999999999999, &r),
          0.0000000000000001859909045360404141006990777005730889214240731552738712999107388172995,
          TEST_TOL6, GSL_SUCCESS);

  /* Taylor / GRAT1 boundary at x = 1.1. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (0.5, 1.0999999999999999, &r),
          0.1380107375686595834948258856203389855742051734445364985421586772368083,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (0.5, 1.1, &r),
          0.1380107375686595437349503278468820426552714475645287367936163742093227,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (1e-15, 1.1, &r),
          0.0000000000000001859909045360403469077642861486071509569826048124350227155455081377325,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (1e-15, 1.1, &r),
          0.0000000000000001859909045360403469077642861486071509569826048124350227155455081377325,
          TEST_TOL6, GSL_SUCCESS);

  /* Small-a dispatcher boundary at a = 2 on the upper side. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (2.0, 2.1, &r),
          0.3796149275842438988375439109755415514383948216212482822343301812350767,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (2.0000000000000004, 2.1, &r),
          0.3796149275842440259474293729505075776885926910374229866180099522629085,
          TEST_TOL6, GSL_SUCCESS);

  /* P-complement / CF crossover Q = 0.34 at a = 3. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (3.0, 3.3986519284532797, &r),
          0.3400000000000000810115389722380825780645659556813111565058763671759007,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (3.0, 3.39865192845328, &r),
          0.3399999999999999953002898130928835147735890967885963382071745255016641,
          TEST_TOL6, GSL_SUCCESS);

  /* P-complement / CF crossover Q = 0.34 at a = 5. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (5.0, 5.614217848114592, &r),
          0.3400000000000001039793996593399417034741005594873160697085494319913876,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (5.0, 5.614217848114593, &r),
          0.3399999999999999699431753704844219265153909706816707814888254906634789,
          TEST_TOL6, GSL_SUCCESS);

  /* P-complement / CF crossover Q = 0.34 at a = 10. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (10.0, 11.004707291055073, &r),
          0.3400000000000001408199862462630347256444676416601113340038612221730333,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (10.0, 11.004707291055075, &r),
          0.3399999999999999482050385099544272742411115826926357989450657182470458,
          TEST_TOL6, GSL_SUCCESS);

  /* CF / asymptotic boundary at x = 40 for a < 12. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (3.63, 39.99999999999999, &r),
          0.00000000000001929066798764417101492936506930472554840166017040062674120267279832098,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (3.63, 40.00000000000001, &r),
          0.0000000000000192906679876439144456394547127513079553205926631808196012128489254564,
          TEST_TOL6, GSL_SUCCESS);

  /* CF / asymptotic boundary at x = 200. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (100.0, 199.99999999999997, &r),
          0.000000000000001843893649711600864181415584100472216951049628983108132053164594526883,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (100.0, 200.00000000000003, &r),
          0.000000000000001843893649711547438542306715110632954436657136349785547851840119662165,
          TEST_TOL6, GSL_SUCCESS);

  /* CF / asymptotic boundary at s = (x-a)/a = 1.5. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (20.0, 49.99999999999999, &r),
          0.0000004791357300338122633149414895742769890730990334786615946410103484564934,
          TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (20.0, 50.00000000000001, &r),
          0.0000004791357300338079656597784055490702223443811616345638865176316554942969,
          TEST_TOL6, GSL_SUCCESS);

  /* Known asymptotic worst case from exhaustive qualification. */
  TEST_SF(s, gsl_new_sf_gamma_inc_Q_e, (592.2792169054019, 1003.6433369303401, &r),
          2.40133782217006808872686439537317301485415191398539929285230021634865e-45,
          TEST_TOL6, GSL_SUCCESS);

  return s;
}
