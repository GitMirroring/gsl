/* specfunc/gamma/test_gamma_inc_Q.c
 *
 * Unit tests for gsl_new_sf_gamma_inc_Q_e().
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

  return s;
}
