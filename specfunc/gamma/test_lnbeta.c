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
#include <gsl/gsl_test.h>

#include "../test_sf.h"

/* log|B(a,b)| and the sign of B(a,b) */
#define TEST_LNBETA_SGN(stat, a, b, expected, expected_sgn, tol)               \
  do {                                                                          \
    double sgn_ = 0.0;                                                          \
    TEST_SF_RESULT(stat, gsl_sf_lnbeta_sgn_e, (a, b, &r, &sgn_), expected, tol, GSL_SUCCESS); \
    gsl_test(sgn_ != (expected_sgn), "gsl_sf_lnbeta_sgn_e(%g, %g) sign %g",      \
             (double) (a), (double) (b), sgn_);                                 \
    stat += (sgn_ != (expected_sgn));                                           \
  } while (0)

int
test_lnbeta(void)
{
  int s = 0;

  /* B1: 0 < a < 1, b < 10. */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.25, 0.75, &r),
          1.4913034761293729, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.75, 1.5, &r),
          -0.042373001096346445, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.99999999999999989, 1.0000000000000002, &r),
          -1.110223024625156e-16, TEST_TOL2, GSL_SUCCESS);

  /* B2: 0 < a < 1, b >= 10. */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.25, 10.0, &r),
          0.72182633330350077, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.75, 100.0, &r),
          -3.2496599722271737, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.007085510039214777, 4.3737812239094046e+295, &r),
          0.12227713696720754, TEST_TOL2, GSL_SUCCESS);

  /* B3: 1 <= a <= 2, b <= 2. */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.0, 1.0, &r),
          0.0, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.25, 1.75, &r),
          -0.87582013800224401, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (2.0, 2.0, &r),
          -1.791759469228055, TEST_TOL2, GSL_SUCCESS);

  /* B4: 1 <= a <= 2, 2 < b < 10. */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.0, 3.0, &r),
          -1.0986122886681098, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.5, 5.0, &r),
          -2.605290467144441, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (2.0, 9.5, &r),
          -4.6026670557699729, TEST_TOL2, GSL_SUCCESS);

  /* B5: 1 <= a <= 2, b >= 10. */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.0, 10.0, &r),
          -2.3025850929940459, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.5, 100.0, &r),
          -7.0322750633365985, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (2.0, 10000000000.0, &r),
          -46.051701859980916, TEST_TOL2, GSL_SUCCESS);

  /* B6: 2 < a < 10, b < 10. */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (2.5, 3.5, &r),
          -3.3018352699620528, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (5.0, 7.0, &r),
          -7.7450028035158391, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (9.0, 9.5, &r),
          -12.6493794533343, TEST_TOL2, GSL_SUCCESS);

  /* B7: 2 < a < 10, 10 <= b <= 1000. */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (2.5, 10.0, &r),
          -5.647837161382057, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (5.0, 100.0, &r),
          -19.946329573136683, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (9.0, 1000.0, &r),
          -51.601093037913266, TEST_TOL2, GSL_SUCCESS);

  /* B8: 2 < a < 10, b > 1000. */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (2.5, 1001.0, &r),
          -16.989075958351606, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (5.0, 10000000000.0, &r),
          -111.95120082035434, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (9.0, 1e+100, &r),
          -2061.7219807918959, TEST_TOL2, GSL_SUCCESS);

  /* B9: a >= 10. */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (10.0, 10.0, &r),
          -13.736229227036555, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (20.0, 100.0, &r),
          -54.550806681721241, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (10000000000.0, 1e+20, &r),
          -240258509310.49857, TEST_TOL2, GSL_SUCCESS);

  /*
   * Near the log(B) = 0 contour (min(a,b) < 1): relative accuracy is kept by
   * the double-double refinement.  Reference values from 300-bit arithmetic.
   */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.3, 38.926161425305125, &r),
          2.0532668681600889154e-17, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.3, 38.9261614, &r),
          1.9678350069260651621e-10, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.0703125, 1e15, &r),
          1.8964563891003876547e-01, TEST_TOL2, GSL_SUCCESS);

  /* local expansion about (1,1) */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.97, 1.0296, &r),
          1.8616358920484127504e-03, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.98, 1.0204081632653061, &r),
          2.6321159339894338624e-04, TEST_TOL2, GSL_SUCCESS);

  /*
   * Negative arguments, by reflection to positive arguments; the last
   * group of points is ill-conditioned for the log Gamma formula
   * (e.g. B(-100.1, 99.3), formerly a failing case of gsl_sf_beta).
   */
  TEST_LNBETA_SGN(s, -0.5, 0.75, 1.807705502178633106e-01, -1.0, TEST_TOL2);
  TEST_LNBETA_SGN(s, -3.5, 4.25, 6.0216929102603405764e-01, 1.0, TEST_TOL3);
  TEST_LNBETA_SGN(s, -0.25, 0.125, 1.4436934211998120839, 1.0, TEST_TOL2);
  TEST_LNBETA_SGN(s, -2.5, -1.25, 2.6274559871856658577, -1.0, TEST_TOL2);
  TEST_LNBETA_SGN(s, -10.3, -0.6, 1.7893765605701333318, -1.0, TEST_TOL2);
  TEST_LNBETA_SGN(s, -1e-5, 2.0, 1.1512935465020228672e+01, -1.0, TEST_TOL2);
  TEST_LNBETA_SGN(s, -100.1, 99.3, -7.7119997994469479387, 1.0, TEST_TOL2);
  TEST_LNBETA_SGN(s, -150.5, 3.25, -1.5681708371596409987e+01, -1.0, TEST_TOL2);

  /* negative arguments near |B| = 1: refined in double-double */
  TEST_LNBETA_SGN(s, -0.0077497, 0.00781, -3.6274355689883591215e-03, -1.0, TEST_TOL1);
  TEST_LNBETA_SGN(s, -3.975, 4.056, -5.9763551730297757358e-02, 1.0, TEST_TOL1);
  TEST_LNBETA_SGN(s, -38.109, 0.094646, -4.3683482011070221539e-02, 1.0, TEST_TOL1);
  TEST_LNBETA_SGN(s, -71.645, -0.3698, -2.9885687940504820207e-02, 1.0, TEST_TOL1);

  /* B(a,b) < 0: log B undefined */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (-0.5, 0.75, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);

  /* poles of B, and zeros (a+b a non-positive integer) */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.0, 1.0, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (-2.0, 3.5, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (-1.5, -0.5, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);


  /* cases carried over from the original specfunc test suite */
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.0e-8, 1.0e-8, &r),  19.113827924512310617 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.0e-8, 0.01, &r),  18.420681743788563403 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.0e-8, 1.0, &r),  18.420680743952365472 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.0e-8, 10.0, &r),  18.420680715662683009 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.0e-8, 1000.0, &r),  18.420680669107656949 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.1, 0.1, &r), 2.9813614810376273949 , TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.1, 1.0, &r),  2.3025850929940456840 , TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.1, 100.0, &r),  1.7926462324527931217 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (0.1, 1000, &r),  1.5619821298353164928 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.0, 1.00025, &r),  -0.0002499687552073570, TEST_TOL4, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.0, 1.01, &r),  -0.009950330853168082848 , TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (1.0, 1000.0, &r),  -6.907755278982137052 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (100.0, 100.0, &r),  -139.66525908670663927 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (100.0, 1000.0, &r),  -336.4348576477366051 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_lnbeta_e, (100.0, 1.0e+8, &r),  -1482.9339185256447309 , TEST_TOL0, GSL_SUCCESS);
  return s;
}
