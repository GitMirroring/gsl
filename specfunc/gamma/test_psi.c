/* specfunc/gamma/test_psi.c
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
#include <gsl/gsl_sf_psi.h>

#include "../test_sf.h"

int
test_psi(void)
{
  int s = 0;

  /*
   * Digamma: x > 0 (series about the zero x0 = 1.4616..., recurrence,
   * asymptotic series).  Reference values from 200-bit arithmetic.
   */
  TEST_SF_RESULT(s, gsl_sf_psi_e, (1e-300, &r), -9.99999999999999974941e+299, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (0.25, &r), -4.22745353337626540809, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (1.0, &r), -5.77215664901532860607e-01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (1.4616321449683622, &r), -9.24126552172942751679e-17, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (2.5, &r), 7.03156640645243187226e-01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (7.5, &r), 1.94675748424608678807, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (10.0, &r), 2.25175258906672110765, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (1e15, &r), 3.45387763949106847603e+01, TEST_TOL1, GSL_SUCCESS);

  /* x < 0: reflection; near a zero of psi (double-double refinement) */
  TEST_SF_RESULT(s, gsl_sf_psi_e, (-0.5, &r), 3.6489973978576520559e-02, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (-2.5, &r), 1.10315664064524318723, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (-3.75, &r), -1.69238859831133034315, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (-100.5, &r), 4.61512460133806411734, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (-1e-5, &r), 9.99994227678856294116e+04, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (-2.6107208684441447, &r), -1.07202759364100020713e-15, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (-10.000000001, &r), 9.9999991961138843275e+08, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (0.0, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (-3.0, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);

  /* trigamma */
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (1e-100, &r), 9.99999999999999960016e+199, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (0.5, &r), 4.93480220054467930942, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (1.0, &r), 1.64493406684822643647, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (3.5, &r), 3.30357756100234864973e-01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (10.0, &r), 1.05166335681685746122e-01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (1e10, &r), 1.00000000005e-10, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (-0.5, &r), 8.93480220054467930942, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (-2.25, &r), 1.93794105118691373626e+01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (-50.000001, &r), 1.00000000505278467017e+12, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (-999.75, &r), 1.97382090521995036942e+01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (1e-160, &r), GSL_POSINF, GSL_POSINF, GSL_EOVRFLW);

  /* Re psi(1+iy): every branch, and its zero y0 = 0.8838... */
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (0.0, &r), -5.77215664901532860607e-01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (0.25, &r), -5.0590581986957357454e-01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (0.5, &r), -3.28886357229459350344e-01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (0.75, &r), -1.13933483172446434034e-01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (0.8838206378941847, &r), 7.76537442298767551232e-19, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (1.0, &r), 9.46503206224769772719e-02, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (1.25, &r), 2.81294825824918783487e-01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (3.0, &r), 1.10798071071015088079, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (4.0, &r), 1.39153628792164620682, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (7.0, &r), 1.94761433458434866918, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (10.0, &r), 2.30341926367141253517, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (1e5, &r), 1.15129254649785617534e+01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (-2.0, &r), 7.14591515373977526657e-01, TEST_TOL1, GSL_SUCCESS);

  /* integer arguments */
  TEST_SF_RESULT(s, gsl_sf_psi_int_e, (1, &r), -5.77215664901532860607e-01, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_int_e, (100, &r), 4.6001618527380874002, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_int_e, (101, &r), 4.6101618527380874002, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_int_e, (1000000, &r), 1.38155100579641907708e+01, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_int_e, (2, &r), 6.44934066848226436472e-01, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_int_e, (101, &r), 9.95016666333357139525e-03, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_int_e, (1000000, &r), 1.00000050000016666667e-06, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_int_e, (0, &r), GSL_NAN, TEST_TOL0, GSL_EDOM);


  /* cases carried over from the original specfunc test suite */
  TEST_SF_RESULT(s, gsl_sf_psi_int_e, (2, &r), 0.42278433509846713939, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_int_e, (3, &r), 0.92278433509846713939, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_int_e, (4, &r), 1.2561176684318004727, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_int_e, (5, &r), 1.5061176684318004727, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_int_e, (110, &r), 4.695928024251535633, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_int_e, (5000, &r), 8.517093188082904107, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (5000.0, &r), 8.517093188082904107, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (5.0, &r), 1.5061176684318004727, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (-10.5, &r),       2.3982391295357816134,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (-1.0e+5-0.5, &r), 11.512935464924395337, 4.0*TEST_TOL4, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_e, (-262144.0-0.5, &r), 12.476653064769611581, 4.0*TEST_TOL4, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (0.8, &r), -0.07088340212750589223, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (5.0, &r),  1.6127848446157465854, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (100.0, &r),  4.605178519404762003, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (2000.0, &r), 7.600902480375416216, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (-0.8, &r), -0.07088340212750589223, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (-1.0, &r),  0.09465032062247697727, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (-5.0, &r),  1.6127848446157465854, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (-100.0, &r),  4.605178519404762003, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1piy_e, (-2000.0, &r), 7.600902480375416216, TEST_TOL0, GSL_SUCCESS);
  /* Additional test values 1-4 computed using gp-pari and
  Abramowitz & Stegun 6.4.6 */
  TEST_SF_RESULT(s, gsl_sf_psi_1_int_e, (1, &r), 1.6449340668482264364,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_int_e, (3, &r), 0.39493406684822643647, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_int_e, (4, &r), 0.28382295573711532536, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_int_e, (1, &r), 1.6449340668482264365,      TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_int_e, (5, &r), 0.22132295573711532536,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_int_e, (100, &r), 0.010050166663333571395,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_int_e, (110, &r), 0.009132356622022545705,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_int_e, (500, &r), 0.0020020013333322666697, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (1.0/32.0, &r), 1025.5728544782377089,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (5.0, &r), 0.22132295573711532536,      TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (100.0, &r), 0.010050166663333571395,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (110.0, &r), 0.009132356622022545705,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (500.0, &r), 0.0020020013333322666697,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (-1.0 - 1.0/128.0, &r), 16386.648472598746587, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (-1.50, &r), 9.3792466449891237539, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (-10.5, &r), 9.7787577398148123845, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (-15.5, &r), 9.8071247184113896201, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (-50.5, &r), 9.8499971860824842274, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_1_e, (-1000.5, &r), 9.8686054001734414233, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (1, 1, &r), 1.6449340668482264364,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (1, 2, &r), 0.64493406684822643647,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (1, 3, &r), 0.39493406684822643647,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (1, 4, &r), 0.28382295573711532536,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (1, 5, &r), 0.22132295573711532536,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (1, 100, &r), 0.010050166663333571395,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (1, 110, &r), 0.009132356622022545705,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (1, 500, &r), 0.0020020013333322666697, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (3, 5.0, &r), 0.021427828192755075022,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (3, 500.0, &r), 1.6048063999872000683e-08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (10, 5.0, &r), -0.08675107579196581317,    TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (10, 50.0, &r), -4.101091112731268288e-12, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (0, -1.5, &r), 0.70315664064524318723,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_psi_n_e, (1, -1.5, &r), 9.3792466449891237539,   TEST_TOL0, GSL_SUCCESS);
  return s;
}
