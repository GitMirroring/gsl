/* specfunc/exp/test_exp_mult.c
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
#include <gsl/gsl_sf_exp.h>

#include "../test_sf.h"
#include "test_all.h"

int
test_exp_mult(void)
{
  gsl_sf_result r;
  gsl_sf_result_e10 re;
  double x;
  int s = 0, sa;

  /* cases carried over from the original specfunc test suite */
  x = 0.8*GSL_LOG_DBL_MAX;
  TEST_SF(s, gsl_sf_exp_mult_e, (-10.0,  1.0e-06, &r), 1.0e-06*exp(-10.0), TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (-10.0,  2.0, &r),     2.0*exp(-10.0),     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (-10.0, -2.0, &r),    -2.0*exp(-10.0),     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, ( 10.0,  1.0e-06, &r), 1.0e-06*exp( 10.0), TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, ( 10.0, -2.0, &r),    -2.0*exp( 10.0),     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, 1.00001, &r),      1.00001*exp(x),     TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, 1.000001, &r),     1.000001*exp(x),    TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, 1.000000001, &r),  1.000000001*exp(x), TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, 100.0, &r),        100.0*exp(x),       TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, 1.0e+20, &r),      1.0e+20*exp(x),     TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_e, (x, exp(-x)*exp(M_LN2), &r),  2.0, TEST_TOL4, GSL_SUCCESS );

  TEST_SF(s, gsl_sf_exp_mult_err_e, (-10.0, TEST_SQRT_TOL0, 2.0, TEST_SQRT_TOL0, &r), 2.0*exp(-10.0), TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_exp_mult_err_e, (x, TEST_SQRT_TOL0*x, exp(-x)*exp(M_LN2), TEST_SQRT_TOL0*exp(-x)*exp(M_LN2), &r),  2.0, TEST_SQRT_TOL0, GSL_SUCCESS );

  sa = 0;
  sa += gsl_sf_exp_mult_e10_e(1.0, 1.0, &re);
  sa += ( test_sf_frac_diff(re.val, M_E ) > TEST_TOL0 );
  sa += ( re.err > TEST_TOL2 );
  sa += ( re.e10 != 0 );
  gsl_test(sa, "gsl_sf_exp_mult_e10_e(1.0, 1.0, &re)");
  TEST_SF_E10(s,  gsl_sf_exp_mult_e10_e, (1.0, 1.0, &re), M_E, 0, TEST_TOL0, GSL_SUCCESS);

  sa = 0;
  sa += gsl_sf_exp_mult_e10_e(1000.0, 1.0e+200, &re);
  sa += ( test_sf_frac_diff(re.val, 1.970071114017046993888879352) > TEST_TOL3 );
  sa += ( re.err > 1.0e-11 );
  sa += ( re.e10 != 634 );
  gsl_test(sa, "gsl_sf_exp_mult_e10_e(1000.0, 1.0e+200, &re)");
  TEST_SF_E10(s,  gsl_sf_exp_mult_e10_e, (1000.0, 1.0e+200, &re), 1.970071114017046993888879352, 634, TEST_TOL3, GSL_SUCCESS);

  sa = 0;
  sa += gsl_sf_exp_mult_err_e10_e(1.0, TEST_TOL0, 1.0, TEST_TOL0, &re);
  sa += ( test_sf_frac_diff(re.val, M_E ) > TEST_TOL0 );
  sa += ( re.err > TEST_TOL2 );
  sa += ( re.e10 != 0 );
  gsl_test(sa, "gsl_sf_exp_mult_err_e10_e(1.0, TEST_TOL0, 1.0, TEST_TOL0, &re)");
  TEST_SF_E10(s,  gsl_sf_exp_mult_err_e10_e, (1.0, TEST_TOL0, 1.0, TEST_TOL0, &re), M_E, 0, TEST_TOL0, GSL_SUCCESS);

  sa = 0;
  sa += gsl_sf_exp_mult_err_e10_e(1000.0, 1.0e-12, 1.0e+200, 1.0e+190, &re);
  sa += ( test_sf_frac_diff(re.val, 1.9700711140165661 ) > TEST_TOL3 );
  sa += ( re.err > 1.0e-09 );
  sa += ( re.e10 != 634 );
  gsl_test(sa, "gsl_sf_exp_mult_err_e10_e(1.0e-12, 1.0e+200, &re)");
  TEST_SF_E10(s,  gsl_sf_exp_mult_err_e10_e, (1000.0, 1.0e-12, 1.0e+200, 1.0e+190, &re), 1.9700711140165661,634, TEST_TOL3, GSL_SUCCESS);

  /* Test cases from Szymon Jaroszewicz */
  TEST_SF_E10(s,  gsl_sf_exp_mult_e10_e, (10000.0, 1.0, &re), 8.806818225662921587261496007, 4342, TEST_TOL5, GSL_SUCCESS);
  TEST_SF_E10(s,  gsl_sf_exp_mult_e10_e, (100.0, 1.0, &re), 2.688117141816135448412625551e43, 0, TEST_TOL2, GSL_SUCCESS);

  /*
   * large |x| and |log|y||, where forming exp(x + log|y|) loses
   * eps |log|y|| unaccounted: full accuracy, and the reported error
   * covers the actual one
   */
  TEST_SF_RESULT(s, gsl_sf_exp_mult_e, (-500.0, 1.0e+200, &r), 7.12457640674128531591e-18, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_exp_mult_e, (500.0, 1.0e-200, &r), 140359221785283738.5616, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_exp_mult_e, (-700.0, 1.0e+300, &r), 0.00009859676543759771374385, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_exp_mult_e, (700.0, 1.0e-300, &r), 10142.32054735004534871, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_exp_mult_e, (-1.0, GSL_DBL_MAX, &r), 6.613343458508869992351e+307, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_exp_mult_e, (709.0, 1.5, &r), 1.232761119233245828386e+308, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_exp_mult_err_e, (-500.0, 0.0, 1.0e+200, 0.0, &r), 7.12457640674128531591e-18, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RESULT(s, gsl_sf_exp_mult_err_e, (700.0, 0.0, 1.0e-300, 0.0, &r), 10142.32054735004534871, TEST_TOL0, GSL_SUCCESS);

  /* range limits */
  TEST_SF_RETURN(s, gsl_sf_exp_mult_e, (710.0, 1.0, &r), GSL_EOVRFLW);
  TEST_SF_RETURN(s, gsl_sf_exp_mult_e, (-744.0, 1.0e+10, &r), GSL_EUNDRFLW);
  TEST_SF_RETURN(s, gsl_sf_exp_mult_e, (1.0e+300, 1.0e-300, &r), GSL_EOVRFLW);
  TEST_SF_RETURN(s, gsl_sf_exp_mult_e, (-1.0e+300, 1.0e+300, &r), GSL_EUNDRFLW);

  TEST_SF_E10(s, gsl_sf_exp_mult_e10_e, (1000.0, 1.0e+200, &re), 1.970071114017046934261, 634, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_E10(s, gsl_sf_exp_mult_e10_e, (10000.0, 1.0, &re), 8.806818225662921587261, 4342, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_E10(s, gsl_sf_exp_mult_e10_e, (-10000.0, 1.0e-300, &re), 1.135483865314736126995, -4643, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_E10(s, gsl_sf_exp_mult_e10_e, (1.0e+09, 1.0, &re), 8.002981770660972533042, 434294481, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_E10(s, gsl_sf_exp_mult_e10_e, (-5000.0, 3.0, &re), 1.010908244492675254335, -2171, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_E10(s, gsl_sf_exp_mult_err_e10_e, (-5000.0, 0.0, 3.0, 0.0, &re), 1.010908244492675254335, -2171, TEST_TOL0, GSL_SUCCESS);

  return s;
}
