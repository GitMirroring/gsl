#include <config.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

#include "../test_sf_new.h"

int
test_lngamma(void)
{
  int s = 0;

  /* Near x = 1 rational region: 0.8 <= x <= 1.23395 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (0.8, &r),
          0.152059678399837588778292602291, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (0.9, &r),
          0.0663762397347429711887167398671, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1.1, &r),
          -0.0498724412598397241482898083159, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1.23395, &r),
          -0.09446606702381254116508568555, TEST_TOL1, GSL_SUCCESS);

  /* Middle rational region: 1.23395 < x <= 1.67895 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (1.25, &r),
          -0.0982718364218131614638538026966, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1.5, &r),
          -0.120782237635245222345518445782, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1.67895, &r),
          -0.1000209420168271613054818956, TEST_TOL1, GSL_SUCCESS);

  /* Near x = 2 rational region: 1.67895 < x < 2.3 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (1.75, &r),
          -0.0844011210204855559577860341314, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (2.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (2.25, &r),
          0.124871714892396594302441287613, TEST_TOL1, GSL_SUCCESS);

  /* Recurrence region: 2.3 <= x < 8 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (2.3, &r),
          0.154189454959630581089917911489, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (2.5, &r),
          0.284682870472919159632494669683, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (2.67895, &r),
          0.418147656004422021021414196095, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (3.0, &r),
          0.693147180559945309417232121458, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (3.3, &r),
          0.987098577894734587878679288615, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (3.67895, &r),
          1.4035725827122472653951624762, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (4.5, &r),
          2.45373657084244222050414250344, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (5.5, &r),
          3.95781396761871629387740085582, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (6.5, &r),
          5.66256205985714152852211231233, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (7.5, &r),
          7.53436423675873295515836763244, TEST_TOL2, GSL_SUCCESS);

  /* Stirling U4 region: x >= 8 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (8.0, &r),
          8.52516136106541430016553103635, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (10.0, &r),
          12.8018274800814696112077178746, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (12.0, &r),
          17.5023078458738858392876529072, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (32.0, &r),
          78.0922235533153106314168080587, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (64.0, &r),
          201.009316399281526679282039157, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (100.0, &r),
          359.13420536957539877604401046, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1000.0, &r),
          5905.22042320918121182607691236, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1000000.0, &r),
          12815504.569147611659976971785, TEST_TOL2, GSL_SUCCESS);

  return s;
}
