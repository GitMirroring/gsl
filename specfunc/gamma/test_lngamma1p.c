/* specfunc/gamma/test_lngamma1p.c
 *
 * Unit tests for gsl_sf_lngamma1p_e().
 */

#include <config.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

#include "../test_sf_new.h"

int
test_lngamma1p(void)
{
  int s = 0;

  /* Domain errors. */
  TEST_SF(s, gsl_sf_lngamma1p_e, (-1.0, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  TEST_SF(s, gsl_sf_lngamma1p_e, (-2.0, &r),
          GSL_NAN, GSL_NAN, GSL_EDOM);

  /* Outer negative recurrence region: -1 < x < -0.5. */
  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.9999999999999999, &r),
          36.736800569677101335, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.999999, &r),
          13.815509980720675988, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.99, &r),
          4.5994798780420208294, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.9, &r),
          2.2527126517342061913, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.75, &r),
          1.2880225246980774574, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.5000000000000001, &r),
          0.5723649429247003051, TEST_TOL1, GSL_SUCCESS);

  /* Central negative recurrence region: -0.5 <= x < -0.1. */
  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.5, &r),
          0.57236494292470008, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.25, &r),
          0.20328095143129535, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.125, &r),
          0.085858707225334327, TEST_TOL1, GSL_SUCCESS);

  /* [5/5] Pade region: |x| <= 0.1. */
  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.1, &r),
          0.066376239734742964, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.01, &r),
          0.0058548067647097761, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (-0.001, &r),
          0.0005780385328913798, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (0.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (0.001, &r),
          -0.0005763935982833696, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (0.01, &r),
          -0.0056903079460696449, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (0.1, &r),
          -0.049872441259839723, TEST_TOL1, GSL_SUCCESS);

  /* Lower rational region: 0.1 < x <= 0.5. */
  TEST_SF(s, gsl_sf_lngamma1p_e, (0.125, &r),
          -0.060023184126039582, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (0.25, &r),
          -0.098271836421813169, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (0.5, &r),
          -0.12078223763524522, TEST_TOL1, GSL_SUCCESS);

  /* Upper rational region: 0.5 < x <= 1. */
  TEST_SF(s, gsl_sf_lngamma1p_e, (0.75, &r),
          -0.084401121020485553, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (1.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  /* Positive outer region: x > 1. */
  TEST_SF(s, gsl_sf_lngamma1p_e, (1.0000000000000002, &r),
          9.387698065543119941e-17, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (1.25, &r),
          0.12487171489239659430, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (1.5, &r),
          0.28468287047291915963, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (2.0, &r),
          0.69314718055994530942, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (5.0, &r),
          4.7874917427820459942, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (10.0, &r),
          15.104412573075515295, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (100.0, &r),
          363.73937555556349014, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_lngamma1p_e, (1.0e6, &r),
          12815518.384658169624, TEST_TOL1, GSL_SUCCESS);

  return s;
}
