/* statistics/test_empty_source.c
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

/* Statistics of an empty dataset (n = 0) are NaN, except sums, which
   are 0. The data pointers are NULL, so a function which reads an
   element when n = 0 crashes the test. */

#if defined(BASE_DOUBLE) || defined(BASE_FLOAT) || defined(BASE_LONG_DOUBLE)
#define TEST_EMPTY_FP 1
#endif

#define TEST_NAN(x, desc) \
  gsl_test (!gsl_isnan (x), NAME(gsl_stats) desc " n=0 returns NaN")

void FUNCTION (test, empty) (void);

void
FUNCTION (test, empty) (void)
{
  const BASE * data = NULL;
  const BASE * data2 = NULL;
  BASE * mdata = NULL;
  double * work = NULL;
  BASE * bwork = NULL;
  int * work_int = NULL;
  const BASE one[1] = { 1 };

  TEST_NAN (FUNCTION (gsl_stats, mean) (data, 1, 0), "_mean");
  TEST_NAN (FUNCTION (gsl_stats, variance) (data, 1, 0), "_variance");
  TEST_NAN (FUNCTION (gsl_stats, variance_m) (data, 1, 0, 0.0), "_variance_m");
  TEST_NAN (FUNCTION (gsl_stats, sd) (data, 1, 0), "_sd");
  TEST_NAN (FUNCTION (gsl_stats, sd_m) (data, 1, 0, 0.0), "_sd_m");
  TEST_NAN (FUNCTION (gsl_stats, variance_with_fixed_mean) (data, 1, 0, 0.0),
            "_variance_with_fixed_mean");
  TEST_NAN (FUNCTION (gsl_stats, sd_with_fixed_mean) (data, 1, 0, 0.0),
            "_sd_with_fixed_mean");
  TEST_NAN (FUNCTION (gsl_stats, absdev) (data, 1, 0), "_absdev");
  TEST_NAN (FUNCTION (gsl_stats, absdev_m) (data, 1, 0, 0.0), "_absdev_m");
  TEST_NAN (FUNCTION (gsl_stats, skew) (data, 1, 0), "_skew");
  TEST_NAN (FUNCTION (gsl_stats, skew_m_sd) (data, 1, 0, 0.0, 1.0), "_skew_m_sd");
  TEST_NAN (FUNCTION (gsl_stats, kurtosis) (data, 1, 0), "_kurtosis");
  TEST_NAN (FUNCTION (gsl_stats, kurtosis_m_sd) (data, 1, 0, 0.0, 1.0),
            "_kurtosis_m_sd");
  TEST_NAN (FUNCTION (gsl_stats, lag1_autocorrelation) (data, 1, 0),
            "_lag1_autocorrelation");
  TEST_NAN (FUNCTION (gsl_stats, lag1_autocorrelation_m) (data, 1, 0, 0.0),
            "_lag1_autocorrelation_m");
  TEST_NAN (FUNCTION (gsl_stats, covariance) (data, 1, data2, 1, 0), "_covariance");
  TEST_NAN (FUNCTION (gsl_stats, covariance_m) (data, 1, data2, 1, 0, 0.0, 0.0),
            "_covariance_m");
  TEST_NAN (FUNCTION (gsl_stats, correlation) (data, 1, data2, 1, 0), "_correlation");
  TEST_NAN (FUNCTION (gsl_stats, spearman) (data, 1, data2, 1, 0, work), "_spearman");
  TEST_NAN (FUNCTION (gsl_stats, pvariance) (data, 1, 0, data2, 1, 0), "_pvariance");
  TEST_NAN (FUNCTION (gsl_stats, pvariance) (data, 1, 0, one, 1, 1),
            "_pvariance (one empty sample)");
  TEST_NAN (FUNCTION (gsl_stats, ttest) (data, 1, 0, data2, 1, 0), "_ttest");
  TEST_NAN (FUNCTION (gsl_stats, median_from_sorted_data) (data, 1, 0),
            "_median_from_sorted_data");
  TEST_NAN (FUNCTION (gsl_stats, median) (mdata, 1, 0), "_median");
  TEST_NAN (FUNCTION (gsl_stats, quantile_from_sorted_data) (data, 1, 0, 0.5),
            "_quantile_from_sorted_data");
  TEST_NAN (FUNCTION (gsl_stats, trmean_from_sorted_data) (0.1, data, 1, 0),
            "_trmean_from_sorted_data");
  TEST_NAN (FUNCTION (gsl_stats, gastwirth_from_sorted_data) (data, 1, 0),
            "_gastwirth_from_sorted_data");
  TEST_NAN (FUNCTION (gsl_stats, mad0) (data, 1, 0, work), "_mad0");
  TEST_NAN (FUNCTION (gsl_stats, mad) (data, 1, 0, work), "_mad");
  TEST_NAN (FUNCTION (gsl_stats, Sn_from_sorted_data) (data, 1, 0, bwork),
            "_Sn_from_sorted_data");
  TEST_NAN (FUNCTION (gsl_stats, Qn_from_sorted_data) (data, 1, 0, bwork, work_int),
            "_Qn_from_sorted_data");

  /* the sum of squares of an empty dataset is 0 */
  gsl_test (FUNCTION (gsl_stats, tss) (data, 1, 0) != 0.0,
            NAME(gsl_stats) "_tss n=0 returns 0");
  gsl_test (FUNCTION (gsl_stats, tss_m) (data, 1, 0, 0.0) != 0.0,
            NAME(gsl_stats) "_tss_m n=0 returns 0");

  /* with a sample of size 1, the other sample gives the pooled variance */
  {
    const BASE b[3] = { 1, 2, 4 };
    double pv = FUNCTION (gsl_stats, pvariance) (one, 1, 1, b, 1, 3);
    double v = FUNCTION (gsl_stats, variance) (b, 1, 3);
    gsl_test_rel (pv, v, 4.0 * GSL_DBL_EPSILON,
                  NAME(gsl_stats) "_pvariance (sample of size 1)");
  }

  /* max and min of an empty dataset: NaN for floating point types, 0
     for integer types; the index functions return 0 */
  {
    BASE min, max;
    size_t imin, imax;

    FUNCTION (gsl_stats, minmax) (&min, &max, data, 1, 0);
    FUNCTION (gsl_stats, minmax_index) (&imin, &imax, data, 1, 0);

#ifdef TEST_EMPTY_FP
    TEST_NAN (FUNCTION (gsl_stats, max) (data, 1, 0), "_max");
    TEST_NAN (FUNCTION (gsl_stats, min) (data, 1, 0), "_min");
    gsl_test (!gsl_isnan (min) || !gsl_isnan (max),
              NAME(gsl_stats) "_minmax n=0 returns NaN");
#else
    gsl_test (FUNCTION (gsl_stats, max) (data, 1, 0) != 0,
              NAME(gsl_stats) "_max n=0 returns 0");
    gsl_test (FUNCTION (gsl_stats, min) (data, 1, 0) != 0,
              NAME(gsl_stats) "_min n=0 returns 0");
    gsl_test (min != 0 || max != 0, NAME(gsl_stats) "_minmax n=0 returns 0");
#endif

    gsl_test (FUNCTION (gsl_stats, max_index) (data, 1, 0) != 0,
              NAME(gsl_stats) "_max_index n=0 returns 0");
    gsl_test (FUNCTION (gsl_stats, min_index) (data, 1, 0) != 0,
              NAME(gsl_stats) "_min_index n=0 returns 0");
    gsl_test (imin != 0 || imax != 0,
              NAME(gsl_stats) "_minmax_index n=0 returns 0");
  }

#ifdef TEST_EMPTY_FP
  /* weighted statistics are undefined without positive weights */
  {
    const BASE * w = NULL;
    const BASE w0[2] = { 0, 0 };
    const BASE x[2] = { 1, 2 };

    TEST_NAN (FUNCTION (gsl_stats, wmean) (w, 1, data, 1, 0), "_wmean");
    TEST_NAN (FUNCTION (gsl_stats, wvariance) (w, 1, data, 1, 0), "_wvariance");
    TEST_NAN (FUNCTION (gsl_stats, wsd) (w, 1, data, 1, 0), "_wsd");
    TEST_NAN (FUNCTION (gsl_stats, wvariance_with_fixed_mean) (w, 1, data, 1, 0, 0.0),
              "_wvariance_with_fixed_mean");
    TEST_NAN (FUNCTION (gsl_stats, wabsdev) (w, 1, data, 1, 0), "_wabsdev");
    TEST_NAN (FUNCTION (gsl_stats, wskew) (w, 1, data, 1, 0), "_wskew");
    TEST_NAN (FUNCTION (gsl_stats, wkurtosis) (w, 1, data, 1, 0), "_wkurtosis");
    TEST_NAN (FUNCTION (gsl_stats, wmean) (w0, 1, x, 1, 2), "_wmean (zero weights)");
    TEST_NAN (FUNCTION (gsl_stats, wvariance_with_fixed_mean) (w0, 1, x, 1, 2, 0.0),
              "_wvariance_with_fixed_mean (zero weights)");

    gsl_test (FUNCTION (gsl_stats, wtss) (w, 1, data, 1, 0) != 0.0,
              NAME(gsl_stats) "_wtss n=0 returns 0");
  }
#endif
}

#undef TEST_NAN
#undef TEST_EMPTY_FP
