/* rstat/test_rstat.c
 * 
 * Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2026 Patrick Alken
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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_ieee_utils.h>

static void
test_basic(const size_t n, const double data[], const double tol, const char * desc)
{
  gsl_rstat_workspace *rstat_workspace_p = gsl_rstat_alloc();
  const double expected_mean = gsl_stats_mean(data, 1, n);
  const double expected_skew = gsl_stats_skew(data, 1, n);
  const double expected_kurtosis = gsl_stats_kurtosis(data, 1, n);
  double expected_rms = 0.0;
  double expected_norm;
  double mean, var, sd, sd_mean, rms, skew, kurtosis, norm;
  size_t i, num;
  int status;

  /* compute expected rms */
  for (i = 0; i < n; ++i)
    expected_rms += data[i] * data[i];

  expected_norm = sqrt(expected_rms);
  expected_rms = sqrt(expected_rms / n);

  /* add data to rstat workspace */
  for (i = 0; i < n; ++i)
    gsl_rstat_add(data[i], rstat_workspace_p);

  mean     = gsl_rstat_mean(rstat_workspace_p);
  rms      = gsl_rstat_rms(rstat_workspace_p);
  norm     = gsl_rstat_norm(rstat_workspace_p);
  skew     = gsl_rstat_skew(rstat_workspace_p);
  kurtosis = gsl_rstat_kurtosis(rstat_workspace_p);
  num      = gsl_rstat_n(rstat_workspace_p);

  gsl_test_int(num, n, "%s n n=%zu", desc, n);
  gsl_test_rel(mean, expected_mean, tol, "%s mean n=%zu", desc, n);
  gsl_test_rel(rms, expected_rms, tol, "%s rms n=%zu", desc, n);
  gsl_test_rel(norm, expected_norm, tol, "%s norm n=%zu", desc, n);
  gsl_test_rel(skew, expected_skew, tol, "%s skew n=%zu", desc, n);
  gsl_test_rel(kurtosis, expected_kurtosis, tol, "%s kurtosis n=%zu", desc, n);
  
  if (n > 1)
    {
      const double expected_var = gsl_stats_variance(data, 1, n);
      const double expected_sd = gsl_stats_sd(data, 1, n);
      const double expected_sd_mean = expected_sd / sqrt((double) n);

      var      = gsl_rstat_variance(rstat_workspace_p);
      sd       = gsl_rstat_sd(rstat_workspace_p);
      sd_mean  = gsl_rstat_sd_mean(rstat_workspace_p);

      gsl_test_rel(var, expected_var, tol, "%s variance n=%zu", desc, n);
      gsl_test_rel(sd, expected_sd, tol, "%s stddev n=%zu", desc, n);
      gsl_test_rel(sd_mean, expected_sd_mean, tol, "%s stddev_mean n=%zu", desc, n);
    }

  if (n <= 200)
    {
      /* median, quantiles and MAD are exact for n <= k = 200 (default) */
      const double p[] = { 0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 };
      double * data_copy = malloc(n * sizeof(double));
      double * work = malloc(n * sizeof(double));
      double expected_median, median;

      memcpy(data_copy, data, n * sizeof(double));
      gsl_sort(data_copy, 1, n);
      expected_median = gsl_stats_median_from_sorted_data(data_copy, 1, n);
      median = gsl_rstat_median(rstat_workspace_p);

      gsl_test_rel(median, expected_median, tol, "%s median n=%zu", desc, n);

      for (i = 0; i < sizeof(p) / sizeof(double); ++i)
        {
          double expected = gsl_stats_quantile_from_sorted_data(data_copy, 1, n, p[i]);
          gsl_test_rel(gsl_rstat_quantile(p[i], rstat_workspace_p), expected, tol,
                       "%s quantile n=%zu p=%g", desc, n, p[i]);
        }

      gsl_test_rel(gsl_rstat_mad0(rstat_workspace_p), gsl_stats_mad0(data, 1, n, work), tol,
                   "%s mad0 n=%zu", desc, n);
      gsl_test_rel(gsl_rstat_mad(rstat_workspace_p), gsl_stats_mad(data, 1, n, work), tol,
                   "%s mad n=%zu", desc, n);

      free(data_copy);
      free(work);
    }

  status = gsl_rstat_reset(rstat_workspace_p);
  gsl_test_int(status, GSL_SUCCESS, "%s: rstat returned success", desc);
  num = gsl_rstat_n(rstat_workspace_p);

  gsl_test_int(num, 0, "n n=%zu" , n);

  gsl_rstat_free(rstat_workspace_p);
}

static void
test_quantile(const double p, const double data[], const size_t n,
              const double expected, const double tol, const char *desc)
{
  gsl_rstat_quantile_workspace *w = gsl_rstat_quantile_alloc(p);
  double result;
  size_t i;

  for (i = 0; i < n; ++i)
    gsl_rstat_quantile_add(data[i], w);

  result = gsl_rstat_quantile_get(w);

  if (fabs(expected) < 1.0e-4)
    gsl_test_abs(result, expected, tol, "%s p=%g", desc, p);
  else
    gsl_test_rel(result, expected, tol, "%s p=%g", desc, p);

  gsl_rstat_quantile_free(w);
}

/* results must be exact for n <= k, with many ties */
static void
test_rstat_exact(const size_t k, gsl_rng * r)
{
  gsl_rstat_workspace *w = gsl_rstat_alloc2(k);
  double *x = malloc(k * sizeof(double));
  double *work = malloc(k * sizeof(double));
  const double p[] = { 0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 };
  size_t n, i;

  for (n = 1; n <= k; ++n)
    {
      double expected;

      gsl_rstat_reset(w);

      for (i = 0; i < n; ++i)
        {
          x[i] = floor(10.0 * gsl_rng_uniform(r)) - 3.0;
          gsl_rstat_add(x[i], w);
        }

      expected = gsl_stats_mad0(x, 1, n, work);
      gsl_test_rel(gsl_rstat_mad0(w), expected, 1.0e-15, "exact k=%zu n=%zu mad0", k, n);

      expected = gsl_stats_mad(x, 1, n, work);
      gsl_test_rel(gsl_rstat_mad(w), expected, 1.0e-15, "exact k=%zu n=%zu mad", k, n);

      gsl_sort(x, 1, n);

      expected = gsl_stats_median_from_sorted_data(x, 1, n);
      gsl_test_rel(gsl_rstat_median(w), expected, 1.0e-15, "exact k=%zu n=%zu median", k, n);

      for (i = 0; i < sizeof(p) / sizeof(double); ++i)
        {
          expected = gsl_stats_quantile_from_sorted_data(x, 1, n, p[i]);
          gsl_test_rel(gsl_rstat_quantile(p[i], w), expected, 1.0e-15,
                       "exact k=%zu n=%zu quantile p=%g", k, n, p[i]);
        }
    }

  gsl_rstat_free(w);
  free(x);
  free(work);
}

/*
test_rstat_accuracy()
  Check normalized rank errors of median, MAD and quantiles against tol
for a data stream of the given ordering. The MAD rank error is the rank
of the estimated MAD among the exact absolute deviations from the exact
median, minus 1/2.
*/
static void
test_rstat_accuracy(const size_t k, const int order, const size_t n, const double tol, gsl_rng * r)
{
  gsl_rstat_workspace *w = gsl_rstat_alloc2(k);
  double *x = malloc(n * sizeof(double));
  double *dev = malloc(n * sizeof(double));
  const double p[] = { 0.01, 0.1, 0.25, 0.75, 0.9, 0.99 };
  const char *oname = order_names[order];
  double median, err;
  size_t i;

  ordered_data(order, n, x, r);

  for (i = 0; i < n; ++i)
    gsl_rstat_add(x[i], w);

  gsl_sort(x, 1, n);

  median = gsl_stats_median_from_sorted_data(x, 1, n);
  for (i = 0; i < n; ++i)
    dev[i] = fabs(x[i] - median);
  gsl_sort(dev, 1, n);

  err = test_rank(gsl_rstat_median(w), x, NULL, n) - 0.5;
  gsl_test(fabs(err) > tol, "accuracy k=%zu %s median rank error %.3e", k, oname, err);

  err = test_rank(gsl_rstat_mad0(w), dev, NULL, n) - 0.5;
  gsl_test(fabs(err) > tol, "accuracy k=%zu %s mad rank error %.3e", k, oname, err);

  for (i = 0; i < sizeof(p) / sizeof(double); ++i)
    {
      err = test_rank(gsl_rstat_quantile(p[i], w), x, NULL, n) - p[i];
      gsl_test(fabs(err) > tol, "accuracy k=%zu %s quantile p=%g rank error %.3e",
               k, oname, p[i], err);
    }

  gsl_test_rel(gsl_rstat_quantile(0.0, w), x[0], 0.0, "accuracy k=%zu %s quantile p=0", k, oname);
  gsl_test_rel(gsl_rstat_quantile(1.0, w), x[n - 1], 0.0, "accuracy k=%zu %s quantile p=1", k, oname);

  for (i = 0; i < sizeof(p) / sizeof(double); ++i)
    {
      const double v = x[(size_t) (p[i] * n)];

      err = gsl_rstat_cdf(v, w) - test_cdf(v, x, n);
      gsl_test(fabs(err) > tol, "accuracy k=%zu %s cdf at p=%g error %.3e", k, oname, p[i], err);
    }

  gsl_test_rel(gsl_rstat_cdf(x[0] - 1.0, w), 0.0, 0.0, "accuracy k=%zu %s cdf below min", k, oname);
  gsl_test_rel(gsl_rstat_cdf(x[n - 1], w), 1.0, 0.0, "accuracy k=%zu %s cdf at max", k, oname);

  gsl_rstat_free(w);
  free(x);
  free(dev);
}

/* the CDF must be exact for n <= k, with many ties */
static void
test_rstat_cdf_exact(const size_t k, gsl_rng * r)
{
  gsl_rstat_workspace *w = gsl_rstat_alloc2(k);
  double *x = malloc(k * sizeof(double));
  size_t n, i;

  for (n = 1; n <= k; ++n)
    {
      double v;

      gsl_rstat_reset(w);

      for (i = 0; i < n; ++i)
        {
          x[i] = floor(10.0 * gsl_rng_uniform(r)) - 3.0;
          gsl_rstat_add(x[i], w);
        }

      /* at the data values, between them, and outside their range */
      for (v = -4.5; v <= 7.5; v += 0.5)
        gsl_test_rel(gsl_rstat_cdf(v, w), test_cdf(v, x, n), 1.0e-15,
                     "cdf exact k=%zu n=%zu x=%g", k, n, v);
    }

  gsl_rstat_free(w);
  free(x);
}

/*
test_rstat_merge_moments()
  Split data into nparts parts of random sizes (some empty), accumulate each
part in its own workspace, merge them into one, and compare the moments
with those of the whole dataset
*/
static void
test_rstat_merge_moments(const size_t n, const size_t nparts, gsl_rng * r)
{
  gsl_rstat_workspace **ws = malloc(nparts * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace *all = gsl_rstat_alloc();
  double *x = malloc(n * sizeof(double));
  const double tol = 1.0e-10;
  size_t i, j;

  for (j = 0; j < nparts; ++j)
    ws[j] = gsl_rstat_alloc();

  for (i = 0; i < n; ++i)
    {
      x[i] = 5.0 + gsl_ran_exponential(r, 2.0);

      /* part j chosen at random, with part 0 left empty */
      j = 1 + (size_t) (gsl_rng_uniform(r) * gsl_rng_uniform(r) * (nparts - 1));
      gsl_rstat_add(x[i], ws[j]);
      gsl_rstat_add(x[i], all);
    }

  /* merge parts 1..nparts-1 into part 0 */
  for (j = 1; j < nparts; ++j)
    gsl_rstat_merge(ws[j], ws[0]);

  gsl_test_int(gsl_rstat_n(ws[0]), n, "merge n=%zu parts=%zu count", n, nparts);
  gsl_test_rel(gsl_rstat_min(ws[0]), gsl_stats_min(x, 1, n), 0.0, "merge n=%zu parts=%zu min", n, nparts);
  gsl_test_rel(gsl_rstat_max(ws[0]), gsl_stats_max(x, 1, n), 0.0, "merge n=%zu parts=%zu max", n, nparts);
  gsl_test_rel(gsl_rstat_mean(ws[0]), gsl_stats_mean(x, 1, n), tol, "merge n=%zu parts=%zu mean", n, nparts);
  gsl_test_rel(gsl_rstat_variance(ws[0]), gsl_stats_variance(x, 1, n), tol,
               "merge n=%zu parts=%zu variance", n, nparts);
  gsl_test_rel(gsl_rstat_skew(ws[0]), gsl_stats_skew(x, 1, n), tol, "merge n=%zu parts=%zu skew", n, nparts);
  gsl_test_rel(gsl_rstat_kurtosis(ws[0]), gsl_stats_kurtosis(x, 1, n), tol,
               "merge n=%zu parts=%zu kurtosis", n, nparts);

  /* same retained data size as a single workspace, within the budget */
  gsl_test(ws[0]->nmax - ws[0]->lev[0] > ws[0]->budget, "merge n=%zu parts=%zu budget", n, nparts);

  for (j = 0; j < nparts; ++j)
    gsl_rstat_free(ws[j]);

  free(ws);
  gsl_rstat_free(all);
  free(x);
}

/* merging parts with at most k data in total must give the exact median, quantiles and MAD */
static void
test_rstat_merge_exact(gsl_rng * r)
{
  const size_t k = 200;
  const size_t n = 190;
  gsl_rstat_workspace *w1 = gsl_rstat_alloc2(k);
  gsl_rstat_workspace *w2 = gsl_rstat_alloc2(k);
  gsl_rstat_workspace *w3 = gsl_rstat_alloc2(k);
  double *x = malloc(n * sizeof(double));
  double *work = malloc(n * sizeof(double));
  double median, mad0;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      x[i] = gsl_ran_gaussian(r, 1.0);
      gsl_rstat_add(x[i], (i < 50) ? w1 : (i < 120) ? w2 : w3);
    }

  median = gsl_rstat_median(w2);
  mad0 = gsl_rstat_mad0(w2);

  gsl_rstat_merge(w2, w1);
  gsl_rstat_merge(w3, w1);

  /* src is unchanged */
  gsl_test_int(gsl_rstat_n(w2), 70, "merge exact src count");
  gsl_test_rel(gsl_rstat_median(w2), median, 0.0, "merge exact src median");
  gsl_test_rel(gsl_rstat_mad0(w2), mad0, 0.0, "merge exact src mad0");

  gsl_test_rel(gsl_rstat_mad0(w1), gsl_stats_mad0(x, 1, n, work), 1.0e-15, "merge exact mad0");

  gsl_sort(x, 1, n);
  gsl_test_rel(gsl_rstat_median(w1), gsl_stats_median_from_sorted_data(x, 1, n), 1.0e-15,
               "merge exact median");
  gsl_test_rel(gsl_rstat_quantile(0.1, w1), gsl_stats_quantile_from_sorted_data(x, 1, n, 0.1), 1.0e-15,
               "merge exact quantile");
  gsl_test_rel(gsl_rstat_cdf(x[37], w1), test_cdf(x[37], x, n), 1.0e-15, "merge exact cdf");

  gsl_rstat_free(w1);
  gsl_rstat_free(w2);
  gsl_rstat_free(w3);
  free(x);
  free(work);
}

/*
test_rstat_merge_accuracy()
  Accumulate a data stream of the given ordering in nparts consecutive
chunks, merge the chunks, and check normalized rank errors against tol
*/
static void
test_rstat_merge_accuracy(const size_t k, const int order, const size_t n, const size_t nparts,
                          const double tol, gsl_rng * r)
{
  gsl_rstat_workspace **ws = malloc(nparts * sizeof(gsl_rstat_workspace *));
  double *x = malloc(n * sizeof(double));
  double *dev = malloc(n * sizeof(double));
  const char *oname = order_names[order];
  double median, err;
  size_t i, j;

  for (j = 0; j < nparts; ++j)
    ws[j] = gsl_rstat_alloc2(k);

  ordered_data(order, n, x, r);

  for (i = 0; i < n; ++i)
    gsl_rstat_add(x[i], ws[i * nparts / n]);

  for (j = 1; j < nparts; ++j)
    gsl_rstat_merge(ws[j], ws[0]);

  gsl_sort(x, 1, n);

  median = gsl_stats_median_from_sorted_data(x, 1, n);
  for (i = 0; i < n; ++i)
    dev[i] = fabs(x[i] - median);
  gsl_sort(dev, 1, n);

  err = test_rank(gsl_rstat_median(ws[0]), x, NULL, n) - 0.5;
  gsl_test(fabs(err) > tol, "merge accuracy k=%zu %s parts=%zu median rank error %.3e",
           k, oname, nparts, err);

  err = test_rank(gsl_rstat_mad0(ws[0]), dev, NULL, n) - 0.5;
  gsl_test(fabs(err) > tol, "merge accuracy k=%zu %s parts=%zu mad rank error %.3e",
           k, oname, nparts, err);

  err = gsl_rstat_cdf(x[n / 10], ws[0]) - test_cdf(x[n / 10], x, n);
  gsl_test(fabs(err) > tol, "merge accuracy k=%zu %s parts=%zu cdf error %.3e",
           k, oname, nparts, err);

  for (j = 0; j < nparts; ++j)
    gsl_rstat_free(ws[j]);

  free(ws);
  free(x);
  free(dev);
}

/* all statistics of an empty workspace are 0 */
static void
test_rstat_empty_check(const gsl_rstat_workspace * w, const char * desc)
{
  gsl_test_int(gsl_rstat_n(w), 0, "empty %s n", desc);
  gsl_test_rel(gsl_rstat_min(w), 0.0, 0.0, "empty %s min", desc);
  gsl_test_rel(gsl_rstat_max(w), 0.0, 0.0, "empty %s max", desc);
  gsl_test_rel(gsl_rstat_mean(w), 0.0, 0.0, "empty %s mean", desc);
  gsl_test_rel(gsl_rstat_variance(w), 0.0, 0.0, "empty %s variance", desc);
  gsl_test_rel(gsl_rstat_sd(w), 0.0, 0.0, "empty %s sd", desc);
  gsl_test_rel(gsl_rstat_rms(w), 0.0, 0.0, "empty %s rms", desc);
  gsl_test_rel(gsl_rstat_norm(w), 0.0, 0.0, "empty %s norm", desc);
  gsl_test_rel(gsl_rstat_sd_mean(w), 0.0, 0.0, "empty %s sd_mean", desc);
  gsl_test_rel(gsl_rstat_skew(w), 0.0, 0.0, "empty %s skew", desc);
  gsl_test_rel(gsl_rstat_kurtosis(w), 0.0, 0.0, "empty %s kurtosis", desc);
  gsl_test_rel(gsl_rstat_median(w), 0.0, 0.0, "empty %s median", desc);
  gsl_test_rel(gsl_rstat_quantile(0.0, w), 0.0, 0.0, "empty %s quantile p=0", desc);
  gsl_test_rel(gsl_rstat_quantile(0.5, w), 0.0, 0.0, "empty %s quantile p=0.5", desc);
  gsl_test_rel(gsl_rstat_quantile(1.0, w), 0.0, 0.0, "empty %s quantile p=1", desc);
  gsl_test_rel(gsl_rstat_mad0(w), 0.0, 0.0, "empty %s mad0", desc);
  gsl_test_rel(gsl_rstat_mad(w), 0.0, 0.0, "empty %s mad", desc);
  gsl_test_rel(gsl_rstat_cdf(1.0, w), 0.0, 0.0, "empty %s cdf", desc);
}

/* statistics of a fresh workspace, a reset workspace, and empty merges */
static void
test_rstat_empty(void)
{
  gsl_rstat_workspace *w = gsl_rstat_alloc2(8);
  gsl_rstat_workspace *w2 = gsl_rstat_alloc2(8);
  size_t i;

  test_rstat_empty_check(w, "fresh");

  for (i = 0; i < 100; ++i)
    gsl_rstat_add((double) i, w);

  gsl_rstat_reset(w);
  test_rstat_empty_check(w, "reset");

  gsl_test_int(gsl_rstat_merge(w2, w), GSL_SUCCESS, "empty merge status");
  test_rstat_empty_check(w, "merge");

  gsl_rstat_free(w);
  gsl_rstat_free(w2);
}

/* bug #68283: skewness and kurtosis are undefined (NaN) for data with zero
   variance, including a single value; n = 2 is well defined */
static void
test_rstat_zero_variance(void)
{
  gsl_rstat_workspace *w = gsl_rstat_alloc();
  size_t n;

  for (n = 1; n <= 4; ++n)
    {
      gsl_rstat_add(2.5, w);
      gsl_test_rel(gsl_rstat_variance(w), 0.0, 0.0, "constant n=%zu variance", n);
      gsl_test_rel(gsl_rstat_skew(w), GSL_NAN, 0.0, "constant n=%zu skew", n);
      gsl_test_rel(gsl_rstat_kurtosis(w), GSL_NAN, 0.0, "constant n=%zu kurtosis", n);
    }

  gsl_rstat_reset(w);
  gsl_rstat_add(1.0, w);
  gsl_rstat_add(3.0, w);
  gsl_test_rel(gsl_rstat_skew(w), 0.0, 0.0, "n=2 skew");
  gsl_test_rel(gsl_rstat_kurtosis(w), -2.75, GSL_DBL_EPSILON, "n=2 kurtosis");

  gsl_rstat_free(w);
}

/*
test_rstat_small()
  For every n = 1, ..., nmax, with k = 8 so that compactions start at
n = k + 1, check properties which hold whether or not the estimates are
exact: exact extremes, ordering of quantiles, range of median, MAD and
CDF, and a loose bound on the normalized rank error of the median; the
results must also not depend on whether queries are made in between.

type: 0 = random, 1 = sorted, 2 = reverse sorted, 3 = many ties
*/
static void
test_rstat_small(const size_t k, const size_t nmax, const int type, gsl_rng * r)
{
  const double p[] = { 0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 };
  const size_t np = sizeof(p) / sizeof(double);
  gsl_rstat_workspace *w = gsl_rstat_alloc2(k);
  gsl_rstat_workspace *w0 = gsl_rstat_alloc2(k); /* no intermediate queries */
  double *x = malloc(nmax * sizeof(double));
  double xmin = 0.0, xmax = 0.0;
  size_t n, i;

  for (n = 1; n <= nmax; ++n)
    {
      double v, q, qprev, median, mad0, err;

      if (type == 1)
        v = (double) n;
      else if (type == 2)
        v = -(double) n;
      else if (type == 3)
        v = floor(4.0 * gsl_rng_uniform(r));
      else
        v = gsl_ran_gaussian(r, 1.0);

      x[n - 1] = v;
      xmin = (n == 1 || v < xmin) ? v : xmin;
      xmax = (n == 1 || v > xmax) ? v : xmax;

      gsl_rstat_add(v, w);
      gsl_rstat_add(v, w0);

      gsl_test_rel(gsl_rstat_quantile(0.0, w), xmin, 0.0, "small type=%d n=%zu quantile p=0", type, n);
      gsl_test_rel(gsl_rstat_quantile(1.0, w), xmax, 0.0, "small type=%d n=%zu quantile p=1", type, n);

      qprev = xmin;
      for (i = 0; i < np; ++i)
        {
          q = gsl_rstat_quantile(p[i], w);
          gsl_test(q < qprev, "small type=%d n=%zu quantile p=%g ordered", type, n, p[i]);
          qprev = q;
        }

      median = gsl_rstat_median(w);
      gsl_test(median < xmin || median > xmax, "small type=%d n=%zu median in range", type, n);

      /* distance of 1/2 from the rank interval [fraction < median, fraction <= median] */
      {
        const double le = test_cdf(median, x, n);
        const double lt = 2.0 * test_rank(median, x, NULL, n) - le;

        err = (0.5 < lt) ? lt - 0.5 : (0.5 > le) ? 0.5 - le : 0.0;
        gsl_test(err > 0.25, "small type=%d n=%zu median rank error %.3e", type, n, err);
      }

      mad0 = gsl_rstat_mad0(w);
      gsl_test(mad0 < 0.0 || mad0 > xmax - xmin, "small type=%d n=%zu mad0 in range", type, n);

      gsl_test_rel(gsl_rstat_cdf(xmax, w), 1.0, 0.0, "small type=%d n=%zu cdf at max", type, n);
      gsl_test_rel(gsl_rstat_cdf(xmin - 1.0, w), 0.0, 0.0, "small type=%d n=%zu cdf below min", type, n);
      {
        /* volatile: compare the results rounded to double (x87 returns
           them in extended precision) */
        volatile double cdf_median = gsl_rstat_cdf(median, w);
        volatile double cdf_min = gsl_rstat_cdf(xmin, w);
        gsl_test(cdf_median < cdf_min, "small type=%d n=%zu cdf ordered", type, n);
      }

      gsl_test_rel(gsl_rstat_median(w0), median, 0.0, "small type=%d n=%zu median without queries", type, n);
    }

  gsl_rstat_free(w);
  gsl_rstat_free(w0);
  free(x);
}

/* merges of small workspaces, crossing n = k, must keep the exact extremes and counts */
static void
test_rstat_small_merge(const size_t k, gsl_rng * r)
{
  size_t n1, n2, i;

  for (n1 = 0; n1 <= 2 * k + 1; ++n1)
    {
      for (n2 = 0; n2 <= 2 * k + 1; n2 += 3)
        {
          gsl_rstat_workspace *w1 = gsl_rstat_alloc2(k);
          gsl_rstat_workspace *w2 = gsl_rstat_alloc2(k);
          double xmin = GSL_POSINF, xmax = GSL_NEGINF;

          for (i = 0; i < n1 + n2; ++i)
            {
              const double v = gsl_ran_gaussian(r, 1.0);

              gsl_rstat_add(v, (i < n1) ? w1 : w2);
              xmin = GSL_MIN(xmin, v);
              xmax = GSL_MAX(xmax, v);
            }

          gsl_rstat_merge(w2, w1);

          gsl_test_int(gsl_rstat_n(w1), n1 + n2, "small merge n1=%zu n2=%zu count", n1, n2);

          if (n1 + n2 > 0)
            {
              const double median = gsl_rstat_median(w1);

              gsl_test_rel(gsl_rstat_quantile(0.0, w1), xmin, 0.0, "small merge n1=%zu n2=%zu min", n1, n2);
              gsl_test_rel(gsl_rstat_quantile(1.0, w1), xmax, 0.0, "small merge n1=%zu n2=%zu max", n1, n2);
              gsl_test(median < xmin || median > xmax, "small merge n1=%zu n2=%zu median in range", n1, n2);
              gsl_test_rel(gsl_rstat_cdf(xmax, w1), 1.0, 0.0, "small merge n1=%zu n2=%zu cdf at max", n1, n2);
            }

          gsl_rstat_free(w1);
          gsl_rstat_free(w2);
        }
    }
}

/* reset must restore the initial state; queries between additions must not change results */
static void
test_rstat_reset(gsl_rng * r)
{
  const size_t n = 100000;
  gsl_rstat_workspace *w = gsl_rstat_alloc2(64);
  double *x = malloc(n * sizeof(double));
  double median, mad0;
  size_t i;

  ordered_data(ORDER_RANDOM, n, x, r);

  for (i = 0; i < n; ++i)
    gsl_rstat_add(x[i], w);

  median = gsl_rstat_median(w);
  mad0 = gsl_rstat_mad0(w);

  gsl_rstat_reset(w);
  gsl_test_rel(gsl_rstat_median(w), 0.0, 0.0, "reset empty median");
  gsl_test_rel(gsl_rstat_mad0(w), 0.0, 0.0, "reset empty mad0");
  gsl_test_rel(gsl_rstat_quantile(0.3, w), 0.0, 0.0, "reset empty quantile");

  for (i = 0; i < n; ++i)
    {
      gsl_rstat_add(x[i], w);

      if (i % 1000 == 0)
        gsl_rstat_mad0(w);
    }

  gsl_test_rel(gsl_rstat_median(w), median, 0.0, "reset repeat median");
  gsl_test_rel(gsl_rstat_mad0(w), mad0, 0.0, "reset repeat mad0");

  gsl_rstat_free(w);
  free(x);
}

/* NaN data are not stored in the quantile estimator; invalid arguments */
static void
test_rstat_errors(void)
{
  gsl_error_handler_t *old = gsl_set_error_handler_off();
  gsl_rstat_workspace *w;

  w = gsl_rstat_alloc2(7);
  gsl_test(w != NULL, "alloc2 k < 8");

  w = gsl_rstat_alloc2(8);
  gsl_rstat_add(1.0, w);
  gsl_rstat_add(GSL_NAN, w);
  gsl_rstat_add(3.0, w);
  gsl_test_int(gsl_rstat_n(w), 3, "NaN count");
  gsl_test_rel(gsl_rstat_median(w), 2.0, 0.0, "NaN median");
  gsl_test_rel(gsl_rstat_mad0(w), 1.0, 0.0, "NaN mad0");
  gsl_test(!gsl_isnan(gsl_rstat_quantile(-0.1, w)), "quantile p < 0");
  gsl_test(!gsl_isnan(gsl_rstat_quantile(1.1, w)), "quantile p > 1");
  gsl_test(!gsl_isnan(gsl_rstat_cdf(GSL_NAN, w)), "cdf NaN");

  {
    gsl_rstat_workspace *w2 = gsl_rstat_alloc2(9);
    gsl_rstat_workspace *w3 = gsl_rstat_alloc2(8);

    gsl_test_int(gsl_rstat_merge(w, w), GSL_EINVAL, "merge src == dest");
    gsl_test_int(gsl_rstat_merge(w2, w), GSL_EBADLEN, "merge different k");

    /* merging an empty workspace changes nothing */
    gsl_test_int(gsl_rstat_merge(w3, w), GSL_SUCCESS, "merge empty status");
    gsl_test_int(gsl_rstat_n(w), 3, "merge empty count");
    gsl_test_rel(gsl_rstat_median(w), 2.0, 0.0, "merge empty median");

    gsl_rstat_free(w2);
    gsl_rstat_free(w3);
  }

  gsl_rstat_free(w);

  gsl_set_error_handler(old);
}

static int
test_rstat(gsl_rng * r)
{
  const double tol1 = 1.0e-8;
  const double tol2 = 1.0e-3;

  {
    int order;

    test_rstat_exact(8, r);
    test_rstat_exact(200, r);
    test_rstat_cdf_exact(8, r);
    test_rstat_cdf_exact(200, r);
    test_rstat_empty();
    test_rstat_zero_variance();
    test_rstat_reset(r);

    {
      int type;

      for (type = 0; type < 4; ++type)
        test_rstat_small(8, 32, type, r);

      test_rstat_small_merge(8, r);
    }
    test_rstat_errors();

    test_rstat_merge_moments(100, 3, r);
    test_rstat_merge_moments(100000, 8, r);
    test_rstat_merge_exact(r);

    for (order = 0; order < ORDER_NUM; ++order)
      {
        test_rstat_accuracy(200, order, 200000, 0.02, r);
        test_rstat_accuracy(1000, order, 200000, 0.005, r);
        test_rstat_merge_accuracy(200, order, 200000, 8, 0.02, r);
      }
  }

  {
    const size_t N = 2000000;
    double *data = malloc(N * sizeof(double));
    double data2[5];
    size_t i, j;
    char buf[64];

    /* test1: test on small datasets n <= 5 (median will be exact in this case) */
    for (i = 0; i < 100; ++i)
      {
        random_data(-1.0, 1.0, 5, data2, r);

        for (j = 1; j <= 5; ++j)
          {
            sprintf(buf, "test1 j=%zu", j);
            test_basic(j, data2, tol1, buf);
          }
      }

    /* test2: test on large datasets */

    random_data(-1.0, 1.0, N, data, r);

    for (i = 1; i <= 10; ++i)
      test_basic(i, data, tol1, "test2");

    test_basic(100, data, tol1, "test2");
    test_basic(1000, data, tol1, "test2");
    test_basic(10000, data, tol1, "test2");
    test_basic(50000, data, tol1, "test2");
    test_basic(80000, data, tol1, "test2");
    test_basic(1500000, data, tol1, "test2");
    test_basic(2000000, data, tol1, "test2");

    /* test3: add large constant */

    for (i = 0; i < 5; ++i)
      data2[i] += 1.0e9;

    test_basic(5, data2, 1.0e-6, "test3");

    free(data);
  }

  {
    /* dataset from Jain and Chlamtac paper */
    const size_t n_jain = 20;
    const double data_jain[] = {  0.02,  0.15,  0.74,  3.39,  0.83,
                                  22.37, 10.15, 15.43, 38.62, 15.92,
                                  34.60, 10.28,  1.47,  0.40,  0.05,
                                  11.39,  0.27,  0.42,  0.09, 11.37 };
    double expected_jain = 4.44063435326;
  
    test_quantile(0.5, data_jain, n_jain, expected_jain, tol1, "jain");
  }

  {
    size_t n = 1000000;
    double *data = malloc(n * sizeof(double));
    double *sorted_data = malloc(n * sizeof(double));
    gsl_rstat_workspace *rstat_workspace_p = gsl_rstat_alloc();
    double p;
    size_t i;

    for (i = 0; i < n; ++i)
      {
        data[i] = gsl_ran_gaussian_tail(r, 1.3, 1.0);
        gsl_rstat_add(data[i], rstat_workspace_p);
      }

    memcpy(sorted_data, data, n * sizeof(double));
    gsl_sort(sorted_data, 1, n);

    /* test quantile calculation */
    for (p = 0.1; p <= 0.9; p += 0.1)
      {
        double expected = gsl_stats_quantile_from_sorted_data(sorted_data, 1, n, p);
        test_quantile(p, data, n, expected, tol2, "gauss");
      }

    /* test mean, variance */
    {
      const double expected_mean = gsl_stats_mean(data, 1, n);
      const double expected_var = gsl_stats_variance(data, 1, n);
      const double expected_sd = gsl_stats_sd(data, 1, n);
      const double expected_skew = gsl_stats_skew(data, 1, n);
      const double expected_kurtosis = gsl_stats_kurtosis(data, 1, n);
      const double mean = gsl_rstat_mean(rstat_workspace_p);
      const double var = gsl_rstat_variance(rstat_workspace_p);
      const double sd = gsl_rstat_sd(rstat_workspace_p);
      const double skew = gsl_rstat_skew(rstat_workspace_p);
      const double kurtosis = gsl_rstat_kurtosis(rstat_workspace_p);
      const double median = gsl_rstat_median(rstat_workspace_p);

      gsl_test_rel(mean, expected_mean, tol1, "mean");
      gsl_test_rel(var, expected_var, tol1, "variance");
      gsl_test_rel(sd, expected_sd, tol1, "stddev");
      gsl_test_rel(skew, expected_skew, tol1, "skew");
      gsl_test_rel(kurtosis, expected_kurtosis, tol1, "kurtosis");

      /* normalized rank error of median is O(1/k), k = 200 (default) */
      gsl_test_abs(test_rank(median, sorted_data, NULL, n), 0.5, 0.02, "median rank");
    }

    free(data);
    free(sorted_data);
    gsl_rstat_free(rstat_workspace_p);
  }

  return GSL_SUCCESS;
}
