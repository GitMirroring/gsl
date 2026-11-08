/* rstat/test_wrstat.c
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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_wrstat.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static void
test_wrstat_basic(const double tol, const size_t n, const double weights[], const double data[],
                  const char * desc)
{
  gsl_wrstat_workspace *wrstat_workspace_p = gsl_wrstat_alloc();
  const double expected_wmean = gsl_stats_wmean(weights, 1, data, 1, n);
  double expected_rms = 0.0;
  double wmean, var, sd, rms;
  double sumW = 0.0;
  size_t i, num;
  int status;

  /* compute expected rms */
  for (i = 0; i < n; ++i)
    {
      expected_rms += weights[i] * data[i] * data[i];
      sumW += weights[i];
    }

  expected_rms = sqrt(expected_rms / sumW);

  /* add data to wrstat workspace */
  for (i = 0; i < n; ++i)
    gsl_wrstat_add(weights[i], data[i], wrstat_workspace_p);

  wmean    = gsl_wrstat_mean(wrstat_workspace_p);
  rms      = gsl_wrstat_rms(wrstat_workspace_p);
  num      = gsl_wrstat_n(wrstat_workspace_p);

  gsl_test_int(num, n, "%s n n=%zu", desc, n);
  gsl_test_rel(wmean, expected_wmean, tol, "%s wmean n=%zu", desc, n);
  gsl_test_rel(rms, expected_rms, tol, "%s wrms n=%zu", desc, n);
  
  {
    double sumW2 = 0.0, sumWx2 = 0.0;

    for (i = 0; i < n; ++i)
      {
        sumW2 += weights[i] * weights[i];
        sumWx2 += weights[i] * data[i] * data[i];
      }

    gsl_test_rel(gsl_wrstat_norm(wrstat_workspace_p), sqrt(sumWx2), tol, "%s wnorm n=%zu", desc, n);
    gsl_test_rel(gsl_wrstat_neff(wrstat_workspace_p), sumW * sumW / sumW2, tol, "%s neff n=%zu", desc, n);

    if (n > 1)
      {
        const double expected_var = gsl_stats_wvariance(weights, 1, data, 1, n);
        const double expected_sd = gsl_stats_wsd(weights, 1, data, 1, n);
        const double expected_skew = gsl_stats_wskew(weights, 1, data, 1, n);
        const double expected_kurtosis = gsl_stats_wkurtosis(weights, 1, data, 1, n);

        var      = gsl_wrstat_variance(wrstat_workspace_p);
        sd       = gsl_wrstat_sd(wrstat_workspace_p);

        gsl_test_rel(var, expected_var, tol, "%s wvariance n=%zu", desc, n);
        gsl_test_rel(sd, expected_sd, tol, "%s wstddev n=%zu", desc, n);
        gsl_test_rel(gsl_wrstat_sd_mean(wrstat_workspace_p), expected_sd * sqrt(sumW2) / sumW, tol,
                     "%s wstddev_mean n=%zu", desc, n);
        gsl_test_rel(gsl_wrstat_skew(wrstat_workspace_p), expected_skew, tol, "%s wskew n=%zu", desc, n);
        gsl_test_rel(gsl_wrstat_kurtosis(wrstat_workspace_p), expected_kurtosis, tol,
                     "%s wkurtosis n=%zu", desc, n);
      }
  }

  status = gsl_wrstat_reset(wrstat_workspace_p);
  gsl_test_int(status, GSL_SUCCESS, "%s: wrstat returned success", desc);
  num = gsl_wrstat_n(wrstat_workspace_p);

  gsl_test_int(num, 0, "n n=%zu" , n);

  gsl_wrstat_free(wrstat_workspace_p);
}

/* weight schemes for weighted quantile tests */
enum
{
  WEIGHT_UNIFORM = 0,  /* U(0.5, 1.5) */
  WEIGHT_LOGNORMAL,    /* exp(N(0,1)) */
  WEIGHT_HUBER,        /* Huber weights min(1, c / |x - center|), correlated with x */
  WEIGHT_WIDE,         /* 10^U(-3,3), six decades */
  WEIGHT_NUM
};

static const char * weight_names[WEIGHT_NUM] = { "uniform", "lognormal", "huber", "wide" };

static void
weight_data(const int scheme, const size_t n, const double x[], double wt[], gsl_rng * r)
{
  size_t i;

  for (i = 0; i < n; ++i)
    {
      if (scheme == WEIGHT_UNIFORM)
        wt[i] = 0.5 + gsl_rng_uniform(r);
      else if (scheme == WEIGHT_LOGNORMAL)
        wt[i] = exp(gsl_ran_gaussian(r, 1.0));
      else if (scheme == WEIGHT_HUBER)
        wt[i] = GSL_MIN(1.0, 3.0 / (fabs(x[i] - 1.0) + 1.0e-12));
      else
        wt[i] = pow(10.0, 6.0 * gsl_rng_uniform(r) - 3.0);
    }
}

/*
brute_wselect()
  Reference lo(t) = min { x_i : C_i >= t } and hi(t) = min { x_i : C_i > t }
for sorted x with weights wt, skipping zero weights
*/
static void
brute_wselect(const double t, const double x[], const double wt[], const size_t n,
              double * lo, double * hi)
{
  double cum = 0.0;
  int found_lo = 0, found_hi = 0;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      if (wt[i] == 0.0)
        continue;

      cum += wt[i];

      if (!found_lo && cum >= t)
        {
          *lo = x[i];
          found_lo = 1;
        }

      if (!found_hi && cum > t)
        {
          *hi = x[i];
          found_hi = 1;
        }

      if (!found_lo)
        *lo = x[i];

      if (!found_hi)
        *hi = x[i];
    }
}

/*
brute_wmedian_mad()
  Reference weighted median and MAD0; sorts x with wt, and stores the sorted
absolute deviations from the median in dev, with weights dwt
*/
static void
brute_wmedian_mad(double x[], double wt[], const size_t n, double dev[], double dwt[],
                  double * median, double * mad0)
{
  double W = 0.0, lo, hi;
  size_t i;

  gsl_sort2(x, 1, wt, 1, n);

  for (i = 0; i < n; ++i)
    W += wt[i];

  brute_wselect(0.5 * W, x, wt, n, &lo, &hi);
  *median = 0.5 * (lo + hi);

  for (i = 0; i < n; ++i)
    {
      dev[i] = fabs(x[i] - *median);
      dwt[i] = wt[i];
    }

  gsl_sort2(dev, 1, dwt, 1, n);

  brute_wselect(0.5 * W, dev, dwt, n, &lo, &hi);
  *mad0 = 0.5 * (lo + hi);
}

/*
test_wrstat_exact()
  Results must equal the reference for at most k data with positive
weight; dyadic weights make the cumulative sums exact, and about 10%
zero weights test that these are ignored
*/
static void
test_wrstat_exact(const size_t k, gsl_rng * r)
{
  gsl_wrstat_workspace *w = gsl_wrstat_alloc2(k);
  const size_t nmax = 2 * k + 100;
  double *x = malloc(nmax * sizeof(double));
  double *wt = malloc(nmax * sizeof(double));
  double *dev = malloc(nmax * sizeof(double));
  double *dwt = malloc(nmax * sizeof(double));
  const double p[] = { 0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 };
  size_t n, i;

  for (n = 1; n <= k; ++n)
    {
      double median, mad0, W = 0.0, lo, hi;
      size_t npos = 0, ntot = 0;

      gsl_wrstat_reset(w);

      while (npos < n && ntot < nmax)
        {
          const double u = gsl_rng_uniform(r);

          x[ntot] = floor(10.0 * gsl_rng_uniform(r)) - 3.0;
          wt[ntot] = (u < 0.1) ? 0.0 : 0.25 * (1.0 + floor(8.0 * gsl_rng_uniform(r)));

          gsl_wrstat_add(wt[ntot], x[ntot], w);

          W += wt[ntot];
          npos += (wt[ntot] > 0.0);
          ++ntot;
        }

      brute_wmedian_mad(x, wt, ntot, dev, dwt, &median, &mad0);

      gsl_test_rel(gsl_wrstat_median(w), median, 1.0e-15, "wrstat exact k=%zu n=%zu median", k, n);
      gsl_test_rel(gsl_wrstat_mad0(w), mad0, 1.0e-15, "wrstat exact k=%zu n=%zu mad0", k, n);
      gsl_test_rel(gsl_wrstat_mad(w), 1.482602218505602 * mad0, 1.0e-15,
                   "wrstat exact k=%zu n=%zu mad", k, n);

      for (i = 0; i < sizeof(p) / sizeof(double); ++i)
        {
          brute_wselect(p[i] * W, x, wt, ntot, &lo, &hi);
          gsl_test_rel(gsl_wrstat_quantile(p[i], w), 0.5 * (lo + hi), 1.0e-15,
                       "wrstat exact k=%zu n=%zu quantile p=%g", k, n, p[i]);
        }
    }

  gsl_wrstat_free(w);
  free(x);
  free(wt);
  free(dev);
  free(dwt);
}

/* with unit weights, wrstat must reproduce rstat exactly, also after many compactions */
static void
test_wrstat_unit(const size_t k, const int order, gsl_rng * r)
{
  const size_t n = 100000;
  gsl_rstat_workspace *uw = gsl_rstat_alloc2(k);
  gsl_wrstat_workspace *w = gsl_wrstat_alloc2(k);
  double *x = malloc(n * sizeof(double));
  const char *oname = order_names[order];
  size_t i;

  ordered_data(order, n, x, r);

  for (i = 0; i < n; ++i)
    {
      gsl_rstat_add(x[i], uw);
      gsl_wrstat_add(1.0, x[i], w);

      if ((i + 1) % 9973 == 0 || i + 1 == k)
        {
          gsl_test_rel(gsl_wrstat_median(w), gsl_rstat_median(uw), 0.0,
                       "wrstat unit k=%zu %s n=%zu median", k, oname, i + 1);
          gsl_test_rel(gsl_wrstat_mad0(w), gsl_rstat_mad0(uw), 0.0,
                       "wrstat unit k=%zu %s n=%zu mad0", k, oname, i + 1);
        }
    }

  gsl_rstat_free(uw);
  gsl_wrstat_free(w);
  free(x);
}

/*
queries between additions convert the weights of level 0 between raw and
cumulative form; with dyadic weights this is exact, so the results must
equal those without intermediate queries
*/
static void
test_wrstat_interleaved(const size_t k, const int order, gsl_rng * r)
{
  const size_t n = 50000;
  gsl_wrstat_workspace *w1 = gsl_wrstat_alloc2(k);
  gsl_wrstat_workspace *w2 = gsl_wrstat_alloc2(k);
  double *x = malloc(n * sizeof(double));
  const char *oname = order_names[order];
  size_t i;

  ordered_data(order, n, x, r);

  for (i = 0; i < n; ++i)
    {
      const double wt = 0.25 * (1.0 + floor(8.0 * gsl_rng_uniform(r)));

      gsl_wrstat_add(wt, x[i], w1);
      gsl_wrstat_add(wt, x[i], w2);

      if (i % 101 == 0)
        {
          gsl_wrstat_median(w2);
          gsl_wrstat_mad0(w2);
        }
    }

  gsl_test_rel(gsl_wrstat_median(w2), gsl_wrstat_median(w1), 0.0,
               "wrstat interleaved k=%zu %s median", k, oname);
  gsl_test_rel(gsl_wrstat_mad0(w2), gsl_wrstat_mad0(w1), 0.0,
               "wrstat interleaved k=%zu %s mad0", k, oname);
  gsl_test_rel(gsl_wrstat_quantile(0.1, w2), gsl_wrstat_quantile(0.1, w1), 0.0,
               "wrstat interleaved k=%zu %s quantile", k, oname);

  gsl_wrstat_free(w1);
  gsl_wrstat_free(w2);
  free(x);
}

/* fraction of the total weight of data in x[0..n-1] which are <= v */
static double
test_wcdf(const double v, const double x[], const double wt[], const size_t n)
{
  double below = 0.0, W = 0.0;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      W += wt[i];
      if (x[i] <= v)
        below += wt[i];
    }

  return below / W;
}

/* the weighted CDF must be exact for at most k data with positive (dyadic) weight */
static void
test_wrstat_cdf_exact(const size_t k, gsl_rng * r)
{
  gsl_wrstat_workspace *w = gsl_wrstat_alloc2(k);
  double *x = malloc(k * sizeof(double));
  double *wt = malloc(k * sizeof(double));
  size_t n, i;

  for (n = 1; n <= k; ++n)
    {
      double v;

      gsl_wrstat_reset(w);

      for (i = 0; i < n; ++i)
        {
          x[i] = floor(10.0 * gsl_rng_uniform(r)) - 3.0;
          wt[i] = 0.25 * (1.0 + floor(8.0 * gsl_rng_uniform(r)));
          gsl_wrstat_add(wt[i], x[i], w);
        }

      for (v = -4.5; v <= 7.5; v += 0.5)
        gsl_test_rel(gsl_wrstat_cdf(v, w), test_wcdf(v, x, wt, n), 1.0e-15,
                     "wrstat cdf exact k=%zu n=%zu x=%g", k, n, v);
    }

  gsl_wrstat_free(w);
  free(x);
  free(wt);
}

/*
test_wrstat_merge_moments()
  Split weighted data into nparts parts of random sizes (part 0 empty, and
some zero weights), merge them, and compare with the whole dataset
*/
static void
test_wrstat_merge_moments(const size_t n, const size_t nparts, gsl_rng * r)
{
  gsl_wrstat_workspace **ws = malloc(nparts * sizeof(gsl_wrstat_workspace *));
  double *x = malloc(n * sizeof(double));
  double *wt = malloc(n * sizeof(double));
  const double tol = 1.0e-10;
  size_t i, j;

  for (j = 0; j < nparts; ++j)
    ws[j] = gsl_wrstat_alloc();

  for (i = 0; i < n; ++i)
    {
      x[i] = 5.0 + gsl_ran_exponential(r, 2.0);
      wt[i] = (gsl_rng_uniform(r) < 0.05) ? 0.0 : exp(gsl_ran_gaussian(r, 1.0));

      j = 1 + (size_t) (gsl_rng_uniform(r) * gsl_rng_uniform(r) * (nparts - 1));
      gsl_wrstat_add(wt[i], x[i], ws[j]);
    }

  for (j = 1; j < nparts; ++j)
    gsl_wrstat_merge(ws[j], ws[0]);

  gsl_test_int(gsl_wrstat_n(ws[0]), n, "wrstat merge n=%zu parts=%zu count", n, nparts);
  gsl_test_rel(gsl_wrstat_min(ws[0]), gsl_stats_min(x, 1, n), 0.0, "wrstat merge n=%zu parts=%zu min", n, nparts);
  gsl_test_rel(gsl_wrstat_max(ws[0]), gsl_stats_max(x, 1, n), 0.0, "wrstat merge n=%zu parts=%zu max", n, nparts);
  gsl_test_rel(gsl_wrstat_mean(ws[0]), gsl_stats_wmean(wt, 1, x, 1, n), tol,
               "wrstat merge n=%zu parts=%zu mean", n, nparts);
  gsl_test_rel(gsl_wrstat_variance(ws[0]), gsl_stats_wvariance(wt, 1, x, 1, n), tol,
               "wrstat merge n=%zu parts=%zu variance", n, nparts);
  gsl_test_rel(gsl_wrstat_skew(ws[0]), gsl_stats_wskew(wt, 1, x, 1, n), tol,
               "wrstat merge n=%zu parts=%zu skew", n, nparts);
  gsl_test_rel(gsl_wrstat_kurtosis(ws[0]), gsl_stats_wkurtosis(wt, 1, x, 1, n), tol,
               "wrstat merge n=%zu parts=%zu kurtosis", n, nparts);

  for (j = 0; j < nparts; ++j)
    gsl_wrstat_free(ws[j]);

  free(ws);
  free(x);
  free(wt);
}

/*
test_wrstat_merge_unit()
  Merging with unit weights must reproduce gsl_rstat_merge() exactly, and
for at most k data in total the merged results must be exact
*/
static void
test_wrstat_merge_unit(const size_t k, const int order, const size_t n, const size_t nparts, gsl_rng * r)
{
  gsl_rstat_workspace **us = malloc(nparts * sizeof(gsl_rstat_workspace *));
  gsl_wrstat_workspace **ws = malloc(nparts * sizeof(gsl_wrstat_workspace *));
  double *x = malloc(n * sizeof(double));
  const char *oname = order_names[order];
  size_t i, j;

  for (j = 0; j < nparts; ++j)
    {
      us[j] = gsl_rstat_alloc2(k);
      ws[j] = gsl_wrstat_alloc2(k);
    }

  ordered_data(order, n, x, r);

  for (i = 0; i < n; ++i)
    {
      gsl_rstat_add(x[i], us[i * nparts / n]);
      gsl_wrstat_add(1.0, x[i], ws[i * nparts / n]);
    }

  for (j = 1; j < nparts; ++j)
    {
      gsl_rstat_merge(us[j], us[0]);
      gsl_wrstat_merge(ws[j], ws[0]);
    }

  gsl_test_rel(gsl_wrstat_median(ws[0]), gsl_rstat_median(us[0]), 0.0,
               "wrstat merge unit k=%zu %s n=%zu median", k, oname, n);
  gsl_test_rel(gsl_wrstat_mad0(ws[0]), gsl_rstat_mad0(us[0]), 0.0,
               "wrstat merge unit k=%zu %s n=%zu mad0", k, oname, n);
  gsl_test_rel(gsl_wrstat_cdf(x[n / 3], ws[0]), gsl_rstat_cdf(x[n / 3], us[0]), 0.0,
               "wrstat merge unit k=%zu %s n=%zu cdf", k, oname, n);

  if (n <= k)
    {
      gsl_sort(x, 1, n);
      gsl_test_rel(gsl_wrstat_median(ws[0]), gsl_stats_median_from_sorted_data(x, 1, n), 1.0e-15,
                   "wrstat merge unit k=%zu %s n=%zu exact median", k, oname, n);
    }

  for (j = 0; j < nparts; ++j)
    {
      gsl_rstat_free(us[j]);
      gsl_wrstat_free(ws[j]);
    }

  free(us);
  free(ws);
  free(x);
}

/* merge weighted chunks of a data stream and check weighted rank errors against tol */
static void
test_wrstat_merge_accuracy(const size_t k, const int order, const int scheme, const size_t n,
                           const size_t nparts, const double tol, gsl_rng * r)
{
  gsl_wrstat_workspace **ws = malloc(nparts * sizeof(gsl_wrstat_workspace *));
  double *x = malloc(n * sizeof(double));
  double *wt = malloc(n * sizeof(double));
  double *dev = malloc(n * sizeof(double));
  double *dwt = malloc(n * sizeof(double));
  const char *oname = order_names[order];
  const char *wname = weight_names[scheme];
  double median, mad0, err;
  size_t i, j;

  for (j = 0; j < nparts; ++j)
    ws[j] = gsl_wrstat_alloc2(k);

  ordered_data(order, n, x, r);
  weight_data(scheme, n, x, wt, r);

  for (i = 0; i < n; ++i)
    gsl_wrstat_add(wt[i], x[i], ws[i * nparts / n]);

  for (j = 1; j < nparts; ++j)
    gsl_wrstat_merge(ws[j], ws[0]);

  brute_wmedian_mad(x, wt, n, dev, dwt, &median, &mad0);

  err = test_rank(gsl_wrstat_median(ws[0]), x, wt, n) - 0.5;
  gsl_test(fabs(err) > tol, "wrstat merge accuracy k=%zu %s %s median rank error %.3e",
           k, oname, wname, err);

  err = test_rank(gsl_wrstat_mad0(ws[0]), dev, dwt, n) - 0.5;
  gsl_test(fabs(err) > tol, "wrstat merge accuracy k=%zu %s %s mad rank error %.3e",
           k, oname, wname, err);

  err = gsl_wrstat_cdf(x[n / 4], ws[0]) - test_wcdf(x[n / 4], x, wt, n);
  gsl_test(fabs(err) > tol, "wrstat merge accuracy k=%zu %s %s cdf error %.3e",
           k, oname, wname, err);

  for (j = 0; j < nparts; ++j)
    gsl_wrstat_free(ws[j]);

  free(ws);
  free(x);
  free(wt);
  free(dev);
  free(dwt);
}

/* statistics other than n, min and max of a workspace without data of positive weight are 0 */
static void
test_wrstat_empty_check(const gsl_wrstat_workspace * w, const char * desc)
{
  gsl_test_rel(gsl_wrstat_mean(w), 0.0, 0.0, "wrstat empty %s mean", desc);
  gsl_test_rel(gsl_wrstat_variance(w), 0.0, 0.0, "wrstat empty %s variance", desc);
  gsl_test_rel(gsl_wrstat_sd(w), 0.0, 0.0, "wrstat empty %s sd", desc);
  gsl_test_rel(gsl_wrstat_rms(w), 0.0, 0.0, "wrstat empty %s rms", desc);
  gsl_test_rel(gsl_wrstat_norm(w), 0.0, 0.0, "wrstat empty %s norm", desc);
  gsl_test_rel(gsl_wrstat_neff(w), 0.0, 0.0, "wrstat empty %s neff", desc);
  gsl_test_rel(gsl_wrstat_sd_mean(w), 0.0, 0.0, "wrstat empty %s sd_mean", desc);
  gsl_test_rel(gsl_wrstat_skew(w), 0.0, 0.0, "wrstat empty %s skew", desc);
  gsl_test_rel(gsl_wrstat_kurtosis(w), 0.0, 0.0, "wrstat empty %s kurtosis", desc);
  gsl_test_rel(gsl_wrstat_median(w), 0.0, 0.0, "wrstat empty %s median", desc);
  gsl_test_rel(gsl_wrstat_quantile(0.0, w), 0.0, 0.0, "wrstat empty %s quantile p=0", desc);
  gsl_test_rel(gsl_wrstat_quantile(0.5, w), 0.0, 0.0, "wrstat empty %s quantile p=0.5", desc);
  gsl_test_rel(gsl_wrstat_quantile(1.0, w), 0.0, 0.0, "wrstat empty %s quantile p=1", desc);
  gsl_test_rel(gsl_wrstat_mad0(w), 0.0, 0.0, "wrstat empty %s mad0", desc);
  gsl_test_rel(gsl_wrstat_mad(w), 0.0, 0.0, "wrstat empty %s mad", desc);
  gsl_test_rel(gsl_wrstat_cdf(1.0, w), 0.0, 0.0, "wrstat empty %s cdf", desc);
}

/* fresh and reset workspaces, data with zero weight only, a single positive weight, empty merges */
static void
test_wrstat_empty(void)
{
  gsl_wrstat_workspace *w = gsl_wrstat_alloc2(8);
  gsl_wrstat_workspace *w2 = gsl_wrstat_alloc2(8);
  size_t i;

  test_wrstat_empty_check(w, "fresh");
  gsl_test_int(gsl_wrstat_n(w), 0, "wrstat empty fresh n");
  gsl_test_rel(gsl_wrstat_min(w), 0.0, 0.0, "wrstat empty fresh min");
  gsl_test_rel(gsl_wrstat_max(w), 0.0, 0.0, "wrstat empty fresh max");

  for (i = 0; i < 100; ++i)
    gsl_wrstat_add(1.0 + i, (double) i, w);

  gsl_wrstat_reset(w);
  test_wrstat_empty_check(w, "reset");
  gsl_test_int(gsl_wrstat_n(w), 0, "wrstat empty reset n");

  gsl_test_int(gsl_wrstat_merge(w2, w), GSL_SUCCESS, "wrstat empty merge status");
  test_wrstat_empty_check(w, "merge");

  /* zero weights only: counted in n, min and max, but no other statistic */
  for (i = 0; i < 20; ++i)
    gsl_wrstat_add(0.0, (double) i - 5.0, w);

  test_wrstat_empty_check(w, "zero weights");
  gsl_test_int(gsl_wrstat_n(w), 20, "wrstat empty zero weights n");
  gsl_test_rel(gsl_wrstat_min(w), -5.0, 0.0, "wrstat empty zero weights min");
  gsl_test_rel(gsl_wrstat_max(w), 14.0, 0.0, "wrstat empty zero weights max");

  /* a single datum of positive weight among zero weights */
  gsl_wrstat_add(2.5, 3.0, w);
  gsl_test_rel(gsl_wrstat_mean(w), 3.0, 0.0, "wrstat single weight mean");
  gsl_test_rel(gsl_wrstat_variance(w), 0.0, 0.0, "wrstat single weight variance");
  gsl_test_rel(gsl_wrstat_rms(w), 3.0, 0.0, "wrstat single weight rms");
  gsl_test_rel(gsl_wrstat_neff(w), 1.0, 0.0, "wrstat single weight neff");
  gsl_test_rel(gsl_wrstat_median(w), 3.0, 0.0, "wrstat single weight median");
  gsl_test_rel(gsl_wrstat_mad0(w), 0.0, 0.0, "wrstat single weight mad0");
  gsl_test_rel(gsl_wrstat_cdf(3.0, w), 1.0, 0.0, "wrstat single weight cdf");

  gsl_wrstat_free(w);
  gsl_wrstat_free(w2);
}

/*
test_wrstat_small()
  For every n = 1, ..., nmax, with k = 8, check properties which hold whether
or not the estimates are exact, for weighted data, and that unit weights
reproduce gsl_rstat exactly
*/
static void
test_wrstat_small(const size_t k, const size_t nmax, const int scheme, gsl_rng * r)
{
  const double p[] = { 0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 };
  const size_t np = sizeof(p) / sizeof(double);
  const char *wname = weight_names[scheme];
  gsl_wrstat_workspace *w = gsl_wrstat_alloc2(k);
  gsl_wrstat_workspace *wu = gsl_wrstat_alloc2(k);
  gsl_rstat_workspace *u = gsl_rstat_alloc2(k);
  double xmin = 0.0, xmax = 0.0;
  size_t n, i;

  for (n = 1; n <= nmax; ++n)
    {
      const double v = gsl_ran_gaussian(r, 1.0);
      double wt, q, qprev, median, mad0;

      weight_data(scheme, 1, &v, &wt, r);

      xmin = (n == 1 || v < xmin) ? v : xmin;
      xmax = (n == 1 || v > xmax) ? v : xmax;

      gsl_wrstat_add(wt, v, w);
      gsl_wrstat_add(1.0, v, wu);
      gsl_rstat_add(v, u);

      qprev = xmin;
      for (i = 0; i < np; ++i)
        {
          q = gsl_wrstat_quantile(p[i], w);
          gsl_test(q < qprev || q > xmax, "wrstat small %s n=%zu quantile p=%g ordered", wname, n, p[i]);
          qprev = q;
        }

      median = gsl_wrstat_median(w);
      gsl_test(median < xmin || median > xmax, "wrstat small %s n=%zu median in range", wname, n);

      mad0 = gsl_wrstat_mad0(w);
      gsl_test(mad0 < 0.0 || mad0 > xmax - xmin, "wrstat small %s n=%zu mad0 in range", wname, n);

      gsl_test_rel(gsl_wrstat_cdf(xmax, w), 1.0, 0.0, "wrstat small %s n=%zu cdf at max", wname, n);
      gsl_test_rel(gsl_wrstat_cdf(xmin - 1.0, w), 0.0, 0.0, "wrstat small %s n=%zu cdf below min", wname, n);
      gsl_test_rel(gsl_wrstat_cdf(gsl_wrstat_quantile(1.0, w), w), 1.0, 0.0,
                   "wrstat small %s n=%zu cdf at quantile p=1", wname, n);

      if (n <= k)
        {
          /* all data retained: the extreme quantiles are the exact extremes */
          gsl_test_rel(gsl_wrstat_quantile(0.0, w), xmin, 0.0, "wrstat small %s n=%zu quantile p=0", wname, n);
          gsl_test_rel(gsl_wrstat_quantile(1.0, w), xmax, 0.0, "wrstat small %s n=%zu quantile p=1", wname, n);
        }

      gsl_test_rel(gsl_wrstat_median(wu), gsl_rstat_median(u), 0.0, "wrstat small unit n=%zu median", n);
      gsl_test_rel(gsl_wrstat_mad0(wu), gsl_rstat_mad0(u), 0.0, "wrstat small unit n=%zu mad0", n);
      gsl_test_rel(gsl_wrstat_cdf(v, wu), gsl_rstat_cdf(v, u), 0.0, "wrstat small unit n=%zu cdf", n);
    }

  gsl_wrstat_free(w);
  gsl_wrstat_free(wu);
  gsl_rstat_free(u);
}

/* check normalized weighted rank errors of median, MAD and quantiles against tol */
static void
test_wrstat_accuracy(const size_t k, const int order, const int scheme, const size_t n,
                     const double tol, gsl_rng * r)
{
  gsl_wrstat_workspace *w = gsl_wrstat_alloc2(k);
  double *x = malloc(n * sizeof(double));
  double *wt = malloc(n * sizeof(double));
  double *dev = malloc(n * sizeof(double));
  double *dwt = malloc(n * sizeof(double));
  const double p[] = { 0.01, 0.1, 0.25, 0.75, 0.9, 0.99 };
  const char *oname = order_names[order];
  const char *wname = weight_names[scheme];
  double median, mad0, err;
  size_t i;

  ordered_data(order, n, x, r);
  weight_data(scheme, n, x, wt, r);

  for (i = 0; i < n; ++i)
    gsl_wrstat_add(wt[i], x[i], w);

  brute_wmedian_mad(x, wt, n, dev, dwt, &median, &mad0);

  err = test_rank(gsl_wrstat_median(w), x, wt, n) - 0.5;
  gsl_test(fabs(err) > tol, "wrstat accuracy k=%zu %s %s median rank error %.3e",
           k, oname, wname, err);

  err = test_rank(gsl_wrstat_mad0(w), dev, dwt, n) - 0.5;
  gsl_test(fabs(err) > tol, "wrstat accuracy k=%zu %s %s mad rank error %.3e",
           k, oname, wname, err);

  for (i = 0; i < sizeof(p) / sizeof(double); ++i)
    {
      err = test_rank(gsl_wrstat_quantile(p[i], w), x, wt, n) - p[i];
      gsl_test(fabs(err) > tol, "wrstat accuracy k=%zu %s %s quantile p=%g rank error %.3e",
               k, oname, wname, p[i], err);
    }

  gsl_wrstat_free(w);
  free(x);
  free(wt);
  free(dev);
  free(dwt);
}

/* data with zero, infinite or NaN weights, and NaN data, are not stored; invalid arguments */
static void
test_wrstat_errors(void)
{
  gsl_error_handler_t *old = gsl_set_error_handler_off();
  gsl_wrstat_workspace *w;

  w = gsl_wrstat_alloc2(7);
  gsl_test(w != NULL, "wrstat alloc2 k < 8");

  w = gsl_wrstat_alloc2(8);
  gsl_test_rel(gsl_wrstat_median(w), 0.0, 0.0, "wrstat empty median");

  gsl_wrstat_add(0.0, 5.0, w);
  gsl_test_rel(gsl_wrstat_median(w), 0.0, 0.0, "wrstat zero weight median");
  gsl_test_int(gsl_wrstat_n(w), 1, "wrstat zero weight count");

  gsl_wrstat_add(2.0, 1.0, w);
  gsl_test_int(gsl_wrstat_add(-1.0, 1.0, w), GSL_EDOM, "wrstat negative weight");
  gsl_wrstat_add(GSL_NAN, 7.0, w);
  gsl_wrstat_add(GSL_POSINF, 7.0, w);
  gsl_wrstat_add(1.0, GSL_NAN, w);
  gsl_test_rel(gsl_wrstat_median(w), 1.0, 0.0, "wrstat skipped data median");
  gsl_test_rel(gsl_wrstat_mad0(w), 0.0, 0.0, "wrstat skipped data mad0");
  gsl_test(!gsl_isnan(gsl_wrstat_quantile(-0.1, w)), "wrstat quantile p < 0");
  gsl_test(!gsl_isnan(gsl_wrstat_quantile(1.1, w)), "wrstat quantile p > 1");
  gsl_test(!gsl_isnan(gsl_wrstat_cdf(GSL_NAN, w)), "wrstat cdf NaN");

  {
    gsl_wrstat_workspace *w2 = gsl_wrstat_alloc2(9);
    gsl_wrstat_workspace *w3 = gsl_wrstat_alloc2(8);

    gsl_test_int(gsl_wrstat_merge(w, w), GSL_EINVAL, "wrstat merge src == dest");
    gsl_test_int(gsl_wrstat_merge(w2, w), GSL_EBADLEN, "wrstat merge different k");
    gsl_test_int(gsl_wrstat_merge(w3, w), GSL_SUCCESS, "wrstat merge empty status");
    gsl_test_rel(gsl_wrstat_median(w), 1.0, 0.0, "wrstat merge empty median");

    /* merging into an empty workspace copies */
    gsl_test_int(gsl_wrstat_merge(w, w3), GSL_SUCCESS, "wrstat merge into empty status");
    gsl_test_int(gsl_wrstat_n(w3), gsl_wrstat_n(w), "wrstat merge into empty count");
    gsl_test_rel(gsl_wrstat_median(w3), 1.0, 0.0, "wrstat merge into empty median");

    gsl_wrstat_free(w2);
    gsl_wrstat_free(w3);
  }

  gsl_wrstat_free(w);

  gsl_set_error_handler(old);
}

static int
test_wrstat(gsl_rng * r)
{
  const double tol1 = 1.0e-8;

  {
    int order, scheme;

    test_wrstat_exact(8, r);
    test_wrstat_exact(200, r);
    test_wrstat_cdf_exact(8, r);
    test_wrstat_cdf_exact(200, r);
    test_wrstat_errors();
    test_wrstat_empty();

    for (scheme = 0; scheme < WEIGHT_NUM; ++scheme)
      test_wrstat_small(8, 32, scheme, r);

    test_wrstat_merge_moments(100, 3, r);
    test_wrstat_merge_moments(100000, 8, r);

    for (order = 0; order < ORDER_NUM; ++order)
      {
        test_wrstat_unit(64, order, r);
        test_wrstat_interleaved(64, order, r);
        test_wrstat_merge_unit(200, order, 190, 3, r);
        test_wrstat_merge_unit(64, order, 100000, 8, r);
      }

    for (scheme = 0; scheme < WEIGHT_NUM; ++scheme)
      {
        test_wrstat_accuracy(200, ORDER_RANDOM, scheme, 200000, 0.02, r);
        test_wrstat_accuracy(200, ORDER_ORBIT, scheme, 200000, 0.02, r);
        test_wrstat_accuracy(1000, ORDER_SORTED, scheme, 200000, 0.005, r);
        test_wrstat_accuracy(1000, ORDER_ORBIT, scheme, 200000, 0.005, r);
        test_wrstat_merge_accuracy(200, ORDER_ORBIT, scheme, 200000, 8, 0.02, r);
      }
  }

  {
    const size_t N = 100000;
    double *data = malloc(N * sizeof(double));
    double *weights = malloc(N * sizeof(double));
    double data2[5], weights2[5];
    size_t i, j;
    char buf[64];

    /* test1: test on small datasets n <= 5 */
    for (i = 0; i < 100; ++i)
      {
        random_data(-1.0, 1.0, 5, data2, r);
        random_data(0.0, 10.0, 5, weights2, r);

        for (j = 1; j <= 5; ++j)
          {
            sprintf(buf, "wrstat test1[%zu,%zu]", i, j);
            test_wrstat_basic(tol1, j, weights2, data2, buf);
          }
      }

    /* test2: test on large datasets */

    random_data(-1.0, 1.0, N, data, r);
    random_data(0.0, 10.0, N, weights, r);

    for (i = 1; i <= 10; ++i)
      test_wrstat_basic(tol1, i, weights, data, "wrstat test2");

    test_wrstat_basic(tol1, 100, weights, data, "wrstat test2");
    test_wrstat_basic(tol1, 1000, weights, data, "wrstat test2");
    test_wrstat_basic(tol1, 10000, weights, data, "wrstat test2");
    test_wrstat_basic(tol1, N, weights, data, "wrstat test2");

    /* test3: add large constant */

    for (i = 0; i < 5; ++i)
      data2[i] += 1.0e9;

    test_wrstat_basic(1.0e-6, 5, weights2, data2, "wrstat test3");

    free(data);
    free(weights);
  }

  return GSL_SUCCESS;
}
