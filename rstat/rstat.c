/* rstat/rstat.c
 *
 * Copyright (C) 2015, 2026 Patrick Alken
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

/*
 * The median, quantiles and median absolute deviation (MAD) are estimated
 * with the KLL quantile algorithm (Karnin, Lang and Liberty, 2016), in its
 * lazy form (Ivkin et al., 2019); see kll.c for the storage layout.
 *
 * An item retained at level h represents 2^h data. New data are inserted
 * into level 0 at their sorted position, so that all levels are always
 * sorted. When the buffer is full, the lowest level h at or above its
 * capacity is compacted: the first item is kept at level h if the number
 * of items is odd, and of each remaining pair (x_{2i}, x_{2i+1})
 * one member is promoted to level h+1 with twice the weight. The same
 * member (first or second, with probability 1/2 each) is chosen for all
 * pairs of a compaction. The random choice makes the rank error of each
 * compaction zero-mean, so the errors of different compactions average
 * out instead of accumulating (always keeping the first member biases
 * the median by about 16% in rank for k = 200, n = 10^6). The random
 * number generator is seeded by gsl_rstat_reset(), so results are
 * reproducible.
 *
 * The number of retained data is about 3k, independent of n, and the
 * normalized rank error of any quantile is O(1/k) for any order of the
 * input data. Queries are exact on the weighted retained data, so the
 * results are exact as long as n <= k.
 */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rstat.h>

#include "kll.c"

static int rstat_insert(const double x, const size_t h, gsl_rstat_workspace *w);
static int rstat_compress(gsl_rstat_workspace *w);
static void rstat_merge_moments(const gsl_rstat_workspace *src, gsl_rstat_workspace *dest);
static size_t rstat_count(const gsl_rstat_workspace *w);

gsl_rstat_workspace *
gsl_rstat_alloc(void)
{
  return gsl_rstat_alloc2(KLL_DEFAULT_K);
} /* gsl_rstat_alloc() */

/*
gsl_rstat_alloc2()
  Allocate workspace with accuracy parameter k >= 8 for the median,
quantiles and MAD: results are exact for up to k data, and the
normalized rank error is O(1/k) beyond; about 3k data are retained
*/

gsl_rstat_workspace *
gsl_rstat_alloc2(const size_t k)
{
  gsl_rstat_workspace *w;

  if (k < KLL_MIN_K)
    {
      GSL_ERROR_NULL ("k must be at least 8", GSL_EINVAL);
    }

  w = calloc(1, sizeof(gsl_rstat_workspace));

  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->cap = malloc(KLL_MAX_LEVELS * sizeof(size_t));
  w->lev = malloc((KLL_MAX_LEVELS + 1) * sizeof(size_t));

  if (w->cap == 0 || w->lev == 0)
    {
      gsl_rstat_free(w);
      GSL_ERROR_NULL ("failed to allocate space for levels", GSL_ENOMEM);
    }

  w->k = k;
  w->nmax = kll_capacity(k, w->cap);

  w->buf = malloc(w->nmax * sizeof(double));

  if (w->buf == 0)
    {
      gsl_rstat_free(w);
      GSL_ERROR_NULL ("failed to allocate space for quantile data", GSL_ENOMEM);
    }

  gsl_rstat_reset(w);

  return w;
} /* gsl_rstat_alloc2() */

void
gsl_rstat_free(gsl_rstat_workspace *w)
{
  if (w->cap)
    free(w->cap);

  if (w->lev)
    free(w->lev);

  if (w->buf)
    free(w->buf);

  free(w);
} /* gsl_rstat_free() */

size_t
gsl_rstat_n(const gsl_rstat_workspace *w)
{
  return w->n;
} /* gsl_rstat_n() */

/* add a data point to the running totals */
int
gsl_rstat_add(const double x, gsl_rstat_workspace *w)
{
  double delta = x - w->mean;
  double delta_n, delta_nsq, term1, n;

  /* update min and max */
  if (w->n == 0)
    {
      w->min = x;
      w->max = x;
    }
  else
    {
      if (x < w->min)
        w->min = x;
      if (x > w->max)
        w->max = x;
    }

  /* update mean and variance */
  n = (double) ++(w->n);
  delta_n = delta / n;
  delta_nsq = delta_n * delta_n;
  term1 = delta * delta_n * (n - 1.0);
  w->mean += delta_n;
  w->M4 += term1 * delta_nsq * (n * n - 3.0 * n + 3.0) +
           6.0 * delta_nsq * w->M2 - 4.0 * delta_n * w->M3;
  w->M3 += term1 * delta_n * (n - 2.0) - 3.0 * delta_n * w->M2;
  w->M2 += term1;

  /* update quantile estimator; NaN has no rank and is not stored */
  if (!gsl_isnan(x))
    return rstat_insert(x, 0, w);

  return GSL_SUCCESS;
} /* gsl_rstat_add() */

/*
gsl_rstat_merge()
  Merge the data of workspace src into workspace dest, as if the data
added to src had been added to dest; src is unchanged. Both workspaces
must have the same accuracy parameter k.

The moments are combined exactly with the formulas of Pebay et al. (2016),
and the retained data of src are inserted into the same levels of dest,
with the usual compactions, so the merged estimates of the median,
quantiles and MAD have the same error bound as for a single workspace.
*/

int
gsl_rstat_merge(const gsl_rstat_workspace *src, gsl_rstat_workspace *dest)
{
  if (src == dest)
    {
      GSL_ERROR ("src and dest must be different workspaces", GSL_EINVAL);
    }
  else if (src->k != dest->k)
    {
      GSL_ERROR ("src and dest must have the same accuracy parameter k", GSL_EBADLEN);
    }
  else if (src->n == 0)
    {
      return GSL_SUCCESS;
    }
  else
    {
      size_t h, i;

      rstat_merge_moments(src, dest);

      /* insert retained data of src, top level first */
      for (h = src->nlevels; h-- > 0; )
        {
          for (i = src->lev[h]; i < src->lev[h + 1]; ++i)
            {
              int status = rstat_insert(src->buf[i], h, dest);
              if (status)
                return status;
            }
        }

      return GSL_SUCCESS;
    }
} /* gsl_rstat_merge() */

double
gsl_rstat_min(const gsl_rstat_workspace *w)
{
  return w->min;
} /* gsl_rstat_min() */

double
gsl_rstat_max(const gsl_rstat_workspace *w)
{
  return w->max;
} /* gsl_rstat_max() */

double
gsl_rstat_mean(const gsl_rstat_workspace *w)
{
  return w->mean;
} /* gsl_rstat_mean() */

double
gsl_rstat_variance(const gsl_rstat_workspace *w)
{
  if (w->n > 1)
    {
      double n = (double) w->n;
      return (w->M2 / (n - 1.0));
    }
  else
    return 0.0;
} /* gsl_rstat_variance() */

double
gsl_rstat_sd(const gsl_rstat_workspace *w)
{
  double var = gsl_rstat_variance(w);

  return (sqrt(var));
} /* gsl_rstat_sd() */

double
gsl_rstat_rms(const gsl_rstat_workspace *w)
{
  double rms = 0.0;

  if (w->n > 0)
    {
      double mean = gsl_rstat_mean(w);
      double sigma = gsl_rstat_sd(w);
      double n = (double) w->n;
      double a = sqrt((n - 1.0) / n);
      rms = gsl_hypot(mean, a * sigma);
    }

  return rms;
}

double
gsl_rstat_norm(const gsl_rstat_workspace *w)
{
  double rms = gsl_rstat_rms(w);
  double norm = sqrt((double) w->n) * rms;
  return norm;
}

/* standard deviation of the mean: sigma / sqrt(n) */
double
gsl_rstat_sd_mean(const gsl_rstat_workspace *w)
{
  if (w->n > 0)
    {
      double sd = gsl_rstat_sd(w);
      return (sd / sqrt((double) w->n));
    }
  else
    return 0.0;
} /* gsl_rstat_sd_mean() */

double
gsl_rstat_median(const gsl_rstat_workspace *w)
{
  if (w->lev[0] == w->nmax)
    {
      return 0.0;
    }
  else
    {
      double lo, hi;

      kll_select(0.5 * rstat_count(w), w->buf, NULL, w->lev, w->nlevels, &lo, &hi);

      return (lo + hi) / 2.0;
    }
} /* gsl_rstat_median() */

/*
gsl_rstat_quantile()
  Return estimate of p-quantile, p in [0,1], with linear interpolation
between order statistics as in gsl_stats_quantile_from_sorted_data();
p = 0 and p = 1 return the exact minimum and maximum
*/

double
gsl_rstat_quantile(const double p, const gsl_rstat_workspace *w)
{
  if (p < 0.0 || p > 1.0)
    {
      GSL_ERROR_VAL ("p must be in [0,1]", GSL_EDOM, GSL_NAN);
    }

  if (w->lev[0] == w->nmax)
    {
      return 0.0;
    }
  else
    {
      const size_t n = rstat_count(w);
      const double index = p * (n - 1);
      const size_t lhs = (size_t) index;
      const double delta = index - lhs;
      double xmin = 0.0, xmax = 0.0, x1, x2;

      kll_extremes(w->buf, w->lev, w->nlevels, &xmin, &xmax);

      /*
       * the extreme order statistics are known exactly, from w->min and
       * w->max; these are NaN if NaN was added first, in which case the
       * comparisons are false and the retained extremes are used
       */
      if (lhs >= n - 1)
        return GSL_MAX(w->max, xmax);

      /* order statistics x_(lhs) = lo(lhs+1), x_(lhs+1) = hi(lhs+1) */
      kll_select(lhs + 1.0, w->buf, NULL, w->lev, w->nlevels, &x1, &x2);

      if (lhs == 0)
        x1 = GSL_MIN(w->min, xmin);

      if (lhs + 1 == n - 1)
        x2 = GSL_MAX(w->max, xmax);

      return (1.0 - delta) * x1 + delta * x2;
    }
} /* gsl_rstat_quantile() */

/*
gsl_rstat_mad0()
  Return estimate of the median absolute deviation,
MAD0 = median_i { | x_i - median(x) | }, as in gsl_stats_mad0()
*/

double
gsl_rstat_mad0(const gsl_rstat_workspace *w)
{
  if (w->lev[0] == w->nmax)
    {
      return 0.0;
    }
  else
    {
      const double t = 0.5 * rstat_count(w);
      double lo, hi;

      kll_select(t, w->buf, NULL, w->lev, w->nlevels, &lo, &hi);
      kll_mad_select((lo + hi) / 2.0, t, w->buf, NULL, w->lev, w->nlevels, &lo, &hi);

      return (lo + hi) / 2.0;
    }
} /* gsl_rstat_mad0() */

/*
gsl_rstat_mad()
  Return estimate of the scaled median absolute deviation, 1.4826 MAD0,
as in gsl_stats_mad()
*/

double
gsl_rstat_mad(const gsl_rstat_workspace *w)
{
  return KLL_MAD_SCALE * gsl_rstat_mad0(w);
} /* gsl_rstat_mad() */

/*
gsl_rstat_cdf()
  Return estimate of the empirical cumulative distribution function at x,
the fraction of data <= x; exact as long as n <= k
*/

double
gsl_rstat_cdf(const double x, const gsl_rstat_workspace *w)
{
  if (w->lev[0] == w->nmax)
    {
      return 0.0;
    }
  else if (gsl_isnan(x))
    {
      return GSL_NAN;
    }
  else
    {
      double lw[KLL_MAX_LEVELS];

      kll_level_weights(lw);

      return kll_rank(x, w->buf, NULL, lw, w->lev, w->nlevels) / (double) rstat_count(w);
    }
} /* gsl_rstat_cdf() */

double
gsl_rstat_skew(const gsl_rstat_workspace *w)
{
  if (w->n > 0)
    {
      double n = (double) w->n;
      double fac = pow(n - 1.0, 1.5) / n;
      return ((fac * w->M3) / pow(w->M2, 1.5));
    }
  else
    return 0.0;
} /* gsl_rstat_skew() */

double
gsl_rstat_kurtosis(const gsl_rstat_workspace *w)
{
  if (w->n > 0)
    {
      double n = (double) w->n;
      double fac = ((n - 1.0) / n) * (n - 1.0);
      return ((fac * w->M4) / (w->M2 * w->M2) - 3.0);
    }
  else
    return 0.0;
} /* gsl_rstat_kurtosis() */

int
gsl_rstat_reset(gsl_rstat_workspace *w)
{
  w->min = 0.0;
  w->max = 0.0;
  w->mean = 0.0;
  w->M2 = 0.0;
  w->M3 = 0.0;
  w->M4 = 0.0;
  w->n = 0;

  w->nlevels = 1;
  w->lev[0] = w->nmax;
  w->lev[1] = w->nmax;
  w->budget = w->cap[0];
  w->rng = KLL_SEED;

  return GSL_SUCCESS;
} /* gsl_rstat_reset() */

/*
rstat_insert()
  Insert x into level h, which is kept sorted, compacting first if the
buffer is full: the items of levels 0..h-1, and the items <= x of level
h, move one slot toward the free space in front of level 0
*/

static int
rstat_insert(const double x, const size_t h, gsl_rstat_workspace *w)
{
  int status;
  size_t p, i;

  status = kll_grow_levels(h + 1, w->cap, w->lev, &(w->nlevels), &(w->budget));
  if (status)
    return status;

  /* make room if the buffer is full */
  if (w->nmax - w->lev[0] >= w->budget)
    {
      status = rstat_compress(w);
      if (status)
        return status;
    }

  p = kll_bound(w->buf, w->lev[h], w->lev[h + 1], x, 1);

  memmove(w->buf + w->lev[0] - 1, w->buf + w->lev[0], (p - w->lev[0]) * sizeof(double));
  w->buf[p - 1] = x;

  for (i = 0; i <= h; ++i)
    --(w->lev[i]);

  return GSL_SUCCESS;
}

/*
rstat_merge_moments()
  Combine min, max, mean and central moments of src into dest, with
dest as dataset A and src as dataset B (Pebay et al., 2016):

M2 = M2_A + M2_B + n_A n_B d^2 / n
M3 = M3_A + M3_B + n_A n_B (n_A - n_B) d^3 / n^2 + 3 d (n_A M2_B - n_B M2_A) / n
M4 = M4_A + M4_B + n_A n_B (n_A^2 - n_A n_B + n_B^2) d^4 / n^3
     + 6 d^2 (n_A^2 M2_B + n_B^2 M2_A) / n^2 + 4 d (n_A M3_B - n_B M3_A) / n

with n = n_A + n_B and d = mean_B - mean_A
*/

static void
rstat_merge_moments(const gsl_rstat_workspace *src, gsl_rstat_workspace *dest)
{
  if (dest->n == 0)
    {
      dest->min = src->min;
      dest->max = src->max;
      dest->mean = src->mean;
      dest->M2 = src->M2;
      dest->M3 = src->M3;
      dest->M4 = src->M4;
    }
  else
    {
      const double na = (double) dest->n;
      const double nb = (double) src->n;
      const double n = na + nb;
      const double d = src->mean - dest->mean;
      const double d2 = d * d;

      dest->M4 += src->M4 + na * nb * (na * na - na * nb + nb * nb) * d2 * d2 / (n * n * n) +
                  6.0 * d2 * (na * na * src->M2 + nb * nb * dest->M2) / (n * n) +
                  4.0 * d * (na * src->M3 - nb * dest->M3) / n;
      dest->M3 += src->M3 + na * nb * (na - nb) * d2 * d / (n * n) +
                  3.0 * d * (na * src->M2 - nb * dest->M2) / n;
      dest->M2 += src->M2 + na * nb * d2 / n;
      dest->mean += nb * d / n;

      if (src->min < dest->min)
        dest->min = src->min;

      if (src->max > dest->max)
        dest->max = src->max;
    }

  dest->n += src->n;
}

/*
rstat_compress()
  Compact the lowest level at or above its capacity, freeing at least one
slot in the buffer. The odd item (if any) stays at level h, the remaining
items are halved in place into the front of level h, merged in place with
level h+1, and levels 0..h are shifted up to close the gap.
*/

static int
rstat_compress(gsl_rstat_workspace *w)
{
  size_t *lev = w->lev;
  double *buf = w->buf;
  size_t h;
  int status;

  status = kll_compact_level(w->cap, lev, &(w->nlevels), &(w->budget), &h);
  if (status)
    return status;

  {
    const size_t a = lev[h];
    const size_t b = lev[h + 1];
    const size_t e = lev[h + 2];
    const size_t s = b - a;
    const size_t a1 = a + (s & 1);          /* odd item at a stays at level h */
    const size_t half = (s - (s & 1)) / 2;
    const size_t r = (kll_uniform(&(w->rng)) < 0.5) ? 0 : 1;
    size_t i, ia, ib, out;

    /* keep member r of each pair, in buf[a1 .. a1+half-1] */
    for (i = 0; i < half; ++i)
      buf[a1 + i] = buf[a1 + r + 2 * i];

    /* merge buf[a1 .. a1+half-1] with level h+1 = buf[b .. e-1] into buf[a1+half .. e-1] */
    ia = a1;
    ib = b;
    out = a1 + half;

    while (ia < a1 + half && ib < e)
      {
        if (buf[ib] < buf[ia])
          buf[out++] = buf[ib++];
        else
          buf[out++] = buf[ia++];
      }

    while (ia < a1 + half)
      buf[out++] = buf[ia++];

    /* shift levels 0..h-1 and the odd item of level h up by half */
    memmove(buf + lev[0] + half, buf + lev[0], (a1 - lev[0]) * sizeof(double));

    for (i = 0; i <= h; ++i)
      lev[i] += half;

    lev[h + 1] = a1 + half;
  }

  return GSL_SUCCESS;
}

/* number of data represented by the retained data (n, excluding NaN) */
static size_t
rstat_count(const gsl_rstat_workspace *w)
{
  size_t count = 0;
  size_t h;

  for (h = 0; h < w->nlevels; ++h)
    count += (w->lev[h + 1] - w->lev[h]) << h;

  return count;
}
