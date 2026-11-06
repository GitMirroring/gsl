/* rstat/wrstat.c
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

/*
 * This module uses the method described in the following paper,
 *
 * [1] Pebay, Terriberry, Kolla, and Bennett, "Numerically stable,
 *     scalable formulas for parallel and online computation of
 *     higher-order multivariate central moments with arbitrary weights",
 *     Comput. Stat., 31, pg. 1305-1325, 2016.
 *
 * [2] Kish, "Survey Sampling", John Wiley & Sons, New York, 1965
 *     (effective sample size).
 *
 * The weighted median, quantiles and MAD are estimated with a weighted
 * form of the KLL quantile algorithm (see rstat.c and kll.c), in which
 * every retained item carries its weight explicitly. When level h is
 * compacted (all levels are kept sorted), the first item is kept at level h if
 * their number is odd, and the remaining 2m items, with total weight
 * W_h, are reduced by systematic resampling: with U = W_h / m and one
 * uniform random number u, the point (j + u) U, j = 0..m-1, selects the
 * item whose cumulative weight interval contains it, and the selected
 * items are promoted to level h+1 with weight U per point (an item
 * selected by several points is promoted once, with the combined weight).
 *
 * This keeps the expected weighted rank of every value unchanged, and the
 * rank error of a compaction is less than U regardless of how the weights
 * are distributed. For unit weights, U = 2 and point j selects item
 * 2j + [u >= 1/2], so the algorithm reduces exactly to the unweighted one.
 *
 * For fast queries (see kll.c), buf_w holds the cumulative weights within
 * each level, so that the weight of the first j items of a level is a
 * single lookup. New data are inserted into level 0 at their sorted
 * position, updating the cumulative weights of the larger items. A
 * compaction converts the two levels involved to raw weights and back.
 *
 * The weighted quantiles are defined as
 *
 * quantile(p) = (lo(pW) + hi(pW)) / 2
 *
 * with lo, hi as in kll.c and W the total weight; the median is
 * quantile(1/2) and MAD0 is the weighted median of | x_i - median |.
 * For unit weights the median and MAD agree with gsl_stats_median()
 * and gsl_stats_mad0(). Results are exact as long as at most k data with
 * positive weight have been added.
 */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_wrstat.h>

#include "kll.c"

static int wrstat_insert(const double wt, const double x, const size_t h, gsl_wrstat_workspace *w);
static int wrstat_compress(gsl_wrstat_workspace *w);
static void wrstat_merge_moments(const gsl_wrstat_workspace *src, gsl_wrstat_workspace *dest);
static double wrstat_total_weight(const gsl_wrstat_workspace *w);
static void wrstat_to_raw(gsl_wrstat_workspace *w, const size_t h);
static void wrstat_to_cum(gsl_wrstat_workspace *w, const size_t h);

gsl_wrstat_workspace *
gsl_wrstat_alloc(void)
{
  return gsl_wrstat_alloc2(KLL_DEFAULT_K);
}

/*
gsl_wrstat_alloc2()
  Allocate workspace with accuracy parameter k >= 8 for the weighted
median, quantiles and MAD: results are exact for up to k data, and the
normalized rank error is O(1/k) beyond; about 3k data are retained
*/

gsl_wrstat_workspace *
gsl_wrstat_alloc2(const size_t k)
{
  gsl_wrstat_workspace *w;

  if (k < KLL_MIN_K)
    {
      GSL_ERROR_NULL ("k must be at least 8", GSL_EINVAL);
    }

  w = calloc(1, sizeof(gsl_wrstat_workspace));

  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->cap = malloc(KLL_MAX_LEVELS * sizeof(size_t));
  w->lev = malloc((KLL_MAX_LEVELS + 1) * sizeof(size_t));

  if (w->cap == 0 || w->lev == 0)
    {
      gsl_wrstat_free(w);
      GSL_ERROR_NULL ("failed to allocate space for levels", GSL_ENOMEM);
    }

  w->k = k;
  w->nmax = kll_capacity(k, w->cap);

  w->buf_x = malloc(w->nmax * sizeof(double));
  w->buf_w = malloc(w->nmax * sizeof(double));

  if (w->buf_x == 0 || w->buf_w == 0)
    {
      gsl_wrstat_free(w);
      GSL_ERROR_NULL ("failed to allocate space for quantile data", GSL_ENOMEM);
    }

  gsl_wrstat_reset(w);

  return w;
}

void
gsl_wrstat_free(gsl_wrstat_workspace *w)
{
  if (w->cap)
    free(w->cap);

  if (w->lev)
    free(w->lev);

  if (w->buf_x)
    free(w->buf_x);

  if (w->buf_w)
    free(w->buf_w);

  free(w);
}

size_t
gsl_wrstat_n(const gsl_wrstat_workspace *w)
{
  return w->n;
}

/* add a data point to the running totals */
int
gsl_wrstat_add(const double wt, const double x, gsl_wrstat_workspace * w)
{
  if (wt < 0.0)
    {
      GSL_ERROR ("weight must be >= 0", GSL_EDOM);
    }
  else
    {
      double delta = x - w->wmean;

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

      /*
       * update mean and central moments, eq. (3.1)-(3.4) of [1] with
       * dataset B a single point of weight wt:
       *
       * M2 += W_A wt d^2 / W
       * M3 += W_A wt (W_A - wt) d^3 / W^2 - 3 wt d M2_A / W
       * M4 += W_A wt (W_A^2 - W_A wt + wt^2) d^4 / W^3
       *       + 6 wt^2 d^2 M2_A / W^2 - 4 wt d M3_A / W
       *
       * with W = W_A + wt and d = x - wmean_A
       */
      if (wt > 0.0)
        {
          double prevW = w->W;
          double ratio, ratio_a, delta2;

          w->W += wt;
          w->W2 += wt * wt;

          ratio = wt / w->W;
          ratio_a = prevW / w->W;
          delta2 = delta * delta;

          w->M4 += prevW * ratio * (ratio_a * ratio_a - ratio_a * ratio + ratio * ratio) * delta2 * delta2 +
                   6.0 * ratio * ratio * delta2 * w->M2 - 4.0 * ratio * delta * w->M3;
          w->M3 += prevW * ratio * (ratio_a - ratio) * delta2 * delta - 3.0 * ratio * delta * w->M2;
          w->wmean += delta * ratio;
          w->M2 += prevW * ratio * delta2;
        }

      ++(w->n);

      /*
       * update quantile estimator; data with zero, infinite or NaN weight,
       * and NaN data, are not stored
       */
      if (wt > 0.0 && gsl_finite(wt) && !gsl_isnan(x))
        return wrstat_insert(wt, x, 0, w);

      return GSL_SUCCESS;
    }
}

/*
gsl_wrstat_merge()
  Merge the data of workspace src into workspace dest, as if the data
added to src had been added to dest; src is unchanged. Both workspaces
must have the same accuracy parameter k. The weighted moments are
combined exactly with the formulas of [1], and the retained data of src
are inserted into the same levels of dest, as in gsl_rstat_merge().
*/

int
gsl_wrstat_merge(const gsl_wrstat_workspace * src, gsl_wrstat_workspace * dest)
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

      wrstat_merge_moments(src, dest);

      /* insert retained data of src with their weights, top level first */
      for (h = src->nlevels; h-- > 0; )
        {
          for (i = src->lev[h]; i < src->lev[h + 1]; ++i)
            {
              const double wt = (i > src->lev[h]) ? src->buf_w[i] - src->buf_w[i - 1] : src->buf_w[i];
              int status = wrstat_insert(wt, src->buf_x[i], h, dest);
              if (status)
                return status;
            }
        }

      return GSL_SUCCESS;
    }
}

double
gsl_wrstat_min(const gsl_wrstat_workspace * w)
{
  return w->min;
}

double
gsl_wrstat_max(const gsl_wrstat_workspace * w)
{
  return w->max;
}

double
gsl_wrstat_mean(const gsl_wrstat_workspace * w)
{
  return w->wmean;
}

double
gsl_wrstat_variance(const gsl_wrstat_workspace * w)
{
  double denom = w->W * w->W - w->W2;

  /*
   * denom = 2 sum_{i<j} w_i w_j is zero if fewer than two data have
   * positive weight, in which case the variance is 0, as for n = 1
   */
  if (w->n > 1 && denom > 0.0)
    return (w->W * w->M2 / denom);
  else
    return 0.0;
}

double
gsl_wrstat_sd(const gsl_wrstat_workspace * w)
{
  double var = gsl_wrstat_variance(w);
  return (sqrt(var));
}

double
gsl_wrstat_rms(const gsl_wrstat_workspace * w)
{
  double rms = 0.0;

  /* data with zero weight do not contribute */
  if (w->W > 0.0)
    {
      double wmean = gsl_wrstat_mean(w);
      rms = gsl_hypot(wmean, sqrt(w->M2 / w->W));
    }

  return rms;
}

/* weighted Euclidean norm, sqrt(sum w_i x_i^2) */
double
gsl_wrstat_norm(const gsl_wrstat_workspace * w)
{
  return sqrt(w->W) * gsl_wrstat_rms(w);
}

/*
Kish's effective sample size [2], n_eff = (sum w_i)^2 / sum w_i^2, which equals
the number of data with positive weight for equal weights
*/
double
gsl_wrstat_neff(const gsl_wrstat_workspace * w)
{
  if (w->W2 > 0.0)
    return (w->W * w->W / w->W2);
  else
    return 0.0;
}

/* standard deviation of the weighted mean, sd / sqrt(n_eff) */
double
gsl_wrstat_sd_mean(const gsl_wrstat_workspace * w)
{
  if (w->W > 0.0)
    return (gsl_wrstat_sd(w) * sqrt(w->W2) / w->W);
  else
    return 0.0;
}

/* weighted skewness, (M3 / W) / sd^3, as in gsl_stats_wskew() */
double
gsl_wrstat_skew(const gsl_wrstat_workspace * w)
{
  if (w->W > 0.0)
    {
      const double var = gsl_wrstat_variance(w);
      return ((w->M3 / w->W) / pow(var, 1.5));
    }
  else
    return 0.0;
}

/* weighted kurtosis, (M4 / W) / sd^4 - 3, as in gsl_stats_wkurtosis() */
double
gsl_wrstat_kurtosis(const gsl_wrstat_workspace * w)
{
  if (w->W > 0.0)
    {
      const double var = gsl_wrstat_variance(w);
      return ((w->M4 / w->W) / (var * var) - 3.0);
    }
  else
    return 0.0;
}

/*
gsl_wrstat_cdf()
  Return estimate of the weighted empirical cumulative distribution
function at x, the fraction of the total weight of data <= x
*/

double
gsl_wrstat_cdf(const double x, const gsl_wrstat_workspace * w)
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

      return kll_rank(x, w->buf_x, w->buf_w, lw, w->lev, w->nlevels) / wrstat_total_weight(w);
    }
}

/*
gsl_wrstat_quantile()
  Return estimate of weighted p-quantile, p in [0,1],
quantile(p) = (lo(pW) + hi(pW)) / 2; for unit weights this is the
average of the two order statistics bracketing rank pn, which differs
from the interpolation of gsl_rstat_quantile() except for the median
*/

double
gsl_wrstat_quantile(const double p, const gsl_wrstat_workspace * w)
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
      const double W = wrstat_total_weight(w);
      double lo, hi;

      kll_select(p * W, w->buf_x, w->buf_w, w->lev, w->nlevels, &lo, &hi);

      return (lo + hi) / 2.0;
    }
}

double
gsl_wrstat_median(const gsl_wrstat_workspace * w)
{
  return gsl_wrstat_quantile(0.5, w);
}

/* estimate of weighted MAD0, the weighted median of | x_i - median | */
double
gsl_wrstat_mad0(const gsl_wrstat_workspace * w)
{
  if (w->lev[0] == w->nmax)
    {
      return 0.0;
    }
  else
    {
      const double t = 0.5 * wrstat_total_weight(w);
      double lo, hi;

      kll_select(t, w->buf_x, w->buf_w, w->lev, w->nlevels, &lo, &hi);
      kll_mad_select((lo + hi) / 2.0, t, w->buf_x, w->buf_w, w->lev, w->nlevels, &lo, &hi);

      return (lo + hi) / 2.0;
    }
}

/* estimate of scaled weighted MAD, 1.4826 MAD0 */
double
gsl_wrstat_mad(const gsl_wrstat_workspace * w)
{
  return KLL_MAD_SCALE * gsl_wrstat_mad0(w);
}

int
gsl_wrstat_reset(gsl_wrstat_workspace * w)
{
  w->W = 0.0;
  w->W2 = 0.0;
  w->min = 0.0;
  w->max = 0.0;
  w->wmean = 0.0;
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
}

/*
wrstat_insert()
  Insert x with weight wt into level h, which is kept sorted with
cumulative weights, compacting first if the buffers are full: the items
of levels 0..h-1, and the items <= x of level h, move one slot toward
the free space in front of level 0, and the cumulative weights of the
items > x of level h increase by wt
*/

static int
wrstat_insert(const double wt, const double x, const size_t h, gsl_wrstat_workspace *w)
{
  double *bx = w->buf_x;
  double *bw = w->buf_w;
  double cprev;
  size_t a, p, i;
  int status;

  status = kll_grow_levels(h + 1, w->cap, w->lev, &(w->nlevels), &(w->budget));
  if (status)
    return status;

  /* make room if the buffers are full */
  if (w->nmax - w->lev[0] >= w->budget)
    {
      status = wrstat_compress(w);
      if (status)
        return status;
    }

  a = w->lev[0];
  p = kll_bound(bx, w->lev[h], w->lev[h + 1], x, 1);
  cprev = (p > w->lev[h]) ? bw[p - 1] : 0.0;

  memmove(bx + a - 1, bx + a, (p - a) * sizeof(double));
  memmove(bw + a - 1, bw + a, (p - a) * sizeof(double));
  bx[p - 1] = x;
  bw[p - 1] = cprev + wt;

  for (i = p; i < w->lev[h + 1]; ++i)
    bw[i] += wt;

  for (i = 0; i <= h; ++i)
    --(w->lev[i]);

  return GSL_SUCCESS;
}

/*
wrstat_merge_moments()
  Combine n, min, max, weights, weighted mean and central moments of src
into dest, with dest as dataset A and src as dataset B, eq. (3.1)-(3.4)
of [1] (see rstat_merge_moments() in rstat.c, with n replaced by W)
*/

static void
wrstat_merge_moments(const gsl_wrstat_workspace *src, gsl_wrstat_workspace *dest)
{
  if (dest->n == 0)
    {
      dest->min = src->min;
      dest->max = src->max;
    }
  else
    {
      if (src->min < dest->min)
        dest->min = src->min;

      if (src->max > dest->max)
        dest->max = src->max;
    }

  if (dest->W == 0.0)
    {
      dest->wmean = src->wmean;
      dest->M2 = src->M2;
      dest->M3 = src->M3;
      dest->M4 = src->M4;
    }
  else if (src->W > 0.0)
    {
      const double wa = dest->W;
      const double wb = src->W;
      const double W = wa + wb;
      const double d = src->wmean - dest->wmean;
      const double d2 = d * d;

      dest->M4 += src->M4 + wa * wb * (wa * wa - wa * wb + wb * wb) * d2 * d2 / (W * W * W) +
                  6.0 * d2 * (wa * wa * src->M2 + wb * wb * dest->M2) / (W * W) +
                  4.0 * d * (wa * src->M3 - wb * dest->M3) / W;
      dest->M3 += src->M3 + wa * wb * (wa - wb) * d2 * d / (W * W) +
                  3.0 * d * (wa * src->M2 - wb * dest->M2) / W;
      dest->M2 += src->M2 + wa * wb * d2 / W;
      dest->wmean += wb * d / W;
    }

  dest->W += src->W;
  dest->W2 += src->W2;
  dest->n += src->n;
}

/*
wrstat_compress()
  Compact the lowest level at or above its capacity, freeing at least one
slot in the buffers. The odd item (if any) stays at level h; the m = s1/2
resampling points produce half <= m items in [a1 .. a1+half-1], which are
merged in place with level h+1 into [b-half .. e-1], and levels 0..h-1 and
the odd item are shifted up by s1 - half >= m.
*/

static int
wrstat_compress(gsl_wrstat_workspace *w)
{
  size_t *lev = w->lev;
  double *bx = w->buf_x;
  double *bw = w->buf_w;
  size_t h;
  int status;

  status = kll_compact_level(w->cap, lev, &(w->nlevels), &(w->budget), &h);
  if (status)
    return status;

  /* compaction works on raw weights */
  wrstat_to_raw(w, h);
  wrstat_to_raw(w, h + 1);

  {
    const size_t a = lev[h];
    const size_t b = lev[h + 1];
    const size_t e = lev[h + 2];
    const size_t s = b - a;
    const size_t a1 = a + (s & 1);          /* odd item at a stays at level h */
    const size_t s1 = b - a1;               /* number of items to resample, even */
    const size_t m = s1 / 2;                /* number of resampling points */
    const double u = kll_uniform(&(w->rng));
    double Wlevel = 0.0, U, cum;
    size_t half = 0;                        /* number of items produced */
    size_t i, j, ia, ib, out, last, shift;

    for (i = a1; i < b; ++i)
      Wlevel += bw[i];

    U = Wlevel / m;

    /*
     * systematic resampling: point j at cumulative weight (j + u) U selects
     * the item whose cumulative weight interval contains it; points
     * selecting the same item are combined. Items are written to
     * [a1 .. a1+half-1], never beyond the item being read
     */
    i = 0;
    cum = bw[a1];
    last = s1;   /* no item selected yet */

    for (j = 0; j < m; ++j)
      {
        const double pos = (j + u) * U;

        while (pos >= cum && i < s1 - 1)
          cum += bw[a1 + (++i)];

        if (i == last)
          {
            bw[a1 + half - 1] += U;
          }
        else
          {
            bx[a1 + half] = bx[a1 + i];
            bw[a1 + half] = U;
            ++half;
            last = i;
          }
      }

    /* merge [a1 .. a1+half-1] with level h+1 = [b .. e-1] into [b-half .. e-1] */
    ia = a1;
    ib = b;
    out = b - half;

    while (ia < a1 + half && ib < e)
      {
        if (bx[ib] < bx[ia])
          {
            bx[out] = bx[ib];
            bw[out++] = bw[ib++];
          }
        else
          {
            bx[out] = bx[ia];
            bw[out++] = bw[ia++];
          }
      }

    while (ia < a1 + half)
      {
        bx[out] = bx[ia];
        bw[out++] = bw[ia++];
      }

    /* shift levels 0..h-1 and the odd item of level h up to close the gap */
    shift = s1 - half;

    memmove(bx + lev[0] + shift, bx + lev[0], (a1 - lev[0]) * sizeof(double));
    memmove(bw + lev[0] + shift, bw + lev[0], (a1 - lev[0]) * sizeof(double));

    for (i = 0; i <= h; ++i)
      lev[i] += shift;

    lev[h + 1] = b - half;
  }

  /*
   * levels h >= 1 hold cumulative weights; the odd item left at level h
   * (if any) is its own cumulative weight
   */
  wrstat_to_cum(w, h + 1);

  return GSL_SUCCESS;
}

/* total weight of the retained data: sum of the last cumulative weight of each level */
static double
wrstat_total_weight(const gsl_wrstat_workspace *w)
{
  double W = 0.0;
  size_t h;

  for (h = 0; h < w->nlevels; ++h)
    {
      if (w->lev[h + 1] > w->lev[h])
        W += w->buf_w[w->lev[h + 1] - 1];
    }

  return W;
}

/* convert cumulative weights of level h to raw weights, in place */
static void
wrstat_to_raw(gsl_wrstat_workspace *w, const size_t h)
{
  size_t i;

  for (i = w->lev[h + 1]; i > w->lev[h] + 1; --i)
    w->buf_w[i - 1] -= w->buf_w[i - 2];
}

/* convert raw weights of level h to cumulative weights, in place */
static void
wrstat_to_cum(gsl_wrstat_workspace *w, const size_t h)
{
  size_t i;

  for (i = w->lev[h] + 1; i < w->lev[h + 1]; ++i)
    w->buf_w[i] += w->buf_w[i - 1];
}
