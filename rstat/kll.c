/* rstat/kll.c
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
 * Internal routines for the KLL quantile estimator (Karnin, Lang and
 * Liberty, 2016), shared by rstat.c and wrstat.c, which include this file.
 *
 * The retained data are stored in a stack of levels h = 0, 1, ..., H-1 in
 * a single buffer. Level h occupies buf[lev[h] .. lev[h+1]-1], level 0
 * first, with the free space in front of level 0. All levels are kept
 * sorted: new data are inserted into level 0 at their sorted position,
 * moving the smaller items of level 0 one slot toward the free space.
 * The capacity of a level at depth d = H-1-h below the top level is
 *
 * cap(d) = max(KLL_MIN_WIDTH, round(k (2/3)^d))
 *
 * When the buffer holds sum_{d<H} cap(d) items, the lowest level at or
 * above its capacity is compacted into the level above it.
 *
 * Queries see the retained data in ascending order, x_0 <= ... <= x_{m-1},
 * with weights w_i > 0 and cumulative weights C_i = w_0 + ... + w_i. For
 * a target weight t, define
 *
 * lo(t) = min { x_i : C_i >= t }
 * hi(t) = min { x_i : C_i >  t }
 *
 * The median is (lo(W/2) + hi(W/2)) / 2 with W = C_{m-1}; for unit
 * weights this is the usual median (average of the two middle values
 * for even m).
 *
 * Since all levels are sorted, the weight of the data <= x is the sum over levels of
 * the weight of the first j_h items of level h, where j_h is found by
 * binary search; for unweighted data this is j_h 2^h, and for weighted
 * data buf_w holds the cumulative weights within each level. A binary
 * search over the items of the largest level, using these ranks,
 * brackets the answer, and the few items in the bracket are then
 * visited in ascending order by merging the levels, with one cursor per
 * level. No sorted copy of the data is stored.
 */

#define KLL_MAX_LEVELS       64                  /* maximum number of levels */
#define KLL_MIN_WIDTH        2                   /* minimum capacity of a level */
#define KLL_MIN_K            8                   /* minimum accuracy parameter k */
#define KLL_DEFAULT_K        200                 /* default accuracy parameter */
#define KLL_SEED             0x2545f491UL        /* random number seed */
#define KLL_MAD_SCALE        1.482602218505602   /* 1 / Phi^{-1}(3/4), as in gsl_stats_mad() */

/*
kll_uniform()
  Return a uniform random number in [0,1), using the 32-bit Mulberry32
generator. Arithmetic is done modulo 2^32 with unsigned long (at least
32 bits in C89), so the sequence is the same on all platforms.
*/

static double
kll_uniform(unsigned long * state)
{
  const unsigned long mask = 0xffffffffUL;
  unsigned long z;

  *state = (*state + 0x6d2b79f5UL) & mask;
  z = *state;
  z = ((z ^ (z >> 15)) * (z | 1UL)) & mask;
  z = (z ^ ((z + (((z ^ (z >> 7)) * (z | 61UL)) & mask)) & mask)) & mask;
  z = (z ^ (z >> 14)) & mask;

  return (double) z / 4294967296.0;
}

/*
kll_capacity()
  Compute capacities cap[d], d = 0, ..., KLL_MAX_LEVELS-1, of the levels
at depth d below the top level, and return their sum, which is the
maximum number of retained items
*/

static size_t
kll_capacity(const size_t k, size_t * cap)
{
  size_t nmax = 0;
  size_t d;

  for (d = 0; d < KLL_MAX_LEVELS; ++d)
    {
      const size_t cd = (size_t) (k * pow(2.0 / 3.0, (double) d) + 0.5);

      cap[d] = GSL_MAX(cd, KLL_MIN_WIDTH);
      nmax += cap[d];
    }

  return nmax;
}

/*
kll_compact_level()
  Find the level h to compact in a full buffer: the lowest level at or
above its capacity, which exists since the number of items equals the
sum of the capacities. If h is the top level, a new empty level is added
above it.
*/

static int
kll_compact_level(const size_t * cap, size_t * lev, size_t * nlevels, size_t * budget, size_t * h)
{
  const size_t H = *nlevels;
  size_t i;

  for (i = 0; i < H - 1; ++i)
    {
      if (lev[i + 1] - lev[i] >= cap[H - 1 - i])
        break;
    }

  if (i == H - 1)
    {
      if (H == KLL_MAX_LEVELS)
        {
          GSL_ERROR ("maximum number of levels exceeded", GSL_EMAXITER);
        }

      lev[H + 1] = lev[H];
      *budget += cap[H];
      *nlevels = H + 1;
    }

  *h = i;

  return GSL_SUCCESS;
}

/* add empty levels on top until there are at least H levels */
static int
kll_grow_levels(const size_t H, const size_t * cap, size_t * lev, size_t * nlevels, size_t * budget)
{
  while (*nlevels < H)
    {
      const size_t n = *nlevels;

      if (n == KLL_MAX_LEVELS)
        {
          GSL_ERROR ("maximum number of levels exceeded", GSL_EMAXITER);
        }

      lev[n + 1] = lev[n];
      *budget += cap[n];
      *nlevels = n + 1;
    }

  return GSL_SUCCESS;
}

/* weights 2^h of unweighted items at levels h = 0, ..., KLL_MAX_LEVELS-1 */
static void
kll_level_weights(double * lw)
{
  size_t h;

  lw[0] = 1.0;

  for (h = 1; h < KLL_MAX_LEVELS; ++h)
    lw[h] = 2.0 * lw[h - 1];
}

/* level with the most items */
static size_t
kll_largest_level(const size_t * lev, const size_t nlevels)
{
  size_t L = 0;
  size_t h;

  for (h = 1; h < nlevels; ++h)
    {
      if (lev[h + 1] - lev[h] > lev[L + 1] - lev[L])
        L = h;
    }

  return L;
}

/* smallest and largest retained items; the levels must be sorted and not all empty */
static void
kll_extremes(const double * buf, const size_t * lev, const size_t nlevels,
             double * xmin, double * xmax)
{
  int found = 0;
  size_t h;

  for (h = 0; h < nlevels; ++h)
    {
      if (lev[h] == lev[h + 1])
        continue;

      if (!found || buf[lev[h]] < *xmin)
        *xmin = buf[lev[h]];

      if (!found || buf[lev[h + 1] - 1] > *xmax)
        *xmax = buf[lev[h + 1] - 1];

      found = 1;
    }
}

/*
kll_bound()
  Return the index of the first item of the sorted range buf[a .. b-1]
which is >= x (upper = 0) or > x (upper = 1), or b if there is none
*/

static size_t
kll_bound(const double * buf, size_t a, size_t b, const double x, const int upper)
{
  while (a < b)
    {
      const size_t mid = a + (b - a) / 2;

      if (buf[mid] < x || (upper && buf[mid] == x))
        a = mid + 1;
      else
        b = mid;
    }

  return a;
}

/*
kll_cumw()
  Return the total weight of the items buf[lev[h] .. i-1] of level h;
buf_w holds cumulative weights within each level, or is NULL for
unweighted data with weights lw[h] = 2^h
*/

static double
kll_cumw(const double * buf_w, const double * lw, const size_t * lev,
         const size_t h, const size_t i)
{
  if (i == lev[h])
    return 0.0;
  else if (buf_w)
    return buf_w[i - 1];
  else
    return (double) (i - lev[h]) * lw[h];
}

/* total weight of the retained items <= x */
static double
kll_rank(const double x, const double * buf, const double * buf_w, const double * lw,
         const size_t * lev, const size_t nlevels)
{
  double r = 0.0;
  size_t h;

  for (h = 0; h < nlevels; ++h)
    r += kll_cumw(buf_w, lw, lev, h, kll_bound(buf, lev[h], lev[h + 1], x, 1));

  return r;
}

/*
kll_select()
  Compute lo(t) and hi(t) of the retained items, which must be sorted
within each level and not all empty; if no item satisfies a condition,
the largest item is returned.

A binary search over the items of the largest level finds the largest
item v of that level whose rank (total weight of items <= v) is < t, so
that lo(t) and hi(t) are > v. The items > v are then visited in
ascending order by merging the levels, starting from the first item > v
in each level, until hi(t) is found.
*/

static void
kll_select(const double t, const double * buf, const double * buf_w,
           const size_t * lev, const size_t nlevels, double * lo, double * hi)
{
  const size_t L = kll_largest_level(lev, nlevels);
  size_t pos[KLL_MAX_LEVELS]; /* next item of each level */
  double lw[KLL_MAX_LEVELS];  /* weights of unweighted items by level */
  double cum = 0.0;
  double x = 0.0;
  int found_lo = 0;
  size_t a = lev[L], b = lev[L + 1];
  size_t h;

  kll_level_weights(lw);

  /* first item of level L with rank >= t */
  while (a < b)
    {
      const size_t mid = a + (b - a) / 2;

      if (kll_rank(buf[mid], buf, buf_w, lw, lev, nlevels) < t)
        a = mid + 1;
      else
        b = mid;
    }

  for (h = 0; h < nlevels; ++h)
    {
      /* skip the items <= v = buf[a-1], whose total weight is < t */
      pos[h] = (a > lev[L]) ? kll_bound(buf, lev[h], lev[h + 1], buf[a - 1], 1) : lev[h];
      cum += kll_cumw(buf_w, lw, lev, h, pos[h]);
    }

  while (1)
    {
      size_t best = nlevels; /* level with smallest next item */

      for (h = 0; h < nlevels; ++h)
        {
          if (pos[h] < lev[h + 1] && (best == nlevels || buf[pos[h]] < buf[pos[best]]))
            best = h;
        }

      if (best == nlevels)
        break;

      x = buf[pos[best]];
      cum += kll_cumw(buf_w, lw, lev, best, pos[best] + 1) - kll_cumw(buf_w, lw, lev, best, pos[best]);
      ++pos[best];

      if (!found_lo && cum >= t)
        {
          *lo = x;
          found_lo = 1;
        }

      if (cum > t)
        {
          *hi = x;
          return;
        }
    }

  /* t >= total weight: return largest item */
  if (!found_lo)
    *lo = x;

  *hi = x;
}

/*
kll_dev_bounds()
  For the items of level h, split at split[h] (first item >= m), return
the range [*l, *r) of items with deviation <= d: m - buf[i] <= d for items
left of m, buf[i] - m <= d for items right of m. The deviations are
computed as in kll_mad_select(), so the comparisons are consistent.
*/

static void
kll_dev_bounds(const double m, const double d, const double * buf, const size_t * lev,
               const size_t * split, const size_t h, size_t * l, size_t * r)
{
  size_t a = lev[h], b = split[h];

  /* first item left of m with m - buf[i] <= d */
  while (a < b)
    {
      const size_t mid = a + (b - a) / 2;

      if (m - buf[mid] <= d)
        b = mid;
      else
        a = mid + 1;
    }

  *l = a;

  /* first item right of m with buf[i] - m > d */
  a = split[h];
  b = lev[h + 1];

  while (a < b)
    {
      const size_t mid = a + (b - a) / 2;

      if (buf[mid] - m > d)
        b = mid;
      else
        a = mid + 1;
    }

  *r = a;
}

/* total weight of the retained items with deviation |x - m| <= d */
static double
kll_dev_rank(const double m, const double d, const double * buf, const double * buf_w,
             const double * lw, const size_t * lev, const size_t * split, const size_t nlevels)
{
  double r = 0.0;
  size_t h, l, u;

  for (h = 0; h < nlevels; ++h)
    {
      kll_dev_bounds(m, d, buf, lev, split, h, &l, &u);
      r += kll_cumw(buf_w, lw, lev, h, u) - kll_cumw(buf_w, lw, lev, h, l);
    }

  return r;
}

/*
kll_mad_select()
  Compute lo(t) and hi(t) of the weighted absolute deviations |x_i - m|
of the retained items, which must be sorted within each level and not
all empty.

The deviations m - x_i of items left of m increase moving left, and
x_i - m of items right of m increase moving right. Binary searches over
the left and right items of the largest level find the largest deviation
d of that level whose rank (total weight of items with deviation <= d)
is < t. The remaining deviations are then visited in ascending order by
walking outward from the items with deviation <= d: in each level a left
cursor moves down and a right cursor moves up, and each step takes the
closer of the nearest remaining items on the left and on the right.
*/

static void
kll_mad_select(const double m, const double t, const double * buf, const double * buf_w,
               const size_t * lev, const size_t nlevels, double * lo, double * hi)
{
  const size_t L = kll_largest_level(lev, nlevels);
  size_t split[KLL_MAX_LEVELS]; /* first item >= m in each level */
  size_t left[KLL_MAX_LEVELS];  /* items buf[lev[h] .. left[h]-1] remain left of m */
  size_t right[KLL_MAX_LEVELS]; /* items buf[right[h] .. lev[h+1]-1] remain right of m */
  double lw[KLL_MAX_LEVELS];    /* weights of unweighted items by level */
  double cum = 0.0;
  double d = 0.0;
  double dstar = 0.0;
  int found_d = 0;
  int found_lo = 0;
  size_t a, b, h;

  kll_level_weights(lw);

  for (h = 0; h < nlevels; ++h)
    split[h] = kll_bound(buf, lev[h], lev[h + 1], m, 0);

  /*
   * left items of level L: deviations decrease with index, so find the
   * first index whose deviation has rank < t
   */
  a = lev[L];
  b = split[L];

  while (a < b)
    {
      const size_t mid = a + (b - a) / 2;

      if (kll_dev_rank(m, m - buf[mid], buf, buf_w, lw, lev, split, nlevels) < t)
        b = mid;
      else
        a = mid + 1;
    }

  if (a < split[L])
    {
      dstar = m - buf[a];
      found_d = 1;
    }

  /*
   * right items of level L: deviations increase with index, so find the
   * last index whose deviation has rank < t
   */
  a = split[L];
  b = lev[L + 1];

  while (a < b)
    {
      const size_t mid = a + (b - a) / 2;

      if (kll_dev_rank(m, buf[mid] - m, buf, buf_w, lw, lev, split, nlevels) < t)
        a = mid + 1;
      else
        b = mid;
    }

  if (a > split[L] && (!found_d || buf[a - 1] - m > dstar))
    {
      dstar = buf[a - 1] - m;
      found_d = 1;
    }

  /* skip the items with deviation <= dstar, whose total weight is < t */
  for (h = 0; h < nlevels; ++h)
    {
      if (found_d)
        {
          kll_dev_bounds(m, dstar, buf, lev, split, h, &left[h], &right[h]);
          cum += kll_cumw(buf_w, lw, lev, h, right[h]) - kll_cumw(buf_w, lw, lev, h, left[h]);
        }
      else
        {
          left[h] = split[h];
          right[h] = split[h];
        }
    }

  while (1)
    {
      size_t bl = nlevels; /* level with largest item left of m */
      size_t br = nlevels; /* level with smallest item right of m */
      size_t i;

      for (h = 0; h < nlevels; ++h)
        {
          if (left[h] > lev[h] && (bl == nlevels || buf[left[h] - 1] > buf[left[bl] - 1]))
            bl = h;

          if (right[h] < lev[h + 1] && (br == nlevels || buf[right[h]] < buf[right[br]]))
            br = h;
        }

      if (bl == nlevels && br == nlevels)
        break;

      if (bl != nlevels && (br == nlevels || m - buf[left[bl] - 1] <= buf[right[br]] - m))
        {
          i = --left[bl];
          d = m - buf[i];
          cum += kll_cumw(buf_w, lw, lev, bl, i + 1) - kll_cumw(buf_w, lw, lev, bl, i);
        }
      else
        {
          i = right[br]++;
          d = buf[i] - m;
          cum += kll_cumw(buf_w, lw, lev, br, i + 1) - kll_cumw(buf_w, lw, lev, br, i);
        }

      if (!found_lo && cum >= t)
        {
          *lo = d;
          found_lo = 1;
        }

      if (cum > t)
        {
          *hi = d;
          return;
        }
    }

  /* t >= total weight: return largest deviation */
  if (!found_lo)
    *lo = d;

  *hi = d;
}
