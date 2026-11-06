/* rstat/test.c
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
#include <gsl/gsl_wrstat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_ieee_utils.h>

/* generate random data on [a,b] */
int
random_data(const double a, const double b, const size_t n, double data[], gsl_rng *r)
{
  size_t i;

  for (i = 0; i < n; ++i)
    data[i] = (b - a) * gsl_rng_uniform(r) + a;

  return 0;
}

/* orderings of test data streams for quantile tests */
enum
{
  ORDER_RANDOM = 0,    /* i.i.d. Gaussian */
  ORDER_SORTED,        /* Gaussian, sorted ascending */
  ORDER_REVERSE,       /* Gaussian, sorted descending */
  ORDER_ORBIT,         /* slowly varying periodic signal plus Cauchy noise */
  ORDER_NUM
};

static const char * order_names[ORDER_NUM] = { "random", "sorted", "reverse", "orbit" };

/* generate data stream of the given ordering */
static void
ordered_data(const int order, const size_t n, double data[], gsl_rng *r)
{
  size_t i;

  for (i = 0; i < n; ++i)
    {
      if (order == ORDER_ORBIT)
        {
          const double t = (double) i / 5400.0;
          data[i] = 20.0 * sin(2.0 * M_PI * t) + 5.0 * sin(2.0 * M_PI * t / 15.0) +
                    gsl_ran_cauchy(r, 1.0);
        }
      else
        {
          data[i] = gsl_ran_gaussian(r, 3.0) + 1.0;
        }
    }

  if (order == ORDER_SORTED || order == ORDER_REVERSE)
    gsl_sort(data, 1, n);

  if (order == ORDER_REVERSE)
    {
      for (i = 0; i < n / 2; ++i)
        {
          double tmp = data[i];
          data[i] = data[n - 1 - i];
          data[n - 1 - i] = tmp;
        }
    }
}

/*
test_rank()
  Normalized weighted rank of v in data x with weights wt (wt = NULL for
unit weights): weight of data < v plus half the weight of data == v,
divided by the total weight
*/
static double
test_rank(const double v, const double x[], const double wt[], const size_t n)
{
  double below = 0.0, equal = 0.0, W = 0.0;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      const double wi = wt ? wt[i] : 1.0;

      W += wi;

      if (x[i] < v)
        below += wi;
      else if (x[i] == v)
        equal += wi;
    }

  return (below + 0.5 * equal) / W;
}

/* fraction of data in x[0..n-1] which are <= v */
static double
test_cdf(const double v, const double x[], const size_t n)
{
  size_t count = 0, i;

  for (i = 0; i < n; ++i)
    count += (x[i] <= v);

  return (double) count / (double) n;
}

#include "test_rstat.c"
#include "test_wrstat.c"

int
main()
{
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  gsl_ieee_env_setup();

  test_rstat(r);
  test_wrstat(r);

  gsl_rng_free(r);

  exit (gsl_test_summary());
}
