#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>

int
main(void)
{
  const size_t N = 10000;
  const double p[3] = { 0.25, 0.5, 0.75 };
  double *data = malloc(N * sizeof(double));
  gsl_rstat_workspace *rstat_random = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_sorted = gsl_rstat_alloc();
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  size_t i, count = 0;

  /* add data to accumulator in random order; also store data for exact
   * comparisons */
  for (i = 0; i < N; ++i)
    {
      data[i] = gsl_ran_rayleigh(r, 1.0);
      gsl_rstat_add(data[i], rstat_random);
    }

  /* add the same data to a second accumulator in ascending order */
  gsl_sort(data, 1, N);

  for (i = 0; i < N; ++i)
    gsl_rstat_add(data[i], rstat_sorted);

  printf ("The dataset is %g, %g, %g, %g, %g, ...\n",
          data[0], data[1], data[2], data[3], data[4]);

  /* quartiles: exact values and estimates from both accumulators */
  for (i = 0; i < 3; ++i)
    {
      double exact = gsl_stats_quantile_from_sorted_data(data, 1, N, p[i]);
      double est_random = gsl_rstat_quantile(p[i], rstat_random);
      double est_sorted = gsl_rstat_quantile(p[i], rstat_sorted);

      printf ("%.2f quantile: exact = %.5f, estimated = %.5f (random order), %.5f (sorted order)\n",
              p[i], exact, est_random, est_sorted);
    }

  /* empirical CDF: fraction of data <= 1 */
  for (i = 0; i < N; ++i)
    count += (data[i] <= 1.0);

  printf ("fraction of data <= 1: exact = %.5f, estimated = %.5f\n",
          (double) count / (double) N, gsl_rstat_cdf(1.0, rstat_random));

  gsl_rstat_free(rstat_random);
  gsl_rstat_free(rstat_sorted);
  gsl_rng_free(r);
  free(data);

  return 0;
}
