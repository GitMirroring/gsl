/* rstat/gsl_wrstat.h
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

#ifndef __GSL_WRSTAT_H__
#define __GSL_WRSTAT_H__

#include <stdlib.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct
{
  double W;        /* sum(w_i) */
  double W2;       /* sum(w_i^2) */
  double min;      /* minimum value added */
  double max;      /* maximum value added */
  double wmean;    /* current weighted mean */
  double M2;       /* M_k = sum_{i=1..n} w_i [ x_i - wmean_n ]^k */
  double M3;
  double M4;
  size_t n;        /* number of data points added */

  /* quantile estimation (weighted KLL algorithm), see wrstat.c */
  size_t k;        /* accuracy parameter (capacity of top level) */
  size_t nmax;     /* maximum number of retained data, length of buf_x, buf_w */
  size_t *cap;     /* capacities of levels by depth below the top level */
  size_t *lev;     /* level h occupies buf_x[lev[h] .. lev[h+1]-1] */
  size_t nlevels;  /* number of levels in use */
  size_t budget;   /* current maximum number of retained data */
  double *buf_x;   /* retained data */
  double *buf_w;   /* cumulative weights of retained data within each level */
  unsigned long rng; /* random number state */
} gsl_wrstat_workspace;

gsl_wrstat_workspace *gsl_wrstat_alloc(void);
gsl_wrstat_workspace *gsl_wrstat_alloc2(const size_t k);
void gsl_wrstat_free(gsl_wrstat_workspace * w);
size_t gsl_wrstat_n(const gsl_wrstat_workspace * w);
int gsl_wrstat_add(const double wt, const double x, gsl_wrstat_workspace * w);
double gsl_wrstat_min(const gsl_wrstat_workspace * w);
double gsl_wrstat_max(const gsl_wrstat_workspace * w);
double gsl_wrstat_mean(const gsl_wrstat_workspace * w);
double gsl_wrstat_variance(const gsl_wrstat_workspace * w);
double gsl_wrstat_sd(const gsl_wrstat_workspace * w);
double gsl_wrstat_rms(const gsl_wrstat_workspace * w);
double gsl_wrstat_norm(const gsl_wrstat_workspace * w);
double gsl_wrstat_sd_mean(const gsl_wrstat_workspace * w);
double gsl_wrstat_skew(const gsl_wrstat_workspace * w);
double gsl_wrstat_kurtosis(const gsl_wrstat_workspace * w);
double gsl_wrstat_neff(const gsl_wrstat_workspace * w);
double gsl_wrstat_median(const gsl_wrstat_workspace * w);
double gsl_wrstat_quantile(const double p, const gsl_wrstat_workspace * w);
double gsl_wrstat_mad0(const gsl_wrstat_workspace * w);
double gsl_wrstat_mad(const gsl_wrstat_workspace * w);
double gsl_wrstat_cdf(const double x, const gsl_wrstat_workspace * w);
int gsl_wrstat_merge(const gsl_wrstat_workspace * src, gsl_wrstat_workspace * dest);
int gsl_wrstat_reset(gsl_wrstat_workspace * w);

__END_DECLS

#endif /* __GSL_WRSTAT_H__ */
