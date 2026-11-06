/* rstat/gsl_rstat.h
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

#ifndef __GSL_RSTAT_H__
#define __GSL_RSTAT_H__

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
  double p;        /* p-quantile */
  double q[5];     /* heights q_i */
  int npos[5];     /* positions n_i */
  double np[5];    /* desired positions n_i' */
  double dnp[5];   /* increments dn_i' */
  size_t n;        /* number of data added */
} gsl_rstat_quantile_workspace;

gsl_rstat_quantile_workspace *gsl_rstat_quantile_alloc(const double p);
void gsl_rstat_quantile_free(gsl_rstat_quantile_workspace *w);
int gsl_rstat_quantile_reset(gsl_rstat_quantile_workspace *w);
int gsl_rstat_quantile_add(const double x, gsl_rstat_quantile_workspace *w);
double gsl_rstat_quantile_get(gsl_rstat_quantile_workspace *w);

typedef struct
{
  double min;      /* minimum value added */
  double max;      /* maximum value added */
  double mean;     /* current mean */
  double M2;       /* M_k = sum_{i=1..n} [ x_i - mean_n ]^k */
  double M3;
  double M4;
  size_t n;        /* number of data points added */

  /* quantile estimation (KLL algorithm), see rstat.c */
  size_t k;        /* accuracy parameter (capacity of top level) */
  size_t nmax;     /* maximum number of retained data, length of buf */
  size_t *cap;     /* capacities of levels by depth below the top level */
  size_t *lev;     /* level h occupies buf[lev[h] .. lev[h+1]-1] */
  size_t nlevels;  /* number of levels in use */
  size_t budget;   /* current maximum number of retained data */
  double *buf;     /* retained data */
  unsigned long rng; /* random number state */
} gsl_rstat_workspace;

gsl_rstat_workspace *gsl_rstat_alloc(void);
gsl_rstat_workspace *gsl_rstat_alloc2(const size_t k);
void gsl_rstat_free(gsl_rstat_workspace *w);
size_t gsl_rstat_n(const gsl_rstat_workspace *w);
int gsl_rstat_add(const double x, gsl_rstat_workspace *w);
double gsl_rstat_min(const gsl_rstat_workspace *w);
double gsl_rstat_max(const gsl_rstat_workspace *w);
double gsl_rstat_mean(const gsl_rstat_workspace *w);
double gsl_rstat_variance(const gsl_rstat_workspace *w);
double gsl_rstat_sd(const gsl_rstat_workspace *w);
double gsl_rstat_rms(const gsl_rstat_workspace *w);
double gsl_rstat_norm(const gsl_rstat_workspace *w);
double gsl_rstat_sd_mean(const gsl_rstat_workspace *w);
double gsl_rstat_median(const gsl_rstat_workspace *w);
double gsl_rstat_quantile(const double p, const gsl_rstat_workspace *w);
double gsl_rstat_mad0(const gsl_rstat_workspace *w);
double gsl_rstat_mad(const gsl_rstat_workspace *w);
double gsl_rstat_cdf(const double x, const gsl_rstat_workspace *w);
int gsl_rstat_merge(const gsl_rstat_workspace *src, gsl_rstat_workspace *dest);
double gsl_rstat_skew(const gsl_rstat_workspace *w);
double gsl_rstat_kurtosis(const gsl_rstat_workspace *w);
int gsl_rstat_reset(gsl_rstat_workspace *w);

__END_DECLS

#endif /* __GSL_RSTAT_H__ */
