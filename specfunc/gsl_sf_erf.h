/* specfunc/gsl_sf_erf.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

#ifndef __GSL_SF_ERF_H__
#define __GSL_SF_ERF_H__

#include <gsl/gsl_sf_result.h>

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

/* erf/erf.c */
int gsl_sf_erf_e(const double x, gsl_sf_result * result);
int gsl_sf_erf_Z_e(const double x, gsl_sf_result * result);
int gsl_sf_erf_Q_e(const double x, gsl_sf_result * result);
int gsl_sf_hazard_e(const double x, gsl_sf_result * result);
double gsl_sf_erf(const double x);
double gsl_sf_erf_Z(const double x);
double gsl_sf_erf_Q(const double x);
double gsl_sf_hazard(const double x);

/* erf/erfc.c */
int gsl_sf_erfc_e(const double x, gsl_sf_result * result);
int gsl_sf_log_erfc_e(const double x, gsl_sf_result * result);
double gsl_sf_erfc(const double x);
double gsl_sf_log_erfc(const double x);

/* erf/erfcx.c */
int gsl_sf_erfcx_e(const double x, gsl_sf_result * result);
double gsl_sf_erfcx(const double x);

__END_DECLS

#endif /* __GSL_SF_ERF_H__ */
