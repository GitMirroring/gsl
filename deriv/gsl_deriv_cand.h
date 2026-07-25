/* deriv/gsl_deriv_cand.h
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

#ifndef __GSL_DERIV_CAND_H__
#define __GSL_DERIV_CAND_H__
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_types.h>

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
  size_t nkeep;                       /* number of derivative candidates to keep */
  size_t index[GSL_DERIV_MAX_STEPS];  /* indices of kept candidates */
} gsl_deriv_cand_t;

typedef struct
{
  const char * name;
  int (*filter)(const gsl_vector * der, gsl_deriv_cand_t * cand);
} gsl_deriv_cand_filter;

typedef struct
{
  const char * name;
  int (*select)(const gsl_deriv_cand_t * cand, const gsl_vector * der,
                const gsl_vector * err, const gsl_vector * steps,
                double * der_est, double * err_est, double * step_est);
} gsl_deriv_cand_select;

/* candidate filter types */
GSL_VAR const gsl_deriv_cand_filter * gsl_deriv_cand_filter_mad;
GSL_VAR const gsl_deriv_cand_filter * gsl_deriv_cand_filter_trim;
GSL_VAR const gsl_deriv_cand_filter * gsl_deriv_cand_filter_none;

/* candidate selector types */
GSL_VAR const gsl_deriv_cand_select * gsl_deriv_cand_select_minerr;

__END_DECLS

#endif /* __GSL_DERIV_CAND_H__ */
