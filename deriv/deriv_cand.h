/* deriv/deriv_cand.h
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

#ifndef __DERIV_CAND_H__
#define __DERIV_CAND_H__
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_deriv.h>

/* candidate filters and selectors used inside the deriv module (not
 * installed: the filters are still experimental) */

typedef struct
{
  const char * name;
  int (*filter)(const gsl_vector * der, const gsl_vector * err, gsl_deriv_cand_t * cand);
} gsl_deriv_cand_filter;

typedef struct
{
  const char * name;
  int (*select)(const gsl_deriv_cand_t * cand, const gsl_vector * der,
                const gsl_vector * err, const gsl_vector * steps,
                double * der_est, double * err_est, double * step_est, size_t * idx_est);
} gsl_deriv_cand_select;

/* candidate filter types */
extern const gsl_deriv_cand_filter * gsl_deriv_cand_filter_mad;
extern const gsl_deriv_cand_filter * gsl_deriv_cand_filter_none;

/* candidate selector types */
extern const gsl_deriv_cand_select * gsl_deriv_cand_select_minerr;

#endif /* __DERIV_CAND_H__ */
