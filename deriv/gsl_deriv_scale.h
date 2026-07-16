/* deriv/gsl_deriv_scale.h
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

#ifndef __GSL_DERIV_SCALE_H__
#define __GSL_DERIV_SCALE_H__
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_deriv_params.h>

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

/* parameters that define the dimensionless scale sequence */
typedef struct
{
  double max_scale; /* maximum scale multiplier */
  double ratio;     /* geometric ratio */
  size_t nscale;    /* number of scale multipliers */
} gsl_deriv_scale_parameters;

typedef struct
{
  const char * name;
  int (*init)(const gsl_deriv_parameters * params, gsl_deriv_scale_parameters * sparams);
} gsl_deriv_scale_generator;

/* scale generator types */
GSL_VAR const gsl_deriv_scale_generator * gsl_deriv_scale_user;

int gsl_deriv_scale_fill(const gsl_deriv_scale_parameters * params, double scales[]);

__END_DECLS

#endif /* __GSL_DERIV_SCALE_H__ */
