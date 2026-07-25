/* deriv/gsl_deriv_params.h
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

#ifndef __GSL_DERIV_PARAMS_H__
#define __GSL_DERIV_PARAMS_H__
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv_stencil.h>
#include <gsl/gsl_deriv_step.h>

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
  unsigned int deriv_order;            /* derivative order (<= GSL_DERIV_MAX_DERIV_ORDER) */
  unsigned int method_order;           /* method order (<= GSL_DERIV_MAX_METHOD_ORDER) */
  gsl_deriv_stencil_t stencil_type;    /* stencil type, GSL_DERIV_STENCIL_xxx */
  unsigned int romberg_terms;          /* number of Romberg terms */

  gsl_deriv_step_policy_t step_policy; /* step policy */
  double min_step;                     /* minimum step size (in units of x) */
  double max_step;                     /* maximum step size (in units of x) */
  size_t nstep;                        /* number of step sizes */
  gsl_vector * scales;                 /* user-supplied scales, same length as x */
} gsl_deriv_parameters;

__END_DECLS

#endif /* __GSL_DERIV_PARAMS_H__ */
