/* deriv/gsl_deriv_jacobian.h
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

#ifndef __GSL_DERIV_JACOBIAN_H__
#define __GSL_DERIV_JACOBIAN_H__
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_deriv_constants.h>
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

/* definition of an arbitrary function F: R^n -> R^m */
typedef struct
{
  int (*f) (const gsl_vector * x, void * params, gsl_vector * f);
  size_t n;      /* number of inputs */
  size_t m;      /* number of outputs */
  void * params;
} gsl_deriv_vec_function;

#define GSL_DERIV_VEC_FN_EVAL(F,x,y) (*((F)->f))(x,(F)->params,y)

typedef struct
{
  gsl_deriv_parameters params;
  
  size_t nrows;                                   /* number of rows of Jacobian matrix */
  size_t ncols;                                   /* number of columns of Jacobian matrix */

  double h[GSL_DERIV_MAX_STEPS];                  /* step sizes (same units as input x) */
  double step_ratio;                              /* step ratio */
  size_t nstep;                                   /* number of steps */

  double stencil[GSL_DERIV_MAX_STENCIL];          /* finite difference stencil */
  size_t stencil_size;                            /* stencil size */

  /* Fornberg finite difference weights for all derivative orders */
  double fornberg_weights[GSL_DERIV_MAX_STENCIL*(GSL_DERIV_MAX_DERIV_ORDER+1)];
  gsl_vector_view unit_weights;                   /* unit Fornberg weights for desired deriv_order */

  gsl_matrix * der_init;                          /* initial derivative estimates, nrows-by-MAX_STEPS */
  gsl_matrix * err_init;                          /* initial derivative roundoff error, nrows-by-MAX_STEPS */
  gsl_vector * work2l;                            /* workspace, length 2*nrows */
  gsl_vector * workn;                             /* workspace, length ncols */
  size_t nromb;                                   /* number of Romberg windows */
  void * romb_state;                              /* Romberg state variables */
  gsl_vector * der_romb;                          /* Romberg derivative estimates, length nromb */
  gsl_vector * err_romb;                          /* Romberg error estimates, length nromb */
} gsl_deriv_jacobian_workspace;

gsl_deriv_jacobian_workspace * gsl_deriv_jacobian_alloc(const gsl_deriv_parameters * params,
                                                        const size_t nrows, const size_t ncols);

void gsl_deriv_jacobian_free(gsl_deriv_jacobian_workspace * w);

int gsl_deriv_jacobian_fixed(const gsl_deriv_vec_function * F, const gsl_vector * x,
                             const gsl_vector * h, gsl_matrix * J, gsl_matrix * esterr,
                             unsigned int * neval, gsl_deriv_jacobian_workspace * w);

int gsl_deriv_jacobian(const gsl_deriv_vec_function * F, const gsl_vector * x,
                       gsl_matrix * J, gsl_matrix * esterr, gsl_matrix * final_step, unsigned int * neval,
                       gsl_deriv_jacobian_workspace * w);

__END_DECLS

#endif /* __GSL_DERIV_JACOBIAN_H__ */
