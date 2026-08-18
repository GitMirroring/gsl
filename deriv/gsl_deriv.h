/* deriv/gsl_deriv.h
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

#ifndef __GSL_DERIV_H__
#define __GSL_DERIV_H__
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_deriv_constants.h>
#include <gsl/gsl_deriv_cand.h>
#include <gsl/gsl_deriv_jacobian.h>
#include <gsl/gsl_deriv_params.h>
#include <gsl/gsl_deriv_stencil.h>

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

/* diagnostic information for adaptive differentiation */
typedef struct
{
  gsl_deriv_parameters params;

  double h[GSL_DERIV_MAX_STEPS];                  /* step sizes (same units as input x) */
  double step_ratio;                              /* step ratio */
  size_t nstep;                                   /* number of steps */

  double stencil[GSL_DERIV_MAX_STENCIL];          /* finite difference stencil */
  size_t stencil_size;                            /* stencil size */

  /* Fornberg finite difference weights for all derivative orders */
  double fornberg_weights[GSL_DERIV_MAX_STENCIL*(GSL_DERIV_MAX_DERIV_ORDER+1)];

  double x;                                       /* point at which to compute derivative */
  double der_init[GSL_DERIV_MAX_STEPS];           /* initial derivative estimates */
  double err_init[GSL_DERIV_MAX_STEPS];           /* initial derivative roundoff error */
  double der_romb[GSL_DERIV_MAX_STEPS];           /* derivative estimates after Romberg extrapolation */
  double err_romb[GSL_DERIV_MAX_STEPS];           /* error estimates for Romberg extrapolated derivatives */
  size_t nromb;                                   /* number of Romberg windows */
  gsl_deriv_cand_t cand;                          /* derivative candidates after filtering */

  double der_select;                              /* final selected derivative */
  double err_select;                              /* final selected error */
  double step_select;                             /* final selected step */
} gsl_deriv_diag;

/* definition of an arbitrary function f: R^n -> R */
typedef struct
{
  double (*f) (const gsl_vector * x, void * params);
  size_t n;
  void * params;
} gsl_deriv_function;

#define GSL_DERIV_FN_EVAL(F,x) (*((F)->f))(x,(F)->params)

gsl_deriv_parameters gsl_deriv_default_parameters(const unsigned int deriv_order, const unsigned int method_order);

int gsl_deriv_fixed(const gsl_deriv_parameters * params, const gsl_function * f, const double x,
                    const double h, double * result, double * esterr, unsigned int * neval);

int gsl_deriv_fixed_v(const gsl_deriv_parameters * params, const gsl_function * f, const gsl_vector * x,
                      const double h, gsl_vector * result, gsl_vector * esterr, unsigned int * neval);

int gsl_deriv_adaptive(const gsl_deriv_parameters * params, const gsl_function * f, const double x,
                       double * result, double * esterr, double * final_step, unsigned int * neval,
                       gsl_deriv_diag * diag);

int gsl_deriv_adaptive_v(const gsl_deriv_parameters * params, const gsl_function * f, const gsl_vector * x,
                         gsl_vector * result, gsl_vector * esterr, gsl_vector * final_step, unsigned int * neval,
                         gsl_deriv_diag * diag);

int gsl_deriv_fornberg_weights(const double nodes[], const unsigned int n,
                               const unsigned int max_deriv, const double x0,
                               double weights[]);

int gsl_deriv_print_diag(FILE * fp, const gsl_deriv_diag * diag);

/* directionaldiff.c */

int gsl_deriv_directional_fixed(const gsl_deriv_parameters * params, const gsl_deriv_function * f,
                                const gsl_vector * x, const gsl_vector * v, const double h,
                                double * result, double * esterr, unsigned int * neval,
                                gsl_vector * work);

int gsl_deriv_directional(const gsl_deriv_parameters * params, const gsl_deriv_function * f,
                          const gsl_vector * x, const gsl_vector * v, double * result,
                          double * esterr, double * final_step, unsigned int * neval,
                          gsl_deriv_diag * diag, gsl_vector * work);

/* gradient.c */

int gsl_deriv_gradient_fixed(const gsl_deriv_parameters * params, const gsl_deriv_function * f,
                             const gsl_vector * x, const gsl_vector * h,
                             gsl_vector * result, gsl_vector * esterr, unsigned int * neval,
                             gsl_vector * work);

int gsl_deriv_gradient(const gsl_deriv_parameters * params, const gsl_deriv_function * f,
                       const gsl_vector * x, gsl_vector * result, gsl_vector * esterr,
                       gsl_vector * final_step, unsigned int * neval, gsl_vector * work);

/* deriv_old.c */

int gsl_deriv_central (const gsl_function *f,
                       double x, double h,
                       double *result, double *abserr);

int gsl_deriv_backward (const gsl_function *f,
                        double x, double h,
                        double *result, double *abserr);

int gsl_deriv_forward (const gsl_function *f,
                       double x, double h,
                       double *result, double *abserr);

__END_DECLS

#endif /* __GSL_DERIV_H__ */
