/* deriv/deriv.c
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>

#define ROMBERG_MAX_TERMS          3
#define ROMBERG_MAX_M              (ROMBERG_MAX_TERMS + 2)
#define ROMBERG_MAX_N              (ROMBERG_MAX_TERMS + 1)

typedef struct
{
  unsigned int expon[ROMBERG_MAX_TERMS];       /* exponent powers appearing in the truncation error expansion, length nexpon */
  size_t nexpon;                               /* length of expon = romberg_terms */
  double QR_data[ROMBERG_MAX_M*ROMBERG_MAX_N]; /* QR factorization of LS matrix, size M-by-N */
  double T_data[ROMBERG_MAX_N*ROMBERG_MAX_N];  /* Householder reflector matrix, size N-by-N */
  double sqrt_cov1;                            /* cov1 = (R^T R)^{-1}]_11 */
} romb_state_t;

static int deriv_init(const gsl_deriv_parameters * params, size_t * stencil_size, double stencil[],
                      double weights[], gsl_vector_view * W);
static int fd_apply(const gsl_function * f, const double x0, const double dx,
                    const gsl_vector * weights, const double stencil[], const size_t stencil_size,
                    const unsigned int deriv_order, double * result, double * abserr_round, unsigned int * neval);
static int vec2mat(const gsl_vector * v, gsl_matrix * A);
static int romb_init(const gsl_deriv_parameters * params, const double scale_ratio, romb_state_t * state);
static int romb_extrap(const gsl_vector * der_init, gsl_vector * der_romb, gsl_vector * errest, romb_state_t * state);

/*
gsl_deriv_default_parameters()
  Return default parameters for adaptive differentiation
*/

gsl_deriv_parameters
gsl_deriv_default_parameters(void)
{
  gsl_deriv_parameters params;

  params.deriv_order = 1;
  params.method_order = 4;
  params.stencil_type = GSL_DERIV_STENCIL_CENTRAL;
  params.romberg_terms = 2;
  params.fixed_step = -1.0;

  /* default scale multipliers - similar to DERIVEST defaults */
  params.max_scale = 1.0e2;
  params.min_scale = 1.0e-6;
  params.nscale = 26;

  return params;
}

/*
gsl_deriv_scalar()
  Computing derivative of function of scalar input with reasonable defaults

Inputs: f           - function f(x)
        x           - point at which to compute derivative
        deriv_order - derivative order
        result      - (output) derivative result
        esterr      - (output) 95% uncertainty estimate of the derivative, such that
                        abs(result - f'(x)) < esterr

Return: success/error
*/

int
gsl_deriv_scalar(const gsl_function * f, const double x,
                 const unsigned int deriv_order, double * result, double * esterr)
{
  return gsl_deriv_scalar_fixed(f, x, -1.0, deriv_order, result, esterr);
}

/*
gsl_deriv_scalar_fixed()
  Differentiate a scalar function with a fixed step size

Inputs: f           - function f(x)
        x           - point at which to compute derivative
        h           - step size
        deriv_order - derivative order
        result      - (output) derivative result
        esterr      - (output) 95% uncertainty estimate of the derivative, such that
                        abs(result - f'(x)) < esterr

Return: success/error
*/

int
gsl_deriv_scalar_fixed(const gsl_function * f, const double x, const double h,
                       const unsigned int deriv_order, double * result, double * esterr)
{
  gsl_deriv_parameters params = gsl_deriv_default_parameters();
  double final_step;
  unsigned int neval;
  gsl_deriv_diag diag;

  params.deriv_order = deriv_order;
  params.fixed_step = h;

  return gsl_deriv_scalar_driver(&params, f, x, result, esterr, &final_step, &neval, &diag);
}

/*
gsl_deriv_scalar_driver()
  Top-level driver for differentiation algorithms for a function
of scalar input

Inputs: params     - derivative parameters
        f          - function to differentiate
        x          - point at which to compute derivative
        result     - (output) derivative result
        esterr     - (output) 95% uncertainty estimate of the derivative(s), such that
                       abs(result - f'(x)) < esterr
        final_step - (output) final overall step size chosen for each derivative
        neval      - (output) total number of function evaluations
        diag       - (output) diagnostic information for adaptive algorithm (may be NULL)

Return: success/error
*/

int
gsl_deriv_scalar_driver(const gsl_deriv_parameters * params, const gsl_function * f, const double x,
                        double * result, double * esterr, double * final_step, unsigned int * neval,
                        gsl_deriv_diag * diag)
{
  gsl_vector_const_view xv = gsl_vector_const_view_array(&x, 1);
  gsl_vector_view resultv = gsl_vector_view_array(result, 1);
  gsl_vector_view esterrv = gsl_vector_view_array(esterr, 1);
  gsl_vector_view fsv = gsl_vector_view_array(final_step, 1);

  return gsl_deriv_driver(params,
                          f,
                          &xv.vector,
                          &resultv.vector,
                          &esterrv.vector,
                          &fsv.vector,
                          neval,
                          diag);
}

/*
gsl_deriv_driver()
  Top-level driver for differentiation algorithms

Inputs: params     - derivative parameters
        f          - function to differentiate
        x          - point(s) at which to compute derivative, length n
        result     - (output) derivative result(s), length n
        esterr     - (output) 95% uncertainty estimate of the derivative(s), such that
                       abs(result - f'(x)) < esterr
                     length n
        final_step - (output) final overall step size chosen for each derivative, length n
        neval      - (output) total number of function evaluations
        diag       - (output) diagnostic information for adaptive algorithm (may be NULL)

Return: success/error
*/

int
gsl_deriv_driver(const gsl_deriv_parameters * params, const gsl_function * f, const gsl_vector * x,
                 gsl_vector * result, gsl_vector * esterr, gsl_vector * final_step, unsigned int * neval,
                 gsl_deriv_diag * diag)
{
  const size_t n = x->size;

  if (result->size != n)
    {
      GSL_ERROR ("result must have same length as x", GSL_EBADLEN);
    }
  else if (esterr->size != n)
    {
      GSL_ERROR ("esterr must have same length as x", GSL_EBADLEN);
    }
  else if (final_step->size != n)
    {
      GSL_ERROR ("final_step must have same length as x", GSL_EBADLEN);
    }
  else
    {
      int status;

      if (params->fixed_step <= 0.0)
        {
          /* apply adaptive algorithm */
          status = gsl_deriv_adaptive(params, f, x, result, esterr, final_step, neval, diag);
        }
      else
        {
          /* apply fixed step algorithm */
          status = gsl_deriv_fixed(params, f, x, result, esterr, final_step, neval);
        }

      return status;
    }
}

/*
gsl_deriv_adaptive()
  Implementation of adaptive differentiation

Inputs: params     - derivative parameters
        f          - function to differentiate
        x          - point(s) at which to compute derivative, length n
        result     - (output) derivative result(s), length n
        esterr     - (output) 95% uncertainty estimate of the derivative(s), such that
                       abs(result - f'(x)) < esterr
                     length n
        final_step - (output) final overall step size chosen for each derivative, length n
        neval      - (output) total number of function evaluations
        diag       - (output) diagnostic information (may be NULL)

Return: success/error
*/

int
gsl_deriv_adaptive(const gsl_deriv_parameters * params, const gsl_function * f, const gsl_vector * x,
                   gsl_vector * result, gsl_vector * esterr, gsl_vector * final_step, unsigned int * neval,
                   gsl_deriv_diag * diag)
{
  const size_t n = x->size;

  if (params->deriv_order > GSL_DERIV_MAX_DERIV_ORDER)
    {
      GSL_ERROR ("derivative order exceeds GSL_DERIV_MAX_DERIV_ORDER", GSL_EDOM);
    }
  else if (params->method_order < 1 || params->method_order > GSL_DERIV_MAX_METHOD_ORDER)
    {
      GSL_ERROR ("method order must be in [1,GSL_DERIV_MAX_METHOD_ORDER]", GSL_EDOM);
    }
  else if (params->romberg_terms > ROMBERG_MAX_TERMS)
    {
      GSL_ERROR ("romberg terms must be <= ROMBERG_MAX_TERMS", GSL_EDOM);
    }
  else if (result->size != n)
    {
      GSL_ERROR ("result must have same length as x", GSL_EBADLEN);
    }
  else if (esterr->size != n)
    {
      GSL_ERROR ("esterr must have same length as x", GSL_EBADLEN);
    }
  else if (final_step->size != n)
    {
      GSL_ERROR ("final_step must have same length as x", GSL_EBADLEN);
    }
  else
    {
      int status;
      const gsl_deriv_scale_generator * scale_gen = gsl_deriv_scale_user;
      const gsl_deriv_cand_filter * cfilter = gsl_deriv_cand_filter_trim;
      const gsl_deriv_cand_select * cselect = gsl_deriv_cand_select_minerr;
      gsl_deriv_scale_parameters scale_params;
      gsl_deriv_diag diag_local;
      gsl_deriv_diag * dptr = (diag == NULL) ? &diag_local : diag;
      gsl_vector_view unit_weights; /* Fornberg weights for deriv_order */
      romb_state_t romb_state;
      size_t i;

      *neval = 0;

      dptr->params = *params;

      /* initialize stencil and weights */
      status = deriv_init(params, &(dptr->stencil_size), dptr->stencil, dptr->fornberg_weights, &unit_weights);
      if (status)
        return status;

      /* initialize scale parameters */
      status = scale_gen->init(params, &scale_params);
      if (status)
        return status;

      dptr->scale_ratio = scale_params.ratio;
      dptr->nscale = scale_params.nscale;

      /* fill scale_multipliers array with geometric scales */
      status = gsl_deriv_scale_fill(&scale_params, dptr->scale_multipliers);
      if (status)
        return status;

      /* number of Romberg windows */
      dptr->nromb = dptr->nscale - params->romberg_terms - 1;

      /* initialize Romberg state */
      romb_init(params, scale_params.ratio, &romb_state);

      for (i = 0; i < n; ++i)
        {
          const double xi = gsl_vector_get(x, i);
          gsl_vector_view der_init = gsl_vector_view_array(dptr->der_init, dptr->nscale);
          gsl_vector_view der_romb = gsl_vector_view_array(dptr->der_romb, dptr->nromb);
          gsl_vector_view err_romb = gsl_vector_view_array(dptr->err_romb, dptr->nromb);
          gsl_vector_view scale_romb = gsl_vector_view_array(dptr->scale_multipliers, dptr->nromb);
          double final_scale;
          size_t j;

          dptr->x = xi;
          dptr->coordinate_scale = GSL_MAX(fabs(xi), 0.02);

          /*
           * approximate derivative of function for each scale multiplier, using
           * the relation:
           *
           * f_j^{(m)}(x_i) = 1/(dx_j)^m * sum_k w_k f(x_i + s_k*dx_j)
           *
           * where:
           *
           * x_i  : current point to evaluate derivative
           * m    : deriv_order
           * j    : scale
           * w_k  : Fornberg weights for deriv_order
           * s_k  : unit stencil point
           * dx_j : coordinate_scale * scale_multiplier[j]
           */

          /* compute derivative estimates at each scale */
          for (j = 0; j < scale_params.nscale; ++j)
            {
              double dx_j = dptr->coordinate_scale * dptr->scale_multipliers[j];

              status = fd_apply(f, xi, dx_j, &unit_weights.vector,
                                dptr->stencil, dptr->stencil_size, params->deriv_order,
                                &(dptr->der_init[j]), &(dptr->err_init[j]), neval);
              if (status)
                return status;
            }

          /* perform Romberg extrapolation */
          romb_extrap(&der_init.vector, &der_romb.vector, &err_romb.vector, &romb_state);

          /* filter out extremely incorrect candidates */
          (cfilter->filter)(&der_romb.vector, &(dptr->cand));

          /* select final candidate */
          (cselect->select)(&(dptr->cand), &der_romb.vector, &err_romb.vector,
                            &scale_romb.vector, &(dptr->der_select), &(dptr->err_select), &final_scale);

          dptr->step_select = dptr->coordinate_scale * final_scale;

          gsl_vector_set(result, i, dptr->der_select);
          gsl_vector_set(esterr, i, dptr->err_select);
          gsl_vector_set(final_step, i, dptr->step_select);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_deriv_fixed()
  Implementation of fixed step differentiation

Inputs: params     - derivative parameters
        f          - function to differentiate
        x          - point(s) at which to compute derivative, length n
        result     - (output) derivative result(s), length n
        esterr     - (output) absolute error estimate of the derivative(s), such that
                       abs(result - f'(x)) < esterr
                     length n
        final_step - (output) final overall step size chosen for each derivative, length n
        neval      - (output) total number of function evaluations

Return: success/error
*/

int
gsl_deriv_fixed(const gsl_deriv_parameters * params, const gsl_function * f, const gsl_vector * x,
                gsl_vector * result, gsl_vector * esterr, gsl_vector * final_step, unsigned int * neval)
{
  const size_t n = x->size;

  if (params->deriv_order > GSL_DERIV_MAX_DERIV_ORDER)
    {
      GSL_ERROR ("derivative order exceeds GSL_DERIV_MAX_DERIV_ORDER", GSL_EDOM);
    }
  else if (params->method_order < 1 || params->method_order > GSL_DERIV_MAX_METHOD_ORDER)
    {
      GSL_ERROR ("method order must be in [1,GSL_DERIV_MAX_METHOD_ORDER]", GSL_EDOM);
    }
  else if (params->fixed_step <= 0.0)
    {
      GSL_ERROR ("fixed_step must be > 0", GSL_EDOM);
    }
  else if (result->size != n)
    {
      GSL_ERROR ("result must have same length as x", GSL_EBADLEN);
    }
  else if (esterr->size != n)
    {
      GSL_ERROR ("esterr must have same length as x", GSL_EBADLEN);
    }
  else if (final_step->size != n)
    {
      GSL_ERROR ("final_step must have same length as x", GSL_EBADLEN);
    }
  else
    {
      int status;
      const double ratio = 2.0;
      double stencil[GSL_DERIV_MAX_STENCIL];
      double unit_weights_all[GSL_DERIV_MAX_STENCIL*(GSL_DERIV_MAX_DERIV_ORDER+1)];
      gsl_vector_view unit_weights; /* Fornberg weights for deriv_order */
      size_t stencil_size;
      size_t i;

      *neval = 0;

      /* initialize stencil and weights */
      status = deriv_init(params, &stencil_size, stencil, unit_weights_all, &unit_weights);
      if (status)
        return status;

      for (i = 0; i < n; ++i)
        {
          const double xi = gsl_vector_get(x, i);
          double D1, D2; /* derivative estimates */
          double round1; /* roundoff error estimate for D1 */

          /*
           * approximate derivative of function at h and h/r:
           * D1 = D(h)   (roundoff error round1)
           * D2 = D(h/r)
           *
           * truncation error = |D(h) - D(h/r)|
           * total error = trunc + round1
           */

          /* compute D(h) and roundoff error */
          status = fd_apply(f, xi, params->fixed_step, &unit_weights.vector,
                            stencil, stencil_size, params->deriv_order,
                            &D1, &round1, neval);
          if (status)
            return status;

          /* compute D(h/r) */
          status = fd_apply(f, xi, params->fixed_step / ratio, &unit_weights.vector,
                            stencil, stencil_size, params->deriv_order,
                            &D2, NULL, neval);
          if (status)
            return status;

          gsl_vector_set(result, i, D1);
          gsl_vector_set(esterr, i, fabs(D1 - D2) + round1);
          gsl_vector_set(final_step, i, params->fixed_step);
        }

      return GSL_SUCCESS;
    }
}

int
gsl_deriv_print_diag(FILE * fp, const gsl_deriv_diag * diag)
{
  const gsl_deriv_parameters * params = &(diag->params);
  gsl_matrix_const_view A = gsl_matrix_const_view_array(diag->fornberg_weights, params->deriv_order + 1, diag->stencil_size);
  gsl_vector_const_view w = gsl_matrix_const_row(&A.matrix, params->deriv_order);
  const gsl_deriv_cand_t * cand = &(diag->cand);
  size_t i;

  fprintf(fp, "# x             = %.12e\n", diag->x);
  fprintf(fp, "# deriv_order   = %u\n", params->deriv_order);
  fprintf(fp, "# method_order  = %u\n", params->method_order);
  fprintf(fp, "# romberg_terms = %u\n", params->romberg_terms);
  fprintf(fp, "# fixed_step    = %.12e\n", params->fixed_step);
  fprintf(fp, "# max_scale     = %.12e\n", params->max_scale);
  fprintf(fp, "# min_scale     = %.12e\n", params->min_scale);
  fprintf(fp, "# scale_ratio   = %.12e\n", diag->scale_ratio);
  fprintf(fp, "# nscale        = %zu\n", params->nscale);

  fprintf(fp, "# coord_scale   = %.12e\n", diag->coordinate_scale);
  fprintf(fp, "# nromb         = %zu\n", diag->nromb);

  fprintf(fp, "# stencil       = [ ");
  for (i = 0; i < diag->stencil_size; ++i)
    fprintf(fp, "%g ", diag->stencil[i]);
  fprintf(fp, "]\n");

  fprintf(fp, "# weights       = [ ");
  for (i = 0; i < diag->stencil_size; ++i)
    fprintf(fp, "%g ", gsl_vector_get(&w.vector, i));
  fprintf(fp, "]\n");

  fprintf(fp, "# der_select    = %.12e\n", diag->der_select);
  fprintf(fp, "# err_select    = %.12e\n", diag->err_select);
  fprintf(fp, "# step_select   = %.12e\n", diag->step_select);

  i = 1;
  fprintf(fp, "# Index 0\n");
  fprintf(fp, "# Field %zu: scale multiplier\n", i++);
  fprintf(fp, "# Field %zu: initial derivative estimate\n", i++);
  fprintf(fp, "# Field %zu: initial derivative roundoff error\n", i++);

  for (i = 0; i < diag->nscale; ++i)
    {
      fprintf(fp, "%20.12e %20.12e %20.12e\n",
              diag->scale_multipliers[i],
              diag->der_init[i],
              diag->err_init[i]);
    }

  i = 1;
  fprintf(fp, "\n\n# Index 1\n");
  fprintf(fp, "# Field %zu: scale multiplier\n", i++);
  fprintf(fp, "# Field %zu: Romberg derivative estimate\n", i++);
  fprintf(fp, "# Field %zu: Romberg error estimate\n", i++);

  for (i = 0; i < diag->nromb; ++i)
    {
      fprintf(fp, "%20.12e %20.12e %20.12e\n",
              diag->scale_multipliers[i],
              diag->der_romb[i],
              diag->err_romb[i]);
    }

  i = 1;
  fprintf(fp, "\n\n# Index 2\n");
  fprintf(fp, "# Field %zu: scale multiplier\n", i++);
  fprintf(fp, "# Field %zu: final candidate derivative estimate\n", i++);
  fprintf(fp, "# Field %zu: final candidate error estimate\n", i++);

  for (i = 0; i < cand->nkeep; ++i)
    {
      fprintf(fp, "%20.12e %20.12e %20.12e\n",
              diag->scale_multipliers[cand->index[i]],
              diag->der_romb[cand->index[i]],
              diag->err_romb[cand->index[i]]);
    }

  return GSL_SUCCESS;
}

/*********************************************
 * INTERNAL ROUTINES                         *
 *********************************************/

/*
deriv_init()
  Common initialization for fixed and adaptive differentiation

Inputs: params       - deriv parameters
        stencil_size - (output) stencil size depending on params
        stencil      - (output) finite difference stencil, length >= GSL_DERIV_MAX_STENCIL
        weights      - (output) Fornberg FD weights, length GSL_DERIV_MAX_STENCIL*(GSL_DERIV_MAX_DERIV_ORDER+1)
        W            - (output) vector view for the weights corresponding to deriv_order, length stencil_size

Return: success/error
*/

static int
deriv_init(const gsl_deriv_parameters * params, size_t * stencil_size, double stencil[],
           double weights[], gsl_vector_view * W)
{
  int status;
  gsl_matrix_view A;

  /* calculate stencil size */
  gsl_deriv_stencil_size(params->deriv_order, params->method_order, params->stencil_type, stencil_size);
  if (*stencil_size > GSL_DERIV_MAX_STENCIL)
    {
      GSL_ERROR("stencil_size exceeds GSL_DERIV_MAX_STENCIL", GSL_EBADLEN);
    }

  /* fill unit scale stencil */
  gsl_deriv_stencil_fill(params->stencil_type, *stencil_size, stencil);

  /* compute Fornberg finite difference weights for unit stencil around point 0 */
  status = gsl_deriv_fornberg_weights(stencil, *stencil_size, params->deriv_order, 0.0, weights);
  if (status)
    return status;

  A = gsl_matrix_view_array(weights, params->deriv_order + 1, *stencil_size);
  *W = gsl_matrix_row(&A.matrix, params->deriv_order);

  return GSL_SUCCESS;
}

/*
fd_apply()
  Apply finite difference weights to a function on a scaled
stencil to approximate the derivative of order m,

f^{(m)}(x0) = 1/(dx)^m * sum_k w_k f(x0 + s_k*dx)

Inputs: f            - function to evaluate
        x0           - point on which to apply stencil
        dx           - stencil scale factor
        weights      - Fornberg weights
        stencil      - FD stencil
        stencil_size - length of stencil
        deriv_order  - derivative order to approximate
        result       - (output) f^{(m)}(x0) approximated by finite differences,
                       where m is deriv_order
        abserr_round - (output) roundoff error estimate (may be NULL)
        neval        - (output) number of function evaluations

Return: success/error

Notes:
1) function evaluations are added to neval as a running sum

2) The roundoff error estimate is:

abserr_round = eps (1/dx)^m sum_k |w_k f(x0 + s_k*dx)|
*/

static int
fd_apply(const gsl_function * f, const double x0, const double dx,
         const gsl_vector * weights, const double stencil[], const size_t stencil_size,
         const unsigned int deriv_order, double * result, double * abserr_round, unsigned int * neval)
{
  if (weights->size != stencil_size)
    {
      GSL_ERROR ("weights does not match stencil_size", GSL_EBADLEN);
    }
  else
    {
      const double fac = gsl_pow_uint(1.0 / dx, deriv_order);
      double sum = 0.0, sumabs = 0.0, c = 0.0;
      size_t k;

      /* compute finite difference approximation using Kahan summation */
      for (k = 0; k < stencil_size; ++k)
        {
          const double wk = gsl_vector_get(weights, k);
          const double x_eval = x0 + dx * stencil[k];
          const double fk = GSL_FN_EVAL(f, x_eval);

          const double term = wk * fk;
          const double y = term - c;
          const double t = sum + y;

          c = (t - sum) - y;
          sum = t;

          /* accumulate absolute magnitude for roundoff estimate */
          sumabs += fabs(term);
        }

      *result = sum * fac;
      *neval += stencil_size;

      if (abserr_round)
        *abserr_round = GSL_DBL_EPSILON * sumabs * fac;

      return GSL_SUCCESS;
    }
}

/* A(i,j) = v(i+j-1) */
static int
vec2mat(const gsl_vector * v, gsl_matrix * A)
{
  const size_t M = A->size1;
  const size_t N = A->size2;
  size_t i, j;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        gsl_matrix_set(A, i, j, gsl_vector_get(v, i + j));
    }

  return GSL_SUCCESS;
}

/*
romb_init()
  Initialize variables for Romberg extrapolation

Inputs: params - deriv parameters
        state  - (output) initialized state

Return: success/error

Notes:
1) state variables initialized by this function are:
expon
nexpon
QR_data
T_data
sqrt_cov1
*/

static int
romb_init(const gsl_deriv_parameters * params, const double scale_ratio, romb_state_t * state)
{
  const size_t nexpon = params->romberg_terms;

  if (nexpon > ROMBERG_MAX_TERMS)
    {
      GSL_ERROR ("nexpon must be <= ROMBERG_MAX_TERMS", GSL_EDOM);
    }
  else
    {
      const size_t M = nexpon + 2;
      const size_t N = nexpon + 1;
      const double srinv = 1.0 / scale_ratio;
      gsl_matrix_view A = gsl_matrix_view_array(state->QR_data, M, N);
      gsl_matrix_view T = gsl_matrix_view_array(state->T_data, N, N);
      gsl_matrix_const_view R = gsl_matrix_const_submatrix(&A.matrix, 0, 0, N, N);
      double y_data[ROMBERG_MAX_N];
      gsl_vector_view y = gsl_vector_view_array(y_data, N);
      size_t i, j;

      state->nexpon = nexpon;

      /* store exponent powers appearing in truncation error expansion (depends on stencil) */
      if (params->stencil_type == GSL_DERIV_STENCIL_CENTRAL)
        {
          for (j = 0; j < nexpon; ++j)
            state->expon[j] = 2*(j+1) + params->method_order - 2;
        }
      else
        {
          for (j = 0; j < nexpon; ++j)
            state->expon[j] = (j+1) + params->method_order - 1;
        }

      /* construct LS matrix for fitting a polynomial to each Romberg window */

      for (i = 0; i < M; ++i)
        gsl_matrix_set(&A.matrix, i, 0, 1.0);

      for (j = 1; j < N; ++j)
        {
          const double base = gsl_pow_uint(srinv, state->expon[j-1]);
          double value = 1.0;

          for (i = 0; i < M; ++i)
            {
              gsl_matrix_set(&A.matrix, i, j, value);
              value *= base;
            }
        }

      /* compute QR decomposition of A matrix */
      gsl_linalg_QR_decomp_r(&A.matrix, &T.matrix);

      /* solve R^T y = e_1 */
      gsl_vector_set_zero(&y.vector);
      y_data[0] = 1.0;
      gsl_blas_dtrsv(CblasUpper, CblasTrans, CblasNonUnit, &R.matrix, &y.vector);

      /* cov1 = (R^T R)^{-1}]_11 */
      state->sqrt_cov1 = gsl_blas_dnrm2(&y.vector);

      return GSL_SUCCESS;
    }
}

/*
romb_extrap()
  Perform Romberg extrapolation

Inputs: der_init - initial derivative estimates
        der_romb - (output) Romberg extrapolated estimates
        errest   - (output) estimated error in extrapolated derivatives
        state    - workspace

Return: success/error

Notes:
1) romb_init() must be called prior to initialize state
*/

static int
romb_extrap(const gsl_vector * der_init, gsl_vector * der_romb, gsl_vector * errest, romb_state_t * state)
{
  const size_t nexpon = state->nexpon;
  const size_t M = nexpon + 2;
  const size_t N = nexpon + 1;
  const size_t nderiv = der_init->size;
  const size_t nromb = nderiv - N;

  if (nexpon > ROMBERG_MAX_TERMS)
    {
      GSL_ERROR ("nexpon must be <= ROMBERG_MAX_TERMS", GSL_EDOM);
    }
  else if (der_romb->size != nromb)
    {
      GSL_ERROR ("der_romb has wrong size", GSL_EBADLEN);
    }
  else if (errest->size != nromb)
    {
      GSL_ERROR ("errest has wrong size", GSL_EBADLEN);
    }
  else
    {
      const double t95 = 12.7062047361747; /* 97.5th percentile Student t-distribution with one dof */
      double rhs_data[ROMBERG_MAX_M*GSL_DERIV_MAX_SCALES];
      double coef_data[ROMBERG_MAX_M*GSL_DERIV_MAX_SCALES], work_data[ROMBERG_MAX_N*GSL_DERIV_MAX_SCALES];
      double s_data[GSL_DERIV_MAX_SCALES];
      gsl_matrix_const_view QR = gsl_matrix_const_view_array(state->QR_data, M, N);
      gsl_matrix_const_view T = gsl_matrix_const_view_array(state->T_data, N, N);
      gsl_matrix_view work = gsl_matrix_view_array(work_data, N, nromb);
      gsl_matrix_view rhs = gsl_matrix_view_array(rhs_data, M, nromb);
      gsl_matrix_view coef = gsl_matrix_view_array(coef_data, M, nromb);
      gsl_vector_view s = gsl_vector_view_array(s_data, nromb);
      gsl_vector_view v;
      size_t j;

      /* generate the Romberg windows and store in rhs matrix */
      vec2mat(der_init, &rhs.matrix);

      /* solve LS system and compute extrapolated coefficients */
      gsl_linalg_QR_lssolvem_r(&QR.matrix, &T.matrix, &rhs.matrix, &coef.matrix, &work.matrix);

      /* store extrapolated derivatives in output vector */
      v = gsl_matrix_row(&coef.matrix, 0);
      gsl_vector_memcpy(der_romb, &v.vector);

      /* compute residual norms of each LS fit */
      for (j = 0; j < nromb; ++j)
        {
          gsl_vector_const_view cj = gsl_matrix_const_subcolumn(&coef.matrix, j, N, M - N);
          s_data[j] = gsl_blas_dnrm2(&cj.vector);
        }

      /* error estimates */
      gsl_vector_axpby(t95 * state->sqrt_cov1, &s.vector, 0.0, errest);

      return GSL_SUCCESS;
    }
}
