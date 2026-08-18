/* deriv/jacobian.c
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

#include "deriv_common.c"

static int fd_apply_vector(const gsl_deriv_vec_function * F, gsl_vector * x0, const size_t j, const double h,
                           const gsl_vector * weights, const double stencil[], const size_t stencil_size,
                           const unsigned int deriv_order, gsl_vector * result, gsl_vector * abserr_round,
                           unsigned int * neval, gsl_vector * work2m);

/*
gsl_deriv_jacobian_alloc()
  Allocate workspace for computing Jacobians

Inputs: params - deriv parameters
        nrows  - number of outputs of F function (rows of Jacobian matrix)
        ncols  - number of inputs of F function (columns of Jacobian matrix)

Return: pointer to workspace
*/

gsl_deriv_jacobian_workspace *
gsl_deriv_jacobian_alloc(const gsl_deriv_parameters * params, const size_t nrows,
                         const size_t ncols)
{
  int status;
  gsl_deriv_jacobian_workspace * w;

  w = calloc(1, sizeof(gsl_deriv_jacobian_workspace));
  if (w == NULL)
    {
      GSL_ERROR_NULL ("unable to allocate workspace", GSL_ENOMEM);
    }

  w->params = *params;
  w->nrows = nrows;
  w->ncols = ncols;

  w->romb_state = malloc(sizeof(romb_state_t));
  if (w->romb_state == NULL)
    {
      gsl_deriv_jacobian_free(w);
      GSL_ERROR_NULL ("unable to allocate romb_state", GSL_ENOMEM);
    }

  /* initialize stencil and weights */
  status = deriv_init(params, &(w->stencil_size), w->stencil, w->fornberg_weights, &(w->unit_weights));
  if (status)
    {
      gsl_deriv_jacobian_free(w);
      GSL_ERROR_NULL("deriv_init failed", GSL_FAILURE);
    }

  w->step_ratio = pow(params->max_step / params->min_step, 1.0 / (params->nstep - 1.0));
  w->nstep = params->nstep;

  /* compute step sizes with geometric ratio */
  status = geometric_fill(params->max_step, w->step_ratio, w->nstep, w->h);
  if (status)
    {
      gsl_deriv_jacobian_free(w);
      GSL_ERROR_NULL("geometric_fill failed", GSL_FAILURE);
    }

  /* number of Romberg windows */
  w->nromb = w->nstep - params->romberg_terms - 1;

  /* initialize Romberg state */
  romb_init(params, w->step_ratio, (romb_state_t *) w->romb_state);

  w->work2l = gsl_vector_alloc(2*nrows);
  if (w->work2l == NULL)
    {
      gsl_deriv_jacobian_free(w);
      GSL_ERROR_NULL("failed to allocate work2l", GSL_ENOMEM);
    }

  w->workn = gsl_vector_alloc(ncols);
  if (w->workn == NULL)
    {
      gsl_deriv_jacobian_free(w);
      GSL_ERROR_NULL("failed to allocate workn", GSL_ENOMEM);
    }

  w->der_init = gsl_matrix_alloc(nrows, GSL_DERIV_MAX_STEPS);
  if (w->der_init == NULL)
    {
      gsl_deriv_jacobian_free(w);
      GSL_ERROR_NULL("failed to allocate der_init", GSL_ENOMEM);
    }

  w->err_init = gsl_matrix_alloc(nrows, GSL_DERIV_MAX_STEPS);
  if (w->err_init == NULL)
    {
      gsl_deriv_jacobian_free(w);
      GSL_ERROR_NULL("failed to allocate err_init", GSL_ENOMEM);
    }

  w->der_romb = gsl_vector_alloc(w->nromb);
  if (w->der_romb == NULL)
    {
      gsl_deriv_jacobian_free(w);
      GSL_ERROR_NULL("failed to allocate der_romb", GSL_ENOMEM);
    }

  w->err_romb = gsl_vector_alloc(w->nromb);
  if (w->err_romb == NULL)
    {
      gsl_deriv_jacobian_free(w);
      GSL_ERROR_NULL("failed to allocate err_romb", GSL_ENOMEM);
    }

  return w;
}

/*
gsl_deriv_jacobian_free()
  Free memory associated with 'w'
*/

void
gsl_deriv_jacobian_free(gsl_deriv_jacobian_workspace * w)
{
  if (w->romb_state)
    free(w->romb_state);

  if (w->work2l)
    gsl_vector_free(w->work2l);

  if (w->workn)
    gsl_vector_free(w->workn);

  if (w->der_init)
    gsl_matrix_free(w->der_init);

  if (w->err_init)
    gsl_matrix_free(w->err_init);

  if (w->der_romb)
    gsl_vector_free(w->der_romb);

  if (w->err_romb)
    gsl_vector_free(w->err_romb);

  free(w);
}

/*
gsl_deriv_jacobian_fixed()
  Compute Jacobian matrix of function:

F: R^n -> R^l

J_ij = dF(x)_i/dx_j

Inputs: F          - function to differentiate, returns l values
        x          - point at which to compute derivative, length n
        h          - vector of absolute step sizes, length n
        J          - (output) J(i,j) = dF(x)_i/dx_j, l-by-n
        esterr     - (output) absolute error estimates, l-by-n
        neval      - (output) total number of function evaluations
        w          - workspace

Return: success/error

Notes:
*/

int
gsl_deriv_jacobian_fixed(const gsl_deriv_vec_function * F, const gsl_vector * x,
                         const gsl_vector * h, gsl_matrix * J, gsl_matrix * esterr,
                         unsigned int * neval, gsl_deriv_jacobian_workspace * w)
{
  const size_t n = x->size;
  const size_t l = F->m;

  if (n != w->ncols)
    {
      GSL_ERROR ("x does not match workspace", GSL_EBADLEN);
    }
  else if (l != w->nrows)
    {
      GSL_ERROR ("F does not match workspace", GSL_EBADLEN);
    }
  else if (n != h->size)
    {
      GSL_ERROR ("h does not match workspace", GSL_EBADLEN);
    }
  else if (J->size1 != l || J->size2 != n)
    {
      GSL_ERROR ("J has wrong size", GSL_EBADLEN);
    }
  else if (esterr->size1 != l || esterr->size2 != n)
    {
      GSL_ERROR ("esterr must have same size as J", GSL_EBADLEN);
    }
  else
    {
      const gsl_deriv_parameters * params = &(w->params);
      const double ratio = 2.0;
      size_t j;

      *neval = 0;

      gsl_vector_memcpy(w->workn, x);

      for (j = 0; j < n; ++j)
        {
          int status;
          double hj = gsl_vector_get(h, j);
          gsl_vector_view D1 = gsl_matrix_column(J, j);
          gsl_vector_view D2 = gsl_matrix_column(w->der_init, 0);
          gsl_vector_view E1 = gsl_matrix_column(esterr, j);
          gsl_vector_view E2 = gsl_matrix_column(w->err_init, 0);
          size_t k;

          /*
           * compute j-th column of Jacobian steps hj and hj/2
           *
           * truncation error = |D(h) - D(h/2)|
           * total error = trunc + round(h)
           */

          status = fd_apply_vector(F, w->workn, j, hj, &(w->unit_weights.vector),
                                   w->stencil, w->stencil_size, params->deriv_order,
                                   &D1.vector, &E1.vector, neval, w->work2l);
          if (status)
            return status;

          status = fd_apply_vector(F, w->workn, j, hj/ratio, &(w->unit_weights.vector),
                                   w->stencil, w->stencil_size, params->deriv_order,
                                   &D2.vector, &E2.vector, neval, w->work2l);
          if (status)
            return status;

          /* add truncation error to esterr */
          for (k = 0; k < l; ++k)
            {
              double D1k = gsl_vector_get(&D1.vector, k);
              double D2k = gsl_vector_get(&D2.vector, k);
              double * ptr = gsl_vector_ptr(&E1.vector, k);
              *ptr += fabs(D1k - D2k);
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_deriv_jacobian()
  Compute Jacobian matrix of function:

F: R^n -> R^l

J_ij = dF(x)_i/dx_j

Inputs: F          - function to differentiate, returns l values
        x          - point at which to compute derivative, length n
        J          - (output) J(i,j) = dF(x)_i/dx_j, l-by-n
        esterr     - (output) 95% uncertainty estimate of the derivative, such that
                       abs(J(i,j) - dF/dx_j) < esterr(i,j)
                     l-by-n
        final_step - (output) final overall step size chosen for derivatives, l-by-n
        neval      - (output) total number of function evaluations
        w          - workspace

Return: success/error

Notes:
*/

int
gsl_deriv_jacobian(const gsl_deriv_vec_function * F, const gsl_vector * x,
                   gsl_matrix * J, gsl_matrix * esterr, gsl_matrix * final_step, unsigned int * neval,
                   gsl_deriv_jacobian_workspace * w)
{
  const size_t n = x->size;
  const size_t l = F->m;

  if (n != w->ncols)
    {
      GSL_ERROR ("x does not match workspace", GSL_EBADLEN);
    }
  else if (l != w->nrows)
    {
      GSL_ERROR ("F does not match workspace", GSL_EBADLEN);
    }
  else if (J->size1 != l || J->size2 != n)
    {
      GSL_ERROR ("J has wrong size", GSL_EBADLEN);
    }
  else if (esterr->size1 != l || esterr->size2 != n)
    {
      GSL_ERROR ("esterr must have same size as J", GSL_EBADLEN);
    }
  else if (final_step->size1 != l || final_step->size2 != n)
    {
      GSL_ERROR ("final_step must have same size as J", GSL_EBADLEN);
    }
  else
    {
      const gsl_deriv_parameters * params = &(w->params);
      const gsl_deriv_cand_filter * cfilter = gsl_deriv_cand_filter_none;
      const gsl_deriv_cand_select * cselect = gsl_deriv_cand_select_minerr;
      gsl_vector_const_view step_romb = gsl_vector_const_view_array(w->h, w->nromb);
      size_t j;

      *neval = 0;

      gsl_vector_memcpy(w->workn, x);

      for (j = 0; j < n; ++j)
        {
          int status;
          double xj = gsl_vector_get(x, j);
          double scale_j;
          size_t k;

          /* get scale factor for steps */
          status = scale_get(params, xj, j, &scale_j);
          if (status)
            return status;

          /*
           * compute approximate j-th column of Jacobian for each step,
           * using the relation:
           *
           * J_k(:,j) = 1/(h_k)^m * sum_l w_l F(x_j + s_l*h_k)
           *
           * where:
           *
           * x_j  : current point to evaluate derivative
           * m    : deriv_order
           * k    : step
           * w_l  : Fornberg weights for deriv_order
           * s_l  : unit stencil point
           * h_k  : step[k]
           * F    : vector-valued function, length l
           */

          /* compute derivative estimates at each step */
          for (k = 0; k < w->nstep; ++k)
            {
              double h_jk = scale_j * w->h[k];
              gsl_vector_view der_init = gsl_matrix_column(w->der_init, k);
              gsl_vector_view err_init = gsl_matrix_column(w->err_init, k);

              status = fd_apply_vector(F, w->workn, j, h_jk, &(w->unit_weights.vector),
                                       w->stencil, w->stencil_size, params->deriv_order,
                                       &der_init.vector, &err_init.vector, neval, w->work2l);
              if (status)
                return status;
            }

          /* loop over rows to finalize estimates for each F_k */
          for (k = 0; k < l; ++k)
            {
              gsl_vector_view der_init = gsl_matrix_subrow(w->der_init, k, 0, w->nstep);
              double der_select, err_select, step_select;
              gsl_deriv_cand_t cand;

              /* perform Romberg extrapolation on each row of der_init */
              romb_extrap(&der_init.vector, w->der_romb, w->err_romb, (romb_state_t *) w->romb_state);

              /* filter out extremely incorrect candidates */
              (cfilter->filter)(w->der_romb, &cand);

              /* select final candidate */
              (cselect->select)(&cand, w->der_romb, w->err_romb,
                                &step_romb.vector, &der_select, &err_select, &step_select);

              /* store outputs */
              gsl_matrix_set(J, k, j, der_select);
              gsl_matrix_set(esterr, k, j, err_select);
              gsl_matrix_set(final_step, k, j, step_select);
            }
        }

      return GSL_SUCCESS;
    }
}

/*********************************************
 * INTERNAL ROUTINES                         *
 *********************************************/

/*
fd_apply_vector()
  Apply finite difference weights to a vector-valued function, given a stencil
and step size, to approximate the derivative,

F: R^n -> R^m

dF(x0)/dx_j = 1/h^m * sum_k w_k F(x0(1),...,x0(j) + s_k*h,...,x0(n))

Inputs: F            - vector-valued function to evaluate
        x0           - point on which to apply stencil, length n
        j            - component of x0 for which we are computing partial derivative
        h            - step size
        weights      - Fornberg weights
        stencil      - FD stencil
        stencil_size - length of stencil
        deriv_order  - derivative order to approximate
        result       - (output) F^{(m)}(x0) approximated by finite differences,
                       where m is deriv_order, length |F|
        abserr_round - (output) roundoff error estimate, length |F|
        neval        - (output) number of function evaluations
        work2m       - workspace, length 2*m

Return: success/error

Notes:
1) Element x0(j) is modified but restored to original value on output

2) function evaluations are added to neval as a running sum

3) The roundoff error estimate is:

abserr_round = eps (1/h)^m sum_k |w_k F(x0 + s_k*h)|
*/

static int
fd_apply_vector(const gsl_deriv_vec_function * F, gsl_vector * x0, const size_t j, const double h,
                const gsl_vector * weights, const double stencil[], const size_t stencil_size,
                const unsigned int deriv_order, gsl_vector * result, gsl_vector * abserr_round,
                unsigned int * neval, gsl_vector * work2m)
{
  const size_t m = F->m;

  if (j >= x0->size)
    {
      GSL_ERROR ("index j larger than x0", GSL_EBADLEN);
    }
  else if (weights->size != stencil_size)
    {
      GSL_ERROR ("weights does not match stencil_size", GSL_EBADLEN);
    }
  else if (result->size != m)
    {
      GSL_ERROR ("result vector does not match F", GSL_EBADLEN);
    }
  else if (abserr_round->size != m)
    {
      GSL_ERROR ("abserr_round vector does not match F", GSL_EBADLEN);
    }
  else if (work2m->size != 2*m)
    {
      GSL_ERROR ("work2m vector has wrong size", GSL_EBADLEN);
    }
  else
    {
      const double x0j = gsl_vector_get(x0, j);
      const double fac = gsl_pow_uint(1.0 / h, deriv_order);
      gsl_vector_view Fv = gsl_vector_subvector(work2m, 0, m);
      gsl_vector_view comp = gsl_vector_subvector(work2m, m, m);
      size_t k;

      gsl_vector_set_zero(result);
      gsl_vector_set_zero(abserr_round);
      gsl_vector_set_zero(&comp.vector);

      /* compute finite difference approximation using Kahan summation */
      for (k = 0; k < stencil_size; ++k)
        {
          int status;
          const double wk = gsl_vector_get(weights, k);
          size_t i;

          /* perturb x0(j) with desired step */
          gsl_vector_set(x0, j, x0j + h * stencil[k]);

          /* evaluate F at modified x0 point */
          status = GSL_DERIV_VEC_FN_EVAL(F, x0, &Fv.vector);
          if (status)
            return status;

          for (i = 0; i < m; ++i)
            {
              double * sum_ptr = gsl_vector_ptr(result, i);
              double * comp_ptr = gsl_vector_ptr(&comp.vector, i);
              double * err_ptr = gsl_vector_ptr(abserr_round, i);
              const double Fi = gsl_vector_get(&Fv.vector, i);
              const double term = wk * Fi;
              const double z = term - (*comp_ptr);
              const double t = *sum_ptr + z;

              /* accumulate roundoff error estimate */
              *err_ptr += fabs(term);

              *comp_ptr = (t - (*sum_ptr)) - z;
              *sum_ptr = t;
            }
        }

      /* finalize result */
      gsl_vector_scale(result, fac);

      /* finalize roundoff error estimate */
      gsl_vector_scale(abserr_round, GSL_DBL_EPSILON * fac);

      *neval += stencil_size;

      /* restore x0(j) */
      gsl_vector_set(x0, j, x0j);

      return GSL_SUCCESS;
    }
}
