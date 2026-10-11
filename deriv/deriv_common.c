/* deriv/deriv_common.c
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
  double g[ROMBERG_MAX_M];                     /* extrapolation weights: der_romb = sum_k g_k der_init_k, length M */
} romb_state_t;

static int geometric_fill(const double max_val, const double ratio, const size_t n, double result[]);
static int deriv_init(const gsl_deriv_parameters * params, size_t * stencil_size, double stencil[],
                      double weights[], gsl_vector_view * W);
static int vec2mat(const gsl_vector * v, gsl_matrix * A);
static int romb_init(const gsl_deriv_parameters * params, const double step_ratio, romb_state_t * state);
static int romb_extrap(const gsl_vector * der_init, const gsl_vector * err_init, gsl_vector * der_romb,
                       gsl_vector * errest, romb_state_t * state);
static int scale_get(const gsl_deriv_parameters * params, const double x, const size_t idx, double * scale);

/*
geometric_fill()
  Fill array with geometric series

result[i] = max_val / ratio^i

Inputs: max_val - maximum value in series
        ratio   - ratio between sequential values
        n       - length of series
        result  - (output) series array, length n
*/

static int
geometric_fill(const double max_val, const double ratio, const size_t n, double result[])
{
  double factor = 1.0; /* ratio^{-i} */
  size_t i;

  for (i = 0; i < n; ++i)
    {
      result[i] = max_val * factor;
      factor /= ratio;
    }

  return GSL_SUCCESS;
}

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

  /* fill unit stencil */
  gsl_deriv_stencil_fill(params->stencil_type, *stencil_size, stencil);

  /* compute Fornberg finite difference weights for unit stencil around point 0 */
  status = gsl_deriv_fornberg_weights(stencil, *stencil_size, params->deriv_order, 0.0, weights);
  if (status)
    return status;

  A = gsl_matrix_view_array(weights, params->deriv_order + 1, *stencil_size);
  *W = gsl_matrix_row(&A.matrix, params->deriv_order);

  return GSL_SUCCESS;
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

Inputs: params     - deriv parameters
        step_ratio - step ratio
        state      - (output) initialized state

Return: success/error

Notes:
1) state variables initialized by this function are:
expon
nexpon
QR_data
T_data
sqrt_cov1
g
*/

static int
romb_init(const gsl_deriv_parameters * params, const double step_ratio, romb_state_t * state)
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
      const double srinv = 1.0 / step_ratio;
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

      /*
       * extrapolation weights g = e_1^T (A^T A)^{-1} A^T: the extrapolated derivative
       * of a window is sum_k g_k der_init_k, so its rounding error is at most
       * sum_k |g_k| err_init_k
       */
      {
        double b_data[ROMBERG_MAX_M], x_data[ROMBERG_MAX_M], work_data[ROMBERG_MAX_N];
        gsl_vector_view b = gsl_vector_view_array(b_data, M);
        gsl_vector_view x = gsl_vector_view_array(x_data, M);
        gsl_vector_view work = gsl_vector_view_array(work_data, N);

        for (i = 0; i < M; ++i)
          {
            gsl_vector_set_basis(&b.vector, i);
            gsl_linalg_QR_lssolve_r(&A.matrix, &T.matrix, &b.vector, &x.vector, &work.vector);
            state->g[i] = x_data[0];
          }
      }

      return GSL_SUCCESS;
    }
}

/*
romb_extrap()
  Perform Romberg extrapolation

Inputs: der_init - initial derivative estimates
        err_init - rounding error estimates of der_init (may be NULL)
        der_romb - (output) Romberg extrapolated estimates
        errest   - (output) estimated error in extrapolated derivatives
        state    - workspace

Return: success/error

Notes:
1) romb_init() must be called prior to initialize state

2) The error estimate of each Romberg window is the largest of:

   a) the 95% uncertainty from the residual of the least-squares fit (as in
      DERIVEST); the fit has a single degree of freedom, so the residual can be
      arbitrarily small by chance
   b) the rounding error propagated through the extrapolation weights,
      sum_k |g_k| err_init_k (if err_init is provided)
   c) the smaller of the differences with the extrapolated estimates of the two
      neighboring windows, which share all but one finite difference estimate

   Without b) and c), selecting the candidate with the smallest error estimate
   favors windows in the rounding-dominated region whose estimates happen to fit
   the truncation error model
*/

static int
romb_extrap(const gsl_vector * der_init, const gsl_vector * err_init, gsl_vector * der_romb,
            gsl_vector * errest, romb_state_t * state)
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
      double rhs_data[ROMBERG_MAX_M*GSL_DERIV_MAX_STEPS];
      double coef_data[ROMBERG_MAX_M*GSL_DERIV_MAX_STEPS], work_data[ROMBERG_MAX_N*GSL_DERIV_MAX_STEPS];
      double s_data[GSL_DERIV_MAX_STEPS];
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

      for (j = 0; j < nromb; ++j)
        {
          double e = gsl_vector_get(errest, j);
          const double dj = gsl_vector_get(der_romb, j);

          /* rounding error propagated through the extrapolation */
          if (err_init != NULL)
            {
              double eround = 0.0;
              size_t k;

              for (k = 0; k < M; ++k)
                eround += fabs(state->g[k]) * gsl_vector_get(err_init, j + k);

              e = GSL_MAX(e, eround);
            }

          /*
           * consistency with the closer of the neighboring windows: a window at
           * the edge of the plateau agrees with its neighbor inside the plateau,
           * while an isolated window in the noise region disagrees with both
           */
          if (nromb > 1)
            {
              double dmin = GSL_POSINF;

              if (j > 0)
                dmin = GSL_MIN(dmin, fabs(dj - gsl_vector_get(der_romb, j - 1)));
              if (j + 1 < nromb)
                dmin = GSL_MIN(dmin, fabs(dj - gsl_vector_get(der_romb, j + 1)));

              e = GSL_MAX(e, dmin);
            }

          gsl_vector_set(errest, j, e);
        }

      return GSL_SUCCESS;
    }
}

/*
scale_get()
  Return scale factor for adaptive step sizes

Inputs: params - deriv parameters
        x      - point at which we are computing deriviatve
        idx    - index of 'x' in the input vector
        scale  - (output) scale factor

Return: success/error
*/

static int
scale_get(const gsl_deriv_parameters * params, const double x, const size_t idx, double * scale)
{
  if (params->step_policy == GSL_DERIV_STEP_RELATIVE)
    {
      *scale = GSL_MAX(fabs(x), 1.0);
    }
  else if (params->step_policy == GSL_DERIV_STEP_USER_SCALE)
    {
      if (params->scales == NULL)
        {
          GSL_ERROR ("user-supplied scales vector is NULL", GSL_EDOM);
        }

      if (idx >= params->scales->size)
        {
          GSL_ERROR ("idx does not exist in user-supplied scales vector", GSL_EBADLEN);
        }

      *scale = gsl_vector_get(params->scales, idx);
    }
  else
    {
      *scale = 1.0; /* absolute step policy */
    }

  return GSL_SUCCESS;
}
