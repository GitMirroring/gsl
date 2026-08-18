/* deriv/gradient.c
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

typedef struct
{
  gsl_vector * x;               /* derivative point */
  const gsl_deriv_function * f; /* original vector-valued function */
  size_t idx;                   /* index for current calculation */
} grad_parameters;

static double grad_func(double x, void * vparams);

/*
gsl_deriv_gradient_fixed()
  Compute vector of partial derivatives of a function

Inputs: params     - derivative parameters
        f          - function to differentiate
        x          - point at which to compute derivative, length n
        h          - vector of absolute fixed steps sizes
        result     - (output) result(i) = d^m f/dx_i^m, length n
        esterr     - (output) absolute error estimate
                       abs(result - f'(x)) < esterr
                     length n
        neval      - (output) total number of function evaluations
        work       - workspace, length n

Return: success/error

Notes:
*/

int
gsl_deriv_gradient_fixed(const gsl_deriv_parameters * params, const gsl_deriv_function * f,
                         const gsl_vector * x, const gsl_vector * h,
                         gsl_vector * result, gsl_vector * esterr, unsigned int * neval,
                         gsl_vector * work)
{
  const size_t n = x->size;

  if (n != h->size)
    {
      GSL_ERROR ("h must have same length as x", GSL_EBADLEN);
    }
  else if (n != result->size)
    {
      GSL_ERROR ("result must have same length as x", GSL_EBADLEN);
    }
  else if (n != esterr->size)
    {
      GSL_ERROR ("esterr must have same length as x", GSL_EBADLEN);
    }
  else if (n != work->size)
    {
      GSL_ERROR ("work must have same length as x", GSL_EBADLEN);
    }
  else
    {
      gsl_function F;
      grad_parameters grad_params;
      size_t i;

      gsl_vector_memcpy(work, x);

      grad_params.x = work;
      grad_params.f = f;

      F.function = grad_func;
      F.params = &grad_params;

      for (i = 0; i < n; ++i)
        {
          int status;
          double xi = gsl_vector_get(x, i);
          double hi = gsl_vector_get(h, i);
          double out_i, err_i;

          grad_params.idx = i;

          /* compute df/dx_i, holding other parameters fixed */
          status = gsl_deriv_fixed(params, &F, xi, hi, &out_i, &err_i, neval);
          if (status)
            return status;

          /* restore x(i) */
          gsl_vector_set(work, i, xi);

          /* store outputs */
          gsl_vector_set(result, i, out_i);
          gsl_vector_set(esterr, i, err_i);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_deriv_gradient()
  Compute vector of partial derivatives of a function

Inputs: params     - derivative parameters
        f          - function to differentiate
        x          - point at which to compute derivative, length n
        result     - (output) result(i) = d^m f/dx_i^m, length n
        esterr     - (output) 95% uncertainty estimate of the derivative, such that
                       abs(result - f'(x)) < esterr
                     length n
        final_step - (output) final overall step size chosen for derivative, length n
        neval      - (output) total number of function evaluations
        work       - workspace, length n

Return: success/error

Notes:
*/

int
gsl_deriv_gradient(const gsl_deriv_parameters * params, const gsl_deriv_function * f,
                   const gsl_vector * x, gsl_vector * result, gsl_vector * esterr,
                   gsl_vector * final_step, unsigned int * neval, gsl_vector * work)
{
  const size_t n = x->size;

  if (n != result->size)
    {
      GSL_ERROR ("result must have same length as x", GSL_EBADLEN);
    }
  else if (n != esterr->size)
    {
      GSL_ERROR ("esterr must have same length as x", GSL_EBADLEN);
    }
  else if (n != final_step->size)
    {
      GSL_ERROR ("final_step must have same length as x", GSL_EBADLEN);
    }
  else if (n != work->size)
    {
      GSL_ERROR ("work must have same length as x", GSL_EBADLEN);
    }
  else
    {
      gsl_function F;
      grad_parameters grad_params;
      size_t i;

      gsl_vector_memcpy(work, x);

      grad_params.x = work;
      grad_params.f = f;

      F.function = grad_func;
      F.params = &grad_params;

      for (i = 0; i < n; ++i)
        {
          int status;
          double xi = gsl_vector_get(x, i);
          double out_i, err_i, step_i;
          gsl_deriv_diag diag;

          grad_params.idx = i;

          /* compute df/dx_i, holding other parameters fixed */
          status = gsl_deriv_adaptive(params, &F, xi, &out_i, &err_i, &step_i, neval, &diag);
          if (status)
            return status;

          /* restore x(i) */
          gsl_vector_set(work, i, xi);

          /* store outputs */
          gsl_vector_set(result, i, out_i);
          gsl_vector_set(esterr, i, err_i);
          gsl_vector_set(final_step, i, step_i);
        }

      return GSL_SUCCESS;
    }
}

/* compute f(x1,...,x(idx)=x_idx,...,xn) */
static double
grad_func(double x_idx, void * vparams)
{
  const grad_parameters * params = (const grad_parameters *) vparams;

  gsl_vector_set(params->x, params->idx, x_idx);

  /* evaluate f with x(idx) = x_idx */
  return GSL_DERIV_FN_EVAL(params->f, params->x);
}
