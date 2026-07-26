/* deriv/directionaldiff.c
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
  const gsl_vector * x;         /* derivative point */
  const gsl_vector * v;         /* direction for derivative */
  const gsl_deriv_function * f; /* original vector-valued function */
  gsl_vector * work;            /* workspace */
} dd_parameters;

static double dd_func(double t, void * vparams);

/*
gsl_deriv_directional()
  Compute directional derivative of a function

Inputs: params     - derivative parameters
        f          - function to differentiate
        x          - point at which to compute derivative, length n
        v          - direction for computing derivative, length n
        result     - (output) derivative result
        esterr     - (output) 95% uncertainty estimate of the derivative, such that
                       abs(result - f'(x)) < esterr
        final_step - (output) final overall step size chosen for derivative
        neval      - (output) total number of function evaluations
        diag       - (output) diagnostic information (may be NULL)
        work       - workspace, length n

Return: success/error
*/

int
gsl_deriv_directional(const gsl_deriv_parameters * params, const gsl_deriv_function * f, const gsl_vector * x,
                      const gsl_vector * v, double * result, double * esterr, double * final_step, unsigned int * neval,
                      gsl_deriv_diag * diag, gsl_vector * work)
{
  const size_t n = x->size;

  if (n != v->size)
    {
      GSL_ERROR ("v must have same length as x", GSL_EBADLEN);
    }
  else if (n != work->size)
    {
      GSL_ERROR ("work must have same length as x", GSL_EBADLEN);
    }
  else
    {
      int status;
      gsl_function F;
      dd_parameters dd_params;

      dd_params.x = x;
      dd_params.v = v;
      dd_params.f = f;
      dd_params.work = work;

      F.function = dd_func;
      F.params = &dd_params;

      status = gsl_deriv_adaptive(params, &F, 0.0, result, esterr, final_step, neval, diag);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

/* compute f(x + t*v) */
static double
dd_func(double t, void * vparams)
{
  const dd_parameters * params = (const dd_parameters *) vparams;

  /* work := t*v */
  gsl_vector_axpby(t, params->v, 0.0, params->work);

  /* work := x + t*v */
  gsl_vector_add(params->work, params->x);

  /* return f(x + t*v) */
  return GSL_DERIV_FN_EVAL(params->f, params->work);
}
