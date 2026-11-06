/* deriv/test_adaptive.c
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
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>

int
test_deriv_adaptive_diag(const deriv_test_t * test_t, const double x, const double a,
                         const unsigned int deriv_order,
                         const double min_step, const double max_step,
                         const gsl_deriv_stencil_t stencil_type,
                         const gsl_deriv_step_policy_t step_policy)
{
  const gsl_function F = { test_t->f, (void *) &a };
  gsl_deriv_parameters params = gsl_deriv_default_parameters(deriv_order, 4);
  double result, abserr, fs;
  unsigned int neval;
  gsl_deriv_diag diag;

  params.stencil_type = stencil_type;
  params.step_policy = step_policy;
  params.deriv_order = deriv_order;

  if (min_step > 0.0)
    params.min_step = min_step;
  if (max_step > 0.0)
    params.max_step = max_step;

  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, &diag);
  gsl_deriv_diag_print(stdout, &diag);

  return 0;
}

int
test_deriv_adaptive_tol(const double tol, const deriv_test_function * TF, const double x,
                        const gsl_deriv_stencil_t stencil_type,
                        const gsl_deriv_step_policy_t step_policy, const char * desc)
{
  const unsigned int method_order = 4;
  const double a = TF->a;
  const deriv_test_t * test_t = TF->test_t;
  const gsl_function F   = { test_t->f,   (void *) &a };
  const gsl_function DF  = { test_t->df,  (void *) &a };
  const gsl_function D2F = { test_t->d2f, (void *) &a };
  const gsl_function D3F = { test_t->d3f, (void *) &a };
  const gsl_function D4F = { test_t->d4f, (void *) &a };
  const double expected_df  = GSL_FN_EVAL (&DF, x);
  const double expected_d2f = GSL_FN_EVAL (&D2F, x);
  const double expected_d3f = GSL_FN_EVAL (&D3F, x);
  const double expected_d4f = GSL_FN_EVAL (&D4F, x);
  double scale = TF->scale;
  gsl_vector_view scale_view = gsl_vector_view_array(&scale, 1);
  const gsl_vector * scales = (step_policy == GSL_DERIV_STEP_USER_SCALE) ? &scale_view.vector : NULL;
  gsl_deriv_parameters params;
  double result, abserr, fs;
  unsigned int neval;
  gsl_deriv_diag diag;

  params = test_parameters(1, method_order, stencil_type, step_policy);
  params.scales = scales;
  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, &diag);
  gsl_test_rel(result, expected_df, tol, "%s df adaptive %s x=%.12e a=%g", TF->desc, desc, x, a);

  params = test_parameters(2, method_order, stencil_type, step_policy);
  params.scales = scales;
  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, &diag);
  gsl_test_rel(result, expected_d2f, 1.0e2*tol, "%s d2f adaptive %s x=%.12e a=%g", TF->desc, desc, x, a);

  params = test_parameters(3, method_order, stencil_type, step_policy);
  params.scales = scales;
  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, &diag);
  gsl_test_rel(result, expected_d3f, 1.0e3*tol, "%s d3f adaptive %s x=%.12e a=%g", TF->desc, desc, x, a);

  params = test_parameters(4, method_order, stencil_type, step_policy);
  params.scales = scales;
  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, &diag);
  gsl_test_rel(result, expected_d4f, 1.0e6*tol, "%s d4f adaptive %s x=%.12e a=%g", TF->desc, desc, x, a);

  return 0;
}

int
test_deriv_adaptive(const double tol, const deriv_test_function * functions, gsl_rng * rng_p)
{
  const size_t N = 20; /* number of random points per function */
  const deriv_test_function * ptr;

  for (ptr = functions; ptr->test_t != NULL; ++ptr)
    {
      size_t j;

      for (j = 0; j < N; ++j)
        {
          double u = gsl_rng_uniform(rng_p);
          double x = ptr->xmin + u * (ptr->xmax - ptr->xmin);

          if (ptr->scale > 0.0)
            {
              /* function with a characteristic length scale: steps scaled by that length */
              test_deriv_adaptive_tol(tol,       ptr, x, GSL_DERIV_STENCIL_CENTRAL,  GSL_DERIV_STEP_USER_SCALE, "central user scale random");
              test_deriv_adaptive_tol(1.0e3*tol, ptr, x, GSL_DERIV_STENCIL_FORWARD,  GSL_DERIV_STEP_USER_SCALE, "forward user scale random");
              test_deriv_adaptive_tol(1.0e2*tol, ptr, x, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_USER_SCALE, "backward user scale random");
            }
          else
            {
              test_deriv_adaptive_tol(tol,       ptr, x, GSL_DERIV_STENCIL_CENTRAL,  GSL_DERIV_STEP_ABSOLUTE, "central absolute random");
              test_deriv_adaptive_tol(tol,       ptr, x, GSL_DERIV_STENCIL_CENTRAL,  GSL_DERIV_STEP_RELATIVE, "central relative random");

              test_deriv_adaptive_tol(1.0e3*tol, ptr, x, GSL_DERIV_STENCIL_FORWARD,  GSL_DERIV_STEP_ABSOLUTE, "forward absolute random");
              test_deriv_adaptive_tol(1.0e3*tol, ptr, x, GSL_DERIV_STENCIL_FORWARD,  GSL_DERIV_STEP_RELATIVE, "forward relative random");

              test_deriv_adaptive_tol(1.0e2*tol, ptr, x, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_ABSOLUTE, "backward absolute random");
              test_deriv_adaptive_tol(1.0e2*tol, ptr, x, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_RELATIVE, "backward relative random");
            }
        }
    }

  return GSL_SUCCESS;
}
