/* deriv/test_directional.c
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

#include "test_d.c"

int
test_directional_tol(const double tol, const directional_test_function * TF, const gsl_vector * x,
                     const gsl_vector * v, const gsl_deriv_stencil_t stencil_type,
                     const gsl_deriv_step_policy_t step_policy, const char * desc)
{
  const unsigned int method_order = 4;
  const size_t n = x->size;
  const gsl_deriv_function F = { TF->f, n, NULL };
  const double expected_df   = (TF->df)(x, v, NULL);
  const double expected_d2f  = (TF->d2f)(x, v, NULL);
  const double expected_d3f  = (TF->d3f)(x, v, NULL);
  gsl_deriv_parameters params;
  double result, abserr, fs;
  unsigned int neval;
  gsl_vector * work = gsl_vector_alloc(n);
  gsl_deriv_diag diag;

  params = test_parameters(1, method_order, stencil_type, step_policy, -1.0, -1.0);
  gsl_deriv_directional(&params, &F, x, v, &result, &abserr, &fs, &neval, &diag, work);
  gsl_test_rel(result, expected_df, tol, "%s df directional %s", TF->desc, desc);

  params = test_parameters(2, method_order, stencil_type, step_policy, -1.0, -1.0);
  gsl_deriv_directional(&params, &F, x, v, &result, &abserr, &fs, &neval, &diag, work);
  gsl_test_rel(result, expected_d2f, tol, "%s d2f directional %s", TF->desc, desc);

  params = test_parameters(3, method_order, stencil_type, step_policy, -1.0, -1.0);
  gsl_deriv_directional(&params, &F, x, v, &result, &abserr, &fs, &neval, &diag, work);
  gsl_test_rel(result, expected_d3f, tol, "%s d3f directional %s", TF->desc, desc);

  gsl_vector_free(work);

  return 0;
}

int
test_directional(const double tol, const directional_test_function * functions, gsl_rng * rng_p)
{
  const size_t M = 10; /* maximum vector size */
  const size_t N = 20; /* number of random points per function */
  const directional_test_function * ptr;
  size_t m;

  for (m = 1; m <= M; ++m)
    {
      gsl_vector * x = gsl_vector_alloc(m);
      gsl_vector * v = gsl_vector_alloc(m);

      for (ptr = functions; ptr->f != NULL; ++ptr)
        {
          size_t j;

          if (ptr->n != m && ptr->n != 0)
            continue;

          for (j = 0; j < N; ++j)
            {
              random_vector(ptr->xmin, ptr->xmax, x, rng_p);
              random_vector(-1.0, 1.0, v, rng_p);

              test_directional_tol(tol, ptr, x, v, GSL_DERIV_STENCIL_CENTRAL, GSL_DERIV_STEP_ABSOLUTE,  "central absolute");
              test_directional_tol(tol, ptr, x, v, GSL_DERIV_STENCIL_CENTRAL, GSL_DERIV_STEP_RELATIVE,  "central relative");

              test_directional_tol(1.0e2*tol, ptr, x, v, GSL_DERIV_STENCIL_FORWARD, GSL_DERIV_STEP_ABSOLUTE,  "forward absolute");
              test_directional_tol(1.0e2*tol, ptr, x, v, GSL_DERIV_STENCIL_FORWARD, GSL_DERIV_STEP_RELATIVE,  "forward relative");

              test_directional_tol(1.0e2*tol, ptr, x, v, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_ABSOLUTE, "backward absolute");
              test_directional_tol(1.0e2*tol, ptr, x, v, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_RELATIVE, "backward relative");
            }
        }

      gsl_vector_free(x);
      gsl_vector_free(v);
    }

  return GSL_SUCCESS;
}
