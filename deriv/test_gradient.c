/* deriv/test_gradient.c
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_test.h>

#include "test_g.c"

int
test_gradient_tol(const double tol, const gradient_test_function * TF, const int fixed_step,
                  const gsl_vector * x, const gsl_matrix * A,
                  const gsl_deriv_stencil_t stencil_type,
                  const gsl_deriv_step_policy_t step_policy, const char * desc)
{
  const unsigned int method_order = 4;
  const size_t n = x->size;
  gsl_vector * g = gsl_vector_alloc(n);
  const gsl_deriv_function F = { TF->f, n, (void *) A };
  gsl_deriv_parameters params;
  gsl_vector * result = gsl_vector_alloc(n);
  gsl_vector * abserr = gsl_vector_alloc(n);
  gsl_vector * fs = gsl_vector_alloc(n);
  gsl_vector * work = gsl_vector_alloc(n);
  unsigned int neval;
  size_t i;

  /* compute gradient */
  (TF->df)(x, (void *) A, g);

  params = test_parameters(1, method_order, stencil_type, step_policy, -1.0, -1.0);

  if (fixed_step)
    {
      gsl_vector * h = gsl_vector_alloc(n);
      double h0 = test_h0(1, method_order);

      gsl_vector_set_all(h, h0);
      gsl_deriv_gradient_fixed(&params, &F, x, h, result, abserr, &neval, work);

      gsl_vector_free(h);
    }
  else
    gsl_deriv_gradient(&params, &F, x, result, abserr, fs, &neval, work);

  for (i = 0; i < n; ++i)
    {
      double result_i = gsl_vector_get(result, i);
      double expected_i = gsl_vector_get(g, i);

      gsl_test_rel(result_i, expected_i, tol, "%s gradient[%zu,%zu] %s fixed=%d", TF->desc, n, i, desc, fixed_step);
    }

  gsl_vector_free(g);
  gsl_vector_free(result);
  gsl_vector_free(abserr);
  gsl_vector_free(fs);
  gsl_vector_free(work);

  return 0;
}

int
test_gradient(const double tol, const gradient_test_function * functions,
              const int fixed_step, gsl_rng * rng_p)
{
  const size_t M = 10; /* maximum vector size */
  const size_t N = 20; /* number of random points per function */
  const gradient_test_function * ptr;
  size_t m;

  for (m = 1; m <= M; ++m)
    {
      gsl_vector * x = gsl_vector_alloc(m);
      gsl_matrix * A = gsl_matrix_alloc(m, m);

      for (ptr = functions; ptr->f != NULL; ++ptr)
        {
          size_t j;

          if (ptr->n != m && ptr->n != 0)
            continue;

          for (j = 0; j < N; ++j)
            {
              random_vector(ptr->xmin, ptr->xmax, x, rng_p);
              random_matrix(-1.0, 1.0, A, rng_p);

              test_gradient_tol(tol, ptr, fixed_step, x, A, GSL_DERIV_STENCIL_CENTRAL, GSL_DERIV_STEP_ABSOLUTE,  "central absolute");
              test_gradient_tol(tol, ptr, fixed_step, x, A, GSL_DERIV_STENCIL_CENTRAL, GSL_DERIV_STEP_RELATIVE,  "central relative");

              test_gradient_tol(tol, ptr, fixed_step, x, A, GSL_DERIV_STENCIL_FORWARD, GSL_DERIV_STEP_ABSOLUTE,  "forward absolute");
              test_gradient_tol(tol, ptr, fixed_step, x, A, GSL_DERIV_STENCIL_FORWARD, GSL_DERIV_STEP_RELATIVE,  "forward relative");

              test_gradient_tol(tol, ptr, fixed_step, x, A, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_ABSOLUTE, "backward absolute");
              test_gradient_tol(tol, ptr, fixed_step, x, A, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_RELATIVE, "backward relative");
            }
        }

      gsl_vector_free(x);
      gsl_matrix_free(A);
    }

  return GSL_SUCCESS;
}
