/* deriv/test_jacobian.c
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

#include "test_j.c"

int
test_jacobian_tol(const double tol, const jacobian_test_function * TF, const gsl_vector * x,
                  const gsl_matrix * A, const gsl_deriv_stencil_t stencil_type,
                  const gsl_deriv_step_policy_t step_policy, const char * desc)
{
  const unsigned int method_order = 4;
  const size_t m = A->size1;
  const size_t n = A->size2;
  gsl_matrix * J = gsl_matrix_alloc(m, n);
  const gsl_deriv_vec_function F = { TF->f, n, m, (void *) A };
  gsl_deriv_parameters params;
  gsl_matrix * result = gsl_matrix_alloc(m, n);
  gsl_matrix * abserr = gsl_matrix_alloc(m, n);
  gsl_matrix * fs = gsl_matrix_alloc(m, n);
  gsl_deriv_jacobian_workspace * w;
  unsigned int neval;
  size_t i, j;

  /* compute Jacobian */
  (TF->df)(x, (void *) A, J);

  params = test_parameters(1, method_order, stencil_type, step_policy, -1.0, -1.0);
  w = gsl_deriv_jacobian_alloc(&params, m);
  gsl_deriv_jacobian(&F, (gsl_vector *) x, result, abserr, fs, &neval, w);

  for (i = 0; i < m; ++i)
    {
      for (j = 0; j < n; ++j)
        {
          double result_ij = gsl_matrix_get(result, i, j);
          double expected_ij = gsl_matrix_get(J, i, j);
          gsl_test_rel(result_ij, expected_ij, tol, "%s jacobian[%zu,%zu] %s", TF->desc, i, j, desc);
        }
    }

  gsl_deriv_jacobian_free(w);
  gsl_matrix_free(J);
  gsl_matrix_free(result);
  gsl_matrix_free(abserr);
  gsl_matrix_free(fs);

  return 0;
}

int
test_jacobian(const double tol, const jacobian_test_function * functions, gsl_rng * rng_p)
{
  const size_t M = 10; /* maximum F size */
  const size_t N = 10; /* maximum x size */
  const size_t K = 20; /* number of random points per function */
  const jacobian_test_function * ptr;
  size_t m, n;

  for (n = 1; n <= N; ++n)
    {
      gsl_vector * x = gsl_vector_alloc(n);

      for (m = 1; m <= M; ++m)
        {
          gsl_matrix * A = gsl_matrix_alloc(m, n);

          for (ptr = functions; ptr->f != NULL; ++ptr)
            {
              size_t k;

              if (ptr->n != n && ptr->n != 0)
                continue;
              if (ptr->m != m && ptr->m != 0)
                continue;

              for (k = 0; k < K; ++k)
                {
                  random_vector(ptr->xmin, ptr->xmax, x, rng_p);
                  random_matrix(-1.0, 1.0, A, rng_p);

                  test_jacobian_tol(tol, ptr, x, A, GSL_DERIV_STENCIL_CENTRAL, GSL_DERIV_STEP_ABSOLUTE,  "central absolute");
                  test_jacobian_tol(tol, ptr, x, A, GSL_DERIV_STENCIL_CENTRAL, GSL_DERIV_STEP_RELATIVE,  "central relative");

                  test_jacobian_tol(1.0e1*tol, ptr, x, A, GSL_DERIV_STENCIL_FORWARD, GSL_DERIV_STEP_ABSOLUTE,  "forward absolute");
                  test_jacobian_tol(1.0e1*tol, ptr, x, A, GSL_DERIV_STENCIL_FORWARD, GSL_DERIV_STEP_RELATIVE,  "forward relative");

                  test_jacobian_tol(1.0e1*tol, ptr, x, A, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_ABSOLUTE,  "backward absolute");
                  test_jacobian_tol(1.0e1*tol, ptr, x, A, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_RELATIVE,  "backward relative");
                }
            }

          gsl_matrix_free(A);
        }

      gsl_vector_free(x);
    }

  return GSL_SUCCESS;
}
