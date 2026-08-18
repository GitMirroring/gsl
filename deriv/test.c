/* deriv/test.c
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

#define U(x) ((void) (x))              /* avoid unused parameter warnings */

/* generate random vector with elements in [a,b] */
static void
random_vector(const double a, const double b, gsl_vector * v, gsl_rng * r)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
    {
      double vi = (b - a) * gsl_rng_uniform(r) + a;
      gsl_vector_set(v, i, vi);
    }
}

/* generate random matrix with elements in [a,b] */
static void
random_matrix(const double a, const double b, gsl_matrix * C, gsl_rng * r)
{
  const size_t M = C->size1;
  const size_t N = C->size2;
  size_t i, j;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double Cij = (b - a) * gsl_rng_uniform(r) + a;
          gsl_matrix_set(C, i, j, Cij);
        }
    }
}

static double
test_h0(const unsigned int deriv_order, const unsigned int method_order)
{
  return pow(GSL_DBL_EPSILON, 1.0 / (double) (deriv_order + method_order));
}

gsl_deriv_parameters
test_parameters(const unsigned int deriv_order, const unsigned int method_order,
                const gsl_deriv_stencil_t stencil_type, const gsl_deriv_step_policy_t step_policy,
                const double min_step_1, const double max_step_1)
{
  /* this factor scales the input step range (tuned for deriv_order=1) to higher deriv orders */
  const double step_fac = pow(GSL_DBL_EPSILON, 1.0/(deriv_order+method_order) - 1.0/(method_order+1.0));
  gsl_deriv_parameters params = gsl_deriv_default_parameters(deriv_order, method_order);

  params.stencil_type = stencil_type;
  params.step_policy = step_policy;

  if (min_step_1 > 0.0)
    params.min_step = min_step_1 * step_fac;
  if (max_step_1 > 0.0)
    params.max_step = max_step_1 * step_fac;

  return params;
}

#include "test_f.c"
#include "test_adaptive.c"
#include "test_directional.c"
#include "test_fixed.c"
#include "test_fornberg.c"
#include "test_gradient.c"
#include "test_jacobian.c"
#include "test_old.c"

int
main ()
{
  const double tol = 1.0e-9;
  gsl_rng * rng_p = gsl_rng_alloc(gsl_rng_default);

  gsl_ieee_env_setup ();

#if 0
  /* diagnostic output */
  test_deriv_adaptive_diag(&sina_test,
                           2.906189932906e+00,
                           1.0,
                           1,
                           1.0e-3,
                           1.0e1,
                           GSL_DERIV_STENCIL_BACKWARD,
                           GSL_DERIV_STEP_ABSOLUTE);
  exit(1);
#endif

  /* test Fornberg weights */
  test_fornberg(tol, rng_p);

  /* test fixed step algorithm */
  test_deriv_fixed(1.0e1*tol, test_normal_functions, rng_p);

  /* test adaptive algorithm */
  test_deriv_adaptive(tol, test_normal_functions, rng_p);
  test_deriv_adaptive(1.0e3*GSL_DBL_EPSILON, test_poly_functions, rng_p); /* tight tolerance for polynomials */
  test_deriv_adaptive(tol, test_scale_functions, rng_p);

  /* test directional derivative */
  test_directional(1.0e2*tol, test_directional_functions, 0, rng_p); /* adaptive */
  test_directional(1.0e5*tol, test_directional_functions, 1, rng_p); /* fixed */

  /* test gradient */
  test_gradient(1.0e1*tol, test_gradient_functions, 0, rng_p); /* adaptive */
  test_gradient(1.0e2*tol, test_gradient_functions, 1, rng_p); /* fixed */

  /* test Jacobian */
  test_jacobian(1.0e2*tol, test_jacobian_functions, 0, rng_p); /* adaptive */
  test_jacobian(1.0e2*tol, test_jacobian_functions, 1, rng_p); /* fixed */

  /* test old/legacy routines */
  test_old();

  gsl_rng_free(rng_p);

  exit (gsl_test_summary ());
}


