/* deriv/test.c
 * 
 * Copyright (C) 2000 David Morrison
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

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

gsl_deriv_parameters
test_parameters(const unsigned int deriv_order, const unsigned int method_order,
                const gsl_deriv_stencil_t stencil_type, const gsl_deriv_step_policy_t step_policy)
{
  gsl_deriv_parameters params = gsl_deriv_default_parameters(deriv_order, method_order);

  params.stencil_type = stencil_type;
  params.step_policy = step_policy;

  return params;
}

#include "test_f.c"
#include "test_adaptive.c"
#include "test_fixed.c"
#include "test_fornberg.c"

typedef int (deriv_fn) (const gsl_function * f, double x, double h, double * res, double *abserr);

void
test (deriv_fn * deriv, gsl_function * f, gsl_function * df, double x, 
      const char * desc)
{
  double result, abserr;
  double expected = GSL_FN_EVAL (df, x);
  (*deriv) (f, x, 1e-4, &result, &abserr);

  gsl_test_abs (result, expected, GSL_MIN(1e-4,fabs(expected)) + GSL_DBL_EPSILON, desc);

  if (abserr < fabs(result-expected)) 
    {
      gsl_test_factor (abserr, fabs(result-expected), 2, "%s error estimate", desc);
    }
  else if (result == expected || expected == 0.0)
    {
      gsl_test_abs (abserr, 0.0, 1e-6, "%s abserr", desc);
    }
  else
    {
      double d = fabs(result - expected);
      gsl_test_abs (abserr, fabs(result-expected), 1e6*d, "%s abserr", desc);
    }
}

int
main ()
{
  const double tol = 1.0e-9;
  gsl_rng * rng_p = gsl_rng_alloc(gsl_rng_default);

  gsl_ieee_env_setup ();

#if 0
  test_deriv_adaptive_diag(&tanh_test,
                           -4.997505021747e+00,
                           0.0,\
                           4,
                           -1.0,\
                           -1.0,
                           GSL_DERIV_STENCIL_BACKWARD,
                           GSL_DERIV_STEP_RELATIVE);
  exit(1);
#endif

  test_fornberg(tol, rng_p);
#if 0
  test_deriv_fixed(1.0e-5, test_normal_functions, rng_p);
#endif
#if 1
  test_deriv_adaptive(tol, test_normal_functions, rng_p);

  /* polynomial tests can use a much tighter tolerance */
  test_deriv_adaptive(1.0e3*GSL_DBL_EPSILON, test_poly_functions, rng_p);
#endif

#if 0
  F2.function = &f2;
  DF2.function = &df2;

  F3.function = &f3;
  DF3.function = &df3;

  F5.function = &f5;
  DF5.function = &df5;

  F6.function = &f6;
  DF6.function = &df6;

  test_de (&F1, &DF1, 1.0, "exp(x), x=1, derivest");
  /*test_de (&F1, &DF1, 2.0, "exp(x), x=2, derivest");*/

#if 0/*XXX*/
  test (&gsl_deriv_central, &F1, &DF1, 1.0, "exp(x), x=1, central deriv");
  test (&gsl_deriv_forward, &F1, &DF1, 1.0, "exp(x), x=1, forward deriv");
  test (&gsl_deriv_backward, &F1, &DF1, 1.0, "exp(x), x=1, backward deriv");

  test (&gsl_deriv_central, &F2, &DF2, 0.1, "x^(3/2), x=0.1, central deriv");
  test (&gsl_deriv_forward, &F2, &DF2, 0.1, "x^(3/2), x=0.1, forward deriv");
  test (&gsl_deriv_backward, &F2, &DF2, 0.1, "x^(3/2), x=0.1, backward deriv");

  test (&gsl_deriv_central, &F3, &DF3, 0.45, "sin(1/x), x=0.45, central deriv");
  test (&gsl_deriv_forward, &F3, &DF3, 0.45, "sin(1/x), x=0.45, forward deriv");
  test (&gsl_deriv_backward, &F3, &DF3, 0.45, "sin(1/x), x=0.45, backward deriv");

  test (&gsl_deriv_central, &F4, &DF4, 0.5, "exp(-x^2), x=0.5, central deriv");
  test (&gsl_deriv_forward, &F4, &DF4, 0.5, "exp(-x^2), x=0.5, forward deriv");
  test (&gsl_deriv_backward, &F4, &DF4, 0.5, "exp(-x^2), x=0.5, backward deriv");

  test (&gsl_deriv_central, &F5, &DF5, 0.0, "x^2, x=0, central deriv");
  test (&gsl_deriv_forward, &F5, &DF5, 0.0, "x^2, x=0, forward deriv");
  test (&gsl_deriv_backward, &F5, &DF5, 0.0, "x^2, x=0, backward deriv");

  test (&gsl_deriv_central, &F6, &DF6, 10.0, "1/x, x=10, central deriv");
  test (&gsl_deriv_forward, &F6, &DF6, 10.0, "1/x, x=10, forward deriv");
  test (&gsl_deriv_backward, &F6, &DF6, 10.0, "1/x, x=10, backward deriv");
#endif
#endif

  gsl_rng_free(rng_p);

  exit (gsl_test_summary ());
}


