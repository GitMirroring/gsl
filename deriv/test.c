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

#include "test_f.c"
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
test_deriv_adaptive_diag(const deriv_test_t * test_t, const double x, const double a,
                         const unsigned int deriv_order,
                         const double min_step, const double max_step,
                         const gsl_deriv_stencil_t stencil_type)
{
  const gsl_function F = { test_t->f, (void *) &a };
  gsl_deriv_parameters params = gsl_deriv_default_parameters();
  double result, abserr, fs;
  unsigned int neval;
  gsl_deriv_diag diag;

  params.stencil_type = stencil_type;
  params.deriv_order = deriv_order;

  if (min_step > 0.0)
    params.min_step = min_step;
  if (max_step > 0.0)
    params.max_step = max_step;

  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, &diag);
  gsl_deriv_print_diag(stdout, &diag);

  return 0;
}

int
test_deriv_adaptive_tol(const double tol, const deriv_test_function * TF, const double x,
                        const gsl_deriv_stencil_t stencil_type,
                        const gsl_deriv_step_policy_t step_policy, const char * desc)
{
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
  gsl_deriv_parameters params = gsl_deriv_default_parameters();
  double result, abserr, fs;
  unsigned int neval;
  gsl_deriv_diag diag;

  params.stencil_type = stencil_type;
  params.step_policy = step_policy;
  if (TF->min_step > 0.0)
    params.min_step = TF->min_step;
  if (TF->max_step > 0.0)
    params.max_step = TF->max_step;

  params.deriv_order = 1;
  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, &diag);
  gsl_test_rel(result, expected_df, tol, "%s df adaptive %s x=%g a=%g", TF->desc, desc, x, a);

  params.deriv_order = 2;
  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, &diag);
  gsl_test_rel(result, expected_d2f, 1.0e2*tol, "%s d2f adaptive %s x=%g a=%g", TF->desc, desc, x, a);

  params.deriv_order = 3;
  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, &diag);
  gsl_test_rel(result, expected_d3f, 1.0e3*tol, "%s d3f adaptive %s x=%g a=%g", TF->desc, desc, x, a);

  params.deriv_order = 4;
  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, &diag);
  gsl_test_rel(result, expected_d4f, 1.0e6*tol, "%s d4f adaptive %s x=%g a=%g", TF->desc, desc, x, a);

  return 0;
}

int
test_deriv_adaptive(const double tol, gsl_rng * rng_p)
{
  const size_t N = 10; /* number of random points per function */
  const deriv_test_function * ptr;

  for (ptr = test_functions; ptr->test_t != NULL; ++ptr)
    {
      size_t j;

      test_deriv_adaptive_tol(tol,       ptr,  1.0, GSL_DERIV_STENCIL_CENTRAL,  GSL_DERIV_STEP_ABSOLUTE, "central absolute");
      test_deriv_adaptive_tol(tol,       ptr,  1.0, GSL_DERIV_STENCIL_CENTRAL,  GSL_DERIV_STEP_RELATIVE, "central relative");

      test_deriv_adaptive_tol(1.0e1*tol, ptr,  1.0, GSL_DERIV_STENCIL_FORWARD,  GSL_DERIV_STEP_ABSOLUTE, "forward absolute");
      test_deriv_adaptive_tol(1.0e1*tol, ptr,  1.0, GSL_DERIV_STENCIL_FORWARD,  GSL_DERIV_STEP_RELATIVE, "forward relative");

      test_deriv_adaptive_tol(1.0e1*tol, ptr,  1.0, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_ABSOLUTE, "backward absolute");
      test_deriv_adaptive_tol(1.0e1*tol, ptr,  1.0, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_RELATIVE, "backward relative");

      test_deriv_adaptive_tol(tol,       ptr, -1.0, GSL_DERIV_STENCIL_CENTRAL,  GSL_DERIV_STEP_ABSOLUTE, "central absolute");
      test_deriv_adaptive_tol(tol,       ptr, -1.0, GSL_DERIV_STENCIL_CENTRAL,  GSL_DERIV_STEP_RELATIVE, "central relative");

      test_deriv_adaptive_tol(1.0e1*tol, ptr, -1.0, GSL_DERIV_STENCIL_FORWARD,  GSL_DERIV_STEP_ABSOLUTE, "forward absolute");
      test_deriv_adaptive_tol(1.0e1*tol, ptr, -1.0, GSL_DERIV_STENCIL_FORWARD,  GSL_DERIV_STEP_RELATIVE, "forward relative");

      test_deriv_adaptive_tol(1.0e1*tol, ptr, -1.0, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_ABSOLUTE, "backward absolute");
      test_deriv_adaptive_tol(1.0e1*tol, ptr, -1.0, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_RELATIVE, "backward relative");

      for (j = 0; j < N; ++j)
        {
          double u = gsl_rng_uniform(rng_p);
          double x = ptr->xmin + u * (ptr->xmax - ptr->xmin);

          test_deriv_adaptive_tol(tol,       ptr, x, GSL_DERIV_STENCIL_CENTRAL,  GSL_DERIV_STEP_ABSOLUTE, "central absolute random");
          test_deriv_adaptive_tol(tol,       ptr, x, GSL_DERIV_STENCIL_CENTRAL,  GSL_DERIV_STEP_RELATIVE, "central relative random");

          test_deriv_adaptive_tol(1.0e3*tol, ptr, x, GSL_DERIV_STENCIL_FORWARD,  GSL_DERIV_STEP_ABSOLUTE, "forward absolute random");
          test_deriv_adaptive_tol(1.0e3*tol, ptr, x, GSL_DERIV_STENCIL_FORWARD,  GSL_DERIV_STEP_RELATIVE, "forward relative random");

          test_deriv_adaptive_tol(1.0e2*tol, ptr, x, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_ABSOLUTE, "backward absolute random");
          test_deriv_adaptive_tol(1.0e2*tol, ptr, x, GSL_DERIV_STENCIL_BACKWARD, GSL_DERIV_STEP_RELATIVE, "backward relative random");
        }
    }

  return GSL_SUCCESS;
}

int
test_deriv_fixed_tol(const double tol, const deriv_test_function * TF, const double x,
                     const gsl_deriv_stencil_t stencil_type, const double h, const char * desc)
{
  const double a = TF->a;
  const deriv_test_t * test_t = TF->test_t;
  const gsl_function F   = { test_t->f,   (void *) &a };
  const gsl_function DF  = { test_t->df,  (void *) &a };
  const gsl_function D2F = { test_t->d2f, (void *) &a };
  const double expected_df  = GSL_FN_EVAL (&DF, x);
  const double expected_d2f = GSL_FN_EVAL (&D2F, x);
  gsl_deriv_parameters params = gsl_deriv_default_parameters();
  double result, abserr;
  unsigned int neval;

  params.stencil_type = stencil_type;

  params.deriv_order = 1;
  gsl_deriv_fixed(&params, &F, x, h, &result, &abserr, &neval);
  gsl_test_rel(result, expected_df, tol, "%s df fixed %s x=%g a=%g h=%g", TF->desc, desc, x, a, h);

  params.deriv_order = 2;
  gsl_deriv_fixed(&params, &F, x, h, &result, &abserr, &neval);
  gsl_test_rel(result, expected_d2f, 1.0e2*tol, "%s d2f fixed %s x=%g a=%g h=%g", TF->desc, desc, x, a, h);

  return 0;
}

int
test_deriv_fixed(const double tol, gsl_rng * rng_p)
{
  const size_t N = 20; /* number of random points per function */
  const deriv_test_function * ptr;

  for (ptr = test_functions; ptr->test_t != NULL; ++ptr)
    {
      double h_fixed;
      size_t j;

      if (ptr->min_step > 0.0 && ptr->max_step > 0.0)
        h_fixed = sqrt(ptr->min_step * ptr->max_step);
      else
        h_fixed = GSL_ROOT3_DBL_EPSILON;

      for (j = 0; j < N; ++j)
        {
          double u = gsl_rng_uniform(rng_p);
          double x = ptr->xmin + u * (ptr->xmax - ptr->xmin);

          test_deriv_fixed_tol(tol, ptr, x, GSL_DERIV_STENCIL_CENTRAL, h_fixed, "central random");
          test_deriv_fixed_tol(1.0e1*tol, ptr, x, GSL_DERIV_STENCIL_FORWARD, h_fixed, "forward random");
          test_deriv_fixed_tol(tol, ptr, x, GSL_DERIV_STENCIL_BACKWARD, h_fixed, "backward random");
        }
    }

  return GSL_SUCCESS;
}

int
main ()
{
  const double tol = 1.0e-11;
  gsl_rng * rng_p = gsl_rng_alloc(gsl_rng_default);

  gsl_ieee_env_setup ();

#if 0
  test_deriv_adaptive_diag(&sina_test, -1.0, 100.0, 2, 1.0e-7, 1.0e-1, GSL_DERIV_STENCIL_CENTRAL);
  exit(1);
#endif

  test_fornberg(tol, rng_p);
  test_deriv_fixed(1.0e-5, rng_p);
  test_deriv_adaptive(tol, rng_p);

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


