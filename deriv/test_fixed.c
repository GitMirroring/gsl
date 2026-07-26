/* deriv/test_fixed.c
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

int
test_deriv_fixed_tol(const double tol, const deriv_test_function * TF, const double x,
                     const gsl_deriv_stencil_t stencil_type, const char * desc)
{
  const unsigned int method_order = 4;
  const double a = TF->a;
  const deriv_test_t * test_t = TF->test_t;
  const gsl_function F   = { test_t->f,   (void *) &a };
  const gsl_function DF  = { test_t->df,  (void *) &a };
  const gsl_function D2F = { test_t->d2f, (void *) &a };
  const double expected_df  = GSL_FN_EVAL (&DF, x);
  const double expected_d2f = GSL_FN_EVAL (&D2F, x);
  gsl_deriv_parameters params;
  double result, abserr, h0;
  unsigned int neval;

  params = test_parameters(1, method_order, stencil_type, GSL_DERIV_STEP_ABSOLUTE, TF->min_step, TF->max_step);
  h0 = pow(GSL_DBL_EPSILON, 1.0 / (double)(1 + method_order));
  gsl_deriv_fixed(&params, &F, x, h0, &result, &abserr, &neval);
  gsl_test_rel(result, expected_df, tol, "%s df fixed %s x=%g a=%g h=%g", TF->desc, desc, x, a, h0);

  params = test_parameters(2, method_order, stencil_type, GSL_DERIV_STEP_ABSOLUTE, TF->min_step, TF->max_step);
  h0 = pow(GSL_DBL_EPSILON, 1.0 / (double)(2 + method_order));
  gsl_deriv_fixed(&params, &F, x, h0, &result, &abserr, &neval);
  gsl_test_rel(result, expected_d2f, 1.0e2*tol, "%s d2f fixed %s x=%g a=%g h=%g", TF->desc, desc, x, a, h0);

  return 0;
}

int
test_deriv_fixed(const double tol, const deriv_test_function * functions, gsl_rng * rng_p)
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

          test_deriv_fixed_tol(tol,       ptr, x, GSL_DERIV_STENCIL_CENTRAL,  "central random");
          test_deriv_fixed_tol(1.0e1*tol, ptr, x, GSL_DERIV_STENCIL_FORWARD,  "forward random");
          test_deriv_fixed_tol(tol,       ptr, x, GSL_DERIV_STENCIL_BACKWARD, "backward random");
        }
    }

  return GSL_SUCCESS;
}
