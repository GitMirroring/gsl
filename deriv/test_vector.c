/* deriv/test_vector.c
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

/*
 * Tests for gsl_deriv_fixed_v() and gsl_deriv_adaptive_v(): each element
 * of the vector result must equal the scalar result at that point, the
 * evaluation counts must add up, user scales must be taken element by
 * element, and vectors of the wrong length must be rejected
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

#define TEST_VECTOR_N 5   /* number of points per call */

static int
test_deriv_vector_one(const deriv_test_function * TF, const unsigned int deriv_order,
                      const gsl_deriv_stencil_t stencil_type,
                      const gsl_deriv_step_policy_t step_policy, gsl_rng * rng_p)
{
  const size_t n = TEST_VECTOR_N;
  const double a = TF->a;
  const gsl_function F = { TF->test_t->f, (void *) &a };
  const double h = test_h0(deriv_order, 4);
  /* strided vectors, to test element access */
  double x_data[2 * TEST_VECTOR_N], r_data[2 * TEST_VECTOR_N], e_data[2 * TEST_VECTOR_N];
  double s_data[2 * TEST_VECTOR_N], scale_data[TEST_VECTOR_N];
  gsl_vector_view x = gsl_vector_view_array_with_stride(x_data, 2, n);
  gsl_vector_view result = gsl_vector_view_array_with_stride(r_data, 2, n);
  gsl_vector_view esterr = gsl_vector_view_array_with_stride(e_data, 2, n);
  gsl_vector_view final_step = gsl_vector_view_array_with_stride(s_data, 2, n);
  gsl_vector_view scales = gsl_vector_view_array(scale_data, n);
  gsl_deriv_parameters params = test_parameters(deriv_order, 4, stencil_type, step_policy);
  unsigned int neval_v, neval_sum;
  char desc[256];
  size_t i;

  sprintf(desc, "%s deriv_order=%u stencil=%d step_policy=%d", TF->desc, deriv_order,
          (int) stencil_type, (int) step_policy);

  for (i = 0; i < n; ++i)
    {
      double u = gsl_rng_uniform(rng_p);

      gsl_vector_set(&x.vector, i, TF->xmin + u * (TF->xmax - TF->xmin));

      /* a different scale for each element */
      scale_data[i] = (TF->scale > 0.0) ? TF->scale * (1.0 + 0.5 * i) : 1.0;
    }

  if (step_policy == GSL_DERIV_STEP_USER_SCALE)
    params.scales = &scales.vector;

  /* fixed step */

  gsl_deriv_fixed_v(&params, &F, &x.vector, h, &result.vector, &esterr.vector, &neval_v);

  neval_sum = 0;
  for (i = 0; i < n; ++i)
    {
      double xi = gsl_vector_get(&x.vector, i);
      double r, e;
      unsigned int neval;

      gsl_deriv_fixed(&params, &F, xi, h, &r, &e, &neval);
      neval_sum += neval;

      gsl_test(gsl_vector_get(&result.vector, i) != r, "%s fixed_v result i=%zu", desc, i);
      gsl_test(gsl_vector_get(&esterr.vector, i) != e, "%s fixed_v esterr i=%zu", desc, i);
    }

  gsl_test(neval_v != neval_sum, "%s fixed_v neval %u, sum of scalar calls %u", desc, neval_v, neval_sum);

  /* adaptive */

  gsl_deriv_adaptive_v(&params, &F, &x.vector, &result.vector, &esterr.vector, &final_step.vector,
                       &neval_v, NULL);

  neval_sum = 0;
  for (i = 0; i < n; ++i)
    {
      double xi = gsl_vector_get(&x.vector, i);
      gsl_vector_view scale_i = gsl_vector_view_array(&scale_data[i], 1);
      gsl_deriv_parameters params_i = params;
      double r, e, fs;
      unsigned int neval;

      if (step_policy == GSL_DERIV_STEP_USER_SCALE)
        params_i.scales = &scale_i.vector;

      gsl_deriv_adaptive(&params_i, &F, xi, &r, &e, &fs, &neval, NULL);
      neval_sum += neval;

      gsl_test(gsl_vector_get(&result.vector, i) != r, "%s adaptive_v result i=%zu", desc, i);
      gsl_test(gsl_vector_get(&esterr.vector, i) != e, "%s adaptive_v esterr i=%zu", desc, i);
      gsl_test(gsl_vector_get(&final_step.vector, i) != fs, "%s adaptive_v final_step i=%zu", desc, i);
    }

  gsl_test(neval_v != neval_sum, "%s adaptive_v neval %u, sum of scalar calls %u", desc, neval_v, neval_sum);

  return GSL_SUCCESS;
}

/* vectors of the wrong length must be rejected */
static int
test_deriv_vector_badlen(void)
{
  const double a = 0.0;
  const gsl_function F = { exp_f, (void *) &a };
  double x_data[3] = { 0.1, 0.2, 0.3 }, r_data[3], e_data[3], s_data[3], scale_data[3] = { 1.0, 1.0, 1.0 };
  gsl_vector_view x = gsl_vector_view_array(x_data, 3);
  gsl_vector_view r = gsl_vector_view_array(r_data, 3), r2 = gsl_vector_view_array(r_data, 2);
  gsl_vector_view e = gsl_vector_view_array(e_data, 3), e2 = gsl_vector_view_array(e_data, 2);
  gsl_vector_view s = gsl_vector_view_array(s_data, 3), s2 = gsl_vector_view_array(s_data, 2);
  gsl_vector_view scales2 = gsl_vector_view_array(scale_data, 2);
  gsl_deriv_parameters params = test_parameters(1, 4, GSL_DERIV_STENCIL_CENTRAL, GSL_DERIV_STEP_ABSOLUTE);
  gsl_error_handler_t * handler = gsl_set_error_handler_off();
  unsigned int neval;
  int status;

  status = gsl_deriv_fixed_v(&params, &F, &x.vector, 1.0e-3, &r2.vector, &e.vector, &neval);
  gsl_test_int(status, GSL_EBADLEN, "fixed_v short result");

  status = gsl_deriv_fixed_v(&params, &F, &x.vector, 1.0e-3, &r.vector, &e2.vector, &neval);
  gsl_test_int(status, GSL_EBADLEN, "fixed_v short esterr");

  status = gsl_deriv_adaptive_v(&params, &F, &x.vector, &r2.vector, &e.vector, &s.vector, &neval, NULL);
  gsl_test_int(status, GSL_EBADLEN, "adaptive_v short result");

  status = gsl_deriv_adaptive_v(&params, &F, &x.vector, &r.vector, &e2.vector, &s.vector, &neval, NULL);
  gsl_test_int(status, GSL_EBADLEN, "adaptive_v short esterr");

  status = gsl_deriv_adaptive_v(&params, &F, &x.vector, &r.vector, &e.vector, &s2.vector, &neval, NULL);
  gsl_test_int(status, GSL_EBADLEN, "adaptive_v short final_step");

  params.step_policy = GSL_DERIV_STEP_USER_SCALE;
  params.scales = &scales2.vector;
  status = gsl_deriv_adaptive_v(&params, &F, &x.vector, &r.vector, &e.vector, &s.vector, &neval, NULL);
  gsl_test_int(status, GSL_EBADLEN, "adaptive_v short scales");

  gsl_set_error_handler(handler);

  return GSL_SUCCESS;
}

int
test_deriv_vector(const deriv_test_function * normal, const deriv_test_function * scaled,
                  gsl_rng * rng_p)
{
  const deriv_test_function * ptr;
  unsigned int deriv_order;

  for (deriv_order = 1; deriv_order <= 2; ++deriv_order)
    {
      for (ptr = normal; ptr->test_t != NULL; ++ptr)
        {
          test_deriv_vector_one(ptr, deriv_order, GSL_DERIV_STENCIL_CENTRAL, GSL_DERIV_STEP_ABSOLUTE, rng_p);
          test_deriv_vector_one(ptr, deriv_order, GSL_DERIV_STENCIL_FORWARD, GSL_DERIV_STEP_RELATIVE, rng_p);
        }

      for (ptr = scaled; ptr->test_t != NULL; ++ptr)
        test_deriv_vector_one(ptr, deriv_order, GSL_DERIV_STENCIL_CENTRAL, GSL_DERIV_STEP_USER_SCALE, rng_p);
    }

  test_deriv_vector_badlen();

  return GSL_SUCCESS;
}
