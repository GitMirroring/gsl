/* deriv/test_esterr.c
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
 * Tests of the error estimates and candidate selection of the adaptive algorithm:
 *
 * 1. the estimated error must bound the actual error (within a safety factor) for
 *    the standard test functions; with the original one degree of freedom error
 *    estimate of each Romberg window, the actual error exceeded the estimate by
 *    factors up to ~1e3 for these functions
 *
 * 2. functions with noise much larger than machine precision must be differentiated
 *    accurately; selecting the candidate with the smallest error estimate must not
 *    pick a window in the rounding-dominated region whose estimate happens to be small
 *
 * 3. the MAD candidate filter must ignore non-finite candidates
 */

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include "deriv_cand.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_test.h>

/*
 * deterministic noise in [-1,1] depending on the bits of x; computed with exact
 * floating point operations (frexp, ldexp, floor) and 32 bit unsigned arithmetic,
 * so that it is identical on all IEEE platforms
 */
static double
test_noise(const double x)
{
  const double two26 = 67108864.0;
  int e;
  double m = frexp(x, &e);
  double t = ldexp(fabs(m), 53);         /* integer < 2^53 */
  double thi = floor(t / two26);
  unsigned long hi = (unsigned long) thi;
  unsigned long lo = (unsigned long) (t - thi * two26);
  unsigned long h;

  h = (hi * 2654435761UL) & 0xffffffffUL;
  h ^= (lo * 2246822519UL) & 0xffffffffUL;
  h ^= ((unsigned long) (e + 2048) * 3266489917UL) & 0xffffffffUL;
  h ^= (x < 0.0) ? 0x9e3779b9UL : 0UL;
  h ^= h >> 16;
  h = (h * 0x45d9f3bUL) & 0xffffffffUL;
  h ^= h >> 16;
  h = (h * 0x45d9f3bUL) & 0xffffffffUL;
  h ^= h >> 16;

  return 2.0 * ((double) h / 4294967295.0) - 1.0;
}

/* sin(x) with deterministic relative noise of size 1e-12 */
static double
noisy_sin(double x, void * params)
{
  U(params);
  return sin(x) + 1.0e-12 * test_noise(x);
}

/* test 1: |result - exact| <= fac * (esterr + 10 eps |exact|) */
static int
test_esterr_bound(const double fac, const deriv_test_function * functions, gsl_rng * rng_p)
{
  const size_t N = 20; /* number of random points per function */
  const gsl_deriv_stencil_t stencils[3] = { GSL_DERIV_STENCIL_CENTRAL, GSL_DERIV_STENCIL_FORWARD,
                                            GSL_DERIV_STENCIL_BACKWARD };
  const char * stencil_names[3] = { "central", "forward", "backward" };
  const deriv_test_function * TF;

  for (TF = functions; TF->test_t != NULL; ++TF)
    {
      size_t j;

      for (j = 0; j < N; ++j)
        {
          const double a = TF->a;
          const double x = TF->xmin + gsl_rng_uniform(rng_p) * (TF->xmax - TF->xmin);
          const gsl_function F = { TF->test_t->f, (void *) &a };
          unsigned int m;

          for (m = 1; m <= 4; ++m)
            {
              const gsl_function DF = { (m == 1) ? TF->test_t->df : (m == 2) ? TF->test_t->d2f :
                                        (m == 3) ? TF->test_t->d3f : TF->test_t->d4f, (void *) &a };
              const double expected = GSL_FN_EVAL(&DF, x);
              size_t k;
              int policy;

              for (k = 0; k < 3; ++k)
                {
                  for (policy = GSL_DERIV_STEP_ABSOLUTE; policy <= GSL_DERIV_STEP_RELATIVE; ++policy)
                    {
                      gsl_deriv_parameters params = gsl_deriv_default_parameters(m, 4);
                      double result, abserr, fs, bound;
                      unsigned int neval;

                      params.stencil_type = stencils[k];
                      params.step_policy = (gsl_deriv_step_policy_t) policy;

                      gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, NULL);

                      bound = fac * (abserr + 10.0 * GSL_DBL_EPSILON * fabs(expected));

                      gsl_test(!(fabs(result - expected) <= bound),
                               "%s d%uf esterr %s %s x=%.12e |error|=%.3e esterr=%.3e",
                               TF->desc, m, stencil_names[k],
                               (policy == GSL_DERIV_STEP_ABSOLUTE) ? "absolute" : "relative",
                               x, fabs(result - expected), abserr);
                    }
                }
            }
        }
    }

  return GSL_SUCCESS;
}

/* test 2: noisy function */
static int
test_esterr_noise(gsl_rng * rng_p)
{
  const size_t N = 400; /* number of random points */
  const gsl_deriv_stencil_t stencils[3] = { GSL_DERIV_STENCIL_CENTRAL, GSL_DERIV_STENCIL_FORWARD,
                                            GSL_DERIV_STENCIL_BACKWARD };
  const char * stencil_names[3] = { "central", "forward", "backward" };
  /*
   * absolute tolerances tol[method_order/2 - 1][deriv_order - 1]; the largest errors of the
   * algorithm are 1.4e-7, 3.0e-6 (method order 2) and 9.0e-9, 3.3e-7 (method order 4); with the
   * original error estimate the first derivatives had errors up to 7.2e-7 and 1.2e-7
   */
  const double tol[2][2] = { { 4.0e-7, 1.5e-5 }, { 5.0e-8, 1.5e-6 } };
  const gsl_function F = { noisy_sin, NULL };
  size_t j;

  for (j = 0; j < N; ++j)
    {
      const double x = -M_PI + 2.0 * M_PI * gsl_rng_uniform(rng_p);
      unsigned int method_order, m;

      for (method_order = 2; method_order <= 4; method_order += 2)
        {
          for (m = 1; m <= 2; ++m)
            {
              const double expected = (m == 1) ? cos(x) : -sin(x);
              size_t k;

              for (k = 0; k < 3; ++k)
                {
                  gsl_deriv_parameters params = gsl_deriv_default_parameters(m, method_order);
                  double result, abserr, fs;
                  unsigned int neval;

                  params.stencil_type = stencils[k];

                  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, NULL);

                  gsl_test_abs(result, expected, tol[method_order / 2 - 1][m - 1], "noisy sin(x) d%uf p=%u %s x=%.12e",
                               m, method_order, stencil_names[k], x);
                }
            }
        }
    }

  return GSL_SUCCESS;
}

/* test 3: MAD filter with non-finite candidates */
static int
test_esterr_filter_mad(void)
{
  double der_data[8] = { 1.0, 1.0 + 1.0e-12, GSL_NAN, 1.0 - 1.0e-12, 5.0, GSL_POSINF, 1.0 + 2.0e-12, 1.0 };
  double err_data[8] = { 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12 };
  gsl_vector_view der = gsl_vector_view_array(der_data, 8);
  gsl_vector_view err = gsl_vector_view_array(err_data, 8);
  gsl_deriv_cand_t cand;
  int status, nonfinite = 0, outlier = 0;
  size_t i;

  status = (gsl_deriv_cand_filter_mad->filter)(&der.vector, &err.vector, &cand);
  gsl_test(status, "filter_mad status with non-finite candidates");
  /* the 5 finite values near 1 are kept; NaN, +Inf and the outlier 5 are discarded */
  gsl_test(cand.nkeep != 5, "filter_mad nkeep with non-finite candidates (%u observed vs 5 expected)",
           (unsigned int) cand.nkeep);

  for (i = 0; i < cand.nkeep; ++i)
    {
      const double di = der_data[cand.index[i]];
      nonfinite += !gsl_finite(di);
      outlier += (di == 5.0);
    }

  gsl_test(nonfinite, "filter_mad keeps non-finite candidates");
  gsl_test(outlier, "filter_mad keeps outlier");

  return GSL_SUCCESS;
}

static double
test_sqrt_f(double x, void * params)
{
  U(params);
  return sqrt(x);           /* NaN for x < 0 */
}

static int
test_sqrt_vec(const gsl_vector * x, void * params, gsl_vector * f)
{
  U(params);
  gsl_vector_set(f, 0, sqrt(gsl_vector_get(x, 0)));
  return GSL_SUCCESS;
}

/* test 4: functions undefined over part of the step range */
static int
test_esterr_domain(void)
{
  gsl_deriv_parameters params = gsl_deriv_default_parameters(1, 4);
  gsl_function F = { test_sqrt_f, NULL };
  gsl_error_handler_t * handler = gsl_set_error_handler_off();
  double result, abserr, fs;
  unsigned int neval;
  gsl_deriv_diag diag;
  int status;

  /* some central stencils reach x < 0: the remaining candidates must suffice */
  status = gsl_deriv_adaptive(&params, &F, 0.01, &result, &abserr, &fs, &neval, &diag);
  gsl_test(status, "sqrt x=0.01 status");
  gsl_test_rel(result, 5.0, 1.0e-12, "sqrt x=0.01 derivative");
  gsl_test(fabs(result - 5.0) > abserr, "sqrt x=0.01 error bound (error %.3e, esterr %.3e)",
           fabs(result - 5.0), abserr);
  gsl_test(diag.cand.nkeep >= diag.nromb, "sqrt x=0.01 NaN candidates discarded");

  /* all central stencils reach x < 0: no finite candidate */
  status = gsl_deriv_adaptive(&params, &F, 1.0e-6, &result, &abserr, &fs, &neval, NULL);
  gsl_test_int(status, GSL_EDOM, "sqrt x=1e-6 status");
  gsl_test(!gsl_isnan(result) || !gsl_isnan(abserr) || !gsl_isnan(fs), "sqrt x=1e-6 outputs NaN");

  {
    double x_data[2] = { 0.01, 1.0e-6 }, r_data[2], e_data[2], s_data[2];
    gsl_vector_view x = gsl_vector_view_array(x_data, 2);
    gsl_vector_view r = gsl_vector_view_array(r_data, 2);
    gsl_vector_view e = gsl_vector_view_array(e_data, 2);
    gsl_vector_view s = gsl_vector_view_array(s_data, 2);

    status = gsl_deriv_adaptive_v(&params, &F, &x.vector, &r.vector, &e.vector, &s.vector, &neval, NULL);
    gsl_test_int(status, GSL_EDOM, "sqrt adaptive_v status");
    gsl_test_rel(r_data[0], 5.0, 1.0e-12, "sqrt adaptive_v derivative before the failing point");
    gsl_test(!gsl_isnan(r_data[1]) || !gsl_isnan(e_data[1]) || !gsl_isnan(s_data[1]),
             "sqrt adaptive_v outputs NaN at the failing point");
  }

  {
    gsl_deriv_vec_function FV = { test_sqrt_vec, 1, 1, NULL };
    gsl_deriv_jacobian_workspace * w = gsl_deriv_jacobian_alloc(&params, 1, 1);
    gsl_matrix * J = gsl_matrix_alloc(1, 1);
    gsl_matrix * E = gsl_matrix_alloc(1, 1);
    gsl_matrix * S = gsl_matrix_alloc(1, 1);
    gsl_vector * x = gsl_vector_alloc(1);

    gsl_vector_set(x, 0, 1.0e-6);
    status = gsl_deriv_jacobian(&FV, x, J, E, S, &neval, w);
    gsl_test_int(status, GSL_EDOM, "sqrt jacobian status");
    gsl_test(!gsl_isnan(gsl_matrix_get(J, 0, 0)) || !gsl_isnan(gsl_matrix_get(E, 0, 0)),
             "sqrt jacobian outputs NaN");

    gsl_deriv_jacobian_free(w);
    gsl_matrix_free(J);
    gsl_matrix_free(E);
    gsl_matrix_free(S);
    gsl_vector_free(x);
  }

  gsl_set_error_handler(handler);

  return GSL_SUCCESS;
}

/* test 5: diag flags for a selection at an end of the step range */
static int
test_esterr_flags(void)
{
  const double a = 0.0;
  gsl_function F = { exp_f, (void *) &a };
  gsl_deriv_parameters params = gsl_deriv_default_parameters(1, 4);
  double result, abserr, fs;
  unsigned int neval;
  gsl_deriv_diag diag;

  /* default range contains the plateau */
  gsl_deriv_adaptive(&params, &F, 0.0, &result, &abserr, &fs, &neval, &diag);
  gsl_test(diag.flags != 0, "exp flags default range (%u)", diag.flags);

  /* steps far too small: roundoff dominates everywhere, the largest steps are best */
  params.min_step = 1.0e-10;
  params.max_step = 1.0e-8;
  gsl_deriv_adaptive(&params, &F, 0.0, &result, &abserr, &fs, &neval, &diag);
  gsl_test(diag.flags != GSL_DERIV_DIAG_LARGEST_STEPS, "exp flags small steps (%u)", diag.flags);

  /* steps far too large: truncation dominates everywhere, the smallest steps are best */
  params.min_step = 0.5;
  params.max_step = 50.0;
  gsl_deriv_adaptive(&params, &F, 0.0, &result, &abserr, &fs, &neval, &diag);
  gsl_test(diag.flags != GSL_DERIV_DIAG_SMALLEST_STEPS, "exp flags large steps (%u)", diag.flags);

  return GSL_SUCCESS;
}

/* test 6: invalid parameters */
static int
test_esterr_params(void)
{
  const double a = 0.0;
  gsl_function F = { exp_f, (void *) &a };
  gsl_deriv_parameters params;
  gsl_error_handler_t * handler = gsl_set_error_handler_off();
  double result, abserr, fs;
  unsigned int neval;
  int status;

  params = gsl_deriv_default_parameters(1, 4);
  params.flags = 1;
  status = gsl_deriv_adaptive(&params, &F, 0.0, &result, &abserr, &fs, &neval, NULL);
  gsl_test_int(status, GSL_EINVAL, "adaptive nonzero flags");
  status = gsl_deriv_fixed(&params, &F, 0.0, 1.0e-3, &result, &abserr, &neval);
  gsl_test_int(status, GSL_EINVAL, "fixed nonzero flags");

  params = gsl_deriv_default_parameters(1, 4);
  params.filter = (gsl_deriv_filter_t) 99;
  status = gsl_deriv_adaptive(&params, &F, 0.0, &result, &abserr, &fs, &neval, NULL);
  gsl_test_int(status, GSL_EINVAL, "adaptive unknown filter");

  params = gsl_deriv_default_parameters(1, 4);
  params.nstep = GSL_DERIV_MAX_STEPS + 1;
  status = gsl_deriv_adaptive(&params, &F, 0.0, &result, &abserr, &fs, &neval, NULL);
  gsl_test_int(status, GSL_EDOM, "adaptive nstep > GSL_DERIV_MAX_STEPS");
  gsl_test(gsl_deriv_jacobian_alloc(&params, 1, 1) != NULL, "jacobian_alloc nstep > GSL_DERIV_MAX_STEPS");

  params = gsl_deriv_default_parameters(1, 4);
  params.nstep = params.romberg_terms + 1;
  status = gsl_deriv_adaptive(&params, &F, 0.0, &result, &abserr, &fs, &neval, NULL);
  gsl_test_int(status, GSL_EDOM, "adaptive nstep without a Romberg window");
  gsl_test(gsl_deriv_jacobian_alloc(&params, 1, 1) != NULL, "jacobian_alloc nstep without a Romberg window");

  params = gsl_deriv_default_parameters(1, 4);
  params.min_step = params.max_step;
  status = gsl_deriv_adaptive(&params, &F, 0.0, &result, &abserr, &fs, &neval, NULL);
  gsl_test_int(status, GSL_EDOM, "adaptive min_step = max_step");

  /* the largest allowed number of steps works */
  params = gsl_deriv_default_parameters(1, 4);
  params.nstep = GSL_DERIV_MAX_STEPS;
  status = gsl_deriv_adaptive(&params, &F, 0.0, &result, &abserr, &fs, &neval, NULL);
  gsl_test(status, "adaptive nstep = GSL_DERIV_MAX_STEPS status");
  gsl_test_rel(result, 1.0, 1.0e-12, "adaptive nstep = GSL_DERIV_MAX_STEPS derivative");

  gsl_set_error_handler(handler);

  return GSL_SUCCESS;
}

int
test_esterr(gsl_rng * rng_p)
{
  test_esterr_bound(20.0, test_normal_functions, rng_p);
  test_esterr_bound(20.0, test_poly_functions, rng_p);
  test_esterr_noise(rng_p);
  test_esterr_filter_mad();
  test_esterr_domain();
  test_esterr_flags();
  test_esterr_params();

  return GSL_SUCCESS;
}
