/* deriv/test_fornberg.c
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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

int
test_fornberg_stencil(const double tol, const double x0, const unsigned int max_deriv,
                      const unsigned int deriv_order, const double stencil[],
                      const unsigned int n, const double expected[], const char * desc)
{
  double * weights = malloc(n*(max_deriv+1) * sizeof(double));
  gsl_matrix_view W = gsl_matrix_view_array(weights, max_deriv + 1, n);
  size_t j;

  gsl_deriv_fornberg_weights(stencil, n, max_deriv, x0, weights);

  for (j = 0; j < n; ++j)
    {
      double wj = gsl_matrix_get(&W.matrix, deriv_order, j);
      gsl_test_rel(wj, expected[j], tol, "fornberg %s j=%zu", desc, j);
    }

  free(weights);

  return 0;
}

/* test weights(stencil) = weights(stencil + a) */
int
test_fornberg_translate(const double tol, const double a, const double x0,
                        const double stencil[], const unsigned int n, const char * desc)
{
  const unsigned int max_deriv = n - 1;
  double * stencil_a = malloc(n * sizeof(double));
  double * weights = malloc(n*(max_deriv+1) * sizeof(double));
  double * weights_a = malloc(n*(max_deriv+1) * sizeof(double));
  size_t i, j;

  for (j = 0; j < n; ++j)
    stencil_a[j] = stencil[j] + a;

  gsl_deriv_fornberg_weights(stencil, n, max_deriv, x0, weights);
  gsl_deriv_fornberg_weights(stencil_a, n, max_deriv, x0 + a, weights_a);

  for (i = 0; i < n*(max_deriv+1); ++i)
    gsl_test_rel(weights_a[i], weights[i], tol, "%s i=%zu", desc, i);

  free(stencil_a);
  free(weights);
  free(weights_a);

  return 0;
}

/* test that Fornberg weights approximate polynomials and derivatives of order < n */
int
test_fornberg_poly(const double tol, const double x0, const unsigned int poly_order,
                   const double stencil[], const unsigned int n,
                   gsl_rng * rng_p, const char * desc)
{
  const unsigned int max_deriv = GSL_MIN(poly_order, n - 1);
  const size_t ncoef = poly_order + 1;
  double * coef = malloc(ncoef * sizeof(double));
  double * f_poly = malloc(n * sizeof(double));
  double * df_poly = malloc((max_deriv+1) * sizeof(double));
  double * weights = malloc(n*(max_deriv+1) * sizeof(double));
  gsl_matrix_view W = gsl_matrix_view_array(weights, max_deriv + 1, n);
  gsl_vector_view c = gsl_vector_view_array(coef, ncoef);
  gsl_vector_view f = gsl_vector_view_array(f_poly, n);
  size_t j;

  /* generate random polynomial coefficients */
  random_vector(-1.0, 1.0, &c.vector, rng_p);

  /* compute Fornberg weights for this stencil */
  gsl_deriv_fornberg_weights(stencil, n, max_deriv, x0, weights);

  /* evaluate polynomial at node locations */
  for (j = 0; j < n; ++j)
    f_poly[j] = gsl_poly_eval(coef, (int) ncoef, x0 + stencil[j]);

  /* evaluate polynomial derivatives at x0 */
  gsl_poly_eval_derivs(coef, ncoef, x0, df_poly, max_deriv+1);

  /* test Fornberg approximation matches exact polynomial derivatives */
  for (j = 0; j <= max_deriv; ++j)
    {
      gsl_vector_const_view wj = gsl_matrix_const_row(&W.matrix, j);
      double result;

      gsl_blas_ddot(&wj.vector, &f.vector, &result);

      gsl_test_rel(result, df_poly[j], tol, "fornberg poly %s order=%zu j=%zu", desc, poly_order, j);
    }

  free(coef);
  free(weights);
  free(f_poly);
  free(df_poly);

  return 0;
}

int
test_fornberg(const double tol, gsl_rng * rng_p)
{
  const double stencil_L2[]     = { -1.0, 1.0 };
  const double expected_L2_D0[] = {  0.5, 0.5 };
  const double expected_L2_D1[] = { -0.5, 0.5 };

  const double stencil_L3[]     = { -1.0,  0.0, 1.0 };
  const double expected_L3_D0[] = {  0.0,  1.0, 0.0 };
  const double expected_L3_D1[] = { -0.5,  0.0, 0.5 };
  const double expected_L3_D2[] = {  1.0, -2.0, 1.0 };

  const double stencil_L4[]     = {  0.0,  1.0, 2.0, 3.0 };
  const double expected_L4_D0[] = {  1.0,  0.0, 0.0, 0.0 };
  const double expected_L4_D1[] = { -11.0/6.0,  3.0, -1.5, 1.0/3.0 };
  const double expected_L4_D2[] = {  2.0, -5.0,  4.0, -1.0 };
  const double expected_L4_D3[] = { -1.0,  3.0, -3.0,  1.0 };

  const double stencil_L5[]     = { -2.0, -1.0, 0.0, 1.0, 2.0 };
  const double expected_L5_D0[] = {  0.0,  0.0, 1.0, 0.0, 0.0 };
  const double expected_L5_D1[] = {  1.0/12.0, -2.0/3.0,  0.0, 2.0/3.0, -1.0/12.0 };
  const double expected_L5_D2[] = { -1.0/12.0,  4.0/3.0, -2.5, 4.0/3.0, -1.0/12.0 };
  const double expected_L5_D3[] = { -0.5,  1.0,  0.0, -1.0, 0.5 };
  const double expected_L5_D4[] = {  1.0, -4.0,  6.0, -4.0, 1.0 };

  const double a = 100.0;
  double x0 = 0.0;
  size_t i, n;

  /* test known stencils */

  test_fornberg_stencil(tol, x0, 1, 0, stencil_L2, 2, expected_L2_D0, "stencil[L2,D0]");
  test_fornberg_stencil(tol, x0, 1, 1, stencil_L2, 2, expected_L2_D1, "stencil[L2,D1]");

  test_fornberg_stencil(tol, x0, 2, 0, stencil_L3, 3, expected_L3_D0, "stencil[L3,D0]");
  test_fornberg_stencil(tol, x0, 2, 1, stencil_L3, 3, expected_L3_D1, "stencil[L3,D1]");
  test_fornberg_stencil(tol, x0, 2, 2, stencil_L3, 3, expected_L3_D2, "stencil[L3,D2]");

  test_fornberg_stencil(tol, x0, 3, 0, stencil_L4, 4, expected_L4_D0, "stencil[L4,D0]");
  test_fornberg_stencil(tol, x0, 3, 1, stencil_L4, 4, expected_L4_D1, "stencil[L4,D1]");
  test_fornberg_stencil(tol, x0, 3, 2, stencil_L4, 4, expected_L4_D2, "stencil[L4,D2]");
  test_fornberg_stencil(tol, x0, 3, 3, stencil_L4, 4, expected_L4_D3, "stencil[L4,D3]");

  test_fornberg_stencil(tol, x0, 4, 0, stencil_L5, 5, expected_L5_D0, "stencil[L5,D0]");
  test_fornberg_stencil(tol, x0, 4, 1, stencil_L5, 5, expected_L5_D1, "stencil[L5,D1]");
  test_fornberg_stencil(tol, x0, 4, 2, stencil_L5, 5, expected_L5_D2, "stencil[L5,D2]");
  test_fornberg_stencil(tol, x0, 4, 3, stencil_L5, 5, expected_L5_D3, "stencil[L5,D3]");
  test_fornberg_stencil(tol, x0, 4, 4, stencil_L5, 5, expected_L5_D4, "stencil[L5,D4]");

  /* test translated stencils */

  test_fornberg_translate(tol, a, x0, stencil_L2, 2, "translated stencil[L2]");
  test_fornberg_translate(tol, a, x0, stencil_L3, 3, "translated stencil[L3]");
  test_fornberg_translate(tol, a, x0, stencil_L4, 4, "translated stencil[L4]");
  test_fornberg_translate(tol, a, x0, stencil_L5, 5, "translated stencil[L5]");
  
  /* polynomial exactness tests, where poly_order < stencil_length */

  for (i = 0; i < 2; ++i)
    test_fornberg_poly(tol, x0, i, stencil_L2, 2, rng_p, "stencil[L2]");

  for (i = 0; i < 3; ++i)
    test_fornberg_poly(tol, x0, i, stencil_L3, 3, rng_p, "stencil[L3]");

  for (i = 0; i < 4; ++i)
    test_fornberg_poly(tol, x0, i, stencil_L4, 4, rng_p, "stencil[L4]");

  for (i = 0; i < 5; ++i)
    test_fornberg_poly(tol, x0, i, stencil_L5, 5, rng_p, "stencil[L5]");

  /* test random stencils */

  for (n = 2; n <= 7; ++n)
    {
      double * stencil = malloc(n * sizeof(double));
      gsl_vector_view sv = gsl_vector_view_array(stencil, n);

      for (i = 0; i < 5; ++i)
        {
          size_t j;

          random_vector(-1.0, 1.0, &sv.vector, rng_p);

          /* test translated stencil */
          test_fornberg_translate(tol, a, x0, stencil, n, "translated random stencil");

          /* test polynomial approximation */
          for (j = 0; j < n; ++j)
            test_fornberg_poly(1.0e3*tol, x0, j, stencil, n, rng_p, "random stencil");
        }

      free(stencil);
    }

  return 0;
}
