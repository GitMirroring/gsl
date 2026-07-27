/* deriv/test_g.c
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

/* gradient test functions */

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

typedef struct
{
  double (*f) (const gsl_vector * x, void * params);
  int (*df)   (const gsl_vector * x, void * params, gsl_vector * g);
  size_t n;
  double xmin;
  double xmax;
  const char * desc;
} gradient_test_function;

double
g_qform_f(const gsl_vector * x, void * p)
{
  const size_t n = x->size;
  const gsl_matrix * A = (const gsl_matrix *) p;
  gsl_vector * work = gsl_vector_alloc(n);
  double val;

  gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, work);
  gsl_blas_ddot(x, work, &val);

  gsl_vector_free(work);

  return val;
}

int
g_qform_df(const gsl_vector * x, void * p, gsl_vector * g)
{
  const size_t n = x->size;
  const gsl_matrix * A = (const gsl_matrix *) p;
  gsl_matrix * W = gsl_matrix_alloc(n, n);

  /* W := A + A^T */
  gsl_matrix_transpose_memcpy(W, A);
  gsl_matrix_add(W, A);

  gsl_blas_dgemv(CblasNoTrans, 1.0, W, x, 0.0, g);

  gsl_matrix_free(W);

  return GSL_SUCCESS;
}

double
g_coup_f(const gsl_vector * x, void * p)
{
  const double x1 = gsl_vector_get(x, 0);
  const double x2 = gsl_vector_get(x, 1);
  const double x3 = gsl_vector_get(x, 2);
  U(p);
  return exp(x1*x2) + sin(x1*x3) + x2*x2*x3;
}

int
g_coup_df(const gsl_vector * x, void * p, gsl_vector * g)
{
  const double x1 = gsl_vector_get(x, 0);
  const double x2 = gsl_vector_get(x, 1);
  const double x3 = gsl_vector_get(x, 2);
  U(p);
  gsl_vector_set(g, 0, x2*exp(x1*x2) + x3*cos(x1*x3));
  gsl_vector_set(g, 1, x1*exp(x1*x2) + 2.0*x2*x3);
  gsl_vector_set(g, 2, x1*cos(x1*x3) + x2*x2);
  return GSL_SUCCESS;
}

static const gradient_test_function test_gradient_functions[] =
{
  { g_qform_f, g_qform_df, 0, -10.0, 10.0, "quadratic form" },
  { g_coup_f,  g_coup_df,  3, -4.0,  4.0,  "coupled" },

  { NULL, NULL, 0, 0.0, 0.0, NULL }
};
