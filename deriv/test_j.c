/* deriv/test_j.c
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

/* Jacobian test functions */

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
  int (*f)  (const gsl_vector * x, void * params, gsl_vector * y);
  int (*df) (const gsl_vector * x, void * params, gsl_matrix * J);
  size_t n;    /* number of inputs */
  size_t m;    /* number of outputs */
  double xmin;
  double xmax;
  const char * desc;
} jacobian_test_function;

int
j_const_f(const gsl_vector * x, void * p, gsl_vector * y)
{
  U(x);
  U(p);
  gsl_vector_set(y, 0,  5.0);
  gsl_vector_set(y, 1, -2.0);
  gsl_vector_set(y, 2,  17.0);
  return GSL_SUCCESS;
}

int
j_const_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
  U(x);
  U(p);
  gsl_matrix_set_zero(J);
  return GSL_SUCCESS;
}

int
j_linear_f(const gsl_vector * x, void * p, gsl_vector * y)
{
  const gsl_matrix * A = (const gsl_matrix *) p;

  U(x);
  gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, y);
  gsl_vector_add_constant(y, 0.2);

  return GSL_SUCCESS;
}

int
j_linear_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
  const gsl_matrix * A = (const gsl_matrix *) p;
  U(x);
  gsl_matrix_memcpy(J, A);
  return GSL_SUCCESS;
}

int
j_quad_f(const gsl_vector * x, void * p, gsl_vector * y)
{
  const double x1 = gsl_vector_get(x, 0);
  const double x2 = gsl_vector_get(x, 1);
  U(p);
  gsl_vector_set(y, 0, x1*x1);
  gsl_vector_set(y, 1, x1*x2);
  gsl_vector_set(y, 2, x2*x2);
  return GSL_SUCCESS;
}

int
j_quad_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
  const double x1 = gsl_vector_get(x, 0);
  const double x2 = gsl_vector_get(x, 1);
  U(p);
  gsl_matrix_set(J, 0, 0, 2.0*x1);
  gsl_matrix_set(J, 0, 1, 0.0);
  gsl_matrix_set(J, 1, 0, x2);
  gsl_matrix_set(J, 1, 1, x1);
  gsl_matrix_set(J, 2, 0, 0.0);
  gsl_matrix_set(J, 2, 1, 2.0*x2);
  return GSL_SUCCESS;
}

int
j_sep_f(const gsl_vector * x, void * p, gsl_vector * y)
{
  const double x1 = gsl_vector_get(x, 0);
  const double x2 = gsl_vector_get(x, 1);
  const double x3 = gsl_vector_get(x, 2);
  U(p);
  gsl_vector_set(y, 0, exp(x1));
  gsl_vector_set(y, 1, sin(x2));
  gsl_vector_set(y, 2, cos(x3));
  return GSL_SUCCESS;
}

int
j_sep_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
  const double x1 = gsl_vector_get(x, 0);
  const double x2 = gsl_vector_get(x, 1);
  const double x3 = gsl_vector_get(x, 2);
  U(p);
  gsl_matrix_set_zero(J);
  gsl_matrix_set(J, 0, 0,  exp(x1));
  gsl_matrix_set(J, 1, 1,  cos(x2));
  gsl_matrix_set(J, 2, 2, -sin(x3));
  return GSL_SUCCESS;
}

int
j_coup_f(const gsl_vector * x, void * p, gsl_vector * y)
{
  const double x1 = gsl_vector_get(x, 0);
  const double x2 = gsl_vector_get(x, 1);
  const double x3 = gsl_vector_get(x, 2);
  U(p);
  gsl_vector_set(y, 0, exp(x1*x2));
  gsl_vector_set(y, 1, sin(x1+x2+x3));
  gsl_vector_set(y, 2, log(1.0+x1*x1+x2*x2+x3*x3));
  return GSL_SUCCESS;
}

int
j_coup_df(const gsl_vector * x, void * p, gsl_matrix * J)
{
  const double x1 = gsl_vector_get(x, 0);
  const double x2 = gsl_vector_get(x, 1);
  const double x3 = gsl_vector_get(x, 2);
  const double denom = 1.0 + x1*x1 + x2*x2 + x3*x3;
  U(p);
  gsl_matrix_set(J, 0, 0,  x2*exp(x1*x2));
  gsl_matrix_set(J, 0, 1,  x1*exp(x1*x2));
  gsl_matrix_set(J, 0, 2,  0.0);
  gsl_matrix_set(J, 1, 0,  cos(x1+x2+x3));
  gsl_matrix_set(J, 1, 1,  cos(x1+x2+x3));
  gsl_matrix_set(J, 1, 2,  cos(x1+x2+x3));
  gsl_matrix_set(J, 2, 0,  2.0*x1/denom);
  gsl_matrix_set(J, 2, 1,  2.0*x2/denom);
  gsl_matrix_set(J, 2, 2,  2.0*x3/denom);
  return GSL_SUCCESS;
}

static const jacobian_test_function test_jacobian_functions[] =
{
  { j_const_f,  j_const_df,  3, 3, -10.0, 10.0, "constant" },
  { j_linear_f, j_linear_df, 0, 0, -10.0, 10.0, "linear" },
  { j_quad_f,   j_quad_df,   2, 3, -10.0, 10.0, "quadratic" },
  { j_sep_f,    j_sep_df,    3, 3, -5.0,  5.0,  "separable" },
  { j_coup_f,   j_coup_df,   3, 3, -3.0,  3.0,  "coupled" },

  { NULL,       NULL,        0, 0,  0.0,  0.0,  NULL }
};
