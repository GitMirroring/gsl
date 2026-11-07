/* linalg/test_common.c
 *
 * Copyright (C) 2017, 2018, 2019, 2020 Patrick Alken
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

static int create_random_vector(gsl_vector * v, gsl_rng * r);
static int create_random_matrix(gsl_matrix * m, gsl_rng * r);
static int create_posdef_matrix(gsl_matrix * m, gsl_rng * r);
static int create_hilbert_matrix2(gsl_matrix * m);

static int
create_random_vector(gsl_vector * v, gsl_rng * r)
{
  const size_t N = v->size;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      double vi = gsl_rng_uniform(r);
      gsl_vector_set(v, i, vi);
    }

  return GSL_SUCCESS;
}

static int
create_random_complex_vector(gsl_vector_complex * v, gsl_rng * r)
{
  const size_t N = v->size;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      gsl_complex vi;
      GSL_REAL(vi) = gsl_rng_uniform(r);
      GSL_IMAG(vi) = gsl_rng_uniform(r);
      gsl_vector_complex_set(v, i, vi);
    }

  return GSL_SUCCESS;
}

static int
create_random_matrix(gsl_matrix * m, gsl_rng * r)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double mij = gsl_rng_uniform(r);
          gsl_matrix_set(m, i, j, mij);
        }
    }

  return GSL_SUCCESS;
}

static int
create_random_complex_matrix(gsl_matrix_complex * m, gsl_rng * r)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          gsl_complex mij;

          GSL_REAL(mij) = gsl_rng_uniform(r);
          GSL_IMAG(mij) = gsl_rng_uniform(r);

          gsl_matrix_complex_set(m, i, j, mij);
        }
    }

  return GSL_SUCCESS;
}

static int
create_symm_matrix(gsl_matrix * m, gsl_rng * r)
{
  const size_t N = m->size1;
  size_t i, j;

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j <= i; ++j)
        {
          double mij = gsl_rng_uniform(r);
          gsl_matrix_set(m, i, j, mij);
        }
    }

  /* copy lower triangle to upper */
  gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, m, m);

  return GSL_SUCCESS;
}

static int
create_herm_matrix(gsl_matrix_complex * m, gsl_rng * r)
{
  const size_t N = m->size1;
  size_t i, j;

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j <= i; ++j)
        {
          double re = gsl_rng_uniform(r);
          double im = (i != j) ? gsl_rng_uniform(r) : 0.0;
          gsl_complex z = gsl_complex_rect(re, im);

          gsl_matrix_complex_set(m, i, j, z);

          if (i != j)
            gsl_matrix_complex_set(m, j, i, gsl_complex_conjugate(z));
        }
    }

  return GSL_SUCCESS;
}

/* create symmetric banded matrix with p sub/super-diagonals */
static int
create_symm_band_matrix(const size_t p, gsl_matrix * m, gsl_rng * r)
{
  size_t i;

  gsl_matrix_set_zero(m);

  for (i = 0; i < p + 1; ++i)
    {
      gsl_vector_view subdiag = gsl_matrix_subdiagonal(m, i);
      create_random_vector(&subdiag.vector, r);

      if (i > 0)
        {
          gsl_vector_view superdiag = gsl_matrix_superdiagonal(m, i);
          gsl_vector_memcpy(&superdiag.vector, &subdiag.vector);
        }
    }

  return GSL_SUCCESS;
}

/* create (p,q) banded matrix */
static int
create_band_matrix(const size_t p, const size_t q, gsl_matrix * m, gsl_rng * r)
{
  size_t i;

  gsl_matrix_set_zero(m);

  for (i = 0; i <= p; ++i)
    {
      gsl_vector_view v = gsl_matrix_subdiagonal(m, i);
      create_random_vector(&v.vector, r);
    }

  for (i = 1; i <= q; ++i)
    {
      gsl_vector_view v = gsl_matrix_superdiagonal(m, i);
      create_random_vector(&v.vector, r);
    }

  return GSL_SUCCESS;
}

static int
create_posdef_matrix(gsl_matrix * m, gsl_rng * r)
{
  const size_t N = m->size1;
  const double alpha = 10.0 * N;
  size_t i;

  /* The idea is to make a symmetric diagonally dominant
   * matrix. Make a symmetric matrix and add alpha*I to
   * its diagonal */

  create_symm_matrix(m, r);

  for (i = 0; i < N; ++i)
    {
      double mii = gsl_matrix_get(m, i, i);
      gsl_matrix_set(m, i, i, mii + alpha);
    }

  return GSL_SUCCESS;
}

static int
create_posdef_complex_matrix(gsl_matrix_complex *m, gsl_rng *r)
{
  const size_t N = m->size1;
  const double alpha = 10.0 * N;
  size_t i;

  create_herm_matrix(m, r);

  for (i = 0; i < N; ++i)
    {
      gsl_complex * mii = gsl_matrix_complex_ptr(m, i, i);
      GSL_REAL(*mii) += alpha;
    }

  return GSL_SUCCESS;
}

static int
create_hilbert_matrix2(gsl_matrix * m)
{
  const size_t N = m->size1;
  size_t i, j;

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          gsl_matrix_set(m, i, j, 1.0/(i+j+1.0));
        }
    }

  return GSL_SUCCESS;
}

/* create a matrix of a given rank */
static int
create_rank_matrix(const size_t rank, gsl_matrix * m, gsl_rng * r)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i;
  gsl_vector *u = gsl_vector_alloc(M);
  gsl_vector *v = gsl_vector_alloc(N);

  gsl_matrix_set_zero(m);

  /* add several rank-1 matrices together */
  for (i = 0; i < rank; ++i)
    {
      create_random_vector(u, r);
      create_random_vector(v, r);
      gsl_blas_dger(1.0, u, v, m);
    }

  gsl_vector_free(u);
  gsl_vector_free(v);

  return GSL_SUCCESS;
}

static int
create_posdef_band_matrix(const size_t p, gsl_matrix * m, gsl_rng * r)
{
  const size_t N = m->size1;
  const double alpha = 10.0 * N;
  size_t i;

  /* The idea is to make a symmetric diagonally dominant
   * matrix. Make a symmetric matrix and add alpha*I to
   * its diagonal */

  create_symm_band_matrix(p, m, r);

  for (i = 0; i < N; ++i)
    {
      double *mii = gsl_matrix_ptr(m, i, i);
      *mii += alpha;
    }

  return GSL_SUCCESS;
}

/* transform dense symmetric banded matrix to compact form, with bandwidth p */
static int
symm2band_matrix(const size_t p, const gsl_matrix * m, gsl_matrix * bm)
{
  const size_t N = m->size1;

  if (bm->size1 != N)
    {
      GSL_ERROR("banded matrix requires N rows", GSL_EBADLEN);
    }
  else if (bm->size2 != p + 1)
    {
      GSL_ERROR("banded matrix requires p + 1 columns", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      gsl_matrix_set_zero(bm);

      for (i = 0; i < p + 1; ++i)
        {
          gsl_vector_const_view diag = gsl_matrix_const_subdiagonal(m, i);
          gsl_vector_view v = gsl_matrix_subcolumn(bm, i, 0, N - i);

          gsl_vector_memcpy(&v.vector, &diag.vector);
        }

      return GSL_SUCCESS;
    }
}

/* transform general dense (p,q) banded matrix to compact banded form */
static int
gen2band_matrix(const size_t p, const size_t q, const gsl_matrix * A, gsl_matrix * AB)
{
  const size_t N = A->size2;

  if (AB->size1 != N)
    {
      GSL_ERROR("banded matrix requires N rows", GSL_EBADLEN);
    }
  else if (AB->size2 != 2*p + q + 1)
    {
      GSL_ERROR("banded matrix requires 2*p + q + 1 columns", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      gsl_matrix_set_zero(AB);

      /* copy diagonal and subdiagonals */
      for (i = 0; i <= p; ++i)
        {
          gsl_vector_const_view v = gsl_matrix_const_subdiagonal(A, i);
          gsl_vector_view w = gsl_matrix_subcolumn(AB, p + q + i, 0, v.vector.size);
          gsl_vector_memcpy(&w.vector, &v.vector);
        }

      /* copy superdiagonals */
      for (i = 1; i <= q; ++i)
        {
          gsl_vector_const_view v = gsl_matrix_const_superdiagonal(A, i);
          gsl_vector_view w = gsl_matrix_subcolumn(AB, p + q - i, i, v.vector.size);
          gsl_vector_memcpy(&w.vector, &v.vector);
        }

      return GSL_SUCCESS;
    }
}

/*
 * Checks of a computed solution x of A x = b against the exact solution
 * sol, where A or sol is random.  A backward stable solver only
 * guarantees the normwise bound
 *
 *   ||x - sol||_inf / ||sol||_inf <= c N eps kappa(A)
 *
 * (Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed.,
 * section 7.1); a component much smaller than ||sol|| can have a much
 * larger relative error, and random matrices are occasionally
 * ill-conditioned, so a fixed component-wise tolerance fails for a small
 * fraction of random systems, depending on rounding details of the
 * platform (bug #67705).  On 250000 random systems up to N = 50 the ratio of the
 * two sides with c = 1 stayed below 0.9 (real) and 3.1 (complex), so
 * the tests use c = TEST_SOLVE_C.
 */

#define TEST_SOLVE_C 64.0

/* kappa_1 estimate of a square matrix, from the R factor of its QR decomposition */
static double
test_cond(const gsl_matrix * A)
{
  const size_t N = A->size1;
  gsl_matrix * QR = gsl_matrix_alloc(N, N);
  gsl_matrix * T = gsl_matrix_alloc(N, N);
  gsl_vector * work = gsl_vector_alloc(3 * N);
  double rcond;

  gsl_matrix_memcpy(QR, A);
  gsl_linalg_QR_decomp_r(QR, T);
  gsl_linalg_tri_upper_rcond(QR, &rcond, work);

  gsl_matrix_free(QR);
  gsl_matrix_free(T);
  gsl_vector_free(work);

  return 1.0 / rcond;
}

/* the same for complex A, through its real form [Re -Im; Im Re] */
static double
test_cond_complex(const gsl_matrix_complex * A)
{
  const size_t N = A->size1;
  gsl_matrix * B = gsl_matrix_alloc(2 * N, 2 * N);
  double kappa;
  size_t i, j;

  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          gsl_complex z = gsl_matrix_complex_get(A, i, j);
          gsl_matrix_set(B, i, j, GSL_REAL(z));
          gsl_matrix_set(B, i, N + j, -GSL_IMAG(z));
          gsl_matrix_set(B, N + i, j, GSL_IMAG(z));
          gsl_matrix_set(B, N + i, N + j, GSL_REAL(z));
        }
    }

  kappa = test_cond(B);
  gsl_matrix_free(B);

  return kappa;
}

/* normwise check of x against sol with relative tolerance tol (what: "" or a
   qualifier for the description, such as " refined") */
static int
test_solution_normwise(const gsl_vector * x, const gsl_vector * sol, const double tol,
                       const char * desc, const char * what)
{
  const size_t N = sol->size;
  double err = 0.0, nrm = 0.0;
  size_t i;
  int s;

  for (i = 0; i < N; ++i)
    {
      double d = fabs(gsl_vector_get(x, i) - gsl_vector_get(sol, i));
      double a = fabs(gsl_vector_get(sol, i));
      if (!(d <= err))
        err = d;              /* also propagates NaN */
      if (a > nrm)
        nrm = a;
    }

  s = !(err <= tol * nrm);
  gsl_test(s, "%s%s (%3lu): normwise relative error %g, tolerance %g",
           desc, what, N, err / nrm, tol);

  return s;
}

static int
test_solution_normwise_complex(const gsl_vector_complex * x, const gsl_vector_complex * sol,
                               const double tol, const char * desc, const char * what)
{
  const size_t N = sol->size;
  double err = 0.0, nrm = 0.0;
  size_t i;
  int s;

  for (i = 0; i < N; ++i)
    {
      gsl_complex d = gsl_complex_sub(gsl_vector_complex_get(x, i), gsl_vector_complex_get(sol, i));
      double e = gsl_complex_abs(d);
      double a = gsl_complex_abs(gsl_vector_complex_get(sol, i));
      if (!(e <= err))
        err = e;
      if (a > nrm)
        nrm = a;
    }

  s = !(err <= tol * nrm);
  gsl_test(s, "%s%s (%3lu): normwise relative error %g, tolerance %g",
           desc, what, N, err / nrm, tol);

  return s;
}
