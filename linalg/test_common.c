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

/*
 * Normwise comparison of a computed matrix X with a reference R:
 *
 *   max_ij |X_ij - R_ij| <= tol max_ij |R_ij|
 *
 * Decompositions are backward stable in this normwise sense (the
 * computed factors reproduce A up to c n eps ||A||, independent of the
 * conditioning of A), while entries much smaller than ||A|| can have a
 * large relative error; an entry-wise relative check with tolerance tol
 * implies this one.
 */
static int
test_matrix_normwise(const gsl_matrix * X, const gsl_matrix * R, const double tol,
                     const char * desc, const char * what)
{
  const size_t M = R->size1;
  const size_t N = R->size2;
  double err = 0.0, nrm = 0.0;
  size_t i, j;
  int s;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double d = fabs(gsl_matrix_get(X, i, j) - gsl_matrix_get(R, i, j));
          double a = fabs(gsl_matrix_get(R, i, j));
          if (!(d <= err))
            err = d;              /* also propagates NaN */
          if (a > nrm)
            nrm = a;
        }
    }

  s = !(err <= tol * (nrm > 0.0 ? nrm : 1.0));
  gsl_test(s, "%s%s (%3lu,%3lu): normwise relative error %g, tolerance %g",
           desc, what, M, N, err / (nrm > 0.0 ? nrm : 1.0), tol);

  return s;
}

static int
test_matrix_normwise_complex(const gsl_matrix_complex * X, const gsl_matrix_complex * R,
                             const double tol, const char * desc, const char * what)
{
  const size_t M = R->size1;
  const size_t N = R->size2;
  double err = 0.0, nrm = 0.0;
  size_t i, j;
  int s;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double d = gsl_complex_abs(gsl_complex_sub(gsl_matrix_complex_get(X, i, j),
                                                     gsl_matrix_complex_get(R, i, j)));
          double a = gsl_complex_abs(gsl_matrix_complex_get(R, i, j));
          if (!(d <= err))
            err = d;
          if (a > nrm)
            nrm = a;
        }
    }

  s = !(err <= tol * (nrm > 0.0 ? nrm : 1.0));
  gsl_test(s, "%s%s (%3lu,%3lu): normwise relative error %g, tolerance %g",
           desc, what, M, N, err / (nrm > 0.0 ? nrm : 1.0), tol);

  return s;
}

/* orthonormal columns: max_ij |(Q^T Q - I)_ij| <= tol */
static int
test_orthonormal(const gsl_matrix * Q, const double tol, const char * desc)
{
  const size_t N = Q->size2;
  gsl_matrix * QTQ = gsl_matrix_alloc(N, N);
  gsl_matrix * I = gsl_matrix_alloc(N, N);
  int s;

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Q, Q, 0.0, QTQ);
  gsl_matrix_set_identity(I);
  s = test_matrix_normwise(QTQ, I, tol, desc, " Q^T Q = I");

  gsl_matrix_free(QTQ);
  gsl_matrix_free(I);

  return s;
}

/* orthonormal columns: max_ij |(Q^H Q - I)_ij| <= tol */
static int
test_orthonormal_complex(const gsl_matrix_complex * Q, const double tol, const char * desc)
{
  const size_t N = Q->size2;
  gsl_matrix_complex * QHQ = gsl_matrix_complex_alloc(N, N);
  gsl_matrix_complex * I = gsl_matrix_complex_alloc(N, N);
  int s;

  gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, GSL_COMPLEX_ONE, Q, Q, GSL_COMPLEX_ZERO, QHQ);
  gsl_matrix_complex_set_identity(I);
  s = test_matrix_normwise_complex(QHQ, I, tol, desc, " Q^H Q = I");

  gsl_matrix_complex_free(QHQ);
  gsl_matrix_complex_free(I);

  return s;
}

/* max_ij |C_ij - delta_ij| <= tol, e.g. for C = A A^{-1}; for an inverse
   computed by a stable algorithm tol = c N eps kappa(A) */
static int
test_matrix_identity(const gsl_matrix * C, const double tol, const char * desc)
{
  gsl_matrix * I = gsl_matrix_alloc(C->size1, C->size2);
  int s;

  gsl_matrix_set_identity(I);
  s = test_matrix_normwise(C, I, tol, desc, " A A^{-1} = I");
  gsl_matrix_free(I);

  return s;
}

static int
test_matrix_identity_complex(const gsl_matrix_complex * C, const double tol, const char * desc)
{
  gsl_matrix_complex * I = gsl_matrix_complex_alloc(C->size1, C->size2);
  int s;

  gsl_matrix_complex_set_identity(I);
  s = test_matrix_normwise_complex(C, I, tol, desc, " A A^{-1} = I");
  gsl_matrix_complex_free(I);

  return s;
}

/*
 * Least squares: x minimizes ||b - A x|| iff the residual is orthogonal
 * to the columns of A.  For a backward stable solver
 *
 *   ||A^T (b - A x)||_inf <= tol ||A||_inf (||A||_inf ||x||_inf + ||b||_inf)
 *
 * with tol = c M eps, independent of the conditioning of A (a
 * comparison of x with the solution of another method, such as the
 * SVD, would have to allow for a sensitivity growing like kappa^2).
 */
static int
test_lssolve_normal(const gsl_matrix * A, const gsl_vector * b, const gsl_vector * x,
                    const double tol, const char * desc)
{
  const size_t M = A->size1;
  const size_t N = A->size2;
  gsl_vector * r = gsl_vector_alloc(M);
  gsl_vector * g = gsl_vector_alloc(N);
  double anorm = 0.0, scale, gnorm;
  size_t i;
  int s;

  for (i = 0; i < M; ++i)
    {
      gsl_vector_const_view row = gsl_matrix_const_row(A, i);
      double rs = gsl_blas_dasum(&row.vector);
      if (rs > anorm)
        anorm = rs;
    }

  gsl_vector_memcpy(r, b);
  gsl_blas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, r);
  gsl_blas_dgemv(CblasTrans, 1.0, A, r, 0.0, g);

  gnorm = fabs(gsl_vector_get(g, gsl_blas_idamax(g)));
  scale = anorm * (anorm * fabs(gsl_vector_get(x, gsl_blas_idamax(x)))
                   + fabs(gsl_vector_get(b, gsl_blas_idamax(b))));

  s = !(gnorm <= tol * scale);
  gsl_test(s, "%s (%3lu,%3lu): ||A^T r|| / (||A|| (||A|| ||x|| + ||b||)) = %g, tolerance %g",
           desc, M, N, gnorm / scale, tol);

  gsl_vector_free(r);
  gsl_vector_free(g);

  return s;
}

/* a residual norm computed by a solver, against ||b - A x|| recomputed
   here: both are accurate to about eps (||A|| ||x|| + ||b||) in absolute
   terms, not relative to the residual, which may be tiny */
static int
test_rnorm(const double norm, const gsl_matrix * A, const gsl_vector * x, const gsl_vector * b,
           const double tol, const char * desc)
{
  const size_t M = A->size1;
  const size_t N = A->size2;
  gsl_vector * r = gsl_vector_alloc(M);
  double anorm = 0.0, scale, norm_expected;
  size_t i;
  int s;

  for (i = 0; i < M; ++i)
    {
      gsl_vector_const_view row = gsl_matrix_const_row(A, i);
      double rs = gsl_blas_dasum(&row.vector);
      if (rs > anorm)
        anorm = rs;
    }

  gsl_vector_memcpy(r, b);
  gsl_blas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, r);
  norm_expected = gsl_blas_dnrm2(r);

  scale = anorm * fabs(gsl_vector_get(x, gsl_blas_idamax(x)))
          + fabs(gsl_vector_get(b, gsl_blas_idamax(b)));

  s = !(fabs(norm - norm_expected) <= tol * scale);
  gsl_test(s, "%s rnorm (%3lu,%3lu): %.18g, expected %.18g, tolerance %g",
           desc, M, N, norm, norm_expected, tol * scale);

  gsl_vector_free(r);

  return s;
}

/* as test_matrix_normwise() for the upper triangle (j >= i) only */
static int
test_matrix_normwise_upper(const gsl_matrix * X, const gsl_matrix * R, const double tol,
                           const char * desc, const char * what)
{
  const size_t M = GSL_MIN(R->size1, R->size2);
  const size_t N = R->size2;
  double err = 0.0, nrm = 0.0;
  size_t i, j;
  int s;

  for (i = 0; i < M; ++i)
    {
      for (j = i; j < N; ++j)
        {
          double d = fabs(gsl_matrix_get(X, i, j) - gsl_matrix_get(R, i, j));
          double a = fabs(gsl_matrix_get(R, i, j));
          if (!(d <= err))
            err = d;
          if (a > nrm)
            nrm = a;
        }
    }

  s = !(err <= tol * (nrm > 0.0 ? nrm : 1.0));
  gsl_test(s, "%s%s (%3lu,%3lu) upper triangle: normwise relative error %g, tolerance %g",
           desc, what, M, N, err / (nrm > 0.0 ? nrm : 1.0), tol);

  return s;
}

/* a residual vector computed by a solver, against b - A x recomputed
   here, on the absolute scale of test_rnorm() */
static int
test_residual(const gsl_vector * res, const gsl_matrix * A, const gsl_vector * x, const gsl_vector * b,
              const double tol, const char * desc)
{
  const size_t M = A->size1;
  const size_t N = A->size2;
  gsl_vector * r = gsl_vector_alloc(M);
  double anorm = 0.0, scale, err = 0.0;
  size_t i;
  int s;

  for (i = 0; i < M; ++i)
    {
      gsl_vector_const_view row = gsl_matrix_const_row(A, i);
      double rs = gsl_blas_dasum(&row.vector);
      if (rs > anorm)
        anorm = rs;
    }

  gsl_vector_memcpy(r, b);
  gsl_blas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, r);

  for (i = 0; i < M; ++i)
    {
      double d = fabs(gsl_vector_get(res, i) - gsl_vector_get(r, i));
      if (!(d <= err))
        err = d;
    }

  scale = anorm * fabs(gsl_vector_get(x, gsl_blas_idamax(x)))
          + fabs(gsl_vector_get(b, gsl_blas_idamax(b)));

  s = !(err <= tol * scale);
  gsl_test(s, "%s residual (%3lu,%3lu): error %g, tolerance %g",
           desc, M, N, err, tol * scale);

  gsl_vector_free(r);

  return s;
}
