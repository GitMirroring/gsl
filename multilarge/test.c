/* multilarge/test.c
 * 
 * Copyright (C) 2015 Patrick Alken
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_ieee_utils.h>

static void test_random_matrix_orth(gsl_matrix *m, const gsl_rng *r);
static void test_random_matrix_ill(gsl_matrix *m, const gsl_rng *r);
static void test_random_vector(gsl_vector *v, const gsl_rng *r,
                               const double lower, const double upper);
static void test_random_matrix(gsl_matrix *m, const gsl_rng *r,
                               const double lower, const double upper);
static void test_random_vector_noise(const gsl_rng *r, gsl_vector *y);
static void test_compare_vectors(const double tol, const gsl_vector * a,
                                 const gsl_vector * b, const char * desc);
static void test_multifit_solve(const double lambda, const gsl_matrix * X,
                                const gsl_vector * y, const gsl_vector * wts,
                                const gsl_vector * diagL, const gsl_matrix * L,
                                double *rnorm, double *snorm, gsl_vector * c);
static void test_multilarge_solve(const gsl_multilarge_linear_type * T, const size_t rows_per_block, const double lambda,
                                  const gsl_matrix * X, const gsl_vector * y, gsl_vector * wts,
                                  const gsl_vector * diagL, const gsl_matrix * L,
                                  double *rnorm, double *snorm, gsl_vector * c);

/* generate random square orthogonal matrix via QR decomposition */
static void
test_random_matrix_orth(gsl_matrix *m, const gsl_rng *r)
{
  const size_t M = m->size1;
  gsl_matrix *A = gsl_matrix_alloc(M, M);
  gsl_vector *tau = gsl_vector_alloc(M);
  gsl_matrix *R = gsl_matrix_alloc(M, M);

  test_random_matrix(A, r, -1.0, 1.0);
  gsl_linalg_QR_decomp(A, tau);
  gsl_linalg_QR_unpack(A, tau, m, R);

  gsl_matrix_free(A);
  gsl_matrix_free(R);
  gsl_vector_free(tau);
}

/* construct ill-conditioned matrix via SVD */
static void
test_random_matrix_ill(gsl_matrix *m, const gsl_rng *r)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  gsl_matrix *U = gsl_matrix_alloc(M, M);
  gsl_matrix *V = gsl_matrix_alloc(N, N);
  gsl_vector *S = gsl_vector_alloc(N);
  gsl_matrix_view Uv = gsl_matrix_submatrix(U, 0, 0, M, N);
  const double smin = 16.0 * GSL_DBL_EPSILON;
  const double smax = 10.0;
  const double ratio = pow(smin / smax, 1.0 / (N - 1.0));
  double s;
  size_t j;

  test_random_matrix_orth(U, r);
  test_random_matrix_orth(V, r);

  /* compute U * S */

  s = smax;
  for (j = 0; j < N; ++j)
    {
      gsl_vector_view uj = gsl_matrix_column(U, j);

      gsl_vector_scale(&uj.vector, s);
      s *= ratio;
    }

  /* compute m = (U * S) * V' */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &Uv.matrix, V, 0.0, m);

  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(S);
}

static void
test_random_vector(gsl_vector *v, const gsl_rng *r,
                   const double lower, const double upper)
{
  size_t i;
  size_t N = v->size;

  for (i = 0; i < N; ++i)
    {
      gsl_vector_set(v, i,
                     gsl_rng_uniform(r) * (upper - lower) + lower);
    }
}

static void
test_random_matrix(gsl_matrix *m, const gsl_rng *r,
                   const double lower, const double upper)
{
  size_t i, j;
  size_t M = m->size1;
  size_t N = m->size2;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
      {
        gsl_matrix_set(m, i, j,
                       gsl_rng_uniform(r) * (upper - lower) + lower);
      }
    }
}

static void
test_random_vector_noise(const gsl_rng *r, gsl_vector *y)
{
  size_t i;

  for (i = 0; i < y->size; ++i)
    {
      double *ptr = gsl_vector_ptr(y, i);
      *ptr += 1.0e-3 * gsl_rng_uniform(r);
    }
}

static void
test_compare_vectors(const double tol, const gsl_vector * a,
                     const gsl_vector * b, const char * desc)
{
  size_t i;

  for (i = 0; i < a->size; ++i)
    {
      double ai = gsl_vector_get(a, i);
      double bi = gsl_vector_get(b, i);

      gsl_test_rel(bi, ai, tol, "%s i=%zu", desc, i);
    }
}

/* solve least squares system with multifit SVD */
static void
test_multifit_solve(const double lambda, const gsl_matrix * X,
                    const gsl_vector * y, const gsl_vector * wts,
                    const gsl_vector * diagL, const gsl_matrix * L,
                    double *rnorm, double *snorm, gsl_vector * c)
{
  const size_t n = X->size1;
  const size_t p = X->size2;
  gsl_multifit_linear_workspace *w =
    gsl_multifit_linear_alloc(n, p);
  gsl_matrix *Xs = gsl_matrix_alloc(n, p);
  gsl_vector *ys = gsl_vector_alloc(n);
  gsl_vector *cs = gsl_vector_alloc(p);
  gsl_matrix *LQR = NULL;
  gsl_vector *Ltau = NULL;
  gsl_matrix *M = NULL;

  /* convert to standard form */
  if (diagL)
    {
      gsl_multifit_linear_wstdform1(diagL, X, wts, y, Xs, ys, w);
    }
  else if (L)
    {
      const size_t m = L->size1;

      LQR = gsl_matrix_alloc(m, p);
      Ltau = gsl_vector_alloc(GSL_MIN(m, p));
      M = (m >= p) ? gsl_matrix_alloc(m, p) : gsl_matrix_alloc(n, p);

      gsl_matrix_memcpy(LQR, L);
      gsl_multifit_linear_L_decomp(LQR, Ltau);
      gsl_multifit_linear_wstdform2(LQR, Ltau, X, wts, y, Xs, ys, M, w);
    }
  else
    {
      gsl_matrix_memcpy(Xs, X);
      gsl_vector_memcpy(ys, y);
    }

  gsl_multifit_linear_svd(Xs, w);
  gsl_multifit_linear_solve(lambda, Xs, ys, cs, rnorm, snorm, w);

  /* convert to general form */
  if (diagL)
    gsl_multifit_linear_genform1(diagL, cs, c, w);
  else if (L)
    gsl_multifit_linear_wgenform2(LQR, Ltau, X, wts, y, cs, M, c, w);
  else
    gsl_vector_memcpy(c, cs);

  gsl_multifit_linear_free(w);
  gsl_matrix_free(Xs);
  gsl_vector_free(ys);
  gsl_vector_free(cs);

  if (LQR)
    gsl_matrix_free(LQR);
  if (Ltau)
    gsl_vector_free(Ltau);
  if (M)
    gsl_matrix_free(M);
}

/* solve least squares system with multilarge */
static void
test_multilarge_solve(const gsl_multilarge_linear_type * T, const size_t rows_per_block, const double lambda,
                      const gsl_matrix * X, const gsl_vector * y, gsl_vector * wts,
                      const gsl_vector * diagL, const gsl_matrix * L,
                      double *rnorm, double *snorm, gsl_vector * c)
{
  const size_t n = X->size1;
  const size_t p = X->size2;
  const size_t nrows = GSL_MIN(rows_per_block, n);
  gsl_multilarge_linear_workspace *w =
    gsl_multilarge_linear_alloc(T, p);
  gsl_matrix *Xs = gsl_matrix_alloc(nrows, p);
  gsl_vector *ys = gsl_vector_alloc(nrows);
  gsl_vector *cs = gsl_vector_alloc(p);
  gsl_matrix *LQR = NULL;
  gsl_vector *Ltau = NULL;
  size_t rowidx = 0;

  if (L)
    {
      const size_t m = L->size1;

      LQR = gsl_matrix_alloc(m, p);
      Ltau = gsl_vector_alloc(p);

      gsl_matrix_memcpy(LQR, L);
      gsl_multilarge_linear_L_decomp(LQR, Ltau);
    }

  while (rowidx < n)
    {
      size_t nleft = n - rowidx;
      size_t nr = GSL_MIN(nrows, nleft);
      gsl_matrix_const_view Xv = gsl_matrix_const_submatrix(X, rowidx, 0, nr, p);
      gsl_vector_const_view yv = gsl_vector_const_subvector(y, rowidx, nr);
      gsl_vector_view wv;
      gsl_matrix_view Xsv = gsl_matrix_submatrix(Xs, 0, 0, nr, p);
      gsl_vector_view ysv = gsl_vector_subvector(ys, 0, nr);

      if (wts)
        wv = gsl_vector_subvector(wts, rowidx, nr);

      /* convert to standard form */
      if (diagL)
        {
          gsl_multilarge_linear_wstdform1(diagL, &Xv.matrix, wts ? &wv.vector : NULL,
                                          &yv.vector, &Xsv.matrix, &ysv.vector, w);
        }
      else if (L)
        {
          gsl_multilarge_linear_wstdform2(LQR, Ltau, &Xv.matrix, wts ? &wv.vector : NULL,
                                          &yv.vector, &Xsv.matrix, &ysv.vector, w);
        }
      else
        {
          gsl_matrix_memcpy(&Xsv.matrix, &Xv.matrix);
          gsl_vector_memcpy(&ysv.vector, &yv.vector);
        }

      gsl_multilarge_linear_accumulate(&Xsv.matrix, &ysv.vector, w);

      rowidx += nr;
    }

  gsl_multilarge_linear_solve(lambda, cs, rnorm, snorm, w);

  if (diagL)
    gsl_multilarge_linear_genform1(diagL, cs, c, w);
  else if (L)
    gsl_multilarge_linear_genform2(LQR, Ltau, cs, c, w);
  else
    gsl_vector_memcpy(c, cs);

  gsl_multilarge_linear_free(w);
  gsl_matrix_free(Xs);
  gsl_vector_free(ys);
  gsl_vector_free(cs);

  if (LQR)
    gsl_matrix_free(LQR);
  if (Ltau)
    gsl_vector_free(Ltau);
}

static void
test_random(const gsl_multilarge_linear_type * T,
            const size_t n, const size_t p,
            const size_t rows_per_block, const double tol,
            const gsl_rng * r)
{
  const double tol1 = 1.0e3 * tol;
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *c = gsl_vector_alloc(p);
  gsl_vector *w = gsl_vector_alloc(n);
  gsl_vector *diagL = gsl_vector_alloc(p);
  gsl_matrix *Lsqr = gsl_matrix_alloc(p, p);
  gsl_matrix *Ltall = gsl_matrix_alloc(5*p, p);
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
  double rnorm0, snorm0;
  double rnorm1, snorm1;
  char str[2048];
  size_t i;

  /* generate LS system */
  test_random_matrix(X, r, -1.0, 1.0);
  test_random_vector(c, r, -1.0, 1.0);

  /* compute y = X c + noise */
  gsl_blas_dgemv(CblasNoTrans, 1.0, X, c, 0.0, y);
  test_random_vector_noise(r, y);

  /* random weights */
  test_random_vector(w, r, 0.0, 1.0);

  /* random diag(L) */
  test_random_vector(diagL, r, 1.0, 5.0);

  /* random square L */
  test_random_matrix(Lsqr, r, -5.0, 5.0);

  /* random tall L */
  test_random_matrix(Ltall, r, -10.0, 10.0);

  for (i = 0; i < 5; ++i)
    {
      double lambda = i == 0 ? 0.0 : pow(10.0, -(double) i);

      /* unweighted with L = I */
      {
        test_multifit_solve(lambda, X, y, NULL, NULL, NULL, &rnorm0, &snorm0, c0);
        test_multilarge_solve(T, rows_per_block, lambda, X, y, NULL, NULL, NULL, &rnorm1, &snorm1, c1);

        sprintf(str, "random %s unweighted stdform n=%zu p=%zu lambda=%g",
                T->name, n, p, lambda);
        test_compare_vectors(tol, c0, c1, str);

        gsl_test_rel(rnorm1, rnorm0, tol1, "rnorm %s", str);
        gsl_test_rel(snorm1, snorm0, tol, "snorm %s", str);
      }

      /* weighted, L = diag(L) */
      {
        test_multifit_solve(lambda, X, y, w, diagL, NULL, &rnorm0, &snorm0, c0);
        test_multilarge_solve(T, rows_per_block, lambda, X, y, w, diagL, NULL, &rnorm1, &snorm1, c1);

        sprintf(str, "random %s weighted diag(L) n=%zu p=%zu lambda=%g",
                T->name, n, p, lambda);
        test_compare_vectors(tol, c0, c1, str);

        gsl_test_rel(rnorm1, rnorm0, tol1, "rnorm %s", str);
        gsl_test_rel(snorm1, snorm0, tol, "snorm %s", str);
      }

      /* unweighted, L = diag(L) */
      {
        test_multifit_solve(lambda, X, y, NULL, diagL, NULL, &rnorm0, &snorm0, c0);
        test_multilarge_solve(T, rows_per_block, lambda, X, y, NULL, diagL, NULL, &rnorm1, &snorm1, c1);

        sprintf(str, "random %s unweighted diag(L) n=%zu p=%zu lambda=%g",
                T->name, n, p, lambda);
        test_compare_vectors(tol, c0, c1, str);

        gsl_test_rel(rnorm1, rnorm0, tol1, "rnorm %s", str);
        gsl_test_rel(snorm1, snorm0, tol, "snorm %s", str);
      }

      /* weighted, L = square */
      {
        test_multifit_solve(lambda, X, y, w, NULL, Lsqr, &rnorm0, &snorm0, c0);
        test_multilarge_solve(T, rows_per_block, lambda, X, y, w, NULL, Lsqr, &rnorm1, &snorm1, c1);

        sprintf(str, "random %s weighted Lsqr n=%zu p=%zu lambda=%g",
                T->name, n, p, lambda);
        test_compare_vectors(tol, c0, c1, str);

        gsl_test_rel(rnorm1, rnorm0, tol1, "rnorm %s", str);
        gsl_test_rel(snorm1, snorm0, tol, "snorm %s", str);
      }

      /* unweighted, L = square */
      {
        test_multifit_solve(lambda, X, y, NULL, NULL, Lsqr, &rnorm0, &snorm0, c0);
        test_multilarge_solve(T, rows_per_block, lambda, X, y, NULL, NULL, Lsqr, &rnorm1, &snorm1, c1);

        sprintf(str, "random %s unweighted Lsqr n=%zu p=%zu lambda=%g",
                T->name, n, p, lambda);
        test_compare_vectors(tol, c0, c1, str);

        gsl_test_rel(rnorm1, rnorm0, tol1, "rnorm %s", str);
        gsl_test_rel(snorm1, snorm0, tol, "snorm %s", str);
      }

      /* weighted, L = tall */
      {
        test_multifit_solve(lambda, X, y, w, NULL, Ltall, &rnorm0, &snorm0, c0);
        test_multilarge_solve(T, rows_per_block, lambda, X, y, w, NULL, Ltall, &rnorm1, &snorm1, c1);

        sprintf(str, "random %s weighted Ltall n=%zu p=%zu lambda=%g",
                T->name, n, p, lambda);
        test_compare_vectors(tol, c0, c1, str);

        gsl_test_rel(rnorm1, rnorm0, tol1, "rnorm %s", str);
        gsl_test_rel(snorm1, snorm0, tol, "snorm %s", str);
      }

      /* unweighted, L = tall */
      {
        test_multifit_solve(lambda, X, y, NULL, NULL, Ltall, &rnorm0, &snorm0, c0);
        test_multilarge_solve(T, rows_per_block, lambda, X, y, NULL, NULL, Ltall, &rnorm1, &snorm1, c1);

        sprintf(str, "random %s unweighted Ltall n=%zu p=%zu lambda=%g",
                T->name, n, p, lambda);
        test_compare_vectors(tol, c0, c1, str);

        gsl_test_rel(rnorm1, rnorm0, tol1, "rnorm %s", str);
        gsl_test_rel(snorm1, snorm0, tol, "snorm %s", str);
      }
    }

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_vector_free(w);
  gsl_vector_free(diagL);
  gsl_matrix_free(Lsqr);
  gsl_matrix_free(Ltall);
  gsl_vector_free(c0);
  gsl_vector_free(c1);
}

/* accumulate rows [first, first + nrows) of (X,y) into w, in blocks of size rows_per_block */
static void
test_accumulate_rows(const size_t first, const size_t nrows, const size_t rows_per_block,
                     const gsl_matrix * X, const gsl_vector * y,
                     gsl_multilarge_linear_workspace * w)
{
  const size_t p = X->size2;
  gsl_matrix *Xs = gsl_matrix_alloc(rows_per_block, p);
  gsl_vector *ys = gsl_vector_alloc(rows_per_block);
  size_t rowidx = first;

  while (rowidx < first + nrows)
    {
      size_t nr = GSL_MIN(rows_per_block, first + nrows - rowidx);
      gsl_matrix_const_view Xv = gsl_matrix_const_submatrix(X, rowidx, 0, nr, p);
      gsl_vector_const_view yv = gsl_vector_const_subvector(y, rowidx, nr);
      gsl_matrix_view Xsv = gsl_matrix_submatrix(Xs, 0, 0, nr, p);
      gsl_vector_view ysv = gsl_vector_subvector(ys, 0, nr);

      /* accumulate() destroys its inputs, so pass copies */
      gsl_matrix_memcpy(&Xsv.matrix, &Xv.matrix);
      gsl_vector_memcpy(&ysv.vector, &yv.vector);
      gsl_multilarge_linear_accumulate(&Xsv.matrix, &ysv.vector, w);

      rowidx += nr;
    }

  gsl_matrix_free(Xs);
  gsl_vector_free(ys);
}

/* check that the strict lower triangle of the TSQR R factor is zero */
static void
test_tsqr_lower_zero(const gsl_multilarge_linear_workspace * w, const char * desc)
{
  const gsl_matrix * R = gsl_multilarge_linear_matrix_ptr(w);
  size_t i, j;
  int nonzero = 0;

  for (i = 1; i < R->size1; ++i)
    {
      for (j = 0; j < i; ++j)
        nonzero += (gsl_matrix_get(R, i, j) != 0.0);
    }

  gsl_test(nonzero, "tsqr %s lower triangle of R is zero", desc);
}

/* solve with w and compare against multifit SVD solution using the first n rows of (X,y) */
static void
test_tsqr_compare(const double lambda, const size_t n,
                  const gsl_matrix * X, const gsl_vector * y,
                  gsl_multilarge_linear_workspace * w,
                  const double tol, const char * desc)
{
  const size_t p = X->size2;
  gsl_matrix_const_view Xv = gsl_matrix_const_submatrix(X, 0, 0, n, p);
  gsl_vector_const_view yv = gsl_vector_const_subvector(y, 0, n);
  gsl_vector *c0 = gsl_vector_alloc(p);
  gsl_vector *c1 = gsl_vector_alloc(p);
  double rnorm0, snorm0, rnorm1, snorm1;
  char str[2048];
  int status;

  sprintf(str, "tsqr %s n=%zu p=%zu lambda=%g", desc, n, p, lambda);

  test_multifit_solve(lambda, &Xv.matrix, &yv.vector, NULL, NULL, NULL, &rnorm0, &snorm0, c0);
  status = gsl_multilarge_linear_solve(lambda, c1, &rnorm1, &snorm1, w);

  gsl_test(status, "%s solve status", str);
  test_compare_vectors(tol, c0, c1, str);
  gsl_test_rel(rnorm1, rnorm0, 1.0e3 * tol, "rnorm %s", str);
  gsl_test_rel(snorm1, snorm0, tol, "snorm %s", str);

  gsl_vector_free(c0);
  gsl_vector_free(c1);
}

/* test TSQR handling of blocks with fewer than p rows, reset, and
 * re-solving after additional rows are accumulated */
static void
test_tsqr_blocks(const size_t n, const size_t p, const double tol, const gsl_rng * r)
{
  const double lambda = 0.1;
  const size_t nlcurve = 50;
  const size_t n1 = n / 2; /* rows accumulated before the first solve */
  gsl_matrix *X = gsl_matrix_alloc(n, p);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *c = gsl_vector_alloc(p);
  gsl_vector *reg0 = gsl_vector_alloc(nlcurve);
  gsl_vector *rho0 = gsl_vector_alloc(nlcurve);
  gsl_vector *eta0 = gsl_vector_alloc(nlcurve);
  gsl_vector *reg1 = gsl_vector_alloc(nlcurve);
  gsl_vector *rho1 = gsl_vector_alloc(nlcurve);
  gsl_vector *eta1 = gsl_vector_alloc(nlcurve);
  gsl_multilarge_linear_workspace *w0 = gsl_multilarge_linear_alloc(gsl_multilarge_linear_tsqr, p);
  gsl_multilarge_linear_workspace *w = gsl_multilarge_linear_alloc(gsl_multilarge_linear_tsqr, p);
  gsl_error_handler_t *old_handler;
  double rnorm, snorm, rcond1;
  char str[2048];
  size_t i;
  int status;

  test_random_matrix(X, r, -1.0, 1.0);
  test_random_vector(y, r, -1.0, 1.0);

  /* reference workspace: all rows accumulated in a single block */
  test_accumulate_rows(0, n, n, X, y, w0);

  old_handler = gsl_set_error_handler_off();

  /* fewer than p rows accumulated: solve, rcond and lcurve must fail */
  for (i = 0; i < p - 1; ++i)
    {
      test_accumulate_rows(i, 1, 1, X, y, w);

      status = gsl_multilarge_linear_solve(lambda, c, &rnorm, &snorm, w);
      gsl_test(status != GSL_EINVAL, "tsqr solve with %zu/%zu rows returns GSL_EINVAL", i + 1, p);

      status = gsl_multilarge_linear_rcond(&rcond1, w);
      gsl_test(status != GSL_EINVAL, "tsqr rcond with %zu/%zu rows returns GSL_EINVAL", i + 1, p);

      status = gsl_multilarge_linear_lcurve(reg1, rho1, eta1, w);
      gsl_test(status != GSL_EINVAL, "tsqr lcurve with %zu/%zu rows returns GSL_EINVAL", i + 1, p);
    }

  /* remaining rows up to n1, one at a time */
  test_accumulate_rows(p - 1, n1 - (p - 1), 1, X, y, w);
  test_tsqr_lower_zero(w, "single rows");
  test_tsqr_compare(0.0, n1, X, y, w, tol, "single rows");
  test_tsqr_compare(lambda, n1, X, y, w, tol, "single rows");

  /* add the rest of the rows; the SVD of R computed above must not be reused */
  test_accumulate_rows(n1, n - n1, 3, X, y, w);
  test_tsqr_compare(lambda, n, X, y, w, tol, "solve after additional rows");

  sprintf(str, "tsqr small blocks n=%zu p=%zu", n, p);

  /* R is unique up to the signs of its rows, so compare |R| with the single-block reference */
  {
    const gsl_matrix * R0 = gsl_multilarge_linear_matrix_ptr(w0);
    const gsl_matrix * R1 = gsl_multilarge_linear_matrix_ptr(w);
    size_t j;

    for (i = 0; i < p; ++i)
      {
        for (j = i; j < p; ++j)
          {
            gsl_test_abs(fabs(gsl_matrix_get(R1, i, j)), fabs(gsl_matrix_get(R0, i, j)),
                         tol * fabs(gsl_matrix_get(R0, i, i)), "%s |R(%zu,%zu)|", str, i, j);
          }
      }
  }

  /* the rcond estimate depends on the row signs of R, so only check it is valid */
  status = gsl_multilarge_linear_rcond(&rcond1, w);
  gsl_test(status, "%s rcond status", str);
  gsl_test(!(rcond1 > 0.0 && rcond1 <= 1.0), "%s rcond in (0,1] (%g)", str, rcond1);

  /* lcurve depends only on the singular values of R and must match the reference */

  status = gsl_multilarge_linear_lcurve(reg0, rho0, eta0, w0);
  gsl_test(status, "%s reference lcurve status", str);
  status = gsl_multilarge_linear_lcurve(reg1, rho1, eta1, w);
  gsl_test(status, "%s lcurve status", str);
  test_compare_vectors(tol, reg0, reg1, "tsqr lcurve reg_param");
  test_compare_vectors(1.0e3 * tol, rho0, rho1, "tsqr lcurve rho");
  test_compare_vectors(tol, eta0, eta1, "tsqr lcurve eta");

  /* reset and accumulate again in small blocks */
  gsl_multilarge_linear_reset(w);

  test_accumulate_rows(0, 1, 1, X, y, w);
  status = gsl_multilarge_linear_solve(lambda, c, &rnorm, &snorm, w);
  gsl_test(status != GSL_EINVAL, "tsqr solve after reset with 1/%zu rows returns GSL_EINVAL", p);

  test_accumulate_rows(1, n - 1, 3, X, y, w);
  test_tsqr_lower_zero(w, "after reset");
  test_tsqr_compare(0.0, n, X, y, w, tol, "after reset");
  test_tsqr_compare(lambda, n, X, y, w, tol, "after reset");

  gsl_set_error_handler(old_handler);

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_vector_free(reg0);
  gsl_vector_free(rho0);
  gsl_vector_free(eta0);
  gsl_vector_free(reg1);
  gsl_vector_free(rho1);
  gsl_vector_free(eta1);
  gsl_multilarge_linear_free(w0);
  gsl_multilarge_linear_free(w);
}

int
main (void)
{
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  gsl_ieee_env_setup();

  {
    const double tol1 = 1.0e-8;
    const double tol2 = 1.0e-11;
    const size_t n_vals[] = { 200, 356, 501 };
    const size_t p_vals[] = { 10, 21, 34 };
    const size_t rows_per_block[] = { 1, 3, 50 }; /* include blocks with fewer than p rows */
    size_t i, j;

    for (i = 0; i < sizeof(n_vals) / sizeof(n_vals[0]); ++i)
      {
        size_t n = n_vals[i];
        size_t p = p_vals[i];

        for (j = 0; j < sizeof(rows_per_block) / sizeof(rows_per_block[0]); ++j)
          {
            /* generate random ill-conditioned LS system and test */
            test_random(gsl_multilarge_linear_normal, n, p, rows_per_block[j], tol1, r);
            test_random(gsl_multilarge_linear_tsqr, n, p, rows_per_block[j], tol2, r);
          }

        test_tsqr_blocks(n, p, tol2, r);
      }
  }

  gsl_rng_free(r);

  exit (gsl_test_summary ());
}
