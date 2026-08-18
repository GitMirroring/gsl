#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

double
f (const gsl_vector * x, void * params)
{
  const size_t n = x->size;
  const gsl_matrix * A = (const gsl_matrix *) params;
  gsl_vector * work = gsl_vector_alloc(n);
  double val;

  gsl_blas_dsymv(CblasLower, 1.0, A, x, 0.0, work);
  gsl_blas_ddot(x, work, &val);

  gsl_vector_free(work);

  return val;
}

double
df (const gsl_vector * x, void * params, gsl_vector * g)
{
  const gsl_matrix * A = (const gsl_matrix *) params;
  gsl_blas_dsymv(CblasLower, 2.0, A, x, 0.0, g);
  return GSL_SUCCESS;
}

int
main (void)
{
  gsl_deriv_parameters params = gsl_deriv_default_parameters(1, 4);
  const size_t n = 3;
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_matrix * A = gsl_matrix_alloc(n, n);
  gsl_vector * g = gsl_vector_alloc(n);
  gsl_vector * result = gsl_vector_alloc(n);
  gsl_vector * abserr = gsl_vector_alloc(n);
  gsl_vector * fs = gsl_vector_alloc(n);
  gsl_vector * work = gsl_vector_alloc(n);
  gsl_deriv_function F;
  unsigned int neval;

  F.f = &f;
  F.n = 3;
  F.params = (void *) A;

  /* fill in x */
  gsl_vector_set(x, 0, -0.1);
  gsl_vector_set(x, 1,  0.9);
  gsl_vector_set(x, 2,  1.7);

  /* fill in A (lower triangle) */
  gsl_matrix_set(A, 0, 0,  1.0);
  gsl_matrix_set(A, 1, 0, -0.2);
  gsl_matrix_set(A, 2, 0,  2.5);
  gsl_matrix_set(A, 1, 1,  2.0);
  gsl_matrix_set(A, 2, 1, -1.1);
  gsl_matrix_set(A, 2, 2,  3.0);

  printf ("f(x) = x^T A x\n");

  gsl_deriv_gradient(&params, &F, x, result, abserr, fs, &neval, work);
  printf ("h_final    = [%.12e,%.12e,%.12e], neval = %u\n",
          gsl_vector_get(fs, 0),
          gsl_vector_get(fs, 1),
          gsl_vector_get(fs, 2),
          neval);
  printf ("grad(f(x)) = [%.12e,%.12e,%.12e] +/- [%.12e,%.12e,%.12e]\n",
          gsl_vector_get(result, 0),
          gsl_vector_get(result, 1),
          gsl_vector_get(result, 2),
          gsl_vector_get(abserr, 0),
          gsl_vector_get(abserr, 1),
          gsl_vector_get(abserr, 2));

  df(x, A, g);
  printf ("exact      = [%.12e,%.12e,%.12e]\n\n",
          gsl_vector_get(g, 0),
          gsl_vector_get(g, 1),
          gsl_vector_get(g, 2));

  gsl_vector_free(x);
  gsl_matrix_free(A);
  gsl_vector_free(g);
  gsl_vector_free(result);
  gsl_vector_free(abserr);
  gsl_vector_free(fs);
  gsl_vector_free(work);

  return 0;
}
