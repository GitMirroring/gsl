#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

int
f(const gsl_vector * x, void * params, gsl_vector * y)
{
  const double x1 = gsl_vector_get(x, 0);
  const double x2 = gsl_vector_get(x, 1);
  (void) params;
  gsl_vector_set(y, 0, x1*x1);
  gsl_vector_set(y, 1, x1*x2);
  gsl_vector_set(y, 2, x2*x2);
  return GSL_SUCCESS;
}

int
df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  const double x1 = gsl_vector_get(x, 0);
  const double x2 = gsl_vector_get(x, 1);
  (void) params;
  gsl_matrix_set(J, 0, 0, 2.0*x1);
  gsl_matrix_set(J, 0, 1, 0.0);
  gsl_matrix_set(J, 1, 0, x2);
  gsl_matrix_set(J, 1, 1, x1);
  gsl_matrix_set(J, 2, 0, 0.0);
  gsl_matrix_set(J, 2, 1, 2.0*x2);
  return GSL_SUCCESS;
}

int
main (void)
{
  gsl_deriv_parameters params = gsl_deriv_default_parameters(1, 4);
  const size_t m = 3;
  const size_t n = 2;
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_matrix * result = gsl_matrix_alloc(m, n);
  gsl_matrix * abserr = gsl_matrix_alloc(m, n);
  gsl_matrix * fs = gsl_matrix_alloc(m, n);
  gsl_matrix * J = gsl_matrix_alloc(m, n);
  gsl_deriv_jacobian_workspace * w = gsl_deriv_jacobian_alloc(&params, m, n);
  gsl_deriv_vec_function F;
  unsigned int neval;

  F.f = &f;
  F.n = n;
  F.m = m;
  F.params = NULL;

  /* fill in x */
  gsl_vector_set(x, 0,  1.5);
  gsl_vector_set(x, 1, -0.3);

  printf ("F(x) = ( x1*x1, x1*x2, x2*x2 )\n");

  gsl_deriv_jacobian(&F, x, result, abserr, fs, &neval, w);
  printf ("neval = %u\n", neval);
  printf ("h_final = [%20.12e,%20.12e]\n"
          "          [%20.12e,%20.12e]\n"
          "          [%20.12e,%20.12e]\n\n",
          gsl_matrix_get(fs, 0, 0), gsl_matrix_get(fs, 0, 1),
          gsl_matrix_get(fs, 1, 0), gsl_matrix_get(fs, 1, 1),
          gsl_matrix_get(fs, 2, 0), gsl_matrix_get(fs, 2, 1));

  printf ("J =     [%20.12e,%20.12e] +/- [%20.12e,%20.12e]\n"
          "        [%20.12e,%20.12e]     [%20.12e,%20.12e]\n"
          "        [%20.12e,%20.12e]     [%20.12e,%20.12e]\n\n",
          gsl_matrix_get(result, 0, 0), gsl_matrix_get(result, 0, 1),
          gsl_matrix_get(abserr, 0, 0), gsl_matrix_get(abserr, 0, 1),
          gsl_matrix_get(result, 1, 0), gsl_matrix_get(result, 1, 1),
          gsl_matrix_get(abserr, 1, 0), gsl_matrix_get(abserr, 1, 1),
          gsl_matrix_get(result, 2, 0), gsl_matrix_get(result, 2, 1),
          gsl_matrix_get(abserr, 2, 0), gsl_matrix_get(abserr, 2, 1));

  df(x, NULL, J);
  printf ("exact = [%20.12e,%20.12e]\n"
          "        [%20.12e,%20.12e]\n"
          "        [%20.12e,%20.12e]\n",
          gsl_matrix_get(J, 0, 0), gsl_matrix_get(J, 0, 1),
          gsl_matrix_get(J, 1, 0), gsl_matrix_get(J, 1, 1),
          gsl_matrix_get(J, 2, 0), gsl_matrix_get(J, 2, 1));

  gsl_vector_free(x);
  gsl_matrix_free(result);
  gsl_matrix_free(abserr);
  gsl_matrix_free(fs);
  gsl_matrix_free(J);
  gsl_deriv_jacobian_free(w);

  return 0;
}
