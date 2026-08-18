#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

double
f (const gsl_vector * x, void * params)
{
  double xnorm = gsl_blas_dnrm2(x);
  (void)(params); /* avoid unused parameter warning */
  return xnorm*xnorm;
}

double
df (const gsl_vector * x, const gsl_vector * v, void * params)
{
  double val;
  (void)(params); /* avoid unused parameter warning */
  gsl_blas_ddot(x, v, &val);
  return 2.0*val;
}

int
main (void)
{
  gsl_deriv_parameters params = gsl_deriv_default_parameters(1, 4);
  double x_data[] = {  0.5, -0.3, 2.0 };
  double v_data[] = { -0.9,  0.1, 0.5 };
  double work_data[3];
  gsl_vector_const_view x = gsl_vector_const_view_array(x_data, 3);
  gsl_vector_const_view v = gsl_vector_const_view_array(v_data, 3);
  gsl_vector_view work = gsl_vector_view_array(work_data, 3);
  gsl_deriv_function F;
  double result, abserr, fs;
  unsigned int neval;

  F.f = &f;
  F.n = 3;
  F.params = NULL;

  printf ("f(x) = |x|^2\n");

  gsl_deriv_directional(&params, &F, &x.vector, &v.vector,
                        &result, &abserr, &fs, &neval, NULL, &work.vector);
  printf ("h_final  = %.12e, neval = %u\n", fs, neval);
  printf ("D_v f(x) = %.12e +/- %.12e\n", result, abserr);
  printf ("exact    = %.12e\n\n", df(&x.vector, &v.vector, NULL));

  return 0;
}
