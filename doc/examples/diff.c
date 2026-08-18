#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

double f (double x, void * params)
{
  (void)(params); /* avoid unused parameter warning */
  return log(1.0 + x);
}

int
main (void)
{
  gsl_deriv_parameters params = gsl_deriv_default_parameters(1, 4);
  gsl_function F;
  double result, abserr, fs;
  unsigned int neval;

  F.function = &f;
  F.params = 0;

  printf ("f(x) = log(1+x)\n");

  gsl_deriv_adaptive(&params, &F, 5.0, &result, &abserr, &fs, &neval, NULL);
  printf ("x = 5.0, h_final = %.12e, neval = %u\n", fs, neval);
  printf ("f'(x) = %.12e +/- %.12e\n", result, abserr);
  printf ("exact = %.12e\n\n", 1.0/6.0);

  params.stencil_type = GSL_DERIV_STENCIL_FORWARD;
  gsl_deriv_adaptive(&params, &F, -0.9, &result, &abserr, &fs, &neval, NULL);
  printf ("x = -0.9, h_final = %.12e, neval = %u\n", fs, neval);
  printf ("f'(x) = %.12e +/- %.12e\n", result, abserr);
  printf ("exact = %.12e\n\n", 10.0);

  return 0;
}
