#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

double f (double x, void * params)
{
  (void)(params); /* avoid unused parameter warning */
  return exp(x);
}

int
main (void)
{
  gsl_deriv_parameters params;
  gsl_function F;
  const double x = 1.0;
  double result, abserr, fs;
  unsigned int neval;

  F.function = &f;
  F.params = NULL;

  printf ("f(x) = exp(x)\n");

  params = gsl_deriv_default_parameters(1, 4);
  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, NULL);
  printf ("x = %g, h_final = %.12e, neval = %u\n", x, fs, neval);
  printf ("f^{(%u)}(x) = %.12e +/- %.12e\n", params.deriv_order, result, abserr);
  printf ("exact      = %.12e\n\n", exp(x));

  params = gsl_deriv_default_parameters(10, 4);
  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, NULL);
  printf ("x = %g, h_final = %.12e, neval = %u\n", x, fs, neval);
  printf ("f^{(%u)}(x) = %.12e +/- %.12e\n", params.deriv_order, result, abserr);
  printf ("exact       = %.12e\n\n", exp(x));

  return 0;
}
