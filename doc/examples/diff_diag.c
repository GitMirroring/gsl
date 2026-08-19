#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

double
f (double x, void * params)
{
  (void)(params); /* avoid unused parameter warning */
  return pow(x, 1.5);
}

double df(double x) { return 1.5*sqrt(x); }

int
main (int argc, char * argv[])
{
  gsl_deriv_parameters params = gsl_deriv_default_parameters(1, 4);
  gsl_function F;
  double result, abserr, fs;
  unsigned int neval;
  double x = 0.0;
  gsl_deriv_diag diag;

  if (argc > 1)
    x = atof(argv[1]);
  if (argc > 2)
    params.min_step = atof(argv[2]);
  if (argc > 3)
    params.max_step = atof(argv[3]);

  /* use forward stencil */
  params.stencil_type = GSL_DERIV_STENCIL_FORWARD;

  F.function = &f;
  F.params = NULL;

  fprintf (stderr, "f(x) = x^(3/2)\n");

  gsl_deriv_adaptive(&params, &F, x, &result, &abserr, &fs, &neval, &diag);
  fprintf (stderr, "x = %g, h_final = %.12e, neval = %u\n", x, fs, neval);
  fprintf (stderr, "f'(x) = %.12e +/- %.12e\n", result, abserr);
  fprintf (stderr, "exact = %.12e\n\n", df(x));

  /* print diagnostic information */
  gsl_deriv_diag_print(stdout, &diag);

  return 0;
}
