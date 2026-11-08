#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

double
f (double x, void * params)
{
  (void)(params); /* avoid unused parameter warning */
  return pow(x, 1.5);
}

double df(double x) { return 1.5*sqrt(x); }

/* print the diagnostic information in a format suitable for gnuplot */
static void
print_diag (FILE * fp, const gsl_deriv_diag * diag)
{
  const gsl_deriv_parameters * params = &(diag->params);
  const gsl_deriv_cand_t * cand = &(diag->cand);
  /* finite difference weights for the requested derivative order */
  const double * w = diag->fornberg_weights + params->deriv_order * diag->stencil_size;
  size_t i;

  fprintf(fp, "# x             = %.12e\n", diag->x);
  fprintf(fp, "# deriv_order   = %u\n", params->deriv_order);
  fprintf(fp, "# method_order  = %u\n", params->method_order);
  fprintf(fp, "# romberg_terms = %u\n", params->romberg_terms);
  fprintf(fp, "# max_step      = %.12e\n", params->max_step);
  fprintf(fp, "# min_step      = %.12e\n", params->min_step);
  fprintf(fp, "# step_ratio    = %.12e\n", diag->step_ratio);
  fprintf(fp, "# nstep         = %lu\n", (unsigned long) params->nstep);
  fprintf(fp, "# nromb         = %lu\n", (unsigned long) diag->nromb);

  fprintf(fp, "# stencil       = [ ");
  for (i = 0; i < diag->stencil_size; ++i)
    fprintf(fp, "%g ", diag->stencil[i]);
  fprintf(fp, "]\n");

  fprintf(fp, "# weights       = [ ");
  for (i = 0; i < diag->stencil_size; ++i)
    fprintf(fp, "%g ", w[i]);
  fprintf(fp, "]\n");

  fprintf(fp, "# der_select    = %.12e\n", diag->der_select);
  fprintf(fp, "# err_select    = %.12e\n", diag->err_select);
  fprintf(fp, "# step_select   = %.12e\n", diag->step_select);

  /* block 0: initial estimates at each step size */
  fprintf(fp, "# Index 0\n");
  fprintf(fp, "# Field 1: step size\n");
  fprintf(fp, "# Field 2: initial derivative estimate\n");
  fprintf(fp, "# Field 3: initial derivative roundoff error\n");
  for (i = 0; i < diag->nstep; ++i)
    fprintf(fp, "%20.12e %20.12e %20.12e\n", diag->h[i], diag->der_init[i], diag->err_init[i]);

  /* block 1: Romberg extrapolated estimates */
  fprintf(fp, "\n\n# Index 1\n");
  fprintf(fp, "# Field 1: step size\n");
  fprintf(fp, "# Field 2: Romberg derivative estimate\n");
  fprintf(fp, "# Field 3: Romberg error estimate\n");
  for (i = 0; i < diag->nromb; ++i)
    fprintf(fp, "%20.12e %20.12e %20.12e\n", diag->h[i], diag->der_romb[i], diag->err_romb[i]);

  /* block 2: candidates remaining after filtering */
  fprintf(fp, "\n\n# Index 2\n");
  fprintf(fp, "# Field 1: step size\n");
  fprintf(fp, "# Field 2: final candidate derivative estimate\n");
  fprintf(fp, "# Field 3: final candidate error estimate\n");
  for (i = 0; i < cand->nkeep; ++i)
    fprintf(fp, "%20.12e %20.12e %20.12e\n", diag->h[cand->index[i]],
            diag->der_romb[cand->index[i]], diag->err_romb[cand->index[i]]);
}

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
  print_diag (stdout, &diag);

  return 0;
}
