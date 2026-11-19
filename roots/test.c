/* roots/test.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gough
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
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_deriv.h>

#define U(x) ((void) (x))              /* avoid unused parameter warnings */

/* stopping parameters */
static const double EPSREL = (10 * GSL_DBL_EPSILON);
static const double EPSABS = (10 * GSL_DBL_EPSILON);
static const unsigned int MAX_ITERATIONS = 150;

#define WITHIN_TOL(a, b, epsrel, epsabs) \
 ((fabs((a) - (b)) < (epsrel) * GSL_MIN(fabs(a), fabs(b)) + (epsabs)))

#include "test_funcs.c"
#include "test_aps.c"

/* the solvers to test */

typedef struct
{
  const gsl_root_fsolver_type * T;
  int bounded;          /* at most about twice the evaluations of bisection */
} fsolver_info;

typedef struct
{
  const gsl_root_fdfsolver_type * T;
  int multiple;         /* converges fast enough at a multiple root */
  int nonsmooth;        /* converges at a root where f is not differentiable */
} fdfsolver_info;

static gsl_function
make_function (const root_test_function * tf)
{
  gsl_function F;
  F.function = tf->test_t->f;
  F.params = (void *) &tf->params;
  return F;
}

/* gsl_function_fdf for a test function, with fdf calling f and df */

static double
fdf_f (double x, void * params)
{
  const root_test_function * tf = (const root_test_function *) params;
  return tf->test_t->f (x, (void *) &tf->params);
}

static double
fdf_df (double x, void * params)
{
  const root_test_function * tf = (const root_test_function *) params;
  return tf->test_t->df (x, (void *) &tf->params);
}

static void
fdf_fdf (double x, void * params, double * y, double * dy)
{
  *y = fdf_f (x, params);
  *dy = fdf_df (x, params);
}

static gsl_function_fdf
make_function_fdf (const root_test_function * tf)
{
  gsl_function_fdf FDF;
  FDF.f = &fdf_f;
  FDF.df = &fdf_df;
  FDF.fdf = &fdf_fdf;
  FDF.params = (void *) tf;
  return FDF;
}

/* A computed root r is correct if it agrees with the reference root,
   or if f(r) is exactly zero: such an r is a root of the function as
   computed in floating point */

static int
root_ok (gsl_function * F, double r, double root)
{
  return WITHIN_TOL (r, root, EPSREL, EPSABS) || GSL_FN_EVAL (F, r) == 0.0;
}

/* Solve with the fsolver T on the bracket of tf: the root estimate
   must stay within the bracket, the final bracket must contain a sign
   change, and the root must be correct */

static void
test_fsolver (const gsl_root_fsolver_type * T, const root_test_function * tf)
{
  gsl_function F = make_function (tf);
  gsl_root_fsolver * s = gsl_root_fsolver_alloc (T);
  int status;
  size_t iterations = 0;
  double r, a, b;

  gsl_root_fsolver_set (s, &F, tf->xmin, tf->xmax);

  do
    {
      iterations++ ;

      gsl_root_fsolver_iterate (s);

      r = gsl_root_fsolver_root (s);
      a = gsl_root_fsolver_x_lower (s);
      b = gsl_root_fsolver_x_upper (s);

      if (a > b)
        gsl_test (GSL_FAILURE, "%s, %s [%g, %g], interval is invalid (%g,%g)",
                  T->name, tf->desc, tf->xmin, tf->xmax, a, b);

      if (r < a || r > b)
        gsl_test (GSL_FAILURE, "%s, %s [%g, %g], r lies outside interval %g (%g,%g)",
                  T->name, tf->desc, tf->xmin, tf->xmax, r, a, b);

      status = gsl_root_test_interval (a, b, EPSABS, EPSREL);
    }
  while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

  gsl_test (status, "%s, %s [%g, %g] (%g obs vs %g expected)",
            T->name, tf->desc, tf->xmin, tf->xmax, r, tf->root);

  {
    const double fa = GSL_FN_EVAL (&F, a);
    const double fb = GSL_FN_EVAL (&F, b);

    gsl_test ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0),
              "%s, %s [%g, %g], final bracket contains a sign change",
              T->name, tf->desc, tf->xmin, tf->xmax);
  }

  if (!root_ok (&F, r, tf->root))
    {
      gsl_test (GSL_FAILURE, "%s, %s [%g, %g], incorrect precision (%g obs vs %g expected)",
                T->name, tf->desc, tf->xmin, tf->xmax, r, tf->root);
    }

  gsl_root_fsolver_free (s);
}

/* count the function evaluations needed by the solver T to reduce the
   bracket [lower_bound, upper_bound] to the test tolerances */

static size_t evals;

static double
count_evals_f (double x, void * params)
{
  gsl_function * f = (gsl_function *) params;
  evals++;
  return GSL_FN_EVAL (f, x);
}

static size_t
count_evals (const gsl_root_fsolver_type * T, gsl_function *f,
             double lower_bound, double upper_bound)
{
  int status;
  size_t iterations = 0;
  gsl_function F;
  gsl_root_fsolver * s = gsl_root_fsolver_alloc (T);

  F.function = &count_evals_f;
  F.params = f;
  evals = 0;

  gsl_root_fsolver_set (s, &F, lower_bound, upper_bound);

  do
    {
      iterations++;
      gsl_root_fsolver_iterate (s);
      status = gsl_root_test_interval (gsl_root_fsolver_x_lower (s),
                                       gsl_root_fsolver_x_upper (s),
                                       EPSABS, EPSREL);
    }
  while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

  gsl_root_fsolver_free (s);

  return (status == GSL_SUCCESS) ? evals : (size_t) -1;
}

/* the solver T must need at most about twice as many function
   evaluations as bisection, and no more than bisection for a smooth
   function with a simple root (bug #40092) */

static void
test_fsolver_evals (const gsl_root_fsolver_type * T, const root_test_function * tf)
{
  gsl_function F = make_function (tf);
  size_t n_bisect = count_evals (gsl_root_fsolver_bisection, &F,
                                 tf->xmin, tf->xmax);
  size_t n = count_evals (T, &F, tf->xmin, tf->xmax);

  gsl_test (n > 2 * n_bisect + 4,
            "%s, %s [%g, %g], %u evaluations vs %u for bisection",
            T->name, tf->desc, tf->xmin, tf->xmax,
            (unsigned int) n, (unsigned int) n_bisect);

  if (!(tf->flags & (ROOT_TEST_MULTIPLE | ROOT_TEST_NONSMOOTH)))
    {
      gsl_test (n > n_bisect,
                "%s, %s [%g, %g], %u evaluations no more than %u for bisection",
                T->name, tf->desc, tf->xmin, tf->xmax,
                (unsigned int) n, (unsigned int) n_bisect);
    }
}

/* the solver T must be invariant under scaling of the function:
   scaling by a power of two is exact, so the number of evaluations
   must not change */

static double scale_factor;

static double
scaled_f (double x, void * params)
{
  gsl_function * f = (gsl_function *) params;
  return scale_factor * GSL_FN_EVAL (f, x);
}

static void
test_fsolver_scale (const gsl_root_fsolver_type * T, const root_test_function * tf)
{
  const int exponents[2] = { -70, 70 };
  gsl_function F = make_function (tf);
  size_t n = count_evals (T, &F, tf->xmin, tf->xmax);
  gsl_function G;
  int i;

  G.function = &scaled_f;
  G.params = &F;

  for (i = 0; i < 2; i++)
    {
      size_t n_scaled;

      scale_factor = ldexp (1.0, exponents[i]);
      n_scaled = count_evals (T, &G, tf->xmin, tf->xmax);

      gsl_test (n_scaled != n,
                "%s, %s [%g, %g], scaled by %g, %u evaluations vs %u unscaled",
                T->name, tf->desc, tf->xmin, tf->xmax, scale_factor,
                (unsigned int) n_scaled, (unsigned int) n);
    }
}

/* the fsolver T must reject invalid brackets */

static void
test_fsolver_invalid (const gsl_root_fsolver_type * T, double lower_bound,
                      double upper_bound)
{
  gsl_function F = make_function (&test_basic_functions[0]);   /* sin(x) */
  int status;
  size_t iterations = 0;
  double x_lower, x_upper;
  gsl_root_fsolver * s = gsl_root_fsolver_alloc (T);

  status = gsl_root_fsolver_set (s, &F, lower_bound, upper_bound);

  gsl_test (status != GSL_EINVAL, "%s (set), invalid range check [%g, %g]",
            T->name, lower_bound, upper_bound);

  if (status == GSL_EINVAL)
    {
      gsl_root_fsolver_free (s);
      return ;
    }

  do
    {
      iterations++ ;
      gsl_root_fsolver_iterate (s);
      x_lower = gsl_root_fsolver_x_lower (s);
      x_upper = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lower, x_upper, EPSABS, EPSREL);
    }
  while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

  gsl_test (!status, "%s, invalid range check [%g, %g]", T->name,
            lower_bound, upper_bound);

  gsl_root_fsolver_free (s);
}

/* Solve with the fdfsolver T from the starting point of tf: it must
   converge to the reference root, or not converge if tf is flagged
   ROOT_TEST_DIVERGES */

static void
test_fdfsolver (const gsl_root_fdfsolver_type * T, const root_test_function * tf)
{
  gsl_function_fdf FDF = make_function_fdf (tf);
  gsl_function F = make_function (tf);
  gsl_root_fdfsolver * s = gsl_root_fdfsolver_alloc (T);
  const int diverges = (tf->flags & ROOT_TEST_DIVERGES);
  int status;
  size_t iterations = 0;
  double prev, r;

  status = gsl_root_fdfsolver_set (s, &FDF, tf->x0);

  gsl_test (status, "%s (set), %s {%g}", T->name, tf->desc, tf->x0);

  do
    {
      iterations++ ;
      prev = gsl_root_fdfsolver_root (s);
      gsl_root_fdfsolver_iterate (s);
      status = gsl_root_test_delta (gsl_root_fdfsolver_root (s), prev,
                                    EPSABS, EPSREL);
    }
  while (status == GSL_CONTINUE && iterations < MAX_ITERATIONS);

  r = gsl_root_fdfsolver_root (s);

  if (diverges)
    {
      gsl_test (!status, "%s, %s {%g}, does not converge", T->name,
                tf->desc, tf->x0);
    }
  else
    {
      gsl_test (status, "%s, %s {%g} (%g obs vs %g expected)",
                T->name, tf->desc, tf->x0, r, tf->root);

      if (!root_ok (&F, r, tf->root))
        {
          gsl_test (GSL_FAILURE, "%s, %s {%g}, incorrect precision (%g obs vs %g expected)",
                    T->name, tf->desc, tf->x0, r, tf->root);
        }
    }

  gsl_root_fdfsolver_free (s);
}

/* Check the derivative df of tf against a numerical derivative at the
   points used by the fdfsolvers: x0, the midpoint between x0 and the
   root, and the root itself unless f is not differentiable there. The
   variable is scaled by L = |x|/10, so that the stencil stays well
   within the region where f is smooth, and a single fourth order
   central difference with step h = eps^(1/5) is used. Its error is
   deterministic, and at most 8e-7 (near the poles of the second APS
   problem), while a wrong derivative typically differs by 1% or more. */

static const double DF_TOL = 1.0e-5;

typedef struct
{
  gsl_function * F;
  double x;
  double L;
} scaled_point;

static double
df_scaled_f (double t, void * params)
{
  const scaled_point * sp = (const scaled_point *) params;
  return GSL_FN_EVAL (sp->F, sp->x + t * sp->L);
}

static void
test_df_point (const root_test_function * tf, double x)
{
  gsl_function F = make_function (tf);
  gsl_deriv_parameters params = gsl_deriv_default_parameters (1, 4);
  scaled_point sp;
  gsl_function G;
  double result, esterr;
  unsigned int neval;

  sp.F = &F;
  sp.x = x;
  sp.L = (x != 0.0) ? 0.1 * fabs (x) : 0.1;

  G.function = &df_scaled_f;
  G.params = &sp;

  gsl_deriv_fixed (&params, &G, 0.0, pow (GSL_DBL_EPSILON, 0.2), &result,
                   &esterr, &neval);

  gsl_test_rel (result / sp.L, tf->test_t->df (x, (void *) &tf->params), DF_TOL,
                "%s, df at x = %g", tf->desc, x);
}

static void
test_df (const root_test_function * tf)
{
  if (!(tf->flags & ROOT_TEST_NO_X0))
    {
      test_df_point (tf, tf->x0);

      if (!(tf->flags & ROOT_TEST_DIVERGES))
        test_df_point (tf, 0.5 * (tf->x0 + tf->root));
    }

  if (!(tf->flags & (ROOT_TEST_NONSMOOTH | ROOT_TEST_DIVERGES)))
    test_df_point (tf, tf->root);
}

/* count calls to f, df and fdf of f(x) = x^2 - 2 */

typedef struct { int nf, ndf, nfdf; } call_counts;

static double
count_f (double x, void * params)
{
  ((call_counts *) params)->nf++;
  return x * x - 2.0;
}

static double
count_df (double x, void * params)
{
  ((call_counts *) params)->ndf++;
  return 2.0 * x;
}

static void
count_fdf (double x, void * params, double * y, double * dy)
{
  ((call_counts *) params)->nfdf++;
  *y = x * x - 2.0;
  *dy = 2.0 * x;
}

/* setting a solver needs f and f' at the starting point and should get
   both from a single fdf call (bugs #50734, #49465) */

static void
test_fdf_set_calls (const gsl_root_fdfsolver_type * T)
{
  call_counts c = { 0, 0, 0 };
  gsl_function_fdf FDF;
  gsl_root_fdfsolver * s = gsl_root_fdfsolver_alloc (T);

  FDF.f = &count_f;
  FDF.df = &count_df;
  FDF.fdf = &count_fdf;
  FDF.params = &c;

  gsl_root_fdfsolver_set (s, &FDF, 1.0);

  gsl_test (c.nf != 0 || c.ndf != 0 || c.nfdf != 1,
            "%s, set uses one fdf call (f %d, df %d, fdf %d calls)",
            gsl_root_fdfsolver_name (s), c.nf, c.ndf, c.nfdf);

  gsl_root_fdfsolver_free (s);
}

/* run all fsolvers and fdfsolvers on the rows of a table */

static void
test_table (const fsolver_info fsolvers[], const fdfsolver_info fdfsolvers[],
            const root_test_function * table)
{
  const root_test_function * tf;
  size_t i;

  for (i = 0; fsolvers[i].T != NULL; i++)
    {
      for (tf = table; tf->test_t != NULL; ++tf)
        {
          if (tf->flags & ROOT_TEST_NO_BRACKET)
            continue;

          test_fsolver (fsolvers[i].T, tf);

          if (fsolvers[i].bounded)
            {
              test_fsolver_evals (fsolvers[i].T, tf);
              test_fsolver_scale (fsolvers[i].T, tf);
            }
        }
    }

  for (tf = table; tf->test_t != NULL; ++tf)
    {
      if (tf->test_t->df != NULL)
        test_df (tf);
    }

  for (i = 0; fdfsolvers[i].T != NULL; i++)
    {
      for (tf = table; tf->test_t != NULL; ++tf)
        {
          if (tf->flags & ROOT_TEST_NO_X0)
            continue;

          /* newton and secant converge only linearly at a multiple
             root and exceed the iteration limit */
          if ((tf->flags & ROOT_TEST_MULTIPLE) && !fdfsolvers[i].multiple)
            continue;

          /* newton oscillates, x -> -x, on sqrt(|x|) sgn(x) */
          if ((tf->flags & ROOT_TEST_NONSMOOTH) && !fdfsolvers[i].nonsmooth)
            continue;

          test_fdfsolver (fdfsolvers[i].T, tf);
        }
    }
}

static void
my_error_handler (const char *reason, const char *file, int line, int err)
{
  if (0)
    printf ("(caught [%s:%d: %s (%d)])\n", file, line, reason, err);
}

int
main (void)
{
  fsolver_info fsolvers[5];
  fdfsolver_info fdfsolvers[4];
  size_t i;

  gsl_ieee_env_setup();

  fsolvers[0].T = gsl_root_fsolver_bisection; fsolvers[0].bounded = 0;
  fsolvers[1].T = gsl_root_fsolver_brent;     fsolvers[1].bounded = 0;
  fsolvers[2].T = gsl_root_fsolver_falsepos;  fsolvers[2].bounded = 1;
  fsolvers[3].T = gsl_root_fsolver_aps;       fsolvers[3].bounded = 1;
  fsolvers[4].T = NULL;                       fsolvers[4].bounded = 0;

  fdfsolvers[0].T = gsl_root_fdfsolver_newton;
  fdfsolvers[0].multiple = 0; fdfsolvers[0].nonsmooth = 0;
  fdfsolvers[1].T = gsl_root_fdfsolver_secant;
  fdfsolvers[1].multiple = 0; fdfsolvers[1].nonsmooth = 1;
  fdfsolvers[2].T = gsl_root_fdfsolver_steffenson;
  fdfsolvers[2].multiple = 1; fdfsolvers[2].nonsmooth = 1;
  fdfsolvers[3].T = NULL;
  fdfsolvers[3].multiple = 0; fdfsolvers[3].nonsmooth = 0;

  gsl_set_error_handler (&my_error_handler);

  test_table (fsolvers, fdfsolvers, test_basic_functions);
  test_table (fsolvers, fdfsolvers, test_aps_functions);

  for (i = 0; fsolvers[i].T != NULL; i++)
    {
      test_fsolver_invalid (fsolvers[i].T, 4.0, 0.0);
      test_fsolver_invalid (fsolvers[i].T, 1.0, 1.0);
      test_fsolver_invalid (fsolvers[i].T, 0.1, 0.2);
    }

  for (i = 0; fdfsolvers[i].T != NULL; i++)
    {
      test_fdf_set_calls (fdfsolvers[i].T);
    }

  exit (gsl_test_summary ());
}
