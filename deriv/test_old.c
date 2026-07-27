/* deriv/test_old.c
 * 
 * Copyright (C) 2000 David Morrison
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
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_test.h>

static double f1  (double x, void *p) { U(p); return exp (x); }
static double df1 (double x, void *p) { U(p); return exp (x); }

static double f2  (double x, void *p) { U(p); return (x >= 0.0) ? x*sqrt(x) : 0.0; }
static double df2 (double x, void *p) { U(p); return (x >= 0.0) ? 1.5*sqrt(x) : 0.0; }

static double f3  (double x, void *p) { U(p); return (x != 0.0) ? sin(1.0/x) : 0.0; }
static double df3 (double x, void *p) { U(p); return (x != 0.0) ? -cos(1.0/x)/(x*x) : 0.0; }

static double f4  (double x, void *p) { U(p); return exp (-x * x); }
static double df4 (double x, void *p) { U(p); return -2.0 * x * exp (-x * x); }

static double f5  (double x, void *p) { U(p); return x * x; }
static double df5 (double x, void *p) { U(p); return 2.0 * x; }

static double f6  (double x, void *p) { U(p); return 1.0 / x; }
static double df6 (double x, void *p) { U(p); return -1.0 / (x * x); }

typedef int (deriv_fn) (const gsl_function * f, double x, double h, double * res, double *abserr);

void
test (deriv_fn * deriv, gsl_function * f, gsl_function * df, double x, 
      const char * desc)
{
  double result, abserr;
  double expected = GSL_FN_EVAL (df, x);
  (*deriv) (f, x, 1e-4, &result, &abserr);

  gsl_test_abs (result, expected, GSL_MIN(1e-4,fabs(expected)) + GSL_DBL_EPSILON, desc);

  if (abserr < fabs(result-expected)) 
    {
      gsl_test_factor (abserr, fabs(result-expected), 2, "%s error estimate", desc);
    }
  else if (result == expected || expected == 0.0)
    {
      gsl_test_abs (abserr, 0.0, 1e-6, "%s abserr", desc);
    }
  else
    {
      double d = fabs(result - expected);
      gsl_test_abs (abserr, fabs(result-expected), 1e6*d, "%s abserr", desc);
    }
}

static int
test_old()
{
  gsl_function F1, DF1, F2, DF2, F3, DF3, F4, DF4, F5, DF5, F6, DF6;

  F1.function = &f1; DF1.function = &df1;
  F2.function = &f2; DF2.function = &df2;
  F3.function = &f3; DF3.function = &df3;
  F4.function = &f4; DF4.function = &df4;
  F5.function = &f5; DF5.function = &df5;
  F6.function = &f6; DF6.function = &df6;
  
  test (&gsl_deriv_central, &F1, &DF1, 1.0, "exp(x), x=1, central deriv");
  test (&gsl_deriv_forward, &F1, &DF1, 1.0, "exp(x), x=1, forward deriv");
  test (&gsl_deriv_backward, &F1, &DF1, 1.0, "exp(x), x=1, backward deriv");

  test (&gsl_deriv_central, &F2, &DF2, 0.1, "x^(3/2), x=0.1, central deriv");
  test (&gsl_deriv_forward, &F2, &DF2, 0.1, "x^(3/2), x=0.1, forward deriv");
  test (&gsl_deriv_backward, &F2, &DF2, 0.1, "x^(3/2), x=0.1, backward deriv");

  test (&gsl_deriv_central, &F3, &DF3, 0.45, "sin(1/x), x=0.45, central deriv");
  test (&gsl_deriv_forward, &F3, &DF3, 0.45, "sin(1/x), x=0.45, forward deriv");
  test (&gsl_deriv_backward, &F3, &DF3, 0.45, "sin(1/x), x=0.45, backward deriv");

  test (&gsl_deriv_central, &F4, &DF4, 0.5, "exp(-x^2), x=0.5, central deriv");
  test (&gsl_deriv_forward, &F4, &DF4, 0.5, "exp(-x^2), x=0.5, forward deriv");
  test (&gsl_deriv_backward, &F4, &DF4, 0.5, "exp(-x^2), x=0.5, backward deriv");

  test (&gsl_deriv_central, &F5, &DF5, 0.0, "x^2, x=0, central deriv");
  test (&gsl_deriv_forward, &F5, &DF5, 0.0, "x^2, x=0, forward deriv");
  test (&gsl_deriv_backward, &F5, &DF5, 0.0, "x^2, x=0, backward deriv");

  test (&gsl_deriv_central, &F6, &DF6, 10.0, "1/x, x=10, central deriv");
  test (&gsl_deriv_forward, &F6, &DF6, 10.0, "1/x, x=10, forward deriv");
  test (&gsl_deriv_backward, &F6, &DF6, 10.0, "1/x, x=10, backward deriv");

  return 0;
}
