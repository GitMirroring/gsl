/* roots/test_funcs.c
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

/* Test functions for the root finding solvers. Each row of a table is
   one test problem: a function and its derivative, the parameters
   passed to them, a bracket for the fsolvers, a starting point for the
   fdfsolvers, the root and flags describing its properties. The flags
   NO_BRACKET and NO_X0 mark rows which are not used for one type of
   solver. This file is included by test.c. */

typedef struct
{
  double (*f)  (double x, void * params);
  double (*df) (double x, void * params);
} root_test_t;

/* parameters of a family of test functions, passed to f and df */
typedef struct
{
  double n;
  double a;
  double b;
} root_test_params;

/* flags */
#define ROOT_TEST_MULTIPLE  1   /* root of multiplicity greater than one */
#define ROOT_TEST_NONSMOOTH 2   /* f is not differentiable at the root */
#define ROOT_TEST_DIVERGES  4   /* the fdfsolvers must not converge from x0 */
#define ROOT_TEST_NO_BRACKET 8  /* no bracket: not used for the fsolvers */
#define ROOT_TEST_NO_X0     16  /* no starting point: not used for the fdfsolvers */

typedef struct
{
  const root_test_t * test_t;
  root_test_params params;
  double xmin, xmax;            /* bracket for the fsolvers */
  double x0;                    /* starting point for the fdfsolvers */
  double root;
  int flags;
  const char * desc;
} root_test_function;

static double sin_f  (double x, void * p) { U(p); return sin (x); }
static double sin_df (double x, void * p) { U(p); return cos (x); }

static double cos_f  (double x, void * p) { U(p); return cos (x); }
static double cos_df (double x, void * p) { U(p); return -sin (x); }

/* x^20 - 1 */
static double pow20_f  (double x, void * p) { U(p); return pow (x, 20.0) - 1.0; }
static double pow20_df (double x, void * p) { U(p); return 20.0 * pow (x, 19.0); }

/* sqrt(|x|) sgn(x), not differentiable at the root 0: a newton step
   maps x to -x */
static double sqrtsgn_f (double x, void * p)
{
  U(p);
  return (x > 0.0) ? sqrt (x) : (x < 0.0) ? -sqrt (-x) : 0.0;
}
static double sqrtsgn_df (double x, void * p) { U(p); return 0.5 / sqrt (fabs (x)); }

/* x^2 - 1e-8 */
static double sq_f  (double x, void * p) { U(p); return x * x - 1.0e-8; }
static double sq_df (double x, void * p) { U(p); return 2.0 * x; }

/* x exp(-x), which tends to 0 as x -> +inf */
static double xexp_f  (double x, void * p) { U(p); return x * exp (-x); }
static double xexp_df (double x, void * p) { U(p); return (1.0 - x) * exp (-x); }

/* 1 / (1 + exp(x)), which has no root */
static double logistic_f  (double x, void * p) { U(p); return 1.0 / (1.0 + exp (x)); }
static double logistic_df (double x, void * p) { U(p); return -exp (x) / gsl_pow_2 (1.0 + exp (x)); }

/* (x - 1)^7 */
static double pow7_f  (double x, void * p) { U(p); return pow (x - 1.0, 7.0); }
static double pow7_df (double x, void * p) { U(p); return 7.0 * pow (x - 1.0, 6.0); }

/* -pi x + e, linear: the solvers must exit correctly at an exact root */
static double lin_f  (double x, void * p) { U(p); return -M_PI * x + M_E; }
static double lin_df (double x, void * p) { U(x); U(p); return -M_PI; }

static const root_test_t sin_test      = { sin_f,      sin_df };
static const root_test_t cos_test      = { cos_f,      cos_df };
static const root_test_t pow20_test    = { pow20_f,    pow20_df };
static const root_test_t sqrtsgn_test  = { sqrtsgn_f,  sqrtsgn_df };
static const root_test_t sq_test       = { sq_f,       sq_df };
static const root_test_t xexp_test     = { xexp_f,     xexp_df };
static const root_test_t logistic_test = { logistic_f, logistic_df };
static const root_test_t pow7_test     = { pow7_f,     pow7_df };
static const root_test_t lin_test      = { lin_f,      lin_df };

#define NOP  { 0.0, 0.0, 0.0 }  /* no parameters */
#define NOBR ROOT_TEST_NO_BRACKET

static const root_test_function test_basic_functions[] =
{
  /* function       params  xmin        xmax    x0     root         flags */
  { &sin_test,      NOP,     3.0,        4.0,    3.4,   M_PI,        0, "sin(x)" },
  { &sin_test,      NOP,    -4.0,       -3.0,   -3.3,  -M_PI,        0, "sin(x)" },
  { &sin_test,      NOP,    -1.0 / 3.0,  1.0,    0.5,   0.0,         0, "sin(x)" },
  { &cos_test,      NOP,     0.0,        3.0,    0.6,   M_PI_2,      0, "cos(x)" },
  { &cos_test,      NOP,    -3.0,        0.0,   -2.5,  -M_PI_2,      0, "cos(x)" },
  { &pow20_test,    NOP,     0.1,        2.0,    0.9,   1.0,         0, "x^20 - 1" },
  { &pow20_test,    NOP,     0.0,        0.0,    1.1,   1.0,         NOBR, "x^20 - 1" },
  { &sqrtsgn_test,  NOP,    -1.0 / 3.0,  1.0,    0.001, 0.0,         ROOT_TEST_NONSMOOTH, "sqrt(|x|) sgn(x)" },
  { &sq_test,       NOP,     0.0,        1.0,    1.0,   1.0e-4,      0, "x^2 - 1e-8" },
  { &xexp_test,     NOP,    -1.0 / 3.0,  2.0,   -2.0,   0.0,         0, "x exp(-x)" },
  { &xexp_test,     NOP,     0.0,        0.0,    2.0,   0.0,         NOBR | ROOT_TEST_DIVERGES, "x exp(-x), x -> +inf" },
  { &logistic_test, NOP,     0.0,        0.0,    0.0,   0.0,         NOBR | ROOT_TEST_DIVERGES, "1/(1 + exp(x)), x -> +inf" },
  { &pow7_test,     NOP,     0.9995,     1.0002, 0.9,   1.0,         ROOT_TEST_MULTIPLE, "(x - 1)^7" },
  { &lin_test,      NOP,     0.0,        0.0,    1.5,   M_E / M_PI,  NOBR, "-pi x + e" },
  { NULL,           NOP,     0.0,        0.0,    0.0,   0.0,         0, NULL }
};
