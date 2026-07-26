/* deriv/test_f.c
 * 
 * Copyright (C) 2026 Patrick Alken
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
#include <gsl/gsl_errno.h>

typedef struct
{
  double (*f)  (double x, void * params);
  double (*df) (double x, void * params);
  double (*d2f)(double x, void * params);
  double (*d3f)(double x, void * params);
  double (*d4f)(double x, void * params);
} deriv_test_t;

typedef struct
{
  const deriv_test_t * test_t;
  double xmin;
  double xmax;
  double a;
  double min_step;
  double max_step;
  const char * desc;
} deriv_test_function;

#define U(x) ((void) (x))              /* avoid unused parameter warnings */
#define A(x)  double a = *(double*)p   /* avoid unused parameter warnings */

double exp_f    (double x, void * p) { U(p); return exp(x); }
double exp_df   (double x, void * p) { U(p); return exp(x); }
double exp_d2f  (double x, void * p) { U(p); return exp(x); }
double exp_d3f  (double x, void * p) { U(p); return exp(x); }
double exp_d4f  (double x, void * p) { U(p); return exp(x); }

double expa_f   (double x, void * p) { A(p); return exp(a*x); }
double expa_df  (double x, void * p) { A(p); return a*exp(a*x); }
double expa_d2f (double x, void * p) { A(p); return a*a*exp(a*x); }
double expa_d3f (double x, void * p) { A(p); return a*a*a*exp(a*x); }
double expa_d4f (double x, void * p) { A(p); return a*a*a*a*exp(a*x); }

double exp2_f   (double x, void * p) { U(p); return exp(-x*x); }
double exp2_df  (double x, void * p) { U(p); return -2.0*x*exp(-x*x); }
double exp2_d2f (double x, void * p) { U(p); return -2.0*(1.0 - 2.0*x*x)*exp(-x*x); }
double exp2_d3f (double x, void * p) { U(p); return -4.0*x*(-3.0 + 2.0*x*x)*exp(-x*x); }
double exp2_d4f (double x, void * p) { U(p); return  4.0*(3.0 + 4.0*x*x*(-3.0+x*x))*exp(-x*x); }

double sin_f    (double x, void * p) { U(p); return  sin(x); }
double sin_df   (double x, void * p) { U(p); return  cos(x); }
double sin_d2f  (double x, void * p) { U(p); return -sin(x); }
double sin_d3f  (double x, void * p) { U(p); return -cos(x); }
double sin_d4f  (double x, void * p) { U(p); return  sin(x); }

double sina_f   (double x, void * p) { A(p); return  sin(a*x); }
double sina_df  (double x, void * p) { A(p); return  a*cos(a*x); }
double sina_d2f (double x, void * p) { A(p); return -a*a*sin(a*x); }
double sina_d3f (double x, void * p) { A(p); return -a*a*a*cos(a*x); }
double sina_d4f (double x, void * p) { A(p); return  a*a*a*a*sin(a*x); }

double cos_f    (double x, void * p) { U(p); return  cos(x); }
double cos_df   (double x, void * p) { U(p); return -sin(x); }
double cos_d2f  (double x, void * p) { U(p); return -cos(x); }
double cos_d3f  (double x, void * p) { U(p); return  sin(x); }
double cos_d4f  (double x, void * p) { U(p); return  cos(x); }

double atan_f   (double x, void * p) { U(p); return  atan(x); }
double atan_df  (double x, void * p) { U(p); return  1.0/(1.0+x*x); }
double atan_d2f (double x, void * p) { U(p); return -2.0*x/gsl_pow_uint(1.0+x*x,2); }
double atan_d3f (double x, void * p) { U(p); return  (-2.0+6.0*x*x)/gsl_pow_uint(1.0+x*x,3); }
double atan_d4f (double x, void * p) { U(p); return  (24.0*x*(1.0-x*x))/gsl_pow_uint(1.0+x*x,4); }

double esin_f   (double x, void * p) { U(p); return  exp(-x)*sin(x); }
double esin_df  (double x, void * p) { U(p); return  exp(-x)*(cos(x)-sin(x)); }
double esin_d2f (double x, void * p) { U(p); return -2.0*exp(-x)*cos(x); }
double esin_d3f (double x, void * p) { U(p); return  2.0*exp(-x)*(cos(x)+sin(x)); }
double esin_d4f (double x, void * p) { U(p); return -4.0*exp(-x)*sin(x); }

double rat1_f   (double x, void * p) { U(p); return  1.0/(1.0+x); }
double rat1_df  (double x, void * p) { U(p); return -1.0/gsl_pow_uint(1.0+x,2); }
double rat1_d2f (double x, void * p) { U(p); return  2.0/gsl_pow_uint(1.0+x,3); }
double rat1_d3f (double x, void * p) { U(p); return -6.0/gsl_pow_uint(1.0+x,4); }
double rat1_d4f (double x, void * p) { U(p); return  24.0/gsl_pow_uint(1.0+x,5); }

double rat2_f   (double x, void * p) { U(p); return  1.0/(1.0+x*x); }
double rat2_df  (double x, void * p) { U(p); return -2.0*x/gsl_pow_uint(1.0+x*x,2); }
double rat2_d2f (double x, void * p) { U(p); return (-2.0+6.0*x*x)/gsl_pow_uint(1.0+x*x,3); }
double rat2_d3f (double x, void * p) { U(p); return (24.0*x*(1.0-x*x))/gsl_pow_uint(1.0+x*x,4); }
double rat2_d4f (double x, void * p) { U(p); return (24.0*(1.0+5.0*x*x*(-2.0+x*x)))/gsl_pow_uint(1.0+x*x,5); }

double rat3_f   (double x, void * p) { U(p); return  x/(1.0+x*x); }
double rat3_df  (double x, void * p) { U(p); return  (1.0+x)*(1.0-x)/gsl_pow_uint(1.0+x*x,2); }
double rat3_d2f (double x, void * p) { U(p); return  2.0*x*(-3.0+x*x)/gsl_pow_uint(1.0+x*x,3); }
double rat3_d3f (double x, void * p) { U(p); return -6.0*(1.0+x*x*(x*x-6.0))/gsl_pow_uint(1.0+x*x,4); }
double rat3_d4f (double x, void * p) { U(p); return  24.0*x*(5.0+x*x*(x*x-10.0))/gsl_pow_uint(1.0+x*x,5); }

double sinh_f   (double x, void * p) { U(p); return  sinh(x); }
double sinh_df  (double x, void * p) { U(p); return  cosh(x); }
double sinh_d2f (double x, void * p) { U(p); return  sinh(x); }
double sinh_d3f (double x, void * p) { U(p); return  cosh(x); }
double sinh_d4f (double x, void * p) { U(p); return  sinh(x); }

double cosh_f   (double x, void * p) { U(p); return  cosh(x); }
double cosh_df  (double x, void * p) { U(p); return  sinh(x); }
double cosh_d2f (double x, void * p) { U(p); return  cosh(x); }
double cosh_d3f (double x, void * p) { U(p); return  sinh(x); }
double cosh_d4f (double x, void * p) { U(p); return  cosh(x); }

double tanh_f   (double x, void * p) { U(p); return  tanh(x); }
double tanh_df  (double x, void * p) { U(p); return  gsl_pow_uint(1.0/cosh(x), 2); }
double tanh_d2f (double x, void * p) { U(p); return -2.0*tanh(x)*gsl_pow_uint(1.0/cosh(x), 2); }
double tanh_d3f (double x, void * p) { U(p); return  2.0*(-2.0+cosh(2.0*x))*gsl_pow_uint(1.0/cosh(x), 4); }
double tanh_d4f (double x, void * p) { U(p); return -4.0*(-5.0+cosh(2.0*x))*tanh(x)*gsl_pow_uint(1.0/cosh(x), 4); }

double p1_f     (double x, void * p) { U(p);       return  x; }
double p1_df    (double x, void * p) { U(x); U(p); return  1.0; }
double p1_d2f   (double x, void * p) { U(x); U(p); return  0.0; }
double p1_d3f   (double x, void * p) { U(x); U(p); return  0.0; }
double p1_d4f   (double x, void * p) { U(x); U(p); return  0.0; }

double p2_f     (double x, void * p) { U(p);       return  x*x; }
double p2_df    (double x, void * p) { U(p);       return  2.0*x; }
double p2_d2f   (double x, void * p) { U(x); U(p); return  2.0; }
double p2_d3f   (double x, void * p) { U(x); U(p); return  0.0; }
double p2_d4f   (double x, void * p) { U(x); U(p); return  0.0; }

double p3_f     (double x, void * p) { U(p);       return  x*x*x; }
double p3_df    (double x, void * p) { U(p);       return  3.0*x*x; }
double p3_d2f   (double x, void * p) { U(p);       return  6.0*x; }
double p3_d3f   (double x, void * p) { U(x); U(p); return  6.0; }
double p3_d4f   (double x, void * p) { U(x); U(p); return  0.0; }

double p4_f     (double x, void * p) { U(p);       return  x*x*x*x; }
double p4_df    (double x, void * p) { U(p);       return  4.0*x*x*x; }
double p4_d2f   (double x, void * p) { U(p);       return  12.0*x*x; }
double p4_d3f   (double x, void * p) { U(p);       return  24.0*x; }
double p4_d4f   (double x, void * p) { U(x); U(p); return  24.0; }

const deriv_test_t exp_test  = {  exp_f,  exp_df,  exp_d2f,  exp_d3f,  exp_d4f };
const deriv_test_t exp2_test = { exp2_f, exp2_df, exp2_d2f, exp2_d3f, exp2_d4f };
const deriv_test_t sin_test  = {  sin_f,  sin_df,  sin_d2f,  sin_d3f,  sin_d4f };
const deriv_test_t cos_test  = {  cos_f,  cos_df,  cos_d2f,  cos_d3f,  cos_d4f };
const deriv_test_t atan_test = { atan_f, atan_df, atan_d2f, atan_d3f, atan_d4f };
const deriv_test_t esin_test = { esin_f, esin_df, esin_d2f, esin_d3f, esin_d4f };
const deriv_test_t expa_test = { expa_f, expa_df, expa_d2f, expa_d3f, expa_d4f };
const deriv_test_t sina_test = { sina_f, sina_df, sina_d2f, sina_d3f, sina_d4f };
const deriv_test_t rat1_test = { rat1_f, rat1_df, rat1_d2f, rat1_d3f, rat1_d4f };
const deriv_test_t rat2_test = { rat2_f, rat2_df, rat2_d2f, rat2_d3f, rat2_d4f };
const deriv_test_t rat3_test = { rat3_f, rat3_df, rat3_d2f, rat3_d3f, rat3_d4f };
const deriv_test_t sinh_test = { sinh_f, sinh_df, sinh_d2f, sinh_d3f, sinh_d4f };
const deriv_test_t cosh_test = { cosh_f, cosh_df, cosh_d2f, cosh_d3f, cosh_d4f };
const deriv_test_t tanh_test = { tanh_f, tanh_df, tanh_d2f, tanh_d3f, tanh_d4f };
const deriv_test_t p1_test   = {   p1_f,   p1_df,   p1_d2f,   p1_d3f,   p1_d4f };
const deriv_test_t p2_test   = {   p2_f,   p2_df,   p2_d2f,   p2_d3f,   p2_d4f };
const deriv_test_t p3_test   = {   p3_f,   p3_df,   p3_d2f,   p3_d3f,   p3_d4f };
const deriv_test_t p4_test   = {   p4_f,   p4_df,   p4_d2f,   p4_d3f,   p4_d4f };

#if 1
static const deriv_test_function test_normal_functions[] =
{
  {  &exp_test,  -5.0,   5.0,   0.0,   -1.0,    -1.0,    "exp(x)" },
  { &exp2_test,  -3.0,   3.0,   0.0,   -1.0,    -1.0,    "exp(-x*x)" },
  {  &sin_test,  -M_PI,  M_PI,  0.0,   -1.0,    -1.0,    "sin(x)" },
  {  &cos_test,  -M_PI,  M_PI,  0.0,   -1.0,    -1.0,    "cos(x)" },
  { &atan_test,  -10.0,  10.0,  0.0,   -1.0,    -1.0,    "atan(x)" },
  { &esin_test,  -3.0,   5.0,   0.0,   -1.0,    -1.0,    "exp(-x)*sin(x)" },

  { &sinh_test,  -5.0,   5.0,   0.0,   -1.0,    -1.0,    "sinh(x)" },
  { &cosh_test,  -5.0,   5.0,   0.0,   -1.0,    -1.0,    "cosh(x)" },
  { &tanh_test,  -4.0,   4.0,   0.0,   -1.0,    -1.0,    "tanh(x)" },

  { &rat1_test,   0.0,   5.0,   0.0,   -1.0,    -1.0,    "1/(1+x)" },
  { &rat2_test,  -5.0,   5.0,   0.0,   -1.0,    -1.0,    "1/(1+x^2)" },
  { &rat3_test,  -10.0,  10.0,  0.0,   -1.0,    -1.0,    "x/(1+x^2)" },

  { NULL, 0.0, 0.0, 0.0, -1.0, -1.0, NULL }
};

static const deriv_test_function test_poly_functions[] =
{
  {  &p1_test,   -5.0,   5.0,   0.0,   -1.0,    -1.0,    "x" },
  {  &p2_test,   -5.0,   5.0,   0.0,   -1.0,    -1.0,    "x^2" },
  {  &p3_test,   -5.0,   5.0,   0.0,   -1.0,    -1.0,    "x^3" },
  {  &p4_test,   -5.0,   5.0,   0.0,   -1.0,    -1.0,    "x^4" },
  { NULL, 0.0, 0.0, 0.0, -1.0, -1.0, NULL }
};
#else
/*
 * Library of test functions
 * The min_step and max_step were hand-tuned by examining
 * the diagnostic information output by gsl_deriv_print_diag()
 */
static const deriv_test_function test_functions[] =
{
  {  &exp_test,  -5.0,   5.0,   0.0,    1.0e-6,  1.0e0,  "exp(x)" },
  { &exp2_test,  -5.0,   5.0,   0.0,    1.0e-6,  1.0e0,  "exp(-x*x)" },
  {  &sin_test,  -M_PI,  M_PI,  0.0,   -1.0,    -1.0,    "sin(x)" },
  {  &cos_test,  -M_PI,  M_PI,  0.0,   -1.0,    -1.0,    "cos(x)" },
  { &esin_test,  -3.0,   5.0,   0.0,    1.0e-4,  1.0e-1, "exp(-x)*sin(x)" },

  { &expa_test,  -5.0,   5.0,   0.01,   1.0e-6,  1.0e2,  "exp(a*x)" },
  { &expa_test,  -5.0,   5.0,   0.1,    1.0e-5,  1.0e2,  "exp(a*x)" },
  { &expa_test,  -5.0,   5.0,   1.0,   -1.0,    -1.0,    "exp(a*x)" },
  { &expa_test,  -3.0,   3.0,   10.0,   1.0e-8,  1.0e0,  "exp(a*x)" },

  { &sina_test,  -M_PI,  M_PI,  0.01,   1.0e-8,  1.0e3,  "sin(a*x)" },
  { &sina_test,  -M_PI,  M_PI,  0.1,    1.0e-8,  1.0e3,  "sin(a*x)" },
  { &sina_test,  -M_PI,  M_PI,  1.0,    1.0e-5,  1.0e1,  "sin(a*x)" },
  { &sina_test,  -M_PI,  M_PI,  10.0,   1.0e-6,  1.0e-1, "sin(a*x)" },
  { &sina_test,  -M_PI,  M_PI,  100.0,  1.0e-7,  1.0e-1, "sin(a*x)" },

  {  &p1_test,   -5.0,   5.0,   0.0,   -1.0,    -1.0,    "x" },
  {  &p2_test,   -5.0,   5.0,   0.0,   -1.0,    -1.0,    "x^2" },
  {  &p3_test,   -5.0,   5.0,   0.0,   -1.0,    -1.0,    "x^3" },
  {  &p4_test,   -5.0,   5.0,   0.0,   -1.0,    -1.0,    "x^4" },

  { NULL, 0.0, 0.0, 0.0, -1.0, -1.0, NULL }
};
#endif

#if 0
double f2 (double x, void *params) { return (x >= 0.0) ? x * sqrt(x) : 0.0; }
double df2 (double x, void *params) { return (x >= 0.0) ? 1.5 * sqrt(x) : 0.0; }

double f3 (double x, void *params) { return (x != 0.0) ? sin(1.0/x) : 0.0; }
double df3 (double x, void *params) { return (x != 0.0) ? -cos(1.0/x)/(x*x) : 0.0; }

double f6 (double x, void *params) { return 1.0 / x; }
double df6 (double x, void *params) { return -1.0 / (x * x); }
#endif
