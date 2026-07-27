/* deriv/test_d.c
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

/* directional derivative test functions */

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

typedef struct
{
  double (*f)   (const gsl_vector * x, void * params);
  double (*df)  (const gsl_vector * x, const gsl_vector * v, void * params);
  double (*d2f) (const gsl_vector * x, const gsl_vector * v, void * params);
  double (*d3f) (const gsl_vector * x, const gsl_vector * v, void * params);
  size_t n;
  double xmin;
  double xmax;
  const char * desc;
} directional_test_function;

#define S(x,i,val)  double val = gsl_vector_get(x,i)

double d_qform_f   (const gsl_vector * x, void * p) { U(p); double xnorm = gsl_blas_dnrm2(x); return xnorm*xnorm; }
double d_qform_df  (const gsl_vector * x, const gsl_vector * v, void * p) { U(p); double val; gsl_blas_ddot(x, v, &val); return 2.0*val; }
double d_qform_d2f (const gsl_vector * x, const gsl_vector * v, void * p) { U(p); U(x); double vnorm = gsl_blas_dnrm2(v); return 2.0*vnorm*vnorm; }
double d_qform_d3f (const gsl_vector * x, const gsl_vector * v, void * p) { U(p); U(x); U(v); return 0.0; }

double d_rot_f     (const gsl_vector * x, void * p) { U(p); S(x,0,x1); S(x,1,x2); return x1*x1 + x2*x2; }
double d_rot_df    (const gsl_vector * x, const gsl_vector * v, void * p) { U(p); S(x,0,x1); S(x,1,x2); S(v,0,v1); S(v,1,v2); return 2.0*(x1*v1+x2*v2); }
double d_rot_d2f   (const gsl_vector * x, const gsl_vector * v, void * p) { U(p); U(x); S(v,0,v1); S(v,1,v2); return 2.0*(v1*v1+v2*v2); }
double d_rot_d3f   (const gsl_vector * x, const gsl_vector * v, void * p) { U(p); U(x); U(v); return 0.0; }

double d_sep_f     (const gsl_vector * x, void * p) { U(p); S(x,0,x1); S(x,1,x2); S(x,2,x3); return exp(x1) + sin(x2) + x3*x3*x3*x3; }
double d_sep_df    (const gsl_vector * x, const gsl_vector * v, void * p) { U(p); S(x,0,x1); S(x,1,x2); S(x,2,x3); S(v,0,v1); S(v,1,v2); S(v,2,v3); return v1*exp(x1) + v2*cos(x2) + 4.0*v3*x3*x3*x3; }
double d_sep_d2f   (const gsl_vector * x, const gsl_vector * v, void * p) { U(p); S(x,0,x1); S(x,1,x2); S(x,2,x3); S(v,0,v1); S(v,1,v2); S(v,2,v3); return v1*v1*exp(x1) - v2*v2*sin(x2) + 12.0*v3*v3*x3*x3; }
double d_sep_d3f   (const gsl_vector * x, const gsl_vector * v, void * p) { U(p); S(x,0,x1); S(x,1,x2); S(x,2,x3); S(v,0,v1); S(v,1,v2); S(v,2,v3); return v1*v1*v1*exp(x1) - v2*v2*v2*cos(x2) + 24.0*v3*v3*v3*x3; }

static const directional_test_function test_directional_functions[] =
{
  { d_qform_f, d_qform_df, d_qform_d2f, d_qform_d3f, 0, -10.0, 10.0, "quadratic form" },
  { d_rot_f,   d_rot_df,   d_rot_d2f,   d_rot_d3f,   2, -5.0,  5.0,  "rotation" },
  { d_sep_f,   d_sep_df,   d_sep_d2f,   d_sep_d3f,   3, -5.0,  5.0,  "separable" },

  { NULL, NULL, NULL, NULL, 0, 0.0, 0.0, NULL }
};
