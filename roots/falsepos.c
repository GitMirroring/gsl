/* roots/falsepos.c
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

/* falsepos.c -- falsepos root finding algorithm 

   The false position algorithm uses bracketing by linear interpolation.

   Plain false position converges from one side only for a convex or
   concave function: one endpoint is retained on every iteration and
   the bracket does not shrink to zero (bug #40092). We use the
   Anderson-Bjorck modification: when the new point replaces the same
   endpoint as on the previous iteration, the function value at the
   retained endpoint is scaled by

     m = 1 - f(x_new) / f(x_prev)

   where x_prev is the endpoint being replaced, or by m = 1/2 (the
   Illinois method) if m <= 0. This moves the next interpolation
   point towards the retained endpoint, so that both endpoints
   converge to the root superlinearly, with one function evaluation
   per iteration. The scaled values are used only for interpolation.

   On its own this has no worst-case bound: for a function which is
   nearly flat on one side of the root, the retained value must be
   halved many times before the other endpoint moves, and the scaling
   starts again whenever it does. As a safeguard the bracket is
   required to shrink at least half as fast as with bisection: if after
   k iterations it is wider than 2 (x_upper - x_lower) 2^(-k/2), measured
   from the initial bracket, the next step is a bisection. This bounds
   the number of function evaluations by about twice that of bisection,
   while for smooth functions with a simple root the safeguard is rarely
   used. (The factor 2 allows the first, one-sided, steps to approach
   the root before the safeguard applies.)

   As in Brent's method, a new point is placed at least a distance
   delta = 2 GSL_DBL_EPSILON max(|x_lower|,|x_upper|) inside the
   bracket. When an endpoint has already converged to the root the
   interpolated point falls on (or next to) it, and a step of delta
   then collapses the bracket instead of reducing it by bisection.
   A bisection step is taken if the bracket is narrower than 2 delta
   or the interpolated point is not a number.

   The root estimate is the endpoint of the bracket with the smaller
   absolute function value.

   References:

   N. Anderson and A. Bjorck, "A new high order method of regula falsi
   type for computing a root of an equation", BIT 13, 253-264 (1973).

   M. Dowell and P. Jarratt, "A modified regula falsi method for
   computing the root of an equation", BIT 11, 168-174 (1971).
*/

#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "roots.h"

typedef struct
  {
    double f_lower, f_upper;    /* function values at the endpoints */
    double g_lower, g_upper;    /* scaled values used for interpolation */
    int side;                   /* endpoint replaced on the previous
                                   iteration: -1 lower, +1 upper, 0 none */
    double w_max;               /* largest bracket width allowed after
                                   this iteration */
    int bisect;                 /* take a bisection step next */
  }
falsepos_state_t;

static int falsepos_init (void * vstate, gsl_function * f, double * root, double x_lower, double x_upper);
static int falsepos_iterate (void * vstate, gsl_function * f, double * root, double * x_lower, double * x_upper);

static int
falsepos_init (void * vstate, gsl_function * f, double * root, double x_lower, double x_upper)
{
  falsepos_state_t * state = (falsepos_state_t *) vstate;

  double f_lower, f_upper ;

  *root = 0.5 * (x_lower + x_upper);

  SAFE_FUNC_CALL (f, x_lower, &f_lower);
  SAFE_FUNC_CALL (f, x_upper, &f_upper);
  
  state->f_lower = f_lower;
  state->f_upper = f_upper;
  state->g_lower = f_lower;
  state->g_upper = f_upper;
  state->side = 0;
  state->w_max = 2.0 * (x_upper - x_lower);
  state->bisect = 0;

  if ((f_lower < 0.0 && f_upper < 0.0) || (f_lower > 0.0 && f_upper > 0.0))
    {
      GSL_ERROR ("endpoints do not straddle y=0", GSL_EINVAL);
    }

  return GSL_SUCCESS;

}

static int
falsepos_iterate (void * vstate, gsl_function * f, double * root, double * x_lower, double * x_upper)
{
  falsepos_state_t * state = (falsepos_state_t *) vstate;

  double x_linear, f_linear;

  double x_left = *x_lower ;
  double x_right = *x_upper ;

  double f_lower = state->f_lower; 
  double f_upper = state->f_upper;

  double g_lower = state->g_lower; 
  double g_upper = state->g_upper;

  double m;

  if (f_lower == 0.0)
    {
      *root = x_left ;
      *x_upper = x_left;
      return GSL_SUCCESS;
    }
  
  if (f_upper == 0.0)
    {
      *root = x_right ;
      *x_lower = x_right;
      return GSL_SUCCESS;
    }
      
  /* Draw a line between the (scaled) values at the endpoints and note
     where it crosses the X axis; that's where we will split the
     interval. Since g_lower and g_upper have opposite signs the
     ratio lies in [0,1]. */
  
  if (state->bisect)
    {
      x_linear = 0.5 * (x_left + x_right);
    }
  else
    {
      const double delta = 2.0 * GSL_DBL_EPSILON * GSL_MAX (fabs (x_left), fabs (x_right));

      x_linear = x_right + (g_upper / (g_upper - g_lower)) * (x_left - x_right);

      if (x_right - x_left <= 2.0 * delta || gsl_isnan (x_linear))
        {
          x_linear = 0.5 * (x_left + x_right);
        }
      else if (x_linear < x_left + delta)
        {
          x_linear = x_left + delta;
        }
      else if (x_linear > x_right - delta)
        {
          x_linear = x_right - delta;
        }
    }

  SAFE_FUNC_CALL (f, x_linear, &f_linear);

  if (f_linear == 0.0)
    {
      *root = x_linear;
      *x_lower = x_linear;
      *x_upper = x_linear;
      return GSL_SUCCESS;
    }
      
  /* Discard the part of the interval which doesn't contain the root.
     If the same endpoint is replaced twice in a row, scale the value
     at the endpoint which is retained. */
  
  if ((f_lower > 0.0 && f_linear < 0.0) || (f_lower < 0.0 && f_linear > 0.0))
    {
      if (state->side == +1)
        {
          m = 1.0 - f_linear / f_upper;
          state->g_lower = g_lower * (m > 0.0 ? m : 0.5);
        }

      *x_upper = x_linear;
      state->f_upper = f_linear;
      state->g_upper = f_linear;
      state->side = +1;
    }
  else
    {
      if (state->side == -1)
        {
          m = 1.0 - f_linear / f_lower;
          state->g_upper = g_upper * (m > 0.0 ? m : 0.5);
        }

      *x_lower = x_linear;
      state->f_lower = f_linear;
      state->g_lower = f_linear;
      state->side = -1;
    }

  /* Bisect next if the bracket is shrinking more slowly than half the
     rate of bisection */

  state->w_max *= M_SQRT1_2;
  state->bisect = (*x_upper - *x_lower > state->w_max);

  /* The best estimate of the root is the endpoint with the smaller
     function value */

  *root = (fabs (state->f_lower) < fabs (state->f_upper)) ? *x_lower : *x_upper;

  return GSL_SUCCESS;
}


static const gsl_root_fsolver_type falsepos_type =
{"falsepos",                            /* name */
 sizeof (falsepos_state_t),
 &falsepos_init,
 &falsepos_iterate};

const gsl_root_fsolver_type  * gsl_root_fsolver_falsepos = &falsepos_type;
