/* roots/aps.c
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

/* aps.c -- Alefeld-Potra-Shi root finding algorithm

   This is based on Algorithm 4.2 of

   G. E. Alefeld, F. A. Potra and Y. Shi, "Algorithm 748: Enclosing
   zeros of continuous functions", ACM Transactions on Mathematical
   Software 21, 327-344 (1995),

   implemented here from the description in the paper. Each iteration
   of the algorithm takes three steps, each of which evaluates the
   function once and replaces one endpoint of the bracket [a,b]:

   1. an interpolation step, using inverse cubic interpolation through
      a, b and the two most recently discarded endpoints d and e when
      their function values are distinct, and otherwise two Newton
      steps on the quadratic through a, b and d;

   2. a second interpolation step, as in 1 but with three Newton steps;

   3. a double-length secant step from the endpoint u with the smaller
      function value.

   The very first step is a secant step, and the first interpolation
   step always uses the quadratic, since e is not yet defined.

   We differ from the paper in three ways, each of which reduced the
   number of function evaluations on the paper's test problems and on
   families of simple and multiple roots:

   - Safeguard. The paper takes a bisection step whenever an iteration
     does not halve the bracket, which allows four evaluations per
     halving. Instead, as in our false position solver, the bracket is
     required to shrink at least half as fast as with bisection: if
     after k steps it is wider than 2 (x_upper - x_lower) 2^(-k/2),
     measured from the initial bracket, the next step is a bisection,
     after which the sequence of steps resumes. This bounds the number
     of evaluations by about twice that of bisection, while fast
     phases of convergence are not interrupted.

   - Distinct function values. Inverse cubic interpolation is used
     only if the four function values differ pairwise by more than
     32 GSL_DBL_EPSILON times the larger function value at the initial
     endpoints, i.e. by more than rounding noise on the scale of the
     problem. Near a multiple root the function values become
     negligible on this scale and the inverse function is not well
     approximated by a polynomial, so the quadratic is used instead.
     (SciPy uses an absolute threshold 32 GSL_DBL_EPSILON here; the
     relative threshold makes the solver invariant under scaling of
     the function.)

   - Stalled secant step. If the double-length secant step does not
     move away from u, and |f(u)| is smaller than the function value at
     the other endpoint by more than a factor 2^50, the step is taken a
     distance 1/32 of the bracket into it instead, as in SciPy.

   Each call to the iterate function performs one step, so that the
   user's convergence test is applied after every function evaluation.
   Every new point is moved, if necessary, to lie at least a few units
   of machine precision inside the bracket (see interior_point in
   roots.h); this plays the role of the tolerance used in the paper's
   bracketing routine.

   The root estimate is the endpoint of the bracket with the smaller
   absolute function value. */

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

enum
{
  STEP_SECANT,                /* initial secant step */
  STEP_INTERP1,               /* first interpolation step */
  STEP_INTERP2,               /* second interpolation step */
  STEP_DSECANT,               /* double-length secant step */
  STEP_BISECT                 /* safeguarding bisection step */
};

typedef struct
{
  double f_lower, f_upper;    /* function values at the endpoints */
  double d, f_d;              /* endpoint discarded by the last step */
  double e, f_e;              /* an earlier discarded endpoint */
  double f_tol;               /* function values closer than this are
                                 not distinct */
  double w_max;               /* largest bracket width allowed after
                                 this step */
  int step;                   /* next step */
  int resume;                 /* step to take after a bisection */
  int first;                  /* first interpolation step */
} aps_state_t;

static int aps_init (void * vstate, gsl_function * f, double * root, double x_lower, double x_upper);
static int aps_iterate (void * vstate, gsl_function * f, double * root, double * x_lower, double * x_upper);

/* Apply k Newton steps to the quadratic

     P(x) = f(a) + f[a,b] (x - a) + f[a,b,d] (x - a)(x - b)

   starting from the endpoint at which P has the same sign as its
   second derivative, so that the iterates are monotone. If the
   quadratic degenerates to a line, or a step fails, the secant point
   of a and b is returned. */

static double
newton_quadratic (const double a, const double b, const double d,
                  const double fa, const double fb, const double fd,
                  const int k)
{
  const double B = (fb - fa) / (b - a);
  const double A = ((fd - fb) / (d - b) - B) / (d - a);
  const double secant = a - fa / B;
  double r;
  int i;

  if (A == 0.0 || !gsl_finite (A))
    {
      return secant;
    }

  r = (A * fa > 0.0) ? a : b;

  for (i = 0; i < k; i++)
    {
      const double p = fa + (B + A * (r - b)) * (r - a);
      const double dp = B + A * (2.0 * r - a - b);

      if (dp == 0.0)
        {
          return secant;
        }

      r -= p / dp;
    }

  if (!(r > a && r < b))
    {
      return secant;
    }

  return r;
}

/* Inverse cubic interpolation: the value at y = 0 of the cubic
   polynomial x(y) through the points (fa,a), (fb,b), (fd,d), (fe,e),
   computed by Neville's scheme. The function values must be
   distinct. */

static double
inverse_cubic (const double a, const double b, const double d, const double e,
               const double fa, const double fb, const double fd, const double fe)
{
  const double q11 = (d - e) * fd / (fe - fd);
  const double q21 = (b - d) * fb / (fd - fb);
  const double q31 = (a - b) * fa / (fb - fa);
  const double d21 = (b - d) * fd / (fd - fb);
  const double d31 = (a - b) * fb / (fb - fa);
  const double q22 = (d21 - q11) * fb / (fe - fb);
  const double q32 = (d31 - q21) * fa / (fd - fa);
  const double d32 = (d31 - q21) * fd / (fd - fa);
  const double q33 = (d32 - q22) * fa / (fe - fa);

  return a + q31 + q32 + q33;
}

/* Interpolation step: inverse cubic interpolation if the four
   function values are distinct (see above) and the result lies in
   (a,b), otherwise k Newton steps on the quadratic through a, b, d */

static int
distinct (const double x, const double y, const double tol)
{
  return fabs (x - y) > tol;
}

static double
interpolate (const aps_state_t * state, const double a, const double b,
             const int k)
{
  const double fa = state->f_lower, fb = state->f_upper;
  const double fd = state->f_d, fe = state->f_e;
  const double tol = state->f_tol;

  if (!state->first
      && distinct (fa, fb, tol) && distinct (fa, fd, tol)
      && distinct (fa, fe, tol) && distinct (fb, fd, tol)
      && distinct (fb, fe, tol) && distinct (fd, fe, tol))
    {
      const double c = inverse_cubic (a, b, state->d, state->e, fa, fb, fd, fe);

      if (c > a && c < b)
        {
          return c;
        }
    }

  return newton_quadratic (a, b, state->d, fa, fb, fd, k);
}

static int
aps_init (void * vstate, gsl_function * f, double * root, double x_lower, double x_upper)
{
  aps_state_t * state = (aps_state_t *) vstate;

  double f_lower, f_upper ;

  *root = 0.5 * (x_lower + x_upper);

  SAFE_FUNC_CALL (f, x_lower, &f_lower);
  SAFE_FUNC_CALL (f, x_upper, &f_upper);

  state->f_lower = f_lower;
  state->f_upper = f_upper;
  state->d = x_upper;
  state->f_d = f_upper;
  state->e = x_upper;
  state->f_e = f_upper;
  state->f_tol = 32.0 * GSL_DBL_EPSILON * GSL_MAX (fabs (f_lower), fabs (f_upper));
  state->w_max = 2.0 * (x_upper - x_lower);
  state->step = STEP_SECANT;
  state->resume = STEP_INTERP1;
  state->first = 1;

  if ((f_lower < 0.0 && f_upper < 0.0) || (f_lower > 0.0 && f_upper > 0.0))
    {
      GSL_ERROR ("endpoints do not straddle y=0", GSL_EINVAL);
    }

  return GSL_SUCCESS;
}

static int
aps_iterate (void * vstate, gsl_function * f, double * root, double * x_lower, double * x_upper)
{
  aps_state_t * state = (aps_state_t *) vstate;

  const double a = *x_lower;
  const double b = *x_upper;
  const double fa = state->f_lower;
  const double fb = state->f_upper;

  double c, fc;
  int next;

  if (fa == 0.0)
    {
      *root = a;
      *x_upper = a;
      return GSL_SUCCESS;
    }

  if (fb == 0.0)
    {
      *root = b;
      *x_lower = b;
      return GSL_SUCCESS;
    }

  switch (state->step)
    {
    case STEP_SECANT:
      c = a - fa * (b - a) / (fb - fa);
      next = STEP_INTERP1;
      break;

    case STEP_INTERP1:
      c = interpolate (state, a, b, 2);
      state->first = 0;
      state->e = state->d;
      state->f_e = state->f_d;
      next = STEP_INTERP2;
      break;

    case STEP_INTERP2:
      c = interpolate (state, a, b, 3);
      next = STEP_DSECANT;
      break;

    case STEP_DSECANT:
      {
        const double u = (fabs (fa) < fabs (fb)) ? a : b;
        const double fu = (fabs (fa) < fabs (fb)) ? fa : fb;

        c = u - 2.0 * fu * (b - a) / (fb - fa);

        if (fabs (c - u) > 0.5 * (b - a))
          {
            c = a + 0.5 * (b - a);
          }
        else if (fabs (c - u) <= GSL_DBL_EPSILON * fabs (u))
          {
            /* the step has stalled at u: if |f(u)| is much smaller
               than the value at the other endpoint v, move 1/32 of
               the way towards v */

            const double v = (u == a) ? b : a;
            const double fv = (u == a) ? fb : fa;
            int exp_u, exp_v;

            frexp (fu, &exp_u);
            frexp (fv, &exp_v);

            if (exp_u < exp_v - 50)
              {
                c = (31.0 * u + v) / 32.0;
              }
          }

        state->e = state->d;
        state->f_e = state->f_d;
        next = STEP_INTERP1;
      }
      break;

    default:                    /* STEP_BISECT */
      c = a + 0.5 * (b - a);
      state->e = state->d;
      state->f_e = state->f_d;
      next = state->resume;
      break;
    }

  c = interior_point (c, a, b);

  SAFE_FUNC_CALL (f, c, &fc);

  if (fc == 0.0)
    {
      *root = c;
      *x_lower = c;
      *x_upper = c;
      return GSL_SUCCESS;
    }

  /* Replace the endpoint with the same sign as f(c), remembering it as
     the most recently discarded point d */

  if ((fa > 0.0 && fc < 0.0) || (fa < 0.0 && fc > 0.0))
    {
      state->d = b;
      state->f_d = fb;
      *x_upper = c;
      state->f_upper = fc;
    }
  else
    {
      state->d = a;
      state->f_d = fa;
      *x_lower = c;
      state->f_lower = fc;
    }

  /* Bisect next if the bracket is shrinking more slowly than half the
     rate of bisection, then resume the sequence of steps */

  state->w_max *= M_SQRT1_2;

  if (state->step != STEP_BISECT && *x_upper - *x_lower > state->w_max)
    {
      state->resume = next;
      next = STEP_BISECT;
    }

  state->step = next;

  /* The best estimate of the root is the endpoint with the smaller
     function value */

  *root = (fabs (state->f_lower) < fabs (state->f_upper)) ? *x_lower : *x_upper;

  return GSL_SUCCESS;
}


static const gsl_root_fsolver_type aps_type =
{
  "aps",                                 /* name */
  sizeof (aps_state_t),
  &aps_init,
  &aps_iterate
};

const gsl_root_fsolver_type  * gsl_root_fsolver_aps = &aps_type;
