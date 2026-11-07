/* dd/internal_trig.c
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

/* common routines for double-double sin(pi x), cos(pi x) */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>

#include "internal.h"

static int trig_pi_reduce(const gsl_dd x, gsl_dd * r);

/*
 * Reduce a finite double-double x for sin(pi x) and cos(pi x).
 *
 * Both functions have period 2, so each component is first reduced
 * modulo 2 with fmod(), which is exact, and the two remainders are
 * added exactly (TwoSum), giving y = x mod 2 with |y| < 4.  Then
 *
 *   y = q/2 + r,    q the integer nearest to 2 y.hi,
 *
 * where y.hi - q/2 is exact (Sterbenz, or q = 0) and the sum with y.lo
 * is a TwoSum, so the double-double r is exact and |r| <= 1/4 + tiny.
 * Since the reduction is exact, sin(pi x) and cos(pi x) keep their
 * relative accuracy near their zeros.
 *
 * Returns the quadrant q mod 4 in {0,1,2,3}.
 */
static int
trig_pi_reduce(const gsl_dd x, gsl_dd * r)
{
  gsl_dd y;
  double q;
  int quad;

  y = two_sum(fmod(x.hi, 2.0), fmod(x.lo, 2.0));

  q = floor(2.0 * y.hi + 0.5);
  *r = two_sum(y.hi - 0.5 * q, y.lo);

  quad = ((int) q) % 4;
  if (quad < 0)
    quad += 4;

  return quad;
}
