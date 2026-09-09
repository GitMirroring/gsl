/* dd/internal_log.c
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

/* common routines for double-double log module */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>

#include "internal.h"

static gsl_dd log_normalize(const gsl_dd z, int * e);
static int log_center_index(const gsl_dd m);

/*
 * Normalize a positive double-double number z as
 *
 *   z = 2^e m,    1 <= m <= 2.
 *
 * Scale both components of z by the same exact binary power so that the
 * double-double value of the significand is preserved.  The initial exponent
 * is obtained from z.hi with frexp(); the result is then corrected when the
 * low component places the complete double-double value just below 1 or
 * just above 2.
 *
 * On return, *e contains the binary exponent and the returned double-double
 * value is the normalized significand m.
 */
static gsl_dd
log_normalize(const gsl_dd z, int * e)
{
  gsl_dd m;

  (void) frexp(z.hi, e);
  m.hi = ldexp(z.hi, 1 - *e);
  m.lo = ldexp(z.lo, 1 - *e);
  --(*e);

  if (DD_LT_d(m, 1.0))
    {
      m.hi *= 2.0;
      m.lo *= 2.0;
      --(*e);
    }
  else if (DD_GT_d(m, 2.0))
    {
      m.hi *= 0.5;
      m.lo *= 0.5;
      ++(*e);
    }

  return m;
}

/*
 * Select the range-reduction center for a normalized double-double
 * significand m, where
 *
 *   1 <= m <= 2.
 *
 * The interval [1,2] is divided into 128 exact dyadic subintervals
 *
 *   [1 + i/128, 1 + (i+1)/128),    i = 0,...,127.
 *
 * Return the index i of the interval containing the complete double-double
 * value m.  Since the initial index is computed from m.hi, explicitly correct
 * the result when m.hi lies exactly on an interval boundary and m.lo < 0,
 * in which case the full value belongs to the preceding interval.
 *
 * The endpoint m = 2 is assigned to the final interval.
 */
static int
log_center_index(const gsl_dd m)
{
  double q;
  int i;

  q = (m.hi - 1.0) * 128.0;
  i = (int) floor(q);

  if (i == 128)
    i = 127;

  if (i > 0 && m.hi == 1.0 + (double) i / 128.0 && m.lo < 0.0)
    --i;

  return i;
}
