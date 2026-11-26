/* sys/frexp_fallback.c
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

#include "ldfrexp_table.h"

/*
 * Fallback for frexp: x = f * 2^e with 1/2 <= |f| < 1, exactly, using
 * only multiplications and divisions by exact powers of two (no libm
 * calls).
 *
 * A subnormal |x| is first scaled by 2^54 into the normal range, and
 * |x| < 1 by 2^512 (at most twice) into [1, 2^1024).  Binary searches in
 * the tables then give 2^(32 j) <= |x| < 2^(32 (j+1)) and 2^i <= |x| / 2^(32 j)
 * < 2^(i+1); every division is exact.  Zero, infinities and NaN are
 * returned unchanged, with e = 0.
 *
 * Used where the system frexp does not handle subnormal arguments
 * (configure: HAVE_BROKEN_FREXP).
 */

static double
frexp_fallback (const double x, int *e)
{
  double y = fabs (x);
  int n = 0, i = 0, j = 0, step;

  if (x == 0.0 || !gsl_finite (x))
    {
      *e = 0;
      return x;
    }

  if (y < GSL_DBL_MIN)
    {
      y *= LDFREXP_2_54;
      n -= 54;
    }

  while (y < 1.0)
    {
      y *= LDFREXP_2_512;
      n -= 512;
    }

  /* largest j with 2^(32 j) <= y, then largest i with 2^i <= y / 2^(32 j) */
  for (step = 16; step > 0; step >>= 1)
    {
      if (ldfrexp_hi[j + step] <= y)
        j += step;
    }

  y /= ldfrexp_hi[j];

  for (step = 16; step > 0; step >>= 1)
    {
      if (ldfrexp_lo[i + step] <= y)
        i += step;
    }

  y /= ldfrexp_lo[i];           /* 1 <= y < 2 */

  *e = n + 32 * j + i + 1;

  return (x < 0.0) ? -0.5 * y : 0.5 * y;
}
