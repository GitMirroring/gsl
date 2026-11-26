/* sys/ldexp_fallback.c
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
 * Fallback for ldexp: x * 2^e, correctly rounded, using only
 * multiplications and divisions by exact powers of two (no libm calls).
 *
 * The exponent is first brought into [-1022, 1023] in steps of 512; the
 * scaled value is then multiplied or divided by 2^|e|, an exact product
 * of two table entries, in a single rounded operation.  A step can round
 * only if its result is subnormal, and then the final result is below
 * 2^-1076 and rounds to zero either way.
 *
 * Used where the system ldexp does not round subnormal results
 * (configure: HAVE_BROKEN_LDEXP).
 */

static double
ldexp_fallback (const double x, const int e)
{
  double y = x, s;
  int n = e, m;

  if (x == 0.0 || !gsl_finite (x))
    {
      return x;
    }

  /* beyond these limits the result is zero or overflows anyway */
  if (n > 2200)
    n = 2200;
  else if (n < -2200)
    n = -2200;

  while (n > 1023)
    {
      y *= LDFREXP_2_512;
      n -= 512;
    }

  while (n < -1022)
    {
      y /= LDFREXP_2_512;
      n += 512;
    }

  m = (n < 0) ? -n : n;
  s = ldfrexp_lo[m & 31] * ldfrexp_hi[m >> 5];

  return (n < 0) ? y / s : y * s;
}
