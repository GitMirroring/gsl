/* dd/internal_sqrt_source.c
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

/*
 * Compute the square root of a positive double-double number.
 *
 * Scale the argument by an even power of two so that
 *
 *     a = m * 2^e,       1/2 <= m < 2,
 *
 * with e even.  Let y0 = sqrt(m.hi) be the binary64 initial
 * approximation.  Form the complete double-double residual
 *
 *     r = m - y0^2
 *
 * using an error-free TwoProd for y0^2, and apply one Newton
 * correction
 *
 *     delta = r / (2*y0).
 *
 * The corrected result y0 + delta is then scaled by 2^(e/2).
 */
static gsl_dd
DD_WRAP(dd_sqrt)(const gsl_dd a)
{
  gsl_dd m, y, yy, r, corr;
  int e;
  double y0;

  if (a.hi == 0.0)
    {
      gsl_dd z;
      z.hi = 0.0;
      z.lo = 0.0;
      return z;
    }

  m = sqrt_prepare(a, &e);

  y0 = sqrt(m.hi);
  y.hi = y0;
  y.lo = 0.0;

  yy = DD_WRAP(two_prod)(y0, y0);
  r = dd_sub(m, yy);

  corr = DD_WRAP(dd_div_d)(r, 2.0 * y0);
  y = dd_add(y, corr);

  return dd_ldexp(y, e / 2);
}
