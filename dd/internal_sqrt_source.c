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
 * Compute the square root of a double-double number.
 *
 * Scale the argument by an even power of two so that
 *
 *     a = m * 2^e,       1/2 <= m < 2,
 *
 * with e even.  Let y0 = sqrt(m.hi), correctly rounded.  The residual
 *
 *     R = m - y0^2 = ((m.hi - p.hi) - p.lo) + m.lo,    p = y0^2 (TwoProd),
 *
 * is exact as a double-double: m.hi - y0^2 is representable for a
 * correctly rounded square root, and the final addition is a TwoSum.
 * With d = 2 y0, expand
 *
 *     sqrt(m) = y0 + R/d - (R/d)^2/d + ...,
 *
 * whose next term is below 2^-156 y0.  The first correction
 * c1 = R.hi/d is rounded, so its remainder (R - c1 d), exact by TwoProd,
 * is folded into the second correction
 *
 *     c2 = ((R - c1 d) - c1^2) / d.
 *
 * y0 + c1 is split exactly by a fast TwoSum, c2 is added to its low part
 * and the result is scaled by 2^(e/2).
 *
 * Zero arguments (of either sign) are returned unchanged, +Inf gives
 * +Inf, and negative arguments or NaN give NaN.
 */
static gsl_dd
DD_WRAP(dd_sqrt)(const gsl_dd a)
{
  gsl_dd m, p, R, y;
  int e;
  double y0, d, c1, c2;

  if (!(a.hi > 0.0 && a.hi < GSL_POSINF))
    {
      gsl_dd z;
      z.hi = (a.hi == 0.0 || a.hi == GSL_POSINF) ? a.hi : GSL_NAN;
      z.lo = 0.0;
      return z;
    }

  m = sqrt_prepare(a, &e);

  y0 = sqrt(m.hi);
  d = 2.0 * y0;

  p = DD_WRAP(two_prod)(y0, y0);
  R = two_sum((m.hi - p.hi) - p.lo, m.lo);

  c1 = R.hi / d;
  p = DD_WRAP(two_prod)(c1, d);
  c2 = ((((R.hi - p.hi) - p.lo) + R.lo) - c1 * c1) / d;

  y = fast_two_sum(y0, c1);
  y = fast_two_sum(y.hi, y.lo + c2);

  return dd_ldexp(y, e / 2);
}
