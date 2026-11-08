/* sys/log1p_fallback.c
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
 * C89 fallback for log(1+x), the same algorithm as the gsl_sf_log1p
 * fallback (specfunc/log/log1p_fallback.c).  It is repeated here so
 * that sys does not depend on specfunc.  The caller must ensure
 * x > -1 and x finite.
 */

#define LOG1P_TINY   5.5511151231257827021e-17   /* 2^-54 */
#define LOG1P_SMALL  1.86264514923095703125e-9   /* 2^-29 */

static double
log1p_fallback (const double x)
{
  const double ax = fabs (x);

  /*
   * For |x| < 2^-54, log(1+x) rounds to x.  Up to 2^-29,
   * x - x^2/2 avoids loss from forming 1+x and is sufficiently
   * accurate for binary64.
   */
  if (ax < LOG1P_TINY)
    {
      return x;
    }
  else if (ax < LOG1P_SMALL)
    {
      return x - 0.5 * x * x;
    }
  else
    {
      /* c is the rounding error in forming u = 1+x */
      const double u = 1.0 + x;
      const double c = (u - 1.0) - x;

      return log (u) - c / u;
    }
}
