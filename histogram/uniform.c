/* histogram/uniform.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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

#include <gsl/gsl_sys.h>

/* Set range[i] = xmin + i (xmax - xmin) / n for i = 0..n, with
   range[0] = xmin and range[n] = xmax exactly. The edges are computed
   as xmin + ((xmax - xmin) i) / n, which is monotone in i and exact
   whenever xmin + i dx is representable (e.g. integer bin edges), so
   that integer data fall into the expected bins (bug #38548). If
   xmax - xmin overflows, fall back to the weighted mean of xmin and
   xmax. */

static void
make_uniform (double range[], size_t n, double xmin, double xmax)
{
  const double w = xmax - xmin;
  size_t i;

  for (i = 0; i < n; i++)
    {
      if (gsl_finite (w))
        {
          range[i] = xmin + (w * (double) i) / (double) n;
        }
      else
        {
          double f1 = ((double) (n-i) / (double) n);
          double f2 = ((double) i / (double) n);
          range[i] = f1 * xmin +  f2 * xmax;
        }
    }

  range[n] = xmax;
}
