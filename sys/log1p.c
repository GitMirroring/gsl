/* sys/log1p.c
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

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>

#ifndef HAVE_LOG1P
#include "log1p_fallback.c"
#endif

/*
 * gsl_log1p()
 *   Compute log(1+x), accurately for small x.  This is a low-level
 *   primitive for use by other modules (gsl_acosh, gsl_asinh,
 *   gsl_atanh, the complex functions); it does not call the error
 *   handler.  Use the system log1p() when available, otherwise the C89
 *   fallback, which compensates for the rounding error in forming 1+x.
 *
 *   x = -1 gives -Inf, x = +Inf gives +Inf and x < -1 gives NaN.
 */
double
gsl_log1p (const double x)
{
  if (x < -1.0)
    return GSL_NAN;
  else if (x == -1.0)
    return GSL_NEGINF;
  else if (x > GSL_DBL_MAX)
    return x;                   /* +Inf, as in C99 (the fallback would give NaN) */

#ifdef HAVE_LOG1P
  return log1p (x);
#else
  return log1p_fallback (x);
#endif
}
