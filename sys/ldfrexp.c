/* sys/ldfrexp.c
 * 
 * Copyright (C) 2002, Gert Van den Eynde
 * Copyright (C) 2007, Brian Gough
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

/* gsl_ldexp and gsl_frexp use the system ldexp and frexp, unless
   configure found them missing or inaccurate for subnormal values (as
   NetBSD's ldexp, which truncates subnormal results); then the exact
   fallbacks are used.  GSL calls these instead of ldexp and frexp. */

#if HAVE_DECL_LDEXP && !defined(HAVE_BROKEN_LDEXP)
#define USE_SYSTEM_LDEXP 1
#else
#include "ldexp_fallback.c"
#endif

#if HAVE_DECL_FREXP && !defined(HAVE_BROKEN_FREXP)
#define USE_SYSTEM_FREXP 1
#else
#include "frexp_fallback.c"
#endif

double
gsl_ldexp (const double x, const int e)
{
#ifdef USE_SYSTEM_LDEXP
  return ldexp (x, e);
#else
  return ldexp_fallback (x, e);
#endif
}

double
gsl_frexp (const double x, int *e)
{
#ifdef USE_SYSTEM_FREXP
  return frexp (x, e);
#else
  return frexp_fallback (x, e);
#endif
}
