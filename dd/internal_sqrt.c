/* dd/internal_sqrt.c
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

/* common routines for double-double dd_sqrt() */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>

#include "internal.h"

/*
 * Normalize a positive double-double number to m * 2^e with e even.
 * The returned mantissa satisfies 1/2 <= m < 2.
 */
static gsl_dd
sqrt_prepare(const gsl_dd a, int * e)
{
  gsl_dd m;
  int k;

  (void) frexp(a.hi, &k);

  if (k & 1)
    --k;

  m = dd_ldexp(a, -k);
  *e = k;

  return m;
}
