/* dd/internal_exp.c
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

/* common routines for double-double dd_exp() */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>

#include "internal.h"

static gsl_dd mul_pwr2_extred(const gsl_dd a, const double c);

/*
 * Multiply a double-double number by a power of two.
 *
 * The multiplier c is assumed to be an exactly representable power of
 * two.  Scaling the two components separately is therefore exact,
 * provided neither component overflows or underflows.
 */
static gsl_dd
mul_pwr2_extred(const gsl_dd a, const double c)
{
  gsl_dd r;
  r.hi = a.hi * c;
  r.lo = a.lo * c;
  return r;
}
