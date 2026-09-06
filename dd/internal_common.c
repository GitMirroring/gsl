/* dd/internal_common.c
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>

#include "internal.h"

static gsl_dd two_sum(const double a, const double b);
static gsl_dd fast_two_sum(const double a, const double b);
static gsl_dd dd_add(const gsl_dd a, const gsl_dd b);
static gsl_dd dd_add_d(const gsl_dd a, const double b);

/*
 * Error-free transformation of a+b into a leading double and its
 * rounding residual.
 */
static gsl_dd
two_sum(const double a, const double b)
{
  gsl_dd r;
  double s = a + b, bb = s - a;
  r.hi = s;
  r.lo = (a - (s - bb)) + (b - bb);
  return r;
}

static gsl_dd
fast_two_sum(const double a, const double b)
{
  gsl_dd r;
  r.hi = a + b;
  r.lo = b - (r.hi - a);
  return r;
}

static gsl_dd
dd_add(const gsl_dd a, const gsl_dd b)
{
  gsl_dd s, t, v;
  double c, w;

  s = two_sum(a.hi, b.hi);
  t = two_sum(a.lo, b.lo);

  c = s.lo + t.hi;
  v = fast_two_sum(s.hi, c);

  w = t.lo + v.lo;

  return fast_two_sum(v.hi, w);
}

static gsl_dd
dd_add_d(const gsl_dd a, const double b)
{
  gsl_dd s;
  s = two_sum(a.hi, b);
  s.lo += a.lo;
  return fast_two_sum(s.hi, s.lo);
}

/*
 * Subtract two double-double numbers using accurate double-word
 * subtraction and renormalize the result.
 */
static gsl_dd
dd_sub(const gsl_dd a, const gsl_dd b)
{
  gsl_dd s, t, v;
  double c, w;

  s = two_sum(a.hi, -b.hi);
  t = two_sum(a.lo, -b.lo);

  c = s.lo + t.hi;
  v = fast_two_sum(s.hi, c);

  w = t.lo + v.lo;

  return fast_two_sum(v.hi, w);
}

/*
 * Subtract a double from a double-double number and renormalize
 * the result.
 */
static gsl_dd
dd_sub_d(const gsl_dd a, const double b)
{
  gsl_dd s;

  s = two_sum(a.hi, -b);
  s.lo += a.lo;

  return fast_two_sum(s.hi, s.lo);
}

#include "internal_exp.c"
