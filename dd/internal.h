/* dd/internal.h
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

#ifndef __INTERNAL_H__
#define __INTERNAL_H__

#define DD_SPLIT                 134217729.0
#define DD_LN2_HI                6.93147180559945286227e-01 /* double-double expansion of log(2) */
#define DD_LN2_LO                2.31904681384629955842e-17
#define DD_LN2_XTRA              5.70770843841621206578e-34
#define DD_E_HI                  2.718281828459045091e+00   /* double-double expansion of e. */
#define DD_E_LO                  1.445646891729250158e-16

#define DD_CAT_(a,b)             a##b
#define DD_CAT(a,b)              DD_CAT_(a,b)

#define DD_GT_d(x, y)            ( (x.hi > y) || (x.hi == y && x.lo > 0.0) )
#define DD_LT_d(x, y)            ( (x.hi < y) || (x.hi == y && x.lo < 0.0) )

/*
 * Short floating-point expansions used for compensated accumulation.
 *
 * A three-component expansion represents the unevaluated sum
 *
 *   x = a + b + c,
 *
 * and a four-component expansion represents
 *
 *   x = a + b + c + d.
 *
 * New binary64 terms are accumulated with cascaded TwoSum operations,
 * propagating each rounding residual into the next component.  This
 * retains information that would be lost by rounding each intermediate
 * sum directly to double-double.
 *
 * These types are temporary accumulation formats, not general
 * triple-double or quad-double arithmetic types.  They are compressed
 * to a normalized double-double result after accumulation is complete.
 */

typedef struct
{
  double a;
  double b;
  double c;
} dd_ex3;

typedef struct
{
  double a;
  double b;
  double c;
  double d;
} dd_ex4;

/* internal_common.c */
static gsl_dd two_sum(const double a, const double b);
static gsl_dd fast_two_sum(const double a, const double b);
static gsl_dd dd_add(const gsl_dd a, const gsl_dd b);
static gsl_dd dd_add_d(const gsl_dd a, const double b);
static void ex3_add_d(dd_ex3 * x, const double q);
static void ex3_renorm(dd_ex3 * x);
static gsl_dd ex3_to_dd(const dd_ex3 x);
static void ex4_add_d(dd_ex4 * x, const double q);
static gsl_dd ex4_to_dd(const dd_ex4 x);

#endif /* __INTERNAL_H__ */
