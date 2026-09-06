/* dd/dd.c
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

#include "internal_common.c"

gsl_dd gsl_dd_add(const gsl_dd a, const gsl_dd b)      { return dd_add(a, b); }
gsl_dd gsl_dd_add_d(const gsl_dd a, const double b)    { return dd_add_d(a, b); }
gsl_dd gsl_dd_sub(const gsl_dd a, const gsl_dd b)      { return dd_sub(a, b); }
gsl_dd gsl_dd_sub_d(const gsl_dd a, const double b)    { return dd_sub_d(a, b); }

#ifdef HAVE_DECL_FMA

#include "internal_fma.c"

gsl_dd gsl_dd_mul(const gsl_dd a, const gsl_dd b)      { return dd_mul_fma(a, b); }
gsl_dd gsl_dd_mul_d(const gsl_dd a, const double b)    { return dd_mul_d_fma(a, b); }
gsl_dd gsl_dd_div(const gsl_dd a, const gsl_dd b)      { return dd_div_fma(a, b); }
gsl_dd gsl_dd_div_d(const gsl_dd a, const double b)    { return dd_div_d_fma(a, b); }
gsl_dd gsl_dd_exp(const gsl_dd a)                      { return dd_exp_fma(a); }

#else

#include "internal_p.c"

gsl_dd gsl_dd_mul(const gsl_dd a, const gsl_dd b)      { return dd_mul_p(a, b); }
gsl_dd gsl_dd_mul_d(const gsl_dd a, const double b)    { return dd_mul_d_p(a, b); }
gsl_dd gsl_dd_div(const gsl_dd a, const gsl_dd b)      { return dd_div_p(a, b); }
gsl_dd gsl_dd_div_d(const gsl_dd a, const double b)    { return dd_div_d_p(a, b); }
gsl_dd gsl_dd_log(gsl_dd y)                            { return dd_log_p(y); }
gsl_dd gsl_dd_sqrt(const gsl_dd a)                     { return dd_sqrt_p(a); }
gsl_dd gsl_dd_exp(const gsl_dd z)                      { return dd_exp_p(z); }
gsl_dd gsl_dd_expm1(const gsl_dd z)                    { return dd_expm1_p(z); }

#endif
