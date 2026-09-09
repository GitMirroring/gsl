/* dd/gsl_dd.h
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

#ifndef __GSL_DD_H__
#define __GSL_DD_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct
{
  double hi;
  double lo;
} gsl_dd;

gsl_dd gsl_dd_add(const gsl_dd a, const gsl_dd b);
gsl_dd gsl_dd_add_d(const gsl_dd a, const double b);
gsl_dd gsl_dd_sub(const gsl_dd a, const gsl_dd b);
gsl_dd gsl_dd_sub_d(const gsl_dd a, const double b);
gsl_dd gsl_dd_mul(const gsl_dd a, const gsl_dd b);
gsl_dd gsl_dd_mul_d(const gsl_dd a, const double b);
gsl_dd gsl_dd_div(const gsl_dd a, const gsl_dd b);
gsl_dd gsl_dd_div_d(const gsl_dd a, const double b);
gsl_dd gsl_dd_exp(const gsl_dd x);
gsl_dd gsl_dd_expm1(const gsl_dd x);
gsl_dd gsl_dd_log(const gsl_dd x);
gsl_dd gsl_dd_log1p(const gsl_dd x);
gsl_dd gsl_dd_log1pmx(const gsl_dd x);
gsl_dd gsl_dd_sqrt(const gsl_dd x);

__END_DECLS

#endif /* __GSL_DD_H__ */
