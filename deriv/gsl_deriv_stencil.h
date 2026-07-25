/* deriv/gsl_deriv_stencil.h
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

#ifndef __GSL_DERIV_STENCIL_H__
#define __GSL_DERIV_STENCIL_H__
#include <gsl/gsl_math.h>

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

typedef enum
{
  GSL_DERIV_STENCIL_CENTRAL  = 0,
  GSL_DERIV_STENCIL_FORWARD  = 1,
  GSL_DERIV_STENCIL_BACKWARD = 2
} gsl_deriv_stencil_t;

int gsl_deriv_stencil_size(const unsigned int deriv_order,
                           const unsigned int method_order,
                           const gsl_deriv_stencil_t stencil_type,
                           size_t * length);

int gsl_deriv_stencil_fill(const gsl_deriv_stencil_t stencil_type,
                           const size_t n, double stencil[]);

__END_DECLS

#endif /* __GSL_DERIV_STENCIL_H__ */
