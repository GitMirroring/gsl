/* deriv/stencil.c
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
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv_stencil.h>

/*
gsl_deriv_stencil_size()
  Compute size of the optimal stencil given derivative order,
method order, and stencil type

Inputs: deriv_order  - derivative order
        method_order - method order
        stencil_type - GSL_DERIV_STENCIL_{CENTRAL,BACKWARD,FORWARD}
*/

int
gsl_deriv_stencil_size(const unsigned int deriv_order,
                       const unsigned int method_order,
                       const gsl_deriv_stencil_t stencil_type,
                       size_t * length)
{
  if (stencil_type == GSL_DERIV_STENCIL_CENTRAL)
    {
      const unsigned int radius = (deriv_order + method_order - 1) / 2;
      *length = 2 * radius + 1;
    }
  else if (stencil_type == GSL_DERIV_STENCIL_FORWARD ||
           stencil_type == GSL_DERIV_STENCIL_BACKWARD)
    {
      *length = deriv_order + method_order;
    }
  else
    {
      GSL_ERROR("invalid stencil_type", GSL_EDOM);
    }

  return GSL_SUCCESS;
}

/*
gsl_deriv_stencil_fill()
  Fill stencil

Inputs: stencil_type - GSL_DERIV_STENCIL_xxx
        n            - length of stencil (output of gsl_deriv_stencil_size)
        stencil      - (output) dimensionless stencil corresponding to unit scale

Return: success/error
*/

int
gsl_deriv_stencil_fill(const gsl_deriv_stencil_t stencil_type,
                       const size_t n,
                       double stencil[])
{
  if (n == 0)
    {
      GSL_ERROR("stencil length must be > 0", GSL_EBADLEN);
    }
  else
    {
      int i;

      switch (stencil_type)
        {
          case GSL_DERIV_STENCIL_CENTRAL:
            {
              const int mid = (int)n / 2;

              if (!(n&1))
                {
                  GSL_ERROR("central stencils require odd length", GSL_EDOM);
                }

              for (i = 1; i <= mid; ++i)
                {
                  stencil[mid - i] = -(double)i;
                  stencil[mid + i] =  (double)i;
                }

              stencil[mid] = 0.0;

              break;
            }

          case GSL_DERIV_STENCIL_FORWARD:
            {
              for (i = 0; i < (int) n; ++i)
                stencil[i] = (double) i;

              break;
            }

          case GSL_DERIV_STENCIL_BACKWARD:
            {
              for (i = (int)n - 1; i > 0; --i)
                stencil[n - i - 1] = (double) (-i);

              stencil[n - 1] = 0.0;

              break;
            }

          default:
            GSL_ERROR("invalid stencil_type", GSL_EDOM);
            break;
        }

      return GSL_SUCCESS;
    }
}
