/* specfunc/gamma/lngamma_sgn_fallback.c
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
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_sincos_pi.h>

#include "error.h"

static int lngamma_sgn_fallback_e(const double x, gsl_sf_result * result, double * sgn);

static int
lngamma_sgn_fallback_e(const double x, gsl_sf_result * result, double * sgn)
{
  if (x > 0.0)
    {
      *sgn = 1.0;
      return lngamma_fallback_e(x, result);
    }
  else
    {
      gsl_sf_result rs;
      int status;

      status = gsl_sf_sin_pi_e(x, &rs);
      if (status)
        {
          result->val = GSL_NAN;
          result->err = GSL_NAN;
          *sgn = GSL_NAN;
          return status;
        }

      if (rs.val == 0.0)
        {
          *sgn = 0.0;
          DOMAIN_ERROR(result);
        }

      *sgn = (rs.val < 0.0) ? -1.0 : 1.0;

      return lngamma_fallback_e(x, result);
    }
}
