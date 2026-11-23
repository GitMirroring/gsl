/* cdf/fdistinv.c
 *
 * Copyright (C) 2002 Przemyslaw Sliwa and Jason H. Stover.
 * Copyright (C) 2006, 2007 Brian Gough
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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include "error.h"

double
gsl_cdf_fdist_Pinv (const double P, const double nu1, const double nu2)
{
  double result;
  double y, w;

  if (!(P >= 0.0 && P <= 1.0))
    {
      CDF_ERROR ("P must be in range 0 <= P <= 1", GSL_EDOM);
    }
  if (!(nu1 > 0.0))
    {
      CDF_ERROR ("nu1 must be positive", GSL_EDOM);
    }
  if (!(nu2 > 0.0))
    {
      CDF_ERROR ("nu2 must be positive", GSL_EDOM);
    }

  /* x = (nu2/nu1) y/(1-y) with y the inverse of the beta distribution
     B(nu1/2, nu2/2) at P; for y > 1/2, 1 - y is computed directly as the
     inverse of B(nu2/2, nu1/2) at Q = 1 - P, to avoid cancellation */
  y = gsl_cdf_beta_Pinv (P, nu1 / 2.0, nu2 / 2.0);

  if (y <= 0.5)
    {
      result = nu2 * y / (nu1 * (1.0 - y));
    }
  else
    {
      w = gsl_cdf_beta_Qinv (P, nu2 / 2.0, nu1 / 2.0);

      result = nu2 * (1.0 - w) / (nu1 * w);
    }

  return result;
}

double
gsl_cdf_fdist_Qinv (const double Q, const double nu1, const double nu2)
{
  double result;
  double y, w;

  if (!(Q >= 0.0 && Q <= 1.0))
    {
      CDF_ERROR ("Q must be in range 0 <= Q <= 1", GSL_EDOM);
    }
  if (!(nu1 > 0.0))
    {
      CDF_ERROR ("nu1 must be positive", GSL_EDOM);
    }
  if (!(nu2 > 0.0))
    {
      CDF_ERROR ("nu2 must be positive", GSL_EDOM);
    }

  /* as for gsl_cdf_fdist_Pinv, with Q in place of P */
  y = gsl_cdf_beta_Qinv (Q, nu1 / 2.0, nu2 / 2.0);

  if (y <= 0.5)
    {
      result = nu2 * y / (nu1 * (1.0 - y));
    }
  else
    {
      w = gsl_cdf_beta_Pinv (Q, nu2 / 2.0, nu1 / 2.0);

      result = nu2 * (1.0 - w) / (nu1 * w);
    }

  return result;
}
