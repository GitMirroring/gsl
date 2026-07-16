/* deriv/scale_generator.c
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
#include <gsl/gsl_deriv_constants.h>
#include <gsl/gsl_deriv_scale.h>

static int user_init(const gsl_deriv_parameters * params, gsl_deriv_scale_parameters * sparams);

/*
gsl_deriv_scale_fill()
  Fill array with nondimensional scale multipliers,

scale_multipliers[i] = max_scale / ratio^i

Inputs: params - scale parameters
        scales - (output) array of scale multipliers, length GSL_DERIV_MAX_SCALES
*/

int
gsl_deriv_scale_fill(const gsl_deriv_scale_parameters * params, double scales[])
{
  if (params->nscale > GSL_DERIV_MAX_SCALES)
    {
      GSL_ERROR ("nscale must be <= GSL_DERIV_MAX_SCALES", GSL_EBADLEN);
    }
  else
    {
      double factor = 1.0; /* ratio^{-i} */
      size_t i;

      for (i = 0; i < params->nscale; ++i)
        {
          scales[i] = params->max_scale * factor;
          factor /= params->ratio;
        }

      return GSL_SUCCESS;
    }
}

/* user specified scale multipliers */
static int
user_init(const gsl_deriv_parameters * params, gsl_deriv_scale_parameters * sparams)
{
  if (params->max_scale <= 0.0)
    {
      GSL_ERROR ("max_scale must be positive", GSL_EDOM);
    }
  else if (params->min_scale <= 0.0)
    {
      GSL_ERROR ("min_scale must be positive", GSL_EDOM);
    }
  else if (params->nscale < 2)
    {
      GSL_ERROR ("nscale must be >= 2", GSL_EDOM);
    }
  else
    {
      sparams->max_scale = params->max_scale;
      sparams->nscale = params->nscale;
      sparams->ratio = pow(params->max_scale / params->min_scale, 1.0 / (params->nscale - 1.0));
      return GSL_SUCCESS;
    }
}

static const gsl_deriv_scale_generator user_type =
{
  "user",
  user_init
};

const gsl_deriv_scale_generator * gsl_deriv_scale_user = &user_type;
