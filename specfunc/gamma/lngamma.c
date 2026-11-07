/* specfunc/gamma/lngamma.c
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
#include "eval.h"

#include "lngamma1p_kernels.c"

#ifdef HAVE_LGAMMA
static int lngamma_system_e(const double x, gsl_sf_result * result);
#endif

#if !defined(HAVE_LGAMMA) || !defined(HAVE_LGAMMA_R)
#include "lngamma_fallback.c"
#endif

#ifdef HAVE_LGAMMA_R
static int lngamma_system_sgn_e(const double x, gsl_sf_result * result, double * sgn);
#else
#include "lngamma_sgn_fallback.c"
#endif

int
gsl_sf_lngamma_e(const double x, gsl_sf_result * result)
{
#ifdef HAVE_LGAMMA

  /* use system lgamma() */
  return lngamma_system_e(x, result);

#else

  /* use local fallback function */
  return lngamma_fallback_e(x, result);

#endif
}

int
gsl_sf_lngamma_sgn_e(const double x, gsl_sf_result * result, double * sgn)
{
#ifdef HAVE_LGAMMA_R

  /* use system lgamma_r() */
  return lngamma_system_sgn_e(x, result, sgn);

#else

  /* use local fallback function */
  return lngamma_sgn_fallback_e(x, result, sgn);

#endif
}

/*
 * Compute log(Gamma(1+x)) for x > -1.
 *
 * Dedicated zero-preserving rational approximations are used around
 * the zeros at x=0 and x=1.  For negative arguments away from x=0,
 * recurrence is used:
 *
 *   log(Gamma(1+x)) = log(Gamma(2+x)) - log(1+x).
 *
 * The regions are:
 *
 *   -1 < x < LNGAMMA_1P_X_RECUR0_RECUR1
 *       recurrence through the near x=0 kernel, with t=1+x
 *
 *   LNGAMMA_1P_X_RECUR0_RECUR1
 *       <= x < LNGAMMA_1P_X_RECUR1_NEAR0
 *       recurrence through the near x=1 kernel, supplying its local
 *       coordinate directly to avoid forming (1+x)-1
 *
 *   LNGAMMA_1P_X_RECUR1_NEAR0
 *       <= x < LNGAMMA_1P_X_NEAR0_NEAR1
 *       near x=0 rational approximation
 *
 *   LNGAMMA_1P_X_NEAR0_NEAR1
 *       <= x <= LNGAMMA_1P_X_NEAR1_LNGAMMA
 *       near x=1 rational approximation
 *
 *   x > LNGAMMA_1P_X_NEAR1_LNGAMMA
 *       evaluate log(Gamma(1+x)) with the general lngamma routine
 *
 * The function has a domain error for x <= -1.
 */
int
gsl_sf_lngamma1p_e(const double x, gsl_sf_result * result)
{
  if (x <= -1.0)
    {
      DOMAIN_ERROR(result);
    }
  else if (x < LNGAMMA_1P_X_RECUR1_NEAR0)
    {
      return lngamma1p_recurrence(x, result);
    }
  else if (x < LNGAMMA_1P_X_NEAR0_NEAR1)
    {
      return lngamma1p_near0(x, result);
    }
  else if (x <= LNGAMMA_1P_X_NEAR1_LNGAMMA)
    {
      return lngamma1p_near1(x, result);
    }
  else
    {
      return gsl_sf_lngamma_e(1.0 + x, result);
    }
}

double gsl_sf_lngamma(const double x) { EVAL_RESULT(gsl_sf_lngamma_e(x, &result)); }
double gsl_sf_lngamma1p(const double x) { EVAL_RESULT(gsl_sf_lngamma1p_e(x, &result)); }

/********************************************
 * INTERNAL ROUTINES FOR LNGAMMA            *
 ********************************************/

#ifdef HAVE_LGAMMA

static int
lngamma_system_e(const double x, gsl_sf_result * result)
{
  if (x <= 0.0 && x == floor(x))
    {
      DOMAIN_ERROR(result);
    }

  result->val = lgamma(x);
  result->err = 8.0 * GSL_DBL_EPSILON * (fabs(result->val) + 1.0);

  return GSL_SUCCESS;
}

#endif /* HAVE_LGAMMA */

#ifdef HAVE_LGAMMA_R

static int
lngamma_system_sgn_e(const double x, gsl_sf_result * result, double * sgn)
{
  int signp;

  if (x <= 0.0 && x == floor(x))
    {
      *sgn = 0.0;
      DOMAIN_ERROR(result);
    }

  result->val = lgamma_r(x, &signp);
  result->err = 8.0 * GSL_DBL_EPSILON * (fabs(result->val) + 1.0);

  *sgn = signp;

  return GSL_SUCCESS;
}

#endif /* HAVE_LGAMMA_R */
