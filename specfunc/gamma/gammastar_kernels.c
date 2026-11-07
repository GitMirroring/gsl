/* specfunc/gamma/gammastar_kernels.c
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
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_gamma.h>

#include "constants.h"
#include "error.h"
#include "stirling.c"

#include "lngamma1p_kernels.c"

#define GAMMASTAR_X_SMALL 0.4
#define GAMMASTAR_X_MID   2.0
#define GAMMASTAR_X_ASYM  16.0

/*
 * Evaluate Gamma*(x) for 0 < x < 0.4.
 *
 * Define
 *
 *   L = log(Gamma(1+x)) + x (1-log(x)).
 *
 * Then
 *
 *   Gamma*(x) = [1+expm1(L)] / [sqrt(2*pi) sqrt(x)].
 *
 * Splitting the square root avoids multiplying 2*pi by very small x
 * before the square root, which is important in the subnormal range.
 */
static int
gammastar_small_e(const double x, gsl_sf_result * result)
{
  gsl_sf_result lg, lx, ex;
  double u, eu, l, el, f, ef, root, eroot;
  int status;

  status = lngamma1p_near0(x, &lg);
  if (status)
    return status;

  status = gsl_sf_log_e(x, &lx);
  if (status)
    return status;

  u = x * (1.0 - lx.val);
  eu = x * lx.err + 3.0 * GSL_DBL_EPSILON * fabs(u);
  l = lg.val + u;
  el = lg.err + eu + 2.0 * GSL_DBL_EPSILON * fabs(l);

  status = gsl_sf_expm1_e(l, &ex);
  if (status)
    return status;

  f = 1.0 + ex.val;
  ef = ex.err + GSL_DBL_EPSILON * fabs(f);

  root = GSL_SF_SQRT_2PI * sqrt(x);
  eroot = 3.0 * GSL_DBL_EPSILON * fabs(root);

  result->val = f / root;
  result->err = fabs(result->val) *
                (el + ef / fabs(f) + eroot / fabs(root) +
                 2.0 * GSL_DBL_EPSILON);

  if (result->err < 12.0 * GSL_DBL_EPSILON * fabs(result->val))
    result->err = 12.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/*
 * Evaluate Gamma*(x) for 0.4 <= x < 2 using a degree-38
 * Chebyshev approximation.
 *
 * The coefficients represent Gamma*(x) directly on [0.4,2]
 *
 *   Gamma*(x) = c[0]/2 + sum_{k=1}^38 c[k] T_k(t),
 *
 * where
 *
 *   t = (2*x - 2.4)/1.6.
 *
 * The approximation and its rounding error are evaluated by
 * cheb_eval_e().
 */

static int
gammastar_mid1_e(const double x, gsl_sf_result * result)
{
  static double coef[39] = { 2.186686780025690613e+00, -7.013869256961335896e-02,
                             2.599843938309164687e-02, -9.558620670091954458e-03,
                             3.499245449414224805e-03, -1.278908632396545545e-03,
                             4.674617312361324249e-04, -1.710631624987272005e-04,
                             6.270876772690214861e-05, -2.303468480452105255e-05,
                             8.479166117853733693e-06, -3.127694049473782844e-06,
                             1.155981962090224919e-06, -4.280315040887567758e-07,
                             1.587582122076981955e-07, -5.897542969353828013e-08,
                             2.193928956304014263e-08, -8.172167173542489818e-09,
                             3.047662324935399382e-09, -1.137805282559345296e-09,
                             4.252082431145215255e-10, -1.590498671308784199e-10,
                             5.954315283179392074e-11, -2.230845389525679870e-11,
                             8.364141595775561608e-12, -3.138088601191735840e-12,
                             1.178096411035275128e-12, -4.425375339935861723e-13,
                             1.663246422352324927e-13, -6.254395947917066545e-14,
                             2.353005025023755814e-14, -8.856388494631339073e-15,
                             3.334843483052709215e-15, -1.256229125992642869e-15,
                             4.733989874623061019e-16, -1.784604843799515136e-16,
                             6.729843708664810407e-17, -2.538683504736716639e-17,
                             9.579566819709661330e-18 };
  static const cheb_series coef_cs = { coef, 38, 0.4, 2.0, 38 };

  return cheb_eval_e(&coef_cs, x, result);
}

/*
 * Evaluate Gamma*(x) for 2 <= x < 16 using a degree-45
 * Chebyshev approximation.
 *
 * The coefficients represent Gamma*(x) directly on [2,16]
 *
 *   Gamma*(x) = c[0]/2 + sum_{k=1}^45 c[k] T_k(t),
 *
 * where
 *
 *   t = (x - 9)/7.
 *
 * The approximation and its rounding error are evaluated by
 * cheb_eval_e().
 */

static int
gammastar_mid2_e(const double x, gsl_sf_result * result)
{
  static double coef[46] = { 2.029711612003090693e+00, -1.425348800128622799e-02,
                             6.831202492875235177e-03, -3.270920477618201577e-03,
                             1.564791556416853836e-03, -7.479531994260010894e-04,
                             3.572252799844768695e-04, -1.704820479771400408e-04,
                             8.130232913782237255e-05, -3.874667605252050885e-05,
                             1.845407430010521547e-05, -8.784063282676991821e-06,
                             4.178907340023351714e-06, -1.987066573049513921e-06,
                             9.444136742982740340e-07, -4.486725298447965377e-07,
                             2.130740676458901484e-07, -1.011536245333627036e-07,
                             4.800617775626012549e-08, -2.277676073763700992e-08,
                             1.080387324005711082e-08, -5.123567199213609583e-09,
                             2.429309146005882619e-09, -1.151653145831406987e-09,
                             5.458830115966487809e-10, -2.587176631528564771e-10,
                             1.226055491236473584e-10, -5.809783442729140457e-11,
                             2.752853415130225007e-11, -1.304326944434600211e-11,
                             6.179828003196381780e-12, -2.927915803254276262e-12,
                             1.387197327946146587e-12, -6.572340433900616555e-13,
                             3.113925483433529784e-13, -1.475388883021507849e-13,
                             6.990661673217008433e-14, -3.312428625933650169e-14,
                             1.569617652642980568e-14, -7.438110591372987534e-15,
                             3.524965065841918769e-15, -1.670598860427790372e-15,
                             7.918013266505274853e-16, -3.753083528685192606e-16,
                             1.779054588170992217e-16, -8.433740846689962961e-17 };
  static const cheb_series coef_cs = { coef, 45, 2.0, 16.0, 45 };

  return cheb_eval_e(&coef_cs, x, result);
}

/*
 * Evaluate Gamma*(x) = exp(Delta(x)) for x >= 16, with Delta the
 * Stirling correction (stirling.c).
 */
static int
gammastar_asym_e(const double x, gsl_sf_result * result)
{
  result->val = exp(gamma_stirling_delta(x));
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}
