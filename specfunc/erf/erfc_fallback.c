/* specfunc/erf/erfc_fallback.c
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

/*
 * Adaptation of the ARM erfc.
 *
 * Copyright (c) 2020-2024, Arm Limited.
 *
 * Source algorithm and coefficients derived from Arm Optimized Routines
 * math/erf.c and math/erf_data.c
 *
 * SPDX-License-Identifier: MIT OR Apache-2.0 WITH LLVM-exception
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>

#include "constants.h"

static int erfc_fallback(const double x, gsl_sf_result * result);
static void split_exp_term(const double x, const double rs, gsl_sf_result * result);

/*
 * Compute the complementary error function erfc(x).
 *
 * For |x| < 0.25, compute erfc(x) from gsl_sf_erf_e().  For larger
 * arguments, use piecewise rational approximations derived from FDLIBM.
 * The exponentially small positive tail is evaluated with a split
 * exponential to avoid loss of accuracy in forming x*x.  For sufficiently
 * large negative arguments the result rounds to 2, while for sufficiently
 * large positive arguments it underflows to 0.
 *
 * The error estimate includes the propagated error from gsl_sf_erf_e()
 * or split_exp_term(), as appropriate, together with rounding error from
 * the final arithmetic.
 */
static int
erfc_fallback(const double x, gsl_sf_result * result)
{
  int status;
  const double ax = fabs(x);
  const double erx = 0.845062911510467529296875;
  static const double pp[5] = { 1.28379167095512558561e-01, -3.25042107247001499370e-01,
                               -2.84817495755985104766e-02, -5.77027029648944159157e-03,
                               -2.37630166566501626084e-05 };
  static const double qq[6] = { 0.0, 3.97917223959155352819e-01, 6.50222499887672944485e-02,
                                5.08130628187576562776e-03, 1.32494738004321644526e-04,
                               -3.96022827877536812320e-06 };
  static const double pa[7] = { -2.36211856075265944077e-03, 4.14856118683748331666e-01,
                                -3.72207876035701323847e-01, 3.18346619901161753674e-01,
                                -1.10894694282396677476e-01, 3.54783043256182359371e-02,
                                -2.16637559486879084300e-03 };
  static const double qa[7] = { 0.0, 1.06420880400844228286e-01, 5.40397917702171048937e-01,
                                7.18286544141962662868e-02, 1.26171219808761642112e-01,
                                1.36370839120290507362e-02, 1.19844998467991074170e-02 };
  static const double ra[8] = { -9.86494403484714822705e-03, -6.93858572707181764372e-01,
                                -1.05586262253232909814e+01, -6.23753324503260060396e+01,
                                -1.62396669462573470355e+02, -1.84605092906711035994e+02,
                                -8.12874355063065934246e+01, -9.81432934416914548592e+00 };
  static const double sa[9] = { 0.0, 1.96512716674392571292e+01, 1.37657754143519042600e+02,
                                4.34565877475229228821e+02, 6.45387271733267880336e+02,
                                4.29008140027567833386e+02, 1.08635005541779435134e+02,
                                6.57024977031928170135e+00, -6.04244152148580987438e-02 };
  static const double rb[7] = { -9.86494292470009928597e-03, -7.99283237680523006574e-01,
                                -1.77579549177547519889e+01, -1.60636384855821916062e+02,
                                -6.37566443368389627722e+02, -1.02509513161107724954e+03,
                                -4.83519191608651397019e+02 };
  static const double sb[8] = { 0.0, 3.03380607434824582924e+01, 3.25792512996573918826e+02,
                                1.53672958608443695994e+03, 3.19985821950859553908e+03,
                                2.55305040643316442583e+03, 4.74528541206955367215e+02,
                               -2.24409524465858183362e+01 };
  gsl_sf_result tmp;

  if (x > GSL_DBL_MAX / 2.0)
    {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  if (x < -GSL_DBL_MAX / 2.0)
    {
      result->val = 2.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  if (ax < 0.25)
    {
      if (x < 0.0)
        {
          status = gsl_sf_erf_e(-x, &tmp);
          if (status)
            return status;

          result->val = 1.0 + tmp.val;
        }
      else
        {
          status = gsl_sf_erf_e(x, &tmp);
          if (status)
            return status;

          result->val = 1.0 - tmp.val;
        }

      result->err = tmp.err + 2.0 * GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }

  if (ax < 0.84375)
    {
      const double z = x * x;
      const double z2 = z * z;
      const double z4 = z2 * z2;
      const double r = (pp[0] + z * pp[1]) +
                       z2 * (pp[2] + z * pp[3]) + z4 * pp[4];
      const double s = (1.0 + z * qq[1]) +
                       z2 * (qq[2] + z * qq[3]) + z4 * (qq[4] + z * qq[5]);
      const double y = r / s;

      if (x < 0.25)
        {
          result->val = 1.0 - (x + x * y);
        }
      else
        {
          const double t = x * y + (x - 0.5);
          result->val = 0.5 - t;
        }

      result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }

  if (ax < 1.25)
    {
      const double s = ax - 1.0;
      const double P = pa[0] + s * (pa[1] + s * (pa[2] +
                        s * (pa[3] + s * (pa[4] + s * (pa[5] + s * pa[6])))));
      const double Q = 1.0 + s * (qa[1] + s * (qa[2] + s * (qa[3] +
                        s * (qa[4] + s * (qa[5] + s * qa[6])))));
      const double v = P / Q;

      if (x > 0.0)
        result->val = (1.0 - erx) - v;
      else
        result->val = 1.0 + (erx + v);

      result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }

  if (ax < 28.0)
    {
      const double s = 1.0 / (x * x);
      const double axv = ax;
      double R, S;

      if (axv < 2.85711669921875)
        {
          R = ra[0] + s * (ra[1] + s * (ra[2] + s * (ra[3] +
              s * (ra[4] + s * (ra[5] + s * (ra[6] + s * ra[7]))))));
          S = 1.0 + s * (sa[1] + s * (sa[2] + s * (sa[3] + s * (sa[4] + s * (sa[5] + s * (sa[6] + s * (sa[7] + s * sa[8])))))));
        }
      else
        {
          if (x < 0.0 && axv >= 6.0)
            {
              result->val = 2.0;
              result->err = GSL_DBL_EPSILON * fabs(result->val);
              return GSL_SUCCESS;
            }

          R = rb[0] + s * (rb[1] + s * (rb[2] + s * (rb[3] +
              s * (rb[4] + s * (rb[5] + s * rb[6])))));
          S = 1.0 + s * (sb[1] + s * (sb[2] + s * (sb[3] +
              s * (sb[4] + s * (sb[5] + s * (sb[6] + s * sb[7]))))));
        }

      split_exp_term(axv, R / S, &tmp);

      {
        const double r = tmp.val / axv;
        result->val = x > 0.0 ? r : 2.0 - r;
      }

      if (x > 0.0)
        result->err = tmp.err / axv + 2.0 * GSL_DBL_EPSILON * fabs(result->val) + GSL_SF_TRUE_MIN;
      else
        result->err = tmp.err / axv + 2.0 * GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }

  result->val = x > 0.0 ? 0.0 : 2.0;
  result->err = 0.0;

  return GSL_SUCCESS;
}

/*********************************************
 * INTERNAL ROUTINES                         *
 *********************************************/

/* Compute exp(-x*x - 9/16 + r/s) with the low bits of x*x controlled. */
static void
split_exp_term(const double x, const double rs, gsl_sf_result * result)
{
  int e;
  const double f = frexp(x, &e);
  const double hi = ldexp(floor(ldexp(f, 20.0)), e - 20);
  const double e1 = -hi * hi - 0.5625;
  const double e2 = (hi - x) * (hi + x) + rs;

  result->val = exp(e1) * exp(e2);
  result->err = (4.0 + x * x) * GSL_DBL_EPSILON * fabs(result->val) + GSL_SF_TRUE_MIN;
}
