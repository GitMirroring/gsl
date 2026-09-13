/* specfunc/erf/erf_fallback.c
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
 * Adaptation of the ARM erf.
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

static int erf_fallback(const double x, gsl_sf_result * result);

#define ERF_LINEAR_COEF              0.128379167095512586316     /* ARM tuned coefficient */
#define ERF_TWO_MINUS_28             3.725290298461914062500e-09 /* 2^{-28} */

static const double PA[10] = {
  0.128379167095512558561, -0.376126389031828711751,
  0.112837916708965924073,  -0.026866170630075902892,
  0.00522397742864975922483,  -0.000854831222998997429866,
  0.000120546545021512873362, -1.49063150674988910036e-05,
  1.61264443490708235249e-06, -1.30742590566799656928e-07
};

static const double NA[5] = {
  0.128379167095512558561, -0.32504210724700149937,
  -0.0284817495755985104766, -0.00577027029648944159157,
  -2.37630166566501626084e-05
};

static const double DA[5] = {
  0.397917223959155352819, 0.0650222499887672944485,
  0.00508130628187576562776, 0.000132494738004321644526,
  -3.9602282787753681232e-06
};

static const double NB[7] = {
  -0.00236211856075265944077, 0.414856118683748331666,
  -0.372207876035701323847, 0.318346619901161753674,
  -0.110894694282396677476, 0.0354783043256182359371,
  -0.002166375594868790843
};

static const double DB[6] = {
  0.106420880400844228286, 0.540397917702171048937,
  0.0718286544141962662868, 0.126171219808761642112,
  0.0136370839120290507362, 0.011984499846799107417
};

static const double PC[16] = {
  0.0770998717435417701971, -0.236521122447293008495,
  0.295651403059370543325, -0.167535795078651794654,
  0.00615940451616571720689,  0.0471810322209430671458,
  -0.0213012326687174943785, -0.00362625151743203779503,
  0.0056989638892621322791, -0.000882133759048259582028,
  -0.000781880568597137135499, 0.000301688204790364776411,
  8.4021280013454630476e-05, -8.73252701598456612057e-05,
  2.46326714413737165731e-05, -2.55964166811578430023e-06
};

static const double PD[18] = {
  0.00467773498104724627655, -0.0413339707081636434371,
  0.165335882829129171556, -0.385783726357600831314,
  0.551119600143394672287, -0.418850694287877134414,
  -0.0293961559930560022147, 0.432530081665319143802,
  -0.407548680724586664947, 0.0270381624052067996711,
  0.263263799297783429587, -0.193642630502174645235,
  -0.0677531396120535506222, 0.218809436648854216401,
  -0.18160468454384678938, 0.0799438600635083845836,
  -0.0184279538809101230135, 0.00162231430381370653193
};

static const double PE[14] = {
  4.30277946368441051327e-06, -2.91890253884814524182e-05,
  9.48643329252567719121e-05, -0.000195809726719565921306,
  0.00028656960349963429965, -0.000313800091194331204521,
  0.000263548541257612410887, -0.000170086616863952561768,
  8.19633045953011311544e-05, -2.66523835297490039719e-05,
  3.54271160510927491568e-06, 1.49667727163638605003e-06,
  -9.12583414628413246765e-07, 1.6434702071850266132e-07
};

static const double PF[17] = {
  1.54172578620972537307e-08, -2.53964670939008739871e-07,
  2.03171534820294374431e-06, -1.04971225284679385109e-05,
  3.92781628171582532639e-05, -0.000113078253657390167359,
  0.000259536459565380423636, -0.000485026022056390189625,
  0.000746082121944444576385, -0.000946698312642611718333,
  0.000984800102486704104221, -0.000827208962269068964063,
  0.000546794702345263622899, -0.000273170178514800477895,
  9.66426095142485354437e-05, -2.15042743701411657445e-05,
  2.25579946402808531391e-06
};

/*
 * Implementation of erf(x), adapted from the Arm Optimized Routines implementation.
 *
 * The approximation is piecewise in |x|.  Very small arguments use the
 * leading linear term of the Maclaurin expansion.  The interval below
 * 0.84375 uses polynomial and rational approximations which preserve the
 * small magnitude of erf(x).  Around x = 1 a rational approximation is
 * evaluated in |x|-1.  For larger arguments, polynomials approximate the
 * complementary tail and erf(x) is reconstructed as +/- (1-r).
 *
 * For |x| >= 5.90625, erf(x) is sufficiently close to +/-1 that its
 * binary64 value is returned directly.
 *
 * The original Arm implementation uses FMA and integer manipulation of
 * the floating-point representation.  This adaptation uses only portable
 * C89 arithmetic.
 */

static int
erf_fallback(const double x, gsl_sf_result * result)
{
  const double split_constant = 0.845062911510467529297;
  const double ax = fabs(x);

  /*
   * Small and intermediate arguments.  In this region erf(x) is evaluated
   * directly rather than as 1-erfc(x), avoiding cancellation near zero.
   */
  if (ax < 0.84375)
    {
      /*
       * For |x| < 2^-28 the cubic and higher terms are negligible at
       * double precision.  Evaluate x * (2/sqrt(pi)).
       */
      if (ax < ERF_TWO_MINUS_28)
        {
          result->val = x + ERF_LINEAR_COEF * x;
          result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
          return GSL_SUCCESS;
        }

      if (ax < 0.5)
        {
          /* polynomial approximation in x^2 for |x| < 1/2. */
          const double x2 = x*x;
          const double r = PA[0] + x2*(PA[1] + x2*(PA[2] + x2*(PA[3] +
                                   x2*(PA[4] + x2*(PA[5] + x2*(PA[6] +
                                   x2*(PA[7] + x2*(PA[8] + x2*(PA[9])))))))));
          result->val = r*x + x;
          result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);
          return GSL_SUCCESS;
        }
      else
        {
          /* rational approximation in x^2 for 1/2 <= |x| < 0.84375. */
          const double x2 = x*x;
          const double P = NA[0] + x2*(NA[1] + x2*(NA[2] + x2*(NA[3] + x2*(NA[4]))));
          const double Q = 1.0 + x2*(DA[0] + x2*(DA[1] + x2*(DA[2] + x2*(DA[3] + x2*(DA[4])))));
          result->val = (P/Q)*x + x;
          result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);
          return GSL_SUCCESS;
        }
    }

  /*
   * Rational approximation centered at |x| = 1.  The constant split_constant contains
   * the leading part of erf(1), while P/Q approximates the correction.
   */
  if (ax < 1.25)
    {
      const double a = ax - 1.0;
      const double P = NB[0] + a*(NB[1] + a*(NB[2] + a*(NB[3] + a*(NB[4] + a*(NB[5] + a*(NB[6]))))));
      const double Q = 1.0 + a*(DB[0] + a*(DB[1] + a*(DB[2] + a*(DB[3] + a*(DB[4] + a*(DB[5]))))));
      const double y = split_constant + P/Q;
      result->val = x < 0.0 ? -y : y;
      result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }

  /*
   * For |x| >= 1.25, approximate the small complementary quantity
   * r = 1-erf(|x|) and reconstruct erf(x) using odd symmetry.
   */

  if (ax < 2.0)
    {
      /* polynomial in |x|-1.25. */
      const double a = ax - 1.25;
      const double r = PC[0] + a*(PC[1] + a*(PC[2] + a*(PC[3] + a*(PC[4] +
                               a*(PC[5] + a*(PC[6] + a*(PC[7] + a*(PC[8] +
                               a*(PC[9] + a*(PC[10] + a*(PC[11] + a*(PC[12] +
                               a*(PC[13] + a*(PC[14] + a*(PC[15])))))))))))))));
      result->val = x < 0.0 ? -1.0 + r : 1.0 - r;
      result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }

  if (ax < 3.25)
    {
      /* polynomial in |x|/2-1 */
      const double a = ax*0.5 - 1.0;
      const double r = PD[0] + a*(PD[1] + a*(PD[2] + a*(PD[3] + a*(PD[4] +
                               a*(PD[5] + a*(PD[6] + a*(PD[7] + a*(PD[8] +
                               a*(PD[9] + a*(PD[10] + a*(PD[11] + a*(PD[12] +
                               a*(PD[13] + a*(PD[14] + a*(PD[15] + a*(PD[16] + a*(PD[17])))))))))))))))));
      result->val = x < 0.0 ? -1.0 + r : 1.0 - r;
      result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }

  if (ax < 4.0)
    {
      /* polynomial in |x|-3.25 */
      const double a = ax - 3.25;
      const double r = PE[0] + a*(PE[1] + a*(PE[2] + a*(PE[3] + a*(PE[4] +
                               a*(PE[5] + a*(PE[6] + a*(PE[7] + a*(PE[8] +
                               a*(PE[9] + a*(PE[10] + a*(PE[11] + a*(PE[12] +
                               a*(PE[13])))))))))))));
      result->val = x < 0.0 ? -1.0 + r : 1.0 - r;
      result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }

  if (ax < 5.90625)
    {
      /* final tail polynomial in |x|/2-2 */
      const double a = ax*0.5 - 2.0;
      const double r = PF[0] + a*(PF[1] + a*(PF[2] + a*(PF[3] + a*(PF[4] +
                               a*(PF[5] + a*(PF[6] + a*(PF[7] + a*(PF[8] +
                               a*(PF[9] + a*(PF[10] + a*(PF[11] + a*(PF[12] +
                               a*(PF[13] + a*(PF[14] + a*(PF[15] + a*(PF[16]))))))))))))))));
      result->val = x < 0.0 ? -1.0 + r : 1.0 - r;
      result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }

  /* NaN check */
  if (x != x)
    {
      result->val = x;
      result->err = GSL_NAN;
      return GSL_SUCCESS;
    }

  /* the result has saturated to +/-1 */
  result->val = x > 0.0 ? 1.0 : -1.0;
  result->err = 0.5 * GSL_DBL_EPSILON;

  return GSL_SUCCESS;
}
