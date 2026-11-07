/* specfunc/gamma/stirling.c
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
 * References: Stirling's asymptotic series, NIST Digital Library of
 * Mathematical Functions, eq. 5.11.1; the coefficients B_{2k}/(2k(2k-1))
 * are mathematical constants.
 */

#ifndef __GSL_SF_GAMMA_STIRLING_C__
#define __GSL_SF_GAMMA_STIRLING_C__

/* smallest argument for which gamma_stirling_delta() is accurate */
#define GAMMA_STIRLING_X_MIN 8.0

/*
 * Stirling correction
 *
 *   Delta(x) = log(Gamma(x)) - [(x - 1/2) log(x) - x + log(sqrt(2 pi))],
 *
 * so that Gamma*(x) = exp(Delta(x)), for x >= GAMMA_STIRLING_X_MIN.
 *
 * Uses the asymptotic series
 *
 *   Delta(x) = sum_{k>=1} B_{2k} / (2k (2k-1) x^{2k-1})
 *            = 1/(12x) - 1/(360x^3) + 1/(1260x^5) - ...,
 *
 * truncated after 11 terms and evaluated by Horner's method in 1/x^2.
 * For x >= 8 the truncation error is below the rounding error, and the
 * result is accurate to about 1 ulp relative to Delta(x).
 */
static double
gamma_stirling_delta(const double x)
{
  /* B_{2k} / (2k (2k-1)), k = 1, ..., 11 */
  static const double c[11] = {  8.33333333333333287e-02,   /*  1/12            */
                                -2.77777777777777788e-03,   /* -1/360           */
                                 7.93650793650793650e-04,   /*  1/1260          */
                                -5.95238095238095292e-04,   /* -1/1680          */
                                 8.41750841750841714e-04,   /*  1/1188          */
                                -1.91752691752691763e-03,   /* -691/360360      */
                                 6.41025641025641003e-03,   /*  1/156           */
                                -2.95506535947712423e-02,   /* -3617/122400     */
                                 1.79644372368830574e-01,   /*  43867/244188    */
                                -1.39243221690590113e+00,   /* -174611/125400   */
                                 1.34028640441683926e+01 }; /*  77683/5796      */
  const double z = 1.0 / (x * x);
  double p = c[10];
  int k;

  for (k = 9; k >= 0; --k)
    p = p * z + c[k];

  return p / x;
}

#endif /* __GSL_SF_GAMMA_STIRLING_C__ */
