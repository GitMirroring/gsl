/* specfunc/gamma/pow_exp.c
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
 * References:
 *
 *   W. J. Cody and W. Waite, Software Manual for the Elementary
 *   Functions, Prentice-Hall (1980) -- the argument reduction of exp();
 *
 *   T. J. Dekker, "A floating-point technique for extending the
 *   available precision", Numer. Math. 18 (1971) 224-242 -- the exact
 *   product (via gsl_dd_two_prod()) and the exact remainder of a
 *   rounded division.
 *
 * No external code was used.
 */

#ifndef __GSL_SF_GAMMA_POW_EXP_C__
#define __GSL_SF_GAMMA_POW_EXP_C__

/* requires <gsl/gsl_dd.h> and "constants.h" */

/*
 * Compute x - q*a exactly, where q = fl(x/a).  The remainder of a
 * rounded division is representable, and x - fl(q*a) is exact by
 * Sterbenz's lemma.
 */
static double
gamma_div_rem(const double x, const double a, const double q)
{
  const gsl_dd p = gsl_dd_two_prod(q, a);
  return (x - p.hi) - p.lo;
}

/*
 * Compute
 *
 *   (x/a)^a exp(d + s) = m (1 + dm) 2^k,
 *
 * for a > 0, x > 0, an exactly known exponent d and a small term s
 * (|s| < 1).  This is the core of the Stirling-type products
 * x^a exp(-x) / Gamma(a) and x^n / n!, in which the large exponents
 * a log(x/a) and d nearly cancel, so neither may be rounded:
 *
 *   - q = fl(x/a) and the exact remainder rho = x - q a give
 *     (x/a)^a = q^a exp(rho/q) to within rounding.  q^a is formed by
 *     pow() with a single rounding error; when it is outside the double
 *     range, for integer a as f^a 2^(e a) with q = f 2^e exactly,
 *     f in [1/sqrt(2), sqrt(2)), and otherwise (or if f^a is still out of
 *     range) as the scaled square of q^(a/2), at the cost of one rounding;
 *
 *   - exp(d) = 2^k exp(r), with the Cody-Waite reduction r = d - k log(2),
 *     |r| <= log(2)/2.
 *
 * The small terms are collected into one exponential, so m carries a few
 * rounding errors independent of the size of a, x and d.  The rounding
 * errors of the products forming m are returned as the relative
 * correction dm, for the caller to apply after its own (compensated)
 * product.  A result found to be below exp(GSL_SF_LOG_TRUE_MIN - 50)
 * before the full evaluation is returned as m = 0; otherwise an
 * underflowing result has k so negative that the caller's ldexp() gives 0.
 *
 * Returns 0 when a >= 1e20, |d| > 1e6, or no scaling of q^a is in range;
 * the caller must then use double-double.
 */
static int
gamma_pow_exp_scaled(const double a, const double x, const double d, const double s,
                     double * m_out, double * dm_out, int * k_out)
{
  const double q = x / a;
  double base = q, m, dm = 0.0, k, t;
  long eb = 0;             /* q^a = base^a 2^eb */
  int e = 0;
  gsl_dd p;

  *m_out = 0.0;
  *dm_out = 0.0;
  *k_out = 0;

  if (!(a < 1.0e20) || !(fabs(d) <= 1.0e6))
    {
      /* with q^a in range, the exponent is below 650 + d */
      return (d < -1.0e6 && a < 1.0e20 && a * log(q) + d < GSL_SF_LOG_TRUE_MIN - 50.0);
    }

  m = pow(base, a);

  if (!(m > 1.0e-300 && m < 1.0e280) && a == floor(a) && a <= 1.0e6)
    {
      base = frexp(q, &e);
      if (base < M_SQRT1_2)
        {
          base *= 2.0;
          --e;
        }

      eb = (long) e * (long) a;
      m = pow(base, a);
    }

  if (!(m > 1.0e-300 && m < 1.0e280))
    {
      /* base^a = 2^e (base^(a/2) 2^(-e/2))^2 */
      const double h = pow(base, 0.5 * a);

      if (h > 1.0e-280 && h < 1.0e280)
        {
          const double f = frexp(h, &e);
          const gsl_dd f2 = gsl_dd_two_prod(f, f);

          m = f2.hi;
          dm = f2.lo / m;
          eb += 2 * (long) e;
        }
      else
        {
          /* a log(q) + d in double is accurate to well within the margin
             for a < 1e20, |d| <= 1e6: skip the slow path on underflow */
          return (a * log(q) + d + s < GSL_SF_LOG_TRUE_MIN - 50.0);
        }
    }

  k = floor(d / M_LN2 + 0.5);    /* |k| < 2^21: k * GSL_SF_LN2_HI is exact */
  t = ((d - k * GSL_SF_LN2_HI) - k * GSL_SF_LN2_LO) + (gamma_div_rem(x, a, q) / q + s);
  p = gsl_dd_two_prod(m, exp(t));

  /* |eb| <= 1075e6 and |k| < 2^21 fit in an int */
  *m_out = p.hi;
  *dm_out = dm + p.lo / p.hi;
  *k_out = (int) (k + (double) eb);

  return 1;
}

#endif /* __GSL_SF_GAMMA_POW_EXP_C__ */
