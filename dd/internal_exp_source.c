/* dd/internal_exp_source.c
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

#define DD_EXP_REDUCTION_FACTOR      (128.0)
#define DD_EXP_ARG_MIN_GUARD         (-746.0)
#define DD_EXP_ARG_MAX_GUARD         (710.0)
#define DD_EXP_MIN_SCALE             (-1075)
#define DD_EXP_MAX_SCALE             (1024)
#define DD_HALF_ULP_ONE              (0.5 * GSL_DBL_EPSILON)
#define DD_QUARTER_ULP_ONE           (0.25 * GSL_DBL_EPSILON)

#define DD_EXPM1_TAYLOR_D10          0
#define DD_EXPM1_TAYLOR_D32          1

/*
 * Evaluate expm1(x) by Taylor expansion with compensated accumulation.
 *
 * The expansion
 *
 *   exp(x) - 1 = x + x^2/2! + x^3/3! + ...
 *
 * is evaluated by forming the powers in double-double arithmetic and
 * accumulating each complete double-double term into a three-component
 * floating-point expansion.  The expansion is compressed to double-double
 * only once, after all terms have been accumulated.
 *
 * Two evaluation modes are provided because this kernel is shared by exp()
 * and expm1():
 *
 *   flag == DD_EXPM1_TAYLOR_D10:
 *     Evaluate through degree 10.  This mode is used by dd_exp() after
 *     128-fold argument reduction.  The x^2/2 term is formed by exact
 *     power-of-two scaling
 *
 *   flag == DD_EXPM1_TAYLOR_D32:
 *     Evaluate through degree 32.  This mode is used directly by dd_expm1()
 *     on its central interval, avoiding both cancellation in exp(x)-1 and
 *     the rounding errors associated with repeated expm1 reconstruction.
 *
 * The coefficients 1/n!, n = 2,...,32, are stored as double-double
 * constants.
 */
static gsl_dd
DD_WRAP(kernel_expm1_taylor)(const gsl_dd x, const int flag)
{
  /* 1/n! table, n = 2,...,32. */
  static const gsl_dd inv_fact_table[31] = {
    { 5.00000000000000000000e-01,  0.00000000000000000000e+00 },
    { 1.66666666666666657415e-01,  9.25185853854297065662e-18 },
    { 4.16666666666666643537e-02,  2.31296463463574266415e-18 },
    { 8.33333333333333321769e-03,  1.15648231731787138000e-19 },
    { 1.38888888888888894189e-03, -5.30054395437357705906e-20 },
    { 1.98412698412698412526e-04,  1.72095582934207052868e-22 },
    { 2.48015873015873015658e-05,  2.15119478667758816085e-23 },
    { 2.75573192239858925110e-06, -1.85839327404647208104e-22 },
    { 2.75573192239858882758e-07,  2.37677146222502973185e-23 },
    { 2.50521083854417202239e-08, -1.44881407093591196603e-24 },
    { 2.08767569878681001866e-09, -1.20734505911325997169e-25 },
    { 1.60590438368216133409e-10,  1.25852945887520980521e-26 },
    { 1.14707455977297245073e-11,  2.06555127528307454245e-28 },
    { 7.64716373181981640551e-13,  7.03872877733453001061e-30 },
    { 4.77947733238738525345e-14,  4.39920548583408125663e-31 },
    { 2.81145725434552059811e-15,  1.65088427308614325994e-31 },
    { 1.56192069685862252711e-16,  1.19106796602737540024e-32 },
    { 8.22063524662432949554e-18,  2.21418941196042653637e-34 },
    { 4.11031762331216484407e-19,  1.44129733786595271498e-36 },
    { 1.95729410633912625952e-20, -1.36435038300879084872e-36 },
    { 8.89679139245057407789e-22, -7.91140261487237621703e-38 },
    { 3.86817017063068412615e-23, -8.84317765548234384789e-40 },
    { 1.61173757109611838590e-24, -3.68465735645097660329e-41 },
    { 6.44695028438447358950e-26, -1.93304042337034648242e-42 },
    { 2.47959626322479758961e-27, -1.29537309647652287758e-43 },
    { 9.18368986379554600539e-29,  1.43031503967873220242e-45 },
    { 3.27988923706983775898e-30,  1.51175427440298786897e-46 },
    { 1.13099628864477158818e-31,  1.04980154129595060239e-47 },
    { 3.76998762881590538515e-33,  2.58703478327503238348e-49 },
    { 1.21612504155351789377e-34,  5.58629056788880576964e-51 },
    { 3.80039075485474341802e-36,  1.74571580246525180301e-52 } };
  gsl_dd p, t;
  dd_ex3 s;
  int i;

  s.a = 0.0;
  s.b = 0.0;
  s.c = 0.0;

  ex3_add_d(&s, x.hi);
  ex3_add_d(&s, x.lo);

  p = DD_WRAP(dd_mul)(x, x);

  if (flag == DD_EXPM1_TAYLOR_D10)
    {
      t = mul_pwr2_extred(p, 0.5);
      ex3_add_d(&s, t.hi);
      ex3_add_d(&s, t.lo);

      for (i = 0; i < 8; ++i)
        {
          p = DD_WRAP(dd_mul)(p, x);
          t = DD_WRAP(dd_mul)(p, inv_fact_table[i + 1]);
          ex3_add_d(&s, t.hi);
          ex3_add_d(&s, t.lo);
        }
    }
  else
    {
      for (i = 0; i < 31; ++i)
        {
          t = DD_WRAP(dd_mul)(p, inv_fact_table[i]);
          ex3_add_d(&s, t.hi);
          ex3_add_d(&s, t.lo);

          if (i != 30)
            p = DD_WRAP(dd_mul)(p, x);
        }
    }

  return ex3_to_dd(s);
}

/*
 * Reduce a double-double argument by an integral multiple of log(2).
 *
 * Compute
 *
 *   r = a - m*log(2)
 *
 * using the three-component expansion
 *
 *   log(2) = DD_LN2_HI + DD_LN2_LO + DD_LN2_XTRA.
 *
 * Each product m*DD_LN2_* is formed with an error-free TwoProd operation.
 * Both components of all three products are then accumulated together with
 * the two components of a in an EX3 expansion.  The expansion is
 * renormalized and compressed only after the complete subtraction has been
 * accumulated.
 *
 * This delayed compression reduces the range-reduction error that would
 * result from performing the three subtractions successively in
 * double-double arithmetic.
 */
static gsl_dd
DD_WRAP(reduce_extred)(const gsl_dd a, const double m)
{
  gsl_dd p0, p1, p2;
  dd_ex3 x;

  p0 = DD_WRAP(two_prod)(DD_LN2_HI, m);
  p1 = DD_WRAP(two_prod)(DD_LN2_LO, m);
  p2 = DD_WRAP(two_prod)(DD_LN2_XTRA, m);

  x.a = a.hi;
  x.b = a.lo;
  x.c = 0.0;

  ex3_add_d(&x, -p0.hi);
  ex3_add_d(&x, -p0.lo);
  ex3_add_d(&x, -p1.hi);
  ex3_add_d(&x, -p1.lo);
  ex3_add_d(&x, -p2.hi);
  ex3_add_d(&x, -p2.lo);
  ex3_renorm(&x);

  return ex3_to_dd(x);
}

/*
 * Compute exp(a) in double-double arithmetic.
 *
 * The argument is first reduced as
 *
 *   a = m*log(2) + r,
 *
 * where m is an integer chosen so that |r| <= log(2)/2.  The remainder is
 * then scaled by 1/128:
 *
 *   y = r/128.
 *
 * On this small interval, expm1(y) is evaluated with the degree-10 mode of
 * kernel_expm1_taylor().  Seven applications of
 *
 *   s <- 2*s + s*s
 *
 * use the identity
 *
 *   expm1(2x) = 2*expm1(x) + expm1(x)^2
 *
 * to undo the factor-of-128 reduction.  Thus, after reconstruction,
 *
 *   s = expm1(r).
 *
 * The reduced exponential is then formed as 1+s and the original binary
 * scale is restored exactly by multiplication by 2^m.
 *
 * Special handling near the binary64 underflow and overflow boundaries
 * avoids premature overflow/underflow during the final scaling and returns
 * the appropriate double-double representation at the ends of the
 * representable range.
 *
 * The logarithmic range reduction uses the three-limb log(2) expansion and
 * EX3 accumulation in reduce_extred().
 */
static gsl_dd
DD_WRAP(dd_exp)(const gsl_dd a)
{
  const double k = DD_EXP_REDUCTION_FACTOR;
  const double inv_k = 1.0 / k;
  double m;
  gsl_dd r, s, one;
  int i;

  if (a.hi < DD_EXP_ARG_MIN_GUARD)
    {
      gsl_dd z;
      z.hi = 0.0;
      z.lo = 0.0;
      return z;
    }

  if (a.hi > DD_EXP_ARG_MAX_GUARD)
    {
      gsl_dd z;
      z.hi = GSL_POSINF;
      z.lo = GSL_POSINF;
      return z;
    }

  if (a.hi == 0.0 && a.lo == 0.0)
    {
      gsl_dd z;
      z.hi = 1.0;
      z.lo = 0.0;
      return z;
    }

  if (a.hi == 1.0 && a.lo == 0.0)
    {
      gsl_dd z;
      z.hi = DD_E_HI;
      z.lo = DD_E_LO;
      return z;
    }

  m = floor(a.hi / DD_LN2_HI + 0.5);

  r = DD_WRAP(reduce_extred)(a, m);
  r = mul_pwr2_extred(r, inv_k);
  s = DD_WRAP(kernel_expm1_taylor)(r, DD_EXPM1_TAYLOR_D10);

  for (i = 0; i < 7; ++i)
    s = dd_add(mul_pwr2_extred(s, 2.0), DD_WRAP(dd_mul)(s, s));

  one.hi = 1.0;
  one.lo = 0.0;
  s = dd_add(s, one);

  if (m > DD_EXP_MAX_SCALE)
    {
      gsl_dd z;
      z.hi = GSL_POSINF;
      z.lo = GSL_POSINF;
      return z;
    }

  if (m < DD_EXP_MIN_SCALE)
    {
      gsl_dd z;
      z.hi = 0.0;
      z.lo = 0.0;
      return z;
    }

  if (m == DD_EXP_MIN_SCALE)
    {
      gsl_dd z;
      if (s.hi > 1.0 || (s.hi == 1.0 && s.lo > 0.0))
        z.hi = ldexp(1.0, -1074);
      else
        z.hi = 0.0;
      z.lo = 0.0;
      return z;
    }

  if (m == DD_EXP_MAX_SCALE)
    {
      if (s.hi < 1.0)
        return dd_ldexp(s, 1024);
      else if (s.hi == 1.0 && s.lo <= -DD_QUARTER_ULP_ONE)
        {
          gsl_dd z;
          double d;
          z.hi = GSL_DBL_MAX;
          d = (s.hi - 1.0) + DD_HALF_ULP_ONE + s.lo;
          z.lo = ldexp(d, 1024);
          return z;
        }
      else
        {
          gsl_dd z;
          z.hi = GSL_POSINF;
          z.lo = GSL_POSINF;
          return z;
        }
    }

  return dd_ldexp(s, (int) m);
}

#define DD_EXPM1_XMIN (-5.75000000000000000e-01)
#define DD_EXPM1_XMAX ( 1.13600000000000000e+00)

/*
 * Compute expm1(x) = exp(x) - 1 in double-double arithmetic.
 *
 * Direct subtraction of one from dd_exp(x) is poorly conditioned near zero,
 * so a dedicated direct Taylor evaluation is used on the asymmetric central
 * interval
 *
 *   -0.575 <= x <= 1.136.
 *
 * There,
 *
 *   expm1(x) = x + x^2/2! + ... + x^32/32!
 *
 * is evaluated by the degree-32 mode of kernel_expm1_taylor().  Complete
 * double-double terms are accumulated in EX3 and compressed only once.
 * This avoids both exp(x)-1 cancellation and the repeated reconstruction
 * roundings of the logarithmically reduced exp algorithm.
 *
 * Outside the central interval, subtraction is sufficiently well
 * conditioned and the result is formed as
 *
 *   dd_exp(x) - 1.
 *
 * The central interval is intentionally asymmetric because cancellation in
 * exp(x)-1 persists to larger positive arguments than negative arguments.
 *
 * Comparisons at the two transition points use both components of the
 * normalized double-double input.  Thus a nonzero low component at an exact
 * high-word threshold selects the branch according to the actual represented
 * double-double value rather than x.hi alone.
 */
static gsl_dd
DD_WRAP(dd_expm1)(const gsl_dd x)
{
  gsl_dd one;
  int central;

  if (x.hi == 0.0 && x.lo == 0.0)
    return x;

  central = 1;

  if (x.hi < DD_EXPM1_XMIN)
    central = 0;
  else if (x.hi == DD_EXPM1_XMIN && x.lo < 0.0)
    central = 0;

  if (x.hi > DD_EXPM1_XMAX)
    central = 0;
  else if (x.hi == DD_EXPM1_XMAX && x.lo > 0.0)
    central = 0;

  if (central)
    return DD_WRAP(kernel_expm1_taylor)(x, DD_EXPM1_TAYLOR_D32);

  one.hi = 1.0;
  one.lo = 0.0;

  return dd_sub(DD_WRAP(dd_exp)(x), one);
}
