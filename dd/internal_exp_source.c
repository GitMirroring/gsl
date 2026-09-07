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

#define DD_EXP_TAYLOR_EPS        4.9303806576313237838233035330174139354575e-32 /* threshold used by iterative series, 2^{-104} */
#define DD_EXP_REDUCTION_FACTOR  1024.0
#define DD_HALF_ULP_ONE          (0.5 * GSL_DBL_EPSILON)                        /* 2^{-53} */
#define DD_QUARTER_ULP_ONE       (0.25 * GSL_DBL_EPSILON)                       /* 2^{-54} */
#define DD_EXP_ARG_MIN_GUARD     (-746.0)
#define DD_EXP_ARG_MAX_GUARD     (710.0)
#define DD_EXP_MIN_SCALE         (-1075)
#define DD_EXP_MAX_SCALE         (1024)

/*
 * Reduce a double-double argument by an integral multiple of log(2).
 *
 * Compute
 *
 *   a - m*log(2)
 *
 * using a three-limb expansion of log(2).  Each product m*log2_i is
 * formed with the backend-specific error-free TwoProd operation, and
 * both the rounded product and its residual are accumulated in a
 * temporary three-term expansion before compression to double-double.
 *
 * The three log(2) limbs provide substantially more range-reduction
 * accuracy than the original two-limb QD constant when |m| is large.
 * The temporary three-term expansion is used only for range reduction;
 * all subsequent exp() arithmetic is ordinary double-double.
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
 * This routine is based on the dd_real exponential algorithm from the
 * QD library, with enhanced range reduction for improved double-double
 * accuracy.
 *
 * First choose an integer m approximately equal to a/log(2), and form
 *
 *   r = (a - m*log(2)) / 1024.
 *
 * The subtraction a - m*log(2) is evaluated by
 * reduce_extred(), which uses three binary64 limbs of log(2)
 * and retains the TwoProd residuals during the reduction.
 *
 * For the small reduced argument r, evaluate exp(r)-1 using
 *
 *   r + r^2/2! + r^3/3! + ...
 *
 * in double-double arithmetic with the QD double-double inverse
 * factorial constants.
 *
 * Recover exp(a - m*log(2)) by applying ten times
 *
 *   s <- 2*s + s^2,
 *
 * since 1 + (2*s + s^2) = (1+s)^2.  The factor 1024 = 2^10 is
 * thereby restored.  Finally add one and scale by 2^m.
 *
 * Relative to the original QD algorithm, this implementation uses a
 * reduction factor of 1024 rather than 512 and a three-limb log(2)
 * range reduction.  The Taylor evaluation and reconstruction use the
 * backend-specific double-double multiplication selected by DD_WRAP().
 */
static gsl_dd
DD_WRAP(dd_exp)(const gsl_dd a)
{
  static const double inv_fact[15][2] = {
    { 1.66666666666666657e-01,  9.25185853854297066e-18 },
    { 4.16666666666666644e-02,  2.31296463463574266e-18 },
    { 8.33333333333333322e-03,  1.15648231731787138e-19 },
    { 1.38888888888888894e-03, -5.30054395437357706e-20 },
    { 1.98412698412698413e-04,  1.72095582934207053e-22 },
    { 2.48015873015873016e-05,  2.15119478667758816e-23 },
    { 2.75573192239858925e-06, -1.85839327404647208e-22 },
    { 2.75573192239858883e-07,  2.37677146222502973e-23 },
    { 2.50521083854417202e-08, -1.44881407093591197e-24 },
    { 2.08767569878681002e-09, -1.20734505911325997e-25 },
    { 1.60590438368216133e-10,  1.25852945887520981e-26 },
    { 1.14707455977297245e-11,  2.06555127528307454e-28 },
    { 7.64716373181981641e-13,  7.03872877733453001e-30 },
    { 4.77947733238738525e-14,  4.39920548583408126e-31 },
    { 2.81145725434552060e-15,  1.65088427308614326e-31 }
  };

  const double k = DD_EXP_REDUCTION_FACTOR;
  const double inv_k = 1.0 / k;
  double m;
  gsl_dd r, s, t, p, one;
  int i;

  /* Protect the integer range reduction for arguments far outside the
     representable exponential range.  The actual lower/upper boundaries
     are determined below from the reduced exponent m. */
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

  /* QD reduction, with the product/subtraction residual retained. */
  r = DD_WRAP(reduce_extred)(a, m);
  r = mul_pwr2_extred(r, inv_k);

  p = DD_WRAP(dd_mul)(r, r);
  s = dd_add(r, mul_pwr2_extred(p, 0.5));

  p = DD_WRAP(dd_mul)(p, r);
  t.hi = inv_fact[0][0];
  t.lo = inv_fact[0][1];
  t = DD_WRAP(dd_mul)(p, t);

  i = 0;
  do
    {
      s = dd_add(s, t);
      p = DD_WRAP(dd_mul)(p, r);
      ++i;
      t.hi = inv_fact[i][0];
      t.lo = inv_fact[i][1];
      t = DD_WRAP(dd_mul)(p, t);
    }
  while (fabs(t.hi) > inv_k * DD_EXP_TAYLOR_EPS && i < 5);

  s = dd_add(s, t);

  /* Ten reconstruction steps: s <- 2*s + s*s. */
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_WRAP(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_WRAP(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_WRAP(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_WRAP(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_WRAP(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_WRAP(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_WRAP(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_WRAP(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_WRAP(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_WRAP(dd_mul)(s, s));

  one.hi = 1.0;
  one.lo = 0.0;
  s = dd_add(s, one);

  /* The reduced exponent identifies the representable DD output range. */
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

      /* The only nonzero normalized DD value in this range is the
         smallest positive binary64 number.  Since
         exp(a) = s * 2^-1075, round to 0 or 2^-1074 according to
         whether the reduced significand is below or above one. */
      if (s.hi > 1.0 || (s.hi == 1.0 && s.lo > 0.0))
        z.hi = ldexp(1.0, -1074); /* smallest positive binary64 subnormal, 2^{-1074} */
      else
        z.hi = 0.0;
      z.lo = 0.0;
      return z;
    }

  if (m == DD_EXP_MAX_SCALE)
    {
      /* A finite result with exponent 1024 is possible only in the
         narrow interval immediately above GSL_DBL_MAX.  Form it with
         GSL_DBL_MAX as the leading word so that no intermediate binary64
         addition overflows. */
      if (s.hi < 1.0)
        return ldexp_extred(s, 1024);
      else if (s.hi == 1.0 && s.lo <= -DD_QUARTER_ULP_ONE)
        {
          gsl_dd z;
          double d;

          /*
           * z.hi + z.lo = 2^1024 (s.hi + s.lo)
           *             = GSL_DBL_MAX + 2^1024 (2^-53 + s.lo).
           */
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

  return ldexp_extred(s, (int) m);
}
