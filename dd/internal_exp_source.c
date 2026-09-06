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
DD_FUNC(reduce_extred)(const gsl_dd a, const double m)
{
  const double log2_hi = 6.931471805599452862e-01;
  const double log2_lo = 2.319046813846299558e-17;
  const double log2_xtra = 5.7077084384162121e-34;
  gsl_dd p0, p1, p2;
  dd_exp3 x;

  p0 = DD_FUNC(two_prod)(log2_hi, m);
  p1 = DD_FUNC(two_prod)(log2_lo, m);
  p2 = DD_FUNC(two_prod)(log2_xtra, m);

  x.a = a.hi;
  x.b = a.lo;
  x.c = 0.0;

  exp3_add_double(&x, -p0.hi);
  exp3_add_double(&x, -p0.lo);
  exp3_add_double(&x, -p1.hi);
  exp3_add_double(&x, -p1.lo);
  exp3_add_double(&x, -p2.hi);
  exp3_add_double(&x, -p2.lo);
  exp3_renorm(&x);

  return exp3_to_dd(x);
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
 * backend-specific double-double multiplication selected by DD_FUNC().
 */
static gsl_dd
DD_FUNC(dd_exp)(const gsl_dd a)
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

  const double k = 1024.0;
  const double inv_k = 1.0 / k;
  const double log2_hi = 6.931471805599452862e-01;
  const double qd_eps = 4.93038065763132e-32;
  const double e_hi = 2.718281828459045091e+00;
  const double e_lo = 1.445646891729250158e-16;
  double m;
  gsl_dd r, s, t, p, one;
  int i;

  if (a.hi <= -709.0)
    {
      gsl_dd z;
      z.hi = 0.0;
      z.lo = 0.0;
      return z;
    }

  if (a.hi >= 709.0)
    {
      gsl_dd z;
      z.hi = HUGE_VAL;
      z.lo = HUGE_VAL;
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
      z.hi = e_hi;
      z.lo = e_lo;
      return z;
    }

  m = floor(a.hi / log2_hi + 0.5);

  /* QD reduction, with the product/subtraction residual retained. */
  r = DD_FUNC(reduce_extred)(a, m);
  r = mul_pwr2_extred(r, inv_k);

  p = DD_FUNC(dd_mul)(r, r);
  s = dd_add(r, mul_pwr2_extred(p, 0.5));

  p = DD_FUNC(dd_mul)(p, r);
  t.hi = inv_fact[0][0];
  t.lo = inv_fact[0][1];
  t = DD_FUNC(dd_mul)(p, t);

  i = 0;
  do
    {
      s = dd_add(s, t);
      p = DD_FUNC(dd_mul)(p, r);
      ++i;
      t.hi = inv_fact[i][0];
      t.lo = inv_fact[i][1];
      t = DD_FUNC(dd_mul)(p, t);
    }
  while (fabs(t.hi) > inv_k * qd_eps && i < 5);

  s = dd_add(s, t);

  /* Ten reconstruction steps: s <- 2*s + s*s. */
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_FUNC(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_FUNC(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_FUNC(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_FUNC(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_FUNC(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_FUNC(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_FUNC(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_FUNC(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_FUNC(dd_mul)(s, s));
  s = dd_add(mul_pwr2_extred(s, 2.0), DD_FUNC(dd_mul)(s, s));

  one.hi = 1.0;
  one.lo = 0.0;
  s = dd_add(s, one);

  return ldexp_extred(s, (int) m);
}
