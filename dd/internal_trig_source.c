/* dd/internal_trig_source.c
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
 * Compute sin(pi x) and/or cos(pi x) for a double-double x.
 *
 * With x reduced exactly to q/2 + r, |r| <= 1/4 (trig_pi_reduce()),
 * form u = pi r to triple precision, u + u3: pi = DD_PI_HI + DD_PI_LO +
 * DD_PI_XTRA and the products are TwoProd terms accumulated in a
 * four-component expansion, from which the rounding error u3 of the
 * double-double u is recovered.  On |u| <= pi/4,
 *
 *   sin(u) = u + u^3 S(u^2),   S(v) = sum_{k=1}^{14} (-1)^k v^(k-1) / (2k+1)!,
 *   cos(u) = 1 - u^2/2 + u^4 C(v),   C(v) = sum_{k=2}^{14} (-1)^k v^(k-2) / (2k)!,
 *
 * with truncation errors below 1e-35.  The leading terms u and u^2/2
 * are added exactly (u^2 split by TwoProd, including the u3
 * correction); the remaining terms have relative weight below 0.11 and
 * are evaluated by double-double Horner for the coefficients k <= 8 and
 * in double precision from v.hi for k >= 9 (relative weight below
 * 2e-18).  Since |v c_(k+1)| < |c_k| in every Horner step, the steps use
 * an addition without the cancellation handling of dd_add().  Each
 * result is rounded to double-double once, and the quadrant q mod 4
 * selects the signs.
 *
 * s or c may be NULL; only the series needed are evaluated.  For
 * infinite or NaN x both results are NaN.
 */
static void
DD_WRAP(dd_sincos_pi)(const gsl_dd x, gsl_dd * s, gsl_dd * c)
{
  static const gsl_dd sc[8] = {
    { -1.66666666666666657415e-01, -9.25185853854297065662e-18 },
    {  8.33333333333333321769e-03,  1.15648231731787138023e-19 },
    { -1.98412698412698412526e-04, -1.72095582934207052868e-22 },
    {  2.75573192239858925110e-06, -1.85839327404647208104e-22 },
    { -2.50521083854417202239e-08,  1.44881407093591196603e-24 },
    {  1.60590438368216133409e-10,  1.25852945887520980521e-26 },
    { -7.64716373181981640551e-13, -7.03872877733453001061e-30 },
    {  2.81145725434552059811e-15,  1.65088427308614325994e-31 } };
  static const gsl_dd cc[7] = {
    {  4.16666666666666643537e-02,  2.31296463463574266415e-18 },
    { -1.38888888888888894189e-03,  5.30054395437357705906e-20 },
    {  2.48015873015873015658e-05,  2.15119478667758816085e-23 },
    { -2.75573192239858882758e-07, -2.37677146222502973185e-23 },
    {  2.08767569878681001866e-09, -1.20734505911325997169e-25 },
    { -1.14707455977297245073e-11, -2.06555127528307454245e-28 },
    {  4.77947733238738525345e-14,  4.39920548583408125663e-31 } };
  gsl_dd r, u, v, p, q, a, su, cu;
  dd_ex4 sum;
  double u3, w, t;
  int quad, need_sin, need_cos, k;

  if (!(x.hi - x.hi == 0.0))
    {
      /* Inf or NaN */
      if (s != NULL)
        {
          s->hi = GSL_NAN;
          s->lo = 0.0;
        }
      if (c != NULL)
        {
          c->hi = GSL_NAN;
          c->lo = 0.0;
        }
      return;
    }

  quad = trig_pi_reduce(x, &r);
  need_sin = (s != NULL && !(quad & 1)) || (c != NULL && (quad & 1));
  need_cos = (s != NULL && (quad & 1)) || (c != NULL && !(quad & 1));

  /* u + u3 = pi r */
  ex4_zero(&sum);
  p = DD_WRAP(two_prod)(DD_PI_HI, r.hi);
  ex4_add_d(&sum, p.hi);
  ex4_add_d(&sum, p.lo);
  p = DD_WRAP(two_prod)(DD_PI_LO, r.hi);
  ex4_add_d(&sum, p.hi);
  ex4_add_d(&sum, p.lo);
  p = DD_WRAP(two_prod)(DD_PI_HI, r.lo);
  ex4_add_d(&sum, p.hi);
  ex4_add_d(&sum, p.lo + (DD_PI_LO * r.lo + DD_PI_XTRA * r.hi));
  u = ex4_to_dd(sum);
  ex4_add_d(&sum, -u.hi);
  ex4_add_d(&sum, -u.lo);
  u3 = ex4_to_dd(sum).hi;

  /* v = u^2: the exact parts p + 2 q and a double-double value */
  p = DD_WRAP(two_prod)(u.hi, u.hi);
  q = DD_WRAP(two_prod)(u.hi, u.lo);
  w = 0.5 * u.lo * u.lo + u.hi * u3;       /* (u^2 - p - 2q)/2 */
  v = fast_two_sum(p.hi, p.lo + 2.0 * q.hi);

  su.hi = su.lo = cu.hi = cu.lo = 0.0;

  if (need_sin)
    {
      t = -8.22063524662432949554e-18 + v.hi * (1.95729410633912625952e-20
          + v.hi * (-3.86817017063068412615e-23 + v.hi * (6.44695028438447358950e-26
          + v.hi * (-9.18368986379554600539e-29 + v.hi * 1.13099628864477158818e-31))));
      a.hi = t;
      a.lo = 0.0;
      for (k = 7; k >= 0; --k)
        {
          gsl_dd m = DD_WRAP(dd_mul)(v, a);
          t = sc[k].hi + m.hi;
          a = fast_two_sum(t, ((sc[k].hi - t) + m.hi) + (sc[k].lo + m.lo));
        }
      a = DD_WRAP(dd_mul)(DD_WRAP(dd_mul)(v, u), a);   /* u^3 S(v) */

      ex4_zero(&sum);
      ex4_add_d(&sum, u.hi);
      ex4_add_d(&sum, u.lo);
      ex4_add_d(&sum, u3);
      ex4_add_d(&sum, a.hi);
      ex4_add_d(&sum, a.lo);
      su = ex4_to_dd(sum);
    }

  if (need_cos)
    {
      gsl_dd q2;

      t = -1.56192069685862252711e-16 + v.hi * (4.11031762331216484407e-19
          + v.hi * (-8.89679139245057407789e-22 + v.hi * (1.61173757109611838590e-24
          + v.hi * (-2.47959626322479758961e-27 + v.hi * 3.27988923706983775898e-30))));
      q2.hi = t;
      q2.lo = 0.0;
      for (k = 6; k >= 0; --k)
        {
          gsl_dd m = DD_WRAP(dd_mul)(v, q2);
          t = cc[k].hi + m.hi;
          q2 = fast_two_sum(t, ((cc[k].hi - t) + m.hi) + (cc[k].lo + m.lo));
        }
      q2 = DD_WRAP(dd_mul)(DD_WRAP(dd_mul)(v, v), q2);  /* u^4 C(v) */

      ex4_zero(&sum);
      ex4_add_d(&sum, 1.0);
      ex4_add_d(&sum, -0.5 * p.hi);
      ex4_add_d(&sum, -0.5 * p.lo);
      ex4_add_d(&sum, -q.hi);
      ex4_add_d(&sum, -q.lo - w);
      ex4_add_d(&sum, q2.hi);
      ex4_add_d(&sum, q2.lo);
      cu = ex4_to_dd(sum);
    }

  switch (quad)
    {
    case 0:
      if (s != NULL) *s = su;
      if (c != NULL) *c = cu;
      break;
    case 1:
      if (s != NULL) *s = cu;
      if (c != NULL) { c->hi = -su.hi; c->lo = -su.lo; }
      break;
    case 2:
      if (s != NULL) { s->hi = -su.hi; s->lo = -su.lo; }
      if (c != NULL) { c->hi = -cu.hi; c->lo = -cu.lo; }
      break;
    default:
      if (s != NULL) { s->hi = -cu.hi; s->lo = -cu.lo; }
      if (c != NULL) *c = su;
      break;
    }
}
