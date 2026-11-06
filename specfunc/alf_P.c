/* specfunc/alf_P.c
 * 
 * Copyright (C) 2009-2026 Patrick Alken
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
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_alf.h>

/*
 * The routines in this module compute associated Legendre functions
 * (ALFs) up to order and degree 2700, using the method described
 * in
 *
 * [1] S. A. Holmes and W. E. Featherstone, A unified approach
 *     to the Clenshaw summation and the recursive computation of very
 *     high degree and order normalised associated Legendre functions,
 *     Journal of Geodesy, 76, pg. 279-299, 2002.
 *
 * Computations of the alternative derivatives, d^k/dtheta^k Plm
 * use the methods described in,
 *
 * [2] W. Bosch, On the computation of derivatives of Legendre functions,
 *     Phys. Chem. Earth (A), 25, 9-11, pg. 655-659, 2000.
 *
 * Further information on ALFs can be found in
 *
 * [3] Abramowitz and Stegun, Handbook of Mathematical Functions,
 *     Chapter 8, 1972.
 */

/*
 * Scaling of the sectoral seeds
 *
 * The seeds P(m,m) contain the factor sin^m(theta), which underflows the
 * double range for large m, although P(l,m) is still significant at
 * higher degrees l. Each seed (and its derivative seeds) is therefore
 * carried as a mantissa and a power-of-two exponent: when |P(m,m)| falls
 * below 2^-960, it is renormalized with frexp(). (One seed update
 * multiplies by sin(theta) >= 2^-27, so seeds above 2^-960 stay well
 * inside the normal double range.) Columns whose seed never
 * needed renormalizing are computed exactly as before. The others are
 * computed by the alf_*_column_scaled() functions in two phases:
 *
 * 1) while the true values are below 2^-1000, the recursion runs on the
 *    scaled values P(l,m) 2^-e; each stored value is multiplied by the
 *    exact factor 2^e, or set to 0 if it is below the smallest normal
 *    double 2^-1022, and the running values are rescaled by 2^-600 when
 *    they grow beyond 2^600;
 *
 * 2) once they are representable, the running values are converted to
 *    true values and the column continues with the ordinary recursion.
 *
 * All scale factors are powers of two, so the scaling introduces no
 * rounding: every value equals what double arithmetic with an unbounded
 * exponent range would give, and values below the normal double range
 * (2.2e-308) are returned as 0. Columns whose seed is still a normal
 * double go straight to phase 2, giving bitwise the same results as the
 * unscaled recursion. This works for every normalization.
 */

#define ALF_SMALL      1.0261342003245941e-289  /* 2^-960: renormalize seeds below this */
#define ALF_BIG        4.1495155688809930e+180  /* 2^600: rescale running values above this */
#define ALF_RESCALE    2.4099198651028841e-181  /* 2^-600 */
#define ALF_RESCALE_E  600

/*
 * Header of the precomputed region
 *
 * gsl_sf_alf_precompute() stores a small header right after the nlm entries
 * reserved for the ALF values, followed by the recurrence factors. The array
 * functions check it, so that an array initialized for a different lmax,
 * mmax (or not at all) is reported instead of silently giving wrong values.
 */

#define ALF_HDR_SIZE   8                /* header length (doubles) */
#define ALF_HDR_MAGIC  0                /* check value ALF_MAGIC */
#define ALF_HDR_LMAX   1                /* lmax */
#define ALF_HDR_MMAX   2                /* mmax */
#define ALF_HDR_NORM   3                /* normalization (gsl_sf_alf_t) */
#define ALF_HDR_FLAGS  4                /* GSL_SF_ALF_FLG_xxx */
#define ALF_HDR_CSFAC  5                /* Condon-Shortley factor, -1 or 1 */
#define ALF_HDR_P00    6                /* P(0,0) */
#define ALF_MAGIC      1196640321.0     /* "GSLA" */

/*
alf_header()
  Header of an array initialized by gsl_sf_alf_precompute() for lmax, mmax

Inputs: array - array
        nlm   - gsl_sf_alf_nlm(lmax, mmax)
        lmax  - maximum degree
        mmax  - maximum order

Return: pointer to the header, or NULL if the array was not initialized for
these lmax, mmax
*/

static inline const double *
alf_header(const double * array, const size_t nlm, const size_t lmax, const size_t mmax)
{
  const double * hdr = array + nlm;

  if (hdr[ALF_HDR_MAGIC] == ALF_MAGIC && hdr[ALF_HDR_LMAX] == (double) lmax &&
      hdr[ALF_HDR_MMAX] == (double) mmax)
    return hdr;

  return NULL;
}

/*
alf_scale_state()
  Factors and thresholds for scaled values v with true value v 2^e:

f1, f2 - factors with f1 * f2 = 2^e, applied as (v * f1) * f2: neither factor
         underflows or overflows prematurely, so the product is exact for
         normal results
tz     - |v| < tz: the true value is below the smallest normal double 2^-1022
         and is returned as 0 (subnormal results would be slow to compute
         on most processors, and are negligible)
tn     - |v| >= tn: the true value is at least 2^-1000 (a normal double with
         margin), so the recursion can continue on true values
*/

static void
alf_scale_state(const int e, double * f1, double * f2, double * tz, double * tn)
{
  if (e < -1000)
    {
      *f1 = ldexp(1.0, -1000);
      *f2 = ldexp(1.0, e + 1000);
    }
  else if (e > 1000)
    {
      *f1 = ldexp(1.0, 1000);
      *f2 = ldexp(1.0, e - 1000);
    }
  else
    {
      *f1 = ldexp(1.0, e);
      *f2 = 1.0;
    }

  *tz = ldexp(1.0, -1022 - e);
  *tn = ldexp(1.0, -1000 - e);
}

/* true value of a scaled value v (see alf_scale_state) */
#define ALF_TRUE(v, f1, f2, tz) ((fabs(v) < (tz)) ? 0.0 : ((v) * (f1)) * (f2))

/*
alf_P_column()
  Column m >= 1 of P(l,m), l = m..lmax, from the seed P(m,m) = pmm 2^e

Inputs: lmax - maximum degree
        m    - order
        x    - cos(theta)
        pmm  - seed P(m,m) 2^-e
        e    - seed exponent (0 for seeds in the double range)
        cm   - recurrence factor c_m
        al   - recurrence factors a(l,m), b(l,m) for l = m+2..lmax
        ptr  - (output) ptr[l] = P(l,m)

Notes:
1) For e != 0 the column starts in scaled form (phase 1, see the notes at
the top of this file); once the values are representable, and from the
start for e = 0, the ordinary recursion runs on true values (phase 2)
*/

static inline void
alf_P_column(const size_t lmax, const size_t m, const double x, const double pmm, int e,
             const double cm, const double * al, double * ptr)
{
  double prev, cur;
  size_t l;

  if (e == 0)
    {
      ptr[m] = pmm;
      if (m == lmax)
        return;

      prev = pmm;
      cur = (cm * x) * pmm;
      ptr[m + 1] = cur;
      l = m + 1;
    }
  else
    {
      double f1, f2, tz, tn;

      alf_scale_state(e, &f1, &f2, &tz, &tn);

      ptr[m] = ALF_TRUE(pmm, f1, f2, tz);
      if (m == lmax)
        return;

      prev = pmm;
      cur = (cm * x) * pmm;
      l = m + 1;

      /* phase 1: true values below 2^-1000, two degrees per check */
      while (fabs(prev) < tn || fabs(cur) < tn)
        {
          ptr[l] = ALF_TRUE(cur, f1, f2, tz);
          if (l + 2 > lmax)
            {
              if (l < lmax)
                {
                  const double next = (al[0] * x) * cur + al[1] * prev;
                  ptr[l + 1] = ALF_TRUE(next, f1, f2, tz);
                }
              return;
            }

          prev = (al[0] * x) * cur + al[1] * prev;
          cur = (al[2] * x) * prev + al[3] * cur;
          ptr[l + 1] = ALF_TRUE(prev, f1, f2, tz);
          al += 4;
          l += 2;

          if (fabs(cur) > ALF_BIG)
            {
              prev *= ALF_RESCALE;
              cur *= ALF_RESCALE;
              e += ALF_RESCALE_E;
              alf_scale_state(e, &f1, &f2, &tz, &tn);
            }
        }

      /* exact conversion to true values (the results are normal) */
      prev = (prev * f1) * f2;
      cur = (cur * f1) * f2;
      ptr[l] = cur;
    }

  /* phase 2: degrees l+1..lmax, two per step */
  for (++l; l < lmax; l += 2)
    {
      prev = (al[0] * x) * cur + al[1] * prev;
      cur = (al[2] * x) * prev + al[3] * cur;
      ptr[l] = prev;
      ptr[l + 1] = cur;
      al += 4;
    }

  if (l == lmax)
    ptr[l] = (al[0] * x) * cur + al[1] * prev;
}

/*
alf_dx_column()
  Column m >= 1 of P(l,m), dP(l,m)/dx and (if d2ptr != NULL) d^2P(l,m)/dx^2,
l = m..lmax, from the seeds (P(m,m), P'(m,m), P''(m,m)) 2^e

Inputs: lmax  - maximum degree
        m     - order
        x     - cos(theta)
        pmm   - seed P(m,m) 2^-e
        dpmm  - seed P'(m,m) 2^-e
        d2pmm - seed P''(m,m) 2^-e (unused if d2ptr == NULL)
        e     - seed exponent (0 for seeds in the double range)
        cm    - recurrence factor c_m
        al    - recurrence factors a(l,m), b(l,m) for l = m+2..lmax
        ptr   - (output) ptr[l] = P(l,m)
        dptr  - (output) dptr[l] = dP(l,m)/dx
        d2ptr - (output) d2ptr[l] = d^2P(l,m)/dx^2, or NULL

Notes:
1) Phases as in alf_P_column(); the switch is decided by P
*/

static inline void
alf_dx_column(const size_t lmax, const size_t m, const double x, const double pmm,
              const double dpmm, const double d2pmm, int e, const double cm,
              const double * al, double * ptr, double * dptr, double * d2ptr)
{
  double prev = pmm, dprev = dpmm, d2prev = d2pmm;
  double cur, dcur, d2cur;
  double f1 = 1.0, f2 = 1.0, tz = 0.0, tn = 0.0;
  size_t l;

  if (e != 0)
    alf_scale_state(e, &f1, &f2, &tz, &tn);

  ptr[m] = (e == 0) ? pmm : ALF_TRUE(pmm, f1, f2, tz);
  dptr[m] = (e == 0) ? dpmm : ALF_TRUE(dpmm, f1, f2, tz);
  if (d2ptr)
    d2ptr[m] = (e == 0) ? d2pmm : ALF_TRUE(d2pmm, f1, f2, tz);

  if (m == lmax)
    return;

  cur = (cm * x) * prev;
  dcur = cm * (prev + x * dprev);
  d2cur = cm * (2.0 * dprev + x * d2prev);
  l = m + 1;

  if (e == 0)
    {
      ptr[l] = cur;
      dptr[l] = dcur;
      if (d2ptr)
        d2ptr[l] = d2cur;
    }
  else
    {
      /* phase 1: true values below 2^-1000, two degrees per check */
      while (fabs(prev) < tn || fabs(cur) < tn)
        {
          ptr[l] = ALF_TRUE(cur, f1, f2, tz);
          dptr[l] = ALF_TRUE(dcur, f1, f2, tz);
          if (d2ptr)
            d2ptr[l] = ALF_TRUE(d2cur, f1, f2, tz);

          if (l + 2 > lmax)
            {
              if (l < lmax)
                {
                  const double next = (al[0] * x) * cur + al[1] * prev;
                  const double dnext = al[0] * (x * dcur + cur) + al[1] * dprev;
                  const double d2next = al[0] * (x * d2cur + 2.0 * dcur) + al[1] * d2prev;

                  ptr[l + 1] = ALF_TRUE(next, f1, f2, tz);
                  dptr[l + 1] = ALF_TRUE(dnext, f1, f2, tz);
                  if (d2ptr)
                    d2ptr[l + 1] = ALF_TRUE(d2next, f1, f2, tz);
                }
              return;
            }

          prev = (al[0] * x) * cur + al[1] * prev;
          dprev = al[0] * (x * dcur + cur) + al[1] * dprev;
          d2prev = al[0] * (x * d2cur + 2.0 * dcur) + al[1] * d2prev;
          cur = (al[2] * x) * prev + al[3] * cur;
          dcur = al[2] * (x * dprev + prev) + al[3] * dcur;
          d2cur = al[2] * (x * d2prev + 2.0 * dprev) + al[3] * d2cur;
          ptr[l + 1] = ALF_TRUE(prev, f1, f2, tz);
          dptr[l + 1] = ALF_TRUE(dprev, f1, f2, tz);
          if (d2ptr)
            d2ptr[l + 1] = ALF_TRUE(d2prev, f1, f2, tz);
          al += 4;
          l += 2;

          if (fabs(cur) > ALF_BIG)
            {
              prev *= ALF_RESCALE;
              cur *= ALF_RESCALE;
              dprev *= ALF_RESCALE;
              dcur *= ALF_RESCALE;
              d2prev *= ALF_RESCALE;
              d2cur *= ALF_RESCALE;
              e += ALF_RESCALE_E;
              alf_scale_state(e, &f1, &f2, &tz, &tn);
            }
        }

      /* exact conversion to true values */
      prev = (prev * f1) * f2;
      cur = (cur * f1) * f2;
      dprev = (dprev * f1) * f2;
      dcur = (dcur * f1) * f2;
      d2prev = (d2prev * f1) * f2;
      d2cur = (d2cur * f1) * f2;
      ptr[l] = cur;
      dptr[l] = dcur;
      if (d2ptr)
        d2ptr[l] = d2cur;
    }

  /* phase 2: degrees l+1..lmax, two per step */
  if (d2ptr)
    {
      for (++l; l < lmax; l += 2)
        {
          prev = (al[0] * x) * cur + al[1] * prev;
          dprev = al[0] * (x * dcur + cur) + al[1] * dprev;
          d2prev = al[0] * (x * d2cur + 2.0 * dcur) + al[1] * d2prev;
          cur = (al[2] * x) * prev + al[3] * cur;
          dcur = al[2] * (x * dprev + prev) + al[3] * dcur;
          d2cur = al[2] * (x * d2prev + 2.0 * dprev) + al[3] * d2cur;
          ptr[l] = prev; ptr[l + 1] = cur;
          dptr[l] = dprev; dptr[l + 1] = dcur;
          d2ptr[l] = d2prev; d2ptr[l + 1] = d2cur;
          al += 4;
        }

      if (l == lmax)
        {
          ptr[l] = (al[0] * x) * cur + al[1] * prev;
          dptr[l] = al[0] * (x * dcur + cur) + al[1] * dprev;
          d2ptr[l] = al[0] * (x * d2cur + 2.0 * dcur) + al[1] * d2prev;
        }
    }
  else
    {
      for (++l; l < lmax; l += 2)
        {
          prev = (al[0] * x) * cur + al[1] * prev;
          dprev = al[0] * (x * dcur + cur) + al[1] * dprev;
          cur = (al[2] * x) * prev + al[3] * cur;
          dcur = al[2] * (x * dprev + prev) + al[3] * dcur;
          ptr[l] = prev; ptr[l + 1] = cur;
          dptr[l] = dprev; dptr[l + 1] = dcur;
          al += 4;
        }

      if (l == lmax)
        {
          ptr[l] = (al[0] * x) * cur + al[1] * prev;
          dptr[l] = al[0] * (x * dcur + cur) + al[1] * dprev;
        }
    }
}

/*
alf_vsh_column()
  Column m >= 1 of P(l,m)/sin(theta) and dP(l,m)/dtheta, l = m..lmax, from
the seeds (P(m,m)/sin(theta), dP(m,m)/dtheta) 2^e

Inputs: lmax - maximum degree
        m    - order
        x    - cos(theta)
        usq  - sin^2(theta)
        pmm  - seed P(m,m)/sin(theta) 2^-e
        dpmm - seed dP(m,m)/dtheta 2^-e
        e    - seed exponent (0 for seeds in the double range)
        cm   - recurrence factor c_m
        al   - recurrence factors a(l,m), b(l,m) for l = m+2..lmax
        ptr  - (output) ptr[l] = P(l,m)/sin(theta)
        dptr - (output) dptr[l] = dP(l,m)/dtheta

Notes:
1) Phases as in alf_P_column(); the switch is decided by P/sin(theta)
*/

static inline void
alf_vsh_column(const size_t lmax, const size_t m, const double x, const double usq,
               const double pmm, const double dpmm, int e, const double cm,
               const double * al, double * ptr, double * dptr)
{
  double prev = pmm, dprev = dpmm, cur, dcur;
  size_t l;

  if (e == 0)
    {
      ptr[m] = pmm;
      dptr[m] = dpmm;
      if (m == lmax)
        return;

      cur = (cm * x) * prev;
      dcur = cm * (-usq * prev + x * dprev);
      l = m + 1;
      ptr[l] = cur;
      dptr[l] = dcur;
    }
  else
    {
      double f1, f2, tz, tn;

      alf_scale_state(e, &f1, &f2, &tz, &tn);

      ptr[m] = ALF_TRUE(pmm, f1, f2, tz);
      dptr[m] = ALF_TRUE(dpmm, f1, f2, tz);
      if (m == lmax)
        return;

      cur = (cm * x) * prev;
      dcur = cm * (-usq * prev + x * dprev);
      l = m + 1;

      /* phase 1: true values below 2^-1000, two degrees per check */
      while (fabs(prev) < tn || fabs(cur) < tn)
        {
          ptr[l] = ALF_TRUE(cur, f1, f2, tz);
          dptr[l] = ALF_TRUE(dcur, f1, f2, tz);

          if (l + 2 > lmax)
            {
              if (l < lmax)
                {
                  const double next = (al[0] * x) * cur + al[1] * prev;
                  const double dnext = al[0] * (x * dcur - usq * cur) + al[1] * dprev;

                  ptr[l + 1] = ALF_TRUE(next, f1, f2, tz);
                  dptr[l + 1] = ALF_TRUE(dnext, f1, f2, tz);
                }
              return;
            }

          prev = (al[0] * x) * cur + al[1] * prev;
          dprev = al[0] * (x * dcur - usq * cur) + al[1] * dprev;
          cur = (al[2] * x) * prev + al[3] * cur;
          dcur = al[2] * (x * dprev - usq * prev) + al[3] * dcur;
          ptr[l + 1] = ALF_TRUE(prev, f1, f2, tz);
          dptr[l + 1] = ALF_TRUE(dprev, f1, f2, tz);
          al += 4;
          l += 2;

          if (fabs(cur) > ALF_BIG)
            {
              prev *= ALF_RESCALE;
              cur *= ALF_RESCALE;
              dprev *= ALF_RESCALE;
              dcur *= ALF_RESCALE;
              e += ALF_RESCALE_E;
              alf_scale_state(e, &f1, &f2, &tz, &tn);
            }
        }

      /* exact conversion to true values */
      prev = (prev * f1) * f2;
      cur = (cur * f1) * f2;
      dprev = (dprev * f1) * f2;
      dcur = (dcur * f1) * f2;
      ptr[l] = cur;
      dptr[l] = dcur;
    }

  /* phase 2: degrees l+1..lmax, two per step */
  for (++l; l < lmax; l += 2)
    {
      prev = (al[0] * x) * cur + al[1] * prev;
      dprev = al[0] * (x * dcur - usq * cur) + al[1] * dprev;
      cur = (al[2] * x) * prev + al[3] * cur;
      dcur = al[2] * (x * dprev - usq * prev) + al[3] * dcur;
      ptr[l] = prev; ptr[l + 1] = cur;
      dptr[l] = dprev; dptr[l + 1] = dcur;
      al += 4;
    }

  if (l == lmax)
    {
      ptr[l] = (al[0] * x) * cur + al[1] * prev;
      dptr[l] = al[0] * (x * dcur - usq * cur) + al[1] * dprev;
    }
}

/*
alf_vsh2_column()
  Column m >= 1 of P(l,m)/sin(theta), dP(l,m)/dtheta and d^2P(l,m)/dtheta^2,
l = m..lmax, from the seeds (P(m,m)/sin(theta), dP(m,m)/dtheta,
d^2P(m,m)/dtheta^2) 2^e

Inputs: lmax  - maximum degree
        m     - order
        x     - cos(theta)
        u     - sin(theta)
        pmm   - seed P(m,m)/sin(theta) 2^-e
        dpmm  - seed dP(m,m)/dtheta 2^-e
        d2pmm - seed d^2P(m,m)/dtheta^2 2^-e
        e     - seed exponent (0 for seeds in the double range)
        cm    - recurrence factor c_m
        al    - recurrence factors a(l,m), b(l,m) for l = m+2..lmax
        ptr   - (output) ptr[l] = P(l,m)/sin(theta)
        dptr  - (output) dptr[l] = dP(l,m)/dtheta
        d2ptr - (output) d2ptr[l] = d^2P(l,m)/dtheta^2

Notes:
1) With q = P/sin(theta), the recurrence P(l) = a x P(l-1) + b P(l-2) and its
theta derivatives (dx/dtheta = -u, d^2x/dtheta^2 = -x) give

q(l)   = a x q(l-1) + b q(l-2)
P'(l)  = a (x P'(l-1) - u^2 q(l-1)) + b P'(l-2)
P''(l) = a (x P''(l-1) - 2u P'(l-1) - x u q(l-1)) + b P''(l-2)

2) Phases as in alf_P_column(); the switch is decided by P/sin(theta)
*/

static inline void
alf_vsh2_column(const size_t lmax, const size_t m, const double x, const double u,
                  const double pmm, const double dpmm, const double d2pmm, int e,
                  const double cm, const double * al, double * ptr, double * dptr, double * d2ptr)
{
  const double usq = u * u;
  const double xu = x * u;
  const double u2 = 2.0 * u;
  double prev = pmm, dprev = dpmm, d2prev = d2pmm, cur, dcur, d2cur;
  size_t l;

  if (e == 0)
    {
      ptr[m] = pmm;
      dptr[m] = dpmm;
      d2ptr[m] = d2pmm;
      if (m == lmax)
        return;

      cur = (cm * x) * prev;
      dcur = cm * (x * dprev - usq * prev);
      d2cur = cm * (x * d2prev - u2 * dprev - xu * prev);
      l = m + 1;
      ptr[l] = cur;
      dptr[l] = dcur;
      d2ptr[l] = d2cur;
    }
  else
    {
      double f1, f2, tz, tn;

      alf_scale_state(e, &f1, &f2, &tz, &tn);

      ptr[m] = ALF_TRUE(pmm, f1, f2, tz);
      dptr[m] = ALF_TRUE(dpmm, f1, f2, tz);
      d2ptr[m] = ALF_TRUE(d2pmm, f1, f2, tz);
      if (m == lmax)
        return;

      cur = (cm * x) * prev;
      dcur = cm * (x * dprev - usq * prev);
      d2cur = cm * (x * d2prev - u2 * dprev - xu * prev);
      l = m + 1;

      /* phase 1: true values below 2^-1000, two degrees per check */
      while (fabs(prev) < tn || fabs(cur) < tn)
        {
          ptr[l] = ALF_TRUE(cur, f1, f2, tz);
          dptr[l] = ALF_TRUE(dcur, f1, f2, tz);
          d2ptr[l] = ALF_TRUE(d2cur, f1, f2, tz);

          if (l + 2 > lmax)
            {
              if (l < lmax)
                {
                  const double next = (al[0] * x) * cur + al[1] * prev;
                  const double dnext = al[0] * (x * dcur - usq * cur) + al[1] * dprev;
                  const double d2next = al[0] * (x * d2cur - u2 * dcur - xu * cur) + al[1] * d2prev;

                  ptr[l + 1] = ALF_TRUE(next, f1, f2, tz);
                  dptr[l + 1] = ALF_TRUE(dnext, f1, f2, tz);
                  d2ptr[l + 1] = ALF_TRUE(d2next, f1, f2, tz);
                }
              return;
            }

          prev = (al[0] * x) * cur + al[1] * prev;
          dprev = al[0] * (x * dcur - usq * cur) + al[1] * dprev;
          d2prev = al[0] * (x * d2cur - u2 * dcur - xu * cur) + al[1] * d2prev;
          cur = (al[2] * x) * prev + al[3] * cur;
          dcur = al[2] * (x * dprev - usq * prev) + al[3] * dcur;
          d2cur = al[2] * (x * d2prev - u2 * dprev - xu * prev) + al[3] * d2cur;
          ptr[l + 1] = ALF_TRUE(prev, f1, f2, tz);
          dptr[l + 1] = ALF_TRUE(dprev, f1, f2, tz);
          d2ptr[l + 1] = ALF_TRUE(d2prev, f1, f2, tz);
          al += 4;
          l += 2;

          if (fabs(cur) > ALF_BIG)
            {
              prev *= ALF_RESCALE;
              cur *= ALF_RESCALE;
              dprev *= ALF_RESCALE;
              dcur *= ALF_RESCALE;
              d2prev *= ALF_RESCALE;
              d2cur *= ALF_RESCALE;
              e += ALF_RESCALE_E;
              alf_scale_state(e, &f1, &f2, &tz, &tn);
            }
        }

      /* exact conversion to true values */
      prev = (prev * f1) * f2;
      cur = (cur * f1) * f2;
      dprev = (dprev * f1) * f2;
      dcur = (dcur * f1) * f2;
      d2prev = (d2prev * f1) * f2;
      d2cur = (d2cur * f1) * f2;
      ptr[l] = cur;
      dptr[l] = dcur;
      d2ptr[l] = d2cur;
    }

  /* phase 2: degrees l+1..lmax, two per step */
  for (++l; l < lmax; l += 2)
    {
      prev = (al[0] * x) * cur + al[1] * prev;
      dprev = al[0] * (x * dcur - usq * cur) + al[1] * dprev;
      d2prev = al[0] * (x * d2cur - u2 * dcur - xu * cur) + al[1] * d2prev;
      cur = (al[2] * x) * prev + al[3] * cur;
      dcur = al[2] * (x * dprev - usq * prev) + al[3] * dcur;
      d2cur = al[2] * (x * d2prev - u2 * dprev - xu * prev) + al[3] * d2cur;
      ptr[l] = prev; ptr[l + 1] = cur;
      dptr[l] = dprev; dptr[l + 1] = dcur;
      d2ptr[l] = d2prev; d2ptr[l + 1] = d2cur;
      al += 4;
    }

  if (l == lmax)
    {
      ptr[l] = (al[0] * x) * cur + al[1] * prev;
      dptr[l] = al[0] * (x * dcur - usq * cur) + al[1] * dprev;
      d2ptr[l] = al[0] * (x * d2cur - u2 * dcur - xu * cur) + al[1] * d2prev;
    }
}

/*
gsl_sf_alf_precompute()
  Precompute recurrence factors for ALFs

Inputs: norm         - ALF normalization
        lmax         - maximum degree
        mmax         - maximum order, must be <= lmax
        flags        - GSL_SF_ALF_FLG_xxx
        output_array - (output) output array, length determined by
                       gsl_sf_alf_array_size()

Notes:
1) Precomputed factors are stored at the back of output_array, skipping
nlm entries which will be used for the ALF values

2) The back of output_array holds:
header: length ALF_HDR_SIZE (lmax, mmax, normalization, flags,
        Condon-Shortley factor, P(0,0), and a check value; the array
        functions use it to verify that output_array was initialized
        for the lmax, mmax they are called with)
alm:    length nlm
blm:    length nlm (interleaved with alm)
cl:     length lmax+1
dl:     length lmax+1
el:     length lmax+1
flm:    length nlm
glm:    length nlm (interleaved with flm)

Unused entries: a(l,m), b(l,m) for (0,0), (1,0), (1,1) and l = m, m+1;
dl[0], el[0]
*/

int
gsl_sf_alf_precompute(const gsl_sf_alf_t norm, const size_t lmax,
                      const size_t mmax, const size_t flags, double output_array[])
{
  if (mmax > lmax)
    {
      GSL_ERROR ("mmax must be <= lmax", GSL_EDOM);
    }
  else
    {
      const double csfac = (flags & GSL_SF_ALF_FLG_CSPHASE) ? -1.0 : 1.0;
      const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
      double * hdr = &output_array[nlm]; /* after the nlm entries for the ALFs */
      double * alm = hdr + ALF_HDR_SIZE;
      double * cl = alm + 2 * nlm;
      double * dl = cl + lmax + 1;
      double * el = dl + lmax + 1;
      double * flm = el + lmax + 1;
      size_t l, m, k;

      for (k = 0; k < ALF_HDR_SIZE; ++k)
        hdr[k] = 0.0;

      hdr[ALF_HDR_MAGIC] = ALF_MAGIC;
      hdr[ALF_HDR_LMAX] = (double) lmax;
      hdr[ALF_HDR_MMAX] = (double) mmax;
      hdr[ALF_HDR_NORM] = (double) norm;
      hdr[ALF_HDR_FLAGS] = (double) flags;
      hdr[ALF_HDR_CSFAC] = csfac;

      /* unused entries */
      alm[0] = alm[1] = 0.0;
      dl[0] = 0.0;
      el[0] = 0.0;

      if (norm == GSL_SF_ALF_SCHMIDT)
        {
          hdr[ALF_HDR_P00] = 1.0; /* S(0,0) */
          cl[0] = 1.0;

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = M_SQRT3;
          dl[1] = csfac;
          el[1] = 1.0;

          /* m = 0 terms */
          k = 2; /* idx(2,0) */
          for (l = 2; l <= lmax; ++l)
            {
              /* a_l^0 */
              alm[2*k] = 2.0 - 1.0 / (double) l;

              /* b_l^0 */
              alm[2*k + 1] = -(1.0 - 1.0 / (double) l);

              ++k;
            }

          /* m > 0 terms */
          k = lmax + 1; /* idx(1,1) */
          for (m = 1; m <= mmax; ++m)
            {
              /* alm and blm are unused for l=m and l=m+1 */
              k += 2;

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = ((2.0*l - 1.0) / sqrt((double) (l + m))) / sqrt((double) (l - m));

                  /* b_l^m */
                  alm[2*k + 1] = -(sqrt((double) (l + m - 1)) / sqrt((double) (l + m))) *
                                  (sqrt((double) (l - m - 1)) / sqrt((double) (l - m)));

                  ++k;
                }
            }

          /* m>0 terms for flm/glm */
          k = lmax + 1; /* idx(1,1) */
          for (m = 1; m <= mmax; ++m)
            {
              double flm_fac = (m == 1) ? M_SQRT2 : 1.0;

              for (l = m; l <= lmax; ++l)
                {
                  /* f_l^m */
                  flm[2*k] = 0.5 * flm_fac * sqrt((l+m)*(l-m+1.0));

                  /* g_l^m */
                  flm[2*k+1] = 0.5 * sqrt((l+m+1.0)*(l-m));

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = sqrt(2.0 * l + 1.0);
              dl[l] = csfac * sqrt(1.0 - 0.5 / l);
              el[l] = sqrt(0.5*l*(l+1.0));
            }
        }
      else if (norm == GSL_SF_ALF_FOURPI)
        {
          hdr[ALF_HDR_P00] = 1.0;    /* R(0,0) */
          cl[0] = M_SQRT3;

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = 2.236067977499789696409174; /* sqrt(5) */
          dl[1] = csfac * M_SQRT3;
          el[1] = 1.0;

          k = 2; /* idx(2,0) */
          for (m = 0; m <= mmax; ++m)
            {
              if (m > 0)
                {
                  /* alm and blm are unused for l=m and l=m+1 */
                  k += 2;
                }

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = (sqrt(2.0 * l + 1.0) / sqrt((double) (l + m))) *
                             (sqrt(2.0 * l - 1.0) / sqrt((double) (l - m)));

                  /* b_l^m */
                  alm[2*k + 1] = -(sqrt(l + m - 1.0) / sqrt((double) (l + m))) *
                                  (sqrt(l - m - 1.0) / sqrt((double) (l - m))) *
                                  (sqrt(2.0*l + 1.0) / sqrt(2.0*l - 3.0));

                  ++k;
                }
            }

          /* m>0 terms for flm/glm */
          k = lmax + 1; /* idx(1,1) */
          for (m = 1; m <= mmax; ++m)
            {
              double flm_fac = (m == 1) ? M_SQRT2 : 1.0;

              for (l = m; l <= lmax; ++l)
                {
                  /* f_l^m */
                  flm[2*k] = 0.5 * flm_fac * sqrt((l+m)*(l-m+1.0));

                  /* g_l^m */
                  flm[2*k+1] = 0.5 * sqrt((l+m+1.0)*(l-m));

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = sqrt(2.0 * l + 3.0);
              dl[l] = csfac * sqrt(1.0 + 0.5 / l);
              el[l] = sqrt(0.5*l*(l+1.0)); /* same as Schmidt: the factor sqrt(2l+1) is common to degree l */
            }
        }
      else if (norm == GSL_SF_ALF_SPHARM)
        {
          hdr[ALF_HDR_P00] = 0.5 / M_SQRTPI;    /* Y(0,0) */
          cl[0] = M_SQRT3;

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = sqrt(5.0);
          dl[1] = csfac * (M_SQRT3 / M_SQRT2);
          el[1] = M_SQRT2;

          k = 2; /* idx(2,0) */
          for (m = 0; m <= mmax; ++m)
            {
              if (m > 0)
                {
                  /* alm and blm are unused for l=m and l=m+1 */
                  k += 2;
                }

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = (sqrt(2.0 * l + 1.0) / sqrt((double) (l + m))) *
                             (sqrt(2.0 * l - 1.0) / sqrt((double) (l - m)));

                  /* b_l^m */
                  alm[2*k + 1] = -(sqrt(l + m - 1.0) / sqrt((double) (l + m))) *
                                  (sqrt(l - m - 1.0) / sqrt((double) (l - m))) *
                                  (sqrt(2.0*l + 1.0) / sqrt(2.0*l - 3.0));

                  ++k;
                }
            }

          /* m>0 terms for flm/glm */
          k = lmax + 1; /* idx(1,1) */
          for (m = 1; m <= mmax; ++m)
            {
              for (l = m; l <= lmax; ++l)
                {
                  /* f_l^m */
                  flm[2*k] = 0.5 * sqrt((l+m)*(l-m+1.0));

                  /* g_l^m */
                  flm[2*k+1] = 0.5 * sqrt((l+m+1.0)*(l-m));

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = sqrt(2.0 * l + 3.0);
              dl[l] = csfac * sqrt(1.0 + 0.5 / l);
              el[l] = sqrt(l*(l+1.0));
            }
        }
      else if (norm == GSL_SF_ALF_FULL)
        {
          cl[0] = M_SQRT3;
          hdr[ALF_HDR_P00] = M_SQRT1_2;         /* N(0,0) */

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = sqrt(5.0);
          dl[1] = csfac * (M_SQRT3 / M_SQRT2);
          el[1] = M_SQRT2;

          k = 2; /* idx(2,0) */
          for (m = 0; m <= mmax; ++m)
            {
              if (m > 0)
                {
                  /* alm and blm are unused for l=m and l=m+1 */
                  k += 2;
                }

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = (sqrt(2.0 * l + 1.0) / sqrt((double) (l + m))) *
                             (sqrt(2.0 * l - 1.0) / sqrt((double) (l - m)));

                  /* b_l^m */
                  alm[2*k + 1] = -(sqrt(l + m - 1.0) / sqrt((double) (l + m))) *
                                  (sqrt(l - m - 1.0) / sqrt((double) (l - m))) *
                                  (sqrt(2.0*l + 1.0) / sqrt(2.0*l - 3.0));

                  ++k;
                }
            }

          /* m>0 terms for flm/glm */
          k = lmax + 1; /* idx(1,1) */
          for (m = 1; m <= mmax; ++m)
            {
              for (l = m; l <= lmax; ++l)
                {
                  /* f_l^m */
                  flm[2*k] = 0.5 * sqrt((l+m)*(l-m+1.0));

                  /* g_l^m */
                  flm[2*k+1] = 0.5 * sqrt((l+m+1.0)*(l-m));

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = sqrt(2.0 * l + 3.0);
              dl[l] = csfac * sqrt(1.0 + 0.5 / l);
              el[l] = sqrt(l*(l+1.0));
            }
        }
      else if (norm == GSL_SF_ALF_NONE)
        {
          cl[0] = 1.0;
          hdr[ALF_HDR_P00] = 1.0;         /* P(0,0) */

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = 3.0;
          dl[1] = csfac;
          el[1] = 1.0;

          k = 2; /* idx(2,0) */
          for (m = 0; m <= mmax; ++m)
            {
              if (m > 0)
                {
                  /* alm and blm are unused for l=m and l=m+1 */
                  k += 2;
                }

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = (2.0 * l - 1.0) / ((double) (l - m));

                  /* b_l^m */
                  alm[2*k + 1] = -(l + m - 1.0) / ((double) (l - m));

                  ++k;
                }
            }

          /* m>0 terms for flm/glm */
          k = lmax + 1; /* idx(1,1) */
          for (m = 1; m <= mmax; ++m)
            {
              for (l = m; l <= lmax; ++l)
                {
                  /* f_l^m */
                  flm[2*k] = 0.5*(l+m)*(l-m+1.0);

                  /* g_l^m */
                  flm[2*k+1] = 0.5;

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = 2.0 * l + 1.0;
              dl[l] = csfac * (2.0 * l - 1.0);
              el[l] = 1.0;
            }
        }
      else
        {
          GSL_ERROR ("unknown normalization", GSL_EDOM);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_alf_array_size()
  Compute total size of array needed for array functions

size = nlm    +  // for ALFs
       8      +  // header
       nlm    +  // for alm factors
       nlm    +  // for blm factors
       lmax+1 +  // for cl factors
       lmax+1 +  // for dl factors
       lmax+1 +  // for el factors
       nlm    +  // for flm factors
       nlm       // for glm factors

Inputs: lmax - maximum degree
*/

size_t
gsl_sf_alf_array_size(const size_t lmax, const size_t mmax)
{
  const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
  const size_t size = 5 * nlm + 3 * (lmax + 1) + ALF_HDR_SIZE;
  return size;
}

/*
gsl_sf_alf_array()
  Compute array of associated Legendre functions at a
given point x

Inputs: lmax         - maximum degree
        mmax         - maximum order (<= lmax)
        x            - input point
        result_array - (input/output) output array of ALFs

Notes:
1) result_array must be of length returned by gsl_sf_alf_array_size()
and be initialized by gsl_sf_alf_precompute()

2) result_array is indexed by gsl_sf_alf_array_index(l,m,lmax),
result_array[index(l,m,lmax)] = Plm(x)
*/

int
gsl_sf_alf_array(const size_t lmax, const size_t mmax, const double x,
                 double result_array[])
{
  if (mmax > lmax)
    {
      GSL_ERROR ("mmax must be <= lmax", GSL_EDOM);
    }
  else if (x < -1.0 || x > 1.0)
    {
      GSL_ERROR ("x is outside [-1,1]", GSL_EDOM);
    }
  else
    {
      /* interior point */

      const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
      const double * hdr = alf_header(result_array, nlm, lmax, mmax);
      const double * alm = &result_array[nlm + ALF_HDR_SIZE];
      const double * cl = alm + 2 * nlm;
      const double * dl = cl + lmax + 1;
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const size_t Lp2 = lmax + 2;
      size_t l, m;
      double plm,      /* P(l,m) */
             pmm;      /* P(m,m) 2^-pexp */
      double plmp1;    /* P(l+1,m) */
      int pexp = 0;    /* exponent of the scaled seed pmm */
      size_t idxmm;    /* idx(m,m) */
      const double *al;

      if (hdr == NULL)
        {
          GSL_ERROR ("result_array not initialized by gsl_sf_alf_precompute() for this lmax, mmax",
                     GSL_EINVAL);
        }

      /* initial values P(0,0) and P(1,0) */

      plm = hdr[ALF_HDR_P00];            /* P(0,0) */
      result_array[0] = plm;

      /* check for quick return */
      if (lmax == 0)
        return GSL_SUCCESS;

      plmp1 = cl[0] * x * plm; /* P(1,0) */
      result_array[1] = plmp1;

      /* Compute P(l,0), l=2:lmax */

      l = 2;        /* idx(2,0) */
      al = &alm[4]; /* a(2,0) */
      while (l < lmax)
        {
          plm   = (al[0] * x) * plmp1 + al[1] * plm;
          plmp1 = (al[2] * x) * plm   + al[3] * plmp1;

          result_array[l] = plm; result_array[l+1] = plmp1;
          l += 2; al += 4;
        }

      if (l == lmax)
        result_array[l] = (al[0] * x) * plmp1 + al[1] * plm;

      if (mmax == 0)
        return GSL_SUCCESS;

      /* columns m >= 1: seed P(m,m) = d_m u P(m-1,m-1), carried as pmm 2^pexp */

      pmm = hdr[ALF_HDR_P00];
      idxmm = 0; /* idx(0,0) */

      for (m = 1; m <= mmax; ++m)
        {
          idxmm += Lp2 - m; /* idx(m,m) = idx(m-1,m-1) + L + 2 - m */
          pmm *= dl[m] * u;

          if (fabs(pmm) < ALF_SMALL && pmm != 0.0)
            {
              int e;

              pmm = frexp(pmm, &e);
              pexp += e;
            }

          alf_P_column(lmax, m, x, pmm, pexp, cl[m], &alm[2 * (idxmm + 2)], result_array + idxmm - m);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_alf_deriv_array()
  Compute array of associated Legendre functions P and derivatives P' at a
given point x

Inputs: lmax               - maximum degree
        mmax               - maximum order (<= lmax)
        x                  - input point in (-1,1)
        result_array       - (input/output) output array of ALFs
        result_deriv_array - (input/output) output array of d/dx Plm(x)

Notes:
1) result_array must be of length returned by gsl_sf_alf_array_size()
and be initialized by gsl_sf_alf_precompute()

2) output arrays are indexed by gsl_sf_alf_array_index(l,m,lmax),
result_array[index(l,m,lmax)] = Plm(x)
result_deriv_array[index(l,m,lmax)] = d/dx Plm(x)
*/

int
gsl_sf_alf_deriv_array(const size_t lmax, const size_t mmax, const double x,
                       double result_array[], double result_deriv_array[])
{
  if (mmax > lmax)
    {
      GSL_ERROR ("mmax must be <= lmax", GSL_EDOM);
    }
  else if (x < -1.0 || x > 1.0)
    {
      GSL_ERROR ("x is outside [-1,1]", GSL_EDOM);
    }
  else if (x == -1.0 || x == 1.0) /* endpoints */
    {
      GSL_ERROR ("x cannot equal endpoints", GSL_EDOM);
    }
  else
    {
      /* interior point */

      const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
      const double * hdr = alf_header(result_array, nlm, lmax, mmax);
      const double * alm = &result_array[nlm + ALF_HDR_SIZE];
      const double * cl = alm + 2 * nlm;
      const double * dl = cl + lmax + 1;
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const double ratio = x / u;
      const size_t Lp2 = lmax + 2;
      size_t l, m, k;
      double plm,      /* P(l,m) */
             pmm;      /* P(m,m) */
      double plmp1;    /* P(l+1,m) */
      double dplm,
             dplmp1,
             dpmm;
      int pexp = 0;    /* exponent of the scaled seeds pmm, dpmm */
      size_t idxmm;    /* idx(m,m) */
      const double *al;

      if (hdr == NULL)
        {
          GSL_ERROR ("result_array not initialized by gsl_sf_alf_precompute() for this lmax, mmax",
                     GSL_EINVAL);
        }

      /* initial values P(0,0) and dP(0,0)/dx */
      plm = hdr[ALF_HDR_P00];                /* P(0,0) */
      dplm = 0.0;                  /* dP(0,0)/dx */
      result_array[0] = plm;
      result_deriv_array[0] = dplm;

      /* check for quick return */
      if (lmax == 0)
        return GSL_SUCCESS;

      /* compute P(1,0) and dP(1,0)/dx */
      plmp1 = cl[0] * x * plm; /* P(1,0) */
      dplmp1 = cl[0] * plm;    /* dP(1,0)/dx */
      result_array[1] = plmp1;
      result_deriv_array[1] = dplmp1;

      /* compute P(l,0) and dP(l,0)/dx, l=2:lmax */

      k = 2; /* idx(2,0) */
      al = &alm[k << 1];
      for (l = 2; l < lmax; l += 2)
        {
          plm    = (al[0] * x) * plmp1 + al[1] * plm;
          dplm   = al[0] * (x * dplmp1 + plmp1) + al[1] * dplm;

          plmp1  = (al[2] * x) * plm + al[3] * plmp1;
          dplmp1 = al[2] * (x * dplm + plm) + al[3] * dplmp1;

          result_array[k] = plm;        result_array[k+1] = plmp1;
          result_deriv_array[k] = dplm; result_deriv_array[k+1] = dplmp1;
          k += 2; al += 4;
        }

      if (l == lmax)
        {
          result_array[k] = (al[0] * x) * plmp1 + al[1] * plm;
          result_deriv_array[k] = al[0] * (x * dplmp1 + plmp1) + al[1] * dplm;
        }

      if (mmax == 0)
        return GSL_SUCCESS;

      /* columns m >= 1: seeds P(m,m) = d_m u P(m-1,m-1) and P'(m,m), carried as (pmm, dpmm) 2^pexp */

      pmm = hdr[ALF_HDR_P00];
      dpmm = 0.0;
      idxmm = 0; /* idx(0,0) */

      for (m = 1; m <= mmax; ++m)
        {
          idxmm += Lp2 - m; /* idx(m,m) = idx(m-1,m-1) + L + 2 - m */
          dpmm = dl[m] * (-ratio * pmm + u * dpmm);
          pmm *= dl[m] * u;

          if (fabs(pmm) < ALF_SMALL && pmm != 0.0)
            {
              int e;

              pmm = frexp(pmm, &e);
              dpmm = ldexp(dpmm, -e);
              pexp += e;
            }

          alf_dx_column(lmax, m, x, pmm, dpmm, 0.0, pexp, cl[m], &alm[2 * (idxmm + 2)],
                        result_array + idxmm - m, result_deriv_array + idxmm - m, NULL);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_alf_deriv2_array()
  Compute array of associated Legendre functions P, 1st derivatives P',
and second derivatives P'' at a given point x

Inputs: lmax                - maximum degree
        mmax                - maximum order (<= lmax)
        x                   - input point in (-1,1)
        result_array        - (input/output) output array of ALFs
        result_deriv_array  - (input/output) output array of d/dx Plm(x), size nlm
        result_deriv2_array - (input/output) output array of d^2/dx^2 Plm(x), size nlm

Notes:
1) result_array must be of length returned by gsl_sf_alf_array_size()
and be initialized by gsl_sf_alf_precompute()

2) output arrays are indexed by gsl_sf_alf_array_index(l,m,lmax),
result_array[index(l,m,lmax)] = Plm(x)
result_deriv_array[index(l,m,lmax)] = d/dx Plm(x)
result_deriv2_array[index(l,m,lmax)] = d^2/dx^2 Plm(x)
*/

int
gsl_sf_alf_deriv2_array(const size_t lmax, const size_t mmax, const double x,
                        double result_array[], double result_deriv_array[],
                        double result_deriv2_array[])
{
  if (mmax > lmax)
    {
      GSL_ERROR ("mmax must be <= lmax", GSL_EDOM);
    }
  else if (x < -1.0 || x > 1.0)
    {
      GSL_ERROR ("x is outside [-1,1]", GSL_EDOM);
    }
  else if (x == -1.0 || x == 1.0) /* endpoints */
    {
      GSL_ERROR ("x cannot equal endpoints", GSL_EDOM);
    }
  else
    {
      /* interior point */

      const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
      const double * hdr = alf_header(result_array, nlm, lmax, mmax);
      const double * alm = &result_array[nlm + ALF_HDR_SIZE];
      const double * cl = alm + 2 * nlm;
      const double * dl = cl + lmax + 1;
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const double ratio = x / u;
      const double uinv3 = 1.0 / u / u / u;
      const size_t Lp2 = lmax + 2;
      size_t l, m, k;
      double plm,      /* P(l,m) */
             pmm;      /* P(m,m) */
      double plmp1;    /* P(l+1,m) */
      double dplm,     /* P'(l,m) */
             dplmp1,   /* P'(l+1,m) */
             dpmm;     /* P'(m,m) */
      double d2plm,    /* P''(l,m) */
             d2plmp1,  /* P''(l+1,m) */
             d2pmm;    /* P''(m,m) */
      int pexp = 0;    /* exponent of the scaled seeds pmm, dpmm, d2pmm */
      size_t idxmm;    /* idx(m,m) */
      const double *al;

      if (hdr == NULL)
        {
          GSL_ERROR ("result_array not initialized by gsl_sf_alf_precompute() for this lmax, mmax",
                     GSL_EINVAL);
        }

      /* initial values P(0,0) and dP(0,0)/dx */
      plm = hdr[ALF_HDR_P00];                /* P(0,0) */
      dplm = 0.0;                  /* P'(0,0) */
      d2plm = 0.0;                 /* P''(0,0) */
      result_array[0] = plm;
      result_deriv_array[0] = dplm;
      result_deriv2_array[0] = d2plm;

      /* check for quick return */
      if (lmax == 0)
        return GSL_SUCCESS;

      /* compute P(1,0), P'(1,0), P''(1,0) */
      plmp1 = cl[0] * x * plm; /* P(1,0) */
      dplmp1 = cl[0] * plm;    /* P'(1,0) */
      d2plmp1 = 0.0;           /* P''(1,0) */
      result_array[1] = plmp1;
      result_deriv_array[1] = dplmp1;
      result_deriv2_array[1] = d2plmp1;

      /* compute P(l,0), P'(l,0), P''(l,0), l=2:lmax */

      k = 2; /* idx(2,0) */
      al = &alm[k << 1];
      for (l = 2; l < lmax; l += 2)
        {
          plm     = (al[0] * x) * plmp1 + al[1] * plm;
          dplm    = al[0] * (x * dplmp1 + plmp1) + al[1] * dplm;
          d2plm   = al[0] * (x * d2plmp1 + 2.0 * dplmp1) + al[1] * d2plm;

          plmp1   = (al[2] * x) * plm + al[3] * plmp1;
          dplmp1  = al[2] * (x * dplm + plm) + al[3] * dplmp1;
          d2plmp1 = al[2] * (x * d2plm + 2.0 * dplm) + al[3] * d2plmp1;

          result_array[k] = plm;          result_array[k+1] = plmp1;
          result_deriv_array[k] = dplm;   result_deriv_array[k+1] = dplmp1;
          result_deriv2_array[k] = d2plm; result_deriv2_array[k+1] = d2plmp1;
          k += 2; al += 4;
        }

      if (l == lmax)
        {
          result_array[k]        = (al[0] * x) * plmp1 + al[1] * plm;
          result_deriv_array[k]  = al[0] * (x * dplmp1 + plmp1) + al[1] * dplm;
          result_deriv2_array[k] = al[0] * (x * d2plmp1 + 2.0 * dplmp1) + al[1] * d2plm;
        }

      if (mmax == 0)
        return GSL_SUCCESS;

      /* columns m >= 1: seeds P(m,m) = d_m u P(m-1,m-1), P'(m,m), P''(m,m),
       * carried as (pmm, dpmm, d2pmm) 2^pexp */

      pmm = hdr[ALF_HDR_P00];
      dpmm = 0.0;
      d2pmm = 0.0;
      idxmm = 0; /* idx(0,0) */

      for (m = 1; m <= mmax; ++m)
        {
          idxmm += Lp2 - m; /* idx(m,m) = idx(m-1,m-1) + L + 2 - m */
          d2pmm = dl[m] * (-uinv3 * pmm - 2.0 * ratio * dpmm + u * d2pmm);
          dpmm = dl[m] * (-ratio * pmm + u * dpmm);
          pmm *= dl[m] * u;

          if (fabs(pmm) < ALF_SMALL && pmm != 0.0)
            {
              int e;

              pmm = frexp(pmm, &e);
              dpmm = ldexp(dpmm, -e);
              d2pmm = ldexp(d2pmm, -e);
              pexp += e;
            }

          alf_dx_column(lmax, m, x, pmm, dpmm, d2pmm, pexp, cl[m], &alm[2 * (idxmm + 2)],
                        result_array + idxmm - m, result_deriv_array + idxmm - m,
                        result_deriv2_array + idxmm - m);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_alf_vsh_array()
  Compute array of associated Legendre functions and their
theta derivatives at a given point x suitable for vector spherical
harmonic expansions

Inputs: lmax               - maximum degree
        mmax               - maximum order
        x                  - input point
        result_array       - (input/output) output array of ALFs
                             Plm(x)            for m = 0
                             Plm(x)/sin(theta) for m > 0
        result_deriv_array - (output) output array of ALF derivatives
                             d/dtheta Plm(x)

Notes:
1) result_array must be of length returned by gsl_sf_alf_array_size()
and be initialized by gsl_sf_alf_precompute()

2) result_array and result_deriv_array are indexed by gsl_sf_alf_array_index(l,m,lmax),
result_array[index(l,0,lmax)]       = Pl0(x)
result_array[index(l,m,lmax)]       = Plm(x)/sin(theta), m > 0
result_deriv_array[index(l,m,lmax)] = dPlm(x)/dtheta
*/

int
gsl_sf_alf_vsh_array(const size_t lmax, const size_t mmax, const double x,
                     double result_array[], double result_deriv_array[])
{
  if (mmax > lmax)
    {
      GSL_ERROR ("mmax must be <= lmax", GSL_EDOM);
    }
  else if (x < -1.0 || x > 1.0)
    {
      GSL_ERROR ("x is outside [-1,1]", GSL_EDOM);
    }
  else
    {
      const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
      const double * hdr = alf_header(result_array, nlm, lmax, mmax);
      const double * alm = &result_array[nlm + ALF_HDR_SIZE];
      const double * cl = alm + 2 * nlm;
      const double * dl = cl + lmax + 1;
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const double usq = u * u;
      const size_t Lp1 = lmax + 1;
      size_t l, m;
      double plm,      /* P(l,m) */
             pmm;      /* P(m,m) */
      double plmp1;    /* P(l+1,m) */
      double dplm,
             dplmp1,
             dpmm;
      int pexp = 0;    /* exponent of the scaled seeds pmm, dpmm */
      size_t idxmm;    /* idx(m,m) */
      const double *al;

      if (hdr == NULL)
        {
          GSL_ERROR ("result_array not initialized by gsl_sf_alf_precompute() for this lmax, mmax",
                     GSL_EINVAL);
        }

      /* initial values P(0,0) and dP(0,0)/dtheta */
      plm = hdr[ALF_HDR_P00];                /* P(0,0) */
      dplm = 0.0;                  /* dP(0,0)/dtheta */
      result_array[0] = plm;
      result_deriv_array[0] = dplm;

      /* check for quick return */
      if (lmax == 0)
        return GSL_SUCCESS;

      /* compute P(1,0) and dP(1,0)/dtheta */
      plmp1 = cl[0] * x * plm;    /* P(1,0) */
      dplmp1 = -u * cl[0] * plm;  /* dP(1,0)/dtheta */
      result_array[1] = plmp1;
      result_deriv_array[1] = dplmp1;

      /* Compute P(l,0), l=2:lmax */

      l = 2;        /* idx(2,0) */
      al = &alm[4]; /* a(2,0) */
      while (l < lmax)
        {
          plm    = (al[0] * x) * plmp1 + al[1] * plm;
          dplm   = al[0] * (x * dplmp1 - u * plmp1) + al[1] * dplm;

          plmp1  = (al[2] * x) * plm + al[3] * plmp1;
          dplmp1 = al[2] * (x * dplm - u * plm) + al[3] * dplmp1;

          result_array[l] = plm; result_array[l+1] = plmp1;
          result_deriv_array[l] = dplm; result_deriv_array[l+1] = dplmp1;
          l += 2; al += 4;
        }

      if (l == lmax)
        {
          result_array[l] = (al[0] * x) * plmp1 + al[1] * plm;
          result_deriv_array[l] = al[0] * (x * dplmp1 - u * plmp1) + al[1] * dplm;
        }

      if (mmax == 0)
        return GSL_SUCCESS;

      /* columns m >= 1: seeds P(m,m)/u and dP(m,m)/dtheta, carried as (pmm, dpmm) 2^pexp;
       * P(1,1)/u = d_1 P(0,0), dP(1,1)/dtheta = d_1 x P(0,0) */

      pmm = dl[1] * hdr[ALF_HDR_P00];
      dpmm = dl[1] * x * hdr[ALF_HDR_P00];
      idxmm = Lp1; /* idx(1,1) */

      for (m = 1; m <= mmax; ++m)
        {
          if (m > 1)
            {
              /* P(m,m)/u = d_m u P(m-1,m-1)/u */
              idxmm += Lp1 - m + 1; /* idx(m,m) = idx(m-1,m-1) + L + 1 - (m-1) */
              dpmm = dl[m] * u * (x * pmm + dpmm);
              pmm *= dl[m] * u;

              if (fabs(pmm) < ALF_SMALL && pmm != 0.0)
                {
                  int e;

                  pmm = frexp(pmm, &e);
                  dpmm = ldexp(dpmm, -e);
                  pexp += e;
                }
            }

          alf_vsh_column(lmax, m, x, usq, pmm, dpmm, pexp, cl[m], &alm[2 * (idxmm + 2)],
                         result_array + idxmm - m, result_deriv_array + idxmm - m);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_alf_vsh2_array()
  Compute array of associated Legendre functions and their first and
second derivatives with respect to theta, for any mmax

Inputs: lmax                - maximum degree
        mmax                - maximum order (<= lmax)
        x                   - input point, cos(theta), -1 <= x <= 1
        result_array        - (input/output) array initialized by gsl_sf_alf_precompute();
                              on output:
                              result_array[index(l,0,lmax)] = P(l,0)
                              result_array[index(l,m,lmax)] = P(l,m)/sin(theta), m > 0
        result_deriv_array  - (output) result_deriv_array[index(l,m,lmax)] = dP(l,m)/dtheta
        result_deriv2_array - (output) result_deriv2_array[index(l,m,lmax)] = d^2P(l,m)/dtheta^2

Notes:
1) The arrays use the conventions of gsl_sf_alf_vsh_array(); the values
are finite at the poles

2) The derivatives are computed by differentiating the recurrence
relations in theta (see alf_vsh2_column()), which needs no other orders,
so any mmax <= lmax is allowed (unlike gsl_sf_alf_theta_derivk_array)
*/

int
gsl_sf_alf_vsh2_array(const size_t lmax, const size_t mmax, const double x,
                              double result_array[], double result_deriv_array[],
                              double result_deriv2_array[])
{
  if (mmax > lmax)
    {
      GSL_ERROR ("mmax must be <= lmax", GSL_EDOM);
    }
  else if (x < -1.0 || x > 1.0)
    {
      GSL_ERROR ("x is outside [-1,1]", GSL_EDOM);
    }
  else
    {
      const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
      const double * hdr = alf_header(result_array, nlm, lmax, mmax);
      const double * alm = &result_array[nlm + ALF_HDR_SIZE];
      const double * cl = alm + 2 * nlm;
      const double * dl = cl + lmax + 1;
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const size_t Lp1 = lmax + 1;
      size_t l, m;
      double plm, plmp1, dplm, dplmp1, d2plm, d2plmp1;
      double pmm, dpmm, d2pmm, pmm_prev; /* seeds 2^-pexp */
      int pexp = 0;
      size_t idxmm;
      const double *al;

      if (hdr == NULL)
        {
          GSL_ERROR ("result_array not initialized by gsl_sf_alf_precompute() for this lmax, mmax",
                     GSL_EINVAL);
        }

      /* m = 0: P(l,0), dP/dtheta, d^2P/dtheta^2 */
      plm = hdr[ALF_HDR_P00];
      dplm = 0.0;
      d2plm = 0.0;
      result_array[0] = plm;
      result_deriv_array[0] = dplm;
      result_deriv2_array[0] = d2plm;

      if (lmax == 0)
        return GSL_SUCCESS;

      plmp1 = cl[0] * x * plm;     /* P(1,0) */
      dplmp1 = -u * cl[0] * plm;   /* dP(1,0)/dtheta */
      d2plmp1 = -x * cl[0] * plm;  /* d^2P(1,0)/dtheta^2 */
      result_array[1] = plmp1;
      result_deriv_array[1] = dplmp1;
      result_deriv2_array[1] = d2plmp1;

      l = 2;
      al = &alm[4]; /* a(2,0) */
      while (l < lmax)
        {
          plm    = (al[0] * x) * plmp1 + al[1] * plm;
          dplm   = al[0] * (x * dplmp1 - u * plmp1) + al[1] * dplm;
          d2plm  = al[0] * (x * d2plmp1 - 2.0 * u * dplmp1 - x * plmp1) + al[1] * d2plm;

          plmp1  = (al[2] * x) * plm + al[3] * plmp1;
          dplmp1 = al[2] * (x * dplm - u * plm) + al[3] * dplmp1;
          d2plmp1 = al[2] * (x * d2plm - 2.0 * u * dplm - x * plm) + al[3] * d2plmp1;

          result_array[l] = plm; result_array[l + 1] = plmp1;
          result_deriv_array[l] = dplm; result_deriv_array[l + 1] = dplmp1;
          result_deriv2_array[l] = d2plm; result_deriv2_array[l + 1] = d2plmp1;
          l += 2; al += 4;
        }

      if (l == lmax)
        {
          result_array[l] = (al[0] * x) * plmp1 + al[1] * plm;
          result_deriv_array[l] = al[0] * (x * dplmp1 - u * plmp1) + al[1] * dplm;
          result_deriv2_array[l] = al[0] * (x * d2plmp1 - 2.0 * u * dplmp1 - x * plmp1) + al[1] * d2plm;
        }

      if (mmax == 0)
        return GSL_SUCCESS;

      /* columns m >= 1: seeds q = P(m,m)/u, dP(m,m)/dtheta, d^2P(m,m)/dtheta^2, carried as
       * (pmm, dpmm, d2pmm) 2^pexp; from P(m,m) = d_m u P(m-1,m-1):
       * q(m)   = d_m u q(m-1)
       * dP(m)  = d_m u (x q(m-1) + dP(m-1))
       * d2P(m) = d_m (-u^2 q(m-1) + 2x dP(m-1) + u d2P(m-1))
       * and for m = 1 (P(0,0) not divided by u): q = d_1 P00, dP = d_1 x P00, d2P = -d_1 u P00 */

      pmm = dl[1] * hdr[ALF_HDR_P00];
      dpmm = dl[1] * x * hdr[ALF_HDR_P00];
      d2pmm = -dl[1] * u * hdr[ALF_HDR_P00];
      idxmm = Lp1; /* idx(1,1) */

      for (m = 1; m <= mmax; ++m)
        {
          if (m > 1)
            {
              idxmm += Lp1 - m + 1; /* idx(m,m) = idx(m-1,m-1) + L + 1 - (m-1) */
              pmm_prev = pmm;
              d2pmm = dl[m] * (-u * u * pmm_prev + 2.0 * x * dpmm + u * d2pmm);
              dpmm = dl[m] * u * (x * pmm_prev + dpmm);
              pmm *= dl[m] * u;

              if (fabs(pmm) < ALF_SMALL && pmm != 0.0)
                {
                  int e;

                  pmm = frexp(pmm, &e);
                  dpmm = ldexp(dpmm, -e);
                  d2pmm = ldexp(d2pmm, -e);
                  pexp += e;
                }
            }

          alf_vsh2_column(lmax, m, x, u, pmm, dpmm, d2pmm, pexp, cl[m], &alm[2 * (idxmm + 2)],
                            result_array + idxmm - m, result_deriv_array + idxmm - m,
                            result_deriv2_array + idxmm - m);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_alf_theta_derivk_array()
  Compute the kth derivative of Plm(x) with respect to theta
using Eq. 14 of [2]

Inputs: lmax         - maximum degree and order
        Plm          - array initialized by gsl_sf_alf_precompute()
        input_array  - (input) array of d^{k-1}/dtheta^{k-1} Plm(x), length nlm
        output_array - (output) d^k/dtheta^k Plm(x), length nlm

Notes:
1) Plm can equal input_array to calculate first derivatives
2) input_array must be computed using lmax == mmax, since all degrees and orders
are required to compute the derivatives
*/

int
gsl_sf_alf_theta_derivk_array(const size_t lmax, const double Plm[],
                              const double input_array[], double output_array[])
{
  if (lmax == 0)
    {
      /* quick return */
      output_array[0] = 0.0;
      return GSL_SUCCESS;
    }
  else
    {
      const size_t mmax = lmax;
      const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
      const double * hdr = alf_header(Plm, nlm, lmax, mmax);
      const double * alm = &Plm[nlm + ALF_HDR_SIZE];
      const double * cl = alm + 2 * nlm;
      const double * dl = cl + lmax + 1;
      const double * el = dl + lmax + 1;
      const double * flm = el + lmax + 1;
      double csfac;
      size_t l, m, p;

      if (hdr == NULL)
        {
          GSL_ERROR ("Plm not initialized by gsl_sf_alf_precompute() with lmax = mmax", GSL_EINVAL);
        }

      csfac = hdr[ALF_HDR_CSFAC];

      /* d^k/dtheta^k P(0,0) = 0 */
      output_array[0] = 0.0;

      /* m=0 terms */
      for (l = 1; l <= lmax; ++l)
        output_array[l] = -csfac * el[l] * input_array[lmax + l];

      /* 0<m<L terms */
      p = lmax + 1; /* idx(l,m) */
      for (m = 1; m <= mmax - 1; ++m)
        {
          for (l = m; l <= lmax; ++l)
            {
              size_t q = p + m - lmax - 1; /* idx(l,m-1) */
              size_t t = p + lmax - m;     /* idx(l,m+1) */

              output_array[p] = csfac * (flm[2*p]   * input_array[q] -
                                         flm[2*p+1] * input_array[t]);
              ++p;
            }
        }

      /* m=L terms */
      p = lmax + 1; /* idx(l,l) starting at l=1 */
      for (l = 1; l <= lmax; ++l)
        {
          size_t q = p + l - 1 - lmax; /* idx(l,l-1) */
          output_array[p] = csfac * flm[2*p] * input_array[q];
          p += lmax + 1 - l;
        }

      return GSL_SUCCESS;
    }
}
