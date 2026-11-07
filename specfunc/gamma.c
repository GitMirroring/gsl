/* specfunc/gamma.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_gamma.h>

#include "error.h"


#define LogRootTwoPi_  0.9189385332046727418


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* coefficients for gamma=7, kmax=8  Lanczos method */
static double lanczos_7_c[9] = {
  0.99999999999980993227684700473478,
  676.520368121885098567009190444019,
 -1259.13921672240287047156078755283,
  771.3234287776530788486528258894,
 -176.61502916214059906584551354,
  12.507343278686904814458936853,
 -0.13857109526572011689554707,
  9.984369578019570859563e-6,
  1.50563273514931155834e-7
};

/* complex version of Lanczos method; this is not safe for export
 * since it becomes bad in the left half-plane
 */
static
int
lngamma_lanczos_complex(double zr, double zi, gsl_sf_result * yr, gsl_sf_result * yi)
{
  int k;
  gsl_sf_result log1_r,    log1_i;
  gsl_sf_result logAg_r,   logAg_i;
  double Ag_r, Ag_i;
  double yi_tmp_val, yi_tmp_err;

  zr -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

  Ag_r = lanczos_7_c[0];
  Ag_i = 0.0;
  for(k=1; k<=8; k++) {
    double R = zr + k;
    double I = zi;
    double a = lanczos_7_c[k] / (R*R + I*I);
    Ag_r +=  a * R;
    Ag_i -=  a * I;
  }

  gsl_sf_complex_log_e(zr + 7.5, zi, &log1_r,  &log1_i);
  gsl_sf_complex_log_e(Ag_r, Ag_i,   &logAg_r, &logAg_i);

  /* (z+0.5)*log(z+7.5) - (z+7.5) + LogRootTwoPi_ + log(Ag(z)) */
  yr->val = (zr+0.5)*log1_r.val - zi*log1_i.val - (zr+7.5) + LogRootTwoPi_ + logAg_r.val;
  yi->val = zi*log1_r.val + (zr+0.5)*log1_i.val - zi + logAg_i.val;
  yr->err = 4.0 * GSL_DBL_EPSILON * fabs(yr->val);
  yi->err = 4.0 * GSL_DBL_EPSILON * fabs(yi->val);
  yi_tmp_val = yi->val;
  yi_tmp_err = yi->err;
  gsl_sf_angle_restrict_symm_err_e(yi_tmp_val, yi);
  yi->err += yi_tmp_err;
  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg)
{
  if(zr <= 0.5) {
    /* Transform to right half plane using reflection;
     * in fact we do a little better by stopping at 1/2.
     */
    double x = 1.0-zr;
    double y = -zi;
    gsl_sf_result a, b;
    gsl_sf_result lnsin_r, lnsin_i;

    int stat_l = lngamma_lanczos_complex(x, y, &a, &b);
    int stat_s = gsl_sf_complex_logsin_e(M_PI*zr, M_PI*zi, &lnsin_r, &lnsin_i);

    if(stat_s == GSL_SUCCESS) {
      int stat_r;
      lnr->val = M_LNPI - lnsin_r.val - a.val;
      lnr->err = lnsin_r.err + a.err + 2.0 * GSL_DBL_EPSILON * fabs(lnr->val);
      arg->val = -lnsin_i.val - b.val;
      arg->err = lnsin_i.err + b.err + 2.0 * GSL_DBL_EPSILON * fabs(arg->val);
      stat_r = gsl_sf_angle_restrict_symm_e(&(arg->val));
      return GSL_ERROR_SELECT_2(stat_r, stat_l);
    }
    else {
      DOMAIN_ERROR_2(lnr,arg);
    }
  }
  else {
    /* otherwise plain vanilla Lanczos */
    return lngamma_lanczos_complex(zr, zi, lnr, arg);
  }
}


int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(m > n) {
    DOMAIN_ERROR(result);
  }
  else if(m == n || m == 0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result nf;
    gsl_sf_result mf;
    gsl_sf_result nmmf;
    if(m*2 > n) m = n-m;
    gsl_sf_lnfact_e(n, &nf);
    gsl_sf_lnfact_e(m, &mf);
    gsl_sf_lnfact_e(n-m, &nmmf);
    result->val  = nf.val - mf.val - nmmf.val;
    result->err  = nf.err + mf.err + nmmf.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result)
{
  if(m > n) {
    DOMAIN_ERROR(result);
  }
  else if(m == n || m == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if (n <= GSL_SF_FACT_NMAX) {
    result->val = (gsl_sf_fact(n) / gsl_sf_fact(m)) / gsl_sf_fact(n-m);
    result->err = 6.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  } else {
    if(m*2 < n) m = n-m;

    if (n - m < 64)  /* compute product for a manageable number of terms */
      {
        double prod = 1.0;
        unsigned int k;
        
        for(k=n; k>=m+1; k--) {
          double tk = (double)k / (double)(k-m);
          if(tk > GSL_DBL_MAX/prod) {
            OVERFLOW_ERROR(result);
          }
          prod *= tk;
        }
        result->val = prod;
        result->err = 2.0 * GSL_DBL_EPSILON * prod * fabs((double) (n-m));
        return GSL_SUCCESS;
      }
    else
      {
        gsl_sf_result lc;
        const int stat_lc = gsl_sf_lnchoose_e (n, m, &lc);
        const int stat_e  = gsl_sf_exp_err_e(lc.val, lc.err, result);
        return GSL_ERROR_SELECT_2(stat_lc, stat_e);
      }
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_choose(unsigned int n, unsigned int m)
{
  EVAL_RESULT(gsl_sf_choose_e(n, m, &result));
}

double gsl_sf_lnchoose(unsigned int n, unsigned int m)
{
  EVAL_RESULT(gsl_sf_lnchoose_e(n, m, &result));
}
