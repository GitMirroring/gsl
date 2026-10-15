/* specfunc/test_alf.c
 * 
 * Copyright (C) 2023, 2024, 2025, 2026 Patrick Alken
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
#include <stdarg.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_sf.h>
#include <assert.h>
#include "test_sf.h"

#include "test_alf_e.c"
#include "test_alf_theta_e.c"

static double
test_alf_dx(const size_t l)
{
  const double dx_max = 0.4;
  double dx;

  if (l < 1000)
    dx = exp((double)l / 1000.0) / exp(2.0);
  else
    dx = dx_max;

  return dx;
}

/*
test_alf_sum()
  This routine computes the sum:

  Sum_{m=0}^l [P(l,m)(x)]^2

This sum should equate to 1.0 for Schmidt semi-normalized
ALFs for all l.
*/

static double
test_alf_sum(const size_t lmax, const size_t l, double *p)
{
  double sum = 0.0;
  size_t m;

  for (m = 0; m <= l; ++m)
    {
      size_t idx = gsl_sf_alf_array_index(l, m, lmax);
      sum += p[idx] * p[idx];
    }

  return sum;
}

/*
test_alf_sum_deriv()
  This routine computes the sum:

  Sum_{m=0}^l P(l,m)(x) * dP(l,m)/dx

which should equal 0 in the case of Schmidt normalized ALFs.
*/

static double
test_alf_sum_deriv(const size_t lmax, const size_t l, double *p, double *dp)
{
  double sum = 0.0;
  size_t m;

  for (m = 0; m <= l; ++m)
    {
      size_t idx = gsl_sf_alf_array_index(l, m, lmax);
      sum += p[idx] * dp[idx];
    }

  return sum;
}

/*
test_alf_sum_deriv2()
  This routine computes the sum:

  Sum_{m=0}^l P(l,m)(x) * dP(l,m)/dx

which should equal 0 in the case of Schmidt normalized ALFs.
*/

static double
test_alf_sum_deriv2(const size_t lmax, const size_t l, double *p, double *dp, double *d2p)
{
  double sum = 0.0;
  size_t m;

  for (m = 0; m <= l; ++m)
    {
      size_t idx = gsl_sf_alf_array_index(l, m, lmax);
      sum += dp[idx] * dp[idx] + p[idx] * d2p[idx];
    }

  return sum;
}

static void
test_value(const size_t lmax, const size_t mmax, const size_t l, const size_t m,
           const double *p, const double expected, const double tol,
           const char *desc, ...)
{
  size_t idx = gsl_sf_alf_array_index(l, m, lmax);
  va_list args;
  char buf[512];

  if (l > lmax || m > mmax)
    return;

  va_start(args, desc);

  vsnprintf(buf, sizeof(buf), desc, args);

  gsl_test_rel(p[idx], expected, tol, "%s lmax=%zu mmax=%zu l=%zu m=%zu", buf, lmax, mmax, l, m);

  va_end(args);
}

/* Y_{lm} = factor * S_{lm} */
static double
test_factor_spharm(const size_t l, const size_t m)
{
  double factor = sqrt( (2.0 * l + 1.0) / 4.0 / M_PI);

  if (m == 0)
    return factor;
  else
    return (factor / sqrt(2.0));
}

/* N_{lm} = factor * S_{lm} */
static double
test_factor_full(const size_t l, const size_t m)
{
  double factor = sqrt(l + 0.5);

  if (m == 0)
    return factor;
  else
    return (factor / sqrt(2.0));
}

/* R_{lm} = factor * S_{lm} */
static double
test_factor_fourpi(const size_t l, const size_t m)
{
  double factor = sqrt(2.0 * l + 1.0);

  (void) m;

  return factor;
}

/* test that p = factor * p_expected */
static int
test_legendre_compare(const double tol, const size_t lmax, const size_t mmax,
                      const double *p_expected, const double *p, const double x,
                      double (*factor)(const size_t l, const size_t m),
                      const char *desc, const char *desc2)
{
  size_t l, m;

  for (m = 0; m <= mmax; ++m)
    {
      for (l = m; l <= lmax; ++l)
        {
          size_t idx = gsl_sf_alf_array_index(l, m, lmax);
          double fac = (*factor)(l, m);

          if (fabs(p_expected[idx]) < GSL_DBL_MIN)
            continue;

          gsl_test_rel(p[idx] / fac, p_expected[idx], tol,
                       "%s %s x=%g l=%zu m=%zu", desc, desc2, x, l, m);
        }
    }

  return 0;
}

/* test gsl_sf_alf_array function */
static int
test_alf_schmidt1(const double tol, const size_t lmax, const size_t mmax, const double x,
                  const size_t expected_L, const double expected_S[], const char * desc)
{
  const size_t plm_size = gsl_sf_alf_array_size(lmax, mmax);
  double * Plm = malloc(plm_size * sizeof(double));
  size_t l, m;

  /* precompute ALF factors */
  gsl_sf_alf_precompute(GSL_SF_ALF_SCHMIDT, lmax, mmax, 0, Plm);

  /* compute ALFs */
  gsl_sf_alf_array(lmax, mmax, x, Plm);

  for (l = 0; l <= GSL_MIN(lmax, expected_L); ++l)
    {
      for (m = 0; m <= GSL_MIN(l, mmax); ++m)
        {
          size_t idx = gsl_sf_alf_array_index(l, m, expected_L);
          test_value(lmax, mmax, l, m, Plm, expected_S[idx], tol, "[schmidt1] 0th %s x=%g", desc, x);
        }
    }

  free(Plm);

  return GSL_SUCCESS;
}

/* test gsl_sf_alf_deriv_array function */
static int
test_alf_schmidt2(const double tol, const size_t lmax, const size_t mmax, const double x,
                  const size_t expected_L, const double expected_S[],
                  const double expected_DS[], const char * desc)
{
  const size_t plm_size = gsl_sf_alf_array_size(lmax, mmax);
  const size_t dplm_size = gsl_sf_alf_nlm(lmax, mmax);
  double * Plm = malloc(plm_size * sizeof(double));
  double * dPlm = malloc(dplm_size * sizeof(double));
  size_t l, m;

  /* precompute ALF factors */
  gsl_sf_alf_precompute(GSL_SF_ALF_SCHMIDT, lmax, mmax, 0, Plm);

  /* compute ALFs */
  gsl_sf_alf_deriv_array(lmax, mmax, x, Plm, dPlm);

  for (l = 0; l <= GSL_MIN(lmax, expected_L); ++l)
    {
      for (m = 0; m <= GSL_MIN(l, mmax); ++m)
        {
          size_t idx = gsl_sf_alf_array_index(l, m, expected_L);
          test_value(lmax, mmax, l, m, Plm, expected_S[idx], tol, "[schmidt2] 0th %s x=%g", desc, x);
          test_value(lmax, mmax, l, m, dPlm, expected_DS[idx], tol, "[schmidt2] 1st %s x=%g", desc, x);
        }
    }

  free(Plm);
  free(dPlm);

  return GSL_SUCCESS;
}

/* test gsl_sf_alf_deriv2_array function */
static int
test_alf_schmidt3(const double tol, const size_t lmax, const size_t mmax, const double x,
                  const size_t expected_L, const double expected_S[],
                  const double expected_DS[], const double expected_DDS[],
                  const char * desc)
{
  const size_t plm_size = gsl_sf_alf_array_size(lmax, mmax);
  const size_t dplm_size = gsl_sf_alf_nlm(lmax, mmax);
  double * Plm = malloc(plm_size * sizeof(double));
  double * dPlm = malloc(dplm_size * sizeof(double));
  double * d2Plm = malloc(dplm_size * sizeof(double));
  size_t l, m;

  /* precompute ALF factors */
  gsl_sf_alf_precompute(GSL_SF_ALF_SCHMIDT, lmax, mmax, 0, Plm);

  /* compute ALFs */
  gsl_sf_alf_deriv2_array(lmax, mmax, x, Plm, dPlm, d2Plm);

  for (l = 0; l <= GSL_MIN(lmax, expected_L); ++l)
    {
      for (m = 0; m <= GSL_MIN(l, mmax); ++m)
        {
          size_t idx = gsl_sf_alf_array_index(l, m, expected_L);
          test_value(lmax, mmax, l, m, Plm, expected_S[idx], tol, "[schmidt3] 0th %s x=%g", desc, x);
          test_value(lmax, mmax, l, m, dPlm, expected_DS[idx], tol, "[schmidt3] 1st %s x=%g", desc, x);
          test_value(lmax, mmax, l, m, d2Plm, expected_DDS[idx], tol, "[schmidt3] 2nd %s x=%g", desc, x);
        }
    }

  free(Plm);
  free(dPlm);
  free(d2Plm);

  return GSL_SUCCESS;
}

/* test gsl_sf_alf_vsh_array function, this function requires sin(theta) != 0 */
static int
test_alf_schmidt4(const double tol, const size_t lmax, const size_t mmax, const double theta,
                  const size_t expected_L, const double * expected_S, const double * expected_DS,
                  const double * expected_DDS, const char * desc)
{
  const double cost = cos(theta);
  const double sint = sin(theta);
  const size_t plm_size = gsl_sf_alf_array_size(lmax, mmax);
  const size_t dplm_size = gsl_sf_alf_nlm(lmax, mmax);
  double * Plm = malloc(plm_size * sizeof(double));
  double * dPlm = malloc(dplm_size * sizeof(double));
  double * d2Plm = malloc(dplm_size * sizeof(double));
  size_t l, m;

  /* precompute ALF factors */
  gsl_sf_alf_precompute(GSL_SF_ALF_SCHMIDT, lmax, mmax, 0, Plm);

  /* compute ALFs */
  gsl_sf_alf_vsh_array(lmax, mmax, cost, Plm, dPlm);

  if (lmax == mmax)
    gsl_sf_alf_theta_derivk_array(lmax, Plm, dPlm, d2Plm);

  for (l = 0; l <= GSL_MIN(lmax, expected_L); ++l)
    {
      for (m = 0; m <= GSL_MIN(l, mmax); ++m)
        {
          size_t idx = gsl_sf_alf_array_index(l, m, expected_L);

          test_value(lmax, mmax, l, m, Plm,
                     (m == 0) ? expected_S[idx] : expected_S[idx]/sint,
                     tol, "[schmidt4] 0th %s th=%g", desc, theta);

          test_value(lmax, mmax, l, m, dPlm, expected_DS[idx], tol, "[schmidt4] 1st %s th=%g", desc, theta);

          if (lmax == mmax)
            test_value(lmax, mmax, l, m, d2Plm, expected_DDS[idx], tol, "[schmidt4] 2nd %s th=%g", desc, theta);
        }
    }

  free(Plm);
  free(dPlm);
  free(d2Plm);

  return GSL_SUCCESS;
}

static int
test_alf_schmidt(const size_t lmax, const size_t mmax, const size_t flags, const char *desc)
{
  int s = 0;
  const double tol = 1.0e-10;
  const size_t nlm = gsl_sf_alf_nlm(lmax, mmax);
  const size_t plm_size = gsl_sf_alf_array_size(lmax, mmax);
  const size_t dplm_size = nlm;
  double * Plm = malloc(plm_size * sizeof(double));
  double * Plm2 = malloc(plm_size * sizeof(double));
  double * dPlm = malloc(dplm_size * sizeof(double));
  double * dPlm2 = malloc(dplm_size * sizeof(double));
  double * d2Plm = malloc(dplm_size * sizeof(double));
  double * Plm_theta = malloc(plm_size * sizeof(double));
  double * dPlm_theta = malloc(dplm_size * sizeof(double));
  double * d2Plm_theta = malloc(dplm_size * sizeof(double));
  double x, dx;

  /* test specific values - the alf.nb Mathematica notebook is used to generated the expected results */

  /* test gsl_sf_alf_array */
  test_alf_schmidt1(tol, lmax, mmax, -1.0,  expected_L_MNB, expected_S_m1,  desc);
  test_alf_schmidt1(tol, lmax, mmax, -0.9,  expected_L_MNB, expected_S_m09, desc);
  test_alf_schmidt1(tol, lmax, mmax, -0.2,  expected_L_MNB, expected_S_m02, desc);
  test_alf_schmidt1(tol, lmax, mmax,  0.15, expected_L_MNB, expected_S_015, desc);
  test_alf_schmidt1(tol, lmax, mmax,  0.23, expected_L_MNB, expected_S_023, desc);
  test_alf_schmidt1(tol, lmax, mmax,  0.6,  expected_L_MNB, expected_S_06,  desc);
  test_alf_schmidt1(tol, lmax, mmax,  0.75, expected_L_MNB, expected_S_075, desc);
  test_alf_schmidt1(tol, lmax, mmax,  1.0,  expected_L_MNB, expected_S_1,   desc);

  /* test gsl_sf_alf_deriv_array */
  test_alf_schmidt2(tol, lmax, mmax, -0.9,  expected_L_MNB, expected_S_m09, expected_DS_m09, desc);
  test_alf_schmidt2(tol, lmax, mmax, -0.2,  expected_L_MNB, expected_S_m02, expected_DS_m02, desc);
  test_alf_schmidt2(tol, lmax, mmax,  0.15, expected_L_MNB, expected_S_015, expected_DS_015, desc);
  test_alf_schmidt2(tol, lmax, mmax,  0.23, expected_L_MNB, expected_S_023, expected_DS_023, desc);
  test_alf_schmidt2(tol, lmax, mmax,  0.6,  expected_L_MNB, expected_S_06,  expected_DS_06,  desc);
  test_alf_schmidt2(tol, lmax, mmax,  0.75, expected_L_MNB, expected_S_075, expected_DS_075, desc);

  /* test gsl_sf_alf_deriv2_array */
  test_alf_schmidt3(tol, lmax, mmax, -0.9,  expected_L_MNB, expected_S_m09, expected_DS_m09, expected_DDS_m09, desc);
  test_alf_schmidt3(tol, lmax, mmax, -0.2,  expected_L_MNB, expected_S_m02, expected_DS_m02, expected_DDS_m02, desc);
  test_alf_schmidt3(tol, lmax, mmax,  0.15, expected_L_MNB, expected_S_015, expected_DS_015, expected_DDS_015, desc);
  test_alf_schmidt3(tol, lmax, mmax,  0.23, expected_L_MNB, expected_S_023, expected_DS_023, expected_DDS_023, desc);
  test_alf_schmidt3(tol, lmax, mmax,  0.6,  expected_L_MNB, expected_S_06,  expected_DS_06,  expected_DDS_06,  desc);
  test_alf_schmidt3(tol, lmax, mmax,  0.75, expected_L_MNB, expected_S_075, expected_DS_075, expected_DDS_075, desc);

  /* test gsl_sf_alf_vsh_array and gsl_sf_alf_theta_derivk_array */
  test_alf_schmidt4(tol, lmax, mmax, 0.1,  expected_t_L_MNB, expected_t_S_01,  expected_t_DS_01,  expected_t_DDS_01,  desc);
  test_alf_schmidt4(tol, lmax, mmax, 0.49, expected_t_L_MNB, expected_t_S_049, expected_t_DS_049, expected_t_DDS_049, desc);
  test_alf_schmidt4(tol, lmax, mmax, 1.2,  expected_t_L_MNB, expected_t_S_12,  expected_t_DS_12,  expected_t_DDS_12,  desc);
  test_alf_schmidt4(tol, lmax, mmax, 2.01, expected_t_L_MNB, expected_t_S_201, expected_t_DS_201, expected_t_DDS_201, desc);
  test_alf_schmidt4(tol, lmax, mmax, 2.99, expected_t_L_MNB, expected_t_S_299, expected_t_DS_299, expected_t_DDS_299, desc);
  test_alf_schmidt4(tol, lmax, mmax, 3.10, expected_t_L_MNB, expected_t_S_310, expected_t_DS_310, expected_t_DDS_310, desc);

  gsl_sf_alf_precompute(GSL_SF_ALF_SCHMIDT, lmax, mmax, 0, Plm);

  x = 0.1;
  gsl_sf_alf_array(lmax, mmax, x, Plm);
  test_value(lmax, mmax, 2700, 500, Plm, -7.421910573369699e-03, tol, "%s x=%g", desc, x);
  test_value(lmax, mmax, 2700, 2500, Plm, 2.717612388452281e-02, tol, "%s x=%g", desc, x);
  test_value(lmax, mmax, 2700, 2700, Plm, 1.887509917445211e-07, tol, "%s x=%g", desc, x);

  x = 1.0;
  gsl_sf_alf_vsh_array(lmax, mmax, x, Plm, dPlm_theta);
  test_value(lmax, mmax, 0, 0, Plm,  1.000000000000000,  tol, "%s theta x=%g", desc, x);
  test_value(lmax, mmax, 1, 0, Plm,  1.000000000000000,  tol, "%s theta x=%g", desc, x);
  test_value(lmax, mmax, 2, 0, Plm,  1.000000000000000,  tol, "%s theta x=%g", desc, x);
  test_value(lmax, mmax, 3, 0, Plm,  1.000000000000000,  tol, "%s theta x=%g", desc, x);

  test_value(lmax, mmax, 0, 0, dPlm_theta, 0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 1, 0, dPlm_theta, 0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 1, 1, dPlm_theta, 1.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 2, 0, dPlm_theta, 0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 2, 1, dPlm_theta, M_SQRT3,            tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 2, 2, dPlm_theta, 0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 3, 0, dPlm_theta, 0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 3, 1, dPlm_theta, 2.44948974278318,   tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 3, 2, dPlm_theta, 0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 3, 3, dPlm_theta, 0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);

  if (lmax == mmax)
    {
      gsl_sf_alf_theta_derivk_array(lmax, Plm, dPlm_theta, d2Plm_theta);
      test_value(lmax, mmax, 0, 0, d2Plm_theta,  0.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 1, 0, d2Plm_theta, -1.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 1, 1, d2Plm_theta,  0.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 2, 0, d2Plm_theta, -3.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 2, 1, d2Plm_theta,  0.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 2, 2, d2Plm_theta,  M_SQRT3,           tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 3, 0, d2Plm_theta, -6.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 3, 1, d2Plm_theta,  0.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 3, 2, d2Plm_theta,  3.872983346207417, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 3, 3, d2Plm_theta,  0.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
    }

  x = -1.0;
  gsl_sf_alf_vsh_array(lmax, mmax, x, Plm, dPlm_theta);
  test_value(lmax, mmax, 0, 0, Plm,   1.000000000000000,  tol, "%s theta x=%g", desc, x);
  test_value(lmax, mmax, 1, 0, Plm,  -1.000000000000000,  tol, "%s theta x=%g", desc, x);
  test_value(lmax, mmax, 2, 0, Plm,   1.000000000000000,  tol, "%s theta x=%g", desc, x);
  test_value(lmax, mmax, 3, 0, Plm,  -1.000000000000000,  tol, "%s theta x=%g", desc, x);

  test_value(lmax, mmax, 0, 0, dPlm_theta,  0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 1, 0, dPlm_theta,  0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 1, 1, dPlm_theta, -1.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 2, 0, dPlm_theta,  0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 2, 1, dPlm_theta,  M_SQRT3,            tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 2, 2, dPlm_theta,  0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 3, 0, dPlm_theta,  0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 3, 1, dPlm_theta, -2.44948974278318,   tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 3, 2, dPlm_theta,  0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);
  test_value(lmax, mmax, 3, 3, dPlm_theta,  0.000000000000000,  tol, "%s deriv theta x=%g", desc, x);

  if (lmax == mmax)
    {
      gsl_sf_alf_theta_derivk_array(lmax, Plm, dPlm_theta, d2Plm_theta);
      test_value(lmax, mmax, 0, 0, d2Plm_theta,  0.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 1, 0, d2Plm_theta,  1.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 1, 1, d2Plm_theta,  0.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 2, 0, d2Plm_theta, -3.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 2, 1, d2Plm_theta,  0.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 2, 2, d2Plm_theta,  M_SQRT3,           tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 3, 0, d2Plm_theta,  6.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 3, 1, d2Plm_theta,  0.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 3, 2, d2Plm_theta, -3.872983346207417, tol, "%s derivk k=2 theta x=%g", desc, x);
      test_value(lmax, mmax, 3, 3, d2Plm_theta,  0.000000000000000, tol, "%s derivk k=2 theta x=%g", desc, x);
    }

  /* test array routines */
  dx = test_alf_dx(lmax);
  gsl_sf_alf_precompute(GSL_SF_ALF_SCHMIDT, lmax, mmax, flags, Plm);
  gsl_sf_alf_precompute(GSL_SF_ALF_SCHMIDT, lmax, mmax, flags, Plm2);
  gsl_sf_alf_precompute(GSL_SF_ALF_SCHMIDT, lmax, mmax, flags, Plm_theta);

  /* test deriv array routines */
  for (x = -1.0; x <= 1.0; x += dx)
    {
      double u = sqrt((1.0 - x) * (1.0 + x));
      size_t idx = 0;
      size_t l, m;

      s += gsl_sf_alf_array(lmax, mmax, x, Plm2);
      s += gsl_sf_alf_vsh_array(lmax, mmax, x, Plm_theta, dPlm_theta);

      if (x < 1.0 && x > -1.0)
        s += gsl_sf_alf_deriv_array(lmax, mmax, x, Plm, dPlm);

      for (m = 0; m <= mmax; ++m)
        {
          for (l = m; l <= lmax; ++l)
            {
              if (fabs(Plm2[idx]) < GSL_DBL_MIN)
                {
                  ++idx;
                  continue;
                }

              assert(idx == gsl_sf_alf_array_index(l, m, lmax));

              if (x < 1.0 && x > -1.0)
                {
                  /* check Plm = Plm2 */
                  gsl_test_rel(Plm[idx], Plm2[idx], tol, "%s deriv P x=%g lmax=%zu mmax=%zu l=%zu m=%zu",
                               desc, x, lmax, mmax, l, m);

                  /* check dPlm = -1/u * dPlm_theta */
                  gsl_test_rel(-u * dPlm[idx], dPlm_theta[idx], tol,
                               "%s deriv_theta x=%g lmax=%zu mmax=%zu l=%zu m=%zu",
                               desc, x, lmax, mmax, l, m);
                }

              if (m == 0)
                {
                  /* check Plm_theta = Plm2 */
                  gsl_test_rel(Plm_theta[idx], Plm2[idx], tol, "%s deriv theta P x=%g lmax=%zu mmax=%zu l=%zu m=%zu",
                               desc, x, lmax, mmax, l, m);
                }
              else
                {
                  /* check Plm_theta*sint = Plm2 */
                  gsl_test_rel(Plm_theta[idx]*u, Plm2[idx], tol, "%s deriv theta P x=%g lmax=%zu mmax=%zu l=%zu m=%zu",
                               desc, x, lmax, mmax, l, m);
                }

              ++idx;
            }
        }

      if (lmax == mmax && x < 1.0 && x > -1.0)
        {
          for (l = 0; l <= lmax; ++l)
            {
              double sum = test_alf_sum_deriv(lmax, l, Plm, dPlm);

              gsl_test_abs(sum, 0.0, tol,
                           "%s deriv lmax=%zu mmax=%zu l=%zu, x=%f, sum=%.12e", desc, lmax, mmax, l, x, sum);
            }
        }
    }

  /* test deriv2 array routines against deriv */
  for (x = -1.0 + dx; x < 1.0 - dx; x += dx)
    {
      size_t i;

      s += gsl_sf_alf_deriv_array(lmax, mmax, x, Plm, dPlm);
      s += gsl_sf_alf_deriv2_array(lmax, mmax, x, Plm2, dPlm2, d2Plm);

      /* check Plm == Plm2 and dPlm == dPlm2 */
      for (i = 0; i < nlm; ++i)
        {
          if (fabs(Plm[i]) < 1.0e3 * GSL_DBL_EPSILON)
            gsl_test_abs(Plm2[i], Plm[i], tol, "%s deriv2 0th i=%zu", desc, i);
          else
            gsl_test_rel(Plm2[i], Plm[i], tol, "%s deriv2 0th i=%zu", desc, i);

          if (fabs(dPlm[i]) < 1.0e3 * GSL_DBL_EPSILON)
            gsl_test_abs(dPlm2[i], dPlm[i], tol, "%s deriv2 1st i=%zu", desc, i);
          else
            gsl_test_rel(dPlm2[i], dPlm[i], tol, "%s deriv2 1st i=%zu", desc, i);
        }

      if (lmax == mmax)
        {
          size_t l;

          for (l = 0; l <= lmax; ++l)
            {
              double sum = test_alf_sum_deriv(lmax, l, Plm, dPlm);
              double sum2 = test_alf_sum_deriv(lmax, l, Plm2, dPlm2);
              double sum3 = test_alf_sum_deriv2(lmax, l, Plm2, dPlm2, d2Plm);

              gsl_test_abs(sum, 0.0, tol,
                           "%s deriv2 lmax=%zu mmax=%zu l=%zu, x=%f, sum=%.12e", desc, lmax, mmax, l, x, sum);
              gsl_test_abs(sum2, 0.0, tol,
                           "%s deriv2 lmax=%zu mmax=%zu l=%zu, x=%f, sum2=%.12e", desc, lmax, mmax, l, x, sum2);
              gsl_test_abs(sum3, 0.0, tol,
                           "%s deriv2 lmax=%zu mmax=%zu l=%zu, x=%f, sum3=%.12e", desc, lmax, mmax, l, x, sum3);
            }
        }
    }

  free(Plm);
  free(Plm2);
  free(dPlm);
  free(dPlm2);
  free(d2Plm);
  free(Plm_theta);
  free(dPlm_theta);
  free(d2Plm_theta);

  return s;
}

/* test other normalizations (other than schmidt) */
static int
test_alf_norm(const gsl_sf_alf_t norm, const size_t lmax,
              const size_t mmax, const size_t flags,
              const char *desc)
{
  int s = 0;
  const double tol = 1.0e-10;
  const size_t plm_size = gsl_sf_alf_array_size(lmax, mmax);
  const size_t dplm_size = gsl_sf_alf_nlm(lmax, mmax);
  double x, dx;
  double (*factor)(const size_t l, const size_t m) = NULL;

  double *p = malloc(sizeof(double) * plm_size);
  double *dp = malloc(sizeof(double) * dplm_size);
  double *p_schmidt = malloc(sizeof(double) * plm_size);
  double *dp_schmidt = malloc(sizeof(double) * dplm_size);

  if (norm == GSL_SF_ALF_SPHARM)
    {
      factor = &test_factor_spharm;
    }
  else if (norm == GSL_SF_ALF_FULL)
    {
      factor = &test_factor_full;

      /* test specific values (computed from GNU octave) */
      gsl_sf_alf_precompute(norm, lmax, mmax, 0, p);
      x = 0.45;
      s += gsl_sf_alf_array(lmax, mmax, x, p);
      test_value(lmax, mmax, 0, 0, p,  0.707106781186548,   tol, desc, "x=0.45");
      test_value(lmax, mmax, 1, 0, p,  0.551135192126215,   tol, desc, "x=0.45");
      test_value(lmax, mmax, 1, 1, p,  0.773385414912901,   tol, desc, "x=0.45");
      test_value(lmax, mmax, 2, 0, p, -0.310298495404022,   tol, desc, "x=0.45");
      test_value(lmax, mmax, 2, 1, p,  0.778204062248457,   tol, desc, "x=0.45");
      test_value(lmax, mmax, 2, 2, p,  0.772176054650104,   tol, desc, "x=0.45");
      test_value(lmax, mmax, 3, 0, p, -0.83661120632398589, tol, desc, "x=0.45");
      test_value(lmax, mmax, 3, 1, p,  0.00904294765791280, tol, desc, "x=0.45");
      test_value(lmax, mmax, 3, 2, p,  0.91934361403343767, tol, desc, "x=0.45");
      test_value(lmax, mmax, 3, 3, p,  0.74482641545541073, tol, desc, "x=0.45");
    }
  else if (norm == GSL_SF_ALF_FOURPI)
    {
      factor = &test_factor_fourpi;
    }

  /*
   * test the scale factors between the Schmidts and these
   * normalized functions
   */
  gsl_sf_alf_precompute(GSL_SF_ALF_SCHMIDT, lmax, mmax, flags, p_schmidt);
  gsl_sf_alf_precompute(norm, lmax, mmax, flags, p);

  dx = test_alf_dx(lmax);
  for (x = -1.0; x <= 1.0; x += dx)
    {
      s += gsl_sf_alf_array(lmax, mmax, x, p_schmidt);
      s += gsl_sf_alf_array(lmax, mmax, x, p);
      test_legendre_compare(tol, lmax, mmax, p_schmidt, p, x, factor, desc, "p");
    }

  /* test derivatives */
  for (x = -1.0 + dx; x < 1.0 - dx; x += dx)
    {
      s += gsl_sf_alf_deriv_array(lmax, mmax, x, p_schmidt, dp_schmidt);
      s += gsl_sf_alf_deriv_array(lmax, mmax, x, p, dp);
      test_legendre_compare(tol, lmax, mmax, p_schmidt, p, x, factor, desc, "deriv p");
      test_legendre_compare(tol, lmax, mmax, dp_schmidt, dp, x, factor, desc, "deriv dp");
    }

  free(p);
  free(dp);
  free(p_schmidt);
  free(dp_schmidt);

  return s;
}

/*
test_alf_unnorm()
  This routine tests the unnormalized ALFs using the relation

S(l,m)(x) = a(l,m) * P(l,m)(x)

where

a(l,0) = 1
a(l,1) = -sqrt(2)/sqrt(l * (l+1))
a(l,m+1) = a(l,m) / sqrt((l+m+1) * (l-m)), m > 1

and

S(l,m) are the Schmidt semi-normalized ALFs
*/

static int
test_alf_unnorm(const size_t lmax_orig, const size_t mmax_orig,
                const size_t flags, const char *desc)
{
  int s = 0;
  const double tol = 1.0e-10;
  const size_t lmax = GSL_MIN(lmax_orig, 140);
  const size_t mmax = GSL_MIN(lmax, mmax_orig);
  const size_t plm_size = gsl_sf_alf_array_size(lmax, mmax);
  const size_t dplm_size = gsl_sf_alf_nlm(lmax, mmax);
  size_t l, m;
  double x, dx;

  double *p = malloc(sizeof(double) * plm_size);
  double *dp = malloc(sizeof(double) * dplm_size);
  double *p2 = malloc(sizeof(double) * plm_size);
  double *p_schmidt = malloc(sizeof(double) * plm_size);
  double *dp_schmidt = malloc(sizeof(double) * dplm_size);

  gsl_sf_alf_precompute(GSL_SF_ALF_SCHMIDT, lmax, mmax, flags, p_schmidt);
  gsl_sf_alf_precompute(GSL_SF_ALF_NONE, lmax, mmax, flags, p);
  gsl_sf_alf_precompute(GSL_SF_ALF_NONE, lmax, mmax, flags, p2);

  dx = test_alf_dx(lmax);
  for (x = -1.0 + dx; x < 1.0 - dx; x += dx)
    {
      gsl_sf_alf_deriv_array(lmax, mmax, x, p_schmidt, dp_schmidt);
      gsl_sf_alf_deriv_array(lmax, mmax, x, p, dp);

      for (l = 0; l <= lmax; ++l)
        {
          size_t M = GSL_MIN(l, mmax);
          double a_lm = sqrt(2.0 / (double)l / (l + 1.0));
          size_t idx;

          /* test S(l,0) = P(l,0) */
          idx = gsl_sf_alf_array_index(l, 0, lmax);
          gsl_test_rel(p[idx], p_schmidt[idx], tol,
                       "%s l=%zu, m=0, x=%f", desc, l, x);
          gsl_test_rel(dp[idx], dp_schmidt[idx], tol,
                       "%s deriv l=%zu, m=0, x=%f", desc, l, x);

          /* test S(l,m) = a_{lm} * P(l,m) for m > 0 */
          for (m = 1; m <= M; ++m)
            {
              idx = gsl_sf_alf_array_index(l, m, lmax);

              gsl_test_rel(a_lm * p[idx], p_schmidt[idx], tol,
                           "%s l=%zu, m=%zu, x=%f", desc, l, m, x);
              gsl_test_abs(a_lm * dp[idx], dp_schmidt[idx], tol,
                           "%s deriv l=%zu, m=%zu, x=%f", desc, l, m, x);

              a_lm /= sqrt((double) (l + m + 1)) *
                      sqrt((double) (l - m));
            }
        }

      gsl_sf_alf_array(lmax, mmax, x, p2);

      /* test if p = p2 */
      for (m = 0; m <= mmax; ++m)
        {
          for (l = m; l <= lmax; ++l)
            {
              size_t idx = gsl_sf_alf_array_index(l, m, lmax);
              gsl_test_rel(p2[idx], p[idx], tol,
                           "%s compare l=%zu, m=%zu, x=%f",
                           desc, l, m, x);
            }
        }
    }

  free(p);
  free(p2);
  free(dp);
  free(p_schmidt);
  free(dp_schmidt);

  return s;
}

static int
test_alf_all(const size_t lmax, const size_t mmax)
{
  int s = 0;

  s += test_alf_schmidt(lmax, mmax, 0, "schmidt nocsphase");
  s += test_alf_schmidt(lmax, mmax, GSL_SF_ALF_FLG_CSPHASE, "schmidt csphase");

  s += test_alf_norm(GSL_SF_ALF_SPHARM, lmax, mmax, 0, "spharm nocsphase");
  s += test_alf_norm(GSL_SF_ALF_SPHARM, lmax, mmax, GSL_SF_ALF_FLG_CSPHASE, "spharm csphase");

  s += test_alf_norm(GSL_SF_ALF_FULL, lmax, mmax, 0, "full nocsphase");
  s += test_alf_norm(GSL_SF_ALF_FULL, lmax, mmax, GSL_SF_ALF_FLG_CSPHASE, "full csphase");

  s += test_alf_norm(GSL_SF_ALF_FOURPI, lmax, mmax, 0, "fourpi nocsphase");
  s += test_alf_norm(GSL_SF_ALF_FOURPI, lmax, mmax, GSL_SF_ALF_FLG_CSPHASE, "fourpi csphase");

  s += test_alf_unnorm(lmax, mmax, 0, "unnorm nocsphase");
  s += test_alf_unnorm(lmax, mmax, GSL_SF_ALF_FLG_CSPHASE, "unnorm csphase");

  return s;
}

/*
test_alf_check()
  The array functions verify that the array was initialized by
gsl_sf_alf_precompute() for the lmax, mmax they are called with
*/

static int
test_alf_check(void)
{
  const size_t lmax = 10;
  const double x = 0.3;
  const size_t asize = gsl_sf_alf_array_size(lmax, lmax);
  double * p = malloc(asize * sizeof(double));
  double * q = calloc(asize, sizeof(double));
  double * d1 = malloc(gsl_sf_alf_nlm(lmax, lmax) * sizeof(double));
  double * d2 = malloc(gsl_sf_alf_nlm(lmax, lmax) * sizeof(double));
  size_t i;

  gsl_sf_alf_precompute(GSL_SF_ALF_FOURPI, lmax, lmax, 0, p);

  for (i = 0; i < 3; ++i)
    {
      /* matching, smaller lmax, smaller mmax */
      const size_t L = (i == 1) ? lmax - 1 : lmax;
      const size_t M = (i == 2) ? lmax / 2 : L;
      const int expected = (i == 0) ? GSL_SUCCESS : GSL_EINVAL;

      gsl_test_int(gsl_sf_alf_array(L, M, x, p), expected, "alf check array lmax=%zu mmax=%zu", L, M);
      gsl_test_int(gsl_sf_alf_vsh_array(L, M, x, p, d1), expected, "alf check vsh_array lmax=%zu mmax=%zu", L, M);
      gsl_test_int(gsl_sf_alf_deriv_array(L, M, x, p, d1), expected, "alf check deriv_array lmax=%zu mmax=%zu", L, M);
      gsl_test_int(gsl_sf_alf_deriv2_array(L, M, x, p, d1, d2), expected, "alf check deriv2_array lmax=%zu mmax=%zu", L, M);
      gsl_test_int(gsl_sf_alf_vsh2_array(L, M, x, p, d1, d2), expected, "alf check vsh2_array lmax=%zu mmax=%zu", L, M);
    }

  /* an array never initialized by gsl_sf_alf_precompute() */
  gsl_test_int(gsl_sf_alf_array(lmax, lmax, x, q), GSL_EINVAL, "alf check uninitialized array");

  /* theta_derivk needs an array precomputed with mmax = lmax */
  gsl_sf_alf_vsh_array(lmax, lmax, x, p, d1);
  gsl_test_int(gsl_sf_alf_theta_derivk_array(lmax, p, d1, d2), GSL_SUCCESS, "alf check derivk mmax = lmax");

  gsl_sf_alf_precompute(GSL_SF_ALF_FOURPI, lmax, lmax / 2, 0, q);
  gsl_sf_alf_vsh_array(lmax, lmax / 2, x, q, d1);
  gsl_test_int(gsl_sf_alf_theta_derivk_array(lmax, q, d1, d2), GSL_EINVAL, "alf check derivk mmax < lmax");

  free(p);
  free(q);
  free(d1);
  free(d2);

  return 0;
}

/*
test_alf_derivk()
  Test gsl_sf_alf_theta_derivk_array against the other functions, for one
normalization, lmax = mmax, at colatitude theta:

k = 1: derivk applied to P(l,m) from gsl_sf_alf_array equals dP/dtheta from
       gsl_sf_alf_vsh_array
k = 2: derivk applied to dP/dtheta equals d^2P/dtheta^2 = (1-x^2) P'' - x P',
       with P' = dP/dx and P'' = d^2P/dx^2 from gsl_sf_alf_deriv2_array

Errors are relative to the largest value of the expected array.
*/

static int
test_alf_derivk(const gsl_sf_alf_t norm, const size_t flags, const size_t lmax,
                const double theta, const double tol, const char * desc)
{
  const double x = cos(theta);
  const size_t nlm = gsl_sf_alf_nlm(lmax, lmax);
  double * p = malloc(gsl_sf_alf_array_size(lmax, lmax) * sizeof(double));
  double * dth = malloc(nlm * sizeof(double));   /* dP/dtheta */
  double * dx = malloc(nlm * sizeof(double));    /* dP/dx */
  double * d2x = malloc(nlm * sizeof(double));   /* d^2P/dx^2 */
  double * d2th = malloc(nlm * sizeof(double));  /* expected d^2P/dtheta^2 */
  double * out = malloc(nlm * sizeof(double));
  double err1 = 0.0, max1 = 0.0, err2 = 0.0, max2 = 0.0;
  size_t i;

  gsl_sf_alf_precompute(norm, lmax, lmax, flags, p);

  /* each function writes only the first nlm entries of p */
  gsl_sf_alf_vsh_array(lmax, lmax, x, p, dth);
  gsl_sf_alf_deriv2_array(lmax, lmax, x, p, dx, d2x);
  gsl_sf_alf_array(lmax, lmax, x, p);

  /* k = 1: input P(l,m) */
  gsl_sf_alf_theta_derivk_array(lmax, p, p, out);
  for (i = 0; i < nlm; ++i)
    {
      err1 = GSL_MAX(err1, fabs(out[i] - dth[i]));
      max1 = GSL_MAX(max1, fabs(dth[i]));
    }

  /* k = 2: input dP/dtheta */
  for (i = 0; i < nlm; ++i)
    d2th[i] = (1.0 - x * x) * d2x[i] - x * dx[i];

  gsl_sf_alf_theta_derivk_array(lmax, p, dth, out);
  for (i = 0; i < nlm; ++i)
    {
      err2 = GSL_MAX(err2, fabs(out[i] - d2th[i]));
      max2 = GSL_MAX(max2, fabs(d2th[i]));
    }

  gsl_test(err1 > tol * max1, "alf derivk k=1 %s lmax=%zu theta=%g (relative error %.2e)",
           desc, lmax, theta, err1 / max1);
  gsl_test(err2 > tol * max2, "alf derivk k=2 %s lmax=%zu theta=%g (relative error %.2e)",
           desc, lmax, theta, err2 / max2);

  free(p);
  free(dth);
  free(dx);
  free(d2x);
  free(d2th);
  free(out);

  return 0;
}

/*
test_alf_vsh2()
  Test gsl_sf_alf_vsh2_array at colatitude theta:

- P and dP/dtheta are bitwise equal to gsl_sf_alf_vsh_array (same recurrences)
- for mmax = lmax, d^2P/dtheta^2 equals gsl_sf_alf_theta_derivk_array applied
  to dP/dtheta
- for mmax < lmax, every column equals the corresponding column of the
  mmax = lmax computation (bitwise)
*/

static int
test_alf_vsh2(const gsl_sf_alf_t norm, const size_t flags, const size_t lmax,
                const double theta, const double tol, const char * desc)
{
  const double x = cos(theta);
  const size_t mmax2 = lmax / 2;
  const size_t nlm = gsl_sf_alf_nlm(lmax, lmax);
  const size_t nlm2 = gsl_sf_alf_nlm(lmax, mmax2);
  double * p = malloc(gsl_sf_alf_array_size(lmax, lmax) * sizeof(double));
  double * p2 = malloc(gsl_sf_alf_array_size(lmax, mmax2) * sizeof(double));
  double * v = malloc(nlm * sizeof(double));
  double * dv = malloc(nlm * sizeof(double));
  double * t = malloc(nlm * sizeof(double));
  double * dt = malloc(nlm * sizeof(double));
  double * d2t = malloc(nlm * sizeof(double));
  double * d2k = malloc(nlm * sizeof(double));
  double * dt2 = malloc(nlm2 * sizeof(double));
  double * d2t2 = malloc(nlm2 * sizeof(double));
  double err = 0.0, vmax = 0.0;
  size_t i;
  int same = 1;

  gsl_sf_alf_precompute(norm, lmax, lmax, flags, p);
  gsl_sf_alf_precompute(norm, lmax, mmax2, flags, p2);

  /* reference: vsh_array (+ derivk) */
  gsl_sf_alf_vsh_array(lmax, lmax, x, p, dv);
  for (i = 0; i < nlm; ++i)
    v[i] = p[i];
  gsl_sf_alf_theta_derivk_array(lmax, p, dv, d2k);

  gsl_sf_alf_vsh2_array(lmax, lmax, x, p, dt, d2t);
  for (i = 0; i < nlm; ++i)
    {
      t[i] = p[i];
      same = same && (t[i] == v[i]) && (dt[i] == dv[i]);
      err = GSL_MAX(err, fabs(d2t[i] - d2k[i]));
      vmax = GSL_MAX(vmax, fabs(d2k[i]));
    }

  gsl_test(!same, "alf vsh2 %s lmax=%zu theta=%g P, dP equal to vsh_array", desc, lmax, theta);
  gsl_test(err > tol * vmax, "alf vsh2 %s lmax=%zu theta=%g d2P vs derivk (relative error %.2e)",
           desc, lmax, theta, err / vmax);

  /* mmax < lmax: columns m <= mmax2 equal (the index does not depend on mmax) */
  gsl_sf_alf_vsh2_array(lmax, mmax2, x, p2, dt2, d2t2);
  same = 1;
  for (i = 0; i < nlm2; ++i)
    same = same && (p2[i] == t[i]) && (dt2[i] == dt[i]) && (d2t2[i] == d2t[i]);
  gsl_test(!same, "alf vsh2 %s lmax=%zu mmax=%zu theta=%g equal to mmax=lmax columns",
           desc, lmax, mmax2, theta);

  free(p); free(p2); free(v); free(dv); free(t); free(dt); free(d2t); free(d2k); free(dt2); free(d2t2);

  return 0;
}

/* gsl_sf_alf_theta_derivk_array in all normalizations, with and without CS phase */
static int
test_alf_derivk_all(void)
{
  const gsl_sf_alf_t norms[] = { GSL_SF_ALF_SCHMIDT, GSL_SF_ALF_SPHARM, GSL_SF_ALF_FULL,
                                 GSL_SF_ALF_FOURPI, GSL_SF_ALF_NONE };
  const char * names[] = { "schmidt", "spharm", "full", "fourpi", "unnorm" };
  const double thetas[] = { 0.1, 0.7, 1.5, 2.3, 3.0 };
  size_t i, j, k;

  for (i = 0; i < 5; ++i)
    {
      for (j = 0; j < 2; ++j)
        {
          const size_t flags = j ? GSL_SF_ALF_FLG_CSPHASE : 0;
          char desc[64];

          sprintf(desc, "%s %s", names[i], j ? "csphase" : "nocsphase");

          for (k = 0; k < 5; ++k)
            {
              test_alf_derivk(norms[i], flags, 12, thetas[k], 1.0e-13, desc);
              test_alf_vsh2(norms[i], flags, 12, thetas[k], 1.0e-13, desc);

              /* unnormalized ALFs overflow at high degree */
              if (norms[i] != GSL_SF_ALF_NONE)
                {
                  test_alf_derivk(norms[i], flags, 600, thetas[k], 1.0e-11, desc);
                  test_alf_vsh2(norms[i], flags, 600, thetas[k], 1.0e-11, desc);
                }
            }

          /* poles (the x derivatives are not defined there) */
          test_alf_vsh2(norms[i], flags, 12, 0.0, 1.0e-13, desc);
          test_alf_vsh2(norms[i], flags, 12, M_PI, 1.0e-13, desc);
        }
    }

  return 0;
}

/*
test_alf_highdeg()
  High degree regression test (lmax = 2190): values in columns whose sectoral
seed P(m,m) is far below the double range, while P(l,m) is significant at
higher degrees. Reference values from 40-digit mpmath (4-pi normalization,
no Condon-Shortley phase); d^2P/dx^2 from the Legendre differential equation.
Tolerances are relative to the size of each quantity at these degrees.
*/

static int
test_alf_highdeg(void)
{
  const size_t lmax = 2190;
  const size_t asize = gsl_sf_alf_array_size(lmax, lmax);
  /* theta, l, m, P, P/sin(theta), dP/dtheta, dP/dx, d^2P/dx^2 */
  const double ref[][8] = {
    { 0.4, 2054, 793, 6.6882690752656492, 17.175023229811219, 9.5458100832329508, -24.512995527229179, -3271312.6559891376 },
    { 0.4, 2188, 845, 6.7589426098981151, 17.356508093142208, 1.5422383112334136, -3.9603638136054859, -3612621.1208961394 },
    { 0.4, 2190, 900, 2.7125821298943312e-5, 6.9657276896945853e-5, 0.020263399998632477, -0.052035042516235241, 96.510447033687591 },
    { 0.25, 2100, 560, 4.5882801888663946e-5, 0.0001854570235129539, 0.039234875199752664, -0.15858628664624085, 528.21820767679027 },
    { 0.25, 2190, 600, 2.256617080656388e-8, 9.1211841857988921e-8, 2.3677627809792891e-5, -9.5704320501330161e-5, 0.39632886968629414 }
  };
  /* typical size of each quantity at these degrees, for absolute tolerances */
  const double scale[5] = { 10.0, 30.0, 3.0e4, 1.0e5, 1.0e8 };
  const double tol = 1.0e-11;
  double * p = malloc(asize * sizeof(double));   /* FOURPI */
  double * ps = malloc(asize * sizeof(double));  /* SCHMIDT */
  double * dp = malloc(asize * sizeof(double));
  double * d2p = malloc(asize * sizeof(double));
  size_t i;

  /* the functions overwrite only the first nlm entries, not the recurrence factors */
  gsl_sf_alf_precompute(GSL_SF_ALF_FOURPI, lmax, lmax, 0, p);
  gsl_sf_alf_precompute(GSL_SF_ALF_SCHMIDT, lmax, lmax, 0, ps);

  for (i = 0; i < sizeof(ref) / sizeof(ref[0]); ++i)
    {
      const double theta = ref[i][0];
      const double x = cos(theta);
      const size_t l = (size_t) ref[i][1], m = (size_t) ref[i][2];
      const size_t idx = gsl_sf_alf_array_index(l, m, lmax);

      gsl_sf_alf_array(lmax, lmax, x, p);
      gsl_test_abs(p[idx], ref[i][3], tol * scale[0], "alf highdeg array theta=%g l=%zu m=%zu", theta, l, m);

      gsl_sf_alf_vsh_array(lmax, lmax, x, p, dp);
      gsl_test_abs(p[idx], ref[i][4], tol * scale[1], "alf highdeg vsh P/sin theta=%g l=%zu m=%zu", theta, l, m);
      gsl_test_abs(dp[idx], ref[i][5], tol * scale[2], "alf highdeg vsh dP/dtheta theta=%g l=%zu m=%zu", theta, l, m);

      gsl_sf_alf_deriv_array(lmax, lmax, x, p, dp);
      gsl_test_abs(p[idx], ref[i][3], tol * scale[0], "alf highdeg deriv P theta=%g l=%zu m=%zu", theta, l, m);
      gsl_test_abs(dp[idx], ref[i][6], tol * scale[3], "alf highdeg deriv dP/dx theta=%g l=%zu m=%zu", theta, l, m);

      gsl_sf_alf_deriv2_array(lmax, lmax, x, p, dp, d2p);
      gsl_test_abs(p[idx], ref[i][3], tol * scale[0], "alf highdeg deriv2 P theta=%g l=%zu m=%zu", theta, l, m);
      gsl_test_abs(dp[idx], ref[i][6], tol * scale[3], "alf highdeg deriv2 dP/dx theta=%g l=%zu m=%zu", theta, l, m);
      gsl_test_abs(d2p[idx], ref[i][7], tol * scale[4], "alf highdeg deriv2 d2P/dx2 theta=%g l=%zu m=%zu", theta, l, m);

      gsl_sf_alf_vsh2_array(lmax, lmax, x, p, dp, d2p);
      gsl_test_abs(p[idx], ref[i][4], tol * scale[1], "alf highdeg vsh2 P/sin theta=%g l=%zu m=%zu", theta, l, m);
      gsl_test_abs(dp[idx], ref[i][5], tol * scale[2], "alf highdeg vsh2 dP/dtheta theta=%g l=%zu m=%zu", theta, l, m);
      gsl_test_abs(d2p[idx], (1.0 - x * x) * ref[i][7] - x * ref[i][6], tol * scale[4],
                   "alf highdeg vsh2 d2P/dtheta2 theta=%g l=%zu m=%zu", theta, l, m);

      /* Schmidt normalization: S(l,m) = P(l,m) / sqrt(2l+1) */
      gsl_sf_alf_array(lmax, lmax, x, ps);
      gsl_test_abs(ps[idx], ref[i][3] / sqrt(2.0 * l + 1.0), tol * scale[0],
                   "alf highdeg schmidt array theta=%g l=%zu m=%zu", theta, l, m);
    }

  /* a value far below the double range (about 1e-522) is exactly 0 */
  gsl_sf_alf_array(lmax, lmax, cos(0.4), p);
  gsl_test_abs(p[gsl_sf_alf_array_index(1600, 1500, lmax)], 0.0, 0.0, "alf highdeg underflow to 0");

  free(p);
  free(ps);
  free(dp);
  free(d2p);

  return 0;
}

/* test associated legendre functions */
int
test_alf(void)
{
  int s = 0;
  size_t l, m;

  for (l = 0; l <= 10; ++l)
    {
      for (m = 0; m <= l; ++m)
        test_alf_all(l, m);
    }

  test_alf_all(140, 140);
  test_alf_highdeg();
  test_alf_derivk_all();
  test_alf_check();

  return s;
}
