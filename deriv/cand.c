/* deriv/cand.c
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
 * This module contains code related to filtering derivative candidates
 * and selecting the final candidate
 *
 * Filtering step: filter out derivative candidates which are extreme outliers
 * Selection step: select final candidate from remaining unfiltered candidates
 */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv_constants.h>
#include <gsl/gsl_deriv_cand.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_vector.h>

static int filter_mad(const gsl_vector * der, gsl_deriv_cand_t * cand);
static int filter_trim(const gsl_vector * der, gsl_deriv_cand_t * cand);
static int filter_none(const gsl_vector * der, gsl_deriv_cand_t * cand);

static int select_minerr(const gsl_deriv_cand_t * cand, const gsl_vector * der,
                         const gsl_vector * err, const gsl_vector * steps,
                         double * der_est, double * err_est, double * step_est);

/* candidate filtering methods */

/*
filter_mad()
  Filter derivative candidates based on MAD

Inputs: der  - computed derivative estimates, max length GSL_DERIV_MAX_STEPS
        cand - (output) structure containing candidates to keep

Return: success/error
*/

static int
filter_mad(const gsl_vector * der, gsl_deriv_cand_t * cand)
{
  if (der->size > GSL_DERIV_MAX_STEPS)
    {
      GSL_ERROR ("der size cannot exceed GSL_DERIV_MAX_STEPS", GSL_EBADLEN);
    }
  else
    {
      const double alpha = 5.0;
      const size_t n = der->size;
      double finite[GSL_DERIV_MAX_STEPS], work[GSL_DERIV_MAX_STEPS];
      double median, MAD, tol;
      size_t nfinite = 0;
      size_t i;

      /*
       * non-finite candidates (e.g. overflow at very small steps) are discarded;
       * the median and MAD are computed from the finite candidates
       */
      for (i = 0; i < n; ++i)
        {
          double di = gsl_vector_get(der, i);
          if (gsl_finite(di))
            finite[nfinite++] = di;
        }

      cand->nkeep = 0;

      if (nfinite > 0)
        {
          /* compute MAD (before median, since the median sorts the array) */
          MAD = gsl_stats_mad(finite, 1, nfinite, work);

          /* compute median */
          median = gsl_stats_median(finite, 1, nfinite);

          tol = GSL_MAX(alpha * MAD, GSL_DBL_EPSILON * fabs(median));

          for (i = 0; i < n; ++i)
            {
              double di = gsl_vector_get(der, i);

              if (gsl_finite(di) && ((fabs(di - median) <= tol) || (MAD == 0.0)))
                cand->index[cand->nkeep++] = i;
            }
        }

      if (cand->nkeep == 0)
        {
          GSL_ERROR("no suitable candidates found", GSL_FAILURE);
        }

      return GSL_SUCCESS;
    }
}

/*
filter_trim()
  Filter derivative candidates based trimming largest and smallest values

Inputs: der  - computed derivative estimates, max length GSL_DERIV_MAX_STEPS
        cand - (output) structure containing candidates to keep

Return: success/error
*/

static int
filter_trim(const gsl_vector * der, gsl_deriv_cand_t * cand)
{
  const size_t n = der->size;

  if (n > GSL_DERIV_MAX_STEPS)
    {
      GSL_ERROR ("der size cannot exceed GSL_DERIV_MAX_STEPS", GSL_EBADLEN);
    }
  else
    {
      const double trim_fraction = 0.15;
      size_t perm_data[GSL_DERIV_MAX_STEPS];
      gsl_permutation perm = gsl_permutation_view_array(perm_data, n);
      size_t ntrim;
      size_t i;

      /* number of candidates to discard from each end */
      ntrim = (size_t) (trim_fraction * n + 0.5);

      /* always leave at least one candidate */
      if (2 * ntrim >= n)
        ntrim = (n - 1) / 2;

      /* sort only to identify outliers */
      gsl_sort_vector_index(&perm, der);

      cand->nkeep = 0;
      for (i = ntrim; i < n - ntrim; ++i)
        cand->index[cand->nkeep++] = gsl_permutation_get(&perm, i);

      if (cand->nkeep == 0)
        {
          GSL_ERROR("no suitable candidates found", GSL_FAILURE);
        }

      return GSL_SUCCESS;
    }
}

/*
filter_none()
  Keep all derivative candidates provided they are finite

Inputs: der  - computed derivative estimates, max length GSL_DERIV_MAX_STEPS
        cand - (output) structure containing candidates to keep

Return: success/error
*/

static int
filter_none(const gsl_vector * der, gsl_deriv_cand_t * cand)
{
  const size_t n = der->size;

  if (n > GSL_DERIV_MAX_STEPS)
    {
      GSL_ERROR ("der size cannot exceed GSL_DERIV_MAX_STEPS", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      cand->nkeep = 0;
      for (i = 0; i < n; ++i)
        {
          double Di = gsl_vector_get(der, i);

          if (gsl_finite(Di))
            cand->index[cand->nkeep++] = i;
        }

      return GSL_SUCCESS;
    }
}

static const gsl_deriv_cand_filter filter_mad_type =
{
  "MAD",
  filter_mad
};

const gsl_deriv_cand_filter * gsl_deriv_cand_filter_mad = &filter_mad_type;

static const gsl_deriv_cand_filter filter_trim_type =
{
  "trim",
  filter_trim
};

const gsl_deriv_cand_filter * gsl_deriv_cand_filter_trim = &filter_trim_type;

static const gsl_deriv_cand_filter filter_none_type =
{
  "none",
  filter_none
};

const gsl_deriv_cand_filter * gsl_deriv_cand_filter_none = &filter_none_type;

/* candidate selection methods */

/*
select_minerr()
  After filtering candidates, select the derivative candidate
with minimal error estimate

Inputs: cand      - candidates selected by candidate_filter()
        der       - derivative estimates
        err       - error estimates
        steps     - step sizes h
        der_est   - (output) estimate of derivative with lowest error
        err_est   - (output) error estimate
        step_est  - (output) step size corresponding to der_est/err_est

Return: success/error
*/

static int
select_minerr(const gsl_deriv_cand_t * cand, const gsl_vector * der,
              const gsl_vector * err, const gsl_vector * steps,
              double * der_est, double * err_est, double * step_est)
{
  const size_t n = der->size;

  if (cand->nkeep == 0)
    {
      GSL_ERROR ("no suitable candidates", GSL_EDOM);
    }
  else if (err->size != n)
    {
      GSL_ERROR ("err has different size than der", GSL_EBADLEN);
    }
  else if (steps->size != n)
    {
      GSL_ERROR ("steps has different size than der", GSL_EBADLEN);
    }
  else
    {
      size_t minidx = cand->index[0];
      double minerr = gsl_vector_get(err, minidx);
      size_t i;

      for (i = 1; i < cand->nkeep; ++i)
        {
          size_t idx = cand->index[i];
          double ei = gsl_vector_get(err, idx);

          if (ei < minerr)
            {
              minidx = idx;
              minerr = ei;
            }
        }

      *der_est = gsl_vector_get(der, minidx);
      *err_est = gsl_vector_get(err, minidx);
      *step_est = gsl_vector_get(steps, minidx);

      return GSL_SUCCESS;
    }
}

static const gsl_deriv_cand_select select_minerr_type =
{
  "minerr",
  select_minerr
};

const gsl_deriv_cand_select * gsl_deriv_cand_select_minerr = &select_minerr_type;
