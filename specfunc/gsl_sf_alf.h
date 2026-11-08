/* specfunc/gsl_sf_alf.h
 * 
 * Copyright (C) 2023 Patrick Alken
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

#ifndef __GSL_SF_ALF_H__
#define __GSL_SF_ALF_H__

#include <stdlib.h>
#include <gsl/gsl_inline.h>
#include <gsl/gsl_sf_result.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/*
 * Normalizations: P(l,m) = K(l,m) P_l^m(x), with the unnormalized P_l^m
 * (without the Condon-Shortley phase, see GSL_SF_ALF_FLG_CSPHASE) and
 *
 *   SCHMIDT     K = sqrt(c_m (l-m)!/(l+m)!)
 *   FOURPI      K = sqrt(c_m (2l+1) (l-m)!/(l+m)!)
 *   ORTHO       K = sqrt(c_m (2l+1)/(4 pi) (l-m)!/(l+m)!)
 *
 * where c_m = 2 - delta_{m0} for the _REAL variants (for real spherical
 * harmonics P(l,m) cos(m phi), P(l,m) sin(m phi)) and c_m = 1 for the
 * _COMPLEX variants (for P(l,m) exp(i m phi)); the spherical harmonics then
 * have squared norm 4 pi/(2l+1), 4 pi and 1 respectively over the sphere.
 *
 *   FULL        K = sqrt((l + 1/2) (l-m)!/(l+m)!)   (int_{-1}^1 P^2 dx = 1)
 *   NONE        K = 1
 *
 * The values of the first five are those of GSL 2.8 (where the interface
 * was present but undocumented) under their previous names.
 */
typedef enum
{
  GSL_SF_ALF_SCHMIDT_REAL = 0,
  GSL_SF_ALF_ORTHO_COMPLEX = 1,
  GSL_SF_ALF_FULL = 2,
  GSL_SF_ALF_FOURPI_REAL = 3,
  GSL_SF_ALF_NONE = 4,
  GSL_SF_ALF_SCHMIDT_COMPLEX = 5,
  GSL_SF_ALF_FOURPI_COMPLEX = 6,
  GSL_SF_ALF_ORTHO_REAL = 7
} gsl_sf_alf_t;

#ifndef GSL_DISABLE_DEPRECATED
/* previous names (GSL 2.8, undocumented) */
#define GSL_SF_ALF_SCHMIDT  GSL_SF_ALF_SCHMIDT_REAL
#define GSL_SF_ALF_SPHARM   GSL_SF_ALF_ORTHO_COMPLEX
#define GSL_SF_ALF_FOURPI   GSL_SF_ALF_FOURPI_REAL
#endif /* !GSL_DISABLE_DEPRECATED */

#define GSL_SF_ALF_FLG_CSPHASE       (1 << 0) /* include Condon-Shortley phase */

int gsl_sf_alf_precompute(const gsl_sf_alf_t norm, const size_t lmax,
                          const size_t mmax, const size_t flags, double result_array[]);
size_t gsl_sf_alf_array_size(const size_t lmax, const size_t mmax);
int gsl_sf_alf_array(const size_t lmax, const size_t mmax, const double x, double result_array[]);
int gsl_sf_alf_deriv_array(const size_t lmax, const size_t mmax, const double x,
                           double result_array[], double result_deriv_array[]);
int gsl_sf_alf_deriv2_array(const size_t lmax, const size_t mmax, const double x,
                            double result_array[], double result_deriv_array[],
                            double result_deriv2_array[]);
int gsl_sf_alf_vsh_array(const size_t lmax, const size_t mmax, const double x,
                         double result_array[], double result_deriv_array[]);
int gsl_sf_alf_vsh2_array(const size_t lmax, const size_t mmax, const double x,
                          double result_array[], double result_deriv_array[],
                          double result_deriv2_array[]);
int gsl_sf_alf_theta_derivk_array(const size_t lmax, const double Plm[],
                                  const double input_array[], double output_array[]);

INLINE_DECL size_t gsl_sf_alf_nlm(const size_t lmax, const size_t mmax);
INLINE_DECL size_t gsl_sf_alf_array_index(const size_t l, const size_t m, const size_t lmax);

#ifdef HAVE_INLINE

/*
gsl_sf_alf_array_index()
This routine computes the index into a result_array[] corresponding
to a given (l,m) using 'M' major indexing:

(l,m) = (0,0) (1,0) (2,0) ... (L,0) (1,1) (2,1) ... (L,1) ... (L,L)

index(l,m,L) = m*L - m(m-1)/2 + l
*/

INLINE_FUN
size_t
gsl_sf_alf_array_index(const size_t l, const size_t m, const size_t lmax)
{
  return (((m * (2 * lmax + 1 - m)) >> 1) + l);
}

/* return number of ALFs for a given lmax = (mmax+1)*(lmax+1) - mmax*(mmax+1)/2 */
INLINE_FUN
size_t
gsl_sf_alf_nlm(const size_t lmax, const size_t mmax)
{
  const size_t M = (mmax > lmax) ? lmax : mmax;
  return ((lmax + 1) * (M + 1) - ((M * (M + 1)) >> 1));
}

#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_SF_ALF_H__ */
