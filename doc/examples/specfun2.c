#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_alf.h>

int
main(void)
{
  const size_t lmax = 5;
  const size_t nlm = gsl_sf_alf_nlm(lmax, lmax);
  const size_t idx21 = gsl_sf_alf_array_index(2, 1, lmax);
  double * Plm = malloc(gsl_sf_alf_array_size(lmax, lmax) * sizeof(double));
  double * dPlm = malloc(nlm * sizeof(double));
  double * d2Plm = malloc(nlm * sizeof(double));
  double x;

  /* spherical harmonic normalization, without the Condon-Shortley phase */
  gsl_sf_alf_precompute(GSL_SF_ALF_ORTHO_COMPLEX, lmax, lmax, 0, Plm);

  for (x = -0.99; x <= 0.99; x += 0.01)
    {
      const double sint = sqrt((1.0 - x) * (1.0 + x));

      /* P(l,m)/sin(theta) for m > 0, and the first two theta derivatives */
      gsl_sf_alf_vsh2_array(lmax, lmax, x, Plm, dPlm, d2Plm);
      printf("%f %e %e %e\n", x, Plm[idx21] * sint, dPlm[idx21], d2Plm[idx21]);
    }

  free(Plm);
  free(dPlm);
  free(d2Plm);

  return 0;
}
