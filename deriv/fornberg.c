/* deriv/fornberg.c
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv_scale.h>

/*
gsl_deriv_fornberg_weights()
  Generate finite difference weights using Fornberg algorithm
for arbitrary nodes,

[1] Fornberg, B. (1998). Calculation of weights in finite difference formulas.
    SIAM review, 40(3), 685-691.

Inputs: nodes     - node locations, length n
        n         - length of nodes[]
        max_deriv - maximum derivative order
        x0        - point for finite difference approximation
        weights   - (output) finite difference weights, length n*(max_deriv+1)

Return: success/error

Notes:
1) On output, weights[j + i*n] contains the jth weight for derivative order i.
In terms of gsl_matrix, one can do:

gsl_matrix_view W = gsl_matrix_view_array(weights, max_deriv+1, n)

and then row i will contain the weights for derivative order i
*/

int
gsl_deriv_fornberg_weights(const double nodes[], const unsigned int n,
                           const unsigned int max_deriv, const double x0,
                           double weights[])
{
  if (n < max_deriv + 1)
    {
      GSL_ERROR ("length of nodes[] insufficient", GSL_EBADLEN);
    }
  else
    {
      double c1, c4, c5;
      unsigned int i;

      c1 = 1;
      c4 = nodes[0] - x0;

      weights[0] = 1.0;
      for (i = 1; i < n*(max_deriv+1); ++i)
        weights[i] = 0.0; /* clear weights */

      for (i = 1; i < n; ++i)
        {
          const unsigned int mn = (int) GSL_MIN(i, max_deriv);
          double c2 = 1.0;
          unsigned int j;

          c5 = c4;
          c4 = nodes[i] - x0;

          for (j = 0; j < i; ++j)
            {
              const double c3 = nodes[i] - nodes[j];
              const double c3_r = 1.0 / c3;
              unsigned int k;

              c2 = c2 * c3;
              if (j == i - 1)
                {
                  const double c2_r = 1.0 / c2;

                  for (k = mn; k >= 1; --k)
                    {
                      const double tmp1 = weights[i - 1 + (k-1)*n];
                      const double tmp2 = weights[i - 1 + k*n];
                      weights[i + k*n] = c1*(k*tmp1 - c5*tmp2)*c2_r;
                    }

                  weights[i] = -c1*c5*weights[i-1]*c2_r;
                }

              for (k = mn; k >= 1; --k)
                {
                  const double tmp1 = weights[j + k*n];
                  const double tmp2 = weights[j + (k-1)*n];
                  weights[j + k*n] = (c4*tmp1 - k*tmp2)*c3_r;
                }

              weights[j] = c4*weights[j]*c3_r;
            }

          c1 = c2;
        }

      return GSL_SUCCESS;
    }
}
