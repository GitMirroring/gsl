/* bspline/augment.c
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

/* Build the augmented knot vector from the breakpoints tau (length nbreak)
 * without checking their order; see gsl_bspline_init_augment().  Used
 * directly by gsl_bspline_init_greville(), whose least squares breakpoints
 * can be out of order by rounding when breakpoints are repeated. */
static void
bspline_augment (const gsl_vector * tau, gsl_bspline_workspace * w)
{
  const double b1 = gsl_vector_get(tau, 0);
  const double bn = gsl_vector_get(tau, w->nbreak - 1);
  size_t i;

  for (i = 0; i < w->spline_order; i++)
    gsl_vector_set (w->knots, i, b1);

  for (i = 1; i < w->nbreak - 1; i++)
    {
      gsl_vector_set (w->knots, w->spline_order - 1 + i,
                      gsl_vector_get (tau, i));
    }

  for (i = w->ncontrol; i < w->ncontrol + w->spline_order; i++)
    gsl_vector_set (w->knots, i, bn);
}
