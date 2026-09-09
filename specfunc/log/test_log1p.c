/* specfunc/log/test_log1p.c
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>

#include "../test_sf_new.h"
#include "log1p_fallback.c"

static int
helper_log1p(int (*func)(const double x, gsl_sf_result * result))
{
  int s = 0;

  /*
   * Negative arguments, including points very close to the singular endpoint.
   */
  TEST_SF(s, func, (-0.9999999999999999, &r),
          -36.7368005696771013991133024372833581080015071211, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-0.999999999999, &r),
          -27.6310432378933585714557033162676284370216543357, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-0.99, &r),
          -4.60517018598809047985756320924389050674947901027, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-0.9, &r),
          -2.30258509299404590606259637971569694423072326342, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-0.75, &r),
          -1.38629436111989061883446424291635313615100026872, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-0.5, &r),
          -0.69314718055994530941723212145817656807550013436, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-0.25, &r),
          -0.287682072451780927439219005993827431503509710898, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-0.1, &r),
          -0.105360515657826307395406673201293597458937732868, TEST_TOL0, GSL_SUCCESS);

  /*
   * Fallback transition at |x| = 2^-29.  Exercise both sides and the boundary.
   */
  TEST_SF(s, func, (-1.8626451492309574e-09, &r),
          -0.00000000186264515096568092297123032283051621948714758097, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-1.862645149230957e-09, &r),
          -0.00000000186264515096568050938092327594469959395874465375, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-1.8626451492309568e-09, &r),
          -0.00000000186264515096568030258576975250179128119460733649, TEST_TOL0, GSL_SUCCESS);

  /*
   * Fallback transition at |x| = 2^-54 and very small negative arguments.
   */
  TEST_SF(s, func, (-5.551115123125784e-17, &r),
          -5.55111512312578408878771829926398988747874521805e-17, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-5.551115123125783e-17, &r),
          -5.55111512312578285619255389143297550887628360348e-17, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-5.551115123125782e-17, &r),
          -5.5511151231257822398949716875174683195750527962e-17, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (-1e-20, &r),
          -9.99999999999999945158271454209571651181069750663e-21, TEST_TOL0, GSL_SUCCESS);

  /*
   * Origin and very small positive arguments, including the 2^-54 transition.
   */
  TEST_SF(s, func, (0.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (1e-20, &r),
          9.99999999999999945148271454209571652278004321579e-21, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (5.551115123125782e-17, &r),
          5.55111512312578193174618058555980025339516034281e-17, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (5.551115123125783e-17, &r),
          5.55111512312578254804376278947523901991981278989e-17, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (5.551115123125784e-17, &r),
          5.55111512312578378063892719730611655296911768403e-17, TEST_TOL0, GSL_SUCCESS);

  /*
   * Positive side of the |x| = 2^-29 transition and the observed stress region.
   */
  TEST_SF(s, func, (1.8626451492309568e-09, &r),
          0.00000000186264514749623335063215632803138899709498293613, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (1.862645149230957e-09, &r),
          0.00000000186264514749623355742730908110231955496481912246, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (1.8626451492309574e-09, &r),
          0.00000000186264514749623397101761458724418067070436320243, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (1e-08, &r),
          0.00000000999999995000000054255893704236201799631283904207, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (1.0458344267283384e-08, &r),
          0.0000000104583442125949024602757830239163138757117029219, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (1e-06, &r),
          0.000000999999500000333287831490611062347623241482110589, TEST_TOL0, GSL_SUCCESS);

  /*
   * Representative moderate positive arguments.
   */
  TEST_SF(s, func, (0.1, &r),
          0.0953101798043248650904204170314766268673281912604, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (0.25, &r),
          0.223143551314209755766295090309834503374601085548, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (0.5, &r),
          0.405465108108164381978013115464349136571990423462, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (1.0, &r),
          0.69314718055994530941723212145817656807550013436, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (2.0, &r),
          1.09861228866810969139524523692252570464749055782, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (10.0, &r),
          2.39789527279837054406194357796512929982170685394, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (50.0, &r),
          3.93182563272432577164477985479565224023569357041, TEST_TOL0, GSL_SUCCESS);

  /*
   * Large positive arguments spanning many orders of magnitude.
   */
  TEST_SF(s, func, (1000000.0, &r),
          13.8155115579637741044412818114397185787732757413, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (1e+100, &r),
          230.258509299404568417702036578196338680495256895, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, func, (1e+300, &r),
          690.775527898213705257902196660513681150659990441, TEST_TOL0, GSL_SUCCESS);

  return s;
}

int
test_log1p(void)
{
  int s = 0;

  /*
   * Apply the same numerical corpus to both the configured public
   * implementation and the portable C89 fallback.
   */
  s += helper_log1p(gsl_sf_log1p_e);
  s += helper_log1p(log1p_fallback);

  /*
   * API/domain behavior belongs to the public function.  The fallback
   * kernel is called directly only on its valid domain x > -1.
   */
  TEST_SF(s, gsl_sf_log1p_e, (-1.0, &r),
          GSL_NEGINF, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_log1p_e, (-1.0000000000000002, &r),
          GSL_NAN, TEST_TOL0, GSL_EDOM);

  return s;
}
