/* sys/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>

#include "log1p_fallback.c"

int
main (void)
{
  double y, y_expected;
  int e, e_expected;

  gsl_ieee_env_setup ();

  /* Test for gsl_hypot */

  y = gsl_hypot (0.0, 0.0);
  y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(0.0, 0.0)");

  y = gsl_hypot (1e-10, 1e-10);
  y_expected = 1.414213562373095048801688e-10;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e-10, 1e-10)");

  y = gsl_hypot (1e-38, 1e-38);
  y_expected = 1.414213562373095048801688e-38;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e-38, 1e-38)");

  y = gsl_hypot (1e-10, -1.0);
  y_expected = 1.000000000000000000005;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e-10, -1)");

  y = gsl_hypot (-1.0, 1e-10);
  y_expected = 1.000000000000000000005;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(-1, 1e-10)");

  y = gsl_hypot (1e307, 1e301);
  y_expected = 1.000000000000499999999999e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e307, 1e301)");

  y = gsl_hypot (1e301, 1e307);
  y_expected = 1.000000000000499999999999e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e301, 1e307)");

  y = gsl_hypot (1e307, 1e307);
  y_expected = 1.414213562373095048801688e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e307, 1e307)");

  /* Test +-Inf, finite */
  
  y = gsl_hypot (GSL_POSINF, 1.2);
  y_expected = GSL_POSINF;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(GSL_POSINF, 1.2)");

  y = gsl_hypot (GSL_NEGINF, 1.2);
  y_expected = GSL_POSINF;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(GSL_NEGINF, 1.2)");

  y = gsl_hypot (1.2, GSL_POSINF);
  y_expected = GSL_POSINF;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1.2, GSL_POSINF)");

  y = gsl_hypot (1.2, GSL_NEGINF);
  y_expected = GSL_POSINF;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1.2, GSL_NEGINF)");

  /* Test NaN, finite */
  
  y = gsl_hypot (GSL_NAN, 1.2);
  y_expected = GSL_NAN;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(GSL_NAN, 1.2)");

  y = gsl_hypot (1.2, GSL_NAN);
  y_expected = GSL_NAN;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1.2, GSL_NAN)");

  /* Test NaN, NaN */

  y = gsl_hypot (GSL_NAN, GSL_NAN);
  y_expected = GSL_NAN;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(GSL_NAN, GSL_NAN)");

  /* Test +Inf, NaN */

  y = gsl_hypot (GSL_POSINF, GSL_NAN);
  y_expected = GSL_POSINF;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(GSL_POSINF, GSL_NAN)");

  /* Test -Inf, NaN */

  y = gsl_hypot (GSL_NEGINF, GSL_NAN);
  y_expected = GSL_POSINF;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(GSL_NEGINF, GSL_NAN)");

  /* Test NaN, +Inf */

  y = gsl_hypot (GSL_NAN, GSL_POSINF);
  y_expected = GSL_POSINF;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(GSL_NAN, GSL_POSINF)");

  /* Test NaN, -Inf */

  y = gsl_hypot (GSL_NAN, GSL_NEGINF);
  y_expected = GSL_POSINF;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(GSL_NAN, GSL_NEGINF)");

  /* Test for gsl_hypot3 */

  y = gsl_hypot3 (0.0, 0.0, 0.0);
  y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot3(0.0, 0.0, 0.0)");

  y = gsl_hypot3 (1e-10, 1e-10, 1e-10);
  y_expected = 1.732050807568877293527446e-10;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot3(1e-10, 1e-10, 1e-10)");

  y = gsl_hypot3 (1e-38, 1e-38, 1e-38);
  y_expected = 1.732050807568877293527446e-38;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot3(1e-38, 1e-38, 1e-38)");

  y = gsl_hypot3 (1e-10, 1e-10, -1.0);
  y_expected = 1.000000000000000000099;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot3(1e-10, 1e-10, -1)");

  y = gsl_hypot3 (1e-10, -1.0, 1e-10);
  y_expected = 1.000000000000000000099;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot3(1e-10, -1, 1e-10)");

  y = gsl_hypot3 (-1.0, 1e-10, 1e-10);
  y_expected = 1.000000000000000000099;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot3(-1, 1e-10, 1e-10)");

  y = gsl_hypot3 (1e307, 1e301, 1e301);
  y_expected = 1.0000000000009999999999995e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot3(1e307, 1e301, 1e301)");

  y = gsl_hypot3 (1e307, 1e307, 1e307);
  y_expected = 1.732050807568877293527446e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot3(1e307, 1e307, 1e307)");

  y = gsl_hypot3 (1e307, 1e-307, 1e-307);
  y_expected = 1.0e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot3(1e307, 1e-307, 1e-307)");

  /* Test for acosh */

  y = gsl_acosh (1.0);
  y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_acosh(1.0)");

  y = gsl_acosh (1.1);
  y_expected = 4.435682543851151891329110663525e-1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_acosh(1.1)");

  y = gsl_acosh (10.0);
  y_expected = 2.9932228461263808979126677137742e0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_acosh(10.0)");

  y = gsl_acosh (1e10);
  y_expected = 2.3718998110500402149594646668302e1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_acosh(1e10)");

  /* Test for asinh */

  y = gsl_asinh (0.0);
  y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(0.0)");

  y = gsl_asinh (1e-10);
  y_expected = 9.9999999999999999999833333333346e-11;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(1e-10)");

  y = gsl_asinh (-1e-10);
  y_expected = -9.9999999999999999999833333333346e-11;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(1e-10)");

  y = gsl_asinh (0.1);
  y_expected = 9.983407889920756332730312470477e-2;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(0.1)");

  y = gsl_asinh (-0.1);
  y_expected = -9.983407889920756332730312470477e-2;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(-0.1)");

  y = gsl_asinh (1.0);
  y_expected = 8.8137358701954302523260932497979e-1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(1.0)");

  y = gsl_asinh (-1.0);
  y_expected = -8.8137358701954302523260932497979e-1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(-1.0)");

  y = gsl_asinh (10.0);
  y_expected = 2.9982229502979697388465955375965e0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(10)");

  y = gsl_asinh (-10.0);
  y_expected = -2.9982229502979697388465955375965e0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(-10)");

  y = gsl_asinh (1e10);
  y_expected = 2.3718998110500402149599646668302e1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(1e10)");

  y = gsl_asinh (-1e10);
  y_expected = -2.3718998110500402149599646668302e1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(-1e10)");

  /* Test for atanh */

  y = gsl_atanh (0.0);
  y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(0.0)");

  y = gsl_atanh (1e-20);
  y_expected = 1e-20;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(1e-20)");

  y = gsl_atanh (-1e-20);
  y_expected = -1e-20;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(-1e-20)");

  y = gsl_atanh (0.1);
  y_expected = 1.0033534773107558063572655206004e-1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(0.1)");

  y = gsl_atanh (-0.1);
  y_expected = -1.0033534773107558063572655206004e-1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(-0.1)");

  y = gsl_atanh (0.9);
  y_expected = 1.4722194895832202300045137159439e0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(0.9)");

  y = gsl_atanh (-0.9);
  y_expected = -1.4722194895832202300045137159439e0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(0.9)");

  /* Test for log1p: gsl_log1p and the C89 fallback (tested directly,
     since gsl_log1p uses the system log1p when available) */

  {
    static const double log1p_tab[][2] = {
      { -0.9999999999999999, -36.7368005696771013991133024372833581080015071211 },
      { -0.999999999999, -27.6310432378933585714557033162676284370216543357 },
      { -0.99, -4.60517018598809047985756320924389050674947901027 },
      { -0.9, -2.30258509299404590606259637971569694423072326342 },
      { -0.75, -1.38629436111989061883446424291635313615100026872 },
      { -0.5, -0.69314718055994530941723212145817656807550013436 },
      { -0.25, -0.287682072451780927439219005993827431503509710898 },
      { -0.1, -0.105360515657826307395406673201293597458937732868 },
      { -1.8626451492309574e-09, -0.00000000186264515096568092297123032283051621948714758097 },
      { -1.862645149230957e-09, -0.00000000186264515096568050938092327594469959395874465375 },
      { -1.8626451492309568e-09, -0.00000000186264515096568030258576975250179128119460733649 },
      { -5.551115123125784e-17, -5.55111512312578408878771829926398988747874521805e-17 },
      { -5.551115123125783e-17, -5.55111512312578285619255389143297550887628360348e-17 },
      { -5.551115123125782e-17, -5.5511151231257822398949716875174683195750527962e-17 },
      { -1e-20, -9.99999999999999945158271454209571651181069750663e-21 },
      { 0.0, 0.0 },
      { 1e-20, 9.99999999999999945148271454209571652278004321579e-21 },
      { 5.551115123125782e-17, 5.55111512312578193174618058555980025339516034281e-17 },
      { 5.551115123125783e-17, 5.55111512312578254804376278947523901991981278989e-17 },
      { 5.551115123125784e-17, 5.55111512312578378063892719730611655296911768403e-17 },
      { 1.8626451492309568e-09, 0.00000000186264514749623335063215632803138899709498293613 },
      { 1.862645149230957e-09, 0.00000000186264514749623355742730908110231955496481912246 },
      { 1.8626451492309574e-09, 0.00000000186264514749623397101761458724418067070436320243 },
      { 1e-08, 0.00000000999999995000000054255893704236201799631283904207 },
      { 1.0458344267283384e-08, 0.0000000104583442125949024602757830239163138757117029219 },
      { 1e-06, 0.000000999999500000333287831490611062347623241482110589 },
      { 0.1, 0.0953101798043248650904204170314766268673281912604 },
      { 0.25, 0.223143551314209755766295090309834503374601085548 },
      { 0.5, 0.405465108108164381978013115464349136571990423462 },
      { 1.0, 0.69314718055994530941723212145817656807550013436 },
      { 2.0, 1.09861228866810969139524523692252570464749055782 },
      { 10.0, 2.39789527279837054406194357796512929982170685394 },
      { 50.0, 3.93182563272432577164477985479565224023569357041 },
      { 1000000.0, 13.8155115579637741044412818114397185787732757413 },
      { 1e+100, 230.258509299404568417702036578196338680495256895 },
      { 1e+300, 690.775527898213705257902196660513681150659990441 }
    };
    size_t i;

    for (i = 0; i < sizeof (log1p_tab) / sizeof (log1p_tab[0]); i++)
      {
        const double x = log1p_tab[i][0];
        const double y_ref = log1p_tab[i][1];

        gsl_test_rel (gsl_log1p (x), y_ref, 2 * GSL_DBL_EPSILON,
                      "gsl_log1p(%.17g)", x);
        gsl_test_rel (log1p_fallback (x), y_ref, 2 * GSL_DBL_EPSILON,
                      "log1p_fallback(%.17g)", x);
      }

    gsl_test_rel (gsl_log1p (-1.0), GSL_NEGINF, 0.0, "gsl_log1p(-1)");
    gsl_test_rel (gsl_log1p (GSL_POSINF), GSL_POSINF, 0.0, "gsl_log1p(+Inf)");
    gsl_test (!gsl_isnan (gsl_log1p (-1.0000000000000002)),
              "gsl_log1p(-1-eps) is NaN");
    gsl_test (!gsl_isnan (gsl_log1p (GSL_NEGINF)), "gsl_log1p(-Inf) is NaN");
    gsl_test (!gsl_isnan (gsl_log1p (GSL_NAN)), "gsl_log1p(NaN) is NaN");
  }

  /* Test for pow_int */

  y = gsl_pow_2 (-3.14);
  y_expected = pow (-3.14, 2.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_2(-3.14)");

  y = gsl_pow_3 (-3.14);
  y_expected = pow (-3.14, 3.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_3(-3.14)");

  y = gsl_pow_4 (-3.14);
  y_expected = pow (-3.14, 4.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_4(-3.14)");

  y = gsl_pow_5 (-3.14);
  y_expected = pow (-3.14, 5.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_5(-3.14)");

  y = gsl_pow_6 (-3.14);
  y_expected = pow (-3.14, 6.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_6(-3.14)");

  y = gsl_pow_7 (-3.14);
  y_expected = pow (-3.14, 7.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_7(-3.14)");

  y = gsl_pow_8 (-3.14);
  y_expected = pow (-3.14, 8.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_8(-3.14)");

  y = gsl_pow_9 (-3.14);
  y_expected = pow (-3.14, 9.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_9(-3.14)");

  {
    int n;
    for (n = -9; n < 10; n++)
      {
        y = gsl_pow_int (-3.14, n);
        y_expected = pow (-3.14, n);
        gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_int(-3.14,%d)", n);
      }
  }


  {
    unsigned int n;
    for (n = 0; n < 10; n++)
      {
        y = gsl_pow_uint (-3.14, n);
        y_expected = pow (-3.14, n);
        gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_uint(-3.14,%d)", n);
      }
  }

  /* Test case for n at INT_MAX, INT_MIN */

  {
    double u = 1.0000001;
    int n = INT_MAX;
    y = gsl_pow_int (u, n);
    y_expected = pow (u, n);
    gsl_test_rel (y, y_expected, 1e-6, "gsl_pow_int(%.7f,%d)", u, n);

    n = INT_MIN;
    y = gsl_pow_int (u, n);
    y_expected = pow (u, n);
    gsl_test_rel (y, y_expected, 1e-6, "gsl_pow_int(%.7f,%d)", u, n);

    /* exact results; INT_MIN is even, INT_MAX odd (bug #68495: -INT_MIN overflowed) */

    gsl_test_rel (gsl_pow_int (-1.0, INT_MIN), 1.0, 0.0, "gsl_pow_int(-1,INT_MIN)");
    gsl_test_rel (gsl_pow_int (-1.0, INT_MAX), -1.0, 0.0, "gsl_pow_int(-1,INT_MAX)");
    gsl_test_rel (gsl_pow_int (2.0, INT_MIN), 0.0, 0.0, "gsl_pow_int(2,INT_MIN)");
  }

  /* Test for ldexp */

  y = gsl_ldexp (M_PI, -2);
  y_expected = M_PI_4;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_ldexp(pi,-2)");

  y = gsl_ldexp (1.0, 2);
  y_expected = 4.000000;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_ldexp(1.0,2)");

  y = gsl_ldexp (0.0, 2);
  y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_ldexp(0.0,2)");

  y = gsl_ldexp (9.999999999999998890e-01, 1024);
  y_expected = GSL_DBL_MAX;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_ldexp DBL_MAX");

  y = gsl_ldexp (1e308, -2000);
  y_expected = 8.7098098162172166755761e-295;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_ldexp(1e308,-2000)");

  y = gsl_ldexp (GSL_DBL_MIN, 2000);
  y_expected = 2.554675596204441378334779940e294;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_ldexp(DBL_MIN,2000)");

  /* Test subnormals */

  {
    int i = 0;
    volatile double x = GSL_DBL_MIN;
    y_expected = 2.554675596204441378334779940e294;
    
    x /= 2;
    while (x > 0)
      {
        i++ ;
        y = gsl_ldexp (x, 2000 + i);
        gsl_test_rel (y, y_expected, 1e-15, "gsl_ldexp(DBL_MIN/2**%d,%d)",i,2000+i);
        x /= 2;
      }
  }


  /* Test for frexp */

  y = gsl_frexp (0.0, &e);
  y_expected = 0;
  e_expected = 0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_frexp(0) fraction");
  gsl_test_int (e, e_expected, "gsl_frexp(0) exponent");

  y = gsl_frexp (M_PI, &e);
  y_expected = M_PI_4;
  e_expected = 2;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_frexp(pi) fraction");
  gsl_test_int (e, e_expected, "gsl_frexp(pi) exponent");

  y = gsl_frexp (2.0, &e);
  y_expected = 0.5;
  e_expected = 2;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_frexp(2.0) fraction");
  gsl_test_int (e, e_expected, "gsl_frexp(2.0) exponent");

  y = gsl_frexp (1.0 / 4.0, &e);
  y_expected = 0.5;
  e_expected = -1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_frexp(0.25) fraction");
  gsl_test_int (e, e_expected, "gsl_frexp(0.25) exponent");

  y = gsl_frexp (1.0 / 4.0 - 4.0 * GSL_DBL_EPSILON, &e);
  y_expected = 0.999999999999996447;
  e_expected = -2;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_frexp(0.25-eps) fraction");
  gsl_test_int (e, e_expected, "gsl_frexp(0.25-eps) exponent");

  y = gsl_frexp (GSL_DBL_MAX, &e);
  y_expected = 9.999999999999998890e-01;
  e_expected = 1024;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_frexp(DBL_MAX) fraction");
  gsl_test_int (e, e_expected, "gsl_frexp(DBL_MAX) exponent");

  y = gsl_frexp (-GSL_DBL_MAX, &e);
  y_expected = -9.999999999999998890e-01;
  e_expected = 1024;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_frexp(-DBL_MAX) fraction");
  gsl_test_int (e, e_expected, "gsl_frexp(-DBL_MAX) exponent");

  y = gsl_frexp (GSL_DBL_MIN, &e);
  y_expected = 0.5;
  e_expected = -1021;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_frexp(DBL_MIN) fraction");
  gsl_test_int (e, e_expected, "gsl_frexp(DBL_MIN) exponent");

  y = gsl_frexp (-GSL_DBL_MIN, &e);
  y_expected = -0.5;
  e_expected = -1021;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_frexp(-DBL_MIN) fraction");
  gsl_test_int (e, e_expected, "gsl_frexp(-DBL_MIN) exponent");

  /* Test subnormals */

  {
    int i = 0;
    volatile double x = GSL_DBL_MIN;
    y_expected = 0.5;
    e_expected = -1021;
    
    x /= 2;

    while (x > 0)
      {
        e_expected--;
        i++ ;
        
        y = gsl_frexp (x, &e);
        gsl_test_rel (y, y_expected, 1e-15, "gsl_frexp(DBL_MIN/2**%d) fraction",i);
        gsl_test_int (e, e_expected, "gsl_frexp(DBL_MIN/2**%d) exponent", i);
        x /= 2;
      }
  }


  /* Test for approximate floating point comparison */
  {
    double x, y;
    int i;

    x = M_PI;
    y = 22.0 / 7.0;

    /* test the basic function */

    for (i = 0; i < 10; i++)
      {
        double tol = pow (10, -i);
        int res = gsl_fcmp (x, y, tol);
        gsl_test_int (res, -(i >= 4), "gsl_fcmp(%.5f,%.5f,%g)", x, y, tol);
      }

    for (i = 0; i < 10; i++)
      {
        double tol = pow (10, -i);
        int res = gsl_fcmp (y, x, tol);
        gsl_test_int (res, (i >= 4), "gsl_fcmp(%.5f,%.5f,%g)", y, x, tol);
      }
  }
    

#if HAVE_IEEE_COMPARISONS
  /* Test for isinf, isnan, finite */

  {
    double zero, one, inf, nan;
    int s;

    zero = 0.0;
    one = 1.0;
    inf = exp (1.0e10);
    nan = inf / inf;

    s = gsl_isinf (zero);
    gsl_test_int (s, 0, "gsl_isinf(0)");

    s = gsl_isinf (one);
    gsl_test_int (s, 0, "gsl_isinf(1)");

    s = gsl_isinf (inf);
    gsl_test_int (s, 1, "gsl_isinf(inf)");

    s = gsl_isinf (-inf);  
    gsl_test_int (s, -1, "gsl_isinf(-inf)");

    s = gsl_isinf (nan);
    gsl_test_int (s, 0, "gsl_isinf(nan)");


    s = gsl_isnan (zero);
    gsl_test_int (s, 0, "gsl_isnan(0)");

    s = gsl_isnan (one);
    gsl_test_int (s, 0, "gsl_isnan(1)");

    s = gsl_isnan (inf);
    gsl_test_int (s, 0, "gsl_isnan(inf)");

    s = gsl_isnan (-inf);
    gsl_test_int (s, 0, "gsl_isnan(-inf)");

    s = gsl_isnan (nan) != 0;
    gsl_test_int (s, 1, "gsl_isnan(nan)");


    s = gsl_finite (zero);
    gsl_test_int (s, 1, "gsl_finite(0)");

    s = gsl_finite (one);
    gsl_test_int (s, 1, "gsl_finite(1)");

    s = gsl_finite (inf);
    gsl_test_int (s, 0, "gsl_finite(inf)");

    s = gsl_finite (-inf);
    gsl_test_int (s, 0, "gsl_finite(-inf)");

    s = gsl_finite (nan);
    gsl_test_int (s, 0, "gsl_finite(nan)");
  }
#endif


  {
    double x = gsl_fdiv (2.0, 3.0);
    gsl_test_rel (x, 2.0 / 3.0, 4 * GSL_DBL_EPSILON, "gsl_fdiv(2,3)");
  }


  /* Test constants in gsl_math.h */

  {
    double x = log(M_E);
    gsl_test_rel (x, 1.0, 4 * GSL_DBL_EPSILON, "ln(M_E)");
  }
  
  {
    double x=pow(2.0,M_LOG2E);
    gsl_test_rel (x, exp(1.0), 4 * GSL_DBL_EPSILON, "2^M_LOG2E");
  }
 
  {
    double x=pow(10.0,M_LOG10E);
    gsl_test_rel (x, exp(1.0), 4 * GSL_DBL_EPSILON, "10^M_LOG10E");
  }

  {
    double x=pow(M_SQRT2, 2.0);
    gsl_test_rel (x, 2.0, 4 * GSL_DBL_EPSILON, "M_SQRT2^2");
  }    

  {
    double x=pow(M_SQRT1_2, 2.0);
    gsl_test_rel (x, 1.0/2.0, 4 * GSL_DBL_EPSILON, "M_SQRT1_2");
  }    

  {
    double x=pow(M_SQRT3, 2.0);
    gsl_test_rel (x, 3.0, 4 * GSL_DBL_EPSILON, "M_SQRT3^2");
  }    

  {
    double x = M_PI;
    gsl_test_rel (x, 3.1415926535897932384626433832795, 4 * GSL_DBL_EPSILON, "M_PI");
  }    

  {
    double x = 2 * M_PI_2;
    gsl_test_rel (x, M_PI, 4 * GSL_DBL_EPSILON, "2*M_PI_2");
  }    

  {
    double x = 4 * M_PI_4;
    gsl_test_rel (x, M_PI, 4 * GSL_DBL_EPSILON, "4*M_PI_4");
  }    

  {
    double x = pow(M_SQRTPI, 2.0);
    gsl_test_rel (x, M_PI, 4 * GSL_DBL_EPSILON, "M_SQRTPI^2");
  }    

  {
    double x = pow(M_2_SQRTPI, 2.0);
    gsl_test_rel (x, 4/M_PI, 4 * GSL_DBL_EPSILON, "M_SQRTPI^2");
  }    

  {
    double x = M_1_PI;
    gsl_test_rel (x, 1/M_PI, 4 * GSL_DBL_EPSILON, "M_1_SQRTPI");
  }    

  {
    double x = M_2_PI;
    gsl_test_rel (x, 2.0/M_PI, 4 * GSL_DBL_EPSILON, "M_2_PI");
  }    

  {
    double x = exp(M_LN10);
    gsl_test_rel (x, 10, 4 * GSL_DBL_EPSILON, "exp(M_LN10)");
  }    

  {
    double x = exp(M_LN2);
    gsl_test_rel (x, 2, 4 * GSL_DBL_EPSILON, "exp(M_LN2)");
  }    

  {
    double x = exp(M_LNPI);
    gsl_test_rel (x, M_PI, 4 * GSL_DBL_EPSILON, "exp(M_LNPI)");
  }    

  {
    double x = M_EULER;
    gsl_test_rel (x, 0.5772156649015328606065120900824, 4 * GSL_DBL_EPSILON, "M_EULER");
  }    

  exit (gsl_test_summary ());
}
