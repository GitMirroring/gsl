/* specfunc/test_gamma.c
 * 
 * Copyright (C) 2007 Brian Gough
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
#include <gsl/gsl_test.h>
#include <gsl/gsl_sf.h>
#include "test_sf.h"

/* functions not (yet) in the gamma/ subdirectory */
int test_gamma_misc(void)
{
  gsl_sf_result r;
  gsl_sf_result r1, r2;
  double sgn;
  int s = 0;

  TEST_SF_2(s, gsl_sf_lngamma_complex_e, (5.0, 2.0, &r1, &r2),
            2.7487017561338026749, TEST_TOL0,
            3.0738434100497007915, TEST_TOL0,
            GSL_SUCCESS);
            
  TEST_SF_2(s, gsl_sf_lngamma_complex_e, (100.0, 100.0, &r1, &r2),
            315.07804459949331323, TEST_TOL1,
            2.0821801804113110099, TEST_TOL3,
            GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_lngamma_complex_e, (100.0, -1000.0, &r1, &r2),
            -882.3920483010362817000, TEST_TOL1,
            -2.1169293725678813270, TEST_TOL3,
            GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_lngamma_complex_e, (-100.0, -1.0, &r1, &r2),
            -365.0362469529239516000, TEST_TOL1,
            -3.0393820262864361140, TEST_TOL1,
            GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_lnchoose_e, (7,3, &r), 3.555348061489413680 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnchoose_e, (5,2, &r), 2.302585092994045684 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_choose_e, (7,3, &r), 35.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_choose_e, (7,4, &r), 35.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_choose_e, (5,2, &r), 10.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_choose_e, (5,3, &r), 10.0 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_choose_e, (500,495, &r), 255244687600.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_choose_e, (500,5, &r), 255244687600.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_choose_e, (500,200, &r), 5.054949849935532221e+144 , TEST_TOL5, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_choose_e, (500,300, &r), 5.054949849935532221e+144 , TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_lnpoch_e, (5, 0.0, &r), 0.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnpoch_e, (5, 1.0/65536.0, &r), 0.000022981557571259389129, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnpoch_e, (5, 1.0/256.0, &r),   0.005884960217985189004,    TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnpoch_e, (7,3, &r), 6.222576268071368616, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_lnpoch_e, (5,2, &r), 3.401197381662155375, TEST_TOL0, GSL_SUCCESS);

  TEST_SF_SGN(s, gsl_sf_lnpoch_sgn_e, (5.0, 0.0, &r, &sgn), 0.0, TEST_TOL1, 1.0, GSL_SUCCESS);
  TEST_SF_SGN(s, gsl_sf_lnpoch_sgn_e, (-4.5, 0.25, &r, &sgn), 0.7430116475119920117, TEST_TOL1, 1.0, GSL_SUCCESS);
  TEST_SF_SGN(s, gsl_sf_lnpoch_sgn_e, (-4.5, 1.25, &r, &sgn), 2.1899306304483174731, TEST_TOL1, -1.0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_poch_e, (5, 0.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (7,3, &r), 504.0 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (5,2, &r), 30.0  , TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (5,1.0/256.0, &r), 1.0059023106151364982 , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_pochrel_e, (5,0, &r), 1.506117668431800472, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (7,3, &r), 503.0/3.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (5,2, &r), 29.0/2.0, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (5,0.01, &r), 1.5186393661368275330, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (-5.5,0.01, &r), 1.8584945633829063516, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (-5.5,-1.0/8.0, &r), 1.0883319303552135488, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (-5.5,-1.0/256.0, &r), 1.7678268037726177453, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_pochrel_e, (-5.5,-11.0, &r), 0.09090909090939652475, TEST_TOL0, GSL_SUCCESS);

  /* Add tests for special cases with negative arguments */
  TEST_SF(s,  gsl_sf_poch_e, (-9.0, -4.0, &r), 1.0/17160.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (-9.0, -3.0, &r), -1.0/1320.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (-9.0, -3.5, &r), 0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_poch_e, (-9.0, 4.0, &r), 3024.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (-9.0, 3.0, &r), -504.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (-9.0, 3.5, &r), 0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_poch_e, (-9.0, 0.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_poch_e, (-8.0, -4.0, &r), 1.0/11880.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (-8.0, -3.0, &r), -1.0/990.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (-8.0, +4.0, &r), 1680.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (-8.0, +3.0, &r), -336.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (-3.0, +4.0, &r), 0.0, TEST_TOL0, GSL_SUCCESS);

  /* FIXME: we should be able to get an exact answer for poch(-a,a) if
     gsl_sf_lngamma functions were fixed to handle integer arguments
     exactly as a special case */
  TEST_SF(s,  gsl_sf_poch_e, (-3.0, +3.0, &r), -6.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_poch_e, (-4.0, +4.0, &r), 24.0, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_poch_e, (-3.0, +100.0, &r), 0.0, TEST_TOL0, GSL_SUCCESS);

  

  TEST_SF(s,  gsl_sf_beta_inc_e, (1.0, 1.0, 0.0, &r), 0.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (1.0, 1.0, 1.0, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (0.1, 0.1, 1.0, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, ( 1.0,  1.0, 0.5, &r), 0.5, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, ( 0.1,  1.0, 0.5, &r), 0.9330329915368074160, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (10.0,  1.0, 0.5, &r), 0.0009765625000000000000, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (50.0,  1.0, 0.5, &r), 8.881784197001252323e-16, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, ( 1.0,  0.1, 0.5, &r), 0.06696700846319258402, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, ( 1.0, 10.0, 0.5, &r), 0.99902343750000000000, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, ( 1.0, 50.0, 0.5, &r), 0.99999999999999911180, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, ( 1.0,  1.0, 0.1, &r), 0.10, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, ( 1.0,  2.0, 0.1, &r), 0.19, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, ( 1.0,  2.0, 0.9, &r), 0.99, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (50.0, 60.0, 0.5, &r), 0.8309072939016694143, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (90.0, 90.0, 0.5, &r), 0.5, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, ( 500.0,  500.0, 0.6, &r), 0.9999999999157549630, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (5000.0, 5000.0, 0.4, &r), 4.518543727260666383e-91, TEST_TOL5, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (5000.0, 5000.0, 0.6, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (5000.0, 2000.0, 0.6, &r), 8.445388773903332659e-89, TEST_TOL5, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_beta_inc_e, (-0.1, -0.1, 1.0, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (-0.1, -0.2, 1.0, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (-0.2, -0.1, 1.0, &r), 1.0, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_beta_inc_e, (-0.1, -0.2, 0.5, &r), 0.675252001958389971991335, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (-0.2, -0.1, 0.5, &r), 0.324747998041610028008665, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_beta_inc_e, (-0.1, -0.1, 0.0, &r), 0.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (-0.1, -0.2, 0.0, &r), 0.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (-0.2, -0.1, 0.0, &r), 0.0, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_beta_inc_e, (-0.1, -0.2, 0.3, &r), 0.7469186777964287252, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_beta_inc_e, (-0.2, -0.1, 0.3, &r), 0.3995299653262016818, TEST_TOL2, GSL_SUCCESS);

  /* Bug report from Thomas Tanner <tanner@gmx.de> */

  TEST_SF(s,  gsl_sf_beta_inc_e, (0.5, 101.5, 0.999457, &r), 1.0, TEST_TOL2, GSL_SUCCESS);

  return s;
}
