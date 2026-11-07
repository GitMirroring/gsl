/* specfunc/test_hyperg.c
 * 
 * Copyright (C) 2007, 2009, 2010 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
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

int test_hyperg(void)
{
  gsl_sf_result r;
  int s = 0;

  /* 0F1 */

  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1, 0.5, &r),     1.5660829297563505373, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (5, 0.5, &r),     1.1042674404828684574, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (100, 30, &r),    1.3492598639485110176, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-0.5, 3, &r),   -39.29137997543434276,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-100.5, 50, &r), 0.6087930289227538496, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1, -5.0, &r),   -0.3268752818235339109, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-0.5, -5.0, &r),-4.581634759005381184,  TEST_TOL1, GSL_SUCCESS);

  /* 1F1: see hyperg/test_1F1.c */

  /* U for integer parameters */

  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 0.0001, &r),  8.634088070212725330, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 0.01, &r),  4.078511443456425847, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 0.5, &r),  0.9229106324837304688, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 2.0, &r),  0.3613286168882225847, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 100, &r),  0.009901942286733018406, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 1000, &r),  0.0009990019940238807150, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 0.01, &r),  7.272361203006010000e+16, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 1, &r),  1957.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 5, &r),  1.042496, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 8, &r),  0.3207168579101562500, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 50, &r),  0.022660399001600000000, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 100, &r),  0.010631236727200000000, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 1000, &r),  0.0010060301203607207200, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 20, 1, &r),  1.7403456103284421000e+16, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 20, 20, &r),  0.22597813610531052969, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 50, 1, &r),  3.374452117521520758e+61, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 50, 50, &r),  0.15394136814987651785, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 0.1, &r),  1.0418325171990852858e+253, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 1, &r),  2.5624945006073464385e+154, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 50, &r),  3.0978624160896431391e+07, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 100, &r),  0.11323192555773717475, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 200, &r),  0.009715680951406713589, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 1000, &r),  0.0011085142546061528661, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1000, 2000, &r),  0.0009970168547036318206, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -1, 1, &r),  0.29817368116159703717, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -1, 10, &r),  0.07816669698940409380, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 1, &r),  0.08271753756946041959, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 5, &r),  0.06127757419425055261, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 10, &r),  0.04656199948873187212, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 20, &r),  0.031606421847946077709, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 0.01, &r),  0.009900000099999796950, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 1, &r),  0.009802970197050404429, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 10, &r),  0.009001648897173103447, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 20, &r),  0.008253126487166557546, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 50, &r),  0.006607993916432051008, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 90, &r),  0.005222713769726871937, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 110, &r),  0.004727658137692606210, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_int_e, (1, -1000, 1, &r),   0.0009980029970019970050, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -1000, 1010, &r),  0.0004971408839859245170, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 0.001, &r),  0.0007505359326875706975, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 0.5, &r),  6.449509938973479986e-06, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 8, &r),  6.190694573035761284e-10, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 20, &r),  3.647213845460374016e-12, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 8, 1, &r),  0.12289755012652317578, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 8, 10, &r),  5.687710359507564272e-09, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 8, 20, &r),  2.8175404594901039724e-11, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 0.01, &r),  1.0099979491941914867e+196, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 0.1, &r),  1.0090713562719862833e+97, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 1, &r),  0.009998990209084729106, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 20, &r),  1.3239363905866130603e-131, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 0.01, &r),  3.274012540759009536e+06, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 1, &r),  1.5202710000000000000e+06, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 10, &r),  1.0154880000000000000e+08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 100, &r),  3.284529863685482880e+19, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 10, 1, &r),  1.1043089864100000000e+11, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 1, &r),  1.3991152402448957897e+20, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 10, &r),  5.364469916567136000e+19, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 100, &r),  3.909797568000000000e+12, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 500, &r),  8.082625576697984130e+25, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 1, 0.01, &r),  1.6973422555823855798e+64, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 1, 1, &r),  7.086160198304780325e+63, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 1, 10, &r),  5.332862895528712200e+65, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 10, 1, &r),  -7.106713471565790573e+71, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 100, 1, &r),  2.4661377199407186476e+104, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 10, 10, &r),  5.687538583671241287e+68, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 100, 10, &r),  1.7880761664553373445e+102, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1, 0.01, &r),  4.185245354032917715e+137, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1, 0.1, &r),  2.4234043408007841358e+137, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1, 10, &r),  -1.8987677149221888807e+139, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 10, 10, &r),  -5.682999988842066677e+143, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 100, 10, &r),  2.3410029853990624280e+189, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1000, 10, &r),  1.9799451517572225316e+271, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -1, 10, &r),  -9.083195466262584149e+64, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -10, 10, &r),  -1.4418257327071634407e+62, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -100, 0.01, &r),  3.0838993811468983931e+93, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -100, 10, &r),  4.014552630378340665e+95, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-100, -100, 10, &r),  2.0556466922347982030e+162, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-100, -200, 10, &r),  1.1778399522973555582e+219, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-100, -200, 100, &r),  9.861313408898201873e+235, TEST_TOL3, GSL_SUCCESS);

  /* U */

  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 0.0001, 0.0001, &r), 1.0000576350699863577, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 0.0001, 1.0, &r), 0.9999403679233247536, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 0.0001, 100.0, &r), 0.9995385992657260887, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1, 0.0001, &r), 1.0009210608660065989, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1.0, 1.0, &r), 0.9999999925484179084, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 10, 1, &r), 13.567851006281412726, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 10, 5, &r), 1.0006265020064596364, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 10, 10, &r), 0.9999244381454633265, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 100, 1, &r),  2.5890615708804247881e+150, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 100, 10, &r),  2.3127845417739661466e+55, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 100, 50, &r), 6402.818715083582554, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 100, 98, &r), 0.9998517867411840044, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 1000, 300, &r),  2.5389557274938010716e+213, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1000, 999, &r), 0.9997195294193261604, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1000, 1100, &r),  0.9995342990014584713, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.5, 1000, 300, &r), 1.1977955438214207486e+217, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.5, 1000, 800, &r), 9.103916020464797207e+08, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.5, 1000, 998, &r), 0.21970269691801966806, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.5, 0.5, 1.0, &r), 0.7578721561413121060, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 0.0001, 0.0001, &r), 0.9992361337764090785, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 0.0001, 1, &r), 0.4036664068111504538, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 0.0001, 100, &r), 0.009805780851264329587, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 1.2, 2.0, &r), 0.3835044780075602550, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, -0.0001, 1, &r), 0.4036388693605999482, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (8, 10.5, 1, &r),  27.981926466707438538, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (8, 10.5, 10, &r),  2.4370135607662056809e-8, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (8, 10.5, 100, &r),  1.1226567526311488330e-16, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (10, -2.5, 10, &r),  6.734690720346560349e-14, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (10, 2.5, 10, &r),  6.787780794037971638e-13, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (10, 2.5, 50, &r),  2.4098720076596087125e-18, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 1, &r),  -3.990841457734147e+6, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 10, &r),  1.307472052129343e+8, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 50, &r),  3.661978424114088e+16, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 90, &r),  8.09469542130868e+19, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 99, &r),  2.546328328942063e+20, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 100, &r),  2.870463201832814e+20, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 200, &r),  8.05143453769373e+23, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 0.1, &r),  -3.043016255306515e+20, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 1, &r),  -3.194745265896115e+12, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 4, &r),  -6.764203430361954e+07, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 10, &r),  -2.067399425480545e+09, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 50, &r),  4.661837330822824e+14, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 10, &r),  -6.805460513724838e+66, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 50, &r),  -2.081052558162805e+18, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 80, &r),  2.034113191014443e+14, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 100, &r),  6.85047268436107e+13, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 200, &r),  1.430815706105649e+20, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-19.5, 82.1, 10, &r),  5.464313196201917432e+60, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 10, &r),  -5.5740216266953e+126, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 40, &r),  5.937463786613894e+91, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 50, &r),  -1.631898534447233e+89, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 70, &r),  3.249026971618851e+84, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 100, &r),  1.003401902126641e+85, TEST_TOL1, GSL_SUCCESS);

  /* Bug report from Stefan Gerlach */

  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.0, 4.0, 1.0, &r),  11.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.0, 0.5, 3.14, &r),  1.1896, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.0, 0.5, 1.13, &r),  -1.3631, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.0, 0.5, 0.0, &r),  0.75, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.0, 0.5, 1e-20, &r),  0.75, TEST_TOL2, GSL_SUCCESS);

  /* U(a,b,x) for x<0 [bug #27859] */

  /* Tests for b >= 0 */

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, 0, -0.1, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 0, -0.1, &r), -0.1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 0, -0.1, &r), 0.21, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 0, -0.1, &r), -0.661, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, 0, -0.1, &r), 2.7721, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, 0, -0.1, &r), -14.52201, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, 0, -0.1, &r), 91.230301, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, 1, -0.1, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 1, -0.1, &r), -1.1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 1, -0.1, &r), 2.41, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 1, -0.1, &r), -7.891, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, 1, -0.1, &r), 34.3361, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, 1, -0.1, &r), -186.20251, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, 1, -0.1, &r), 1208.445361, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 1, 2, -0.1, &r), -10.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, 2, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 2, -0.1, &r), -2.1, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 2, -0.1, &r), 6.61, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 2, -0.1, &r), -27.721, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, 2, -0.1, &r), 145.2201, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, 2, -0.1, &r), -912.30301, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, 2, -0.1, &r), 6682.263421, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 2, 3, -0.1, &r), 100.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 1, 3, -0.1, &r), 90.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, 3, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 3, -0.1, &r), -3.10, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 3, -0.1, &r), 12.81, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 3, -0.1, &r), -66.151, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, 3, -0.1, &r), 409.8241, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, 3, -0.1, &r), -2961.42351, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, 3, -0.1, &r), 24450.804481, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 3, 4, -0.1, &r), -1000.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 2, 4, -0.1, &r), -1900.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 1, 4, -0.1, &r), -1810.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, 4, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 4, -0.1, &r), -4.10, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 4, -0.1, &r), 21.01, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 4, -0.1, &r), -129.181, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, 4, -0.1, &r), 926.5481, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, 4, -0.1, &r), -7594.16401, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, 4, -0.1, &r), 70015.788541, TEST_TOL2, GSL_SUCCESS);

  /* Tests for b < 0 */

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, -1, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, -1, -0.1, &r), 0.9, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, -1, -0.1, &r), 0.01, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, -1, -0.1, &r), -0.031, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, -1, -0.1, &r), 0.1281, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, -1, -0.1, &r), -0.66151000000000007694, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, -1, -0.1, &r), 4.098241, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, -2, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, -2, -0.1, &r), 1.9, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, -2, -0.1, &r), 1.81, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, -2, -0.1, &r), -0.001, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, -2, -0.1, &r), 0.0041, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, -2, -0.1, &r), -0.02101, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, -2, -0.1, &r), 0.129181, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, -3, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, -3, -0.1, &r), 2.9, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, -3, -0.1, &r), 5.61, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, -3, -0.1, &r), 5.429, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, -3, -0.1, &r), 0.0001, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, -3, -0.1, &r), -0.00051, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, -3, -0.1, &r), 0.003121, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, -4, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, -4, -0.1, &r), 3.9, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, -4, -0.1, &r), 11.41, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, -4, -0.1, &r), 22.259, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, -4, -0.1, &r), 21.7161, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, -4, -0.1, &r), -1e-5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, -4, -0.1, &r), 0.000061, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-7, -4, -0.1, &r), -0.0004341, TEST_TOL0, GSL_SUCCESS);

  /* Tests for integer a */

  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 0.5, -0.5, &r), -9.5, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-8, 0.5, -0.5, &r), 180495.0625, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-8, 1.5, -0.5, &r), 827341.0625, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-8, 1.5, -10, &r), 7.162987810253906e9, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (3, 6, -0.5, &r), -296.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (3, 7, -0.5, &r), 2824, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (5, 12, -1.7, &r), -153.262676210016018065768591104, TEST_TOL2, GSL_SUCCESS);

  /* A few random tests */

  TEST_SF(s, gsl_sf_hyperg_U_e, (0, 0, -0.5, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0, 1, -0.5, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0, 1, -0.001, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 0.99, -0.1, &r), -1.09, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 0, -0.5, &r), -0.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 0, -0.5, &r), 1.25, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-7, 0, -0.1, &r), -668.2263421, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_int_e, (3, 6, -0.5, &r), -296.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (3, 7, -0.5, &r), 2824, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (5, 12, -1.7, &r), -153.262676210016018065768591104, TEST_TOL2, GSL_SUCCESS);

  /* Bug report from Raymond Rogers */

  TEST_SF(s, gsl_sf_hyperg_U_e, (4.11, 0.11, 6.4, &r), 6.422378238765078623739153038e-5, TEST_TOL2, GSL_SUCCESS);
  /* Addition tests from Raymond Rogers */
TEST_SF(s, gsl_sf_hyperg_U_e, (5, 4, 6.4, &r), 3.2586223825343211136628535e-05, TEST_TOL2, GSL_SUCCESS);
TEST_SF(s, gsl_sf_hyperg_U_e, (2.2,1.2 , 8.7, &r), 5.7250017539318661177749625e-03, TEST_TOL2, GSL_SUCCESS);
TEST_SF(s, gsl_sf_hyperg_U_e, (2, -6.4, 1, &r),1.2141502795806162484648638e-02
, TEST_TOL2, GSL_SUCCESS);

  /* added special handling with x == 1.0 , Richard J. Mathar, 2008-01-09 */

  /* Taylor Binnington a = 0 */

  /* Andrew Benson <abenson@caltech.edu> bug #24812
     in Pari:
     poch(a,x) = { gamma(a+x)/gamma(a) }
     t(a,b,c,x,k) = { (poch(a,k)*poch(b,k)/poch(c,k)) * (x^k)/(k!) }
     suminf(k=0,t(-10.34, 2.05, 3.05, 0.1725,k))  */

  /* Didier Pinchon also bug #24812 */

  /* Bill Maier - bug #45926 */

#if 0 /* XXX - bug #39056 */
  /* Test case from Hatef Monajemi <monajemi@stanford.edu> */

  /* Test case from Robert L Wolpert <Wolpert@stat.duke.edu> */

  /* Test case from ldnlwm@163.com */

#endif

  /* updated correct values, testing enabled, Richard J. Mathar, 2008-01-09 */

  TEST_SF(s, gsl_sf_hyperg_2F0_e, (0.01, 1.0, -0.02, &r), .99980388665511730901180717   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (0.1,  0.5, -0.02, &r), .99901595171179281891589794   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (1,   1, -0.02, &r), .98075549650574351826538049000    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (8,   8, -0.02, &r), .32990592849626965538692141   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (50, 50, -0.02, &r), .2688995263772964415245902e-12 , TEST_TOL0, GSL_SUCCESS);

  

  return s;
}
