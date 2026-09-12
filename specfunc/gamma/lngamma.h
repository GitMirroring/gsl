/* specfunc/gamma/lngamma.h
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

#ifndef __LNGAMMA_H__
#define __LNGAMMA_H__

/* lngamma1p region boundaries */
#define LNGAMMA_1P_X_RECUR0_RECUR1         (-0.6)
#define LNGAMMA_1P_X_RECUR1_NEAR0          (-0.4)
#define LNGAMMA_1P_X_NEAR0_NEAR1           ( 0.4)
#define LNGAMMA_1P_X_NEAR1_LNGAMMA         ( 2.0)

/* lngamma region boundaries */
#define LNGAMMA_X_NEG0_MIN                 (-28.0)
#define LNGAMMA_X_NEG0_MAX                 (-3.0)
#define LNGAMMA_X_NEG1_MAX                 (-2.0)
#define LNGAMMA_X_NEG2_MAX                 (-1.0)
#define LNGAMMA_X_SMALL_RECUR0_RECUR1      (LNGAMMA_1P_X_NEAR0_NEAR1)
#define LNGAMMA_X_SMALL_RECUR1_NEAR1       (0.8)
#define LNGAMMA_X_NEAR1_MIDDLE             (1.23395)
#define LNGAMMA_X_MIDDLE_NEAR2             (1.67895)
#define LNGAMMA_X_NEAR2_RECUR              (2.3)
#define LNGAMMA_X_RECUR_STIRLING           (8.0)

/* mathematical constants */
#define LNGAMMA_M_PI2_12_MINUS_HALF        0.3224670334241132182362075833230125946094749506033992188677791146850037 /* pi^2 / 12 - 1/2 */
#define LNGAMMA_M_1_MINUS_EULER            0.4227843350984671393934879099175975689578406640600764011942327651151323 /* 1 - M_EULER */
#define LNGAMMA_M_LNGAMMA_3_2             -0.1207822376352452223455184457816472122518527279025994683638684737573247 /* log(Gamma(3/2)) */
#define LNGAMMA_M_ABS_LNGAMMA_3_2          0.1207822376352452223455184457816472122518527279025994683638684737573247 /* |log(Gamma(3/2))| */

#define LNGAMMA_STIRLING_SCALE             (2.0 * LNGAMMA_X_RECUR_STIRLING * LNGAMMA_X_RECUR_STIRLING)

#endif /* __LNGAMMA_H__ */
