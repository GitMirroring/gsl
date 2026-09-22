/* specfunc/constants.h
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

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

/* internal constants needed in special functions routines */

#define GSL_SF_LN2_HI              6.93147180369123816490e-01     /* high part of log(2) */
#define GSL_SF_LN2_LO              1.90821492927058770002e-10     /* low part of log(2) */
#define GSL_SF_DBL_EPSILON_4       5.55111512312578270212e-17     /* GSL_DBL_EPSILON / 4 = 2^{-54} */
#define GSL_SF_LOG_SQRT_2PI        9.18938533204672741780e-01     /* log(sqrt(2*pi)) */
#define GSL_SF_SQRT_2PI            2.50662827463100050242e+00     /* sqrt(2*pi) */
#define GSL_SF_INV_SQRT_2PI        3.98942280401432677940e-01     /* 1/sqrt(2*pi) */
#define GSL_SF_SQRT2_PI            7.97884560802865355880e-01     /* sqrt(2/pi) */
#define GSL_SF_1_MINUS_EULER       4.22784335098467139393e-01     /* 1 - M_EULER */
#define GSL_SF_PI2_12              8.22467033424113218236e-01     /* pi^2 / 12 */

#define GSL_SF_TRUE_MIN            4.94065645841246544177e-324    /* GSL_DBL_MIN * GSL_DBL_EPSILON = 2^{-1074} */
#define GSL_SF_LOG_TRUE_MIN       -7.44440071921381262314e+02     /* log(GSL_SF_TRUE_MIN) */
#define GSL_SF_1_DBL_EPSILON       4.50359962737049600000e+15     /* 1 / GSL_DBL_EPSILON = 2^{52} */

#endif /* !__CONSTANTS_H__ */
