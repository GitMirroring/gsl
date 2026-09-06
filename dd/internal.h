/* dd/internal.h
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

#ifndef __INTERNAL_H__
#define __INTERNAL_H__

#define DD_SPLIT                 134217729.0
#define DD_LN2_HI                6.93147180369123816490e-01
#define DD_LN2_LO                1.90821492927058770002e-10

#define DD_CAT_(a,b)             a##b
#define DD_CAT(a,b)              DD_CAT_(a,b)

#endif /* __INTERNAL_H__ */
