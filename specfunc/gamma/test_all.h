/* specfunc/gamma/test_all.h
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

#ifndef TEST_GAMMA_ALL_H
#define TEST_GAMMA_ALL_H

int test_gamma_all(void);

int test_gamma(void);
int test_gammainv(void);

int test_lngamma(void);
int test_lngamma1p(void);

int test_gamma_inc_P(void);
int test_gamma_inc_Q(void);

#endif /* TEST_GAMMA_ALL_H */
