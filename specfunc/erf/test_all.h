/* specfunc/erf/test_all.h
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

#ifndef TEST_ERF_ALL_H
#define TEST_ERF_ALL_H

int test_erf_all(void);
int test_erf(void);
int test_erf_Z(void);
int test_erf_Q(void);
int test_erfc(void);
int test_erfcx(void);
int test_log_erfc(void);

#endif /* TEST_ERF_ALL_H */
