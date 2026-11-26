/* sys/ldfrexp_table.h
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

#ifndef __LDFREXP_TABLE_H__
#define __LDFREXP_TABLE_H__

/* 2^i and 2^(32 i), i = 0, ..., 31, for ldexp_fallback and
   frexp_fallback: 2^n = ldfrexp_lo[n % 32] * ldfrexp_hi[n / 32] exactly
   for 0 <= n <= 1023.  17 significant digits determine these exactly
   (checked in sys/test.c). */

static const double ldfrexp_lo[32] = {
  1.0, 2.0, 4.0, 8.0,
  16.0, 32.0, 64.0, 128.0,
  256.0, 512.0, 1024.0, 2048.0,
  4096.0, 8192.0, 16384.0, 32768.0,
  65536.0, 131072.0, 262144.0, 524288.0,
  1048576.0, 2097152.0, 4194304.0, 8388608.0,
  16777216.0, 33554432.0, 67108864.0, 134217728.0,
  268435456.0, 536870912.0, 1073741824.0, 2147483648.0
};

static const double ldfrexp_hi[32] = {
  1.0, 4294967296.0, 1.8446744073709552e+19,
  7.9228162514264338e+28, 3.4028236692093846e+38, 1.4615016373309029e+48,
  6.2771017353866808e+57, 2.695994666715064e+67, 1.157920892373162e+77,
  4.9732323640978664e+86, 2.1359870359209101e+96, 9.173994463960286e+105,
  3.9402006196394479e+115, 1.6923032801030364e+125, 7.2683872429560689e+134,
  3.1217485503159922e+144, 1.3407807929942597e+154, 5.7586096570152914e+163,
  2.4733040147310453e+173, 1.0622759856335342e+183, 4.5624406176221952e+192,
  1.959553324262937e+202, 8.4162174424773976e+211, 3.6147378671465184e+221,
  1.5525180923007089e+231, 6.6680144328798543e+240, 2.8638903918474961e+250,
  1.2300315572313621e+260, 5.2829453113566525e+269, 2.269007733883336e+279,
  9.7453140113999991e+288, 4.1855804968213567e+298
};

#define LDFREXP_2_512 (ldfrexp_hi[16])            /* 2^512 */
#define LDFREXP_2_54  (ldfrexp_hi[1] * ldfrexp_lo[22])   /* 2^54 */

#endif /* __LDFREXP_TABLE_H__ */
