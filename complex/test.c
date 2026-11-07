/* complex/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * Copyright (C) 2020, 2021 Patrick Alken
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

#if defined(HAVE_COMPLEX_H) && defined(TEST_C11)
#include <complex.h>
#endif

#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>

#define BASE_DOUBLE
#include "templates_on.h"
#include "test_source.c"
#include "templates_off.h"
#undef BASE_DOUBLE

/*
 * Library functions that take and return complex numbers by value, for
 * all three types.  test_c11 compiles this with <complex.h>, where
 * gsl_complex may be double _Complex while the library is built with
 * the struct; the calls go through volatile pointers so that inline
 * versions cannot hide a calling-convention mismatch (bug #61342).
 */

#define TEST_BYVALUE(VEC, Z, R, ALLOC, FREE, GET, SET, SETALL)                \
  do {                                                                        \
    VEC * v = ALLOC (3);                                                      \
    Z (* volatile get) (const VEC *, const size_t) = GET;                     \
    void (* volatile set) (VEC *, const size_t, Z) = SET;                     \
    void (* volatile setall) (VEC *, Z) = SETALL;                             \
    Z z, w;                                                                   \
    GSL_SET_COMPLEX (&z, (R) 1.5, (R) -2.25);                                 \
    setall (v, z);                                                            \
    gsl_test (v->data[2 * 2] != (R) 1.5 || v->data[2 * 2 + 1] != (R) -2.25,   \
              #SETALL " passes the argument by value");                       \
    GSL_SET_COMPLEX (&z, (R) -0.5, (R) 4.0);                                  \
    set (v, 1, z);                                                            \
    gsl_test (v->data[2] != (R) -0.5 || v->data[3] != (R) 4.0,                \
              #SET " passes the argument by value");                          \
    w = get (v, 1);                                                           \
    gsl_test (GSL_REAL (w) != (R) -0.5 || GSL_IMAG (w) != (R) 4.0,            \
              #GET " returns the value");                                     \
    FREE (v);                                                                 \
  } while (0)

static void
test_byvalue (void)
{
  TEST_BYVALUE (gsl_vector_complex, gsl_complex, double,
                gsl_vector_complex_alloc, gsl_vector_complex_free,
                gsl_vector_complex_get, gsl_vector_complex_set,
                gsl_vector_complex_set_all);
  TEST_BYVALUE (gsl_vector_complex_float, gsl_complex_float, float,
                gsl_vector_complex_float_alloc, gsl_vector_complex_float_free,
                gsl_vector_complex_float_get, gsl_vector_complex_float_set,
                gsl_vector_complex_float_set_all);
  TEST_BYVALUE (gsl_vector_complex_long_double, gsl_complex_long_double, long double,
                gsl_vector_complex_long_double_alloc, gsl_vector_complex_long_double_free,
                gsl_vector_complex_long_double_get, gsl_vector_complex_long_double_set,
                gsl_vector_complex_long_double_set_all);

  /* GSL_SET_COMPLEX sets each part exactly (bug #65912) */
  {
    gsl_complex z;
    GSL_SET_COMPLEX (&z, 1.0, -0.0);
    gsl_test (GSL_REAL (z) != 1.0 || GSL_IMAG (z) != 0.0 || !(1.0 / GSL_IMAG (z) < 0.0),
              "GSL_SET_COMPLEX keeps the sign of a zero imaginary part");
    GSL_SET_COMPLEX (&z, 1.0, GSL_POSINF);
    gsl_test (GSL_REAL (z) != 1.0 || GSL_IMAG (z) != GSL_POSINF,
              "GSL_SET_COMPLEX with an infinite imaginary part");
    GSL_SET_COMPLEX (&z, GSL_NAN, 2.0);
    gsl_test (GSL_REAL (z) == GSL_REAL (z) || GSL_IMAG (z) != 2.0,
              "GSL_SET_COMPLEX with a NaN real part");
  }
}

int
main (void)
{
  test_all();
  test_byvalue();

  exit (gsl_test_summary ());
}
