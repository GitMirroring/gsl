/*
 * specfunc/gamma/test_all.c
 *
 * Gamma-family unit-test aggregator.
 */

#include <config.h>
#include "test_all.h"

int
test_gamma_all(void)
{
  int status = 0;

  status += test_gamma_inc_P();
  status += test_gamma_inc_Q();
  status += test_lngamma();
  status += test_lngamma1p();

  return status;
}
