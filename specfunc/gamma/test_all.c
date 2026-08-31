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

  status += test_lngamma();
  status += test_lngamma_1p();

  return status;
}
