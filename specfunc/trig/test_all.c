/*
 * specfunc/trig/test_all.c
 *
 * Unit-test aggregator.
 */

#include <config.h>

#include "test_all.h"

int
test_trig_all(void)
{
  int status = 0;

  status += test_sinc();

  return status;
}
