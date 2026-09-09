/*
 * specfunc/log/test_all.c
 *
 * Unit-test aggregator.
 */

#include <config.h>

#include "test_all.h"

int
test_log_all(void)
{
  int status = 0;

  status += test_log1p();
  status += test_log1pmx();

  return status;
}
