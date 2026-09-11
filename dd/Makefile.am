noinst_LTLIBRARIES = libgsldd.la 

check_PROGRAMS = test

libgsldd_la_SOURCES = dd.c

AM_CPPFLAGS = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgsldd.la ../ieee-utils/libgslieeeutils.la ../test/libgsltest.la ../sys/libgslsys.la ../err/libgslerr.la

pkginclude_HEADERS = gsl_dd.h

noinst_HEADERS =         \
  internal.h             \
  internal_common.c      \
  internal_exp.c         \
  internal_exp_source.c  \
  internal_fma.c         \
  internal_log.c         \
  internal_log_source.c  \
  internal_p.c           \
  internal_sqrt.c        \
  internal_sqrt_source.c
