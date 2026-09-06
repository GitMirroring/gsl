noinst_LTLIBRARIES = libgsldd.la 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_dd.h

libgsldd_la_SOURCES = dd.c

AM_CPPFLAGS = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgsldd.la ../ieee-utils/libgslieeeutils.la ../test/libgsltest.la ../sys/libgslsys.la ../err/libgslerr.la
