noinst_LTLIBRARIES = libgslfit.la 

pkginclude_HEADERS = gsl_fit.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslfit_la_SOURCES = linear.c

check_PROGRAMS = test

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslfit.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la

