noinst_LTLIBRARIES = libgslntuple.la

pkginclude_HEADERS = gsl_ntuple.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslntuple_la_SOURCES = ntuple.c

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslntuple.la ../histogram/libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

CLEANFILES = test.dat
