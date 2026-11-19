# -*-makefile-*-

noinst_LTLIBRARIES = libgslroots.la 

pkginclude_HEADERS = gsl_roots.h

noinst_HEADERS = roots.h test_funcs.c test_aps.c

AM_CPPFLAGS = -I$(top_srcdir)

libgslroots_la_SOURCES = bisection.c brent.c falsepos.c aps.c	newton.c  secant.c steffenson.c	convergence.c fsolver.c fdfsolver.c 

check_PROGRAMS = test

TESTS = $(check_PROGRAMS)

test_SOURCES = test.c
test_LDADD = libgslroots.la ../deriv/libgslderiv.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../sort/libgslsort.la ../permutation/libgslpermutation.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../statistics/libgslstatistics.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

