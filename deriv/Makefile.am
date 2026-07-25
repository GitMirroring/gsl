noinst_LTLIBRARIES = libgslderiv.la

AM_CPPFLAGS = -I$(top_srcdir)

libgslderiv_la_SOURCES = \
  cand.c                 \
  deriv.c                \
  deriv_old.c            \
  fornberg.c             \
  stencil.c

TESTS = $(check_PROGRAMS)
check_PROGRAMS = test #demo

test_SOURCES = test.c
test_LDADD =  libgslderiv.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../sort/libgslsort.la ../permutation/libgslpermutation.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../rng/libgslrng.la ../poly/libgslpoly.la ../statistics/libgslstatistics.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la

#demo_SOURCES = demo.c
#demo_LDADD =  libgslderiv.la  ../vector/libgslvector.la ../block/libgslblock.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la

pkginclude_HEADERS = gsl_deriv.h
noinst_HEADERS = test_f.c
