noinst_LTLIBRARIES = libgslderiv.la

AM_CPPFLAGS = -I$(top_srcdir)

libgslderiv_la_SOURCES = \
  cand.c                 \
  deriv.c                \
  deriv_old.c            \
  directionaldiff.c      \
  fornberg.c             \
  gradient.c             \
  jacobian.c             \
  stencil.c

TESTS = $(check_PROGRAMS)
check_PROGRAMS = test #demo

test_SOURCES = test.c
test_LDADD =  libgslderiv.la ../linalg/libgsllinalg.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../sort/libgslsort.la ../permutation/libgslpermutation.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../rng/libgslrng.la ../poly/libgslpoly.la ../statistics/libgslstatistics.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la

#demo_SOURCES = demo.c
#demo_LDADD =  libgslderiv.la  ../vector/libgslvector.la ../block/libgslblock.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la

pkginclude_HEADERS =    \
  gsl_deriv.h           \
  gsl_deriv_cand.h      \
  gsl_deriv_constants.h \
  gsl_deriv_jacobian.h  \
  gsl_deriv_params.h    \
  gsl_deriv_stencil.h   \
  gsl_deriv_step.h

noinst_HEADERS =     \
  deriv_common.c     \
  test_adaptive.c    \
  test_d.c           \
  test_directional.c \
  test_f.c           \
  test_fixed.c       \
  test_fornberg.c    \
  test_g.c           \
  test_gradient.c    \
  test_j.c           \
  test_jacobian.c    \
  test_old.c
