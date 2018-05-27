#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT(gemstat, 1.6)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/seq2expr.cpp])
AM_INIT_AUTOMAKE([subdir-objects])

AC_ARG_ENABLE([debug],
[  --enable-debug	Turn on debugging],
[case "${enableval}" in
yes) debug=true ;;
no) debug=false ;;
*) AC_MSG_ERROR([bad value ${enablevar} for --enable-debug]) ;;
esac], [debug=false])
AM_CONDITIONAL([DEBUG], [test x${debug} = xtrue])

#beta optimization
AC_ARG_WITH([beta-optimize],
[  --with-beta-optimize	How to optimize beta, one of 'broken' (traditional), 'separate', 'together'],
[case "${withval}" in
broken) betaopt="broken" ;;
separate) betaopt="separate" ;;
together) betaopt="together" ;;
*) AC_MSG_ERROR([bad value ${withvar} for --with-beta-optimize])  ;;
esac], [betaopt="together"])
AM_CONDITIONAL([BETAOPTBROKEN], [test "x${betaopt}" = "xbroken"])
AM_CONDITIONAL([BETAOPTSEPARATE], [test "x${betaopt}" = "xseparate"])
AM_CONDITIONAL([BETAOPTTOGETHER], [test "x${betaopt}" = "xtogether"])

AC_ARG_ENABLE([gtest],
[  --enable-gtest decides to use google test or not.],
[to_use_gtest="yes"],
[to_use_gtest="no"])
AM_CONDITIONAL([USE_GTEST], [test "x${to_use_gtest}" = "xyes"])


# Checks for libraries.
AX_PATH_GSL([1.13], [], [AC_MSG_ERROR(Could not find the required version of GSL)])

AX_PATH_NLOPT([2.4], [], [AC_MSG_ERROR(Could not find the required version of NLOPT)])

#
#turn NLOpt on or off.
#
AC_ARG_ENABLE([nlopt],
[  --disable-nlopt  Disable using the NLOpt optimizer, which enables using the gsl optimizer.],
[with_nlopt="no"],
[with_nlopt="yes"])
AM_CONDITIONAL([DISABLE_NLOPT], [test "x${with_nlopt}" = "xno"])

LT_INIT


# Checks for programs.
AC_PROG_CXX
AC_LANG(C++)

# Order matters here.
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([gslcblas], [cblas_dgemm])
AC_CHECK_LIB([gsl], [gsl_blas_dgemm])

# Checks for header files.
AC_CHECK_HEADERS([gsl/gsl_linalg.h gsl/gsl_math.h gsl/gsl_multimin.h gsl/gsl_randist.h])


# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt strchr])

AC_CONFIG_FILES([Makefile
		 src/Makefile])


# Setup the test suite.
AC_CONFIG_TESTDIR(tests)
AC_CONFIG_FILES([tests/Makefile tests/atlocal])
AC_CONFIG_FILES([tests/testexec], [chmod +x tests/testexec])
AM_MISSING_PROG([AUTOM4TE], [autom4te])

AM_COND_IF([USE_GTEST], [

#Some of these variables are copied or inspired by: https://github.com/eiichiroi/autotools-unittest
# Define variables for unit test using gtest
GTEST_VERSION="1.8.0"
GTEST_CPPFLAGS="-DGTEST_HAS_PTHREAD=0 -I\$(top_srcdir)/tests/gtest/googletest/include/"
GTEST_CXXFLAGS=
GTEST_LDFLAGS="-L\$(top_builddir)/tests/gtest/googletest/lib/"
GTEST_LIBS=
AC_SUBST([GTEST_VERSION])
AC_SUBST([GTEST_CPPFLAGS])
AC_SUBST([GTEST_CXXFLAGS])
AC_SUBST([GTEST_LDFLAGS])
AC_SUBST([GTEST_LIBS])


	AC_CONFIG_SUBDIRS([tests/gtest/])
	AC_CONFIG_FILES([tests/gtest/Makefile])
	AC_CONFIG_SRCDIR([tests/unittests/src/dummy.cpp])
	AC_CONFIG_FILES([tests/unittests/Makefile])
])

AC_OUTPUT
