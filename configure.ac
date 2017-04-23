#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT(gemstat, 1.1)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/seq2expr.cpp])
AM_INIT_AUTOMAKE

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
esac], [betaopt="broken"])
AM_CONDITIONAL([BETAOPTBROKEN], [test "x${betaopt}" = "xbroken"])
AM_CONDITIONAL([BETAOPTSEPARATE], [test "x${betaopt}" = "xseparate"])
AM_CONDITIONAL([BETAOPTTOGETHER], [test "x${betaopt}" = "xtogether"])


# Checks for libraries.
AX_PATH_GSL([1.13], [], [AC_MSG_ERROR(Could not find the required version of GSL)])

AX_PATH_NLOPT()

#
#turn NLOpt on or off.
#
AC_ARG_ENABLE([nlopt],
[  --disable-nlopt  Disable using the NLOpt optimizer, which enables using the gsl optimizer.],
[with_nlopt="no"],
[with_nlopt="yes"])
AM_CONDITIONAL([DISABLE_NLOPT], [test "x${with_nlopt}" = "xno"])


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

AC_OUTPUT
