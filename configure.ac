#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT(gemstat, 0.2)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/seq2expr.cpp])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_LANG(C++)

# Checks for libraries.
AX_PATH_GSL([1.14], [], [AC_MSG_ERROR(Could not find the required version of GSL)])

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
