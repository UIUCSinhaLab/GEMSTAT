#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT(gemstat, 0.2)
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/seq2expr.cpp])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_LANG(C++)

# Checks for libraries.
# Order matters here.
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([gslcblas], [cblas_dgemm])
AC_CHECK_LIB([gsl], [gsl_blas_dgemm])

# Checks for header files.
AC_CHECK_HEADERS([gsl/gsl_linalg.h gsl/gsl_math.h gsl/gsl_multimin.h gsl/gsl_randist.h], [GSLINC=-I/usr/include])

AC_SUBST([GSLINC])



# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt strchr])

AC_CONFIG_FILES([Makefile
		 src/Makefile])
AC_OUTPUT
