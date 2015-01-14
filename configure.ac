#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(gemstat, 0.2)
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/seq2expr.cpp])
#AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
#AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgsl':
AC_CHECK_LIB([gsl], [main])
# FIXME: Replace `main' with a function in `-lgslcblas':
AC_CHECK_LIB([gslcblas], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt strchr])

AC_CONFIG_FILES([Makefile
		 src/Makefile])
AC_OUTPUT
