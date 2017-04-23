


AC_DEFUN([AX_PATH_NLOPT],
[
AC_ARG_WITH(nlopt-prefix,[  --with-nlopt-prefix=PFX   Prefix where NLOpt is installed (optional)],
            nlopt_prefix="$withval", nlopt_prefix="")
AC_ARG_WITH(nlopt-exec-prefix,[  --with-nlopt-exec-prefix=PFX Exec prefix where NLOpt is installed (optional)],
            nlopt_exec_prefix="$withval", nlopt_exec_prefix="")
AC_ARG_ENABLE(nlopttest, [  --disable-nlopttest       Do not try to compile and run a test NLOpt program],
		    , enable_nlopttest=yes)

have_nlopt="yes"

PKG_CHECK_MODULES([NLOPT], [nlopt],
[
  have_nlopt="yes"
  NLOPT_LIBS="-lnlopt_cxx -lm"
],
[
  have_nlopt="no"
])

AC_SUBST(NLOPT_CFLAGS)
AC_SUBST(NLOPT_LIBS)
])

AU_ALIAS([AM_PATH_NLOPT], [AX_PATH_NLOPT])
