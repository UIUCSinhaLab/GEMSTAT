


AC_DEFUN([AX_PATH_NLOPT],
[
AC_ARG_WITH(nlopt-prefix,[  --with-nlopt-prefix=PFX   Prefix where NLOpt is installed (optional)],
            nlopt_prefix="$withval", nlopt_prefix="")
AC_ARG_WITH(nlopt-exec-prefix,[  --with-nlopt-exec-prefix=PFX Exec prefix where NLOpt is installed (optional)],
            nlopt_exec_prefix="$withval", nlopt_exec_prefix="")
AC_ARG_ENABLE(nlopttest, [  --disable-nlopttest       Do not try to compile and run a test NLOpt program],
		    , enable_nlopttest=yes)

have_nlopt="yes"

ac_save_PKG_CONFIG_PATH="$PKG_CONFIG_PATH"
if test "x${nlopt_prefix}" != x ; then
	export PKG_CONFIG_PATH="${nlopt_prefix}/lib/pkgconfig"
fi


PKG_CHECK_MODULES([NLOPT], [nlopt],
[
  have_nlopt="yes"
  NLOPT_LIBS=${NLOPT_LIBS/-lnlopt /-lnlopt_cxx }
],
[
  have_nlopt="no"
])

AC_MSG_NOTICE("NLOPT_LIBS : $NLOPT_LIBS")
AC_MSG_NOTICE("NLOPT_CFLAGS : $NLOPT_CFLAGS")

export PKG_CONFIG_PATH=${ac_save_PKG_CONFIG_PATH}

if test "x$have_nlopt" = xyes ; then
	ifelse([$2], , :, [$2])
else
	ifelse([$3], , :, [$3])
fi


AC_SUBST(NLOPT_CFLAGS)
AC_SUBST(NLOPT_LIBS)
])

AU_ALIAS([AM_PATH_NLOPT], [AX_PATH_NLOPT])
