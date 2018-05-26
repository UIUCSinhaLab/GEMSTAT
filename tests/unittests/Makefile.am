ACLOCAL_AMFLAGS= -I m4
AUTOMAKE_OPTIONS = foreign

#optional_nlopt
if DISABLE_NLOPT
	NLOPT_CPPFLAGS = -DDISABLE_NLOPT
else
	NLOPT_CPPFLAGS =
endif

#check_LTLIBRARIES = $(top_builddir)/tests/gtest/googletest/lib/libgtest.la $(top_builddir)/tests/gtest/googletest/lib/libgtest_main.la

AM_CPPFLAGS = -I$(top_srcdir)/src $(GSL_CFLAGS) $(BETAOPTFLAGS) $(NLOPT_CPPFLAGS) $(NLOPT_CFLAGS)

$(top_srcdir)/tests/gtest/googletest/lib/libgtest.la:
	$(MAKE) -C $(top_builddir)/tests/gtest/googletest/ lib/libgtest.la

UNIT_TEST=
check_PROGRAMS=

check_PROGRAMS += bin/dummy
bin_dummy_SOURCES = src/dummy.cpp
bin_dummy_CPPFLAGS = -I$(top_srcdir)/src $(GTEST_CPPFLAGS) $(GSL_CFLAGS) $(BETAOPTFLAGS) $(NLOPT_CPPFLAGS) $(NLOPT_CFLAGS)
bin_dummy_CXXFLAGS = $(GTEST_CXXFLAGS)
bin_dummy_LDFLAGS = $(GTEST_LDFLAGS) 
bin_dummy_LDADD = $(GSL_LIBS) $(NLOPT_LIBS) \
		  $(top_builddir)/tests/gtest/googletest/lib/libgtest.la $(top_builddir)/tests/gtest/googletest/lib/libgtest_main.la
