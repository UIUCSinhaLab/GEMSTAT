AUTOMAKE_OPTIONS = foreign

if BETAOPTTOGETHER
BETAOPTFLAGS = -DBETAOPTTOGETHER
endif

if BETAOPTSEPARATE
BETAOPTFLAGS = -DBETAOPTSEPARATE
endif

if BETAOPTBROKEN
BETAOPTFLAGS = -DBETAOPTBROKEN
endif

AM_CPPFLAGS = $(GSL_CFLAGS) $(BETAOPTFLAGS)
bin_PROGRAMS = seq2expr
seq2expr_SOURCES = seq2expr.cpp INDExprPredictor.cpp INDExprPredictor.h ExprPredictor.cpp ExprPredictor.h SeqAnnotator.cpp SeqAnnotator.h Tools.cpp Tools.h Utils.h IO.cpp IO.h
seq2expr_LDADD = $(GSL_LIBS)
