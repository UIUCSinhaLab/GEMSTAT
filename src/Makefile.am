AUTOMAKE_OPTIONS = foreign

#optional_nlopt
if DISABLE_NLOPT
NLOPT_LIBS =
NLOPT_CPPFLAGS = -DDISABLE_NLOPT
else
NLOPT_LIBS = -lnlopt_cxx
NLOPT_CPPFLAGS =
endif

#beta optimization stuff
if BETAOPTTOGETHER
BETAOPTFLAGS = -DBETAOPTTOGETHER
endif

if BETAOPTSEPARATE
BETAOPTFLAGS = -DBETAOPTSEPARATE
endif

if BETAOPTBROKEN
BETAOPTFLAGS = -DBETAOPTBROKEN
endif

AM_CPPFLAGS = $(GSL_CFLAGS) $(BETAOPTFLAGS) $(NLOPT_CPPFLAGS)
bin_PROGRAMS = seq2expr seqannot
seq2expr_SOURCES = seq2expr.cpp FactorIntFunc.cpp FactorIntFunc.h ExprModel.cpp ExprPar.cpp ExprPar.h ExprModel.h PredictorTrainer.h PredictorTrainer.cpp ExprPredictor.cpp ExprPredictor.h SeqAnnotator.cpp SeqAnnotator.h Tools.cpp Tools.h Utils.h IO.cpp IO.h ObjFunc.cpp ObjFunc.h
seq2expr_LDADD = $(GSL_LIBS) $(NLOPT_LIBS)
seqannot_SOURCES = seqannot.cpp SeqAnnotator.cpp SeqAnnotator.h Tools.cpp Tools.h Utils.h
seqannot_LDADD = $(GSL_LIBS) $(NLOPT_LIBS)
