AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = $(GSL_CFLAGS)
bin_PROGRAMS = seq2expr seqannot
seq2expr_SOURCES = seq2expr.cpp ExprPredictor.cpp ExprPredictor.h SeqAnnotator.cpp SeqAnnotator.h Tools.cpp Tools.h Utils.h IO.cpp IO.h
seq2expr_LDADD = $(GSL_LIBS)
seqannot_SOURCES = seqannot.cpp SeqAnnotator.cpp SeqAnnotator.h Tools.cpp Tools.h Utils.h IO.cpp IO.h
