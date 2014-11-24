CC = g++

GSL_DIR = /shared-mounts/hunchback/gsl-1.8

INC = -I$(GSL_DIR)
CFLAGS = -O3 $(INC) 

LIB_DIRS = -L$(GSL_DIR)/lib
LIBS = -lgsl -lgslcblas -lm
LFLAGS = -O3 $(LIB_DIRS) $(LIBS)

BIN = seq2expr 

all: $(BIN)

clean:
	rm -f $(BIN)

Tools.o : Tools.h Tools.cpp
	$(CC) $(CFLAGS) -c Tools.cpp
SeqAnnotator.o : Tools.h SeqAnnotator.h SeqAnnotator.cpp
	$(CC) $(CFLAGS) -c SeqAnnotator.cpp
ExprPredictor.o : Tools.h SeqAnnotator.h ExprPredictor.h ExprPredictor.cpp
	$(CC) $(CFLAGS) -c ExprPredictor.cpp
seq2expr.o : Tools.h SeqAnnotator.h ExprPredictor.h seq2expr.cpp
	$(CC) $(CFLAGS) -c seq2expr.cpp

seq2expr : Tools.o SeqAnnotator.o ExprPredictor.o seq2expr.o 
	$(CC) $(LFLAGS) -o $@ Tools.o SeqAnnotator.o ExprPredictor.o seq2expr.o

