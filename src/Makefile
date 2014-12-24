CC = g++

GSL_DIR = /home/grad/samee1/packages/gsl-1.14

INC = -I$(GSL_DIR)
CFLAGS = -O3 $(INC) 

LIB_DIRS = -L/usr/lib64
LIBS = -lgsl -lgslcblas -lm
LFLAGS = -O3 $(LIB_DIRS) $(LIBS)

BIN = seq2expr 

all: $(BIN)

clean:
	rm -f $(BIN)

Tools.o : Tools.h Tools.cpp
	$(CC) $(CFLAGS) -g -c Tools.cpp
SeqAnnotator.o : Tools.h SeqAnnotator.h SeqAnnotator.cpp
	$(CC) $(CFLAGS) -g -c SeqAnnotator.cpp
ExprPredictor.o : Tools.h SeqAnnotator.h ExprPredictor.h ExprPredictor.cpp
	$(CC) $(CFLAGS) -g -c ExprPredictor.cpp
seq2expr.o : Tools.h SeqAnnotator.h ExprPredictor.h seq2expr.cpp
	$(CC) $(CFLAGS) -g -c seq2expr.cpp

seq2expr : Tools.o SeqAnnotator.o ExprPredictor.o seq2expr.o 
	$(CC) $(LFLAGS) -g -o $@ Tools.o SeqAnnotator.o ExprPredictor.o seq2expr.o

