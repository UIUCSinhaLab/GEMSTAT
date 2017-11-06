all: snot.a

snot.tab.o snot.lex.o:  snot.tab.cpp snot.tab.hpp snot.lex.cpp snot.yy.hpp param_storage.h
	g++ -g -c snot.tab.cpp snot.lex.cpp

snot.a: snot.tab.o snot.lex.o
	ar -rcs ./libsnot.a snot.tab.o snot.lex.o

