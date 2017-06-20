#!/bin/bash

./src/seq2expr -s data/seqs.fa -e data/expr.tab -m data/factors.wtmx -f data/factor_expr.tab -c data/coop.txt -i data/factor_info.txt -o Markov -oo SSE -na 1 -fo foo.out -po foo.par
