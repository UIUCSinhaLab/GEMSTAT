#!/bin/bash
../../src/seqannot -m factors.wtmx -s seqs.fa > seq.fa.annotations

../../src/seq2expr -m factors.wtmx -s seqs.fa -e expr.1 -i factor_info.txt  -f factor_expr.tab -c coop.1 -o Direct -oo SSE -fo out.1 -po par.1 -na 1 -no_gt_out

../../src/seq2expr -m factors.wtmx -s seqs.fa -e expr.1 -i factor_info.txt  -f factor_expr.tab -c coop.2 -o Direct -oo SSE -fo out.2 -po par.2 -na 0 -p par.coop.in -no_gt_out

../../src/seq2expr -m factors.wtmx -s seqs.fa -e expr.1 -i factor_info.txt  -f factor_expr.tab -c coop.3 -o Direct -oo SSE -fo out.3 -po par.3 -na 0 -p par.coop.in -no_gt_out

../../src/seq2expr -m factors.wtmx -s seqs.4.fa -e expr.4 -i factor_info.txt  -f factor_expr.tab -c coop.4 -o Direct -oo SSE -fo out.4 -po par.4 -na 0 -p par.coop.in.4 -no_gt_out
