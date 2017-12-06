#!/bin/bash

## \
#Snippet from \
# \
#http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in \
# \
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" \
# \
#end snippet \
##

BASE_DIR=${SCRIPT_DIR}/../../

DATA_DIR=${BASE_DIR}/data/

${BASE_DIR}/src/seq2expr -s ${DATA_DIR}/seqs.fa -e ${DATA_DIR}/expr.tab -m ${DATA_DIR}/factors.wtmx -f ${DATA_DIR}/factor_expr.tab -c ${DATA_DIR}/coop.txt -i ${DATA_DIR}/factor_info.txt -o Direct -oo SSE -na 1 -fo example_01.out -po example_01.par
