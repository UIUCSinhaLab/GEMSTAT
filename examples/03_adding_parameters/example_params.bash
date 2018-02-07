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

DATA_ARGUMENTS="-s ${DATA_DIR}/seqs.fa -e ${DATA_DIR}/expr.tab -m ${DATA_DIR}/factors.wtmx -f ${DATA_DIR}/factor_expr.tab -c ${DATA_DIR}/coop.txt -i ${DATA_DIR}/factor_info.txt"

MODEL_ARGUMENTS="-o Direct -ct 50 -rt 0"

TRAINING_ARGUMENTS="-oo SSE -na 1"

${BASE_DIR}/src/seq2expr ${DATA_ARGUMENTS} ${MODEL_ARGUMENTS} ${TRAINING_ARGUMENTS} -fo example_03.out -po example_03.par \
	-onebeta \
	-ff ${SCRIPT_DIR}/ff.txt \
	-lower_bound ${SCRIPT_DIR}/lower.par \
	-upper_bound ${SCRIPT_DIR}/upper.par \
	-p ${SCRIPT_DIR}/start.par
