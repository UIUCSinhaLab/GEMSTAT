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

cd ${SCRIPT_DIR}

mkdir -p usr

curl -O "http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz"
tar -zxvof nlopt-2.4.2.tar.gz
(cd nlopt-2.4.2 ; ./configure --without-guile --without-python --without-octave --without-matlab --with-cxx --prefix=${SCRIPT_DIR}/usr; make ; make install)
