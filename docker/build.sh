#!/bin/bash

###Snippet from http://stackoverflow.com/questions/59895/
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
###end snippet

CONTAINER_VERSION="v0.1"

docker build -t luntlab/gemstat_build:${CONTAINER_VERSION} --target gemstat_build ${SCRIPT_DIR} || exit 1
docker tag luntlab/gemstat_build:${CONTAINER_VERSION} luntlab/gemstat_build:latest

#docker push luntlab/gemstat_build:${CONTAINER_VERSION}
#docker push luntlab/gemstat_build:latest
