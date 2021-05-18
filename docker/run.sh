#!/bin/sh

###Snippet from http://stackoverflow.com/questions/59895/
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
###end snippet

IMAGE_NAME="luntlab/gemstat_build:latest"

#See this note about IPC https://discuss.pytorch.org/t/unable-to-write-to-file-torch-18692-1954506624/9990

docker run -it --rm -v${SCRIPT_DIR}/..:/workspace/src --ipc=host ${IMAGE_NAME}
