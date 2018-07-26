#! /bin/bash

# Linux lib path
export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH}
# OSX, Darwin lib path
export DYLD_LIBRARY_PATH=${PWD}/lib:${DYLD_LIBRARY_PATH}
