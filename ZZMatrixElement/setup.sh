#!/bin/sh


#echo "$@"
pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/
./setup.sh "$@" || return 1 >& /dev/null || exit 1 #return only works when sourced, exit will exit your whole session if sourced
popd
scramv1 b "$@"
