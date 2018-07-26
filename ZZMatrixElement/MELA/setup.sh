#!/bin/sh

#echo "$@"
if [[ "$1" == *"clean"* ]];then
	scramv1 b "$@"
	pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/fortran/
	make clean
	rm -f ../data/${SCRAM_ARCH}/libjhugenmela.so
	popd
else
	pushd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/fortran/
	make all
	if mv libjhugenmela.so ../data/${SCRAM_ARCH}/; then
		echo
		echo "...and you are running setup.sh, so this was just done."
		echo
		popd
		scramv1 b "$@"
	else
		echo
		echo "ERROR: something went wrong in mv, see ^ error message"
		echo
		popd
		return 1 >& /dev/null || exit 1 #return only works when sourced, exit will exit your whole session if sourced
	fi
fi
