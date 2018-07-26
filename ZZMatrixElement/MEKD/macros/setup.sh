export MEKD_COMPILE_WITH_ROOT=No
export MEKD_STANDALONE=Yes

if [ -f ${ROOTSYS}/bin/root-config ]; then
	export MEKD_COMPILE_WITH_ROOT=Yes
fi

echo "User arguments that are being passed to the compiler: " $@
echo "Compiling..."

#make clean
make $@
. setLocalLibrary.sh