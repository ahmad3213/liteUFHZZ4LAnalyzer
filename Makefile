CXX=c++
CC=gcc
CFLAGS= -std=c++14 -O2 -Wall
INSS=-I./include
INCMSSW=-I${CMSSW_BASE}/src

CFLAGS += `root-config --cflags`
LIBS += `root-config --ldflags --glibs`
LIBS += -L$(ROOFITSYS)/lib -lRooFit -lRooFitCore -lRooStats
#LIBS += -L${CMSSW_BASE}/lib/slc6_amd64_gcc530/ -lZZMatrixElementMEKD -lZZMatrixElementMELA -lZZMatrixElementMEMCalculators -lZZMatrixElementPythonWrapper

zz4lOBJ=ZZ4L_Ana.o
#kinZfitterOBJ=KinZfitter/KinZfitter.o
#helperFunctionOBJ=KinZfitter/HelperFunction.o

.PHONY: clean all main test

all: ZZ4L_Ana

ZZ4L_Ana: ZZ4L_Ana.o 
		$(CXX) -o ZZ4L_Ana.exe $(zz4lOBJ) $(helperFunctionOBJ) $(LIBS)

clean:
	@rm *.o *.exe

##############RULES##############
.cc.o:
	$(CXX) $(CFLAGS) $(INSS) ${INCMSSW} -I. -c $<
.cpp.o:
	$(CXX) $(CFLAGS) $(INSS) ${INCMSSW} -I. -c $<

