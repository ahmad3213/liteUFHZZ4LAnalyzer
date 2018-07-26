c++ -std=c++14 -O2 -Wall `root-config --cflags --ldflags --glibs` -L/cvmfs/cms.cern.ch/sl6_amd64_gcc530/lcg/root/6.06.00-ikhhed4//lib -lRooFit -lRooFitCore -lRooStats  -c KinZfitter.cpp
c++ -std=c++14 -O2 -Wall `root-config --cflags --ldflags --glibs` -L/cvmfs/cms.cern.ch/sl6_amd64_gcc530/lcg/root/6.06.00-ikhhed4//lib -lRooFit -lRooFitCore -lRooStats  -c HelperFunction.cc
