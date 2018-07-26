LIBDIR=lib
MEKDDIR=..
INCDIR1=$(MEKDDIR)/src
INCDIR2=$(MEKDDIR)/src/Extra_code
INCDIR3=$(MEKDDIR)/src/higgs_properties
INCDIR4=$(MEKDDIR)/src/MadGraphSrc
INCDIR5=$(MEKDDIR)/interface
PDF_DIR=$(MEKDDIR)/src/PDFTables
CXXFLAGS= -pedantic -Wunused-variable -Werror=sign-compare -Wunused-but-set-variable -Werror=return-type -Werror=missing-braces -fPIC -O -I $(INCDIR1) -I $(INCDIR2) -I $(INCDIR3) -I $(INCDIR4) -I $(INCDIR5) -I.
CFLAGS= -pedantic -Wunused-variable -Werror=sign-compare -Wunused-but-set-variable -Werror=return-type -Werror=missing-braces -fPIC -std=c99
LIBFLAGS= -L$(LIBDIR)
UNAME := $(shell uname)

########## ROOT PART ##########
ROOTLIB=$(ROOTSYS)/lib
ROOTINCLUDES=$(ROOTSYS)/include
###############################




#############################
### Amplitudes for models ###
#############################

objects1=$(INCDIR4)/HelAmps_HEF_UFO_bkg.o $(INCDIR4)/HelAmps_HEF_MEKD.o $(INCDIR4)/HelAmps_HEF_MEKD2_1.o $(INCDIR4)/HelAmps_HZZ_Unitary_bkgpA.o $(INCDIR4)/Parameters_MEKD.o $(INCDIR4)/read_slha.o



############################
### SM Z-based processes ###
############################

objects2_Z4l=$(INCDIR4)/Z4l/qq_Z4l_BKG_DN_OF.o $(INCDIR4)/Z4l/qq_Z4l_BKG_DN_OFpA.o $(INCDIR4)/Z4l/qq_Z4l_BKG_DN_SF.o $(INCDIR4)/Z4l/qq_Z4l_BKG_DN_SFpA.o $(INCDIR4)/Z4l/qq_Z4l_BKG_UP_OF.o $(INCDIR4)/Z4l/qq_Z4l_BKG_UP_OFpA.o $(INCDIR4)/Z4l/qq_Z4l_BKG_UP_SF.o $(INCDIR4)/Z4l/qq_Z4l_BKG_UP_SFpA.o $(INCDIR4)/Z4l/qq_Z4l_SIG_DN_OF.o $(INCDIR4)/Z4l/qq_Z4l_SIG_DN_OFpA.o $(INCDIR4)/Z4l/qq_Z4l_SIG_DN_SF.o $(INCDIR4)/Z4l/qq_Z4l_SIG_DN_SFpA.o $(INCDIR4)/Z4l/qq_Z4l_SIG_UP_OF.o $(INCDIR4)/Z4l/qq_Z4l_SIG_UP_OFpA.o $(INCDIR4)/Z4l/qq_Z4l_SIG_UP_SF.o $(INCDIR4)/Z4l/qq_Z4l_SIG_UP_SFpA.o

objects2_ZZ=$(INCDIR4)/ZZ/qq_ZZ_DN_OF.o $(INCDIR4)/ZZ/qq_ZZ_UP_OF.o $(INCDIR4)/ZZ/qq_ZZ_DN_SF.o $(INCDIR4)/ZZ/qq_ZZ_UP_SF.o $(INCDIR4)/ZZ/qq_ZZ_DN_OFpA.o $(INCDIR4)/ZZ/qq_ZZ_UP_OFpA.o $(INCDIR4)/ZZ/qq_ZZ_DN_SFpA.o $(INCDIR4)/ZZ/qq_ZZ_UP_SFpA.o

objects2_DY=$(INCDIR4)/DY/qq_DY_DN_2l.o $(INCDIR4)/DY/qq_DY_UP_2l.o $(INCDIR4)/DY/qq_DY_DN_2lpA.o $(INCDIR4)/DY/qq_DY_UP_2lpA.o $(INCDIR4)/DY/DY_2l.o $(INCDIR4)/DY/DY_2lpA.o



##################################
### Spin-0 resonance processes ###
##################################

objects2_spin0=$(INCDIR4)/Spin0/gg_Spin0_OF.o $(INCDIR4)/Spin0/gg_Spin0_SF.o $(INCDIR4)/Spin0/gg_Spin0_OFpA.o $(INCDIR4)/Spin0/gg_Spin0_SFpA.o $(INCDIR4)/Spin0/Spin0_OF.o $(INCDIR4)/Spin0/Spin0_SF.o $(INCDIR4)/Spin0/Spin0_OFpA.o $(INCDIR4)/Spin0/Spin0_SFpA.o
# $(INCDIR4)/Spin0/qq_Spin0_DN_OF.o $(INCDIR4)/Spin0/qq_Spin0_UP_OF.o $(INCDIR4)/Spin0/qq_Spin0_DN_SF.o $(INCDIR4)/Spin0/qq_Spin0_UP_SF.o $(INCDIR4)/Spin0/qq_Spin0_DN_OFpA.o $(INCDIR4)/Spin0/qq_Spin0_UP_OFpA.o $(INCDIR4)/Spin0/qq_Spin0_DN_SFpA.o $(INCDIR4)/Spin0/qq_Spin0_UP_SFpA.o

objects2_spin0_2=$(INCDIR4)/Spin0/gg_Spin0_2l.o $(INCDIR4)/Spin0/gg_Spin0_2lpA.o $(INCDIR4)/Spin0/Spin0_2l.o $(INCDIR4)/Spin0/Spin0_2lpA.o
# $(INCDIR4)/Spin0/qq_Spin0_DN_2l.o $(INCDIR4)/Spin0/qq_Spin0_DN_2lpA.o $(INCDIR4)/Spin0/qq_Spin0_UP_2l.o $(INCDIR4)/Spin0/qq_Spin0_UP_2lpA.o

objects2_spin0_3=$(INCDIR4)/Spin0/gg_Spin0_2f_OF.o $(INCDIR4)/Spin0/gg_Spin0_2f_SF.o $(INCDIR4)/Spin0/gg_Spin0_2f_OFpA.o $(INCDIR4)/Spin0/gg_Spin0_2f_SFpA.o $(INCDIR4)/Spin0/Spin0_2f_OF.o $(INCDIR4)/Spin0/Spin0_2f_SF.o $(INCDIR4)/Spin0/Spin0_2f_OFpA.o $(INCDIR4)/Spin0/Spin0_2f_SFpA.o



##################################
### Spin-1 resonance processes ###
##################################

objects2_spin1=$(INCDIR4)/Spin1/qq_Spin1_DN_OF.o $(INCDIR4)/Spin1/qq_Spin1_UP_OF.o $(INCDIR4)/Spin1/qq_Spin1_DN_SF.o $(INCDIR4)/Spin1/qq_Spin1_UP_SF.o $(INCDIR4)/Spin1/qq_Spin1_DN_OFpA.o $(INCDIR4)/Spin1/qq_Spin1_UP_OFpA.o $(INCDIR4)/Spin1/qq_Spin1_DN_SFpA.o $(INCDIR4)/Spin1/qq_Spin1_UP_SFpA.o $(INCDIR4)/Spin1/Spin1_OF.o $(INCDIR4)/Spin1/Spin1_SF.o $(INCDIR4)/Spin1/Spin1_OFpA.o $(INCDIR4)/Spin1/Spin1_SFpA.o

objects2_spin1_2=$(INCDIR4)/Spin1/qq_Spin1_DN_2l.o $(INCDIR4)/Spin1/qq_Spin1_UP_2l.o $(INCDIR4)/Spin1/qq_Spin1_DN_2lpA.o $(INCDIR4)/Spin1/qq_Spin1_UP_2lpA.o  $(INCDIR4)/Spin1/Spin1_2l.o $(INCDIR4)/Spin1/Spin1_2lpA.o

objects2_spin1_3=$(INCDIR4)/Spin1/qq_Spin1_2f_DN_OF.o $(INCDIR4)/Spin1/qq_Spin1_2f_UP_OF.o $(INCDIR4)/Spin1/qq_Spin1_2f_DN_SF.o $(INCDIR4)/Spin1/qq_Spin1_2f_UP_SF.o $(INCDIR4)/Spin1/qq_Spin1_2f_DN_OFpA.o $(INCDIR4)/Spin1/qq_Spin1_2f_UP_OFpA.o $(INCDIR4)/Spin1/qq_Spin1_2f_DN_SFpA.o $(INCDIR4)/Spin1/qq_Spin1_2f_UP_SFpA.o $(INCDIR4)/Spin1/Spin1_2f_OF.o $(INCDIR4)/Spin1/Spin1_2f_SF.o $(INCDIR4)/Spin1/Spin1_2f_OFpA.o $(INCDIR4)/Spin1/Spin1_2f_SFpA.o



##################################
### Spin-2 resonance processes ###
##################################

objects2_spin2=$(INCDIR4)/Spin2/gg_Spin2_OF.o $(INCDIR4)/Spin2/gg_Spin2_SF.o $(INCDIR4)/Spin2/gg_Spin2_OFpA.o $(INCDIR4)/Spin2/gg_Spin2_SFpA.o $(INCDIR4)/Spin2/qq_Spin2_DN_OF.o $(INCDIR4)/Spin2/qq_Spin2_UP_OF.o $(INCDIR4)/Spin2/qq_Spin2_DN_SF.o $(INCDIR4)/Spin2/qq_Spin2_UP_SF.o $(INCDIR4)/Spin2/qq_Spin2_DN_OFpA.o $(INCDIR4)/Spin2/qq_Spin2_UP_OFpA.o $(INCDIR4)/Spin2/qq_Spin2_DN_SFpA.o $(INCDIR4)/Spin2/qq_Spin2_UP_SFpA.o $(INCDIR4)/Spin2/Spin2_OF.o $(INCDIR4)/Spin2/Spin2_SF.o $(INCDIR4)/Spin2/Spin2_OFpA.o $(INCDIR4)/Spin2/Spin2_SFpA.o

objects2_spin2_2=$(INCDIR4)/Spin2/gg_Spin2_2l.o $(INCDIR4)/Spin2/qq_Spin2_DN_2l.o $(INCDIR4)/Spin2/qq_Spin2_UP_2l.o $(INCDIR4)/Spin2/gg_Spin2_2lpA.o $(INCDIR4)/Spin2/qq_Spin2_DN_2lpA.o $(INCDIR4)/Spin2/qq_Spin2_UP_2lpA.o $(INCDIR4)/Spin2/Spin2_2l.o $(INCDIR4)/Spin2/Spin2_2lpA.o 

objects2_spin2_3=$(INCDIR4)/Spin2/gg_Spin2_2f_OF.o $(INCDIR4)/Spin2/gg_Spin2_2f_SF.o $(INCDIR4)/Spin2/gg_Spin2_2f_OFpA.o $(INCDIR4)/Spin2/gg_Spin2_2f_SFpA.o $(INCDIR4)/Spin2/qq_Spin2_2f_DN_OF.o $(INCDIR4)/Spin2/qq_Spin2_2f_UP_OF.o $(INCDIR4)/Spin2/qq_Spin2_2f_DN_SF.o $(INCDIR4)/Spin2/qq_Spin2_2f_UP_SF.o $(INCDIR4)/Spin2/qq_Spin2_2f_DN_OFpA.o $(INCDIR4)/Spin2/qq_Spin2_2f_UP_OFpA.o $(INCDIR4)/Spin2/qq_Spin2_2f_DN_SFpA.o $(INCDIR4)/Spin2/qq_Spin2_2f_UP_SFpA.o $(INCDIR4)/Spin2/Spin2_2f_OF.o $(INCDIR4)/Spin2/Spin2_2f_SF.o $(INCDIR4)/Spin2/Spin2_2f_OFpA.o $(INCDIR4)/Spin2/Spin2_2f_SFpA.o



#####################
### Misc. objects ###
#####################

objects3=$(INCDIR2)/MEKD_CalcHEP_PDF.o $(INCDIR2)/MEKD_CalcHEP_Extra_functions.o $(INCDIR2)/MEKD_MG_Boosts.o $(INCDIR3)/hggeffective.o $(PDF_DIR)/pdt.o

objects4=$(INCDIR1)/MEKD_MG.o

## RAW ME from MG5_aMC
# objects5=$(INCDIR4)/RAW/CPPProcess.o $(INCDIR4)/RAW/HelAmps_sm.o $(INCDIR4)/RAW/Parameters_sm.o $(INCDIR4)/RAW/read_slha.o



########################
### Combined objects ###
########################

objects2=$(objects2_ZZ) $(objects2_Z4l) $(objects2_DY) $(objects2_spin0) $(objects2_spin0_2) $(objects2_spin0_3) $(objects2_spin1) $(objects2_spin1_2) $(objects2_spin1_3) $(objects2_spin2) $(objects2_spin2_2) $(objects2_spin2_3)

objects=$(objects1) $(objects2) $(objects3) $(objects4) $(objects5)



###################
### All targets ###
###################

targetlib=$(LIBDIR)/libMEKD_MG.so
target=runKD_MAD
test=MEKD_Test


### The uncommented target that is going to be built
test: $(test)

# all: $(target)




############################
### Linux building block ###
############################

ifeq ($(UNAME), Linux)
# do something Linux-y
ifeq ($(MEKD_STANDALONE), Yes)
CXXFLAGS+= -D MEKD_STANDALONE
endif
$(targetlib): $(objects)
	$(CXX) $(CXXFLAGS) -fPIC -shared -Wl,-soname,libMEKD_MG.so -o $(targetlib) $(objects) $(LIBFLAGS)


ifeq ($(MEKD_COMPILE_WITH_ROOT), Yes)
CXXFLAGS+= -D MEKD_with_ROOT -std=c++0x

$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) `$(ROOTSYS)/bin/root-config --ldflags` -I $(ROOTINCLUDES) $(target).cpp -o $(target) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)

$(test): $(targetlib)
	$(CXX) $(CXXFLAGS) `$(ROOTSYS)/bin/root-config --ldflags` -I $(ROOTINCLUDES) $(test).cpp -o $(test) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)
else
$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) $(target).cpp -o $(target) -L. $(LIBFLAGS) $(targetlib)

$(test): $(targetlib)
	$(CXX) $(CXXFLAGS) $(test).cpp -o $(test) -L. $(LIBFLAGS) $(targetlib)
endif

endif




#############################
### Mac OS building block ###
#############################

ifeq ($(UNAME), $(filter $(UNAME),OSX Darwin))
# do something OSX-y
ifeq ($(MEKD_STANDALONE), Yes)
CXXFLAGS+= -D MEKD_STANDALONE
endif
$(targetlib): $(objects)
	$(CXX) $(CXXFLAGS) -shared -Wl,-install_name,libMEKD_MG.so -o $(targetlib) $(objects) $(LIBFLAGS)


ifeq ($(MEKD_COMPILE_WITH_ROOT), Yes)
CXXFLAGS+= -D MEKD_with_ROOT -std=c++0x
$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) `$(ROOTSYS)/bin/root-config --ldflags` -I $(ROOTINCLUDES) $(target).cpp -o $(target) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)

$(test): $(targetlib)
	$(CXX) $(CXXFLAGS) `$(ROOTSYS)/bin/root-config --ldflags` -I $(ROOTINCLUDES) $(test).cpp -o $(test) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)
else
$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) $(target).cpp -o $(target) -L. $(LIBFLAGS) $(targetlib)

$(test): $(targetlib)
	$(CXX) $(CXXFLAGS) $(test).cpp -o $(test) -L. $(LIBFLAGS) $(targetlib)
endif

endif



.PHONY: clean

clean:
	rm -f $(targetlib); \
	rm -f $(target); \
	rm -f $(test); \
	rm -f $(objects)
