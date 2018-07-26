#ifndef MEKD_MG_cpp
#define MEKD_MG_cpp

/// C++ libraries
#include <iostream>
// #include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>	// for sorting
#include <cmath>


/// CMSSW includes
#ifndef MEKD_STANDALONE
#include "FWCore/ParameterSet/interface/FileInPath.h"
#endif


/// MEs.
#include "MadGraphSrc/DY/MEKD_MEs_DY.h"
#include "MadGraphSrc/Z4l/MEKD_MEs_Z4l.h"
#include "MadGraphSrc/ZZ/MEKD_MEs_ZZ.h"

#include "MadGraphSrc/Spin0/MEKD_MEs_Spin0.h"
#include "MadGraphSrc/Spin1/MEKD_MEs_Spin1.h"
#include "MadGraphSrc/Spin2/MEKD_MEs_Spin2.h"

/// RAW MG5_aMC ME. For testing purposes only.
// #include "MadGraphSrc/RAW/CPPProcess.h"	// ME_RAW




extern "C"
{
#include "Extra_code/MEKD_CalcHEP_PDF.h"
#include "PDFTables/pdt.h"
}

#include "Extra_code/MEKD_CalcHEP_Extra_functions.h"
#include "Extra_code/MEKD_MG_Boosts.h"
#include "higgs_properties/hggeffective.h"
#include "MadGraphSrc/read_slha.h"

#include "../interface/MEKD_MG.h"

using namespace std;


/// Part of pdfreader
extern "C" pdtStr pdtSg, pdtSd, pdtSu, pdtSs, pdtSc,
	pdtSad, pdtSau, pdtSas, pdtSac;

// #define PDTFILE "PDFTables/cteq6l.pdt" // CalCHEP reads a table for CTEQ6L. You can change PDF set as you want.


/// 4l final state (+photon)
qq_Z4l_SIG_DN_OF ME_qq_Z4l_SIG_DownType_OF;
qq_Z4l_SIG_DN_SF ME_qq_Z4l_SIG_DownType_SF;
qq_Z4l_SIG_UP_OF ME_qq_Z4l_SIG_UpType_OF;
qq_Z4l_SIG_UP_SF ME_qq_Z4l_SIG_UpType_SF;
qq_Z4l_SIG_DN_OFpA ME_qq_Z4l_SIG_DownType_OFpA;
qq_Z4l_SIG_DN_SFpA ME_qq_Z4l_SIG_DownType_SFpA;
qq_Z4l_SIG_UP_OFpA ME_qq_Z4l_SIG_UpType_OFpA;
qq_Z4l_SIG_UP_SFpA ME_qq_Z4l_SIG_UpType_SFpA;

qq_Z4l_BKG_DN_OF ME_qq_Z4l_BKG_DownType_OF;
qq_Z4l_BKG_DN_SF ME_qq_Z4l_BKG_DownType_SF;
qq_Z4l_BKG_UP_OF ME_qq_Z4l_BKG_UpType_OF;
qq_Z4l_BKG_UP_SF ME_qq_Z4l_BKG_UpType_SF;
qq_Z4l_BKG_DN_OFpA ME_qq_Z4l_BKG_DownType_OFpA;
qq_Z4l_BKG_DN_SFpA ME_qq_Z4l_BKG_DownType_SFpA;
qq_Z4l_BKG_UP_OFpA ME_qq_Z4l_BKG_UpType_OFpA;
qq_Z4l_BKG_UP_SFpA ME_qq_Z4l_BKG_UpType_SFpA;

qq_ZZ_DN_OF ME_qq_ZZ_DownType_OF;
qq_ZZ_DN_SF ME_qq_ZZ_DownType_SF;
qq_ZZ_UP_OF ME_qq_ZZ_UpType_OF;
qq_ZZ_UP_SF ME_qq_ZZ_UpType_SF;
qq_ZZ_DN_OFpA ME_qq_ZZ_DownType_OFpA;
qq_ZZ_DN_SFpA ME_qq_ZZ_DownType_SFpA;
qq_ZZ_UP_OFpA ME_qq_ZZ_UpType_OFpA;
qq_ZZ_UP_SFpA ME_qq_ZZ_UpType_SFpA;

gg_Spin0_OF ME_Signal_gg_Spin0_OF;
gg_Spin0_SF ME_Signal_gg_Spin0_SF;
gg_Spin0_OFpA ME_Signal_gg_Spin0_OFpA;
gg_Spin0_SFpA ME_Signal_gg_Spin0_SFpA;

// qq_Spin0_DN_OF ME_Signal_qq_Spin0_DownType_OF;
// qq_Spin0_DN_SF ME_Signal_qq_Spin0_DownType_SF;
// qq_Spin0_UP_OF ME_Signal_qq_Spin0_UpType_OF;
// qq_Spin0_UP_SF ME_Signal_qq_Spin0_UpType_SF;
// qq_Spin0_DN_OFpA ME_Signal_qq_Spin0_DownType_OFpA;
// qq_Spin0_DN_SFpA ME_Signal_qq_Spin0_DownType_SFpA;
// qq_Spin0_UP_OFpA ME_Signal_qq_Spin0_UpType_OFpA;
// qq_Spin0_UP_SFpA ME_Signal_qq_Spin0_UpType_SFpA;

qq_Spin1_DN_OF ME_Signal_qq_Spin1_DownType_OF;
qq_Spin1_DN_SF ME_Signal_qq_Spin1_DownType_SF;
qq_Spin1_UP_OF ME_Signal_qq_Spin1_UpType_OF;
qq_Spin1_UP_SF ME_Signal_qq_Spin1_UpType_SF;
qq_Spin1_DN_OFpA ME_Signal_qq_Spin1_DownType_OFpA;
qq_Spin1_DN_SFpA ME_Signal_qq_Spin1_DownType_SFpA;
qq_Spin1_UP_OFpA ME_Signal_qq_Spin1_UpType_OFpA;
qq_Spin1_UP_SFpA ME_Signal_qq_Spin1_UpType_SFpA;

gg_Spin2_OF ME_Signal_gg_Spin2_OF;
gg_Spin2_SF ME_Signal_gg_Spin2_SF;
gg_Spin2_OFpA ME_Signal_gg_Spin2_OFpA;
gg_Spin2_SFpA ME_Signal_gg_Spin2_SFpA;

qq_Spin2_DN_OF ME_Signal_qq_Spin2_DownType_OF;
qq_Spin2_DN_SF ME_Signal_qq_Spin2_DownType_SF;
qq_Spin2_UP_OF ME_Signal_qq_Spin2_UpType_OF;
qq_Spin2_UP_SF ME_Signal_qq_Spin2_UpType_SF;
qq_Spin2_DN_OFpA ME_Signal_qq_Spin2_DownType_OFpA;
qq_Spin2_DN_SFpA ME_Signal_qq_Spin2_DownType_SFpA;
qq_Spin2_UP_OFpA ME_Signal_qq_Spin2_UpType_OFpA;
qq_Spin2_UP_SFpA ME_Signal_qq_Spin2_UpType_SFpA;


/// ZZ and two-fermion couplings. 4l final state (+photon)
gg_Spin0_2f_OF ME_Signal_gg_Spin0_2f_OF;
gg_Spin0_2f_SF ME_Signal_gg_Spin0_2f_SF;
gg_Spin0_2f_OFpA ME_Signal_gg_Spin0_2f_OFpA;
gg_Spin0_2f_SFpA ME_Signal_gg_Spin0_2f_SFpA;

qq_Spin1_2f_DN_OF ME_Signal_qq_Spin1_2f_DownType_OF;
qq_Spin1_2f_DN_SF ME_Signal_qq_Spin1_2f_DownType_SF;
qq_Spin1_2f_UP_OF ME_Signal_qq_Spin1_2f_UpType_OF;
qq_Spin1_2f_UP_SF ME_Signal_qq_Spin1_2f_UpType_SF;
qq_Spin1_2f_DN_OFpA ME_Signal_qq_Spin1_2f_DownType_OFpA;
qq_Spin1_2f_DN_SFpA ME_Signal_qq_Spin1_2f_DownType_SFpA;
qq_Spin1_2f_UP_OFpA ME_Signal_qq_Spin1_2f_UpType_OFpA;
qq_Spin1_2f_UP_SFpA ME_Signal_qq_Spin1_2f_UpType_SFpA;

gg_Spin2_2f_OF ME_Signal_gg_Spin2_2f_OF;
gg_Spin2_2f_SF ME_Signal_gg_Spin2_2f_SF;
gg_Spin2_2f_OFpA ME_Signal_gg_Spin2_2f_OFpA;
gg_Spin2_2f_SFpA ME_Signal_gg_Spin2_2f_SFpA;

qq_Spin2_2f_DN_OF ME_Signal_qq_Spin2_2f_DownType_OF;
qq_Spin2_2f_DN_SF ME_Signal_qq_Spin2_2f_DownType_SF;
qq_Spin2_2f_UP_OF ME_Signal_qq_Spin2_2f_UpType_OF;
qq_Spin2_2f_UP_SF ME_Signal_qq_Spin2_2f_UpType_SF;
qq_Spin2_2f_DN_OFpA ME_Signal_qq_Spin2_2f_DownType_OFpA;
qq_Spin2_2f_DN_SFpA ME_Signal_qq_Spin2_2f_DownType_SFpA;
qq_Spin2_2f_UP_OFpA ME_Signal_qq_Spin2_2f_UpType_OFpA;
qq_Spin2_2f_UP_SFpA ME_Signal_qq_Spin2_2f_UpType_SFpA;


/// 1->4l final state (+photon)
DY_2l ME_DY_2l;
DY_2lpA ME_DY_2lpA;

Spin0_2l ME_Signal_Spin0_2l;
Spin0_OF ME_Signal_Spin0_OF;
Spin0_SF ME_Signal_Spin0_SF;
Spin0_2f_OF ME_Signal_Spin0_2f_OF;
Spin0_2f_SF ME_Signal_Spin0_2f_SF;
Spin0_2lpA ME_Signal_Spin0_2lpA;
Spin0_OFpA ME_Signal_Spin0_OFpA;
Spin0_SFpA ME_Signal_Spin0_SFpA;
Spin0_2f_OFpA ME_Signal_Spin0_2f_OFpA;
Spin0_2f_SFpA ME_Signal_Spin0_2f_SFpA;

Spin1_2l ME_Signal_Spin1_2l;
Spin1_OF ME_Signal_Spin1_OF;
Spin1_SF ME_Signal_Spin1_SF;
Spin1_2f_OF ME_Signal_Spin1_2f_OF;
Spin1_2f_SF ME_Signal_Spin1_2f_SF;
Spin1_2lpA ME_Signal_Spin1_2lpA;
Spin1_OFpA ME_Signal_Spin1_OFpA;
Spin1_SFpA ME_Signal_Spin1_SFpA;
Spin1_2f_OFpA ME_Signal_Spin1_2f_OFpA;
Spin1_2f_SFpA ME_Signal_Spin1_2f_SFpA;

Spin2_2l ME_Signal_Spin2_2l;
Spin2_OF ME_Signal_Spin2_OF;
Spin2_SF ME_Signal_Spin2_SF;
Spin2_2f_OF ME_Signal_Spin2_2f_OF;
Spin2_2f_SF ME_Signal_Spin2_2f_SF;
Spin2_2lpA ME_Signal_Spin2_2lpA;
Spin2_OFpA ME_Signal_Spin2_OFpA;
Spin2_SFpA ME_Signal_Spin2_SFpA;
Spin2_2f_OFpA ME_Signal_Spin2_2f_OFpA;
Spin2_2f_SFpA ME_Signal_Spin2_2f_SFpA;


/// 2mu final state (+photon)
qq_DY_DN_2l ME_qq_DY_DownType_2l;
qq_DY_UP_2l ME_qq_DY_UpType_2l;
qq_DY_DN_2lpA ME_qq_DY_DownType_2lpA;
qq_DY_UP_2lpA ME_qq_DY_UpType_2lpA;

gg_Spin0_2l ME_Signal_gg_Spin0_2l;
gg_Spin0_2lpA ME_Signal_gg_Spin0_2lpA;

// qq_Spin0_DN_2l ME_Signal_qq_Spin0_DownType_2l;
// qq_Spin0_UP_2l ME_Signal_qq_Spin0_UpType_2l;
// qq_Spin0_DN_2lpA ME_Signal_qq_Spin0_DownType_2lpA;
// qq_Spin0_UP_2lpA ME_Signal_qq_Spin0_UpType_2lpA;

qq_Spin1_DN_2l ME_Signal_qq_Spin1_DownType_2l;
qq_Spin1_UP_2l ME_Signal_qq_Spin1_UpType_2l;
qq_Spin1_DN_2lpA ME_Signal_qq_Spin1_DownType_2lpA;
qq_Spin1_UP_2lpA ME_Signal_qq_Spin1_UpType_2lpA;

gg_Spin2_2l ME_Signal_gg_Spin2_2l;
gg_Spin2_2lpA ME_Signal_gg_Spin2_2lpA;

qq_Spin2_DN_2l ME_Signal_qq_Spin2_DownType_2l;
qq_Spin2_UP_2l ME_Signal_qq_Spin2_UpType_2l;
qq_Spin2_DN_2lpA ME_Signal_qq_Spin2_DownType_2lpA;
qq_Spin2_UP_2lpA ME_Signal_qq_Spin2_UpType_2lpA;


/// RAW MG5_aMC ME. For testing purposes only.
// CPPProcess ME_RAW;




MEKD_MG::MEKD_MG()
{
	Mixing_Coefficients_Spin0 = new complex<double>[4];
	Mixing_Coefficients_Spin0_internal = new complex<double>[4];
	Mixing_Coefficients_Spin1 = new complex<double>[8];
	Mixing_Coefficients_Spin1_internal = new complex<double>[8];
	Mixing_Coefficients_Spin2 = new complex<double>[20];
	Mixing_Coefficients_Spin2_internal = new complex<double>[20];
	
	Set_Default_MEKD_MG_Parameters();
	
	
	/// Cross-cheking MEs for consistency. ZZ
	if( ME_qq_ZZ_DownType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_ZZ_DownType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_ZZ_UpType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_ZZ_UpType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_ZZ_DownType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_ZZ_DownType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_ZZ_UpType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_ZZ_UpType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_qq_Z4l_BKG_DownType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_BKG_DownType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_BKG_UpType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_BKG_UpType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_BKG_DownType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_BKG_DownType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_BKG_UpType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_BKG_UpType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_qq_Z4l_SIG_DownType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_SIG_DownType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_SIG_UpType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_SIG_UpType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_SIG_DownType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_SIG_DownType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_SIG_UpType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_Z4l_SIG_UpType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_gg_Spin0_OF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin0_SF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin0_OFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin0_SFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
// 	if( ME_Signal_qq_Spin0_DownType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
// 	if( ME_Signal_qq_Spin0_DownType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
// 	if( ME_Signal_qq_Spin0_UpType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
// 	if( ME_Signal_qq_Spin0_UpType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
// 	if( ME_Signal_qq_Spin0_DownType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
// 	if( ME_Signal_qq_Spin0_DownType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
// 	if( ME_Signal_qq_Spin0_UpType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
// 	if( ME_Signal_qq_Spin0_UpType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_qq_Spin1_DownType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_DownType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_UpType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_UpType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_DownType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_DownType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_UpType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_UpType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_gg_Spin2_OF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin2_SF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin2_OFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin2_SFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_qq_Spin2_DownType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_DownType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_UpType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_UpType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_DownType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_DownType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_UpType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_UpType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	
	/// Cross-cheking MEs for consistency. ZZ and two-fermion couplings
	if( ME_Signal_gg_Spin0_2f_OF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin0_2f_SF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin0_2f_OFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin0_2f_SFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_qq_Spin1_2f_DownType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_2f_DownType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_2f_UpType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_2f_UpType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_2f_DownType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_2f_DownType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_2f_UpType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_2f_UpType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_gg_Spin2_2f_OF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin2_2f_SF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin2_2f_OFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin2_2f_SFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_qq_Spin2_2f_DownType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_2f_DownType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_2f_UpType_OF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_2f_UpType_SF.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_2f_DownType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_2f_DownType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_2f_UpType_OFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_2f_UpType_SFpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	
	///  Cross-cheking MEs for consistency. ZZ, no initial state
	if( ME_DY_2l.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_DY_2lpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_Spin0_2l.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin0_OF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin0_SF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin0_2f_OF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin0_2f_SF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin0_2lpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin0_OFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin0_SFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin0_2f_OFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin0_2f_SFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_Spin1_OF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin1_SF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin1_2f_OF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin1_2f_SF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin1_OFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin1_SFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin1_2f_OFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin1_2f_SFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_Spin2_OF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin2_SF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin2_2f_OF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin2_2f_SF.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin2_OFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin2_SFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin2_2f_OFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_Spin2_2f_SFpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	
	/// Cross-cheking MEs for consistency. 2mu
	if( ME_qq_DY_DownType_2l.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_DY_UpType_2l.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_DY_DownType_2lpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_qq_DY_UpType_2lpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_gg_Spin0_2l.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin0_2lpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
// 	if( ME_Signal_qq_Spin0_DownType_2l.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
// 	if( ME_Signal_qq_Spin0_UpType_2l.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
// 	if( ME_Signal_qq_Spin0_DownType_2lpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
// 	if( ME_Signal_qq_Spin0_UpType_2lpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_qq_Spin1_DownType_2l.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_UpType_2l.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_DownType_2lpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin1_UpType_2lpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_gg_Spin2_2l.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_gg_Spin2_2lpA.nprocesses!=1 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	if( ME_Signal_qq_Spin2_DownType_2l.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_UpType_2l.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_DownType_2lpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	if( ME_Signal_qq_Spin2_UpType_2lpA.nprocesses!=2 ) { cerr << "Problem in ME class detected. Exiting.\n"; exit(1); }
	
	
	p_set.push_back( new double[4] );
	p_set.push_back( new double[4] );
	p_set.push_back( new double[4] );
	p_set.push_back( new double[4] );
	p_set.push_back( new double[4] );
	p_set.push_back( new double[4] );
	p_set.push_back( new double[4] );	// a photon comes here, otherwise, unused
	
	p1 = new double[4];
	p2 = new double[4];
	p3 = new double[4];
	p4 = new double[4];
	p5 = new double[4];
	
	id1 = 10000;
	id2 = 10000;
	id3 = 10000;
	id4 = 10000;
	id5 = 10000;
	
	id_set.push_back( id1 );
	id_set.push_back( id2 );
	id_set.push_back( id3 );
	id_set.push_back( id4 );
	id_set.push_back( id5 );
	
	pl1_internal = NULL;
	pl2_internal = NULL;
	pl3_internal = NULL;
	pl4_internal = NULL;
	pA1_internal = NULL;
	
	Parameters_Are_Loaded = false;
	Predefined_Model = false;
}



MEKD_MG::~MEKD_MG()
{
	delete Mixing_Coefficients_Spin0;
	delete Mixing_Coefficients_Spin0_internal;
	delete Mixing_Coefficients_Spin1;
	delete Mixing_Coefficients_Spin1_internal;
	delete Mixing_Coefficients_Spin2;
	delete Mixing_Coefficients_Spin2_internal;
	
	if( Parameters_Are_Loaded ) Unload_pdfreader();
	
	p_set.clear();
	id_set.clear();
}



int MEKD_MG::Load_Parameters()
{
	Set_Of_Model_Parameters.read_slha_file( Parameter_file );
	
	/// Initializing parameters
	// ZZ
	ME_qq_ZZ_UpType_SF.initProc( Parameter_file );
	ME_qq_ZZ_UpType_OF.initProc( Parameter_file );
	ME_qq_ZZ_DownType_SF.initProc( Parameter_file );
	ME_qq_ZZ_DownType_OF.initProc( Parameter_file );
	ME_qq_ZZ_UpType_SFpA.initProc( Parameter_file );
	ME_qq_ZZ_UpType_OFpA.initProc( Parameter_file );
	ME_qq_ZZ_DownType_SFpA.initProc( Parameter_file );
	ME_qq_ZZ_DownType_OFpA.initProc( Parameter_file );
	
	ME_qq_Z4l_BKG_UpType_SF.initProc( Parameter_file );
	ME_qq_Z4l_BKG_UpType_OF.initProc( Parameter_file );
	ME_qq_Z4l_BKG_DownType_SF.initProc( Parameter_file );
	ME_qq_Z4l_BKG_DownType_OF.initProc( Parameter_file );
	ME_qq_Z4l_BKG_UpType_SFpA.initProc( Parameter_file );
	ME_qq_Z4l_BKG_UpType_OFpA.initProc( Parameter_file );
	ME_qq_Z4l_BKG_DownType_SFpA.initProc( Parameter_file );
	ME_qq_Z4l_BKG_DownType_OFpA.initProc( Parameter_file );
	
	ME_qq_Z4l_SIG_UpType_SF.initProc( Parameter_file );
	ME_qq_Z4l_SIG_UpType_OF.initProc( Parameter_file );
	ME_qq_Z4l_SIG_DownType_SF.initProc( Parameter_file );
	ME_qq_Z4l_SIG_DownType_OF.initProc( Parameter_file );
	ME_qq_Z4l_SIG_UpType_SFpA.initProc( Parameter_file );
	ME_qq_Z4l_SIG_UpType_OFpA.initProc( Parameter_file );
	ME_qq_Z4l_SIG_DownType_SFpA.initProc( Parameter_file );
	ME_qq_Z4l_SIG_DownType_OFpA.initProc( Parameter_file );
	
	ME_Signal_gg_Spin0_SF.initProc( Parameter_file );
	ME_Signal_gg_Spin0_OF.initProc( Parameter_file );
	ME_Signal_gg_Spin0_SFpA.initProc( Parameter_file );
	ME_Signal_gg_Spin0_OFpA.initProc( Parameter_file );
	
// 	ME_Signal_qq_Spin0_DownType_SF.initProc( Parameter_file );
// 	ME_Signal_qq_Spin0_DownType_OF.initProc( Parameter_file );
// 	ME_Signal_qq_Spin0_UpType_SF.initProc( Parameter_file );
// 	ME_Signal_qq_Spin0_UpType_OF.initProc( Parameter_file );
// 	ME_Signal_qq_Spin0_DownType_SFpA.initProc( Parameter_file );
// 	ME_Signal_qq_Spin0_DownType_OFpA.initProc( Parameter_file );
// 	ME_Signal_qq_Spin0_UpType_SFpA.initProc( Parameter_file );
// 	ME_Signal_qq_Spin0_UpType_OFpA.initProc( Parameter_file );
	
	ME_Signal_qq_Spin1_DownType_SF.initProc( Parameter_file );
	ME_Signal_qq_Spin1_DownType_OF.initProc( Parameter_file );
	ME_Signal_qq_Spin1_UpType_SF.initProc( Parameter_file );
	ME_Signal_qq_Spin1_UpType_OF.initProc( Parameter_file );
	ME_Signal_qq_Spin1_DownType_SFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin1_DownType_OFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin1_UpType_SFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin1_UpType_OFpA.initProc( Parameter_file );
	
	ME_Signal_gg_Spin2_SF.initProc( Parameter_file );
	ME_Signal_gg_Spin2_OF.initProc( Parameter_file );
	ME_Signal_gg_Spin2_SFpA.initProc( Parameter_file );
	ME_Signal_gg_Spin2_OFpA.initProc( Parameter_file );
	
	ME_Signal_qq_Spin2_DownType_SF.initProc( Parameter_file );
	ME_Signal_qq_Spin2_DownType_OF.initProc( Parameter_file );
	ME_Signal_qq_Spin2_UpType_SF.initProc( Parameter_file );
	ME_Signal_qq_Spin2_UpType_OF.initProc( Parameter_file );
	ME_Signal_qq_Spin2_DownType_SFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin2_DownType_OFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin2_UpType_SFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin2_UpType_OFpA.initProc( Parameter_file );
	
	
	// ZZ and two-fermion couplings
	ME_Signal_gg_Spin0_2f_SF.initProc( Parameter_file );
	ME_Signal_gg_Spin0_2f_OF.initProc( Parameter_file );
	ME_Signal_gg_Spin0_2f_SFpA.initProc( Parameter_file );
	ME_Signal_gg_Spin0_2f_OFpA.initProc( Parameter_file );
	
	ME_Signal_qq_Spin1_2f_DownType_SF.initProc( Parameter_file );
	ME_Signal_qq_Spin1_2f_DownType_OF.initProc( Parameter_file );
	ME_Signal_qq_Spin1_2f_UpType_SF.initProc( Parameter_file );
	ME_Signal_qq_Spin1_2f_UpType_OF.initProc( Parameter_file );
	ME_Signal_qq_Spin1_2f_DownType_SFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin1_2f_DownType_OFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin1_2f_UpType_SFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin1_2f_UpType_OFpA.initProc( Parameter_file );
	
	ME_Signal_gg_Spin2_2f_SF.initProc( Parameter_file );
	ME_Signal_gg_Spin2_2f_OF.initProc( Parameter_file );
	ME_Signal_gg_Spin2_2f_SFpA.initProc( Parameter_file );
	ME_Signal_gg_Spin2_2f_OFpA.initProc( Parameter_file );
	
	ME_Signal_qq_Spin2_2f_DownType_SF.initProc( Parameter_file );
	ME_Signal_qq_Spin2_2f_DownType_OF.initProc( Parameter_file );
	ME_Signal_qq_Spin2_2f_UpType_SF.initProc( Parameter_file );
	ME_Signal_qq_Spin2_2f_UpType_OF.initProc( Parameter_file );
	ME_Signal_qq_Spin2_2f_DownType_SFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin2_2f_DownType_OFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin2_2f_UpType_SFpA.initProc( Parameter_file );
	ME_Signal_qq_Spin2_2f_UpType_OFpA.initProc( Parameter_file );
	
	
	// ZZ, no itial state
	ME_DY_2l.initProc( Parameter_file );
	ME_DY_2lpA.initProc( Parameter_file );
	
	ME_Signal_Spin0_2l.initProc( Parameter_file );
	ME_Signal_Spin0_SF.initProc( Parameter_file );
	ME_Signal_Spin0_OF.initProc( Parameter_file );
	ME_Signal_Spin0_2f_SF.initProc( Parameter_file );
	ME_Signal_Spin0_2f_OF.initProc( Parameter_file );
	ME_Signal_Spin0_2lpA.initProc( Parameter_file );
	ME_Signal_Spin0_SFpA.initProc( Parameter_file );
	ME_Signal_Spin0_OFpA.initProc( Parameter_file );
	ME_Signal_Spin0_2f_SFpA.initProc( Parameter_file );
	ME_Signal_Spin0_2f_OFpA.initProc( Parameter_file );
	
	ME_Signal_Spin1_2l.initProc( Parameter_file );
	ME_Signal_Spin1_SF.initProc( Parameter_file );
	ME_Signal_Spin1_OF.initProc( Parameter_file );
	ME_Signal_Spin1_2f_SF.initProc( Parameter_file );
	ME_Signal_Spin1_2f_OF.initProc( Parameter_file );
	ME_Signal_Spin1_2lpA.initProc( Parameter_file );
	ME_Signal_Spin1_SFpA.initProc( Parameter_file );
	ME_Signal_Spin1_OFpA.initProc( Parameter_file );
	ME_Signal_Spin1_2f_SFpA.initProc( Parameter_file );
	ME_Signal_Spin1_2f_OFpA.initProc( Parameter_file );
	
	ME_Signal_Spin2_2l.initProc( Parameter_file );
	ME_Signal_Spin2_SF.initProc( Parameter_file );
	ME_Signal_Spin2_OF.initProc( Parameter_file );
	ME_Signal_Spin2_2f_SF.initProc( Parameter_file );
	ME_Signal_Spin2_2f_OF.initProc( Parameter_file );
	ME_Signal_Spin2_2lpA.initProc( Parameter_file );
	ME_Signal_Spin2_SFpA.initProc( Parameter_file );
	ME_Signal_Spin2_OFpA.initProc( Parameter_file );
	ME_Signal_Spin2_2f_SFpA.initProc( Parameter_file );
	ME_Signal_Spin2_2f_OFpA.initProc( Parameter_file );
	
	
	// 2mu
	ME_qq_DY_UpType_2l.initProc( Parameter_file );
	ME_qq_DY_DownType_2l.initProc( Parameter_file );
	ME_qq_DY_UpType_2lpA.initProc( Parameter_file );
	ME_qq_DY_DownType_2lpA.initProc( Parameter_file );
	
	ME_Signal_gg_Spin0_2l.initProc( Parameter_file );
	ME_Signal_gg_Spin0_2lpA.initProc( Parameter_file );
	
// 	ME_Signal_qq_Spin0_DownType_2l.initProc( Parameter_file );
// 	ME_Signal_qq_Spin0_UpType_2l.initProc( Parameter_file );
// 	ME_Signal_qq_Spin0_DownType_2lpA.initProc( Parameter_file );
// 	ME_Signal_qq_Spin0_UpType_2lpA.initProc( Parameter_file );
	
	ME_Signal_qq_Spin1_DownType_2l.initProc( Parameter_file );
	ME_Signal_qq_Spin1_UpType_2l.initProc( Parameter_file );
	ME_Signal_qq_Spin1_DownType_2lpA.initProc( Parameter_file );
	ME_Signal_qq_Spin1_UpType_2lpA.initProc( Parameter_file );
	
	ME_Signal_gg_Spin2_2l.initProc( Parameter_file );
	ME_Signal_gg_Spin2_2lpA.initProc( Parameter_file );
	
	ME_Signal_qq_Spin2_DownType_2l.initProc( Parameter_file );
	ME_Signal_qq_Spin2_UpType_2l.initProc( Parameter_file );
	ME_Signal_qq_Spin2_DownType_2lpA.initProc( Parameter_file );
	ME_Signal_qq_Spin2_UpType_2lpA.initProc( Parameter_file );
	
	
	/// RAW MG5_aMC ME. For testing purposes only.
// 	ME_RAW.initProc( "../src/MadGraphSrc/RAW/param_card.dat" );
	
	
	params_m_d = Set_Of_Model_Parameters.get_block_entry( "mass", 1, 0 ).real();
	params_m_u = Set_Of_Model_Parameters.get_block_entry( "mass", 2, 0 ).real();
	params_m_s = Set_Of_Model_Parameters.get_block_entry( "mass", 3, 0 ).real();
	params_m_c = Set_Of_Model_Parameters.get_block_entry( "mass", 4, 0 ).real();
	params_m_e = Set_Of_Model_Parameters.get_block_entry( "mass", 11, 0 ).real();
	params_m_mu = Set_Of_Model_Parameters.get_block_entry( "mass", 13, 0 ).real();
	params_m_Z = Set_Of_Model_Parameters.get_block_entry( "mass", 23, 9.11876e+01 ).real();
	
	params_rhou11 = Set_Of_Model_Parameters.get_block_entry( "heff", 9, 0 );
	params_rhou12 = Set_Of_Model_Parameters.get_block_entry( "heff", 10, 0 );
	params_rhoc11 = Set_Of_Model_Parameters.get_block_entry( "heff", 11, 0 );
	params_rhoc12 = Set_Of_Model_Parameters.get_block_entry( "heff", 12, 0 );
	params_rhod11 = Set_Of_Model_Parameters.get_block_entry( "heff", 13, 0 );
	params_rhod12 = Set_Of_Model_Parameters.get_block_entry( "heff", 14, 0 );
	params_rhos11 = Set_Of_Model_Parameters.get_block_entry( "heff", 15, 0 );
	params_rhos12 = Set_Of_Model_Parameters.get_block_entry( "heff", 16, 0 );
	params_rhob11 = Set_Of_Model_Parameters.get_block_entry( "heff", 17, 0 );
	params_rhob12 = Set_Of_Model_Parameters.get_block_entry( "heff", 18, 0 );
	
	params_rhou11 = Set_Of_Model_Parameters.get_block_entry( "vec", 3, 0 );
	params_rhou12 = Set_Of_Model_Parameters.get_block_entry( "vec", 4, 0 );
	params_rhou13 = Set_Of_Model_Parameters.get_block_entry( "vec", 5, 0 );
	params_rhou14 = Set_Of_Model_Parameters.get_block_entry( "vec", 6, 0 );
	params_rhoc11 = Set_Of_Model_Parameters.get_block_entry( "vec", 7, 0 );
	params_rhoc12 = Set_Of_Model_Parameters.get_block_entry( "vec", 8, 0 );
	params_rhoc13 = Set_Of_Model_Parameters.get_block_entry( "vec", 9, 0 );
	params_rhoc14 = Set_Of_Model_Parameters.get_block_entry( "vec", 10, 0 );
	params_rhod11 = Set_Of_Model_Parameters.get_block_entry( "vec", 11, 0 );
	params_rhod12 = Set_Of_Model_Parameters.get_block_entry( "vec", 12, 0 );
	params_rhod13 = Set_Of_Model_Parameters.get_block_entry( "vec", 13, 0 );
	params_rhod14 = Set_Of_Model_Parameters.get_block_entry( "vec", 14, 0 );
	params_rhos11 = Set_Of_Model_Parameters.get_block_entry( "vec", 15, 0 );
	params_rhos12 = Set_Of_Model_Parameters.get_block_entry( "vec", 16, 0 );
	params_rhos13 = Set_Of_Model_Parameters.get_block_entry( "vec", 17, 0 );
	params_rhos14 = Set_Of_Model_Parameters.get_block_entry( "vec", 18, 0 );
	params_rhob11 = Set_Of_Model_Parameters.get_block_entry( "vec", 19, 0 );
	params_rhob12 = Set_Of_Model_Parameters.get_block_entry( "vec", 20, 0 );
	params_rhob13 = Set_Of_Model_Parameters.get_block_entry( "vec", 21, 0 );
	params_rhob14 = Set_Of_Model_Parameters.get_block_entry( "vec", 22, 0 );
	
	params_rhou21 = Set_Of_Model_Parameters.get_block_entry( "gravity", 21, 0 );
	params_rhou22 = Set_Of_Model_Parameters.get_block_entry( "gravity", 22, 0 );
	params_rhou23 = Set_Of_Model_Parameters.get_block_entry( "gravity", 23, 0 );
	params_rhou24 = Set_Of_Model_Parameters.get_block_entry( "gravity", 24, 0 );
	params_rhoc21 = Set_Of_Model_Parameters.get_block_entry( "gravity", 25, 0 );
	params_rhoc22 = Set_Of_Model_Parameters.get_block_entry( "gravity", 26, 0 );
	params_rhoc23 = Set_Of_Model_Parameters.get_block_entry( "gravity", 27, 0 );
	params_rhoc24 = Set_Of_Model_Parameters.get_block_entry( "gravity", 28, 0 );
	params_rhod21 = Set_Of_Model_Parameters.get_block_entry( "gravity", 29, 0 );
	params_rhod22 = Set_Of_Model_Parameters.get_block_entry( "gravity", 30, 0 );
	params_rhod23 = Set_Of_Model_Parameters.get_block_entry( "gravity", 31, 0 );
	params_rhod24 = Set_Of_Model_Parameters.get_block_entry( "gravity", 32, 0 );
	params_rhos21 = Set_Of_Model_Parameters.get_block_entry( "gravity", 33, 0 );
	params_rhos22 = Set_Of_Model_Parameters.get_block_entry( "gravity", 34, 0 );
	params_rhos23 = Set_Of_Model_Parameters.get_block_entry( "gravity", 35, 0 );
	params_rhos24 = Set_Of_Model_Parameters.get_block_entry( "gravity", 36, 0 );
	params_rhob21 = Set_Of_Model_Parameters.get_block_entry( "gravity", 37, 0 );
	params_rhob22 = Set_Of_Model_Parameters.get_block_entry( "gravity", 38, 0 );
	params_rhob23 = Set_Of_Model_Parameters.get_block_entry( "gravity", 39, 0 );
	params_rhob24 = Set_Of_Model_Parameters.get_block_entry( "gravity", 40, 0 );
	
	
	v_expectation = 1.0/sqrt( sqrt(2)*Set_Of_Model_Parameters.get_block_entry( "sminputs", 2, 1.166370e-05 ).real() );
	hZZ_coupling = 2.0*params_m_Z*params_m_Z/v_expectation;
	
	
	Load_pdfreader( const_cast<char*>(PDF_file.c_str()) );
	
	Parameters_Are_Loaded = true;
	return 0;
}



int MEKD_MG::Reload_Parameters()
{
	if( !Parameters_Are_Loaded ) return 1;
	
	Set_Of_Model_Parameters.read_slha_file( static_cast<string>(Parameter_file) );
	
	params_m_d = Set_Of_Model_Parameters.get_block_entry( "mass", 1, 0 ).real();
	params_m_u = Set_Of_Model_Parameters.get_block_entry( "mass", 2, 0 ).real();
	params_m_s = Set_Of_Model_Parameters.get_block_entry( "mass", 3, 0 ).real();
	params_m_c = Set_Of_Model_Parameters.get_block_entry( "mass", 4, 0 ).real();
	params_m_e = Set_Of_Model_Parameters.get_block_entry( "mass", 11, 0 ).real();
	params_m_mu = Set_Of_Model_Parameters.get_block_entry( "mass", 13, 0 ).real();
	params_m_Z = Set_Of_Model_Parameters.get_block_entry( "mass", 23, 9.11876e+01 ).real();
	
	params_rhou11 = Set_Of_Model_Parameters.get_block_entry( "vec", 3, 0 );
	params_rhou12 = Set_Of_Model_Parameters.get_block_entry( "vec", 4, 0 );
	params_rhou13 = Set_Of_Model_Parameters.get_block_entry( "vec", 5, 0 );
	params_rhou14 = Set_Of_Model_Parameters.get_block_entry( "vec", 6, 0 );
	params_rhoc11 = Set_Of_Model_Parameters.get_block_entry( "vec", 7, 0 );
	params_rhoc12 = Set_Of_Model_Parameters.get_block_entry( "vec", 8, 0 );
	params_rhoc13 = Set_Of_Model_Parameters.get_block_entry( "vec", 9, 0 );
	params_rhoc14 = Set_Of_Model_Parameters.get_block_entry( "vec", 10, 0 );
	params_rhod11 = Set_Of_Model_Parameters.get_block_entry( "vec", 11, 0 );
	params_rhod12 = Set_Of_Model_Parameters.get_block_entry( "vec", 12, 0 );
	params_rhod13 = Set_Of_Model_Parameters.get_block_entry( "vec", 13, 0 );
	params_rhod14 = Set_Of_Model_Parameters.get_block_entry( "vec", 14, 0 );
	params_rhos11 = Set_Of_Model_Parameters.get_block_entry( "vec", 15, 0 );
	params_rhos12 = Set_Of_Model_Parameters.get_block_entry( "vec", 16, 0 );
	params_rhos13 = Set_Of_Model_Parameters.get_block_entry( "vec", 17, 0 );
	params_rhos14 = Set_Of_Model_Parameters.get_block_entry( "vec", 18, 0 );
	params_rhob11 = Set_Of_Model_Parameters.get_block_entry( "vec", 19, 0 );
	params_rhob12 = Set_Of_Model_Parameters.get_block_entry( "vec", 20, 0 );
	params_rhob13 = Set_Of_Model_Parameters.get_block_entry( "vec", 21, 0 );
	params_rhob14 = Set_Of_Model_Parameters.get_block_entry( "vec", 22, 0 );
	
	params_rhou11 = Set_Of_Model_Parameters.get_block_entry( "heff", 9, 0 );
	params_rhou12 = Set_Of_Model_Parameters.get_block_entry( "heff", 10, 0 );
	params_rhoc11 = Set_Of_Model_Parameters.get_block_entry( "heff", 11, 0 );
	params_rhoc12 = Set_Of_Model_Parameters.get_block_entry( "heff", 12, 0 );
	params_rhod11 = Set_Of_Model_Parameters.get_block_entry( "heff", 13, 0 );
	params_rhod12 = Set_Of_Model_Parameters.get_block_entry( "heff", 14, 0 );
	params_rhos11 = Set_Of_Model_Parameters.get_block_entry( "heff", 15, 0 );
	params_rhos12 = Set_Of_Model_Parameters.get_block_entry( "heff", 16, 0 );
	params_rhob11 = Set_Of_Model_Parameters.get_block_entry( "heff", 17, 0 );
	params_rhob12 = Set_Of_Model_Parameters.get_block_entry( "heff", 18, 0 );
	
	params_rhou21 = Set_Of_Model_Parameters.get_block_entry( "gravity", 21, 0 );
	params_rhou22 = Set_Of_Model_Parameters.get_block_entry( "gravity", 22, 0 );
	params_rhou23 = Set_Of_Model_Parameters.get_block_entry( "gravity", 23, 0 );
	params_rhou24 = Set_Of_Model_Parameters.get_block_entry( "gravity", 24, 0 );
	params_rhoc21 = Set_Of_Model_Parameters.get_block_entry( "gravity", 25, 0 );
	params_rhoc22 = Set_Of_Model_Parameters.get_block_entry( "gravity", 26, 0 );
	params_rhoc23 = Set_Of_Model_Parameters.get_block_entry( "gravity", 27, 0 );
	params_rhoc24 = Set_Of_Model_Parameters.get_block_entry( "gravity", 28, 0 );
	params_rhod21 = Set_Of_Model_Parameters.get_block_entry( "gravity", 29, 0 );
	params_rhod22 = Set_Of_Model_Parameters.get_block_entry( "gravity", 30, 0 );
	params_rhod23 = Set_Of_Model_Parameters.get_block_entry( "gravity", 31, 0 );
	params_rhod24 = Set_Of_Model_Parameters.get_block_entry( "gravity", 32, 0 );
	params_rhos21 = Set_Of_Model_Parameters.get_block_entry( "gravity", 33, 0 );
	params_rhos22 = Set_Of_Model_Parameters.get_block_entry( "gravity", 34, 0 );
	params_rhos23 = Set_Of_Model_Parameters.get_block_entry( "gravity", 35, 0 );
	params_rhos24 = Set_Of_Model_Parameters.get_block_entry( "gravity", 36, 0 );
	params_rhob21 = Set_Of_Model_Parameters.get_block_entry( "gravity", 37, 0 );
	params_rhob22 = Set_Of_Model_Parameters.get_block_entry( "gravity", 38, 0 );
	params_rhob23 = Set_Of_Model_Parameters.get_block_entry( "gravity", 39, 0 );
	params_rhob24 = Set_Of_Model_Parameters.get_block_entry( "gravity", 40, 0 );
	
	v_expectation = 1.0/sqrt( sqrt(2)*Set_Of_Model_Parameters.get_block_entry( "sminputs", 2, 1.166370e-05 ).real() );
	hZZ_coupling = 2.0*params_m_Z*params_m_Z/v_expectation;
	
	
	Unload_pdfreader();
	Load_pdfreader( const_cast<char*>(PDF_file.c_str()) );
	
	return 0;
}



void MEKD_MG::Set_Default_MEKD_MG_Parameters()
{
	Boost_To_CM = true;	// for a boosted data
	Debug_Mode = false;	// Enable debugging mode
	Fix_Spin0_Production = true;	// use the SM Higgs production mechanism
	Fix_Spin1_Production = false;	// use the a hybrid production mechanism
// 	Force_g3_running = false;	// unused. At some point was included for alpha_QCD
	Overwrite_e_and_mu_masses = false;	// switch for manual m_e, m_mu masses
	Use_mh_eq_m4l = true;	// Set mh to m4l for every event
	Use_mZ4l_eq_m4l = true;	// Set m_Z to m4l for Z4l events
	Use_Higgs_width = true;	//	if false, width is fixed to =1
	Use_PDF_w_pT0 = false;	// Use PDFs in the pT=0 frame. If true, Boost_To_CM is ignored
	Vary_resonance_width = true;	// Allow width to be varied with mass
	Vary_signal_couplings = true;	// Allow couplings to change with mass
	Warning_Mode = true;	// Print warnings
	
	ContributionCoeff_d = 0;	//42	/// the value has no effect if PDF is used but the variable is always used
	ContributionCoeff_u = 1;	//217
	ContributionCoeff_s = 0;	//5
	ContributionCoeff_c = 0;	//3
// 	GG=0;	// Assign QCD coupling, force g3 running if needed
	Sqrt_s = 8000;	//Max energy, collision energy
	
	Mixing_Coefficients_Spin0[0] = complex<double>(1/sqrt(4), 0);
	Mixing_Coefficients_Spin0[1] = complex<double>(1/sqrt(4), 0);
	Mixing_Coefficients_Spin0[2] = complex<double>(1/sqrt(4), 0);
	Mixing_Coefficients_Spin0[3] = complex<double>(1/sqrt(4), 0);
	
	Mixing_Coefficients_Spin1[0] = complex<double>(1/sqrt(4), 0);
	Mixing_Coefficients_Spin1[1] = complex<double>(1/sqrt(4), 0);
	Mixing_Coefficients_Spin1[2] = complex<double>(1/sqrt(4), 0);
	Mixing_Coefficients_Spin1[3] = complex<double>(1/sqrt(4), 0);
	Mixing_Coefficients_Spin1[4] = complex<double>(1/sqrt(2), 0);
	Mixing_Coefficients_Spin1[5] = complex<double>(1/sqrt(2), 0);
	
	Electron_mass = 0;	//0.0005109989, for enabled overwriting
	Higgs_mass = 126;	// Works only if Use_mh_eq_m4l=false
	Higgs_width = 5.753088e-03;	// Practically not used, for future implementations
	Muon_mass = 0;	//0.10565837, for enabled overwriting
	Proton_mass = 0.93827205;	// Always used if needed
	
	Final_state = "2e2m";	// Final state, for the moment: 4e, 4mu, 2e2mu
	Resonance_decay_mode = "ZZ";	// default: ZZ. Alternatives: 2l, 2l_s
	Test_Model = "ggSpin0Pm";	// Models: ZZ, DY, Custom, ggSpin0Pm, ggSpin0M, ggSpin0Ph, qqSpin1P, qqSpin1M, ggSpin2Pm, ggSpin2Ph, ggSpin2Mh, ggSpin2Pb, qqSpin2Pm, qqSpin2Ph, qqSpin2Mh, qqSpin2Pb, Spin0Pm, Spin0M, Spin0Ph, Spin1P, Spin1M, Spin2Pm, Spin2Ph, Spin2Mh, Spin2Pb
	
#ifndef MEKD_STANDALONE
	string inputParameterFile = "ZZMatrixElement/MEKD/src/Cards/param_card.dat";	// Location where a parameter card is stored
	string inputPDFFile = "ZZMatrixElement/MEKD/src/PDFTables/cteq6l.pdt";	// PDF/PDT table file
	edm::FileInPath parameterFileWithFullPath(inputParameterFile);
	edm::FileInPath pdfFileWithFullPath(inputPDFFile);
	Parameter_file = parameterFileWithFullPath.fullPath();
	PDF_file = pdfFileWithFullPath.fullPath();
#else
	Parameter_file = "../src/Cards/param_card.dat";	// Location where a parameter card is stored
	PDF_file = "../src/PDFTables/cteq6l.pdt";	// PDF/PDT table file
#endif
}



int MEKD_MG::Run_MEKD_MG()
{
	if( !Parameters_Are_Loaded ) Load_Parameters();
	if( Arrange_Internal_pls() == 1 ) { cerr << "Particle id error. Exiting.\n"; exit(1); }
	
	double CollisionE;
	
	PDFx1 = 0;
	PDFx2 = 0;
	Background_ME = 0;
	Signal_ME = 0;
	
	if( Overwrite_e_and_mu_masses )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 11, Electron_mass );
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, Muon_mass );
		params_m_e = Electron_mass;
		params_m_mu = Muon_mass;
	}
	
	if( Final_state == "4e" || Final_state == "4eA" )
		{ ml1=Set_Of_Model_Parameters.get_block_entry( "mass", 11, Electron_mass ).real(); ml2=ml1; ml3=ml1; ml4=ml1; }
	if( Final_state == "4m" || Final_state == "4mu" || Final_state == "4mA" || Final_state == "4muA" )
		{ ml1=Set_Of_Model_Parameters.get_block_entry( "mass", 13, Muon_mass ).real(); ml2=ml1; ml3=ml1; ml4=ml1; }
	if( Final_state == "2e2m" || Final_state == "2e2mu" || Final_state == "2e2mA" || Final_state == "2e2muA" )
		{ ml1=Set_Of_Model_Parameters.get_block_entry( "mass", 11, Electron_mass ).real(); ml2=ml1; ml3=Set_Of_Model_Parameters.get_block_entry( "mass", 13, Muon_mass ).real(); ml4=ml3; }
		
	if( Final_state == "2m" || Final_state == "2mu" || Final_state == "2mA" || Final_state == "2muA" )
		{ ml1=Set_Of_Model_Parameters.get_block_entry( "mass", 13, Muon_mass ).real(); ml2=ml1; ml3=0; ml4=0; }
	
	
	/// No boosting setup for initial partons
	for( int i=0; i<4; i++ )
	{
		p_set[0][i] = 0;
		p_set[1][i] = 0;
		
		if( pl1_internal == NULL ) p_set[2][i] = 0;
		else p_set[2][i] = pl1_internal[i];
		if( pl2_internal == NULL ) p_set[3][i] = 0;
		else p_set[3][i] = pl2_internal[i];
		if( pl3_internal == NULL ) p_set[4][i] = 0;
		else p_set[4][i] = pl3_internal[i];
		if( pl4_internal == NULL ) p_set[5][i] = 0;
		else p_set[5][i] = pl4_internal[i];
		
		// Adaptive photon handling
		if( pA1_internal == NULL ) p_set[6][i] = 0;
		else
		{
			if( Final_state == "2mA" || Final_state == "2muA" )
			{
				p_set[4][i] = pA1_internal[i];
				p_set[6][i] = 0;
			}
			else p_set[6][i] = pA1_internal[i];
		}
	}
	
	PDFx1 = ( (p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]+p_set[6][0]) + (p_set[2][3]+p_set[3][3]+p_set[4][3]+p_set[5][3]+p_set[6][3]) )/Sqrt_s;
	PDFx2 = ( (p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]+p_set[6][0]) - (p_set[2][3]+p_set[3][3]+p_set[4][3]+p_set[5][3]+p_set[6][3]) )/Sqrt_s;
	
	// 0 mass approximation
	p_set[0][0] = 0.5*PDFx1*Sqrt_s;
	p_set[0][1] = 0;
	p_set[0][2] = 0;
	p_set[0][3] = 0.5*PDFx1*Sqrt_s;	// to be updated

	p_set[1][0] = 0.5*PDFx2*Sqrt_s;
	p_set[1][1] = 0;
	p_set[1][2] = 0;
	p_set[1][3] = -0.5*PDFx2*Sqrt_s;	// to be updated
	
	
	/// Calculate values needed for the PDF in the pT=0 frame
	if( Use_PDF_w_pT0 )
	{
		p_set[0][0] = Sqrt_s/2;
		p_set[1][0] = Sqrt_s/2;
// 		p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - Proton_mass*Proton_mass );
// 		p_set[1][3] = -p_set[0][3];
// 		/// Using forward-beam approximation
// 		Boost_4p_and_2p_2_pT0( ml1, p_set[2], ml2, p_set[3], ml3, p_set[4], ml4, p_set[5], Proton_mass, p_set[0], Proton_mass, p_set[1] );
// 		Boost_4p_2_pT0( ml1, p_set[2], ml2, p_set[3], ml3, p_set[4], ml4, p_set[5] );
		Boost_5p_2_pT0( ml1, p_set[2], ml2, p_set[3], ml3, p_set[4], ml4, p_set[5], 0, p_set[6] );
		
		PDFx1 = ( (p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]+p_set[6][0]) + (p_set[2][3]+p_set[3][3]+p_set[4][3]+p_set[5][3]+p_set[6][3]) )/( 2*p_set[0][0] );
		PDFx2 = ( (p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]+p_set[6][0]) - (p_set[2][3]+p_set[3][3]+p_set[4][3]+p_set[5][3]+p_set[6][3]) )/( 2*p_set[0][0] );
		/// Setting up partons
		p_set[0][0]*= PDFx1;
		p_set[0][1] = 0;
		p_set[0][2] = 0;
		p_set[0][3] = 0;	// to be updated
		p_set[1][0]*= PDFx2;
		p_set[1][1] = 0;
		p_set[1][2] = 0;
		p_set[1][3] = 0;	// to be updated
		if( Debug_Mode ) { printf( "Coefficients for PDF ( x1, x2 ): ( %.10E, %.10E )\n", PDFx1, PDFx2 ); }
	}
	
	
	/// If flag is true, boost to CM frame iff PDF is NOT included.
	if( Boost_To_CM && !Use_PDF_w_pT0 ) 
	{
// 		Boost2CM( ml1, p_set[2], ml2, p_set[3], ml3, p_set[4], ml4, p_set[5] );
		Boost2CM( ml1, p_set[2], ml2, p_set[3], ml3, p_set[4], ml4, p_set[5], 0, p_set[6] );
		CollisionE = p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]+p_set[6][0];
		p_set[0][0] = CollisionE/2;
		p_set[1][0] = CollisionE/2;
	}
	
	
	Mass_4l = sqrt( (p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]+p_set[6][0])*(p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]+p_set[6][0])
		- (p_set[2][1]+p_set[3][1]+p_set[4][1]+p_set[5][1]+p_set[6][1])*(p_set[2][1]+p_set[3][1]+p_set[4][1]+p_set[5][1]+p_set[6][1])
		- (p_set[2][2]+p_set[3][2]+p_set[4][2]+p_set[5][2]+p_set[6][2])*(p_set[2][2]+p_set[3][2]+p_set[4][2]+p_set[5][2]+p_set[6][2])
		- (p_set[2][3]+p_set[3][3]+p_set[4][3]+p_set[5][3]+p_set[6][3])*(p_set[2][3]+p_set[3][3]+p_set[4][3]+p_set[5][3]+p_set[6][3]) );
	
	
	/// Pick quark flavors to use if PDFs are not set. Normalizing coefficients here.
	if( !Use_PDF_w_pT0 )
	{
		buffer = new double;
		(*buffer) = ( ContributionCoeff_d+ContributionCoeff_u+ContributionCoeff_s+ContributionCoeff_c );
		ContributionCoeff_d = (ContributionCoeff_d)/(*buffer);
		ContributionCoeff_u = (ContributionCoeff_u)/(*buffer);
		ContributionCoeff_c = (ContributionCoeff_c)/(*buffer);
		ContributionCoeff_s = (ContributionCoeff_s)/(*buffer);
		delete buffer;
	}
	
	
	if( Debug_Mode )
	{
		printf( "Energy of Parton 1: %.10E\nEnergy of Parton 2: %.10E\n", p_set[0][0], p_set[1][0] );
		printf( "Final-state four-momenta entering ME (E px py px):\n" );
		printf( "%.10E %.10E %.10E %.10E\n", p_set[2][0], p_set[2][1], p_set[2][2], p_set[2][3] );
		printf( "%.10E %.10E %.10E %.10E\n", p_set[3][0], p_set[3][1], p_set[3][2], p_set[3][3] );
		printf( "%.10E %.10E %.10E %.10E\n", p_set[4][0], p_set[4][1], p_set[4][2], p_set[4][3] );
		printf( "%.10E %.10E %.10E %.10E\n", p_set[5][0], p_set[5][1], p_set[5][2], p_set[5][3] );
		printf( "%.10E %.10E %.10E %.10E\n", p_set[6][0], p_set[6][1], p_set[6][2], p_set[6][3] );
		printf( "Sum px=%.10E\n", (p_set[2][1]+p_set[3][1]+p_set[4][1]+p_set[5][1]+p_set[6][1]) );
		printf( "Sum py=%.10E\n", (p_set[2][2]+p_set[3][2]+p_set[4][2]+p_set[5][2]+p_set[6][2]) );
		printf( "Sum pz=%.10E\n", (p_set[2][3]+p_set[3][3]+p_set[4][3]+p_set[5][3]+p_set[6][3]) );
		printf( "Sum E=%.10E\n", (p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]+p_set[6][0]) );
	}
	
	
	/// Background is interesting in any case, except for the Signal Runs or '!' is indicated in the first model to save CPU
	if( Test_Model[0]!='!' && Test_Models.size() == 0 ) { Run_MEKD_MG_ME_Configurator_BKG_ZZ( "qq" ); Background_ME=Signal_ME; }
	else if( Test_Models.size() > 0 )
	{
		if( Test_Models[0][0]!='!' ) { Run_MEKD_MG_ME_Configurator_BKG_ZZ( "qq" ); Background_ME=Signal_ME; }
	}
	
	
	/// Signal ME(s) is(are) chosen here
	if( Test_Models.size() > 0 && Test_Model[0]!='!' )
	{
		Signal_MEs.clear();
		Test_Model_buffer = &(Test_Models[0]);	// Should be NULL or undefined before this point; works as counter=0
	}
	else
		Test_Model_buffer = &Test_Model;	// Should be NULL or undefined before this point
	
	counter = 1;
	while( Test_Model_buffer!=NULL )
	{
		// Is it a parameter card defined?
		if( (*Test_Model_buffer)=="Custom" || (*Test_Model_buffer)=="!Custom" )
			Run_MEKD_MG_ME_Configurator_Custom();
		
		
		// Is it a "background"?
		else if( (*Test_Model_buffer)=="qqZZ" || (*Test_Model_buffer)=="!qqZZ" )
			Run_MEKD_MG_ME_Configurator_BKG_ZZ( "qq" );
		else if( (*Test_Model_buffer)==  "ZZ" || (*Test_Model_buffer)==  "!ZZ" )
			Run_MEKD_MG_ME_Configurator_BKG_ZZ( "qq" );
			
		else if( (*Test_Model_buffer)=="qqDY" || (*Test_Model_buffer)=="!qqDY" )
			Run_MEKD_MG_ME_Configurator_BKG_ZZ( "qq" );
		else if( (*Test_Model_buffer)==  "DY" || (*Test_Model_buffer)==  "!DY" )
			Run_MEKD_MG_ME_Configurator_BKG_ZZ( "NO" );
		
		
		// Is it a Z boson resonance?
		else if( (*Test_Model_buffer)=="qqZ4l_Background" || (*Test_Model_buffer)=="!qqZ4l_Background" )
			Run_MEKD_MG_ME_Configurator_Z4l_BKG( "qq" );
		else if( (*Test_Model_buffer)=="qqZ4l_Signal" || (*Test_Model_buffer)=="!qqZ4l_Signal" )
			Run_MEKD_MG_ME_Configurator_Z4l_SIG( "qq" );
		
		
		/// Resonance to ZZ decay modes (with exceptions for 2l, 2l+A states). Final states: 4 leptons (+photon) also 2 muons (+photon)
		Resonance_decay_mode = "ZZ";
		
		// Is it a spin-0 resonance?
		if(      (*Test_Model_buffer)=="ggSpin0Pm" || (*Test_Model_buffer)=="!ggSpin0Pm" )	// SM Higgs
			Run_MEKD_MG_ME_Configurator_Spin0Pm( "gg" );
		else if( (*Test_Model_buffer)==  "Spin0Pm" || (*Test_Model_buffer)==  "!Spin0Pm" )
			Run_MEKD_MG_ME_Configurator_Spin0Pm( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin0M" || (*Test_Model_buffer)=="!ggSpin0M" )
			Run_MEKD_MG_ME_Configurator_Spin0M( "gg" );
		else if( (*Test_Model_buffer)==  "Spin0M" || (*Test_Model_buffer)==  "!Spin0M" )
			Run_MEKD_MG_ME_Configurator_Spin0M( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin0Ph" || (*Test_Model_buffer)=="!ggSpin0Ph" )
			Run_MEKD_MG_ME_Configurator_Spin0Ph( "gg" );
		else if( (*Test_Model_buffer)==  "Spin0Ph" || (*Test_Model_buffer)==  "!Spin0Ph" )
			Run_MEKD_MG_ME_Configurator_Spin0Ph( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin0" || (*Test_Model_buffer)=="!ggSpin0" )
			Run_MEKD_MG_ME_Configurator_Spin0( "gg" );
		else if( (*Test_Model_buffer)==  "Spin0" || (*Test_Model_buffer)==  "!Spin0" )
			Run_MEKD_MG_ME_Configurator_Spin0( "NO" );
		
		
		// Is it a spin-1 resonance?
		else if( (*Test_Model_buffer)=="qqSpin1M" || (*Test_Model_buffer)=="!qqSpin1M" )
			Run_MEKD_MG_ME_Configurator_Spin1M( "qq" );
		else if( (*Test_Model_buffer)==  "Spin1M" || (*Test_Model_buffer)==  "!Spin1M" )
			Run_MEKD_MG_ME_Configurator_Spin1M( "NO" );
		
		else if( (*Test_Model_buffer)=="qqSpin1P" || (*Test_Model_buffer)=="!qqSpin1P" )
			Run_MEKD_MG_ME_Configurator_Spin1P( "qq" );
		else if( (*Test_Model_buffer)==  "Spin1P" || (*Test_Model_buffer)==  "!Spin1P" )
			Run_MEKD_MG_ME_Configurator_Spin1P( "NO" );
		
		else if( (*Test_Model_buffer)=="qqSpin1" || (*Test_Model_buffer)=="!qqSpin1" )
			Run_MEKD_MG_ME_Configurator_Spin1( "qq" );
		else if( (*Test_Model_buffer)==  "Spin1" || (*Test_Model_buffer)==  "!Spin1" )
			Run_MEKD_MG_ME_Configurator_Spin1( "NO" );
		
		
		// Is it a spin-2 resonance?
		else if( (*Test_Model_buffer)=="ggSpin2Pm" || (*Test_Model_buffer)=="!ggSpin2Pm" )
			Run_MEKD_MG_ME_Configurator_Spin2Pm( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2Pm" || (*Test_Model_buffer)=="!qqSpin2Pm" )
			Run_MEKD_MG_ME_Configurator_Spin2Pm( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2Pm" || (*Test_Model_buffer)==  "!Spin2Pm" )
			Run_MEKD_MG_ME_Configurator_Spin2Pm( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin2Ph" || (*Test_Model_buffer)=="!ggSpin2Ph" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2Ph" || (*Test_Model_buffer)=="!qqSpin2Ph" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2Ph" || (*Test_Model_buffer)==  "!Spin2Ph" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin2Mh" || (*Test_Model_buffer)=="!ggSpin2Mh" )
			Run_MEKD_MG_ME_Configurator_Spin2Mh( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2Mh" || (*Test_Model_buffer)=="!qqSpin2Mh" )
			Run_MEKD_MG_ME_Configurator_Spin2Mh( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2Mh" || (*Test_Model_buffer)==  "!Spin2Mh" )
			Run_MEKD_MG_ME_Configurator_Spin2Mh( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin2Pb" || (*Test_Model_buffer)=="!ggSpin2Pb" )
			Run_MEKD_MG_ME_Configurator_Spin2Pb( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2Pb" || (*Test_Model_buffer)=="!qqSpin2Pb" )
			Run_MEKD_MG_ME_Configurator_Spin2Pb( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2Pb" || (*Test_Model_buffer)==  "!Spin2Pb" )
			Run_MEKD_MG_ME_Configurator_Spin2Pb( "NO" );
		
		
		else if( (*Test_Model_buffer)=="ggSpin2Ph2" || (*Test_Model_buffer)=="!ggSpin2Ph2" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph2( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2Ph2" || (*Test_Model_buffer)=="!qqSpin2Ph2" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph2( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2Ph2" || (*Test_Model_buffer)==  "!Spin2Ph2" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph2( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin2Ph3" || (*Test_Model_buffer)=="!ggSpin2Ph3" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph3( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2Ph3" || (*Test_Model_buffer)=="!qqSpin2Ph3" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph3( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2Ph3" || (*Test_Model_buffer)==  "!Spin2Ph3" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph3( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin2Ph6" || (*Test_Model_buffer)=="!ggSpin2Ph6" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph6( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2Ph6" || (*Test_Model_buffer)=="!qqSpin2Ph6" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph6( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2Ph6" || (*Test_Model_buffer)==  "!Spin2Ph6" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph6( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin2Ph7" || (*Test_Model_buffer)=="!ggSpin2Ph7" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph7( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2Ph7" || (*Test_Model_buffer)=="!qqSpin2Ph7" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph7( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2Ph7" || (*Test_Model_buffer)==  "!Spin2Ph7" )
			Run_MEKD_MG_ME_Configurator_Spin2Ph7( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin2Mh9" || (*Test_Model_buffer)=="!ggSpin2Mh9" )
			Run_MEKD_MG_ME_Configurator_Spin2Mh9( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2Mh9" || (*Test_Model_buffer)=="!qqSpin2Mh9" )
			Run_MEKD_MG_ME_Configurator_Spin2Mh9( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2Mh9" || (*Test_Model_buffer)==  "!Spin2Mh9" )
			Run_MEKD_MG_ME_Configurator_Spin2Mh9( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin2Mh10" || (*Test_Model_buffer)=="!ggSpin2Mh10" )
			Run_MEKD_MG_ME_Configurator_Spin2Mh10( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2Mh10" || (*Test_Model_buffer)=="!qqSpin2Mh10" )
			Run_MEKD_MG_ME_Configurator_Spin2Mh10( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2Mh10" || (*Test_Model_buffer)==  "!Spin2Mh10" )
			Run_MEKD_MG_ME_Configurator_Spin2Mh10( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin2" || (*Test_Model_buffer)=="!ggSpin2" )
			Run_MEKD_MG_ME_Configurator_Spin2( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2" || (*Test_Model_buffer)=="!qqSpin2" )
			Run_MEKD_MG_ME_Configurator_Spin2( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2" || (*Test_Model_buffer)==  "!Spin2" )
			Run_MEKD_MG_ME_Configurator_Spin2( "NO" );
		
		
		
		/// Resonance to 2l decay modes. Final states: 4 leptons (+photon)
		Resonance_decay_mode = "2l";
		
		// Is it a spin-0 resonance?
		if(      (*Test_Model_buffer)=="ggSpin0Pm_2f" || (*Test_Model_buffer)=="!ggSpin0Pm_2f" )	// SM Higgs
			Run_MEKD_MG_ME_Configurator_Spin0Pm( "gg" );
		else if( (*Test_Model_buffer)==  "Spin0Pm_2f" || (*Test_Model_buffer)==  "!Spin0Pm_2f" )
			Run_MEKD_MG_ME_Configurator_Spin0Pm( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin0M_2f" || (*Test_Model_buffer)=="!ggSpin0M_2f" )
			Run_MEKD_MG_ME_Configurator_Spin0M( "gg" );
		else if( (*Test_Model_buffer)==  "Spin0M_2f" || (*Test_Model_buffer)==  "!Spin0M_2f" )
			Run_MEKD_MG_ME_Configurator_Spin0M( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin0_2f" || (*Test_Model_buffer)=="!ggSpin0_2f" )
			Run_MEKD_MG_ME_Configurator_Spin0( "gg" );
		else if( (*Test_Model_buffer)==  "Spin0_2f" || (*Test_Model_buffer)==  "!Spin0_2f" )
			Run_MEKD_MG_ME_Configurator_Spin0( "NO" );
		
		
		// Is it a spin-1 resonance?
		else if( (*Test_Model_buffer)=="qqSpin1M_2f" || (*Test_Model_buffer)=="!qqSpin1M_2f" )
			Run_MEKD_MG_ME_Configurator_Spin1M( "qq" );
		else if( (*Test_Model_buffer)==  "Spin1M" || (*Test_Model_buffer)==  "!Spin1M_2f" )
			Run_MEKD_MG_ME_Configurator_Spin1M( "NO" );
		
		else if( (*Test_Model_buffer)=="qqSpin1P_2f" || (*Test_Model_buffer)=="!qqSpin1P_2f" )
			Run_MEKD_MG_ME_Configurator_Spin1P( "qq" );
		else if( (*Test_Model_buffer)==  "Spin1P_2f" || (*Test_Model_buffer)==  "!Spin1P_2f" )
			Run_MEKD_MG_ME_Configurator_Spin1P( "NO" );
		
		else if( (*Test_Model_buffer)=="qqSpin1_2f" || (*Test_Model_buffer)=="!qqSpin1_2f" )
			Run_MEKD_MG_ME_Configurator_Spin1( "qq" );
		else if( (*Test_Model_buffer)==  "Spin1_2f" || (*Test_Model_buffer)==  "!Spin1_2f" )
			Run_MEKD_MG_ME_Configurator_Spin1( "NO" );
		
		
		// Is it a spin-2 resonance?
		else if( (*Test_Model_buffer)=="ggSpin2Pm_2f" || (*Test_Model_buffer)=="!ggSpin2Pm_2f" )
			Run_MEKD_MG_ME_Configurator_Spin2Pm( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2Pm_2f" || (*Test_Model_buffer)=="!qqSpin2Pm_2f" )
			Run_MEKD_MG_ME_Configurator_Spin2Pm( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2Pm_2f" || (*Test_Model_buffer)==  "!Spin2Pm_2f" )
			Run_MEKD_MG_ME_Configurator_Spin2Pm( "NO" );
		
		else if( (*Test_Model_buffer)=="ggSpin2_2f" || (*Test_Model_buffer)=="!ggSpin2_2f" )
			Run_MEKD_MG_ME_Configurator_Spin2Pm( "gg" );
		else if( (*Test_Model_buffer)=="qqSpin2_2f" || (*Test_Model_buffer)=="!qqSpin2_2f" )
			Run_MEKD_MG_ME_Configurator_Spin2Pm( "qq" );
		else if( (*Test_Model_buffer)==  "Spin2_2f" || (*Test_Model_buffer)==  "!Spin2_2f" )
			Run_MEKD_MG_ME_Configurator_Spin2Pm( "NO" );
		
		
		// Is it a RAW MG5_aMC ME?
		if(      (*Test_Model_buffer)=="ggCPPProcess" || (*Test_Model_buffer)=="!ggCPPProcess" )	// ME_RAW
			Run_MEKD_MG_ME_Configurator_CPPProcess( "gg" );
		else if( (*Test_Model_buffer)=="qqCPPProcess" || (*Test_Model_buffer)=="!qqCPPProcess" )
			Run_MEKD_MG_ME_Configurator_CPPProcess( "qq" );
		else if( (*Test_Model_buffer)==  "CPPProcess" || (*Test_Model_buffer)==  "!CPPProcess" )
			Run_MEKD_MG_ME_Configurator_CPPProcess( "NO" );
		
		
		if( Debug_Mode ) cout << "Evaluated model: " << (*Test_Model_buffer) << "; calculated ME: " << Signal_ME << endl;
		
		if( Test_Models.size() > 0 && Test_Model[0]!='!' )
		{
			Signal_MEs.push_back( Signal_ME );
			
			if( counter < Test_Models.size() )
				Test_Model_buffer = &(Test_Models[counter]);
			else
				Test_Model_buffer = NULL;
		}
		else
			Test_Model_buffer = NULL;
		
		counter++;
	}
	
	
	if( Debug_Mode )
	{
		printf( "Status after ME calculations are over:\n" );
		printf( "Energy of Parton 1: %.10E\nEnergy of Parton 2: %.10E\n", p_set[0][0], p_set[1][0] );
		printf( "Final-state four-momenta exiting ME(s) (E px py px):\n" );
		printf( "%.10E %.10E %.10E %.10E\n", p_set[2][0], p_set[2][1], p_set[2][2], p_set[2][3] );
		printf( "%.10E %.10E %.10E %.10E\n", p_set[3][0], p_set[3][1], p_set[3][2], p_set[3][3] );
		printf( "%.10E %.10E %.10E %.10E\n", p_set[4][0], p_set[4][1], p_set[4][2], p_set[4][3] );
		printf( "%.10E %.10E %.10E %.10E\n", p_set[5][0], p_set[5][1], p_set[5][2], p_set[5][3] );
		printf( "%.10E %.10E %.10E %.10E\n", p_set[6][0], p_set[6][1], p_set[6][2], p_set[6][3] );
		printf( "Sum px=%.10E\n", (p_set[2][1]+p_set[3][1]+p_set[4][1]+p_set[5][1]+p_set[6][1]) );
		printf( "Sum py=%.10E\n", (p_set[2][2]+p_set[3][2]+p_set[4][2]+p_set[5][2]+p_set[6][2]) );
		printf( "Sum pz=%.10E\n", (p_set[2][3]+p_set[3][3]+p_set[4][3]+p_set[5][3]+p_set[6][3]) );
		printf( "Sum E=%.10E\n", (p_set[2][0]+p_set[3][0]+p_set[4][0]+p_set[5][0]+p_set[6][0]) );
	}
	
	
	if( Test_Model[0]!='!' ) KD = log( Signal_ME/Background_ME );
	
	return 0;
}



int MEKD_MG::Run_MEKD_MG(string Input_Model)
{
	buffer_string = Test_Model;
	Test_Model = "!";
	Test_Model += Input_Model;
	
	error_value = Run_MEKD_MG();
	
	Test_Model = buffer_string;
	return error_value;
}



///#include "MEKD_MG_RunMEs.cpp"
///////////////////////////////////
/// INCLUDED MEKD_MG_RunMEs.cpp ///
/// code follows below          ///
///                             ///
///  Part responsible for ME    ///
///  calculations               ///
///////////////////////////////////




/// ZZ and DY processes
int MEKD_MG::Run_MEKD_MG_ME_Configurator_BKG_ZZ(string initial_state)
{
	return Run_MEKD_MG_ME_Dispatcher_BKG_ZZ( initial_state );
}



/// Z4l Background processes (t channel)
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Z4l_BKG(string initial_state)
{
	if( Use_mZ4l_eq_m4l ) Set_Of_Model_Parameters.set_block_entry( "mass", 23, Mass_4l );
	return Run_MEKD_MG_ME_Dispatcher_Z4l_BKG( initial_state );
}



/// Z4l Signal processes (s channel)
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Z4l_SIG(string initial_state)
{
	if( Use_mZ4l_eq_m4l ) Set_Of_Model_Parameters.set_block_entry( "mass", 23, Mass_4l );
	return Run_MEKD_MG_ME_Dispatcher_Z4l_SIG( initial_state );
}



int MEKD_MG::Run_MEKD_MG_ME_Configurator_Custom()
{
	if( (error_value=Run_MEKD_MG_ME_Dispatcher_SIG_Spin0( "gg" ))!=0 ) return error_value;
	buffer_Custom = Signal_ME;
	if( (error_value=Run_MEKD_MG_ME_Dispatcher_SIG_Spin1( "qq" ))!=0 ) return error_value;
	buffer_Custom += Signal_ME;
	if( (error_value=Run_MEKD_MG_ME_Dispatcher_SIG_Spin2( "gg" ))!=0 ) return error_value;
	buffer_Custom += Signal_ME;
	if( (error_value=Run_MEKD_MG_ME_Dispatcher_SIG_Spin2( "qq" ))!=0 ) return error_value;
	Signal_ME += buffer_Custom;
	
	return 0;
}



int MEKD_MG::Run_MEKD_MG_ME_Configurator_CPPProcess(string initial_state)
{
	return Run_MEKD_MG_ME_Dispatcher_CPPProcess( initial_state );
}



////////////////////////////////////
/// Spin-0 resonances start here ///
////////////////////////////////////



/// A SM Higgs
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin0Pm(string initial_state)
{
	Predefined_Model = true;
	Mixing_Coefficients_Spin0_internal[0] = complex<double>(1, 0);	// fits 2l case
	Mixing_Coefficients_Spin0_internal[1] = complex<double>(0, 0);	// fits 2l case
	Mixing_Coefficients_Spin0_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin0_internal[3] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin0_internal[0] = complex<double>(1, 0);
		Mixing_Coefficients_Spin0_internal[1] = complex<double>(0, 0);
	}
	
	return Run_MEKD_MG_ME_Configurator_Spin0( initial_state );
}



/// A pseudoscalar
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin0M(string initial_state)
{
	Predefined_Model = true;
	Mixing_Coefficients_Spin0_internal[0] = complex<double>(0, 0);
	Mixing_Coefficients_Spin0_internal[1] = complex<double>(0, 0);
	Mixing_Coefficients_Spin0_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin0_internal[3] = complex<double>(1, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin0_internal[0] = complex<double>(0, 0);
		Mixing_Coefficients_Spin0_internal[1] = complex<double>(1, 0);
	}
	
	return Run_MEKD_MG_ME_Configurator_Spin0( initial_state );
}



/// A scalar with higher-order couplings
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin0Ph(string initial_state)
{
	Predefined_Model = true;
	Mixing_Coefficients_Spin0_internal[0] = complex<double>(0, 0);
	Mixing_Coefficients_Spin0_internal[1] = complex<double>(1, 0);
	Mixing_Coefficients_Spin0_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin0_internal[3] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin0_internal[0] = complex<double>(0, 0);
		Mixing_Coefficients_Spin0_internal[1] = complex<double>(0, 0);
	}
	
	return Run_MEKD_MG_ME_Configurator_Spin0( initial_state );
}



////////////////////////////////////
/// Spin-1 resonances start here ///
////////////////////////////////////



/// A vector default configuration
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin1M(string initial_state)
{
	Predefined_Model = true;
	Mixing_Coefficients_Spin1_internal[0] = complex<double>(1, 0);	// Production, same for 2l
	Mixing_Coefficients_Spin1_internal[1] = complex<double>(0, 0);
	Mixing_Coefficients_Spin1_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin1_internal[3] = complex<double>(0, 0);
	Mixing_Coefficients_Spin1_internal[4] = complex<double>(1, 0);	// Decay
	Mixing_Coefficients_Spin1_internal[5] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin1_internal[6] = complex<double>(0, 0);
		Mixing_Coefficients_Spin1_internal[7] = complex<double>(0, 0);
	}
	
	return Run_MEKD_MG_ME_Configurator_Spin1( initial_state );
}



/// A vector default configuration
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin1P(string initial_state)
{
	Predefined_Model = true;
	Mixing_Coefficients_Spin1_internal[0] = complex<double>(0, 0);	// Production
	Mixing_Coefficients_Spin1_internal[1] = complex<double>(1, 0);	// same for 2l
	Mixing_Coefficients_Spin1_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin1_internal[3] = complex<double>(0, 0);
	Mixing_Coefficients_Spin1_internal[4] = complex<double>(0, 0);	// Decay
	Mixing_Coefficients_Spin1_internal[5] = complex<double>(1, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin1_internal[6] = complex<double>(0, 0);
		Mixing_Coefficients_Spin1_internal[7] = complex<double>(0, 0);
	}
	
	return Run_MEKD_MG_ME_Configurator_Spin1( initial_state );
}



////////////////////////////////////
/// Spin-2 resonances start here ///
////////////////////////////////////



/// A minimal-coupling KK graviton
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin2Pm(string initial_state)
{
	Predefined_Model = true;
	if( initial_state=="gg" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);	// Production
	Mixing_Coefficients_Spin2_internal[1] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[3] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[4] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[5] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[6] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[7] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[8] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[9] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[10] = complex<double>(-1, 0);	// Decay
	Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[14] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[15] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[16] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[17] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[18] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[19] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin2_internal[10] = complex<double>(1, 0);
		Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	}
	
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[1] = complex<double>(1, 0);
	
	return Run_MEKD_MG_ME_Configurator_Spin2( initial_state );
}



///
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin2Ph(string initial_state)
{
	Predefined_Model = true;
	Mixing_Coefficients_Spin2_internal[0] = complex<double>(0, 0);	// Production
	Mixing_Coefficients_Spin2_internal[1] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[2] = complex<double>(0, 0);
	if( initial_state=="gg" ) Mixing_Coefficients_Spin2_internal[3] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[4] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[5] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[6] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[7] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[8] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[9] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);	// Decay
	Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[13] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[14] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[15] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[16] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[17] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[18] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[19] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	}
	
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[1] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[3] = complex<double>(0, 0);
	
	return Run_MEKD_MG_ME_Configurator_Spin2( initial_state );
}



///
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin2Mh(string initial_state)
{
	Predefined_Model = true;
	Mixing_Coefficients_Spin2_internal[0] = complex<double>(0, 0);	// Production
	Mixing_Coefficients_Spin2_internal[1] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[3] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[4] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[5] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[6] = complex<double>(0, 0);
	if( initial_state=="gg" ) Mixing_Coefficients_Spin2_internal[7] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[8] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[9] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);	// Decay
	Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[14] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[15] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[16] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[17] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[18] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[19] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	}
	
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[1] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[7] = complex<double>(0, 0);
	
	return Run_MEKD_MG_ME_Configurator_Spin2( initial_state );
}



///
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin2Pb(string initial_state)
{
	Predefined_Model = true;
	if( initial_state=="gg" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);	// Production
	Mixing_Coefficients_Spin2_internal[1] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[3] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[4] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[5] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[6] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[7] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[8] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[9] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);	// Decay
	Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[14] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[15] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[16] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[17] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[18] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[19] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	}
	
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[1] = complex<double>(1, 0);
	
	return Run_MEKD_MG_ME_Configurator_Spin2( initial_state );
}



///
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin2Ph2(string initial_state)
{
	Predefined_Model = true;
	Mixing_Coefficients_Spin2_internal[0] = complex<double>(0, 0);	// Production
	if( initial_state=="gg" ) Mixing_Coefficients_Spin2_internal[1] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[3] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[4] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[5] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[6] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[7] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[8] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[9] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);	// Decay
	Mixing_Coefficients_Spin2_internal[11] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[14] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[15] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[16] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[17] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[18] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[19] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	}
	
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[1] = complex<double>(1, 0);
	
	return Run_MEKD_MG_ME_Configurator_Spin2( initial_state );
}



///
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin2Ph3(string initial_state)
{
	Predefined_Model = true;
	Mixing_Coefficients_Spin2_internal[0] = complex<double>(0, 0);	// Production
	Mixing_Coefficients_Spin2_internal[1] = complex<double>(0, 0);
	if( initial_state=="gg" ) Mixing_Coefficients_Spin2_internal[2] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[3] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[4] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[5] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[6] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[7] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[8] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[9] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);	// Decay
	Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[12] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[14] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[15] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[16] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[17] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[18] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[19] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	}
	
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[1] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[2] = complex<double>(0, 0);
	
	return Run_MEKD_MG_ME_Configurator_Spin2( initial_state );
}



///
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin2Ph6(string initial_state)
{
	Predefined_Model = true;
	if( initial_state=="gg" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);	// Production
	Mixing_Coefficients_Spin2_internal[1] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[3] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[4] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[5] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[6] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[7] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[8] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[9] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);	// Decay
	Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[14] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[15] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[16] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[17] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[18] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[19] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	}
	
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[1] = complex<double>(1, 0);
	
	return Run_MEKD_MG_ME_Configurator_Spin2( initial_state );
}



///
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin2Ph7(string initial_state)
{
	Predefined_Model = true;
	if( initial_state=="gg" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);	// Production
	Mixing_Coefficients_Spin2_internal[1] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[3] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[4] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[5] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[6] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[7] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[8] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[9] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);	// Decay
	Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[14] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[15] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[16] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[17] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[18] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[19] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	}
	
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[1] = complex<double>(1, 0);
	
	return Run_MEKD_MG_ME_Configurator_Spin2( initial_state );
}



///
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin2Mh9(string initial_state)
{
	Predefined_Model = true;
	Mixing_Coefficients_Spin2_internal[0] = complex<double>(0, 0);	// Production
	Mixing_Coefficients_Spin2_internal[1] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[3] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[4] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[5] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[6] = complex<double>(0, 0);
	if( initial_state=="gg" ) Mixing_Coefficients_Spin2_internal[7] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[8] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[9] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);	// Decay
	Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[14] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[15] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[16] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[17] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[18] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[19] = complex<double>(0, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	}
	
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[1] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[7] = complex<double>(0, 0);
	
	return Run_MEKD_MG_ME_Configurator_Spin2( initial_state );
}



///
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin2Mh10(string initial_state)
{
	Predefined_Model = true;
	Mixing_Coefficients_Spin2_internal[0] = complex<double>(0, 0);	// Production
	Mixing_Coefficients_Spin2_internal[1] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[2] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[3] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[4] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[5] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[6] = complex<double>(0, 0);
	if( initial_state=="gg" ) Mixing_Coefficients_Spin2_internal[7] = complex<double>(1, 0);
	Mixing_Coefficients_Spin2_internal[8] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[9] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);	// Decay
	Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[14] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[15] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[16] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[17] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[18] = complex<double>(0, 0);
	Mixing_Coefficients_Spin2_internal[19] = complex<double>(1, 0);
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" || Resonance_decay_mode=="2l" )
	{
		Mixing_Coefficients_Spin2_internal[10] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[11] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[12] = complex<double>(0, 0);
		Mixing_Coefficients_Spin2_internal[13] = complex<double>(0, 0);
	}
	
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[0] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[1] = complex<double>(1, 0);
	if( initial_state=="qq" ) Mixing_Coefficients_Spin2_internal[7] = complex<double>(0, 0);
	
	return Run_MEKD_MG_ME_Configurator_Spin2( initial_state );
}



/// A generic spin-0 resonance handler
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin0(string initial_state)
{
	if( Predefined_Model )
	{
		buffer_complex = Mixing_Coefficients_Spin0_internal;
		Predefined_Model = false;
	}
	else buffer_complex = Mixing_Coefficients_Spin0;
	
	if( Use_mh_eq_m4l )
	{	
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000006, Mass_4l );
		
		if( Use_Higgs_width )
		{
			if( Vary_resonance_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, static_cast<double>( MEKD_CalcHEP_Extra::Higgs_width(Mass_4l) ) );
			else Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, Higgs_width );
		}
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, 1 );
		
		
		LmbdGG_calculated = LmbdGG(Mass_4l);
	}
	else
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000006, Higgs_mass );
		
		if( Use_Higgs_width )
		{
			if( Vary_resonance_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, static_cast<double>( MEKD_CalcHEP_Extra::Higgs_width(Mass_4l) ) );
			else Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, Higgs_width );
		}
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000006, 1 );
		
		
		LmbdGG_calculated = LmbdGG(Higgs_mass);
	}
	
	if( Vary_signal_couplings )
	{
		//gg
		if( Fix_Spin0_Production )
		{
			Set_Of_Model_Parameters.set_block_entry( "heff", 1, complex<double>(0, 0) );
			Set_Of_Model_Parameters.set_block_entry( "heff", 2, complex<double>(4*LmbdGG_calculated, 0) );	// Spin0Pm
			Set_Of_Model_Parameters.set_block_entry( "heff", 3, complex<double>(0, 0) );
			Set_Of_Model_Parameters.set_block_entry( "heff", 4, complex<double>(0, 0) );	// Spin0M
		}
		else
		{
			Set_Of_Model_Parameters.set_block_entry( "heff", 1, complex<double>(0, 0) );
			Set_Of_Model_Parameters.set_block_entry( "heff", 2, (buffer_complex[0]+buffer_complex[1]+buffer_complex[2])*complex<double>(4*LmbdGG_calculated, 0) );	// Spin0P
			Set_Of_Model_Parameters.set_block_entry( "heff", 3, complex<double>(0, 0) );
			Set_Of_Model_Parameters.set_block_entry( "heff", 4, buffer_complex[3]*complex<double>(4*LmbdGG_calculated, 0) );	// Spin0M
		}
		
		// Decay to ZZ
		if( Use_mh_eq_m4l )
		{
			Set_Of_Model_Parameters.set_block_entry( "heff", 5, buffer_complex[0]*complex<double>(hZZ_coupling, 0) );	// Spin0Pm
			Set_Of_Model_Parameters.set_block_entry( "heff", 6, buffer_complex[1]*complex<double>(hZZ_coupling/params_m_Z/params_m_Z, 0) );	// Spin0Ph
			Set_Of_Model_Parameters.set_block_entry( "heff", 7, buffer_complex[2]*complex<double>(hZZ_coupling/params_m_Z/params_m_Z/Mass_4l/Mass_4l, 0) );	// Spin0Ph+
			Set_Of_Model_Parameters.set_block_entry( "heff", 8, buffer_complex[3]*complex<double>(hZZ_coupling/params_m_Z/params_m_Z, 0) );	// Spin0M
		}
		else
		{
			Set_Of_Model_Parameters.set_block_entry( "heff", 5, buffer_complex[0]*complex<double>(hZZ_coupling, 0) );	// Spin0Pm
			Set_Of_Model_Parameters.set_block_entry( "heff", 6, buffer_complex[1]*complex<double>(hZZ_coupling/params_m_Z/params_m_Z, 0) );	// Spin0Ph
			Set_Of_Model_Parameters.set_block_entry( "heff", 7, buffer_complex[2]*complex<double>(hZZ_coupling/params_m_Z/params_m_Z/Higgs_mass/Higgs_mass, 0) );	// Spin0Ph+
			Set_Of_Model_Parameters.set_block_entry( "heff", 8, buffer_complex[3]*complex<double>(hZZ_coupling/params_m_Z/params_m_Z, 0) );	// Spin0M
		}
		
		// Decay to 2l (or 2mu in 2f)
		Set_Of_Model_Parameters.set_block_entry( "heff", 19, buffer_complex[0]*complex<double>(4.291210e-04, 0) );
		Set_Of_Model_Parameters.set_block_entry( "heff", 20, buffer_complex[1]*complex<double>(4.291210e-04, 0) );
		
		// Decay to 2e (2f)
		Set_Of_Model_Parameters.set_block_entry( "heff", 21, buffer_complex[0]*complex<double>(4.291210e-04, 0) );	// for Hee should be 2.075371e-06
		Set_Of_Model_Parameters.set_block_entry( "heff", 22, buffer_complex[1]*complex<double>(4.291210e-04, 0) );
	}
	
	//qq
	params_rhou01 = complex<double>(0, 0);
	params_rhou02 = complex<double>(0, 0);
	params_rhoc01 = complex<double>(0, 0);
	params_rhoc02 = complex<double>(0, 0);
	params_rhod01 = complex<double>(0, 0);
	params_rhod02 = complex<double>(0, 0);
	params_rhos01 = complex<double>(0, 0);
	params_rhos02 = complex<double>(0, 0);
	params_rhob01 = complex<double>(0, 0);
	params_rhob02 = complex<double>(0, 0);
	
	buffer_complex = NULL;
	
	return Run_MEKD_MG_ME_Dispatcher_SIG_Spin0( initial_state );
}



/// A generic spin-1 resonance handler
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin1(string initial_state)
{
	if( Predefined_Model )
	{
		buffer_complex = Mixing_Coefficients_Spin1_internal;
		Predefined_Model = false;
	}
	else buffer_complex = Mixing_Coefficients_Spin1;
	
	if( Use_mh_eq_m4l )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 300, Mass_4l );
		
		if( Use_Higgs_width )
		{
			if( Vary_resonance_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 300, static_cast<double>( MEKD_CalcHEP_Extra::Higgs_width(Mass_4l) ) );
			else Set_Of_Model_Parameters.set_block_entry( "decay", 300, Higgs_width );
		}
		else Set_Of_Model_Parameters.set_block_entry( "decay", 300, 1 );
		
		
		LmbdGG_calculated = LmbdGG(Mass_4l);
	}
	else
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 300, Higgs_mass );
		
		if( Use_Higgs_width )
		{
			if( Vary_resonance_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 300, static_cast<double>( MEKD_CalcHEP_Extra::Higgs_width(Mass_4l) ) );
			else Set_Of_Model_Parameters.set_block_entry( "decay", 300, Higgs_width );
		}
		else Set_Of_Model_Parameters.set_block_entry( "decay", 300, 1 );
		
		
		LmbdGG_calculated = LmbdGG(Higgs_mass);
	}
	
	if( Vary_signal_couplings )
	{
		//qq
		if( Fix_Spin1_Production )
		{
			params_rhod11 = complex<double>(sqrt(1/2)*LmbdGG_calculated*v_expectation, 0);
			params_rhos11 = complex<double>(sqrt(1/2)*LmbdGG_calculated*v_expectation, 0);
			params_rhob11 = complex<double>(sqrt(1/2)*LmbdGG_calculated*v_expectation, 0);
			params_rhou11 = complex<double>(sqrt(1/2)*LmbdGG_calculated*v_expectation, 0);
			params_rhoc11 = complex<double>(sqrt(1/2)*LmbdGG_calculated*v_expectation, 0);
			params_rhod12 = complex<double>(0, 0);
			params_rhos12 = complex<double>(0, 0);
			params_rhob12 = complex<double>(0, 0);
			params_rhou12 = complex<double>(0, 0);
			params_rhoc12 = complex<double>(0, 0);
			params_rhod13 = complex<double>(0, 0);
			params_rhos13 = complex<double>(0, 0);
			params_rhob13 = complex<double>(0, 0);
			params_rhou13 = complex<double>(0, 0);
			params_rhoc13 = complex<double>(0, 0);
			params_rhod14 = complex<double>(0, 0);
			params_rhos14 = complex<double>(0, 0);
			params_rhob14 = complex<double>(0, 0);
			params_rhou14 = complex<double>(0, 0);
			params_rhoc14 = complex<double>(0, 0);
		}
		else
		{
			params_rhod11 = buffer_complex[0]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhos11 = buffer_complex[0]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhob11 = buffer_complex[0]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhou11 = buffer_complex[0]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhoc11 = buffer_complex[0]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhod12 = buffer_complex[1]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhos12 = buffer_complex[1]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhob12 = buffer_complex[1]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhou12 = buffer_complex[1]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhoc12 = buffer_complex[1]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhod13 = buffer_complex[2]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhos13 = buffer_complex[2]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhob13 = buffer_complex[2]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhou13 = buffer_complex[2]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhoc13 = buffer_complex[2]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhod14 = buffer_complex[3]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhos14 = buffer_complex[3]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhob14 = buffer_complex[3]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhou14 = buffer_complex[3]*complex<double>(LmbdGG_calculated*v_expectation, 0);
			params_rhoc14 = buffer_complex[3]*complex<double>(LmbdGG_calculated*v_expectation, 0);
		}
		
		// Decay to ZZ
		Set_Of_Model_Parameters.set_block_entry( "vec", 1, buffer_complex[4]*complex<double>(hZZ_coupling/2/params_m_Z, 0) );
		Set_Of_Model_Parameters.set_block_entry( "vec", 2, buffer_complex[5]*complex<double>(hZZ_coupling/4/params_m_Z, 0) );
		
		// Decay to 2l (or 2mu in 2f)
		Set_Of_Model_Parameters.set_block_entry( "vec", 23, buffer_complex[4]*complex<double>(4.291210e-04, 0) );
		Set_Of_Model_Parameters.set_block_entry( "vec", 24, buffer_complex[5]*complex<double>(4.291210e-04, 0) );
		Set_Of_Model_Parameters.set_block_entry( "vec", 25, buffer_complex[6]*complex<double>(4.291210e-04, 0) );
		Set_Of_Model_Parameters.set_block_entry( "vec", 26, buffer_complex[7]*complex<double>(4.291210e-04, 0) );
		
		// Decay to 2e (2f)
		Set_Of_Model_Parameters.set_block_entry( "vec", 27, buffer_complex[4]*complex<double>(4.291210e-04, 0) );	// for Hee should be 2.075371e-06
		Set_Of_Model_Parameters.set_block_entry( "vec", 28, buffer_complex[5]*complex<double>(4.291210e-04, 0) );
		Set_Of_Model_Parameters.set_block_entry( "vec", 29, buffer_complex[6]*complex<double>(4.291210e-04, 0) );
		Set_Of_Model_Parameters.set_block_entry( "vec", 30, buffer_complex[7]*complex<double>(4.291210e-04, 0) );
	}
	
	buffer_complex = NULL;
	
	return Run_MEKD_MG_ME_Dispatcher_SIG_Spin1( initial_state );
}



/// A generic spin-2 resonance handler
int MEKD_MG::Run_MEKD_MG_ME_Configurator_Spin2(string initial_state)
{
	if( Predefined_Model )
	{
		buffer_complex = Mixing_Coefficients_Spin2_internal;
		Predefined_Model = false;
	}
	else buffer_complex = Mixing_Coefficients_Spin2;
	
	if( Use_mh_eq_m4l )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000007, Mass_4l );
		
		if( Use_Higgs_width )
		{
			if( Vary_resonance_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000007, static_cast<double>( MEKD_CalcHEP_Extra::Higgs_width(Mass_4l) ) );
			else Set_Of_Model_Parameters.set_block_entry( "decay", 9000007, Higgs_width );
		}
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000007, 1 );
		
		
		LmbdGG_calculated = LmbdGG(Mass_4l);
	}
	else
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 9000007, Higgs_mass );
		
		if( Use_Higgs_width )
		{
			if( Vary_resonance_width ) Set_Of_Model_Parameters.set_block_entry( "decay", 9000007, static_cast<double>( MEKD_CalcHEP_Extra::Higgs_width(Mass_4l) ) );
			else Set_Of_Model_Parameters.set_block_entry( "decay", 9000007, Higgs_width );
		}
		else Set_Of_Model_Parameters.set_block_entry( "decay", 9000007, 1 );
		
		
		LmbdGG_calculated = LmbdGG(Higgs_mass);
	}
	
	
	if( Vary_signal_couplings )
	{
		//gg
		if( Use_mh_eq_m4l )
		{
			Set_Of_Model_Parameters.set_block_entry( "gravity", 1, buffer_complex[0]*complex<double>(8.0*LmbdGG_calculated, 0) );	// 8 flavors
			Set_Of_Model_Parameters.set_block_entry( "gravity", 2, buffer_complex[1]*complex<double>(8.0*LmbdGG_calculated/Mass_4l/Mass_4l, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 3, buffer_complex[2]*complex<double>(8.0*LmbdGG_calculated/Mass_4l/Mass_4l, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 4, buffer_complex[3]*complex<double>(8.0*LmbdGG_calculated/Mass_4l/Mass_4l, 0) );	// old dummy scale factor a_s(100 GeV)^2/2Pi/m4l
			Set_Of_Model_Parameters.set_block_entry( "gravity", 5, buffer_complex[4]*complex<double>(8.0*LmbdGG_calculated*Mass_4l*Mass_4l, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 6, buffer_complex[5]*complex<double>(8.0*LmbdGG_calculated/Mass_4l/Mass_4l, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 7, buffer_complex[6]*complex<double>(8.0*LmbdGG_calculated, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 8, buffer_complex[7]*complex<double>(8.0*LmbdGG_calculated/Mass_4l/Mass_4l, 0) );	// old dummy scale factor a_s(100 GeV)^2/2Pi/m4l
			Set_Of_Model_Parameters.set_block_entry( "gravity", 9, buffer_complex[8]*complex<double>(8.0*LmbdGG_calculated, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity",10, buffer_complex[9]*complex<double>(8.0*LmbdGG_calculated/Mass_4l/Mass_4l, 0) );
		}
		else
		{
			Set_Of_Model_Parameters.set_block_entry( "gravity", 1, buffer_complex[0]*complex<double>(8.0*LmbdGG_calculated, 0) );	// 8 flavors
			Set_Of_Model_Parameters.set_block_entry( "gravity", 2, buffer_complex[1]*complex<double>(8.0*LmbdGG_calculated/Higgs_mass/Higgs_mass, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 3, buffer_complex[2]*complex<double>(8.0*LmbdGG_calculated/Higgs_mass/Higgs_mass, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 4, buffer_complex[3]*complex<double>(8.0*LmbdGG_calculated/Higgs_mass/Higgs_mass, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 5, buffer_complex[4]*complex<double>(8.0*LmbdGG_calculated*Higgs_mass*Higgs_mass, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 6, buffer_complex[5]*complex<double>(8.0*LmbdGG_calculated/Higgs_mass/Higgs_mass, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 7, buffer_complex[6]*complex<double>(8.0*LmbdGG_calculated, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 8, buffer_complex[7]*complex<double>(8.0*LmbdGG_calculated/Higgs_mass/Higgs_mass, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 9, buffer_complex[8]*complex<double>(8.0*LmbdGG_calculated, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity",10, buffer_complex[9]*complex<double>(8.0*LmbdGG_calculated/Higgs_mass/Higgs_mass, 0) );
			
		}
			
		//qq
		params_rhod21 = buffer_complex[0]*complex<double>(LmbdGG_calculated, 0);
		params_rhos21 = buffer_complex[0]*complex<double>(LmbdGG_calculated, 0);
		params_rhob21 = buffer_complex[0]*complex<double>(LmbdGG_calculated, 0);
		params_rhou21 = buffer_complex[0]*complex<double>(LmbdGG_calculated, 0);
		params_rhoc21 = buffer_complex[0]*complex<double>(LmbdGG_calculated, 0);
		params_rhod22 = buffer_complex[1]*complex<double>(LmbdGG_calculated, 0);
		params_rhos22 = buffer_complex[1]*complex<double>(LmbdGG_calculated, 0);
		params_rhob22 = buffer_complex[1]*complex<double>(LmbdGG_calculated, 0);
		params_rhou22 = buffer_complex[1]*complex<double>(LmbdGG_calculated, 0);
		params_rhoc22 = buffer_complex[1]*complex<double>(LmbdGG_calculated, 0);
		params_rhod23 = buffer_complex[2]*complex<double>(LmbdGG_calculated, 0);
		params_rhos23 = buffer_complex[2]*complex<double>(LmbdGG_calculated, 0);
		params_rhob23 = buffer_complex[2]*complex<double>(LmbdGG_calculated, 0);
		params_rhou23 = buffer_complex[2]*complex<double>(LmbdGG_calculated, 0);
		params_rhoc23 = buffer_complex[2]*complex<double>(LmbdGG_calculated, 0);
		params_rhod24 = buffer_complex[3]*complex<double>(LmbdGG_calculated, 0);
		params_rhos24 = buffer_complex[3]*complex<double>(LmbdGG_calculated, 0);
		params_rhob24 = buffer_complex[3]*complex<double>(LmbdGG_calculated, 0);
		params_rhou24 = buffer_complex[3]*complex<double>(LmbdGG_calculated, 0);
		params_rhoc24 = buffer_complex[3]*complex<double>(LmbdGG_calculated, 0);
		
		// Decay to ZZ
		if( Use_mh_eq_m4l )
		{
			Set_Of_Model_Parameters.set_block_entry( "gravity", 11, buffer_complex[10]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/sqrt(2), 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 12, buffer_complex[11]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/Mass_4l/Mass_4l, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 13, buffer_complex[12]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/params_m_Z, 0) );	// or /mZ^4?
			Set_Of_Model_Parameters.set_block_entry( "gravity", 14, buffer_complex[13]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/params_m_Z, 0) );	// or /mZ^4?
			Set_Of_Model_Parameters.set_block_entry( "gravity", 15, buffer_complex[14]*complex<double>(hZZ_coupling/2/sqrt(2), 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 16, buffer_complex[15]*complex<double>(hZZ_coupling/2/params_m_Z/Mass_4l, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 17, buffer_complex[16]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 18, buffer_complex[17]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/params_m_Z, 0) );	// or /mZ^4?
			Set_Of_Model_Parameters.set_block_entry( "gravity", 19, buffer_complex[18]*complex<double>(hZZ_coupling/2/params_m_Z/Mass_4l, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 20, buffer_complex[19]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/Mass_4l/Mass_4l, 0) );
			
		}
		else
		{
			Set_Of_Model_Parameters.set_block_entry( "gravity", 11, buffer_complex[10]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/sqrt(2), 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 12, buffer_complex[11]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/Higgs_mass/Higgs_mass, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 13, buffer_complex[12]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/params_m_Z, 0) );	// or /mZ^4?
			Set_Of_Model_Parameters.set_block_entry( "gravity", 14, buffer_complex[13]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/params_m_Z, 0) );	// or /mZ^4?
			Set_Of_Model_Parameters.set_block_entry( "gravity", 15, buffer_complex[14]*complex<double>(hZZ_coupling/2/sqrt(2), 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 16, buffer_complex[15]*complex<double>(hZZ_coupling/2/params_m_Z/Higgs_mass, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 17, buffer_complex[16]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 18, buffer_complex[17]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/params_m_Z, 0) );	// or /mZ^4?
			Set_Of_Model_Parameters.set_block_entry( "gravity", 19, buffer_complex[18]*complex<double>(hZZ_coupling/2/params_m_Z/Higgs_mass, 0) );
			Set_Of_Model_Parameters.set_block_entry( "gravity", 20, buffer_complex[19]*complex<double>(hZZ_coupling/2/params_m_Z/params_m_Z/Higgs_mass/Higgs_mass, 0) );
		}
		
		// Decay to 2l (or 2mu in 2f)
		Set_Of_Model_Parameters.set_block_entry( "gravity", 41, buffer_complex[10]*complex<double>(4.291210e-04, 0) );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 42, buffer_complex[11]*complex<double>(4.291210e-04, 0) );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 43, buffer_complex[12]*complex<double>(4.291210e-04, 0) );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 44, buffer_complex[13]*complex<double>(4.291210e-04, 0) );
		
		// Decay to 2e (2f)
		Set_Of_Model_Parameters.set_block_entry( "gravity", 45, buffer_complex[10]*complex<double>(4.291210e-04, 0) );	// for Hee should be 2.075371e-06
		Set_Of_Model_Parameters.set_block_entry( "gravity", 46, buffer_complex[11]*complex<double>(4.291210e-04, 0) );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 47, buffer_complex[12]*complex<double>(4.291210e-04, 0) );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 48, buffer_complex[13]*complex<double>(4.291210e-04, 0) );
	}
	
	buffer_complex = NULL;
	
	return Run_MEKD_MG_ME_Dispatcher_SIG_Spin2( initial_state );
}



/// ME_RAW (RAW MG5_aMC ME) dispatcher
int MEKD_MG::Run_MEKD_MG_ME_Dispatcher_CPPProcess(string initial_state)
{
// 	if( Resonance_decay_mode=="ZZ" )
// 	{
// 		if( Final_state=="4e" || Final_state=="4eA" )
// 		{
// // 			/// Common mass for the same-flavor leptons
// // 			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_e );
// 			
// 			if( initial_state=="NO" && Final_state=="4e" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_RAW );
// 			
// 			if( initial_state=="NO" && Final_state=="4eA" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_RAW );
// 			
// 			if( initial_state=="gg" && Final_state=="4e" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_RAW );
// 			
// 			if( initial_state=="gg" && Final_state=="4eA" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_RAW );
// 			
// 			if( initial_state=="qq" && Final_state=="4e" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_RAW, ME_RAW );
// 			
// 			if( initial_state=="qq" && Final_state=="4eA" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_RAW, ME_RAW );
// 		}
// 		
// 		if( Final_state=="2e2m" || Final_state=="2e2mu" || Final_state=="2e2mA" || Final_state=="2e2muA" )
// 		{
// // 			/// Common mass for the opposite-flavor leptons
// // 			Set_Of_Model_Parameters.set_block_entry( "mass", 11, params_m_e );
// // 			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
// 			
// 			if( initial_state=="NO" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_RAW );
// 			
// 			if( initial_state=="NO" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_RAW );
// 			
// 			if( initial_state=="gg" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_RAW );
// 			
// 			if( initial_state=="gg" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_RAW );
// 			
// 			if( initial_state=="qq" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_RAW, ME_RAW );
// 			
// 			if( initial_state=="qq" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_RAW, ME_RAW );
// 		}
// 		
// 		if( Final_state=="4m" || Final_state=="4mu" || Final_state=="4mA" || Final_state=="4muA" )
// 		{
// // 			/// Common mass for the same-flavor leptons
// // 			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
// 			
// 			if( initial_state=="NO" && (Final_state=="4m" || Final_state=="4mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_RAW );
// 			
// 			if( initial_state=="NO" && (Final_state=="4mA" || Final_state=="4muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_RAW );
// 			
// 			if( initial_state=="gg" && (Final_state=="4m" || Final_state=="4mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_RAW );
// 			
// 			if( initial_state=="gg" && (Final_state=="4mA" || Final_state=="4muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_RAW );
// 			
// 			if( initial_state=="qq" && (Final_state=="4m" || Final_state=="4mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_RAW, ME_RAW );
// 			
// 			if( initial_state=="qq" && (Final_state=="4mA" || Final_state=="4muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_RAW, ME_RAW );
// 		}
// 	}
// 	
// 	if( Resonance_decay_mode=="2l" )
// 	{
// 		if( Final_state=="4e" || Final_state=="4eA" )
// 		{
// // 			/// Common mass for the same-flavor leptons
// // 			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_e );
// 			
// 			if( initial_state=="NO" && Final_state=="4e" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_RAW );
// 			
// 			if( initial_state=="NO" && Final_state=="4eA" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_RAW );
// 			
// 			if( initial_state=="gg" && Final_state=="4e" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_RAW );
// 			
// 			if( initial_state=="gg" && Final_state=="4eA" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_RAW );
// 			
// 			if( initial_state=="qq" && Final_state=="4e" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_RAW, ME_RAW );
// 			
// 			if( initial_state=="qq" && Final_state=="4eA" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_RAW, ME_RAW );
// 		}
// 		
// 		if( Final_state=="2e2m" || Final_state=="2e2mu" || Final_state=="2e2mA" || Final_state=="2e2muA" )
// 		{
// // 			/// Common mass for the opposite-flavor leptons
// // 			Set_Of_Model_Parameters.set_block_entry( "mass", 11, params_m_e );
// // 			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
// 			
// 			if( initial_state=="NO" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_RAW );
// 			
// 			if( initial_state=="NO" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_RAW );
// 			
// 			if( initial_state=="gg" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_RAW );
// 			
// 			if( initial_state=="gg" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_RAW );
// 			
// 			if( initial_state=="qq" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_RAW, ME_RAW );
// 			
// 			if( initial_state=="qq" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_RAW, ME_RAW );
// 		}
// 		
// 		if( Final_state=="4m" || Final_state=="4mu" || Final_state=="4mA" || Final_state=="4muA" )
// 		{
// // 			/// Common mass for the same-flavor leptons
// // 			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
// 			
// 			if( initial_state=="NO" && (Final_state=="4m" || Final_state=="4mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_RAW );
// 			
// 			if( initial_state=="NO" && (Final_state=="4mA" || Final_state=="4muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_RAW );
// 			
// 			if( initial_state=="gg" && (Final_state=="4m" || Final_state=="4mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_RAW );
// 			
// 			if( initial_state=="gg" && (Final_state=="4mA" || Final_state=="4muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_RAW );
// 			
// 			if( initial_state=="qq" && (Final_state=="4m" || Final_state=="4mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_RAW, ME_RAW );
// 			
// 			if( initial_state=="qq" && (Final_state=="4mA" || Final_state=="4muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_RAW, ME_RAW );
// 		}
// 	}
// 	
// 	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" )
// 	{
// // 		/// Mass for the muons
// // 		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
// 		
// 		if( initial_state=="NO" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_RAW );
// 		
// 		if( initial_state=="NO" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_RAW );
// 		
// 		if( initial_state=="gg" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_RAW );
// 		
// 		if( initial_state=="gg" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_RAW );
// 		
// 		if( initial_state=="qq" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_RAW, ME_RAW );
// 		
// 		if( initial_state=="qq" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_RAW, ME_RAW );
// 	}
	
	return 1;
}



/// ZZ and DY ME dispatcher
int MEKD_MG::Run_MEKD_MG_ME_Dispatcher_BKG_ZZ(string initial_state)
{
	if( initial_state=="NO" && !(Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA") )
	{
		cerr << "ZZ gg initial state is not supported.\n";
		return 1;
	}
	if( initial_state=="gg" )
	{
		cerr << "ZZ productionless state is not supported.\n";
		return 1;
	}
	
	if( Final_state=="4e" || Final_state=="4eA" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_e );
		
// 		if( initial_state=="NO" && Final_state=="4e" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_ZZ_SF );
// 		
// 		if( initial_state=="NO" && Final_state=="4eA" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_ZZ_SFpA );
// 		
// 		if( initial_state=="gg" && Final_state=="4e" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_ZZ_SF );
// 		
// 		if( initial_state=="gg" && Final_state=="4eA" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_ZZ_SFpA );
		
		if( initial_state=="qq" && Final_state=="4e" )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_ZZ_DownType_SF, ME_qq_ZZ_UpType_SF );
		
		if( initial_state=="qq" && Final_state=="4eA" )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_ZZ_DownType_SFpA, ME_qq_ZZ_UpType_SFpA );
	}
	
	if( Final_state=="2e2m" || Final_state=="2e2mu" || Final_state=="2e2mA" || Final_state=="2e2muA" )
	{
		/// Common mass for the opposite-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 11, params_m_e );
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
		
// 		if( initial_state=="NO" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_ZZ_OF );
// 		
// 		if( initial_state=="NO" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_ZZ_OFpA );
// 		
// 		if( initial_state=="gg" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_ZZ_OF );
// 		
// 		if( initial_state=="gg" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_ZZ_OFpA );
		
		if( initial_state=="qq" && (Final_state=="2e2m" || Final_state=="2e2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_ZZ_DownType_OF, ME_qq_ZZ_UpType_OF );
		
		if( initial_state=="qq" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_ZZ_DownType_OFpA, ME_qq_ZZ_UpType_OFpA );
	}
	
	if( Final_state=="4m" || Final_state=="4mu" || Final_state=="4mA" || Final_state=="4muA" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
		
// 		if( initial_state=="NO" && (Final_state=="4m" || Final_state=="4mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_ZZ_SF );
// 		
// 		if( initial_state=="NO" && (Final_state=="4mA" || Final_state=="4muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_ZZ_SFpA );
// 		
// 		if( initial_state=="gg" && (Final_state=="4m" || Final_state=="4mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_ZZ_SF );
// 		
// 		if( initial_state=="gg" && (Final_state=="4mA" || Final_state=="4muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_ZZ_SFpA );
		
		if( initial_state=="qq" && (Final_state=="4m" || Final_state=="4mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_ZZ_DownType_SF, ME_qq_ZZ_UpType_SF );
		
		if( initial_state=="qq" && (Final_state=="4mA" || Final_state=="4muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_ZZ_DownType_SFpA, ME_qq_ZZ_UpType_SFpA );
	}
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" )
	{
		/// Mass for the muons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
		
		if( initial_state=="NO" && (Final_state=="2m" || Final_state=="2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_DY_2l );
		
		if( initial_state=="NO" && (Final_state=="2mA" || Final_state=="2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_DY_2lpA );
// 		
// 		if( initial_state=="gg" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_DY_2l );
// 		
// 		if( initial_state=="gg" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_DY_2lpA );
		
		if( initial_state=="qq" && (Final_state=="2m" || Final_state=="2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_DY_DownType_2l, ME_qq_DY_UpType_2l );
		
		if( initial_state=="qq" && (Final_state=="2mA" || Final_state=="2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_DY_DownType_2lpA, ME_qq_DY_UpType_2lpA );
	}
	
	return 1;
}



/// Z4l Background ME dispatcher
int MEKD_MG::Run_MEKD_MG_ME_Dispatcher_Z4l_BKG(string initial_state)
{
	if( initial_state=="NO" )
	{
		cerr << "Z -> 4l gg initial state is not supported.\n";
		return 1;
	}
	if( initial_state=="gg" )
	{
		cerr << "Z -> 4l productionless state is not supported.\n";
		return 1;
	}
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" )
	{
		cerr << "Z -> 4l to 2 lepton is pointless. Run Drell-Yan (DY) instead (ZZ-based).\n";
		return 1;
	}
	
	if( Final_state=="4e" || Final_state=="4eA" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_e );
		
// 		if( initial_state=="NO" && Final_state=="4e" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Z4l_BKG_SF );
// 		
// 		if( initial_state=="NO" && Final_state=="4eA" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Z4l_BKG_SFpA );
// 		
// 		if( initial_state=="gg" && Final_state=="4e" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_Z4l_BKG_SF );
// 		
// 		if( initial_state=="gg" && Final_state=="4eA" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_Z4l_BKG_SFpA );
		
		if( initial_state=="qq" && Final_state=="4e" )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_Z4l_BKG_DownType_SF, ME_qq_Z4l_BKG_UpType_SF );
		
		if( initial_state=="qq" && Final_state=="4eA" )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_Z4l_BKG_DownType_SFpA, ME_qq_Z4l_BKG_UpType_SFpA );
	}
	
	if( Final_state=="2e2m" || Final_state=="2e2mu" || Final_state=="2e2mA" || Final_state=="2e2muA" )
	{
		/// Common mass for the opposite-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 11, params_m_e );
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
		
// 		if( initial_state=="NO" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Z4l_BKG_OF );
// 		
// 		if( initial_state=="NO" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Z4l_BKG_OFpA );
// 		
// 		if( initial_state=="gg" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_Z4l_BKG_OF );
// 		
// 		if( initial_state=="gg" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_Z4l_BKG_OFpA );
		
		if( initial_state=="qq" && (Final_state=="2e2m" || Final_state=="2e2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_Z4l_BKG_DownType_OF, ME_qq_Z4l_BKG_UpType_OF );
		
		if( initial_state=="qq" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_Z4l_BKG_DownType_OFpA, ME_qq_Z4l_BKG_UpType_OFpA );
	}
	
	if( Final_state=="4m" || Final_state=="4mu" || Final_state=="4mA" || Final_state=="4muA" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
		
// 		if( initial_state=="NO" && (Final_state=="4m" || Final_state=="4mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Z4l_BKG_SF );
// 		
// 		if( initial_state=="NO" && (Final_state=="4mA" || Final_state=="4muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Z4l_BKG_SFpA );
// 		
// 		if( initial_state=="gg" && (Final_state=="4m" || Final_state=="4mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_Z4l_BKG_SF );
// 		
// 		if( initial_state=="gg" && (Final_state=="4mA" || Final_state=="4muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_Z4l_BKG_SFpA );
		
		if( initial_state=="qq" && (Final_state=="4m" || Final_state=="4mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_Z4l_BKG_DownType_SF, ME_qq_Z4l_BKG_UpType_SF );
		
		if( initial_state=="qq" && (Final_state=="4mA" || Final_state=="4muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_Z4l_BKG_DownType_SFpA, ME_qq_Z4l_BKG_UpType_SFpA );
	}
	
// 	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" )
// 	{
// 		/// Mass for the muons
// 		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
// 		
// 		if( initial_state=="NO" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Z4l_BKG_2l );
// 		
// 		if( initial_state=="NO" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Z4l_BKG_2lpA );
// 		
// 		if( initial_state=="gg" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_Z4l_BKG__2l );
// 		
// 		if( initial_state=="gg" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_Z4l_BKG__2lpA );
// 		
// 		if( initial_state=="qq" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_Z4l_BKG_DownType_2l, ME_qq_Z4l_BKG_UpType_2l );
// 		
// 		if( initial_state=="qq" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_Z4l_BKG_DownType_2lpA, ME_qq_Z4l_BKG_UpType_2lpA );
// 	}
	
	return 1;
}



/// Z4l Signal ME dispatcher
int MEKD_MG::Run_MEKD_MG_ME_Dispatcher_Z4l_SIG(string initial_state)
{
	if( initial_state=="NO" )
	{
		cerr << "Z -> 4l gg initial state is not supported.\n";
		return 1;
	}
	if( initial_state=="gg" )
	{
		cerr << "Z -> 4l productionless state is not supported.\n";
		return 1;
	}
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" )
	{
		cerr << "Z -> 4l to 2 lepton is pointless. Run Drell-Yan (DY) instead (ZZ-based).\n";
		return 1;
	}
	
	if( Final_state=="4e" || Final_state=="4eA" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_e );
		
// 		if( initial_state=="NO" && Final_state=="4e" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Z4l_SIG_SF );
// 		
// 		if( initial_state=="NO" && Final_state=="4eA" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Z4l_SIG_SFpA );
// 		
// 		if( initial_state=="gg" && Final_state=="4e" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_Z4l_SIG_SF );
// 		
// 		if( initial_state=="gg" && Final_state=="4eA" )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_Z4l_SIG_SFpA );
		
		if( initial_state=="qq" && Final_state=="4e" )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_Z4l_SIG_DownType_SF, ME_qq_Z4l_SIG_UpType_SF );
		
		if( initial_state=="qq" && Final_state=="4eA" )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_Z4l_SIG_DownType_SFpA, ME_qq_Z4l_SIG_UpType_SFpA );
	}
	
	if( Final_state=="2e2m" || Final_state=="2e2mu" || Final_state=="2e2mA" || Final_state=="2e2muA" )
	{
		/// Common mass for the opposite-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 11, params_m_e );
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
		
// 		if( initial_state=="NO" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Z4l_SIG_OF );
// 		
// 		if( initial_state=="NO" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Z4l_SIG_OFpA );
// 		
// 		if( initial_state=="gg" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_Z4l_SIG_OF );
// 		
// 		if( initial_state=="gg" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_Z4l_SIG_OFpA );
		
		if( initial_state=="qq" && (Final_state=="2e2m" || Final_state=="2e2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_Z4l_SIG_DownType_OF, ME_qq_Z4l_SIG_UpType_OF );
		
		if( initial_state=="qq" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_Z4l_SIG_DownType_OFpA, ME_qq_Z4l_SIG_UpType_OFpA );
	}
	
	if( Final_state=="4m" || Final_state=="4mu" || Final_state=="4mA" || Final_state=="4muA" )
	{
		/// Common mass for the same-flavor leptons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
		
// 		if( initial_state=="NO" && (Final_state=="4m" || Final_state=="4mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Z4l_SIG_SF );
// 		
// 		if( initial_state=="NO" && (Final_state=="4mA" || Final_state=="4muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Z4l_SIG_SFpA );
// 		
// 		if( initial_state=="gg" && (Final_state=="4m" || Final_state=="4mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_Z4l_SIG_SF );
// 		
// 		if( initial_state=="gg" && (Final_state=="4mA" || Final_state=="4muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_Z4l_SIG_SFpA );
		
		if( initial_state=="qq" && (Final_state=="4m" || Final_state=="4mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_Z4l_SIG_DownType_SF, ME_qq_Z4l_SIG_UpType_SF );
		
		if( initial_state=="qq" && (Final_state=="4mA" || Final_state=="4muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_Z4l_SIG_DownType_SFpA, ME_qq_Z4l_SIG_UpType_SFpA );
	}
	
// 	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" )
// 	{
// 		/// Mass for the muons
// 		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
// 		
// 		if( initial_state=="NO" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Z4l_SIG_2l );
// 		
// 		if( initial_state=="NO" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Z4l_SIG_2lpA );
// 		
// 		if( initial_state=="gg" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_gg_Z4l_SIG__2l );
// 		
// 		if( initial_state=="gg" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_gg_Z4l_SIG__2lpA );
// 		
// 		if( initial_state=="qq" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_qq_Z4l_SIG_DownType_2l, ME_qq_Z4l_SIG_UpType_2l );
// 		
// 		if( initial_state=="qq" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_qq_Z4l_SIG_DownType_2lpA, ME_qq_Z4l_SIG_UpType_2lpA );
// 	}
	
	return 1;
}



/// Spin-0 ME dispatcher
int MEKD_MG::Run_MEKD_MG_ME_Dispatcher_SIG_Spin0(string initial_state)
{
	if( initial_state=="qq" )
	{
		cerr << "Spin-0 qqbar initial state is redundant, thus not provided.\n";
		return 1;
	}
	
	if( Resonance_decay_mode=="ZZ" )
	{
		if( Final_state=="4e" || Final_state=="4eA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_e );
			
			if( initial_state=="NO" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin0_SF );
			
			if( initial_state=="NO" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin0_SFpA );
			
			if( initial_state=="gg" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin0_SF );
			
			if( initial_state=="gg" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin0_SFpA );
			
// 			if( initial_state=="qq" && Final_state=="4e" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin0_DownType_SF, ME_Signal_qq_Spin0_UpType_SF );
// 			
// 			if( initial_state=="qq" && Final_state=="4eA" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin0_DownType_SFpA, ME_Signal_qq_Spin0_UpType_SFpA );
		}
		
		if( Final_state=="2e2m" || Final_state=="2e2mu" || Final_state=="2e2mA" || Final_state=="2e2muA" )
		{
			/// Common mass for the opposite-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 11, params_m_e );
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin0_OF );
			
			if( initial_state=="NO" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin0_OFpA );
			
			if( initial_state=="gg" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin0_OF );
			
			if( initial_state=="gg" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin0_OFpA );
			
// 			if( initial_state=="qq" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin0_DownType_OF, ME_Signal_qq_Spin0_UpType_OF );
// 			
// 			if( initial_state=="qq" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin0_DownType_OFpA, ME_Signal_qq_Spin0_UpType_OFpA );
		}
		
		if( Final_state=="4m" || Final_state=="4mu" || Final_state=="4mA" || Final_state=="4muA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin0_SF );
			
			if( initial_state=="NO" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin0_SFpA );
			
			if( initial_state=="gg" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin0_SF );
			
			if( initial_state=="gg" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin0_SFpA );
			
// 			if( initial_state=="qq" && (Final_state=="4m" || Final_state=="4mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin0_DownType_SF, ME_Signal_qq_Spin0_UpType_SF );
// 			
// 			if( initial_state=="qq" && (Final_state=="4mA" || Final_state=="4muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin0_DownType_SFpA, ME_Signal_qq_Spin0_UpType_SFpA );
		}
	}
	
	if( Resonance_decay_mode=="2l" )
	{
		if( Final_state=="4e" || Final_state=="4eA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_e );
			
			if( initial_state=="NO" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin0_2f_SF );
			
			if( initial_state=="NO" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin0_2f_SFpA );
			
			if( initial_state=="gg" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin0_2f_SF );
			
			if( initial_state=="gg" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin0_2f_SFpA );
			
// 			if( initial_state=="qq" && Final_state=="4e" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin0_2f_DownType_SF, ME_Signal_qq_Spin0_2f_UpType_SF );
// 			
// 			if( initial_state=="qq" && Final_state=="4eA" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin0_2f_DownType_SFpA, ME_Signal_qq_Spin0_2f_UpType_SFpA );
		}
		
		if( Final_state=="2e2m" || Final_state=="2e2mu" || Final_state=="2e2mA" || Final_state=="2e2muA" )
		{
			/// Common mass for the opposite-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 11, params_m_e );
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin0_2f_OF );
			
			if( initial_state=="NO" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin0_2f_OFpA );
			
			if( initial_state=="gg" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin0_2f_OF );
			
			if( initial_state=="gg" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin0_2f_OFpA );
			
// 			if( initial_state=="qq" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin0_2f_DownType_OF, ME_Signal_qq_Spin0_2f_UpType_OF );
// 			
// 			if( initial_state=="qq" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin0_2f_DownType_OFpA, ME_Signal_qq_Spin0_2f_UpType_OFpA );
		}
		
		if( Final_state=="4m" || Final_state=="4mu" || Final_state=="4mA" || Final_state=="4muA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin0_2f_SF );
			
			if( initial_state=="NO" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin0_2f_SFpA );
			
			if( initial_state=="gg" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin0_2f_SF );
			
			if( initial_state=="gg" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin0_2f_SFpA );
			
// 			if( initial_state=="qq" && (Final_state=="4m" || Final_state=="4mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin0_2f_DownType_SF, ME_Signal_qq_Spin0_2f_UpType_SF );
// 			
// 			if( initial_state=="qq" && (Final_state=="4mA" || Final_state=="4muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin0_2f_DownType_SFpA, ME_Signal_qq_Spin0_2f_UpType_SFpA );
		}
	}
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" )
	{
		/// Mass for the muons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
		
		if( initial_state=="NO" && (Final_state=="2m" || Final_state=="2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin0_2l );
		
		if( initial_state=="NO" && (Final_state=="2mA" || Final_state=="2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin0_2lpA );
		
		if( initial_state=="gg" && (Final_state=="2m" || Final_state=="2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin0_2l );
		
		if( initial_state=="gg" && (Final_state=="2mA" || Final_state=="2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin0_2lpA );
		
// 		if( initial_state=="qq" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin0_DownType_2l, ME_Signal_qq_Spin0_UpType_2l );
		
// 		if( initial_state=="qq" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin0_DownType_2lpA, ME_Signal_qq_Spin0_UpType_2lpA );
	}
	
	return 1;
}



/// Spin-1 ME dispatcher
int MEKD_MG::Run_MEKD_MG_ME_Dispatcher_SIG_Spin1(string initial_state)
{
	if( initial_state=="gg" )
	{
		cerr << "Spin-1 gg initial state is not possible, thus not provided.\n";
		return 1;
	}
	
	if( Resonance_decay_mode=="ZZ" )
	{
		if( Final_state=="4e" || Final_state=="4eA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_e );
			
			if( initial_state=="NO" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin1_SF );
			
			if( initial_state=="NO" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin1_SFpA );
			
// 			if( initial_state=="gg" && Final_state=="4e" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin1_SF );
// 			
// 			if( initial_state=="gg" && Final_state=="4eA" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin1_SFpA );
			
			if( initial_state=="qq" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin1_DownType_SF, ME_Signal_qq_Spin1_UpType_SF );
			
			if( initial_state=="qq" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin1_DownType_SFpA, ME_Signal_qq_Spin1_UpType_SFpA );
		}
		
		if( Final_state=="2e2m" || Final_state=="2e2mu" || Final_state=="2e2mA" || Final_state=="2e2muA" )
		{
			/// Common mass for the opposite-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 11, params_m_e );
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin1_OF );
			
			if( initial_state=="NO" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin1_OFpA );
			
//	 		if( initial_state=="gg" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin1_OF );
// 			
// 			if( initial_state=="gg" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin1_OFpA );
			
			if( initial_state=="qq" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin1_DownType_OF, ME_Signal_qq_Spin1_UpType_OF );
			
			if( initial_state=="qq" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin1_DownType_OFpA, ME_Signal_qq_Spin1_UpType_OFpA );
		}
		
		if( Final_state=="4m" || Final_state=="4mu" || Final_state=="4mA" || Final_state=="4muA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin1_SF );
			
			if( initial_state=="NO" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin1_SFpA );
			
// 			if( initial_state=="gg" && (Final_state=="4m" || Final_state=="4mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin1_SF );
// 			
// 			if( initial_state=="gg" && (Final_state=="4mA" || Final_state=="4muA") )
//	 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin1_SFpA );
			
			if( initial_state=="qq" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin1_DownType_SF, ME_Signal_qq_Spin1_UpType_SF );
			
			if( initial_state=="qq" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin1_DownType_SFpA, ME_Signal_qq_Spin1_UpType_SFpA );
		}
	}
	
	if( Resonance_decay_mode=="2l" )
	{
		if( Final_state=="4e" || Final_state=="4eA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_e );
			
			if( initial_state=="NO" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin1_2f_SF );
			
			if( initial_state=="NO" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin1_2f_SFpA );
			
// 			if( initial_state=="gg" && Final_state=="4e" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin1_2f_SF );
// 			
// 			if( initial_state=="gg" && Final_state=="4eA" )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin1_2f_SFpA );
			
			if( initial_state=="qq" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin1_2f_DownType_SF, ME_Signal_qq_Spin1_2f_UpType_SF );
			
			if( initial_state=="qq" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin1_2f_DownType_SFpA, ME_Signal_qq_Spin1_2f_UpType_SFpA );
		}
		
		if( Final_state=="2e2m" || Final_state=="2e2mu" || Final_state=="2e2mA" || Final_state=="2e2muA" )
		{
			/// Common mass for the opposite-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 11, params_m_e );
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin1_2f_OF );
			
			if( initial_state=="NO" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin1_2f_OFpA );
			
//	 		if( initial_state=="gg" && (Final_state=="2e2m" || Final_state=="2e2mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin1_2f_OF );
// 			
// 			if( initial_state=="gg" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin1_2f_OFpA );
			
			if( initial_state=="qq" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin1_2f_DownType_OF, ME_Signal_qq_Spin1_2f_UpType_OF );
			
			if( initial_state=="qq" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin1_2f_DownType_OFpA, ME_Signal_qq_Spin1_2f_UpType_OFpA );
		}
		
		if( Final_state=="4m" || Final_state=="4mu" || Final_state=="4mA" || Final_state=="4muA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin1_2f_SF );
			
			if( initial_state=="NO" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin1_2f_SFpA );
			
// 			if( initial_state=="gg" && (Final_state=="4m" || Final_state=="4mu") )
// 				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin1_2f_SF );
// 			
// 			if( initial_state=="gg" && (Final_state=="4mA" || Final_state=="4muA") )
//	 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin1_2f_SFpA );
			
			if( initial_state=="qq" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin1_2f_DownType_SF, ME_Signal_qq_Spin1_2f_UpType_SF );
			
			if( initial_state=="qq" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin1_2f_DownType_SFpA, ME_Signal_qq_Spin1_2f_UpType_SFpA );
		}
	}
		
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" )
	{
		/// Mass for the muons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
		
		if( initial_state=="NO" && (Final_state=="2m" || Final_state=="2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin1_2l );
		
		if( initial_state=="NO" && (Final_state=="2mA" || Final_state=="2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin1_2lpA );
		
// 		if( initial_state=="gg" && (Final_state=="2m" || Final_state=="2mu") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin1_2l );
// 		
// 		if( initial_state=="gg" && (Final_state=="2mA" || Final_state=="2muA") )
// 			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin1_2lpA );
		
		if( initial_state=="qq" && (Final_state=="2m" || Final_state=="2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin1_DownType_2l, ME_Signal_qq_Spin1_UpType_2l );
		
		if( initial_state=="qq" && (Final_state=="2mA" || Final_state=="2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin1_DownType_2lpA, ME_Signal_qq_Spin1_UpType_2lpA );
	}
	
	return 1;
}



/// Spin-2 ME dispatcher
int MEKD_MG::Run_MEKD_MG_ME_Dispatcher_SIG_Spin2(string initial_state)
{
	if( Resonance_decay_mode=="ZZ" )
	{
		if( Final_state=="4e" || Final_state=="4eA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_e );
			
			if( initial_state=="NO" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin2_SF );
			
			if( initial_state=="NO" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin2_SFpA );
			
			if( initial_state=="gg" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin2_SF );
			
			if( initial_state=="gg" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin2_SFpA );
			
			if( initial_state=="qq" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin2_DownType_SF, ME_Signal_qq_Spin2_UpType_SF );
			
			if( initial_state=="qq" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin2_DownType_SFpA, ME_Signal_qq_Spin2_UpType_SFpA );
		}
		
		if( Final_state=="2e2m" || Final_state=="2e2mu" || Final_state=="2e2mA" || Final_state=="2e2muA" )
		{
			/// Common mass for the opposite-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 11, params_m_e );
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin2_OF );
			
			if( initial_state=="NO" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin2_OFpA );
			
			if( initial_state=="gg" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin2_OF );
			
			if( initial_state=="gg" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin2_OFpA );
			
			if( initial_state=="qq" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin2_DownType_OF, ME_Signal_qq_Spin2_UpType_OF );
			
			if( initial_state=="qq" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin2_DownType_OFpA, ME_Signal_qq_Spin2_UpType_OFpA );
		}
		
		if( Final_state=="4m" || Final_state=="4mu" || Final_state=="4mA" || Final_state=="4muA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin2_SF );
			
			if( initial_state=="NO" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin2_SFpA );
			
			if( initial_state=="gg" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin2_SF );
			
			if( initial_state=="gg" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin2_SFpA );
			
			if( initial_state=="qq" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin2_DownType_SF, ME_Signal_qq_Spin2_UpType_SF );
			
			if( initial_state=="qq" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin2_DownType_SFpA, ME_Signal_qq_Spin2_UpType_SFpA );
		}
	}
	
	if( Resonance_decay_mode=="2l" )
	{
		if( Final_state=="4e" || Final_state=="4eA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_e );
			
			if( initial_state=="NO" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin2_2f_SF );
			
			if( initial_state=="NO" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin2_2f_SFpA );
			
			if( initial_state=="gg" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin2_2f_SF );
			
			if( initial_state=="gg" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin2_2f_SFpA );
			
			if( initial_state=="qq" && Final_state=="4e" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin2_2f_DownType_SF, ME_Signal_qq_Spin2_2f_UpType_SF );
			
			if( initial_state=="qq" && Final_state=="4eA" )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin2_2f_DownType_SFpA, ME_Signal_qq_Spin2_2f_UpType_SFpA );
		}
		
		if( Final_state=="2e2m" || Final_state=="2e2mu" || Final_state=="2e2mA" || Final_state=="2e2muA" )
		{
			/// Common mass for the opposite-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 11, params_m_e );
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin2_2f_OF );
			
			if( initial_state=="NO" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin2_2f_OFpA );
			
			if( initial_state=="gg" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin2_2f_OF );
			
			if( initial_state=="gg" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin2_2f_OFpA );
			
			if( initial_state=="qq" && (Final_state=="2e2m" || Final_state=="2e2mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin2_2f_DownType_OF, ME_Signal_qq_Spin2_2f_UpType_OF );
			
			if( initial_state=="qq" && (Final_state=="2e2mA" || Final_state=="2e2muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin2_2f_DownType_OFpA, ME_Signal_qq_Spin2_2f_UpType_OFpA );
		}
		
		if( Final_state=="4m" || Final_state=="4mu" || Final_state=="4mA" || Final_state=="4muA" )
		{
			/// Common mass for the same-flavor leptons
			Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
			
			if( initial_state=="NO" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin2_2f_SF );
			
			if( initial_state=="NO" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin2_2f_SFpA );
			
			if( initial_state=="gg" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin2_2f_SF );
			
			if( initial_state=="gg" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin2_2f_SFpA );
			
			if( initial_state=="qq" && (Final_state=="4m" || Final_state=="4mu") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin2_2f_DownType_SF, ME_Signal_qq_Spin2_2f_UpType_SF );
			
			if( initial_state=="qq" && (Final_state=="4mA" || Final_state=="4muA") )
				return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin2_2f_DownType_SFpA, ME_Signal_qq_Spin2_2f_UpType_SFpA );
		}
	}
	
	if( Final_state=="2m" || Final_state=="2mu" || Final_state=="2mA" || Final_state=="2muA" )
	{
		/// Mass for the muons
		Set_Of_Model_Parameters.set_block_entry( "mass", 13, params_m_mu );
		
		if( initial_state=="NO" && (Final_state=="2m" || Final_state=="2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( false, ME_Signal_Spin2_2l );
		
		if( initial_state=="NO" && (Final_state=="2mA" || Final_state=="2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_NO( true, ME_Signal_Spin2_2lpA );
		
		if( initial_state=="gg" && (Final_state=="2m" || Final_state=="2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( false, ME_Signal_gg_Spin2_2l );
		
		if( initial_state=="gg" && (Final_state=="2mA" || Final_state=="2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_gg( true, ME_Signal_gg_Spin2_2lpA );
		
		if( initial_state=="qq" && (Final_state=="2m" || Final_state=="2mu") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( false, ME_Signal_qq_Spin2_DownType_2l, ME_Signal_qq_Spin2_UpType_2l );
		
		if( initial_state=="qq" && (Final_state=="2mA" || Final_state=="2muA") )
			return Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar( true, ME_Signal_qq_Spin2_DownType_2lpA, ME_Signal_qq_Spin2_UpType_2lpA );
	}
	
	return 1;
}



template<class Generic_MEKD_MG_ME>
int MEKD_MG::Run_MEKD_MG_MEs_Evaluator_Initial_State_NO(bool photon, Generic_MEKD_MG_ME &Generic_ME)
{
	Signal_ME = 0;
	
	/// No initial state block
	for( unsigned int count=0; count < 4; count++ ) buffer_p[count] = p_set[0][count];
	
	if( !photon )
	{
		p_set[0][0] = p_set[2][0] + p_set[3][0] + p_set[4][0] + p_set[5][0];
		p_set[0][1] = p_set[2][1] + p_set[3][1] + p_set[4][1] + p_set[5][1];
		p_set[0][2] = p_set[2][2] + p_set[3][2] + p_set[4][2] + p_set[5][2];
		p_set[0][3] = p_set[2][3] + p_set[3][3] + p_set[4][3] + p_set[5][3];
		
		buffer = p_set[1];
		for( unsigned int count=0; count < 4; count++ ) p_set[count+1] = p_set[count+2];
		p_set[5] = buffer;
	}
	else
	{
		p_set[0][0] = p_set[2][0] + p_set[3][0] + p_set[4][0] + p_set[5][0] + p_set[6][0];
		p_set[0][1] = p_set[2][1] + p_set[3][1] + p_set[4][1] + p_set[5][1] + p_set[6][1];
		p_set[0][2] = p_set[2][2] + p_set[3][2] + p_set[4][2] + p_set[5][2] + p_set[6][2];
		p_set[0][3] = p_set[2][3] + p_set[3][3] + p_set[4][3] + p_set[5][3] + p_set[6][3];
		
		buffer = p_set[1];
		for( unsigned int count=0; count < 5; count++ ) p_set[count+1] = p_set[count+2];
		p_set[6] = buffer;
	}
	
	
	Generic_ME.updateProc( Set_Of_Model_Parameters );
	Generic_ME.setMomenta( p_set );
	Generic_ME.sigmaKin();
	buffer = const_cast<double*>( Generic_ME.getMatrixElements() );
	
	Signal_ME = buffer[0];
	
	
	if( !photon )
	{
		buffer = p_set[5];
		for( unsigned int count=4; count > 0; count-- ) p_set[count+1] = p_set[count];
		p_set[1] = buffer;
	}
	else
	{
		buffer = p_set[6];
		for( unsigned int count=5; count > 0; count-- ) p_set[count+1] = p_set[count];
		p_set[1] = buffer;
	}
	
	for( unsigned int count=0; count < 4; count++ ) p_set[0][count] = buffer_p[count];
	
	return 0;
}



template<class Generic_MEKD_MG_ME>
int MEKD_MG::Run_MEKD_MG_MEs_Evaluator_Initial_State_gg(bool photon, Generic_MEKD_MG_ME &Generic_ME)
{
	Signal_ME = 0;
	
	/// gg block
	p_set[0][3] = p_set[0][0];
	p_set[1][3] = -p_set[1][0];
	
	Generic_ME.updateProc( Set_Of_Model_Parameters );
	Generic_ME.setMomenta( p_set );
	Generic_ME.sigmaKin();
	buffer = const_cast<double*>( Generic_ME.getMatrixElements() );
	
	if( Use_PDF_w_pT0 ) { Signal_ME = pdfreader( 21, PDFx1, Mass_4l )*pdfreader( 21, PDFx2, Mass_4l )*buffer[0]; }
	else Signal_ME = buffer[0];
	
	return 0;
}



template<class Generic_MEKD_MG_ME_s, class Generic_MEKD_MG_ME_c>
int MEKD_MG::Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar(bool photon, Generic_MEKD_MG_ME_s &Generic_ME_s, Generic_MEKD_MG_ME_c &Generic_ME_c)
{
	Signal_ME = 0;
	
	/// Down quark block. Down type (s-like)
	if( ContributionCoeff_d != 0 )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 3, params_m_d );
// 		Set_Of_Model_Parameters.set_block_entry( "heff", 15, params_rhod01 );
// 		Set_Of_Model_Parameters.set_block_entry( "heff", 16, params_rhod02 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 15, params_rhod11 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 16, params_rhod12 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 17, params_rhod13 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 18, params_rhod14 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 33, params_rhod21 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 34, params_rhod22 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 35, params_rhod23 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 36, params_rhod24 );
		p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - params_m_d*params_m_d );
		p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - params_m_d*params_m_d );
		
		
		Generic_ME_s.updateProc( Set_Of_Model_Parameters );
		Generic_ME_s.setMomenta( p_set );
		Generic_ME_s.sigmaKin();
		buffer = const_cast<double*>( Generic_ME_s.getMatrixElements() );
		
		if( Use_PDF_w_pT0 )
		{
			ContributionCoeff_d = pdfreader( 1, PDFx1, Mass_4l )*pdfreader( -1, PDFx2, Mass_4l );
			Signal_ME = ContributionCoeff_d*buffer[0];
			ContributionCoeff_d = pdfreader( -1, PDFx1, Mass_4l )*pdfreader( 1, PDFx2, Mass_4l );
			Signal_ME += ContributionCoeff_d*buffer[1];
		}
		else Signal_ME = ContributionCoeff_d*(buffer[0]+buffer[1]);
	}
		
		
		
	/// Strange quark block. Down type (s-like)
	if( ContributionCoeff_s != 0 )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 3, params_m_s );
// 		Set_Of_Model_Parameters.set_block_entry( "heff", 15, params_rhos01 );
// 		Set_Of_Model_Parameters.set_block_entry( "heff", 16, params_rhos02 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 15, params_rhos11 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 16, params_rhos12 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 17, params_rhos13 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 18, params_rhos14 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 33, params_rhos21 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 34, params_rhos22 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 35, params_rhos23 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 36, params_rhos24 );
		p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - params_m_s*params_m_s );
		p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - params_m_s*params_m_s );
		
		Generic_ME_s.updateProc( Set_Of_Model_Parameters );
		Generic_ME_s.setMomenta( p_set );
		Generic_ME_s.sigmaKin();
		buffer = const_cast<double*>( Generic_ME_s.getMatrixElements() );
		
		if( Use_PDF_w_pT0 )
		{
			ContributionCoeff_s= pdfreader( 3, PDFx1, Mass_4l )*pdfreader( -3, PDFx2, Mass_4l );
			Signal_ME += ContributionCoeff_s*buffer[0];
			ContributionCoeff_s = pdfreader( -3, PDFx1, Mass_4l )*pdfreader( 3, PDFx2, Mass_4l );
			Signal_ME += ContributionCoeff_s*buffer[1];
		}
		else Signal_ME += ContributionCoeff_s*(buffer[0]+buffer[1]);
	}
		
		
		
	/// Up quark block. Up type (c-like)
	if( ContributionCoeff_u != 0 )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 4, params_m_u );
// 		Set_Of_Model_Parameters.set_block_entry( "heff", 11, params_rhou01 );
// 		Set_Of_Model_Parameters.set_block_entry( "heff", 12, params_rhou02 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 7, params_rhou11 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 8, params_rhou12 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 9, params_rhou13 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 10, params_rhou14 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 25, params_rhou21 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 26, params_rhou22 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 27, params_rhou23 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 28, params_rhou24 );
		p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - params_m_u*params_m_u );
		p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - params_m_u*params_m_u );
		
		
		Generic_ME_c.updateProc( Set_Of_Model_Parameters );
		Generic_ME_c.setMomenta( p_set );
		Generic_ME_c.sigmaKin();
		buffer = const_cast<double*>( Generic_ME_c.getMatrixElements() );
		
		if( Use_PDF_w_pT0 )
		{
			ContributionCoeff_u = pdfreader( 2, PDFx1, Mass_4l )*pdfreader( -2, PDFx2, Mass_4l );
			Signal_ME += ContributionCoeff_u*buffer[0];
			ContributionCoeff_u = pdfreader( -2, PDFx1, Mass_4l )*pdfreader( 2, PDFx2, Mass_4l );
			Signal_ME += ContributionCoeff_u*buffer[1];
		}
		else Signal_ME += ContributionCoeff_u*(buffer[0]+buffer[1]);
	}
		
		
		
	/// Charm quark block. Up type (c-like)
	if( ContributionCoeff_c != 0 )
	{
		Set_Of_Model_Parameters.set_block_entry( "mass", 4, params_m_c );
// 		Set_Of_Model_Parameters.set_block_entry( "heff", 11, params_rhoc01 );
// 		Set_Of_Model_Parameters.set_block_entry( "heff", 12, params_rhoc02 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 7, params_rhoc11 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 8, params_rhoc12 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 9, params_rhoc13 );
		Set_Of_Model_Parameters.set_block_entry( "vec", 10, params_rhoc14 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 25, params_rhoc21 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 26, params_rhoc22 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 27, params_rhoc23 );
		Set_Of_Model_Parameters.set_block_entry( "gravity", 28, params_rhoc24 );
		p_set[0][3] = sqrt( p_set[0][0]*p_set[0][0] - params_m_c*params_m_c );
		p_set[1][3] = -sqrt( p_set[1][0]*p_set[1][0] - params_m_c*params_m_c );
		
		
		Generic_ME_c.updateProc( Set_Of_Model_Parameters );
		Generic_ME_c.setMomenta( p_set );
		Generic_ME_c.sigmaKin();
		buffer = const_cast<double*>( Generic_ME_c.getMatrixElements() );
		
		if( Use_PDF_w_pT0 )
		{
			ContributionCoeff_c = pdfreader( 4, PDFx1, Mass_4l )*pdfreader( -4, PDFx2, Mass_4l );
			Signal_ME += ContributionCoeff_c*buffer[0];
			ContributionCoeff_c = pdfreader( -4, PDFx1, Mass_4l )*pdfreader( 4, PDFx2, Mass_4l );
			Signal_ME += ContributionCoeff_c*buffer[1];
		}
		else Signal_ME += ContributionCoeff_c*(buffer[0]+buffer[1]);
	}
	
	if( Use_mZ4l_eq_m4l ) Set_Of_Model_Parameters.set_block_entry( "mass", 23, params_m_Z );	//return to real mass. Used in Z -> 4l
	
	return 0;
}



///////////////////////////////////
/// END OF MEKD_MG_RunMEs.cpp   ///
///////////////////////////////////





///#include "MEKD_MG_Sorter.cpp"
///////////////////////////////////
/// INCLUDED MEKD_MG_Sorter.cpp ///
/// code follows below          ///
///                             ///
/// Part responsible for        ///
/// momenta rearrangement       ///
///////////////////////////////////



int MEKD_MG::Arrange_Internal_pls()
{
	id_set[0]=id1; id_set[1]=id2; id_set[2]=id3; id_set[3]=id4; id_set[4]=id5;
	
	if( id_set[0] == 0 ) id_set[0]=10000;
	if( id_set[1] == 0 ) id_set[1]=10000;
	if( id_set[2] == 0 ) id_set[2]=10000;
	if( id_set[3] == 0 ) id_set[3]=10000;
	if( id_set[4] == 0 ) id_set[4]=10000;
	sort( id_set.begin(), id_set.end() );
	
	
	////////////////////////////
	/// 2mu-decay-mode block ///
	////////////////////////////
	
	/// Two-lepton final state block
	if( id_set[0] == -13 && id_set[1] == 13 && id_set[2] == 10000 && id_set[3] == 10000 && id_set[4] == 10000 )
	{
		if( id1 == 13 ) pl1_internal = p1;
		if( id2 == 13 ) pl1_internal = p2;
		
		if( id1 == -13 ) pl2_internal = p1;
		if( id2 == -13 ) pl2_internal = p2;
		
		pA1_internal = NULL;
		Final_state = "2m";
		
		return 0;
	}
	
	/// Two-lepton + photon final state block
	if( id_set[0] == -13 && id_set[1] == 13 && id_set[2] == 22 && id_set[3] == 10000 && id_set[4] == 10000 )
	{
		if( id1 == 13 ) pl1_internal = p1;
		if( id2 == 13 ) pl1_internal = p2;
		if( id3 == 13 ) pl1_internal = p3;
		
		if( id1 == -13 ) pl2_internal = p1;
		if( id2 == -13 ) pl2_internal = p2;
		if( id3 == -13 ) pl2_internal = p3;
		
		if( id1 == 22 ) pA1_internal = p1;
		if( id2 == 22 ) pA1_internal = p2;
		if( id3 == 22 ) pA1_internal = p3;
		Final_state = "2mA";
		
		return 0;
	}
	
	
	///////////////////////////
	/// ZZ-decay-mode block ///
	///////////////////////////
	
	/// Four-lepton final state block
	if( id_set[0] == -13 && id_set[1] == -11 && id_set[2] == 11 && id_set[3] == 13 && id_set[4] == 10000 )
	{
		if( id1 == 11 ) pl1_internal = p1;
		if( id2 == 11 ) pl1_internal = p2;
		if( id3 == 11 ) pl1_internal = p3;
		if( id4 == 11 ) pl1_internal = p4;
		
		if( id1 == -11 ) pl2_internal = p1;
		if( id2 == -11 ) pl2_internal = p2;
		if( id3 == -11 ) pl2_internal = p3;
		if( id4 == -11 ) pl2_internal = p4;
		
		if( id1 == 13 ) pl3_internal = p1;
		if( id2 == 13 ) pl3_internal = p2;
		if( id3 == 13 ) pl3_internal = p3;
		if( id4 == 13 ) pl3_internal = p4;
		
		if( id1 == -13 ) pl4_internal = p1;
		if( id2 == -13 ) pl4_internal = p2;
		if( id3 == -13 ) pl4_internal = p3;
		if( id4 == -13 ) pl4_internal = p4;
		
		pA1_internal = NULL;
		Final_state = "2e2m";
		
		return 0;
	}
	
	if( id_set[0] == -13 && id_set[1] == -13 && id_set[2] == 13 && id_set[3] == 13 && id_set[4] == 10000 )
	{
		buffer_bool = false;	//first muon has beed caught
		if( id1 == 13 && !buffer_bool ) { pl1_internal=p1; buffer_bool=true; }
		if( id2 == 13 && buffer_bool ) pl3_internal = p2;
		if( id2 == 13 && !buffer_bool ) { pl1_internal=p2; buffer_bool=true; }
		if( id3 == 13 && buffer_bool ) pl3_internal = p3;
		if( id3 == 13 && !buffer_bool ) { pl1_internal=p3; buffer_bool=true; }
		if( id4 == 13 && buffer_bool ) pl3_internal = p4;
		
		buffer_bool = false;	//first antimuon has beed caught
		if( id1 == -13 && !buffer_bool ) { pl2_internal=p1; buffer_bool=true; }
		if( id2 == -13 && buffer_bool ) pl4_internal = p2;
		if( id2 == -13 && !buffer_bool ) { pl2_internal=p2; buffer_bool=true; }
		if( id3 == -13 && buffer_bool ) pl4_internal = p3;
		if( id3 == -13 && !buffer_bool ) { pl2_internal=p3; buffer_bool=true; }
		if( id4 == -13 && buffer_bool ) pl4_internal = p4;
		
		pA1_internal = NULL;
		Final_state = "4mu";
		
		return 0;
	}
	
	if( id_set[0] == -11 && id_set[1] == -11 && id_set[2] == 11 && id_set[3] == 11 && id_set[4] == 10000 )
	{
		buffer_bool = false;	//first electron has beed caught
		if( id1 == 11 && !buffer_bool ) { pl1_internal=p1; buffer_bool=true; }
		if( id2 == 11 && buffer_bool ) pl3_internal = p2;
		if( id2 == 11 && !buffer_bool ) { pl1_internal=p2; buffer_bool=true; }
		if( id3 == 11 && buffer_bool ) pl3_internal = p3;
		if( id3 == 11 && !buffer_bool ) { pl1_internal=p3; buffer_bool=true; }
		if( id4 == 11 && buffer_bool ) pl3_internal = p4;
		
		buffer_bool = false;	//first positron has beed caught
		if( id1 == -11 && !buffer_bool ) { pl2_internal=p1; buffer_bool=true; }
		if( id2 == -11 && buffer_bool ) pl4_internal = p2;
		if( id2 == -11 && !buffer_bool ) { pl2_internal=p2; buffer_bool=true; }
		if( id3 == -11 && buffer_bool ) pl4_internal = p3;
		if( id3 == -11 && !buffer_bool ) { pl2_internal=p3; buffer_bool=true; }
		if( id4 == -11 && buffer_bool ) pl4_internal = p4;
		
		pA1_internal = NULL;
		Final_state = "4e";
		
		return 0;
	}
	
	/// Four-lepton + photon final state block
	if( id_set[0] == -13 && id_set[1] == -11 && id_set[2] == 11 && id_set[3] == 13 && id_set[4] == 22 )
	{
		if( id1 == 11 ) pl1_internal = p1;
		if( id2 == 11 ) pl1_internal = p2;
		if( id3 == 11 ) pl1_internal = p3;
		if( id4 == 11 ) pl1_internal = p4;
		if( id5 == 11 ) pl1_internal = p5;
		
		if( id1 == -11 ) pl2_internal = p1;
		if( id2 == -11 ) pl2_internal = p2;
		if( id3 == -11 ) pl2_internal = p3;
		if( id4 == -11 ) pl2_internal = p4;
		if( id5 == -11 ) pl1_internal = p5;
		
		if( id1 == 13 ) pl3_internal = p1;
		if( id2 == 13 ) pl3_internal = p2;
		if( id3 == 13 ) pl3_internal = p3;
		if( id4 == 13 ) pl3_internal = p4;
		if( id5 == 13 ) pl3_internal = p5;
		
		if( id1 == -13 ) pl4_internal = p1;
		if( id2 == -13 ) pl4_internal = p2;
		if( id3 == -13 ) pl4_internal = p3;
		if( id4 == -13 ) pl4_internal = p4;
		if( id5 == -13 ) pl4_internal = p5;
		
		if( id1 == 22 ) pA1_internal = p1;
		if( id2 == 22 ) pA1_internal = p2;
		if( id3 == 22 ) pA1_internal = p3;
		if( id4 == 22 ) pA1_internal = p4;
		if( id5 == 22 ) pA1_internal = p5;
		Final_state = "2e2mA";
		
		return 0;
	}
	
	if( id_set[0] == -13 && id_set[1] == -13 && id_set[2] == 13 && id_set[3] == 13 && id_set[4] == 22 )
	{
		buffer_bool = false;	//first muon has beed caught
		if( id1 == 13 && !buffer_bool ) { pl1_internal=p1; buffer_bool=true; }
		if( id2 == 13 && buffer_bool ) pl3_internal = p2;
		if( id2 == 13 && !buffer_bool ) { pl1_internal=p2; buffer_bool=true; }
		if( id3 == 13 && buffer_bool ) pl3_internal = p3;
		if( id3 == 13 && !buffer_bool ) { pl1_internal=p3; buffer_bool=true; }
		if( id4 == 13 && buffer_bool ) pl3_internal = p4;
		if( id4 == 13 && !buffer_bool ) { pl1_internal=p4; buffer_bool=true; }
		if( id5 == 13 && buffer_bool ) pl3_internal = p5;
		
		buffer_bool = false;	//first antimuon has beed caught
		if( id1 == -13 && !buffer_bool ) { pl2_internal=p1; buffer_bool=true; }
		if( id2 == -13 && buffer_bool ) pl4_internal = p2;
		if( id2 == -13 && !buffer_bool ) { pl2_internal=p2; buffer_bool=true; }
		if( id3 == -13 && buffer_bool ) pl4_internal = p3;
		if( id3 == -13 && !buffer_bool ) { pl2_internal=p3; buffer_bool=true; }
		if( id4 == -13 && buffer_bool ) pl4_internal = p4;
		if( id4 == -13 && !buffer_bool ) { pl2_internal=p4; buffer_bool=true; }
		if( id5 == -13 && buffer_bool ) pl4_internal = p5;
		
		if( id1 == 22 ) pA1_internal = p1;
		if( id2 == 22 ) pA1_internal = p2;
		if( id3 == 22 ) pA1_internal = p3;
		if( id4 == 22 ) pA1_internal = p4;
		if( id5 == 22 ) pA1_internal = p5;
		Final_state = "4muA";
		
		return 0;
	}
	
	if( id_set[0] == -11 && id_set[1] == -11 && id_set[2] == 11 && id_set[3] == 11 && id_set[4] == 22 )
	{
		buffer_bool = false;	//first electron has beed caught
		if( id1 == 11 && !buffer_bool ) { pl1_internal=p1; buffer_bool=true; }
		if( id2 == 11 && buffer_bool ) pl3_internal = p2;
		if( id2 == 11 && !buffer_bool ) { pl1_internal=p2; buffer_bool=true; }
		if( id3 == 11 && buffer_bool ) pl3_internal = p3;
		if( id3 == 11 && !buffer_bool ) { pl1_internal=p3; buffer_bool=true; }
		if( id4 == 11 && buffer_bool ) pl3_internal = p4;
		if( id4 == 11 && !buffer_bool ) { pl1_internal=p4; buffer_bool=true; }
		if( id5 == 11 && buffer_bool ) pl3_internal = p5;
		
		buffer_bool = false;	//first positron has beed caught
		if( id1 == -11 && !buffer_bool ) { pl2_internal=p1; buffer_bool=true; }
		if( id2 == -11 && buffer_bool ) pl4_internal = p2;
		if( id2 == -11 && !buffer_bool ) { pl2_internal=p2; buffer_bool=true; }
		if( id3 == -11 && buffer_bool ) pl4_internal = p3;
		if( id3 == -11 && !buffer_bool ) { pl2_internal=p3; buffer_bool=true; }
		if( id4 == -11 && buffer_bool ) pl4_internal = p4;
		if( id4 == -11 && !buffer_bool ) { pl2_internal=p4; buffer_bool=true; }
		if( id5 == -11 && buffer_bool ) pl4_internal = p5;
		
		if( id1 == 22 ) pA1_internal = p1;
		if( id2 == 22 ) pA1_internal = p2;
		if( id3 == 22 ) pA1_internal = p3;
		if( id4 == 22 ) pA1_internal = p4;
		if( id5 == 22 ) pA1_internal = p5;
		Final_state = "4eA";
		
		return 0;
	}
	
	
	
	if( id_set[0] == 10000 && id_set[1] == 10000 && id_set[2] == 10000 && id_set[3] == 10000 && id_set[4] == 10000 )
	{
		if( Warning_Mode ) cout << "Warning. Particle ids are not set. Assuming a proper input-particle configuration.\n";
		if( Warning_Mode ) cout << "Proceeding according to a specified final state (" << Final_state  << ").\n";
		pl1_internal=p1; pl2_internal=p2; pl3_internal=p3; pl4_internal=p4; pA1_internal=p5;
		
		return 0;
	}
	
	return 1;
}


///////////////////////////////////
/// END OF MEKD_MG_Sorter.cpp   ///
///////////////////////////////////





///#include "MEKD_MG_2Model_Mixer.cpp"
/////////////////////////////////////////
/// INCLUDED MEKD_MG_2Model_Mixer.cpp ///
/// code follows below                ///
///                                   ///
/// Part responsible for              ///
/// holding mixed couplings           ///
/////////////////////////////////////////






///////////////////////////////////////
/// END OF MEKD_MG_2Model_Mixer.cpp ///
///////////////////////////////////////



#endif