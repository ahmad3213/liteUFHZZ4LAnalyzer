/*************************************************************************
*  Authors:   MEKD fans
*  More info: http://mekd.ihepa.ufl.edu
*  Contact:   mekd@phys.ufl.edu
*************************************************************************/
#ifndef MEKD_MEKD_h
#define MEKD_MEKD_h

// C++ includes
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#if (defined MEKD_STANDALONE && defined MEKD_with_ROOT) || !defined MEKD_STANDALONE
// ROOT includes
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TTreeIndex.h"
#endif

// MadGraph-based ME calculator
#include "MEKD_MG.h"


//////////////////////////////////////////////////////////////////////////
///  MEKD 2 interface class.
///
///  Provides neessary interface to the MadGraph-based ME calculator
///  and computes MEs and KDs for the process specified by the user.
///
//////////////////////////////////////////////////////////////////////////

using namespace std;


class MEKD {
public:
	///
	/// Constructor.
	///
	/// \param collisionEnergy							The sqrt(s) value in TeV (double, DEFAULT = 8).
	/// \param PDFName									The name of the parton density functions to be used (string, DEFAULT = "CTEQ6L", NONE = "").
	///
	MEKD(double collisionEnergy = 8, string PDFName = "CTEQ6L");
	
	
	///
	/// Compute and extract individual KD and MEs. Works after running int computeMEs(...).
	///
	/// Supported process names: "ZZ", "ggSpin0Pm", "ggSpin0M", "ggSpin0Ph", "qqSpin1P", "qqSpin1M", "ggSpin2Pm", "qqSpin2Pm".
	///
	/// \param[in]  processA, processB					The names of the processes X = A, B for which the KDs and MEs are computed (string, REQUIRED).
	/// \param[out] kd									The computed KD value for discrimination of processes A and B (double).
	/// \param[out] me2processA							The computed |ME|^2 for process A (double).
	/// \param[out] me2processB							The computed |ME|^2 for process B (double).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeKD( string processA, string processB,	// names of the processes
					double& kd,							// return KD
					double& me2processA,				// return |ME|^2 for process A
					double& me2processB );				// return |ME|^2 for process B
	
	
	///
	/// Compute KDs and MEs for process A and process B out of the 4-momenta of 4 leptons (lepton ordering does not matter).
	///
	/// Supported process names: "ZZ", "DY", "Custom", "ggSpin0Pm", "ggSpin0Ph", "ggSpin0M", "qqSpin1P", "qqSpin1M", "ggSpin2Pm", "ggSpin2Ph", "ggSpin2Mh", "ggSpin2Pb", "qqSpin2Pm", "qqSpin2Ph", "qqSpin2Mh", "qqSpin2Pb", "Spin0Pm", "Spin0Ph", "Spin0M" "Spin1P", "Spin1M", "Spin2Pm", "Spin2Ph", "Spin2Mh", "Spin2Pb", "ggSpin0", "qqSpin1", "ggSpin2", "qqSpin2", "Spin0", "Spin1", "Spin2", "qqZ4l_Background", "qqZ4l_Signal".
	///
	/// \param[in]  processA, processB					The names of the processes X = A, B for which the KDs and MEs are computed (string, REQUIRED).
	/// \param[in]  lept1P, lept2P, lept3P, lept4P		The input arrays with 4-momentum (E,px,py,pz) values of leptons N=1..4 (double*, REQUIRED).
	/// \param[in]  lept1Id, lept2Id, lept3Id, lept4Id	The input IDs (PDG) of leptons N=1..4 (int, REQUIRED).
	/// \param[out] kd									The computed KD value for discrimination of processes A and B (double).
	/// \param[out] me2processA							The computed |ME|^2 for process A (double).
	/// \param[out] me2processB							The computed |ME|^2 for process B (double).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeKD( string processA, string processB,	// names of the processes
					double lept1P[], int lept1Id,		// 4-momentum (E,px,py,pz) and id (PDG) of lepton 1
					double lept2P[], int lept2Id,		// 4-momentum (E,px,py,pz) and id (PDG) of lepton 2
					double lept3P[], int lept3Id,		// 4-momentum (E,px,py,pz) and id (PDG) of lepton 3
					double lept4P[], int lept4Id,		// 4-momentum (E,px,py,pz) and id (PDG) of lepton 4
					double& kd,							// return KD
					double& me2processA,				// return |ME|^2 for process A
					double& me2processB );				// return |ME|^2 for process B
	
	
	///
	/// Compute KDs and MEs for process A and process B out of the 4-momenta of the input particles (ordering does not matter).
	///
	/// Supported process names: "ZZ", "DY", "Custom", "ggSpin0Pm", "ggSpin0Ph", "ggSpin0M", "qqSpin1P", "qqSpin1M", "ggSpin2Pm", "ggSpin2Ph", "ggSpin2Mh", "ggSpin2Pb", "qqSpin2Pm", "qqSpin2Ph", "qqSpin2Mh", "qqSpin2Pb", "Spin0Pm", "Spin0Ph", "Spin0M" "Spin1P", "Spin1M", "Spin2Pm", "Spin2Ph", "Spin2Mh", "Spin2Pb", "ggSpin0", "qqSpin1", "ggSpin2", "qqSpin2", "Spin0", "Spin1", "Spin2", "qqZ4l_Background", "qqZ4l_Signal".
	///
	/// \param[in]  processA, processB					The names of the processes X = A, B for which the KDs and MEs are computed (string, REQUIRED).
	/// \param[in]  input_Ps							The input vector of arrays with 4-momentum (E,px,py,pz) values of particles N=1..5 (vector<double*>, REQUIRED).
	/// \param[in]  input_IDs							The input vector of IDs (PDG) of particles N=1..5 (vector<int>, REQUIRED).
	/// \param[out] kd									The computed KD value for discrimination of processes A and B (double).
	/// \param[out] me2processA							The computed |ME|^2 for process A (double).
	/// \param[out] me2processB							The computed |ME|^2 for process B (double).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeKD( string processA, string processB,					// names of the processes
					vector<double*> input_Ps, vector<int> input_IDs,	// 4-momenta (E,px,py,pz) and IDs (PDG) of input particles
					double& kd,											// return KD
					double& me2processA,								// return |ME|^2 for process A
					double& me2processB );								// return |ME|^2 for process B
	
	
	///
	/// Compute ME for a processName out of the 4-momenta of the input particles (ordering does not matter).
	///
	/// Supported process names: "ZZ", "DY", "Custom", "ggSpin0Pm", "ggSpin0Ph", "ggSpin0M", "qqSpin1P", "qqSpin1M", "ggSpin2Pm", "ggSpin2Ph", "ggSpin2Mh", "ggSpin2Pb", "qqSpin2Pm", "qqSpin2Ph", "qqSpin2Mh", "qqSpin2Pb", "Spin0Pm", "Spin0Ph", "Spin0M" "Spin1P", "Spin1M", "Spin2Pm", "Spin2Ph", "Spin2Mh", "Spin2Pb", "ggSpin0", "qqSpin1", "ggSpin2", "qqSpin2", "Spin0", "Spin1", "Spin2", "qqZ4l_Background", "qqZ4l_Signal".
	///
	/// \param[in]  processName							The Name of the process for which the ME is to be computed (string, REQUIRED).
	/// \param[in]  input_Ps							The input vector of arrays with 4-momentum (E,px,py,pz) values of particles N=1..5 (vector<double*>, REQUIRED).
	/// \param[in]  input_IDs							The input vector of IDs (PDG) of particles N=1..5 (vector<int>, REQUIRED).
	/// \param[out] me2process							The computed |ME|^2 for process of interest (double).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeME( string processName,
					vector<double*> input_Ps, vector<int> input_IDs,
					double& me2process );
	
	
	///
	/// Compute all important/default MEs first for the use with computeKD(string, string, double&, double&, double&). Option Custom is outside this scope.
	///
	/// \param[in]  lept1P, lept2P, lept3P, lept4P		The input arrays with 4-momentum (E,px,py,pz) values of leptons N=1..4 (double*, REQUIRED).
	/// \param[in]  lept1Id, lept2Id, lept3Id, lept4Id	The input IDs (PDG) of leptons N=1..4 (int, REQUIRED).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeMEs( double lept1P[], int lept1Id,		// 4-momentum (E,px,py,pz) and id (PDG) of lepton 1
					double lept2P[], int lept2Id,		// 4-momentum (E,px,py,pz) and id (PDG) of lepton 2
					double lept3P[], int lept3Id,		// 4-momentum (E,px,py,pz) and id (PDG) of lepton 3
					double lept4P[], int lept4Id );		// 4-momentum (E,px,py,pz) and id (PDG) of lepton 4
	
	
	///
	/// Compute all important/default MEs first for the use with computeKD(string, string, double&, double&, double&). Option Custom is outside this scope.
	///
	/// \param[in]  input_Ps							The input vector of arrays with 4-momentum (E,px,py,pz) values of particles N=1..5 (vector<double*>, REQUIRED).
	/// \param[in]  input_IDs							The input vector of IDs (PDG) of particles N=1..5 (vector<int>, REQUIRED).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeMEs( vector<double*> input_Ps, vector<int> input_IDs );	// 4-momenta (E,px,py,pz) and IDs (PDG) of input particles
	
	
	///
	/// Mixed-state ME mixer of (gg)Spin0Pm, (gg)Spin0Ph, (gg)Spin0Phexotic, and (gg)Spin0M states, corresponding couplings kappa_1/rhomu01, kappa_2/rhomu02, kappa_3/rhomu03, kappa_4/rhomu04.
	///
	/// \param[in]	Spin0Pm_relamp						The relative complex amplitude of the Spin0Pm state
	/// \param[in]	Spin0Ph_relamp						The relative complex amplitude of the Spin0Ph state
	/// \param[in]	Spin0Phexo_relamp					The relative complex amplitude of the Spin0PhExotic state
	/// \param[in]	Spin0M_relamp						The relative complex amplitude of the Spin0M state
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int Mix_Spin0( complex<double> Spin0Pm_relamp, complex<double> Spin0Ph_relamp, complex<double> Spin0Phexo_relamp, complex<double> Spin0M_relamp );
	
	
	///
	/// Mixed-state ME mixer of production like qqSpin1M, qqSpin1P, exotic, exotic pseudo, and decay like Spin1M, Spin1P states, corresponding couplings rhoQ11, rhoQ12, rhoQ13, rhoQ14, b1z/rhomu11, b2z/rhomu12, rhomu13, rhomu14.
	///
	/// \param[in]	prod_Spin1M_relamp					The relative complex amplitude for the Spin1M-like production
	/// \param[in]	prod_Spin1P_relamp					The relative complex amplitude for the Spin1P-like production
	/// \param[in]	prod_Spin1Mexo_relamp				The relative complex amplitude for the Spin1Mexotic-like production
	/// \param[in]	prod_Spin1Pexo_relamp				The relative complex amplitude for the Spin1Pexotic-like production
	/// \param[in]	dec_Spin1M_relamp					The relative complex amplitude for the Spin1M-like decay
	/// \param[in]	dec_Spin1P_relamp					The relative complex amplitude for the Spin1P-like decay
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int Mix_Spin1( complex<double> prod_Spin1M_relamp, complex<double> prod_Spin1P_relamp, complex<double> prod_Spin1Mexo_relamp, complex<double> prod_Spin1Pexo_relamp, complex<double> dec_Spin1M_relamp, complex<double> dec_Spin1P_relamp, complex<double> dec_Spin1_rhomu13_relamp, complex<double> dec_Spin1_rhomu14_relamp );
	
	
	///
	/// Mixed-state ME mixer of Spin-2 states. Production couplings: k1g/rhoQ21, ..., k4g/rhoQ24, ..., k10g, decay couplings: k1z/rhomu21, k4z/rhomu24, ..., k10z.
	///
	/// \param[in]	*prod_Spin2_relamp					The relative complex amplitudes for the Spin-2 state production. An array of size 10
	/// \param[in]	*dec_Spin2_relamp					The relative complex amplitudes for the Spin-2 state decay. An array of size 10
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int Mix_Spin2( complex<double> *prod_Spin2_relamp, complex<double> *dec_Spin2_relamp );
	
	
	/// Mixed-state ME mixer. Variables
	complex<double> *m_Mixing_Coefficients_Spin0;
	complex<double> *m_Mixing_Coefficients_Spin1;
	complex<double> *m_Mixing_Coefficients_Spin2;
	
	
	/// Madgraph-based ME calculator
	MEKD_MG MEKD_MG_Calc;
	
private:
	/// Properties. Variables.
	unsigned int buffer_uint;		// For counters as may be constantly created and destroyed
	int buffer_int;					// For internal collection of return values
	double m_collisionEnergy;		// c.m. collision energy sqrt(s) in TeV
	double ME_ZZ, ME_Spin0PSMH, ME_Spin0Ph, ME_Spin0M;	// computeMEs(...) results
	double ME_Spin1P, ME_Spin1M, ME_ggSpin2Pm, ME_qqSpin2Pm;	// computeMEs(...) results
	string m_PDFName;				// Name of the parton density functions to be used. Supported: CTEQ6l;
	string m_process;				// Name of the process (background, signal hypotheses, etc.). Supported: Custom, SMHiggs, CPoddScalar, CPevenScalar, Spin2particle, ZZ
	string m_processA;				// Name of the process A (background, signal hypotheses, etc.). Supported: Custom, SMHiggs, CPoddScalar, CPevenScalar, Spin2particle, ZZ
	string m_processB;				// Name of the process B (background, signal hypotheses, etc.). Supported: Custom, SMHiggs, CPoddScalar, CPevenScalar, Spin2particle, ZZ
	bool m_usePDF;					// flag to use PDFs (true) or not (false)
	bool m_runBackgroundME;			// flat to run the ME for ZZ process (true) or not (false)
	enum ERRCodes	{NO_ERR, ERR_SQRT, ERR_PDFS, ERR_PROCESS, NUM_ERRORS, ERR_INPUT, ERR_OTHER};
	
	vector<int> four_particle_IDs_i;		// For the storage of the four IDs
	vector<double*> four_particle_Ps_i;		// For the storage of the four four-momenta
	
	/// Methods
	int setProcessName(string process);		// sanity check for input process name, translation to the the names supported by MEKD_MG
	int setProcessNames(string processA, string processB);	// sanity check for input process names, translation to the the names supported by MEKD_MG
	int processParameters();				// sanity check for internal paramters
	
#if (defined MEKD_STANDALONE && defined MEKD_with_ROOT) || !defined MEKD_STANDALONE
//------------------------------------------------------------------------
// ROOT-compatabile members
//------------------------------------------------------------------------
public:
	///
	/// Compute KDs and MEs for process A and process B. Works after running int computeMEs(...).
	/// The overloaded method that supports input parameters of ROOT types TString and TLorentzVector.
	///
	/// Supported process names: "ZZ", "ggSpin0Pm", "ggSpin0M", "ggSpin0Ph", "qqSpin1P", "qqSpin1M", "ggSpin2Pm", "qqSpin2Pm".
	///
	/// \param[in]  processA, processB					The names of the processes X = A, B for which the KDs and MEs are computed (TString, REQUIRED).
	/// \param[out] kd									The computed KD value for discrimination of processes A and B (double).
	/// \param[out] me2processA							The computed |ME|^2 for process A (double).
	/// \param[out] me2processB							The computed |ME|^2 for process B (double).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeKD( TString processA, TString processB,	// names of the processes
					double& kd,							// return KD
					double& me2processA,				// return |ME|^2 for process A
					double& me2processB );				// return |ME|^2 for process B
	
	
	///
	/// Compute KDs and MEs for process A and process B out of the 4-momenta of 4 leptons (lepton ordering does not matter).
	/// The overloaded method that supports input parameters of ROOT types TString and TLorentzVector.
	///
	/// Supported process names: "ZZ", "DY", "Custom", "ggSpin0Pm", "ggSpin0Ph", "ggSpin0M", "qqSpin1P", "qqSpin1M", "ggSpin2Pm", "ggSpin2Ph", "ggSpin2Mh", "ggSpin2Pb", "qqSpin2Pm", "qqSpin2Ph", "qqSpin2Mh", "qqSpin2Pb", "Spin0Pm", "Spin0Ph", "Spin0M" "Spin1P", "Spin1M", "Spin2Pm", "Spin2Ph", "Spin2Mh", "Spin2Pb", "ggSpin0", "qqSpin1", "ggSpin2", "qqSpin2", "Spin0", "Spin1", "Spin2", "qqZ4l_Background", "qqZ4l_Signal".
	///
	/// \param[in]  processA, processB					The names of the processes X = A, B for which the KDs and MEs are computed (TString, REQUIRED).
	/// \param[in]  lept1P, lept2P, lept3P, lept4P		The input arrays with 4-momentum (E,px,py,pz) values of leptons N=1..4 (TLorentzVector, REQUIRED).
	/// \param[in]  lept1Id, lept2Id, lept3Id, lept4Id	The input IDs (PDG) of leptons N=1..4 (int, REQUIRED).
	/// \param[out] kd									The computed KD value for discrimination of processes A and B (double).
	/// \param[out] me2processA							The computed |ME|^2 for process A (double).
	/// \param[out] me2processB							The computed |ME|^2 for process B (double).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeKD( TString processA, TString processB,	// names of the processes
					TLorentzVector lept1P, int lept1Id,	// 4-momentum and id (PDG) of lepton 1
					TLorentzVector lept2P, int lept2Id,	// 4-momentum and id (PDG) of lepton 2
					TLorentzVector lept3P, int lept3Id,	// 4-momentum and id (PDG) of lepton 3
					TLorentzVector lept4P, int lept4Id,	// 4-momentum and id (PDG) of lepton 4
					double& kd,							// return KD
					double& me2processA,				// return |ME|^2 for process A
					double& me2processB );				// return |ME|^2 for process B
	
	
	///
	/// Compute KDs and MEs for process A and process B out of the 4-momenta of the input particles (ordering does not matter).
	/// The overloaded method that supports input parameters of ROOT types TString and TLorentzVector.
	///
	/// Supported process names: "ZZ", "DY", "Custom", "ggSpin0Pm", "ggSpin0Ph", "ggSpin0M", "qqSpin1P", "qqSpin1M", "ggSpin2Pm", "ggSpin2Ph", "ggSpin2Mh", "ggSpin2Pb", "qqSpin2Pm", "qqSpin2Ph", "qqSpin2Mh", "qqSpin2Pb", "Spin0Pm", "Spin0Ph", "Spin0M" "Spin1P", "Spin1M", "Spin2Pm", "Spin2Ph", "Spin2Mh", "Spin2Pb", "ggSpin0", "qqSpin1", "ggSpin2", "qqSpin2", "Spin0", "Spin1", "Spin2", "qqZ4l_Background", "qqZ4l_Signal".
	///
	/// \param[in]  processA, processB					The names of the processes X = A, B for which the KDs and MEs are computed (Tstring, REQUIRED).
	/// \param[in]  input_Ps							The input vector of TLorentzVectors with 4-momentum (E,px,py,pz) values of particles N=1..5 (vector<TLorentzVector>, REQUIRED).
	/// \param[in]  input_IDs							The input vector of IDs (PDG) of particles N=1..5 (int, REQUIRED).
	/// \param[out] kd									The computed KD value for discrimination of processes A and B (double).
	/// \param[out] me2processA							The computed |ME|^2 for process A (double).
	/// \param[out] me2processB							The computed |ME|^2 for process B (double).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeKD( TString processA, TString processB,						// names of the processes
					vector<TLorentzVector> input_Ps, vector<int> input_IDs,	// 4-momenta (E,px,py,pz) and IDs (PDG) of input particles
					double& kd,												// return KD
					double& me2processA,									// return |ME|^2 for process A
					double& me2processB );									// return |ME|^2 for process B
	
	
	///
	/// Compute ME for a processName out of the 4-momenta of the input particles (ordering does not matter).
	///
	/// Supported process names: "ZZ", "DY", "Custom", "ggSpin0Pm", "ggSpin0Ph", "ggSpin0M", "qqSpin1P", "qqSpin1M", "ggSpin2Pm", "ggSpin2Ph", "ggSpin2Mh", "ggSpin2Pb", "qqSpin2Pm", "qqSpin2Ph", "qqSpin2Mh", "qqSpin2Pb", "Spin0Pm", "Spin0Ph", "Spin0M" "Spin1P", "Spin1M", "Spin2Pm", "Spin2Ph", "Spin2Mh", "Spin2Pb", "ggSpin0", "qqSpin1", "ggSpin2", "qqSpin2", "Spin0", "Spin1", "Spin2", "qqZ4l_Background", "qqZ4l_Signal".
	///
	/// \param[in]  processName							The name of the process for which the ME is to be computed (TString, REQUIRED).
	/// \param[in]  input_Ps							The input vector of TLorentzVectors with 4-momentum (E,px,py,pz) values of particles N=1..5 (vector<TLorentzVector>, REQUIRED).
	/// \param[in]  input_IDs							The input vector of IDs (PDG) of particles N=1..5 (vector<int>, REQUIRED).
	/// \param[out] me2process							The computed |ME|^2 for process of interest (double).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeME( TString processName,
					vector<TLorentzVector> input_Ps, vector<int> input_IDs,
					double& me2process );
	
	
	///
	/// Compute all important/default MEs first for the use with computeKD(TString, TString, double&, double&, double&). Option Custom is outside this scope.
	///
	/// \param[in]  lept1P, lept2P, lept3P, lept4P		The input arrays with 4-momentum (E,px,py,pz) values of leptons N=1..4 (TLorentzVector, REQUIRED).
	/// \param[in]  lept1Id, lept2Id, lept3Id, lept4Id	The input IDs (PDG) of leptons N=1..4 (int, REQUIRED).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeMEs( TLorentzVector lept1P, int lept1Id,		// 4-momentum and id (PDG) of lepton 1
					TLorentzVector lept2P, int lept2Id,		// 4-momentum and id (PDG) of lepton 2
					TLorentzVector lept3P, int lept3Id,		// 4-momentum and id (PDG) of lepton 3
					TLorentzVector lept4P, int lept4Id );	// 4-momentum and id (PDG) of lepton 4
	
	
	///
	/// Compute all important/default MEs first for the use with computeKD(string/TString, string/TString, double&, double&, double&). Option Custom is outside this scope.
	///
	/// \param[in]  input_Ps							The input vector of TLorentzVectors with 4-momentum (E,px,py,pz) values of particles N=1..5 (vector<TLorentzVector>, REQUIRED).
	/// \param[in]  input_IDs							The input vector of IDs (PDG) of particles N=1..5 (vector<int>, REQUIRED).
	/// \return											The error code of the computation: 0 = NO_ERR, 1 = ERR_SQRTS, 2 = ERR_PDFS, 3 = ERR_PROCESS, 4 = NUM_ERRORS, 5 = ERR_INPUT, 6 = ERR_OTHER
	///
	int computeMEs( vector<TLorentzVector> input_Ps, vector<int> input_IDs );	// 4-momenta (E,px,py,pz) and IDs (PDG) of input particles
	
	
private:
	double lept1P_i[4], lept2P_i[4], lept3P_i[4], lept4P_i[4];	// For storing TLorentzVectors for internal use
	vector<double*> input_Ps_i;									// For storing vector of TLorentzVectors for internal use
	
#endif
};

//////////////////////////////////////////////////////////////////////////

#endif

