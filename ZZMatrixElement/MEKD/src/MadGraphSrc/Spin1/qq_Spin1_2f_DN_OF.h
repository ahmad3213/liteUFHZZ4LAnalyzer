//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.2, 2014-02-07
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_HEF_MEKD2_1_ssx_zp_emepmummup_H
#define MG5_Sigma_HEF_MEKD2_1_ssx_zp_emepmummup_H

#include <complex> 
#include <vector> 

#include "../Parameters_MEKD.h"	// Changed by Convert_source 0.2
#include "../read_slha.h"	// Added by Convert_source 0.2

using namespace std;

//==========================================================================
// A class for calculating the matrix elements for
// Process: s s~ > zp > e- e+ mu- mu+ S1VV=1 QED=2 S1QQ=2 / h xg
//--------------------------------------------------------------------------

class qq_Spin1_2f_DN_OF
{
  public:

    // Constructor.
    qq_Spin1_2f_DN_OF() {}

    // Initialize process.
	virtual void initProc(string param_card_name);
	
	// Update process.
	virtual void updateProc(SLHAReader_MEKD &slha);

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin();

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat();

    // Info on the subprocess.
    virtual string name() const {return "s s~ > e- e+ mu- mu+ (HEF_MEKD2_1)";}

    virtual int code() const {return 0;}

    const vector<double> & getMasses() const {return mME;}

    // Get and set momenta for matrix element evaluation
    vector < double * > getMomenta(){return p;}
    void setMomenta(vector < double * > & momenta){p = momenta;}
    void setInitial(int inid1, int inid2){id1 = inid1;id2 = inid2;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
    static const int ninitial = 2;
    static const int nexternal = 6;
    static const int nprocesses = 2;

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]);
    static const int nwavefuncs = 33;
    std::complex<double> w[nwavefuncs][18];
    static const int namplitudes = 27;
    std::complex<double> amp[namplitudes];
	int ntry, sum_hel, ngood;	// Moved here by Convert_source 0.2 
    double matrix_ssx_zp_emepmummup_no_hxg();

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses];

    // Color flows, used when selecting color
    double * jamp2[nprocesses];

    // Pointer to the model parameters
	Parameters_MEKD * pars;	// Changed by Convert_source 0.2

    // vector with external particle masses
    vector<double> mME;

    // vector with momenta (to be changed each event)
    vector < double * > p;
    // Initial particle ids
    int id1, id2;

};


#endif  // MG5_Sigma_HEF_MEKD2_1_ssx_zp_emepmummup_H
