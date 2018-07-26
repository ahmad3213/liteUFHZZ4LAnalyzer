/*
 * MEKD Test Suite 1
 * 
 * It is not the technology show-off and the class is NOT used in the optimal way
 * 
 */
#ifndef MEKD_TEST_H
#define MEKD_TEST_H


// C++ includes
#include <iostream>
#include <vector>
#include <cmath>

/// ME calculator
#include "../src/MEKD.cpp"

/// ROOT includes
#if (defined MEKD_STANDALONE && defined MEKD_with_ROOT) || !defined MEKD_STANDALONE
// ROOT includes
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"	// For the workaround for >> "Warning in <TClass::TClass>: no dictionary for class...
#include "TTree.h"
#include "TTreeIndex.h"
#endif




const bool Show_Description = true;
const bool Show_Basic_Data = true;
const bool Show_Advanced_Data = true;
const bool Show_Debug = true;

const unsigned int correlation_points = 500;
const unsigned int interference_gen_points = 200;	// number of fake events for averaging
const unsigned int interference_steps = 45;	// steps for the mixing angle from 0 to pi/2
const unsigned int scale_check_points = 100;	// must be smaller than correlation_points number
const unsigned int shuffles_for_ghosts = 100;	// number of permutations to search for remnant ghost effects

const double Precision_of_interest = 1E-12;
const double Match_between_ZZ_Full_and_Z4l_rel = 0.2;	// fraction of expected interference between Z propagator (s ch.) and the rest in the qq -> ZZ -> 4l.

const double m_e = 0.0005109989;
const double m_mu = 0.10565837;


int error_value, id1, id2, id3, id4;

double p1[4], p2[4], p3[4], p4[4];

TLorentzVector Lp1, Lp2, Lp3, Lp4;

vector<double*> Set_Of_Arrays;
vector<int> Set_Of_IDs;



#endif