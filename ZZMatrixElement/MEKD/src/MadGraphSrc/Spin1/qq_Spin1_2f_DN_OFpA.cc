//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.2, 2014-02-07
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "qq_Spin1_2f_DN_OFpA.h"
#include "../HelAmps_HEF_MEKD2_1.h"	// Changed by Convert_source 0.2

using namespace MG5_HEF_MEKD2_1;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: s s~ > zp > e- e+ mu- mu+ a S1VV=1 QED=3 S1QQ=2 / h xg

//--------------------------------------------------------------------------
// Initialize process.

void qq_Spin1_2f_DN_OFpA::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MEKD::getInstance();	// Changed by Convert_source 0.2 
  SLHAReader_MEKD slha(param_card_name);
  pars->setIndependentParameters(slha);
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// Moved here by Convert_source 0.2
  // Set external particle masses for this matrix element
  mME.push_back(pars->MS);
  mME.push_back(pars->MS);
  mME.push_back(pars->Me);
  mME.push_back(pars->Me);
  mME.push_back(pars->MM);
  mME.push_back(pars->MM);
  mME.push_back(pars->ZERO);
  jamp2[0] = new double[1];
	for( int count=0; count<namplitudes; count++ ) amp[count] = 0;
}

//--------------------------------------------------------------------------
// Update process.	// Created here by Convert_source 0.2

void qq_Spin1_2f_DN_OFpA::updateProc(SLHAReader_MEKD &slha)
{
	pars->setIndependentParameters(slha);
	pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// needed when altering couplings
	
	// Set external particle masses for this matrix element
	// Should correspond to initProc
	mME[0]=(pars->MS);
	mME[1]=(pars->MS);
	mME[2]=(pars->Me);
	mME[3]=(pars->Me);
	mME[4]=(pars->MM);
	mME[5]=(pars->MM);
	mME[6]=(pars->ZERO);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void qq_Spin1_2f_DN_OFpA::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters();
  pars->setDependentCouplings();
	// Deleted by Convert_source 0.2
	
  // Reset color flows
  for(int i = 0;i < 1;i++ )
    jamp2[0][i] = 0.;

  // Local variables and constants
  const int ncomb = 128;
  static bool goodhel[ncomb] = {ncomb * false};
//	static int ntry = 0, sum_hel = 0, ngood = 0;	// Moved by Convert_source 0.2
  static int igood[ncomb];
  static int jhel;
//	std::complex<double> * * wfs;	// Changed by Convert_source 0.2
  double t[nprocesses];
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1,
      -1}, {-1, -1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, -1, 1, -1}, {-1, -1,
      -1, -1, -1, 1, 1}, {-1, -1, -1, -1, 1, -1, -1}, {-1, -1, -1, -1, 1, -1,
      1}, {-1, -1, -1, -1, 1, 1, -1}, {-1, -1, -1, -1, 1, 1, 1}, {-1, -1, -1,
      1, -1, -1, -1}, {-1, -1, -1, 1, -1, -1, 1}, {-1, -1, -1, 1, -1, 1, -1},
      {-1, -1, -1, 1, -1, 1, 1}, {-1, -1, -1, 1, 1, -1, -1}, {-1, -1, -1, 1, 1,
      -1, 1}, {-1, -1, -1, 1, 1, 1, -1}, {-1, -1, -1, 1, 1, 1, 1}, {-1, -1, 1,
      -1, -1, -1, -1}, {-1, -1, 1, -1, -1, -1, 1}, {-1, -1, 1, -1, -1, 1, -1},
      {-1, -1, 1, -1, -1, 1, 1}, {-1, -1, 1, -1, 1, -1, -1}, {-1, -1, 1, -1, 1,
      -1, 1}, {-1, -1, 1, -1, 1, 1, -1}, {-1, -1, 1, -1, 1, 1, 1}, {-1, -1, 1,
      1, -1, -1, -1}, {-1, -1, 1, 1, -1, -1, 1}, {-1, -1, 1, 1, -1, 1, -1},
      {-1, -1, 1, 1, -1, 1, 1}, {-1, -1, 1, 1, 1, -1, -1}, {-1, -1, 1, 1, 1,
      -1, 1}, {-1, -1, 1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1, 1}, {-1, 1, -1,
      -1, -1, -1, -1}, {-1, 1, -1, -1, -1, -1, 1}, {-1, 1, -1, -1, -1, 1, -1},
      {-1, 1, -1, -1, -1, 1, 1}, {-1, 1, -1, -1, 1, -1, -1}, {-1, 1, -1, -1, 1,
      -1, 1}, {-1, 1, -1, -1, 1, 1, -1}, {-1, 1, -1, -1, 1, 1, 1}, {-1, 1, -1,
      1, -1, -1, -1}, {-1, 1, -1, 1, -1, -1, 1}, {-1, 1, -1, 1, -1, 1, -1},
      {-1, 1, -1, 1, -1, 1, 1}, {-1, 1, -1, 1, 1, -1, -1}, {-1, 1, -1, 1, 1,
      -1, 1}, {-1, 1, -1, 1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1, 1}, {-1, 1, 1, -1,
      -1, -1, -1}, {-1, 1, 1, -1, -1, -1, 1}, {-1, 1, 1, -1, -1, 1, -1}, {-1,
      1, 1, -1, -1, 1, 1}, {-1, 1, 1, -1, 1, -1, -1}, {-1, 1, 1, -1, 1, -1, 1},
      {-1, 1, 1, -1, 1, 1, -1}, {-1, 1, 1, -1, 1, 1, 1}, {-1, 1, 1, 1, -1, -1,
      -1}, {-1, 1, 1, 1, -1, -1, 1}, {-1, 1, 1, 1, -1, 1, -1}, {-1, 1, 1, 1,
      -1, 1, 1}, {-1, 1, 1, 1, 1, -1, -1}, {-1, 1, 1, 1, 1, -1, 1}, {-1, 1, 1,
      1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1, 1}, {1, -1, -1, -1, -1, -1, -1}, {1,
      -1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, -1, 1, -1}, {1, -1, -1, -1, -1,
      1, 1}, {1, -1, -1, -1, 1, -1, -1}, {1, -1, -1, -1, 1, -1, 1}, {1, -1, -1,
      -1, 1, 1, -1}, {1, -1, -1, -1, 1, 1, 1}, {1, -1, -1, 1, -1, -1, -1}, {1,
      -1, -1, 1, -1, -1, 1}, {1, -1, -1, 1, -1, 1, -1}, {1, -1, -1, 1, -1, 1,
      1}, {1, -1, -1, 1, 1, -1, -1}, {1, -1, -1, 1, 1, -1, 1}, {1, -1, -1, 1,
      1, 1, -1}, {1, -1, -1, 1, 1, 1, 1}, {1, -1, 1, -1, -1, -1, -1}, {1, -1,
      1, -1, -1, -1, 1}, {1, -1, 1, -1, -1, 1, -1}, {1, -1, 1, -1, -1, 1, 1},
      {1, -1, 1, -1, 1, -1, -1}, {1, -1, 1, -1, 1, -1, 1}, {1, -1, 1, -1, 1, 1,
      -1}, {1, -1, 1, -1, 1, 1, 1}, {1, -1, 1, 1, -1, -1, -1}, {1, -1, 1, 1,
      -1, -1, 1}, {1, -1, 1, 1, -1, 1, -1}, {1, -1, 1, 1, -1, 1, 1}, {1, -1, 1,
      1, 1, -1, -1}, {1, -1, 1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, 1, -1}, {1, -1,
      1, 1, 1, 1, 1}, {1, 1, -1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, -1, 1},
      {1, 1, -1, -1, -1, 1, -1}, {1, 1, -1, -1, -1, 1, 1}, {1, 1, -1, -1, 1,
      -1, -1}, {1, 1, -1, -1, 1, -1, 1}, {1, 1, -1, -1, 1, 1, -1}, {1, 1, -1,
      -1, 1, 1, 1}, {1, 1, -1, 1, -1, -1, -1}, {1, 1, -1, 1, -1, -1, 1}, {1, 1,
      -1, 1, -1, 1, -1}, {1, 1, -1, 1, -1, 1, 1}, {1, 1, -1, 1, 1, -1, -1}, {1,
      1, -1, 1, 1, -1, 1}, {1, 1, -1, 1, 1, 1, -1}, {1, 1, -1, 1, 1, 1, 1}, {1,
      1, 1, -1, -1, -1, -1}, {1, 1, 1, -1, -1, -1, 1}, {1, 1, 1, -1, -1, 1,
      -1}, {1, 1, 1, -1, -1, 1, 1}, {1, 1, 1, -1, 1, -1, -1}, {1, 1, 1, -1, 1,
      -1, 1}, {1, 1, 1, -1, 1, 1, -1}, {1, 1, 1, -1, 1, 1, 1}, {1, 1, 1, 1, -1,
      -1, -1}, {1, 1, 1, 1, -1, -1, 1}, {1, 1, 1, 1, -1, 1, -1}, {1, 1, 1, 1,
      -1, 1, 1}, {1, 1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, 1, -1, 1}, {1, 1, 1, 1,
      1, 1, -1}, {1, 1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {36, 36};

  ntry = ntry + 1;

  // Reset the matrix elements
  for(int i = 0;i < nprocesses;i++ )
  {
    matrix_element[i] = 0.;
  }
  // Define permutation
  int perm[nexternal];
  for(int i = 0;i < nexternal;i++ )
  {
    perm[i] = i;
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0;ihel < ncomb;ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]);
        t[0] = matrix_ssx_zp_emepmummupa_no_hxg();
        // Mirror initial state momenta for mirror process
        perm[0] = 1;
        perm[1] = 0;
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]);
        // Mirror back
        perm[0] = 0;
        perm[1] = 1;
        // Calculate matrix elements
        t[1] = matrix_ssx_zp_emepmummupa_no_hxg();
        double tsum = 0;
        for(int iproc = 0;iproc < nprocesses;iproc++ )
        {
          matrix_element[iproc] += t[iproc];
          tsum += t[iproc];
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true;
          ngood++;
          igood[ngood] = ihel;
        }
      }
    }
    jhel = 0;
    sum_hel = min(sum_hel, ngood);
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0;j < sum_hel;j++ )
    {
      jhel++;
      if (jhel >= ngood)
        jhel = 0;
      double hwgt = double(ngood)/double(sum_hel);
      int ihel = igood[jhel];
      calculate_wavefunctions(perm, helicities[ihel]);
      t[0] = matrix_ssx_zp_emepmummupa_no_hxg();
      // Mirror initial state momenta for mirror process
      perm[0] = 1;
      perm[1] = 0;
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]);
      // Mirror back
      perm[0] = 0;
      perm[1] = 1;
      // Calculate matrix elements
      t[1] = matrix_ssx_zp_emepmummupa_no_hxg();
      for(int iproc = 0;iproc < nprocesses;iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt;
      }
    }
  }

  for (int i = 0;i < nprocesses;i++ )
    matrix_element[i] /= denominators[i];



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double qq_Spin1_2f_DN_OFpA::sigmaHat() 
{
  // Select between the different processes
  if(id1 == -3 && id2 == 3)
  {
    // Add matrix elements for processes with beams (-3, 3)
    return matrix_element[1];
  }
  else if(id1 == 3 && id2 == -3)
  {
    // Add matrix elements for processes with beams (3, -3)
    return matrix_element[0];
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.;
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void qq_Spin1_2f_DN_OFpA::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
//	int i, j;	// Changed by Convert_source 0.2

  // Calculate all wavefunctions
  ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]);
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]);
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]);
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]);
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]);
  ixxxxx(p[perm[5]], mME[5], hel[5], -1, w[5]);
  vxxxxx(p[perm[6]], mME[6], hel[6], +1, w[6]);
  FFV1_2_3_4_3(w[0], w[1], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
      pars->HEF_MEKD2_1_GC_150, pars->MZp, pars->WZp, w[7]);
  FFV2P0_3(w[3], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[8]);
  FFV1_2_3_4_1(w[4], w[7], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, pars->MM, pars->ZERO, w[9]);
  FFV2_2(w[5], w[8], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[10]);
  FFV1_2_3_4_2(w[5], w[7], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, pars->MM, pars->ZERO, w[11]);
  FFV2_1(w[4], w[8], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[12]);
  FFV5_7_3(w[3], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[13]);
  FFV5_7_2(w[5], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[14]);
  FFV5_7_1(w[4], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[15]);
//   FFV2P0_3(w[0], w[1], pars->HEF_MEKD2_1_GC_3, pars->ZERO, pars->ZERO, w[16]);
//   FFV1_2_3_4_3(w[3], w[2], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, pars->MZp, pars->WZp, w[17]);
//   FFV2_1(w[4], w[16], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[18]);
//   FFV1_2_3_4_2(w[5], w[17], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, pars->MM, pars->ZERO, w[19]);
//   FFV2_2(w[5], w[16], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[20]);
//   FFV1_2_3_4_1(w[4], w[17], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, pars->MM, pars->ZERO, w[21]);
//   FFV5_6_3(w[0], w[1], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, pars->MZ, pars->WZ, w[22]);
//   FFV5_7_1(w[4], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
//       w[23]);
//   FFV5_7_2(w[5], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
//       w[24]);
  FFV2_1(w[4], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[25]);
  FFV1_2_3_4_1(w[25], w[7], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, pars->MM, pars->ZERO, w[26]);
//   VVV1_2_1(w[13], w[7], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, pars->MZ, pars->WZ, w[27]);
//   FFV2_1(w[25], w[16], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[28]);
//   VVV1_2_1(w[22], w[17], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, pars->MZ, pars->WZ, w[29]);
//   FFV5_7_1(w[25], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
//       w[30]);
//   VVV1_2_3(w[22], w[13], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, pars->MZp, pars->WZp, w[31]);
  FFV2_2(w[5], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[32]);
  FFV1_2_3_4_2(w[32], w[7], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, pars->MM, pars->ZERO, w[33]);
//   FFV2_2(w[32], w[16], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[34]);
//   FFV5_7_2(w[32], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
//       w[35]);
  FFV2_1(w[2], w[6], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[36]);
  FFV2P0_3(w[3], w[36], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[37]);
  FFV5_7_3(w[3], w[36], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[38]);
//   FFV1_2_3_4_3(w[3], w[36], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, pars->MZp, pars->WZp, w[39]);
  FFV2P0_3(w[5], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[40]);
  FFV1_2_3_4_1(w[36], w[7], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
      pars->HEF_MEKD2_1_GC_130, pars->Me, pars->ZERO, w[41]);
  FFV1_2_3_4_2(w[3], w[7], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
      pars->HEF_MEKD2_1_GC_130, pars->Me, pars->ZERO, w[42]);
  FFV5_7_3(w[5], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[43]);
//   VVV1_2_1(w[43], w[7], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, pars->MZ, pars->WZ, w[44]);
//   FFV1_2_3_4_3(w[5], w[4], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, pars->MZp, pars->WZp, w[45]);
//   FFV2_1(w[36], w[16], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[46]);
//   FFV2_2(w[3], w[16], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[47]);
//   FFV5_7_1(w[36], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->Me, pars->ZERO,
//       w[48]);
//   FFV5_7_2(w[3], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->Me, pars->ZERO,
//       w[49]);
//   VVV1_2_1(w[22], w[45], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, pars->MZ, pars->WZ, w[50]);
//   VVV1_2_3(w[22], w[43], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, pars->MZp, pars->WZp, w[51]);
  FFV2_2(w[3], w[6], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[52]);
  FFV2P0_3(w[52], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[53]);
  FFV5_7_3(w[52], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[54]);
//   FFV1_2_3_4_3(w[52], w[2], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, pars->MZp, pars->WZp, w[55]);
  FFV1_2_3_4_1(w[2], w[7], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
      pars->HEF_MEKD2_1_GC_130, pars->Me, pars->ZERO, w[56]);
  FFV1_2_3_4_2(w[52], w[7], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
      pars->HEF_MEKD2_1_GC_130, pars->Me, pars->ZERO, w[57]);
//   FFV2_1(w[2], w[16], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[58]);
//   FFV2_2(w[52], w[16], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[59]);
//   FFV5_7_1(w[2], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->Me, pars->ZERO,
//       w[60]);
//   FFV5_7_2(w[52], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->Me, pars->ZERO,
//       w[61]);
  FFV2_2(w[3], w[40], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[62]);
  FFV2_1(w[2], w[40], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[63]);
  FFV5_7_2(w[3], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->Me, pars->ZERO,
      w[64]);
  FFV5_7_1(w[2], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->Me, pars->ZERO,
      w[65]);
//   FFV1_2_3_4_2(w[3], w[45], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, pars->Me, pars->ZERO, w[66]);
//   FFV1_2_3_4_1(w[2], w[45], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, pars->Me, pars->ZERO, w[67]);
  FFV2P0_3(w[5], w[25], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[68]);
  FFV5_7_3(w[5], w[25], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[69]);
//   FFV1_2_3_4_3(w[5], w[25], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, pars->MZp, pars->WZp, w[70]);
  FFV2P0_3(w[32], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[71]);
  FFV5_7_3(w[32], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[72]);
//   FFV1_2_3_4_3(w[32], w[4], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, pars->MZp, pars->WZp, w[73]);
//   FFV2_2(w[0], w[6], pars->HEF_MEKD2_1_GC_3, pars->MS, pars->ZERO, w[74]);
//   FFV2P0_3(w[74], w[1], pars->HEF_MEKD2_1_GC_3, pars->ZERO, pars->ZERO, w[75]);
//   FFV5_6_3(w[74], w[1], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, pars->MZ, pars->WZ, w[76]);
//   FFV1_2_3_4_3(w[74], w[1], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, pars->MZp, pars->WZp, w[77]);
//   FFV1_2_3_4_2(w[74], w[17], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, pars->MS, pars->ZERO, w[78]);
//   FFV2_2(w[74], w[40], pars->HEF_MEKD2_1_GC_3, pars->MS, pars->ZERO, w[79]);
//   FFV5_6_2(w[74], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, pars->MS, pars->ZERO,
//       w[80]);
//   FFV2_2(w[74], w[8], pars->HEF_MEKD2_1_GC_3, pars->MS, pars->ZERO, w[81]);
//   FFV1_2_3_4_2(w[74], w[45], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, pars->MS, pars->ZERO, w[82]);
//   FFV5_6_2(w[74], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, pars->MS, pars->ZERO,
//       w[83]);
//   FFV2_1(w[1], w[6], pars->HEF_MEKD2_1_GC_3, pars->MS, pars->ZERO, w[84]);
//   FFV2P0_3(w[0], w[84], pars->HEF_MEKD2_1_GC_3, pars->ZERO, pars->ZERO, w[85]);
//   FFV5_6_3(w[0], w[84], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, pars->MZ, pars->WZ, w[86]);
//   FFV1_2_3_4_3(w[0], w[84], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, pars->MZp, pars->WZp, w[87]);
//   FFV1_2_3_4_2(w[0], w[17], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, pars->MS, pars->ZERO, w[88]);
//   FFV2_2(w[0], w[40], pars->HEF_MEKD2_1_GC_3, pars->MS, pars->ZERO, w[89]);
//   FFV5_6_2(w[0], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, pars->MS, pars->ZERO,
//       w[90]);
//   FFV2_2(w[0], w[8], pars->HEF_MEKD2_1_GC_3, pars->MS, pars->ZERO, w[91]);
//   FFV1_2_3_4_2(w[0], w[45], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, pars->MS, pars->ZERO, w[92]);
//   FFV5_6_2(w[0], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, pars->MS, pars->ZERO,
//       w[93]);
//   FFV2_1(w[1], w[40], pars->HEF_MEKD2_1_GC_3, pars->MS, pars->ZERO, w[94]);
//   FFV1_2_3_4_1(w[1], w[17], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, pars->MS, pars->ZERO, w[95]);
//   FFV5_6_1(w[1], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, pars->MS, pars->ZERO,
//       w[96]);
//   FFV1_2_3_4_1(w[1], w[45], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, pars->MS, pars->ZERO, w[97]);
//   FFV2_1(w[1], w[8], pars->HEF_MEKD2_1_GC_3, pars->MS, pars->ZERO, w[98]);
//   FFV5_6_1(w[1], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, pars->MS, pars->ZERO,
//       w[99]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV2_0(w[10], w[9], w[6], pars->HEF_MEKD2_1_GC_5, amp[0]);
  FFV2_0(w[11], w[12], w[6], pars->HEF_MEKD2_1_GC_5, amp[1]);
  FFV2_0(w[14], w[9], w[6], pars->HEF_MEKD2_1_GC_5, amp[2]);
  FFV2_0(w[11], w[15], w[6], pars->HEF_MEKD2_1_GC_5, amp[3]);
//   FFV2_0(w[19], w[18], w[6], pars->HEF_MEKD2_1_GC_5, amp[4]);	// s ch. A
//   FFV2_0(w[20], w[21], w[6], pars->HEF_MEKD2_1_GC_5, amp[5]);	// s ch. A
//   FFV2_0(w[19], w[23], w[6], pars->HEF_MEKD2_1_GC_5, amp[6]);	// s ch. Z
//   FFV2_0(w[24], w[21], w[6], pars->HEF_MEKD2_1_GC_5, amp[7]);	// s ch. Z
  FFV2_0(w[5], w[26], w[8], pars->HEF_MEKD2_1_GC_5, amp[8]);
  FFV2_0(w[11], w[25], w[8], pars->HEF_MEKD2_1_GC_5, amp[9]);
//   FFV5_7_0(w[5], w[25], w[27], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[10]);	// ZZ
  FFV5_7_0(w[5], w[26], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[11]);
  FFV5_7_0(w[11], w[25], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[12]);
//   FFV1_2_3_4_0(w[5], w[28], w[17], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[13]);	// s ch. A
//   FFV1_2_3_4_0(w[20], w[25], w[17], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[14]);	// s ch. A
//   FFV5_7_0(w[5], w[25], w[29], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[15]);	// s ch. Z, VBF
//   FFV1_2_3_4_0(w[5], w[30], w[17], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[16]);	// s ch. Z
//   FFV1_2_3_4_0(w[24], w[25], w[17], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[17]);	// s ch. Z
//   FFV1_2_3_4_0(w[5], w[25], w[31], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[18]);	// s ch. Z, VBF
  FFV2_0(w[32], w[9], w[8], pars->HEF_MEKD2_1_GC_5, amp[19]);
  FFV2_0(w[33], w[4], w[8], pars->HEF_MEKD2_1_GC_5, amp[20]);
//   FFV5_7_0(w[32], w[4], w[27], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[21]);	// ZZ
  FFV5_7_0(w[32], w[9], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[22]);
  FFV5_7_0(w[33], w[4], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[23]);
//   FFV1_2_3_4_0(w[32], w[18], w[17], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[24]);	// s ch. A
//   FFV1_2_3_4_0(w[34], w[4], w[17], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[25]);	// s ch. A
//   FFV5_7_0(w[32], w[4], w[29], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[26]);	// s ch. Z, VBF
//   FFV1_2_3_4_0(w[32], w[23], w[17], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[27]);	// s ch. Z
//   FFV1_2_3_4_0(w[35], w[4], w[17], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[28]);	// s ch. Z
//   FFV1_2_3_4_0(w[32], w[4], w[31], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[29]);	// s ch. Z, VBF
  FFV2_0(w[5], w[9], w[37], pars->HEF_MEKD2_1_GC_5, amp[30]);
  FFV5_7_0(w[5], w[9], w[38], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[31]);
  FFV2_0(w[11], w[4], w[37], pars->HEF_MEKD2_1_GC_5, amp[32]);
  FFV5_7_0(w[11], w[4], w[38], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[33]);
//   FFV1_2_3_4_0(w[5], w[18], w[39], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[34]);	// s ch. A
//   FFV1_2_3_4_0(w[20], w[4], w[39], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[35]);	// s ch. A
//   FFV1_2_3_4_0(w[5], w[23], w[39], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[36]);	// s ch. Z
//   FFV1_2_3_4_0(w[24], w[4], w[39], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[37]);	// s ch. Z
  FFV2_0(w[3], w[41], w[40], pars->HEF_MEKD2_1_GC_5, amp[38]);
  FFV2_0(w[42], w[36], w[40], pars->HEF_MEKD2_1_GC_5, amp[39]);
  FFV5_7_0(w[3], w[41], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[40]);
  FFV5_7_0(w[42], w[36], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[41]);
//   FFV5_7_0(w[3], w[36], w[44], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[42]);	// ZZ
//   FFV1_2_3_4_0(w[3], w[46], w[45], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[43]);	// s ch. A
//   FFV1_2_3_4_0(w[47], w[36], w[45], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[44]);	// s ch. A
//   FFV1_2_3_4_0(w[3], w[48], w[45], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[45]);	// s ch. Z
//   FFV1_2_3_4_0(w[49], w[36], w[45], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[46]);	// s ch. Z
//   FFV5_7_0(w[3], w[36], w[50], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[47]);	// s ch. Z, VBF
//   FFV1_2_3_4_0(w[3], w[36], w[51], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[48]);	// s ch. Z, VBF
  FFV2_0(w[5], w[9], w[53], pars->HEF_MEKD2_1_GC_5, amp[49]);
  FFV5_7_0(w[5], w[9], w[54], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[50]);
  FFV2_0(w[11], w[4], w[53], pars->HEF_MEKD2_1_GC_5, amp[51]);
  FFV5_7_0(w[11], w[4], w[54], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[52]);
//   FFV1_2_3_4_0(w[5], w[18], w[55], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[53]);	// s ch. A
//   FFV1_2_3_4_0(w[20], w[4], w[55], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[54]);	// s ch. A
//   FFV1_2_3_4_0(w[5], w[23], w[55], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[55]);	// s ch. Z
//   FFV1_2_3_4_0(w[24], w[4], w[55], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[56]);	// s ch. Z
  FFV2_0(w[52], w[56], w[40], pars->HEF_MEKD2_1_GC_5, amp[57]);
  FFV2_0(w[57], w[2], w[40], pars->HEF_MEKD2_1_GC_5, amp[58]);
  FFV5_7_0(w[52], w[56], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[59]);
  FFV5_7_0(w[57], w[2], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[60]);
//   FFV5_7_0(w[52], w[2], w[44], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[61]);	// ZZ
//   FFV1_2_3_4_0(w[52], w[58], w[45], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[62]);	// s ch. A
//   FFV1_2_3_4_0(w[59], w[2], w[45], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[63]);	// s ch. A
//   FFV1_2_3_4_0(w[52], w[60], w[45], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[64]);	// s ch. Z
//   FFV1_2_3_4_0(w[61], w[2], w[45], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[65]);	// s ch. Z
//   FFV5_7_0(w[52], w[2], w[50], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[66]);	// s ch. Z, VBF
//   FFV1_2_3_4_0(w[52], w[2], w[51], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[67]);	// s ch. Z, VBF
  FFV2_0(w[62], w[56], w[6], pars->HEF_MEKD2_1_GC_5, amp[68]);
  FFV2_0(w[42], w[63], w[6], pars->HEF_MEKD2_1_GC_5, amp[69]);
  FFV2_0(w[64], w[56], w[6], pars->HEF_MEKD2_1_GC_5, amp[70]);
  FFV2_0(w[42], w[65], w[6], pars->HEF_MEKD2_1_GC_5, amp[71]);
//   FFV2_0(w[66], w[58], w[6], pars->HEF_MEKD2_1_GC_5, amp[72]);	// s ch. A
//   FFV2_0(w[47], w[67], w[6], pars->HEF_MEKD2_1_GC_5, amp[73]);	// s ch. A
//   FFV2_0(w[66], w[60], w[6], pars->HEF_MEKD2_1_GC_5, amp[74]);	// s ch. Z
//   FFV2_0(w[49], w[67], w[6], pars->HEF_MEKD2_1_GC_5, amp[75]);	// s ch. Z
  FFV2_0(w[3], w[56], w[68], pars->HEF_MEKD2_1_GC_5, amp[76]);
  FFV5_7_0(w[3], w[56], w[69], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[77]);
  FFV2_0(w[42], w[2], w[68], pars->HEF_MEKD2_1_GC_5, amp[78]);
  FFV5_7_0(w[42], w[2], w[69], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[79]);
//   FFV1_2_3_4_0(w[3], w[58], w[70], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[80]);	// s ch. A
//   FFV1_2_3_4_0(w[47], w[2], w[70], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[81]);	// s ch. A
//   FFV1_2_3_4_0(w[3], w[60], w[70], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[82]);	// s ch. Z
//   FFV1_2_3_4_0(w[49], w[2], w[70], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[83]);	// s ch. Z
  FFV2_0(w[3], w[56], w[71], pars->HEF_MEKD2_1_GC_5, amp[84]);
  FFV5_7_0(w[3], w[56], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[85]);
  FFV2_0(w[42], w[2], w[71], pars->HEF_MEKD2_1_GC_5, amp[86]);
  FFV5_7_0(w[42], w[2], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[87]);
//   FFV1_2_3_4_0(w[3], w[58], w[73], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[88]);	// s ch. A
//   FFV1_2_3_4_0(w[47], w[2], w[73], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[89]);	// s ch. A
//   FFV1_2_3_4_0(w[3], w[60], w[73], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[90]);	// s ch. Z
//   FFV1_2_3_4_0(w[49], w[2], w[73], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[91]);	// s ch. Z
//   FFV2_0(w[5], w[21], w[75], pars->HEF_MEKD2_1_GC_5, amp[92]);	// s ch. A
//   FFV5_7_0(w[5], w[21], w[76], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[93]);	// s ch. Z
//   FFV2_0(w[19], w[4], w[75], pars->HEF_MEKD2_1_GC_5, amp[94]);	// s ch. A
//   FFV5_7_0(w[19], w[4], w[76], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[95]);	// s ch. Z
//   FFV1_2_3_4_0(w[5], w[12], w[77], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[96]);	// ISR
//   FFV1_2_3_4_0(w[10], w[4], w[77], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[97]);	// ISR
//   FFV1_2_3_4_0(w[5], w[15], w[77], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[98]);	// ISR
//   FFV1_2_3_4_0(w[14], w[4], w[77], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[99]);	// ISR
//   FFV2_0(w[78], w[1], w[40], pars->HEF_MEKD2_1_GC_3, amp[100]);	// t ch.
//   FFV1_2_3_4_0(w[79], w[1], w[17], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[101]);	// t ch.
//   VVV1_2_0(w[43], w[76], w[17], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[102]);	// t ch.
//   FFV5_6_0(w[78], w[1], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[103]);	// t ch.
//   FFV1_2_3_4_0(w[80], w[1], w[17], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[104]);	// t ch.
//   FFV1_2_3_4_0(w[81], w[1], w[45], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[105]);	// t ch.
//   FFV2_0(w[82], w[1], w[8], pars->HEF_MEKD2_1_GC_3, amp[106]);	// t ch.
//   VVV1_2_0(w[13], w[76], w[45], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[107]);	// t ch.
//   FFV1_2_3_4_0(w[83], w[1], w[45], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[108]);	// t ch.
//   FFV5_6_0(w[82], w[1], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[109]);	// t ch.
//   VVV1_2_0(w[13], w[43], w[77], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[110]);	// s ch. Z, VBF, ISR
//   FFV2_0(w[3], w[67], w[75], pars->HEF_MEKD2_1_GC_5, amp[111]);
//   FFV5_7_0(w[3], w[67], w[76], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[112]);
//   FFV2_0(w[66], w[2], w[75], pars->HEF_MEKD2_1_GC_5, amp[113]);
//   FFV5_7_0(w[66], w[2], w[76], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[114]);
//   FFV1_2_3_4_0(w[3], w[63], w[77], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[115]);	// ISR
//   FFV1_2_3_4_0(w[62], w[2], w[77], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[116]);	// ISR
//   FFV1_2_3_4_0(w[3], w[65], w[77], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[117]);	// ISR
//   FFV1_2_3_4_0(w[64], w[2], w[77], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[118]);	// ISR
//   FFV2_0(w[5], w[21], w[85], pars->HEF_MEKD2_1_GC_5, amp[119]);
//   FFV5_7_0(w[5], w[21], w[86], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[120]);
//   FFV2_0(w[19], w[4], w[85], pars->HEF_MEKD2_1_GC_5, amp[121]);
//   FFV5_7_0(w[19], w[4], w[86], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[122]);
//   FFV1_2_3_4_0(w[5], w[12], w[87], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[123]);	// ISR
//   FFV1_2_3_4_0(w[10], w[4], w[87], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[124]);	// ISR
//   FFV1_2_3_4_0(w[5], w[15], w[87], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[125]);	// ISR
//   FFV1_2_3_4_0(w[14], w[4], w[87], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[126]);	// ISR
//   FFV2_0(w[88], w[84], w[40], pars->HEF_MEKD2_1_GC_3, amp[127]);
//   FFV1_2_3_4_0(w[89], w[84], w[17], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[128]);
//   VVV1_2_0(w[43], w[86], w[17], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[129]);
//   FFV5_6_0(w[88], w[84], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[130]);
//   FFV1_2_3_4_0(w[90], w[84], w[17], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[131]);
//   FFV1_2_3_4_0(w[91], w[84], w[45], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[132]);
//   FFV2_0(w[92], w[84], w[8], pars->HEF_MEKD2_1_GC_3, amp[133]);
//   VVV1_2_0(w[13], w[86], w[45], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[134]);
//   FFV1_2_3_4_0(w[93], w[84], w[45], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[135]);
//   FFV5_6_0(w[92], w[84], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[136]);
//   VVV1_2_0(w[13], w[43], w[87], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[137]);
//   FFV2_0(w[3], w[67], w[85], pars->HEF_MEKD2_1_GC_5, amp[138]);
//   FFV5_7_0(w[3], w[67], w[86], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[139]);
//   FFV2_0(w[66], w[2], w[85], pars->HEF_MEKD2_1_GC_5, amp[140]);
//   FFV5_7_0(w[66], w[2], w[86], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[141]);
//   FFV1_2_3_4_0(w[3], w[63], w[87], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[142]);	// ISR
//   FFV1_2_3_4_0(w[62], w[2], w[87], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[143]);	// ISR
//   FFV1_2_3_4_0(w[3], w[65], w[87], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[144]);	// ISR
//   FFV1_2_3_4_0(w[64], w[2], w[87], pars->HEF_MEKD2_1_GC_131, pars->HEF_MEKD2_1_GC_128, pars->HEF_MEKD2_1_GC_129,
//       pars->HEF_MEKD2_1_GC_130, amp[145]);	// ISR
//   FFV2_0(w[88], w[94], w[6], pars->HEF_MEKD2_1_GC_3, amp[146]);
//   FFV2_0(w[89], w[95], w[6], pars->HEF_MEKD2_1_GC_3, amp[147]);
//   FFV2_0(w[88], w[96], w[6], pars->HEF_MEKD2_1_GC_3, amp[148]);
//   FFV2_0(w[90], w[95], w[6], pars->HEF_MEKD2_1_GC_3, amp[149]);
//   FFV2_0(w[91], w[97], w[6], pars->HEF_MEKD2_1_GC_3, amp[150]);
//   FFV2_0(w[92], w[98], w[6], pars->HEF_MEKD2_1_GC_3, amp[151]);
//   FFV2_0(w[93], w[97], w[6], pars->HEF_MEKD2_1_GC_3, amp[152]);
//   FFV2_0(w[92], w[99], w[6], pars->HEF_MEKD2_1_GC_3, amp[153]);
//   FFV2_0(w[88], w[1], w[68], pars->HEF_MEKD2_1_GC_3, amp[154]);
//   FFV5_6_0(w[88], w[1], w[69], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[155]);
//   FFV2_0(w[0], w[95], w[68], pars->HEF_MEKD2_1_GC_3, amp[156]);
//   FFV5_6_0(w[0], w[95], w[69], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[157]);
//   FFV1_2_3_4_0(w[91], w[1], w[70], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[158]);
//   FFV1_2_3_4_0(w[0], w[98], w[70], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[159]);
//   FFV1_2_3_4_0(w[93], w[1], w[70], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[160]);
//   FFV1_2_3_4_0(w[0], w[99], w[70], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[161]);
//   FFV2_0(w[88], w[1], w[71], pars->HEF_MEKD2_1_GC_3, amp[162]);
//   FFV5_6_0(w[88], w[1], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[163]);
//   FFV2_0(w[0], w[95], w[71], pars->HEF_MEKD2_1_GC_3, amp[164]);
//   FFV5_6_0(w[0], w[95], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[165]);
//   FFV1_2_3_4_0(w[91], w[1], w[73], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[166]);
//   FFV1_2_3_4_0(w[0], w[98], w[73], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[167]);
//   FFV1_2_3_4_0(w[93], w[1], w[73], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[168]);
//   FFV1_2_3_4_0(w[0], w[99], w[73], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[169]);
//   FFV2_0(w[92], w[1], w[37], pars->HEF_MEKD2_1_GC_3, amp[170]);
//   FFV5_6_0(w[92], w[1], w[38], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[171]);
//   FFV2_0(w[0], w[97], w[37], pars->HEF_MEKD2_1_GC_3, amp[172]);
//   FFV5_6_0(w[0], w[97], w[38], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[173]);
//   FFV1_2_3_4_0(w[89], w[1], w[39], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[174]);
//   FFV1_2_3_4_0(w[0], w[94], w[39], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[175]);
//   FFV1_2_3_4_0(w[90], w[1], w[39], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[176]);
//   FFV1_2_3_4_0(w[0], w[96], w[39], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[177]);
//   FFV2_0(w[92], w[1], w[53], pars->HEF_MEKD2_1_GC_3, amp[178]);
//   FFV5_6_0(w[92], w[1], w[54], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[179]);
//   FFV2_0(w[0], w[97], w[53], pars->HEF_MEKD2_1_GC_3, amp[180]);
//   FFV5_6_0(w[0], w[97], w[54], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_187, amp[181]);
//   FFV1_2_3_4_0(w[89], w[1], w[55], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[182]);
//   FFV1_2_3_4_0(w[0], w[94], w[55], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[183]);
//   FFV1_2_3_4_0(w[90], w[1], w[55], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[184]);
//   FFV1_2_3_4_0(w[0], w[96], w[55], pars->HEF_MEKD2_1_GC_151, pars->HEF_MEKD2_1_GC_148, pars->HEF_MEKD2_1_GC_149,
//       pars->HEF_MEKD2_1_GC_150, amp[185]);

}
double qq_Spin1_2f_DN_OFpA::matrix_ssx_zp_emepmummupa_no_hxg() 
{
  int i, j;
  // Local variables
	// Commented out by Convert_source 0.2
//  const int ngraphs = 186;
  const int ncolor = 1;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[ncolor] = {1};
  static const double cf[ncolor][ncolor] = {{3}};

  // Calculate color flows
  jamp[0] = +amp[0] + amp[1] + amp[2] + amp[3] + amp[4] + amp[5] + amp[6] +
      amp[7] + amp[8] + amp[9] + amp[10] + amp[11] + amp[12] + amp[13] +
      amp[14] + amp[15] + amp[16] + amp[17] + amp[18] + amp[19] + amp[20] +
      amp[21] + amp[22] + amp[23] + amp[24] + amp[25] + amp[26] + amp[27] +
      amp[28] + amp[29] + amp[30] + amp[31] + amp[32] + amp[33] + amp[34] +
      amp[35] + amp[36] + amp[37] + amp[38] + amp[39] + amp[40] + amp[41] +
      amp[42] + amp[43] + amp[44] + amp[45] + amp[46] + amp[47] + amp[48] +
      amp[49] + amp[50] + amp[51] + amp[52] + amp[53] + amp[54] + amp[55] +
      amp[56] + amp[57] + amp[58] + amp[59] + amp[60] + amp[61] + amp[62] +
      amp[63] + amp[64] + amp[65] + amp[66] + amp[67] + amp[68] + amp[69] +
      amp[70] + amp[71] + amp[72] + amp[73] + amp[74] + amp[75] + amp[76] +
      amp[77] + amp[78] + amp[79] + amp[80] + amp[81] + amp[82] + amp[83] +
      amp[84] + amp[85] + amp[86] + amp[87] + amp[88] + amp[89] + amp[90] +
      amp[91] + amp[92] + amp[93] + amp[94] + amp[95] + amp[96] + amp[97] +
      amp[98] + amp[99] + amp[100] + amp[101] + amp[102] + amp[103] + amp[104]
      + amp[105] + amp[106] + amp[107] + amp[108] + amp[109] + amp[110] +
      amp[111] + amp[112] + amp[113] + amp[114] + amp[115] + amp[116] +
      amp[117] + amp[118] + amp[119] + amp[120] + amp[121] + amp[122] +
      amp[123] + amp[124] + amp[125] + amp[126] + amp[127] + amp[128] +
      amp[129] + amp[130] + amp[131] + amp[132] + amp[133] + amp[134] +
      amp[135] + amp[136] + amp[137] + amp[138] + amp[139] + amp[140] +
      amp[141] + amp[142] + amp[143] + amp[144] + amp[145] + amp[146] +
      amp[147] + amp[148] + amp[149] + amp[150] + amp[151] + amp[152] +
      amp[153] + amp[154] + amp[155] + amp[156] + amp[157] + amp[158] +
      amp[159] + amp[160] + amp[161] + amp[162] + amp[163] + amp[164] +
      amp[165] + amp[166] + amp[167] + amp[168] + amp[169] + amp[170] +
      amp[171] + amp[172] + amp[173] + amp[174] + amp[175] + amp[176] +
      amp[177] + amp[178] + amp[179] + amp[180] + amp[181] + amp[182] +
      amp[183] + amp[184] + amp[185];

  // Sum and square the color flows to get the matrix element
  double matrix = 0;
  for(i = 0;i < ncolor;i++ )
  {
    ztemp = 0.;
    for(j = 0;j < ncolor;j++ )
      ztemp = ztemp + cf[i][j] * jamp[j];
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i];
  }

  // Store the leading color flows for choice of color
  for(i = 0;i < ncolor;i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i]));

  return matrix;
}



