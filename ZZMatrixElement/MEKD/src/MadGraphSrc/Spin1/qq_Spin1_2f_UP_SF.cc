//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.2, 2014-02-07
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "qq_Spin1_2f_UP_SF.h"
#include "../HelAmps_HEF_MEKD2_1.h"	// Changed by Convert_source 0.2

using namespace MG5_HEF_MEKD2_1;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: c c~ > zp > mu- mu+ mu- mu+ S1VV=1 QED=2 S1QQ=2 / h xg

//--------------------------------------------------------------------------
// Initialize process.

void qq_Spin1_2f_UP_SF::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MEKD::getInstance();	// Changed by Convert_source 0.2 
  SLHAReader_MEKD slha(param_card_name);
  pars->setIndependentParameters(slha);
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// Moved here by Convert_source 0.2
  // Set external particle masses for this matrix element
  mME.push_back(pars->MC);
  mME.push_back(pars->MC);
  mME.push_back(pars->MM);
  mME.push_back(pars->MM);
  mME.push_back(pars->MM);
  mME.push_back(pars->MM);
  jamp2[0] = new double[1];
	for( int count=0; count<namplitudes; count++ ) amp[count] = 0;
}

//--------------------------------------------------------------------------
// Update process.	// Created here by Convert_source 0.2

void qq_Spin1_2f_UP_SF::updateProc(SLHAReader_MEKD &slha)
{
	pars->setIndependentParameters(slha);
	pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// needed when altering couplings
	
	// Set external particle masses for this matrix element
	// Should correspond to initProc
	mME[0]=(pars->MC);
	mME[1]=(pars->MC);
	mME[2]=(pars->MM);
	mME[3]=(pars->MM);
	mME[4]=(pars->MM);
	mME[5]=(pars->MM);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void qq_Spin1_2f_UP_SF::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters();
  pars->setDependentCouplings();
	// Deleted by Convert_source 0.2
	
  // Reset color flows
  for(int i = 0;i < 1;i++ )
    jamp2[0][i] = 0.;

  // Local variables and constants
  const int ncomb = 64;
  static bool goodhel[ncomb] = {ncomb * false};
//	static int ntry = 0, sum_hel = 0, ngood = 0;	// Moved by Convert_source 0.2
  static int igood[ncomb];
  static int jhel;
//	std::complex<double> * * wfs;	// Changed by Convert_source 0.2
  double t[nprocesses];
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1},
      {-1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1},
      {-1, -1, -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1},
      {-1, -1, -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1},
      {-1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1},
      {-1, -1, 1, 1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1}, {-1,
      1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1, -1, 1}, {-1, 1,
      -1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, 1}, {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, 1}, {-1, 1, 1, 1,
      -1, -1}, {-1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1},
      {1, -1, -1, -1, -1, -1}, {1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, 1, -1},
      {1, -1, -1, -1, 1, 1}, {1, -1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1}, {1,
      -1, -1, 1, 1, -1}, {1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1,
      1, -1, -1, 1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1,
      -1, -1}, {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1},
      {1, 1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1,
      1, -1, -1, 1, 1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1,
      1, 1, -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1,
      1}, {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1}, {1,
      1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {144, 144};

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
        t[0] = matrix_ccx_zp_mummupmummup_no_hxg();
        // Mirror initial state momenta for mirror process
        perm[0] = 1;
        perm[1] = 0;
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]);
        // Mirror back
        perm[0] = 0;
        perm[1] = 1;
        // Calculate matrix elements
        t[1] = matrix_ccx_zp_mummupmummup_no_hxg();
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
      t[0] = matrix_ccx_zp_mummupmummup_no_hxg();
      // Mirror initial state momenta for mirror process
      perm[0] = 1;
      perm[1] = 0;
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]);
      // Mirror back
      perm[0] = 0;
      perm[1] = 1;
      // Calculate matrix elements
      t[1] = matrix_ccx_zp_mummupmummup_no_hxg();
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

double qq_Spin1_2f_UP_SF::sigmaHat() 
{
  // Select between the different processes
  if(id1 == -4 && id2 == 4)
  {
    // Add matrix elements for processes with beams (-4, 4)
    return matrix_element[1];
  }
  else if(id1 == 4 && id2 == -4)
  {
    // Add matrix elements for processes with beams (4, -4)
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

void qq_Spin1_2f_UP_SF::calculate_wavefunctions(const int perm[], const int hel[])
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
  FFV1_2_3_4_3(w[0], w[1], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
      pars->HEF_MEKD2_1_GC_110, pars->MZp, pars->WZp, w[6]);
  FFV2P0_3(w[3], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[7]);
  FFV1_2_3_4_1(w[4], w[6], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, pars->MM, pars->ZERO, w[8]);
  FFV1_2_3_4_2(w[5], w[6], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, pars->MM, pars->ZERO, w[9]);
  FFV5_7_3(w[3], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[10]);
//   FFV2P0_3(w[0], w[1], pars->HEF_MEKD2_1_GC_4, pars->ZERO, pars->ZERO, w[11]);
//   FFV1_2_3_4_3(w[3], w[2], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, pars->MZp, pars->WZp, w[12]);
//   FFV2_1(w[4], w[11], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[13]);
//   FFV2_2(w[5], w[11], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[14]);
//   FFV5_8_3(w[0], w[1], pars->HEF_MEKD2_1_GC_182, pars->HEF_MEKD2_1_GC_187, pars->MZ, pars->WZ, w[15]);
//   FFV5_7_1(w[4], w[15], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
//       w[16]);
//   FFV5_7_2(w[5], w[15], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
//       w[17]);
  FFV5_7_3(w[5], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[18]);
//   FFV1_2_3_4_3(w[5], w[4], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, pars->MZp, pars->WZp, w[19]);
  FFV2P0_3(w[5], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[20]);
  FFV1_2_3_4_2(w[3], w[6], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, pars->MM, pars->ZERO, w[21]);
  FFV5_7_3(w[5], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[22]);
//   FFV1_2_3_4_3(w[5], w[2], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, pars->MZp, pars->WZp, w[23]);
//   FFV2_2(w[3], w[11], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[24]);
//   FFV5_7_2(w[3], w[15], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
//       w[25]);
  FFV5_7_3(w[3], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[26]);
//   FFV1_2_3_4_3(w[3], w[4], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, pars->MZp, pars->WZp, w[27]);
  FFV2P0_3(w[3], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[28]);
  FFV1_2_3_4_1(w[2], w[6], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, pars->MM, pars->ZERO, w[29]);
//   FFV2_1(w[2], w[11], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[30]);
//   FFV5_7_1(w[2], w[15], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
//       w[31]);
  FFV2P0_3(w[5], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[32]);
//   FFV1_2_3_4_2(w[0], w[12], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, pars->MC, pars->ZERO, w[33]);
//   FFV2_2(w[0], w[32], pars->HEF_MEKD2_1_GC_4, pars->MC, pars->ZERO, w[34]);
//   FFV5_8_2(w[0], w[18], pars->HEF_MEKD2_1_GC_182, pars->HEF_MEKD2_1_GC_187, pars->MC, pars->ZERO,
//       w[35]);
//   FFV2_2(w[0], w[7], pars->HEF_MEKD2_1_GC_4, pars->MC, pars->ZERO, w[36]);
//   FFV1_2_3_4_2(w[0], w[19], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, pars->MC, pars->ZERO, w[37]);
//   FFV5_8_2(w[0], w[10], pars->HEF_MEKD2_1_GC_182, pars->HEF_MEKD2_1_GC_187, pars->MC, pars->ZERO,
//       w[38]);
//   FFV1_2_3_4_2(w[0], w[23], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, pars->MC, pars->ZERO, w[39]);
//   FFV2_2(w[0], w[28], pars->HEF_MEKD2_1_GC_4, pars->MC, pars->ZERO, w[40]);
//   FFV5_8_2(w[0], w[26], pars->HEF_MEKD2_1_GC_182, pars->HEF_MEKD2_1_GC_187, pars->MC, pars->ZERO,
//       w[41]);
//   FFV2_2(w[0], w[20], pars->HEF_MEKD2_1_GC_4, pars->MC, pars->ZERO, w[42]);
//   FFV1_2_3_4_2(w[0], w[27], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, pars->MC, pars->ZERO, w[43]);
//   FFV5_8_2(w[0], w[22], pars->HEF_MEKD2_1_GC_182, pars->HEF_MEKD2_1_GC_187, pars->MC, pars->ZERO,
//       w[44]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV2_0(w[5], w[8], w[7], pars->HEF_MEKD2_1_GC_5, amp[0]);
  FFV2_0(w[9], w[4], w[7], pars->HEF_MEKD2_1_GC_5, amp[1]);
  FFV5_7_0(w[5], w[8], w[10], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[2]);
  FFV5_7_0(w[9], w[4], w[10], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[3]);
//   FFV1_2_3_4_0(w[5], w[13], w[12], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[4]);	// s ch. A
//   FFV1_2_3_4_0(w[14], w[4], w[12], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[5]);	// s ch. A
//   FFV1_2_3_4_0(w[5], w[16], w[12], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[6]);	// s ch. Z
//   FFV1_2_3_4_0(w[17], w[4], w[12], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[7]);	// s ch. Z
//   VVV1_2_0(w[10], w[18], w[6], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[8]);	// ZZ
//   VVV1_2_0(w[15], w[18], w[12], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[9]);	// s ch. Z, VBF
//   VVV1_2_0(w[15], w[10], w[19], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[10]);	// s ch. Z, VBF
  FFV2_0(w[21], w[4], w[20], pars->HEF_MEKD2_1_GC_5, amp[11]);
  FFV2_0(w[3], w[8], w[20], pars->HEF_MEKD2_1_GC_5, amp[12]);
  FFV5_7_0(w[21], w[4], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[13]);
  FFV5_7_0(w[3], w[8], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[14]);
//   FFV1_2_3_4_0(w[24], w[4], w[23], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[15]);	// s ch. A
//   FFV1_2_3_4_0(w[3], w[13], w[23], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[16]);	// s ch. A
//   FFV1_2_3_4_0(w[25], w[4], w[23], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[17]);	// s ch. Z
//   FFV1_2_3_4_0(w[3], w[16], w[23], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[18]);	// s ch. Z
//   VVV1_2_0(w[22], w[26], w[6], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[19]);	// ZZ
//   VVV1_2_0(w[15], w[26], w[23], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[20]);	// s ch. Z, VBF
//   VVV1_2_0(w[15], w[22], w[27], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[21]);	// s ch. Z, VBF
  FFV2_0(w[5], w[29], w[28], pars->HEF_MEKD2_1_GC_5, amp[22]);
  FFV2_0(w[9], w[2], w[28], pars->HEF_MEKD2_1_GC_5, amp[23]);
  FFV5_7_0(w[5], w[29], w[26], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[24]);
  FFV5_7_0(w[9], w[2], w[26], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[25]);
//   FFV1_2_3_4_0(w[5], w[30], w[27], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[26]);	// s ch. A
//   FFV1_2_3_4_0(w[14], w[2], w[27], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[27]);	// s ch. A
//   FFV1_2_3_4_0(w[5], w[31], w[27], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[28]);	// s ch. Z
//   FFV1_2_3_4_0(w[17], w[2], w[27], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[29]);	// s ch. Z
  FFV2_0(w[3], w[29], w[32], pars->HEF_MEKD2_1_GC_5, amp[30]);
  FFV2_0(w[21], w[2], w[32], pars->HEF_MEKD2_1_GC_5, amp[31]);
  FFV5_7_0(w[3], w[29], w[18], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[32]);
  FFV5_7_0(w[21], w[2], w[18], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[33]);
//   FFV1_2_3_4_0(w[3], w[30], w[19], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[34]);	// s ch. A
//   FFV1_2_3_4_0(w[24], w[2], w[19], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[35]);	// s ch. A
//   FFV1_2_3_4_0(w[3], w[31], w[19], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[36]);	// s ch. Z
//   FFV1_2_3_4_0(w[25], w[2], w[19], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
//       pars->HEF_MEKD2_1_GC_140, amp[37]);	// s ch. Z
//   FFV2_0(w[33], w[1], w[32], pars->HEF_MEKD2_1_GC_4, amp[38]);
//   FFV1_2_3_4_0(w[34], w[1], w[12], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, amp[39]);
//   FFV5_8_0(w[33], w[1], w[18], pars->HEF_MEKD2_1_GC_182, pars->HEF_MEKD2_1_GC_187, amp[40]);
//   FFV1_2_3_4_0(w[35], w[1], w[12], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, amp[41]);
//   FFV1_2_3_4_0(w[36], w[1], w[19], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, amp[42]);
//   FFV2_0(w[37], w[1], w[7], pars->HEF_MEKD2_1_GC_4, amp[43]);
//   FFV1_2_3_4_0(w[38], w[1], w[19], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, amp[44]);
//   FFV5_8_0(w[37], w[1], w[10], pars->HEF_MEKD2_1_GC_182, pars->HEF_MEKD2_1_GC_187, amp[45]);
//   FFV2_0(w[39], w[1], w[28], pars->HEF_MEKD2_1_GC_4, amp[46]);
//   FFV1_2_3_4_0(w[40], w[1], w[23], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, amp[47]);
//   FFV5_8_0(w[39], w[1], w[26], pars->HEF_MEKD2_1_GC_182, pars->HEF_MEKD2_1_GC_187, amp[48]);
//   FFV1_2_3_4_0(w[41], w[1], w[23], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, amp[49]);
//   FFV1_2_3_4_0(w[42], w[1], w[27], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, amp[50]);
//   FFV2_0(w[43], w[1], w[20], pars->HEF_MEKD2_1_GC_4, amp[51]);
//   FFV1_2_3_4_0(w[44], w[1], w[27], pars->HEF_MEKD2_1_GC_111, pars->HEF_MEKD2_1_GC_108, pars->HEF_MEKD2_1_GC_109,
//       pars->HEF_MEKD2_1_GC_110, amp[52]);
//   FFV5_8_0(w[43], w[1], w[22], pars->HEF_MEKD2_1_GC_182, pars->HEF_MEKD2_1_GC_187, amp[53]);

}
double qq_Spin1_2f_UP_SF::matrix_ccx_zp_mummupmummup_no_hxg() 
{
  int i, j;
  // Local variables
	// Commented out by Convert_source 0.2
//  const int ngraphs = 54;
  const int ncolor = 1;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[ncolor] = {1};
  static const double cf[ncolor][ncolor] = {{3}};

  // Calculate color flows
  jamp[0] = +amp[0] + amp[1] + amp[2] + amp[3] + amp[4] + amp[5] + amp[6] +
      amp[7] + amp[8] + amp[9] + amp[10] - amp[11] - amp[12] - amp[13] -
      amp[14] - amp[15] - amp[16] - amp[17] - amp[18] - amp[19] - amp[20] -
      amp[21] - amp[22] - amp[23] - amp[24] - amp[25] - amp[26] - amp[27] -
      amp[28] - amp[29] + amp[30] + amp[31] + amp[32] + amp[33] + amp[34] +
      amp[35] + amp[36] + amp[37] + amp[38] + amp[39] + amp[40] + amp[41] +
      amp[42] + amp[43] + amp[44] + amp[45] - amp[46] - amp[47] - amp[48] -
      amp[49] - amp[50] - amp[51] - amp[52] - amp[53];

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



