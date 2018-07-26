//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.2, 2014-02-07
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "gg_Spin0_2f_OFpA.h"
#include "../HelAmps_HEF_MEKD2_1.h"	// Changed by Convert_source 0.2

using namespace MG5_HEF_MEKD2_1;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > h > e- e+ mu- mu+ a S0QQ=1 QED=3 HIW=1 HIGS=1 HIG=1 HIWS=1 /
// zp xg

//--------------------------------------------------------------------------
// Initialize process.

void gg_Spin0_2f_OFpA::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MEKD::getInstance();	// Changed by Convert_source 0.2 
  SLHAReader_MEKD slha(param_card_name);
  pars->setIndependentParameters(slha);
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// Moved here by Convert_source 0.2
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO);
  mME.push_back(pars->ZERO);
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

void gg_Spin0_2f_OFpA::updateProc(SLHAReader_MEKD &slha)
{
	pars->setIndependentParameters(slha);
	pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// needed when altering couplings
	
	// Set external particle masses for this matrix element
	// Should correspond to initProc
	mME[0]=(pars->ZERO);
	mME[1]=(pars->ZERO);
	mME[2]=(pars->Me);
	mME[3]=(pars->Me);
	mME[4]=(pars->MM);
	mME[5]=(pars->MM);
	mME[6]=(pars->ZERO);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void gg_Spin0_2f_OFpA::sigmaKin() 
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
  const int denominators[nprocesses] = {256};

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
        t[0] = matrix_gg_h_emepmummupa_no_zpxg();

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
      t[0] = matrix_gg_h_emepmummupa_no_zpxg();

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

double gg_Spin0_2f_OFpA::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 21 && id2 == 21)
  {
    // Add matrix elements for processes with beams (21, 21)
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

void gg_Spin0_2f_OFpA::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
//	int i, j;	// Changed by Convert_source 0.2

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
  vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]);
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]);
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]);
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]);
  ixxxxx(p[perm[5]], mME[5], hel[5], -1, w[5]);
  vxxxxx(p[perm[6]], mME[6], hel[6], +1, w[6]);
  VVS3_4_5_3(w[0], w[1], pars->HEF_MEKD2_1_GC_13, pars->HEF_MEKD2_1_GC_15, pars->HEF_MEKD2_1_GC_19, pars->MH,
      pars->WH, w[7]);
  FFV2P0_3(w[3], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[8]);
  FFS1_2_1(w[4], w[7], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO, w[9]);
  FFV2_2(w[5], w[8], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[10]);
  FFS1_2_2(w[5], w[7], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[11]);
  FFV2_1(w[4], w[8], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[12]);
  FFV5_7_3(w[3], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[13]);
  FFV5_7_2(w[5], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[14]);
  FFV5_7_1(w[4], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[15]);
  VVS1_3(w[0], w[1], pars->HEF_MEKD2_1_GC_23, pars->MH, pars->WH, w[16]);
  FFS1_2_1(w[4], w[16], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[17]);
  FFS1_2_2(w[5], w[16], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[18]);
  FFV2_1(w[4], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[19]);
  FFS1_2_1(w[19], w[7], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[20]);
//   VVS3_4_5_1(w[13], w[7], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22, pars->MZ,
//       pars->WZ, w[21]);
//   VVS2_1(w[13], w[7], pars->HEF_MEKD2_1_GC_25, pars->MZ, pars->WZ, w[22]);
  FFS1_2_1(w[19], w[16], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[23]);
//   VVS3_4_5_1(w[13], w[16], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22, pars->MZ,
//       pars->WZ, w[24]);
//   VVS2_1(w[13], w[16], pars->HEF_MEKD2_1_GC_25, pars->MZ, pars->WZ, w[25]);
  FFV2_2(w[5], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[26]);
  FFS1_2_2(w[26], w[7], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[27]);
  FFS1_2_2(w[26], w[16], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[28]);
  FFV2_1(w[2], w[6], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[29]);
  FFV2P0_3(w[3], w[29], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[30]);
  FFV5_7_3(w[3], w[29], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[31]);
  FFV2P0_3(w[5], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[32]);
  FFS1_2_1(w[29], w[7], pars->HEF_MEKD2_1_GC_127, pars->HEF_MEKD2_1_GC_126, pars->Me, pars->ZERO,
      w[33]);
  FFS1_2_2(w[3], w[7], pars->HEF_MEKD2_1_GC_127, pars->HEF_MEKD2_1_GC_126, pars->Me, pars->ZERO,
      w[34]);
  FFV5_7_3(w[5], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[35]);
//   VVS3_4_5_1(w[35], w[7], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22, pars->MZ,
//       pars->WZ, w[36]);
//   VVS2_1(w[35], w[7], pars->HEF_MEKD2_1_GC_25, pars->MZ, pars->WZ, w[37]);
  FFS1_2_1(w[29], w[16], pars->HEF_MEKD2_1_GC_127, pars->HEF_MEKD2_1_GC_126, pars->Me, pars->ZERO,
      w[38]);
  FFS1_2_2(w[3], w[16], pars->HEF_MEKD2_1_GC_127, pars->HEF_MEKD2_1_GC_126, pars->Me, pars->ZERO,
      w[39]);
//   VVS3_4_5_1(w[35], w[16], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22, pars->MZ,
//       pars->WZ, w[40]);
//   VVS2_1(w[35], w[16], pars->HEF_MEKD2_1_GC_25, pars->MZ, pars->WZ, w[41]);
  FFV2_2(w[3], w[6], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[42]);
  FFV2P0_3(w[42], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[43]);
  FFV5_7_3(w[42], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[44]);
  FFS1_2_1(w[2], w[7], pars->HEF_MEKD2_1_GC_127, pars->HEF_MEKD2_1_GC_126, pars->Me, pars->ZERO,
      w[45]);
  FFS1_2_2(w[42], w[7], pars->HEF_MEKD2_1_GC_127, pars->HEF_MEKD2_1_GC_126, pars->Me, pars->ZERO,
      w[46]);
  FFS1_2_1(w[2], w[16], pars->HEF_MEKD2_1_GC_127, pars->HEF_MEKD2_1_GC_126, pars->Me, pars->ZERO,
      w[47]);
  FFS1_2_2(w[42], w[16], pars->HEF_MEKD2_1_GC_127, pars->HEF_MEKD2_1_GC_126, pars->Me, pars->ZERO,
      w[48]);
  FFV2_2(w[3], w[32], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[49]);
  FFV2_1(w[2], w[32], pars->HEF_MEKD2_1_GC_5, pars->Me, pars->ZERO, w[50]);
  FFV5_7_2(w[3], w[35], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->Me, pars->ZERO,
      w[51]);
  FFV5_7_1(w[2], w[35], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->Me, pars->ZERO,
      w[52]);
  FFV2P0_3(w[5], w[19], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[53]);
  FFV5_7_3(w[5], w[19], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[54]);
  FFV2P0_3(w[26], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[55]);
  FFV5_7_3(w[26], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[56]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV2_0(w[10], w[9], w[6], pars->HEF_MEKD2_1_GC_5, amp[0]);
  FFV2_0(w[11], w[12], w[6], pars->HEF_MEKD2_1_GC_5, amp[1]);
  FFV2_0(w[14], w[9], w[6], pars->HEF_MEKD2_1_GC_5, amp[2]);
  FFV2_0(w[11], w[15], w[6], pars->HEF_MEKD2_1_GC_5, amp[3]);
  FFV2_0(w[10], w[17], w[6], pars->HEF_MEKD2_1_GC_5, amp[4]);
  FFV2_0(w[18], w[12], w[6], pars->HEF_MEKD2_1_GC_5, amp[5]);
  FFV2_0(w[14], w[17], w[6], pars->HEF_MEKD2_1_GC_5, amp[6]);
  FFV2_0(w[18], w[15], w[6], pars->HEF_MEKD2_1_GC_5, amp[7]);
  FFV2_0(w[5], w[20], w[8], pars->HEF_MEKD2_1_GC_5, amp[8]);
  FFV2_0(w[11], w[19], w[8], pars->HEF_MEKD2_1_GC_5, amp[9]);
//   FFV5_7_0(w[5], w[19], w[21], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[10]);	// ZZ
//   FFV5_7_0(w[5], w[19], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[11]);	// ZZ
  FFV5_7_0(w[5], w[20], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[12]);
  FFV5_7_0(w[11], w[19], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[13]);
  FFV2_0(w[5], w[23], w[8], pars->HEF_MEKD2_1_GC_5, amp[14]);
  FFV2_0(w[18], w[19], w[8], pars->HEF_MEKD2_1_GC_5, amp[15]);
//   FFV5_7_0(w[5], w[19], w[24], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[16]);	// ZZ
//   FFV5_7_0(w[5], w[19], w[25], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[17]);	// ZZ
  FFV5_7_0(w[5], w[23], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[18]);
  FFV5_7_0(w[18], w[19], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[19]);
  FFV2_0(w[26], w[9], w[8], pars->HEF_MEKD2_1_GC_5, amp[20]);
  FFV2_0(w[27], w[4], w[8], pars->HEF_MEKD2_1_GC_5, amp[21]);
//   FFV5_7_0(w[26], w[4], w[21], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[22]);	// ZZ
//   FFV5_7_0(w[26], w[4], w[22], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[23]);	// ZZ
  FFV5_7_0(w[26], w[9], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[24]);
  FFV5_7_0(w[27], w[4], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[25]);
  FFV2_0(w[26], w[17], w[8], pars->HEF_MEKD2_1_GC_5, amp[26]);
  FFV2_0(w[28], w[4], w[8], pars->HEF_MEKD2_1_GC_5, amp[27]);
//   FFV5_7_0(w[26], w[4], w[24], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[28]);	// ZZ
//   FFV5_7_0(w[26], w[4], w[25], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[29]);	// ZZ
  FFV5_7_0(w[26], w[17], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[30]);
  FFV5_7_0(w[28], w[4], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[31]);
  FFV2_0(w[5], w[9], w[30], pars->HEF_MEKD2_1_GC_5, amp[32]);
  FFV5_7_0(w[5], w[9], w[31], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[33]);
  FFV2_0(w[11], w[4], w[30], pars->HEF_MEKD2_1_GC_5, amp[34]);
  FFV5_7_0(w[11], w[4], w[31], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[35]);
  FFV2_0(w[5], w[17], w[30], pars->HEF_MEKD2_1_GC_5, amp[36]);
  FFV5_7_0(w[5], w[17], w[31], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[37]);
  FFV2_0(w[18], w[4], w[30], pars->HEF_MEKD2_1_GC_5, amp[38]);
  FFV5_7_0(w[18], w[4], w[31], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[39]);
  FFV2_0(w[3], w[33], w[32], pars->HEF_MEKD2_1_GC_5, amp[40]);
  FFV2_0(w[34], w[29], w[32], pars->HEF_MEKD2_1_GC_5, amp[41]);
  FFV5_7_0(w[3], w[33], w[35], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[42]);
  FFV5_7_0(w[34], w[29], w[35], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[43]);
//   FFV5_7_0(w[3], w[29], w[36], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[44]);	// ZZ
//   FFV5_7_0(w[3], w[29], w[37], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[45]);	// ZZ
  FFV2_0(w[3], w[38], w[32], pars->HEF_MEKD2_1_GC_5, amp[46]);
  FFV2_0(w[39], w[29], w[32], pars->HEF_MEKD2_1_GC_5, amp[47]);
  FFV5_7_0(w[3], w[38], w[35], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[48]);
  FFV5_7_0(w[39], w[29], w[35], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[49]);
//   FFV5_7_0(w[3], w[29], w[40], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[50]);	// ZZ
//   FFV5_7_0(w[3], w[29], w[41], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[51]);	// ZZ
  FFV2_0(w[5], w[9], w[43], pars->HEF_MEKD2_1_GC_5, amp[52]);
  FFV5_7_0(w[5], w[9], w[44], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[53]);
  FFV2_0(w[11], w[4], w[43], pars->HEF_MEKD2_1_GC_5, amp[54]);
  FFV5_7_0(w[11], w[4], w[44], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[55]);
  FFV2_0(w[5], w[17], w[43], pars->HEF_MEKD2_1_GC_5, amp[56]);
  FFV5_7_0(w[5], w[17], w[44], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[57]);
  FFV2_0(w[18], w[4], w[43], pars->HEF_MEKD2_1_GC_5, amp[58]);
  FFV5_7_0(w[18], w[4], w[44], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[59]);
  FFV2_0(w[42], w[45], w[32], pars->HEF_MEKD2_1_GC_5, amp[60]);
  FFV2_0(w[46], w[2], w[32], pars->HEF_MEKD2_1_GC_5, amp[61]);
  FFV5_7_0(w[42], w[45], w[35], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[62]);
  FFV5_7_0(w[46], w[2], w[35], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[63]);
//   FFV5_7_0(w[42], w[2], w[36], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[64]);	// ZZ
//   FFV5_7_0(w[42], w[2], w[37], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[65]);	// ZZ
  FFV2_0(w[42], w[47], w[32], pars->HEF_MEKD2_1_GC_5, amp[66]);
  FFV2_0(w[48], w[2], w[32], pars->HEF_MEKD2_1_GC_5, amp[67]);
  FFV5_7_0(w[42], w[47], w[35], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[68]);
  FFV5_7_0(w[48], w[2], w[35], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[69]);
//   FFV5_7_0(w[42], w[2], w[40], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[70]);	// ZZ
//   FFV5_7_0(w[42], w[2], w[41], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[71]);	// ZZ
  FFV2_0(w[49], w[45], w[6], pars->HEF_MEKD2_1_GC_5, amp[72]);
  FFV2_0(w[34], w[50], w[6], pars->HEF_MEKD2_1_GC_5, amp[73]);
  FFV2_0(w[51], w[45], w[6], pars->HEF_MEKD2_1_GC_5, amp[74]);
  FFV2_0(w[34], w[52], w[6], pars->HEF_MEKD2_1_GC_5, amp[75]);
  FFV2_0(w[49], w[47], w[6], pars->HEF_MEKD2_1_GC_5, amp[76]);
  FFV2_0(w[39], w[50], w[6], pars->HEF_MEKD2_1_GC_5, amp[77]);
  FFV2_0(w[51], w[47], w[6], pars->HEF_MEKD2_1_GC_5, amp[78]);
  FFV2_0(w[39], w[52], w[6], pars->HEF_MEKD2_1_GC_5, amp[79]);
  FFV2_0(w[3], w[45], w[53], pars->HEF_MEKD2_1_GC_5, amp[80]);
  FFV5_7_0(w[3], w[45], w[54], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[81]);
  FFV2_0(w[34], w[2], w[53], pars->HEF_MEKD2_1_GC_5, amp[82]);
  FFV5_7_0(w[34], w[2], w[54], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[83]);
  FFV2_0(w[3], w[47], w[53], pars->HEF_MEKD2_1_GC_5, amp[84]);
  FFV5_7_0(w[3], w[47], w[54], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[85]);
  FFV2_0(w[39], w[2], w[53], pars->HEF_MEKD2_1_GC_5, amp[86]);
  FFV5_7_0(w[39], w[2], w[54], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[87]);
  FFV2_0(w[3], w[45], w[55], pars->HEF_MEKD2_1_GC_5, amp[88]);
  FFV5_7_0(w[3], w[45], w[56], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[89]);
  FFV2_0(w[34], w[2], w[55], pars->HEF_MEKD2_1_GC_5, amp[90]);
  FFV5_7_0(w[34], w[2], w[56], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[91]);
  FFV2_0(w[3], w[47], w[55], pars->HEF_MEKD2_1_GC_5, amp[92]);
  FFV5_7_0(w[3], w[47], w[56], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[93]);
  FFV2_0(w[39], w[2], w[55], pars->HEF_MEKD2_1_GC_5, amp[94]);
  FFV5_7_0(w[39], w[2], w[56], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[95]);

}
double gg_Spin0_2f_OFpA::matrix_gg_h_emepmummupa_no_zpxg() 
{
  int i, j;
  // Local variables
	// Commented out by Convert_source 0.2
//  const int ngraphs = 96;
  const int ncolor = 1;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[ncolor] = {1};
  static const double cf[ncolor][ncolor] = {{2}};

  // Calculate color flows
  jamp[0] = +2. * (+amp[0] + amp[1] + amp[2] + amp[3] + amp[4] + amp[5] +
      amp[6] + amp[7] + amp[8] + amp[9] + amp[10] + amp[11] + amp[12] + amp[13]
      + amp[14] + amp[15] + amp[16] + amp[17] + amp[18] + amp[19] + amp[20] +
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
      amp[91] + amp[92] + amp[93] + amp[94] + amp[95]);

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



