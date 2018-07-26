//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.2, 2014-02-07
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "gg_Spin0_2f_SFpA.h"
#include "../HelAmps_HEF_MEKD2_1.h"	// Changed by Convert_source 0.2

using namespace MG5_HEF_MEKD2_1;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > h > mu- mu+ mu- mu+ a S0QQ=1 QED=3 HIW=1 HIGS=1 HIG=1 HIWS=1
// / zp xg

//--------------------------------------------------------------------------
// Initialize process.

void gg_Spin0_2f_SFpA::initProc(string param_card_name) 
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
  mME.push_back(pars->MM);
  mME.push_back(pars->MM);
  mME.push_back(pars->MM);
  mME.push_back(pars->MM);
  mME.push_back(pars->ZERO);
  jamp2[0] = new double[1];
	for( int count=0; count<namplitudes; count++ ) amp[count] = 0;
}

//--------------------------------------------------------------------------
// Update process.	// Created here by Convert_source 0.2

void gg_Spin0_2f_SFpA::updateProc(SLHAReader_MEKD &slha)
{
	pars->setIndependentParameters(slha);
	pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// needed when altering couplings
	
	// Set external particle masses for this matrix element
	// Should correspond to initProc
	mME[0]=(pars->ZERO);
	mME[1]=(pars->ZERO);
	mME[2]=(pars->MM);
	mME[3]=(pars->MM);
	mME[4]=(pars->MM);
	mME[5]=(pars->MM);
	mME[6]=(pars->ZERO);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void gg_Spin0_2f_SFpA::sigmaKin() 
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
  const int denominators[nprocesses] = {1024};

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
        t[0] = matrix_gg_h_mummupmummupa_no_zpxg();

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
      t[0] = matrix_gg_h_mummupmummupa_no_zpxg();

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

double gg_Spin0_2f_SFpA::sigmaHat() 
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

void gg_Spin0_2f_SFpA::calculate_wavefunctions(const int perm[], const int hel[])
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
  FFV2P0_3(w[5], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[29]);
  FFS1_2_2(w[3], w[7], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[30]);
  FFV2_1(w[4], w[29], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[31]);
  FFV2_2(w[3], w[29], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[32]);
  FFV5_7_3(w[5], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[33]);
  FFV5_7_1(w[4], w[33], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[34]);
  FFV5_7_2(w[3], w[33], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[35]);
  FFS1_2_2(w[3], w[16], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[36]);
  FFV2_2(w[3], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[37]);
  FFS1_2_2(w[37], w[7], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[38]);
//   VVS3_4_5_1(w[33], w[7], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22, pars->MZ,
//       pars->WZ, w[39]);
//   VVS2_1(w[33], w[7], pars->HEF_MEKD2_1_GC_25, pars->MZ, pars->WZ, w[40]);
  FFS1_2_2(w[37], w[16], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[41]);
//   VVS3_4_5_1(w[33], w[16], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22, pars->MZ,
//       pars->WZ, w[42]);
//   VVS2_1(w[33], w[16], pars->HEF_MEKD2_1_GC_25, pars->MZ, pars->WZ, w[43]);
  FFV2_1(w[2], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[44]);
  FFV2P0_3(w[5], w[44], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[45]);
  FFV5_7_3(w[5], w[44], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[46]);
  FFV2P0_3(w[3], w[44], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[47]);
  FFV5_7_3(w[3], w[44], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[48]);
  FFV2P0_3(w[3], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[49]);
  FFS1_2_1(w[44], w[7], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[50]);
  FFV5_7_3(w[3], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[51]);
//   VVS3_4_5_1(w[51], w[7], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22, pars->MZ,
//       pars->WZ, w[52]);
//   VVS2_1(w[51], w[7], pars->HEF_MEKD2_1_GC_25, pars->MZ, pars->WZ, w[53]);
  FFS1_2_1(w[44], w[16], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[54]);
//   VVS3_4_5_1(w[51], w[16], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22, pars->MZ,
//       pars->WZ, w[55]);
//   VVS2_1(w[51], w[16], pars->HEF_MEKD2_1_GC_25, pars->MZ, pars->WZ, w[56]);
  FFV2P0_3(w[5], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[57]);
  FFV5_7_3(w[5], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[58]);
//   VVS3_4_5_1(w[58], w[7], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22, pars->MZ,
//       pars->WZ, w[59]);
//   VVS2_1(w[58], w[7], pars->HEF_MEKD2_1_GC_25, pars->MZ, pars->WZ, w[60]);
//   VVS3_4_5_1(w[58], w[16], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22, pars->MZ,
//       pars->WZ, w[61]);
//   VVS2_1(w[58], w[16], pars->HEF_MEKD2_1_GC_25, pars->MZ, pars->WZ, w[62]);
  FFS1_2_1(w[2], w[7], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[63]);
  FFV2_2(w[5], w[49], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[64]);
  FFV2_1(w[2], w[49], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[65]);
  FFV5_7_2(w[5], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[66]);
  FFV5_7_1(w[2], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[67]);
  FFS1_2_1(w[2], w[16], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[68]);
  FFV2P0_3(w[37], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[69]);
  FFV5_7_3(w[37], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[70]);
  FFV2P0_3(w[37], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[71]);
  FFV5_7_3(w[37], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[72]);
  FFV2_2(w[3], w[57], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[73]);
  FFV2_1(w[2], w[57], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[74]);
  FFV5_7_2(w[3], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[75]);
  FFV5_7_1(w[2], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[76]);
  FFV2P0_3(w[3], w[19], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[77]);
  FFV5_7_3(w[3], w[19], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[78]);
  FFV2P0_3(w[5], w[19], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[79]);
  FFV5_7_3(w[5], w[19], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[80]);
  FFV2P0_3(w[26], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[81]);
  FFV5_7_3(w[26], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[82]);
  FFV2P0_3(w[26], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[83]);
  FFV5_7_3(w[26], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[84]);

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
  FFV2_0(w[30], w[31], w[6], pars->HEF_MEKD2_1_GC_5, amp[32]);
  FFV2_0(w[32], w[9], w[6], pars->HEF_MEKD2_1_GC_5, amp[33]);
  FFV2_0(w[30], w[34], w[6], pars->HEF_MEKD2_1_GC_5, amp[34]);
  FFV2_0(w[35], w[9], w[6], pars->HEF_MEKD2_1_GC_5, amp[35]);
  FFV2_0(w[36], w[31], w[6], pars->HEF_MEKD2_1_GC_5, amp[36]);
  FFV2_0(w[32], w[17], w[6], pars->HEF_MEKD2_1_GC_5, amp[37]);
  FFV2_0(w[36], w[34], w[6], pars->HEF_MEKD2_1_GC_5, amp[38]);
  FFV2_0(w[35], w[17], w[6], pars->HEF_MEKD2_1_GC_5, amp[39]);
  FFV2_0(w[38], w[4], w[29], pars->HEF_MEKD2_1_GC_5, amp[40]);
  FFV2_0(w[37], w[9], w[29], pars->HEF_MEKD2_1_GC_5, amp[41]);
//   FFV5_7_0(w[37], w[4], w[39], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[42]);	// ZZ
//   FFV5_7_0(w[37], w[4], w[40], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[43]);	// ZZ
  FFV5_7_0(w[38], w[4], w[33], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[44]);
  FFV5_7_0(w[37], w[9], w[33], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[45]);
  FFV2_0(w[41], w[4], w[29], pars->HEF_MEKD2_1_GC_5, amp[46]);
  FFV2_0(w[37], w[17], w[29], pars->HEF_MEKD2_1_GC_5, amp[47]);
//   FFV5_7_0(w[37], w[4], w[42], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[48]);	// ZZ
//   FFV5_7_0(w[37], w[4], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[49]);	// ZZ
  FFV5_7_0(w[41], w[4], w[33], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[50]);
  FFV5_7_0(w[37], w[17], w[33], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[51]);
  FFV2_0(w[30], w[19], w[29], pars->HEF_MEKD2_1_GC_5, amp[52]);
  FFV2_0(w[3], w[20], w[29], pars->HEF_MEKD2_1_GC_5, amp[53]);
//   FFV5_7_0(w[3], w[19], w[39], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[54]);	// ZZ
//   FFV5_7_0(w[3], w[19], w[40], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[55]);	// ZZ
  FFV5_7_0(w[30], w[19], w[33], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[56]);
  FFV5_7_0(w[3], w[20], w[33], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[57]);
  FFV2_0(w[36], w[19], w[29], pars->HEF_MEKD2_1_GC_5, amp[58]);
  FFV2_0(w[3], w[23], w[29], pars->HEF_MEKD2_1_GC_5, amp[59]);
//   FFV5_7_0(w[3], w[19], w[42], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[60]);	// ZZ
//   FFV5_7_0(w[3], w[19], w[43], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[61]);	// ZZ
  FFV5_7_0(w[36], w[19], w[33], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[62]);
//   FFV5_7_0(w[3], w[23], w[33], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[63]);	// ZZ
  FFV2_0(w[30], w[4], w[45], pars->HEF_MEKD2_1_GC_5, amp[64]);	// ZZ
  FFV5_7_0(w[30], w[4], w[46], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[65]);
  FFV2_0(w[5], w[9], w[47], pars->HEF_MEKD2_1_GC_5, amp[66]);
  FFV5_7_0(w[5], w[9], w[48], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[67]);
  FFV2_0(w[3], w[9], w[45], pars->HEF_MEKD2_1_GC_5, amp[68]);
  FFV5_7_0(w[3], w[9], w[46], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[69]);
  FFV2_0(w[11], w[4], w[47], pars->HEF_MEKD2_1_GC_5, amp[70]);
  FFV5_7_0(w[11], w[4], w[48], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[71]);
  FFV2_0(w[36], w[4], w[45], pars->HEF_MEKD2_1_GC_5, amp[72]);
  FFV5_7_0(w[36], w[4], w[46], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[73]);
  FFV2_0(w[5], w[17], w[47], pars->HEF_MEKD2_1_GC_5, amp[74]);
  FFV5_7_0(w[5], w[17], w[48], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[75]);
  FFV2_0(w[3], w[17], w[45], pars->HEF_MEKD2_1_GC_5, amp[76]);
  FFV5_7_0(w[3], w[17], w[46], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[77]);
  FFV2_0(w[18], w[4], w[47], pars->HEF_MEKD2_1_GC_5, amp[78]);
  FFV5_7_0(w[18], w[4], w[48], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[79]);
  FFV2_0(w[5], w[50], w[49], pars->HEF_MEKD2_1_GC_5, amp[80]);
  FFV2_0(w[11], w[44], w[49], pars->HEF_MEKD2_1_GC_5, amp[81]);
  FFV5_7_0(w[5], w[50], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[82]);
//   FFV5_7_0(w[5], w[44], w[52], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[83]);	// ZZ
//   FFV5_7_0(w[5], w[44], w[53], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[84]);	// ZZ
  FFV5_7_0(w[11], w[44], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[85]);
  FFV2_0(w[5], w[54], w[49], pars->HEF_MEKD2_1_GC_5, amp[86]);
  FFV2_0(w[18], w[44], w[49], pars->HEF_MEKD2_1_GC_5, amp[87]);
  FFV5_7_0(w[5], w[54], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[88]);
//   FFV5_7_0(w[5], w[44], w[55], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[89]);	// ZZ
//   FFV5_7_0(w[5], w[44], w[56], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[90]);	// ZZ
  FFV5_7_0(w[18], w[44], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[91]);
  FFV2_0(w[3], w[50], w[57], pars->HEF_MEKD2_1_GC_5, amp[92]);
  FFV2_0(w[30], w[44], w[57], pars->HEF_MEKD2_1_GC_5, amp[93]);
  FFV5_7_0(w[3], w[50], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[94]);
  FFV5_7_0(w[30], w[44], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[95]);
//   FFV5_7_0(w[3], w[44], w[59], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[96]);	// ZZ
//   FFV5_7_0(w[3], w[44], w[60], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[97]);	// ZZ
  FFV2_0(w[3], w[54], w[57], pars->HEF_MEKD2_1_GC_5, amp[98]);
  FFV2_0(w[36], w[44], w[57], pars->HEF_MEKD2_1_GC_5, amp[99]);
  FFV5_7_0(w[3], w[54], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[100]);
  FFV5_7_0(w[36], w[44], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[101]);
//   FFV5_7_0(w[3], w[44], w[61], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[102]);	// ZZ
//   FFV5_7_0(w[3], w[44], w[62], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[103]);	// ZZ
  FFV2_0(w[64], w[63], w[6], pars->HEF_MEKD2_1_GC_5, amp[104]);
  FFV2_0(w[11], w[65], w[6], pars->HEF_MEKD2_1_GC_5, amp[105]);
  FFV2_0(w[66], w[63], w[6], pars->HEF_MEKD2_1_GC_5, amp[106]);
  FFV2_0(w[11], w[67], w[6], pars->HEF_MEKD2_1_GC_5, amp[107]);
  FFV2_0(w[64], w[68], w[6], pars->HEF_MEKD2_1_GC_5, amp[108]);
  FFV2_0(w[18], w[65], w[6], pars->HEF_MEKD2_1_GC_5, amp[109]);
  FFV2_0(w[66], w[68], w[6], pars->HEF_MEKD2_1_GC_5, amp[110]);
  FFV2_0(w[18], w[67], w[6], pars->HEF_MEKD2_1_GC_5, amp[111]);
  FFV2_0(w[26], w[63], w[49], pars->HEF_MEKD2_1_GC_5, amp[112]);
  FFV2_0(w[27], w[2], w[49], pars->HEF_MEKD2_1_GC_5, amp[113]);
  FFV5_7_0(w[26], w[63], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[114]);
//   FFV5_7_0(w[26], w[2], w[52], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[115]);	// ZZ
//   FFV5_7_0(w[26], w[2], w[53], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[116]);	// ZZ
  FFV5_7_0(w[27], w[2], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[117]);
  FFV2_0(w[26], w[68], w[49], pars->HEF_MEKD2_1_GC_5, amp[118]);
  FFV2_0(w[28], w[2], w[49], pars->HEF_MEKD2_1_GC_5, amp[119]);
  FFV5_7_0(w[26], w[68], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[120]);
//   FFV5_7_0(w[26], w[2], w[55], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[121]);	// ZZ
//   FFV5_7_0(w[26], w[2], w[56], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[122]);	// ZZ
  FFV5_7_0(w[28], w[2], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[123]);
  FFV2_0(w[5], w[63], w[69], pars->HEF_MEKD2_1_GC_5, amp[124]);
  FFV5_7_0(w[5], w[63], w[70], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[125]);
  FFV2_0(w[5], w[9], w[71], pars->HEF_MEKD2_1_GC_5, amp[126]);
  FFV5_7_0(w[5], w[9], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[127]);
  FFV2_0(w[11], w[4], w[71], pars->HEF_MEKD2_1_GC_5, amp[128]);
  FFV5_7_0(w[11], w[4], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[129]);
  FFV2_0(w[11], w[2], w[69], pars->HEF_MEKD2_1_GC_5, amp[130]);
  FFV5_7_0(w[11], w[2], w[70], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[131]);
  FFV2_0(w[5], w[68], w[69], pars->HEF_MEKD2_1_GC_5, amp[132]);
  FFV5_7_0(w[5], w[68], w[70], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[133]);
  FFV2_0(w[5], w[17], w[71], pars->HEF_MEKD2_1_GC_5, amp[134]);
  FFV5_7_0(w[5], w[17], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[135]);
  FFV2_0(w[18], w[4], w[71], pars->HEF_MEKD2_1_GC_5, amp[136]);
  FFV5_7_0(w[18], w[4], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[137]);
  FFV2_0(w[18], w[2], w[69], pars->HEF_MEKD2_1_GC_5, amp[138]);
  FFV5_7_0(w[18], w[2], w[70], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[139]);
  FFV2_0(w[37], w[63], w[57], pars->HEF_MEKD2_1_GC_5, amp[140]);
  FFV2_0(w[38], w[2], w[57], pars->HEF_MEKD2_1_GC_5, amp[141]);
  FFV5_7_0(w[37], w[63], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[142]);
  FFV5_7_0(w[38], w[2], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[143]);
//   FFV5_7_0(w[37], w[2], w[59], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[144]);	// ZZ
//   FFV5_7_0(w[37], w[2], w[60], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[145]);	// ZZ
  FFV2_0(w[37], w[68], w[57], pars->HEF_MEKD2_1_GC_5, amp[146]);
  FFV2_0(w[41], w[2], w[57], pars->HEF_MEKD2_1_GC_5, amp[147]);
  FFV5_7_0(w[37], w[68], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[148]);
  FFV5_7_0(w[41], w[2], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[149]);
//   FFV5_7_0(w[37], w[2], w[61], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[150]);	// ZZ
//   FFV5_7_0(w[37], w[2], w[62], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[151]);	// ZZ
  FFV2_0(w[73], w[63], w[6], pars->HEF_MEKD2_1_GC_5, amp[152]);
  FFV2_0(w[30], w[74], w[6], pars->HEF_MEKD2_1_GC_5, amp[153]);
  FFV2_0(w[75], w[63], w[6], pars->HEF_MEKD2_1_GC_5, amp[154]);
  FFV2_0(w[30], w[76], w[6], pars->HEF_MEKD2_1_GC_5, amp[155]);
  FFV2_0(w[73], w[68], w[6], pars->HEF_MEKD2_1_GC_5, amp[156]);
  FFV2_0(w[36], w[74], w[6], pars->HEF_MEKD2_1_GC_5, amp[157]);
  FFV2_0(w[75], w[68], w[6], pars->HEF_MEKD2_1_GC_5, amp[158]);
  FFV2_0(w[36], w[76], w[6], pars->HEF_MEKD2_1_GC_5, amp[159]);
  FFV2_0(w[5], w[63], w[77], pars->HEF_MEKD2_1_GC_5, amp[160]);
  FFV5_7_0(w[5], w[63], w[78], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[161]);
  FFV2_0(w[3], w[63], w[79], pars->HEF_MEKD2_1_GC_5, amp[162]);
  FFV5_7_0(w[3], w[63], w[80], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[163]);
  FFV2_0(w[30], w[2], w[79], pars->HEF_MEKD2_1_GC_5, amp[164]);
  FFV5_7_0(w[30], w[2], w[80], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[165]);
  FFV2_0(w[11], w[2], w[77], pars->HEF_MEKD2_1_GC_5, amp[166]);
  FFV5_7_0(w[11], w[2], w[78], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[167]);
  FFV2_0(w[5], w[68], w[77], pars->HEF_MEKD2_1_GC_5, amp[168]);
  FFV5_7_0(w[5], w[68], w[78], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[169]);
  FFV2_0(w[3], w[68], w[79], pars->HEF_MEKD2_1_GC_5, amp[170]);
  FFV5_7_0(w[3], w[68], w[80], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[171]);
  FFV2_0(w[36], w[2], w[79], pars->HEF_MEKD2_1_GC_5, amp[172]);
  FFV5_7_0(w[36], w[2], w[80], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[173]);
  FFV2_0(w[18], w[2], w[77], pars->HEF_MEKD2_1_GC_5, amp[174]);
  FFV5_7_0(w[18], w[2], w[78], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[175]);
  FFV2_0(w[3], w[63], w[81], pars->HEF_MEKD2_1_GC_5, amp[176]);
  FFV5_7_0(w[3], w[63], w[82], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[177]);
  FFV2_0(w[30], w[4], w[83], pars->HEF_MEKD2_1_GC_5, amp[178]);
  FFV5_7_0(w[30], w[4], w[84], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[179]);
  FFV2_0(w[30], w[2], w[81], pars->HEF_MEKD2_1_GC_5, amp[180]);
  FFV5_7_0(w[30], w[2], w[82], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[181]);
  FFV2_0(w[3], w[9], w[83], pars->HEF_MEKD2_1_GC_5, amp[182]);
  FFV5_7_0(w[3], w[9], w[84], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[183]);
  FFV2_0(w[3], w[68], w[81], pars->HEF_MEKD2_1_GC_5, amp[184]);
  FFV5_7_0(w[3], w[68], w[82], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[185]);
  FFV2_0(w[36], w[4], w[83], pars->HEF_MEKD2_1_GC_5, amp[186]);
  FFV5_7_0(w[36], w[4], w[84], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[187]);
  FFV2_0(w[36], w[2], w[81], pars->HEF_MEKD2_1_GC_5, amp[188]);
  FFV5_7_0(w[36], w[2], w[82], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[189]);
  FFV2_0(w[3], w[17], w[83], pars->HEF_MEKD2_1_GC_5, amp[190]);
  FFV5_7_0(w[3], w[17], w[84], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[191]);

}

double gg_Spin0_2f_SFpA::matrix_gg_h_mummupmummupa_no_zpxg() 
{
  int i, j;
  // Local variables
	// Commented out by Convert_source 0.2
//  const int ngraphs = 192;
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
      amp[28] + amp[29] + amp[30] + amp[31] - amp[32] - amp[33] - amp[34] -
      amp[35] - amp[36] - amp[37] - amp[38] - amp[39] - amp[40] - amp[41] -
      amp[42] - amp[43] - amp[44] - amp[45] - amp[46] - amp[47] - amp[48] -
      amp[49] - amp[50] - amp[51] - amp[52] - amp[53] - amp[54] - amp[55] -
      amp[56] - amp[57] - amp[58] - amp[59] - amp[60] - amp[61] - amp[62] -
      amp[63] - amp[64] - amp[65] + amp[66] + amp[67] - amp[68] - amp[69] +
      amp[70] + amp[71] - amp[72] - amp[73] + amp[74] + amp[75] - amp[76] -
      amp[77] + amp[78] + amp[79] - amp[80] - amp[81] - amp[82] - amp[83] -
      amp[84] - amp[85] - amp[86] - amp[87] - amp[88] - amp[89] - amp[90] -
      amp[91] + amp[92] + amp[93] + amp[94] + amp[95] + amp[96] + amp[97] +
      amp[98] + amp[99] + amp[100] + amp[101] + amp[102] + amp[103] - amp[104]
      - amp[105] - amp[106] - amp[107] - amp[108] - amp[109] - amp[110] -
      amp[111] - amp[112] - amp[113] - amp[114] - amp[115] - amp[116] -
      amp[117] - amp[118] - amp[119] - amp[120] - amp[121] - amp[122] -
      amp[123] - amp[124] - amp[125] + amp[126] + amp[127] + amp[128] +
      amp[129] - amp[130] - amp[131] - amp[132] - amp[133] + amp[134] +
      amp[135] + amp[136] + amp[137] - amp[138] - amp[139] + amp[140] +
      amp[141] + amp[142] + amp[143] + amp[144] + amp[145] + amp[146] +
      amp[147] + amp[148] + amp[149] + amp[150] + amp[151] + amp[152] +
      amp[153] + amp[154] + amp[155] + amp[156] + amp[157] + amp[158] +
      amp[159] - amp[160] - amp[161] + amp[162] + amp[163] + amp[164] +
      amp[165] - amp[166] - amp[167] - amp[168] - amp[169] + amp[170] +
      amp[171] + amp[172] + amp[173] - amp[174] - amp[175] + amp[176] +
      amp[177] - amp[178] - amp[179] + amp[180] + amp[181] - amp[182] -
      amp[183] + amp[184] + amp[185] - amp[186] - amp[187] + amp[188] +
      amp[189] - amp[190] - amp[191]);

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



