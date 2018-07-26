//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.2, 2014-02-07
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "gg_Spin2_2f_SFpA.h"
#include "../HelAmps_HEF_MEKD2_1.h"	// Changed by Convert_source 0.2

using namespace MG5_HEF_MEKD2_1;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > xg > mu- mu+ mu- mu+ a GIG=1 QED=3 S2QQ=1 GIZ=1 / zp h

//--------------------------------------------------------------------------
// Initialize process.

void gg_Spin2_2f_SFpA::initProc(string param_card_name) 
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

void gg_Spin2_2f_SFpA::updateProc(SLHAReader_MEKD &slha)
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

void gg_Spin2_2f_SFpA::sigmaKin() 
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
        t[0] = matrix_gg_xg_mummupmummupa_no_zph();

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
      t[0] = matrix_gg_xg_mummupmummupa_no_zph();

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

double gg_Spin2_2f_SFpA::sigmaHat() 
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

void gg_Spin2_2f_SFpA::calculate_wavefunctions(const int perm[], const int hel[])
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
  VVT1_10_11_12_13_3_5_7_8_9_3(w[0], w[1], pars->HEF_MEKD2_1_GC_86, pars->HEF_MEKD2_1_GC_72,
      pars->HEF_MEKD2_1_GC_80, pars->HEF_MEKD2_1_GC_64, pars->HEF_MEKD2_1_GC_68, pars->HEF_MEKD2_1_GC_91, pars->HEF_MEKD2_1_GC_62,
      pars->HEF_MEKD2_1_GC_84, pars->HEF_MEKD2_1_GC_76, pars->HEF_MEKD2_1_GC_82, pars->MXG, pars->WXG, w[7]);
  FFV2P0_3(w[3], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[8]);
  FFT1_2_3_5_1(w[4], w[7], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[9]);
  FFV2_2(w[5], w[8], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[10]);
  FFT1_2_3_5_2(w[5], w[7], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[11]);
  FFV2_1(w[4], w[8], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[12]);
  FFV5_7_3(w[3], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[13]);
  FFV5_7_2(w[5], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[14]);
  FFV5_7_1(w[4], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[15]);
  VVT4_3(w[0], w[1], pars->HEF_MEKD2_1_GC_63, pars->MXG, pars->WXG, w[16]);
  FFT1_2_3_5_1(w[4], w[16], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[17]);
  FFT1_2_3_5_2(w[5], w[16], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[18]);
  FFV2_1(w[4], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[19]);
  FFT1_2_3_5_1(w[19], w[7], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[20]);
//   VVT10_11_12_13_2_3_6_7_8_9_1(w[13], w[7], pars->HEF_MEKD2_1_GC_75, pars->HEF_MEKD2_1_GC_81,
//       pars->HEF_MEKD2_1_GC_67, pars->HEF_MEKD2_1_GC_71, pars->HEF_MEKD2_1_GC_90, pars->HEF_MEKD2_1_GC_92, pars->HEF_MEKD2_1_GC_63,
//       pars->HEF_MEKD2_1_GC_85, pars->HEF_MEKD2_1_GC_79, pars->HEF_MEKD2_1_GC_83, pars->MZ, pars->WZ, w[21]);
  FFT1_2_3_5_1(w[19], w[16], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[22]);
  FFV2_2(w[5], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[23]);
  FFT1_2_3_5_2(w[23], w[7], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[24]);
  FFT1_2_3_5_2(w[23], w[16], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[25]);
  FFV2P0_3(w[5], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[26]);
  FFT1_2_3_5_2(w[3], w[7], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[27]);
  FFV2_1(w[4], w[26], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[28]);
  FFV2_2(w[3], w[26], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[29]);
  FFV5_7_3(w[5], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[30]);
  FFV5_7_1(w[4], w[30], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[31]);
  FFV5_7_2(w[3], w[30], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[32]);
  FFT1_2_3_5_2(w[3], w[16], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[33]);
  FFV2_2(w[3], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[34]);
  FFT1_2_3_5_2(w[34], w[7], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[35]);
//   VVT10_11_12_13_2_3_6_7_8_9_1(w[30], w[7], pars->HEF_MEKD2_1_GC_75, pars->HEF_MEKD2_1_GC_81,
//       pars->HEF_MEKD2_1_GC_67, pars->HEF_MEKD2_1_GC_71, pars->HEF_MEKD2_1_GC_90, pars->HEF_MEKD2_1_GC_92, pars->HEF_MEKD2_1_GC_63,
//       pars->HEF_MEKD2_1_GC_85, pars->HEF_MEKD2_1_GC_79, pars->HEF_MEKD2_1_GC_83, pars->MZ, pars->WZ, w[36]);
  FFT1_2_3_5_2(w[34], w[16], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[37]);
  FFV2_1(w[2], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[38]);
  FFV2P0_3(w[5], w[38], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[39]);
  FFV5_7_3(w[5], w[38], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[40]);
  FFV2P0_3(w[3], w[38], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[41]);
  FFV5_7_3(w[3], w[38], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[42]);
  FFV2P0_3(w[3], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[43]);
  FFT1_2_3_5_1(w[38], w[7], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[44]);
  FFV5_7_3(w[3], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[45]);
  VVT10_11_12_13_2_3_6_7_8_9_1(w[45], w[7], pars->HEF_MEKD2_1_GC_75, pars->HEF_MEKD2_1_GC_81,
      pars->HEF_MEKD2_1_GC_67, pars->HEF_MEKD2_1_GC_71, pars->HEF_MEKD2_1_GC_90, pars->HEF_MEKD2_1_GC_92, pars->HEF_MEKD2_1_GC_63,
      pars->HEF_MEKD2_1_GC_85, pars->HEF_MEKD2_1_GC_79, pars->HEF_MEKD2_1_GC_83, pars->MZ, pars->WZ, w[46]);
  FFT1_2_3_5_1(w[38], w[16], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[47]);
  FFV2P0_3(w[5], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[48]);
  FFV5_7_3(w[5], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[49]);
//   VVT10_11_12_13_2_3_6_7_8_9_1(w[49], w[7], pars->HEF_MEKD2_1_GC_75, pars->HEF_MEKD2_1_GC_81,
//       pars->HEF_MEKD2_1_GC_67, pars->HEF_MEKD2_1_GC_71, pars->HEF_MEKD2_1_GC_90, pars->HEF_MEKD2_1_GC_92, pars->HEF_MEKD2_1_GC_63,
//       pars->HEF_MEKD2_1_GC_85, pars->HEF_MEKD2_1_GC_79, pars->HEF_MEKD2_1_GC_83, pars->MZ, pars->WZ, w[50]);
  FFT1_2_3_5_1(w[2], w[7], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[51]);
  FFV2_2(w[5], w[43], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[52]);
  FFV2_1(w[2], w[43], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[53]);
  FFV5_7_2(w[5], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[54]);
  FFV5_7_1(w[2], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[55]);
  FFT1_2_3_5_1(w[2], w[16], pars->HEF_MEKD2_1_GC_145, pars->HEF_MEKD2_1_GC_142, pars->HEF_MEKD2_1_GC_143,
      pars->HEF_MEKD2_1_GC_144, pars->MM, pars->ZERO, w[56]);
  FFV2P0_3(w[34], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[57]);
  FFV5_7_3(w[34], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[58]);
  FFV2P0_3(w[34], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[59]);
  FFV5_7_3(w[34], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[60]);
  FFV2_2(w[3], w[48], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[61]);
  FFV2_1(w[2], w[48], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[62]);
  FFV5_7_2(w[3], w[49], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[63]);
  FFV5_7_1(w[2], w[49], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[64]);
  FFV2P0_3(w[3], w[19], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[65]);
  FFV5_7_3(w[3], w[19], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[66]);
  FFV2P0_3(w[5], w[19], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[67]);
  FFV5_7_3(w[5], w[19], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[68]);
  FFV2P0_3(w[23], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[69]);
  FFV5_7_3(w[23], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[70]);
  FFV2P0_3(w[23], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[71]);
  FFV5_7_3(w[23], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[72]);

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
  FFV5_7_0(w[5], w[20], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[11]);
  FFV5_7_0(w[11], w[19], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[12]);
  FFV2_0(w[5], w[22], w[8], pars->HEF_MEKD2_1_GC_5, amp[13]);
  FFV2_0(w[18], w[19], w[8], pars->HEF_MEKD2_1_GC_5, amp[14]);
  FFV5_7_0(w[5], w[22], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[15]);
  FFV5_7_0(w[18], w[19], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[16]);
  FFV2_0(w[23], w[9], w[8], pars->HEF_MEKD2_1_GC_5, amp[17]);
  FFV2_0(w[24], w[4], w[8], pars->HEF_MEKD2_1_GC_5, amp[18]);
//   FFV5_7_0(w[23], w[4], w[21], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[19]);	// ZZ
  FFV5_7_0(w[23], w[9], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[20]);
  FFV5_7_0(w[24], w[4], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[21]);
  FFV2_0(w[23], w[17], w[8], pars->HEF_MEKD2_1_GC_5, amp[22]);
  FFV2_0(w[25], w[4], w[8], pars->HEF_MEKD2_1_GC_5, amp[23]);
  FFV5_7_0(w[23], w[17], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[24]);
  FFV5_7_0(w[25], w[4], w[13], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[25]);
  FFV2_0(w[27], w[28], w[6], pars->HEF_MEKD2_1_GC_5, amp[26]);
  FFV2_0(w[29], w[9], w[6], pars->HEF_MEKD2_1_GC_5, amp[27]);
  FFV2_0(w[27], w[31], w[6], pars->HEF_MEKD2_1_GC_5, amp[28]);
  FFV2_0(w[32], w[9], w[6], pars->HEF_MEKD2_1_GC_5, amp[29]);
  FFV2_0(w[33], w[28], w[6], pars->HEF_MEKD2_1_GC_5, amp[30]);
  FFV2_0(w[29], w[17], w[6], pars->HEF_MEKD2_1_GC_5, amp[31]);
  FFV2_0(w[33], w[31], w[6], pars->HEF_MEKD2_1_GC_5, amp[32]);
  FFV2_0(w[32], w[17], w[6], pars->HEF_MEKD2_1_GC_5, amp[33]);
  FFV2_0(w[35], w[4], w[26], pars->HEF_MEKD2_1_GC_5, amp[34]);
  FFV2_0(w[34], w[9], w[26], pars->HEF_MEKD2_1_GC_5, amp[35]);
//   FFV5_7_0(w[34], w[4], w[36], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[36]);	// ZZ
  FFV5_7_0(w[35], w[4], w[30], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[37]);
  FFV5_7_0(w[34], w[9], w[30], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[38]);
  FFV2_0(w[37], w[4], w[26], pars->HEF_MEKD2_1_GC_5, amp[39]);
  FFV2_0(w[34], w[17], w[26], pars->HEF_MEKD2_1_GC_5, amp[40]);
  FFV5_7_0(w[37], w[4], w[30], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[41]);
  FFV5_7_0(w[34], w[17], w[30], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[42]);
  FFV2_0(w[27], w[19], w[26], pars->HEF_MEKD2_1_GC_5, amp[43]);
  FFV2_0(w[3], w[20], w[26], pars->HEF_MEKD2_1_GC_5, amp[44]);
//   FFV5_7_0(w[3], w[19], w[36], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[45]);	// ZZ
  FFV5_7_0(w[27], w[19], w[30], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[46]);
  FFV5_7_0(w[3], w[20], w[30], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[47]);
  FFV2_0(w[33], w[19], w[26], pars->HEF_MEKD2_1_GC_5, amp[48]);
  FFV2_0(w[3], w[22], w[26], pars->HEF_MEKD2_1_GC_5, amp[49]);
  FFV5_7_0(w[33], w[19], w[30], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[50]);
  FFV5_7_0(w[3], w[22], w[30], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[51]);
  FFV2_0(w[27], w[4], w[39], pars->HEF_MEKD2_1_GC_5, amp[52]);
  FFV5_7_0(w[27], w[4], w[40], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[53]);
  FFV2_0(w[5], w[9], w[41], pars->HEF_MEKD2_1_GC_5, amp[54]);
  FFV5_7_0(w[5], w[9], w[42], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[55]);
  FFV2_0(w[3], w[9], w[39], pars->HEF_MEKD2_1_GC_5, amp[56]);
  FFV5_7_0(w[3], w[9], w[40], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[57]);
  FFV2_0(w[11], w[4], w[41], pars->HEF_MEKD2_1_GC_5, amp[58]);
  FFV5_7_0(w[11], w[4], w[42], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[59]);
  FFV2_0(w[33], w[4], w[39], pars->HEF_MEKD2_1_GC_5, amp[60]);
  FFV5_7_0(w[33], w[4], w[40], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[61]);
  FFV2_0(w[5], w[17], w[41], pars->HEF_MEKD2_1_GC_5, amp[62]);
  FFV5_7_0(w[5], w[17], w[42], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[63]);
  FFV2_0(w[3], w[17], w[39], pars->HEF_MEKD2_1_GC_5, amp[64]);
  FFV5_7_0(w[3], w[17], w[40], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[65]);
  FFV2_0(w[18], w[4], w[41], pars->HEF_MEKD2_1_GC_5, amp[66]);
  FFV5_7_0(w[18], w[4], w[42], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[67]);
  FFV2_0(w[5], w[44], w[43], pars->HEF_MEKD2_1_GC_5, amp[68]);
  FFV2_0(w[11], w[38], w[43], pars->HEF_MEKD2_1_GC_5, amp[69]);
  FFV5_7_0(w[5], w[44], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[70]);
//   FFV5_7_0(w[5], w[38], w[46], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[71]);	// ZZ
  FFV5_7_0(w[11], w[38], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[72]);
  FFV2_0(w[5], w[47], w[43], pars->HEF_MEKD2_1_GC_5, amp[73]);
  FFV2_0(w[18], w[38], w[43], pars->HEF_MEKD2_1_GC_5, amp[74]);
  FFV5_7_0(w[5], w[47], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[75]);
  FFV5_7_0(w[18], w[38], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[76]);
  FFV2_0(w[3], w[44], w[48], pars->HEF_MEKD2_1_GC_5, amp[77]);
  FFV2_0(w[27], w[38], w[48], pars->HEF_MEKD2_1_GC_5, amp[78]);
  FFV5_7_0(w[3], w[44], w[49], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[79]);
  FFV5_7_0(w[27], w[38], w[49], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[80]);
//   FFV5_7_0(w[3], w[38], w[50], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[81]);	// ZZ
  FFV2_0(w[3], w[47], w[48], pars->HEF_MEKD2_1_GC_5, amp[82]);
  FFV2_0(w[33], w[38], w[48], pars->HEF_MEKD2_1_GC_5, amp[83]);
  FFV5_7_0(w[3], w[47], w[49], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[84]);
  FFV5_7_0(w[33], w[38], w[49], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[85]);
  FFV2_0(w[52], w[51], w[6], pars->HEF_MEKD2_1_GC_5, amp[86]);
  FFV2_0(w[11], w[53], w[6], pars->HEF_MEKD2_1_GC_5, amp[87]);
  FFV2_0(w[54], w[51], w[6], pars->HEF_MEKD2_1_GC_5, amp[88]);
  FFV2_0(w[11], w[55], w[6], pars->HEF_MEKD2_1_GC_5, amp[89]);
  FFV2_0(w[52], w[56], w[6], pars->HEF_MEKD2_1_GC_5, amp[90]);
  FFV2_0(w[18], w[53], w[6], pars->HEF_MEKD2_1_GC_5, amp[91]);
  FFV2_0(w[54], w[56], w[6], pars->HEF_MEKD2_1_GC_5, amp[92]);
  FFV2_0(w[18], w[55], w[6], pars->HEF_MEKD2_1_GC_5, amp[93]);
  FFV2_0(w[23], w[51], w[43], pars->HEF_MEKD2_1_GC_5, amp[94]);
  FFV2_0(w[24], w[2], w[43], pars->HEF_MEKD2_1_GC_5, amp[95]);
  FFV5_7_0(w[23], w[51], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[96]);
//   FFV5_7_0(w[23], w[2], w[46], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[97]);	// ZZ
  FFV5_7_0(w[24], w[2], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[98]);
  FFV2_0(w[23], w[56], w[43], pars->HEF_MEKD2_1_GC_5, amp[99]);
  FFV2_0(w[25], w[2], w[43], pars->HEF_MEKD2_1_GC_5, amp[100]);
  FFV5_7_0(w[23], w[56], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[101]);
  FFV5_7_0(w[25], w[2], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[102]);
  FFV2_0(w[5], w[51], w[57], pars->HEF_MEKD2_1_GC_5, amp[103]);
  FFV5_7_0(w[5], w[51], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[104]);
  FFV2_0(w[5], w[9], w[59], pars->HEF_MEKD2_1_GC_5, amp[105]);
  FFV5_7_0(w[5], w[9], w[60], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[106]);
  FFV2_0(w[11], w[4], w[59], pars->HEF_MEKD2_1_GC_5, amp[107]);
  FFV5_7_0(w[11], w[4], w[60], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[108]);
  FFV2_0(w[11], w[2], w[57], pars->HEF_MEKD2_1_GC_5, amp[109]);
  FFV5_7_0(w[11], w[2], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[110]);
  FFV2_0(w[5], w[56], w[57], pars->HEF_MEKD2_1_GC_5, amp[111]);
  FFV5_7_0(w[5], w[56], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[112]);
  FFV2_0(w[5], w[17], w[59], pars->HEF_MEKD2_1_GC_5, amp[113]);
  FFV5_7_0(w[5], w[17], w[60], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[114]);
  FFV2_0(w[18], w[4], w[59], pars->HEF_MEKD2_1_GC_5, amp[115]);
  FFV5_7_0(w[18], w[4], w[60], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[116]);
  FFV2_0(w[18], w[2], w[57], pars->HEF_MEKD2_1_GC_5, amp[117]);
  FFV5_7_0(w[18], w[2], w[58], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[118]);
  FFV2_0(w[34], w[51], w[48], pars->HEF_MEKD2_1_GC_5, amp[119]);
  FFV2_0(w[35], w[2], w[48], pars->HEF_MEKD2_1_GC_5, amp[120]);
  FFV5_7_0(w[34], w[51], w[49], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[121]);
  FFV5_7_0(w[35], w[2], w[49], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[122]);
//   FFV5_7_0(w[34], w[2], w[50], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[123]);	// ZZ
  FFV2_0(w[34], w[56], w[48], pars->HEF_MEKD2_1_GC_5, amp[124]);
  FFV2_0(w[37], w[2], w[48], pars->HEF_MEKD2_1_GC_5, amp[125]);
  FFV5_7_0(w[34], w[56], w[49], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[126]);
  FFV5_7_0(w[37], w[2], w[49], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[127]);
  FFV2_0(w[61], w[51], w[6], pars->HEF_MEKD2_1_GC_5, amp[128]);
  FFV2_0(w[27], w[62], w[6], pars->HEF_MEKD2_1_GC_5, amp[129]);
  FFV2_0(w[63], w[51], w[6], pars->HEF_MEKD2_1_GC_5, amp[130]);
  FFV2_0(w[27], w[64], w[6], pars->HEF_MEKD2_1_GC_5, amp[131]);
  FFV2_0(w[61], w[56], w[6], pars->HEF_MEKD2_1_GC_5, amp[132]);
  FFV2_0(w[33], w[62], w[6], pars->HEF_MEKD2_1_GC_5, amp[133]);
  FFV2_0(w[63], w[56], w[6], pars->HEF_MEKD2_1_GC_5, amp[134]);
  FFV2_0(w[33], w[64], w[6], pars->HEF_MEKD2_1_GC_5, amp[135]);
  FFV2_0(w[5], w[51], w[65], pars->HEF_MEKD2_1_GC_5, amp[136]);
  FFV5_7_0(w[5], w[51], w[66], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[137]);
  FFV2_0(w[3], w[51], w[67], pars->HEF_MEKD2_1_GC_5, amp[138]);
  FFV5_7_0(w[3], w[51], w[68], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[139]);
  FFV2_0(w[27], w[2], w[67], pars->HEF_MEKD2_1_GC_5, amp[140]);
  FFV5_7_0(w[27], w[2], w[68], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[141]);
  FFV2_0(w[11], w[2], w[65], pars->HEF_MEKD2_1_GC_5, amp[142]);
  FFV5_7_0(w[11], w[2], w[66], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[143]);
  FFV2_0(w[5], w[56], w[65], pars->HEF_MEKD2_1_GC_5, amp[144]);
  FFV5_7_0(w[5], w[56], w[66], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[145]);
  FFV2_0(w[3], w[56], w[67], pars->HEF_MEKD2_1_GC_5, amp[146]);
  FFV5_7_0(w[3], w[56], w[68], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[147]);
  FFV2_0(w[33], w[2], w[67], pars->HEF_MEKD2_1_GC_5, amp[148]);
  FFV5_7_0(w[33], w[2], w[68], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[149]);
  FFV2_0(w[18], w[2], w[65], pars->HEF_MEKD2_1_GC_5, amp[150]);
  FFV5_7_0(w[18], w[2], w[66], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[151]);
  FFV2_0(w[3], w[51], w[69], pars->HEF_MEKD2_1_GC_5, amp[152]);
  FFV5_7_0(w[3], w[51], w[70], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[153]);
  FFV2_0(w[27], w[4], w[71], pars->HEF_MEKD2_1_GC_5, amp[154]);
  FFV5_7_0(w[27], w[4], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[155]);
  FFV2_0(w[27], w[2], w[69], pars->HEF_MEKD2_1_GC_5, amp[156]);
  FFV5_7_0(w[27], w[2], w[70], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[157]);
  FFV2_0(w[3], w[9], w[71], pars->HEF_MEKD2_1_GC_5, amp[158]);
  FFV5_7_0(w[3], w[9], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[159]);
  FFV2_0(w[3], w[56], w[69], pars->HEF_MEKD2_1_GC_5, amp[160]);
  FFV5_7_0(w[3], w[56], w[70], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[161]);
  FFV2_0(w[33], w[4], w[71], pars->HEF_MEKD2_1_GC_5, amp[162]);
  FFV5_7_0(w[33], w[4], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[163]);
  FFV2_0(w[33], w[2], w[69], pars->HEF_MEKD2_1_GC_5, amp[164]);
  FFV5_7_0(w[33], w[2], w[70], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[165]);
  FFV2_0(w[3], w[17], w[71], pars->HEF_MEKD2_1_GC_5, amp[166]);
  FFV5_7_0(w[3], w[17], w[72], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[167]);

}
double gg_Spin2_2f_SFpA::matrix_gg_xg_mummupmummupa_no_zph() 
{
  int i, j;
  // Local variables
	// Commented out by Convert_source 0.2
//  const int ngraphs = 168;
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
      amp[21] + amp[22] + amp[23] + amp[24] + amp[25] - amp[26] - amp[27] -
      amp[28] - amp[29] - amp[30] - amp[31] - amp[32] - amp[33] - amp[34] -
      amp[35] - amp[36] - amp[37] - amp[38] - amp[39] - amp[40] - amp[41] -
      amp[42] - amp[43] - amp[44] - amp[45] - amp[46] - amp[47] - amp[48] -
      amp[49] - amp[50] - amp[51] - amp[52] - amp[53] + amp[54] + amp[55] -
      amp[56] - amp[57] + amp[58] + amp[59] - amp[60] - amp[61] + amp[62] +
      amp[63] - amp[64] - amp[65] + amp[66] + amp[67] - amp[68] - amp[69] -
      amp[70] - amp[71] - amp[72] - amp[73] - amp[74] - amp[75] - amp[76] +
      amp[77] + amp[78] + amp[79] + amp[80] + amp[81] + amp[82] + amp[83] +
      amp[84] + amp[85] - amp[86] - amp[87] - amp[88] - amp[89] - amp[90] -
      amp[91] - amp[92] - amp[93] - amp[94] - amp[95] - amp[96] - amp[97] -
      amp[98] - amp[99] - amp[100] - amp[101] - amp[102] - amp[103] - amp[104]
      + amp[105] + amp[106] + amp[107] + amp[108] - amp[109] - amp[110] -
      amp[111] - amp[112] + amp[113] + amp[114] + amp[115] + amp[116] -
      amp[117] - amp[118] + amp[119] + amp[120] + amp[121] + amp[122] +
      amp[123] + amp[124] + amp[125] + amp[126] + amp[127] + amp[128] +
      amp[129] + amp[130] + amp[131] + amp[132] + amp[133] + amp[134] +
      amp[135] - amp[136] - amp[137] + amp[138] + amp[139] + amp[140] +
      amp[141] - amp[142] - amp[143] - amp[144] - amp[145] + amp[146] +
      amp[147] + amp[148] + amp[149] - amp[150] - amp[151] + amp[152] +
      amp[153] - amp[154] - amp[155] + amp[156] + amp[157] - amp[158] -
      amp[159] + amp[160] + amp[161] - amp[162] - amp[163] + amp[164] +
      amp[165] - amp[166] - amp[167]);

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



