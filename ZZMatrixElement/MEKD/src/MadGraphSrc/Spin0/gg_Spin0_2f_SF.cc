//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.2, 2014-02-07
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "gg_Spin0_2f_SF.h"
#include "../HelAmps_HEF_MEKD2_1.h"	// Changed by Convert_source 0.2

using namespace MG5_HEF_MEKD2_1;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > h > mu- mu+ mu- mu+ S0QQ=1 QED=2 HIW=1 HIGS=1 HIG=1 HIWS=1 /
// zp xg

//--------------------------------------------------------------------------
// Initialize process.

void gg_Spin0_2f_SF::initProc(string param_card_name) 
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
  jamp2[0] = new double[1];
	for( int count=0; count<namplitudes; count++ ) amp[count] = 0;
}

//--------------------------------------------------------------------------
// Update process.	// Created here by Convert_source 0.2

void gg_Spin0_2f_SF::updateProc(SLHAReader_MEKD &slha)
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
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void gg_Spin0_2f_SF::sigmaKin() 
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
        t[0] = matrix_gg_h_mummupmummup_no_zpxg();

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
      t[0] = matrix_gg_h_mummupmummup_no_zpxg();

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

double gg_Spin0_2f_SF::sigmaHat() 
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

void gg_Spin0_2f_SF::calculate_wavefunctions(const int perm[], const int hel[])
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
  VVS3_4_5_3(w[0], w[1], pars->HEF_MEKD2_1_GC_13, pars->HEF_MEKD2_1_GC_15, pars->HEF_MEKD2_1_GC_19, pars->MH,
      pars->WH, w[6]);
  FFV2P0_3(w[3], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[7]);
  FFS1_2_1(w[4], w[6], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO, w[8]);
  FFS1_2_2(w[5], w[6], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO, w[9]);
  FFV5_7_3(w[3], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[10]);
  VVS1_3(w[0], w[1], pars->HEF_MEKD2_1_GC_23, pars->MH, pars->WH, w[11]);
  FFS1_2_1(w[4], w[11], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[12]);
  FFS1_2_2(w[5], w[11], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[13]);
  FFV5_7_3(w[5], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[14]);
  FFV2P0_3(w[5], w[2], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[15]);
  FFS1_2_2(w[3], w[6], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[16]);
  FFV5_7_3(w[5], w[2], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[17]);
  FFS1_2_2(w[3], w[11], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[18]);
  FFV5_7_3(w[3], w[4], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[19]);
  FFV2P0_3(w[3], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[20]);
  FFS1_2_1(w[2], w[6], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[21]);
  FFS1_2_1(w[2], w[11], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, pars->MM, pars->ZERO,
      w[22]);
  FFV2P0_3(w[5], w[4], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[23]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV2_0(w[5], w[8], w[7], pars->HEF_MEKD2_1_GC_5, amp[0]);
  FFV2_0(w[9], w[4], w[7], pars->HEF_MEKD2_1_GC_5, amp[1]);
  FFV5_7_0(w[5], w[8], w[10], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[2]);
  FFV5_7_0(w[9], w[4], w[10], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[3]);
  FFV2_0(w[5], w[12], w[7], pars->HEF_MEKD2_1_GC_5, amp[4]);
  FFV2_0(w[13], w[4], w[7], pars->HEF_MEKD2_1_GC_5, amp[5]);
  FFV5_7_0(w[5], w[12], w[10], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[6]);
  FFV5_7_0(w[13], w[4], w[10], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[7]);
//   VVS3_4_5_0(w[10], w[14], w[6], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[8]);	// ZZ
//   VVS2_0(w[10], w[14], w[6], pars->HEF_MEKD2_1_GC_25, amp[9]);	// ZZ
//   VVS3_4_5_0(w[10], w[14], w[11], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[10]);	// ZZ
//   VVS2_0(w[10], w[14], w[11], pars->HEF_MEKD2_1_GC_25, amp[11]);	// ZZ
  FFV2_0(w[16], w[4], w[15], pars->HEF_MEKD2_1_GC_5, amp[12]);
  FFV2_0(w[3], w[8], w[15], pars->HEF_MEKD2_1_GC_5, amp[13]);
  FFV5_7_0(w[16], w[4], w[17], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[14]);
  FFV5_7_0(w[3], w[8], w[17], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[15]);
  FFV2_0(w[18], w[4], w[15], pars->HEF_MEKD2_1_GC_5, amp[16]);
  FFV2_0(w[3], w[12], w[15], pars->HEF_MEKD2_1_GC_5, amp[17]);
  FFV5_7_0(w[18], w[4], w[17], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[18]);
  FFV5_7_0(w[3], w[12], w[17], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[19]);
//   VVS3_4_5_0(w[17], w[19], w[6], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[20]);	// ZZ
//   VVS2_0(w[17], w[19], w[6], pars->HEF_MEKD2_1_GC_25, amp[21]);	// ZZ
//   VVS3_4_5_0(w[17], w[19], w[11], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[22]);	// ZZ
//   VVS2_0(w[17], w[19], w[11], pars->HEF_MEKD2_1_GC_25, amp[23]);	// ZZ
  FFV2_0(w[5], w[21], w[20], pars->HEF_MEKD2_1_GC_5, amp[24]);
  FFV2_0(w[9], w[2], w[20], pars->HEF_MEKD2_1_GC_5, amp[25]);
  FFV5_7_0(w[5], w[21], w[19], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[26]);
  FFV5_7_0(w[9], w[2], w[19], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[27]);
  FFV2_0(w[5], w[22], w[20], pars->HEF_MEKD2_1_GC_5, amp[28]);
  FFV2_0(w[13], w[2], w[20], pars->HEF_MEKD2_1_GC_5, amp[29]);
  FFV5_7_0(w[5], w[22], w[19], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[30]);
  FFV5_7_0(w[13], w[2], w[19], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[31]);
  FFV2_0(w[3], w[21], w[23], pars->HEF_MEKD2_1_GC_5, amp[32]);
  FFV2_0(w[16], w[2], w[23], pars->HEF_MEKD2_1_GC_5, amp[33]);
  FFV5_7_0(w[3], w[21], w[14], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[34]);
  FFV5_7_0(w[16], w[2], w[14], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[35]);
  FFV2_0(w[3], w[22], w[23], pars->HEF_MEKD2_1_GC_5, amp[36]);
  FFV2_0(w[18], w[2], w[23], pars->HEF_MEKD2_1_GC_5, amp[37]);
  FFV5_7_0(w[3], w[22], w[14], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[38]);
  FFV5_7_0(w[18], w[2], w[14], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, amp[39]);

}
double gg_Spin0_2f_SF::matrix_gg_h_mummupmummup_no_zpxg() 
{
  int i, j;
  // Local variables
	// Commented out by Convert_source 0.2
//  const int ngraphs = 40;
  const int ncolor = 1;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[ncolor] = {1};
  static const double cf[ncolor][ncolor] = {{2}};

  // Calculate color flows
  jamp[0] = +2. * (+amp[0] + amp[1] + amp[2] + amp[3] + amp[4] + amp[5] +
      amp[6] + amp[7] + amp[8] + amp[9] + amp[10] + amp[11] - amp[12] - amp[13]
      - amp[14] - amp[15] - amp[16] - amp[17] - amp[18] - amp[19] - amp[20] -
      amp[21] - amp[22] - amp[23] - amp[24] - amp[25] - amp[26] - amp[27] -
      amp[28] - amp[29] - amp[30] - amp[31] + amp[32] + amp[33] + amp[34] +
      amp[35] + amp[36] + amp[37] + amp[38] + amp[39]);

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



