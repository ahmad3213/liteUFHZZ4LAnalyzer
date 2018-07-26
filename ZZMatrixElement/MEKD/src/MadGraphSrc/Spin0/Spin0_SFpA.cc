//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.9, 2013-04-01
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "Spin0_SFpA.h"
#include "../HelAmps_HEF_MEKD.h"

using namespace MG5_HEF_MEKD; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: h > mu- mu+ mu- mu+ a

//--------------------------------------------------------------------------
// Initialize process.

void Spin0_SFpA::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MEKD::getInstance(); 
  SLHAReader_MEKD slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// moved here by Ghost remover v. 0.1 


  // Set external particle masses for this matrix element
  mME.push_back(pars->MH); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[1]; 
}

//--------------------------------------------------------------------------
// Update process.

void Spin0_SFpA::updateProc(SLHAReader_MEKD &slha) 
{
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// moved here by Ghost remover v. 0.1 
  
  // Set external particle masses for this matrix element
  mME[0]=(pars->MH);
  mME[1]=(pars->MM);
  mME[2]=(pars->MM);
  mME[3]=(pars->MM);
  mME[4]=(pars->MM);
  mME[5]=(pars->ZERO);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void Spin0_SFpA::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
//  static bool firsttime = true; 
//   if (firsttime)
//   {
//     pars->printDependentParameters(); 
//     pars->printDependentCouplings(); 
//     firsttime = false; 
//   }

  // Reset color flows
  for(int i = 0; i < 1; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 32; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  //static int ntry = 0, sum_hel = 0, ngood = 0;	// picked out by Ghost remover v. 0.1 
  static int igood[ncomb]; 
  static int jhel; 
//  std::complex<double> * * wfs;
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{0, -1, -1, -1, -1, -1}, {0,
      -1, -1, -1, -1, 1}, {0, -1, -1, -1, 1, -1}, {0, -1, -1, -1, 1, 1}, {0,
      -1, -1, 1, -1, -1}, {0, -1, -1, 1, -1, 1}, {0, -1, -1, 1, 1, -1}, {0, -1,
      -1, 1, 1, 1}, {0, -1, 1, -1, -1, -1}, {0, -1, 1, -1, -1, 1}, {0, -1, 1,
      -1, 1, -1}, {0, -1, 1, -1, 1, 1}, {0, -1, 1, 1, -1, -1}, {0, -1, 1, 1,
      -1, 1}, {0, -1, 1, 1, 1, -1}, {0, -1, 1, 1, 1, 1}, {0, 1, -1, -1, -1,
      -1}, {0, 1, -1, -1, -1, 1}, {0, 1, -1, -1, 1, -1}, {0, 1, -1, -1, 1, 1},
      {0, 1, -1, 1, -1, -1}, {0, 1, -1, 1, -1, 1}, {0, 1, -1, 1, 1, -1}, {0, 1,
      -1, 1, 1, 1}, {0, 1, 1, -1, -1, -1}, {0, 1, 1, -1, -1, 1}, {0, 1, 1, -1,
      1, -1}, {0, 1, 1, -1, 1, 1}, {0, 1, 1, 1, -1, -1}, {0, 1, 1, 1, -1, 1},
      {0, 1, 1, 1, 1, -1}, {0, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {4}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_h_mummupmummupa(); 

        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
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
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_h_mummupmummupa(); 

      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double Spin0_SFpA::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 9000006 && id2 == 13)
  {
    // Add matrix elements for processes with beams (9000006, 13)
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

void Spin0_SFpA::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
//   int i, j; 

  // Calculate all wavefunctions
  sxxxxx(p[perm[0]], -1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], +1, w[1]); 
  ixxxxx(p[perm[2]], mME[2], hel[2], -1, w[2]); 
  oxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]); 
  ixxxxx(p[perm[4]], mME[4], hel[4], -1, w[4]); 
  vxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]); 
  FFV5_7_3(w[2], w[1], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[6]); 
  FFV2_1(w[3], w[5], pars->HEF_MEKD_GC_5, pars->MM, pars->ZERO, w[7]); 
  FFV5_7_3(w[4], w[7], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[8]); 
  FFV2_2(w[4], w[5], pars->HEF_MEKD_GC_5, pars->MM, pars->ZERO, w[9]); 
  FFV5_7_3(w[9], w[3], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[10]); 
  FFV5_7_3(w[4], w[1], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[11]); 
  FFV2_2(w[2], w[5], pars->HEF_MEKD_GC_5, pars->MM, pars->ZERO, w[12]); 
  FFV5_7_3(w[12], w[3], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[13]); 
  FFV5_7_3(w[2], w[7], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[14]); 
  FFV2_1(w[1], w[5], pars->HEF_MEKD_GC_5, pars->MM, pars->ZERO, w[15]); 
  FFV5_7_3(w[2], w[3], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[16]); 
  FFV5_7_3(w[4], w[15], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[17]); 
  FFV5_7_3(w[4], w[3], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[18]); 
  FFV5_7_3(w[2], w[15], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[19]); 
  FFV5_7_3(w[9], w[1], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[20]); 
  FFV5_7_3(w[12], w[1], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[21]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  VVS3_4_5_0(w[6], w[8], w[0], pars->HEF_MEKD_GC_14, pars->HEF_MEKD_GC_18, pars->HEF_MEKD_GC_22, amp[0]); 
  VVS2_0(w[6], w[8], w[0], pars->HEF_MEKD_GC_25, amp[1]); 
  VVS3_4_5_0(w[6], w[10], w[0], pars->HEF_MEKD_GC_14, pars->HEF_MEKD_GC_18, pars->HEF_MEKD_GC_22, amp[2]); 
  VVS2_0(w[6], w[10], w[0], pars->HEF_MEKD_GC_25, amp[3]); 
  VVS3_4_5_0(w[11], w[13], w[0], pars->HEF_MEKD_GC_14, pars->HEF_MEKD_GC_18, pars->HEF_MEKD_GC_22,
      amp[4]);
  VVS2_0(w[11], w[13], w[0], pars->HEF_MEKD_GC_25, amp[5]); 
  VVS3_4_5_0(w[11], w[14], w[0], pars->HEF_MEKD_GC_14, pars->HEF_MEKD_GC_18, pars->HEF_MEKD_GC_22,
      amp[6]);
  VVS2_0(w[11], w[14], w[0], pars->HEF_MEKD_GC_25, amp[7]); 
  VVS3_4_5_0(w[17], w[16], w[0], pars->HEF_MEKD_GC_14, pars->HEF_MEKD_GC_18, pars->HEF_MEKD_GC_22,
      amp[8]);
  VVS2_0(w[17], w[16], w[0], pars->HEF_MEKD_GC_25, amp[9]); 
  VVS3_4_5_0(w[19], w[18], w[0], pars->HEF_MEKD_GC_14, pars->HEF_MEKD_GC_18, pars->HEF_MEKD_GC_22,
      amp[10]);
  VVS2_0(w[19], w[18], w[0], pars->HEF_MEKD_GC_25, amp[11]); 
  VVS3_4_5_0(w[20], w[16], w[0], pars->HEF_MEKD_GC_14, pars->HEF_MEKD_GC_18, pars->HEF_MEKD_GC_22,
      amp[12]);
  VVS2_0(w[20], w[16], w[0], pars->HEF_MEKD_GC_25, amp[13]); 
  VVS3_4_5_0(w[21], w[18], w[0], pars->HEF_MEKD_GC_14, pars->HEF_MEKD_GC_18, pars->HEF_MEKD_GC_22,
      amp[14]);
  VVS2_0(w[21], w[18], w[0], pars->HEF_MEKD_GC_25, amp[15]); 

}
double Spin0_SFpA::matrix_h_mummupmummupa() 
{
  int i, j; 
  // Local variables
//   const int ngraphs = 16; 
  const int ncolor = 1; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[1] = {1.}; 
  static const double cf[1][1] = {{1.}}; 

  // Calculate color flows
  jamp[0] = +amp[0] + amp[1] + amp[2] + amp[3] - amp[4] - amp[5] - amp[6] -
      amp[7] - amp[8] - amp[9] + amp[10] + amp[11] - amp[12] - amp[13] +
      amp[14] + amp[15];

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



