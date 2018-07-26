//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.5, 2012-11-18
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "qq_Z4l_SIG_DN_OFpA.h"
#include "../HelAmps_HZZ_Unitary_bkgpA.h"
#include "../read_slha.h"

using namespace MG5_HZZ_Unitary_bkgpA; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: s s~ > e- e+ mu- mu+ a / h zp xg WEIGHTED=10

//--------------------------------------------------------------------------
// Initialize process.

void qq_Z4l_SIG_DN_OFpA::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MEKD::getInstance(); 
  SLHAReader_MEKD slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// moved here by Ghost remover v. 0.1 
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
// Update process.

void qq_Z4l_SIG_DN_OFpA::updateProc(SLHAReader_MEKD &slha) 
{
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// moved here by Ghost remover v. 0.1 
  
  // Set external particle masses for this matrix element
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

void qq_Z4l_SIG_DN_OFpA::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 

  // Reset color flows
  for(int i = 0; i < 1; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 128; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  //static int ntry = 0, sum_hel = 0, ngood = 0;	// picked out by Ghost remover v. 0.1 
  static int igood[ncomb]; 
  static int jhel; 
//   std::complex<double> * * wfs; 
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
        t[0] = matrix_ssx_emepmummupa_no_hzpxg(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_ssx_emepmummupa_no_hzpxg(); 
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
      t[0] = matrix_ssx_emepmummupa_no_hzpxg(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_ssx_emepmummupa_no_hzpxg(); 
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

double qq_Z4l_SIG_DN_OFpA::sigmaHat() 
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

void qq_Z4l_SIG_DN_OFpA::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
//   int i, j; 

  // Calculate all wavefunctions
  ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  ixxxxx(p[perm[5]], mME[5], hel[5], -1, w[5]); 
  vxxxxx(p[perm[6]], mME[6], hel[6], +1, w[6]); 
//   FFV1_3(w[0], w[1], pars->Unitary_GC_5, pars->ZERO, pars->ZERO, w[7]); 
  FFV1_3(w[3], w[2], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[8]); 
//   FFV1_1(w[4], w[7], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[9]); 
  FFV1_2(w[5], w[8], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[10]); 
//   FFV1_2(w[5], w[7], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[11]); 
  FFV1_1(w[4], w[8], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[12]); 
  FFV2_4_3(w[3], w[2], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[13]); 
  FFV2_4_2(w[5], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[14]); 
  FFV2_4_1(w[4], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[15]); 
  FFV2_3_3(w[0], w[1], pars->Unitary_GC_70, pars->Unitary_GC_74, pars->MZ, pars->WZ, w[16]); 
  FFV2_4_1(w[4], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[17]); 
  FFV2_4_2(w[5], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[18]); 
  FFV1_1(w[4], w[6], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[19]); 
//   FFV1_1(w[19], w[7], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[20]); 
  FFV2_4_1(w[19], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO,
      w[21]);
  FFV1_2(w[5], w[6], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[22]); 
//   FFV1_2(w[22], w[7], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[23]); 
  FFV2_4_2(w[22], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO,
      w[24]);
  FFV1_1(w[2], w[6], pars->Unitary_GC_7, pars->Me, pars->ZERO, w[25]); 
  FFV1_3(w[3], w[25], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[26]); 
  FFV2_4_3(w[3], w[25], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[27]); 
  FFV1_3(w[5], w[4], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[28]); 
//   FFV1_1(w[25], w[7], pars->Unitary_GC_7, pars->Me, pars->ZERO, w[29]); 
//   FFV1_2(w[3], w[7], pars->Unitary_GC_7, pars->Me, pars->ZERO, w[30]); 
  FFV2_4_3(w[5], w[4], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[31]); 
  FFV2_4_1(w[25], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->Me, pars->ZERO,
      w[32]);
  FFV2_4_2(w[3], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->Me, pars->ZERO, w[33]); 
  FFV1_2(w[3], w[6], pars->Unitary_GC_7, pars->Me, pars->ZERO, w[34]); 
  FFV1_3(w[34], w[2], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[35]); 
  FFV2_4_3(w[34], w[2], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[36]); 
//   FFV1_1(w[2], w[7], pars->Unitary_GC_7, pars->Me, pars->ZERO, w[37]); 
//   FFV1_2(w[34], w[7], pars->Unitary_GC_7, pars->Me, pars->ZERO, w[38]); 
  FFV2_4_1(w[2], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->Me, pars->ZERO, w[39]); 
  FFV2_4_2(w[34], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->Me, pars->ZERO,
      w[40]);
  FFV1_2(w[3], w[28], pars->Unitary_GC_7, pars->Me, pars->ZERO, w[41]); 
  FFV1_1(w[2], w[28], pars->Unitary_GC_7, pars->Me, pars->ZERO, w[42]); 
  FFV2_4_2(w[3], w[31], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->Me, pars->ZERO, w[43]); 
  FFV2_4_1(w[2], w[31], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->Me, pars->ZERO, w[44]); 
  FFV1_3(w[5], w[19], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[45]); 
  FFV2_4_3(w[5], w[19], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[46]); 
  FFV1_3(w[22], w[4], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[47]); 
  FFV2_4_3(w[22], w[4], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[48]); 
//   FFV1_2(w[0], w[6], pars->Unitary_GC_5, pars->MS, pars->ZERO, w[49]); 
//   FFV1_3(w[49], w[1], pars->Unitary_GC_5, pars->ZERO, pars->ZERO, w[50]); 
//   FFV2_3_3(w[49], w[1], pars->Unitary_GC_70, pars->Unitary_GC_74, pars->MZ, pars->WZ, w[51]); 
//   FFV1_2(w[49], w[8], pars->Unitary_GC_5, pars->MS, pars->ZERO, w[52]); 
//   FFV1_2(w[49], w[28], pars->Unitary_GC_5, pars->MS, pars->ZERO, w[53]); 
//   FFV2_3_2(w[49], w[31], pars->Unitary_GC_70, pars->Unitary_GC_74, pars->MS, pars->ZERO,
//       w[54]);
//   FFV2_3_2(w[49], w[13], pars->Unitary_GC_70, pars->Unitary_GC_74, pars->MS, pars->ZERO,
//       w[55]);
//   FFV1_1(w[1], w[6], pars->Unitary_GC_5, pars->MS, pars->ZERO, w[56]); 
//   FFV1_3(w[0], w[56], pars->Unitary_GC_5, pars->ZERO, pars->ZERO, w[57]); 
//   FFV2_3_3(w[0], w[56], pars->Unitary_GC_70, pars->Unitary_GC_74, pars->MZ, pars->WZ, w[58]); 
//   FFV1_2(w[0], w[8], pars->Unitary_GC_5, pars->MS, pars->ZERO, w[59]); 
//   FFV1_2(w[0], w[28], pars->Unitary_GC_5, pars->MS, pars->ZERO, w[60]); 
//   FFV2_3_2(w[0], w[31], pars->Unitary_GC_70, pars->Unitary_GC_74, pars->MS, pars->ZERO, w[61]); 
//   FFV2_3_2(w[0], w[13], pars->Unitary_GC_70, pars->Unitary_GC_74, pars->MS, pars->ZERO, w[62]); 
//   FFV1_1(w[1], w[28], pars->Unitary_GC_5, pars->MS, pars->ZERO, w[63]); 
//   FFV1_1(w[1], w[8], pars->Unitary_GC_5, pars->MS, pars->ZERO, w[64]); 
//   FFV2_3_1(w[1], w[31], pars->Unitary_GC_70, pars->Unitary_GC_74, pars->MS, pars->ZERO, w[65]); 
//   FFV2_3_1(w[1], w[13], pars->Unitary_GC_70, pars->Unitary_GC_74, pars->MS, pars->ZERO, w[66]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
//   FFV1_0(w[10], w[9], w[6], pars->Unitary_GC_7, amp[0]); 
//   FFV1_0(w[11], w[12], w[6], pars->Unitary_GC_7, amp[1]); 
//   FFV1_0(w[14], w[9], w[6], pars->Unitary_GC_7, amp[2]); 
//   FFV1_0(w[11], w[15], w[6], pars->Unitary_GC_7, amp[3]); 
  FFV1_0(w[10], w[17], w[6], pars->Unitary_GC_7, amp[4]); 
  FFV1_0(w[18], w[12], w[6], pars->Unitary_GC_7, amp[5]); 
  FFV1_0(w[14], w[17], w[6], pars->Unitary_GC_7, amp[6]); 
  FFV1_0(w[18], w[15], w[6], pars->Unitary_GC_7, amp[7]); 
//   FFV1_0(w[5], w[20], w[8], pars->Unitary_GC_7, amp[8]); 
//   FFV1_0(w[11], w[19], w[8], pars->Unitary_GC_7, amp[9]); 
//   FFV2_4_0(w[5], w[20], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[10]); 
//   FFV2_4_0(w[11], w[19], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[11]); 
  FFV1_0(w[5], w[21], w[8], pars->Unitary_GC_7, amp[12]); 
  FFV1_0(w[18], w[19], w[8], pars->Unitary_GC_7, amp[13]); 
  FFV2_4_0(w[5], w[21], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[14]); 
  FFV2_4_0(w[18], w[19], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[15]); 
//   FFV1_0(w[22], w[9], w[8], pars->Unitary_GC_7, amp[16]); 
//   FFV1_0(w[23], w[4], w[8], pars->Unitary_GC_7, amp[17]); 
//   FFV2_4_0(w[22], w[9], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[18]); 
//   FFV2_4_0(w[23], w[4], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[19]); 
  FFV1_0(w[22], w[17], w[8], pars->Unitary_GC_7, amp[20]); 
  FFV1_0(w[24], w[4], w[8], pars->Unitary_GC_7, amp[21]); 
  FFV2_4_0(w[22], w[17], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[22]); 
  FFV2_4_0(w[24], w[4], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[23]); 
//   FFV1_0(w[5], w[9], w[26], pars->Unitary_GC_7, amp[24]); 
//   FFV2_4_0(w[5], w[9], w[27], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[25]); 
//   FFV1_0(w[11], w[4], w[26], pars->Unitary_GC_7, amp[26]); 
//   FFV2_4_0(w[11], w[4], w[27], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[27]); 
  FFV1_0(w[5], w[17], w[26], pars->Unitary_GC_7, amp[28]); 
  FFV2_4_0(w[5], w[17], w[27], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[29]); 
  FFV1_0(w[18], w[4], w[26], pars->Unitary_GC_7, amp[30]); 
  FFV2_4_0(w[18], w[4], w[27], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[31]); 
//   FFV1_0(w[3], w[29], w[28], pars->Unitary_GC_7, amp[32]); 
//   FFV1_0(w[30], w[25], w[28], pars->Unitary_GC_7, amp[33]); 
//   FFV2_4_0(w[3], w[29], w[31], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[34]); 
//   FFV2_4_0(w[30], w[25], w[31], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[35]); 
  FFV1_0(w[3], w[32], w[28], pars->Unitary_GC_7, amp[36]); 
  FFV1_0(w[33], w[25], w[28], pars->Unitary_GC_7, amp[37]); 
  FFV2_4_0(w[3], w[32], w[31], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[38]); 
  FFV2_4_0(w[33], w[25], w[31], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[39]); 
//   FFV1_0(w[5], w[9], w[35], pars->Unitary_GC_7, amp[40]); 
//   FFV2_4_0(w[5], w[9], w[36], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[41]); 
//   FFV1_0(w[11], w[4], w[35], pars->Unitary_GC_7, amp[42]); 
//   FFV2_4_0(w[11], w[4], w[36], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[43]); 
  FFV1_0(w[5], w[17], w[35], pars->Unitary_GC_7, amp[44]); 
  FFV2_4_0(w[5], w[17], w[36], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[45]); 
  FFV1_0(w[18], w[4], w[35], pars->Unitary_GC_7, amp[46]); 
  FFV2_4_0(w[18], w[4], w[36], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[47]); 
//   FFV1_0(w[34], w[37], w[28], pars->Unitary_GC_7, amp[48]); 
//   FFV1_0(w[38], w[2], w[28], pars->Unitary_GC_7, amp[49]); 
//   FFV2_4_0(w[34], w[37], w[31], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[50]); 
//   FFV2_4_0(w[38], w[2], w[31], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[51]); 
  FFV1_0(w[34], w[39], w[28], pars->Unitary_GC_7, amp[52]); 
  FFV1_0(w[40], w[2], w[28], pars->Unitary_GC_7, amp[53]); 
  FFV2_4_0(w[34], w[39], w[31], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[54]); 
  FFV2_4_0(w[40], w[2], w[31], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[55]); 
//   FFV1_0(w[41], w[37], w[6], pars->Unitary_GC_7, amp[56]); 
//   FFV1_0(w[30], w[42], w[6], pars->Unitary_GC_7, amp[57]); 
//   FFV1_0(w[43], w[37], w[6], pars->Unitary_GC_7, amp[58]); 
//   FFV1_0(w[30], w[44], w[6], pars->Unitary_GC_7, amp[59]); 
//   FFV1_0(w[41], w[39], w[6], pars->Unitary_GC_7, amp[60]); 
//   FFV1_0(w[33], w[42], w[6], pars->Unitary_GC_7, amp[61]); 
//   FFV1_0(w[43], w[39], w[6], pars->Unitary_GC_7, amp[62]); 
//   FFV1_0(w[33], w[44], w[6], pars->Unitary_GC_7, amp[63]); 
//   FFV1_0(w[3], w[37], w[45], pars->Unitary_GC_7, amp[64]); 
//   FFV2_4_0(w[3], w[37], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[65]); 
//   FFV1_0(w[30], w[2], w[45], pars->Unitary_GC_7, amp[66]); 
//   FFV2_4_0(w[30], w[2], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[67]); 
  FFV1_0(w[3], w[39], w[45], pars->Unitary_GC_7, amp[68]); 
  FFV2_4_0(w[3], w[39], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[69]); 
  FFV1_0(w[33], w[2], w[45], pars->Unitary_GC_7, amp[70]); 
  FFV2_4_0(w[33], w[2], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[71]); 
//   FFV1_0(w[3], w[37], w[47], pars->Unitary_GC_7, amp[72]); 
//   FFV2_4_0(w[3], w[37], w[48], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[73]); 
//   FFV1_0(w[30], w[2], w[47], pars->Unitary_GC_7, amp[74]); 
//   FFV2_4_0(w[30], w[2], w[48], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[75]); 
  FFV1_0(w[3], w[39], w[47], pars->Unitary_GC_7, amp[76]); 
  FFV2_4_0(w[3], w[39], w[48], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[77]); 
  FFV1_0(w[33], w[2], w[47], pars->Unitary_GC_7, amp[78]); 
  FFV2_4_0(w[33], w[2], w[48], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[79]); 
//   FFV1_0(w[5], w[12], w[50], pars->Unitary_GC_7, amp[80]); 
//   FFV2_4_0(w[5], w[12], w[51], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[81]); 
//   FFV1_0(w[10], w[4], w[50], pars->Unitary_GC_7, amp[82]); 
//   FFV2_4_0(w[10], w[4], w[51], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[83]); 
//   FFV1_0(w[5], w[15], w[50], pars->Unitary_GC_7, amp[84]); 
//   FFV2_4_0(w[5], w[15], w[51], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[85]); 
//   FFV1_0(w[14], w[4], w[50], pars->Unitary_GC_7, amp[86]); 
//   FFV2_4_0(w[14], w[4], w[51], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[87]);
//   FFV1_0(w[52], w[1], w[28], pars->Unitary_GC_5, amp[88]); 
//   FFV1_0(w[53], w[1], w[8], pars->Unitary_GC_5, amp[89]); 
//   FFV2_3_0(w[52], w[1], w[31], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[90]); 
//   FFV1_0(w[54], w[1], w[8], pars->Unitary_GC_5, amp[91]); 
//   FFV1_0(w[55], w[1], w[28], pars->Unitary_GC_5, amp[92]); 
//   FFV2_3_0(w[53], w[1], w[13], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[93]); 
//   FFV2_3_0(w[55], w[1], w[31], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[94]); 
//   FFV2_3_0(w[54], w[1], w[13], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[95]); 
//   FFV1_0(w[3], w[42], w[50], pars->Unitary_GC_7, amp[96]); 
//   FFV2_4_0(w[3], w[42], w[51], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[97]); 
//   FFV1_0(w[41], w[2], w[50], pars->Unitary_GC_7, amp[98]); 
//   FFV2_4_0(w[41], w[2], w[51], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[99]); 
//   FFV1_0(w[3], w[44], w[50], pars->Unitary_GC_7, amp[100]); 
//   FFV2_4_0(w[3], w[44], w[51], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[101]); 
//   FFV1_0(w[43], w[2], w[50], pars->Unitary_GC_7, amp[102]); 
//   FFV2_4_0(w[43], w[2], w[51], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[103]); 
//   FFV1_0(w[5], w[12], w[57], pars->Unitary_GC_7, amp[104]); 
//   FFV2_4_0(w[5], w[12], w[58], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[105]); 
//   FFV1_0(w[10], w[4], w[57], pars->Unitary_GC_7, amp[106]); 
//   FFV2_4_0(w[10], w[4], w[58], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[107]); 
//   FFV1_0(w[5], w[15], w[57], pars->Unitary_GC_7, amp[108]); 
//   FFV2_4_0(w[5], w[15], w[58], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[109]); 
//   FFV1_0(w[14], w[4], w[57], pars->Unitary_GC_7, amp[110]); 
//   FFV2_4_0(w[14], w[4], w[58], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[111]); 
//   FFV1_0(w[59], w[56], w[28], pars->Unitary_GC_5, amp[112]); 
//   FFV1_0(w[60], w[56], w[8], pars->Unitary_GC_5, amp[113]); 
//   FFV2_3_0(w[59], w[56], w[31], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[114]); 
//   FFV1_0(w[61], w[56], w[8], pars->Unitary_GC_5, amp[115]); 
//   FFV1_0(w[62], w[56], w[28], pars->Unitary_GC_5, amp[116]); 
//   FFV2_3_0(w[60], w[56], w[13], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[117]); 
//   FFV2_3_0(w[62], w[56], w[31], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[118]); 
//   FFV2_3_0(w[61], w[56], w[13], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[119]); 
//   FFV1_0(w[3], w[42], w[57], pars->Unitary_GC_7, amp[120]); 
//   FFV2_4_0(w[3], w[42], w[58], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[121]); 
//   FFV1_0(w[41], w[2], w[57], pars->Unitary_GC_7, amp[122]); 
//   FFV2_4_0(w[41], w[2], w[58], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[123]); 
//   FFV1_0(w[3], w[44], w[57], pars->Unitary_GC_7, amp[124]); 
//   FFV2_4_0(w[3], w[44], w[58], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[125]); 
//   FFV1_0(w[43], w[2], w[57], pars->Unitary_GC_7, amp[126]); 
//   FFV2_4_0(w[43], w[2], w[58], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[127]); 
//   FFV1_0(w[59], w[63], w[6], pars->Unitary_GC_5, amp[128]); 
//   FFV1_0(w[60], w[64], w[6], pars->Unitary_GC_5, amp[129]); 
//   FFV1_0(w[59], w[65], w[6], pars->Unitary_GC_5, amp[130]); 
//   FFV1_0(w[61], w[64], w[6], pars->Unitary_GC_5, amp[131]); 
//   FFV1_0(w[62], w[63], w[6], pars->Unitary_GC_5, amp[132]); 
//   FFV1_0(w[60], w[66], w[6], pars->Unitary_GC_5, amp[133]); 
//   FFV1_0(w[62], w[65], w[6], pars->Unitary_GC_5, amp[134]); 
//   FFV1_0(w[61], w[66], w[6], pars->Unitary_GC_5, amp[135]); 
//   FFV1_0(w[59], w[1], w[45], pars->Unitary_GC_5, amp[136]); 
//   FFV2_3_0(w[59], w[1], w[46], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[137]); 
//   FFV1_0(w[0], w[64], w[45], pars->Unitary_GC_5, amp[138]); 
//   FFV2_3_0(w[0], w[64], w[46], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[139]); 
//   FFV1_0(w[62], w[1], w[45], pars->Unitary_GC_5, amp[140]); 
//   FFV2_3_0(w[62], w[1], w[46], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[141]); 
//   FFV1_0(w[0], w[66], w[45], pars->Unitary_GC_5, amp[142]); 
//   FFV2_3_0(w[0], w[66], w[46], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[143]); 
//   FFV1_0(w[59], w[1], w[47], pars->Unitary_GC_5, amp[144]); 
//   FFV2_3_0(w[59], w[1], w[48], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[145]); 
//   FFV1_0(w[0], w[64], w[47], pars->Unitary_GC_5, amp[146]); 
//   FFV2_3_0(w[0], w[64], w[48], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[147]); 
//   FFV1_0(w[62], w[1], w[47], pars->Unitary_GC_5, amp[148]); 
//   FFV2_3_0(w[62], w[1], w[48], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[149]); 
//   FFV1_0(w[0], w[66], w[47], pars->Unitary_GC_5, amp[150]); 
//   FFV2_3_0(w[0], w[66], w[48], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[151]); 
//   FFV1_0(w[60], w[1], w[26], pars->Unitary_GC_5, amp[152]); 
//   FFV2_3_0(w[60], w[1], w[27], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[153]); 
//   FFV1_0(w[0], w[63], w[26], pars->Unitary_GC_5, amp[154]); 
//   FFV2_3_0(w[0], w[63], w[27], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[155]); 
//   FFV1_0(w[61], w[1], w[26], pars->Unitary_GC_5, amp[156]); 
//   FFV2_3_0(w[61], w[1], w[27], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[157]); 
//   FFV1_0(w[0], w[65], w[26], pars->Unitary_GC_5, amp[158]); 
//   FFV2_3_0(w[0], w[65], w[27], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[159]); 
//   FFV1_0(w[60], w[1], w[35], pars->Unitary_GC_5, amp[160]); 
//   FFV2_3_0(w[60], w[1], w[36], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[161]); 
//   FFV1_0(w[0], w[63], w[35], pars->Unitary_GC_5, amp[162]); 
//   FFV2_3_0(w[0], w[63], w[36], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[163]); 
//   FFV1_0(w[61], w[1], w[35], pars->Unitary_GC_5, amp[164]); 
//   FFV2_3_0(w[61], w[1], w[36], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[165]); 
//   FFV1_0(w[0], w[65], w[35], pars->Unitary_GC_5, amp[166]); 
//   FFV2_3_0(w[0], w[65], w[36], pars->Unitary_GC_70, pars->Unitary_GC_74, amp[167]); 

}


double qq_Z4l_SIG_DN_OFpA::matrix_ssx_emepmummupa_no_hzpxg() 
{
  int i, j; 
  // Local variables
//   const int ngraphs = 168; 
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
      amp[165] + amp[166] + amp[167];

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



