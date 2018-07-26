//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.5, 2012-11-18
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "qq_Z4l_SIG_UP_SFpA.h"
#include "../HelAmps_HZZ_Unitary_bkgpA.h"
#include "../read_slha.h"

using namespace MG5_HZZ_Unitary_bkgpA; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: c c~ > mu- mu+ mu- mu+ a / h zp xg WEIGHTED=10

//--------------------------------------------------------------------------
// Initialize process.

void qq_Z4l_SIG_UP_SFpA::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MEKD::getInstance(); 
  SLHAReader_MEKD slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// moved here by Ghost remover v. 0.1 
  // Set external particle masses for this matrix element
  mME.push_back(pars->MC); 
  mME.push_back(pars->MC); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[1];
	for( int count=0; count<namplitudes; count++ ) amp[count] = 0;
}

//--------------------------------------------------------------------------
// Update process.

void qq_Z4l_SIG_UP_SFpA::updateProc(SLHAReader_MEKD &slha) 
{
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// moved here by Ghost remover v. 0.1 
  
  // Set external particle masses for this matrix element
  mME[0]=(pars->MC);
  mME[1]=(pars->MC);
  mME[2]=(pars->MM);
  mME[3]=(pars->MM);
  mME[4]=(pars->MM);
  mME[5]=(pars->MM);
  mME[6]=(pars->ZERO);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void qq_Z4l_SIG_UP_SFpA::sigmaKin() 
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
  const int denominators[nprocesses] = {144, 144}; 

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
        t[0] = matrix_ccx_mummupmummupa_no_hzpxg(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_ccx_mummupmummupa_no_hzpxg(); 
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
      t[0] = matrix_ccx_mummupmummupa_no_hzpxg(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_ccx_mummupmummupa_no_hzpxg(); 
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

double qq_Z4l_SIG_UP_SFpA::sigmaHat() 
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

void qq_Z4l_SIG_UP_SFpA::calculate_wavefunctions(const int perm[], const int hel[])
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
//   FFV1_3(w[0], w[1], pars->Unitary_GC_6, pars->ZERO, pars->ZERO, w[7]); 
  FFV1_3(w[3], w[2], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[8]); 
//   FFV1_1(w[4], w[7], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[9]); 
  FFV1_2(w[5], w[8], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[10]); 
//   FFV1_2(w[5], w[7], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[11]); 
  FFV1_1(w[4], w[8], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[12]); 
  FFV2_4_3(w[3], w[2], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[13]); 
  FFV2_4_2(w[5], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[14]); 
  FFV2_4_1(w[4], w[13], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[15]); 
  FFV2_5_3(w[0], w[1], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MZ, pars->WZ, w[16]); 
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
  FFV1_3(w[5], w[2], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[25]); 
//   FFV1_2(w[3], w[7], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[26]); 
  FFV1_1(w[4], w[25], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[27]); 
  FFV1_2(w[3], w[25], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[28]); 
  FFV2_4_3(w[5], w[2], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[29]); 
  FFV2_4_1(w[4], w[29], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[30]); 
  FFV2_4_2(w[3], w[29], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[31]); 
  FFV2_4_2(w[3], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[32]); 
  FFV1_2(w[3], w[6], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[33]); 
//   FFV1_2(w[33], w[7], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[34]); 
  FFV2_4_2(w[33], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO,
      w[35]);
  FFV1_1(w[2], w[6], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[36]); 
  FFV1_3(w[5], w[36], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[37]); 
  FFV2_4_3(w[5], w[36], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[38]); 
  FFV1_3(w[3], w[36], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[39]); 
  FFV2_4_3(w[3], w[36], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[40]); 
  FFV1_3(w[3], w[4], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[41]); 
//   FFV1_1(w[36], w[7], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[42]); 
  FFV2_4_3(w[3], w[4], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[43]); 
  FFV2_4_1(w[36], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO,
      w[44]);
  FFV1_3(w[5], w[4], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[45]); 
  FFV2_4_3(w[5], w[4], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[46]); 
//   FFV1_1(w[2], w[7], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[47]); 
  FFV1_2(w[5], w[41], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[48]); 
  FFV1_1(w[2], w[41], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[49]); 
  FFV2_4_2(w[5], w[43], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[50]); 
  FFV2_4_1(w[2], w[43], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[51]); 
  FFV2_4_1(w[2], w[16], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[52]); 
  FFV1_3(w[33], w[4], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[53]); 
  FFV2_4_3(w[33], w[4], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[54]); 
  FFV1_3(w[33], w[2], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[55]); 
  FFV2_4_3(w[33], w[2], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[56]); 
  FFV1_2(w[3], w[45], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[57]); 
  FFV1_1(w[2], w[45], pars->Unitary_GC_7, pars->MM, pars->ZERO, w[58]); 
  FFV2_4_2(w[3], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[59]); 
  FFV2_4_1(w[2], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MM, pars->ZERO, w[60]); 
  FFV1_3(w[3], w[19], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[61]); 
  FFV2_4_3(w[3], w[19], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[62]); 
  FFV1_3(w[5], w[19], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[63]); 
  FFV2_4_3(w[5], w[19], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[64]); 
  FFV1_3(w[22], w[4], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[65]); 
  FFV2_4_3(w[22], w[4], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[66]); 
  FFV1_3(w[22], w[2], pars->Unitary_GC_7, pars->ZERO, pars->ZERO, w[67]); 
  FFV2_4_3(w[22], w[2], pars->Unitary_GC_70, pars->Unitary_GC_75, pars->MZ, pars->WZ, w[68]); 
//   FFV1_2(w[0], w[6], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[69]); 
  FFV1_3(w[69], w[1], pars->Unitary_GC_6, pars->ZERO, pars->ZERO, w[70]); 
  FFV2_5_3(w[69], w[1], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MZ, pars->WZ, w[71]); 
  FFV1_2(w[69], w[8], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[72]); 
  FFV1_2(w[69], w[45], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[73]); 
  FFV2_5_2(w[69], w[46], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO,
      w[74]);
  FFV2_5_2(w[69], w[13], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO,
      w[75]);
  FFV1_2(w[69], w[25], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[76]); 
  FFV1_2(w[69], w[41], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[77]); 
  FFV2_5_2(w[69], w[43], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO,
      w[78]);
  FFV2_5_2(w[69], w[29], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO,
      w[79]);
//   FFV1_1(w[1], w[6], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[80]); 
//   FFV1_3(w[0], w[80], pars->Unitary_GC_6, pars->ZERO, pars->ZERO, w[81]); 
  FFV2_5_3(w[0], w[80], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MZ, pars->WZ, w[82]); 
//   FFV1_2(w[0], w[8], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[83]); 
//   FFV1_2(w[0], w[45], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[84]); 
//   FFV2_5_2(w[0], w[46], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO, w[85]); 
//   FFV2_5_2(w[0], w[13], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO, w[86]); 
//   FFV1_2(w[0], w[25], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[87]); 
//   FFV1_2(w[0], w[41], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[88]); 
//   FFV2_5_2(w[0], w[43], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO, w[89]); 
//   FFV2_5_2(w[0], w[29], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO, w[90]); 
//   FFV1_1(w[1], w[45], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[91]); 
//   FFV1_1(w[1], w[8], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[92]); 
//   FFV2_5_1(w[1], w[46], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO, w[93]); 
//   FFV2_5_1(w[1], w[13], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO, w[94]); 
//   FFV1_1(w[1], w[41], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[95]); 
//   FFV1_1(w[1], w[25], pars->Unitary_GC_6, pars->MC, pars->ZERO, w[96]); 
//   FFV2_5_1(w[1], w[43], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO, w[97]); 
//   FFV2_5_1(w[1], w[29], pars->Unitary_GC_71, pars->Unitary_GC_74, pars->MC, pars->ZERO, w[98]); 

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
//   FFV1_0(w[26], w[27], w[6], pars->Unitary_GC_7, amp[24]); 
//   FFV1_0(w[28], w[9], w[6], pars->Unitary_GC_7, amp[25]); 
//   FFV1_0(w[26], w[30], w[6], pars->Unitary_GC_7, amp[26]); 
//   FFV1_0(w[31], w[9], w[6], pars->Unitary_GC_7, amp[27]); 
  FFV1_0(w[32], w[27], w[6], pars->Unitary_GC_7, amp[28]); 
  FFV1_0(w[28], w[17], w[6], pars->Unitary_GC_7, amp[29]); 
  FFV1_0(w[32], w[30], w[6], pars->Unitary_GC_7, amp[30]); 
  FFV1_0(w[31], w[17], w[6], pars->Unitary_GC_7, amp[31]); 
//   FFV1_0(w[34], w[4], w[25], pars->Unitary_GC_7, amp[32]); 
//   FFV1_0(w[33], w[9], w[25], pars->Unitary_GC_7, amp[33]); 
//   FFV2_4_0(w[34], w[4], w[29], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[34]); 
//   FFV2_4_0(w[33], w[9], w[29], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[35]); 
  FFV1_0(w[35], w[4], w[25], pars->Unitary_GC_7, amp[36]); 
  FFV1_0(w[33], w[17], w[25], pars->Unitary_GC_7, amp[37]); 
  FFV2_4_0(w[35], w[4], w[29], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[38]); 
  FFV2_4_0(w[33], w[17], w[29], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[39]); 
//   FFV1_0(w[26], w[19], w[25], pars->Unitary_GC_7, amp[40]); 
//   FFV1_0(w[3], w[20], w[25], pars->Unitary_GC_7, amp[41]); 
//   FFV2_4_0(w[26], w[19], w[29], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[42]); 
//   FFV2_4_0(w[3], w[20], w[29], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[43]); 
  FFV1_0(w[32], w[19], w[25], pars->Unitary_GC_7, amp[44]); 
  FFV1_0(w[3], w[21], w[25], pars->Unitary_GC_7, amp[45]); 
  FFV2_4_0(w[32], w[19], w[29], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[46]); 
  FFV2_4_0(w[3], w[21], w[29], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[47]); 
//   FFV1_0(w[26], w[4], w[37], pars->Unitary_GC_7, amp[48]); 
//   FFV2_4_0(w[26], w[4], w[38], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[49]); 
//   FFV1_0(w[5], w[9], w[39], pars->Unitary_GC_7, amp[50]); 
//   FFV2_4_0(w[5], w[9], w[40], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[51]); 
//   FFV1_0(w[3], w[9], w[37], pars->Unitary_GC_7, amp[52]); 
//   FFV2_4_0(w[3], w[9], w[38], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[53]); 
//   FFV1_0(w[11], w[4], w[39], pars->Unitary_GC_7, amp[54]); 
//   FFV2_4_0(w[11], w[4], w[40], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[55]); 
  FFV1_0(w[32], w[4], w[37], pars->Unitary_GC_7, amp[56]); 
  FFV2_4_0(w[32], w[4], w[38], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[57]); 
  FFV1_0(w[5], w[17], w[39], pars->Unitary_GC_7, amp[58]); 
  FFV2_4_0(w[5], w[17], w[40], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[59]); 
  FFV1_0(w[3], w[17], w[37], pars->Unitary_GC_7, amp[60]); 
  FFV2_4_0(w[3], w[17], w[38], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[61]); 
  FFV1_0(w[18], w[4], w[39], pars->Unitary_GC_7, amp[62]); 
  FFV2_4_0(w[18], w[4], w[40], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[63]); 
//   FFV1_0(w[5], w[42], w[41], pars->Unitary_GC_7, amp[64]); 
//   FFV1_0(w[11], w[36], w[41], pars->Unitary_GC_7, amp[65]); 
//   FFV2_4_0(w[5], w[42], w[43], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[66]); 
//   FFV2_4_0(w[11], w[36], w[43], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[67]); 
  FFV1_0(w[5], w[44], w[41], pars->Unitary_GC_7, amp[68]); 
  FFV1_0(w[18], w[36], w[41], pars->Unitary_GC_7, amp[69]); 
  FFV2_4_0(w[5], w[44], w[43], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[70]); 
  FFV2_4_0(w[18], w[36], w[43], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[71]); 
//   FFV1_0(w[3], w[42], w[45], pars->Unitary_GC_7, amp[72]); 
//   FFV1_0(w[26], w[36], w[45], pars->Unitary_GC_7, amp[73]); 
//   FFV2_4_0(w[3], w[42], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[74]); 
//   FFV2_4_0(w[26], w[36], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[75]); 
  FFV1_0(w[3], w[44], w[45], pars->Unitary_GC_7, amp[76]); 
  FFV1_0(w[32], w[36], w[45], pars->Unitary_GC_7, amp[77]); 
  FFV2_4_0(w[3], w[44], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[78]); 
  FFV2_4_0(w[32], w[36], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[79]); 
//   FFV1_0(w[48], w[47], w[6], pars->Unitary_GC_7, amp[80]); 
//   FFV1_0(w[11], w[49], w[6], pars->Unitary_GC_7, amp[81]); 
//   FFV1_0(w[50], w[47], w[6], pars->Unitary_GC_7, amp[82]); 
//   FFV1_0(w[11], w[51], w[6], pars->Unitary_GC_7, amp[83]); 
  FFV1_0(w[48], w[52], w[6], pars->Unitary_GC_7, amp[84]); 
  FFV1_0(w[18], w[49], w[6], pars->Unitary_GC_7, amp[85]); 
  FFV1_0(w[50], w[52], w[6], pars->Unitary_GC_7, amp[86]); 
  FFV1_0(w[18], w[51], w[6], pars->Unitary_GC_7, amp[87]); 
//   FFV1_0(w[22], w[47], w[41], pars->Unitary_GC_7, amp[88]); 
//   FFV1_0(w[23], w[2], w[41], pars->Unitary_GC_7, amp[89]); 
//   FFV2_4_0(w[22], w[47], w[43], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[90]); 
//   FFV2_4_0(w[23], w[2], w[43], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[91]); 
  FFV1_0(w[22], w[52], w[41], pars->Unitary_GC_7, amp[92]); 
  FFV1_0(w[24], w[2], w[41], pars->Unitary_GC_7, amp[93]); 
  FFV2_4_0(w[22], w[52], w[43], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[94]); 
  FFV2_4_0(w[24], w[2], w[43], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[95]); 
//   FFV1_0(w[5], w[47], w[53], pars->Unitary_GC_7, amp[96]); 
//   FFV2_4_0(w[5], w[47], w[54], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[97]); 
//   FFV1_0(w[5], w[9], w[55], pars->Unitary_GC_7, amp[98]); 
//   FFV2_4_0(w[5], w[9], w[56], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[99]); 
//   FFV1_0(w[11], w[4], w[55], pars->Unitary_GC_7, amp[100]); 
//   FFV2_4_0(w[11], w[4], w[56], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[101]); 
//   FFV1_0(w[11], w[2], w[53], pars->Unitary_GC_7, amp[102]); 
//   FFV2_4_0(w[11], w[2], w[54], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[103]); 
  FFV1_0(w[5], w[52], w[53], pars->Unitary_GC_7, amp[104]); 
  FFV2_4_0(w[5], w[52], w[54], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[105]); 
  FFV1_0(w[5], w[17], w[55], pars->Unitary_GC_7, amp[106]); 
  FFV2_4_0(w[5], w[17], w[56], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[107]); 
  FFV1_0(w[18], w[4], w[55], pars->Unitary_GC_7, amp[108]); 
  FFV2_4_0(w[18], w[4], w[56], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[109]); 
  FFV1_0(w[18], w[2], w[53], pars->Unitary_GC_7, amp[110]); 
  FFV2_4_0(w[18], w[2], w[54], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[111]); 
//   FFV1_0(w[33], w[47], w[45], pars->Unitary_GC_7, amp[112]); 
//   FFV1_0(w[34], w[2], w[45], pars->Unitary_GC_7, amp[113]); 
//   FFV2_4_0(w[33], w[47], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[114]); 
//   FFV2_4_0(w[34], w[2], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[115]); 
  FFV1_0(w[33], w[52], w[45], pars->Unitary_GC_7, amp[116]); 
  FFV1_0(w[35], w[2], w[45], pars->Unitary_GC_7, amp[117]); 
  FFV2_4_0(w[33], w[52], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[118]); 
  FFV2_4_0(w[35], w[2], w[46], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[119]); 
//   FFV1_0(w[57], w[47], w[6], pars->Unitary_GC_7, amp[120]); 
//   FFV1_0(w[26], w[58], w[6], pars->Unitary_GC_7, amp[121]); 
//   FFV1_0(w[59], w[47], w[6], pars->Unitary_GC_7, amp[122]); 
//   FFV1_0(w[26], w[60], w[6], pars->Unitary_GC_7, amp[123]); 
  FFV1_0(w[57], w[52], w[6], pars->Unitary_GC_7, amp[124]); 
  FFV1_0(w[32], w[58], w[6], pars->Unitary_GC_7, amp[125]); 
  FFV1_0(w[59], w[52], w[6], pars->Unitary_GC_7, amp[126]); 
  FFV1_0(w[32], w[60], w[6], pars->Unitary_GC_7, amp[127]); 
//   FFV1_0(w[5], w[47], w[61], pars->Unitary_GC_7, amp[128]); 
//   FFV2_4_0(w[5], w[47], w[62], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[129]); 
//   FFV1_0(w[3], w[47], w[63], pars->Unitary_GC_7, amp[130]); 
//   FFV2_4_0(w[3], w[47], w[64], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[131]); 
//   FFV1_0(w[26], w[2], w[63], pars->Unitary_GC_7, amp[132]); 
//   FFV2_4_0(w[26], w[2], w[64], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[133]); 
//   FFV1_0(w[11], w[2], w[61], pars->Unitary_GC_7, amp[134]); 
//   FFV2_4_0(w[11], w[2], w[62], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[135]); 
  FFV1_0(w[5], w[52], w[61], pars->Unitary_GC_7, amp[136]); 
  FFV2_4_0(w[5], w[52], w[62], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[137]); 
  FFV1_0(w[3], w[52], w[63], pars->Unitary_GC_7, amp[138]); 
  FFV2_4_0(w[3], w[52], w[64], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[139]); 
  FFV1_0(w[32], w[2], w[63], pars->Unitary_GC_7, amp[140]); 
  FFV2_4_0(w[32], w[2], w[64], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[141]); 
  FFV1_0(w[18], w[2], w[61], pars->Unitary_GC_7, amp[142]); 
  FFV2_4_0(w[18], w[2], w[62], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[143]); 
//   FFV1_0(w[3], w[47], w[65], pars->Unitary_GC_7, amp[144]); 
//   FFV2_4_0(w[3], w[47], w[66], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[145]); 
//   FFV1_0(w[26], w[4], w[67], pars->Unitary_GC_7, amp[146]); 
//   FFV2_4_0(w[26], w[4], w[68], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[147]); 
//   FFV1_0(w[26], w[2], w[65], pars->Unitary_GC_7, amp[148]); 
//   FFV2_4_0(w[26], w[2], w[66], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[149]); 
//   FFV1_0(w[3], w[9], w[67], pars->Unitary_GC_7, amp[150]); 
//   FFV2_4_0(w[3], w[9], w[68], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[151]); 
  FFV1_0(w[3], w[52], w[65], pars->Unitary_GC_7, amp[152]); 
  FFV2_4_0(w[3], w[52], w[66], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[153]); 
  FFV1_0(w[32], w[4], w[67], pars->Unitary_GC_7, amp[154]); 
  FFV2_4_0(w[32], w[4], w[68], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[155]); 
  FFV1_0(w[32], w[2], w[65], pars->Unitary_GC_7, amp[156]); 
  FFV2_4_0(w[32], w[2], w[66], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[157]); 
  FFV1_0(w[3], w[17], w[67], pars->Unitary_GC_7, amp[158]); 
  FFV2_4_0(w[3], w[17], w[68], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[159]); 
//   FFV1_0(w[5], w[12], w[70], pars->Unitary_GC_7, amp[160]); 
//   FFV2_4_0(w[5], w[12], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[161]); 
//   FFV1_0(w[10], w[4], w[70], pars->Unitary_GC_7, amp[162]); 
//   FFV2_4_0(w[10], w[4], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[163]); 
//   FFV1_0(w[5], w[15], w[70], pars->Unitary_GC_7, amp[164]); 
//   FFV2_4_0(w[5], w[15], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[165]); 
//   FFV1_0(w[14], w[4], w[70], pars->Unitary_GC_7, amp[166]); 
//   FFV2_4_0(w[14], w[4], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[167]); 
//   FFV1_0(w[72], w[1], w[45], pars->Unitary_GC_6, amp[168]); 
//   FFV1_0(w[73], w[1], w[8], pars->Unitary_GC_6, amp[169]); 
//   FFV2_5_0(w[72], w[1], w[46], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[170]); 
//   FFV1_0(w[74], w[1], w[8], pars->Unitary_GC_6, amp[171]); 
//   FFV1_0(w[75], w[1], w[45], pars->Unitary_GC_6, amp[172]); 
//   FFV2_5_0(w[73], w[1], w[13], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[173]); 
//   FFV2_5_0(w[75], w[1], w[46], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[174]); 
//   FFV2_5_0(w[74], w[1], w[13], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[175]); 
//   FFV1_0(w[28], w[4], w[70], pars->Unitary_GC_7, amp[176]); 
//   FFV2_4_0(w[28], w[4], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[177]); 
//   FFV1_0(w[3], w[27], w[70], pars->Unitary_GC_7, amp[178]); 
//   FFV2_4_0(w[3], w[27], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[179]); 
//   FFV1_0(w[31], w[4], w[70], pars->Unitary_GC_7, amp[180]); 
//   FFV2_4_0(w[31], w[4], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[181]); 
//   FFV1_0(w[3], w[30], w[70], pars->Unitary_GC_7, amp[182]); 
//   FFV2_4_0(w[3], w[30], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[183]); 
//   FFV1_0(w[76], w[1], w[41], pars->Unitary_GC_6, amp[184]); 
//   FFV1_0(w[77], w[1], w[25], pars->Unitary_GC_6, amp[185]); 
//   FFV2_5_0(w[76], w[1], w[43], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[186]); 
//   FFV1_0(w[78], w[1], w[25], pars->Unitary_GC_6, amp[187]); 
//   FFV1_0(w[79], w[1], w[41], pars->Unitary_GC_6, amp[188]); 
//   FFV2_5_0(w[77], w[1], w[29], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[189]); 
//   FFV2_5_0(w[79], w[1], w[43], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[190]); 
//   FFV2_5_0(w[78], w[1], w[29], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[191]); 
//   FFV1_0(w[5], w[49], w[70], pars->Unitary_GC_7, amp[192]); 
//   FFV2_4_0(w[5], w[49], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[193]); 
//   FFV1_0(w[48], w[2], w[70], pars->Unitary_GC_7, amp[194]); 
//   FFV2_4_0(w[48], w[2], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[195]); 
//   FFV1_0(w[5], w[51], w[70], pars->Unitary_GC_7, amp[196]); 
//   FFV2_4_0(w[5], w[51], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[197]); 
//   FFV1_0(w[50], w[2], w[70], pars->Unitary_GC_7, amp[198]); 
//   FFV2_4_0(w[50], w[2], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[199]); 
//   FFV1_0(w[3], w[58], w[70], pars->Unitary_GC_7, amp[200]); 
//   FFV2_4_0(w[3], w[58], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[201]); 
//   FFV1_0(w[57], w[2], w[70], pars->Unitary_GC_7, amp[202]); 
//   FFV2_4_0(w[57], w[2], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[203]); 
//   FFV1_0(w[3], w[60], w[70], pars->Unitary_GC_7, amp[204]); 
//   FFV2_4_0(w[3], w[60], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[205]); 
//   FFV1_0(w[59], w[2], w[70], pars->Unitary_GC_7, amp[206]); 
//   FFV2_4_0(w[59], w[2], w[71], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[207]); 
//   FFV1_0(w[5], w[12], w[81], pars->Unitary_GC_7, amp[208]); 
//   FFV2_4_0(w[5], w[12], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[209]); 
//   FFV1_0(w[10], w[4], w[81], pars->Unitary_GC_7, amp[210]); 
//   FFV2_4_0(w[10], w[4], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[211]); 
//   FFV1_0(w[5], w[15], w[81], pars->Unitary_GC_7, amp[212]); 
//   FFV2_4_0(w[5], w[15], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[213]); 
//   FFV1_0(w[14], w[4], w[81], pars->Unitary_GC_7, amp[214]); 
//   FFV2_4_0(w[14], w[4], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[215]); 
//   FFV1_0(w[83], w[80], w[45], pars->Unitary_GC_6, amp[216]); 
//   FFV1_0(w[84], w[80], w[8], pars->Unitary_GC_6, amp[217]); 
//   FFV2_5_0(w[83], w[80], w[46], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[218]); 
//   FFV1_0(w[85], w[80], w[8], pars->Unitary_GC_6, amp[219]); 
//   FFV1_0(w[86], w[80], w[45], pars->Unitary_GC_6, amp[220]); 
//   FFV2_5_0(w[84], w[80], w[13], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[221]); 
//   FFV2_5_0(w[86], w[80], w[46], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[222]); 
//   FFV2_5_0(w[85], w[80], w[13], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[223]); 
//   FFV1_0(w[28], w[4], w[81], pars->Unitary_GC_7, amp[224]); 
//   FFV2_4_0(w[28], w[4], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[225]); 
//   FFV1_0(w[3], w[27], w[81], pars->Unitary_GC_7, amp[226]); 
//   FFV2_4_0(w[3], w[27], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[227]); 
//   FFV1_0(w[31], w[4], w[81], pars->Unitary_GC_7, amp[228]); 
//   FFV2_4_0(w[31], w[4], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[229]); 
//   FFV1_0(w[3], w[30], w[81], pars->Unitary_GC_7, amp[230]); 
//   FFV2_4_0(w[3], w[30], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[231]); 
//   FFV1_0(w[87], w[80], w[41], pars->Unitary_GC_6, amp[232]); 
//   FFV1_0(w[88], w[80], w[25], pars->Unitary_GC_6, amp[233]); 
//   FFV2_5_0(w[87], w[80], w[43], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[234]); 
//   FFV1_0(w[89], w[80], w[25], pars->Unitary_GC_6, amp[235]); 
//   FFV1_0(w[90], w[80], w[41], pars->Unitary_GC_6, amp[236]); 
//   FFV2_5_0(w[88], w[80], w[29], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[237]); 
//   FFV2_5_0(w[90], w[80], w[43], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[238]); 
//   FFV2_5_0(w[89], w[80], w[29], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[239]); 
//   FFV1_0(w[5], w[49], w[81], pars->Unitary_GC_7, amp[240]); 
//   FFV2_4_0(w[5], w[49], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[241]); 
//   FFV1_0(w[48], w[2], w[81], pars->Unitary_GC_7, amp[242]); 
//   FFV2_4_0(w[48], w[2], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[243]); 
//   FFV1_0(w[5], w[51], w[81], pars->Unitary_GC_7, amp[244]); 
//   FFV2_4_0(w[5], w[51], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[245]); 
//   FFV1_0(w[50], w[2], w[81], pars->Unitary_GC_7, amp[246]); 
//   FFV2_4_0(w[50], w[2], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[247]); 
//   FFV1_0(w[3], w[58], w[81], pars->Unitary_GC_7, amp[248]); 
//   FFV2_4_0(w[3], w[58], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[249]); 
//   FFV1_0(w[57], w[2], w[81], pars->Unitary_GC_7, amp[250]); 
//   FFV2_4_0(w[57], w[2], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[251]); 
//   FFV1_0(w[3], w[60], w[81], pars->Unitary_GC_7, amp[252]); 
//   FFV2_4_0(w[3], w[60], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[253]); 
//   FFV1_0(w[59], w[2], w[81], pars->Unitary_GC_7, amp[254]); 
//   FFV2_4_0(w[59], w[2], w[82], pars->Unitary_GC_70, pars->Unitary_GC_75, amp[255]); 
//   FFV1_0(w[83], w[91], w[6], pars->Unitary_GC_6, amp[256]); 
//   FFV1_0(w[84], w[92], w[6], pars->Unitary_GC_6, amp[257]); 
//   FFV1_0(w[83], w[93], w[6], pars->Unitary_GC_6, amp[258]); 
//   FFV1_0(w[85], w[92], w[6], pars->Unitary_GC_6, amp[259]); 
//   FFV1_0(w[86], w[91], w[6], pars->Unitary_GC_6, amp[260]); 
//   FFV1_0(w[84], w[94], w[6], pars->Unitary_GC_6, amp[261]); 
//   FFV1_0(w[86], w[93], w[6], pars->Unitary_GC_6, amp[262]); 
//   FFV1_0(w[85], w[94], w[6], pars->Unitary_GC_6, amp[263]); 
//   FFV1_0(w[83], w[1], w[63], pars->Unitary_GC_6, amp[264]); 
//   FFV2_5_0(w[83], w[1], w[64], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[265]); 
//   FFV1_0(w[0], w[92], w[63], pars->Unitary_GC_6, amp[266]); 
//   FFV2_5_0(w[0], w[92], w[64], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[267]); 
//   FFV1_0(w[86], w[1], w[63], pars->Unitary_GC_6, amp[268]); 
//   FFV2_5_0(w[86], w[1], w[64], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[269]); 
//   FFV1_0(w[0], w[94], w[63], pars->Unitary_GC_6, amp[270]); 
//   FFV2_5_0(w[0], w[94], w[64], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[271]); 
//   FFV1_0(w[83], w[1], w[65], pars->Unitary_GC_6, amp[272]); 
//   FFV2_5_0(w[83], w[1], w[66], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[273]); 
//   FFV1_0(w[0], w[92], w[65], pars->Unitary_GC_6, amp[274]); 
//   FFV2_5_0(w[0], w[92], w[66], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[275]); 
//   FFV1_0(w[86], w[1], w[65], pars->Unitary_GC_6, amp[276]); 
//   FFV2_5_0(w[86], w[1], w[66], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[277]); 
//   FFV1_0(w[0], w[94], w[65], pars->Unitary_GC_6, amp[278]); 
//   FFV2_5_0(w[0], w[94], w[66], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[279]); 
//   FFV1_0(w[87], w[95], w[6], pars->Unitary_GC_6, amp[280]); 
//   FFV1_0(w[88], w[96], w[6], pars->Unitary_GC_6, amp[281]); 
//   FFV1_0(w[87], w[97], w[6], pars->Unitary_GC_6, amp[282]); 
//   FFV1_0(w[89], w[96], w[6], pars->Unitary_GC_6, amp[283]); 
//   FFV1_0(w[90], w[95], w[6], pars->Unitary_GC_6, amp[284]); 
//   FFV1_0(w[88], w[98], w[6], pars->Unitary_GC_6, amp[285]); 
//   FFV1_0(w[90], w[97], w[6], pars->Unitary_GC_6, amp[286]); 
//   FFV1_0(w[89], w[98], w[6], pars->Unitary_GC_6, amp[287]); 
//   FFV1_0(w[87], w[1], w[53], pars->Unitary_GC_6, amp[288]); 
//   FFV2_5_0(w[87], w[1], w[54], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[289]); 
//   FFV1_0(w[0], w[96], w[53], pars->Unitary_GC_6, amp[290]); 
//   FFV2_5_0(w[0], w[96], w[54], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[291]); 
//   FFV1_0(w[90], w[1], w[53], pars->Unitary_GC_6, amp[292]); 
//   FFV2_5_0(w[90], w[1], w[54], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[293]); 
//   FFV1_0(w[0], w[98], w[53], pars->Unitary_GC_6, amp[294]); 
//   FFV2_5_0(w[0], w[98], w[54], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[295]); 
//   FFV1_0(w[87], w[1], w[61], pars->Unitary_GC_6, amp[296]); 
//   FFV2_5_0(w[87], w[1], w[62], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[297]); 
//   FFV1_0(w[0], w[96], w[61], pars->Unitary_GC_6, amp[298]); 
//   FFV2_5_0(w[0], w[96], w[62], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[299]); 
//   FFV1_0(w[90], w[1], w[61], pars->Unitary_GC_6, amp[300]); 
//   FFV2_5_0(w[90], w[1], w[62], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[301]); 
//   FFV1_0(w[0], w[98], w[61], pars->Unitary_GC_6, amp[302]); 
//   FFV2_5_0(w[0], w[98], w[62], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[303]); 
//   FFV1_0(w[88], w[1], w[37], pars->Unitary_GC_6, amp[304]); 
//   FFV2_5_0(w[88], w[1], w[38], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[305]); 
//   FFV1_0(w[0], w[95], w[37], pars->Unitary_GC_6, amp[306]); 
//   FFV2_5_0(w[0], w[95], w[38], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[307]); 
//   FFV1_0(w[89], w[1], w[37], pars->Unitary_GC_6, amp[308]); 
//   FFV2_5_0(w[89], w[1], w[38], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[309]); 
//   FFV1_0(w[0], w[97], w[37], pars->Unitary_GC_6, amp[310]); 
//   FFV2_5_0(w[0], w[97], w[38], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[311]); 
//   FFV1_0(w[84], w[1], w[39], pars->Unitary_GC_6, amp[312]); 
//   FFV2_5_0(w[84], w[1], w[40], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[313]); 
//   FFV1_0(w[0], w[91], w[39], pars->Unitary_GC_6, amp[314]); 
//   FFV2_5_0(w[0], w[91], w[40], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[315]); 
//   FFV1_0(w[85], w[1], w[39], pars->Unitary_GC_6, amp[316]); 
//   FFV2_5_0(w[85], w[1], w[40], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[317]); 
//   FFV1_0(w[0], w[93], w[39], pars->Unitary_GC_6, amp[318]); 
//   FFV2_5_0(w[0], w[93], w[40], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[319]); 
//   FFV1_0(w[88], w[1], w[67], pars->Unitary_GC_6, amp[320]); 
//   FFV2_5_0(w[88], w[1], w[68], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[321]); 
//   FFV1_0(w[0], w[95], w[67], pars->Unitary_GC_6, amp[322]); 
//   FFV2_5_0(w[0], w[95], w[68], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[323]); 
//   FFV1_0(w[89], w[1], w[67], pars->Unitary_GC_6, amp[324]); 
//   FFV2_5_0(w[89], w[1], w[68], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[325]); 
//   FFV1_0(w[0], w[97], w[67], pars->Unitary_GC_6, amp[326]); 
//   FFV2_5_0(w[0], w[97], w[68], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[327]); 
//   FFV1_0(w[84], w[1], w[55], pars->Unitary_GC_6, amp[328]); 
//   FFV2_5_0(w[84], w[1], w[56], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[329]); 
//   FFV1_0(w[0], w[91], w[55], pars->Unitary_GC_6, amp[330]); 
//   FFV2_5_0(w[0], w[91], w[56], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[331]); 
//   FFV1_0(w[85], w[1], w[55], pars->Unitary_GC_6, amp[332]); 
//   FFV2_5_0(w[85], w[1], w[56], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[333]); 
//   FFV1_0(w[0], w[93], w[55], pars->Unitary_GC_6, amp[334]); 
//   FFV2_5_0(w[0], w[93], w[56], pars->Unitary_GC_71, pars->Unitary_GC_74, amp[335]); 

}


double qq_Z4l_SIG_UP_SFpA::matrix_ccx_mummupmummupa_no_hzpxg() 
{
  int i, j; 
  // Local variables
//   const int ngraphs = 336; 
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
      amp[21] + amp[22] + amp[23] - amp[24] - amp[25] - amp[26] - amp[27] -
      amp[28] - amp[29] - amp[30] - amp[31] - amp[32] - amp[33] - amp[34] -
      amp[35] - amp[36] - amp[37] - amp[38] - amp[39] - amp[40] - amp[41] -
      amp[42] - amp[43] - amp[44] - amp[45] - amp[46] - amp[47] - amp[48] -
      amp[49] + amp[50] + amp[51] - amp[52] - amp[53] + amp[54] + amp[55] -
      amp[56] - amp[57] + amp[58] + amp[59] - amp[60] - amp[61] + amp[62] +
      amp[63] - amp[64] - amp[65] - amp[66] - amp[67] - amp[68] - amp[69] -
      amp[70] - amp[71] + amp[72] + amp[73] + amp[74] + amp[75] + amp[76] +
      amp[77] + amp[78] + amp[79] - amp[80] - amp[81] - amp[82] - amp[83] -
      amp[84] - amp[85] - amp[86] - amp[87] - amp[88] - amp[89] - amp[90] -
      amp[91] - amp[92] - amp[93] - amp[94] - amp[95] - amp[96] - amp[97] +
      amp[98] + amp[99] + amp[100] + amp[101] - amp[102] - amp[103] - amp[104]
      - amp[105] + amp[106] + amp[107] + amp[108] + amp[109] - amp[110] -
      amp[111] + amp[112] + amp[113] + amp[114] + amp[115] + amp[116] +
      amp[117] + amp[118] + amp[119] + amp[120] + amp[121] + amp[122] +
      amp[123] + amp[124] + amp[125] + amp[126] + amp[127] - amp[128] -
      amp[129] + amp[130] + amp[131] + amp[132] + amp[133] - amp[134] -
      amp[135] - amp[136] - amp[137] + amp[138] + amp[139] + amp[140] +
      amp[141] - amp[142] - amp[143] + amp[144] + amp[145] - amp[146] -
      amp[147] + amp[148] + amp[149] - amp[150] - amp[151] + amp[152] +
      amp[153] - amp[154] - amp[155] + amp[156] + amp[157] - amp[158] -
      amp[159] + amp[160] + amp[161] + amp[162] + amp[163] + amp[164] +
      amp[165] + amp[166] + amp[167] + amp[168] + amp[169] + amp[170] +
      amp[171] + amp[172] + amp[173] + amp[174] + amp[175] - amp[176] -
      amp[177] - amp[178] - amp[179] - amp[180] - amp[181] - amp[182] -
      amp[183] - amp[184] - amp[185] - amp[186] - amp[187] - amp[188] -
      amp[189] - amp[190] - amp[191] - amp[192] - amp[193] - amp[194] -
      amp[195] - amp[196] - amp[197] - amp[198] - amp[199] + amp[200] +
      amp[201] + amp[202] + amp[203] + amp[204] + amp[205] + amp[206] +
      amp[207] + amp[208] + amp[209] + amp[210] + amp[211] + amp[212] +
      amp[213] + amp[214] + amp[215] + amp[216] + amp[217] + amp[218] +
      amp[219] + amp[220] + amp[221] + amp[222] + amp[223] - amp[224] -
      amp[225] - amp[226] - amp[227] - amp[228] - amp[229] - amp[230] -
      amp[231] - amp[232] - amp[233] - amp[234] - amp[235] - amp[236] -
      amp[237] - amp[238] - amp[239] - amp[240] - amp[241] - amp[242] -
      amp[243] - amp[244] - amp[245] - amp[246] - amp[247] + amp[248] +
      amp[249] + amp[250] + amp[251] + amp[252] + amp[253] + amp[254] +
      amp[255] + amp[256] + amp[257] + amp[258] + amp[259] + amp[260] +
      amp[261] + amp[262] + amp[263] + amp[264] + amp[265] + amp[266] +
      amp[267] + amp[268] + amp[269] + amp[270] + amp[271] + amp[272] +
      amp[273] + amp[274] + amp[275] + amp[276] + amp[277] + amp[278] +
      amp[279] - amp[280] - amp[281] - amp[282] - amp[283] - amp[284] -
      amp[285] - amp[286] - amp[287] - amp[288] - amp[289] - amp[290] -
      amp[291] - amp[292] - amp[293] - amp[294] - amp[295] - amp[296] -
      amp[297] - amp[298] - amp[299] - amp[300] - amp[301] - amp[302] -
      amp[303] - amp[304] - amp[305] - amp[306] - amp[307] - amp[308] -
      amp[309] - amp[310] - amp[311] + amp[312] + amp[313] + amp[314] +
      amp[315] + amp[316] + amp[317] + amp[318] + amp[319] - amp[320] -
      amp[321] - amp[322] - amp[323] - amp[324] - amp[325] - amp[326] -
      amp[327] + amp[328] + amp[329] + amp[330] + amp[331] + amp[332] +
      amp[333] + amp[334] + amp[335];

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



