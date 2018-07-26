//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.9, 2013-04-01
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "Spin2_OFpA.h"
#include "../HelAmps_HEF_MEKD.h"

using namespace MG5_HEF_MEKD; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: xg > e- e+ mu- mu+ a

//--------------------------------------------------------------------------
// Initialize process.

void Spin2_OFpA::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MEKD::getInstance(); 
  SLHAReader_MEKD slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// moved here by Ghost remover v. 0.1 


  // Set external particle masses for this matrix element
  mME.push_back(pars->MXG); 
  mME.push_back(pars->Me); 
  mME.push_back(pars->Me); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[1]; 
}

//--------------------------------------------------------------------------
// Update process.

void Spin2_OFpA::updateProc(SLHAReader_MEKD &slha) 
{
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// moved here by Ghost remover v. 0.1 
  
  // Set external particle masses for this matrix element
  mME[0]=(pars->MXG);
  mME[1]=(pars->Me);
  mME[2]=(pars->Me);
  mME[3]=(pars->MM);
  mME[4]=(pars->MM);
  mME[5]=(pars->ZERO);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void Spin2_OFpA::sigmaKin() 
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
  const int ncomb = 160; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  //static int ntry = 0, sum_hel = 0, ngood = 0;	// picked out by Ghost remover v. 0.1 
  static int igood[ncomb]; 
  static int jhel; 
//  std::complex<double> * * wfs;
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-2, -1, -1, -1, -1, -1},
      {-2, -1, -1, -1, -1, 1}, {-2, -1, -1, -1, 1, -1}, {-2, -1, -1, -1, 1, 1},
      {-2, -1, -1, 1, -1, -1}, {-2, -1, -1, 1, -1, 1}, {-2, -1, -1, 1, 1, -1},
      {-2, -1, -1, 1, 1, 1}, {-2, -1, 1, -1, -1, -1}, {-2, -1, 1, -1, -1, 1},
      {-2, -1, 1, -1, 1, -1}, {-2, -1, 1, -1, 1, 1}, {-2, -1, 1, 1, -1, -1},
      {-2, -1, 1, 1, -1, 1}, {-2, -1, 1, 1, 1, -1}, {-2, -1, 1, 1, 1, 1}, {-2,
      1, -1, -1, -1, -1}, {-2, 1, -1, -1, -1, 1}, {-2, 1, -1, -1, 1, -1}, {-2,
      1, -1, -1, 1, 1}, {-2, 1, -1, 1, -1, -1}, {-2, 1, -1, 1, -1, 1}, {-2, 1,
      -1, 1, 1, -1}, {-2, 1, -1, 1, 1, 1}, {-2, 1, 1, -1, -1, -1}, {-2, 1, 1,
      -1, -1, 1}, {-2, 1, 1, -1, 1, -1}, {-2, 1, 1, -1, 1, 1}, {-2, 1, 1, 1,
      -1, -1}, {-2, 1, 1, 1, -1, 1}, {-2, 1, 1, 1, 1, -1}, {-2, 1, 1, 1, 1, 1},
      {-1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, 1,
      -1}, {-1, -1, -1, -1, 1, 1}, {-1, -1, -1, 1, -1, -1}, {-1, -1, -1, 1, -1,
      1}, {-1, -1, -1, 1, 1, -1}, {-1, -1, -1, 1, 1, 1}, {-1, -1, 1, -1, -1,
      -1}, {-1, -1, 1, -1, -1, 1}, {-1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1,
      1}, {-1, -1, 1, 1, -1, -1}, {-1, -1, 1, 1, -1, 1}, {-1, -1, 1, 1, 1, -1},
      {-1, -1, 1, 1, 1, 1}, {-1, 1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1},
      {-1, 1, -1, -1, 1, -1}, {-1, 1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, -1},
      {-1, 1, -1, 1, -1, 1}, {-1, 1, -1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1}, {-1,
      1, 1, -1, -1, -1}, {-1, 1, 1, -1, -1, 1}, {-1, 1, 1, -1, 1, -1}, {-1, 1,
      1, -1, 1, 1}, {-1, 1, 1, 1, -1, -1}, {-1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1,
      1, -1}, {-1, 1, 1, 1, 1, 1}, {0, -1, -1, -1, -1, -1}, {0, -1, -1, -1, -1,
      1}, {0, -1, -1, -1, 1, -1}, {0, -1, -1, -1, 1, 1}, {0, -1, -1, 1, -1,
      -1}, {0, -1, -1, 1, -1, 1}, {0, -1, -1, 1, 1, -1}, {0, -1, -1, 1, 1, 1},
      {0, -1, 1, -1, -1, -1}, {0, -1, 1, -1, -1, 1}, {0, -1, 1, -1, 1, -1}, {0,
      -1, 1, -1, 1, 1}, {0, -1, 1, 1, -1, -1}, {0, -1, 1, 1, -1, 1}, {0, -1, 1,
      1, 1, -1}, {0, -1, 1, 1, 1, 1}, {0, 1, -1, -1, -1, -1}, {0, 1, -1, -1,
      -1, 1}, {0, 1, -1, -1, 1, -1}, {0, 1, -1, -1, 1, 1}, {0, 1, -1, 1, -1,
      -1}, {0, 1, -1, 1, -1, 1}, {0, 1, -1, 1, 1, -1}, {0, 1, -1, 1, 1, 1}, {0,
      1, 1, -1, -1, -1}, {0, 1, 1, -1, -1, 1}, {0, 1, 1, -1, 1, -1}, {0, 1, 1,
      -1, 1, 1}, {0, 1, 1, 1, -1, -1}, {0, 1, 1, 1, -1, 1}, {0, 1, 1, 1, 1,
      -1}, {0, 1, 1, 1, 1, 1}, {1, -1, -1, -1, -1, -1}, {1, -1, -1, -1, -1, 1},
      {1, -1, -1, -1, 1, -1}, {1, -1, -1, -1, 1, 1}, {1, -1, -1, 1, -1, -1},
      {1, -1, -1, 1, -1, 1}, {1, -1, -1, 1, 1, -1}, {1, -1, -1, 1, 1, 1}, {1,
      -1, 1, -1, -1, -1}, {1, -1, 1, -1, -1, 1}, {1, -1, 1, -1, 1, -1}, {1, -1,
      1, -1, 1, 1}, {1, -1, 1, 1, -1, -1}, {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1,
      1, -1}, {1, -1, 1, 1, 1, 1}, {1, 1, -1, -1, -1, -1}, {1, 1, -1, -1, -1,
      1}, {1, 1, -1, -1, 1, -1}, {1, 1, -1, -1, 1, 1}, {1, 1, -1, 1, -1, -1},
      {1, 1, -1, 1, -1, 1}, {1, 1, -1, 1, 1, -1}, {1, 1, -1, 1, 1, 1}, {1, 1,
      1, -1, -1, -1}, {1, 1, 1, -1, -1, 1}, {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1,
      1, 1}, {1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, -1},
      {1, 1, 1, 1, 1, 1}, {2, -1, -1, -1, -1, -1}, {2, -1, -1, -1, -1, 1}, {2,
      -1, -1, -1, 1, -1}, {2, -1, -1, -1, 1, 1}, {2, -1, -1, 1, -1, -1}, {2,
      -1, -1, 1, -1, 1}, {2, -1, -1, 1, 1, -1}, {2, -1, -1, 1, 1, 1}, {2, -1,
      1, -1, -1, -1}, {2, -1, 1, -1, -1, 1}, {2, -1, 1, -1, 1, -1}, {2, -1, 1,
      -1, 1, 1}, {2, -1, 1, 1, -1, -1}, {2, -1, 1, 1, -1, 1}, {2, -1, 1, 1, 1,
      -1}, {2, -1, 1, 1, 1, 1}, {2, 1, -1, -1, -1, -1}, {2, 1, -1, -1, -1, 1},
      {2, 1, -1, -1, 1, -1}, {2, 1, -1, -1, 1, 1}, {2, 1, -1, 1, -1, -1}, {2,
      1, -1, 1, -1, 1}, {2, 1, -1, 1, 1, -1}, {2, 1, -1, 1, 1, 1}, {2, 1, 1,
      -1, -1, -1}, {2, 1, 1, -1, -1, 1}, {2, 1, 1, -1, 1, -1}, {2, 1, 1, -1, 1,
      1}, {2, 1, 1, 1, -1, -1}, {2, 1, 1, 1, -1, 1}, {2, 1, 1, 1, 1, -1}, {2,
      1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {5}; 

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
        t[0] = matrix_xg_emepmummupa(); 

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
      t[0] = matrix_xg_emepmummupa(); 

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

double Spin2_OFpA::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 9000007 && id2 == 11)
  {
    // Add matrix elements for processes with beams (9000007, 11)
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

void Spin2_OFpA::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
//   int i, j; 

  // Calculate all wavefunctions
  txxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
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
  FFV2_1(w[1], w[5], pars->HEF_MEKD_GC_5, pars->Me, pars->ZERO, w[11]); 
  FFV5_7_3(w[4], w[3], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[12]); 
  FFV5_7_3(w[2], w[11], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[13]); 
  FFV2_2(w[2], w[5], pars->HEF_MEKD_GC_5, pars->Me, pars->ZERO, w[14]); 
  FFV5_7_3(w[14], w[1], pars->HEF_MEKD_GC_161, pars->HEF_MEKD_GC_168, pars->MZ, pars->WZ, w[15]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  VVT10_11_12_13_2_3_6_7_8_9_0(w[6], w[8], w[0], pars->HEF_MEKD_GC_75, pars->HEF_MEKD_GC_81,
      pars->HEF_MEKD_GC_67, pars->HEF_MEKD_GC_71, pars->HEF_MEKD_GC_90, pars->HEF_MEKD_GC_92, pars->HEF_MEKD_GC_63,
      pars->HEF_MEKD_GC_85, pars->HEF_MEKD_GC_79, pars->HEF_MEKD_GC_83, amp[0]);
  VVT10_11_12_13_2_3_6_7_8_9_0(w[6], w[10], w[0], pars->HEF_MEKD_GC_75, pars->HEF_MEKD_GC_81,
      pars->HEF_MEKD_GC_67, pars->HEF_MEKD_GC_71, pars->HEF_MEKD_GC_90, pars->HEF_MEKD_GC_92, pars->HEF_MEKD_GC_63,
      pars->HEF_MEKD_GC_85, pars->HEF_MEKD_GC_79, pars->HEF_MEKD_GC_83, amp[1]);
  VVT10_11_12_13_2_3_6_7_8_9_0(w[13], w[12], w[0], pars->HEF_MEKD_GC_75, pars->HEF_MEKD_GC_81,
      pars->HEF_MEKD_GC_67, pars->HEF_MEKD_GC_71, pars->HEF_MEKD_GC_90, pars->HEF_MEKD_GC_92, pars->HEF_MEKD_GC_63,
      pars->HEF_MEKD_GC_85, pars->HEF_MEKD_GC_79, pars->HEF_MEKD_GC_83, amp[2]);
  VVT10_11_12_13_2_3_6_7_8_9_0(w[15], w[12], w[0], pars->HEF_MEKD_GC_75, pars->HEF_MEKD_GC_81,
      pars->HEF_MEKD_GC_67, pars->HEF_MEKD_GC_71, pars->HEF_MEKD_GC_90, pars->HEF_MEKD_GC_92, pars->HEF_MEKD_GC_63,
      pars->HEF_MEKD_GC_85, pars->HEF_MEKD_GC_79, pars->HEF_MEKD_GC_83, amp[3]);

}
double Spin2_OFpA::matrix_xg_emepmummupa() 
{
  int i, j; 
  // Local variables
//   const int ngraphs = 4; 
  const int ncolor = 1; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[1] = {1.}; 
  static const double cf[1][1] = {{1.}}; 

  // Calculate color flows
  jamp[0] = +amp[0] + amp[1] + amp[2] + amp[3]; 

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



