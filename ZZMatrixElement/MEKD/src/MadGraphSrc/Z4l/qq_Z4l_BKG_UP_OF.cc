//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.3, 2012-11-01
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "qq_Z4l_BKG_UP_OF.h"
#include "../HelAmps_HEF_UFO_bkg.h"
#include "../read_slha.h"

using namespace MG5_HEF_UFO_bkg;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: c c~ > e- e+ mu- mu+ / h xg WEIGHTED=8

//--------------------------------------------------------------------------
// Initialize process.

void qq_Z4l_BKG_UP_OF::initProc(string param_card_name) 
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
  mME.push_back(pars->Me); 
  mME.push_back(pars->Me); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->MM); 
  jamp2[0] = new double[1];
	for( int count=0; count<namplitudes; count++ ) amp[count] = 0;
}

//--------------------------------------------------------------------------
// Update process.

void qq_Z4l_BKG_UP_OF::updateProc(SLHAReader_MEKD &slha) 
{
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// moved here by Ghost remover v. 0.1 
  
  // Set external particle masses for this matrix element
  mME[0]=(pars->MC);
  mME[1]=(pars->MC);
  mME[2]=(pars->Me);
  mME[3]=(pars->Me);
  mME[4]=(pars->MM);
  mME[5]=(pars->MM);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void qq_Z4l_BKG_UP_OF::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings();
    
  // Reset color flows
  for(int i = 0; i < 1; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 64; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  //static int ntry = 0, sum_hel = 0, ngood = 0;	// picked out by Ghost remover v. 0.1 
  static int igood[ncomb]; 
  static int jhel; 
//   std::complex<double> * * wfs; 
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
        t[0] = matrix_ccx_emepmummup_no_hxg(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_ccx_emepmummup_no_hxg(); 
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
      t[0] = matrix_ccx_emepmummup_no_hxg(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_ccx_emepmummup_no_hxg(); 
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

double qq_Z4l_BKG_UP_OF::sigmaHat() 
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

void qq_Z4l_BKG_UP_OF::calculate_wavefunctions(const int perm[], const int hel[])
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
  FFV41_3(w[0], w[1], pars->GC_2, pars->ZERO, pars->ZERO, w[6]); 
  FFV41_3(w[3], w[2], pars->GC_3, pars->ZERO, pars->ZERO, w[7]); 
  FFV41_1(w[4], w[6], pars->GC_3, pars->MM, pars->ZERO, w[8]); 
  FFV41_2(w[5], w[6], pars->GC_3, pars->MM, pars->ZERO, w[9]); 
  FFV42_44_3(w[3], w[2], pars->GC_109, pars->GC_116, pars->MZ, pars->WZ,
      w[10]);
//   FFV42_45_3(w[0], w[1], pars->GC_110, pars->GC_115, pars->MZ, pars->WZ,
//       w[11]);
//   FFV42_44_1(w[4], w[11], pars->GC_109, pars->GC_116, pars->MM, pars->ZERO,
//       w[12]);
//   FFV42_44_2(w[5], w[11], pars->GC_109, pars->GC_116, pars->MM, pars->ZERO,
//       w[13]);
  FFV41_3(w[5], w[4], pars->GC_3, pars->ZERO, pars->ZERO, w[14]); 
  FFV41_1(w[2], w[6], pars->GC_3, pars->Me, pars->ZERO, w[15]); 
  FFV41_2(w[3], w[6], pars->GC_3, pars->Me, pars->ZERO, w[16]); 
  FFV42_44_3(w[5], w[4], pars->GC_109, pars->GC_116, pars->MZ, pars->WZ,
      w[17]);
//   FFV42_44_1(w[2], w[11], pars->GC_109, pars->GC_116, pars->Me, pars->ZERO,
//       w[18]);
//   FFV42_44_2(w[3], w[11], pars->GC_109, pars->GC_116, pars->Me, pars->ZERO,
//       w[19]);
  FFV41_2(w[0], w[7], pars->GC_2, pars->MC, pars->ZERO, w[20]); 
  FFV41_2(w[0], w[14], pars->GC_2, pars->MC, pars->ZERO, w[21]); 
  FFV42_45_2(w[0], w[17], pars->GC_110, pars->GC_115, pars->MC, pars->ZERO,
      w[22]);
  FFV42_45_2(w[0], w[10], pars->GC_110, pars->GC_115, pars->MC, pars->ZERO,
      w[23]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV41_0(w[5], w[8], w[7], pars->GC_3, amp[0]); 
  FFV41_0(w[9], w[4], w[7], pars->GC_3, amp[1]); 
  FFV42_44_0(w[5], w[8], w[10], pars->GC_109, pars->GC_116, amp[2]); 
  FFV42_44_0(w[9], w[4], w[10], pars->GC_109, pars->GC_116, amp[3]); 
//   FFV41_0(w[5], w[12], w[7], pars->GC_3, amp[4]); 
//   FFV41_0(w[13], w[4], w[7], pars->GC_3, amp[5]); 
//   FFV42_44_0(w[5], w[12], w[10], pars->GC_109, pars->GC_116, amp[6]); 
//   FFV42_44_0(w[13], w[4], w[10], pars->GC_109, pars->GC_116, amp[7]); 
  FFV41_0(w[3], w[15], w[14], pars->GC_3, amp[8]); 
  FFV41_0(w[16], w[2], w[14], pars->GC_3, amp[9]); 
  FFV42_44_0(w[3], w[15], w[17], pars->GC_109, pars->GC_116, amp[10]); 
  FFV42_44_0(w[16], w[2], w[17], pars->GC_109, pars->GC_116, amp[11]); 
//   FFV41_0(w[3], w[18], w[14], pars->GC_3, amp[12]); 
//   FFV41_0(w[19], w[2], w[14], pars->GC_3, amp[13]); 
//   FFV42_44_0(w[3], w[18], w[17], pars->GC_109, pars->GC_116, amp[14]); 
//   FFV42_44_0(w[19], w[2], w[17], pars->GC_109, pars->GC_116, amp[15]); 
  FFV41_0(w[20], w[1], w[14], pars->GC_2, amp[16]); 
  FFV41_0(w[21], w[1], w[7], pars->GC_2, amp[17]); 
  FFV42_45_0(w[20], w[1], w[17], pars->GC_110, pars->GC_115, amp[18]); 
  FFV41_0(w[22], w[1], w[7], pars->GC_2, amp[19]); 
  FFV41_0(w[23], w[1], w[14], pars->GC_2, amp[20]); 
  FFV42_45_0(w[21], w[1], w[10], pars->GC_110, pars->GC_115, amp[21]); 
  FFV42_45_0(w[23], w[1], w[17], pars->GC_110, pars->GC_115, amp[22]); 
  FFV42_45_0(w[22], w[1], w[10], pars->GC_110, pars->GC_115, amp[23]); 

}

double qq_Z4l_BKG_UP_OF::matrix_ccx_emepmummup_no_hxg() 
{
  int i, j; 
  // Local variables
//   const int ngraphs = 24; 
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
      amp[21] + amp[22] + amp[23];

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



