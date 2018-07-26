//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.2, 2014-02-07
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "Spin1_2f_SF.h"
#include "../HelAmps_HEF_MEKD2_1.h"	// Changed by Convert_source 0.2

using namespace MG5_HEF_MEKD2_1;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: zp > mu- mu+ mu- mu+ S1QQ=1 / h xg

//--------------------------------------------------------------------------
// Initialize process.

void Spin1_2f_SF::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MEKD::getInstance();	// Changed by Convert_source 0.2 
  SLHAReader_MEKD slha(param_card_name);
  pars->setIndependentParameters(slha);
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// Moved here by Convert_source 0.2
  // Set external particle masses for this matrix element
  mME.push_back(pars->MZp);
  mME.push_back(pars->MM);
  mME.push_back(pars->MM);
  mME.push_back(pars->MM);
  mME.push_back(pars->MM);
  jamp2[0] = new double[1];
	for( int count=0; count<namplitudes; count++ ) amp[count] = 0;
}

//--------------------------------------------------------------------------
// Update process.	// Created here by Convert_source 0.2

void Spin1_2f_SF::updateProc(SLHAReader_MEKD &slha)
{
	pars->setIndependentParameters(slha);
	pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// needed when altering couplings
	
	// Set external particle masses for this matrix element
	// Should correspond to initProc
	mME[0]=(pars->MZp);
	mME[1]=(pars->MM);
	mME[2]=(pars->MM);
	mME[3]=(pars->MM);
	mME[4]=(pars->MM);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void Spin1_2f_SF::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters();
  pars->setDependentCouplings();
	// Deleted by Convert_source 0.2
	
  // Reset color flows
  for(int i = 0;i < 1;i++ )
    jamp2[0][i] = 0.;

  // Local variables and constants
  const int ncomb = 48;
  static bool goodhel[ncomb] = {ncomb * false};
//	static int ntry = 0, sum_hel = 0, ngood = 0;	// Moved by Convert_source 0.2
  static int igood[ncomb];
  static int jhel;
//	std::complex<double> * * wfs;	// Changed by Convert_source 0.2
  double t[nprocesses];
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1}, {-1,
      -1, -1, -1, 1}, {-1, -1, -1, 1, -1}, {-1, -1, -1, 1, 1}, {-1, -1, 1, -1,
      -1}, {-1, -1, 1, -1, 1}, {-1, -1, 1, 1, -1}, {-1, -1, 1, 1, 1}, {-1, 1,
      -1, -1, -1}, {-1, 1, -1, -1, 1}, {-1, 1, -1, 1, -1}, {-1, 1, -1, 1, 1},
      {-1, 1, 1, -1, -1}, {-1, 1, 1, -1, 1}, {-1, 1, 1, 1, -1}, {-1, 1, 1, 1,
      1}, {0, -1, -1, -1, -1}, {0, -1, -1, -1, 1}, {0, -1, -1, 1, -1}, {0, -1,
      -1, 1, 1}, {0, -1, 1, -1, -1}, {0, -1, 1, -1, 1}, {0, -1, 1, 1, -1}, {0,
      -1, 1, 1, 1}, {0, 1, -1, -1, -1}, {0, 1, -1, -1, 1}, {0, 1, -1, 1, -1},
      {0, 1, -1, 1, 1}, {0, 1, 1, -1, -1}, {0, 1, 1, -1, 1}, {0, 1, 1, 1, -1},
      {0, 1, 1, 1, 1}, {1, -1, -1, -1, -1}, {1, -1, -1, -1, 1}, {1, -1, -1, 1,
      -1}, {1, -1, -1, 1, 1}, {1, -1, 1, -1, -1}, {1, -1, 1, -1, 1}, {1, -1, 1,
      1, -1}, {1, -1, 1, 1, 1}, {1, 1, -1, -1, -1}, {1, 1, -1, -1, 1}, {1, 1,
      -1, 1, -1}, {1, 1, -1, 1, 1}, {1, 1, 1, -1, -1}, {1, 1, 1, -1, 1}, {1, 1,
      1, 1, -1}, {1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {12};

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
        t[0] = matrix_zp_mummupmummup_no_hxg();

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
      t[0] = matrix_zp_mummupmummup_no_hxg();

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

double Spin1_2f_SF::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 300 && id2 == 13)
  {
    // Add matrix elements for processes with beams (300, 13)
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

void Spin1_2f_SF::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
//	int i, j;	// Changed by Convert_source 0.2

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
  oxxxxx(p[perm[1]], mME[1], hel[1], +1, w[1]);
  ixxxxx(p[perm[2]], mME[2], hel[2], -1, w[2]);
  oxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]);
  ixxxxx(p[perm[4]], mME[4], hel[4], -1, w[4]);
  FFV2P0_3(w[2], w[1], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[5]);
  FFV2_1(w[3], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[6]);
  FFV2_2(w[4], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[7]);
  FFV5_7_3(w[2], w[1], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[8]);
  FFV5_7_1(w[3], w[8], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO, w[9]);
  FFV5_7_2(w[4], w[8], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[10]);
  FFV5_7_3(w[4], w[3], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[11]);
  FFV2P0_3(w[4], w[1], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[12]);
  FFV2_2(w[2], w[12], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[13]);
  FFV2_1(w[3], w[12], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[14]);
  FFV5_7_3(w[4], w[1], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[15]);
  FFV5_7_2(w[2], w[15], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[16]);
  FFV5_7_1(w[3], w[15], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[17]);
  FFV5_7_3(w[2], w[3], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[18]);
  FFV2P0_3(w[2], w[3], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[19]);
  FFV2_1(w[1], w[19], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[20]);
  FFV2_2(w[4], w[19], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[21]);
  FFV5_7_1(w[1], w[18], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[22]);
  FFV5_7_2(w[4], w[18], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[23]);
  FFV2P0_3(w[4], w[3], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[24]);
  FFV2_1(w[1], w[24], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[25]);
  FFV2_2(w[2], w[24], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[26]);
  FFV5_7_1(w[1], w[11], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[27]);
  FFV5_7_2(w[2], w[11], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[28]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_2_3_4_0(w[4], w[6], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[0]);
  FFV1_2_3_4_0(w[7], w[3], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[1]);
  FFV1_2_3_4_0(w[4], w[9], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[2]);
  FFV1_2_3_4_0(w[10], w[3], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[3]);
//   VVV1_2_0(w[8], w[11], w[0], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[4]);	// ZZ
  FFV1_2_3_4_0(w[13], w[3], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[5]);
  FFV1_2_3_4_0(w[2], w[14], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[6]);
  FFV1_2_3_4_0(w[16], w[3], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[7]);
  FFV1_2_3_4_0(w[2], w[17], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[8]);
//   VVV1_2_0(w[15], w[18], w[0], pars->HEF_MEKD2_1_GC_2, pars->HEF_MEKD2_1_GC_1, amp[9]);	// ZZ
  FFV1_2_3_4_0(w[4], w[20], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[10]);
  FFV1_2_3_4_0(w[21], w[1], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[11]);
  FFV1_2_3_4_0(w[4], w[22], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[12]);
  FFV1_2_3_4_0(w[23], w[1], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[13]);
  FFV1_2_3_4_0(w[2], w[25], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[14]);
  FFV1_2_3_4_0(w[26], w[1], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[15]);
  FFV1_2_3_4_0(w[2], w[27], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[16]);
  FFV1_2_3_4_0(w[28], w[1], w[0], pars->HEF_MEKD2_1_GC_141, pars->HEF_MEKD2_1_GC_138, pars->HEF_MEKD2_1_GC_139,
      pars->HEF_MEKD2_1_GC_140, amp[17]);

}
double Spin1_2f_SF::matrix_zp_mummupmummup_no_hxg() 
{
  int i, j;
  // Local variables
	// Commented out by Convert_source 0.2
//  const int ngraphs = 18;
  const int ncolor = 1;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[1] = {1.};
  static const double cf[1][1] = {{1.}};

  // Calculate color flows
  jamp[0] = +amp[0] + amp[1] + amp[2] + amp[3] + amp[4] - amp[5] - amp[6] -
      amp[7] - amp[8] - amp[9] - amp[10] - amp[11] - amp[12] - amp[13] +
      amp[14] + amp[15] + amp[16] + amp[17];

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



