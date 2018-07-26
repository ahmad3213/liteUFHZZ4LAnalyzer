//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.2, 2014-02-07
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "Spin0_2f_SFpA.h"
#include "../HelAmps_HEF_MEKD2_1.h"	// Changed by Convert_source 0.2

using namespace MG5_HEF_MEKD2_1;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: h > mu- mu+ mu- mu+ a S0QQ=1 / zp xg

//--------------------------------------------------------------------------
// Initialize process.

void Spin0_2f_SFpA::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MEKD::getInstance();	// Changed by Convert_source 0.2 
  SLHAReader_MEKD slha(param_card_name);
  pars->setIndependentParameters(slha);
  pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// Moved here by Convert_source 0.2
  // Set external particle masses for this matrix element
  mME.push_back(pars->MH);
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

void Spin0_2f_SFpA::updateProc(SLHAReader_MEKD &slha)
{
	pars->setIndependentParameters(slha);
	pars->setIndependentCouplings();
	ntry = 0, sum_hel = 0, ngood = 0;	// needed when altering couplings
	
	// Set external particle masses for this matrix element
	// Should correspond to initProc
	mME[0]=(pars->MH);
	mME[1]=(pars->MM);
	mME[2]=(pars->MM);
	mME[3]=(pars->MM);
	mME[4]=(pars->MM);
	mME[5]=(pars->ZERO);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void Spin0_2f_SFpA::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters();
  pars->setDependentCouplings();
	// Deleted by Convert_source 0.2
	
  // Reset color flows
  for(int i = 0;i < 1;i++ )
    jamp2[0][i] = 0.;

  // Local variables and constants
  const int ncomb = 32;
  static bool goodhel[ncomb] = {ncomb * false};
//	static int ntry = 0, sum_hel = 0, ngood = 0;	// Moved by Convert_source 0.2
  static int igood[ncomb];
  static int jhel;
//	std::complex<double> * * wfs;	// Changed by Convert_source 0.2
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
        t[0] = matrix_h_mummupmummupa_no_zpxg();

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
      t[0] = matrix_h_mummupmummupa_no_zpxg();

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

double Spin0_2f_SFpA::sigmaHat() 
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

void Spin0_2f_SFpA::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
//	int i, j;	// Changed by Convert_source 0.2

  // Calculate all wavefunctions
  sxxxxx(p[perm[0]], -1, w[0]);
  oxxxxx(p[perm[1]], mME[1], hel[1], +1, w[1]);
  ixxxxx(p[perm[2]], mME[2], hel[2], -1, w[2]);
  oxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]);
  ixxxxx(p[perm[4]], mME[4], hel[4], -1, w[4]);
  vxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]);
  FFV2P0_3(w[2], w[1], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[6]);
  FFV2_1(w[3], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[7]);
  FFV2_1(w[7], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[8]);
  FFV2_2(w[4], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[9]);
  FFV2_2(w[9], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[10]);
  FFV5_7_3(w[2], w[1], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[11]);
  FFV5_7_1(w[3], w[11], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[12]);
  FFV2_1(w[12], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[13]);
  FFV5_7_2(w[4], w[11], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[14]);
  FFV2_2(w[14], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[15]);
  FFV2_1(w[3], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[16]);
  FFV2_1(w[16], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[17]);
  FFV5_7_1(w[16], w[11], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[18]);
  FFV5_7_3(w[4], w[16], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[19]);
  FFV2_2(w[4], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[20]);
  FFV2_2(w[20], w[6], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[21]);
  FFV5_7_2(w[20], w[11], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[22]);
  FFV5_7_3(w[20], w[3], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[23]);
  FFV2P0_3(w[4], w[1], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[24]);
  FFV2_2(w[2], w[24], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[25]);
  FFV2_2(w[25], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[26]);
  FFV2_1(w[3], w[24], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[27]);
  FFV2_1(w[27], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[28]);
  FFV5_7_3(w[4], w[1], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[29]);
  FFV5_7_2(w[2], w[29], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[30]);
  FFV2_2(w[30], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[31]);
  FFV5_7_1(w[3], w[29], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[32]);
  FFV2_1(w[32], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[33]);
  FFV2_2(w[2], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[34]);
  FFV2_2(w[34], w[24], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[35]);
  FFV5_7_2(w[34], w[29], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[36]);
  FFV5_7_3(w[34], w[3], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[37]);
  FFV2_1(w[16], w[24], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[38]);
  FFV5_7_1(w[16], w[29], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[39]);
  FFV5_7_3(w[2], w[16], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[40]);
  FFV2_1(w[1], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[41]);
  FFV2P0_3(w[2], w[41], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[42]);
  FFV2_1(w[3], w[42], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[43]);
  FFV2_2(w[4], w[42], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[44]);
  FFV5_7_3(w[2], w[41], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[45]);
  FFV5_7_1(w[3], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[46]);
  FFV5_7_2(w[4], w[45], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[47]);
  FFV2P0_3(w[4], w[41], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[48]);
  FFV2_2(w[2], w[48], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[49]);
  FFV2_1(w[3], w[48], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[50]);
  FFV5_7_3(w[4], w[41], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[51]);
  FFV5_7_2(w[2], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[52]);
  FFV5_7_1(w[3], w[51], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[53]);
  FFV2P0_3(w[2], w[3], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[54]);
  FFV2_1(w[41], w[54], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[55]);
  FFV2_2(w[4], w[54], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[56]);
  FFV5_7_3(w[2], w[3], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[57]);
  FFV5_7_1(w[41], w[57], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[58]);
  FFV5_7_2(w[4], w[57], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[59]);
  FFV2P0_3(w[4], w[3], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[60]);
  FFV2_1(w[41], w[60], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[61]);
  FFV2_2(w[2], w[60], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[62]);
  FFV5_7_3(w[4], w[3], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[63]);
  FFV5_7_1(w[41], w[63], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[64]);
  FFV5_7_2(w[2], w[63], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[65]);
  FFV2_1(w[1], w[54], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[66]);
  FFV2_1(w[66], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[67]);
  FFV2_2(w[56], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[68]);
  FFV5_7_1(w[1], w[57], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[69]);
  FFV2_1(w[69], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[70]);
  FFV2_2(w[59], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[71]);
  FFV2_2(w[20], w[54], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[72]);
  FFV5_7_3(w[20], w[1], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[73]);
  FFV5_7_2(w[20], w[57], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[74]);
  FFV2P0_3(w[34], w[1], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[75]);
  FFV2_1(w[3], w[75], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[76]);
  FFV2_2(w[4], w[75], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[77]);
  FFV5_7_3(w[34], w[1], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MZ, pars->WZ, w[78]);
  FFV5_7_1(w[3], w[78], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[79]);
  FFV5_7_2(w[4], w[78], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[80]);
  FFV2P0_3(w[34], w[3], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[81]);
  FFV2_1(w[1], w[81], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[82]);
  FFV2_2(w[4], w[81], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[83]);
  FFV5_7_1(w[1], w[37], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[84]);
  FFV5_7_2(w[4], w[37], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[85]);
  FFV2_1(w[1], w[60], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[86]);
  FFV2_2(w[34], w[60], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[87]);
  FFV5_7_1(w[1], w[63], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[88]);
  FFV5_7_2(w[34], w[63], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[89]);
  FFV2_1(w[86], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[90]);
  FFV2_2(w[62], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[91]);
  FFV2_1(w[88], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[92]);
  FFV2_2(w[65], w[5], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[93]);
  FFV2P0_3(w[2], w[16], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[94]);
  FFV2_1(w[1], w[94], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[95]);
  FFV2_2(w[4], w[94], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[96]);
  FFV5_7_1(w[1], w[40], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[97]);
  FFV5_7_2(w[4], w[40], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[98]);
  FFV2P0_3(w[4], w[16], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[99]);
  FFV2_1(w[1], w[99], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[100]);
  FFV2_2(w[2], w[99], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[101]);
  FFV5_7_1(w[1], w[19], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[102]);
  FFV5_7_2(w[2], w[19], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[103]);
  FFV2P0_3(w[20], w[1], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[104]);
  FFV2_2(w[2], w[104], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[105]);
  FFV2_1(w[3], w[104], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[106]);
  FFV5_7_2(w[2], w[73], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[107]);
  FFV5_7_1(w[3], w[73], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[108]);
  FFV2P0_3(w[20], w[3], pars->HEF_MEKD2_1_GC_5, pars->ZERO, pars->ZERO, w[109]);
  FFV2_1(w[1], w[109], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[110]);
  FFV2_2(w[2], w[109], pars->HEF_MEKD2_1_GC_5, pars->MM, pars->ZERO, w[111]);
  FFV5_7_1(w[1], w[23], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[112]);
  FFV5_7_2(w[2], w[23], pars->HEF_MEKD2_1_GC_181, pars->HEF_MEKD2_1_GC_188, pars->MM, pars->ZERO,
      w[113]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFS1_2_0(w[4], w[8], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[0]);
  FFS1_2_0(w[10], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[1]);
  FFS1_2_0(w[4], w[13], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[2]);
  FFS1_2_0(w[15], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[3]);
  FFS1_2_0(w[4], w[17], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[4]);
  FFS1_2_0(w[9], w[16], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[5]);
  FFS1_2_0(w[4], w[18], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[6]);
  FFS1_2_0(w[14], w[16], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[7]);
//   VVS3_4_5_0(w[11], w[19], w[0], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[8]);	// ZZ
//   VVS2_0(w[11], w[19], w[0], pars->HEF_MEKD2_1_GC_25, amp[9]);	// ZZ
  FFS1_2_0(w[20], w[7], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[10]);
  FFS1_2_0(w[21], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[11]);
  FFS1_2_0(w[20], w[12], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[12]);
  FFS1_2_0(w[22], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[13]);
//   VVS3_4_5_0(w[11], w[23], w[0], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[14]);	// ZZ
//   VVS2_0(w[11], w[23], w[0], pars->HEF_MEKD2_1_GC_25, amp[15]);	// ZZ
  FFS1_2_0(w[26], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[16]);
  FFS1_2_0(w[2], w[28], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[17]);
  FFS1_2_0(w[31], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[18]);
  FFS1_2_0(w[2], w[33], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[19]);
  FFS1_2_0(w[35], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[20]);
  FFS1_2_0(w[34], w[27], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[21]);
  FFS1_2_0(w[36], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[22]);
  FFS1_2_0(w[34], w[32], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[23]);
//   VVS3_4_5_0(w[29], w[37], w[0], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[24]);	// ZZ
//   VVS2_0(w[29], w[37], w[0], pars->HEF_MEKD2_1_GC_25, amp[25]);	// ZZ
  FFS1_2_0(w[25], w[16], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[26]);
  FFS1_2_0(w[2], w[38], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[27]);
  FFS1_2_0(w[30], w[16], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[28]);
  FFS1_2_0(w[2], w[39], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[29]);
//   VVS3_4_5_0(w[29], w[40], w[0], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[30]);	// ZZ
//   VVS2_0(w[29], w[40], w[0], pars->HEF_MEKD2_1_GC_25, amp[31]);	// ZZ
  FFS1_2_0(w[4], w[43], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[32]);
  FFS1_2_0(w[44], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[33]);
  FFS1_2_0(w[4], w[46], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[34]);
  FFS1_2_0(w[47], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[35]);
  FFS1_2_0(w[49], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[36]);
  FFS1_2_0(w[2], w[50], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[37]);
  FFS1_2_0(w[52], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[38]);
  FFS1_2_0(w[2], w[53], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[39]);
  FFS1_2_0(w[4], w[55], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[40]);
  FFS1_2_0(w[56], w[41], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[41]);
  FFS1_2_0(w[4], w[58], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[42]);
//   VVS3_4_5_0(w[51], w[57], w[0], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[43]);	// ZZ
//   VVS2_0(w[51], w[57], w[0], pars->HEF_MEKD2_1_GC_25, amp[44]);	// ZZ
  FFS1_2_0(w[59], w[41], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[45]);
  FFS1_2_0(w[2], w[61], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[46]);
  FFS1_2_0(w[62], w[41], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[47]);
//   VVS3_4_5_0(w[45], w[63], w[0], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[48]);	// ZZ
//   VVS2_0(w[45], w[63], w[0], pars->HEF_MEKD2_1_GC_25, amp[49]);	// ZZ
  FFS1_2_0(w[2], w[64], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[50]);
  FFS1_2_0(w[65], w[41], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[51]);
  FFS1_2_0(w[4], w[67], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[52]);
  FFS1_2_0(w[68], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[53]);
  FFS1_2_0(w[4], w[70], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[54]);
  FFS1_2_0(w[71], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[55]);
  FFS1_2_0(w[20], w[66], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[56]);
  FFS1_2_0(w[72], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[57]);
  FFS1_2_0(w[20], w[69], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[58]);
//   VVS3_4_5_0(w[73], w[57], w[0], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[59]);	// ZZ
//   VVS2_0(w[73], w[57], w[0], pars->HEF_MEKD2_1_GC_25, amp[60]);	// ZZ
  FFS1_2_0(w[74], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[61]);
  FFS1_2_0(w[4], w[76], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[62]);
  FFS1_2_0(w[77], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[63]);
  FFS1_2_0(w[4], w[79], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[64]);
  FFS1_2_0(w[80], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[65]);
  FFS1_2_0(w[4], w[82], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[66]);
  FFS1_2_0(w[83], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[67]);
  FFS1_2_0(w[4], w[84], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[68]);
  FFS1_2_0(w[85], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[69]);
  FFS1_2_0(w[34], w[86], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[70]);
  FFS1_2_0(w[87], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[71]);
//   VVS3_4_5_0(w[78], w[63], w[0], pars->HEF_MEKD2_1_GC_14, pars->HEF_MEKD2_1_GC_18, pars->HEF_MEKD2_1_GC_22,
//       amp[72]);	// ZZ
//   VVS2_0(w[78], w[63], w[0], pars->HEF_MEKD2_1_GC_25, amp[73]);	// ZZ
  FFS1_2_0(w[34], w[88], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[74]);
  FFS1_2_0(w[89], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[75]);
  FFS1_2_0(w[2], w[90], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[76]);
  FFS1_2_0(w[91], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[77]);
  FFS1_2_0(w[2], w[92], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[78]);
  FFS1_2_0(w[93], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[79]);
  FFS1_2_0(w[4], w[95], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[80]);
  FFS1_2_0(w[96], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[81]);
  FFS1_2_0(w[4], w[97], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[82]);
  FFS1_2_0(w[98], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[83]);
  FFS1_2_0(w[2], w[100], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[84]);
  FFS1_2_0(w[101], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[85]);
  FFS1_2_0(w[2], w[102], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[86]);
  FFS1_2_0(w[103], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[87]);
  FFS1_2_0(w[105], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[88]);
  FFS1_2_0(w[2], w[106], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[89]);
  FFS1_2_0(w[107], w[3], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[90]);
  FFS1_2_0(w[2], w[108], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[91]);
  FFS1_2_0(w[2], w[110], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[92]);
  FFS1_2_0(w[111], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[93]);
  FFS1_2_0(w[2], w[112], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[94]);
  FFS1_2_0(w[113], w[1], w[0], pars->HEF_MEKD2_1_GC_137, pars->HEF_MEKD2_1_GC_136, amp[95]);

}
double Spin0_2f_SFpA::matrix_h_mummupmummupa_no_zpxg() 
{
  int i, j;
  // Local variables
	// Commented out by Convert_source 0.2
//  const int ngraphs = 96;
  const int ncolor = 1;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[1] = {1.};
  static const double cf[1][1] = {{1.}};

  // Calculate color flows
  jamp[0] = +amp[0] + amp[1] + amp[2] + amp[3] + amp[4] + amp[5] + amp[6] +
      amp[7] + amp[8] + amp[9] + amp[10] + amp[11] + amp[12] + amp[13] +
      amp[14] + amp[15] - amp[16] - amp[17] - amp[18] - amp[19] - amp[20] -
      amp[21] - amp[22] - amp[23] - amp[24] - amp[25] - amp[26] - amp[27] -
      amp[28] - amp[29] - amp[30] - amp[31] + amp[32] + amp[33] + amp[34] +
      amp[35] - amp[36] - amp[37] - amp[38] - amp[39] - amp[40] - amp[41] -
      amp[42] - amp[43] - amp[44] - amp[45] + amp[46] + amp[47] + amp[48] +
      amp[49] + amp[50] + amp[51] - amp[52] - amp[53] - amp[54] - amp[55] -
      amp[56] - amp[57] - amp[58] - amp[59] - amp[60] - amp[61] + amp[62] +
      amp[63] + amp[64] + amp[65] - amp[66] - amp[67] - amp[68] - amp[69] +
      amp[70] + amp[71] + amp[72] + amp[73] + amp[74] + amp[75] + amp[76] +
      amp[77] + amp[78] + amp[79] - amp[80] - amp[81] - amp[82] - amp[83] +
      amp[84] + amp[85] + amp[86] + amp[87] - amp[88] - amp[89] - amp[90] -
      amp[91] + amp[92] + amp[93] + amp[94] + amp[95];

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



