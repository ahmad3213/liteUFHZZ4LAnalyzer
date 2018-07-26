//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph 5 v. 1.5.2, 2012-10-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef Parameters_MEKD_CC
#define Parameters_MEKD_CC

#include <iostream> 
#include <iomanip> 
#include "Parameters_MEKD.h"

// Initialize static instance
Parameters_MEKD * Parameters_MEKD::instance = 0; 

// Function to get static instance - only one instance per program
Parameters_MEKD * Parameters_MEKD::getInstance()
{
  if (instance == 0)
    instance = new Parameters_MEKD(); 

  return instance; 
}

void Parameters_MEKD::setIndependentParameters(SLHAReader_MEKD& slha)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  WXG = slha.get_block_entry("decay", 9000007, 5.753088e-03).real(); 
  WH = slha.get_block_entry("decay", 9000006, 5.753088e-03).real();
  WZp = slha.get_block_entry("decay", 300, 5.753088e-03).real();
  WW = slha.get_block_entry("decay", 24, 2.085000e+00).real(); 
  WZ = slha.get_block_entry("decay", 23, 2.495200e+00).real(); 
  WT = slha.get_block_entry("decay", 6, 1.508336e+00).real(); 
  ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00).real(); 
  ymm = slha.get_block_entry("yukawa", 13, 1.056600e-01).real(); 
  yme = slha.get_block_entry("yukawa", 11, 5.110000e-04).real(); 
  ymt = slha.get_block_entry("yukawa", 6, 1.720000e+02).real(); 
  ymb = slha.get_block_entry("yukawa", 5, 4.700000e+00).real(); 
  ymc = slha.get_block_entry("yukawa", 4, 1.270000e+00).real(); 
  yms = slha.get_block_entry("yukawa", 3, 1.010000e-01).real(); 
  ymup = slha.get_block_entry("yukawa", 2, 2.550000e-03).real(); 
  ymdo = slha.get_block_entry("yukawa", 1, 5.040000e-03).real(); 
  aS = slha.get_block_entry("sminputs", 3, 1.184000e-01).real(); 
  Gf = slha.get_block_entry("sminputs", 2, 1.166370e-05).real(); 
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.279000e+02).real(); 
  MXG = slha.get_block_entry("mass", 9000007, 1.250000e+02).real(); 
  MH = slha.get_block_entry("mass", 9000006, 1.250000e+02).real();
  MZp = slha.get_block_entry("mass", 300, 1.250000e+02).real();
  MZ = slha.get_block_entry("mass", 23, 9.118760e+01).real(); 
  MTA = slha.get_block_entry("mass", 15, 1.777000e+00).real(); 
  MM = slha.get_block_entry("mass", 13, 1.056600e-01).real(); 
  Me = slha.get_block_entry("mass", 11, 5.110000e-04).real(); 
  MT = slha.get_block_entry("mass", 6, 1.720000e+02).real(); 
  MB = slha.get_block_entry("mass", 5, 4.700000e+00).real(); 
  MC = slha.get_block_entry("mass", 4, 1.270000e+00).real(); 
  MS = slha.get_block_entry("mass", 3, 1.010000e-01).real(); 
  MU = slha.get_block_entry("mass", 2, 2.550000e-03).real(); 
  MD = slha.get_block_entry("mass", 1, 5.040000e-03).real();
  rhoe02 = slha.get_block_entry("heff", 22, 0.000000e+00); 
  rhoe01 = slha.get_block_entry("heff", 21, 2.075371e-06);
  rhomu02 = slha.get_block_entry("heff", 20, 0.000000e+00); 
  rhomu01 = slha.get_block_entry("heff", 19, 4.291210e-04); 
  rhob02 = slha.get_block_entry("heff", 18, 1.000000e-01); 
  rhob01 = slha.get_block_entry("heff", 17, 1.000000e-01); 
  rhos02 = slha.get_block_entry("heff", 16, 1.000000e-01); 
  rhos01 = slha.get_block_entry("heff", 15, 1.000000e-01); 
  rhod02 = slha.get_block_entry("heff", 14, 1.000000e-01); 
  rhod01 = slha.get_block_entry("heff", 13, 1.000000e-01); 
  rhoc02 = slha.get_block_entry("heff", 12, 1.000000e-01); 
  rhoc01 = slha.get_block_entry("heff", 11, 1.000000e-01); 
  rhou02 = slha.get_block_entry("heff", 10, 1.000000e-01); 
  rhou01 = slha.get_block_entry("heff", 9, 1.000000e-01); 
  g4z = slha.get_block_entry("heff", 8, 1.000000e-01); 
  g3z = slha.get_block_entry("heff", 7, 1.000000e-01); 
  g2z = slha.get_block_entry("heff", 6, 1.000000e-01); 
  g1z = slha.get_block_entry("heff", 5, 1.000000e-01); 
  g4g = slha.get_block_entry("heff", 4, 1.000000e-01); 
  g3g = slha.get_block_entry("heff", 3, 1.000000e-01); 
  g2g = slha.get_block_entry("heff", 2, 1.000000e-01); 
  g1g = slha.get_block_entry("heff", 1, 1.000000e-01);
  rhoe14 = slha.get_block_entry("vec", 30, 0.000000e+00);
  rhoe13 = slha.get_block_entry("vec", 29, 0.000000e+00);
  rhoe12 = slha.get_block_entry("vec", 28, 0.000000e+00);
  rhoe11 = slha.get_block_entry("vec", 27, 2.075371e-06);
  rhomu14 = slha.get_block_entry("vec", 26, 0.000000e+00); 
  rhomu13 = slha.get_block_entry("vec", 25, 0.000000e+00); 
  rhomu12 = slha.get_block_entry("vec", 24, 0.000000e+00); 
  rhomu11 = slha.get_block_entry("vec", 23, 4.291210e-04);
  rhob14 = slha.get_block_entry("vec", 22, 1.000000e-01); 
  rhob13 = slha.get_block_entry("vec", 21, 1.000000e-01); 
  rhob12 = slha.get_block_entry("vec", 20, 1.000000e-01); 
  rhob11 = slha.get_block_entry("vec", 19, 1.000000e-01); 
  rhos14 = slha.get_block_entry("vec", 18, 1.000000e-01); 
  rhos13 = slha.get_block_entry("vec", 17, 1.000000e-01); 
  rhos12 = slha.get_block_entry("vec", 16, 1.000000e-01); 
  rhos11 = slha.get_block_entry("vec", 15, 1.000000e-01); 
  rhod14 = slha.get_block_entry("vec", 14, 1.000000e-01); 
  rhod13 = slha.get_block_entry("vec", 13, 1.000000e-01); 
  rhod12 = slha.get_block_entry("vec", 12, 1.000000e-01); 
  rhod11 = slha.get_block_entry("vec", 11, 1.000000e-01); 
  rhoc14 = slha.get_block_entry("vec", 10, 1.000000e-01); 
  rhoc13 = slha.get_block_entry("vec", 9, 1.000000e-01); 
  rhoc12 = slha.get_block_entry("vec", 8, 1.000000e-01); 
  rhoc11 = slha.get_block_entry("vec", 7, 1.000000e-01); 
  rhou14 = slha.get_block_entry("vec", 6, 1.000000e-01); 
  rhou13 = slha.get_block_entry("vec", 5, 1.000000e-01); 
  rhou12 = slha.get_block_entry("vec", 4, 1.000000e-01); 
  rhou11 = slha.get_block_entry("vec", 3, 1.000000e-01); 
  b2z = slha.get_block_entry("vec", 2, 1.000000e-01); 
  b1z = slha.get_block_entry("vec", 1, 1.000000e-01);
  rhoe24 = slha.get_block_entry("gravity", 48, 0.000000e+00);
  rhoe23 = slha.get_block_entry("gravity", 47, 0.000000e+00);
  rhoe22 = slha.get_block_entry("gravity", 46, 0.000000e+00);
  rhoe21 = slha.get_block_entry("gravity", 45, 2.075371e-06);
  rhomu24 = slha.get_block_entry("gravity", 44, 0.000000e+00); 
  rhomu23 = slha.get_block_entry("gravity", 43, 0.000000e+00); 
  rhomu22 = slha.get_block_entry("gravity", 42, 0.000000e+00); 
  rhomu21 = slha.get_block_entry("gravity", 41, 4.291210e-04); 
  rhob24 = slha.get_block_entry("gravity", 40, 1.000000e-01); 
  rhob23 = slha.get_block_entry("gravity", 39, 1.000000e-01); 
  rhob22 = slha.get_block_entry("gravity", 38, 1.000000e-01); 
  rhob21 = slha.get_block_entry("gravity", 37, 1.000000e-01); 
  rhos24 = slha.get_block_entry("gravity", 36, 1.000000e-01); 
  rhos23 = slha.get_block_entry("gravity", 35, 1.000000e-01); 
  rhos22 = slha.get_block_entry("gravity", 34, 1.000000e-01); 
  rhos21 = slha.get_block_entry("gravity", 33, 1.000000e-01); 
  rhod24 = slha.get_block_entry("gravity", 32, 1.000000e-01); 
  rhod23 = slha.get_block_entry("gravity", 31, 1.000000e-01); 
  rhod22 = slha.get_block_entry("gravity", 30, 1.000000e-01); 
  rhod21 = slha.get_block_entry("gravity", 29, 1.000000e-01); 
  rhoc24 = slha.get_block_entry("gravity", 28, 1.000000e-01); 
  rhoc23 = slha.get_block_entry("gravity", 27, 1.000000e-01); 
  rhoc22 = slha.get_block_entry("gravity", 26, 1.000000e-01); 
  rhoc21 = slha.get_block_entry("gravity", 25, 1.000000e-01); 
  rhou24 = slha.get_block_entry("gravity", 24, 1.000000e-01); 
  rhou23 = slha.get_block_entry("gravity", 23, 1.000000e-01); 
  rhou22 = slha.get_block_entry("gravity", 22, 1.000000e-01); 
  rhou21 = slha.get_block_entry("gravity", 21, 1.000000e-01); 
  k10z = slha.get_block_entry("gravity", 20, 1.000000e-01); 
  k9z = slha.get_block_entry("gravity", 19, 1.000000e-01); 
  k8z = slha.get_block_entry("gravity", 18, 1.000000e-01); 
  k7z = slha.get_block_entry("gravity", 17, 1.000000e-01); 
  k6z = slha.get_block_entry("gravity", 16, 1.000000e-01); 
  k5z = slha.get_block_entry("gravity", 15, 1.000000e-01); 
  k4z = slha.get_block_entry("gravity", 14, 1.000000e-01); 
  k3z = slha.get_block_entry("gravity", 13, 1.000000e-01); 
  k2z = slha.get_block_entry("gravity", 12, 1.000000e-01); 
  k1z = slha.get_block_entry("gravity", 11, 1.000000e-01);
  k10g = slha.get_block_entry("gravity", 10, 1.000000e-01); 
  k9g = slha.get_block_entry("gravity", 9, 1.000000e-01); 
  k8g = slha.get_block_entry("gravity", 8, 1.000000e-01); 
  k7g = slha.get_block_entry("gravity", 7, 1.000000e-01); 
  k6g = slha.get_block_entry("gravity", 6, 1.000000e-01); 
  k5g = slha.get_block_entry("gravity", 5, 1.000000e-01); 
  k4g = slha.get_block_entry("gravity", 4, 1.000000e-01); 
  k3g = slha.get_block_entry("gravity", 3, 1.000000e-01); 
  k2g = slha.get_block_entry("gravity", 2, 1.000000e-01); 
  k1g = slha.get_block_entry("gravity", 1, 1.000000e-01);
  cabi = slha.get_block_entry("ckmblock", 1, 2.277360e-01).real();
  gw = 1.; 
  g1 = 1.; 
  cos__cabi = cos(cabi); 
  CKM11 = cos__cabi; 
  sin__cabi = sin(cabi); 
  CKM12 = sin__cabi; 
  CKM13 = 0.; 
  CKM21 = -sin__cabi; 
  CKM22 = cos__cabi; 
  CKM23 = 0.; 
  CKM31 = 0.; 
  CKM32 = 0.; 
  CKM33 = 1.; 
  MZ__exp__2 = pow(MZ, 2.); 
  MZ__exp__4 = pow(MZ, 4.); 
  sqrt__2 = sqrt(2.); 
  MH__exp__2 = pow(MH, 2.); 
  conjg__CKM11 = conj(CKM11); 
  conjg__CKM21 = conj(CKM21); 
  conjg__CKM31 = conj(CKM31); 
  conjg__CKM12 = conj(CKM12); 
  conjg__CKM22 = conj(CKM22); 
  conjg__CKM32 = conj(CKM32); 
  conjg__CKM13 = conj(CKM13); 
  conjg__CKM23 = conj(CKM23); 
  conjg__CKM33 = conj(CKM33); 
  complexi = std::complex<double> (0., 1.); 
  aEW = 1./aEWM1; 
  MW = sqrt(MZ__exp__2/2. + sqrt(MZ__exp__4/4. - (aEW * M_PI * MZ__exp__2)/(Gf
      * sqrt__2)));
  sqrt__aEW = sqrt(aEW); 
  ee = 2. * sqrt__aEW * sqrt(M_PI); 
  MW__exp__2 = pow(MW, 2.); 
  sw2 = 1. - MW__exp__2/MZ__exp__2; 
  cw = sqrt(1. - sw2); 
  sqrt__sw2 = sqrt(sw2); 
  sw = sqrt__sw2; 
  vev = (2. * MW * sw)/ee; 
  vev__exp__2 = pow(vev, 2.); 
  lam = MH__exp__2/(2. * vev__exp__2); 
  yb = (ymb * sqrt__2)/vev; 
  yc = (ymc * sqrt__2)/vev; 
  ydo = (ymdo * sqrt__2)/vev; 
  ye = (yme * sqrt__2)/vev; 
  ym = (ymm * sqrt__2)/vev; 
  ys = (yms * sqrt__2)/vev; 
  yt = (ymt * sqrt__2)/vev; 
  ytau = (ymtau * sqrt__2)/vev; 
  yup = (ymup * sqrt__2)/vev; 
  muH = sqrt(lam * vev__exp__2); 
  I1x11 = ydo * conjg__CKM11; 
  I1x12 = ydo * conjg__CKM21; 
  I1x13 = ydo * conjg__CKM31; 
  I1x21 = ys * conjg__CKM12; 
  I1x22 = ys * conjg__CKM22; 
  I1x23 = ys * conjg__CKM32; 
  I1x31 = yb * conjg__CKM13; 
  I1x32 = yb * conjg__CKM23; 
  I1x33 = yb * conjg__CKM33; 
  I2x11 = yup * conjg__CKM11; 
  I2x12 = yc * conjg__CKM21; 
  I2x13 = yt * conjg__CKM31; 
  I2x21 = yup * conjg__CKM12; 
  I2x22 = yc * conjg__CKM22; 
  I2x23 = yt * conjg__CKM32; 
  I2x31 = yup * conjg__CKM13; 
  I2x32 = yc * conjg__CKM23; 
  I2x33 = yt * conjg__CKM33; 
  I3x11 = CKM11 * yup; 
  I3x12 = CKM21 * yc; 
  I3x13 = CKM31 * yt; 
  I3x21 = CKM12 * yup; 
  I3x22 = CKM22 * yc; 
  I3x23 = CKM32 * yt; 
  I3x31 = CKM13 * yup; 
  I3x32 = CKM23 * yc; 
  I3x33 = CKM33 * yt; 
  I4x11 = CKM11 * ydo; 
  I4x12 = CKM21 * ydo; 
  I4x13 = CKM31 * ydo; 
  I4x21 = CKM12 * ys; 
  I4x22 = CKM22 * ys; 
  I4x23 = CKM32 * ys; 
  I4x31 = CKM13 * yb; 
  I4x32 = CKM23 * yb; 
  I4x33 = CKM33 * yb; 
  ee__exp__2 = pow(ee, 2.); 
  sw__exp__2 = pow(sw, 2.); 
  cw__exp__2 = pow(cw, 2.); 
}

void Parameters_MEKD::setIndependentCouplings()
{
	GC_1 = -(ee * complexi)/3.;
	GC_2 = (2. * ee * complexi)/3.;
	GC_3 = -(ee * complexi);
	GC_109 = -(cw * ee * complexi)/(2. * sw);
	GC_110 = (cw * ee * complexi)/(2. * sw);
	GC_115 = -(ee * complexi * sw)/(6. * cw);  
	GC_116 = (ee * complexi * sw)/(2. * cw);
	
	Unitary_GC_5 = -(ee * complexi)/3.; 
	Unitary_GC_6 = (2. * ee * complexi)/3.; 
	Unitary_GC_7 = -(ee * complexi); 
	Unitary_GC_70 = -(cw * ee * complexi)/(2. * sw); 
	Unitary_GC_71 = (cw * ee * complexi)/(2. * sw); 
	Unitary_GC_74 = -(ee * complexi * sw)/(6. * cw); 
	Unitary_GC_75 = (ee * complexi * sw)/(2. * cw);
	
	// Coming from spin 0
	HEF_MEKD_GC_3 = -(ee * complexi)/3.; 
	HEF_MEKD_GC_4 = (2. * ee * complexi)/3.; 
	HEF_MEKD_GC_5 = -(ee * complexi);
	HEF_MEKD_GC_13 = -(complexi * g1g); 
	HEF_MEKD_GC_14 = -(complexi * g1z); 
	HEF_MEKD_GC_15 = -(complexi * g2g); 
	HEF_MEKD_GC_18 = -(complexi * g2z); 
	HEF_MEKD_GC_19 = -(complexi * g3g); 
	HEF_MEKD_GC_22 = -(complexi * g3z); 
	HEF_MEKD_GC_23 = (complexi * g4g)/8.; 
	HEF_MEKD_GC_25 = (complexi * g4z)/2.; 
	HEF_MEKD_GC_106 = complexi * rhoc01; 
	HEF_MEKD_GC_107 = -rhoc02; 
	HEF_MEKD_GC_116 = complexi * rhod01; 
	HEF_MEKD_GC_117 = -rhod02; 
	HEF_MEKD_GC_126 = complexi * rhos01; 
	HEF_MEKD_GC_127 = -rhos02; 
	HEF_MEKD_GC_136 = complexi * rhou01; 
	HEF_MEKD_GC_137 = -rhou02;
	HEF_MEKD_GC_161 = -(cw * ee * complexi)/(2. * sw); 
	HEF_MEKD_GC_168 = (ee * complexi * sw)/(2. * cw); 
	
	// Coming extra from spin 1
	HEF_MEKD_GC_1 = -b1z; 
	HEF_MEKD_GC_2 = -2. * b2z; 
// 	HEF_MEKD_GC_3 = -(ee * complexi)/3.; 
// 	HEF_MEKD_GC_4 = (2. * ee * complexi)/3.; 
// 	HEF_MEKD_GC_5 = -(ee * complexi); 
	HEF_MEKD_GC_108 = complexi * rhoc11; 
	HEF_MEKD_GC_109 = complexi * rhoc12; 
	HEF_MEKD_GC_110 = complexi * rhoc13; 
	HEF_MEKD_GC_111 = rhoc14; 
	HEF_MEKD_GC_118 = complexi * rhod11; 
	HEF_MEKD_GC_119 = complexi * rhod12; 
	HEF_MEKD_GC_120 = complexi * rhod13; 
	HEF_MEKD_GC_121 = rhod14; 
	HEF_MEKD_GC_128 = complexi * rhos11; 
	HEF_MEKD_GC_129 = complexi * rhos12; 
	HEF_MEKD_GC_130 = complexi * rhos13; 
	HEF_MEKD_GC_131 = rhos14; 
	HEF_MEKD_GC_138 = complexi * rhou11; 
	HEF_MEKD_GC_139 = complexi * rhou12; 
	HEF_MEKD_GC_140 = complexi * rhou13; 
	HEF_MEKD_GC_141 = rhou14;
// 	HEF_MEKD_GC_161 = -(cw * ee * complexi)/(2. * sw); 
// 	HEF_MEKD_GC_168 = (ee * complexi * sw)/(2. * cw); 
	
	// Coming extra from spin 2
// 	HEF_MEKD_GC_3 = -(ee * complexi)/3.; 
// 	HEF_MEKD_GC_4 = (2. * ee * complexi)/3.; 
// 	HEF_MEKD_GC_5 = -(ee * complexi); 
	HEF_MEKD_GC_62 = -(complexi * k10g)/2.; 
	HEF_MEKD_GC_63 = -(complexi * k10z)/2.; 
	HEF_MEKD_GC_64 = -(complexi * k1g); 
	HEF_MEKD_GC_67 = -(complexi * k1z); 
	HEF_MEKD_GC_68 = complexi * k2g; 
	HEF_MEKD_GC_71 = complexi * k2z; 
	HEF_MEKD_GC_72 = complexi * k3g; 
	HEF_MEKD_GC_75 = complexi * k3z; 
	HEF_MEKD_GC_76 = -2. * complexi * k4g; 
	HEF_MEKD_GC_79 = -2. * complexi * k4z; 
	HEF_MEKD_GC_80 = complexi * k5g; 
	HEF_MEKD_GC_81 = complexi * k5z; 
	HEF_MEKD_GC_82 = -(complexi * k6g)/2.; 
	HEF_MEKD_GC_83 = -(complexi * k6z)/2.; 
	HEF_MEKD_GC_84 = -(complexi * k7g); 
	HEF_MEKD_GC_85 = -(complexi * k7z); 
	HEF_MEKD_GC_86 = (complexi * k8g)/4.; 
	HEF_MEKD_GC_90 = complexi * k8z; 
	HEF_MEKD_GC_91 = -(complexi * k9g)/2.; 
	HEF_MEKD_GC_92 = -(complexi * k9z)/2.; 
	HEF_MEKD_GC_112 = (complexi * rhoc21)/2.; 
	HEF_MEKD_GC_113 = (complexi * rhoc22)/2.; 
	HEF_MEKD_GC_114 = -(complexi * rhoc23)/2.; 
	HEF_MEKD_GC_115 = rhoc24/2.; 
	HEF_MEKD_GC_122 = (complexi * rhod21)/2.; 
	HEF_MEKD_GC_123 = (complexi * rhod22)/2.; 
	HEF_MEKD_GC_124 = -(complexi * rhod23)/2.; 
	HEF_MEKD_GC_125 = rhod24/2.; 
	HEF_MEKD_GC_132 = (complexi * rhos21)/2.; 
	HEF_MEKD_GC_133 = (complexi * rhos22)/2.; 
	HEF_MEKD_GC_134 = -(complexi * rhos23)/2.; 
	HEF_MEKD_GC_135 = rhos24/2.; 
	HEF_MEKD_GC_142 = (complexi * rhou21)/2.; 
	HEF_MEKD_GC_143 = (complexi * rhou22)/2.; 
	HEF_MEKD_GC_144 = -(complexi * rhou23)/2.; 
	HEF_MEKD_GC_145 = rhou24/2.;
// 	HEF_MEKD_GC_161 = -(cw * ee * complexi)/(2. * sw); 
// 	HEF_MEKD_GC_168 = (ee * complexi * sw)/(2. * cw);
	
	
	/// Model HEF_MEKD2_1 with 2l couplings
	// Coming extra from DY
	HEF_MEKD2_1_GC_3 = -(ee * complexi)/3.;
	HEF_MEKD2_1_GC_4 = (2. * ee * complexi)/3.;
	HEF_MEKD2_1_GC_5 = -(ee * complexi);
	HEF_MEKD2_1_GC_181 = -(cw * ee * complexi)/(2. * sw);
	HEF_MEKD2_1_GC_182 = (cw * ee * complexi)/(2. * sw);
	HEF_MEKD2_1_GC_187 = -(ee * complexi * sw)/(6. * cw);
	HEF_MEKD2_1_GC_188 = (ee * complexi * sw)/(2. * cw);
	
	// Coming extra from spin 0
// 	HEF_MEKD2_1_GC_5 = -(ee * complexi);
	HEF_MEKD2_1_GC_13 = -(complexi * g1g);
	HEF_MEKD2_1_GC_14 = -(complexi * g1z);
	HEF_MEKD2_1_GC_15 = -(complexi * g2g);
	HEF_MEKD2_1_GC_18 = -(complexi * g2z);
	HEF_MEKD2_1_GC_19 = -(complexi * g3g);
	HEF_MEKD2_1_GC_22 = -(complexi * g3z);
	HEF_MEKD2_1_GC_23 = (complexi * g4g)/8.;
	HEF_MEKD2_1_GC_25 = (complexi * g4z)/2.;
	HEF_MEKD2_1_GC_126 = complexi * rhoe01;
	HEF_MEKD2_1_GC_127 = -rhoe02;
	HEF_MEKD2_1_GC_136 = complexi * rhomu01;
	HEF_MEKD2_1_GC_137 = -rhomu02;
// 	HEF_MEKD2_1_GC_181 = -(cw * ee * complexi)/(2. * sw);
// 	HEF_MEKD2_1_GC_188 = (ee * complexi * sw)/(2. * cw);
	
	// Coming extra from spin 1
	HEF_MEKD2_1_GC_1 = -b1z;
	HEF_MEKD2_1_GC_2 = -2. * b2z;
// 	HEF_MEKD2_1_GC_3 = -(ee * complexi)/3.;
// 	HEF_MEKD2_1_GC_4 = (2. * ee * complexi)/3.;
// 	HEF_MEKD2_1_GC_5 = -(ee * complexi);
	HEF_MEKD2_1_GC_108 = complexi * rhoc11;
	HEF_MEKD2_1_GC_109 = complexi * rhoc12;
	HEF_MEKD2_1_GC_110 = complexi * rhoc13;
	HEF_MEKD2_1_GC_111 = rhoc14;
	HEF_MEKD2_1_GC_128 = complexi * rhoe11;
	HEF_MEKD2_1_GC_129 = complexi * rhoe12;
	HEF_MEKD2_1_GC_130 = complexi * rhoe13;
	HEF_MEKD2_1_GC_131 = rhoe14;
	HEF_MEKD2_1_GC_138 = complexi * rhomu11;
	HEF_MEKD2_1_GC_139 = complexi * rhomu12;
	HEF_MEKD2_1_GC_140 = complexi * rhomu13;
	HEF_MEKD2_1_GC_141 = rhomu14;
	HEF_MEKD2_1_GC_148 = complexi * rhos11;
	HEF_MEKD2_1_GC_149 = complexi * rhos12;
	HEF_MEKD2_1_GC_150 = complexi * rhos13;
	HEF_MEKD2_1_GC_151 = rhos14;
	
	// Coming extra from spin 2
// 	HEF_MEKD2_1_GC_3 = -(ee * complexi)/3.;
// 	HEF_MEKD2_1_GC_4 = (2. * ee * complexi)/3.;
// 	HEF_MEKD2_1_GC_5 = -(ee * complexi);
	HEF_MEKD2_1_GC_62 = -(complexi * k10g)/2.;
	HEF_MEKD2_1_GC_63 = -(complexi * k10z)/2.;
	HEF_MEKD2_1_GC_64 = -(complexi * k1g);
	HEF_MEKD2_1_GC_67 = -(complexi * k1z);
	HEF_MEKD2_1_GC_68 = complexi * k2g;
	HEF_MEKD2_1_GC_71 = complexi * k2z;
	HEF_MEKD2_1_GC_72 = complexi * k3g;
	HEF_MEKD2_1_GC_75 = complexi * k3z;
	HEF_MEKD2_1_GC_76 = -2. * complexi * k4g;
	HEF_MEKD2_1_GC_79 = -2. * complexi * k4z;
	HEF_MEKD2_1_GC_80 = complexi * k5g;
	HEF_MEKD2_1_GC_81 = complexi * k5z;
	HEF_MEKD2_1_GC_82 = -(complexi * k6g)/2.;
	HEF_MEKD2_1_GC_83 = -(complexi * k6z)/2.;
	HEF_MEKD2_1_GC_84 = -(complexi * k7g);
	HEF_MEKD2_1_GC_85 = -(complexi * k7z);
	HEF_MEKD2_1_GC_86 = (complexi * k8g)/4.;
	HEF_MEKD2_1_GC_90 = complexi * k8z;
	HEF_MEKD2_1_GC_91 = -(complexi * k9g)/2.;
	HEF_MEKD2_1_GC_92 = -(complexi * k9z)/2.;
	HEF_MEKD2_1_GC_112 = (complexi * rhoc21)/2.;
	HEF_MEKD2_1_GC_113 = (complexi * rhoc22)/2.;
	HEF_MEKD2_1_GC_114 = -(complexi * rhoc23)/2.;
	HEF_MEKD2_1_GC_115 = rhoc24/2.;
	HEF_MEKD2_1_GC_132 = (complexi * rhoe21)/2.;
	HEF_MEKD2_1_GC_133 = (complexi * rhoe22)/2.;
	HEF_MEKD2_1_GC_134 = -(complexi * rhoe23)/2.;
	HEF_MEKD2_1_GC_135 = rhoe24/2.;
	HEF_MEKD2_1_GC_142 = (complexi * rhomu21)/2.;
	HEF_MEKD2_1_GC_143 = (complexi * rhomu22)/2.;
	HEF_MEKD2_1_GC_144 = -(complexi * rhomu23)/2.;
	HEF_MEKD2_1_GC_145 = rhomu24/2.;
	HEF_MEKD2_1_GC_152 = (complexi * rhos21)/2.;
	HEF_MEKD2_1_GC_153 = (complexi * rhos22)/2.;
	HEF_MEKD2_1_GC_154 = -(complexi * rhos23)/2.;
	HEF_MEKD2_1_GC_155 = rhos24/2.;
// 	HEF_MEKD2_1_GC_181 = -(cw * ee * complexi)/(2. * sw);
// 	HEF_MEKD2_1_GC_182 = (cw * ee * complexi)/(2. * sw);
// 	HEF_MEKD2_1_GC_187 = -(ee * complexi * sw)/(6. * cw);
// 	HEF_MEKD2_1_GC_188 = (ee * complexi * sw)/(2. * cw);
}

void Parameters_MEKD::setDependentParameters()
{
  sqrt__aS = sqrt(aS); 
  G = 2. * sqrt__aS * sqrt(M_PI); 
  G__exp__2 = pow(G, 2.); 
}

void Parameters_MEKD::setDependentCouplings()
{

}



#endif
