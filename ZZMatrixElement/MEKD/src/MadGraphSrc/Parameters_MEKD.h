//==========================================================================
// This file has been automatically generated for C++
// MadGraph 5 v. 1.5.2, 2012-10-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef Parameters_MEKD_H
#define Parameters_MEKD_H

#include <complex> 

#include "read_slha.h"
using namespace std;

 
class Parameters_MEKD
{
public:
	static Parameters_MEKD * getInstance();
	
	// Define "zero"
	double zero, ZERO;
	// Model parameters independent of aS
	double WXG, WH, WZp, WW, WZ, WT, ymtau, ymm, yme, ymt, ymb, ymc, yms, ymup,
		ymdo, aS, Gf, aEWM1, MXG, MH, MZp, MZ, MTA, MM, Me, MT, MB, MC, MS, MU, MD,
		cabi, gw, g1, cos__cabi, sin__cabi, MZ__exp__2,
		MZ__exp__4, sqrt__2, MH__exp__2, aEW, MW, sqrt__aEW, ee, MW__exp__2,
		sw2, cw, sqrt__sw2, sw, vev, vev__exp__2, lam, yb, yc, ydo, ye, ym,
		ys, yt, ytau, yup, muH, ee__exp__2, sw__exp__2, cw__exp__2;
	// Spin 0
	std::complex<double> rhob02, rhob01, rhos02, rhos01, rhod02, rhod01, rhoc02, rhoc01, rhou02, rhou01, g4z, g3z, g2z, g1z, g4g, g3g, g2g, g1g;
	std::complex<double> rhomu02, rhomu01, rhoe02, rhoe01;	// 2l
	
	// Spin 1
	std::complex<double> rhob14, rhob13, rhob12, rhob11, rhos14, rhos13, rhos12, rhos11, rhod14, rhod13, rhod12, rhod11, rhoc14, rhoc13, rhoc12, rhoc11, rhou14, rhou13, rhou12, rhou11, b2z, b1z;
	std::complex<double> rhomu14, rhomu13, rhomu12, rhomu11, rhoe14, rhoe13, rhoe12, rhoe11;	// 2l
	
	// Spin 2
	std::complex<double> rhob24, rhob23, rhob22, rhob21, rhos24, rhos23, rhos22, rhos21, rhod24, rhod23, rhod22, rhod21, rhoc24, rhoc23, rhoc22, rhoc21, rhou24, rhou23, rhou22, rhou21, k10g, k9g, k8g, k7g, k6g, k5g, k4g, k3g, k2g, k1g, k10z, k9z, k8z, k7z, k6z, k5z, k4z, k3z, k2z, k1z;
	std::complex<double> rhomu24, rhomu23, rhomu22, rhomu21, rhoe24, rhoe23, rhoe22, rhoe21;	// 2l
	
	std::complex<double> CKM11, CKM12, CKM13, CKM21, CKM22, CKM23, CKM31,
        CKM32, CKM33, conjg__CKM11, conjg__CKM21, conjg__CKM31, conjg__CKM12,
        conjg__CKM22, conjg__CKM32, conjg__CKM13, conjg__CKM23, conjg__CKM33,
        complexi, I1x11, I1x12, I1x13, I1x21, I1x22, I1x23, I1x31, I1x32,
        I1x33, I2x11, I2x12, I2x13, I2x21, I2x22, I2x23, I2x31, I2x32, I2x33,
        I3x11, I3x12, I3x13, I3x21, I3x22, I3x23, I3x31, I3x32, I3x33, I4x11,
        I4x12, I4x13, I4x21, I4x22, I4x23, I4x31, I4x32, I4x33;
    // Model parameters dependent on aS
    double sqrt__aS, G, G__exp__2; 
    // Model couplings independent of aS, ZZ part
    std::complex<double> GC_1, GC_2, GC_3, GC_109, GC_110, GC_115, GC_116,
		Unitary_GC_5, Unitary_GC_6, Unitary_GC_7, Unitary_GC_70, Unitary_GC_71, Unitary_GC_74, Unitary_GC_75;
		
	// Coming extra from spin 0
	std::complex<double> HEF_MEKD_GC_3, HEF_MEKD_GC_4, HEF_MEKD_GC_5, HEF_MEKD_GC_13, HEF_MEKD_GC_14, HEF_MEKD_GC_15, HEF_MEKD_GC_18, HEF_MEKD_GC_19, HEF_MEKD_GC_22, HEF_MEKD_GC_23, HEF_MEKD_GC_25, HEF_MEKD_GC_106, HEF_MEKD_GC_107, HEF_MEKD_GC_116, HEF_MEKD_GC_117, HEF_MEKD_GC_126, HEF_MEKD_GC_127, HEF_MEKD_GC_136, HEF_MEKD_GC_137, HEF_MEKD_GC_161, HEF_MEKD_GC_168;
		
	// Coming extra from spin 1
	std::complex<double> HEF_MEKD_GC_1, HEF_MEKD_GC_2, HEF_MEKD_GC_108, HEF_MEKD_GC_109, HEF_MEKD_GC_110, HEF_MEKD_GC_111, HEF_MEKD_GC_118, HEF_MEKD_GC_119, HEF_MEKD_GC_120, HEF_MEKD_GC_121, HEF_MEKD_GC_128, HEF_MEKD_GC_129, HEF_MEKD_GC_130, HEF_MEKD_GC_131, HEF_MEKD_GC_138, HEF_MEKD_GC_139, HEF_MEKD_GC_140, HEF_MEKD_GC_141;
		
	// Coming extra from spin 2
	std::complex<double> HEF_MEKD_GC_62, HEF_MEKD_GC_63, HEF_MEKD_GC_64, HEF_MEKD_GC_67, HEF_MEKD_GC_68, HEF_MEKD_GC_71, HEF_MEKD_GC_72, HEF_MEKD_GC_75 , HEF_MEKD_GC_76, HEF_MEKD_GC_79, HEF_MEKD_GC_80, HEF_MEKD_GC_81, HEF_MEKD_GC_82, HEF_MEKD_GC_83, HEF_MEKD_GC_84, HEF_MEKD_GC_85, HEF_MEKD_GC_86, HEF_MEKD_GC_90, HEF_MEKD_GC_91, HEF_MEKD_GC_92, HEF_MEKD_GC_112, HEF_MEKD_GC_113, HEF_MEKD_GC_114, HEF_MEKD_GC_115, HEF_MEKD_GC_122, HEF_MEKD_GC_123, HEF_MEKD_GC_124, HEF_MEKD_GC_125, HEF_MEKD_GC_132, HEF_MEKD_GC_133, HEF_MEKD_GC_134, HEF_MEKD_GC_135, HEF_MEKD_GC_142, HEF_MEKD_GC_143, HEF_MEKD_GC_144, HEF_MEKD_GC_145;
	
	
	/// HEF_MEKD2_1 model
	// Coming from DY
	std::complex<double> HEF_MEKD2_1_GC_3, HEF_MEKD2_1_GC_4, HEF_MEKD2_1_GC_5, HEF_MEKD2_1_GC_181, HEF_MEKD2_1_GC_182, HEF_MEKD2_1_GC_187, HEF_MEKD2_1_GC_188;
	
	// Coming extra from spin 0
	std::complex<double> HEF_MEKD2_1_GC_13, HEF_MEKD2_1_GC_14, HEF_MEKD2_1_GC_15, HEF_MEKD2_1_GC_18, HEF_MEKD2_1_GC_19, HEF_MEKD2_1_GC_22, HEF_MEKD2_1_GC_23, HEF_MEKD2_1_GC_25, HEF_MEKD2_1_GC_126, HEF_MEKD2_1_GC_127, HEF_MEKD2_1_GC_136, HEF_MEKD2_1_GC_137;
	
	// Coming extra from spin 1
	std::complex<double> HEF_MEKD2_1_GC_1, HEF_MEKD2_1_GC_2, HEF_MEKD2_1_GC_108, HEF_MEKD2_1_GC_109, HEF_MEKD2_1_GC_110, HEF_MEKD2_1_GC_111, HEF_MEKD2_1_GC_128, HEF_MEKD2_1_GC_129, HEF_MEKD2_1_GC_130, HEF_MEKD2_1_GC_131, HEF_MEKD2_1_GC_138, HEF_MEKD2_1_GC_139, HEF_MEKD2_1_GC_140, HEF_MEKD2_1_GC_141, HEF_MEKD2_1_GC_148, HEF_MEKD2_1_GC_149, HEF_MEKD2_1_GC_150, HEF_MEKD2_1_GC_151;
	
	// Coming extra from spin 2
	std::complex<double> HEF_MEKD2_1_GC_62, HEF_MEKD2_1_GC_63, HEF_MEKD2_1_GC_64, HEF_MEKD2_1_GC_67, HEF_MEKD2_1_GC_68, HEF_MEKD2_1_GC_71, HEF_MEKD2_1_GC_72, HEF_MEKD2_1_GC_75, HEF_MEKD2_1_GC_76, HEF_MEKD2_1_GC_79, HEF_MEKD2_1_GC_80, HEF_MEKD2_1_GC_81, HEF_MEKD2_1_GC_82, HEF_MEKD2_1_GC_83, HEF_MEKD2_1_GC_84, HEF_MEKD2_1_GC_85, HEF_MEKD2_1_GC_86, HEF_MEKD2_1_GC_90, HEF_MEKD2_1_GC_91, HEF_MEKD2_1_GC_92, HEF_MEKD2_1_GC_112, HEF_MEKD2_1_GC_113, HEF_MEKD2_1_GC_114, HEF_MEKD2_1_GC_115, HEF_MEKD2_1_GC_132, HEF_MEKD2_1_GC_133, HEF_MEKD2_1_GC_134, HEF_MEKD2_1_GC_135, HEF_MEKD2_1_GC_142, HEF_MEKD2_1_GC_143, HEF_MEKD2_1_GC_144, HEF_MEKD2_1_GC_145, HEF_MEKD2_1_GC_152, HEF_MEKD2_1_GC_153, HEF_MEKD2_1_GC_154, HEF_MEKD2_1_GC_155;
	
	
    // Model couplings dependent on aS


    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader_MEKD& slha); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 


  private:
    static Parameters_MEKD * instance; 
}; 

#endif  // Parameters_MEKD_H

