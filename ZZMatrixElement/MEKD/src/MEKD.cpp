/*************************************************************************
*  Authors:   MEKD fans
*  More info: http://mekd.ihepa.ufl.edu
*  Contact:   mekd@phys.ufl.edu
*************************************************************************/
#ifndef MEKD_MEKD_cpp
#define MEKD_MEKD_cpp

/// MEKD header
#include "../interface/MEKD.h"

using namespace std;


//////////////////////////////////////////////////////////////////////////
///  MEKD class member implementation
///
///  Provides necessary interface to the MadGraph-based ME calculator
///  and computes MEs and KDs for the process specified by the user.
///
//////////////////////////////////////////////////////////////////////////


///------------------------------------------------------------------------
/// MEKD::MEKD - a default constructor
///------------------------------------------------------------------------
MEKD::MEKD(double collisionEnergy, string PDFName)
{
	m_collisionEnergy = collisionEnergy;
	m_PDFName = PDFName;
	
	MEKD_MG_Calc.Sqrt_s = m_collisionEnergy*1000;	// translate TeV to GeV
	
	m_Mixing_Coefficients_Spin0 = MEKD_MG_Calc.Mixing_Coefficients_Spin0;
	m_Mixing_Coefficients_Spin1 = MEKD_MG_Calc.Mixing_Coefficients_Spin1;
	m_Mixing_Coefficients_Spin2 = MEKD_MG_Calc.Mixing_Coefficients_Spin2;
	
	ME_ZZ = 0;
	ME_Spin0PSMH = 0;
	ME_Spin0Ph = 0;
	ME_Spin0M = 0;
	ME_Spin1P = 0;
	ME_Spin1M = 0;
	ME_ggSpin2Pm = 0;
	ME_qqSpin2Pm = 0;
	
	four_particle_IDs_i.resize( 4, 0 );
	four_particle_Ps_i.resize( 4, NULL );
}


///------------------------------------------------------------------------
/// MEKD::processParameters - sanity check for internal parameters
///------------------------------------------------------------------------
/// Might be merged to the constructor if these flags are not used as is(was)
///------------------------------------------------------------------------
int MEKD::processParameters()
{
    /// Check if the PDF name is supported and set PDF flag
    if (m_PDFName!="CTEQ6L" && m_PDFName!="" && m_PDFName!="no PDFs") return ERR_PDFS;
    m_usePDF = (m_PDFName=="CTEQ6L");
	
	MEKD_MG_Calc.Use_PDF_w_pT0 = m_usePDF;
	
    /// Check if sqrt(s) is 7 or 8 TeV
	if (m_collisionEnergy!=7 && m_collisionEnergy!=8 && m_collisionEnergy!=13) cerr << "WARNING! You have set energy to be " << m_collisionEnergy << " TeV\n";
	
    return 0;
}


///------------------------------------------------------------------------
/// MEKD::setProcessName - sanity check and setting of process names
///------------------------------------------------------------------------
int MEKD::setProcessName(string process)
{
/// Check if process is supported, translation of namings
	if( process=="Custom" )														{m_process="Custom"; }	// Parameter-card-defined model
	
	else if( process=="qqZZ" || process=="ZZ" )									{m_process="qqZZ"; }	// "Background"
	else if( process=="qqDY" )													{m_process="qqDY"; }	// "Background"
	else if( process=="DY" )													{m_process="DY"; }
	
	else if( process=="Higgs" || process=="SMHiggs" || process=="ggSpin0Pm" )	{m_process="ggSpin0Pm"; }	// Spin-0 models
	else if( process==  "Spin0Pm" )												{m_process=  "Spin0Pm"; }
	else if( process=="CP-odd" || process=="Higgs0M" || process=="ggSpin0M" )	{m_process="ggSpin0M"; }
	else if( process==  "Spin0M" )												{m_process=  "Spin0M"; }
	else if( process=="ggSpin0PH" || process=="ggSpin0Ph" )						{m_process="ggSpin0Ph"; }
	else if( process==  "Spin0Ph" )												{m_process=  "Spin0Ph";}
	else if( process=="ggSpin0" )												{m_process="ggSpin0"; }
	else if( process==  "Spin0" )												{m_process=  "Spin0"; }
	
	else if( process=="qqSpin1P" )												{m_process="qqSpin1P"; }	// Spin-1 models
	else if( process==  "Spin1P" )												{m_process=  "Spin1P"; }
	else if( process=="qqSpin1M" )												{m_process="qqSpin1M"; }
	else if( process==  "Spin1M" )												{m_process=  "Spin1M"; }
	else if( process=="qqSpin1" )												{m_process="qqSpin1"; }
	else if( process==  "Spin1" )												{m_process=  "Spin1"; }
	
	else if( process=="ggSpin2PM" || process=="Graviton2PM" ||
		process=="ggSpin2Pm" )													{m_process="ggSpin2Pm"; }	// Spin-2 models
	else if( process=="qqSpin2PM" || process=="qqGraviton2PM" ||
		process=="qqSpin2Pm" )													{m_process="qqSpin2Pm"; }
	else if( process==  "Spin2Pm" )												{m_process=  "Spin2Pm"; }
	else if( process=="ggSpin2Ph" )												{m_process="ggSpin2Ph"; }
	else if( process=="qqSpin2Ph" )												{m_process="qqSpin2Ph"; }
	else if( process==  "Spin2Ph" )												{m_process=  "Spin2Ph"; }
	else if( process=="ggSpin2Mh" )												{m_process="ggSpin2Mh"; }
	else if( process=="qqSpin2Mh" )												{m_process="qqSpin2Mh"; }
	else if( process==  "Spin2Mh" )												{m_process=  "Spin2Mh"; }
	else if( process=="ggSpin2Pb" )												{m_process="ggSpin2Pb"; }
	else if( process=="qqSpin2Pb" )												{m_process="qqSpin2Pb"; }
	else if( process==  "Spin2Pb" )												{m_process=  "Spin2Pb"; }
	else if( process=="ggSpin2Ph2" )											{m_process="ggSpin2Ph2"; }
	else if( process=="qqSpin2Ph2" )											{m_process="qqSpin2Ph2"; }
	else if( process==  "Spin2Ph2" )											{m_process=  "Spin2Ph2"; }
	else if( process=="ggSpin2Ph3" )											{m_process="ggSpin2Ph3"; }
	else if( process=="qqSpin2Ph3" )											{m_process="qqSpin2Ph3"; }
	else if( process==  "Spin2Ph3" )											{m_process=  "Spin2Ph3"; }
	else if( process=="ggSpin2Ph6" )											{m_process="ggSpin2Ph6"; }
	else if( process=="qqSpin2Ph6" )											{m_process="qqSpin2Ph6"; }
	else if( process==  "Spin2Ph6" )											{m_process=  "Spin2Ph6"; }
	else if( process=="ggSpin2Ph7" )											{m_process="ggSpin2Ph7"; }
	else if( process=="qqSpin2Ph7" )											{m_process="qqSpin2Ph7"; }
	else if( process==  "Spin2Ph7" )											{m_process=  "Spin2Ph7"; }
	else if( process=="ggSpin2Mh9" )											{m_process="ggSpin2Mh9"; }
	else if( process=="qqSpin2Mh9" )											{m_process="qqSpin2Mh9"; }
	else if( process==  "Spin2Mh9" )											{m_process=  "Spin2Mh9"; }
	else if( process=="ggSpin2Mh10" )											{m_process="ggSpin2Mh10"; }
	else if( process=="qqSpin2Mh10" )											{m_process="qqSpin2Mh10"; }
	else if( process==  "Spin2Mh10" )											{m_process=  "Spin2Mh10"; }
	else if( process=="ggSpin2" )												{m_process="ggSpin2"; }
	else if( process=="qqSpin2" )												{m_process="qqSpin2"; }
	else if( process==  "Spin2" )												{m_process=  "Spin2"; }
	
	else if( process=="ggSpin0Pm_2f" )											{m_process="ggSpin0Pm_2f"; }	// Spin-0-to-2-leptons models
	else if( process==  "Spin0Pm_2f" )											{m_process=  "Spin0Pm_2f"; }
	else if( process=="ggSpin0M_2f" )											{m_process="ggSpin0M_2f"; }
	else if( process==  "Spin0M_2f" )											{m_process=  "Spin0M_2f"; }
// 	else if( process=="ggSpin0_2f" )											{m_process="ggSpin0_2f"; }
// 	else if( process==  "Spin0_2f" )											{m_process=  "Spin0_2f"; }
	
	else if( process=="qqSpin1P_2f" )											{m_process="qqSpin1P_2f"; }	// Spin-1-to-2-leptons models
	else if( process==  "Spin1P_2f" )											{m_process=  "Spin1P_2f"; }
	else if( process=="qqSpin1M_2f" )											{m_process="qqSpin1M_2f"; }
	else if( process==  "Spin1M_2f" )											{m_process=  "Spin1M_2f"; }
// 	else if( process=="qqSpin1_2f" )											{m_process="qqSpin1_2f"; }
// 	else if( process==  "Spin1_2f" )											{m_process=  "Spin1_2f"; }
	
	else if( process=="ggSpin2Pm_2f" )											{m_process="ggSpin2Pm_2f"; }	// Spin-2-to-2-leptons models
	else if( process=="qqSpin2Pm_2f" )											{m_process="qqSpin2Pm_2f"; }
	else if( process==  "Spin2Pm_2f" )											{m_process=  "Spin2Pm_2f"; }
// 	else if( process=="ggSpin2_2f" )											{m_process="ggSpin2_2f"; }
// 	else if( process=="qqSpin2_2f" )											{m_process="qqSpin2_2f"; }
// 	else if( process==  "Spin2_2f" )											{m_process=  "Spin2_2f"; }
	
	else if( process=="qqZ4l_Background" )										{m_process="qqZ4l_Background"; }	// Z->4l models
	else if( process=="qqZ4l_Signal" )											{m_process="qqZ4l_Signal"; }
	else return ERR_PROCESS;
	
	return 0;
}


///------------------------------------------------------------------------
/// MEKD::setProcessNames - sanity check and setting of process names
///------------------------------------------------------------------------
int MEKD::setProcessNames(string processA, string processB) {
	/// processes A and B should be different
	if( processA == processB ) return ERR_PROCESS;
	/// check if processA is supported, translation of namings
	if( processA=="Custom" )													{m_processA="Custom"; }	// Parameter-card-defined model
	
	else if( processA=="qqZZ" || processA=="ZZ" )								{m_processA="qqZZ"; }	// "Background"
	else if( processA=="qqDY" )													{m_processA="qqDY"; }	// "Background"
	else if( processA=="DY" )													{m_processA="DY"; }
	
	else if( processA=="Higgs" || processA=="SMHiggs" || processA=="ggSpin0Pm" )	{m_processA="ggSpin0Pm"; }	// Spin-0 models
	else if( processA==  "Spin0Pm" )												{m_processA=  "Spin0Pm"; }
	else if( processA=="CP-odd" || processA=="Higgs0M" || processA=="ggSpin0M" )	{m_processA="ggSpin0M"; }
	else if( processA==  "Spin0M" )												{m_processA=  "Spin0M"; }
	else if( processA=="ggSpin0PH" || processA=="ggSpin0Ph" )					{m_processA="ggSpin0Ph"; }
	else if( processA==  "Spin0Ph" )											{m_processA=  "Spin0Ph";}
	else if( processA=="ggSpin0" )												{m_processA="ggSpin0"; }
	else if( processA==  "Spin0" )												{m_processA=  "Spin0"; }
	
	else if( processA=="qqSpin1P" )												{m_processA="qqSpin1P"; }	// Spin-1 models
	else if( processA==  "Spin1P" )												{m_processA=  "Spin1P"; }
	else if( processA=="qqSpin1M" )												{m_processA="qqSpin1M"; }
	else if( processA==  "Spin1M" )												{m_processA=  "Spin1M"; }
	else if( processA=="qqSpin1" )												{m_processA="qqSpin1"; }
	else if( processA==  "Spin1" )												{m_processA=  "Spin1"; }
	
	else if( processA=="ggSpin2PM" || processA=="Graviton2PM" ||
		processA=="ggSpin2Pm" )													{m_processA="ggSpin2Pm"; }	// Spin-2 models
	else if( processA=="qqSpin2PM" || processA=="qqGraviton2PM" ||
		processA=="qqSpin2Pm" )													{m_processA="qqSpin2Pm"; }
	else if( processA==  "Spin2Pm" )											{m_processA=  "Spin2Pm"; }
	else if( processA=="ggSpin2Ph" )											{m_processA="ggSpin2Ph"; }
	else if( processA=="qqSpin2Ph" )											{m_processA="qqSpin2Ph"; }
	else if( processA==  "Spin2Ph" )											{m_processA=  "Spin2Ph"; }
	else if( processA=="ggSpin2Mh" )											{m_processA="ggSpin2Mh"; }
	else if( processA=="qqSpin2Mh" )											{m_processA="qqSpin2Mh"; }
	else if( processA==  "Spin2Mh" )											{m_processA=  "Spin2Mh"; }
	else if( processA=="ggSpin2Pb" )											{m_processA="ggSpin2Pb"; }
	else if( processA=="qqSpin2Pb" )											{m_processA="qqSpin2Pb"; }
	else if( processA==  "Spin2Pb" )											{m_processA=  "Spin2Pb"; }
	else if( processA=="ggSpin2Ph2" )											{m_processA="ggSpin2Ph2"; }
	else if( processA=="qqSpin2Ph2" )											{m_processA="qqSpin2Ph2"; }
	else if( processA==  "Spin2Ph2" )											{m_processA=  "Spin2Ph2"; }
	else if( processA=="ggSpin2Ph3" )											{m_processA="ggSpin2Ph3"; }
	else if( processA=="qqSpin2Ph3" )											{m_processA="qqSpin2Ph3"; }
	else if( processA==  "Spin2Ph3" )											{m_processA=  "Spin2Ph3"; }
	else if( processA=="ggSpin2Ph6" )											{m_processA="ggSpin2Ph6"; }
	else if( processA=="qqSpin2Ph6" )											{m_processA="qqSpin2Ph6"; }
	else if( processA==  "Spin2Ph6" )											{m_processA=  "Spin2Ph6"; }
	else if( processA=="ggSpin2Ph7" )											{m_processA="ggSpin2Ph7"; }
	else if( processA=="qqSpin2Ph7" )											{m_processA="qqSpin2Ph7"; }
	else if( processA==  "Spin2Ph7" )											{m_processA=  "Spin2Ph7"; }
	else if( processA=="ggSpin2Mh9" )											{m_processA="ggSpin2Mh9"; }
	else if( processA=="qqSpin2Mh9" )											{m_processA="qqSpin2Mh9"; }
	else if( processA==  "Spin2Mh9" )											{m_processA=  "Spin2Mh9"; }
	else if( processA=="ggSpin2Mh10" )											{m_processA="ggSpin2Mh10"; }
	else if( processA=="qqSpin2Mh10" )											{m_processA="qqSpin2Mh10"; }
	else if( processA==  "Spin2Mh10" )											{m_processA=  "Spin2Mh10"; }
	else if( processA=="ggSpin2" )												{m_processA="ggSpin2"; }
	else if( processA=="qqSpin2" )												{m_processA="qqSpin2"; }
	else if( processA==  "Spin2" )												{m_processA=  "Spin2"; }
	
	else if( processA=="ggSpin0Pm_2f" )											{m_processA="ggSpin0Pm_2f"; }	// Spin-0-to-2-leptons models
	else if( processA==  "Spin0Pm_2f" )											{m_processA=  "Spin0Pm_2f"; }
	else if( processA=="ggSpin0M_2f" )											{m_processA="ggSpin0M_2f"; }
	else if( processA==  "Spin0M_2f" )											{m_processA=  "Spin0M_2f"; }
// 	else if( processA=="ggSpin0_2f" )											{m_processA="ggSpin0_2f"; }
// 	else if( processA==  "Spin0_2f" )											{m_processA=  "Spin0_2f"; }
	
	else if( processA=="qqSpin1P_2f" )											{m_processA="qqSpin1P_2f"; }	// Spin-1-to-2-leptons models
	else if( processA==  "Spin1P_2f" )											{m_processA=  "Spin1P_2f"; }
	else if( processA=="qqSpin1M_2f" )											{m_processA="qqSpin1M_2f"; }
	else if( processA==  "Spin1M_2f" )											{m_processA=  "Spin1M_2f"; }
// 	else if( processA=="qqSpin1_2f" )											{m_processA="qqSpin1_2f"; }
// 	else if( processA==  "Spin1_2f" )											{m_processA=  "Spin1_2f"; }
	
	else if( processA=="ggSpin2Pm_2f" )											{m_processA="ggSpin2Pm_2f"; }	// Spin-2-to-2-leptons models
	else if( processA=="qqSpin2Pm_2f" )											{m_processA="qqSpin2Pm_2f"; }
	else if( processA==  "Spin2Pm_2f" )											{m_processA=  "Spin2Pm_2f"; }
// 	else if( processA=="ggSpin2_2f" )											{m_processA="ggSpin2_2f"; }
// 	else if( processA=="qqSpin2_2f" )											{m_processA="qqSpin2_2f"; }
// 	else if( processA==  "Spin2_2f" )											{m_processA=  "Spin2_2f"; }
	
	else if( processA=="qqZ4l_Background" )										{m_processA="qqZ4l_Background"; }	// Z->4l models
	else if( processA=="qqZ4l_Signal" )											{m_processA="qqZ4l_Signal"; }
	else return ERR_PROCESS;
	
	/// check if processB is supported, translation of namings
	if( processB=="Custom" )													{m_processB="Custom"; }	// Parameter-card-defined model
	
	else if( processB=="qqZZ" || processB=="ZZ" )								{m_processB="qqZZ"; }	// "Background"
	else if( processB=="qqDY" )													{m_processB="qqDY"; }	// "Background"
	else if( processB=="DY" )													{m_processB="DY"; }
	
	else if( processB=="Higgs" || processB=="SMHiggs" || processB=="ggSpin0Pm" )	{m_processB="ggSpin0Pm"; }	// Spin-0 models
	else if( processB==  "Spin0Pm" )												{m_processB=  "Spin0Pm"; }
	else if( processB=="CP-odd" || processB=="Higgs0M" || processB=="ggSpin0M" )	{m_processB="ggSpin0M"; }
	else if( processB==  "Spin0M" )												{m_processB=  "Spin0M"; }
	else if( processB=="ggSpin0PH" || processB=="ggSpin0Ph" )					{m_processB="ggSpin0Ph"; }
	else if( processB==  "Spin0Ph" )											{m_processB=  "Spin0Ph";}
	else if( processB=="ggSpin0" )												{m_processB="ggSpin0"; }
	else if( processB==  "Spin0" )												{m_processB=  "Spin0"; }
	
	else if( processB=="qqSpin1P" )												{m_processB="qqSpin1P"; }	// Spin-1 models
	else if( processB==  "Spin1P" )												{m_processB=  "Spin1P"; }
	else if( processB=="qqSpin1M" )												{m_processB="qqSpin1M"; }
	else if( processB==  "Spin1M" )												{m_processB=  "Spin1M"; }
	else if( processB=="qqSpin1" )												{m_processB="qqSpin1"; }
	else if( processB==  "Spin1" )												{m_processB=  "Spin1"; }
	
	else if( processB=="ggSpin2PM" || processB=="Graviton2PM" ||
		processB=="ggSpin2Pm" )													{m_processB="ggSpin2Pm"; }	// Spin-2 models
	else if( processB=="qqSpin2PM" || processB=="qqGraviton2PM" ||
		processB=="qqSpin2Pm" )													{m_processB="qqSpin2Pm"; }
	else if( processB==  "Spin2Pm" )											{m_processB=  "Spin2Pm"; }
	else if( processB=="ggSpin2Ph" )											{m_processB="ggSpin2Ph"; }
	else if( processB=="qqSpin2Ph" )											{m_processB="qqSpin2Ph"; }
	else if( processB==  "Spin2Ph" )											{m_processB=  "Spin2Ph"; }
	else if( processB=="ggSpin2Mh" )											{m_processB="ggSpin2Mh"; }
	else if( processB=="qqSpin2Mh" )											{m_processB="qqSpin2Mh"; }
	else if( processB==  "Spin2Mh" )											{m_processB=  "Spin2Mh"; }
	else if( processB=="ggSpin2Pb" )											{m_processB="ggSpin2Pb"; }
	else if( processB=="qqSpin2Pb" )											{m_processB="qqSpin2Pb"; }
	else if( processB==  "Spin2Pb" )											{m_processB=  "Spin2Pb"; }
	else if( processB=="ggSpin2Ph2" )											{m_processB="ggSpin2Ph2"; }
	else if( processB=="qqSpin2Ph2" )											{m_processB="qqSpin2Ph2"; }
	else if( processB==  "Spin2Ph2" )											{m_processB=  "Spin2Ph2"; }
	else if( processB=="ggSpin2Ph3" )											{m_processB="ggSpin2Ph3"; }
	else if( processB=="qqSpin2Ph3" )											{m_processB="qqSpin2Ph3"; }
	else if( processB==  "Spin2Ph3" )											{m_processB=  "Spin2Ph3"; }
	else if( processB=="ggSpin2Ph6" )											{m_processB="ggSpin2Ph6"; }
	else if( processB=="qqSpin2Ph6" )											{m_processB="qqSpin2Ph6"; }
	else if( processB==  "Spin2Ph6" )											{m_processB=  "Spin2Ph6"; }
	else if( processB=="ggSpin2Ph7" )											{m_processB="ggSpin2Ph7"; }
	else if( processB=="qqSpin2Ph7" )											{m_processB="qqSpin2Ph7"; }
	else if( processB==  "Spin2Ph7" )											{m_processB=  "Spin2Ph7"; }
	else if( processB=="ggSpin2Mh9" )											{m_processB="ggSpin2Mh9"; }
	else if( processB=="qqSpin2Mh9" )											{m_processB="qqSpin2Mh9"; }
	else if( processB==  "Spin2Mh9" )											{m_processB=  "Spin2Mh9"; }
	else if( processB=="ggSpin2Mh10" )											{m_processB="ggSpin2Mh10"; }
	else if( processB=="qqSpin2Mh10" )											{m_processB="qqSpin2Mh10"; }
	else if( processB==  "Spin2Mh10" )											{m_processB=  "Spin2Mh10"; }
	else if( processB=="ggSpin2" )												{m_processB="ggSpin2"; }
	else if( processB=="qqSpin2" )												{m_processB="qqSpin2"; }
	else if( processB==  "Spin2" )												{m_processB=  "Spin2"; }
	
	else if( processB=="ggSpin0Pm_2f" )											{m_processB="ggSpin0Pm_2f"; }	// Spin-0-to-2-leptons models
	else if( processB==  "Spin0Pm_2f" )											{m_processB=  "Spin0Pm_2f"; }
	else if( processB=="ggSpin0M_2f" )											{m_processB="ggSpin0M_2f"; }
	else if( processB==  "Spin0M_2f" )											{m_processB=  "Spin0M_2f"; }
// 	else if( processB=="ggSpin0_2f" )											{m_processB="ggSpin0_2f"; }
// 	else if( processB==  "Spin0_2f" )											{m_processB=  "Spin0_2f"; }
	
	else if( processB=="qqSpin1P_2f" )											{m_processB="qqSpin1P_2f"; }	// Spin-1-to-2-leptons models
	else if( processB==  "Spin1P_2f" )											{m_processB=  "Spin1P_2f"; }
	else if( processB=="qqSpin1M_2f" )											{m_processB="qqSpin1M_2f"; }
	else if( processB==  "Spin1M_2f" )											{m_processB=  "Spin1M_2f"; }
// 	else if( processB=="qqSpin1_2f" )											{m_processB="qqSpin1_2f"; }
// 	else if( processB==  "Spin1_2f" )											{m_processB=  "Spin1_2f"; }
	
	else if( processB=="ggSpin2Pm_2f" )											{m_processB="ggSpin2Pm_2f"; }	// Spin-2-to-2-leptons models
	else if( processB=="qqSpin2Pm_2f" )											{m_processB="qqSpin2Pm_2f"; }
	else if( processB==  "Spin2Pm_2f" )											{m_processB=  "Spin2Pm_2f"; }
// 	else if( processB=="ggSpin2_2f" )											{m_processB="ggSpin2_2f"; }
// 	else if( processB=="qqSpin2_2f" )											{m_processB="qqSpin2_2f"; }
// 	else if( processB==  "Spin2_2f" )											{m_processB=  "Spin2_2f"; }
	
	else if( processB=="qqZ4l_Background" )										{m_processB="qqZ4l_Background"; }	// Z->4l models
	else if( processB=="qqZ4l_Signal" )											{m_processB="qqZ4l_Signal"; }
    else return ERR_PROCESS;
	
    return 0;
}


///------------------------------------------------------------------------
/// MEKD::computeKD - compute KD from precalculated MEs for input prcesses A and B
///------------------------------------------------------------------------
int MEKD::computeKD( string processA, string processB,
					double& kd,
					double& me2processA,
					double& me2processB )
{
	/// Sanity check for input process names
	if( (buffer_int=setProcessNames(processA, processB)) != 0 ) return buffer_int;
	
	/// Looking for the precalculated MEs
	if( ME_ZZ == 0 ) { cerr << "ERROR! The requested process has not been precalculated.\n" ; return ERR_PROCESS; }
	else if( m_processA=="ZZ" )			me2processA = ME_ZZ;
	else if( m_processA=="ggSpin0Pm" )	me2processA = ME_Spin0PSMH;
	else if( m_processA=="ggSpin0M" )	me2processA = ME_Spin0M;
	else if( m_processA=="ggSpin0Ph" )	me2processA = ME_Spin0Ph;
	else if( m_processA=="qqSpin1P" )	me2processA = ME_Spin1P;
	else if( m_processA=="qqSpin1M" )	me2processA = ME_Spin1M;
	else if( m_processA=="ggSpin2Pm" )	me2processA = ME_ggSpin2Pm;
	else if( m_processA=="qqSpin2Pm" )	me2processA = ME_qqSpin2Pm;
	else { cerr << "ERROR! The requested process has not been precalculated.\n" ; return ERR_PROCESS; }
	
	
	/// Looking for the precalculated MEs
	if( ME_ZZ == 0 ) { cerr << "ERROR! The requested process has not been precalculated.\n" ; return ERR_PROCESS; }
	else if( m_processB=="ZZ" )			me2processB = ME_ZZ;
	else if( m_processB=="ggSpin0Pm" )	me2processB = ME_Spin0PSMH;
	else if( m_processB=="ggSpin0M" )	me2processB = ME_Spin0M;
	else if( m_processB=="ggSpin0Ph" )	me2processB = ME_Spin0Ph;
	else if( m_processB=="qqSpin1P" )	me2processB = ME_Spin1P;
	else if( m_processB=="qqSpin1M" )	me2processB = ME_Spin1M;
	else if( m_processB=="ggSpin2Pm" )	me2processB = ME_ggSpin2Pm;
	else if( m_processB=="qqSpin2Pm" )	me2processB = ME_qqSpin2Pm;
	else { cerr << "ERROR! The requested process has not been precalculated.\n" ; return ERR_PROCESS; }
	
	
	/// Build Kinematic Discriminant (KD) as a ratio of logs of MEs
	kd = log(me2processA / me2processB);
	
	return 0;
}


///------------------------------------------------------------------------
/// MEKD::computeKD - compute KD and MEs for input prcesses A and B
///------------------------------------------------------------------------
int MEKD::computeKD( string processA, string processB,
                     double lept1P[], int lept1Id, double lept2P[], int lept2Id,
                     double lept3P[], int lept3Id, double lept4P[], int lept4Id,
                     double& kd, double& me2processA, double& me2processB )
{
	/// Load internal containers
	four_particle_Ps_i[0] = lept1P;
	four_particle_Ps_i[1] = lept2P;
	four_particle_Ps_i[2] = lept3P;
	four_particle_Ps_i[3] = lept4P;
	
	four_particle_IDs_i[0] = lept1Id;
	four_particle_IDs_i[1] = lept2Id;
	four_particle_IDs_i[2] = lept3Id;
	four_particle_IDs_i[3] = lept4Id;
	
	return computeKD( processA, processB, four_particle_Ps_i, four_particle_IDs_i, kd, me2processA, me2processB );
}


///------------------------------------------------------------------------
/// MEKD::computeKD - compute KD and MEs for the input processes A and B
///------------------------------------------------------------------------
int MEKD::computeKD( string processA, string processB,
                     vector<double*> input_Ps, vector<int> input_IDs,
                     double& kd, double& me2processA, double& me2processB )
{
	/// Checks input for compatibility
	if( input_Ps.size() != input_IDs.size() ) return ERR_INPUT;
	if( input_Ps.size() < 2 || input_Ps.size() > 5 ) return ERR_INPUT;
	
	/// Set an expected resonance decay mode
	if( input_Ps.size() == 2 || input_Ps.size() == 3 ) MEKD_MG_Calc.Resonance_decay_mode="2mu";
	else MEKD_MG_Calc.Resonance_decay_mode="ZZ";
	
    /// Sanity check for input process names and internal parameters
    if( (buffer_int=setProcessNames(processA, processB)) != 0 ) return buffer_int;
    if( (buffer_int=processParameters()) != 0 ) return buffer_int;

	/// Set input-particle kinematics
	MEKD_MG_Calc.p1 = input_Ps[0];
	MEKD_MG_Calc.id1 = input_IDs[0];
	MEKD_MG_Calc.p2 = input_Ps[1];
	MEKD_MG_Calc.id2 = input_IDs[1];
	if( input_IDs.size() >= 3 )
	{
		MEKD_MG_Calc.p3 = input_Ps[2];
		MEKD_MG_Calc.id3 = input_IDs[2];
	}
	else MEKD_MG_Calc.id3 = 0;
	if( input_IDs.size() >= 4 )
	{
		MEKD_MG_Calc.p4 = input_Ps[3];
		MEKD_MG_Calc.id4 = input_IDs[3];
	}
	else MEKD_MG_Calc.id4 = 0;
	if( input_IDs.size() >= 5 )
	{
		MEKD_MG_Calc.p5 = input_Ps[4];
		MEKD_MG_Calc.id5 = input_IDs[4];
	}
	else MEKD_MG_Calc.id5 = 0;

    /// Compute ME for process A only (e.g. signal 1)
    buffer_int = MEKD_MG_Calc.Run_MEKD_MG( m_processA );
    /// Get ME for process A
    me2processA = MEKD_MG_Calc.Signal_ME;
	
    /// Compute ME for process B only (e.g. signal 2 or background)
    buffer_int = MEKD_MG_Calc.Run_MEKD_MG( m_processB );
    /// Get ME for process B
    me2processB = MEKD_MG_Calc.Signal_ME;
    /// Build Kinematic Discriminant (KD) as a ratio of logs of MEs
    kd = log(me2processA / me2processB);
	
	return buffer_int;
}


///------------------------------------------------------------------------
/// MEKD::computeME - compute ME for the input process
///------------------------------------------------------------------------
int MEKD::computeME( string processName,
					vector<double*> input_Ps, vector<int> input_IDs,
					double& me2process )
{
	/// Checks input for compatibility
	if( input_Ps.size() != input_IDs.size() ) return ERR_INPUT;
	if( input_Ps.size() < 2 || input_Ps.size() > 5 ) return ERR_INPUT;
	
	/// Set an expected resonance decay mode
	if( input_Ps.size() == 2 || input_Ps.size() == 3 ) MEKD_MG_Calc.Resonance_decay_mode="2mu";
	else MEKD_MG_Calc.Resonance_decay_mode="ZZ";
	
    /// Sanity check for input process names and internal parameters
    if( (buffer_int=setProcessName(processName)) != 0 ) return buffer_int;
    if( (buffer_int=processParameters()) != 0 ) return buffer_int;

	/// Set input-particle kinematics
	MEKD_MG_Calc.p1 = input_Ps[0];
	MEKD_MG_Calc.id1 = input_IDs[0];
	MEKD_MG_Calc.p2 = input_Ps[1];
	MEKD_MG_Calc.id2 = input_IDs[1];
	if( input_IDs.size() >= 3 )
	{
		MEKD_MG_Calc.p3 = input_Ps[2];
		MEKD_MG_Calc.id3 = input_IDs[2];
	}
	else MEKD_MG_Calc.id3 = 0;
	if( input_IDs.size() >= 4 )
	{
		MEKD_MG_Calc.p4 = input_Ps[3];
		MEKD_MG_Calc.id4 = input_IDs[3];
	}
	else MEKD_MG_Calc.id4 = 0;
	if( input_IDs.size() >= 5 )
	{
		MEKD_MG_Calc.p5 = input_Ps[4];
		MEKD_MG_Calc.id5 = input_IDs[4];
	}
	else MEKD_MG_Calc.id5 = 0;
	
	/// Compute ME for the process (e.g. signal 1)
	buffer_int = MEKD_MG_Calc.Run_MEKD_MG( m_process );
	/// Get ME for the process
	me2process = MEKD_MG_Calc.Signal_ME;
	
	return buffer_int;
}


///------------------------------------------------------------------------
/// MEKD::computeMEs - compute MEs for a multiple reuse
///------------------------------------------------------------------------
int MEKD::computeMEs( double lept1P[], int lept1Id, double lept2P[], int lept2Id,
					double lept3P[], int lept3Id, double lept4P[], int lept4Id )
{
	/// Load internal containers
	four_particle_Ps_i[0] = lept1P;
	four_particle_Ps_i[1] = lept2P;
	four_particle_Ps_i[2] = lept3P;
	four_particle_Ps_i[3] = lept4P;
	
	four_particle_IDs_i[0] = lept1Id;
	four_particle_IDs_i[1] = lept2Id;
	four_particle_IDs_i[2] = lept3Id;
	four_particle_IDs_i[3] = lept4Id;
	
	return computeMEs( four_particle_Ps_i, four_particle_IDs_i );
}


///------------------------------------------------------------------------
/// MEKD::computeMEs - compute MEs for a multiple reuse
///------------------------------------------------------------------------
int MEKD::computeMEs( vector<double*> input_Ps, vector<int> input_IDs )
{
	/// Checks input for compatibility
	if( input_Ps.size() != input_IDs.size() ) return ERR_INPUT;
	if( input_Ps.size() < 2 || input_Ps.size() > 5 ) return ERR_INPUT;
	
	/// Set an expected resonance decay mode
	if( input_Ps.size() == 2 || input_Ps.size() == 3 ) MEKD_MG_Calc.Resonance_decay_mode="2mu";
	else MEKD_MG_Calc.Resonance_decay_mode="ZZ";
	
	/// Sanity check for internal parameters
	if( (buffer_int=processParameters()) != 0 ) return buffer_int;
	
	/// Parameterize MEKD_MG_Calc
	if( MEKD_MG_Calc.Test_Models.size()==0 )	// Fills in intersting models to compute only once
	{
		MEKD_MG_Calc.Test_Models.push_back( "ZZ" );
		MEKD_MG_Calc.Test_Models.push_back( "ggSpin0Pm" );
		MEKD_MG_Calc.Test_Models.push_back( "ggSpin0M" );
		MEKD_MG_Calc.Test_Models.push_back( "ggSpin0Ph" );
		MEKD_MG_Calc.Test_Models.push_back( "qqSpin1P" );
		MEKD_MG_Calc.Test_Models.push_back( "qqSpin1M" );
		MEKD_MG_Calc.Test_Models.push_back( "ggSpin2Pm" );
		MEKD_MG_Calc.Test_Models.push_back( "qqSpin2Pm" );
	}
	else if( MEKD_MG_Calc.Test_Models.size()!=4 ) return ERR_PROCESS;
	
	/// Set input-particle kinematics
	MEKD_MG_Calc.p1 = input_Ps[0];
	MEKD_MG_Calc.id1 = input_IDs[0];
	MEKD_MG_Calc.p2 = input_Ps[1];
	MEKD_MG_Calc.id2 = input_IDs[1];
	if( input_IDs.size() >= 3 )
	{
		MEKD_MG_Calc.p3 = input_Ps[2];
		MEKD_MG_Calc.id3 = input_IDs[2];
	}
	else MEKD_MG_Calc.id3 = 0;
	if( input_IDs.size() >= 4 )
	{
		MEKD_MG_Calc.p4 = input_Ps[3];
		MEKD_MG_Calc.id4 = input_IDs[3];
	}
	else MEKD_MG_Calc.id4 = 0;
	if( input_IDs.size() >= 5 )
	{
		MEKD_MG_Calc.p5 = input_Ps[4];
		MEKD_MG_Calc.id5 = input_IDs[4];
	}
	else MEKD_MG_Calc.id5 = 0;
	
	/// Compute MEs
	buffer_int = MEKD_MG_Calc.Run_MEKD_MG();
	
	/// ME value readouts
	ME_ZZ = MEKD_MG_Calc.Signal_MEs[0];
	ME_Spin0PSMH = MEKD_MG_Calc.Signal_MEs[1];
	ME_Spin0M = MEKD_MG_Calc.Signal_MEs[2];
	ME_Spin0Ph = MEKD_MG_Calc.Signal_MEs[3];
	ME_Spin1P = MEKD_MG_Calc.Signal_MEs[4];
	ME_Spin1M = MEKD_MG_Calc.Signal_MEs[5];
	ME_ggSpin2Pm = MEKD_MG_Calc.Signal_MEs[6];
	ME_qqSpin2Pm = MEKD_MG_Calc.Signal_MEs[7];
	
	return buffer_int;
}


///------------------------------------------------------------------------
/// Mixed-state ME mixer of production like qqSpin1M, qqSpin1P, exotic, exotic pseudo, and decay like Spin1M, Spin1P states, corresponding couplings rhoQ11, rhoQ12, rhoQ13, rhoQ14, b1z/rhomu11, b2z/rhomu12, rhomu13, rhomu14.
///------------------------------------------------------------------------
int MEKD::Mix_Spin0( complex<double> Spin0Pm_relamp, complex<double> Spin0Ph_relamp, complex<double> Spin0Phexo_relamp, complex<double> Spin0M_relamp )
{
	m_Mixing_Coefficients_Spin0[0] = Spin0Pm_relamp;
	m_Mixing_Coefficients_Spin0[1] = Spin0Ph_relamp;
	m_Mixing_Coefficients_Spin0[2] = Spin0Phexo_relamp;
	m_Mixing_Coefficients_Spin0[3] = Spin0M_relamp;
	
	if( m_Mixing_Coefficients_Spin0 == NULL ) return ERR_OTHER;
	return 0;
}


///------------------------------------------------------------------------
/// Mixed-state ME mixer of production like qqSpin1M, qqSpin1P, exotic, exotic pseudo, and decay like Spin1M, Spin1P states, corresponding couplings rhoQ11, rhoQ12, rhoQ13, rhoQ14, b1z/rhomu11, b2z/rhomu12, rhomu13, rhomu14.
///------------------------------------------------------------------------
int MEKD::Mix_Spin1( complex<double> prod_Spin1M_relamp, complex<double> prod_Spin1P_relamp, complex<double> prod_Spin1Mexo_relamp, complex<double> prod_Spin1Pexo_relamp, complex<double> dec_Spin1M_relamp, complex<double> dec_Spin1P_relamp, complex<double> dec_Spin1_rhomu13_relamp, complex<double> dec_Spin1_rhomu14_relamp )
{
	m_Mixing_Coefficients_Spin1[0] = prod_Spin1M_relamp;
	m_Mixing_Coefficients_Spin1[1] = prod_Spin1P_relamp;
	m_Mixing_Coefficients_Spin1[2] = prod_Spin1Mexo_relamp;
	m_Mixing_Coefficients_Spin1[3] = prod_Spin1Pexo_relamp;
	m_Mixing_Coefficients_Spin1[4] = dec_Spin1M_relamp;
	m_Mixing_Coefficients_Spin1[5] = dec_Spin1M_relamp;
	m_Mixing_Coefficients_Spin1[6] = dec_Spin1_rhomu13_relamp;
	m_Mixing_Coefficients_Spin1[7] = dec_Spin1_rhomu14_relamp;
	
	if( m_Mixing_Coefficients_Spin1 == NULL ) return ERR_OTHER;
	return 0;
}


///------------------------------------------------------------------------
/// Mixed-state ME mixer of Spin-2 states. Production couplings: k1g/rhoQ21, ..., k4g/rhoQ24, ..., k10g, decay couplings: k1z/rhomu21, k4z/rhomu24, ..., k10z.
///------------------------------------------------------------------------
int MEKD::Mix_Spin2( complex<double> *prod_Spin2_relamp, complex<double> *dec_Spin2_relamp )
{
	for( buffer_uint=0; buffer_uint<10; buffer_uint++)
		m_Mixing_Coefficients_Spin2[buffer_uint] = prod_Spin2_relamp[buffer_uint];
	for( buffer_uint=10; buffer_uint<20; buffer_uint++)
		m_Mixing_Coefficients_Spin2[buffer_uint] = dec_Spin2_relamp[buffer_uint];
	
	if( m_Mixing_Coefficients_Spin2 == NULL ) return ERR_OTHER;
	return 0;
}



#if (defined(MEKD_STANDALONE) && defined(MEKD_with_ROOT)) || !(defined(MEKD_STANDALONE))
///------------------------------------------------------------------------
/// (ROOT-compatabile overloads)
///------------------------------------------------------------------------


///------------------------------------------------------------------------
/// MEKD::computeKD - compute KD and MEs for the input processes A and B
///------------------------------------------------------------------------
int MEKD::computeKD( TString processA, TString processB,
                     TLorentzVector lept1P, int lept1Id, TLorentzVector lept2P, int lept2Id,
                     TLorentzVector lept3P, int lept3Id, TLorentzVector lept4P, int lept4Id,
                     double& kd, double& me2processA, double& me2processB )
{
	/// Prepare 4-momenta in the required format
	lept1P_i[0] = lept1P.E();
	lept1P_i[1] = lept1P.Px();
	lept1P_i[2] = lept1P.Py();
	lept1P_i[3] = lept1P.Pz();
	
	lept2P_i[0] = lept2P.E();
	lept2P_i[1] = lept2P.Px();
	lept2P_i[2] = lept2P.Py();
	lept2P_i[3] = lept2P.Pz();
	
	lept3P_i[0] = lept3P.E();
	lept3P_i[1] = lept3P.Px();
	lept3P_i[2] = lept3P.Py();
	lept3P_i[3] = lept3P.Pz();
	
	lept4P_i[0] = lept4P.E();
	lept4P_i[1] = lept4P.Px();
	lept4P_i[2] = lept4P.Py();
	lept4P_i[3] = lept4P.Pz();
	
	/// Load internal containers
	four_particle_Ps_i[0] = lept1P_i;
	four_particle_Ps_i[1] = lept2P_i;
	four_particle_Ps_i[2] = lept3P_i;
	four_particle_Ps_i[3] = lept4P_i;
	
	four_particle_IDs_i[0] = lept1Id;
	four_particle_IDs_i[1] = lept2Id;
	four_particle_IDs_i[2] = lept3Id;
	four_particle_IDs_i[3] = lept4Id;
	
	return computeKD( (string) processA.Data(), (string) processB.Data(), four_particle_Ps_i, four_particle_IDs_i, kd, me2processA, me2processB );
}


///------------------------------------------------------------------------
/// MEKD::computeKD - compute KD and MEs for the input processes A and B
///------------------------------------------------------------------------
int MEKD::computeKD( TString processA, TString processB,
					vector<TLorentzVector> input_Ps, vector<int> input_IDs,
					double& kd, double& me2processA, double& me2processB )
{
	/// Resize internal vector<double*> if needed
	if( input_Ps_i.size() != input_Ps.size() )
	{
		for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ ) { delete input_Ps_i[buffer_uint]; input_Ps_i[buffer_uint]=NULL; }
		input_Ps_i.resize( input_Ps.size(), NULL );
		for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ ) { input_Ps_i[buffer_uint]=new double[4]; }
	}
	
	/// Put vector<TLorentzVector> into internal containers
	for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ )
	{
		input_Ps_i[buffer_uint][0] = input_Ps[buffer_uint].E();
		input_Ps_i[buffer_uint][1] = input_Ps[buffer_uint].Px();
		input_Ps_i[buffer_uint][2] = input_Ps[buffer_uint].Py();
		input_Ps_i[buffer_uint][3] = input_Ps[buffer_uint].Pz();
	}
	
	return computeKD( (string) processA.Data(), (string) processB.Data(), input_Ps_i, input_IDs, kd, me2processA, me2processB );
}


///------------------------------------------------------------------------
/// MEKD::computeME - compute ME for the input process
///------------------------------------------------------------------------
int MEKD::computeME( TString processName,
					vector<TLorentzVector> input_Ps, vector<int> input_IDs,
					double& me2process )
{
	/// Resize internal vector<double*> if needed
	if( input_Ps_i.size() != input_Ps.size() )
	{
		for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ ) { delete input_Ps_i[buffer_uint]; input_Ps_i[buffer_uint]=NULL; }
		input_Ps_i.resize( input_Ps.size(), NULL );
		for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ ) { input_Ps_i[buffer_uint]=new double[4]; }
	}
	
	/// Put vector<TLorentzVector> into internal containers
	for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ )
	{
		input_Ps_i[buffer_uint][0] = input_Ps[buffer_uint].E();
		input_Ps_i[buffer_uint][1] = input_Ps[buffer_uint].Px();
		input_Ps_i[buffer_uint][2] = input_Ps[buffer_uint].Py();
		input_Ps_i[buffer_uint][3] = input_Ps[buffer_uint].Pz();
	}
	
	return computeME( (string) processName.Data(), input_Ps_i, input_IDs, me2process );
}


///------------------------------------------------------------------------
/// MEKD::computeMEs - compute MEs for a multiple reuse
///------------------------------------------------------------------------
int MEKD::computeMEs( TLorentzVector lept1P, int lept1Id, TLorentzVector lept2P, int lept2Id,
					TLorentzVector lept3P, int lept3Id, TLorentzVector lept4P, int lept4Id )
{
	/// Prepare 4-momenta in the required format
	lept1P_i[0] = lept1P.E();
	lept1P_i[1] = lept1P.Px();
	lept1P_i[2] = lept1P.Py();
	lept1P_i[3] = lept1P.Pz();
	
	lept2P_i[0] = lept2P.E();
	lept2P_i[1] = lept2P.Px();
	lept2P_i[2] = lept2P.Py();
	lept2P_i[3] = lept2P.Pz();
	
	lept3P_i[0] = lept3P.E();
	lept3P_i[1] = lept3P.Px();
	lept3P_i[2] = lept3P.Py();
	lept3P_i[3] = lept3P.Pz();
	
	lept4P_i[0] = lept4P.E();
	lept4P_i[1] = lept4P.Px();
	lept4P_i[2] = lept4P.Py();
	lept4P_i[3] = lept4P.Pz();
	
	/// Load internal containers
	four_particle_Ps_i[0] = lept1P_i;
	four_particle_Ps_i[1] = lept2P_i;
	four_particle_Ps_i[2] = lept3P_i;
	four_particle_Ps_i[3] = lept4P_i;
	
	four_particle_IDs_i[0] = lept1Id;
	four_particle_IDs_i[1] = lept2Id;
	four_particle_IDs_i[2] = lept3Id;
	four_particle_IDs_i[3] = lept4Id;
	
	return computeMEs( four_particle_Ps_i, four_particle_IDs_i );
}


///------------------------------------------------------------------------
/// MEKD::computeMEs - compute MEs for a multiple reuse
///------------------------------------------------------------------------
int MEKD::computeMEs( vector<TLorentzVector> input_Ps, vector<int> input_IDs )
{
	/// Resize internal vector<double*> if needed
	if( input_Ps_i.size() != input_Ps.size() )
	{
		for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ ) { delete input_Ps_i[buffer_uint]; input_Ps_i[buffer_uint]=NULL; }
		input_Ps_i.resize( input_Ps.size(), NULL );
		for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ ) { input_Ps_i[buffer_uint]=new double[4]; }
	}
	
	/// Put vector<TLorentzVector> into internal containers
	for( buffer_uint=0; buffer_uint < input_Ps_i.size(); buffer_uint++ )
	{
		input_Ps_i[buffer_uint][0] = input_Ps[buffer_uint].E();
		input_Ps_i[buffer_uint][1] = input_Ps[buffer_uint].Px();
		input_Ps_i[buffer_uint][2] = input_Ps[buffer_uint].Py();
		input_Ps_i[buffer_uint][3] = input_Ps[buffer_uint].Pz();
	}
	
	return computeMEs( input_Ps_i, input_IDs );
}

//////////////////////////////////////////////////////////////////////////
#endif	// end of (defined(MEKD_STANDALONE) && defined(MEKD_with_ROOT)) || !(defined(MEKD_STANDALONE))


#endif	// end of MEKD_MEKD_cpp
