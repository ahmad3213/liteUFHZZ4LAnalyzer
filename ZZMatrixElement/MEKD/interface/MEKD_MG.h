#ifndef MEKD_MG_h
#define MEKD_MG_h

#include <string>
#include <vector>

#include <sstream>
#include <complex>

#include "../src/MadGraphSrc/read_slha.h"

using namespace std;



class MEKD_MG
{
public:
	/// Flags
	bool Boost_To_CM;	// for a boosted data
	bool Debug_Mode;	// Enable debugging mode
	bool Fix_Spin0_Production;	// use the SM Higgs production mechanism
	bool Fix_Spin1_Production;	// use the a hybrid production mechanism
// 	bool Force_g3_running;	// unused. At some point was included for alpha_QCD
	bool Overwrite_e_and_mu_masses;	// switch for manual m_e, m_mu masses
	bool Use_Higgs_width;	//	if false, width is fixed to =1
	bool Use_mh_eq_m4l;	// Set mh to m4l for every event
	bool Use_mZ4l_eq_m4l;	// Set m_Z to m4l for Z4l events
	bool Use_PDF_w_pT0;	// Use PDFs in the pT=0 frame. If true, Boost_To_CM is ignored
	bool Vary_resonance_width;	// Allow width to be varied with mass
	bool Vary_signal_couplings;	// Allow couplings to change with mass
	bool Warning_Mode;	// Print warnings
	
	/// General parameters
	double ContributionCoeff_d;	//42	/// the value has no effect if PDF is used but the variable is always used
	double ContributionCoeff_u;	//217
	double ContributionCoeff_s;	//5
	double ContributionCoeff_c;	//3
// 	double GG;	// Assign QCD coupling, force g3 running if needed
	double Sqrt_s;	//Max energy, collision energy
	
	/// State mixing
	complex<double> *Mixing_Coefficients_Spin0, *Mixing_Coefficients_Spin1, *Mixing_Coefficients_Spin2;
	
	/// Physical parameters
	double Electron_mass;	//0.0005109989, for enabled overwriting
	double Higgs_mass;	// Works only if Use_mh_eq_m4l=false
	double Higgs_width;	// Practically not used, for future implementations
	double Muon_mass;	//0.10565837, for enabled overwriting
	double Proton_mass;	// Always used if needed
	
	/// Final-state lepton/photon information
	double *p1, *p2, *p3, *p4, *p5;
	double id1, id2, id3, id4, id5;
	
	/// String flags and file locations
	string Final_state;	// Final state, for the moment: 4e, 4mu, 2e2mu
	string Resonance_decay_mode;	// default: ZZ. Alternatives: 2l, 2l_s
	string Test_Model;	// Models: ZZ, DY, Custom, CPevenScalar, ggSpin0Pm, ggSpin0M, ggSpin0Ph, qqSpin1P, qqSpin1M, ggSpin2Pm, ggSpin2Ph, ggSpin2Mh, ggSpin2Pb, qqSpin2Pm, qqSpin2Ph, qqSpin2Mh, qqSpin2Pb, Spin0Pm, Spin0M, Spin0Ph, Spin1P, Spin1M, Spin2Pm, Spin2Ph, Spin2Mh, Spin2Pb, qqZ4l_Signal, qqZ4l_Background
	vector<string> Test_Models;	// same names as for the Test_Model
	string Parameter_file;	// Location where a parameter card is stored
	string PDF_file;	// PDF/PDT table file
	
	/// Calculation results
	double Mass_4l;	//is filled after running RUN_XXXX(...). Invariant mass of the final-state system
	double Background_ME;	//may not be used if running RUN_MEKD_MG( string ) is chosen
	double Signal_ME;	//is filled after running RUN_XXXX(...)
	vector<double> Signal_MEs;	//is filled if Test_Models are set after running RUN_XXXX(...)
	double KD;	//is not filled with RUN_MEKD_MG( string )
	
	/// Parameter container. For experts only
	SLHAReader_MEKD Set_Of_Model_Parameters;
	
	/// Functions
	void Set_Default_MEKD_MG_Parameters();
	
	int Reload_Parameters();	// reloads parameter set and updates PDF file reader
	int Run_MEKD_MG();	// main routine to evaluate matrix elements; updates "Calculation results"
	int Run_MEKD_MG(string Input_Model);	// Calculates a ME ONLY for a chosen model; ignores automatic background calculation. Updates Signal_ME
	
	/// Constructors, destructors
	MEKD_MG();
	~MEKD_MG();
	
private:
	bool Parameters_Are_Loaded, buffer_bool, Predefined_Model;
	
	int error_value;
	unsigned int counter;
	
	double *buffer, buffer_p[4], buffer_Custom, ml1, ml2, ml3, ml4, PDFx1, PDFx2, LmbdGG_calculated;
	double *pl1_internal, *pl2_internal, *pl3_internal, *pl4_internal, *pA1_internal;
	
	complex<double> *buffer_complex, *Mixing_Coefficients_Spin0_internal, *Mixing_Coefficients_Spin1_internal, *Mixing_Coefficients_Spin2_internal;
	
	string *Test_Model_buffer;
	
	// Parameters
	double v_expectation;	// Vacuum expectation value
	double hZZ_coupling;
	double params_m_d, params_m_u, params_m_s, params_m_c, params_m_e, params_m_mu, params_m_Z;
	complex<double> params_rhou01, params_rhou02, params_rhoc01, params_rhoc02,
		params_rhod01, params_rhod02, params_rhos01, params_rhos02,
		params_rhob01, params_rhob02;
	complex<double> params_rhou11, params_rhou12, params_rhou13, params_rhou14,
		params_rhoc11, params_rhoc12, params_rhoc13, params_rhoc14,
		params_rhod11, params_rhod12, params_rhod13, params_rhod14,
		params_rhos11, params_rhos12, params_rhos13, params_rhos14,
		params_rhob11, params_rhob12, params_rhob13, params_rhob14;
	complex<double> params_rhou21, params_rhou22, params_rhou23, params_rhou24,
		params_rhoc21, params_rhoc22, params_rhoc23, params_rhoc24,
		params_rhod21, params_rhod22, params_rhod23, params_rhod24,
		params_rhos21, params_rhos22, params_rhos23, params_rhos24,
		params_rhob21, params_rhob22, params_rhob23, params_rhob24;
	
	
	string buffer_string;
	
	vector<double> id_set;
	vector<double*> p_set;
	
	/// Internal functions ///
	int Load_Parameters();
	
	int Arrange_Internal_pls();
	
	/// Sets up particular choices. Tier 3
	int Run_MEKD_MG_ME_Configurator_BKG_ZZ(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Custom();
	int Run_MEKD_MG_ME_Configurator_CPPProcess(string initial_state);	// RAW MG5_aMC ME
	int Run_MEKD_MG_ME_Configurator_Spin0(string initial_state);	// A general mixed spin-0 state
	int Run_MEKD_MG_ME_Configurator_Spin1(string initial_state);	// A general mixed spin-1 state
	int Run_MEKD_MG_ME_Configurator_Spin2(string initial_state);	// A general mixed spin-2 state
	int Run_MEKD_MG_ME_Configurator_Spin0Pm(string initial_state);	// SM Higgs
	int Run_MEKD_MG_ME_Configurator_Spin0M(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin0Ph(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin1P(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin1M(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin2Pm(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin2Ph(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin2Mh(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin2Pb(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin2Ph2(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin2Ph3(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin2Ph6(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin2Ph7(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin2Mh9(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin2Mh10(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin0Pm_Spin0M(string initial_state);	// A mixed state of two contributions
	int Run_MEKD_MG_ME_Configurator_Spin0Pm_Spin0Ph(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Spin0M_Spin0Ph(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Z4l_BKG(string initial_state);
	int Run_MEKD_MG_ME_Configurator_Z4l_SIG(string initial_state);
	
	/// Dispatches MEs that have correct parameters. Tier 2
	int Run_MEKD_MG_ME_Dispatcher_CPPProcess(string initial_state);	// RAW MG5_aMC ME
	int Run_MEKD_MG_ME_Dispatcher_BKG_ZZ(string initial_state);
	int Run_MEKD_MG_ME_Dispatcher_Z4l_BKG(string initial_state);
	int Run_MEKD_MG_ME_Dispatcher_Z4l_SIG(string initial_state);
	int Run_MEKD_MG_ME_Dispatcher_SIG_Spin0(string initial_state);
	int Run_MEKD_MG_ME_Dispatcher_SIG_Spin1(string initial_state);
	int Run_MEKD_MG_ME_Dispatcher_SIG_Spin2(string initial_state);
	
	/// Evaluators. Blind-calculation functions. Handles MEs from Dispatchers. Tier 1
	template<class Generic_MEKD_MG_ME>
	int Run_MEKD_MG_MEs_Evaluator_Initial_State_NO(bool photon, Generic_MEKD_MG_ME &Generic_ME);
	
	template<class Generic_MEKD_MG_ME>
	int Run_MEKD_MG_MEs_Evaluator_Initial_State_gg(bool photon, Generic_MEKD_MG_ME &Generic_ME);
	
	template<class Generic_MEKD_MG_ME_s, class Generic_MEKD_MG_ME_c>
	int Run_MEKD_MG_MEs_Evaluator_Initial_State_qqbar(bool photon, Generic_MEKD_MG_ME_s &Generic_ME_s, Generic_MEKD_MG_ME_c &Generic_ME_c);
};


#endif