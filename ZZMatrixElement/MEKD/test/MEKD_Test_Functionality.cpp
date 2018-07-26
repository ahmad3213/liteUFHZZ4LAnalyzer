#ifndef MEKD_TEST_FUNCTIONALITY_CPP
#define MEKD_TEST_FUNCTIONALITY_CPP

#include "MEKD_Test.h"


bool Reference_is_calculated = false;

double KD1, ME11, ME21;




///////////////////////////
/// TEST Block "on the fly"
///////////////////////////


int MEKD_Test_Functionality_Test1()
{
	/// TEST 1
	if( Show_Description ) cout << "\n -- STARTING TEST 1 -- \n";
// 	double KD1, ME11, ME21;
	MEKD test1(8.0, "");
	
	
	if( Show_Description ) cout << "Testing int computeKD( string, string, double[], int, double[], int, double[], int, double[], int, double&, double&, double& );\n";
	
	test1.computeKD( (string) "ggSpin0Pm", (string) "ZZ", p1, id1, p2, id2, p3, id3, p4, id4, KD1, ME11, ME21 );
	
	if( Show_Basic_Data )
	{
		cout << "Output: ME1 ME2 KD\n";
		cout << ME11 << " " << ME21 << " " << KD1 << endl;
	}
	
	Reference_is_calculated = true;
	
	return 0;
}



int MEKD_Test_Functionality_Test2()
{
	if( !Reference_is_calculated ) MEKD_Test_Functionality_Test1();
	/// TEST 2
	if( Show_Description ) cout << "\n -- STARTING TEST 2 -- \n";
	double KD2, ME12, ME22;
	MEKD test2(8.0, "");
	
	
	if( Show_Description ) cout << "Testing int computeKD( string, string, vector<double*>, vector<int>, double&, double&, double& );\n";
	
	test2.computeKD( (string) "ggSpin0Pm", (string) "ZZ", Set_Of_Arrays, Set_Of_IDs, KD2, ME12, ME22 );
	
	if( Show_Basic_Data )
	{
		cout << "Output: ME1 ME2 KD\n";
		cout << ME12 << " " << ME22 << " " << KD2 << endl;
	}
	
	if( fabs(KD1-KD2) < Precision_of_interest && fabs(ME11-ME12) < Precision_of_interest && fabs(ME21-ME22) < Precision_of_interest ) return 0;
	return 1;
}



int MEKD_Test_Functionality_Test3()
{
	if( !Reference_is_calculated ) MEKD_Test_Functionality_Test1();
	/// TEST 3
	if( Show_Description ) cout << "\n -- STARTING TEST 3 -- \n";
	double KD3, ME13, ME23;
	MEKD test3(8.0, "");
	
	
	if( Show_Description ) cout << "Testing int computeKD( TString, TString, TLorentzVector, int, TLorentzVector, int, TLorentzVector, int, TLorentzVector, int, double&, double&, double& );\n";
	
	test3.computeKD( (TString) "ggSpin0Pm", (TString) "ZZ", Lp1, id1, Lp2, id2, Lp3, id3, Lp4, id4, KD3, ME13, ME23 );
	
	if( Show_Basic_Data )
	{
		cout << "Output: ME1 ME2 KD\n";
		cout << ME13 << " " << ME23 << " " << KD3 << endl;
	}
	
	if( fabs(KD1-KD3) < Precision_of_interest && fabs(ME11-ME13) < Precision_of_interest && fabs(ME21-ME23) < Precision_of_interest ) return 0;
	return 1;
}



int MEKD_Test_Functionality_Test4()
{
	if( !Reference_is_calculated ) MEKD_Test_Functionality_Test1();
	/// TEST 4
	if( Show_Description ) cout << "\n -- STARTING TEST 4 -- \n";
	double KD4, ME14, ME24;
	MEKD test4(8.0, "");
	
	vector<TLorentzVector> Set_Of_TLorentzVectors;
	
	Set_Of_TLorentzVectors.push_back( Lp1 );
	Set_Of_TLorentzVectors.push_back( Lp2 );
	Set_Of_TLorentzVectors.push_back( Lp3 );
	Set_Of_TLorentzVectors.push_back( Lp4 );
	
	
	if( Show_Description ) cout << "Testing int computeKD( TString, TString, vector<TLorentzVector>, vector<int>, double&, double&, double& );\n";
	
	test4.computeKD( (TString) "ggSpin0Pm", (TString) "ZZ", Set_Of_TLorentzVectors, Set_Of_IDs, KD4, ME14, ME24 );
	
	if( Show_Basic_Data )
	{
		cout << "Output: ME1 ME2 KD\n";
		cout << ME14 << " " << ME24 << " " << KD4 << endl;
	}
	
	if( fabs(KD1-KD4) < Precision_of_interest && fabs(ME11-ME14) < Precision_of_interest && fabs(ME21-ME24) < Precision_of_interest ) return 0;
	return 1;
}



//////////////////////////////
/// TEST Block "precalculated"
//////////////////////////////


int MEKD_Test_Functionality_Test5()
{
	if( !Reference_is_calculated ) MEKD_Test_Functionality_Test1();	
	/// TEST 1
	if( Show_Description ) cout << "\n -- STARTING TEST 5 -- \n";
	double KD1_p2, ME11_p2, ME21_p2;
	MEKD test1_p2(8.0, "");
	
	
	if( Show_Description ) cout << "Testing int computeMEs( double[], int, double[], int, double[], int, double[], int );\n";
	
	test1_p2.computeMEs( p1, id1, p2, id2, p3, id3, p4, id4 );
	test1_p2.computeKD( (string) "ggSpin0Pm", (string) "ZZ", KD1_p2, ME11_p2, ME21_p2 );
	
	if( Show_Basic_Data )
	{
		cout << "Output: ME1 ME2 KD\n";
		cout << ME11_p2 << " " << ME21_p2 << " " << KD1_p2 << endl;
	}
	
	if( fabs(KD1-KD1_p2) < Precision_of_interest && fabs(ME11-ME11_p2) < Precision_of_interest && fabs(ME21-ME21_p2) < Precision_of_interest ) return 0;
	return 1;
}
	
	
	
int MEKD_Test_Functionality_Test6()
{
	if( !Reference_is_calculated ) MEKD_Test_Functionality_Test1();	
	/// TEST 2
	if( Show_Description ) cout << "\n -- STARTING TEST 6 -- \n";
	double KD2_p2, ME12_p2, ME22_p2;
	MEKD test2_p2(8.0, "");
	
	
	if( Show_Description ) cout << "Testing int computeMEs( vector<double*>, vector<int> );\n";
	
	test2_p2.computeMEs( Set_Of_Arrays, Set_Of_IDs );
	test2_p2.computeKD( (string) "ggSpin0Pm", (string) "ZZ", KD2_p2, ME12_p2, ME22_p2 );
	
	if( Show_Basic_Data )
	{
		cout << "Output: ME1 ME2 KD\n";
		cout << ME12_p2 << " " << ME22_p2 << " " << KD2_p2 << endl;
	}
	
	if( fabs(KD1-KD2_p2) < Precision_of_interest && fabs(ME11-ME12_p2) < Precision_of_interest && fabs(ME21-ME22_p2) < Precision_of_interest ) return 0;
	return 1;
}
	
	
	
int MEKD_Test_Functionality_Test7()
{
	if( !Reference_is_calculated ) MEKD_Test_Functionality_Test1();	
	/// TEST 3
	if( Show_Description ) cout << "\n -- STARTING TEST 7 -- \n";
	double KD3_p2, ME13_p2, ME23_p2;
	MEKD test3_p2(8.0, "");
	
	
	if( Show_Description ) cout << "Testing int computeMEs( TLorentzVector, int, TLorentzVector, int, TLorentzVector, int, TLorentzVector, int );\n";
	
	test3_p2.computeMEs( Lp1, id1, Lp2, id2, Lp3, id3, Lp4, id4 );
	test3_p2.computeKD( (string) "ggSpin0Pm", (string) "ZZ", KD3_p2, ME13_p2, ME23_p2 );
	
	if( Show_Basic_Data )
	{
		cout << "Output: ME1 ME2 KD\n";
		cout << ME13_p2 << " " << ME23_p2 << " " << KD3_p2 << endl;
	}
	
	if( fabs(KD1-KD3_p2) < Precision_of_interest && fabs(ME11-ME13_p2) < Precision_of_interest && fabs(ME21-ME23_p2) < Precision_of_interest ) return 0;
	return 1;
}
	
	
	
int MEKD_Test_Functionality_Test8()
{
	if( !Reference_is_calculated ) MEKD_Test_Functionality_Test1();	
	/// TEST 4
	if( Show_Description ) cout << "\n -- STARTING TEST 8 -- \n";
	double KD4_p2, ME14_p2, ME24_p2;
	MEKD test4_p2(8.0, "");
	
	vector<TLorentzVector> Set_Of_TLorentzVectors;
	
	Set_Of_TLorentzVectors.push_back( Lp1 );
	Set_Of_TLorentzVectors.push_back( Lp2 );
	Set_Of_TLorentzVectors.push_back( Lp3 );
	Set_Of_TLorentzVectors.push_back( Lp4 );
	
	
	if( Show_Description ) cout << "Testing int computeMEs( vector<TLorentzVector>, vector<int> );\n";
	
	test4_p2.computeMEs( Set_Of_TLorentzVectors, Set_Of_IDs );
	test4_p2.computeKD( (string) "ggSpin0Pm", (string) "ZZ", KD4_p2, ME14_p2, ME24_p2 );
	
	if( Show_Basic_Data )
	{
		cout << "Output: ME1 ME2 KD\n";
		cout << ME14_p2 << " " << ME24_p2 << " " << KD4_p2 << endl;
	}
	
	if( fabs(KD1-KD4_p2) < Precision_of_interest && fabs(ME11-ME14_p2) < Precision_of_interest && fabs(ME21-ME24_p2) < Precision_of_interest ) return 0;
	return 1;
}
	
	
#endif