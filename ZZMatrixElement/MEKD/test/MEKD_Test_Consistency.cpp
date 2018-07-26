#ifndef MEKD_TEST_CONSISTENCY_CPP
#define MEKD_TEST_CONSISTENCY_CPP

#include "MEKD_Test.h"


//////////////////////////////////
/// TEST Block "Consistency Check"
//////////////////////////////////



/// Consistency between qq -> Z -> 4l plus background channels versus the full qq -> ZZ -> 4l. Z -> 4l is a subset of ZZ -> 4l.
int MEKD_Test_Consistency_Test1()
{
	/// TEST 1
	if( Show_Description ) cout << "\n -- STARTING TEST 1 -- \n";
	double ME11_p3, ME21_p3, ME31_p3;
	MEKD test1_p3(8.0, "");
	
	test1_p3.MEKD_MG_Calc.Use_mZ4l_eq_m4l = false;	// use false
// 	test1_p3.MEKD_MG_Calc.Debug_Mode = true;
	
	if( Show_Description ) cout << "Testing Z4l using int computeME( string, vector<double*>, vector<int>, double& );\n";
	
	if( (error_value=test1_p3.computeME( (string) "ZZ", Set_Of_Arrays, Set_Of_IDs, ME11_p3 )) != 0 ) cout << "ERROR CODE in ME1: " << error_value << endl;
	if( (error_value=test1_p3.computeME( (string) "qqZ4l_Background", Set_Of_Arrays, Set_Of_IDs, ME21_p3 )) !=0 ) cout << "ERROR CODE in ME2: " << error_value << endl;
	if( (error_value=test1_p3.computeME( (string) "qqZ4l_Signal", Set_Of_Arrays, Set_Of_IDs, ME31_p3 )) !=0 ) cout << "ERROR CODE in ME3: " << error_value << endl;
	
	if( Show_Basic_Data )
	{
		cout << "Output: ME1 ME2 ME3\n";
		cout << ME11_p3 << " " << ME21_p3 << " " << ME31_p3 << endl;
	}
	
	if( fabs(ME11_p3-ME21_p3-ME31_p3)/fabs(ME11_p3) < Match_between_ZZ_Full_and_Z4l_rel ) return 0;
	return 1;
}
	
	
#endif