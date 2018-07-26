#ifndef MEKD_TEST_CHECK_MODELS_CPP
#define MEKD_TEST_CHECK_MODELS_CPP

#include "MEKD_Test.h"
#include <iomanip>


//////////////////////////////////
/// TEST Block "Models checking"
//////////////////////////////////



/// Check 2l (plus photon) final-state models.
int MEKD_Test_Check_Models_Test1()
{
	/// TEST 1
	if( Show_Description ) cout << "\n -- STARTING TEST 1 -- \n";
	MEKD test1(8.0, "");
	
	
	int id2l1, id2l2, id2l3;
	
	double p2l1[4]={0}, p2l2[4]={0}, p2lA1[4]={0}, ME_value;
	
	vector<string> models_to_test;
	vector<double*> Set_Of_Arrays_2l, Set_Of_Arrays_2lpA;
	vector<int> Set_Of_IDs_2l, Set_Of_IDs_2lpA;
	
	
	if( Show_Description )
	{
		cout << "\n ------------------------------------------------ \n";
		cout << " ------------------------------------------------ \n";
		cout << " -- STARTING MODEL TESTs, using 2mu (+A) event -- \n";
		cout << " ------------------------------------------------ \n";
		cout << " ------------------------------------------------ \n";
	}
	
	id2l1=13;
	id2l2=-13;
	id2l3=22;
	
	p2l1[0] = 62;
	p2l1[1] = 50;
	p2l1[2] = -17;
	p2l1[3] = -32.4805915;
	
	p2l2[0] = 62;
	p2l2[1] = -50;
	p2l2[2] = 17;
	p2l2[3] = 32.4805915;
	
	// setting collinear to the first lepton
	p2lA1[0] = 4;
	p2lA1[1] = 3.23;
	p2lA1[2] = -1.1;
	p2lA1[3] = -2.087367;
	
	
	Set_Of_Arrays_2l.push_back( p2l1 );
	Set_Of_Arrays_2l.push_back( p2l2 );
	
	Set_Of_IDs_2l.push_back( id2l1 );
	Set_Of_IDs_2l.push_back( id2l2 );
	
	Set_Of_Arrays_2lpA.push_back( p2l1 );
	Set_Of_Arrays_2lpA.push_back( p2l2 );
	Set_Of_Arrays_2lpA.push_back( p2lA1 );
	
	Set_Of_IDs_2lpA.push_back( id2l1 );
	Set_Of_IDs_2lpA.push_back( id2l2 );
	Set_Of_IDs_2lpA.push_back( id2l3 );
	
	
	models_to_test.push_back( "qqDY" );
	models_to_test.push_back( "ggSpin0Pm" );
	models_to_test.push_back( "ggSpin0M" );
	models_to_test.push_back( "qqSpin1M" );
	models_to_test.push_back( "qqSpin1P" );
	models_to_test.push_back( "ggSpin2Pm" );
	models_to_test.push_back( "qqSpin2Pm" );
	
	models_to_test.push_back( "DY" );
	models_to_test.push_back( "Spin0Pm" );
	models_to_test.push_back( "Spin0M" );
	models_to_test.push_back( "Spin1M" );
	models_to_test.push_back( "Spin1P" );
	models_to_test.push_back( "Spin2Pm" );
	
	
	
	
	/// RUNING CALCULATIONS BELOW
	cout << setw(30) << left << "Description " << setw(15) << left << "Model name" << ": " << "ME value" << endl;
	cout << "-------------------------------------------------------------\n";
	for( unsigned int count_model=0; count_model<models_to_test.size(); count_model++ )
	{
		if( (error_value=test1.computeME( models_to_test[count_model], Set_Of_Arrays_2l, Set_Of_IDs_2l, ME_value )) != 0 ) cout << "ERROR CODE in ME for " << models_to_test[count_model]  << ": " << error_value << endl;
		else cout << setw(30) << left << "Result for " << setw(15) << left << models_to_test[count_model] << ": " << ME_value << endl;
		
		if( (error_value=test1.computeME( models_to_test[count_model], Set_Of_Arrays_2lpA, Set_Of_IDs_2lpA, ME_value )) != 0 ) cout << "ERROR CODE in ME for " << models_to_test[count_model]  << ": " << error_value << endl;
		else cout << setw(30) << left << "Result for (with a photon)" << setw(15) << left << models_to_test[count_model] << ": " << ME_value << endl;
	}
	
	
	
	return 1;
}



/// Check 2l (plus photon) final-state models.
int MEKD_Test_Check_Models_Test2()
{
	/// TEST 2
	if( Show_Description ) cout << "\n -- STARTING TEST 2 -- \n";
	
	if( Show_Description )
	{
		cout << "\n -------------------------------------------------- \n";
		cout << " -------------------------------------------------- \n";
		cout << " -- STARTING MODEL TESTs, using 2e2mu (+A) event -- \n";
		cout << " -------------------------------------------------- \n";
		cout << " -------------------------------------------------- \n";
	}
	
	MEKD test2(8.0, "");
	
	
	double ME_value, pA1[4]={0};
	
	vector<string> models_to_test;
	vector<double*> Set_Of_Arrays_pA;
	vector<int> Set_Of_IDs_pA;
	
	
	// setting collinear to the first lepton
	pA1[0] = 0.4;
	pA1[1] = 0.23;
	pA1[2] = 0.221;
	pA1[3] = -0.24136902;
	
	
	Set_Of_Arrays_pA.push_back( Set_Of_Arrays[0] );
	Set_Of_Arrays_pA.push_back( Set_Of_Arrays[1] );
	Set_Of_Arrays_pA.push_back( Set_Of_Arrays[2] );
	Set_Of_Arrays_pA.push_back( Set_Of_Arrays[3] );
	Set_Of_Arrays_pA.push_back( pA1 );
	
	Set_Of_IDs_pA.push_back( Set_Of_IDs[0] );
	Set_Of_IDs_pA.push_back( Set_Of_IDs[1] );
	Set_Of_IDs_pA.push_back( Set_Of_IDs[2] );
	Set_Of_IDs_pA.push_back( Set_Of_IDs[3] );
	Set_Of_IDs_pA.push_back( 22 );
	
	
	models_to_test.push_back( "qqZZ" );
	models_to_test.push_back( "ggSpin0Pm" );
	models_to_test.push_back( "ggSpin0Ph" );
	models_to_test.push_back( "ggSpin0M" );
// 	models_to_test.push_back( "ggSpin0Pm_2f" );
// 	models_to_test.push_back( "ggSpin0M_2f" );
	models_to_test.push_back( "qqSpin1M" );
	models_to_test.push_back( "qqSpin1P" );
	models_to_test.push_back( "ggSpin2Pm" );
	models_to_test.push_back( "qqSpin2Pm" );
	models_to_test.push_back( "ggSpin2Pm_2f" );
	models_to_test.push_back( "qqSpin2Pm_2f" );
	
	models_to_test.push_back( "Spin0Pm" );
	models_to_test.push_back( "Spin0Ph" );
	models_to_test.push_back( "Spin0M" );
// 	models_to_test.push_back( "Spin0Pm_2f" );
// 	models_to_test.push_back( "Spin0M_2f" );
	models_to_test.push_back( "Spin1M" );
	models_to_test.push_back( "Spin1P" );
	models_to_test.push_back( "Spin2Pm" );
	models_to_test.push_back( "Spin2Pm_2f" );
	
	
	
	
	/// RUNING CALCULATIONS BELOW
	cout << setw(30) << left << "Description " << setw(15) << left << "Model name" << ": " << "ME value" << endl;
	cout << "-------------------------------------------------------------\n";
	for( unsigned int count_model=0; count_model<models_to_test.size(); count_model++ )
	{
		if( (error_value=test2.computeME( models_to_test[count_model], Set_Of_Arrays, Set_Of_IDs, ME_value )) != 0 ) cout << "ERROR CODE in ME for " << models_to_test[count_model]  << ": " << error_value << endl;
		else cout << setw(30) << left << "Result for " << setw(15) << left << models_to_test[count_model] << ": " << ME_value << endl;
		
		if( (error_value=test2.computeME( models_to_test[count_model], Set_Of_Arrays_pA, Set_Of_IDs_pA, ME_value )) != 0 ) cout << "ERROR CODE in ME for " << models_to_test[count_model]  << ": " << error_value << endl;
		else cout << setw(30) << left << "Result for (with a photon)" << setw(15) << left << models_to_test[count_model] << ": " << ME_value << endl;
	}
	
	
	
	if( Show_Description )
	{
		cout << "\n ------------------------------------------------ \n";
		cout << " ------------------------------------------------ \n";
		cout << " -- STARTING MODEL TESTs, using 4mu (+A) event -- \n";
		cout << " ------------------------------------------------ \n";
		cout << " ------------------------------------------------ \n";
	}
	
	
	Set_Of_IDs[0] = 13;
	Set_Of_IDs[1] = -13;
	
	Set_Of_IDs_pA[0] = 13;
	Set_Of_IDs_pA[1] = -13;
	
	
	/// RUNING CALCULATIONS BELOW
	cout << setw(30) << left << "Description " << setw(15) << left << "Model name" << ": " << "ME value" << endl;
	cout << "-------------------------------------------------------------\n";
	for( unsigned int count_model=0; count_model<models_to_test.size(); count_model++ )
	{
		if( (error_value=test2.computeME( models_to_test[count_model], Set_Of_Arrays, Set_Of_IDs, ME_value )) != 0 ) cout << "ERROR CODE in ME for " << models_to_test[count_model]  << ": " << error_value << endl;
		else cout << setw(30) << left << "Result for " << setw(15) << left << models_to_test[count_model] << ": " << ME_value << endl;
		
		if( (error_value=test2.computeME( models_to_test[count_model], Set_Of_Arrays_pA, Set_Of_IDs_pA, ME_value )) != 0 ) cout << "ERROR CODE in ME for " << models_to_test[count_model]  << ": " << error_value << endl;
		else cout << setw(30) << left << "Result for (with a photon)" << setw(15) << left << models_to_test[count_model] << ": " << ME_value << endl;
	}
	
	Set_Of_IDs[0] = 11;
	Set_Of_IDs[1] = -11;
	
	
	
	return 1;
}



#endif