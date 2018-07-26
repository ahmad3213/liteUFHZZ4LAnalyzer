#ifndef MEKD_TEST_DEBUG_GHOSTS_CPP
#define MEKD_TEST_DEBUG_GHOSTS_CPP

#include "MEKD_Test.h"


int MEKD_Test_Debug_Ghosts_Test0()
{
	/// TEST 0
	if( Show_Description ) cout << "\n -- STARTING TEST 0 -- \n";
	const unsigned int Nr_of_models = 4;
	double ME[Nr_of_models];
	MEKD test01(8.0, ""), test02(8.0, "");
	
	string model[Nr_of_models];
	model[0] = "ggSpin2Mh";	//8
	model[1] = "ggSpin2Pm";	//10
	model[2] = "ggSpin2Pm";	//8
	model[3] = "ggSpin2Mh";	//10
	
	if( Show_Description ) cout << "Testing ME ordering using int computeME( string, vector<double*>, vector<int>, double& );\n";
	
	if( (error_value=test01.computeME( model[0], Set_Of_Arrays, Set_Of_IDs, ME[0] )) != 0 ) cout << "ERROR CODE in ME for " << model[0] << "; Err: " << error_value << endl;
	if( (error_value=test01.computeME( model[1], Set_Of_Arrays, Set_Of_IDs, ME[1] )) != 0 ) cout << "ERROR CODE in ME for " << model[1] << "; Err: " << error_value << endl;
	if( (error_value=test02.computeME( model[2], Set_Of_Arrays, Set_Of_IDs, ME[2] )) != 0 ) cout << "ERROR CODE in ME for " << model[2] << "; Err: " << error_value << endl;
	if( (error_value=test02.computeME( model[3], Set_Of_Arrays, Set_Of_IDs, ME[3] )) != 0 ) cout << "ERROR CODE in ME for " << model[3] << "; Err: " << error_value << endl;

	if( Show_Basic_Data )
	{
		if( Show_Description ) cout << "Initial results, before any permutation.\n";
		cout.width( 20 );
		cout << std::left << "Model" << std::right << "ME value\n\n";

		for( unsigned int count=0; count < Nr_of_models; count++ )
		{
			cout.width( 20 );
			cout << std::left << model[count] << std::right << ME[count] << endl;
		}
	}
	
	
	return 0;
}



int MEKD_Test_Debug_Ghosts_Test1()
{
	/// TEST 1
	if( Show_Description ) cout << "\n -- STARTING TEST 1 -- \n";
	const unsigned int Nr_of_models = 21;
	double ME[Nr_of_models];
	double ME_shuffled[Nr_of_models];
	MEKD test1(8.0, "");
	
	string model[Nr_of_models];
	model[0] = "ZZ";
	model[1] = "ggSpin0Pm";
	model[2] = "ggSpin0Ph";
	model[3] = "ggSpin0M";
	model[4] = "qqZ4l_Background";
	model[5] = "qqZ4l_Signal";
	model[6] = "qqSpin1M";
	model[7] = "qqSpin1P";
	model[8] = "ggSpin2Pm";	//8
	model[9] = "ggSpin2Ph";	//9
	model[10] = "ggSpin2Mh";	//10
	model[11] = "ggSpin2Pb";	//11
	model[12] = "Spin0Pm";
	model[13] = "Spin0Ph";
	model[14] = "Spin0M";
	model[15] = "Spin1M";
	model[16] = "Spin1P";
	model[17] = "Spin2Pm";
	model[18] = "Spin2Ph";
	model[19] = "Spin2Mh";
	model[20] = "Spin2Pb";
	
	if( Show_Description ) cout << "Testing ME ordering using int computeME( string, vector<double*>, vector<int>, double& );\n";
	
	for( unsigned int count=0; count < Nr_of_models; count++ )
	{
		if( (error_value=test1.computeME( model[count], Set_Of_Arrays, Set_Of_IDs, ME[count] )) != 0 ) cout << "ERROR CODE in ME for " << model[count] << "; Err: " << error_value << endl;
	}

	if( Show_Basic_Data )
	{
		if( Show_Description ) cout << "Initial results, before any permutation.\n";
		cout.width( 20 );
		cout << std::left << "Model" << std::right << "ME value\n\n";

		for( unsigned int count=0; count < Nr_of_models; count++ )
		{
			cout.width( 20 );
			cout << std::left << model[count] << std::right << ME[count] << endl;
		}
	}
	
	
	unsigned int order[Nr_of_models];
	for( unsigned int count=0; count < Nr_of_models; count++ ) order[count] = count;
	std::sort( order, order+Nr_of_models );
	
	
	for( unsigned int perm=0; perm < shuffles_for_ghosts; perm++ )
	{
		std::random_shuffle( order, order+Nr_of_models );
		for( unsigned int count=0; count < Nr_of_models; count++ )
		{
			if( (error_value=test1.computeME( model[order[count]], Set_Of_Arrays, Set_Of_IDs, ME_shuffled[order[count]] )) != 0 ) cout << "ERROR CODE in ME for " << model[order[count]] << "; Err: " << error_value << endl;
			
			if( (ME[order[count]] - ME_shuffled[order[count]]) != 0 )
			{
				if( Show_Basic_Data )
				{
					cout << "Problem in moving order of calculation for " << model[order[count]] << " detected!\n";
					cout << "Previously calculated models in this shuffle: ";
					for( unsigned int count2=0; count2 < count; count2++ ) cout << model[order[count2]] << " ";
				}
				return 1;
			}
		}
		
		if( Show_Debug )
		{
			if( Show_Description ) cout << "Results after a permutation.\n";
			cout.width( 20 );
			cout << std::left << "Model" << std::right << "ME value\n\n";

			for( unsigned int count=0; count < Nr_of_models; count++ )
			{
				cout.width( 20 );
				cout << std::left << model[count] << std::right << ME_shuffled[count] << endl;
			}
		}
	}
	
	
	return 0;
}




int MEKD_Test_Debug_Ghosts_Test2()
{
	/// TEST 2
	if( Show_Description ) cout << "\n -- STARTING TEST 2 -- \n";
	const unsigned int Nr_of_models = 21;
	double ME[Nr_of_models];
	double ME_shuffled[Nr_of_models];
	MEKD test2(8.0, "");
	MEKD *test2_perm;
	
	string model[Nr_of_models];
	model[0] = "ZZ";
	model[1] = "ggSpin0Pm";
	model[2] = "ggSpin0Ph";
	model[3] = "ggSpin0M";
	model[4] = "qqZ4l_Background";
	model[5] = "qqZ4l_Signal";
	model[6] = "qqSpin1M";
	model[7] = "qqSpin1P";
	model[8] = "ggSpin2Pm";	//8
	model[9] = "ggSpin2Ph";	//9
	model[10] = "ggSpin2Mh";	//10
	model[11] = "ggSpin2Pb";	//11
	model[12] = "Spin0Pm";
	model[13] = "Spin0Ph";
	model[14] = "Spin0M";
	model[15] = "Spin1M";
	model[16] = "Spin1P";
	model[17] = "Spin2Pm";
	model[18] = "Spin2Ph";
	model[19] = "Spin2Mh";
	model[20] = "Spin2Pb";
	
	if( Show_Description ) cout << "Testing ME ordering using int MEKD.MEKD_MG_Calc.Run_MEKD_MG();\n";
	
	test2.MEKD_MG_Calc.p1 = p1;
	test2.MEKD_MG_Calc.p2 = p2;
	test2.MEKD_MG_Calc.p3 = p3;
	test2.MEKD_MG_Calc.p4 = p4;
	
	test2.MEKD_MG_Calc.id1 = id1;
	test2.MEKD_MG_Calc.id2 = id2;
	test2.MEKD_MG_Calc.id3 = id3;
	test2.MEKD_MG_Calc.id4 = id4;
	
	for( unsigned int count=0; count < Nr_of_models; count++ )
		test2.MEKD_MG_Calc.Test_Models.push_back( model[count] );
	
	if( (error_value=test2.MEKD_MG_Calc.Run_MEKD_MG()) != 0 ) cout << "ERROR CODE in MEKD.MEKD_MG_Calc.Run_MEKD_MG(); Err: " << error_value << endl;
	
	for( unsigned int count=0; count < Nr_of_models; count++ )
		ME[count] = test2.MEKD_MG_Calc.Signal_MEs[count];

	if( Show_Basic_Data )
	{
		if( Show_Description ) cout << "Initial results, before any permutation.\n";
		cout.width( 20 );
		cout << std::left << "Model" << std::right << "ME value\n\n";

		for( unsigned int count=0; count < Nr_of_models; count++ )
		{
			cout.width( 20 );
			cout << std::left << model[count] << std::right << ME[count] << endl;
		}
	}
	
	
	unsigned int order[Nr_of_models];
	for( unsigned int count=0; count < Nr_of_models; count++ ) order[count] = count;
	std::sort( order, order+Nr_of_models );
	
	
	for( unsigned int perm=0; perm < shuffles_for_ghosts; perm++ )
	{
		test2_perm = new MEKD(8.0, "");
		
		test2_perm->MEKD_MG_Calc.p1 = p1;
		test2_perm->MEKD_MG_Calc.p2 = p2;
		test2_perm->MEKD_MG_Calc.p3 = p3;
		test2_perm->MEKD_MG_Calc.p4 = p4;
		
		test2_perm->MEKD_MG_Calc.id1 = id1;
		test2_perm->MEKD_MG_Calc.id2 = id2;
		test2_perm->MEKD_MG_Calc.id3 = id3;
		test2_perm->MEKD_MG_Calc.id4 = id4;
		
		std::random_shuffle( order, order+Nr_of_models );
		
		for( unsigned int count=0; count < Nr_of_models; count++ )
			test2_perm->MEKD_MG_Calc.Test_Models.push_back( model[order[count]] );
		
		if( (error_value=test2_perm->MEKD_MG_Calc.Run_MEKD_MG()) != 0 ) cout << "ERROR CODE in MEKD.MEKD_MG_Calc.Run_MEKD_MG(); Err: " << error_value << endl;
		
		for( unsigned int count=0; count < Nr_of_models; count++ )
		ME_shuffled[order[count]] = test2_perm->MEKD_MG_Calc.Signal_MEs[count];
		
		for( unsigned int count=0; count < Nr_of_models; count++ )
		{
			if( (ME[order[count]] - ME_shuffled[order[count]]) != 0 )
			{
				if( Show_Basic_Data )
				{
					cout << "Problem in moving order of calculation for " << model[order[count]] << " detected!\n";
					cout << "Number of shuffle iterations: " << perm << endl;
					cout << "Previously calculated models in this shuffle: ";
					for( unsigned int count2=0; count2 < count; count2++ ) cout << model[order[count2]] << " ";
					cout << endl;
				}
				return 1;
			}
			
		}
		
		if( Show_Debug )
		{
			if( Show_Description ) cout << "Results after a permutation.\n";
			cout.width( 20 );
			cout << std::left << "Model" << std::right << "ME value\n\n";

			for( unsigned int count=0; count < Nr_of_models; count++ )
			{
				cout.width( 20 );
				cout << std::left << model[count] << std::right << ME_shuffled[count] << endl;
			}
		}
		
		delete test2_perm;
	}
	
	
	return 0;
}




#endif
