#ifndef MEKD_TEST_GEN_CPP
#define MEKD_TEST_GEN_CPP

#include "MEKD_Test.h"




void MEKD_Test_Gen_Test1()
{
	//////////////////////////////////
	/// TEST Block "Correlation Check"
	//////////////////////////////////
	
	if( Show_Description )
	{
		cout << "\n -------------------------------------------------------------------- \n";
		cout << " -- Testing ME correlations for consistency. Generating 4e events. -- \n";
		cout << " -------------------------------------------------------------------- \n";
	}
	
	if( Show_Description ) cout << "Testing Z4l correlations (summed up signal and background) to full qq -> ZZ using int computeME( string, vector<double*>, vector<int>, double& );\n";
	
	double ME1_p4[correlation_points], ME2_p4[correlation_points], ME3_p4[correlation_points];
	double val, invariant_mass, p_lepton1[3], p_lepton2[3], p_lepton3[3], p_lepton4[3];
	TRandom3 r1;
	
	/// TEST 1
	if( Show_Description ) cout << "\n -- STARTING TEST 1. Invariant mass 80--100 GeV. m4l != mZ -- \n";
	MEKD test1_p4(8.0, "");
	double mean_ME11_p4=0, mean_ME21_p4=0, mean_ME31_p4=0, corr_nom1=0, corr_denom1_1=0, corr_denom2_1=0;
	
	test1_p4.MEKD_MG_Calc.Use_mZ4l_eq_m4l = false;
	
	
	Set_Of_IDs[0] = 11;
	Set_Of_IDs[1] = -11;
	Set_Of_IDs[2] = 11;
	Set_Of_IDs[3] = -11;
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		invariant_mass = (100-80)*(r1.Uniform())+80;
		
		for( unsigned int count2=0; count2<2; count2++)
		{
			p_lepton1[count2] = 2*r1.Uniform()-1;
			p_lepton2[count2] = 2*r1.Uniform()-1;
			p_lepton3[count2] = 2*r1.Uniform()-1;
		}
		
		p_lepton4[0] = -(p_lepton1[0]+p_lepton2[0]+p_lepton3[0]);
		p_lepton4[1] = -(p_lepton1[1]+p_lepton2[1]+p_lepton3[1]);
		p_lepton4[2] = -(p_lepton1[2]+p_lepton2[2]+p_lepton3[2]);
		
		val = sqrt( p_lepton1[0]*p_lepton1[0]+p_lepton1[1]*p_lepton1[1]+p_lepton1[2]*p_lepton1[2]+m_e*m_e );
		val += sqrt( p_lepton2[0]*p_lepton2[0]+p_lepton2[1]*p_lepton2[1]+p_lepton2[2]*p_lepton2[2]+m_e*m_e );
		val += sqrt( p_lepton3[0]*p_lepton3[0]+p_lepton3[1]*p_lepton3[1]+p_lepton3[2]*p_lepton3[2]+m_e*m_e );
		val += sqrt( p_lepton4[0]*p_lepton4[0]+p_lepton4[1]*p_lepton4[1]+p_lepton4[2]*p_lepton4[2]+m_e*m_e );
		
		val = invariant_mass/val;
		
		for( unsigned int count2=1; count2<3; count2++ )
		{
			Set_Of_Arrays[0][count2] = p_lepton1[count2-1]*val;
			Set_Of_Arrays[1][count2] = p_lepton2[count2-1]*val;
			Set_Of_Arrays[2][count2] = p_lepton3[count2-1]*val;
			Set_Of_Arrays[3][count2] = p_lepton4[count2-1]*val;
		}
		Set_Of_Arrays[0][0] = sqrt( p_lepton1[0]*p_lepton1[0] +p_lepton1[1]*p_lepton1[1] +p_lepton1[2]*p_lepton1[2] + m_e*m_e );
		Set_Of_Arrays[1][0] = sqrt( p_lepton2[0]*p_lepton2[0] +p_lepton2[1]*p_lepton2[1] +p_lepton2[2]*p_lepton2[2] + m_e*m_e );
		Set_Of_Arrays[2][0] = sqrt( p_lepton3[0]*p_lepton3[0] +p_lepton3[1]*p_lepton3[1] +p_lepton3[2]*p_lepton3[2] + m_e*m_e );
		Set_Of_Arrays[3][0] = sqrt( p_lepton4[0]*p_lepton4[0] +p_lepton4[1]*p_lepton4[1] +p_lepton4[2]*p_lepton4[2] + m_e*m_e );
		
		if( (error_value=test1_p4.computeME( (string) "ZZ", Set_Of_Arrays, Set_Of_IDs, ME1_p4[count] )) != 0 ) cout << "ERROR CODE in ME1: " << error_value << endl;
		if( (error_value=test1_p4.computeME( (string) "qqZ4l_Background", Set_Of_Arrays, Set_Of_IDs, ME2_p4[count] )) !=0 ) cout << "ERROR CODE in ME2: " << error_value << endl;
		if( (error_value=test1_p4.computeME( (string) "qqZ4l_Signal", Set_Of_Arrays, Set_Of_IDs, ME3_p4[count] )) !=0 ) cout << "ERROR CODE in ME3: " << error_value << endl;
		
// 		cout << ME1_p4[count]) << " " << ME2_p4[count]+ME3_p4[count] << endl;
	}
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		mean_ME11_p4 += ME1_p4[count];
		mean_ME21_p4 += ME2_p4[count];
		mean_ME31_p4 += ME3_p4[count];
	}
	mean_ME11_p4 /= correlation_points;
	mean_ME21_p4 /= correlation_points;
	mean_ME31_p4 /= correlation_points;
	
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		corr_nom1 += (ME2_p4[count]+ME3_p4[count] - mean_ME21_p4-mean_ME31_p4)*(ME1_p4[count] - mean_ME11_p4);
		corr_denom1_1 += (ME2_p4[count]+ME3_p4[count] - mean_ME21_p4-mean_ME31_p4)*(ME2_p4[count]+ME3_p4[count] - mean_ME21_p4-mean_ME31_p4);
		corr_denom2_1 += (ME1_p4[count] - mean_ME11_p4)*(ME1_p4[count] - mean_ME11_p4);
	}
	
	if( Show_Basic_Data )
	{
		cout << "Output: correlation_coefficient\n";
		cout << (corr_nom1/sqrt(corr_denom1_1)/sqrt(corr_denom2_1)) << endl;
		cout << "Rel. Difference between means ME1 and ME2+ME3 wrt ME1: " << (mean_ME11_p4-mean_ME21_p4-mean_ME31_p4)/mean_ME11_p4 << endl;
	}
	
	
	
	/// TEST 2
	if( Show_Description ) cout << "\n -- STARTING TEST 2. Invariant mass 110--130 GeV. m4l != mZ -- \n";
	MEKD test2_p4(8.0, "");
	double mean_ME12_p4=0, mean_ME22_p4=0, mean_ME32_p4=0, corr_nom2=0, corr_denom1_2=0, corr_denom2_2=0;
	
	test2_p4.MEKD_MG_Calc.Use_mZ4l_eq_m4l = false;
	
	
	Set_Of_IDs[0] = 11;
	Set_Of_IDs[1] = -11;
	Set_Of_IDs[2] = 11;
	Set_Of_IDs[3] = -11;
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		invariant_mass = (130-110)*(r1.Uniform())+110;
		
		for( unsigned int count2=0; count2<2; count2++)
		{
			p_lepton1[count2] = 2*r1.Uniform()-1;
			p_lepton2[count2] = 2*r1.Uniform()-1;
			p_lepton3[count2] = 2*r1.Uniform()-1;
		}
		
		p_lepton4[0] = -(p_lepton1[0]+p_lepton2[0]+p_lepton3[0]);
		p_lepton4[1] = -(p_lepton1[1]+p_lepton2[1]+p_lepton3[1]);
		p_lepton4[2] = -(p_lepton1[2]+p_lepton2[2]+p_lepton3[2]);
		
		val = sqrt( p_lepton1[0]*p_lepton1[0]+p_lepton1[1]*p_lepton1[1]+p_lepton1[2]*p_lepton1[2]+m_e*m_e );
		val += sqrt( p_lepton2[0]*p_lepton2[0]+p_lepton2[1]*p_lepton2[1]+p_lepton2[2]*p_lepton2[2]+m_e*m_e );
		val += sqrt( p_lepton3[0]*p_lepton3[0]+p_lepton3[1]*p_lepton3[1]+p_lepton3[2]*p_lepton3[2]+m_e*m_e );
		val += sqrt( p_lepton4[0]*p_lepton4[0]+p_lepton4[1]*p_lepton4[1]+p_lepton4[2]*p_lepton4[2]+m_e*m_e );
		
		val = invariant_mass/val;
		
		for( unsigned int count2=1; count2<3; count2++ )
		{
			Set_Of_Arrays[0][count2] = p_lepton1[count2-1]*val;
			Set_Of_Arrays[1][count2] = p_lepton2[count2-1]*val;
			Set_Of_Arrays[2][count2] = p_lepton3[count2-1]*val;
			Set_Of_Arrays[3][count2] = p_lepton4[count2-1]*val;
		}
		Set_Of_Arrays[0][0] = sqrt( p_lepton1[0]*p_lepton1[0] +p_lepton1[1]*p_lepton1[1] +p_lepton1[2]*p_lepton1[2] + m_e*m_e );
		Set_Of_Arrays[1][0] = sqrt( p_lepton2[0]*p_lepton2[0] +p_lepton2[1]*p_lepton2[1] +p_lepton2[2]*p_lepton2[2] + m_e*m_e );
		Set_Of_Arrays[2][0] = sqrt( p_lepton3[0]*p_lepton3[0] +p_lepton3[1]*p_lepton3[1] +p_lepton3[2]*p_lepton3[2] + m_e*m_e );
		Set_Of_Arrays[3][0] = sqrt( p_lepton4[0]*p_lepton4[0] +p_lepton4[1]*p_lepton4[1] +p_lepton4[2]*p_lepton4[2] + m_e*m_e );
		
		if( (error_value=test2_p4.computeME( (string) "ZZ", Set_Of_Arrays, Set_Of_IDs, ME1_p4[count] )) != 0 ) cout << "ERROR CODE in ME1: " << error_value << endl;
		if( (error_value=test2_p4.computeME( (string) "qqZ4l_Background", Set_Of_Arrays, Set_Of_IDs, ME2_p4[count] )) !=0 ) cout << "ERROR CODE in ME2: " << error_value << endl;
		if( (error_value=test2_p4.computeME( (string) "qqZ4l_Signal", Set_Of_Arrays, Set_Of_IDs, ME3_p4[count] )) !=0 ) cout << "ERROR CODE in ME3: " << error_value << endl;
		
// 		cout << ME1_p4[count]) << " " << ME2_p4[count]+ME3_p4[count] << endl;
	}
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		mean_ME12_p4 += ME1_p4[count];
		mean_ME22_p4 += ME2_p4[count];
		mean_ME32_p4 += ME3_p4[count];
	}
	mean_ME12_p4 /= correlation_points;
	mean_ME22_p4 /= correlation_points;
	mean_ME32_p4 /= correlation_points;
	
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		corr_nom2 += (ME2_p4[count]+ME3_p4[count] - mean_ME22_p4-mean_ME32_p4)*(ME1_p4[count] - mean_ME12_p4);
		corr_denom1_2 += (ME2_p4[count]+ME3_p4[count] - mean_ME22_p4-mean_ME32_p4)*(ME2_p4[count]+ME3_p4[count] - mean_ME22_p4-mean_ME32_p4);
		corr_denom2_2 += (ME1_p4[count] - mean_ME12_p4)*(ME1_p4[count] - mean_ME12_p4);
	}
	
	if( Show_Basic_Data )
	{
		cout << "Output: correlation_coefficient\n";
		cout << (corr_nom2/sqrt(corr_denom1_2)/sqrt(corr_denom2_2)) << endl;
		cout << "Rel. Difference between means ME1 and ME2+ME3 wrt ME1: " << (mean_ME12_p4-mean_ME22_p4-mean_ME32_p4)/mean_ME12_p4 << endl;
	}
	
	
	
	/// TEST 3
	if( Show_Description ) cout << "\n -- STARTING TEST 3. Invariant mass 80--100 GeV. m4l = mZ -- \n";
	MEKD test3_p4(8.0, "");
	double mean_ME13_p4=0, mean_ME23_p4=0, mean_ME33_p4=0, corr_nom3=0, corr_denom1_3=0, corr_denom2_3=0;
	
	test3_p4.MEKD_MG_Calc.Use_mZ4l_eq_m4l = true;
	
	
	Set_Of_IDs[0] = 11;
	Set_Of_IDs[1] = -11;
	Set_Of_IDs[2] = 11;
	Set_Of_IDs[3] = -11;
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		invariant_mass = (100-80)*(r1.Uniform())+80;
		
		for( unsigned int count2=0; count2<2; count2++)
		{
			p_lepton1[count2] = 2*r1.Uniform()-1;
			p_lepton2[count2] = 2*r1.Uniform()-1;
			p_lepton3[count2] = 2*r1.Uniform()-1;
		}
		
		p_lepton4[0] = -(p_lepton1[0]+p_lepton2[0]+p_lepton3[0]);
		p_lepton4[1] = -(p_lepton1[1]+p_lepton2[1]+p_lepton3[1]);
		p_lepton4[2] = -(p_lepton1[2]+p_lepton2[2]+p_lepton3[2]);
		
		val = sqrt( p_lepton1[0]*p_lepton1[0]+p_lepton1[1]*p_lepton1[1]+p_lepton1[2]*p_lepton1[2]+m_e*m_e );
		val += sqrt( p_lepton2[0]*p_lepton2[0]+p_lepton2[1]*p_lepton2[1]+p_lepton2[2]*p_lepton2[2]+m_e*m_e );
		val += sqrt( p_lepton3[0]*p_lepton3[0]+p_lepton3[1]*p_lepton3[1]+p_lepton3[2]*p_lepton3[2]+m_e*m_e );
		val += sqrt( p_lepton4[0]*p_lepton4[0]+p_lepton4[1]*p_lepton4[1]+p_lepton4[2]*p_lepton4[2]+m_e*m_e );
		
		val = invariant_mass/val;
		
		for( unsigned int count2=1; count2<3; count2++ )
		{
			Set_Of_Arrays[0][count2] = p_lepton1[count2-1]*val;
			Set_Of_Arrays[1][count2] = p_lepton2[count2-1]*val;
			Set_Of_Arrays[2][count2] = p_lepton3[count2-1]*val;
			Set_Of_Arrays[3][count2] = p_lepton4[count2-1]*val;
		}
		Set_Of_Arrays[0][0] = sqrt( p_lepton1[0]*p_lepton1[0] +p_lepton1[1]*p_lepton1[1] +p_lepton1[2]*p_lepton1[2] + m_e*m_e );
		Set_Of_Arrays[1][0] = sqrt( p_lepton2[0]*p_lepton2[0] +p_lepton2[1]*p_lepton2[1] +p_lepton2[2]*p_lepton2[2] + m_e*m_e );
		Set_Of_Arrays[2][0] = sqrt( p_lepton3[0]*p_lepton3[0] +p_lepton3[1]*p_lepton3[1] +p_lepton3[2]*p_lepton3[2] + m_e*m_e );
		Set_Of_Arrays[3][0] = sqrt( p_lepton4[0]*p_lepton4[0] +p_lepton4[1]*p_lepton4[1] +p_lepton4[2]*p_lepton4[2] + m_e*m_e );
		
		if( (error_value=test3_p4.computeME( (string) "ZZ", Set_Of_Arrays, Set_Of_IDs, ME1_p4[count] )) != 0 ) cout << "ERROR CODE in ME1: " << error_value << endl;
		if( (error_value=test3_p4.computeME( (string) "qqZ4l_Background", Set_Of_Arrays, Set_Of_IDs, ME2_p4[count] )) !=0 ) cout << "ERROR CODE in ME2: " << error_value << endl;
		if( (error_value=test3_p4.computeME( (string) "qqZ4l_Signal", Set_Of_Arrays, Set_Of_IDs, ME3_p4[count] )) !=0 ) cout << "ERROR CODE in ME3: " << error_value << endl;
		
// 		cout << ME1_p4[count]) << " " << ME2_p4[count]+ME3_p4[count] << endl;
	}
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		mean_ME13_p4 += ME1_p4[count];
		mean_ME23_p4 += ME2_p4[count];
		mean_ME33_p4 += ME3_p4[count];
	}
	mean_ME13_p4 /= correlation_points;
	mean_ME23_p4 /= correlation_points;
	mean_ME33_p4 /= correlation_points;
	
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		corr_nom3 += (ME2_p4[count]+ME3_p4[count] - mean_ME23_p4-mean_ME33_p4)*(ME1_p4[count] - mean_ME13_p4);
		corr_denom1_3 += (ME2_p4[count]+ME3_p4[count] - mean_ME23_p4-mean_ME33_p4)*(ME2_p4[count]+ME3_p4[count] - mean_ME23_p4-mean_ME33_p4);
		corr_denom2_3 += (ME1_p4[count] - mean_ME13_p4)*(ME1_p4[count] - mean_ME13_p4);
	}
	
	if( Show_Basic_Data )
	{
		cout << "Output: correlation_coefficient\n";
		cout << (corr_nom3/sqrt(corr_denom1_3)/sqrt(corr_denom2_3)) << endl;
		cout << "Rel. Difference between means ME1 and ME2+ME3 wrt ME1: " << (mean_ME13_p4-mean_ME23_p4-mean_ME33_p4)/mean_ME13_p4 << endl;
	}
	
	
	
	/// TEST 4
	if( Show_Description ) cout << "\n -- STARTING TEST 4. Invariant mass 110--130 GeV. m4l = mZ -- \n";
	MEKD test4_p4(8.0, "");
	double mean_ME14_p4=0, mean_ME24_p4=0, mean_ME34_p4=0, corr_nom4=0, corr_denom1_4=0, corr_denom2_4=0;
	
	test4_p4.MEKD_MG_Calc.Use_mZ4l_eq_m4l = true;
	
	
	Set_Of_IDs[0] = 11;
	Set_Of_IDs[1] = -11;
	Set_Of_IDs[2] = 11;
	Set_Of_IDs[3] = -11;
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		invariant_mass = (130-110)*(r1.Uniform())+110;
		
		for( unsigned int count2=0; count2<2; count2++)
		{
			p_lepton1[count2] = 2*r1.Uniform()-1;
			p_lepton2[count2] = 2*r1.Uniform()-1;
			p_lepton3[count2] = 2*r1.Uniform()-1;
		}
		
		p_lepton4[0] = -(p_lepton1[0]+p_lepton2[0]+p_lepton3[0]);
		p_lepton4[1] = -(p_lepton1[1]+p_lepton2[1]+p_lepton3[1]);
		p_lepton4[2] = -(p_lepton1[2]+p_lepton2[2]+p_lepton3[2]);
		
		val = sqrt( p_lepton1[0]*p_lepton1[0]+p_lepton1[1]*p_lepton1[1]+p_lepton1[2]*p_lepton1[2]+m_e*m_e );
		val += sqrt( p_lepton2[0]*p_lepton2[0]+p_lepton2[1]*p_lepton2[1]+p_lepton2[2]*p_lepton2[2]+m_e*m_e );
		val += sqrt( p_lepton3[0]*p_lepton3[0]+p_lepton3[1]*p_lepton3[1]+p_lepton3[2]*p_lepton3[2]+m_e*m_e );
		val += sqrt( p_lepton4[0]*p_lepton4[0]+p_lepton4[1]*p_lepton4[1]+p_lepton4[2]*p_lepton4[2]+m_e*m_e );
		
		val = invariant_mass/val;
		
		for( unsigned int count2=1; count2<3; count2++ )
		{
			Set_Of_Arrays[0][count2] = p_lepton1[count2-1]*val;
			Set_Of_Arrays[1][count2] = p_lepton2[count2-1]*val;
			Set_Of_Arrays[2][count2] = p_lepton3[count2-1]*val;
			Set_Of_Arrays[3][count2] = p_lepton4[count2-1]*val;
		}
		Set_Of_Arrays[0][0] = sqrt( p_lepton1[0]*p_lepton1[0] +p_lepton1[1]*p_lepton1[1] +p_lepton1[2]*p_lepton1[2] + m_e*m_e );
		Set_Of_Arrays[1][0] = sqrt( p_lepton2[0]*p_lepton2[0] +p_lepton2[1]*p_lepton2[1] +p_lepton2[2]*p_lepton2[2] + m_e*m_e );
		Set_Of_Arrays[2][0] = sqrt( p_lepton3[0]*p_lepton3[0] +p_lepton3[1]*p_lepton3[1] +p_lepton3[2]*p_lepton3[2] + m_e*m_e );
		Set_Of_Arrays[3][0] = sqrt( p_lepton4[0]*p_lepton4[0] +p_lepton4[1]*p_lepton4[1] +p_lepton4[2]*p_lepton4[2] + m_e*m_e );
		
		if( (error_value=test4_p4.computeME( (string) "ZZ", Set_Of_Arrays, Set_Of_IDs, ME1_p4[count] )) != 0 ) cout << "ERROR CODE in ME1: " << error_value << endl;
		if( (error_value=test4_p4.computeME( (string) "qqZ4l_Background", Set_Of_Arrays, Set_Of_IDs, ME2_p4[count] )) !=0 ) cout << "ERROR CODE in ME2: " << error_value << endl;
		if( (error_value=test4_p4.computeME( (string) "qqZ4l_Signal", Set_Of_Arrays, Set_Of_IDs, ME3_p4[count] )) !=0 ) cout << "ERROR CODE in ME3: " << error_value << endl;
		
// 		cout << ME1_p4[count]) << " " << ME2_p4[count]+ME3_p4[count] << endl;
	}
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		mean_ME14_p4 += ME1_p4[count];
		mean_ME24_p4 += ME2_p4[count];
		mean_ME34_p4 += ME3_p4[count];
	}
	mean_ME14_p4 /= correlation_points;
	mean_ME24_p4 /= correlation_points;
	mean_ME34_p4 /= correlation_points;
	
	
	for( unsigned int count=0; count < correlation_points; count++ )
	{
		corr_nom4 += (ME2_p4[count]+ME3_p4[count] - mean_ME24_p4-mean_ME34_p4)*(ME1_p4[count] - mean_ME14_p4);
		corr_denom1_4 += (ME2_p4[count]+ME3_p4[count] - mean_ME24_p4-mean_ME34_p4)*(ME2_p4[count]+ME3_p4[count] - mean_ME24_p4-mean_ME34_p4);
		corr_denom2_4 += (ME1_p4[count] - mean_ME14_p4)*(ME1_p4[count] - mean_ME14_p4);
	}
	
	if( Show_Basic_Data )
	{
		cout << "Output: correlation_coefficient\n";
		cout << (corr_nom4/sqrt(corr_denom1_4)/sqrt(corr_denom2_4)) << endl;
		cout << "Rel. Difference between means ME1 and ME2+ME3 wrt ME1: " << (mean_ME14_p4-mean_ME24_p4-mean_ME34_p4)/mean_ME14_p4 << endl;
	}
	
	
	
	/// TEST 5
	if( Show_Description ) cout << "\n -- STARTING TEST 5. Invariant mass 110--130 GeV. Determining average scales for MEs -- \n";
	MEKD test5_p4(8.0, "");
	double mean_ME15_p4=0, mean_ME25_p4=0, mean_ME35_p4=0;
	
	test5_p4.MEKD_MG_Calc.Use_mZ4l_eq_m4l = true;
	
	
	Set_Of_IDs[0] = 11;
	Set_Of_IDs[1] = -11;
	Set_Of_IDs[2] = 11;
	Set_Of_IDs[3] = -11;
	
	for( unsigned int count=0; count < scale_check_points; count++ )
	{
		invariant_mass = (130-110)*(r1.Uniform())+110;
		
		for( unsigned int count2=0; count2<2; count2++)
		{
			p_lepton1[count2] = 2*r1.Uniform()-1;
			p_lepton2[count2] = 2*r1.Uniform()-1;
			p_lepton3[count2] = 2*r1.Uniform()-1;
		}
		
		p_lepton4[0] = -(p_lepton1[0]+p_lepton2[0]+p_lepton3[0]);
		p_lepton4[1] = -(p_lepton1[1]+p_lepton2[1]+p_lepton3[1]);
		p_lepton4[2] = -(p_lepton1[2]+p_lepton2[2]+p_lepton3[2]);
		
		val = sqrt( p_lepton1[0]*p_lepton1[0]+p_lepton1[1]*p_lepton1[1]+p_lepton1[2]*p_lepton1[2]+m_e*m_e );
		val += sqrt( p_lepton2[0]*p_lepton2[0]+p_lepton2[1]*p_lepton2[1]+p_lepton2[2]*p_lepton2[2]+m_e*m_e );
		val += sqrt( p_lepton3[0]*p_lepton3[0]+p_lepton3[1]*p_lepton3[1]+p_lepton3[2]*p_lepton3[2]+m_e*m_e );
		val += sqrt( p_lepton4[0]*p_lepton4[0]+p_lepton4[1]*p_lepton4[1]+p_lepton4[2]*p_lepton4[2]+m_e*m_e );
		
		val = invariant_mass/val;
		
		for( unsigned int count2=1; count2<3; count2++ )
		{
			Set_Of_Arrays[0][count2] = p_lepton1[count2-1]*val;
			Set_Of_Arrays[1][count2] = p_lepton2[count2-1]*val;
			Set_Of_Arrays[2][count2] = p_lepton3[count2-1]*val;
			Set_Of_Arrays[3][count2] = p_lepton4[count2-1]*val;
		}
		Set_Of_Arrays[0][0] = sqrt( p_lepton1[0]*p_lepton1[0] +p_lepton1[1]*p_lepton1[1] +p_lepton1[2]*p_lepton1[2] + m_e*m_e );
		Set_Of_Arrays[1][0] = sqrt( p_lepton2[0]*p_lepton2[0] +p_lepton2[1]*p_lepton2[1] +p_lepton2[2]*p_lepton2[2] + m_e*m_e );
		Set_Of_Arrays[2][0] = sqrt( p_lepton3[0]*p_lepton3[0] +p_lepton3[1]*p_lepton3[1] +p_lepton3[2]*p_lepton3[2] + m_e*m_e );
		Set_Of_Arrays[3][0] = sqrt( p_lepton4[0]*p_lepton4[0] +p_lepton4[1]*p_lepton4[1] +p_lepton4[2]*p_lepton4[2] + m_e*m_e );
		
		if( (error_value=test5_p4.computeME( (string) "ggSpin0Pm", Set_Of_Arrays, Set_Of_IDs, ME1_p4[count] )) != 0 ) cout << "ERROR CODE in ME1: " << error_value << endl;
		if( (error_value=test5_p4.computeME( (string) "ggSpin0Ph", Set_Of_Arrays, Set_Of_IDs, ME2_p4[count] )) !=0 ) cout << "ERROR CODE in ME2: " << error_value << endl;
		if( (error_value=test5_p4.computeME( (string) "ggSpin0M", Set_Of_Arrays, Set_Of_IDs, ME3_p4[count] )) !=0 ) cout << "ERROR CODE in ME3: " << error_value << endl;
		
// 		cout << ME1_p4[count] << " " << ME2_p4[count] << " " << ME3_p4[count] << endl;
	}
	
	for( unsigned int count=0; count < scale_check_points; count++ )
	{
		mean_ME15_p4 += ME1_p4[count];
		mean_ME25_p4 += ME2_p4[count];
		mean_ME35_p4 += ME3_p4[count];
	}
	mean_ME15_p4 /= scale_check_points;
	mean_ME25_p4 /= scale_check_points;
	mean_ME35_p4 /= scale_check_points;
	
	if( Show_Basic_Data )
	{
		cout << "ME | Scale value\n";
		cout << "----------------------------\n";
		cout << "ggSpin0Pm | " << mean_ME15_p4 << endl;
		cout << "ggSpin0Ph | " << mean_ME25_p4 << endl;
		cout << "ggSpin0M | " << mean_ME35_p4 << endl;
	}
	
	
	
	/// TESTs Results
	cout << "\nTEST 1: ";
	if( fabs(corr_nom1/sqrt(corr_denom1_1)/sqrt(corr_denom2_1)) > 0.9 && fabs(mean_ME11_p4-mean_ME21_p4-mean_ME31_p4)/fabs(mean_ME11_p4) < 0.01 )
		cout << "PASSED\n";
	else cout << "FAILED\n";
	
	cout << "TEST 2: ";
	if( fabs(corr_nom2/sqrt(corr_denom1_2)/sqrt(corr_denom2_2)) > 0.9 && fabs(mean_ME12_p4-mean_ME22_p4-mean_ME32_p4)/fabs(mean_ME12_p4) < 0.03 )
		cout << "PASSED\n";
	else cout << "FAILED\n";
	
	cout << "TEST 3: ";
	if( fabs(corr_nom3/sqrt(corr_denom1_3)/sqrt(corr_denom2_3)) > 0.33 && fabs(mean_ME13_p4-mean_ME23_p4-mean_ME33_p4)/fabs(mean_ME13_p4) > 30 )
		cout << "PASSED\n";
	else cout << "FAILED\n";
	
	cout << "TEST 4: ";
	if( fabs(corr_nom4/sqrt(corr_denom1_4)/sqrt(corr_denom2_4)) > 0.60 && fabs(mean_ME14_p4-mean_ME24_p4-mean_ME34_p4)/fabs(mean_ME14_p4) > 200 )
		cout << "PASSED\n";
	else cout << "FAILED\n";

}


#endif