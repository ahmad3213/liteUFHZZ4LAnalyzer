#ifndef MEKD_TEST_GEN_INTERFERENCE_CPP
#define MEKD_TEST_GEN_INTERFERENCE_CPP

#include "MEKD_Test.h"

double PI=2*asin(1.0);


//////////////////////////////////
/// TEST Block Interference Check
//////////////////////////////////


int MEKD_Test_Gen_Interference_Test1()
{	
	if( Show_Description )
	{
		cout << "\n -------------------------------------------------------------------- \n";
		cout << " -- Testing ME Interference effects. ggSpin0Pm and ggSpin0Ph. Generating 4e events. -- \n";
		cout << " -------------------------------------------------------------------- \n";
	}
	
	if( Show_Description ) cout << "Testing intermodel interference using int computeME( string, vector<double*>, vector<int>, double& );\n";
	
	double ME1[interference_gen_points], ME2[interference_gen_points], ME3[interference_gen_points], ME_interf[interference_steps][interference_gen_points];
	double val, invariant_mass, p_lepton1[3], p_lepton2[3], p_lepton3[3], p_lepton4[3], angle;
	TRandom3 r1;
	
	/// TEST 1
	if( Show_Description ) cout << "\n -- STARTING TEST 1. Invariant mass 120--130 GeV. -- \n";
	MEKD test1(8.0, "");
	double mean_interf[interference_steps];
	for( unsigned int count_step=0; count_step<interference_steps; count_step++ ) mean_interf[count_step] = 0;
	
	
	
	
	Set_Of_IDs[0] = 11;
	Set_Of_IDs[1] = -11;
	Set_Of_IDs[2] = 11;
	Set_Of_IDs[3] = -11;
	
	for( unsigned int count=0; count < interference_gen_points; count++ )
	{
		invariant_mass = (130-120)*(r1.Uniform())+120;
		
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
		
		if( (error_value=test1.computeME( (string) "ggSpin0Pm", Set_Of_Arrays, Set_Of_IDs, ME1[count] )) != 0 ) cout << "ERROR CODE in ME1: " << error_value << endl;
		if( (error_value=test1.computeME( (string) "ggSpin0Ph", Set_Of_Arrays, Set_Of_IDs, ME2[count] )) !=0 ) cout << "ERROR CODE in ME2: " << error_value << endl;
		
		if( Show_Debug ) cout << ME1[count] << " " << ME2[count] << endl;
		
		for( unsigned int count_step=0; count_step<interference_steps; count_step++ )
		{
			angle = static_cast<double>((PI/2.0)*count_step/interference_steps);
			test1.Mix_Spin0( complex<double>( cos(angle), 0.0), complex<double>(sin(angle), 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0) );
			if( (error_value=test1.computeME( (string) "ggSpin0", Set_Of_Arrays, Set_Of_IDs, ME3[count] )) !=0 ) cout << "ERROR CODE in ME3: " << error_value << endl;
		
			ME_interf[count_step][count] = 2.0*( ME3[count]
				- ME1[count]*cos(angle)*cos(angle)
				- ME2[count]*sin(angle)*sin(angle) )/( ME3[count] + ME1[count]*cos(angle)*cos(angle) + ME2[count]*sin(angle)*sin(angle) );
				
			if( Show_Advanced_Data ) cout << "angle cos2*ggSpin0Pm sin2*ggSpin0Ph Intf: " << angle << " " <<
				ME1[count]*cos(angle)*cos(angle) << " " <<
				ME2[count]*sin(angle)*sin(angle) << " " <<
				ME3[count] << " " <<
				ME_interf[count_step][count] << " " << endl;
		}
	}
	
	for( unsigned int count_step=0; count_step<interference_steps; count_step++ )
	{
		for( unsigned int count=0; count < interference_gen_points; count++ )
			mean_interf[count_step] += ME_interf[count_step][count];
		mean_interf[count_step] /= interference_gen_points;
	}
	
	if( Show_Basic_Data )
	{
		cout << "Angle PI-1_atan_ratio Interference\n";
		for( unsigned int count_step=0; count_step<interference_steps; count_step++ )
		{
			angle = static_cast<double>((PI/2.0)*count_step/interference_steps);
			cout << angle << " " << (1/PI)*angle<< " " << mean_interf[count_step] << " " << endl;
		}
	}
	
	if( interference_steps > 1 && fabs(mean_interf[1]) > 0 ) return 0;
	return 1;
}



int MEKD_Test_Gen_Interference_Test2()
{	
	if( Show_Description )
	{
		cout << "\n -------------------------------------------------------------------- \n";
		cout << " -- Testing ME Interference effects. ggSpin0Pm and ggSpin0M. Generating 4e events. -- \n";
		cout << " -------------------------------------------------------------------- \n";
	}
	
	if( Show_Description ) cout << "Testing intermodel interference using int computeME( string, vector<double*>, vector<int>, double& );\n";
	
	double ME1[interference_gen_points], ME2[interference_gen_points], ME3[interference_gen_points], ME_interf[interference_steps][interference_gen_points];
	double val, invariant_mass, p_lepton1[3], p_lepton2[3], p_lepton3[3], p_lepton4[3], angle;
	TRandom3 r1;
	
	/// TEST 1
	if( Show_Description ) cout << "\n -- STARTING TEST 1. Invariant mass 120--130 GeV. -- \n";
	MEKD test(8.0, "");
	double mean_interf[interference_steps];
	for( unsigned int count_step=0; count_step<interference_steps; count_step++ ) mean_interf[count_step] = 0;
	
	
	
	
	Set_Of_IDs[0] = 11;
	Set_Of_IDs[1] = -11;
	Set_Of_IDs[2] = 11;
	Set_Of_IDs[3] = -11;
	
	for( unsigned int count=0; count < interference_gen_points; count++ )
	{
		invariant_mass = (130-120)*(r1.Uniform())+120;
		
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
		
		if( (error_value=test.computeME( (string) "ggSpin0Pm", Set_Of_Arrays, Set_Of_IDs, ME1[count] )) != 0 ) cout << "ERROR CODE in ME1: " << error_value << endl;
		if( (error_value=test.computeME( (string) "ggSpin0M", Set_Of_Arrays, Set_Of_IDs, ME2[count] )) !=0 ) cout << "ERROR CODE in ME2: " << error_value << endl;
		
		if( Show_Debug ) cout << ME1[count] << " " << ME2[count] << endl;
		
		for( unsigned int count_step=0; count_step<interference_steps; count_step++ )
		{
			angle = static_cast<double>((PI/2.0)*count_step/interference_steps);
			test.Mix_Spin0( complex<double>( cos(angle), 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(sin(angle), 0.0) );
			if( (error_value=test.computeME( (string) "ggSpin0", Set_Of_Arrays, Set_Of_IDs, ME3[count] )) !=0 ) cout << "ERROR CODE in ME3: " << error_value << endl;
		
			ME_interf[count_step][count] = 2.0*( ME3[count]
				- ME1[count]*cos(angle)*cos(angle)
				- ME2[count]*sin(angle)*sin(angle) )/( ME3[count] + ME1[count]*cos(angle)*cos(angle) + ME2[count]*sin(angle)*sin(angle) );
				
			if( Show_Advanced_Data ) cout << "angle cos2*ggSpin0Pm sin2*ggSpin0M Intf: " << angle << " " <<
				ME1[count]*cos(angle)*cos(angle) << " " <<
				ME2[count]*sin(angle)*sin(angle) << " " <<
				ME3[count] << " " <<
				ME_interf[count_step][count] << endl;
		}
	}
	
	for( unsigned int count_step=0; count_step<interference_steps; count_step++ )
	{
		for( unsigned int count=0; count < interference_gen_points; count++ )
			mean_interf[count_step] += ME_interf[count_step][count];
		mean_interf[count_step] /= interference_gen_points;
	}
	
	if( Show_Basic_Data )
	{
		cout << "Angle PI-1_atan_ratio Interference\n";
		for( unsigned int count_step=0; count_step<interference_steps; count_step++ )
		{
			angle = static_cast<double>((PI/2.0)*count_step/interference_steps);
			cout << angle << " " << (1/PI)*angle << " " << mean_interf[count_step] << endl;
		}
	}
	
	if( interference_steps > 1 && fabs(mean_interf[1]) > 0 ) return 0;
	return 1;
}


#endif