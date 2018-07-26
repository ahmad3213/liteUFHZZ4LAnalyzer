//
//  hggeffective.cpp
//  
//
//  Created by Myeonghun Park on 9/2/12.
//
//  Modified by Aurelijus Rinkevicius on 2013-05-06
//
#ifndef _hggeffective_cpp
#define _hggeffective_cpp

#include "hggeffective.h"


using namespace std;



/// Variables

double EE= 3.12300000000e-1,
	SW= 4.80800000000e-1,
	Me= 5.109989e-4,
	Mm= 1.056584e-1,
	Mcp= 1.40000000000e0,
	Mtp= 1.77500000000e2,
	MZ= 9.118760e1,
	wZ= 2.495200e0;
double Mbp = 4.891000e0;
double CW = sqrt(1-SW*SW);
double MW = MZ*sqrt(1-SW*SW);
double alphaQh = 0.1184;
double nfh = 5.0;



/// Functions

double fiRe(double tau)
{
	double x;
	
	if(tau<1)
	{
		x=asin(sqrt(tau));
		
		return x*x;
	}else if(tau==1) return 0;
	else
	{
		x=sqrt(1-1/tau);
		x=log((1+x)/(1-x));
		
		return -0.25*(x*x - M_NPI*M_NPI);
	}
}


double fiIm(double tau)
{
	double x;
	
	if(tau<=1) return 0;
	else
	{
		x = sqrt(1-1/tau);
		x = log((1+x)/(1-x));
		
		return ( 0.5*x*M_NPI );
	}
}


double HggFr(double tau) { return  2.0*(tau+(tau-1)*fiRe(tau))/(tau*tau); }


double HggFi(double tau) { return  2.0*(tau-1)*fiIm(tau)/(tau*tau); }


double LmbdGG(double Mh)
{
	double FFeven=1+alphaQh/M_NPI*(95.0/4.0-7.0/6.0*nfh)+alphaQh*alphaQh/pow(M_NPI,2)
		*(370.0-47.0*nfh+0.9*pow(nfh,2)-(19.0/8.0+2.0*nfh/3.0)*2.0*log(Mtp/Mh));
	
	complex<double> loops=HggF(pow(Mh/2.0/Mcp,2))+HggF(pow(Mh/2.0/Mbp,2))+HggF(pow(Mh/2.0/Mtp,2));
	
	return ( alphaQh/(8*M_NPI)*1/2*(-EE)/(2*MW*SW)*abs(loops)*sqrt(FFeven) );
}


complex<double> HggF(double tau)
{
	complex<double> value(HggFr(tau),HggFi(tau));
	
	return value;
}


#endif
