//
//  hggeffective.cpp
//  
//
//  Created by Myeonghun Park on 9/2/12.
//
//  Modified by Aurelijus Rinkevicius on 2013-05-06
//
#ifndef _hggeffective_h
#define _hggeffective_h

#include <cstdio>
#include <cmath>
#include <complex>

using namespace std;


/// Definition
#define M_NPI 3.14159265358979323846


/// Functions
double fiRe(double tau);
double fiIm(double tau);
double HggFr(double tau);
double HggFi(double tau);

double LmbdGG(double Mh);

complex<double> HggF(double tau);


#endif
