#ifndef MEKD_MG_Boosts_h
#define MEKD_MG_Boosts_h

#include <cmath>
#include <cstdio>


void Boost_2p_2_pT0(double mass0, double *pi0, double mass1, double *pi1);
void Boost_3p_2_pT0(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2);
void Boost_4p_2_pT0(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3);
void Boost_5p_2_pT0(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3, double mass4, double *pi4);
void Boost_2p_and_2p_2_pT0(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3);
void Boost_3p_and_2p_2_pT0(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3, double mass4, double *pi4);
void Boost_4p_and_2p_2_pT0(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3, double mass4, double *pi4, double mass5, double *pi5);
void Boost_5p_and_2p_2_pT0(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3, double mass4, double *pi4, double mass5, double *pi5, double mass6, double *pi6);
void Boost2CM(double mass0, double *pi0, double mass1, double *pi1);	// 2 particles
void Boost2CM(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2);	// 3 particles
void Boost2CM(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3);	// 4 particles
void Boost2CM(double mass0, double *pi0, double mass1, double *pi1, double mass2, double *pi2, double mass3, double *pi3, double mass4, double *pi4);	// 5 particles
void Boost(double *vector, double *boost);
void Boost_long(long double *vector, long double *boost);



#endif