#ifndef MEKD_CalcHEP_Extra_H
#define MEKD_CalcHEP_Extra_H

namespace MEKD_CalcHEP_Extra
{

/// Functions
void Flip_1_and_2_Six_Fourmomenta(double*);
void Reorder_for_4e_4m(double*);

long double Higgs_width(double);	//wrapper
long double Higgs_width_Poly_Fit_Estm(double);
double Higgs_width_Poly_Fit_Estm_Old(double);

}

#endif