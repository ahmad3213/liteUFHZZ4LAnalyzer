#ifndef MEKD_CalcHEP_PDF_h
#define MEKD_CalcHEP_PDF_h


/// To translate KF code (PYTHIA) into a number defined in CalCHEP
int ConvertID_2_CalcID(int pNum);

double pdfreader(long pNum, double x, double q);

void Load_pdfreader(char *file);
void Unload_pdfreader();


#endif