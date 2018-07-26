#ifndef MEKD_CalcHEP_PDF_C
#define MEKD_CalcHEP_PDF_C

#include "MEKD_CalcHEP_PDF.h"

/// For PDFs
#include "../PDFTables/pdt.h"
// #define PDTFILE "../PDFTables/cteq6l.pdt" // CalCHEP reads a table for CTEQ6L. You can change PDF set as you want.



/// Part of pdfreader
pdtStr pdtSg, pdtSd, pdtSu, pdtSs, pdtSc,
	pdtSad, pdtSau, pdtSas, pdtSac;



int ConvertID_2_CalcID(int pNum) // To translate KF code (PYTHIA) into a number defined in CalCHEP
{
	switch(pNum)
	{
		case   5: case -5: return 1;
		case   4: case -4: return 2;
		case   3: case -3: return 3;
        case  -1:          return 4;
        case  -2:          return 5;
        case  21:          return 6;
        case   2:          return 7;
        case   1:          return 8;
        case  81:          return 9;
        case -81:          return 10;
        case  83:          return 11;
        case -83:          return 12;
        default:           return 0;
    }
}


double pdfreader(long pNum, double x, double q)
{
	if( pNum==1 ) return interFunc( x, q, &pdtSd );
	if( pNum==2 ) return interFunc( x, q, &pdtSu );
	if( pNum==3 ) return interFunc( x, q, &pdtSs );
	if( pNum==4 ) return interFunc( x, q, &pdtSc );
	if( pNum==-1 ) return interFunc( x, q, &pdtSad );
	if( pNum==-2 ) return interFunc( x, q, &pdtSau );
	if( pNum==-3 ) return interFunc( x, q, &pdtSas );
	if( pNum==-4 ) return interFunc( x, q, &pdtSac );
	if( pNum==21 ) return interFunc( x, q, &pdtSg );
	
	return 0;
}


// double pdfreader(long pNum, double x, double q)
// {
// 	if( pNum==1 ) return 0;
// 	if( pNum==2 ) return 1;
// 	if( pNum==3 ) return 0;
// 	if( pNum==4 ) return 0;
// 	if( pNum==-1 ) return 0;
// 	if( pNum==-2 ) return 1;
// 	if( pNum==-3 ) return 0;
// 	if( pNum==-4 ) return 0;
// 	if( pNum==21 ) return 1;
// 	
// 	return 0;
// }


void Load_pdfreader(char* file)
{
	int err;
	
	err=getPdtData( file, ConvertID_2_CalcID(21), &pdtSg);
	if( err!=0 ) { printf( "Setting-up pdfreader has failed! Error number: %d. Exiting.\n", err ); exit(1); }
	err=getPdtData( file, ConvertID_2_CalcID(1), &pdtSd );
	if( err!=0 ) { printf( "Setting-up pdfreader has failed! Error number: %d. Exiting.\n", err ); exit(1); }
	err=getPdtData( file, ConvertID_2_CalcID(2), &pdtSu );
	if( err!=0 ) { printf( "Setting-up pdfreader has failed! Error number: %d. Exiting.\n", err ); exit(1); }
	err=getPdtData( file, ConvertID_2_CalcID(3), &pdtSs );
	if( err!=0 ) { printf( "Setting-up pdfreader has failed! Error number: %d. Exiting.\n", err ); exit(1); }
	err=getPdtData( file, ConvertID_2_CalcID(4), &pdtSc );
	if( err!=0 ) { printf( "Setting-up pdfreader has failed! Error number: %d. Exiting.\n", err ); exit(1); }
	err=getPdtData( file, ConvertID_2_CalcID(-1), &pdtSad);
	if( err!=0 ) { printf( "Setting-up pdfreader has failed! Error number: %d. Exiting.\n", err ); exit(1); }
	err=getPdtData( file, ConvertID_2_CalcID(-2), &pdtSau);
	if( err!=0 ) { printf( "Setting-up pdfreader has failed! Error number: %d. Exiting.\n", err ); exit(1); }
	err=getPdtData( file, ConvertID_2_CalcID(-3), &pdtSas);
	if( err!=0 ) { printf( "Setting-up pdfreader has failed! Error number: %d. Exiting.\n", err ); exit(1); }
	err=getPdtData( file, ConvertID_2_CalcID(-4), &pdtSac);
	if( err!=0 ) { printf( "Setting-up pdfreader has failed! Error number: %d. Exiting.\n", err ); exit(1); }
}


void Unload_pdfreader()
{
	freePdtData( &pdtSg );
	freePdtData( &pdtSd );
	freePdtData( &pdtSu );
	freePdtData( &pdtSs );
	freePdtData( &pdtSc );
	freePdtData( &pdtSad);
	freePdtData( &pdtSau);
	freePdtData( &pdtSas);
	freePdtData( &pdtSac);
}



#endif