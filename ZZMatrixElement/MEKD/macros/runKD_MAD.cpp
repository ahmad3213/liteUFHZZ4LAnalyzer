// C includes
#include <unistd.h>	// neded for getopt, optarg,

// C++ includes
#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

/// ME calculator
#include "../src/MEKD.cpp"

/// ROOT includes
#ifdef MEKD_with_ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TLorentzVector.h"
#include "TSystem.h"	// For the workaround for >> "Warning in <TClass::TClass>: no dictionary for class...
#endif

////////////////////////////////////////////////////////////////////////////////

using namespace std;

// Default parameters
double ppEnergy         = 8; // TeV
string inputFile        = "";
string inputTree        = "passedEvents";
string resonanceName    = "SMHiggs";
string backgroundName   = "ZZ";
string pdfInclude       = "CTEQ6L";
string logName          = "";

////////////////////////////////////////////////////////////////////////////////
// Methods
void usage( int status = 0 );
int calculateKD(string inputFile, string tree, string resonanceName, string pdfInclude, string logName);
int calculateKD_Madgraph(string inputFile, string resonanceName, string pdfInclude, string logName);
int calculateKD_Madgraph(string file, string tree, string resonanceName, string pdfInclude, string logName);
void Change_Momentum_Coord_Convention(double*);
////////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{   
	gSystem->Load("libTree");	// A workaround for >> "Warning in <TClass::TClass>: no dictionary for class ...
    // Parse options
	char ch;
	while ((ch = getopt(argc, argv, "f:t:s:x:p:l:?:h")) != -1 )
	{
		switch (ch)
		{
			case 'f': inputFile = string(optarg); break;
			case 't': inputTree = string(optarg); break;
			case 's': ppEnergy = atof(optarg); break;
			case 'x': resonanceName = string(optarg); break;
			case 'p': pdfInclude = string(optarg); break;
			case 'l': logName = string(optarg); break;
			case '?':
			case 'h': usage(0); break;
			default:
				cerr << "\n[ MEKD - ERROR: Unknown option '" << optarg << "' ]" <<  endl  << endl;
				usage(-1);
		}
	}
    
    // Basic sanity checks
    // Basic sanity checks
    if(inputFile == "" || resonanceName == "") usage(-1);
    if (pdfInclude != "CTEQ6L") pdfInclude = "no PDFs"; // only CTEQ6L supported at the moment
    
    // Print chosen options
	cout << endl;
	cout << "----------------------------------------------------------" << endl;
	cout << " Matrix Element (ME) Kinematic Discriminant (KD) producer " << endl;
	cout << "----------------------------------------------------------" << endl;
	cout << "  Input event file:   " << inputFile << endl;
	cout << "  sqrt(s):            " << ppEnergy  << " TeV" << endl;
	cout << "  Resonance type:     " << resonanceName << endl;
	cout << "  PDFs included:      " << pdfInclude << endl;
	cout << "----------------------------------------------------------" << endl;

    // Run ME KD with chosen options
    int err = calculateKD(inputFile, inputTree, resonanceName, pdfInclude, logName);
	
    return err;
}


//------------------------------------------------------------------------------------------------------------
// Print short user instructions
void usage( int status) {
	cout << "----------------------------------------------------------------------------------------------------------" << endl;
	cout << "                        Matrix Element (ME) Kinematic Discriminant (KD) producer                          " << endl;
	cout << "----------------------------------------------------------------------------------------------------------" << endl;
	cout << "  Usage:                                                                                                  " << endl;
#ifdef MEKD_with_ROOT
	cout << "     ./runKD_MAD [-f input_file] [-t root_tree] [-s sqrts] [-x x_resonance] [-p pdf_include] [-l log_file]" << endl;
#else
	cout << "     ./runKD_MAD [-f input_file] [-s sqrts] [-x x_resonance] [-p pdf_include] [-l log_file]               " << endl;
#endif
	cout << endl;
	cout << "  where:" << endl;
	cout << "     input_file    name of the input file (string, REQUIRED)                                    " << endl << endl;
#ifdef MEKD_with_ROOT
	cout << "     root_tree     name of the input root tree (string, DEFAULT = 'passedEvents')               " << endl << endl;
	cout << "                   Available only if compiled with ROOT.                                        " << endl;
#endif
	cout << "     sqrts         pp collsiion energy in TeV (double, DEFAULT = 8)                             " << endl;
	cout << "                   Available options: 7, 8.                                                     " << endl << endl;
	cout << "     x_resonance   choice of the signal spin-0 resonance (string, DEFAULT = 'SMHiggs')          " << endl;
	cout << "                   Available options: 'SMHiggs', 'Higgs0M', 'Graviton2PM' and 'Custom'.         " << endl << endl;
	cout << "     pdf_include   Name of PDFs, PDFs not used if pdf_include = '' (string, DEFAULT= 'CTEQ6L')  " << endl;
	cout << "                   Available PDFs: 'CTEQ6L' (NONE = '').                                        " << endl << endl;
	cout << "     log_file      name of the log file, no logging if log_file = '' (string, DEFAULT = '')     " << endl;
	cout << endl;
	cout << "----------------------------------------------------------------------------------------------------------" << endl;
	cout << endl;
	exit(status);
}

//------------------------------------------------------------------------------------------------------------
// Save position of current standard output (and redirect to log file)
int redirectAndSaveStdout(fpos_t &pos, string logName) {
    fgetpos(stdout, &pos);
    int fd = dup(fileno(stdout));
    const char* nullDevName = "/dev/null"; // set the name of nul device
    if (fopen (nullDevName, "w") == NULL) nullDevName = "nul"; // change the name if WIN (just simple check)
    if (logName == "") freopen(nullDevName, "w", stdout);
    else freopen((const_cast<char*>(logName.c_str())), "w", stdout);
    return fd;
}

//------------------------------------------------------------------------------------------------------------
// Flush stdout, close file and restore standard output to stdout
void closeAndRestoreStdout(int fd, fpos_t pos) {
    fflush(stdout);
    dup2(fd, fileno(stdout));
    close(fd);
    clearerr(stdout);
    fsetpos(stdout, &pos);
}

//------------------------------------------------------------------------------------------------------------
// Import/read 4-momentum components for one lepton from the input file
void import4Momenta(istringstream &data, double *carray, int n){
    char temp[100];
    double dtemp;
    for (int i=0; i<n; i++) {
        data>>temp;
        dtemp=atof(temp);
        carray[i]=dtemp;
    }
}

//------------------------------------------------------------------------------------------------------------
// Import/read 4 lepton PDG IDs from the input file
void import4IDs(istringstream &data, int &idL1, int &idL2, int &idL3, int &idL4){
    char temp[100];
    data>>temp; idL1=atoi(temp);
    data>>temp; idL2=atoi(temp);
    data>>temp; idL3=atoi(temp);
    data>>temp; idL4=atoi(temp);
}

//------------------------------------------------------------------------------------------------------------
// Import/read lwpton IDs and 4-momenta from the input file
void readLeptonIDsAnd4Momenta(ifstream &data, int &idL1, int &idL2, int &idL3, int &idL4, double p3in[], double p4in[], double p5in[], double p6in[]) {
//     double p1in[4],p2in[4];
    string buff;
    getline(data,buff);
    istringstream iss(buff);
    
    // read lepton PDG IDs
    import4IDs(iss, idL1, idL2, idL3, idL4);
    // read lepton 4-momenta
    import4Momenta(iss,p3in,4);
    import4Momenta(iss,p4in,4);
    import4Momenta(iss,p5in,4);
    import4Momenta(iss,p6in,4);
}

//------------------------------------------------------------------------------------------------------------
// Run the selected ME method (only Madgraph supported in this macro)
int calculateKD(string inputFile, string inputTree, string resonanceName, string pdfInclude, string logName) {
	// Check if a supported resonance name is provided and run KD calculator based on the library meMadgraph.h (which uses some Madgraph libraries)
	if( resonanceName == "SMHiggs" || resonanceName == "Higgs0M" || resonanceName == "Graviton2PM" || resonanceName == "Custom" ){
		if (inputFile.substr(inputFile.length()-5,inputFile.length()-1) == ".root"){
#ifdef MEKD_with_ROOT
			return calculateKD_Madgraph(inputFile, inputTree, resonanceName, pdfInclude, logName);
#else
			cerr << "\n[ MEKD - ERROR: ROOT is not included!. Exiting. ]\n\n";
			exit(EXIT_FAILURE);
#endif
		}else
			return calculateKD_Madgraph(inputFile, resonanceName, pdfInclude, logName);
	} else {
		cerr << "\n[ MEKD - ERROR: Resonance of type '"+resonanceName+"' not supported ]" << endl << endl;
		return -1;
	}
}

//------------------------------------------------------------------------------------------------------------
// Read the kinematic information for sample of events from the input file and calculate MEs and KD
int calculateKD_Madgraph(string inputFile, string resonanceName, string pdfInclude, string logName)
{
	double p1[4], p2[4], p3[4], p4[4];
	int idL1, idL2, idL3, idL4;
    double ME_XtoZZ, ME_ZZ, KD;
	
    //Save position of current standard output and redirect it to user log
    fpos_t pos;
    int fd = redirectAndSaveStdout(pos,logName);
	
	// open the input and output files
    ifstream data;
    ofstream dataout;
	
    data.open(inputFile.c_str());
    // cross check if input file is readable
    if (!data.good()) {
        cerr << "\n[ MEKD - ERROR: File '"+inputFile+"' is not readable. ]" << endl << endl;
        return -1;
    }
    
    string ouputFile = inputFile.substr(0,inputFile.length()-4);
	ouputFile += "_withDiscriminant.dat";
	
    dataout.open( ouputFile.c_str() );
	dataout.setf( ios_base::scientific );
	dataout.precision( 17 );
	
    // Initialise MEKD producer
	MEKD* MEKDwithPDFs = new MEKD(ppEnergy,pdfInclude);
	
    // scan the input file and calculate MEs and KD for each event
	while( !(data.eof()) ){
        // Read in lepton IDs and 4-momenta from the input file
        readLeptonIDsAnd4Momenta(data, idL1, idL2, idL3, idL4, p1, p2, p3, p4);
		Change_Momentum_Coord_Convention( p1 );
		Change_Momentum_Coord_Convention( p2 );
		Change_Momentum_Coord_Convention( p3 );
		Change_Momentum_Coord_Convention( p4 );
        
        cout << "resonanceName: " << resonanceName << ", backgroundName: " << backgroundName << endl;
        
        // Calclualte MEs and KD
        int err = MEKDwithPDFs->computeKD(resonanceName, backgroundName, p1, idL1, p2, idL2, p3, idL3, p4, idL4, KD, ME_XtoZZ, ME_ZZ);

		// store MEs and KD in the output file
		dataout << ME_XtoZZ << " " << ME_ZZ <<  " " << KD << endl;
		
	} // scan the input file
	
	// Close inout and output files
	data.close(); 
	dataout.close();
	
	//Close log and restore standard output
	closeAndRestoreStdout(fd, pos);
	
	cout << "----------------------------------------------------------" << endl;
	cout << "  Output file: " << ouputFile << endl;
	cout << "----------------------------------------------------------" << endl << endl;
	
	return 0;
}

////////////////////////////////////////////////////////////////////////////////

#ifdef MEKD_with_ROOT

//------------------------------------------------------------------------------------------------------------
// Read the input tree with 4-momenta and calculate ME
int calculateKD_Madgraph(string file, string tree, string resonanceName, string pdfInclude, string logName)
{
	// Variables
	double e1, e2, e3, e4, px1, px2, px3, px4;
	double py1, py2, py3, py4, pz1, pz2, pz3, pz4;
	double m4l;
	int idL1, idL2, idL3, idL4;
    double ME_XtoZZ, ME_ZZ, KD;
    
    TString inputFile(file);
    TString treeName(tree);
    
    // intitialise variables and perform sanity checks
    vector<double*> p(1, new double[4]);
    
    //Save position of current standard output and redirect to log
    fpos_t pos;
    int fd = redirectAndSaveStdout(pos,logName);
    
    TFile f(inputFile);
    if (f.IsZombie()) {
        cerr << "\n [ MEKD - ERROR: File '"+inputFile+"' does not exist ]" << endl << endl;
        return -1;
    }
    
    TFile* newFile = new TFile( inputFile(0,inputFile.Length()-5) + "_withDiscriminant.root","RECREATE"); //replace '.root' with '_withDiscriminants.root'
    
    /////////////////////////////////////////////////////////////
    
    TFile* sigFile = new TFile(inputFile,"UPDATE");
    TTree* sigTree = NULL;
    if(sigFile) sigTree = (TTree*) sigFile->Get(treeName);
    if(!sigTree){
        cerr << "\n [ MEKD - ERROR: Can not find the tree '" + treeName + "' ]" << endl << endl;
        return -1;
    }

    sigTree->SetBranchAddress("mass4l",&m4l);
    sigTree->SetBranchAddress("EL1",&e1);
    sigTree->SetBranchAddress("EL2",&e2);
    sigTree->SetBranchAddress("EL3",&e3);
    sigTree->SetBranchAddress("EL4",&e4);
    sigTree->SetBranchAddress("pXL1",&px1);
    sigTree->SetBranchAddress("pXL2",&px2);
    sigTree->SetBranchAddress("pXL3",&px3);
    sigTree->SetBranchAddress("pXL4",&px4);
    sigTree->SetBranchAddress("pYL1",&py1);
    sigTree->SetBranchAddress("pYL2",&py2);
    sigTree->SetBranchAddress("pYL3",&py3);
    sigTree->SetBranchAddress("pYL4",&py4);
    sigTree->SetBranchAddress("pZL1",&pz1);
    sigTree->SetBranchAddress("pZL2",&pz2);
    sigTree->SetBranchAddress("pZL3",&pz3);
    sigTree->SetBranchAddress("pZL4",&pz4);
    sigTree->SetBranchAddress("idL1",&idL1);
    sigTree->SetBranchAddress("idL2",&idL2);
    sigTree->SetBranchAddress("idL3",&idL3);
    sigTree->SetBranchAddress("idL4",&idL4);
    
    TBranch *me_XtoZZ = sigTree->Branch("ME_HiggsSM",&ME_XtoZZ,"ME_HiggsSM/D");
    TBranch *me_ZZ = sigTree->Branch("ME_ZZ",&ME_ZZ,"ME_ZZ/D");
    TBranch *me_KD = sigTree->Branch("MEKD",&KD,"MEKD/D");
    
    /////////////////////////////////////////////////////////////

    // Initialise MEKD producer
	MEKD* MEKDwithPDFs = new MEKD(ppEnergy,pdfInclude);
	
    Long64_t nentries = sigTree->GetEntries();
    for(int iEvt=0; iEvt < nentries; iEvt++)
	{
        sigTree->GetEntry(iEvt);
        
        // input 4-momenta
        TLorentzVector p1, p2, p3, p4;
        p1.SetPxPyPzE(px1,py1,pz1,e1);
        p2.SetPxPyPzE(px2,py2,pz2,e2);
        p3.SetPxPyPzE(px3,py3,pz3,e3);
        p4.SetPxPyPzE(px4,py4,pz4,e4);
        
        // compute MEs and KD
        int err = MEKDwithPDFs->computeKD(resonanceName, backgroundName, p1, idL1, p2, idL2, p3, idL3, p4, idL4, KD, ME_XtoZZ, ME_ZZ);

        // fill the tree
        me_XtoZZ->Fill();
        me_ZZ->Fill();
        me_KD->Fill();
    }
    
    /////////////////////////////////////////////////////////////
    // write the new tree into a new file
    if (treeName.Contains("/")) { // works only for depth 1
        TString newRootDir = treeName(0,treeName.Index("/"));
        newFile->mkdir(newRootDir.Data());
        newFile->cd(newRootDir.Data());
    } else {
        newFile->cd();
    }
    TString newTreeName = treeName(treeName.Index("/")+1,treeName.Length());
    sigTree->CloneTree()->Write(newTreeName, TObject::kOverwrite);
    newFile->Close();
    
    //Close log and restore standard output
    closeAndRestoreStdout(fd, pos);

    cout << "----------------------------------------------------------" << endl;
    cout << "  Output file: " << newFile->GetName() << endl;
    cout << "----------------------------------------------------------" << endl << endl;

    return 0;
    
}

#endif

void Change_Momentum_Coord_Convention(double *input_p)
{
	double buffer[4];
	
	buffer[0] = input_p[3];
	buffer[1] = input_p[0];
	buffer[2] = input_p[1];
	buffer[3] = input_p[2];
	
	for( int i=0; i<4; i++ ) input_p[i]=buffer[i];
}

