#ifndef ZZ4LAnalysisTree_h
#define ZZ4LAnalysisTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <map>
#include <utility>
#include <iterator>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"

#include "TSystem.h"
#include "TStyle.h"
#include "TPaveText.h"

#include "TPaveLabel.h"
#include "TLegend.h"

#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>
#include <fstream>
#include "TRandom3.h"

#include <algorithm>

#include "TColor.h"

using namespace std;  

////
// output tree
bool passedFullSelection, passedZ4lSelection;
bool passedZXCRSelection, passedZ4lZXCRSelection;
bool passSmartCut;
int nZXCRFailedLeptons;
int finalState;
int TnP_l1_pdgId; //[nTnPfloat TnP_l1_pt; //[nTnP]
float TnP_l1_pt; //[nTnP]
float TnP_l1_eta; //[nTnP]
float TnP_l1_phi; //[nTnP]
float TnP_l1_mass; //[nTnP]
int TnP_l1_tightId; //[nTnP]
float TnP_l1_sip3d; //[nTnP]
float TnP_l1_p4WithFSR_pt; //[nTnP]
float TnP_l1_p4WithFSR_eta; //[nTnP]
float TnP_l1_p4WithFSR_phi; //[nTnP]
float TnP_l1_p4WithFSR_mass; //[nTnP] 
int TnP_l1_looseId;
float TnP_l1_relIsoAfterFSR;
float TnP_l1_relIsobeforeFSR;
float TnP_l1_chargedHadIso03;
float met;
float D_bkg_kin;
int njets_pt30_eta4p7;
int njets_pt30_eta2p5;
int nbjets_pt30_eta4p7;
int EventCat;
bool passTrig;
// input tree variables
std::string *triggersPassed;
ULong64_t Run, LumiSect, Event;
bool passedTrig;
bool passedFiducialSelection;
float dataMCWeight, genWeight, pileupWeight, crossSection, sumweight ;
float k_qqZZ_qcd_M,k_qqZZ_ewk,k_ggZZ;
float sumW;
int nVtx, nInt; //nPV
float me_qqZZ_MCFM;
std::vector<float>* lep_mass;
std::vector<float> *lep_pt; std::vector<float> *lep_eta; std::vector<float> *lep_phi;
std::vector<float>* lepFSR_mass;
std::vector<float> *lepFSR_pt; std::vector<float> *lepFSR_eta; std::vector<float> *lepFSR_phi;

std::vector<int> *lep_tightId;
std::vector<int> *lep_ecalDriven;
std::vector<int> *lep_id;
std::vector<int> *lep_Sip;
std::vector<float> *lep_dxy;
std::vector<float> *lep_dz;
int lep_Hindex[4];
std::vector<float> *lep_RelIso;
std::vector<float> *lep_isoCH;
std::vector<float> *lep_RelIsoNoFSR;

std::vector<float> *lep_pterr;
std::vector<float> *lep_dataMC;

std::vector<float> *jet_mass;
std::vector<float> *jet_pt; std::vector<float> *jet_eta; std::vector<float> *jet_phi;
std::vector<int> *jet_iscleanH4l; std::vector<float> *jet_QGTagger; std::vector<float> *jet_csvv2;

std::vector<int> *fsrPhotons_lepindex;
std::vector<float> *fsrPhotons_pt; std::vector<float> *fsrPhotons_eta; std::vector<float> *fsrPhotons_phi;
std::vector<float> *fsrPhotons_pterr;

namespace ZZ4LAnalysisTree {

    void setAddresses(TTree* tree, TString filename){
        
        tree->SetBranchStatus("*",0);

        tree->SetBranchStatus("Run",1);
        tree->SetBranchStatus("LumiSect",1);
        tree->SetBranchStatus("Event",1);
        tree->SetBranchStatus("nVtx",1);
        tree->SetBranchStatus("nInt",1);
        tree->SetBranchStatus("genWeight",1);
        tree->SetBranchStatus("crossSection",1);
        tree->SetBranchStatus("pileupWeight",1);
        tree->SetBranchStatus("passedTrig",1);
        tree->SetBranchStatus("triggersPassed",1);
        tree->SetBranchStatus("passedFiducialSelection",1);
        tree->SetBranchStatus("lep_id",1);
        tree->SetBranchStatus("lep_Sip",1);
        tree->SetBranchStatus("lep_dxy",1);
        tree->SetBranchStatus("lep_dz",1);
        tree->SetBranchStatus("lep_tightId",1);
        tree->SetBranchStatus("lep_ecalDriven",1);
        tree->SetBranchStatus("lep_pt",1);
        tree->SetBranchStatus("lep_pterr",1);
        tree->SetBranchStatus("lep_eta",1);
        tree->SetBranchStatus("lep_phi",1);
        tree->SetBranchStatus("lep_mass",1);
        tree->SetBranchStatus("lep_RelIso",1);
        tree->SetBranchStatus("lep_isoCH",1);
        tree->SetBranchStatus("lep_RelIsoNoFSR",1);
        tree->SetBranchStatus("lepFSR_pt",1);
        tree->SetBranchStatus("lepFSR_eta",1);
        tree->SetBranchStatus("lepFSR_phi",1);
        tree->SetBranchStatus("lepFSR_mass",1);        
        tree->SetBranchStatus("jet_pt",1);
        tree->SetBranchStatus("jet_eta",1);
        tree->SetBranchStatus("jet_phi",1);
        tree->SetBranchStatus("jet_mass",1);
        tree->SetBranchStatus("jet_iscleanH4l",1);
        tree->SetBranchStatus("jet_QGTagger",1);
        tree->SetBranchStatus("jet_csvv2",1);
        tree->SetBranchStatus("fsrPhotons_pt",1);
        tree->SetBranchStatus("fsrPhotons_pterr",1);
        tree->SetBranchStatus("fsrPhotons_eta",1);
        tree->SetBranchStatus("fsrPhotons_phi",1);
        tree->SetBranchStatus("fsrPhotons_lepindex",1);
        tree->SetBranchStatus("met",1);
        tree->SetBranchStatus("me_qqZZ_MCFM",1);
        tree->SetBranchAddress("Run",&Run);
        tree->SetBranchAddress("LumiSect",&LumiSect);
        tree->SetBranchAddress("Event",&Event);
        tree->SetBranchAddress("nVtx",&nVtx);
        tree->SetBranchAddress("nInt",&nInt);
        tree->SetBranchAddress("genWeight",&genWeight);
        tree->SetBranchAddress("crossSection",&crossSection);
        tree->SetBranchAddress("pileupWeight",&pileupWeight);
        tree->SetBranchAddress("passedTrig",&passedTrig);
        tree->SetBranchAddress("triggersPassed",&triggersPassed);       
        tree->SetBranchAddress("lep_tightId", &lep_tightId);
        tree->SetBranchAddress("lep_ecalDriven", &lep_ecalDriven);
        tree->SetBranchAddress("passedFiducialSelection",&passedFiducialSelection);
        tree->SetBranchAddress("lep_id", &lep_id);
        tree->SetBranchAddress("lep_Sip",&lep_Sip);
        tree->SetBranchAddress("lep_dxy",&lep_dxy);  
        tree->SetBranchAddress("lep_dz",&lep_dz);
        tree->SetBranchAddress("lep_pt",&lep_pt);
        tree->SetBranchAddress("lep_pterr",&lep_pterr);
        tree->SetBranchAddress("lep_eta",&lep_eta);
        tree->SetBranchAddress("lep_phi",&lep_phi);
        tree->SetBranchAddress("lep_mass",&lep_mass);
        tree->SetBranchAddress("lep_isoCH",&lep_isoCH);
        tree->SetBranchAddress("lep_RelIso",&lep_RelIso);
        tree->SetBranchAddress("lep_RelIsoNoFSR",&lep_RelIsoNoFSR);
        tree->SetBranchAddress("lepFSR_pt",&lepFSR_pt);
        tree->SetBranchAddress("lepFSR_eta",&lepFSR_eta);
        tree->SetBranchAddress("lepFSR_phi",&lepFSR_phi);
        tree->SetBranchAddress("lepFSR_mass",&lepFSR_mass);
        tree->SetBranchAddress("jet_pt",&jet_pt);
        tree->SetBranchAddress("jet_eta",&jet_eta);
        tree->SetBranchAddress("jet_phi",&jet_phi);
        tree->SetBranchAddress("jet_mass",&jet_mass);
        tree->SetBranchAddress("jet_iscleanH4l",&jet_iscleanH4l);
        tree->SetBranchAddress("jet_QGTagger",&jet_QGTagger);
        tree->SetBranchAddress("jet_csvv2",&jet_csvv2);
        tree->SetBranchAddress("fsrPhotons_pt",&fsrPhotons_pt);
        tree->SetBranchAddress("fsrPhotons_pterr",&fsrPhotons_pterr);
        tree->SetBranchAddress("fsrPhotons_eta",&fsrPhotons_eta);
        tree->SetBranchAddress("fsrPhotons_phi",&fsrPhotons_phi);
        tree->SetBranchAddress("fsrPhotons_lepindex",&fsrPhotons_lepindex);
        tree->SetBranchAddress("met",&met);
        tree->SetBranchAddress("me_qqZZ_MCFM",&me_qqZZ_MCFM);
        // Event Selection

        tree->SetBranchStatus("lep_Hindex",1);
        tree->SetBranchStatus("passedZ4lSelection",1);
        tree->SetBranchStatus("passedFullSelection",1);
        tree->SetBranchStatus("passedZXCRSelection",1);
        tree->SetBranchStatus("nZXCRFailedLeptons",1);
        tree->SetBranchStatus("finalState",1);
        tree->SetBranchStatus("dataMCWeight",1);
        tree->SetBranchStatus("k_qqZZ_qcd_M",1);
        tree->SetBranchStatus("k_qqZZ_ewk",1);
        tree->SetBranchStatus("k_ggZZ",1);
        tree->SetBranchStatus("D_bkg_kin",1);
        tree->SetBranchStatus("njets_pt30_eta4p7",1);
        tree->SetBranchStatus("njets_pt30_eta2p5",1);
        tree->SetBranchStatus("nbjets_pt30_eta4p7",1);
        tree->SetBranchStatus("EventCat",1); 

       
        tree->SetBranchAddress("lep_Hindex",&lep_Hindex);
        tree->SetBranchAddress("passedZ4lSelection",&passedZ4lSelection);
        tree->SetBranchAddress("passedFullSelection",&passedFullSelection);
        tree->SetBranchAddress("passedZXCRSelection",&passedZXCRSelection);
        tree->SetBranchAddress("nZXCRFailedLeptons",&nZXCRFailedLeptons);
        tree->SetBranchAddress("finalState",&finalState);
        tree->SetBranchAddress("dataMCWeight",&dataMCWeight);
        tree->SetBranchAddress("k_qqZZ_qcd_M",&k_qqZZ_qcd_M);
        tree->SetBranchAddress("k_qqZZ_ewk",&k_qqZZ_ewk);
        tree->SetBranchAddress("k_ggZZ",&k_ggZZ);
        tree->SetBranchAddress("D_bkg_kin",&D_bkg_kin);
        tree->SetBranchAddress("njets_pt30_eta4p7",&njets_pt30_eta4p7);
        tree->SetBranchAddress("njets_pt30_eta2p5",&njets_pt30_eta2p5);
        tree->SetBranchAddress("nbjets_pt30_eta4p7",&nbjets_pt30_eta4p7);
        tree->SetBranchAddress("EventCat",&EventCat);

        cout <<"filename "<<filename<<endl;
        cout<<"end loading the tree"<<endl;
        
    }

}

#endif
    
