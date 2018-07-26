#ifndef AnalysisTree_h
#define AnalysisTree_h

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

  // input tree variables
  std::string *triggersPassed;
  ULong64_t Run, LumiSect, Event;

  float dataMCWeight, genWeight;
  float sumW;
  int nInt; //nPV

  std::vector<float>* GENZ_pt;

  std::vector<float>* lep_mass;
  std::vector<float> *lep_pt; std::vector<float> *lep_eta; std::vector<float> *lep_phi;
  std::vector<float>* lepFSR_mass;
  std::vector<float> *lepFSR_pt; std::vector<float> *lepFSR_eta; std::vector<float> *lepFSR_phi;

  std::vector<int> *lep_tightId;
  std::vector<int> *lep_tightIdSUS;
  std::vector<int> *lep_tightIdHiPt;
  std::vector<int> *lep_id;
  std::vector<int> *lep_Sip;
  std::vector<float> *lep_dxy
  std::vector<float> *lep_dz
  std::vector<float> *lep_RelIso;
  std::vector<float> *lep_RelIsoNoFSR;
  std::vector<float> *lep_pterr;
  std::vector<float> *lep_dataMC;
  std::vector<string> *lep_filtersMatched;

  std::vector<float> *fsrPhotons_pt; std::vector<float> *fsrPhotons_eta; std::vector<float> *fsrPhotons_phi;
  std::vector<int> *fsrPhotons_lepindex;

  /////

  std::vector<float> *jet_mass;
  std::vector<float> *jet_pt; std::vector<float> *jet_eta; std::vector<float> *jet_phi;
  std::vector<float> *jet_csvv2;
  std::vector<float> *jet_relpterr; std::vector<float> *jet_phierr;
  std::vector<float> *jet_QGTagger;
  std::vector<int> *jet_partonFlavour;
  std::vector<int> *jet_iscleanH4l;

  /////

  std::vector<float> *mergedjet_mass;
  std::vector<float> *mergedjet_pt; std::vector<float> *mergedjet_eta; std::vector<float> *mergedjet_phi;
  std::vector<float> *mergedjet_prunedmass; std::vector<float> *mergedjet_softdropmass;
  std::vector<float> *mergedjet_tau1; std::vector<float> *mergedjet_tau2;
  std::vector<float> *mergedjet_L1; // correction from current JEC to L1

  std::vector<int> *mergedjet_iscleanH4l;

  // softdrop subjet
  vector<int>* mergedjet_nsubjet;
  vector<vector<float> >* mergedjet_subjet_pt;
  vector<vector<float> >* mergedjet_subjet_eta;
  vector<vector<float> >* mergedjet_subjet_phi;
  vector<vector<float> >* mergedjet_subjet_mass;
  vector<vector<float> >* mergedjet_subjet_btag;
  vector<vector<int> >* mergedjet_subjet_partonFlavour;

  ////
  float met, ht;
  int lheNb, lheNj, nGenStatus2bHad;

  ////
  // out
  
  // lepton, trigger eff.....
  bool passTrig;

  double GENzPt;
  int nGENz;

  int VlepId;
  int L1, L2;
  int nl;
  double mlljj, mll, mjj, pT_lljj;
  double pT_ll;
  double pTl1, etal1, phil1, ml1, pTerrl1;
  double pTl2, etal2, phil2, ml2, pTerrl2;

  double dataeff, mceff;
  bool hlt_se_L1, hlt_sm_L1, hlt_de1_L1, hlt_de2_L1, hlt_dm1_L1, hlt_dm2_L1;
  bool hlt_se_L2, hlt_sm_L2, hlt_de1_L2, hlt_de2_L2, hlt_dm1_L2, hlt_dm2_L2;

  // merged jet
  int nJ;
  double mllJ, mJ_prune, mJ_softdrop, mJ_jec, pT_llJ;
  double tau21, mJ_L1;
  double pT_J, pT_jj;
  int nsubjet;
  int partonFlavour_subjet1, partonFlavour_subjet2;
  double btag_subjet1, btag_subjet2;
  double btagSF_subjet1, btagSF_subjet2;
  double btagSFup_subjet1, btagSFup_subjet2, btagSFdn_subjet1, btagSFdn_subjet2;   
  double btagEff_subjet1, btagEff_subjet2;
  double btagWeight_J, btagWeightUp_J, btagWeightDn_J; int btagCat_J;
   
  // resolved jets
  //
  int nj;
  int nbjl, nbjm, nbjt; 
  double pTj1, etaj1, phij1, mj1, pTerrj1, phierrj1;
  double pTj2, etaj2, phij2, mj2, pTerrj2, phierrj2;
  double qgTag1, qgTag2;
  double btagj1, btagj2, btagave;
  int partonFlavour_jet1, partonFlavour_jet2;
  double btagSFj1, btagSFj2;
  double btagSFupj1, btagSFupj2, btagSFdnj1, btagSFdnj2;
  double btagEffj1, btagEffj2;
  double btagWeight_jj, btagWeightUp_jj, btagWeightDn_jj; int btagCat_jj;

  double pTj1Refit, etaj1Refit, phij1Refit, mj1Refit;
  double pTj2Refit, etaj2Refit, phij2Refit, mj2Refit;  
  double mjjRefit, mlljjRefit;

  double kinFitchi2; 

  float cosThetaStar, cosTheta1, cosTheta2;
  float Phi, Phi1; 

  // VBF
  int njVBF; double detajjVBF, mjjVBF;
  double pTVBFj1, etaVBFj1, phiVBFj1, mVBFj1;
  double pTVBFj2, etaVBFj2, phiVBFj2, mVBFj2;
  double qgTagVBF1, qgTagVBF2;
  double eta1x2VBF;

  /////////
  double weight;
  float xsec,crossSection;

  bool findXbb, findZjj, findZJ;

  bool passedTrig, passedFullSelection;

  double ZPtWeight, Kzpt;

  /////////////

  bool isFromDavid;
  bool isHiPtId;
  bool isData;

namespace AnalysisTree {

  void setAddresses(TTree* tree, TString filename){

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("passedFullSelection",1);
    tree->SetBranchStatus("passedTrig",1);
    tree->SetBranchStatus("triggersPassed",1);
    tree->SetBranchStatus("GENZ_pt",1);
    tree->SetBranchStatus("lep_id",1);
    tree->SetBranchStatus("lep_tightId",1);
    tree->SetBranchStatus("lep_tightIdSUS",1);
    if(isHiPtId) tree->SetBranchStatus("lep_tightIdHiPt",1);
    tree->SetBranchStatus("lep_Sip",1);
    tree->SetBranchStatus("lep_dxy",1);
     tree->SetBranchStatus("lep_dz",1); 

    tree->SetBranchAddress("lep_tightId", &lep_tightId);
    tree->SetBranchAddress("lep_tightIdSUS", &lep_tightIdSUS);
    if(isHiPtId) tree->SetBranchAddress("lep_tightIdHiPt", &lep_tightIdHiPt);
    tree->SetBranchAddress("lep_Sip", &lep_Sip);
    tree->SetBranchAddress("lep_dxy", &lep_dxy);
    tree->SetBranchAddress("lep_dz", &lep_dz);
    tree->SetBranchAddress("lep_id", &lep_id);
    tree->SetBranchAddress("passedFullSelection",&passedFullSelection);
    tree->SetBranchAddress("passedTrig",&passedTrig);
    tree->SetBranchAddress("triggersPassed",&triggersPassed);
    tree->SetBranchAddress("GENZ_pt",&GENZ_pt);

    tree->SetBranchStatus("Run",1);
    tree->SetBranchStatus("LumiSect",1);
    tree->SetBranchStatus("Event",1);
    tree->SetBranchAddress("Run",&Run);
    tree->SetBranchAddress("LumiSect",&LumiSect);
    tree->SetBranchAddress("Event",&Event);

    tree->SetBranchStatus("jet_pt",1);
    tree->SetBranchStatus("jet_eta",1);
    tree->SetBranchStatus("jet_phi",1);
    tree->SetBranchStatus("jet_mass",1);
    tree->SetBranchStatus("jet_relpterr",1);
    tree->SetBranchStatus("jet_phierr",1);
    tree->SetBranchStatus("jet_QGTagger",1);
    tree->SetBranchStatus("jet_csvv2",1);
    tree->SetBranchStatus("jet_partonFlavour",1);

    tree->SetBranchAddress("jet_pt",&jet_pt);
    tree->SetBranchAddress("jet_eta",&jet_eta);
    tree->SetBranchAddress("jet_phi",&jet_phi);
    tree->SetBranchAddress("jet_mass",&jet_mass);
    tree->SetBranchAddress("jet_relpterr",&jet_relpterr);
    tree->SetBranchAddress("jet_phierr",&jet_phierr);
    tree->SetBranchAddress("jet_QGTagger",&jet_QGTagger);
    tree->SetBranchAddress("jet_csvv2",&jet_csvv2);
    tree->SetBranchAddress("jet_partonFlavour",&jet_partonFlavour);

    tree->SetBranchStatus("mergedjet_pt",1);
    tree->SetBranchStatus("mergedjet_eta",1);
    tree->SetBranchStatus("mergedjet_phi",1);
    tree->SetBranchStatus("mergedjet_mass",1);
    tree->SetBranchStatus("mergedjet_prunedmass",1);
    tree->SetBranchStatus("mergedjet_softdropmass",1);
    tree->SetBranchStatus("mergedjet_tau1",1);
    tree->SetBranchStatus("mergedjet_tau2",1);
    tree->SetBranchStatus("mergedjet_L1",1);
    tree->SetBranchStatus("mergedjet_subjet_pt",1);
    tree->SetBranchStatus("mergedjet_subjet_eta",1);
    tree->SetBranchStatus("mergedjet_subjet_phi",1);
    tree->SetBranchStatus("mergedjet_subjet_mass",1);
    tree->SetBranchStatus("mergedjet_subjet_btag",1);
    tree->SetBranchStatus("mergedjet_nsubjet",1);
    tree->SetBranchStatus("mergedjet_subjet_partonFlavour",1);

    tree->SetBranchAddress("mergedjet_pt",&mergedjet_pt);
    tree->SetBranchAddress("mergedjet_eta",&mergedjet_eta);
    tree->SetBranchAddress("mergedjet_phi",&mergedjet_phi);
    tree->SetBranchAddress("mergedjet_mass",&mergedjet_mass);
    tree->SetBranchAddress("mergedjet_tau1",&mergedjet_tau1);
    tree->SetBranchAddress("mergedjet_tau2",&mergedjet_tau2);
    tree->SetBranchAddress("mergedjet_L1",&mergedjet_L1);
    tree->SetBranchAddress("mergedjet_prunedmass",&mergedjet_prunedmass);
    tree->SetBranchAddress("mergedjet_softdropmass",&mergedjet_softdropmass);

    tree->SetBranchAddress("mergedjet_nsubjet",&mergedjet_nsubjet);
    tree->SetBranchAddress("mergedjet_subjet_pt",&mergedjet_subjet_pt);
    tree->SetBranchAddress("mergedjet_subjet_eta",&mergedjet_subjet_eta);
    tree->SetBranchAddress("mergedjet_subjet_phi",&mergedjet_subjet_phi);
    tree->SetBranchAddress("mergedjet_subjet_mass",&mergedjet_subjet_mass);
    tree->SetBranchAddress("mergedjet_subjet_btag",&mergedjet_subjet_btag);
    tree->SetBranchAddress("mergedjet_subjet_partonFlavour",&mergedjet_subjet_partonFlavour);

    tree->SetBranchStatus("lep_pt",1);
    tree->SetBranchStatus("lep_pterr",1);
    tree->SetBranchStatus("lep_eta",1);
    tree->SetBranchStatus("lep_phi",1);
    tree->SetBranchStatus("lep_mass",1);
    tree->SetBranchStatus("lep_RelIso",1);
    tree->SetBranchStatus("lep_RelIsoNoFSR",1);
    tree->SetBranchStatus("lepFSR_pt",1);
    tree->SetBranchStatus("lepFSR_eta",1);
    tree->SetBranchStatus("lepFSR_phi",1);
    tree->SetBranchStatus("lepFSR_mass",1);
    tree->SetBranchAddress("lep_pt",&lep_pt);
    tree->SetBranchAddress("lep_pterr",&lep_pterr);
    tree->SetBranchAddress("lep_eta",&lep_eta);
    tree->SetBranchAddress("lep_phi",&lep_phi);
    tree->SetBranchAddress("lep_mass",&lep_mass);
    tree->SetBranchAddress("lep_RelIso",&lep_RelIso);
    tree->SetBranchAddress("lep_RelIsoNoFSR",&lep_RelIsoNoFSR);
    tree->SetBranchAddress("lepFSR_pt",&lepFSR_pt);
    tree->SetBranchAddress("lepFSR_eta",&lepFSR_eta);
    tree->SetBranchAddress("lepFSR_phi",&lepFSR_phi);
    tree->SetBranchAddress("lepFSR_mass",&lepFSR_mass);

    tree->SetBranchStatus("fsrPhotons_pt",1);
    tree->SetBranchStatus("fsrPhotons_eta",1);
    tree->SetBranchStatus("fsrPhotons_phi",1);
    tree->SetBranchStatus("fsrPhotons_lepindex",1);
    tree->SetBranchAddress("fsrPhotons_pt",&fsrPhotons_pt);
    tree->SetBranchAddress("fsrPhotons_eta",&fsrPhotons_eta);
    tree->SetBranchAddress("fsrPhotons_phi",&fsrPhotons_phi);
    tree->SetBranchAddress("fsrPhotons_lepindex",&fsrPhotons_lepindex);

    tree->SetBranchStatus("dataMCWeight",1);
    tree->SetBranchStatus("genWeight",1);
    tree->SetBranchStatus("crossSection",1);
    tree->SetBranchStatus("lep_dataMC",1);
    tree->SetBranchStatus("lep_filtersMatched",1);
    tree->SetBranchStatus("nInt",1);

    tree->SetBranchAddress("lep_dataMC",&lep_dataMC);
    tree->SetBranchAddress("lep_filtersMatched",&lep_filtersMatched);
    tree->SetBranchAddress("dataMCWeight", &dataMCWeight);
    tree->SetBranchAddress("genWeight", &genWeight);
    tree->SetBranchAddress("crossSection", &crossSection);
    tree->SetBranchAddress("nInt",&nInt);

    tree->SetBranchStatus("met",1);
    tree->SetBranchAddress("met", &met);

    cout <<"filename "<<filename<<endl;
    cout<<"strstr(filename,DY) "<<strstr(filename,"DY")<<endl;

    if (strstr(filename,"DY")) {
        cout <<"here"<<endl;
        tree->SetBranchStatus("lheNb",1);
        tree->SetBranchAddress("lheNb", &lheNb);
        tree->SetBranchStatus("lheNj",1);
        tree->SetBranchAddress("lheNj", &lheNj);
        tree->SetBranchStatus("nGenStatus2bHad",1);
        tree->SetBranchAddress("nGenStatus2bHad", &nGenStatus2bHad);
    }

    cout<<"end loading the tree"<<endl;


 }




 void setCMS(TH1F* Graph1)
 {

    Graph1->GetXaxis()->SetLabelFont(42);
    Graph1->GetXaxis()->SetLabelOffset(0.003);

    Graph1->GetXaxis()->SetLabelSize(0.035);

    Graph1->GetXaxis()->SetTitleSize(0.04);
    Graph1->GetXaxis()->SetTitleOffset(0.8);
    Graph1->GetXaxis()->SetTitleFont(42);

    Graph1->GetYaxis()->SetLabelFont(42);

    Graph1->GetYaxis()->SetLabelOffset(0.004);
    Graph1->GetYaxis()->SetLabelSize(0.035);
    Graph1->GetYaxis()->SetTitleSize(0.04);
    Graph1->GetYaxis()->SetTitleOffset(1.1);
    Graph1->GetYaxis()->SetTitleFont(42);

 }

 void setCMS(TGraph* Graph1)
 {

    Graph1->GetXaxis()->SetLabelFont(42);
    Graph1->GetXaxis()->SetLabelOffset(0.003);

    Graph1->GetXaxis()->SetLabelSize(0.035);

    Graph1->GetXaxis()->SetTitleSize(0.04);
    Graph1->GetXaxis()->SetTitleOffset(0.8);
    Graph1->GetXaxis()->SetTitleFont(42);

    Graph1->GetYaxis()->SetLabelFont(42);

    Graph1->GetYaxis()->SetLabelOffset(0.004);
    Graph1->GetYaxis()->SetLabelSize(0.035);
    Graph1->GetYaxis()->SetTitleSize(0.04);
    Graph1->GetYaxis()->SetTitleOffset(1.0);
    Graph1->GetYaxis()->SetTitleFont(42);

 }

 void setCMS(TH2F* Graph1)
 {

    Graph1->GetXaxis()->SetLabelFont(42);
    Graph1->GetXaxis()->SetLabelOffset(0.003);

    Graph1->GetXaxis()->SetLabelSize(0.035);

    Graph1->GetXaxis()->SetTitleSize(0.04);
    Graph1->GetXaxis()->SetTitleOffset(0.8);
    Graph1->GetXaxis()->SetTitleFont(42);

    Graph1->GetYaxis()->SetLabelFont(42);

    Graph1->GetYaxis()->SetLabelOffset(0.004);
    Graph1->GetYaxis()->SetLabelSize(0.035);
    Graph1->GetYaxis()->SetTitleSize(0.04);
    Graph1->GetYaxis()->SetTitleOffset(1.1);
    Graph1->GetYaxis()->SetTitleFont(42);
 }

}

#endif
