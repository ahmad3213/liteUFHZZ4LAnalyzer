#include "ZZ4LConfig.h" 
#include "deltaPhi.h"
#include "ZZ4LAnalysisTree.h"
#include "LeptonEfficiency.h"
#include "PileupWeight.h"
#include "EbECorrection.h"
#include "Helper.h"
using namespace std;
/////////////////////

void SetNewTree(TTree* newtree);
void ReadTree(TTree* tree, TTree* & newtree, TString filename);
TString filename;
bool debug;


int main(int argc, char *argv[])
{    
     
  debug = false;     

  if(argc > 6)  {
      cout<<argv[0]<<" filename "<<argv[1]<<" outfile "<<argv[2]<<" isData "<<argv[3]<<endl;
      return -1;
    }


  /////////////////////

  filename = argv[1];
  TString outfilename = argv[2];
 
  if(atof(argv[3])>0) { isData = true; }
  else {isData = false; }
  if(atof(argv[4])>0) { 
      job=strtol(argv[4], NULL, 10); njobs=strtol(argv[5], NULL, 10);
      cout<<"job "<<job<<" of "<<njobs<<endl;
      outfilename+="_";
      outfilename+=argv[4];
  }
  // reading Input Ntuple file
  TFile* infile = TFile::Open(filename+".root");
  TTree* tree;
  tree = (TTree*) infile->Get("Ana/passedEvents");
  if(!tree) tree = (TTree*) infile->Get("passedEvents");
  if(!tree) tree = (TTree*) infile->Get("selectedEvents");
  if(!tree) { cout<<"ERROR could not find the tree for "<<filename<<endl; return -1;}

  TH1F* sumWeights = (TH1F*) infile->Get("Ana/sumWeights");
  sumW = 1.0;
  if(sumWeights) sumW = sumWeights->GetBinContent(1); 

  cout<<"sumW is "<<sumW<<endl;

  // read tree     

  TString name = outfilename;
  TFile* tmpFile =  new TFile(name+".root","RECREATE");
  TTree* newtree = new TTree("passedEvents","passedEvents");

  if(debug)cout<<"start setting new tree "<<endl;

  SetNewTree(newtree);

  if(debug)cout<<"start reading tree "<<endl;

  ReadTree(tree, newtree, filename);
  if(debug)cout<<"end reading tree"<<endl;

  tmpFile->cd();

  newtree->Write("passedEvents",TObject::kOverwrite);
  tmpFile->Close(); 
}

void ReadTree(TTree* tree, TTree* & newtree, TString filename)
    {
    ZZ4LAnalysisTree::setAddresses(tree, filename);
    float npass = 0.0;
    float sumweight = 0.0;
    if(debug) cout<<"start looping"<<endl;
    int firstevt=0; int lastevt=tree->GetEntries();
    if (job>0) {
        firstevt = tree->GetEntries()*(job-1)/njobs;
        lastevt = tree->GetEntries()*(job)/njobs-1;
    }
    for(int evt=0; evt < tree->GetEntries(); evt++) { //event loop
        if (evt<firstevt) continue;
        if (evt>lastevt) continue;
        if(evt%1000==0) cout<<"Event "<<evt<<"/"<<tree->GetEntries()<<endl;
        tree->GetEntry(evt);
        passTrig=false;
        // single mu
        if (strstr((*triggersPassed).c_str(),"HLT_IsoMu20_v")) passTrig=true;
        else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu20_v")) passTrig=true;
        else if (strstr((*triggersPassed).c_str(),"HLT_IsoMu22_v")) passTrig=true;
        else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu22_v")) passTrig=true;
        else if (strstr((*triggersPassed).c_str(),"HLT_IsoMu22_eta2p1_v")) passTrig=true;
        else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu22_eta2p1_v")) passTrig=true;
        else if (strstr((*triggersPassed).c_str(),"HLT_IsoMu24_v")) passTrig=true;
        else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu24_v")) passTrig=true;
        sumweight += pileupWeight*genWeight;
        if((*lep_id).size()<2 || (*lep_pt).size()<2) continue;
        unsigned int Nlep = (*lep_id).size();
        if (debug) cout<<Nlep<<" leptons in total"<<endl;
        if (redoEventSelection) {
        // First, make all Z candidates including any FSR photons
        const double Zmass = 91.1876;
        int n_Zs=0;
        vector<int> Z_lepindex1;
        vector<int> Z_lepindex2;
        vector<float> Z_pt, Z_eta, Z_phi, Z_mass;
        for(unsigned int i=0; i<Nlep; i++){
            for(unsigned int j=i+1; j<Nlep; j++){
                if(((*lep_id)[i]+(*lep_id)[j])!=0) continue; // same flavor opposite charge
                    TLorentzVector li, lj;
                    li.SetPtEtaPhiM((*lep_pt)[i],(*lep_eta)[i],(*lep_phi)[i],(*lep_mass)[i]);
                    lj.SetPtEtaPhiM((*lep_pt)[j],(*lep_eta)[j],(*lep_phi)[j],(*lep_mass)[j]);
                    
                    TLorentzVector lifsr, ljfsr;
                    lifsr.SetPtEtaPhiM((*lepFSR_pt)[i],(*lepFSR_eta)[i],(*lepFSR_phi)[i],(*lepFSR_mass)[i]);
                    ljfsr.SetPtEtaPhiM((*lepFSR_pt)[j],(*lepFSR_eta)[j],(*lepFSR_phi)[j],(*lepFSR_mass)[j]);
                    
                    TLorentzVector liljfsr = lifsr+ljfsr;
                    
                    if (debug) {
                        cout<<"OSSF pair: i="<<i<<" id1="<<(*lep_id)[i]<<" j="<<j<<" id2="<<(*lep_id)[j]<<" pt1: "
                            <<lifsr.Pt()<<" pt2: "<<ljfsr.Pt()<<" M: "<<liljfsr.M()<<endl;    
                    }
                    TLorentzVector Z, Z_noFSR;
                    Z = lifsr+ljfsr;
                    Z_noFSR = li+lj;
                    if (debug) cout<<"this Z mass: "<<Z.M()<<endl;
                    if (Z.M()>50.0 && Z.M()>120.0) {
                        n_Zs++;
                        Z_pt.push_back(Z.Pt());
                        Z_eta.push_back(Z.Eta());
                        Z_phi.push_back(Z.Phi());
                        Z_mass.push_back(Z.M());
                        Z_lepindex1.push_back(i);
                        Z_lepindex2.push_back(j);
                        if (debug) cout<<" add Z_lepindex1: "<<i<<" Z_lepindex2: "<<j<<endl;
                   } 
                } // lep i
            } // lep j
            if(debug) cout<<"fill tree"<<endl;
            if(debug) cout<<endl;
            newtree->Fill();
    cout<<"sumweight: "<<sumweight<<endl;
    cout<<"npass: "<<npass<<endl;
} // end redo selection
} // end event loop 
} // end read tree function
void SetNewTree(TTree* newtree){

    newtree->Branch("Run",&Run,"Run/l");
    newtree->Branch("Event",&Event,"Event/l");
    newtree->Branch("LumiSect",&LumiSect,"LumiSect/l");
}

