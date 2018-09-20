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
        if (isData) {
            passTrig = passedTrig;
            /*
            // double ele
            if (strstr((*triggersPassed).c_str(),"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v")) passTrig=true;
            // double mu
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) passTrig=true;
            // single ele
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele25_eta2p1_WPTight_Gsf_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele27_WPTight_Gsf_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele27_eta2p1_WPLoose_Gsf_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele32_eta2p1_WPTight_Gsf_v")) passTrig=true;
            // single mu
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoMu20_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu20_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoMu22_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu22_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoMu22_eta2p1_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu22_eta2p1_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoMu24_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_IsoTkMu24_v")) passTrig=true;
            // multi lepton
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_TripleMu_12_10_5_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v")) passTrig=true;
            else if (strstr((*triggersPassed).c_str(),"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v")) passTrig=true;
            */
        } else {
            passTrig = true;
            sumweight += pileupWeight*genWeight;
        }

        if((*lep_id).size()<4 || (*lep_pt).size()<4) continue;
        unsigned int Nlep = (*lep_id).size();
        if (debug) cout<<Nlep<<" leptons in total"<<endl;


        if (redoEventSelection) {
            // 2 OSSF Pairs
            bool properLep_ID = false; int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
            for(unsigned int i =0; i<Nlep; i++) {
                if((*lep_id)[i]==-13) Nmm = Nmm+1;
                if((*lep_id)[i]==13) Nmp = Nmp+1;
            }
            for(unsigned int i =0; i<Nlep; i++) {
                if((*lep_id)[i]==-11) Nem = Nem+1;
                if((*lep_id)[i]==11) Nep = Nep+1;
            }
            
            if(Nmm>=2 && Nmp>=2) properLep_ID = true; //4mu
            if(Nem>=2 && Nep>=2) properLep_ID = true; //4e
            if(Nmm>0 && Nmp>0 && Nem>0 && Nep>0) properLep_ID = true; //2e2mu
            
            if(!properLep_ID) continue;
            // First, make all Z candidates including any FSR photons
            const double Zmass = 91.1876;
            int n_Zs=0;
            vector<int> Z_lepindex1;
            vector<int> Z_lepindex2;
            vector<float> Z_pt, Z_eta, Z_phi, Z_mass;
            for(unsigned int i=0; i<Nlep; i++){
                for(unsigned int j=i+1; j<Nlep; j++){
                    // same flavor opposite charge
                     if(((*lep_id)[i]+(*lep_id)[j])!=0) continue;
                   // 
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
                    
                    if (debug) cout<<"this Z mass: "<<Z.M()<<" mZ2Low: "<<mZ2Low<<endl;
                    
                    if (Z.M()>0.0) {
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
}
}
}
void SetNewTree(TTree* newtree){

    newtree->Branch("Run",&Run,"Run/l");
    newtree->Branch("Event",&Event,"Event/l");
    newtree->Branch("LumiSect",&LumiSect,"LumiSect/l");
}

