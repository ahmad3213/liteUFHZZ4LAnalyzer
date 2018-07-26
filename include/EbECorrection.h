#ifndef EbECorrection_h
#define EbECorrection_h

// Data
TFile* f_corr_e_1 = new TFile("KinZfitter/DoubleLepton_m2eLUT_m2e_1.root","READ");
TFile* f_corr_e_2 = new TFile("KinZfitter/DoubleLepton_m2eLUT_m2e_2.root","READ");
TFile* f_corr_e_3 = new TFile("KinZfitter/DoubleLepton_m2eLUT_m2e_3.root","READ");
TFile* f_corr_mu = new TFile("KinZfitter/DoubleLepton_m2muLUT_m2mu.root","READ");

TH2F* el_corr_1 = (TH2F*)f_corr_e_1->Get("2e");
TH2F* el_corr_2 = (TH2F*)f_corr_e_2->Get("2e");
TH2F* el_corr_3 = (TH2F*)f_corr_e_3->Get("2e");

TH2F* mu_corr = (TH2F*)f_corr_mu->Get("2mu");

TAxis* x_elpTaxis_1 = el_corr_1->GetXaxis(); TAxis* y_eletaaxis_1 = el_corr_1->GetYaxis();
double maxPtEl_1 = x_elpTaxis_1->GetXmax(); double minPtEl_1 = x_elpTaxis_1->GetXmin();

TAxis* x_elpTaxis_2 = el_corr_2->GetXaxis(); TAxis* y_eletaaxis_2 = el_corr_2->GetYaxis();
double maxPtEl_2 = x_elpTaxis_2->GetXmax(); double minPtEl_2 = x_elpTaxis_2->GetXmin();

TAxis* x_elpTaxis_3 = el_corr_3->GetXaxis(); TAxis* y_eletaaxis_3 = el_corr_3->GetYaxis();
double maxPtEl_3 = x_elpTaxis_3->GetXmax(); double minPtEl_3 = x_elpTaxis_3->GetXmin();

TAxis* x_mupTaxis = mu_corr->GetXaxis(); TAxis* y_muetaaxis = mu_corr->GetYaxis();
double maxPtMu = x_mupTaxis->GetXmax(); double minPtMu = x_mupTaxis->GetXmin();

// MC
TFile* f_corr_e_1_mc = new TFile("KinZfitter/DYJetsToLL_M-50_m2eLUT_m2e_1.root","READ");
TFile* f_corr_e_2_mc = new TFile("KinZfitter/DYJetsToLL_M-50_m2eLUT_m2e_2.root","READ");
TFile* f_corr_e_3_mc = new TFile("KinZfitter/DYJetsToLL_M-50_m2eLUT_m2e_3.root","READ");
TFile* f_corr_mu_mc = new TFile("KinZfitter/DYJetsToLL_M-50_m2muLUT_m2mu.root","READ");

TH2F* el_corr_1_mc = (TH2F*)f_corr_e_1_mc->Get("2e");
TH2F* el_corr_2_mc = (TH2F*)f_corr_e_2_mc->Get("2e");
TH2F* el_corr_3_mc = (TH2F*)f_corr_e_3_mc->Get("2e");

TH2F* mu_corr_mc = (TH2F*)f_corr_mu_mc->Get("2mu");

TAxis* x_elpTaxis_1_mc = el_corr_1_mc->GetXaxis(); TAxis* y_eletaaxis_1_mc = el_corr_1_mc->GetYaxis();
double maxPtEl_1_mc = x_elpTaxis_1_mc->GetXmax(); double minPtEl_1_mc = x_elpTaxis_1_mc->GetXmin();

TAxis* x_elpTaxis_2_mc = el_corr_2_mc->GetXaxis(); TAxis* y_eletaaxis_2_mc = el_corr_2_mc->GetYaxis();
double maxPtEl_2_mc = x_elpTaxis_2_mc->GetXmax(); double minPtEl_2_mc = x_elpTaxis_2_mc->GetXmin();

TAxis* x_elpTaxis_3_mc = el_corr_3_mc->GetXaxis(); TAxis* y_eletaaxis_3_mc = el_corr_3_mc->GetYaxis();
double maxPtEl_3_mc = x_elpTaxis_3_mc->GetXmax(); double minPtEl_3_mc = x_elpTaxis_3_mc->GetXmin();

TAxis* x_mupTaxis_mc = mu_corr_mc->GetXaxis(); TAxis* y_muetaaxis_mc = mu_corr_mc->GetYaxis();
double maxPtMu_mc = x_mupTaxis_mc->GetXmax(); double minPtMu_mc = x_mupTaxis_mc->GetXmin();


float EbeCorrection(float absID, float pt, float pterr, float abseta, bool ecalDriven, bool isData) {

    float corr=1.0;
    if (isData) {
        // Data
        if (absID==13) {
            // Muons
            int xbin = x_mupTaxis->FindBin(pt);
            int ybin = y_muetaaxis->FindBin(abseta);
            if(pt>minPtMu && pt<maxPtMu ){
                corr= mu_corr->GetBinContent(xbin,ybin);
            }
        }
        else if (absID==11) {
            // Electrons
            if (ecalDriven) {
                // Ecal Driven               
                if (abseta < 1) {
                    if (pterr/pt < 0.03) {                        
                        int xbin = x_elpTaxis_1->FindBin(pt);
                        int ybin = y_eletaaxis_1->FindBin(abseta);
                        if(pt > minPtEl_1 && pt < maxPtEl_1 ){
                            corr = el_corr_1->GetBinContent(xbin,ybin);
                        } else {
                            corr = 1.0;
                        }
                    } else {corr = 1.187;}
                } else if (abseta > 1 && abseta < 2.5) {                    
                    if (pterr/pt < 0.07) {                        
                        int xbin = x_elpTaxis_2->FindBin(pt);
                        int ybin = y_eletaaxis_2->FindBin(abseta);
                        if(pt > minPtEl_2 && pt < maxPtEl_2 ){
                            corr = el_corr_2->GetBinContent(xbin,ybin);
                        } else {
                            corr = 1.0;
                        }
                    } else {corr = 0.815;}
                }
            } else {
                // Non-Ecal Driven                
                int xbin = x_elpTaxis_3->FindBin(pt);
                int ybin = y_eletaaxis_3->FindBin(abseta);
                if(pt > minPtEl_3 && pt < maxPtEl_3 ){
                    corr = el_corr_3->GetBinContent(xbin,ybin);
                } else {
                    corr = 1.0;
                }
            }    
        }
    } else {
        // MC
        if (absID==13) {
            // Muons
            int xbin = x_mupTaxis_mc->FindBin(pt);
            int ybin = y_muetaaxis_mc->FindBin(abseta);
            if(pt>minPtMu_mc && pt<maxPtMu_mc ){
                corr= mu_corr_mc->GetBinContent(xbin,ybin);
            }
        }
        else if (absID==11) {
            // Electrons
            if (ecalDriven) {
                // Ecal Driven                
                if (abseta < 1) {
                    if (pterr/pt < 0.03) {                        
                        int xbin = x_elpTaxis_1_mc->FindBin(pt);
                        int ybin = y_eletaaxis_1_mc->FindBin(abseta);
                        if(pt > minPtEl_1_mc && pt < maxPtEl_1_mc ){
                            corr = el_corr_1_mc->GetBinContent(xbin,ybin);
                        } else {
                            corr = 1.0;
                        }
                    } else {corr = 1.224;}
                } else if (abseta > 1 && abseta < 2.5) {                    
                    if (pterr/pt < 0.07) {                        
                        int xbin = x_elpTaxis_2_mc->FindBin(pt);
                        int ybin = y_eletaaxis_2_mc->FindBin(abseta);
                        if(pt > minPtEl_2_mc && pt < maxPtEl_2_mc ){
                            corr = el_corr_2_mc->GetBinContent(xbin,ybin);
                        } else {
                            corr = 1.0;
                        }
                    } else {corr = 0.786;}
                }
            } else {
                // Non-Ecal Driven
                int xbin = x_elpTaxis_3_mc->FindBin(pt);
                int ybin = y_eletaaxis_3_mc->FindBin(abseta);
                if(pt > minPtEl_3_mc && pt < maxPtEl_3_mc ){
                    corr = el_corr_3_mc->GetBinContent(xbin,ybin);
                } else {
                    corr = 1.0;
                }
            }    
        }

    }

    //cout<<"absID: "<<" pt: "<<pt<<" abseta: "<<abseta<<" corr: "<<corr;
    return corr;
    
}

#endif
