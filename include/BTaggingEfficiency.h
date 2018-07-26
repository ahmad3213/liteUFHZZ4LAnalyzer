#ifndef BTaggingEfficiency_h
#define BTaggingEfficiency_h

#include "BTagCalibrationStandalone.h"


TFile *f_LWP  = new TFile("include/ZVTo2L2Q_80X_bTaggingLWPEffAnalyzerAK4PF_bTaggingEffMap.root","READ");
TFile *f_MWP = new TFile("include/ZVTo2L2Q_80X_bTaggingMWPEffAnalyzerAK4PF_bTaggingEffMap.root","READ");

namespace BTaggingEff {

double bTaggingEff(double pt, double eta, int partonFlavour){

  TString pID = "udsg";
  if(abs(partonFlavour)==4) pID="c";
  if(abs(partonFlavour)==5) pID="b";

  TH2D* eff = (TH2D*) f_LWP->Get("efficiency_"+pID);
  int nX = eff->GetXaxis()->FindBin(min(pt,670.0));
  int nY = eff->GetYaxis()->FindBin(min(fabs(eta),2.5));

  return eff->GetBinContent(nX,nY);

}

}

#endif
