#ifndef PileupWeight_h
#define PileupWeight_h
TFile *f_pilup  = new TFile("include/pileup_MC_80x_2016B_Run_271036-274443_mb_69735.root","READ");
TH1D *puweight_dtmc = (TH1D*) f_pilup->Get("puweight_dtmc");
double puweight(int nPV){

 return puweight_dtmc->GetBinContent(nPV);
//HZZ4LPileUp pileUp;
//TFile *f_pileup  = new TFile("include/puWeightsMoriond17_v2.root","READ");
//TH1D *h_pileup = (TH1D*) f_pileup->Get("weights");
// double puweight(int nPV){
//        pileupWeight = pileUp.getPUWeight(h_pileup,nPV);

}

#endif
// return puweight_dtmc->GetBinContent(nPV);
//   return pileupWeight;
