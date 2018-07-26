#ifndef LeptonEfficiency_h
#define LeptonEfficiency_h

#include "TEfficiency.h"

TFile *f_data_se  = new TFile("include/eff_data_2d_se.root","READ");
TFile *f_data_de1 = new TFile("include/eff_data_2d_de1.root","READ");
TFile *f_data_de2 = new TFile("include/eff_data_2d_de2.root","READ");
TFile *f_data_sm  = new TFile("include/eff_data_2d_sm.root","READ");
TFile *f_data_dm1 = new TFile("include/eff_data_2d_dm1.root","READ");
TFile *f_data_dm2 = new TFile("include/eff_data_2d_dm2.root","READ");

TFile *f_mc_se  = new TFile("include/eff_mc_2d_se.root","READ");
TFile *f_mc_de1 = new TFile("include/eff_mc_2d_de1.root","READ");
TFile *f_mc_de2 = new TFile("include/eff_mc_2d_de2.root","READ");
TFile *f_mc_sm  = new TFile("include/eff_mc_2d_sm.root","READ");
TFile *f_mc_dm1 = new TFile("include/eff_mc_2d_dm1.root","READ");
TFile *f_mc_dm2 = new TFile("include/eff_mc_2d_dm2.root","READ");

TCanvas *c_data_se = (TCanvas*)f_data_se->Get("c1");
TEfficiency *eff_data_se = (TEfficiency*)c_data_se->GetPrimitive("den_se_2d_clone");
TCanvas *c_data_de1 = (TCanvas*)f_data_de1->Get("c1");
TEfficiency *eff_data_de1 = (TEfficiency*)c_data_de1->GetPrimitive("den_de1_2d_clone");
TCanvas *c_data_de2 = (TCanvas*)f_data_de2->Get("c1");
TEfficiency *eff_data_de2 = (TEfficiency*)c_data_de2->GetPrimitive("den_de2_2d_clone");
TCanvas *c_data_sm = (TCanvas*)f_data_sm->Get("c1");
TEfficiency *eff_data_sm = (TEfficiency*)c_data_sm->GetPrimitive("den_sm_2d_clone");
TCanvas *c_data_dm1 = (TCanvas*)f_data_dm1->Get("c1");
TEfficiency *eff_data_dm1 = (TEfficiency*)c_data_dm1->GetPrimitive("den_dm1_2d_clone");
TCanvas *c_data_dm2 = (TCanvas*)f_data_dm2->Get("c1");
TEfficiency *eff_data_dm2 = (TEfficiency*)c_data_dm2->GetPrimitive("den_dm2_2d_clone");

TCanvas *c_mc_se = (TCanvas*)f_mc_se->Get("c2");
TEfficiency *eff_mc_se = (TEfficiency*)c_mc_se->GetPrimitive("den_se_2d_clone");
TCanvas *c_mc_de1 = (TCanvas*)f_mc_de1->Get("c2");
TEfficiency *eff_mc_de1 = (TEfficiency*)c_mc_de1->GetPrimitive("den_de1_2d_clone");
TCanvas *c_mc_de2 = (TCanvas*)f_mc_de2->Get("c2");
TEfficiency *eff_mc_de2 = (TEfficiency*)c_mc_de2->GetPrimitive("den_de2_2d_clone");
TCanvas *c_mc_sm = (TCanvas*)f_mc_sm->Get("c2");
TEfficiency *eff_mc_sm = (TEfficiency*)c_mc_sm->GetPrimitive("den_sm_2d_clone");
TCanvas *c_mc_dm1 = (TCanvas*)f_mc_dm1->Get("c2");
TEfficiency *eff_mc_dm1 = (TEfficiency*)c_mc_dm1->GetPrimitive("den_dm1_2d_clone");
TCanvas *c_mc_dm2 = (TCanvas*)f_mc_dm2->Get("c2");
TEfficiency *eff_mc_dm2 = (TEfficiency*)c_mc_dm2->GetPrimitive("den_dm2_2d_clone");

TFile *f_egrecosf = new TFile("include/eleRECO.txt.egamma_SF2D.root","READ");
TH2D *egrecosf = (TH2D*)f_egrecosf->Get("EGamma_SF2D");

#endif
