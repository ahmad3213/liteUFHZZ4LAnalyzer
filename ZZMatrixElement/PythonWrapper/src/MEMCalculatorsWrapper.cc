#include "ZZMatrixElement/PythonWrapper/interface/MEMCalculatorsWrapper.h"
#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h"


MEMCalculatorsWrapper::MEMCalculatorsWrapper(double collisionEnergy, double sKD_mass) {
    mem_ = new MEMs(collisionEnergy,sKD_mass);
}

MEMCalculatorsWrapper::~MEMCalculatorsWrapper() {
    if(mem_ !=0) delete mem_;
}

MEMCalculatorsWrapper::Angles 
MEMCalculatorsWrapper::computeAngles(TLorentzVector Z1_lept1, int Z1_lept1Id,
		   TLorentzVector Z1_lept2, int Z1_lept2Id,
		   TLorentzVector Z2_lept1, int Z2_lept1Id,
		   TLorentzVector Z2_lept2, int Z2_lept2Id) {
    Angles ret;
    TUtil::computeAngles(Z1_lept1,Z1_lept1Id,
                        Z1_lept2,Z1_lept2Id,
                        Z2_lept1,Z2_lept1Id,
                        Z2_lept2,Z2_lept2Id,
                    ret.costhetastar,ret.costheta1,ret.costheta2,ret.phi,ret.phistar1);
    return ret;
  }

MEMCalculatorsWrapper::Angles 
MEMCalculatorsWrapper::computeAngles(const math::XYZTLorentzVector & Z1_lept1, int Z1_lept1Id,
		   const math::XYZTLorentzVector & Z1_lept2, int Z1_lept2Id,
		   const math::XYZTLorentzVector & Z2_lept1, int Z2_lept1Id,
		   const math::XYZTLorentzVector & Z2_lept2, int Z2_lept2Id) {
    return computeAngles(TLorentzVector(Z1_lept1.Px(),Z1_lept1.Py(),Z1_lept1.Pz(),Z1_lept1.E()), Z1_lept1Id,
                         TLorentzVector(Z1_lept2.Px(),Z1_lept2.Py(),Z1_lept2.Pz(),Z1_lept2.E()), Z1_lept2Id,
                         TLorentzVector(Z2_lept1.Px(),Z2_lept1.Py(),Z2_lept1.Pz(),Z2_lept1.E()), Z2_lept1Id,
                         TLorentzVector(Z2_lept2.Px(),Z2_lept2.Py(),Z2_lept2.Pz(),Z2_lept2.E()), Z2_lept2Id);
  }


void  
MEMCalculatorsWrapper::computeAll(TLorentzVector Z1_lept1, int Z1_lept1Id,
		   TLorentzVector Z1_lept2, int Z1_lept2Id,
		   TLorentzVector Z2_lept1, int Z2_lept1Id,
		   TLorentzVector Z2_lept2, int Z2_lept2Id) {

    std::vector<TLorentzVector> ps;
    ps.push_back(Z1_lept1);
    ps.push_back(Z1_lept2);
    ps.push_back(Z2_lept1);
    ps.push_back(Z2_lept2);


    std::vector<int> id;
    id.push_back(Z1_lept1Id);
    id.push_back(Z1_lept2Id);
    id.push_back(Z2_lept1Id);
    id.push_back(Z2_lept2Id);

    mem_->computeMEs(ps,id);

    //Now the SuperKD part
    pm4l_sig_=0.0;
    pm4l_bkg_=0.0;
    mem_->computePm4l(ps, id, MEMNames::kNone, pm4l_sig_, pm4l_bkg_);
}

std::vector<std::pair<std::string,float>>  
MEMCalculatorsWrapper::computeNew(
        const math::XYZTLorentzVector & Z1_lept1, int Z1_lept1Id,
        const math::XYZTLorentzVector & Z1_lept2, int Z1_lept2Id,
        const math::XYZTLorentzVector & Z2_lept1, int Z2_lept1Id,
        const math::XYZTLorentzVector & Z2_lept2, int Z2_lept2Id,
        const std::vector<math::XYZTLorentzVector> & jets)
{
    std::vector<TLorentzVector> partP;
    partP.emplace_back(Z1_lept1.Px(),Z1_lept1.Py(),Z1_lept1.Pz(),Z1_lept1.E());
    partP.emplace_back(Z1_lept2.Px(),Z1_lept2.Py(),Z1_lept2.Pz(),Z1_lept2.E());
    partP.emplace_back(Z2_lept1.Px(),Z2_lept1.Py(),Z2_lept1.Pz(),Z2_lept1.E());
    partP.emplace_back(Z2_lept2.Px(),Z2_lept2.Py(),Z2_lept2.Pz(),Z2_lept2.E());

    std::vector<int> partId;
    partId.push_back(Z1_lept1Id);
    partId.push_back(Z1_lept2Id);
    partId.push_back(Z2_lept1Id);
    partId.push_back(Z2_lept2Id);

    double p0plus_VAJHU, bkg_VAMCFM, p0plus_m4l, bkg_m4l;
    double p0minus_VAJHU, Dgg10_VAMCFM;

    mem_->computeME(MEMNames::kSMHiggs, MEMNames::kJHUGen, partP, partId, p0plus_VAJHU); // Calculation of SM gg->H->4l JHUGen ME
    mem_->computeME(MEMNames::k0minus, MEMNames::kJHUGen, partP, partId, p0minus_VAJHU); // Calculation of PS (0-, fa3=1) gg->H->4l JHUGen ME 
    mem_->computeME(MEMNames::kggHZZ_10, MEMNames::kMCFM, partP, partId, Dgg10_VAMCFM); // Direct calculation of Dgg (D^kin for off-shell) from MCFM MEs
    mem_->computeME(MEMNames::kqqZZ, MEMNames::kMCFM, partP, partId, bkg_VAMCFM); // qq->4l background calculation from MCFM
    mem_->computePm4l(partP,partId, MEMNames::kNone, p0plus_m4l, bkg_m4l); // m4l probabilities for signal and background, nominal resolution

    double D_bkg_kin = p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM ); // D^kin_bkg
    double D_bkg = p0plus_VAJHU * p0plus_m4l / ( p0plus_VAJHU * p0plus_m4l + bkg_VAMCFM * bkg_m4l ); // D^kin including superMELA
    double D_g4 = p0plus_VAJHU / ( p0plus_VAJHU + p0minus_VAJHU ); // D_0-
    std::vector<std::pair<std::string,float>> ret;
    ret.emplace_back("D_bkg^kin",D_bkg_kin);
    ret.emplace_back("D_bkg",D_bkg);
    ret.emplace_back("D_gg",Dgg10_VAMCFM);
    ret.emplace_back("D_0-", D_g4);

    if (jets.size() >= 2) {
        std::vector<TLorentzVector> partPprod(partP);
        std::vector<int>            partIdprod(partId);
        for (const auto &p4 : jets) {
            partPprod.emplace_back(p4.Px(), p4.Py(), p4.Pz(), p4.E());
            partIdprod.push_back(0);
            if (partPprod.size() == 6) break;
        }
        double phjj_VAJHU, pvbf_VAJHU;
        mem_->computeME(MEMNames::kJJ_SMHiggs_GG, MEMNames::kJHUGen,  partPprod, partIdprod, phjj_VAJHU); // SM gg->H+2j
        mem_->computeME(MEMNames::kJJ_SMHiggs_VBF, MEMNames::kJHUGen, partPprod, partIdprod, pvbf_VAJHU);  // SM VBF->H
        double Djet_VAJHU = pvbf_VAJHU / ( pvbf_VAJHU + phjj_VAJHU ); // D^VBF_HJJ
        ret.emplace_back("D_HJJ^VBF", Djet_VAJHU);
    } else {
        ret.emplace_back("D_HJJ^VBF", -1.0);
    }

    return ret;
}


float
MEMCalculatorsWrapper:: getKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_qqZZ;
    mem_->computeKD(kSMHiggs, kJHUGen, kqqZZ, kMCFM, &MEMs::probRatio, KD, ME_ggHiggs, ME_qqZZ);
    return KD;
}

float
MEMCalculatorsWrapper:: getSuperKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_qqZZ;
    mem_->computeKD(kSMHiggs, kJHUGen, kqqZZ, kMCFM, &MEMs::probRatio, KD, ME_ggHiggs, ME_qqZZ);
    return pm4l_sig_/(pm4l_sig_+pm4l_bkg_*(1./KD-1));
}

float
MEMCalculatorsWrapper:: getGG0KD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_gg0Minus;
    mem_->computeKD(kSMHiggs, kJHUGen, k0minus, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_gg0Minus);
    return KD;
}

float
MEMCalculatorsWrapper:: getGG0HKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_gg0hPlus;
    mem_->computeKD(kSMHiggs, kJHUGen, k0hplus, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_gg0hPlus);
    return KD;
}

float
MEMCalculatorsWrapper:: getQQ1MinusKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_qq1Minus;
    mem_->computeKD(kSMHiggs, kJHUGen, k1minus, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_qq1Minus);
    return KD;
}

float
MEMCalculatorsWrapper:: getQQ1PlusKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_qq2Plus;
    mem_->computeKD(kSMHiggs, kJHUGen, k1plus, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_qq2Plus);
    return KD;
}

float
MEMCalculatorsWrapper:: getGG2PlusKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_gg2Plus;
    mem_->computeKD(kSMHiggs, kJHUGen, k2mplus_gg, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_gg2Plus);
    return KD;
}

float
MEMCalculatorsWrapper:: getQQ2PlusKD() {
    using namespace MEMNames;
    double KD,ME_ggHiggs,ME_qq2Plus;
    mem_->computeKD(kSMHiggs, kJHUGen, k2mplus_qqbar, kJHUGen, &MEMs::probRatio, KD, ME_ggHiggs, ME_qq2Plus);
    return KD;
}



float
MEMCalculatorsWrapper:: getInterferenceWeight() {
    return mem_->getMELAWeight();
}




