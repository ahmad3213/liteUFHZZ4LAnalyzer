#include "TMCFMUtils.hh"
#include "TMath.h"

using namespace std;
using namespace PDGHelpers;


namespace TMCFMUtils{
  const std::vector<intQuad_t> MCFMHash_QQVVQQAny = Hash_QQVVQQAny();
}

void TMCFMUtils::AssociatedParticleOrdering_QQVVQQAny(int iSel, int jSel, int rSel, int sSel, int order[2]){
  const std::vector<TMCFMUtils::intQuad_t>& hash = MCFMHash_QQVVQQAny;
  bool outFound=false;
  for (unsigned int ih=0; ih<hash.size(); ih++){
    if (
      !(
      (PDGHelpers::isAnUnknownJet(iSel) || iSel==hash.at(ih)[0])
      &&
      (PDGHelpers::isAnUnknownJet(jSel) || jSel==hash.at(ih)[1])
      )
      ) continue;
    // Final particles are q
    if (PDGHelpers::isAJet(rSel) && PDGHelpers::isAJet(sSel)){
      if (
        (PDGHelpers::isAnUnknownJet(rSel) || rSel==hash.at(ih)[2])
        &&
        (PDGHelpers::isAnUnknownJet(sSel) || sSel==hash.at(ih)[3])
        ){
        order[0]=0;
        order[1]=1;
        outFound=true;
        //cout << "Hash requested outgoing "<< hash.at(ih)[2] << " " << hash.at(ih)[3] << ", unswapped r, s = " << rSel << " " << sSel << endl;
      }
      else if (
        (PDGHelpers::isAnUnknownJet(rSel) || rSel==hash.at(ih)[3])
        &&
        (PDGHelpers::isAnUnknownJet(sSel) || sSel==hash.at(ih)[2])
        ){
        order[0]=1;
        order[1]=0;
        outFound=true;
        //cout << "Hash requested outgoing "<< hash.at(ih)[2] << " " << hash.at(ih)[3] << ", swapped r, s = " << rSel << " " << sSel << endl;
      }
    }
    // Final particles l/nu
    else if ((PDGHelpers::isALepton(rSel) || PDGHelpers::isANeutrino(rSel)) && (PDGHelpers::isALepton(sSel) || PDGHelpers::isANeutrino(sSel))){
      if (abs(hash.at(ih)[0])==abs(hash.at(ih)[1]) && abs(hash.at(ih)[0])==abs(hash.at(ih)[2]) && abs(hash.at(ih)[0])==abs(hash.at(ih)[3])) continue; // Do not consider the ordering in uquq_uquq or dqdq_dqdq

      if (
        (
        TMath::Sign(1, rSel)==TMath::Sign(1, hash.at(ih)[2]) &&
        ((PDGHelpers::isALepton(rSel) && PDGHelpers::isDownTypeQuark(hash.at(ih)[2])) || (PDGHelpers::isANeutrino(rSel) && PDGHelpers::isUpTypeQuark(hash.at(ih)[2])))
        )
        &&
        (
        TMath::Sign(1, sSel)==TMath::Sign(1, hash.at(ih)[3]) &&
        ((PDGHelpers::isALepton(sSel) && PDGHelpers::isDownTypeQuark(hash.at(ih)[3])) || (PDGHelpers::isANeutrino(sSel) && PDGHelpers::isUpTypeQuark(hash.at(ih)[3])))
        )
        ){
        order[0]=0;
        order[1]=1;
        outFound=true;
        //cout << "Hash requested outgoing "<< hash.at(ih)[2] << " " << hash.at(ih)[3] << ", unswapped r, s = " << rSel << " " << sSel << endl;
      }
      else if (
        (
        TMath::Sign(1, rSel)==TMath::Sign(1, hash.at(ih)[3]) &&
        ((PDGHelpers::isALepton(rSel) && PDGHelpers::isDownTypeQuark(hash.at(ih)[3])) || (PDGHelpers::isANeutrino(rSel) && PDGHelpers::isUpTypeQuark(hash.at(ih)[3])))
        )
        &&
        (
        TMath::Sign(1, sSel)==TMath::Sign(1, hash.at(ih)[2]) &&
        ((PDGHelpers::isALepton(sSel) && PDGHelpers::isDownTypeQuark(hash.at(ih)[2])) || (PDGHelpers::isANeutrino(sSel) && PDGHelpers::isUpTypeQuark(hash.at(ih)[2])))
        )
        ){
        order[0]=1;
        order[1]=0;
        outFound=true;
        //cout << "Hash requested outgoing "<< hash.at(ih)[2] << " " << hash.at(ih)[3] << ", swapped r, s = " << rSel << " " << sSel << endl;
      }
      if (PDGHelpers::getCoupledVertex(rSel, sSel)!=PDGHelpers::getCoupledVertex(hash.at(ih)[2], hash.at(ih)[3])) outFound=false;
    }
    if (outFound) break;
  }
  if (!outFound){ for (unsigned int ip=0; ip<2; ip++) order[ip]=-1; }
}
std::vector<TMCFMUtils::intQuad_t> TMCFMUtils::Hash_QQVVQQAny(){
  std::vector<TMCFMUtils::intQuad_t> pcfg;
  std::vector<TMCFMUtils::intQuad_t> hash_qqvvqq = TMCFMUtils::Hash_QQVVQQ();
  std::vector<TMCFMUtils::intQuad_t> hash_qqvvqqstrong = TMCFMUtils::Hash_QQVVQQStrong();
  for (unsigned int c=0; c<hash_qqvvqq.size(); c++) pcfg.push_back(hash_qqvvqq.at(c));
  for (unsigned int c=0; c<hash_qqvvqqstrong.size(); c++) pcfg.push_back(hash_qqvvqqstrong.at(c));
  return pcfg;
}
std::vector<TMCFMUtils::intQuad_t> TMCFMUtils::Hash_QQVVQQ(){
  /*
  Based on the following cases in MCFM:
  parameter(
  & uqcq_uqcq=1,uquq_uquq=2,dqsq_dqsq=3,
  & dqdq_dqdq=4,uqbq_uqbq=5,dqcq_dqcq=6,
  & dquq_dquq=7,dqcq_uqsq=8,uqsq_dqcq=9
  & jmax=12)
  integer,parameter:: j1(jmax)=(/1,2,8,8, 7,2,7,1, 1,7,2,7/)
  integer,parameter:: j2(jmax)=(/2,1,7,7, 2,7,1,7, 7,1,7,2/)
  integer,parameter:: j7(jmax)=(/7,7,2,1, 1,8,2,8, 2,8,1,8/)
  integer,parameter:: j8(jmax)=(/8,8,1,2, 8,1,8,2, 8,2,8,1/)
  */
  std::vector<TMCFMUtils::intQuad_t> base_cfg;
  // uc_uc
  base_cfg.push_back(TMCFMUtils::intQuad_t(2, 4, 2, 4));
  // ds_ds
  base_cfg.push_back(TMCFMUtils::intQuad_t(1, 3, 1, 3));
  base_cfg.push_back(TMCFMUtils::intQuad_t(1, 5, 1, 5));
  base_cfg.push_back(TMCFMUtils::intQuad_t(3, 5, 3, 5));
  // ub_ub
  base_cfg.push_back(TMCFMUtils::intQuad_t(2, 3, 2, 3));
  base_cfg.push_back(TMCFMUtils::intQuad_t(2, 5, 2, 5));
  base_cfg.push_back(TMCFMUtils::intQuad_t(4, 5, 4, 5));
  // dc_dc
  base_cfg.push_back(TMCFMUtils::intQuad_t(1, 4, 1, 4));
  // du_du
  base_cfg.push_back(TMCFMUtils::intQuad_t(1, 2, 1, 2));
  base_cfg.push_back(TMCFMUtils::intQuad_t(3, 4, 3, 4));
  // dc_us
  base_cfg.push_back(TMCFMUtils::intQuad_t(1, 4, 2, 3));
  // us_dc
  base_cfg.push_back(TMCFMUtils::intQuad_t(2, 3, 1, 4));
  // uu_uu
  base_cfg.push_back(TMCFMUtils::intQuad_t(2));
  base_cfg.push_back(TMCFMUtils::intQuad_t(4));
  // dd_dd
  base_cfg.push_back(TMCFMUtils::intQuad_t(1));
  base_cfg.push_back(TMCFMUtils::intQuad_t(3));
  base_cfg.push_back(TMCFMUtils::intQuad_t(5));

  std::vector<TMCFMUtils::intQuad_t> jcfg;
  jcfg.push_back(TMCFMUtils::intQuad_t(0, 1, 2, 3));
  jcfg.push_back(TMCFMUtils::intQuad_t(1, 0, 2, 3));
  jcfg.push_back(TMCFMUtils::intQuad_t(3, 2, 1, 0));
  jcfg.push_back(TMCFMUtils::intQuad_t(3, 2, 0, 1));
  jcfg.push_back(TMCFMUtils::intQuad_t(2, 1, 0, 3));
  jcfg.push_back(TMCFMUtils::intQuad_t(1, 2, 3, 0));
  jcfg.push_back(TMCFMUtils::intQuad_t(2, 0, 1, 3));
  jcfg.push_back(TMCFMUtils::intQuad_t(0, 2, 3, 1));
  jcfg.push_back(TMCFMUtils::intQuad_t(0, 2, 1, 3));
  jcfg.push_back(TMCFMUtils::intQuad_t(2, 0, 3, 1));
  jcfg.push_back(TMCFMUtils::intQuad_t(1, 2, 0, 3));
  jcfg.push_back(TMCFMUtils::intQuad_t(2, 1, 3, 0));

  std::vector<TMCFMUtils::intQuad_t> pcfg;
  for (unsigned int j=0; j<jcfg.size(); j++){
    for (unsigned int p=0; p<base_cfg.size(); p++){
      TMCFMUtils::intQuad_t cfg;
      for (unsigned int ipos=0; ipos<4; ipos++){
        int idpos = jcfg.at(j)[ipos];
        int idAssigned = base_cfg.at(p)[ipos];
        if ((idpos<2 && ipos>=2) || (idpos>=2 && ipos<2)) idAssigned = -idAssigned;
        cfg[jcfg.at(j)[ipos]] = idAssigned;
      }
      pcfg.push_back(cfg);
    }
  }
  // Uncommenting the lines below prints out the hash when the library is loaded.
  /*
  for (unsigned int ic=0; ic<pcfg.size(); ic++) std::cout
  << "TMCFMUtils::Hash_QQVVQQAny: Hash configuration " << ic << " requests ids=( "
  << pcfg.at(ic)[0] << ", " << pcfg.at(ic)[1] << ", " << pcfg.at(ic)[2] << ", " << pcfg.at(ic)[3]
  << ")" << std::endl;
  */
  return pcfg;
}
std::vector<TMCFMUtils::intQuad_t> TMCFMUtils::Hash_QQVVQQStrong(){
  /*
  Based on the following cases in MCFM:
  call qq4lggampf(1,2,3,4,5,6,7,8,3,4,za,zb,msqgg)
  msq(1,-1)=msq(1,-1)+stat*aveqq*msqgg(1)
  call qq4lggampf(2,1,3,4,5,6,7,8,3,4,za,zb,msqgg)
  msq(-1,1)=msq(-1,1)+stat*aveqq*msqgg(1)
  call qq4lggampf(7,8,3,4,5,6,1,2,3,4,za,zb,msqgg)
  msq(0,0)=msq(0,0)+avegg*(3d0*msqgg(1)+2d0*msqgg(2))
  call qq4lggampf(7,2,3,4,5,6,1,8,3,4,za,zb,msqgg)
  msq(0,-1)=msq(0,-1)+aveqg*msqgg(1)
  call qq4lggampf(7,1,3,4,5,6,2,8,3,4,za,zb,msqgg)
  msq(-1,0)=msq(-1,0)+aveqg*msqgg(1)
  call qq4lggampf(2,8,3,4,5,6,1,7,3,4,za,zb,msqgg)
  msq(0,1)=msq(0,1)+aveqg*msqgg(1)
  call qq4lggampf(1,8,3,4,5,6,2,7,3,4,za,zb,msqgg)
  msq(1,0)=msq(1,0)+aveqg*msqgg(1)
  */
  std::vector<TMCFMUtils::intQuad_t> base_cfg;
  // Start with qqb_gg
  for (int iq=1; iq<=5; iq++) base_cfg.push_back(TMCFMUtils::intQuad_t(iq, -iq, 21, 21));

  std::vector<TMCFMUtils::intQuad_t> jcfg;
  jcfg.push_back(TMCFMUtils::intQuad_t(0, 1, 2, 3)); // qqb->gg
  jcfg.push_back(TMCFMUtils::intQuad_t(1, 0, 2, 3)); // qbq->gg
  jcfg.push_back(TMCFMUtils::intQuad_t(2, 3, 0, 1)); // gg->qbq
  jcfg.push_back(TMCFMUtils::intQuad_t(2, 1, 0, 3)); // gqb->qbg
  jcfg.push_back(TMCFMUtils::intQuad_t(2, 0, 1, 3)); // qbg->qbg
  jcfg.push_back(TMCFMUtils::intQuad_t(1, 3, 2, 0)); // gq->gq
  jcfg.push_back(TMCFMUtils::intQuad_t(0, 3, 2, 1)); // qg->gq

  std::vector<TMCFMUtils::intQuad_t> pcfg;
  for (unsigned int j=0; j<jcfg.size(); j++){
    for (unsigned int p=0; p<base_cfg.size(); p++){
      TMCFMUtils::intQuad_t cfg;
      for (unsigned int ipos=0; ipos<4; ipos++){
        int idpos = jcfg.at(j)[ipos];
        int idAssigned = base_cfg.at(p)[ipos];
        if (((idpos<2 && ipos>=2) || (idpos>=2 && ipos<2)) && !PDGHelpers::isAGluon(idAssigned)) idAssigned = -idAssigned;
        cfg[jcfg.at(j)[ipos]] = idAssigned;
      }
      pcfg.push_back(cfg);
    }
  }
  // Uncommenting the lines below prints out the hash when the library is loaded.
  /*
  for (unsigned int ic=0; ic<pcfg.size(); ic++) std::cout
  << "TMCFMUtils::Hash_QQVVQQAny: Hash configuration " << ic << " requests ids=( "
  << pcfg.at(ic)[0] << ", " << pcfg.at(ic)[1] << ", " << pcfg.at(ic)[2] << ", " << pcfg.at(ic)[3]
  << ")" << std::endl;
  */
  return pcfg;
}
