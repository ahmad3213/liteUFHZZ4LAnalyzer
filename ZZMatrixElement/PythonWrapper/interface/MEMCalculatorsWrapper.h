#include <TLorentzVector.h>
#include <vector>
#include <string>
#include <DataFormats/Math/interface/LorentzVector.h>
struct MEMs;

class MEMCalculatorsWrapper 
{
 public:
  MEMCalculatorsWrapper(double collisionEnergy = 8, double sKD_mass = 125.6) ;

  ~MEMCalculatorsWrapper() ;

  struct Angles {
    float costhetastar,costheta1,costheta2,phi,phistar1;
  };
  Angles computeAngles(TLorentzVector Z1_lept1, int Z1_lept1Id,
		   TLorentzVector Z1_lept2, int Z1_lept2Id,
		   TLorentzVector Z2_lept1, int Z2_lept1Id,
		   TLorentzVector Z2_lept2, int Z2_lept2Id) ;
  Angles computeAngles(const math::XYZTLorentzVector & Z1_lept1, int Z1_lept1Id,
		   const math::XYZTLorentzVector & Z1_lept2, int Z1_lept2Id,
		   const math::XYZTLorentzVector & Z2_lept1, int Z2_lept1Id,
		   const math::XYZTLorentzVector & Z2_lept2, int Z2_lept2Id) ;


  void  computeAll(TLorentzVector Z1_lept1, int Z1_lept1Id,
		   TLorentzVector Z1_lept2, int Z1_lept2Id,
		   TLorentzVector Z2_lept1, int Z2_lept1Id,
		   TLorentzVector Z2_lept2, int Z2_lept2Id) ;

  std::vector<std::pair<std::string,float> > computeNew(
              const math::XYZTLorentzVector & Z1_lept1, int Z1_lept1Id,
              const math::XYZTLorentzVector & Z1_lept2, int Z1_lept2Id,
              const math::XYZTLorentzVector & Z2_lept1, int Z2_lept1Id,
              const math::XYZTLorentzVector & Z2_lept2, int Z2_lept2Id,
              const std::vector<math::XYZTLorentzVector> & jetP4s) ;


  float getKD() ;
  float getSuperKD() ;
  float getGG0KD() ;
  float getGG0HKD() ;
  float getQQ1MinusKD() ;
  float getQQ1PlusKD() ;
  float getGG2PlusKD() ;
  float getQQ2PlusKD() ;
  float getInterferenceWeight() ;

 private:

  MEMs * mem_;
  double pm4l_sig_;
  double pm4l_bkg_;
};
