#ifndef STRINGPROCESS_H
#define STRINGPROCESS_H

#include <vector>
#include "Pythia8/Pythia.h"
#include "particledata.h"

namespace Smash {

class StringProcess {
 private:
  //int PDGidA, PDGidB;
  int baryonA, baryonB;
  //int chargeA, chargeB;
  //std::array<int,4> idqsetA;
  //std::array<int,4> idqsetB;
  double PPosA, PNegA, PPosB, PNegB;
  double massA, massB;
  double sqrtsAB, pabscomAB;
  PdgCode PDGcodeA, PDGcodeB;
  FourVector plabA, plabB;
  FourVector pcomA, pcomB;
  FourVector ucomAB;
  ThreeVector vcomAB;
  std::array<ThreeVector,4> evecBasisAB;

  int NpartFinal;
  int NpartString1;
  int NpartString2;

  double pLightConeMin;
  double xfracMin;

  double betapowS;

  double alphapowV;
  double betapowV;

  double sigmaQperp;

  double kappaString;

  Pythia8::Pythia *pythia;

 public:
  StringProcess();
  ~StringProcess();

  // final state array
  std::array<std::vector<int>,2> final_PDGid;
  // final_PDGid[0] : PDGid
  // final_PDGid[1] : if it is leading hadron (1 = Yes, 0 = No)
  std::array<std::vector<double>,5> final_pvec;
  // final_pvec[0] : energy
  // final_pvec[1] : px
  // final_pvec[2] : py
  // final_pvec[3] : pz
  // final_pvec[4] : mass
  std::array<std::vector<double>,2> final_tform;
  // final_tform[0] : formation time (fm)
  // final_tform[1] : suppression factor for the cross section (xtotfac)

  void set_pythia(Pythia8::Pythia *pythiaIn);

  void set_pLightConeMin(double pLightConeMinIn) {
    pLightConeMin = pLightConeMinIn;
  }

  void set_betapowS(double betapowSIn) { betapowS = betapowSIn; }

  void set_alphapowV(double alphapowVIn) { alphapowV = alphapowVIn; }
  void set_betapowV(double betapowVIn) { betapowV = betapowVIn; }

  void set_sigmaQperp(double sigmaQperpIn) { sigmaQperp = sigmaQperpIn; }

  void set_kappaString(double kappaStringIn) { kappaString = kappaStringIn; }

  bool init(const ParticleList &incomingList);
  bool init_lab(PdgCode &idAIn, PdgCode &idBIn, double massAIn, double massBIn,
                Pythia8::Vec4 pvecAIn, Pythia8::Vec4 pvecBIn);
  bool init_com(PdgCode &idAIn, PdgCode &idBIn, double massAIn, double massBIn,
                double sqrtsABIn);

  void make_incoming_com_momenta();
  void make_orthonormal_basis();
  void compute_incoming_lightcone_momenta();

  bool next_SDiff(int channel);  // single-diffractive : A + B -> A + X or X + B
  //bool next_SDiff_AX();  // single-diffractive : A + B -> A + X
  //bool next_SDiff_XB();  // single-diffractive : A + B -> X + B
  bool next_DDiff();             // double-diffractive : A + B -> X + X
  bool next_NDiff();             // non-diffractive
  bool next_BBbarAnn();          // baryon-antibaryon annihilation

  void reset_finalArray();
  int append_finalArray(FourVector &uString, ThreeVector &evecLong);

  void makeStringEnds(PdgCode &pdgcodeIn, int &idq1, int &idq2);
  //void makeStringEnds(std::array<int,4> &idqset, int &idq1, int &idq2);
  int fragmentString(int idq1, int idq2, double mString, ThreeVector &evecLong,
                     bool random_rotation);
};

}  // namespace Smash

#endif
