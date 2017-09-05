#ifndef STRINGPROCESS_H
#define STRINGPROCESS_H

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <vector>
#include "Pythia8/Pythia.h"
#include "particledata.h"

//using namespace Pythia8;
using namespace std;

namespace Smash {

class StringProcess {
 private:
  int PDGidA, PDGidB;
  int baryonA, baryonB;
  int chargeA, chargeB;
  int *idqsetA;
  int *idqsetB;
  double massA, massB;
  double sqrtsAB, pabscomAB;
  FourVector plabA;
  FourVector plabB;
  FourVector pcomA;
  FourVector pcomB;
  FourVector ucomAB;
  ThreeVector vcomAB;
  ThreeVector *evecBasisAB;
  //double **evecBasisAB;

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

  int CINI, CFIN;
  int BINI, BFIN;
  double EINI, EFIN;
  double pxINI, pxFIN;
  double pyINI, pyFIN;
  double pzINI, pzFIN;

  Pythia8::Pythia *pythia;

  //double XSecTot, XSecEl, XSecInel, XSecAX, XSecXB, XSecXX, XSecAXB, XSecND;
  //double *XSecSummed;
  //SigmaTotal sigmaTot;

 public:
  StringProcess();
  ~StringProcess();

  // final state array
  vector<int> *final_PDGid;
  // final_PDGid[0] : PDGid
  // final_PDGid[1] : if it is leading hadron (1 = Yes, 0 = No)
  vector<double> *final_pvec;
  // final_pvec[0] : energy
  // final_pvec[1] : px
  // final_pvec[2] : py
  // final_pvec[3] : pz
  // final_pvec[4] : mass
  vector<double> *final_tform;
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
  bool init_lab(int idAIn, int idBIn, double massAIn, double massBIn,
                Pythia8::Vec4 pvecAIn, Pythia8::Vec4 pvecBIn);
  bool init_com(int idAIn, int idBIn, double massAIn, double massBIn,
                double sqrtsABIn);

  //bool next_Inel();
  bool next_SDiff_AX();  // single-diffractive : A + B -> A + X
  bool next_SDiff_XB();  // single-diffractive : A + B -> X + B
  bool next_DDiff_XX();  // double-diffractive : A + B -> X + X
  bool next_NDiff();     // non-diffractive
  bool next_BBbarAnn(); // baryon-antibaryon annihilation

  //double get_XSecTot() { return XSecTot; }
  //double get_XSecEl() { return XSecEl; }
  //double get_XSecInel() { return XSecInel; }
  //double get_XSecAX() { return XSecAX; }
  //double get_XSecXB() { return XSecXB; }
  //double get_XSecXX() { return XSecXX; }
  //double get_XSecAXB() { return XSecAXB; }
  //double get_XSecND() { return XSecND; }

  void reset_finalArray();
  int append_finalArray(FourVector &uString, ThreeVector &evecLong);
  bool check_conservation();

  void PDGid2idqset(int pdgid, int *idqset);
  void makeStringEnds(int *idqset, int *idq1, int *idq2);
  int fragmentString(int idq1, int idq2, double mString, ThreeVector &evecLong,
                     bool ranrot);

  double sample_XSDIS(double xmin, double b);
  double sample_XVDIS(double xmin, double a, double b);
  double sample_Qperp(double sigQ);
};

}  // namespace Smash

#endif
