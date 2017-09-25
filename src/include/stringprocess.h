#ifndef STRINGPROCESS_H
#define STRINGPROCESS_H

#include <vector>
#include "Pythia8/Pythia.h"
#include "particledata.h"

namespace Smash {

/**
 * \brief String excitation processes used in SMASH
 * This class implements string excitation processes based on the UrQMD model
 * \iref{Bass:1998ca,Bleicher:1999xi} and subsequent fragmentation
 * according to the LUND/PYTHIA fragmentation scheme \iref{Andersson:1983ia,Sjostrand:2014zea}.
 * There are single and double diffractive and non-diffractive excitation subprocesses
 * as described in SMASH wiki.
 */
class StringProcess {
 private:
  int baryonA, baryonB;
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

  double pmin_gluon_lightcone;
  double xmin_gluon_fraction;

  double pow_fgluon_beta;

  double pow_fquark_alpha;
  double pow_fquark_beta;

  double sigmaQperp;

  double kappa_tension_string;

  double time_collision;

  double gamma_factor_com;

  Pythia8::Pythia *pythia;

 public:
  /** constructor */
  StringProcess();
  /** destructor */
  ~StringProcess();
  /**
   * final state array
   * which must be accessed after the collision
   */
  ParticleList final_state;
  /**
   * set the pointer of PYTHIA object
   * \param pythiaIn is a pointer of object created and initialized in somewhere else.
   */
  void set_pythia(Pythia8::Pythia *pythiaIn);
  /**
   * set the minimum lightcone momentum scale carried by gluon.
   * This is relevant for the double-diffractive process.
   * The minimum lightcone momentum fraction is set to be pmin_gluon_lightcone/sqrtsAB.
   * \param pLightConeMinIn is a value that we want to use for pmin_gluon_lightcone.
   */
  void set_pmin_gluon_lightcone(double pLightConeMinIn) {
    pmin_gluon_lightcone = pLightConeMinIn;
  }
  /**
   * lightcone momentum fraction of gluon is sampled
   * according to probability distribution P(x) = 1/x * (1 - x)^{1 + betapowS}
   * in double-diffractive processes.
   * \param betapowSIn is a value that we want to use for pow_fgluon_beta.
   */
  void set_pow_fgluon(double betapowSIn) { pow_fgluon_beta = betapowSIn; }
  /**
   * lightcone momentum fraction of quark is sampled
   * according to probability distribution P(x) = x^{alphapowV - 1} * (1 - x)^{betapowV - 1}
   * in non-diffractive processes.
   * \param alphapowVIn is a value that we want to use for pow_fquark_alpha.
   * \param betapowVIn is a value that we want to use for pow_fquark_beta.
   */
  void set_pow_fquark(double alphapowVIn, double betapowVIn){
    pow_fquark_alpha = alphapowVIn;
    pow_fquark_beta = betapowVIn;
  }
  /**
   * set the average amount of transverse momentum transfer sigmaQperp.
   * \param sigmaQperpIn is a value that we want to use for sigmaQperp.
   */
  void set_sigma_Qperp(double sigmaQperpIn) { sigmaQperp = sigmaQperpIn; }
  /**
   * set the string tension kappaString which is used in append_final_state.
   * \param kappaStringIn is a value that we want to use for kappaString.
   */
  void set_tension_string(double kappaStringIn) { kappa_tension_string = kappaStringIn; }
  /**
   * initialization
   * feed intial particles, time of collision and gamma factor of the center of mass.
   * \param incomingList is the list of initial state particles.
   * \param tcollIn is time of collision.
   * \param gammaFacIn gamma factor of the center of mass.
   * \return if initialization is successful.
   */
  bool init(const ParticleList &incomingList, double tcollIn, double gammaFacIn);
  /** boost the momenta of incoming particles into the center of mass frame. */
  void make_incoming_com_momenta();
  /**
   * compute three orthonormal basis vectors in the center of mass frame
   * such that one vector is along with the three-momentum of particle A.
   */
  void make_orthonormal_basis();
  /**
   * compute the lightcone momenta of incoming particles
   * where the longitudinal direction is set to be same
   * as that of the three-momentum of particle A.
   */
  void compute_incoming_lightcone_momenta();
  /**
   * Single-diffractive process
   * is based on single pomeron exchange described in \iref{Ingelman:1984ns}.
   * \param channel specifies which hadron to excite into a string.
   * channel = 1 : A + B -> A + X
   * channel = 2 : A + B -> X + B
   * \return whether the process is successfully implemented.
   */
  bool next_SDiff(int channel);
  /**
   * Double-diffractive process ( A + B -> X + X )
   * is similar to the single-diffractive process,
   * but lightcone momenta of gluons are sampled
   * in the same was as the UrQMD model \iref{Bass:1998ca,Bleicher:1999xi}.
   * String masses are computed after pomeron exchange
   * aquiring transverse momentum transfer.
   * \return whether the process is successfully implemented.
   */
  bool next_DDiff();
  /**
   * Non-diffractive process
   * is modelled in accordance with dual-topological approach \iref{Capella:1978ig}.
   * This involves a parton exchange in conjunction with momentum transfer.
   * Probability distribution function of the lightcone momentum fraction
   * carried by quark is based on the UrQMD model \iref{Bass:1998ca,Bleicher:1999xi}.
   * \return whether the process is successfully implemented.
   */
  bool next_NDiff();
  /**
   * Baryon-antibaryon annihilation process
   * Based on what UrQMD \iref{Bass:1998ca,Bleicher:1999xi} does,
   * it create two mesonic strings after annihilating one quark-antiquark pair.
   * Each string has mass equal to half of sqrts.
   * \return whether the process is successfully implemented.
   */
  bool next_BBbarAnn();
  /**
   * compute the formation time and fill the arrays with final-state particles
   * as described in \iref{Andersson:1983ia}.
   * \param uString is velocity four vector of the string.
   * \param evecLong is unit 3-vector in which string is stretched.
   * \return number of hadrons fragmented out of string.
   */
  int append_final_state(FourVector &uString, ThreeVector &evecLong);
  /**
   * make a random selection to determine partonic contents at the string ends.
   * \param pdgcodeIn is PdgCode of hadron which transforms into a string.
   * \param idq1 is PDG id of quark or anti-diquark.
   * \param idq2 is PDG id of anti-quark or diquark.
   */
  void makeStringEnds(PdgCode &pdgcodeIn, int &idq1, int &idq2);
  //void makeStringEnds(std::array<int,4> &idqset, int &idq1, int &idq2);
  /**
   * perform string fragmentation to determine species and momenta of hadrons
   * by implementing PYTHIA 8.2 \iref{Andersson:1983ia,Sjostrand:2014zea}.
   * \param idq1 is PDG id of quark or anti-diquark (carrying color index).
   * \param idq2 is PDG id of diquark or anti-quark (carrying anti-color index).
   * \param mString is the string mass.
   * \param evecLong is unit 3-vector specifying the direction of diquark or anti-diquark.
   * \param random_rotation is whether or not we randomly rotate the orientation.
   * \return number of hadrons fragmented out of string.
   */
  int fragmentString(int idq1, int idq2, double mString, ThreeVector &evecLong,
                     bool random_rotation);
};

}  // namespace Smash

#endif
