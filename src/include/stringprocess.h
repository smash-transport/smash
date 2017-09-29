#ifndef STRINGPROCESS_H
#define STRINGPROCESS_H

#include <vector>
#include "Pythia8/Pythia.h"
#include "particledata.h"

namespace Smash {

/**
 * \brief String excitation processes used in SMASH
 *
 * Only one instance of this class should be created.
 *
 * This class implements string excitation processes based on the UrQMD model
 * \iref{Bass:1998ca,Bleicher:1999xi} and subsequent fragmentation
 * according to the LUND/PYTHIA fragmentation scheme
 * \iref{Andersson:1983ia,Sjostrand:2014zea}.
 *
 * The class implemets the following functionality:
 * - given two colliding initial state particles it provides hadronic final
 *   state after single diffractive, double diffractive and non-diffractive
 *   string excitation
 * - owns a Pythia8::SigmaTotal object to compute corresponding cross-sections
 * - owns a Pythia object, that allows to fragment strings
 */
class StringProcess {
 private:
  double PPosA; //!< forward lightcone momentum p^{+} of incoming particle A
                //!< in the center of mass frame
  double PPosB; //!< forward lightcone momentum p^{+} of incoming particle B
                //!< in the center of mass frame
  double PNegA; //!< backward lightcone momentum p^{-} of incoming particle A
                //!< in the center of mass frame
  double PNegB; //!< backward lightcone momentum p^{-} of incoming particle B
                //!< in the center of mass frame
  double massA; //!< mass of incoming particle A
  double massB; //!< mass of incoming particle A
  double sqrtsAB; //!< sqrt of Mandelstam variable s
  double pabscomAB; //!< magnitude of three momentum of incoming particles
                    //!< in the center of mass frame
  std::array<PdgCode, 2> PDGcodes_; //!< PdgCodes of incoming particles
  std::array<FourVector, 2> plab_; //!< momenta of incoming particles in the lab frame
  std::array<FourVector, 2> pcom_; //!< momentum of incoming particles in the center of mass frame
  FourVector ucomAB; //!< velocity four vector of the center of mass in the lab frame
  ThreeVector vcomAB; //!< velocity three vector of the center of mass in the lab frame
  std::array<ThreeVector,4> evecBasisAB; //!< orthonormal basis vectors in the center of mass frame
                                         //!< where the 3rd one is parallel to momtentum
                                         //!< of incoming particle A
  int NpartFinal; //!< total number of final state particles
  std::array<int, 2> NpartString; //!< number of particles fragmented from strings
  double pmin_gluon_lightcone_; //!< the minimum lightcone momentum scale carried by gluon
  double pow_fgluon_beta_; //!< parameter for the gluon distribution function
                          //!< P(x) = 1/x * (1 - x)^{1 + pow_fgluon_beta_}
  double pow_fquark_alpha_; //!< parameter for the quark distribution function
                          //!< P(x) = x^{pow_fquark_alpha_ - 1} * (1 - x)^{pow_fquark_beta_ - 1}
  double pow_fquark_beta_; //!< parameter for the quark distribution function
                          //!< P(x) = x^{pow_fquark_alpha_ - 1} * (1 - x)^{pow_fquark_beta_ - 1}
  double sigma_qperp_; //!< transverse momentum spread of the excited strings
                      //!< transverse momenta of strings are sampled
                      //!< according to gaussian distribution with width sigma_qperp_
  double kappa_tension_string_; //!< string tension
  double time_collision; //!< time of collision in the computational frame
  double gamma_factor_com; //!< Lorentz gamma factor of center of mass
                           //!< in the computational frame

  /// PYTHIA object used in fragmentation
  std::unique_ptr<Pythia8::Pythia> pythia_;

  /// An object to compute cross-sections
  Pythia8::SigmaTotal pythia_sigmatot_;

 public:
  /** Constructor, initializes pythia. Should only be called once. */
  StringProcess();

  /**
   * Interface to pythia_sigmatot_ to compute cross-sections of A+B->
   * different final states.
   * \param pdg_a pdg code of incoming particle A
   * \param pdg_b pdg code of incoming particle B
   * \param sqrt_s collision energy in the center of mass frame [GeV]
   * \return array with single diffractive cross-sections AB->AX, AB->XB and
   * double diffractive AB->XX.
   */
  std::array<double, 3> cross_sections(int pdg_a, int pdg_b, double sqrt_s) {
    // This threshold magic is following Pythia. Todo(ryu): take care of this.
    double sqrts_threshold = 2. * (1. + 1.0e-6);
    pdg_a = std::abs(pdg_a);
    pdg_b = std::abs(pdg_b);
    pdg_a = (pdg_a > 1000) ? pdg_a : 10 * (pdg_a / 10) + 3;
    pdg_b = (pdg_b > 1000) ? pdg_b : 10 * (pdg_b / 10) + 3;
    sqrts_threshold += pythia_->particleData.m0(pdg_a) +
                       pythia_->particleData.m0(pdg_b);
    // Constant cross-section for sub-processes below threshold equal to
    // cross-section at the threshold.
    if (sqrt_s < sqrts_threshold) {
      sqrt_s = sqrts_threshold;
    }
    pythia_sigmatot_.calc(pdg_a, pdg_b, sqrt_s);
    return {pythia_sigmatot_.sigmaAX(),
            pythia_sigmatot_.sigmaXB(),
            pythia_sigmatot_.sigmaXX()};
  }
  /**
   * final state array
   * which must be accessed after the collision
   */
  ParticleList final_state;
  /**
   * set the minimum lightcone momentum scale carried by gluon.
   * This is relevant for the double-diffractive process.
   * The minimum lightcone momentum fraction is set to be pmin_gluon_lightcone_/sqrtsAB.
   * \param pLightConeMinIn is a value that we want to use for pmin_gluon_lightcone_.
   */
  void set_pmin_gluon_lightcone(double pLightConeMinIn) {
    pmin_gluon_lightcone_ = pLightConeMinIn;
  }
  /**
   * lightcone momentum fraction of gluon is sampled
   * according to probability distribution P(x) = 1/x * (1 - x)^{1 + pow_fgluon_beta_}
   * in double-diffractive processes.
   * \param betapowSIn is a value that we want to use for pow_fgluon_beta_.
   */
  void set_pow_fgluon(double betapowSIn) { pow_fgluon_beta_ = betapowSIn; }
  /**
   * lightcone momentum fraction of quark is sampled
   * according to probability distribution P(x) = x^{pow_fquark_alpha_ - 1} * (1 - x)^{pow_fquark_beta_ - 1}
   * in non-diffractive processes.
   * \param alphapowVIn is a value that we want to use for pow_fquark_alpha_.
   * \param betapowVIn is a value that we want to use for pow_fquark_beta_.
   */
  void set_pow_fquark(double alphapowVIn, double betapowVIn){
    pow_fquark_alpha_ = alphapowVIn;
    pow_fquark_beta_ = betapowVIn;
  }
  /**
   * set the average amount of transverse momentum transfer sigma_qperp_.
   * \param sigmaQperpIn is a value that we want to use for sigma_qperp_.
   */
  void set_sigma_qperp_(double sigmaQperpIn) { sigma_qperp_ = sigmaQperpIn; }
  /**
   * set the string tension kappaString which is used in append_final_state.
   * \param kappaStringIn is a value that we want to use for kappaString.
   */
  void set_tension_string(double kappaStringIn) { kappa_tension_string_ = kappaStringIn; }
  /**
   * initialization
   * feed intial particles, time of collision and gamma factor of the center of mass.
   * \param incomingList is the list of initial state particles.
   * \param tcollIn is time of collision.
   * \param gammaFacIn gamma factor of the center of mass.
   */
  void init(const ParticleList &incomingList, double tcollIn, double gammaFacIn);
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
  int append_final_state(const FourVector &uString,
                         const ThreeVector &evecLong);
  /**
   * Construct diquark from two quarks. Order does not matter.
   * \param q1 PDG code of quark 1
   * \param q2 PDG code of quark 2
   * \return PDG code of diquark composed of q1 and q2
   */
  int diquark_from_quarks(int q1, int q2);

  /**
   * make a random selection to determine partonic contents at the string ends.
   * \param pdgcode_in is PdgCode of hadron which transforms into a string.
   * \param idq1 is PDG id of quark or anti-diquark.
   * \param idq2 is PDG id of anti-quark or diquark.
   */
  void make_string_ends(const PdgCode &pdgcode_in, int &idq1, int &idq2);

  /// Easy setter of Pythia Vec4 from SMASH
  Pythia8::Vec4 set_Vec4(double energy, const ThreeVector &mom) {
    return Pythia8::Vec4(energy, mom.x1(), mom.x2(), mom.x3());
  }

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
  int fragment_string(int idq1, int idq2, double mString, ThreeVector &evecLong,
                     bool random_rotation);
};

}  // namespace Smash

#endif
