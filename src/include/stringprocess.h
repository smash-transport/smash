/*
 *
 *    Copyright (c) 2017-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_STRINGPROCESS_H_
#define SRC_INCLUDE_STRINGPROCESS_H_

#include <memory>
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
  // The following 4 variables are in the center of mass frame
  /// forward lightcone momentum p^{+} of incoming particle A
  double PPosA_;
  /// forward lightcone momentum p^{+} of incoming particle B
  double PPosB_;
  /// backward lightcone momentum p^{-} of incoming particle A
  double PNegA_;
  /// backward lightcone momentum p^{-} of incoming particle B
  double PNegB_;
  /// masses of incoming particles
  double massA_, massB_;
  /// sqrt of Mandelstam variable s of collision
  double sqrtsAB_;
  /// PdgCodes of incoming particles
  std::array<PdgCode, 2> PDGcodes_;
  /// momenta of incoming particles in the lab frame
  std::array<FourVector, 2> plab_;
  /// momenta of incoming particles in the center of mass frame
  std::array<FourVector, 2> pcom_;
  /// velocity four vector of the center of mass in the lab frame
  FourVector ucomAB_;
  /// velocity three vector of the center of mass in the lab frame
  ThreeVector vcomAB_;
  /**
   * Orthonormal basis vectors in the center of mass frame,
   * where the 3rd one is parallel to momentum of incoming particle A
   */
  std::array<ThreeVector, 4> evecBasisAB_;
  /// total number of final state particles
  int NpartFinal_;
  /// number of particles fragmented from strings
  std::array<int, 2> NpartString_;
  /// the minimum lightcone momentum scale carried by gluon
  double pmin_gluon_lightcone_;
  /** parameter for the gluon distribution function
   *  P(x) = 1/x * (1 - x)^{1 + pow_fgluon_beta_}
   */
  double pow_fgluon_beta_;
  /** parameter for the quark distribution function
   * P(x) = x^{pow_fquark_alpha_ - 1} * (1 - x)^{pow_fquark_beta_ - 1}
   */
  double pow_fquark_alpha_;
  /** parameter for the quark distribution function
   * P(x) = x^{pow_fquark_alpha_ - 1} * (1 - x)^{pow_fquark_beta_ - 1}
   */
  double pow_fquark_beta_;
  /** Transverse momentum spread of the excited strings.
   * Transverse momenta of strings are sampled according to gaussian
   *  distribution with width sigma_qperp_
   */
  double sigma_qperp_;
  /// string tension
  double kappa_tension_string_;
  /// time of collision in the computational frame
  double time_collision_;
  /// Lorentz gamma factor of center of mass in the computational frame
  double gamma_factor_com_;

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
  std::array<double, 3> cross_sections_diffractive(int pdg_a, int pdg_b,
                                                   double sqrt_s) {
    // This threshold magic is following Pythia. Todo(ryu): take care of this.
    double sqrts_threshold = 2. * (1. + 1.0e-6);
    pdg_a = std::abs(pdg_a);
    pdg_b = std::abs(pdg_b);
    // In the case of mesons, the corresponding vector meson masses
    // are used to evaluate the energy threshold.
    const int pdg_a_mod = (pdg_a > 1000) ? pdg_a : 10 * (pdg_a / 10) + 3;
    const int pdg_b_mod = (pdg_b > 1000) ? pdg_b : 10 * (pdg_b / 10) + 3;
    sqrts_threshold +=
        pythia_->particleData.m0(pdg_a_mod) +
        pythia_->particleData.m0(pdg_b_mod);
    // Constant cross-section for sub-processes below threshold equal to
    // cross-section at the threshold.
    if (sqrt_s < sqrts_threshold) {
      sqrt_s = sqrts_threshold;
    }
    pythia_sigmatot_.calc(pdg_a, pdg_b, sqrt_s);
    return {pythia_sigmatot_.sigmaAX(), pythia_sigmatot_.sigmaXB(),
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
   * The minimum lightcone momentum fraction is set to be
   * pmin_gluon_lightcone_/sqrtsAB.
   * \param pLightConeMinIn is a value that we want to use for
   * pmin_gluon_lightcone_.
   */
  void set_pmin_gluon_lightcone(double p_light_cone_min) {
    pmin_gluon_lightcone_ = p_light_cone_min;
  }
  /**
   * lightcone momentum fraction of gluon is sampled
   * according to probability distribution P(x) = 1/x * (1 - x)^{1 +
   * pow_fgluon_beta_}
   * in double-diffractive processes.
   * \param betapow is a value that we want to use for pow_fgluon_beta_.
   */
  void set_pow_fgluon(double betapow) { pow_fgluon_beta_ = betapow; }
  /**
   * lightcone momentum fraction of quark is sampled
   * according to probability distribution
   * \f$ P(x) = x^{pow_fquark_alpha_ - 1} * (1 - x)^{pow_fquark_beta_ - 1} \f$
   * in non-diffractive processes.
   * \param alphapow is a value that we want to use for pow_fquark_alpha_.
   * \param betapow is a value that we want to use for pow_fquark_beta_.
   */
  void set_pow_fquark(double alphapow, double betapow) {
    pow_fquark_alpha_ = alphapow;
    pow_fquark_beta_ = betapow;
  }
  /**
   * set the average amount of transverse momentum transfer sigma_qperp_.
   * \param sigma_qperp is a value that we want to use for sigma_qperp_.
   */
  void set_sigma_qperp_(double sigma_qperp) { sigma_qperp_ = sigma_qperp; }
  /**
   * set the string tension, which is used in append_final_state.
   * \param kappa_string is a value that we want to use for string tension.
   */
  void set_tension_string(double kappa_string) {
    kappa_tension_string_ = kappa_string;
  }
  /**
   * initialization
   * feed intial particles, time of collision and gamma factor of the center of
   * mass.
   * \param incoming is the list of initial state particles.
   * \param tcoll is time of collision.
   * \param gamma gamma factor of the center of mass.
   */
  void init(const ParticleList &incoming, double tcoll, double gamma);
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
   * Prepare kinematics of two strings, fragment them and append to final_state
   * \param quarks pdg ids of string ends
   * \param pstr_com 4-momenta of strings in the C.o.m. frame
   * \param m_str masses of strings
   * \return whether fragmentations and final state creation was successful
   */
  bool make_final_state_2strings(
      const std::array<std::array<int, 2>, 2> &quarks,
      const std::array<FourVector, 2> &pstr_com,
      const std::array<double, 2> &m_str);

  /**
   * Single-diffractive process
   * is based on single pomeron exchange described in \iref{Ingelman:1984ns}.
   * \param is_AB_to_AX specifies which hadron to excite into a string.
   * true : A + B -> A + X
   * false : A + B -> X + B
   * \return whether the process is successfully implemented.
   */
  bool next_SDiff(bool is_AB_to_AX);
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
   * Soft Non-diffractive process
   * is modelled in accordance with dual-topological approach
   * \iref{Capella:1978ig}.
   * This involves a parton exchange in conjunction with momentum transfer.
   * Probability distribution function of the lightcone momentum fraction
   * carried by quark is based on the UrQMD model
   * \iref{Bass:1998ca,Bleicher:1999xi}.
   * \return whether the process is successfully implemented.
   */
  bool next_NDiffSoft();
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
    return Pythia8::Vec4(mom.x1(), mom.x2(), mom.x3(), energy);
  }

  /**
   * perform string fragmentation to determine species and momenta of hadrons
   * by implementing PYTHIA 8.2 \iref{Andersson:1983ia,Sjostrand:2014zea}.
   * \param idq1 is PDG id of quark or anti-diquark (carrying color index).
   * \param idq2 is PDG id of diquark or anti-quark (carrying anti-color index).
   * \param mString is the string mass.
   * \param evecLong is unit 3-vector specifying the direction of diquark or
   * anti-diquark.
   * \param random_rotation is whether or not we randomly rotate the
   * orientation.
   * \return number of hadrons fragmented out of string.
   */
  int fragment_string(int idq1, int idq2, double mString, ThreeVector &evecLong,
                      bool random_rotation);
};

}  // namespace Smash

#endif  // SRC_INCLUDE_STRINGPROCESS_H_
