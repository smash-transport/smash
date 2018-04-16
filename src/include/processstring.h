/*
 *
 *    Copyright (c) 2017-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_PROCESSSTRING_H_
#define SRC_INCLUDE_PROCESSSTRING_H_

#include <memory>
#include <utility>
#include <vector>

#include "particledata.h"
#include "Pythia8/Pythia.h"

namespace smash {

/**
 * StringSoft Types are used to identify the type of the soft string subprocess
 * (-1) nothing (None)
 * (0)  single diffractive A+B to A+X (SingleDiffAX)
 * (1)  single diffractive A+B to X+B (SingleDiffXB)
 * (2)  double diffractive (DoubleDiff)
 * (3)  soft non-diffractive (NonDiff)
 */
enum class StringSoftType {
  None = -1,
  SingleDiffAX = 0,
  SingleDiffXB = 1,
  DoubleDiff = 2,
  NonDiff = 3
};

// \todo Sangwook: make file (processstring) and class (StringProcess) name
// consistent

/**
 * \brief String excitation processes used in SMASH
 *
 * Only one instance of this class should be created.
 *
 * This class implements string excitation processes based on the UrQMD model
 * \iref{Bass:1998ca}, \iref{Bleicher:1999xi} and subsequent fragmentation
 * according to the LUND/PYTHIA fragmentation scheme
 * \iref{Andersson:1983ia}, \iref{Sjostrand:2014zea}.
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
  /// forward lightcone momentum p^{+} of incoming particle A [GeV]
  double PPosA_;
  /// forward lightcone momentum p^{+} of incoming particle B [GeV]
  double PPosB_;
  /// backward lightcone momentum p^{-} of incoming particle A [GeV]
  double PNegA_;
  /// backward lightcone momentum p^{-} of incoming particle B [GeV]
  double PNegB_;
  /// masses of incoming particles [GeV]
  double massA_, massB_;
  /// sqrt of Mandelstam variable s of collision [GeV]
  double sqrtsAB_;
  /// PdgCodes of incoming particles
  std::array<PdgCode, 2> PDGcodes_;
  /// momenta of incoming particles in the lab frame [GeV]
  std::array<FourVector, 2> plab_;
  /// momenta of incoming particles in the center of mass frame [GeV]
  std::array<FourVector, 2> pcom_;
  /// velocity four vector of the center of mass in the lab frame
  FourVector ucomAB_;
  /// velocity three vector of the center of mass in the lab frame
  ThreeVector vcomAB_;
  /**
   * Orthonormal basis vectors in the center of mass frame,
   * where the 0th one is parallel to momentum of incoming particle A
   */
  std::array<ThreeVector, 3> evecBasisAB_;
  /// total number of final state particles
  int NpartFinal_;
  /// soft subprocess identifier
  StringSoftType subproc_;
  /// number of particles fragmented from strings
  std::array<int, 2> NpartString_;
  /// the minimum lightcone momentum scale carried by gluon [GeV]
  double pmin_gluon_lightcone_;
  /**
   * parameter \f$\beta\f$ for the gluon distribution function
   * \f$ P(x) = x^{-1} (1 - x)^{1 + \beta} \f$
   */
  double pow_fgluon_beta_;
  /**
   * parameter \f$\alpha\f$ for the quark distribution function
   * \f$ P(x) = x^{\alpha - 1} (1 - x)^{\beta - 1} \f$
   */
  double pow_fquark_alpha_;
  /**
   * parameter \f$\beta\f$ for the quark distribution function
   * \f$ P(x) = x^{\alpha - 1} (1 - x)^{\beta - 1} \f$
   */
  double pow_fquark_beta_;
  /**
   * Transverse momentum spread of the excited strings. [GeV]
   * Transverse momenta of strings are sampled according to gaussian
   * distribution with width sigma_qperp_
   */
  double sigma_qperp_;
  /// string tension [GeV/fm]
  double kappa_tension_string_;
  /// constant proper time in the case of constant formation time [fm]
  double time_formation_const_;
  /// time of collision in the computational frame [fm]
  double time_collision_;
  /// Lorentz gamma factor of center of mass in the computational frame
  double gamma_factor_com_;
  /// square root of 2 (\f$\sqrt{2}\f$)
  double sqrt2_;

  /**
   * final state array
   * which must be accessed after the collision
   */
  ParticleList final_state_;

  /// PYTHIA object used in hard string routine
  std::unique_ptr<Pythia8::Pythia> pythia_parton_;

  /// PYTHIA object used in fragmentation
  std::unique_ptr<Pythia8::Pythia> pythia_hadron_;

  /// An object to compute cross-sections
  Pythia8::SigmaTotal pythia_sigmatot_;

 public:
  // clang-format off

  /**
   * Constructor, initializes pythia. Should only be called once.
   * \param[in] string_tension value of kappa_tension_string_ [GeV/fm]
   * \param[in] time_formation value of time_formation_const_ [fm]
   * \param[in] gluon_beta value of pow_fgluon_beta_
   * \param[in] gluon_pmin value of pmin_gluon_lightcone_
   * \param[in] quark_alpha value of pow_fquark_alpha_
   * \param[in] quark_beta value of pow_fquark_beta_
   * \param[in] strange_supp strangeness suppression factor
   *        (StringFlav:probStoUD) in fragmentation
   * \param[in] diquark_supp diquark suppression factor
   *        (StringFlav:probQQtoQ) in fragmentation
   * \param[in] sigma_perp value of sigma_qperp_ [GeV]
   * \param[in] stringz_a parameter (StringZ:aLund)
   *        for the fragmentation function
   * \param[in] stringz_b parameter (StringZ:bLund)
   *        for the fragmentation function
   * \param[in] string_sigma_T transverse momentum spread (StringPT:sigma)
   *        in fragmentation [GeV]
   *
   * \see StringProcess::common_setup_pythia(Pythia8::Pythia *,
   *                     double, double, double, double, double)
   * \see StringProcess::kappa_tension_string_
   * \see StringProcess::time_formation_const_
   * \see StringProcess::pow_fgluon_beta_
   * \see StringProcess::pmin_gluon_lightcone_
   * \see StringProcess::pow_fquark_alpha_
   * \see StringProcess::pow_fquark_beta_
   * \see StringProcess::sigma_qperp_
   * \see 3rdparty/pythia8230/share/Pythia8/xmldoc/FlavourSelection.xml
   * \see 3rdparty/pythia8230/share/Pythia8/xmldoc/Fragmentation.xml
   * \see 3rdparty/pythia8230/share/Pythia8/xmldoc/MasterSwitches.xml
   * \see 3rdparty/pythia8230/share/Pythia8/xmldoc/MultipartonInteractions.xml
   */
  StringProcess(double string_tension, double time_formation,
                double gluon_beta, double gluon_pmin,
                double quark_alpha, double quark_beta,
                double strange_supp, double diquark_supp,
                double sigma_perp,
                double stringz_a, double stringz_b,
                double string_sigma_T);

  /**
   * Common setup of PYTHIA objects for soft and hard string routines
   * \param[out] pythia_in pointer to the PYTHIA object
   * \param[in] strange_supp strangeness suppression factor
   *        (StringFlav:probStoUD) in fragmentation
   * \param[in] diquark_supp diquark suppression factor
   *        (StringFlav:probQQtoQ) in fragmentation
   * \param[in] stringz_a parameter (StringZ:aLund)
   *        for the fragmentation function
   * \param[in] stringz_b parameter (StringZ:bLund)
   *        for the fragmentation function
   * \param[in] string_sigma_T transverse momentum spread (StringPT:sigma)
   *        in fragmentation [GeV]
   *
   * \see 3rdparty/pythia8230/share/Pythia8/xmldoc/FlavourSelection.xml
   * \see 3rdparty/pythia8230/share/Pythia8/xmldoc/Fragmentation.xml
   */
  void common_setup_pythia(Pythia8::Pythia *pythia_in, double strange_supp,
                           double diquark_supp, double stringz_a,
                           double stringz_b, double string_sigma_T);

  // clang-format on

  /**
   * Function to get the PYTHIA object for hard string routine
   * \return pointer to the PYTHIA object used in hard string routine
   */
  Pythia8::Pythia *get_ptr_pythia_parton() { return pythia_parton_.get(); }

  /**
   * Interface to pythia_sigmatot_ to compute cross-sections of A+B->
   * different final states \iref{Schuler:1993wr}.
   * \param[in] pdg_a pdg code of incoming particle A
   * \param[in] pdg_b pdg code of incoming particle B
   * \param[in] sqrt_s collision energy in the center of mass frame [GeV]
   * \return array with single diffractive cross-sections AB->AX, AB->XB and
   * double diffractive AB->XX.
   */
  std::array<double, 3> cross_sections_diffractive(int pdg_a, int pdg_b,
                                                   double sqrt_s) {
    // This threshold magic is following Pythia. Todo(ryu): take care of this.
    double sqrts_threshold = 2. * (1. + 1.0e-6);
    pdg_a = std::abs(pdg_a);
    pdg_b = std::abs(pdg_b);
    /* In the case of mesons, the corresponding vector meson masses
     * are used to evaluate the energy threshold. */
    const int pdg_a_mod = (pdg_a > 1000) ? pdg_a : 10 * (pdg_a / 10) + 3;
    const int pdg_b_mod = (pdg_b > 1000) ? pdg_b : 10 * (pdg_b / 10) + 3;
    sqrts_threshold += pythia_hadron_->particleData.m0(pdg_a_mod) +
                       pythia_hadron_->particleData.m0(pdg_b_mod);
    /* Constant cross-section for sub-processes below threshold equal to
     * cross-section at the threshold. */
    if (sqrt_s < sqrts_threshold) {
      sqrt_s = sqrts_threshold;
    }
    pythia_sigmatot_.calc(pdg_a, pdg_b, sqrt_s);
    return {pythia_sigmatot_.sigmaAX(), pythia_sigmatot_.sigmaXB(),
            pythia_sigmatot_.sigmaXX()};
  }

  /**
   * \todo The following set_ functions are replaced with
   * constructor with arguments.
   * Must be cleaned up if necessary.
   */

  /**
   * set the minimum lightcone momentum scale carried by gluon.
   * This is relevant for the double-diffractive process.
   * The minimum lightcone momentum fraction is set to be
   * pmin_gluon_lightcone_/sqrtsAB.
   * \param p_light_cone_min a value that we want to use for
   *                         pmin_gluon_lightcone_.
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
   * Set the soft subprocess identifier
   * \param[in] iproc soft string subprocess that will be implemented
   */
  void set_subproc(StringSoftType iproc) { subproc_ = iproc; }
  /// \return the soft subprocess identifier
  StringSoftType get_subproc() { return subproc_; }
  /**
   * initialization
   * feed intial particles, time of collision and gamma factor of the center of
   * mass.
   * \param[in] incoming is the list of initial state particles.
   * \param[in] tcoll is time of collision.
   * \param[in] gamma gamma factor of the center of mass.
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
   * Determine string masses and directions in which strings are stretched
   * \param[in] quarks pdg ids of string ends
   * \param[in] pstr_com 4-momenta of strings in the C.o.m. frame
   * \param[out] m_str masses of strings
   * \param[out] evec_str are directions in which strings are stretched.
   * \return whether masses are above the threshold
   */
  bool set_mass_and_direction_2strings(
      const std::array<std::array<int, 2>, 2> &quarks,
      const std::array<FourVector, 2> &pstr_com,
      std::array<double, 2> &m_str,
      std::array<ThreeVector, 2> &evec_str);
  /**
   * Prepare kinematics of two strings, fragment them and append to final_state
   * \param[in] quarks pdg ids of string ends
   * \param[in] pstr_com 4-momenta of strings in the C.o.m. frame
   * \param[in] m_str masses of strings
   * \param[out] evec_str are directions in which strings are stretched.
   * \param[in] flip_string_ends is whether or not we randomly switch string ends.
   * \return whether fragmentations and final state creation was successful
   */
  bool make_final_state_2strings(
      const std::array<std::array<int, 2>, 2> &quarks,
      const std::array<FourVector, 2> &pstr_com,
      const std::array<double, 2> &m_str,
      const std::array<ThreeVector, 2> &evec_str,
      bool flip_string_ends);

  /**
   * Single-diffractive process
   * is based on single pomeron exchange described in \iref{Ingelman:1984ns}.
   * \param[in] is_AB_to_AX specifies which hadron to excite into a string.
   *                    true : A + B -> A + X
   *                    false : A + B -> X + B
   * \return whether the process is successfully implemented.
   */
  bool next_SDiff(bool is_AB_to_AX);
  /**
   * Double-diffractive process ( A + B -> X + X )
   * is similar to the single-diffractive process,
   * but lightcone momenta of gluons are sampled
   * in the same was as the UrQMD model \iref{Bass:1998ca},
   * \iref{Bleicher:1999xi}.
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
   * \iref{Bass:1998ca}, \iref{Bleicher:1999xi}.
   * \return whether the process is successfully implemented.
   */
  bool next_NDiffSoft();
  /**
   * Hard Non-diffractive process
   * is based on PYTHIA 8 with partonic showers and interactions.
   * \return whether the process is successfully implemented.
   */
  bool next_NDiffHard();
  /**
   * Baryon-antibaryon annihilation process
   * Based on what UrQMD \iref{Bass:1998ca}, \iref{Bleicher:1999xi} does,
   * it create two mesonic strings after annihilating one quark-antiquark pair.
   * Each string has mass equal to half of sqrts.
   * \return whether the process is successfully implemented.
   */
  bool next_BBbarAnn();

  /**
   * a function to get the final state particle list
   * which is called after the collision
   * \return ParticleList filled with the final state particles.
   */
  ParticleList get_final_state() { return final_state_; }

  /**
   * compute the formation time and fill the arrays with final-state particles
   * as described in \iref{Andersson:1983ia}.
   * \param[in] uString is velocity four vector of the string.
   * \param[in] evecLong is unit 3-vector in which string is stretched.
   * \return number of hadrons fragmented out of string.
   */
  int append_final_state(const FourVector &uString,
                         const ThreeVector &evecLong);

  /**
   * convert Kaon-L or Kaon-S into K0 or Anti-K0
   * \param pythia_id is PDG id to be converted.
   */
  static void convert_KaonLS(int &pythia_id) {
    if (pythia_id == 310 || pythia_id == 130) {
      pythia_id = (Random::uniform_int(0, 1) == 0) ? 311 : -311;
    }
  }

  /**
   * Construct diquark from two quarks. Order does not matter.
   * \param[in] q1 PDG code of quark 1
   * \param[in] q2 PDG code of quark 2
   * \return PDG code of diquark composed of q1 and q2
   */
  static int diquark_from_quarks(int q1, int q2);

  /**
   * make a random selection to determine partonic contents at the string ends.
   * \param[in] pdgcode_in is PdgCode of hadron which transforms into a string.
   * \param[out] idq1 is PDG id of quark or anti-diquark.
   * \param[out] idq2 is PDG id of anti-quark or diquark.
   */
  static void make_string_ends(const PdgCode &pdgcode_in, int &idq1, int &idq2);

  /**
   * Easy setter of Pythia Vec4 from SMASH
   * \param[in] energy time component
   * \param[in] mom spatial three-vector
   * \return Pythia Vec4 from energy and ThreeVector
   */
  static Pythia8::Vec4 set_Vec4(double energy, const ThreeVector &mom) {
    return Pythia8::Vec4(mom.x1(), mom.x2(), mom.x3(), energy);
  }

  // clang-format off

  /**
   * perform string fragmentation to determine species and momenta of hadrons
   * by implementing PYTHIA 8.2 \iref{Andersson:1983ia}, \iref{Sjostrand:2014zea}.
   * \param[in] idq1 PDG id of quark or anti-diquark (carrying color index).
   * \param[in] idq2 PDG id of diquark or anti-quark (carrying anti-color index).
   * \param[in] mString the string mass.
   * \param[out] evecLong unit 3-vector specifying the direction of diquark or
   *                      anti-diquark.
   * \param[in] flip_string_ends is whether or not we randomly switch string ends.
   * \return number of hadrons fragmented out of string.
   */
  int fragment_string(int idq1, int idq2, double mString,
                      ThreeVector &evecLong, bool flip_string_ends);

  /**
   * Assign a cross section scaling factor to all outgoing particles.
   * Factor is only non-zero, when the outgoing particle carries
   * a valence quark from the excited hadron.
   */
  static void assign_all_scaling_factors(int baryon_string,
                                         ParticleList& outgoing_particles,
                                         ThreeVector &evec_coll,
                                         double suppression_factor);

  /**
   * Find the first particle, which can carry nq1, and the last particle,
   * which can carry nq2 valence quarks and return their indices in
   * the given list.
   */
  static std::pair<int, int> find_leading(int nq1, int nq2, ParticleList& list);

  /**
   * Assign a cross section scaling factor to the given particle.
   * The scaling factor is the number of quarks from the excited hadron,
   * that the fragment carries devided by the total number of quarks in
   * this fragment.
   */
  static void assign_scaling_factor(int nquark, ParticleData& data,
                                    double suppression_factor);

  // clang-format on
};

}  // namespace smash

#endif  // SRC_INCLUDE_PROCESSSTRING_H_
