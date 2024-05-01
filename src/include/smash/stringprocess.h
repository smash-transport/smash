/*
 *
 *    Copyright (c) 2017-2020,2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_STRINGPROCESS_H_
#define SRC_INCLUDE_SMASH_STRINGPROCESS_H_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Pythia8/Pythia.h"

#include "constants.h"
#include "logging.h"
#include "particledata.h"

namespace smash {
static constexpr int LPythia = LogArea::Pythia::id;

/**
 * \brief String excitation processes used in SMASH
 *
 * Only one instance of this class should be created.
 *
 * This class implements string excitation processes based on the UrQMD model
 * \iref{Bass:1998ca}, \iref{Bleicher:1999xi} and subsequent fragmentation
 * according to the LUND/PYTHIA fragmentation scheme
 * \iref{Andersson:1983ia}, \iref{Sjostrand:2014zea}, \iref{Bierlich:2022pfr}.
 *
 * The class implements the following functionality:
 * - given two colliding initial state particles it provides hadronic final
 *   state after single diffractive, double diffractive and non-diffractive
 *   string excitation
 * - owns a Pythia8::SigmaTotal object to compute corresponding cross-sections
 * - owns a Pythia object, that allows to fragment strings
 */
class StringProcess {
 private:
  // The following 4 variables are in the center of mass frame
  /// forward lightcone momentum p^{+} of incoming particle A in CM-frame [GeV]
  double PPosA_;
  /// forward lightcone momentum p^{+} of incoming particle B in CM-frame [GeV]
  double PPosB_;
  /// backward lightcone momentum p^{-} of incoming particle A in CM-frame [GeV]
  double PNegA_;
  /// backward lightcone momentum p^{-} of incoming particle B in CM-frame [GeV]
  double PNegB_;
  /// mass of incoming particle A [GeV]
  double massA_;
  /// mass of incoming particle B [GeV]
  double massB_;
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
  /// number of particles fragmented from strings
  std::array<int, 2> NpartString_;
  /// the minimum lightcone momentum scale carried by a gluon [GeV]
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
  /**
   * parameter (StringZ:aLund) for the fragmentation function
   * of leading baryon in soft non-diffractive string processes
   */
  double stringz_a_leading_;
  /**
   * parameter (StringZ:bLund) for the fragmentation function
   * of leading baryon in soft non-diffractive string processes
   */
  double stringz_b_leading_;
  /**
   * parameter (StringZ:aLund) for the fragmentation function
   * of other (produced) hadrons in soft non-diffractive string processes
   */
  double stringz_a_produce_;
  /**
   * parameter (StringZ:bLund) for the fragmentation function
   * of other (produced) hadrons in soft non-diffractive string processes
   */
  double stringz_b_produce_;
  /// strange quark suppression factor
  double strange_supp_;
  /// diquark suppression factor
  double diquark_supp_;
  /// popcorn rate
  double popcorn_rate_;
  /// transverse momentum spread in string fragmentation
  double string_sigma_T_;
  /// string tension [GeV/fm]
  double kappa_tension_string_;
  /**
   * additional cross-section suppression factor
   * to take coherence effect into account.
   */
  double additional_xsec_supp_;
  /// constant proper time in the case of constant formation time [fm]
  double time_formation_const_;
  /// factor to be multiplied to formation times in soft strings
  double soft_t_form_;
  /// time of collision in the computational frame [fm]
  double time_collision_;
  /**
   * Whether the formation time should depend on the mass of the fragment
   * according to \iref{Andersson:1983ia} eq. 2.45:
   *
   * \f$ \tau = \sqrt{2}\frac{m}{\kappa} \f$
   *
   * The formation time and position is not calculated directly using the yoyo
   * model because the spacetime rapidity where a string fragment forms is not
   * equal to the fragment's momentum space rapidity. This cannot be easily
   * combined with possible interactions before the formation time.
   **/
  bool mass_dependent_formation_times_;
  /**
   * Probability of splitting a nucleon into the quark flavour it has only
   * once and a diquark it has twice.
   */
  double prob_proton_to_d_uu_;

  /// Whether to use a separate fragmentation function for leading baryons.
  bool separate_fragment_baryon_;

  /**
   * Whether to use the monash tune \iref{Skands:2014pea} for all string
   * processes.
   */
  bool use_monash_tune_;

  /**
   * final state array
   * which must be accessed after the collision
   */
  ParticleList final_state_;

  /**
   * Map containing PYTHIA objects for hard string routines.
   * Particle IDs are used as the keys to obtain the respective object.
   * This was introduced to reduce the amount of Pythia init() calls.
   */
  typedef std::map<std::pair<int, int>, std::unique_ptr<Pythia8::Pythia>>
      pythia_map;

  /// Map object to contain the different pythia objects
  pythia_map hard_map_;

  /// PYTHIA object used in fragmentation
  std::unique_ptr<Pythia8::Pythia> pythia_hadron_;

  /// An object to compute cross-sections
  Pythia8::SigmaTotal pythia_sigmatot_;

  /**
   * An object for the flavor selection in string fragmentation
   * in the case of separate fragmentation function for leading baryon
   */
  Pythia8::StringFlav pythia_stringflav_;

  /**
   * event record for intermediate partonic state
   * in the hard string routine
   */
  Pythia8::Event event_intermediate_;

 public:
  // clang-format off

  /**
   * Constructor, initializes pythia. Should only be called once.
   * \param[in] string_tension value of #kappa_tension_string_ [GeV/fm]
   * \param[in] time_formation value of #time_formation_const_ [fm]
   * \param[in] gluon_beta value of #pow_fgluon_beta_
   * \param[in] gluon_pmin value of #pmin_gluon_lightcone_
   * \param[in] quark_alpha value of #pow_fquark_alpha_
   * \param[in] quark_beta value of #pow_fquark_beta_
   * \param[in] strange_supp strangeness suppression factor
   *        (StringFlav:probStoUD) in fragmentation
   * \param[in] diquark_supp diquark suppression factor
   *        (StringFlav:probQQtoQ) in fragmentation
   * \param[in] sigma_perp value of #sigma_qperp_ [GeV]
   * \param[in] stringz_a_leading Parameter a in Lund fragmentation function
   *            for leading baryons.
   * \param[in] stringz_b_leading Parameter b in Lund fragmentation function
   *            for leading baryons.
   * \param[in] stringz_a parameter (StringZ:aLund)
   *        for the fragmentation function
   * \param[in] stringz_b parameter (StringZ:bLund)
   *        for the fragmentation function [GeV^-2]
   * \param[in] string_sigma_T transverse momentum spread (StringPT:sigma)
   *        in fragmentation [GeV]
   * \param[in] factor_t_form to be multiplied to soft string formation times
   * \param[in] mass_dependent_formation_times Whether the formation times of
   *            string fragments should depend on their mass.
   * \param[in] prob_proton_to_d_uu Probability of a nucleon to be split into
   *            the quark it contains once and a diquark another flavour.
   * \param[in] separate_fragment_baryon whether to use a separate
   *            fragmentation function for leading baryons in non-diffractive
   *            string processes.
   * \param[in] popcorn_rate parameter (StringFlav:popcornRate)
   *        to determine the production rate of popcorn mesons from
   *        the diquark end of a string.
   * \param[in] use_monash_tune whether to use the monash tune for all string
   *            processes. This is recommended if one runs smash at LHC
   *            energies
   *
   * \see StringProcess::common_setup_pythia(Pythia8::Pythia *,
   *                     double, double, double, double, double)
   * \see pythia8302/share/Pythia8/xmldoc/FlavourSelection.xml
   * \see pythia8302/share/Pythia8/xmldoc/Fragmentation.xml
   * \see pythia8302/share/Pythia8/xmldoc/MasterSwitches.xml
   * \see pythia8302/share/Pythia8/xmldoc/MultipartonInteractions.xml
   */
  StringProcess(double string_tension, double time_formation,
                double gluon_beta, double gluon_pmin,
                double quark_alpha, double quark_beta,
                double strange_supp, double diquark_supp,
                double sigma_perp, double stringz_a_leading,
                double stringz_b_leading, double stringz_a,
                double stringz_b,  double string_sigma_T,
                double factor_t_form,
                bool mass_dependent_formation_times,
                double prob_proton_to_d_uu,
                bool separate_fragment_baryon, double popcorn_rate,
                bool use_monash_tune);

  /**
   * Common setup of PYTHIA objects for soft and hard string routines
   * \param[out] pythia_in pointer to the PYTHIA object
   * \param[in] strange_supp strangeness suppression factor
   *        (StringFlav:probStoUD) in fragmentation
   * \param[in] diquark_supp diquark suppression factor
   *        (StringFlav:probQQtoQ) in fragmentation
   * \param[in] popcorn_rate parameter (StringFlav:popcornRate)
   *        to determine the production rate of popcorn mesons from
   *        the diquark end of a string.
   * \param[in] stringz_a parameter (StringZ:aLund)
   *        for the fragmentation function
   * \param[in] stringz_b parameter (StringZ:bLund)
   *        for the fragmentation function
   * \param[in] string_sigma_T transverse momentum spread (StringPT:sigma)
   *        in fragmentation [GeV]
   *
   * \see pythia8302/share/Pythia8/xmldoc/FlavourSelection.xml
   * \see pythia8302/share/Pythia8/xmldoc/Fragmentation.xml
   */
  void common_setup_pythia(Pythia8::Pythia *pythia_in, double strange_supp,
                           double diquark_supp, double popcorn_rate,
                           double stringz_a, double stringz_b,
                           double string_sigma_T);

  /**
   * Set PYTHIA random seeds to be desired values.
   * The value is recalculated such that it is allowed by PYTHIA.
   *
   * \see smash::maximum_rndm_seed_in_pythia
   */
  void init_pythia_hadron_rndm() {
    const int seed_new =
        random::uniform_int(1, maximum_rndm_seed_in_pythia);

    pythia_hadron_->rndm.init(seed_new);
    logg[LPythia].debug("pythia_hadron_ : rndm is initialized with seed ",
                        seed_new);
  }

  // clang-format on

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
    /* In the case of mesons, the corresponding vector meson masses
     * are used to evaluate the energy threshold. */
    const int pdg_a_mod =
        (std::abs(pdg_a) > 1000) ? pdg_a : 10 * (std::abs(pdg_a) / 10) + 3;
    const int pdg_b_mod =
        (std::abs(pdg_b) > 1000) ? pdg_b : 10 * (std::abs(pdg_b) / 10) + 3;
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

  // clang-format off

  /**
   * initialization
   * feed intial particles, time of collision and gamma factor of the center of
   * mass.
   * \param[in] incoming is the list of initial state particles.
   * \param[in] tcoll is time of collision.
   */
  void init(const ParticleList &incoming, double tcoll);
  /**
   * compute three orthonormal basis vectors from unit vector
   * in the longitudinal direction
   * \param[in] evec_polar unit three-vector in the longitudinal direction
   * \param[out] evec_basis orthonormal basis vectors of which
   *             evec_basis[0] is in the longitudinal direction while
   *             evec_basis[1] and evec_basis[2] span the transverse plane.
   */
  static void make_orthonormal_basis(ThreeVector &evec_polar,
                                     std::array<ThreeVector, 3> &evec_basis);
  /**
   * compute the lightcone momenta of incoming particles
   * where the longitudinal direction is set to be same
   * as that of the three-momentum of particle A.
   */
  void compute_incoming_lightcone_momenta();
  /**
   * Determine string masses and directions in which strings are stretched
   * \param[in] quarks pdg ids of string ends
   * \param[in] pstr_com 4-momenta of strings in the C.o.m. frame [GeV]
   * \param[out] m_str masses of strings [GeV]
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
   * \param[in] pstr_com 4-momenta of strings in the C.o.m. frame [GeV]
   * \param[in] m_str masses of strings [GeV]
   * \param[out] evec_str are directions in which strings are stretched.
   * \param[in] flip_string_ends is whether or not we randomly switch string ends.
   * \param[in] separate_fragment_baryon is whether to fragment leading baryons
   *            (or anti-baryons) with a separate fragmentation function.
   * \return whether fragmentations and final state creation was successful
   */
  bool make_final_state_2strings(
      const std::array<std::array<int, 2>, 2> &quarks,
      const std::array<FourVector, 2> &pstr_com,
      const std::array<double, 2> &m_str,
      const std::array<ThreeVector, 2> &evec_str,
      bool flip_string_ends, bool separate_fragment_baryon);

  /**
   * Single-diffractive process
   * is based on single pomeron exchange described in \iref{Ingelman:1984ns}.
   * \param[in] is_AB_to_AX specifies which hadron to excite into a string.
   *            true : A + B -> A + X,
   *            false : A + B -> X + B
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
   *
   * \throw std::runtime_error
   *        if incoming particles are neither mesonic nor baryonic
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
   *
   * \throw std::invalid_argument
   *        if incoming particles are not baryon-antibaryon pair
   */
  bool next_BBbarAnn();

  /**
   * Compare the valence quark contents of the actual and mapped hadrons and
   * evaluate how many more constituents the actual hadron has compared to the
   * mapped one.
   * excess_quark[i - 1] is how many more quarks with flavor i (PDG id i)
   * pdg_actual has compared to pdg_mapped.
   * excess_antiq[i - 1] is how many more antiquarks with flavor i (PDG id -i)
   * pdg_actual has compared to pdg_mapped.
   *
   * \param[in] pdg_actual PDG code of actual incoming particle.
   * \param[in] pdg_mapped PDG code of mapped particles used in PYTHIA
   *            event generation.
   * \param[out] excess_quark excess of quarks.
   * \param[out] excess_antiq excess of anti-quarks.
   */
  static void find_excess_constituent(PdgCode &pdg_actual, PdgCode &pdg_mapped,
                                      std::array<int, 5> &excess_quark,
                                      std::array<int, 5> &excess_antiq);
  /**
   * Convert a partonic PYTHIA particle into the desired species
   * and update the excess of constituents.
   * If the quark flavor i is converted into another flavor j,
   * excess_constituent[i - 1] increases by 1 and
   * excess_constituent[j - 1] decreases by 1.
   * Note that this happens only if
   * excess_constituent[i - 1] < 0 and excess_constituent[j - 1] > 0
   * (i.e., the incoming hadron has more constituents with flavor j
   *  and less constituents with flavor i, compared to the mapped hadron),
   * so they get closer to 0 after the function call.
   *
   * \param[out] particle PYTHIA particle object to be converted.
   * \param[out] excess_constituent excess in the number of quark constituents.
   *             If the particle has positive (negative) quark number,
   *             excess of quarks (anti-quarks) should be used.
   *
   * \see StringProcess::restore_constituent(Pythia8::Event &,
   *                         std::array<std::array<int, 5>, 2> &,
   *                         std::array<std::array<int, 5>, 2> &)
   */
  void replace_constituent(Pythia8::Particle &particle,
                           std::array<int, 5> &excess_constituent);

  /**
   * Compute how many quarks and antiquarks we have in the system,
   * and update the correspoing arrays with size 5.
   * Note that elements of the array (0, 1, 2, 3, 4) correspond
   * to d, u, s, c, b flavors.
   *
   * \param[in] event_intermediate PYTHIA partonic event record
   *            which contains output from PYTHIA (hard) event generation.
   * \param[out] nquark_total total number of quarks in the system.
   *             This is computed based on event_intermediate.
   * \param[out] nantiq_total total number of antiquarks in the system.
   *             This is computed based on event_intermediate.
   */
  void find_total_number_constituent(Pythia8::Event &event_intermediate,
                                     std::array<int, 5> &nquark_total,
                                     std::array<int, 5> &nantiq_total);

  /**
   * Take total number of quarks and check if the system has
   * enough constituents that need to be converted into other flavors.
   * If that is not the case, a gluon is splitted into a quark-antiquark pair
   * with desired flavor, so that their flavor can be changed afterwards.
   * For example, if there is no antiquark in the system and we have
   * excess_antiq = (1, -1, 0, 0, 0)
   * (i.e., one ubar has to be converted into dbar),
   * a gluon will be splitted into u-ubar pair.
   *
   * \param[out] event_intermediate PYTHIA partonic event record to be updated
   *             when a gluon happens to split into a qqbar pair.
   * \param[out] nquark_total total number of quarks in the system.
   *             This is computed based on event_intermediate.
   * \param[out] nantiq_total total number of antiquarks in the system.
   *             This is computed based on event_intermediate.
   * \param[in] sign_constituent true (false)
   *            if want to check quarks (antiquarks) and their excesses.
   * \param[in] excess_constituent excess in the number of quark constituents.
   *            If sign_constituent is true (false),
   *            excess of quarks (anti-quarks) should be used.
   * \return false if there are not enough constituents and there is no gluon
   *         to split into desired quark-antiquark pair.
   *         Otherwise, it gives true.
   */
  bool splitting_gluon_qqbar(Pythia8::Event &event_intermediate,
      std::array<int, 5> &nquark_total, std::array<int, 5> &nantiq_total,
      bool sign_constituent,
      std::array<std::array<int, 5>, 2> &excess_constituent);

  /**
   * Take total number of quarks and check if the system has
   * enough constituents that need to be converted into other flavors.
   * If that is not the case, excesses of quarks and antiquarks are
   * modified such that the net quark number of each flavor is
   * conserved.
   * For example, if there is no antiquark in the system and we have
   * excess_antiq = (1, -1, 0, 0, 0)
   * (i.e., one ubar has to be converted into dbar),
   * excess_antiq will be changed into (0, 0, 0, 0, 0) and
   * (-1, 1, 0, 0, 0) will be added to excess_quark
   * (i.e., one d quark has to be converted into u quark instead).
   *
   * Number of quarks is checked if the first argument is
   * the total number of quarks, and the second and third arguments are
   * respectively excesses of quarks and antiquarks.
   * Number of antiquarks is checked if the first argument is
   * the total number of antiquarks, and the second and third arguments are
   * respectively excesses of antiquarks and quarks.
   *
   * \param[in] nquark_total total number of quarks (antiquarks)
   *            in the system.
   * \param[out] excess_quark excess of quarks (antiquarks)
   *             in incoming particles, compared to the mapped ones.
   * \param[out] excess_antiq excess of anti-quarks (quarks)
   *             in incoming particles, compared to the mapped ones.
   */
  void rearrange_excess(std::array<int, 5> &nquark_total,
                        std::array<std::array<int, 5>, 2> &excess_quark,
                        std::array<std::array<int, 5>, 2> &excess_antiq);

  /**
   * Take the intermediate partonic state from PYTHIA event with mapped hadrons
   * and convert constituents into the desired ones according to the excess of
   * quarks and anti-quarks.
   * Quark (antiquark) flavor is changed and excess of quark (antiquark)
   * is also updated by calling StringProcess::replace_constituent.
   * Beginning with the most forward (or backward) constituent,
   * conversion is done until the total net quark number of each flavor
   * is same with that of incoming hadrons.
   * (i.e., excess_quark minus excess_antiq of incoming hadrons becomes zero.)
   *
   * However, note that there are some circumstances where
   * this procedure is not directly carried out.
   * For example, a proton-kaon(+) collision mapped onto a proton-pion(+)
   * might be an issue if it involves d + dbar -> g g partonic interaction,
   * given that we anticipate to change dbar to sbar.
   * If such case occurs, we first try to split gluon into
   * quark-antiquark pair with desired flavor.
   * If there are not enough gluons to split, we try to modify the excesses
   * of constituents such that the net quark number is conserved.
   *
   * \param[out] event_intermediate PYTHIA partonic event record to be updated
   *             according to the valence quark contents of incoming hadrons.
   * \param[out] excess_quark excess of quarks
   *             in incoming particles, compared to the mapped ones.
   * \param[out] excess_antiq excess of anti-quarks
   *             in incoming particles, compared to the mapped ones.
   *
   * \see StringProcess::replace_constituent(Pythia8::Particle &,
   *                                         std::array<int, 5> &)
   * \see StringProcess::splitting_gluon_qqbar(Pythia8::Event &,
   *     std::array<int, 5> &, std::array<int, 5> &,
   *     bool, std::array<std::array<int, 5>, 2> &)
   * \see StringProcess::rearrange_excess(std::array<int, 5> &,
   *                                      std::array<std::array<int, 5>, 2> &,
   *                                      std::array<std::array<int, 5>, 2> &)
   */
  bool restore_constituent(Pythia8::Event &event_intermediate,
                           std::array<std::array<int, 5>, 2> &excess_quark,
                           std::array<std::array<int, 5>, 2> &excess_antiq);

  /**
   * Identify a set of partons, which are connected
   * to form a color-neutral string, from a given PYTHIA event record.
   * All partons found are moved into a new event record for the further
   * hadronization process.
   * Note that col and acol of Pythia8::Particle contain information
   * on the color flow.
   * This function begins with the most forward (or backward) parton.
   *
   * For example,
   * quark (col = 1, acol = 0), gluon (col = 2, acol = 1)
   * and antiquark (col = 0, acol = 2) correspond to
   * a \f$ \bar{q} \, g \, q \f$ mesonic string.
   * quark (col = 1, acol = 0) and diquark (col = 0, acol = 1) correspond to
   * a \f$ qq \, q\f$ baryonic string.
   *
   * \param[in] find_forward_string If it is set to be true (false),
   *                                it begins with forward (backward) parton.
   * \param[out] event_intermediate PYTHIA event record
   *                                from which a string is identified.
   *                                All partons found here are removed.
   * \param[out] event_hadronize PYTHIA event record
   *                             to which partons in a string are added.
   */
  void compose_string_parton(bool find_forward_string,
                             Pythia8::Event &event_intermediate,
                             Pythia8::Event &event_hadronize);
  /**
   * Identify a set of partons and junction(s), which are connected
   * to form a color-neutral string, from a given PYTHIA event record.
   * All partons found are moved into a new event record for the further
   * hadronization process.
   * Junction topology in PYTHIA combines three quarks (antiquarks)
   * to make a color-neutral baryonic (anti-baryonic) configuration.
   * A junction (anti-junction) carries three color (anti-color) indices
   * which are connected with quarks (antiquarks).
   * This function begins with the first junction.
   *
   * For example,
   * if there is a kind-1 junction with legs (col = 1, 2 and 3),
   * it will first look for three partons with color indices col = 1, 2 and 3
   * and trace color indices until each leg is ``closed'' with quark.
   * If there is no quark in the end, there should be an anti-junction
   * and its legs are connected to partons with corresponding anti-colors.
   *
   * \param[out] find_forward_string If it is set to be true (false),
   *                                 it is a string in the forward (backward)
   *                                 direction.
   * \param[out] event_intermediate PYTHIA event record
   *                                from which a string is identified.
   *                                All partons and junction(s) found here
   *                                are removed.
   * \param[out] event_hadronize PYTHIA event record
   *                             to which partons in a string are added.
   *
   * \see StringProcess::find_junction_leg(bool, std::vector<int> &,
   *                                       Pythia8::Event &, Pythia8::Event &)
   */
  void compose_string_junction(bool &find_forward_string,
                               Pythia8::Event &event_intermediate,
                               Pythia8::Event &event_hadronize);

  /**
   * Identify partons, which are associated with junction legs,
   * from a given PYTHIA event record.
   * All partons found are moved into a new event record for the further
   * hadronization process.
   * \param[in] sign_color true (false) if the junction is associated with
   *                       color (anti-color) indices, corresponding
   *                       baryonic (anti-baryonic) string
   * \param[out] col set of color indices that need to be found.
   *                 The value is set to be zero
   *                 once the corresponding partons are found.
   * \param[out] event_intermediate PYTHIA event record
   *                                from which a string is identified.
   *                                All partons and junction(s) found here
   *                                are removed.
   * \param[out] event_hadronize PYTHIA event record
   *                             to which partons in a string are added.
   *
   * \see StringProcess::compose_string_junction(bool &,
   *                         Pythia8::Event &, Pythia8::Event &)
   */
  void find_junction_leg(bool sign_color, std::vector<int> &col,
                         Pythia8::Event &event_intermediate,
                         Pythia8::Event &event_hadronize);

  /**
   * Obtain index of the most forward or backward particle
   * in a given PYTHIA event record.
   * \param[in] find_forward if it looks for the most forward
   *                         or backward particle.
   * \param[in] np_end number of the last particle entries to be excluded
   *                   in lookup. In other words, it finds the most forward
   *                   (or backward) particle among
   *                   event[1, ... , event.size() - 1 - np_end].
   * \param[in] event PYTHIA event record which contains particle entries.
   *                  Note that event[0] is reserved for information
   *                  on the entire system.
   * \return index of the selected particle,
   *         which is used to access the specific particle entry
   *         in the event record.
   */
  int get_index_forward(bool find_forward, int np_end,
                        Pythia8::Event &event) {
    int iforward = 1;
    for (int ip = 2; ip < event.size() - np_end; ip++) {
      const double y_quark_current = event[ip].y();
      const double y_quark_forward = event[iforward].y();
      if ((find_forward && y_quark_current > y_quark_forward) ||
          (!find_forward && y_quark_current < y_quark_forward)) {
        iforward = ip;
      }
    }
    return iforward;
  }

  /**
   * a function to get the final state particle list
   * which is called after the collision
   * \return ParticleList filled with the final state particles.
   */
  ParticleList get_final_state() { return final_state_; }

  /**
   * a function that clears the final state particle list
   * which is used for testing mainly
   */
  void clear_final_state() { final_state_.clear(); }

  /**
   * compute the formation time and fill the arrays with final-state particles
   * as described in \iref{Andersson:1983ia}.
   * \param[out] intermediate_particles list of fragmented particles
                 to be appended
   * \param[in] uString is velocity four vector of the string.
   * \param[in] evecLong is unit 3-vector in which string is stretched.
   * \return number of hadrons fragmented out of string.
   *
   * \throw std::invalid_argument if fragmented particle is not hadron
   * \throw std::invalid_argument if string is neither mesonic nor baryonic
   */
  int append_final_state(ParticleList &intermediate_particles,
                         const FourVector &uString,
                         const ThreeVector &evecLong);

  /**
   * append new particle from PYTHIA to a specific particle list
   * \param[in] pdgid PDG id of particle
   * \param[in] momentum four-momentum of particle
   * \param[out] intermediate_particles particle list to which
   *             the new particle is added.
   * \return whether PDG id exists in ParticleType table.
   */
  static bool append_intermediate_list(int pdgid, FourVector momentum,
                                       ParticleList &intermediate_particles) {
    const std::string s = std::to_string(pdgid);
    PdgCode pythia_code(s);
    ParticleTypePtr new_type = ParticleType::try_find(pythia_code);
    if (new_type) {
      ParticleData new_particle(ParticleType::find(pythia_code));
      new_particle.set_4momentum(momentum);
      new_particle.set_unpolarized_spin_vector();
      intermediate_particles.push_back(new_particle);
      return true;
    } else {
      // if the particle does not exist in SMASH the pythia event is rerun
      return false;
    }
  }

  /**
   * convert Kaon-L or Kaon-S into K0 or Anti-K0
   * \param[out] pythia_id is PDG id to be converted.
   */
  static void convert_KaonLS(int &pythia_id) {
    if (pythia_id == 310 || pythia_id == 130) {
      pythia_id = (random::uniform_int(0, 1) == 0) ? 311 : -311;
    }
  }

  /**
   * find two quarks from a diquark. Order does not matter.
   * \param[in] diquark PDG id of diquark
   * \param[out] q1 PDG id of quark 1
   * \param[out] q2 PDG id of quark 2
   * \param[out] deg_spin spin degeneracy
   */
  static void quarks_from_diquark(int diquark, int &q1, int &q2, int &deg_spin);

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
   * \param[in] xi probability to split a nucleon into the quark it has only
   *            once and a diquark of another flavour.
   */
  static void make_string_ends(const PdgCode &pdgcode_in, int &idq1, int &idq2,
                               double xi);

  /**
   * Easy setter of Pythia Vec4 from SMASH
   * \param[in] energy time component
   * \param[in] mom spatial three-vector
   * \return Pythia Vec4 from energy and ThreeVector
   */
  static Pythia8::Vec4 set_Vec4(double energy, const ThreeVector &mom) {
    return Pythia8::Vec4(mom.x1(), mom.x2(), mom.x3(), energy);
  }

  /**
   * compute the four-momentum properly oriented in the lab frame.
   * While PYTHIA assumes that the collision axis is in z-direction,
   * this is not necessarly the case in SMASH.
   * \param[in] particle particle object from PYTHIA event generation
   *            where z-axis is set to be the collision axis
   * \param[in] evec_basis three basis vectors in the lab frame
   *            evec_basis[0] is unit vector in the collision axis and
   *            other two span the transverse plane
   * \return four-momentum of particle in the lab frame
   */
  static FourVector reorient(Pythia8::Particle &particle,
                             std::array<ThreeVector, 3> &evec_basis) {
    ThreeVector three_momentum = evec_basis[0] * particle.pz() +
                                 evec_basis[1] * particle.px() +
                                 evec_basis[2] * particle.py();
    return FourVector(particle.e(), three_momentum);
  }

  /**
   * perform string fragmentation to determine species and momenta of hadrons
   * by exploiting PYTHIA 8.3 \iref{Andersson:1983ia}, \iref{Sjostrand:2014zea},
   * \iref{Bierlich:2022pfr}.
   * \param[in] idq1 PDG id of quark or anti-diquark (carrying color index).
   * \param[in] idq2 PDG id of diquark or anti-quark (carrying anti-color index).
   * \param[in] mString the string mass. [GeV]
   * \param[out] evecLong unit 3-vector specifying the direction of diquark or
   *                      anti-diquark.
   * \param[in] flip_string_ends is whether or not we randomly switch string ends.
   * \param[in] separate_fragment_baryon is whether fragment leading baryon
   *            (or anti-baryon) with separate fragmentation function.
   * \param[out] intermediate_particles list of fragmented particles
   * \return number of hadrons fragmented out of string.
   *
   * \throw std::runtime_error
   *        if string mass is lower than threshold set by PYTHIA
   */
  int fragment_string(int idq1, int idq2, double mString,
                      ThreeVector &evecLong, bool flip_string_ends,
                      bool separate_fragment_baryon,
                      ParticleList &intermediate_particles);

  /**
   * Fragment one hadron from a given string configuration if the string mass
   * is above threshold (given by the consituent masses).
   * Otherwise the entire string breaks down into (final) two hadrons.
   * \param[in] from_forward whether a hadron is fragmented from the forward end
   *            of the string
   * \param[in] separate_fragment_baryon whether a separate fragmentation function
   *            is used for the leading baryon from the diquark end of string.
   * \param[in] evec_basis three orthonormal basis vectors of which
   *            evec_basis[0] is in the longitudinal direction while
   *            evec_basis[1] and evec_basis[2] span the transverse plane.
   * \param[out] ppos_string lightcone momentum p^+ of the string.
   *             This will be changed according to that of the remaining string.
   * \param[out] pneg_string lightcone momentum p^- of the string.
   *             This will be changed according to that of the remaining string.
   * \param[out] QTrx_string_pos transverse momentum px carried by
   *             the forward end of the string.
   *             This will be changed according to that of the remaining string.
   * \param[out] QTrx_string_neg transverse momentum px carried by
   *             the backward end of the string.
   *             This will be changed according to that of the remaining string.
   * \param[out] QTry_string_pos transverse momentum py carried by
   *             the forward end of the string.
   *             This will be changed according to that of the remaining string.
   * \param[out] QTry_string_neg transverse momentum py carried by
   *             the backward end of the string.
   *             This will be changed according to that of the remaining string.
   * \param[out] flav_string_pos constituent flavor at the forward end of the string.
   *             This will be changed according to that of the remaining string.
   * \param[out] flav_string_neg constituent flavor at the backward end of the string.
   *             This will be changed according to that of the remaining string.
   * \param[out] pdgid_frag PDG id of fragmented hadron(s)
   * \param[out] momentum_frag four-momenta of fragmented hadrons(s)
   * \return number of fragmented hadron(s). It can be 1 or 2 depending on
   *         the string mass. If it fails to fragment, returns 0.
   */
  int fragment_off_hadron(bool from_forward,
                          bool separate_fragment_baryon,
                          std::array<ThreeVector, 3> &evec_basis,
                          double &ppos_string, double &pneg_string,
                          double &QTrx_string_pos, double &QTrx_string_neg,
                          double &QTry_string_pos, double &QTry_string_neg,
                          Pythia8::FlavContainer &flav_string_pos,
                          Pythia8::FlavContainer &flav_string_neg,
                          std::vector<int> &pdgid_frag,
                          std::vector<FourVector> &momentum_frag);

  /**
   * Determines hadron type from valence quark constituents.
   * First, try with PYTHIA routine.
   * If it does not work, select a resonance with the same quantum numbers.
   * The probability to pick each resonance in this case is proportional to
   * spin degeneracy / mass, which is inspired by UrQMD.
   * \param[in] idq1 PDG id of a valence quark constituent.
   * \param[in] idq2 PDG id of another valence quark constituent.
   * \return PDG id of selected hadronic species.
   */
  int get_hadrontype_from_quark(int idq1, int idq2);

  /**
   * \param[in] idq1 id of the valence quark (anti-diquark)
   * \param[in] idq2 id of the valence anti-quark (diquark)
   * \param[in] mass mass of the resonance.
   * \return PDG id of resonance with the closest mass.
   *         If the mass is below the threshold or the input PDG id is invalid,
   *         0 is returned.
   */
  int get_resonance_from_quark(int idq1, int idq2, double mass);

  /**
   * Determines lightcone momenta of two final hadrons fragmented from a string
   * in the same way as StringFragmentation::finalTwo in StringFragmentation.cc
   * of PYTHIA 8.
   * \param[in] separate_fragment_hadron whether a separate fragmentation function
   *            is used for the forward hadron
   * \param[in] ppos_string lightcone momentum p^+ of the string
   * \param[in] pneg_string ligntcone momentum p^- of the string
   * \param[in] mTrn_had_forward transverse mass of the forward hadron
   * \param[in] mTrn_had_backward transverse mass of the backward hadron
   * \param[out] ppos_had_forward lightcone momentum p^+ of the forward hadron
   * \param[out] ppos_had_backward lightcone momentum p^+ of the backward hadron
   * \param[out] pneg_had_forward lightcone momentum p^- of the forward hadron
   * \param[out] pneg_had_backward lightcone momentum p^- of the backward hadron
   * \return if the lightcone momenta are properly determined such that
   *         energy and momentum are conserved.
   */
  bool make_lightcone_final_two(
      bool separate_fragment_hadron,
      double ppos_string, double pneg_string,
      double mTrn_had_forward, double mTrn_had_backward,
      double &ppos_had_forward, double &ppos_had_backward,
      double &pneg_had_forward, double &pneg_had_backward);

  /**
   * Sample lightcone momentum fraction according to
   * the LUND fragmentation function.
   * \f$ f(z) = \frac{1}{z} (1 - z)^a \exp{ \left(- \frac{b m_T^2}{z} \right) } \f$
   * \param[in] a parameter for the fragmentation function
   * \param[in] b parameter for the fragmentation function
   * \param[in] mTrn transverse mass of the fragmented hadron
   * \return sampled lightcone momentum fraction
   */
  static double sample_zLund(double a, double b, double mTrn);

  /**
   * \param[out] event_fragments event record which contains information of particles
   * \param[in] evec_basis three orthonormal basis vectors of which
   *            evec_basis[0] is in the longitudinal direction while
   *            evec_basis[1] and evec_basis[2] span the transverse plane.
   * \param[in] ppos_string lightcone momentum p^+ of the string
   * \param[in] pneg_string ligntcone momentum p^- of the string
   * \param[in] QTrx_string transverse momentum px of the string
   * \param[in] QTry_string transverse momentum py of the string
   * \param[in] QTrx_add_pos transverse momentum px to be added
   *            to the most forward hadron
   * \param[in] QTry_add_pos transverse momentum py to be added
   *            to the most forward hadron
   * \param[in] QTrx_add_neg transverse momentum px to be added
   *            to the most backward hadron
   * \param[in] QTry_add_neg transverse momentum py to be added
   *            to the most backward hadron
   */
  bool remake_kinematics_fragments(
      Pythia8::Event &event_fragments, std::array<ThreeVector, 3> &evec_basis,
      double ppos_string, double pneg_string,
      double QTrx_string, double QTry_string,
      double QTrx_add_pos, double QTry_add_pos,
      double QTrx_add_neg, double QTry_add_neg);

  /**
   * Shift the momentum rapidity of all particles in a given event record.
   * y to factor_yrapid * y + diff_yrapid
   * \param[out] event event record which contains information of particles
   * \param[in] evec_basis three orthonormal basis vectors of which
   *            evec_basis[0] is in the longitudinal direction while
   *            evec_basis[1] and evec_basis[2] span the transverse plane.
   * \param[in] factor_yrapid factor multiplied to the old rapidity
   * \param[in] diff_yrapid rapidity difference added to the old one
   */
  void shift_rapidity_event(
      Pythia8::Event &event, std::array<ThreeVector, 3> &evec_basis,
      double factor_yrapid, double diff_yrapid) {
    Pythia8::Vec4 pvec_string_now = Pythia8::Vec4(0., 0., 0., 0.);
    // loop over all particles in the record
    for (int ipyth = 1; ipyth < event.size(); ipyth++) {
      if (!event[ipyth].isFinal()) {
        continue;
      }

      FourVector p_frag = FourVector(
          event[ipyth].e(), event[ipyth].px(),
          event[ipyth].py(), event[ipyth].pz());
      const double E_frag = p_frag.x0();
      const double pl_frag = p_frag.threevec() * evec_basis[0];
      const double ppos_frag = (E_frag + pl_frag) * M_SQRT1_2;
      const double pneg_frag = (E_frag - pl_frag) * M_SQRT1_2;
      const double mTrn_frag = std::sqrt(2. * ppos_frag * pneg_frag);
      // evaluate the old momentum rapidity
      const double y_frag = 0.5 * std::log(ppos_frag / pneg_frag);

      // obtain the new momentum rapidity
      const double y_new_frag = factor_yrapid * y_frag + diff_yrapid;
      // compute the new four momentum
      const double E_new_frag = mTrn_frag * std::cosh(y_new_frag);
      const double pl_new_frag = mTrn_frag * std::sinh(y_new_frag);
      ThreeVector mom_new_frag =
          p_frag.threevec() + (pl_new_frag - pl_frag) * evec_basis[0];
      Pythia8::Vec4 pvec_new_frag = set_Vec4(E_new_frag, mom_new_frag);
      event[ipyth].p(pvec_new_frag);
      pvec_string_now += pvec_new_frag;
    }
    event[0].p(pvec_string_now);
    event[0].m(pvec_string_now.mCalc());
  }

  /**
   * Assign a cross section scaling factor to all outgoing particles.
   *
   * The factor is only non-zero, when the outgoing particle carries
   * a valence quark from the excited hadron. The assigned cross section
   * scaling factor is equal to the number of the valence quarks from the
   * fragmented hadron contained in the fragment divided by the total number
   * of valence quarks of that fragment multiplied by a coherence factor
   * \param[in] baryon_string baryon number of the string
   * \param[out] outgoing_particles list of string fragments to which scaling
   *             factors are assigned
   * \param[in] evecLong direction in which the string is stretched
   * \param[in] suppression_factor additional coherence factor to be
   *            multiplied with scaling factor
   */
  static void assign_all_scaling_factors(int baryon_string,
                                         ParticleList& outgoing_particles,
                                         const ThreeVector &evecLong,
                                         double suppression_factor);

  /**
   * Find the leading string fragments
   *
   * Find the first particle, which can carry nq1, and the last particle,
   * which can carry nq2 valence quarks and return their indices in
   * the given list.
   *
   * \param[in] nq1 number of valence quarks from excited hadron at forward
   *                end of the string
   * \param[in] nq2 number of valance quarks from excited hadron at backward
   *                end of the string
   * \param[in] list list of string fragments
   * \return indices of the leading hadrons in \p list
   */
  static std::pair<int, int> find_leading(int nq1, int nq2, ParticleList& list);

  /**
   * Assign a cross section scaling factor to the given particle.
   *
   * The scaling factor is the number of quarks from the excited hadron,
   * that the fragment carries devided by the total number of quarks in
   * this fragment multiplied by coherence factor.
   *
   * \param[in] nquark number of valence quarks from the excited hadron
   *            contained in the given string fragment \p data
   * \param[out] data particle to assign a scaling factor to
   * \param[in] suppression_factor coherence factor to decrease scaling factor
   */
  static void assign_scaling_factor(int nquark, ParticleData& data,
                                    double suppression_factor);

  /**
   * Take pdg code and map onto particle specie
   * which can be handled by PYTHIA.
   * Positively charged baryons are mapped onto proton and other baryons are
   * mapped onto neutrons. Same rule applies for anti-baryons.
   * Positively (negatively) charged mesons are mapped onto pi+ (pi-).
   * Negatively and positively charged leptons are mapped respectivly onto
   * electron and positron.
   * Currently, we do not have cross sections for leptons and photons
   * with high energy, so such collisions should not happen.
   *
   * \param[in] pdg PdgCode that will be mapped
   * \return mapped PDG id to be used in PYTHIA
   *
   * \throw std::runtime_error
   *        if the incoming particle is neither hadron nor lepton.
   */
  static int pdg_map_for_pythia(PdgCode &pdg);


  /**
   * \return forward lightcone momentum incoming particle A in CM-frame [GeV]
   * \see PPosA_
   */
  double getPPosA() { return PPosA_;}

  /**
   * \return backward lightcone momentum incoming particle Af in CM-frame [GeV]
   * \see PNegA_
   */
  double getPNegA() { return PNegA_;}

  /**
   * \return forward lightcone momentum incoming particle B in CM-frame [GeV]
   * \see PPosB_
   */
  double getPPosB() { return PPosB_;}

  /**
   * \return backward lightcone momentum incoming particle B in CM-frame [GeV]
   * \see PNegB_
   */
  double getPnegB() { return PNegB_;}

  /**
   * \return mass of incoming particle A [GeV]
   * \see massA_
   */
  double get_massA() { return massA_;}

  /**
   * \return mass of incoming particle B [GeV]
   * \see massB_
   */
  double get_massB() { return massB_;}

  /**
   * \return sqrt of mandelstam s [GeV]
   * \see sqrtsAB_
   */
  double get_sqrts() { return sqrtsAB_;}

  /**
   * \return array with PDG Codes of incoming particles
   * \see PDGcodes_
   */
  std::array<PdgCode, 2> get_PDGs() { return PDGcodes_;}

  /**
   * \return momenta of incoming particles in lab frame [GeV]
   * \see plab_
   */
  std::array<FourVector, 2> get_plab() { return plab_;}

  /**
   * \return momenta of incoming particles in center of mass frame [GeV]
   * \see pcom_
   */
  std::array<FourVector, 2> get_pcom() { return pcom_;}

  /**
   * \return velocity four vector of the COM in the lab frame
   * \see ucomAB_
   */
  FourVector get_ucom() { return ucomAB_;}

  /**
   * \return velocity three vector of the COM in the lab frame
   * \see vcomAB_
   */
  ThreeVector get_vcom() { return vcomAB_;}

  /**
   * \return collision time
   * \see time_collision_
   */
  double get_tcoll() { return time_collision_;}

  // clang-format on
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_STRINGPROCESS_H_
