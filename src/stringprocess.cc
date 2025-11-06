/*
 *
 *    Copyright (c) 2017-2020,2022,2024-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/stringprocess.h"

#include <array>

#include "smash/angles.h"
#include "smash/kinematics.h"
#include "smash/pow.h"
#include "smash/random.h"

namespace smash {
static constexpr int LOutput = LogArea::Output::id;

StringProcess::StringProcess(
    double string_tension, double time_formation, double gluon_beta,
    double gluon_pmin, double quark_alpha, double quark_beta,
    double strange_supp, double diquark_supp, double sigma_perp,
    double stringz_a_leading, double stringz_b_leading, double stringz_a,
    double stringz_b, double string_sigma_T, double factor_t_form,
    bool mass_dependent_formation_times, double prob_proton_to_d_uu,
    bool separate_fragment_baryon, double popcorn_rate, bool use_monash_tune,
    SpinInteractionType spin_interaction_type)
    : pmin_gluon_lightcone_(gluon_pmin),
      pow_fgluon_beta_(gluon_beta),
      pow_fquark_alpha_(quark_alpha),
      pow_fquark_beta_(quark_beta),
      sigma_qperp_(sigma_perp),
      stringz_a_leading_(stringz_a_leading),
      stringz_b_leading_(stringz_b_leading),
      stringz_a_produce_(stringz_a),
      stringz_b_produce_(stringz_b),
      strange_supp_(strange_supp),
      diquark_supp_(diquark_supp),
      popcorn_rate_(popcorn_rate),
      string_sigma_T_(string_sigma_T),
      kappa_tension_string_(string_tension),
      additional_xsec_supp_(0.7),
      time_formation_const_(time_formation),
      soft_t_form_(factor_t_form),
      time_collision_(0.),
      mass_dependent_formation_times_(mass_dependent_formation_times),
      prob_proton_to_d_uu_(prob_proton_to_d_uu),
      separate_fragment_baryon_(separate_fragment_baryon),
      use_monash_tune_(use_monash_tune),
      spin_interaction_type_(spin_interaction_type) {
  // setup and initialize pythia for fragmentation
  pythia_hadron_ = std::make_unique<Pythia8::Pythia>(PYTHIA_XML_DIR, false);
  /* turn off all parton-level processes to implement only hadronization */
  pythia_hadron_->readString("ProcessLevel:all = off");
  common_setup_pythia(pythia_hadron_.get(), strange_supp, diquark_supp,
                      popcorn_rate, stringz_a, stringz_b, string_sigma_T);

  /* initialize PYTHIA */
  pythia_hadron_->init();

  /*
   * The const_cast<type>() function is used to obtain the reference of the
   * PrivateInfo object in the pythia_hadron_.
   * This cast is needed since Pythia 8.302 which included a major architecture
   * change. The Info object of a Pythia object is now private, only a const
   * reference can be obtained.
   * In order to reference the PrivateInfo object during initialization, we
   * cast the const reference to obtain the stored address.
   */
  pythia_sigmatot_.initInfoPtr(
      const_cast<Pythia8::Info &>(pythia_hadron_->info));
  pythia_sigmatot_.init();

  pythia_stringflav_.initInfoPtr(
      const_cast<Pythia8::Info &>(pythia_hadron_->info));
  pythia_stringflav_.init();

  event_intermediate_.init("intermediate partons",
                           &pythia_hadron_->particleData);

  for (int imu = 0; imu < 3; imu++) {
    evecBasisAB_[imu] = ThreeVector(0., 0., 0.);
  }

  final_state_.clear();
}

void StringProcess::common_setup_pythia(Pythia8::Pythia *pythia_in,
                                        double strange_supp,
                                        double diquark_supp,
                                        double popcorn_rate, double stringz_a,
                                        double stringz_b,
                                        double string_sigma_T) {
  // choose parametrization for mass-dependent width
  pythia_in->readString("ParticleData:modeBreitWigner = 4");
  /* choose minimum transverse momentum scale
   * involved in partonic interactions */
  pythia_in->readString("MultipartonInteractions:pTmin = 1.5");
  pythia_in->readString("MultipartonInteractions:nSample = 100000");
  // transverse momentum spread in string fragmentation
  pythia_in->readString("StringPT:sigma = " + std::to_string(string_sigma_T));
  // diquark suppression factor in string fragmentation
  pythia_in->readString("StringFlav:probQQtoQ = " +
                        std::to_string(diquark_supp));
  // strangeness suppression factor in string fragmentation
  pythia_in->readString("StringFlav:probStoUD = " +
                        std::to_string(strange_supp));
  pythia_in->readString("StringFlav:popcornRate = " +
                        std::to_string(popcorn_rate));
  // parameters for the fragmentation function
  pythia_in->readString("StringZ:aLund = " + std::to_string(stringz_a));
  pythia_in->readString("StringZ:bLund = " + std::to_string(stringz_b));

  // manually set the parton distribution function
  pythia_in->readString("PDF:pSet = 13");
  pythia_in->readString("PDF:pSetB = 13");
  pythia_in->readString("PDF:piSet = 1");
  pythia_in->readString("PDF:piSetB = 1");
  pythia_in->readString("Beams:idA = 2212");
  pythia_in->readString("Beams:idB = 2212");
  pythia_in->readString("Beams:eCM = 10.");

  // set PYTHIA random seed from outside
  pythia_in->readString("Random:setSeed = on");
  // suppress unnecessary output
  pythia_in->readString("Print:quiet = on");
  // No resonance decays, since the resonances will be handled by SMASH
  pythia_in->readString("HadronLevel:Decay = off");
  // set particle masses and widths in PYTHIA to be same with those in SMASH
  for (auto &ptype : ParticleType::list_all()) {
    int pdgid = ptype.pdgcode().get_decimal();
    double mass_pole = ptype.mass();
    double width_pole = ptype.width_at_pole();
    // check if the particle species is in PYTHIA
    if (pythia_in->particleData.isParticle(pdgid)) {
      // set mass and width in PYTHIA
      pythia_in->particleData.m0(pdgid, mass_pole);
      pythia_in->particleData.mWidth(pdgid, width_pole);
    } else if (pdgid == 310 || pdgid == 130) {
      // set mass and width of Kaon-L and Kaon-S
      pythia_in->particleData.m0(pdgid, kaon_mass);
      pythia_in->particleData.mWidth(pdgid, 0.);
    }
  }

  // make energy-momentum conservation in PYTHIA more precise
  pythia_in->readString("Check:epTolErr = 1e-6");
  pythia_in->readString("Check:epTolWarn = 1e-8");
  if (use_monash_tune_) {
    pythia_in->readString("Tune:ee = 7");
    pythia_in->readString("Tune:pp = 14");
  }
}

// compute the formation time and fill the arrays with final-state particles
int StringProcess::append_final_state(ParticleList &intermediate_particles,
                                      const FourVector &uString,
                                      const ThreeVector &evecLong) {
  int nfrag = 0;
  int bstring = 0;

  for (ParticleData &data : intermediate_particles) {
    nfrag += 1;
    bstring += data.pdgcode().baryon_number();
  }
  assert(nfrag > 0);

  /* compute the cross section scaling factor for leading hadrons
   * based on the number of valence quarks. */
  assign_all_scaling_factors(bstring, intermediate_particles, evecLong,
                             additional_xsec_supp_);

  // Velocity three-vector to perform Lorentz boost.
  const ThreeVector vstring = uString.velocity();

  // compute the formation times of hadrons
  for (int i = 0; i < nfrag; i++) {
    ThreeVector velocity = intermediate_particles[i].momentum().velocity();
    double gamma = 1. / intermediate_particles[i].inverse_gamma();
    // boost 4-momentum into the center of mass frame
    FourVector momentum =
        intermediate_particles[i].momentum().lorentz_boost(-vstring);
    intermediate_particles[i].set_4momentum(momentum);

    if (mass_dependent_formation_times_) {
      // set the formation time and position in the rest frame of string
      double tau_prod = M_SQRT2 * intermediate_particles[i].effective_mass() /
                        kappa_tension_string_;
      double t_prod = tau_prod * gamma;
      FourVector fragment_position = FourVector(t_prod, t_prod * velocity);
      /* boost formation position into the center of mass frame
       * and then into the lab frame */
      fragment_position = fragment_position.lorentz_boost(-vstring);
      fragment_position = fragment_position.lorentz_boost(-vcomAB_);
      intermediate_particles[i].set_slow_formation_times(
          time_collision_,
          soft_t_form_ * fragment_position.x0() + time_collision_);
    } else {
      ThreeVector v_calc = momentum.lorentz_boost(-vcomAB_).velocity();
      double gamma_factor = 1.0 / std::sqrt(1 - (v_calc).sqr());
      intermediate_particles[i].set_slow_formation_times(
          time_collision_,
          time_formation_const_ * gamma_factor + time_collision_);
    }

    final_state_.push_back(intermediate_particles[i]);
  }

  return nfrag;
}

void StringProcess::init(const ParticleList &incoming, double tcoll) {
  PDGcodes_[0] = incoming[0].pdgcode();
  PDGcodes_[1] = incoming[1].pdgcode();
  massA_ = incoming[0].effective_mass();
  massB_ = incoming[1].effective_mass();

  plab_[0] = incoming[0].momentum();
  plab_[1] = incoming[1].momentum();

  // compute sqrts and velocity of the center of mass.
  sqrtsAB_ = (plab_[0] + plab_[1]).abs();
  ucomAB_ = (plab_[0] + plab_[1]) / sqrtsAB_;
  vcomAB_ = ucomAB_.velocity();

  pcom_[0] = plab_[0].lorentz_boost(vcomAB_);
  pcom_[1] = plab_[1].lorentz_boost(vcomAB_);

  const double pabscomAB = pCM(sqrtsAB_, massA_, massB_);
  ThreeVector evec_coll = pcom_[0].threevec() / pabscomAB;
  make_orthonormal_basis(evec_coll, evecBasisAB_);

  compute_incoming_lightcone_momenta();

  time_collision_ = tcoll;
}

/* single diffractive
 * is_AB_to_AX = true  : A + B -> A + X
 * is_AB_to_AX = false : A + B -> X + B */
bool StringProcess::next_SDiff(bool is_AB_to_AX) {
  NpartFinal_ = 0;
  NpartString_[0] = 0;
  NpartString_[1] = 0;
  final_state_.clear();

  double massH = is_AB_to_AX ? massA_ : massB_;
  double mstrMin = is_AB_to_AX ? massB_ : massA_;
  double mstrMax = sqrtsAB_ - massH;

  int idqX1, idqX2;
  double QTrn, QTrx, QTry;
  double pabscomHX_sqr, massX;

  // decompose hadron into quarks
  make_string_ends(is_AB_to_AX ? PDGcodes_[1] : PDGcodes_[0], idqX1, idqX2,
                   prob_proton_to_d_uu_);
  // string mass must be larger than threshold set by PYTHIA.
  mstrMin = pythia_hadron_->particleData.m0(idqX1) +
            pythia_hadron_->particleData.m0(idqX2);
  // this threshold cannot be larger than maximum of allowed string mass.
  if (mstrMin > mstrMax) {
    return false;
  }
  // sample the transverse momentum transfer
  QTrx = random::normal(0., sigma_qperp_ * M_SQRT1_2);
  QTry = random::normal(0., sigma_qperp_ * M_SQRT1_2);
  QTrn = std::sqrt(QTrx * QTrx + QTry * QTry);
  /* sample the string mass and
   * evaluate the three-momenta of hadron and string. */
  massX = random::power(-1.0, mstrMin, mstrMax);
  pabscomHX_sqr = pCM_sqr(sqrtsAB_, massH, massX);
  /* magnitude of the three momentum must be larger
   * than the transverse momentum. */
  const bool foundPabsX = pabscomHX_sqr > QTrn * QTrn;

  if (!foundPabsX) {
    return false;
  }
  double sign_direction = is_AB_to_AX ? 1. : -1.;
  // determine three momentum of the final state hadron
  const ThreeVector cm_momentum =
      sign_direction *
      (evecBasisAB_[0] * std::sqrt(pabscomHX_sqr - QTrn * QTrn) +
       evecBasisAB_[1] * QTrx + evecBasisAB_[2] * QTry);
  const FourVector pstrHcom(std::sqrt(pabscomHX_sqr + massH * massH),
                            cm_momentum);
  const FourVector pstrXcom(std::sqrt(pabscomHX_sqr + massX * massX),
                            -cm_momentum);

  const FourVector ustrXcom = pstrXcom / massX;
  /* determine direction in which the string is stretched.
   * this is set to be same with the the collision axis
   * in the center of mass frame. */
  const ThreeVector threeMomentum =
      is_AB_to_AX ? pcom_[1].threevec() : pcom_[0].threevec();
  const FourVector pnull = FourVector(threeMomentum.abs(), threeMomentum);
  const FourVector prs = pnull.lorentz_boost(ustrXcom.velocity());
  ThreeVector evec = prs.threevec() / prs.threevec().abs();
  // perform fragmentation and add particles to final_state.
  ParticleList new_intermediate_particles;
  int nfrag = fragment_string(idqX1, idqX2, massX, evec, true, false,
                              new_intermediate_particles);
  if (nfrag < 1) {
    NpartString_[0] = 0;
    return false;
  }

  // Set an unpolarized spin vector for the new intermediate particles
  if (spin_interaction_type_ != SpinInteractionType::Off) {
    for (ParticleData &new_particle : new_intermediate_particles) {
      new_particle.set_unpolarized_spin_vector();
    }
  }

  NpartString_[0] =
      append_final_state(new_intermediate_particles, ustrXcom, evec);

  NpartString_[1] = 1;
  PdgCode hadron_code = is_AB_to_AX ? PDGcodes_[0] : PDGcodes_[1];
  ParticleData new_particle(ParticleType::find(hadron_code));
  new_particle.set_4momentum(pstrHcom);
  new_particle.set_cross_section_scaling_factor(1.);
  new_particle.set_formation_time(time_collision_);
  if (spin_interaction_type_ != SpinInteractionType::Off) {
    new_particle.set_unpolarized_spin_vector();
  }
  final_state_.push_back(new_particle);

  NpartFinal_ = NpartString_[0] + NpartString_[1];
  return true;
}

bool StringProcess::set_mass_and_direction_2strings(
    const std::array<std::array<int, 2>, 2> &quarks,
    const std::array<FourVector, 2> &pstr_com, std::array<double, 2> &m_str,
    std::array<ThreeVector, 2> &evec_str) {
  std::array<bool, 2> found_mass;
  for (int i = 0; i < 2; i++) {
    found_mass[i] = false;

    m_str[i] = pstr_com[i].sqr();
    m_str[i] = (m_str[i] > 0.) ? std::sqrt(m_str[i]) : 0.;
    const double threshold = pythia_hadron_->particleData.m0(quarks[i][0]) +
                             pythia_hadron_->particleData.m0(quarks[i][1]);
    // string mass must be larger than threshold set by PYTHIA.
    if (m_str[i] > threshold) {
      found_mass[i] = true;
      /* Determine direction in which string i is stretched.
       * This is set to be same with the collision axis
       * in the center of mass frame.
       * Initial state partons inside incoming hadrons are
       * moving along the collision axis,
       * which is parallel to three momenta of incoming hadrons
       * in the center of mass frame.
       * Given that partons are assumed to be massless,
       * their four momenta are null vectors and parallel to pnull.
       * If we take unit three-vector of prs,
       * which is pnull in the rest frame of string,
       * it would be the direction in which string ends are moving. */
      const ThreeVector mom = pcom_[i].threevec();
      const FourVector pnull(mom.abs(), mom);
      const FourVector prs = pnull.lorentz_boost(pstr_com[i].velocity());
      evec_str[i] = prs.threevec() / prs.threevec().abs();
    }
  }

  return found_mass[0] && found_mass[1];
}

bool StringProcess::make_final_state_2strings(
    const std::array<std::array<int, 2>, 2> &quarks,
    const std::array<FourVector, 2> &pstr_com,
    const std::array<double, 2> &m_str,
    const std::array<ThreeVector, 2> &evec_str, bool flip_string_ends,
    bool separate_fragment_baryon) {
  const std::array<FourVector, 2> ustr_com = {pstr_com[0] / m_str[0],
                                              pstr_com[1] / m_str[1]};
  for (int i = 0; i < 2; i++) {
    ParticleList new_intermediate_particles;

    // determine direction in which string i is stretched.
    ThreeVector evec = evec_str[i];
    // perform fragmentation and add particles to final_state.
    int nfrag = fragment_string(quarks[i][0], quarks[i][1], m_str[i], evec,
                                flip_string_ends, separate_fragment_baryon,
                                new_intermediate_particles);
    if (nfrag <= 0) {
      NpartString_[i] = 0;
      return false;
    }

    // Set an unpolarized spin vector for the new intermediate particles
    if (spin_interaction_type_ != SpinInteractionType::Off) {
      for (ParticleData &new_particle : new_intermediate_particles) {
        new_particle.set_unpolarized_spin_vector();
      }
    }

    NpartString_[i] =
        append_final_state(new_intermediate_particles, ustr_com[i], evec);
    assert(nfrag == NpartString_[i]);
  }
  if ((NpartString_[0] > 0) && (NpartString_[1] > 0)) {
    NpartFinal_ = NpartString_[0] + NpartString_[1];
    return true;
  }
  return false;
}

// double-diffractive : A + B -> X + X
bool StringProcess::next_DDiff() {
  NpartFinal_ = 0;
  NpartString_[0] = 0;
  NpartString_[1] = 0;
  final_state_.clear();

  std::array<std::array<int, 2>, 2> quarks;
  std::array<FourVector, 2> pstr_com;
  std::array<double, 2> m_str;
  std::array<ThreeVector, 2> evec_str;
  ThreeVector threeMomentum;

  // decompose hadron into quark (and diquark) contents
  make_string_ends(PDGcodes_[0], quarks[0][0], quarks[0][1],
                   prob_proton_to_d_uu_);
  make_string_ends(PDGcodes_[1], quarks[1][0], quarks[1][1],
                   prob_proton_to_d_uu_);
  // sample the lightcone momentum fraction carried by gluons
  const double xmin_gluon_fraction = pmin_gluon_lightcone_ / sqrtsAB_;
  const double xfracA =
      random::beta_a0(xmin_gluon_fraction, pow_fgluon_beta_ + 1.);
  const double xfracB =
      random::beta_a0(xmin_gluon_fraction, pow_fgluon_beta_ + 1.);
  // sample the transverse momentum transfer
  const double QTrx = random::normal(0., sigma_qperp_ * M_SQRT1_2);
  const double QTry = random::normal(0., sigma_qperp_ * M_SQRT1_2);
  const double QTrn = std::sqrt(QTrx * QTrx + QTry * QTry);
  // evaluate the lightcone momentum transfer
  const double QPos = -QTrn * QTrn / (2. * xfracB * PNegB_);
  const double QNeg = QTrn * QTrn / (2. * xfracA * PPosA_);
  // compute four-momentum of string 1
  threeMomentum =
      evecBasisAB_[0] * (PPosA_ + QPos - PNegA_ - QNeg) * M_SQRT1_2 +
      evecBasisAB_[1] * QTrx + evecBasisAB_[2] * QTry;
  pstr_com[0] =
      FourVector((PPosA_ + QPos + PNegA_ + QNeg) * M_SQRT1_2, threeMomentum);
  // compute four-momentum of string 2
  threeMomentum =
      evecBasisAB_[0] * (PPosB_ - QPos - PNegB_ + QNeg) * M_SQRT1_2 -
      evecBasisAB_[1] * QTrx - evecBasisAB_[2] * QTry;
  pstr_com[1] =
      FourVector((PPosB_ - QPos + PNegB_ - QNeg) * M_SQRT1_2, threeMomentum);

  const bool found_masses =
      set_mass_and_direction_2strings(quarks, pstr_com, m_str, evec_str);
  if (!found_masses) {
    return false;
  }
  const bool flip_string_ends = true;
  const bool success = make_final_state_2strings(
      quarks, pstr_com, m_str, evec_str, flip_string_ends, false);
  return success;
}

// soft non-diffractive
bool StringProcess::next_NDiffSoft() {
  NpartFinal_ = 0;
  NpartString_[0] = 0;
  NpartString_[1] = 0;
  final_state_.clear();

  std::array<std::array<int, 2>, 2> quarks;
  std::array<FourVector, 2> pstr_com;
  std::array<double, 2> m_str;
  std::array<ThreeVector, 2> evec_str;

  // decompose hadron into quark (and diquark) contents
  int idqA1, idqA2, idqB1, idqB2;
  make_string_ends(PDGcodes_[0], idqA1, idqA2, prob_proton_to_d_uu_);
  make_string_ends(PDGcodes_[1], idqB1, idqB2, prob_proton_to_d_uu_);

  const int bar_a = PDGcodes_[0].baryon_number(),
            bar_b = PDGcodes_[1].baryon_number();
  if (bar_a == 1 ||  // baryon-baryon, baryon-meson, baryon-antibaryon
      (bar_a == 0 && bar_b == 1) ||  // meson-baryon
      (bar_a == 0 && bar_b == 0)) {  // meson-meson
    quarks[0][0] = idqB1;
    quarks[0][1] = idqA2;
    quarks[1][0] = idqA1;
    quarks[1][1] = idqB2;
  } else if (((bar_a == 0) && (bar_b == -1)) ||  // meson-antibaryon
             (bar_a == -1)) {  // antibaryon-baryon, antibaryon-meson,
                               // antibaryon-antibaryon
    quarks[0][0] = idqA1;
    quarks[0][1] = idqB2;
    quarks[1][0] = idqB1;
    quarks[1][1] = idqA2;
  } else {
    std::stringstream ss;
    ss << "  StringProcess::next_NDiff : baryonA = " << bar_a
       << ", baryonB = " << bar_b;
    throw std::runtime_error(ss.str());
  }
  // sample the lightcone momentum fraction carried by quarks
  const double xfracA = random::beta(pow_fquark_alpha_, pow_fquark_beta_);
  const double xfracB = random::beta(pow_fquark_alpha_, pow_fquark_beta_);
  // sample the transverse momentum transfer
  const double QTrx = random::normal(0., sigma_qperp_ * M_SQRT1_2);
  const double QTry = random::normal(0., sigma_qperp_ * M_SQRT1_2);
  const double QTrn = std::sqrt(QTrx * QTrx + QTry * QTry);
  // evaluate the lightcone momentum transfer
  const double QPos = -QTrn * QTrn / (2. * xfracB * PNegB_);
  const double QNeg = QTrn * QTrn / (2. * xfracA * PPosA_);
  const double dPPos = -xfracA * PPosA_ - QPos;
  const double dPNeg = xfracB * PNegB_ - QNeg;
  // compute four-momentum of string 1
  ThreeVector threeMomentum =
      evecBasisAB_[0] * (PPosA_ + dPPos - PNegA_ - dPNeg) * M_SQRT1_2 +
      evecBasisAB_[1] * QTrx + evecBasisAB_[2] * QTry;
  pstr_com[0] =
      FourVector((PPosA_ + dPPos + PNegA_ + dPNeg) * M_SQRT1_2, threeMomentum);
  m_str[0] = pstr_com[0].sqr();
  // compute four-momentum of string 2
  threeMomentum =
      evecBasisAB_[0] * (PPosB_ - dPPos - PNegB_ + dPNeg) * M_SQRT1_2 -
      evecBasisAB_[1] * QTrx - evecBasisAB_[2] * QTry;
  pstr_com[1] =
      FourVector((PPosB_ - dPPos + PNegB_ - dPNeg) * M_SQRT1_2, threeMomentum);

  const bool found_masses =
      set_mass_and_direction_2strings(quarks, pstr_com, m_str, evec_str);
  if (!found_masses) {
    return false;
  }
  const bool flip_string_ends = false;
  const bool success =
      make_final_state_2strings(quarks, pstr_com, m_str, evec_str,
                                flip_string_ends, separate_fragment_baryon_);
  return success;
}

// hard non-diffractive
bool StringProcess::next_NDiffHard() {
  NpartFinal_ = 0;
  final_state_.clear();

  logg[LPythia].debug("Hard non-diff. with ", PDGcodes_[0], " + ", PDGcodes_[1],
                      " at CM energy [GeV] ", sqrtsAB_);

  std::array<int, 2> pdg_for_pythia;
  std::array<std::array<int, 5>, 2> excess_quark;
  std::array<std::array<int, 5>, 2> excess_antiq;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 5; j++) {
      excess_quark[i][j] = 0;
      excess_antiq[i][j] = 0;
    }

    // get PDG id used in PYTHIA event generation
    pdg_for_pythia[i] = pdg_map_for_pythia(PDGcodes_[i]);
    logg[LPythia].debug("  incoming particle ", i, " : ", PDGcodes_[i],
                        " is mapped onto ", pdg_for_pythia[i]);

    PdgCode pdgcode_for_pythia(std::to_string(pdg_for_pythia[i]));
    /* evaluate how many more constituents incoming hadron has
     * compared to the mapped one. */
    find_excess_constituent(PDGcodes_[i], pdgcode_for_pythia, excess_quark[i],
                            excess_antiq[i]);
    logg[LPythia].debug("    excess_quark[", i, "] = (", excess_quark[i][0],
                        ", ", excess_quark[i][1], ", ", excess_quark[i][2],
                        ", ", excess_quark[i][3], ", ", excess_quark[i][4],
                        ")");
    logg[LPythia].debug("    excess_antiq[", i, "] = (", excess_antiq[i][0],
                        ", ", excess_antiq[i][1], ", ", excess_antiq[i][2],
                        ", ", excess_antiq[i][3], ", ", excess_antiq[i][4],
                        ")");
  }

  std::pair<int, int> idAB{pdg_for_pythia[0], pdg_for_pythia[1]};

  // If an entry for the calculated particle IDs does not exist, create one and
  // initialize it accordingly
  if (hard_map_.count(idAB) == 0) {
    hard_map_[idAB] = std::make_unique<Pythia8::Pythia>(PYTHIA_XML_DIR, false);
    hard_map_[idAB]->readString("SoftQCD:nonDiffractive = on");
    hard_map_[idAB]->readString("MultipartonInteractions:pTmin = 1.5");
    hard_map_[idAB]->readString("HadronLevel:all = off");

    common_setup_pythia(hard_map_[idAB].get(), strange_supp_, diquark_supp_,
                        popcorn_rate_, stringz_a_produce_, stringz_b_produce_,
                        string_sigma_T_);

    hard_map_[idAB]->settings.flag("Beams:allowVariableEnergy", true);

    hard_map_[idAB]->settings.mode("Beams:idA", idAB.first);
    hard_map_[idAB]->settings.mode("Beams:idB", idAB.second);
    hard_map_[idAB]->settings.parm("Beams:eCM", sqrtsAB_);

    logg[LPythia].debug("Pythia object initialized with ", pdg_for_pythia[0],
                        " + ", pdg_for_pythia[1], " at CM energy [GeV] ",
                        sqrtsAB_);

    if (!hard_map_[idAB]->init()) {
      throw std::runtime_error("Pythia failed to initialize.");
    }
  }

  const int seed_new = random::uniform_int(1, maximum_rndm_seed_in_pythia);
  hard_map_[idAB]->rndm.init(seed_new);
  logg[LPythia].debug("hard_map_[", idAB.first, "][", idAB.second,
                      "] : rndm is initialized with seed ", seed_new);

  // Change the energy using the Pythia 8.302+ feature

  // Short notation for Pythia event
  Pythia8::Event &event_hadron = pythia_hadron_->event;
  logg[LPythia].debug("Pythia hard event created");
  // we update the collision energy in the CM frame
  hard_map_[idAB]->setKinematics(sqrtsAB_);
  bool final_state_success = hard_map_[idAB]->next();
  logg[LPythia].debug("Pythia final state computed, success = ",
                      final_state_success);
  if (!final_state_success) {
    return false;
  }

  ParticleList new_intermediate_particles;
  ParticleList new_non_hadron_particles;

  Pythia8::Vec4 pSum = 0.;
  event_intermediate_.reset();
  /* Update the partonic intermediate state from PYTHIA output.
   * Note that hadronization will be performed separately,
   * after identification of strings and replacement of constituents. */
  for (int i = 0; i < hard_map_[idAB]->event.size(); i++) {
    if (hard_map_[idAB]->event[i].isFinal()) {
      const int pdgid = hard_map_[idAB]->event[i].id();
      Pythia8::Vec4 pquark = hard_map_[idAB]->event[i].p();
      const double mass = hard_map_[idAB]->particleData.m0(pdgid);

      const int status = hard_map_[idAB]->event[i].status();
      const int color = hard_map_[idAB]->event[i].col();
      const int anticolor = hard_map_[idAB]->event[i].acol();

      pSum += pquark;
      event_intermediate_.append(pdgid, status, color, anticolor, pquark, mass);
    }
  }
  // add junctions to the intermediate state if there is any.
  event_intermediate_.clearJunctions();
  for (int i = 0; i < hard_map_[idAB]->event.sizeJunction(); i++) {
    const int kind = hard_map_[idAB]->event.kindJunction(i);
    std::array<int, 3> col;
    for (int j = 0; j < 3; j++) {
      col[j] = hard_map_[idAB]->event.colJunction(i, j);
    }
    event_intermediate_.appendJunction(kind, col[0], col[1], col[2]);
  }
  /* The zeroth entry of event record is supposed to have the information
   * on the whole system. Specify the total momentum and invariant mass. */
  event_intermediate_[0].p(pSum);
  event_intermediate_[0].m(pSum.mCalc());

  /* Replace quark constituents according to the excess of valence quarks
   * and then rescale momenta of partons by constant factor
   * to fulfill the energy-momentum conservation. */
  bool correct_constituents =
      restore_constituent(event_intermediate_, excess_quark, excess_antiq);
  if (!correct_constituents) {
    logg[LPythia].debug("failed to find correct partonic constituents.");
    return false;
  }

  int npart = event_intermediate_.size();
  int ipart = 0;
  while (ipart < npart) {
    const int pdgid = event_intermediate_[ipart].id();
    if (event_intermediate_[ipart].isFinal() &&
        !event_intermediate_[ipart].isParton() &&
        !hard_map_[idAB]->particleData.isOctetHadron(pdgid)) {
      logg[LPythia].debug("PDG ID from Pythia: ", pdgid);
      FourVector momentum = reorient(event_intermediate_[ipart], evecBasisAB_);
      logg[LPythia].debug("4-momentum from Pythia: ", momentum);
      bool found_ptype =
          append_intermediate_list(pdgid, momentum, new_non_hadron_particles);
      if (!found_ptype) {
        logg[LPythia].warn("PDG ID ", pdgid,
                           " does not exist in ParticleType - start over.");
        final_state_success = false;
      }
      event_intermediate_.remove(ipart, ipart);
      npart -= 1;
    } else {
      ipart += 1;
    }
  }

  bool hadronize_success = false;
  bool find_forward_string = true;
  logg[LPythia].debug("Hard non-diff: partonic process gives ",
                      event_intermediate_.size(), " partons.");
  // identify and fragment strings until there is no parton left.
  while (event_intermediate_.size() > 1) {
    // dummy event to initialize the internal variables of PYTHIA.
    pythia_hadron_->event.reset();
    if (!pythia_hadron_->next()) {
      logg[LPythia].debug("  Dummy event in hard string routine failed.");
      hadronize_success = false;
      break;
    }

    if (event_intermediate_.sizeJunction() > 0) {
      // identify string from a junction if there is any.
      compose_string_junction(find_forward_string, event_intermediate_,
                              pythia_hadron_->event);
    } else {
      /* identify string from a most forward or backward parton.
       * if there is no junction. */
      compose_string_parton(find_forward_string, event_intermediate_,
                            pythia_hadron_->event);
    }

    // fragment the (identified) string into hadrons.
    hadronize_success = pythia_hadron_->forceHadronLevel(false);
    logg[LPythia].debug("Pythia hadronized, success = ", hadronize_success);

    new_intermediate_particles.clear();
    if (hadronize_success) {
      for (int i = 0; i < event_hadron.size(); i++) {
        if (event_hadron[i].isFinal()) {
          int pythia_id = event_hadron[i].id();
          logg[LPythia].debug("PDG ID from Pythia: ", pythia_id);
          /* K_short and K_long need to be converted to K0
           * since SMASH only knows K0 */
          convert_KaonLS(pythia_id);

          /* evecBasisAB_[0] is a unit 3-vector in the collision axis,
           * while evecBasisAB_[1] and evecBasisAB_[2] spans the transverse
           * plane. Given that PYTHIA assumes z-direction to be the collision
           * axis, pz from PYTHIA should be the momentum compoment in
           * evecBasisAB_[0]. px and py are respectively the momentum components
           * in two transverse directions evecBasisAB_[1] and evecBasisAB_[2].
           */
          FourVector momentum = reorient(event_hadron[i], evecBasisAB_);
          logg[LPythia].debug("4-momentum from Pythia: ", momentum);
          logg[LPythia].debug("appending the particle ", pythia_id,
                              " to the intermediate particle list.");
          bool found_ptype = false;
          if (event_hadron[i].isHadron()) {
            found_ptype = append_intermediate_list(pythia_id, momentum,
                                                   new_intermediate_particles);
          } else {
            found_ptype = append_intermediate_list(pythia_id, momentum,
                                                   new_non_hadron_particles);
          }
          if (!found_ptype) {
            logg[LPythia].warn("PDG ID ", pythia_id,
                               " does not exist in ParticleType - start over.");
            hadronize_success = false;
          }
        }
      }
    }

    /* if hadronization is not successful,
     * reset the event records, return false and then start over. */
    if (!hadronize_success) {
      break;
    }

    // Set an unpolarized spin vector for the new intermediate particles
    if (spin_interaction_type_ != SpinInteractionType::Off) {
      for (ParticleData &new_particle : new_intermediate_particles) {
        new_particle.set_unpolarized_spin_vector();
      }
    }

    FourVector uString = FourVector(1., 0., 0., 0.);
    ThreeVector evec = find_forward_string ? evecBasisAB_[0] : -evecBasisAB_[0];
    int nfrag = append_final_state(new_intermediate_particles, uString, evec);
    NpartFinal_ += nfrag;

    find_forward_string = !find_forward_string;
  }

  if (hadronize_success) {
    // add the final state particles, which are not hadron.
    for (ParticleData data : new_non_hadron_particles) {
      data.set_cross_section_scaling_factor(1.);
      data.set_formation_time(time_collision_);
      if (spin_interaction_type_ != SpinInteractionType::Off) {
        data.set_unpolarized_spin_vector();
      }
      final_state_.push_back(data);
    }
  } else {
    final_state_.clear();
  }

  return hadronize_success;
}

void StringProcess::find_excess_constituent(PdgCode &pdg_actual,
                                            PdgCode &pdg_mapped,
                                            std::array<int, 5> &excess_quark,
                                            std::array<int, 5> &excess_antiq) {
  /* decompose PDG id of the actual hadron and mapped one
   * to get the valence quark constituents */
  std::array<int, 3> qcontent_actual = pdg_actual.quark_content();
  std::array<int, 3> qcontent_mapped = pdg_mapped.quark_content();

  excess_quark = {0, 0, 0, 0, 0};
  excess_antiq = {0, 0, 0, 0, 0};
  for (int i = 0; i < 3; i++) {
    if (qcontent_actual[i] > 0) {  // quark content of the actual hadron
      int j = qcontent_actual[i] - 1;
      excess_quark[j] += 1;
    }

    if (qcontent_mapped[i] > 0) {  // quark content of the mapped hadron
      int j = qcontent_mapped[i] - 1;
      excess_quark[j] -= 1;
    }

    if (qcontent_actual[i] < 0) {  // antiquark content of the actual hadron
      int j = std::abs(qcontent_actual[i]) - 1;
      excess_antiq[j] += 1;
    }

    if (qcontent_mapped[i] < 0) {  // antiquark content of the mapped hadron
      int j = std::abs(qcontent_mapped[i]) - 1;
      excess_antiq[j] -= 1;
    }
  }
}

void StringProcess::replace_constituent(
    Pythia8::Particle &particle, std::array<int, 5> &excess_constituent) {
  // If the particle is neither quark nor diquark, nothing to do.
  if (!particle.isQuark() && !particle.isDiquark()) {
    return;
  }

  // If there is no excess of constituents, nothing to do.
  const std::array<int, 5> excess_null = {0, 0, 0, 0, 0};
  if (excess_constituent == excess_null) {
    return;
  }

  int nq = 0;
  std::array<int, 2> pdgid = {0, 0};
  int spin_deg = 0;
  int pdgid_new = 0;
  if (particle.isQuark()) {
    nq = 1;
    pdgid[0] = particle.id();
  } else if (particle.isDiquark()) {
    nq = 2;
    quarks_from_diquark(particle.id(), pdgid[0], pdgid[1], spin_deg);
  }

  for (int iq = 0; iq < nq; iq++) {
    int jq = std::abs(pdgid[iq]) - 1;
    int k_select = 0;
    std::vector<int> k_found;
    k_found.clear();
    // check if the constituent needs to be converted.
    if (excess_constituent[jq] < 0) {
      for (int k = 0; k < 5; k++) {
        // check which specie it can be converted into.
        if (k != jq && excess_constituent[k] > 0) {
          k_found.push_back(k);
        }
      }
    }

    // make a random selection of specie and update the excess of constituent.
    if (k_found.size() > 0) {
      const int l =
          random::uniform_int(0, static_cast<int>(k_found.size()) - 1);
      k_select = k_found[l];
      /* flavor jq + 1 is converted into k_select + 1
       * and excess_constituent is updated. */
      pdgid[iq] = pdgid[iq] > 0 ? k_select + 1 : -(k_select + 1);
      excess_constituent[jq] += 1;
      excess_constituent[k_select] -= 1;
    }
  }

  // determine PDG id of the converted parton.
  if (particle.isQuark()) {
    pdgid_new = pdgid[0];
  } else if (particle.isDiquark()) {
    if (std::abs(pdgid[0]) < std::abs(pdgid[1])) {
      std::swap(pdgid[0], pdgid[1]);
    }

    pdgid_new = std::abs(pdgid[0]) * 1000 + std::abs(pdgid[1]) * 100;
    if (std::abs(pdgid[0]) == std::abs(pdgid[1])) {
      pdgid_new += 3;
    } else {
      pdgid_new += spin_deg;
    }

    if (particle.id() < 0) {
      pdgid_new *= -1;
    }
  }
  logg[LPythia].debug("  parton id = ", particle.id(), " is converted to ",
                      pdgid_new);

  // update the constituent mass and energy.
  Pythia8::Vec4 pquark = particle.p();
  double mass_new = pythia_hadron_->particleData.m0(pdgid_new);
  double e_new = std::sqrt(mass_new * mass_new + pquark.pAbs() * pquark.pAbs());
  // update the particle object.
  particle.id(pdgid_new);
  particle.e(e_new);
  particle.m(mass_new);
}

void StringProcess::find_total_number_constituent(
    Pythia8::Event &event_intermediate, std::array<int, 5> &nquark_total,
    std::array<int, 5> &nantiq_total) {
  for (int iflav = 0; iflav < 5; iflav++) {
    nquark_total[iflav] = 0;
    nantiq_total[iflav] = 0;
  }

  for (int ip = 1; ip < event_intermediate.size(); ip++) {
    if (!event_intermediate[ip].isFinal()) {
      continue;
    }
    const int pdgid = event_intermediate[ip].id();
    if (pdgid > 0) {
      // quarks
      for (int iflav = 0; iflav < 5; iflav++) {
        nquark_total[iflav] +=
            pythia_hadron_->particleData.nQuarksInCode(pdgid, iflav + 1);
      }
    } else {
      // antiquarks
      for (int iflav = 0; iflav < 5; iflav++) {
        nantiq_total[iflav] += pythia_hadron_->particleData.nQuarksInCode(
            std::abs(pdgid), iflav + 1);
      }
    }
  }
}

bool StringProcess::splitting_gluon_qqbar(
    Pythia8::Event &event_intermediate, std::array<int, 5> &nquark_total,
    std::array<int, 5> &nantiq_total, bool sign_constituent,
    std::array<std::array<int, 5>, 2> &excess_constituent) {
  Pythia8::Vec4 pSum = event_intermediate[0].p();

  /* compute total number of quark and antiquark constituents
   * in the whole system. */
  find_total_number_constituent(event_intermediate, nquark_total, nantiq_total);

  for (int iflav = 0; iflav < 5; iflav++) {
    /* Find how many constituent will be in the system after
     * changing the flavors.
     * Note that nquark_total is number of constituent right after
     * the pythia event (with mapped incoming hadrons), while the excess
     * shows how many constituents we have more or less that nquark_total. */
    int nquark_final =
        excess_constituent[0][iflav] + excess_constituent[1][iflav];
    if (sign_constituent) {
      nquark_final += nquark_total[iflav];
    } else {
      nquark_final += nantiq_total[iflav];
    }
    /* Therefore, nquark_final should not be negative.
     * negative nquark_final means that it will not be possible to
     * find a constituent to change the flavor. */
    bool enough_quark = nquark_final >= 0;
    /* If that is the case, a gluon will be splitted into
     * a quark-antiquark pair with the desired flavor. */
    if (!enough_quark) {
      logg[LPythia].debug("  not enough constituents with flavor ", iflav + 1,
                          " : try to split a gluon to qqbar.");
      for (int ic = 0; ic < std::abs(nquark_final); ic++) {
        /* Since each incoming hadron has its own count of the excess,
         * it is necessary to find which one is problematic. */
        int ih_mod = -1;
        if (excess_constituent[0][iflav] < 0) {
          ih_mod = 0;
        } else {
          ih_mod = 1;
        }

        /* find the most forward or backward gluon
         * depending on which incoming hadron is found to be an issue. */
        int iforward = 1;
        for (int ip = 2; ip < event_intermediate.size(); ip++) {
          if (!event_intermediate[ip].isFinal() ||
              !event_intermediate[ip].isGluon()) {
            continue;
          }

          const double y_gluon_current = event_intermediate[ip].y();
          const double y_gluon_forward = event_intermediate[iforward].y();
          if ((ih_mod == 0 && y_gluon_current > y_gluon_forward) ||
              (ih_mod == 1 && y_gluon_current < y_gluon_forward)) {
            iforward = ip;
          }
        }

        if (!event_intermediate[iforward].isGluon()) {
          logg[LPythia].debug("There is no gluon to split into qqbar.");
          return false;
        }

        // four momentum of the original gluon
        Pythia8::Vec4 pgluon = event_intermediate[iforward].p();

        const int pdgid = iflav + 1;
        const double mass = pythia_hadron_->particleData.m0(pdgid);
        const int status = event_intermediate[iforward].status();
        /* color and anticolor indices.
         * the color index of gluon goes to the quark, while
         * the anticolor index goes to the antiquark */
        const int col = event_intermediate[iforward].col();
        const int acol = event_intermediate[iforward].acol();

        // three momenta of quark and antiquark
        std::array<double, 2> px_quark;
        std::array<double, 2> py_quark;
        std::array<double, 2> pz_quark;
        // energies of quark and antiquark
        std::array<double, 2> e_quark;
        // four momenta of quark and antiquark
        std::array<Pythia8::Vec4, 2> pquark;
        // transverse momentum scale of string fragmentation
        const double sigma_qt_frag = pythia_hadron_->parm("StringPT:sigma");
        // sample relative transverse momentum between quark and antiquark
        const double qx = random::normal(0., sigma_qt_frag * M_SQRT1_2);
        const double qy = random::normal(0., sigma_qt_frag * M_SQRT1_2);
        // setup kinematics
        for (int isign = 0; isign < 2; isign++) {
          /* As the simplest assumption, the three momentum of gluon
           * is equally distributed to quark and antiquark.
           * Then, they have a relative transverse momentum. */
          px_quark[isign] = 0.5 * pgluon.px() + (isign == 0 ? 1. : -1.) * qx;
          py_quark[isign] = 0.5 * pgluon.py() + (isign == 0 ? 1. : -1.) * qy;
          pz_quark[isign] = 0.5 * pgluon.pz();
          e_quark[isign] =
              std::sqrt(mass * mass + px_quark[isign] * px_quark[isign] +
                        py_quark[isign] * py_quark[isign] +
                        pz_quark[isign] * pz_quark[isign]);
          pquark[isign] = Pythia8::Vec4(px_quark[isign], py_quark[isign],
                                        pz_quark[isign], e_quark[isign]);
        }

        /* Total energy is not conserved at this point,
         * but this will be cured later. */
        pSum += pquark[0] + pquark[1] - pgluon;
        // add quark and antiquark to the event record
        event_intermediate.append(pdgid, status, col, 0, pquark[0], mass);
        event_intermediate.append(-pdgid, status, 0, acol, pquark[1], mass);
        // then remove the gluon from the record
        event_intermediate.remove(iforward, iforward);

        logg[LPythia].debug("  gluon at iforward = ", iforward,
                            " is splitted into ", pdgid, ",", -pdgid,
                            " qqbar pair.");
        /* Increase the total number of quarks and antiquarks by 1,
         * as we have extra ones from a gluon. */
        nquark_total[iflav] += 1;
        nantiq_total[iflav] += 1;
      }
    }
  }

  /* The zeroth entry of event record is supposed to have the information
   * on the whole system. Specify the total momentum and invariant mass. */
  event_intermediate[0].p(pSum);
  event_intermediate[0].m(pSum.mCalc());

  return true;
}

void StringProcess::rearrange_excess(
    std::array<int, 5> &nquark_total,
    std::array<std::array<int, 5>, 2> &excess_quark,
    std::array<std::array<int, 5>, 2> &excess_antiq) {
  for (int iflav = 0; iflav < 5; iflav++) {
    /* Find how many constituent will be in the system after
     * changing the flavors.
     * Note that nquark_total is number of constituent right after
     * the pythia event (with mapped incoming hadrons), while the excess
     * shows how many constituents we have more or less that nquark_total. */
    int nquark_final =
        nquark_total[iflav] + excess_quark[0][iflav] + excess_quark[1][iflav];
    /* Therefore, nquark_final should not be negative.
     * negative nquark_final means that it will not be possible to
     * find a constituent to change the flavor. */
    bool enough_quark = nquark_final >= 0;
    // If that is the case, excess of constituents will be modified
    if (!enough_quark) {
      logg[LPythia].debug("  not enough constituents with flavor ", iflav + 1,
                          " : try to modify excess of constituents.");
      for (int ic = 0; ic < std::abs(nquark_final); ic++) {
        /* Since each incoming hadron has its own count of the excess,
         * it is necessary to find which one is problematic. */
        int ih_mod = -1;
        if (excess_quark[0][iflav] < 0) {
          ih_mod = 0;
        } else {
          ih_mod = 1;
        }
        /* Increase the excess of both quark and antiquark
         * with corresponding flavor (iflav + 1) by 1.
         * This is for conservation of the net quark number. */
        excess_quark[ih_mod][iflav] += 1;
        excess_antiq[ih_mod][iflav] += 1;

        /* Since incoming hadrons are mapped onto ones with
         * the same baryon number (or quark number),
         * summation of the excesses over all flavors should be zero.
         * Therefore, we need to find another flavor which has
         * a positive excess and subtract by 1. */
        for (int jflav = 0; jflav < 5; jflav++) {
          // another flavor with positive excess of constituents
          if (jflav != iflav && excess_quark[ih_mod][jflav] > 0) {
            /* Decrease the excess of both quark and antiquark
             * with corresponding flavor (jflav + 1) by 1. */
            excess_quark[ih_mod][jflav] -= 1;
            excess_antiq[ih_mod][jflav] -= 1;
            /* We only need to find one (another) flavor to subtract.
             * No more! */
            break;
          }
        }
      }
    }
  }
}

bool StringProcess::restore_constituent(
    Pythia8::Event &event_intermediate,
    std::array<std::array<int, 5>, 2> &excess_quark,
    std::array<std::array<int, 5>, 2> &excess_antiq) {
  Pythia8::Vec4 pSum = event_intermediate[0].p();
  const double energy_init = pSum.e();
  logg[LPythia].debug("  initial total energy [GeV] : ", energy_init);

  // Total number of quarks and antiquarks, respectively.
  std::array<int, 5> nquark_total;
  std::array<int, 5> nantiq_total;

  /* Split a gluon into qqbar if we do not have enough constituents
   * to be converted in the system. */
  bool split_for_quark = splitting_gluon_qqbar(
      event_intermediate, nquark_total, nantiq_total, true, excess_quark);
  bool split_for_antiq = splitting_gluon_qqbar(
      event_intermediate, nquark_total, nantiq_total, false, excess_antiq);

  /* Modify excess_quark and excess_antiq if we do not have enough constituents
   * to be converted in the system. */
  if (!split_for_quark || !split_for_antiq) {
    rearrange_excess(nquark_total, excess_quark, excess_antiq);
    rearrange_excess(nantiq_total, excess_antiq, excess_quark);
  }

  // Final check if there are enough constituents.
  for (int iflav = 0; iflav < 5; iflav++) {
    if (nquark_total[iflav] + excess_quark[0][iflav] + excess_quark[1][iflav] <
        0) {
      logg[LPythia].debug("Not enough quark constituents of flavor ",
                          iflav + 1);
      return false;
    }

    if (nantiq_total[iflav] + excess_antiq[0][iflav] + excess_antiq[1][iflav] <
        0) {
      logg[LPythia].debug("Not enough antiquark constituents of flavor ",
                          -(iflav + 1));
      return false;
    }
  }

  for (int ih = 0; ih < 2; ih++) {
    logg[LPythia].debug("  initial excess_quark[", ih, "] = (",
                        excess_quark[ih][0], ", ", excess_quark[ih][1], ", ",
                        excess_quark[ih][2], ", ", excess_quark[ih][3], ", ",
                        excess_quark[ih][4], ")");
    logg[LPythia].debug("  initial excess_antiq[", ih, "] = (",
                        excess_antiq[ih][0], ", ", excess_antiq[ih][1], ", ",
                        excess_antiq[ih][2], ", ", excess_antiq[ih][3], ", ",
                        excess_antiq[ih][4], ")");
  }

  bool recovered_quarks = false;
  while (!recovered_quarks) {
    /* Flavor conversion begins with the most forward and backward parton
     * respectively for incoming_particles_[0] and incoming_particles_[1]. */
    std::array<bool, 2> find_forward = {true, false};
    const std::array<int, 5> excess_null = {0, 0, 0, 0, 0};
    std::array<int, 5> excess_total = excess_null;

    for (int ih = 0; ih < 2; ih++) {  // loop over incoming hadrons
      int nfrag = event_intermediate.size();
      for (int np_end = 0; np_end < nfrag - 1; np_end++) {  // constituent loop
        /* select the np_end-th most forward or backward parton and
         * change its specie.
         * np_end = 0 corresponds to the most forward,
         * np_end = 1 corresponds to the second most forward and so on. */
        int iforward =
            get_index_forward(find_forward[ih], np_end, event_intermediate);
        pSum -= event_intermediate[iforward].p();

        if (event_intermediate[iforward].id() > 0) {  // quark and diquark
          replace_constituent(event_intermediate[iforward], excess_quark[ih]);
          logg[LPythia].debug(
              "    excess_quark[", ih, "] = (", excess_quark[ih][0], ", ",
              excess_quark[ih][1], ", ", excess_quark[ih][2], ", ",
              excess_quark[ih][3], ", ", excess_quark[ih][4], ")");
        } else {  // antiquark and anti-diquark
          replace_constituent(event_intermediate[iforward], excess_antiq[ih]);
          logg[LPythia].debug(
              "    excess_antiq[", ih, "] = (", excess_antiq[ih][0], ", ",
              excess_antiq[ih][1], ", ", excess_antiq[ih][2], ", ",
              excess_antiq[ih][3], ", ", excess_antiq[ih][4], ")");
        }

        const int pdgid = event_intermediate[iforward].id();
        Pythia8::Vec4 pquark = event_intermediate[iforward].p();
        const double mass = pythia_hadron_->particleData.m0(pdgid);

        const int status = event_intermediate[iforward].status();
        const int color = event_intermediate[iforward].col();
        const int anticolor = event_intermediate[iforward].acol();

        pSum += pquark;
        event_intermediate.append(pdgid, status, color, anticolor, pquark,
                                  mass);

        event_intermediate.remove(iforward, iforward);
        /* Now the last np_end + 1 entries in event_intermediate
         * are np_end + 1 most forward (or backward) partons. */
      }

      // Compute the excess of net quark numbers.
      for (int j = 0; j < 5; j++) {
        excess_total[j] += (excess_quark[ih][j] - excess_antiq[ih][j]);
      }
    }

    /* If there is no excess of net quark numbers,
     * quark content is considered to be correct. */
    recovered_quarks = excess_total == excess_null;
  }
  logg[LPythia].debug("  valence quark contents of hadons are recovered.");

  logg[LPythia].debug("  current total energy [GeV] : ", pSum.e());
  /* rescale momenta of all partons by a constant factor
   * to conserve the total energy. */
  while (true) {
    if (std::abs(pSum.e() - energy_init) <=
        std::abs(really_small * energy_init)) {
      break;
    }

    double energy_current = pSum.e();
    double slope = 0.;
    for (int i = 1; i < event_intermediate.size(); i++) {
      slope += event_intermediate[i].pAbs2() / event_intermediate[i].e();
    }

    const double rescale_factor = 1. + (energy_init - energy_current) / slope;
    pSum = 0.;
    for (int i = 1; i < event_intermediate.size(); i++) {
      const double px = rescale_factor * event_intermediate[i].px();
      const double py = rescale_factor * event_intermediate[i].py();
      const double pz = rescale_factor * event_intermediate[i].pz();
      const double pabs = rescale_factor * event_intermediate[i].pAbs();
      const double mass = event_intermediate[i].m();

      event_intermediate[i].px(px);
      event_intermediate[i].py(py);
      event_intermediate[i].pz(pz);
      event_intermediate[i].e(std::sqrt(mass * mass + pabs * pabs));
      pSum += event_intermediate[i].p();
    }
    logg[LPythia].debug("  parton momenta are rescaled by factor of ",
                        rescale_factor);
  }

  logg[LPythia].debug("  final total energy [GeV] : ", pSum.e());
  /* The zeroth entry of event record is supposed to have the information
   * on the whole system. Specify the total momentum and invariant mass. */
  event_intermediate[0].p(pSum);
  event_intermediate[0].m(pSum.mCalc());

  return true;
}

void StringProcess::compose_string_parton(bool find_forward_string,
                                          Pythia8::Event &event_intermediate,
                                          Pythia8::Event &event_hadronize) {
  Pythia8::Vec4 pSum = 0.;
  event_hadronize.reset();

  // select the most forward or backward parton.
  int iforward = get_index_forward(find_forward_string, 0, event_intermediate);
  logg[LPythia].debug("Hard non-diff: iforward = ", iforward, "(",
                      event_intermediate[iforward].id(), ")");

  pSum += event_intermediate[iforward].p();
  event_hadronize.append(event_intermediate[iforward]);

  int col_to_find = event_intermediate[iforward].acol();
  int acol_to_find = event_intermediate[iforward].col();
  event_intermediate.remove(iforward, iforward);
  logg[LPythia].debug("Hard non-diff: event_intermediate reduces in size to ",
                      event_intermediate.size());

  // trace color and anti-color indices and find corresponding partons.
  while (col_to_find != 0 || acol_to_find != 0) {
    logg[LPythia].debug("  col_to_find = ", col_to_find,
                        ", acol_to_find = ", acol_to_find);

    int ifound = -1;
    for (int i = 1; i < event_intermediate.size(); i++) {
      const int pdgid = event_intermediate[i].id();
      bool found_col =
          col_to_find != 0 && col_to_find == event_intermediate[i].col();
      bool found_acol =
          acol_to_find != 0 && acol_to_find == event_intermediate[i].acol();
      if (found_col) {
        logg[LPythia].debug("  col_to_find ", col_to_find, " from i ", i, "(",
                            pdgid, ") found");
      }
      if (found_acol) {
        logg[LPythia].debug("  acol_to_find ", acol_to_find, " from i ", i, "(",
                            pdgid, ") found");
      }

      if (found_col && !found_acol) {
        ifound = i;
        col_to_find = event_intermediate[i].acol();
        break;
      } else if (!found_col && found_acol) {
        ifound = i;
        acol_to_find = event_intermediate[i].col();
        break;
      } else if (found_col && found_acol) {
        ifound = i;
        col_to_find = 0;
        acol_to_find = 0;
        break;
      }
    }

    if (ifound < 0) {
      event_intermediate.list();
      event_intermediate.listJunctions();
      event_hadronize.list();
      event_hadronize.listJunctions();
      if (col_to_find != 0) {
        logg[LPythia].error("No parton with col = ", col_to_find);
      }
      if (acol_to_find != 0) {
        logg[LPythia].error("No parton with acol = ", acol_to_find);
      }
      throw std::runtime_error("Hard string could not be identified.");
    } else {
      pSum += event_intermediate[ifound].p();
      // add a parton to the new event record.
      event_hadronize.append(event_intermediate[ifound]);
      // then remove from the original event record.
      event_intermediate.remove(ifound, ifound);
      logg[LPythia].debug(
          "Hard non-diff: event_intermediate reduces in size to ",
          event_intermediate.size());
    }
  }

  /* The zeroth entry of event record is supposed to have the information
   * on the whole system. Specify the total momentum and invariant mass. */
  event_hadronize[0].p(pSum);
  event_hadronize[0].m(pSum.mCalc());
}

void StringProcess::compose_string_junction(bool &find_forward_string,
                                            Pythia8::Event &event_intermediate,
                                            Pythia8::Event &event_hadronize) {
  event_hadronize.reset();

  /* Move the first junction to the event record for hadronization
   * and specify color or anti-color indices to be found.
   * If junction kind is an odd number, it connects three quarks
   * to make a color-neutral baryonic configuration.
   * Otherwise, it connects three antiquarks
   * to make a color-neutral anti-baryonic configuration. */
  const int kind = event_intermediate.kindJunction(0);
  bool sign_color = kind % 2 == 1;
  std::vector<int> col;  // color or anti-color indices of the junction legs
  for (int j = 0; j < 3; j++) {
    col.push_back(event_intermediate.colJunction(0, j));
  }
  event_hadronize.appendJunction(kind, col[0], col[1], col[2]);
  event_intermediate.eraseJunction(0);
  logg[LPythia].debug("junction (", col[0], ", ", col[1], ", ", col[2],
                      ") with kind ", kind, " will be handled.");

  bool found_string = false;
  while (!found_string) {
    // trace color or anti-color indices and find corresponding partons.
    find_junction_leg(sign_color, col, event_intermediate, event_hadronize);
    found_string = true;
    for (unsigned int j = 0; j < col.size(); j++) {
      found_string = found_string && col[j] == 0;
    }
    if (!found_string) {
      /* if there is any leg which is not closed with parton,
       * look over junctions and find connected ones. */
      logg[LPythia].debug("  still has leg(s) unfinished.");
      sign_color = !sign_color;
      std::vector<int> junction_to_move;
      for (int i = 0; i < event_intermediate.sizeJunction(); i++) {
        const int kind_new = event_intermediate.kindJunction(i);
        /* If the original junction is associated with positive baryon number,
         * it looks for anti-junctions whose legs are connected with
         * anti-quarks (anti-colors in general). */
        if (sign_color != (kind_new % 2 == 1)) {
          continue;
        }

        std::array<int, 3> col_new;
        for (int k = 0; k < 3; k++) {
          col_new[k] = event_intermediate.colJunction(i, k);
        }

        int n_legs_connected = 0;
        // loop over remaining legs
        for (unsigned int j = 0; j < col.size(); j++) {
          if (col[j] == 0) {
            continue;
          }
          for (int k = 0; k < 3; k++) {
            if (col[j] == col_new[k]) {
              n_legs_connected += 1;
              col[j] = 0;
              col_new[k] = 0;
            }
          }
        }

        // specify which junction is connected to the original one.
        if (n_legs_connected > 0) {
          for (int k = 0; k < 3; k++) {
            if (col_new[k] != 0) {
              col.push_back(col_new[k]);
            }
          }
          logg[LPythia].debug("  junction ", i, " (",
                              event_intermediate.colJunction(i, 0), ", ",
                              event_intermediate.colJunction(i, 1), ", ",
                              event_intermediate.colJunction(i, 2),
                              ") with kind ", kind_new, " will be added.");
          junction_to_move.push_back(i);
        }
      }

      /* If there is any connected junction,
       * move it to the event record for hadronization. */
      for (unsigned int i = 0; i < junction_to_move.size(); i++) {
        unsigned int imove = junction_to_move[i] - i;
        const int kind_add = event_intermediate.kindJunction(imove);
        std::array<int, 3> col_add;
        for (int k = 0; k < 3; k++) {
          col_add[k] = event_intermediate.colJunction(imove, k);
        }
        // add a junction to the new event record.
        event_hadronize.appendJunction(kind_add, col_add[0], col_add[1],
                                       col_add[2]);
        // then remove from the original event record.
        event_intermediate.eraseJunction(imove);
      }
    }
  }

  Pythia8::Vec4 pSum = event_hadronize[0].p();
  find_forward_string = pSum.pz() > 0.;
}

void StringProcess::find_junction_leg(bool sign_color, std::vector<int> &col,
                                      Pythia8::Event &event_intermediate,
                                      Pythia8::Event &event_hadronize) {
  Pythia8::Vec4 pSum = event_hadronize[0].p();
  for (unsigned int j = 0; j < col.size(); j++) {
    if (col[j] == 0) {
      continue;
    }
    bool found_leg = false;
    while (!found_leg) {
      int ifound = -1;
      for (int i = 1; i < event_intermediate.size(); i++) {
        const int pdgid = event_intermediate[i].id();
        if (sign_color && col[j] == event_intermediate[i].col()) {
          logg[LPythia].debug("  col[", j, "] = ", col[j], " from i ", i, "(",
                              pdgid, ") found");
          ifound = i;
          col[j] = event_intermediate[i].acol();
          break;
        } else if (!sign_color && col[j] == event_intermediate[i].acol()) {
          logg[LPythia].debug("  acol[", j, "] = ", col[j], " from i ", i, "(",
                              pdgid, ") found");
          ifound = i;
          col[j] = event_intermediate[i].col();
          break;
        }
      }

      if (ifound < 0) {
        found_leg = true;
        if (event_intermediate.sizeJunction() == 0) {
          event_intermediate.list();
          event_intermediate.listJunctions();
          event_hadronize.list();
          event_hadronize.listJunctions();
          logg[LPythia].error("No parton with col = ", col[j],
                              " connected with junction leg ", j);
          throw std::runtime_error("Hard string could not be identified.");
        }
      } else {
        pSum += event_intermediate[ifound].p();
        // add a parton to the new event record.
        event_hadronize.append(event_intermediate[ifound]);
        // then remove from the original event record.
        event_intermediate.remove(ifound, ifound);
        logg[LPythia].debug(
            "Hard non-diff: event_intermediate reduces in size to ",
            event_intermediate.size());
        if (col[j] == 0) {
          found_leg = true;
        }
      }
    }
  }

  /* The zeroth entry of event record is supposed to have the information
   * on the whole system. Specify the total momentum and invariant mass. */
  event_hadronize[0].p(pSum);
  event_hadronize[0].m(pSum.mCalc());
}

// baryon-antibaryon annihilation
bool StringProcess::next_BBbarAnn() {
  const std::array<FourVector, 2> ustrcom = {FourVector(1., 0., 0., 0.),
                                             FourVector(1., 0., 0., 0.)};

  NpartFinal_ = 0;
  NpartString_[0] = 0;
  NpartString_[1] = 0;
  final_state_.clear();

  logg[LPythia].debug("Annihilation occurs between ", PDGcodes_[0], "+",
                      PDGcodes_[1], " at CM energy [GeV] ", sqrtsAB_);

  // check if the initial state is baryon-antibaryon pair.
  PdgCode baryon = PDGcodes_[0], antibaryon = PDGcodes_[1];
  if (baryon.baryon_number() == -1) {
    std::swap(baryon, antibaryon);
  }
  if (baryon.baryon_number() != 1 || antibaryon.baryon_number() != -1) {
    throw std::invalid_argument("Expected baryon-antibaryon pair.");
  }

  // Count how many qqbar combinations are possible for each quark type
  constexpr int n_q_types = 5;  // u, d, s, c, b
  std::vector<int> qcount_bar, qcount_antibar;
  std::vector<int> n_combinations;
  bool no_combinations = true;
  for (int i = 0; i < n_q_types; i++) {
    qcount_bar.push_back(baryon.net_quark_number(i + 1));
    qcount_antibar.push_back(-antibaryon.net_quark_number(i + 1));
    const int n_i = qcount_bar[i] * qcount_antibar[i];
    n_combinations.push_back(n_i);
    if (n_i > 0) {
      no_combinations = false;
    }
  }

  /* if it is a BBbar pair but there is no qqbar pair to annihilate,
   * nothing happens */
  if (no_combinations) {
    for (int i = 0; i < 2; i++) {
      NpartString_[i] = 1;
      ParticleData new_particle(ParticleType::find(PDGcodes_[i]));
      new_particle.set_4momentum(pcom_[i]);
      new_particle.set_cross_section_scaling_factor(1.);
      new_particle.set_formation_time(time_collision_);
      if (spin_interaction_type_ != SpinInteractionType::Off) {
        new_particle.set_unpolarized_spin_vector();
      }
      final_state_.push_back(new_particle);
    }
    NpartFinal_ = NpartString_[0] + NpartString_[1];
    return true;
  }

  // Select qqbar pair to annihilate and remove it away
  auto discrete_distr = random::discrete_dist<int>(n_combinations);
  const int q_annihilate = discrete_distr() + 1;
  qcount_bar[q_annihilate - 1]--;
  qcount_antibar[q_annihilate - 1]--;

  // Get the remaining quarks and antiquarks
  std::vector<int> remaining_quarks, remaining_antiquarks;
  for (int i = 0; i < n_q_types; i++) {
    for (int j = 0; j < qcount_bar[i]; j++) {
      remaining_quarks.push_back(i + 1);
    }
    for (int j = 0; j < qcount_antibar[i]; j++) {
      remaining_antiquarks.push_back(-(i + 1));
    }
  }
  assert(remaining_quarks.size() == 2);
  assert(remaining_antiquarks.size() == 2);

  const std::array<double, 2> mstr = {0.5 * sqrtsAB_, 0.5 * sqrtsAB_};

  // randomly select two quark-antiquark pairs
  if (random::uniform_int(0, 1) == 0) {
    std::swap(remaining_quarks[0], remaining_quarks[1]);
  }
  if (random::uniform_int(0, 1) == 0) {
    std::swap(remaining_antiquarks[0], remaining_antiquarks[1]);
  }
  // Make sure it satisfies kinematical threshold constraint
  bool kin_threshold_satisfied = true;
  for (int i = 0; i < 2; i++) {
    const double mstr_min =
        pythia_hadron_->particleData.m0(remaining_quarks[i]) +
        pythia_hadron_->particleData.m0(remaining_antiquarks[i]);
    if (mstr_min > mstr[i]) {
      kin_threshold_satisfied = false;
    }
  }
  if (!kin_threshold_satisfied) {
    return false;
  }
  // Fragment two strings
  for (int i = 0; i < 2; i++) {
    ParticleList new_intermediate_particles;

    ThreeVector evec = pcom_[i].threevec() / pcom_[i].threevec().abs();
    const int nfrag =
        fragment_string(remaining_quarks[i], remaining_antiquarks[i], mstr[i],
                        evec, true, false, new_intermediate_particles);
    if (nfrag <= 0) {
      NpartString_[i] = 0;
      return false;
    }
    // Set an unpolarized spin vector for the new intermediate particles
    if (spin_interaction_type_ != SpinInteractionType::Off) {
      for (ParticleData &intermediate_particle : new_intermediate_particles) {
        intermediate_particle.set_unpolarized_spin_vector();
      }
    }

    NpartString_[i] =
        append_final_state(new_intermediate_particles, ustrcom[i], evec);
  }
  NpartFinal_ = NpartString_[0] + NpartString_[1];
  return true;
}

void StringProcess::make_orthonormal_basis(
    ThreeVector &evec_polar, std::array<ThreeVector, 3> &evec_basis) {
  assert(std::fabs(evec_polar.sqr() - 1.) < really_small);

  if (std::abs(evec_polar.x3()) < (1. - 1.0e-8)) {
    double ex, ey, et;
    double theta, phi;

    // evec_basis[0] is set to be longitudinal direction
    evec_basis[0] = evec_polar;

    theta = std::acos(evec_basis[0].x3());

    ex = evec_basis[0].x1();
    ey = evec_basis[0].x2();
    et = std::sqrt(ex * ex + ey * ey);
    if (ey > 0.) {
      phi = std::acos(ex / et);
    } else {
      phi = -std::acos(ex / et);
    }

    /* The transverse plane is spanned
     * by evec_basis[1] and evec_basis[2]. */
    evec_basis[1].set_x1(std::cos(theta) * std::cos(phi));
    evec_basis[1].set_x2(std::cos(theta) * std::sin(phi));
    evec_basis[1].set_x3(-std::sin(theta));

    evec_basis[2].set_x1(-std::sin(phi));
    evec_basis[2].set_x2(std::cos(phi));
    evec_basis[2].set_x3(0.);
  } else {
    // if evec_polar is very close to the z axis
    if (evec_polar.x3() > 0.) {
      evec_basis[1] = ThreeVector(1., 0., 0.);
      evec_basis[2] = ThreeVector(0., 1., 0.);
      evec_basis[0] = ThreeVector(0., 0., 1.);
    } else {
      evec_basis[1] = ThreeVector(0., 1., 0.);
      evec_basis[2] = ThreeVector(1., 0., 0.);
      evec_basis[0] = ThreeVector(0., 0., -1.);
    }
  }

  assert(std::fabs(evec_basis[1] * evec_basis[2]) < really_small);
  assert(std::fabs(evec_basis[2] * evec_basis[0]) < really_small);
  assert(std::fabs(evec_basis[0] * evec_basis[1]) < really_small);
}

void StringProcess::compute_incoming_lightcone_momenta() {
  PPosA_ = (pcom_[0].x0() + evecBasisAB_[0] * pcom_[0].threevec()) * M_SQRT1_2;
  PNegA_ = (pcom_[0].x0() - evecBasisAB_[0] * pcom_[0].threevec()) * M_SQRT1_2;
  PPosB_ = (pcom_[1].x0() + evecBasisAB_[0] * pcom_[1].threevec()) * M_SQRT1_2;
  PNegB_ = (pcom_[1].x0() - evecBasisAB_[0] * pcom_[1].threevec()) * M_SQRT1_2;
}

void StringProcess::quarks_from_diquark(int diquark, int &q1, int &q2,
                                        int &deg_spin) {
  // The 4-digit pdg id should be diquark.
  assert((std::abs(diquark) > 1000) && (std::abs(diquark) < 5510) &&
         (std::abs(diquark) % 100 < 10));

  // The fourth digit corresponds to the spin degeneracy.
  deg_spin = std::abs(diquark) % 10;
  // Diquark (anti-diquark) is decomposed into two quarks (antiquarks).
  const int sign_anti = diquark > 0 ? 1 : -1;

  // Obtain two quarks (or antiquarks) from the first and second digit.
  q1 = sign_anti * (std::abs(diquark) - (std::abs(diquark) % 1000)) / 1000;
  q2 = sign_anti * (std::abs(diquark) % 1000 - deg_spin) / 100;
}

int StringProcess::diquark_from_quarks(int q1, int q2) {
  assert((q1 > 0 && q2 > 0) || (q1 < 0 && q2 < 0));
  if (std::abs(q1) < std::abs(q2)) {
    std::swap(q1, q2);
  }
  int diquark = std::abs(q1 * 1000 + q2 * 100);
  /* Adding spin degeneracy = 2S+1. For identical quarks spin cannot be 0
   * because of Pauli exclusion principle, so spin 1 is assumed. Otherwise
   * S = 0 with probability 1/4 and S = 1 with probability 3/4. */
  diquark += (q1 != q2 && random::uniform_int(0, 3) == 0) ? 1 : 3;
  return (q1 < 0) ? -diquark : diquark;
}

void StringProcess::make_string_ends(const PdgCode &pdg, int &idq1, int &idq2,
                                     double xi) {
  std::array<int, 3> quarks = pdg.quark_content();
  if (pdg.is_nucleon()) {
    // protons and neutrons treated seperately since single quarks is at a
    // different position in the PDG code
    if (pdg.charge() == 0) {  // (anti)neutron
      if (random::uniform(0., 1.) < xi) {
        idq1 = quarks[0];
        idq2 = diquark_from_quarks(quarks[1], quarks[2]);
      } else {
        idq1 = quarks[1];
        idq2 = diquark_from_quarks(quarks[0], quarks[2]);
      }
    } else {  // (anti)proton
      if (random::uniform(0., 1.) < xi) {
        idq1 = quarks[2];
        idq2 = diquark_from_quarks(quarks[0], quarks[1]);
      } else {
        idq1 = quarks[0];
        idq2 = diquark_from_quarks(quarks[1], quarks[2]);
      }
    }
  } else {
    if (pdg.is_meson()) {
      idq1 = quarks[1];
      idq2 = quarks[2];
      /* Some mesons with PDG id 11X are actually mixed state of uubar and
       * ddbar. have a random selection whether we have uubar or ddbar in this
       * case. */
      if (idq1 == 1 && idq2 == -1 && random::uniform_int(0, 1) == 0) {
        idq1 = 2;
        idq2 = -2;
      }
    } else {
      assert(pdg.is_baryon());
      // Get random quark to position 0
      std::swap(quarks[random::uniform_int(0, 2)], quarks[0]);
      idq1 = quarks[0];
      idq2 = diquark_from_quarks(quarks[1], quarks[2]);
    }
  }
  // Fulfil the convention: idq1 should be quark or anti-diquark
  if (idq1 < 0) {
    std::swap(idq1, idq2);
  }
}

int StringProcess::fragment_string(int idq1, int idq2, double mString,
                                   ThreeVector &evecLong, bool flip_string_ends,
                                   bool separate_fragment_baryon,
                                   ParticleList &intermediate_particles) {
  pythia_hadron_->event.reset();
  intermediate_particles.clear();

  logg[LPythia].debug("initial quark content for fragment_string : ", idq1,
                      ", ", idq2);
  logg[LPythia].debug("initial string mass (GeV) for fragment_string : ",
                      mString);
  // PDG id of quark constituents of string ends
  std::array<int, 2> idqIn;
  idqIn[0] = idq1;
  idqIn[1] = idq2;

  int bstring = 0;
  // constituent masses of the string
  std::array<double, 2> m_const;

  for (int i = 0; i < 2; i++) {
    // evaluate total baryon number of the string times 3
    bstring += pythia_hadron_->particleData.baryonNumberType(idqIn[i]);

    m_const[i] = pythia_hadron_->particleData.m0(idqIn[i]);
  }
  logg[LPythia].debug("baryon number of string times 3 : ", bstring);

  if (flip_string_ends && random::uniform_int(0, 1) == 0) {
    /* in the case where we flip the string ends,
     * we need to flip the longitudinal unit vector itself
     * since it is set to be direction of diquark (anti-quark)
     * or anti-diquark. */
    evecLong = -evecLong;
  }

  if (m_const[0] + m_const[1] > mString) {
    throw std::runtime_error("String fragmentation: m1 + m2 > mString");
  }
  Pythia8::Vec4 pSum = 0.;

  int number_of_fragments = 0;
  bool do_string_fragmentation = false;

  /* lightcone momenta p^+ and p^- of the string
   * p^{\pm} is defined as (E \pm p_{longitudinal}) / sqrts{2}. */
  double ppos_string_new, pneg_string_new;
  /* transverse momentum (and magnitude) acquired
   * by the the string to be fragmented */
  double QTrx_string_new, QTry_string_new, QTrn_string_new;
  // transverse mass of the string to be fragmented
  double mTrn_string_new;
  // mass of the string to be fragmented
  double mass_string_new;

  /* Transverse momentum to be added to the most forward hadron
   * from PYTHIA fragmentation */
  double QTrx_add_pos, QTry_add_pos;
  /* Transverse momentum to be added to the most backward hadron
   * from PYTHIA fragmentation */
  double QTrx_add_neg, QTry_add_neg;

  /* Set those transverse momenta to be zero.
   * This is the case when we solely rely on the PYTHIA fragmentation
   * procedure without separate fragmentation function.
   * In the case of separate fragmentation function for the leading baryon,
   * appropriate values will be assigned later. */
  QTrx_add_pos = 0.;
  QTry_add_pos = 0.;
  QTrx_add_neg = 0.;
  QTry_add_neg = 0.;

  std::array<ThreeVector, 3> evec_basis;
  make_orthonormal_basis(evecLong, evec_basis);

  if (separate_fragment_baryon && (std::abs(bstring) == 3) &&
      (mString > m_const[0] + m_const[1] + 1.)) {
    /* A separate fragmentation function will be used to the leading baryon,
     * if we have a baryonic string and the corresponding option is turned on.
     */
    int n_frag_prior;
    /* PDG id of fragmented hadrons
     * before switching to the PYTHIA fragmentation */
    std::vector<int> pdgid_frag_prior;
    /* four-momenta of fragmented hadrons
     * before switching to the PYTHIA fragmentation */
    std::vector<FourVector> momentum_frag_prior;

    // Transverse momentum px of the forward end of the string
    double QTrx_string_pos;
    // Transverse momentum px of the backward end of the string
    double QTrx_string_neg;
    // Transverse momentum py of the forward end of the string
    double QTry_string_pos;
    // Transverse momentum py of the backward end of the string
    double QTry_string_neg;
    /* Absolute value of the transverse momentum of
     * the forward end of the string */
    double QTrn_string_pos;
    /* Absolute value of the transverse momentum of
     * the backward end of the string */
    double QTrn_string_neg;

    std::array<double, 2> m_trans;

    // How many times we try to fragment leading baryon.
    const int niter_max = 10000;
    bool found_leading_baryon = false;
    for (int iiter = 0; iiter < niter_max; iiter++) {
      n_frag_prior = 0;
      pdgid_frag_prior.clear();
      momentum_frag_prior.clear();
      int n_frag = 0;
      std::vector<int> pdgid_frag;
      std::vector<FourVector> momentum_frag;
      // The original string is aligned in the logitudinal direction.
      ppos_string_new = mString * M_SQRT1_2;
      pneg_string_new = mString * M_SQRT1_2;
      // There is no transverse momentum at the original string ends.
      QTrx_string_pos = 0.;
      QTrx_string_neg = 0.;
      QTrx_string_new = 0.;
      QTry_string_pos = 0.;
      QTry_string_neg = 0.;
      QTry_string_new = 0.;
      // Constituent flavor at the forward (diquark) end of the string
      Pythia8::FlavContainer flav_string_pos =
          bstring > 0 ? Pythia8::FlavContainer(idq2)
                      : Pythia8::FlavContainer(idq1);
      // Constituent flavor at the backward (quark) end of the string
      Pythia8::FlavContainer flav_string_neg =
          bstring > 0 ? Pythia8::FlavContainer(idq1)
                      : Pythia8::FlavContainer(idq2);
      // Whether the first baryon from the forward end is fragmented
      bool found_forward_baryon = false;
      // Whether the first hadron from the forward end is fragmented
      bool done_forward_end = false;
      /* Whether energy of the string is depleted and the string
       * breaks into final two hadrons. */
      bool energy_used_up = false;
      while (!found_forward_baryon && !energy_used_up) {
        /* Keep fragmenting hadrons until
         * the first baryon is fragmented from the forward (diquark) end or
         * energy of the string is used up. */
        // Randomly select the string end from which the hadron is fragmented.
        bool from_forward = random::uniform_int(0, 1) == 0;
        /* The separate fragmentation function for the leading baryon kicks in
         * only if the first fragmented hadron from the forward (diquark) end
         * is a baryon. */
        n_frag = fragment_off_hadron(
            from_forward,
            separate_fragment_baryon && from_forward && !done_forward_end,
            evec_basis, ppos_string_new, pneg_string_new, QTrx_string_pos,
            QTrx_string_neg, QTry_string_pos, QTry_string_neg, flav_string_pos,
            flav_string_neg, pdgid_frag, momentum_frag);
        if (n_frag == 0) {
          /* If it fails to fragment hadron, start over from
           * the initial (baryonic) string configuration. */
          break;
        } else {
          QTrx_string_new = QTrx_string_pos + QTrx_string_neg;
          QTry_string_new = QTry_string_pos + QTry_string_neg;
          /* Quark (antiquark) constituents of the remaining string are
           * different from those of the original string.
           * Therefore, the constituent masses have to be updated. */
          idqIn[0] = bstring > 0 ? flav_string_neg.id : flav_string_pos.id;
          idqIn[1] = bstring > 0 ? flav_string_pos.id : flav_string_neg.id;
          for (int i = 0; i < 2; i++) {
            m_const[i] = pythia_hadron_->particleData.m0(idqIn[i]);
          }
          QTrn_string_pos = std::sqrt(QTrx_string_pos * QTrx_string_pos +
                                      QTry_string_pos * QTry_string_pos);
          QTrn_string_neg = std::sqrt(QTrx_string_neg * QTrx_string_neg +
                                      QTry_string_neg * QTry_string_neg);
          if (bstring > 0) {  // in the case of baryonic string
            /* Quark is coming from the newly produced backward qqbar pair
             * and therefore has transverse momentum, which is opposite to
             * that of the fragmented (backward) meson. */
            m_trans[0] = std::sqrt(m_const[0] * m_const[0] +
                                   QTrn_string_neg * QTrn_string_neg);
            /* Antiquark is coming from the newly produced forward qqbar pair
             * and therefore has transverse momentum, which is opposite to
             * that of the fragmented (leading) baryon. */
            m_trans[1] = std::sqrt(m_const[1] * m_const[1] +
                                   QTrn_string_pos * QTrn_string_pos);
          } else {  // in the case of anti-baryonic string
            /* Quark is coming from the newly produced forward qqbar pair
             * and therefore has transverse momentum, which is opposite to
             * that of the fragmented (leading) antibaryon. */
            m_trans[0] = std::sqrt(m_const[0] * m_const[0] +
                                   QTrn_string_pos * QTrn_string_pos);
            /* Antiquark is coming from the newly produced backward qqbar pair
             * and therefore has transverse momentum, which is opposite to
             * that of the fragmented (backward) meson. */
            m_trans[1] = std::sqrt(m_const[1] * m_const[1] +
                                   QTrn_string_neg * QTrn_string_neg);
          }
          done_forward_end = done_forward_end || from_forward;
          found_forward_baryon =
              found_forward_baryon ||
              (from_forward &&
               pythia_hadron_->particleData.isBaryon(pdgid_frag[0]));
        }
        if (n_frag == 2) {
          energy_used_up = true;
        }
        /* Add PDG id and four-momenta of fragmented hadrons
         * to the list if fragmentation is successful. */
        n_frag_prior += n_frag;
        for (int i_frag = 0; i_frag < n_frag; i_frag++) {
          pdgid_frag_prior.push_back(pdgid_frag[i_frag]);
          momentum_frag_prior.push_back(momentum_frag[i_frag]);
        }
      }
      if (n_frag == 0) {
        continue;
      } else {
        if (n_frag == 1) {
          // Compute transverse mass and momentum of the remaining string.
          double mTsqr_string = 2. * ppos_string_new * pneg_string_new;
          mTrn_string_new = std::sqrt(mTsqr_string);
          QTrn_string_new = std::sqrt(QTrx_string_new * QTrx_string_new +
                                      QTry_string_new * QTry_string_new);
          if (mTrn_string_new < QTrn_string_new) {
            /* If transverse mass is lower than transverse momentum,
             * start over. */
            found_leading_baryon = false;
          } else {
            // Compute mass of the remaining string.
            mass_string_new =
                std::sqrt(mTsqr_string - QTrn_string_new * QTrn_string_new);
            /* Proceed only if the string mass is large enough to call
             * PYTHIA fragmentation routine.
             * Otherwise, start over. */
            if (mass_string_new > m_const[0] + m_const[1]) {
              do_string_fragmentation = true;
              found_leading_baryon = true;
              QTrx_add_pos = QTrx_string_pos;
              QTry_add_pos = QTry_string_pos;
              QTrx_add_neg = QTrx_string_neg;
              QTry_add_neg = QTry_string_neg;
            } else {
              found_leading_baryon = false;
            }
          }
        } else if (n_frag == 2) {
          /* If the string ended up breaking into final two hadrons,
           * there is no need to perform PYTHIA fragmentation. */
          do_string_fragmentation = false;
          found_leading_baryon = true;
        }
      }

      if (found_leading_baryon) {
        break;
      }
    }
    if (found_leading_baryon) {
      /* If the kinematics makes sense, add fragmented hadrons so far
       * to the intermediate particle list. */
      for (int i_frag = 0; i_frag < n_frag_prior; i_frag++) {
        logg[LPythia].debug("appending the the fragmented hadron ",
                            pdgid_frag_prior[i_frag],
                            " to the intermediate particle list.");

        bool found_ptype = append_intermediate_list(pdgid_frag_prior[i_frag],
                                                    momentum_frag_prior[i_frag],
                                                    intermediate_particles);
        if (!found_ptype) {
          logg[LPythia].error("PDG ID ", pdgid_frag_prior[i_frag],
                              " should exist in ParticleType.");
          throw std::runtime_error("string fragmentation failed.");
        }
        number_of_fragments++;
      }
    } else {
      /* If it is not possible to find the leading baryon with appropriate
       * kinematics after trying many times, return failure (no hadron). */
      return 0;
    }

    if (do_string_fragmentation) {
      mTrn_string_new = std::sqrt(2. * ppos_string_new * pneg_string_new);
      // lightcone momentum p^+ of the quark constituents on the string ends
      std::array<double, 2> ppos_parton;
      // lightcone momentum p^- of the quark constituents on the string ends
      std::array<double, 2> pneg_parton;

      /* lightcone momenta of the string ends (quark and antiquark)
       * First, obtain ppos_parton[0] and ppos_parton[1]
       * (p^+ of quark and antiquark) by solving the following equations.
       * ppos_string_new = ppos_parton[0] + ppos_parton[1]
       * pneg_string_new = 0.5 * m_trans[0] * m_trans[0] / ppos_parton[0] +
       *                   0.5 * m_trans[1] * m_trans[1] / ppos_parton[1] */
      const double pb_const =
          (mTrn_string_new * mTrn_string_new + m_trans[0] * m_trans[0] -
           m_trans[1] * m_trans[1]) /
          (4. * pneg_string_new);
      const double pc_const =
          0.5 * m_trans[0] * m_trans[0] * ppos_string_new / pneg_string_new;
      ppos_parton[0] = pb_const + (bstring > 0 ? -1. : 1.) *
                                      std::sqrt(pb_const * pb_const - pc_const);
      ppos_parton[1] = ppos_string_new - ppos_parton[0];
      /* Then, compute pneg_parton[0] and pneg_parton[1]
       * (p^- of quark and antiquark) from the dispersion relation.
       * 2 p^+ p^- = m_transverse^2 */
      for (int i = 0; i < 2; i++) {
        pneg_parton[i] = 0.5 * m_trans[i] * m_trans[i] / ppos_parton[i];
      }

      const int status = 1;
      int color, anticolor;
      ThreeVector three_mom;
      ThreeVector transverse_mom;
      Pythia8::Vec4 pquark;

      // quark end of the remaining (mesonic) string
      color = 1;
      anticolor = 0;
      /* In the case of baryonic string,
       * Quark is coming from the backward end of the remaining string.
       * Transverse momentum of the backward end is subtracted
       * at this point to keep the remaining string aligned in
       * the original longitudinal direction.
       * It will be added to the most backward hadron from
       * PYTHIA fragmentation.
       * In the case of anti-baryonic string,
       * Quark is coming from the forward end of the remaining string.
       * Transverse momentum of the forward end is subtracted
       * at this point to keep the remaining string aligned in
       * the original longitudinal direction.
       * It will be added to the most forward hadron from
       * PYTHIA fragmentation. */
      transverse_mom =
          bstring > 0 ? evec_basis[1] * (QTrx_string_neg - QTrx_add_neg) +
                            evec_basis[2] * (QTry_string_neg - QTry_add_neg)
                      : evec_basis[1] * (QTrx_string_pos - QTrx_add_pos) +
                            evec_basis[2] * (QTry_string_pos - QTry_add_pos);
      three_mom =
          evec_basis[0] * (ppos_parton[0] - pneg_parton[0]) * M_SQRT1_2 +
          transverse_mom;
      const double E_quark =
          std::sqrt(m_const[0] * m_const[0] + three_mom.sqr());
      pquark = set_Vec4(E_quark, three_mom);
      pSum += pquark;
      pythia_hadron_->event.append(idqIn[0], status, color, anticolor, pquark,
                                   m_const[0]);

      // antiquark end of the remaining (mesonic) string
      color = 0;
      anticolor = 1;
      /* In the case of baryonic string,
       * Antiquark is coming from the forward end of the remaining string.
       * Transverse momentum of the forward end is subtracted
       * at this point to keep the remaining string aligned in
       * the original longitudinal direction.
       * It will be added to the most forward hadron from
       * PYTHIA fragmentation.
       * In the case of anti-baryonic string,
       * Antiquark is coming from the backward end of the remaining string.
       * Transverse momentum of the backward end is subtracted
       * at this point to keep the remaining string aligned in
       * the original longitudinal direction.
       * It will be added to the most backward hadron from
       * PYTHIA fragmentation. */
      transverse_mom =
          bstring > 0 ? evec_basis[1] * (QTrx_string_pos - QTrx_add_pos) +
                            evec_basis[2] * (QTry_string_pos - QTry_add_pos)
                      : evec_basis[1] * (QTrx_string_neg - QTrx_add_neg) +
                            evec_basis[2] * (QTry_string_neg - QTry_add_neg);
      three_mom =
          evec_basis[0] * (ppos_parton[1] - pneg_parton[1]) * M_SQRT1_2 +
          transverse_mom;
      const double E_antiq =
          std::sqrt(m_const[1] * m_const[1] + three_mom.sqr());
      pquark = set_Vec4(E_antiq, three_mom);
      pSum += pquark;
      pythia_hadron_->event.append(idqIn[1], status, color, anticolor, pquark,
                                   m_const[1]);
    }
  } else {
    do_string_fragmentation = true;

    ppos_string_new = mString * M_SQRT1_2;
    pneg_string_new = mString * M_SQRT1_2;
    QTrx_string_new = 0.;
    QTry_string_new = 0.;
    QTrn_string_new = 0.;
    mTrn_string_new = mString;
    mass_string_new = mString;

    /* diquark (anti-quark) with PDG id idq2 is going in the direction of
     * evecLong.
     * quark with PDG id idq1 is going in the direction opposite to evecLong. */
    double sign_direction = 1.;
    if (bstring == -3) {  // anti-baryonic string
      /* anti-diquark with PDG id idq1 is going in the direction of evecLong.
       * anti-quark with PDG id idq2 is going in the direction
       * opposite to evecLong. */
      sign_direction = -1;
    }

    // evaluate momenta of quarks
    const double pCMquark = pCM(mString, m_const[0], m_const[1]);
    const double E1 = std::sqrt(m_const[0] * m_const[0] + pCMquark * pCMquark);
    const double E2 = std::sqrt(m_const[1] * m_const[1] + pCMquark * pCMquark);

    ThreeVector direction = sign_direction * evecLong;

    // For status and (anti)color see \iref{Sjostrand:2007gs}.
    const int status1 = 1, color1 = 1, anticolor1 = 0;
    Pythia8::Vec4 pquark = set_Vec4(E1, -direction * pCMquark);
    pSum += pquark;
    pythia_hadron_->event.append(idqIn[0], status1, color1, anticolor1, pquark,
                                 m_const[0]);

    const int status2 = 1, color2 = 0, anticolor2 = 1;
    pquark = set_Vec4(E2, direction * pCMquark);
    pSum += pquark;
    pythia_hadron_->event.append(idqIn[1], status2, color2, anticolor2, pquark,
                                 m_const[1]);
  }

  if (do_string_fragmentation) {
    logg[LPythia].debug("fragmenting a string with ", idqIn[0], ", ", idqIn[1]);
    // implement PYTHIA fragmentation
    pythia_hadron_->event[0].p(pSum);
    pythia_hadron_->event[0].m(pSum.mCalc());
    bool successful_hadronization = pythia_hadron_->next();
    // update_info();
    if (!successful_hadronization) {
      return 0;
    }

    /* Add transverse momenta of string ends to the most forward and
     * backward hadrons from PYTHIA fragmentation. */
    bool successful_kinematics = remake_kinematics_fragments(
        pythia_hadron_->event, evec_basis, ppos_string_new, pneg_string_new,
        QTrx_string_new, QTry_string_new, QTrx_add_pos, QTry_add_pos,
        QTrx_add_neg, QTry_add_neg);
    if (!successful_kinematics) {
      return 0;
    }

    for (int ipyth = 0; ipyth < pythia_hadron_->event.size(); ipyth++) {
      if (!pythia_hadron_->event[ipyth].isFinal()) {
        continue;
      }
      int pythia_id = pythia_hadron_->event[ipyth].id();
      /* K_short and K_long need are converted to K0
       * since SMASH only knows K0 */
      convert_KaonLS(pythia_id);
      FourVector momentum(
          pythia_hadron_->event[ipyth].e(), pythia_hadron_->event[ipyth].px(),
          pythia_hadron_->event[ipyth].py(), pythia_hadron_->event[ipyth].pz());
      logg[LPythia].debug("appending the fragmented hadron ", pythia_id,
                          " to the intermediate particle list.");
      bool found_ptype =
          append_intermediate_list(pythia_id, momentum, intermediate_particles);
      if (!found_ptype) {
        logg[LPythia].warn("PDG ID ", pythia_id,
                           " does not exist in ParticleType - start over.");
        intermediate_particles.clear();
        return 0;
      }

      number_of_fragments++;
    }
  }
  return number_of_fragments;
}

int StringProcess::fragment_off_hadron(
    bool from_forward, bool separate_fragment_baryon,
    std::array<ThreeVector, 3> &evec_basis, double &ppos_string,
    double &pneg_string, double &QTrx_string_pos, double &QTrx_string_neg,
    double &QTry_string_pos, double &QTry_string_neg,
    Pythia8::FlavContainer &flav_string_pos,
    Pythia8::FlavContainer &flav_string_neg, std::vector<int> &pdgid_frag,
    std::vector<FourVector> &momentum_frag) {
  /* How many times we try to find flavor of qqbar pair and corresponding
   * hadronic species */
  const int n_try = 10;
  pdgid_frag.clear();
  momentum_frag.clear();

  if (ppos_string < 0. || pneg_string < 0.) {
    throw std::runtime_error("string has a negative lightcone momentum.");
  }
  double mTsqr_string = 2. * ppos_string * pneg_string;
  // Transverse mass of the original string
  double mTrn_string = std::sqrt(mTsqr_string);
  // Total transverse momentum of the original string
  double QTrx_string_tot = QTrx_string_pos + QTrx_string_neg;
  double QTry_string_tot = QTry_string_pos + QTry_string_neg;
  double QTsqr_string_tot = std::fabs(QTrx_string_tot * QTrx_string_tot) +
                            std::fabs(QTry_string_tot * QTry_string_tot);
  double QTrn_string_tot = std::sqrt(QTsqr_string_tot);
  if (mTrn_string < QTrn_string_tot) {
    return 0;
  }
  // Mass of the original string
  double mass_string = std::sqrt(mTsqr_string - QTsqr_string_tot);
  logg[LPythia].debug("  Fragment off one hadron from a string ( ",
                      flav_string_pos.id, " , ", flav_string_neg.id,
                      " ) with mass ", mass_string, " GeV.");

  // Take relevant parameters from PYTHIA.
  const double sigma_qt_frag = pythia_hadron_->parm("StringPT:sigma");
  const double stop_string_mass =
      pythia_hadron_->parm("StringFragmentation:stopMass");
  const double stop_string_smear =
      pythia_hadron_->parm("StringFragmentation:stopSmear");

  // Enhance the width of transverse momentum with certain probability
  const double prob_enhance_qt =
      pythia_hadron_->parm("StringPT:enhancedFraction");
  double fac_enhance_qt;
  if (random::uniform(0., 1.) < prob_enhance_qt) {
    fac_enhance_qt = pythia_hadron_->parm("StringPT:enhancedWidth");
  } else {
    fac_enhance_qt = 1.;
  }

  /* Sample the transverse momentum of newly created quark-antiquark
   * (or diquark-antidiquark) pair.
   * Note that one constituent carries QT_new while -QT_new is carried by
   * another.
   * The former one is taken by the (first) fragmented hadron and
   * the later one will be assigned to the remaining string or
   * taken by the second fragmented hadron. */
  double QTrx_new =
      random::normal(0., fac_enhance_qt * sigma_qt_frag * M_SQRT1_2);
  double QTry_new =
      random::normal(0., fac_enhance_qt * sigma_qt_frag * M_SQRT1_2);
  logg[LPythia].debug("  Transverse momentum (", QTrx_new, ", ", QTry_new,
                      ") GeV selected for the new qqbar pair.");

  /* Determine the transverse momentum of the (first) fragmented hadron.
   * Transverse momentum of hadron =
   * QT_string (of the string end) +
   * QT_new (of one of the quark-antiquark pair).
   * If the first hadron is fragmented from the forward (backward) end
   * of a string, then transverse momentum carried by the forward (backward)
   * end is taken. */
  double QTrx_had_1st =
      from_forward ? QTrx_string_pos + QTrx_new : QTrx_string_neg + QTrx_new;
  double QTry_had_1st =
      from_forward ? QTry_string_pos + QTry_new : QTry_string_neg + QTry_new;
  double QTrn_had_1st =
      std::sqrt(QTrx_had_1st * QTrx_had_1st + QTry_had_1st * QTry_had_1st);

  // PDG id of the (first) fragmented hadron
  int pdgid_had_1st = 0;
  // Mass of the (first) fragmented hadron
  double mass_had_1st = 0.;
  /* Constituent flavor of the original string end,
   * at which the (first) hadron is fragmented. */
  Pythia8::FlavContainer flav_old =
      from_forward ? flav_string_pos : flav_string_neg;
  /* Constituent flavor of newly created quark-antiquark pair,
   * which is taken by the (first) fragmented hadron.
   * Antiparticle of this flavor will be assigned to the remaining string,
   * or taken by the second fragmented hadron. */
  Pythia8::FlavContainer flav_new = Pythia8::FlavContainer(0);
  /* Sample flavor of the quark-antiquark (or diquark-antidiquark) pair
   * and combine with that of the original string end to find the hadronic
   * species. */
  for (int i_try = 0; i_try < n_try; i_try++) {
    // Sample the new flavor.
    flav_new = pythia_stringflav_.pick(flav_old);
    // Combine to get the PDG id of hadron.
    pdgid_had_1st = pythia_stringflav_.combine(flav_old, flav_new);
    if (pdgid_had_1st != 0) {
      // If the PDG id is found, determine mass.
      mass_had_1st = pythia_hadron_->particleData.mSel(pdgid_had_1st);
      logg[LPythia].debug("    number of tries of flavor selection : ",
                          i_try + 1, " in StringProcess::fragment_off_hadron.");
      break;
    }
  }
  if (pdgid_had_1st == 0) {
    return 0;
  }
  logg[LPythia].debug("  New flavor ", flav_new.id,
                      " selected for the string end with ", flav_old.id);
  logg[LPythia].debug("  PDG id ", pdgid_had_1st,
                      " selected for the (first) fragmented hadron.");
  bool had_1st_baryon = pythia_hadron_->particleData.isBaryon(pdgid_had_1st);
  // Transverse mass of the (first) fragmented hadron
  double mTrn_had_1st =
      std::sqrt(mass_had_1st * mass_had_1st + QTrn_had_1st * QTrn_had_1st);
  logg[LPythia].debug("  Transverse momentum (", QTrx_had_1st, ", ",
                      QTry_had_1st,
                      ") GeV selected for the (first) fragmented hadron.");

  /* Compute the mass threshold to continue string fragmentation.
   * This formula is taken from StringFragmentation::energyUsedUp
   * in StringFragmentation.cc of PYTHIA 8. */
  const double mass_min_to_continue =
      (stop_string_mass + pythia_hadron_->particleData.m0(flav_new.id) +
       pythia_hadron_->particleData.m0(flav_string_pos.id) +
       pythia_hadron_->particleData.m0(flav_string_neg.id)) *
      (1. + (2. * random::uniform(0., 1.) - 1.) * stop_string_smear);
  /* If the string mass is lower than that threshold,
   * the string breaks into the last two hadrons. */
  bool string_into_final_two = mass_string < mass_min_to_continue;
  if (string_into_final_two) {
    logg[LPythia].debug("  The string mass is below the mass threshold ",
                        mass_min_to_continue,
                        " GeV : finishing with two hadrons.");
  }

  // Lightcone momentum of the (first) fragmented hadron
  double ppos_had_1st = 0.;
  double pneg_had_1st = 0.;

  /* Whether the string end, at which the (first) hadron is fragmented,
   * had a diquark or antidiquark */
  bool from_diquark_end =
      from_forward ? pythia_hadron_->particleData.isDiquark(flav_string_pos.id)
                   : pythia_hadron_->particleData.isDiquark(flav_string_neg.id);
  // Whether the forward end of the string has a diquark
  bool has_diquark_pos =
      pythia_hadron_->particleData.isDiquark(flav_string_pos.id);

  int n_frag = 0;
  if (string_into_final_two) {
    /* In the case of a string breaking into the last two hadrons,
     * we determine species, mass and four-momentum of the second hadron. */
    // PDG id of the second fragmented hadron
    int pdgid_had_2nd = 0.;
    // Mass of the second fragmented hadron
    double mass_had_2nd = 0.;
    /* Constituent flavor of newly created quark-antiquark pair,
     * which is taken by the second fragmented hadron. */
    Pythia8::FlavContainer flav_new2 = Pythia8::FlavContainer(0);
    flav_new2.anti(flav_new);
    /* Getting a hadron from diquark and antidiquark does not always work.
     * So, if this is the case, start over. */
    if (pythia_hadron_->particleData.isDiquark(flav_string_neg.id) &&
        pythia_hadron_->particleData.isDiquark(flav_new2.id) && from_forward) {
      return 0;
    }
    if (pythia_hadron_->particleData.isDiquark(flav_string_pos.id) &&
        pythia_hadron_->particleData.isDiquark(flav_new2.id) && !from_forward) {
      return 0;
    }
    for (int i_try = 0; i_try < n_try; i_try++) {
      // Combine to get the PDG id of the second hadron.
      pdgid_had_2nd =
          from_forward ? pythia_stringflav_.combine(flav_string_neg, flav_new2)
                       : pythia_stringflav_.combine(flav_string_pos, flav_new2);
      if (pdgid_had_2nd != 0) {
        // If the PDG id is found, determine mass.
        mass_had_2nd = pythia_hadron_->particleData.mSel(pdgid_had_2nd);
        break;
      }
    }
    if (pdgid_had_2nd == 0) {
      return 0;
    }
    logg[LPythia].debug("  PDG id ", pdgid_had_2nd,
                        " selected for the (second) fragmented hadron.");
    bool had_2nd_baryon = pythia_hadron_->particleData.isBaryon(pdgid_had_2nd);

    /* Determine transverse momentum carried by the second hadron.
     * If the first hadron fragmented from the forward (backward) end
     * of a string, transvere momentum at the backward (forward) end will
     * contribute.
     * Transverse momentum of newly created constituent
     * must be added as well. */
    double QTrx_had_2nd =
        from_forward ? QTrx_string_neg - QTrx_new : QTrx_string_pos - QTrx_new;
    double QTry_had_2nd =
        from_forward ? QTry_string_neg - QTry_new : QTry_string_pos - QTry_new;
    double QTrn_had_2nd =
        std::sqrt(QTrx_had_2nd * QTrx_had_2nd + QTry_had_2nd * QTry_had_2nd);
    double mTrn_had_2nd =
        std::sqrt(mass_had_2nd * mass_had_2nd + QTrn_had_2nd * QTrn_had_2nd);
    logg[LPythia].debug("  Transverse momentum (", QTrx_had_2nd, ", ",
                        QTry_had_2nd,
                        ") GeV selected for the (second) fragmented hadron.");

    double ppos_had_2nd = 0.;
    double pneg_had_2nd = 0.;

    /* Compute lightcone momenta of the final two hadrons.
     * If the fragmentation begins at the forward (backward) end of a string,
     * the first (second) hadron is the forward one and the second (first)
     * hadron is the backward one. */
    bool found_kinematics =
        from_forward
            ? make_lightcone_final_two(
                  separate_fragment_baryon && has_diquark_pos && had_1st_baryon,
                  ppos_string, pneg_string, mTrn_had_1st, mTrn_had_2nd,
                  ppos_had_1st, ppos_had_2nd, pneg_had_1st, pneg_had_2nd)
            : make_lightcone_final_two(
                  separate_fragment_baryon && has_diquark_pos && had_2nd_baryon,
                  ppos_string, pneg_string, mTrn_had_2nd, mTrn_had_1st,
                  ppos_had_2nd, ppos_had_1st, pneg_had_2nd, pneg_had_1st);
    if (!found_kinematics) {
      return 0;
    }

    // The entire string breaks into hadrons, so there is no momentum left.
    ppos_string = 0.;
    pneg_string = 0.;
    QTrx_string_pos = 0.;
    QTry_string_pos = 0.;

    // Add the first hadron to the list.
    pdgid_frag.push_back(pdgid_had_1st);
    FourVector mom_had_1st = FourVector(
        (ppos_had_1st + pneg_had_1st) * M_SQRT1_2,
        evec_basis[0] * (ppos_had_1st - pneg_had_1st) * M_SQRT1_2 +
            evec_basis[1] * QTrx_had_1st + evec_basis[2] * QTry_had_1st);
    momentum_frag.push_back(mom_had_1st);

    // Add the second hadron to the list.
    pdgid_frag.push_back(pdgid_had_2nd);
    FourVector mom_had_2nd = FourVector(
        (ppos_had_2nd + pneg_had_2nd) * M_SQRT1_2,
        evec_basis[0] * (ppos_had_2nd - pneg_had_2nd) * M_SQRT1_2 +
            evec_basis[1] * QTrx_had_2nd + evec_basis[2] * QTry_had_2nd);
    momentum_frag.push_back(mom_had_2nd);

    n_frag += 2;
  } else {
    /* If the string mass is large enough (larger than the threshold),
     * perform the normal fragmentation.
     * Different sets of parameters for the LUND fragmentation function
     * are used, depending on whether the (first) fragmented hadrons is
     * the leading baryon. */
    double stringz_a_use, stringz_b_use;
    /* If the firstly fragmented hadron from the diquark end is a baryon,
     * it can be considered to be the leading baryon. */
    if (separate_fragment_baryon && from_diquark_end && had_1st_baryon) {
      stringz_a_use = stringz_a_leading_;
      stringz_b_use = stringz_b_leading_;
    } else {
      stringz_a_use = stringz_a_produce_;
      stringz_b_use = stringz_b_produce_;
    }

    /* Sample the lightcone momentum fraction from
     * the LUND fragmentation function. */
    double xfrac = sample_zLund(stringz_a_use, stringz_b_use, mTrn_had_1st);
    if (from_forward) {
      /* If the (first) hadron is fragmented from the forward end,
       * it is the lightcone momentum fraction of p^+ */
      ppos_had_1st = xfrac * ppos_string;
      pneg_had_1st = 0.5 * mTrn_had_1st * mTrn_had_1st / ppos_had_1st;
      if (pneg_had_1st > pneg_string) {
        return 0;
      }
    } else {
      /* If the (first) hadron is fragmented from the backward end,
       * it is the lightcone momentum fraction of p^- */
      pneg_had_1st = xfrac * pneg_string;
      ppos_had_1st = 0.5 * mTrn_had_1st * mTrn_had_1st / pneg_had_1st;
      if (ppos_had_1st > ppos_string) {
        return 0;
      }
    }

    // Add PDG id and four-momentum of the (first) hadron to the list.
    pdgid_frag.push_back(pdgid_had_1st);
    FourVector mom_had_1st = FourVector(
        (ppos_had_1st + pneg_had_1st) * M_SQRT1_2,
        evec_basis[0] * (ppos_had_1st - pneg_had_1st) * M_SQRT1_2 +
            evec_basis[1] * QTrx_had_1st + evec_basis[2] * QTry_had_1st);
    momentum_frag.push_back(mom_had_1st);

    // Update lightcone momentum of the string.
    ppos_string -= ppos_had_1st;
    pneg_string -= pneg_had_1st;
    /* Update flavor and transverse momentum of the string end,
     * from which the (first) hadron is fragmented.
     * Flavor of the new string end is antiparticle of
     * the constituent taken by the hadron.
     * Transverse momentum of the new string end is opposite to that
     * of the constituent taken by the hadron. */
    if (from_forward) {
      /* Update the forward end of the string
       * if hadron is fragmented from there. */
      flav_string_pos.anti(flav_new);
      QTrx_string_pos = -QTrx_new;
      QTry_string_pos = -QTry_new;
    } else {
      /* Update the backward end of the string
       * if hadron is fragmented from there. */
      flav_string_neg.anti(flav_new);
      QTrx_string_neg = -QTrx_new;
      QTry_string_neg = -QTry_new;
    }

    n_frag += 1;
  }

  return n_frag;
}

int StringProcess::get_hadrontype_from_quark(int idq1, int idq2) {
  const int baryon_number =
      pythia_hadron_->particleData.baryonNumberType(idq1) +
      pythia_hadron_->particleData.baryonNumberType(idq2);

  int pdgid_hadron = 0;
  /* PDG id of the leading baryon from valence quark constituent.
   * First, try with the PYTHIA machinary. */
  Pythia8::FlavContainer flav1 = Pythia8::FlavContainer(idq1);
  Pythia8::FlavContainer flav2 = Pythia8::FlavContainer(idq2);
  const int n_try = 10;
  for (int i_try = 0; i_try < n_try; i_try++) {
    pdgid_hadron = pythia_stringflav_.combine(flav1, flav2);
    if (pdgid_hadron != 0) {
      return pdgid_hadron;
    }
  }

  /* If PYTHIA machinary does not work, determine type of the leading baryon
   * based on the quantum numbers and mass. */

  // net quark number of d, u, s, c and b flavors
  std::array<int, 5> frag_net_q;
  /* Evaluate total net quark number of baryon (antibaryon)
   * from the valence quark constituents. */
  for (int iq = 0; iq < 5; iq++) {
    int nq1 =
        pythia_hadron_->particleData.nQuarksInCode(std::abs(idq1), iq + 1);
    int nq2 =
        pythia_hadron_->particleData.nQuarksInCode(std::abs(idq2), iq + 1);
    nq1 = idq1 > 0 ? nq1 : -nq1;
    nq2 = idq2 > 0 ? nq2 : -nq2;
    frag_net_q[iq] = nq1 + nq2;
  }
  const int frag_iso3 = frag_net_q[1] - frag_net_q[0];
  const int frag_strange = -frag_net_q[2];
  const int frag_charm = frag_net_q[3];
  const int frag_bottom = -frag_net_q[4];
  logg[LPythia].debug("  conserved charges : iso3 = ", frag_iso3,
                      ", strangeness = ", frag_strange,
                      ", charmness = ", frag_charm,
                      ", bottomness = ", frag_bottom);

  std::vector<int> pdgid_possible;
  std::vector<double> weight_possible;
  std::vector<double> weight_summed;
  /* loop over hadronic species
   * Any hadron with the same valence quark contents is allowed and
   * the probability goes like spin degeneracy over mass. */
  for (auto &ptype : ParticleType::list_all()) {
    if (!ptype.is_hadron()) {
      continue;
    }
    const int pdgid = ptype.pdgcode().get_decimal();
    if ((pythia_hadron_->particleData.isParticle(pdgid)) &&
        (baryon_number == 3 * ptype.pdgcode().baryon_number()) &&
        (frag_iso3 == ptype.pdgcode().isospin3()) &&
        (frag_strange == ptype.pdgcode().strangeness()) &&
        (frag_charm == ptype.pdgcode().charmness()) &&
        (frag_bottom == ptype.pdgcode().bottomness())) {
      const int spin_degeneracy = ptype.pdgcode().spin_degeneracy();
      const double mass_pole = ptype.mass();
      const double weight = static_cast<double>(spin_degeneracy) / mass_pole;
      pdgid_possible.push_back(pdgid);
      weight_possible.push_back(weight);

      logg[LPythia].debug("  PDG id ", pdgid, " is possible with weight ",
                          weight);
    }
  }
  const int n_possible = pdgid_possible.size();
  weight_summed.push_back(0.);
  for (int i = 0; i < n_possible; i++) {
    weight_summed.push_back(weight_summed[i] + weight_possible[i]);
  }

  /* Sample baryon (antibaryon) specie,
   * which is fragmented from the leading diquark (anti-diquark). */
  const double uspc = random::uniform(0., weight_summed[n_possible]);
  for (int i = 0; i < n_possible; i++) {
    if ((uspc >= weight_summed[i]) && (uspc < weight_summed[i + 1])) {
      return pdgid_possible[i];
    }
  }

  return 0;
}

int StringProcess::get_resonance_from_quark(int idq1, int idq2, double mass) {
  // if the mass is too low, return 0 (failure).
  if (mass < pion_mass) {
    return 0;
  }

  /* It checks whether one has a valid input
   * the string ends. */

  // idq1 is supposed to be a quark or anti-diquark.
  bool end1_is_quark = idq1 > 0 && pythia_hadron_->particleData.isQuark(idq1);
  bool end1_is_antidiq =
      idq1 < 0 && pythia_hadron_->particleData.isDiquark(idq1);
  // idq2 is supposed to be a anti-quark or diquark.
  bool end2_is_antiq = idq2 < 0 && pythia_hadron_->particleData.isQuark(idq2);
  bool end2_is_diquark =
      idq2 > 0 && pythia_hadron_->particleData.isDiquark(idq2);

  int baryon;
  if (end1_is_quark) {
    if (end2_is_antiq) {
      // we have a mesonic resonance from a quark-antiquark pair.
      baryon = 0;
    } else if (end2_is_diquark) {
      // we have a baryonic resonance from a quark-diquark pair.
      baryon = 1;
    } else {
      return 0;
    }
  } else if (end1_is_antidiq) {
    if (end2_is_antiq) {
      // we have a antibaryonic resonance from a antiquark-antidiquark pair.
      baryon = -1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }

  /* array for the net quark numbers of the constituents.
   * net_qnumber[0, 1, 2, 3 and 4] correspond respectively to
   * d, u, s, c, b quark flavors. */
  std::array<int, 5> net_qnumber;
  for (int iflav = 0; iflav < 5; iflav++) {
    net_qnumber[iflav] = 0;

    int qnumber1 =
        pythia_hadron_->particleData.nQuarksInCode(std::abs(idq1), iflav + 1);
    if (idq1 < 0) {
      // anti-diquark gets an extra minus sign.
      qnumber1 = -qnumber1;
    }
    net_qnumber[iflav] += qnumber1;

    int qnumber2 =
        pythia_hadron_->particleData.nQuarksInCode(std::abs(idq2), iflav + 1);
    if (idq2 < 0) {
      // anti-quark gets an extra minus sign.
      qnumber2 = -qnumber2;
    }
    net_qnumber[iflav] += qnumber2;
  }

  // List of PDG ids of resonances with the same quantum number.
  std::vector<int> pdgid_possible;
  // Corresponding mass differences.
  std::vector<double> mass_diff;
  for (auto &ptype : ParticleType::list_all()) {
    if (!ptype.is_hadron() || ptype.is_stable() ||
        ptype.pdgcode().baryon_number() != baryon) {
      // Only resonances with the same baryon number are considered.
      continue;
    }
    const int pdgid = ptype.pdgcode().get_decimal();
    const double mass_min = ptype.min_mass_spectral();
    if (mass < mass_min) {
      // A resoance with mass lower than its minimum threshold is not allowed.
      continue;
    }

    if (ptype.pdgcode().isospin3() != net_qnumber[1] - net_qnumber[0]) {
      // check isospin3.
      continue;
    }
    if (ptype.pdgcode().strangeness() != -net_qnumber[2]) {
      // check strangeness.
      continue;
    }
    if (ptype.pdgcode().charmness() != net_qnumber[3]) {
      // check charmness.
      continue;
    }
    if (ptype.pdgcode().bottomness() != -net_qnumber[4]) {
      // check bottomness.
      continue;
    }

    const double mass_pole = ptype.mass();
    // Add the PDG id and mass difference to the vector array.
    pdgid_possible.push_back(pdgid);
    mass_diff.push_back(mass - mass_pole);
  }

  const int n_res = pdgid_possible.size();
  if (n_res == 0) {
    // If there is no possible resonance found, return 0 (failure).
    return 0;
  }

  int ires_closest = 0;
  double mass_diff_min = std::fabs(mass_diff[0]);
  /* Find a resonance whose pole mass is closest to
   * the input mass. */
  for (int ires = 1; ires < n_res; ires++) {
    if (std::fabs(mass_diff[ires]) < mass_diff_min) {
      ires_closest = ires;
      mass_diff_min = mass_diff[ires];
    }
  }
  logg[LPythia].debug("Quark constituents ", idq1, " and ", idq2, " with mass ",
                      mass, " (GeV) turned into a resonance ",
                      pdgid_possible[ires_closest]);
  return pdgid_possible[ires_closest];
}

bool StringProcess::make_lightcone_final_two(
    bool separate_fragment_hadron, double ppos_string, double pneg_string,
    double mTrn_had_forward, double mTrn_had_backward, double &ppos_had_forward,
    double &ppos_had_backward, double &pneg_had_forward,
    double &pneg_had_backward) {
  const double mTsqr_string = 2. * ppos_string * pneg_string;
  if (mTsqr_string < 0.) {
    return false;
  }
  const double mTrn_string = std::sqrt(mTsqr_string);
  if (mTrn_string < mTrn_had_forward + mTrn_had_backward) {
    return false;
  }

  // square of transvere mass of the forward hadron
  const double mTsqr_had_forward = mTrn_had_forward * mTrn_had_forward;
  // square of transvere mass of the backward hadron
  const double mTsqr_had_backward = mTrn_had_backward * mTrn_had_backward;

  /* The following part determines lightcone momentum fraction of p^+
   * carried by each hadron.
   * Lightcone momenta of the forward and backward hadrons are
   * p^+ forward  = (xe_pos + xpz_pos) * p^+ string,
   * p^- forward  = (xe_pos - xpz_pos) * p^- string,
   * p^+ backward = (xe_neg - xpz_pos) * p^+ string and
   * p^- backward = (xe_neg + xpz_pos) * p^- string.
   * where xe_pos and xe_neg satisfy xe_pos + xe_neg = 1.
   * Then evaluate xe_pos, xe_neg and xpz_pos in terms of
   * the transverse masses of hadrons and string. */

  // Express xe_pos and xe_neg in terms of the transverse masses.
  const double xm_diff =
      (mTsqr_had_forward - mTsqr_had_backward) / mTsqr_string;
  const double xe_pos = 0.5 * (1. + xm_diff);
  const double xe_neg = 0.5 * (1. - xm_diff);

  // Express xpz_pos in terms of the transverse masses.
  const double lambda_sqr =
      pow_int(mTsqr_string - mTsqr_had_forward - mTsqr_had_backward, 2) -
      4. * mTsqr_had_forward * mTsqr_had_backward;
  if (lambda_sqr <= 0.) {
    return false;
  }
  const double lambda = std::sqrt(lambda_sqr);
  const double b_lund =
      separate_fragment_hadron ? stringz_b_leading_ : stringz_b_produce_;
  /* The probability to flip sign of xpz_pos is taken from
   * StringFragmentation::finalTwo in StringFragmentation.cc
   * of PYTHIA 8. */
  const double prob_reverse =
      std::exp(-b_lund * lambda) / (1. + std::exp(-b_lund * lambda));
  double xpz_pos = 0.5 * lambda / mTsqr_string;
  if (random::uniform(0., 1.) < prob_reverse) {
    xpz_pos = -xpz_pos;
  }

  ppos_had_forward = (xe_pos + xpz_pos) * ppos_string;
  ppos_had_backward = (xe_neg - xpz_pos) * ppos_string;

  pneg_had_forward = 0.5 * mTsqr_had_forward / ppos_had_forward;
  pneg_had_backward = 0.5 * mTsqr_had_backward / ppos_had_backward;

  return true;
}

double StringProcess::sample_zLund(double a, double b, double mTrn) {
  // the lightcone momentum fraction x
  double xfrac = 0.;
  bool xfrac_accepted = false;
  /* First sample the inverse 1/x of the lightcone momentum fraction.
   * Then, obtain x itself.
   * The probability distribution function for the inverse of x is
   * PDF(u = 1/x) = (1/u) * (1 - 1/u)^a * exp(-b * mTrn^2 * u)
   * with 1 < u < infinity.
   * The rejection method can be used with an envelop function
   * ENV(u) = exp(-b * mTrn^2 * u). */
  while (!xfrac_accepted) {
    const double fac_env = b * mTrn * mTrn;
    const double u_aux = random::uniform(0., 1.);
    /* Sample u = 1/x according to the envelop function
     * ENV(u) = exp(-b * mTrn^2 * u). */
    const double xfrac_inv = 1. - std::log(u_aux) / fac_env;
    assert(xfrac_inv >= 1.);
    /* Evaluate the ratio of the real probability distribution function to
     * the envelop function. */
    const double xf_ratio = std::pow(1. - 1. / xfrac_inv, a) / xfrac_inv;
    // Determine whether the sampled value will be accepted.
    if (random::uniform(0., 1.) <= xf_ratio) {
      /* If the sampled value of 1/x is accepted,
       * obtain the value of x. */
      xfrac = 1. / xfrac_inv;
      xfrac_accepted = true;
    }
  }
  return xfrac;
}

bool StringProcess::remake_kinematics_fragments(
    Pythia8::Event &event_fragments, std::array<ThreeVector, 3> &evec_basis,
    double ppos_string, double pneg_string, double QTrx_string,
    double QTry_string, double QTrx_add_pos, double QTry_add_pos,
    double QTrx_add_neg, double QTry_add_neg) {
  logg[LPythia].debug("Correcting the kinematics of fragmented hadrons...");

  if (ppos_string < 0. || pneg_string < 0.) {
    logg[LPythia].debug(
        "  wrong lightcone momenta of string : ppos_string (GeV) = ",
        ppos_string, " pneg_string (GeV) = ", pneg_string);
    return false;
  }
  // Momentum rapidity of the final string
  const double yrapid_string = 0.5 * std::log(ppos_string / pneg_string);
  logg[LOutput].debug("Momentum-space rapidity of the string should be ",
                      yrapid_string);

  // Transverse mass of the final string
  const double mTrn_string = std::sqrt(2. * ppos_string * pneg_string);
  logg[LOutput].debug("Transvere mass (GeV) of the string should be ",
                      mTrn_string);
  // Transverse momentum of the final string
  const double QTrn_string =
      std::sqrt(QTrx_string * QTrx_string + QTry_string * QTry_string);
  if (mTrn_string < QTrn_string) {
    logg[LOutput].debug(
        "  wrong transverse mass of string : mTrn_string (GeV) = ", mTrn_string,
        " QTrn_string (GeV) = ", QTrn_string);
    return false;
  }
  const double msqr_string =
      mTrn_string * mTrn_string - QTrn_string * QTrn_string;
  // Mass of the final string
  const double mass_string = std::sqrt(msqr_string);
  logg[LOutput].debug("The string mass (GeV) should be ", mass_string);

  /* If there is no transverse momentum to be added to the string ends,
   * skip the entire procedure and return. */
  if (std::fabs(QTrx_add_pos) < small_number * mass_string &&
      std::fabs(QTry_add_pos) < small_number * mass_string &&
      std::fabs(QTrx_add_neg) < small_number * mass_string &&
      std::fabs(QTry_add_neg) < small_number * mass_string) {
    logg[LOutput].debug("  no need to add transverse momenta - skipping.");
    return true;
  }

  FourVector ptot_string_ini = FourVector(0., 0., 0., 0.);
  // Compute total four-momentum of the initial string.
  for (int ipyth = 1; ipyth < event_fragments.size(); ipyth++) {
    if (!event_fragments[ipyth].isFinal()) {
      continue;
    }

    FourVector p_frag =
        FourVector(event_fragments[ipyth].e(), event_fragments[ipyth].px(),
                   event_fragments[ipyth].py(), event_fragments[ipyth].pz());
    ptot_string_ini += p_frag;
  }
  const double E_string_ini = ptot_string_ini.x0();
  const double pz_string_ini = ptot_string_ini.threevec() * evec_basis[0];
  const double ppos_string_ini = (E_string_ini + pz_string_ini) * M_SQRT1_2;
  const double pneg_string_ini = (E_string_ini - pz_string_ini) * M_SQRT1_2;
  // Compute the momentum rapidity of the initial string.
  const double yrapid_string_ini =
      0.5 * std::log(ppos_string_ini / pneg_string_ini);
  /* Then, boost into the frame in which string is at rest in the
   * longitudinal direction. */
  shift_rapidity_event(event_fragments, evec_basis, 1., -yrapid_string_ini);

  int ip_forward = 0;
  int ip_backward = 0;
  double y_forward = 0.;
  double y_backward = 0.;
  ptot_string_ini = FourVector(0., 0., 0., 0.);
  // Find the most forward and backward hadrons based on the momentum rapidity.
  for (int ipyth = 1; ipyth < event_fragments.size(); ipyth++) {
    if (!event_fragments[ipyth].isFinal()) {
      continue;
    }

    FourVector p_frag =
        FourVector(event_fragments[ipyth].e(), event_fragments[ipyth].px(),
                   event_fragments[ipyth].py(), event_fragments[ipyth].pz());
    ptot_string_ini += p_frag;

    const double E_frag = p_frag.x0();
    const double pl_frag = p_frag.threevec() * evec_basis[0];
    double y_current = 0.5 * std::log((E_frag + pl_frag) / (E_frag - pl_frag));
    if (y_current > y_forward) {
      ip_forward = ipyth;
      y_forward = y_current;
    }
    if (y_current < y_backward) {
      ip_backward = ipyth;
      y_backward = y_current;
    }
  }
  logg[LOutput].debug("  The most forward hadron is ip_forward = ", ip_forward,
                      " with rapidity ", y_forward);
  logg[LOutput].debug("  The most backward hadron is ip_backward = ",
                      ip_backward, " with rapidity ", y_backward);

  const double px_string_ini = ptot_string_ini.threevec() * evec_basis[1];
  const double py_string_ini = ptot_string_ini.threevec() * evec_basis[2];

  /* Check if the transverse momentum px is conserved i.e.,
   * px of the initial string + px to be added = px of the final string */
  bool correct_px = std::fabs(px_string_ini + QTrx_add_pos + QTrx_add_neg -
                              QTrx_string) < small_number * mass_string;
  if (!correct_px) {
    logg[LOutput].debug(
        "  input transverse momenta in x-axis are not consistent.");
    return false;
  }
  /* Check if the transverse momentum py is conserved i.e.,
   * py of the initial string + py to be added = py of the final string */
  bool correct_py = std::fabs(py_string_ini + QTry_add_pos + QTry_add_neg -
                              QTry_string) < small_number * mass_string;
  if (!correct_py) {
    logg[LOutput].debug(
        "  input transverse momenta in y-axis are not consistent.");
    return false;
  }

  Pythia8::Vec4 pvec_string_now =
      set_Vec4(ptot_string_ini.x0(), ptot_string_ini.threevec());

  logg[LOutput].debug(
      "  Adding transverse momentum to the most forward hadron.");
  pvec_string_now -= event_fragments[ip_forward].p();
  const double mass_frag_pos = event_fragments[ip_forward].p().mCalc();
  // Four-momentum of the most forward hadron
  FourVector p_old_frag_pos = FourVector(
      event_fragments[ip_forward].e(), event_fragments[ip_forward].px(),
      event_fragments[ip_forward].py(), event_fragments[ip_forward].pz());
  // Add transverse momentum to it.
  ThreeVector mom_new_frag_pos = p_old_frag_pos.threevec() +
                                 QTrx_add_pos * evec_basis[1] +
                                 QTry_add_pos * evec_basis[2];
  // Re-calculate the energy.
  double E_new_frag_pos =
      std::sqrt(mom_new_frag_pos.sqr() + mass_frag_pos * mass_frag_pos);
  Pythia8::Vec4 pvec_new_frag_pos = set_Vec4(E_new_frag_pos, mom_new_frag_pos);
  pvec_string_now += pvec_new_frag_pos;
  // Update the event record.
  event_fragments[ip_forward].p(pvec_new_frag_pos);

  logg[LOutput].debug(
      "  Adding transverse momentum to the most backward hadron.");
  pvec_string_now -= event_fragments[ip_backward].p();
  const double mass_frag_neg = event_fragments[ip_backward].p().mCalc();
  // Four-momentum of the most backward hadron
  FourVector p_old_frag_neg = FourVector(
      event_fragments[ip_backward].e(), event_fragments[ip_backward].px(),
      event_fragments[ip_backward].py(), event_fragments[ip_backward].pz());
  // Add transverse momentum to it.
  ThreeVector mom_new_frag_neg = p_old_frag_neg.threevec() +
                                 QTrx_add_neg * evec_basis[1] +
                                 QTry_add_neg * evec_basis[2];
  // Re-calculate the energy.
  double E_new_frag_neg =
      std::sqrt(mom_new_frag_neg.sqr() + mass_frag_neg * mass_frag_neg);
  Pythia8::Vec4 pvec_new_frag_neg = set_Vec4(E_new_frag_neg, mom_new_frag_neg);
  pvec_string_now += pvec_new_frag_neg;
  // Update the event record.
  event_fragments[ip_backward].p(pvec_new_frag_neg);

  // Update the event record with total four-momentum of the current string.
  event_fragments[0].p(pvec_string_now);
  event_fragments[0].m(pvec_string_now.mCalc());

  // Sum of transverse masses of all fragmented hadrons.
  double mTrn_frag_all = 0.;
  for (int ipyth = 1; ipyth < event_fragments.size(); ipyth++) {
    if (!event_fragments[ipyth].isFinal()) {
      continue;
    }

    FourVector p_frag =
        FourVector(event_fragments[ipyth].e(), event_fragments[ipyth].px(),
                   event_fragments[ipyth].py(), event_fragments[ipyth].pz());
    ptot_string_ini += p_frag;

    const double E_frag = p_frag.x0();
    const double pl_frag = p_frag.threevec() * evec_basis[0];
    const double ppos_frag = (E_frag + pl_frag) * M_SQRT1_2;
    const double pneg_frag = (E_frag - pl_frag) * M_SQRT1_2;
    const double mTrn_frag = std::sqrt(2. * ppos_frag * pneg_frag);
    mTrn_frag_all += mTrn_frag;
  }
  logg[LOutput].debug(
      "Sum of transverse masses (GeV) of all fragmented hadrons : ",
      mTrn_frag_all);
  /* If the transverse mass of the (final) string is smaller than
   * the sum of transverse masses, kinematics cannot be determined. */
  if (mTrn_string < mTrn_frag_all) {
    logg[LOutput].debug("  which is larger than mT of the actual string ",
                        mTrn_string);
    return false;
  }

  double mass_string_now = pvec_string_now.mCalc();
  double msqr_string_now = mass_string_now * mass_string_now;
  // Total four-momentum of the current string
  FourVector p_string_now =
      FourVector(pvec_string_now.e(), pvec_string_now.px(),
                 pvec_string_now.py(), pvec_string_now.pz());
  double E_string_now = p_string_now.x0();
  double pz_string_now = p_string_now.threevec() * evec_basis[0];
  logg[LOutput].debug("The string mass (GeV) at this point : ",
                      mass_string_now);
  double ppos_string_now = (E_string_now + pz_string_now) * M_SQRT1_2;
  double pneg_string_now = (E_string_now - pz_string_now) * M_SQRT1_2;
  // Momentum rapidity of the current string
  double yrapid_string_now = 0.5 * std::log(ppos_string_now / pneg_string_now);
  logg[LOutput].debug("The momentum-space rapidity of string at this point : ",
                      yrapid_string_now);
  logg[LOutput].debug(
      "The momentum-space rapidities of hadrons will be changed.");
  const int niter_max = 10000;
  bool accepted = false;
  double fac_all_yrapid = 1.;
  /* Rescale momentum rapidities of hadrons by replacing
   * y_hadron with y_string_now + fac_yrapid * (y_hadron - y_string_now).
   * This is done iteratively by finding the value of fac_yrapid which gives
   * the correct string mass. */
  for (int iiter = 0; iiter < niter_max; iiter++) {
    if (std::fabs(mass_string_now - mass_string) < really_small * mass_string) {
      accepted = true;
      break;
    }
    double E_deriv = 0.;
    double pz_deriv = 0.;

    /* Have a Taylor series of mass square as a linear function
     * of fac_yrapid and find a trial value of fac_yrapid. */
    for (int ipyth = 1; ipyth < event_fragments.size(); ipyth++) {
      if (!event_fragments[ipyth].isFinal()) {
        continue;
      }

      FourVector p_frag =
          FourVector(event_fragments[ipyth].e(), event_fragments[ipyth].px(),
                     event_fragments[ipyth].py(), event_fragments[ipyth].pz());
      const double E_frag = p_frag.x0();
      const double pl_frag = p_frag.threevec() * evec_basis[0];
      const double ppos_frag = (E_frag + pl_frag) * M_SQRT1_2;
      const double pneg_frag = (E_frag - pl_frag) * M_SQRT1_2;
      const double mTrn_frag = std::sqrt(2. * ppos_frag * pneg_frag);
      const double y_frag = 0.5 * std::log(ppos_frag / pneg_frag);

      E_deriv += mTrn_frag * (y_frag - yrapid_string_now) * std::sinh(y_frag);
      pz_deriv += mTrn_frag * (y_frag - yrapid_string_now) * std::cosh(y_frag);
    }
    double slope = 2. * (E_string_now * E_deriv - pz_string_now * pz_deriv);
    double fac_yrapid = 1. + std::tanh((msqr_string - msqr_string_now) / slope);
    fac_all_yrapid *= fac_yrapid;

    // Replace momentum rapidities of hadrons.
    shift_rapidity_event(event_fragments, evec_basis, fac_yrapid,
                         (1. - fac_yrapid) * yrapid_string_now);
    // Update the four-momentum and mass of the current string
    pvec_string_now = event_fragments[0].p();
    mass_string_now = pvec_string_now.mCalc();
    msqr_string_now = mass_string_now * mass_string_now;
    p_string_now = FourVector(pvec_string_now.e(), pvec_string_now.px(),
                              pvec_string_now.py(), pvec_string_now.pz());
    E_string_now = p_string_now.x0();
    pz_string_now = p_string_now.threevec() * evec_basis[0];
    ppos_string_now = (E_string_now + pz_string_now) * M_SQRT1_2;
    pneg_string_now = (E_string_now - pz_string_now) * M_SQRT1_2;
    yrapid_string_now = 0.5 * std::log(ppos_string_now / pneg_string_now);
    logg[LOutput].debug("  step ", iiter + 1, " : fac_yrapid = ", fac_yrapid,
                        " , string mass (GeV) = ", mass_string_now,
                        " , string rapidity = ", yrapid_string_now);
  }

  if (!accepted) {
    logg[LOutput].debug("  Too many iterations in rapidity rescaling.");
    return false;
  }
  logg[LOutput].debug(
      "The overall factor multiplied to the rapidities of hadrons = ",
      fac_all_yrapid);
  logg[LOutput].debug("The momentum-space rapidity of string at this point : ",
                      yrapid_string_now);
  const double y_diff = yrapid_string - yrapid_string_now;
  logg[LOutput].debug("The hadrons will be boosted by rapidity ", y_diff,
                      " for the longitudinal momentum conservation.");

  // Boost the hadrons back into the original frame.
  shift_rapidity_event(event_fragments, evec_basis, 1., y_diff);

  return true;
}

void StringProcess::assign_scaling_factor(int nquark, ParticleData &data,
                                          double suppression_factor) {
  int nbaryon = data.pdgcode().baryon_number();
  if (nbaryon == 0) {
    // Mesons always get a scaling factor of 1/2 since there is never
    // a q-qbar pair at the end of a string so nquark is always 1
    data.set_cross_section_scaling_factor(0.5 * suppression_factor);
  } else if (data.is_baryon()) {
    // Leading baryons get a factor of 2/3 if they carry 2
    // and 1/3 if they carry 1 of the strings valence quarks
    data.set_cross_section_scaling_factor(suppression_factor * nquark /
                                          (3.0 * nbaryon));
  }
}

std::pair<int, int> StringProcess::find_leading(int nq1, int nq2,
                                                ParticleList &list) {
  assert(list.size() >= 2);
  int end = list.size() - 1;
  int i1, i2;
  for (i1 = 0;
       i1 <= end && !list[i1].pdgcode().contains_enough_valence_quarks(nq1);
       i1++) {
  }
  for (i2 = end;
       i2 >= 0 && !list[i2].pdgcode().contains_enough_valence_quarks(nq2);
       i2--) {
  }
  std::pair<int, int> indices(i1, i2);
  return indices;
}

void StringProcess::assign_all_scaling_factors(int baryon_string,
                                               ParticleList &outgoing_particles,
                                               const ThreeVector &evecLong,
                                               double suppression_factor) {
  // Set each particle's cross section scaling factor to 0 first
  for (ParticleData &data : outgoing_particles) {
    data.set_cross_section_scaling_factor(0.0);
  }
  // sort outgoing particles according to the longitudinal velocity
  std::sort(outgoing_particles.begin(), outgoing_particles.end(),
            [&](ParticleData i, ParticleData j) {
              return i.momentum().velocity() * evecLong >
                     j.momentum().velocity() * evecLong;
            });
  int nq1, nq2;  // number of quarks at both ends of the string
  switch (baryon_string) {
    case 0:
      nq1 = -1;
      nq2 = 1;
      break;
    case 1:
      nq1 = 2;
      nq2 = 1;
      break;
    case -1:
      nq1 = -2;
      nq2 = -1;
      break;
    default:
      throw std::runtime_error("string is neither mesonic nor baryonic");
  }
  // Try to find nq1 on one string end and nq2 on the other string end and the
  // other way around. When the leading particles are close to the string ends,
  // the quarks are assumed to be distributed this way.
  std::pair<int, int> i = find_leading(nq1, nq2, outgoing_particles);
  std::pair<int, int> j = find_leading(nq2, nq1, outgoing_particles);
  if (baryon_string == 0 && i.second - i.first < j.second - j.first) {
    assign_scaling_factor(nq2, outgoing_particles[j.first], suppression_factor);
    assign_scaling_factor(nq1, outgoing_particles[j.second],
                          suppression_factor);
  } else {
    assign_scaling_factor(nq1, outgoing_particles[i.first], suppression_factor);
    assign_scaling_factor(nq2, outgoing_particles[i.second],
                          suppression_factor);
  }
}

int StringProcess::pdg_map_for_pythia(PdgCode &pdg) {
  PdgCode pdg_mapped(0x0);

  if (pdg.baryon_number() == 1) {  // baryon
    pdg_mapped = pdg.charge() > 0 ? PdgCode(pdg::p) : PdgCode(pdg::n);
  } else if (pdg.baryon_number() == -1) {  // antibaryon
    pdg_mapped = pdg.charge() < 0 ? PdgCode(-pdg::p) : PdgCode(-pdg::n);
  } else if (pdg.is_hadron()) {  // meson
    if (pdg.charge() >= 0) {
      pdg_mapped = PdgCode(pdg::pi_p);
    } else {
      pdg_mapped = PdgCode(pdg::pi_m);
    }
  } else if (pdg.is_lepton()) {  // lepton
    pdg_mapped = pdg.charge() < 0 ? PdgCode(0x11) : PdgCode(-0x11);
  } else {
    throw std::runtime_error("StringProcess::pdg_map_for_pythia failed.");
  }

  return pdg_mapped.get_decimal();
}

}  // namespace smash
