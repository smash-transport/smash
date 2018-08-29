/*
 *
 *    Copyright (c) 2017-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <array>

#include "smash/angles.h"
#include "smash/kinematics.h"
#include "smash/logging.h"
#include "smash/processstring.h"
#include "smash/random.h"

namespace smash {

StringProcess::StringProcess(double string_tension, double time_formation,
                             double gluon_beta, double gluon_pmin,
                             double quark_alpha, double quark_beta,
                             double strange_supp, double diquark_supp,
                             double sigma_perp, double leading_frag_mean,
                             double leading_frag_width, double stringz_a,
                             double stringz_b, double string_sigma_T,
                             double factor_t_form, bool use_yoyo_model)
    : pmin_gluon_lightcone_(gluon_pmin),
      pow_fgluon_beta_(gluon_beta),
      pow_fquark_alpha_(quark_alpha),
      pow_fquark_beta_(quark_beta),
      sigma_qperp_(sigma_perp),
      leading_frag_mean_(leading_frag_mean),
      leading_frag_width_(leading_frag_width),
      kappa_tension_string_(string_tension),
      additional_xsec_supp_(0.7),
      time_formation_const_(time_formation),
      soft_t_form_(factor_t_form),
      time_collision_(0.),
      use_yoyo_model_(use_yoyo_model) {
  // setup and initialize pythia for hard string process
  pythia_parton_ = make_unique<Pythia8::Pythia>(PYTHIA_XML_DIR, false);
  /* select only non-diffractive events
   * diffractive ones are implemented in a separate routine */
  pythia_parton_->readString("SoftQCD:nonDiffractive = on");
  pythia_parton_->readString("MultipartonInteractions:pTmin = 1.5");
  pythia_parton_->readString("HadronLevel:all = off");
  common_setup_pythia(pythia_parton_.get(), strange_supp, diquark_supp,
                      stringz_a, stringz_b, string_sigma_T);

  // setup and initialize pythia for fragmentation
  pythia_hadron_ = make_unique<Pythia8::Pythia>(PYTHIA_XML_DIR, false);
  /* turn off all parton-level processes to implement only hadronization */
  pythia_hadron_->readString("ProcessLevel:all = off");
  common_setup_pythia(pythia_hadron_.get(), strange_supp, diquark_supp,
                      stringz_a, stringz_b, string_sigma_T);

  /* initialize PYTHIA */
  pythia_hadron_->init();
  pythia_sigmatot_.init(&pythia_hadron_->info, pythia_hadron_->settings,
                        &pythia_hadron_->particleData);
  event_intermediate_.init("intermediate partons",
                           &pythia_hadron_->particleData);

  sqrt2_ = std::sqrt(2.);

  for (int imu = 0; imu < 3; imu++) {
    evecBasisAB_[imu] = ThreeVector(0., 0., 0.);
  }

  final_state_.clear();
}

void StringProcess::common_setup_pythia(Pythia8::Pythia *pythia_in,
                                        double strange_supp,
                                        double diquark_supp, double stringz_a,
                                        double stringz_b,
                                        double string_sigma_T) {
  // choose parametrization for mass-dependent width
  pythia_in->readString("ParticleData:modeBreitWigner = 4");
  /* choose minimum transverse momentum scale
   * involved in partonic interactions */
  pythia_in->readString("MultipartonInteractions:pTmin = 1.5");
  pythia_in->readString("MultipartonInteractions:nSample = 10000");
  // transverse momentum spread in string fragmentation
  pythia_in->readString("StringPT:sigma = " + std::to_string(string_sigma_T));
  // diquark suppression factor in string fragmentation
  pythia_in->readString("StringFlav:probQQtoQ = " +
                        std::to_string(diquark_supp));
  // strangeness suppression factor in string fragmentation
  pythia_in->readString("StringFlav:probStoUD = " +
                        std::to_string(strange_supp));
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
}

// compute the formation time and fill the arrays with final-state particles
int StringProcess::append_final_state(ParticleList &intermediate_particles,
                                      const FourVector &uString,
                                      const ThreeVector &evecLong) {
  int nfrag = 0;
  int bstring = 0;
  double p_pos_tot = 0.0, p_neg_tot = 0.0;

  for (ParticleData &data : intermediate_particles) {
    nfrag += 1;
    bstring += data.pdgcode().baryon_number();

    const double pparallel = data.momentum().threevec() * evecLong;
    p_pos_tot += (data.momentum().x0() + pparallel) / sqrt2_;
    p_neg_tot += (data.momentum().x0() - pparallel) / sqrt2_;
  }
  assert(nfrag > 0);

  /* compute the cross section scaling factor for leading hadrons
   * based on the number of valence quarks. */
  assign_all_scaling_factors(bstring, intermediate_particles,
                             evecLong, additional_xsec_supp_);

  std::vector<double> xvertex_pos, xvertex_neg;
  xvertex_pos.resize(nfrag + 1);
  xvertex_neg.resize(nfrag + 1);
  // x^{+} coordinates of the forward end
  xvertex_pos[0] = p_pos_tot / kappa_tension_string_;
  // x^{-} coordinates of the backward end
  xvertex_neg[nfrag] = p_neg_tot / kappa_tension_string_;

  for (int i = 0; i < nfrag; i++) {
    const double pparallel =
        intermediate_particles[i].momentum().threevec() * evecLong;

    // recursively compute x^{+} coordinates of q-qbar formation vertex
    xvertex_pos[i + 1] = xvertex_pos[i] -
                         (intermediate_particles[i].momentum().x0() +
                          pparallel) /
                         (kappa_tension_string_ * sqrt2_);

    // recursively compute x^{-} coordinates of q-qbar formation vertex
    const int ineg = nfrag - i - 1;
    xvertex_neg[ineg] = xvertex_neg[ineg + 1] -
                        (intermediate_particles[ineg].momentum().x0() -
                         pparallel) /
                        (kappa_tension_string_ * sqrt2_);
  }

  // Velocity three-vector to perform Lorentz boost.
  const ThreeVector vstring = uString.velocity();

  /* compute the formation times of hadrons
   * from the lightcone coordinates of q-qbar formation vertices. */
  for (int i = 0; i < nfrag; i++) {
    // boost 4-momentum into the center of mass frame
    FourVector momentum =
        intermediate_particles[i].momentum().LorentzBoost(-vstring);
    intermediate_particles[i].set_4momentum(momentum);

    if (use_yoyo_model_) {
      // set the formation time and position in the rest frame of string
      double t_prod = (xvertex_pos[i] + xvertex_neg[i + 1]) / sqrt2_;
      FourVector fragment_position = FourVector(
          t_prod, evecLong * (xvertex_pos[i] - xvertex_neg[i + 1]) / sqrt2_);
      /* boost formation vertex into the center of mass frame
       * and then into the lab frame */
      fragment_position = fragment_position.LorentzBoost(-vstring);
      fragment_position = fragment_position.LorentzBoost(-vcomAB_);
      intermediate_particles[i].set_formation_time(
          soft_t_form_ * fragment_position.x0() + time_collision_);
    } else {
      ThreeVector v_calc =
          momentum.LorentzBoost(-vcomAB_).velocity();
      double gamma_factor = 1.0 / std::sqrt(1 - (v_calc).sqr());
      intermediate_particles[i].set_slow_formation_times(time_collision_,
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

  pcom_[0] = plab_[0].LorentzBoost(vcomAB_);
  pcom_[1] = plab_[1].LorentzBoost(vcomAB_);

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
  make_string_ends(is_AB_to_AX ? PDGcodes_[1] : PDGcodes_[0], idqX1, idqX2);
  // string mass must be larger than threshold set by PYTHIA.
  mstrMin = pythia_hadron_->particleData.m0(idqX1) +
            pythia_hadron_->particleData.m0(idqX2);
  // this threshold cannot be larger than maximum of allowed string mass.
  if (mstrMin > mstrMax) {
    return false;
  }
  // sample the transverse momentum transfer
  QTrx = random::normal(0., sigma_qperp_ / sqrt2_);
  QTry = random::normal(0., sigma_qperp_ / sqrt2_);
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
  const FourVector prs = pnull.LorentzBoost(ustrXcom.velocity());
  ThreeVector evec = prs.threevec() / prs.threevec().abs();
  // perform fragmentation and add particles to final_state.
  ParticleList new_intermediate_particles;
  bool separate_fragment_baryon = false;
  int nfrag = fragment_string(idqX1, idqX2, massX, evec, true,
                              separate_fragment_baryon,
                              new_intermediate_particles);
  if (nfrag < 1) {
    NpartString_[0] = 0;
    return false;
  }
  NpartString_[0] = append_final_state(new_intermediate_particles,
                                       ustrXcom, evec);

  NpartString_[1] = 1;
  PdgCode hadron_code = is_AB_to_AX ? PDGcodes_[0] : PDGcodes_[1];
  ParticleData new_particle(ParticleType::find(hadron_code));
  new_particle.set_4momentum(pstrHcom);
  new_particle.set_cross_section_scaling_factor(1.);
  new_particle.set_formation_time(time_collision_);
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
      const FourVector prs = pnull.LorentzBoost(pstr_com[i].velocity());
      evec_str[i] = prs.threevec() / prs.threevec().abs();
    }
  }

  return found_mass[0] && found_mass[1];
}

bool StringProcess::make_final_state_2strings(
    const std::array<std::array<int, 2>, 2> &quarks,
    const std::array<FourVector, 2> &pstr_com,
    const std::array<double, 2> &m_str,
    const std::array<ThreeVector, 2> &evec_str,
    bool flip_string_ends, bool separate_fragment_baryon) {
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
    NpartString_[i] = append_final_state(new_intermediate_particles,
                                         ustr_com[i], evec);
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
  make_string_ends(PDGcodes_[0], quarks[0][0], quarks[0][1]);
  make_string_ends(PDGcodes_[1], quarks[1][0], quarks[1][1]);
  // sample the lightcone momentum fraction carried by gluons
  const double xmin_gluon_fraction = pmin_gluon_lightcone_ / sqrtsAB_;
  const double xfracA =
      random::beta_a0(xmin_gluon_fraction, pow_fgluon_beta_ + 1.);
  const double xfracB =
      random::beta_a0(xmin_gluon_fraction, pow_fgluon_beta_ + 1.);
  // sample the transverse momentum transfer
  const double QTrx = random::normal(0., sigma_qperp_ / sqrt2_);
  const double QTry = random::normal(0., sigma_qperp_ / sqrt2_);
  const double QTrn = std::sqrt(QTrx * QTrx + QTry * QTry);
  // evaluate the lightcone momentum transfer
  const double QPos = -QTrn * QTrn / (2. * xfracB * PNegB_);
  const double QNeg = QTrn * QTrn / (2. * xfracA * PPosA_);
  // compute four-momentum of string 1
  threeMomentum = evecBasisAB_[0] * (PPosA_ + QPos - PNegA_ - QNeg) / sqrt2_ +
                  evecBasisAB_[1] * QTrx + evecBasisAB_[2] * QTry;
  pstr_com[0] =
      FourVector((PPosA_ + QPos + PNegA_ + QNeg) / sqrt2_, threeMomentum);
  // compute four-momentum of string 2
  threeMomentum = evecBasisAB_[0] * (PPosB_ - QPos - PNegB_ + QNeg) / sqrt2_ -
                  evecBasisAB_[1] * QTrx - evecBasisAB_[2] * QTry;
  pstr_com[1] =
      FourVector((PPosB_ - QPos + PNegB_ - QNeg) / sqrt2_, threeMomentum);

  const bool found_masses =
      set_mass_and_direction_2strings(quarks, pstr_com, m_str, evec_str);
  if (!found_masses) {
    return false;
  }
  const bool flip_string_ends = true;
  const bool separate_fragment_baryon = false;
  const bool success = make_final_state_2strings(quarks, pstr_com, m_str,
                                                 evec_str, flip_string_ends,
                                                 separate_fragment_baryon);
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
  make_string_ends(PDGcodes_[0], idqA1, idqA2);
  make_string_ends(PDGcodes_[1], idqB1, idqB2);

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
  const double QTrx = random::normal(0., sigma_qperp_ / sqrt2_);
  const double QTry = random::normal(0., sigma_qperp_ / sqrt2_);
  const double QTrn = std::sqrt(QTrx * QTrx + QTry * QTry);
  // evaluate the lightcone momentum transfer
  const double QPos = -QTrn * QTrn / (2. * xfracB * PNegB_);
  const double QNeg = QTrn * QTrn / (2. * xfracA * PPosA_);
  const double dPPos = -xfracA * PPosA_ - QPos;
  const double dPNeg = xfracB * PNegB_ - QNeg;
  // compute four-momentum of string 1
  ThreeVector threeMomentum =
      evecBasisAB_[0] * (PPosA_ + dPPos - PNegA_ - dPNeg) / sqrt2_ +
      evecBasisAB_[1] * QTrx + evecBasisAB_[2] * QTry;
  pstr_com[0] =
      FourVector((PPosA_ + dPPos + PNegA_ + dPNeg) / sqrt2_, threeMomentum);
  m_str[0] = pstr_com[0].sqr();
  // compute four-momentum of string 2
  threeMomentum = evecBasisAB_[0] * (PPosB_ - dPPos - PNegB_ + dPNeg) / sqrt2_ -
                  evecBasisAB_[1] * QTrx - evecBasisAB_[2] * QTry;
  pstr_com[1] =
      FourVector((PPosB_ - dPPos + PNegB_ - dPNeg) / sqrt2_, threeMomentum);

  const bool found_masses =
      set_mass_and_direction_2strings(quarks, pstr_com, m_str, evec_str);
  if (!found_masses) {
    return false;
  }
  const bool flip_string_ends = false;
  const bool separate_fragment_baryon = true;
  const bool success = make_final_state_2strings(quarks, pstr_com, m_str,
                                                 evec_str, flip_string_ends,
                                                 separate_fragment_baryon);
  return success;
}

// hard non-diffractive
bool StringProcess::next_NDiffHard() {
  const auto &log = logger<LogArea::Pythia>();
  NpartFinal_ = 0;
  final_state_.clear();

  log.debug("Hard non-diff. with ", PDGcodes_[0], " + ",
            PDGcodes_[1], " at CM energy [GeV] ", sqrtsAB_);

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
    log.debug("  incoming particle ", i, " : ",
              PDGcodes_[i], " is mapped onto ", pdg_for_pythia[i]);

    PdgCode pdgcode_for_pythia(std::to_string(pdg_for_pythia[i]));
    /* evaluate how many more constituents incoming hadron has
     * compared to the mapped one. */
    find_excess_constituent(PDGcodes_[i], pdgcode_for_pythia,
                            excess_quark[i], excess_antiq[i]);
  }

  int previous_idA = pythia_parton_->mode("Beams:idA"),
      previous_idB = pythia_parton_->mode("Beams:idB");
  double previous_eCM = pythia_parton_->parm("Beams:eCM");
  // check if the initial state for PYTHIA remains same.
  bool same_initial_state = previous_idA == pdg_for_pythia[0] &&
                            previous_idB == pdg_for_pythia[1] &&
                            std::abs(previous_eCM - sqrtsAB_) < really_small;

  /* Perform PYTHIA initialization if it was not previously initialized
   * or the initial state changed. */
  if (!pythia_parton_initialized_ || !same_initial_state) {
    pythia_parton_->settings.mode("Beams:idA", pdg_for_pythia[0]);
    pythia_parton_->settings.mode("Beams:idB", pdg_for_pythia[1]);
    pythia_parton_->settings.parm("Beams:eCM", sqrtsAB_);

    pythia_parton_initialized_ = pythia_parton_->init();
    log.info("Pythia initialized with ", pdg_for_pythia[0], " + ",
             pdg_for_pythia[1], " at CM energy [GeV] ", sqrtsAB_);
    if (!pythia_parton_initialized_) {
      throw std::runtime_error("Pythia failed to initialize.");
    }
  }
  /* Set the random seed of the Pythia random Number Generator.
   * Pythia's random is controlled by SMASH in every single collision.
   * In this way we ensure that the results are reproducible
   * for every event if one knows SMASH random seed. */
  pythia_parton_->rndm.init(random::uniform_int(1,
                                std::numeric_limits<int>::max()));

  // Short notation for Pythia event
  Pythia8::Event &event_hadron = pythia_hadron_->event;
  log.debug("Pythia hard event created");
  bool final_state_success = pythia_parton_->next();
  log.debug("Pythia final state computed, success = ", final_state_success);
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
  for (int i = 0; i < pythia_parton_->event.size(); i++) {
    if (pythia_parton_->event[i].isFinal()) {
      const int pdgid = pythia_parton_->event[i].id();
      Pythia8::Vec4 pquark = pythia_parton_->event[i].p();
      const double mass = pythia_parton_->particleData.m0(pdgid);

      const int status = pythia_parton_->event[i].status();
      const int color = pythia_parton_->event[i].col();
      const int anticolor = pythia_parton_->event[i].acol();

      pSum += pquark;
      event_intermediate_.append(pdgid, status,
                                 color, anticolor, pquark, mass);
    }
  }
  if (!final_state_success) {
    return false;
  }
  // add junctions to the intermediate state if there is any.
  event_intermediate_.clearJunctions();
  for (int i = 0; i < pythia_parton_->event.sizeJunction(); i++) {
    const int kind = pythia_parton_->event.kindJunction(i);
    std::array<int, 3> col;
    for (int j = 0; j < 3; j++) {
      col[j] = pythia_parton_->event.colJunction(i, j);
    }
    event_intermediate_.appendJunction(kind, col[0], col[1], col[2]);
  }
  event_intermediate_[0].p(pSum);
  event_intermediate_[0].m(pSum.mCalc());

  /* Replace quark constituents according to the excess of valence quarks
   * and then rescale momenta of partons by constant factor
   * to fulfill the energy-momentum conservation. */
  restore_constituent(event_intermediate_, excess_quark, excess_antiq);

  int npart = event_intermediate_.size();
  int ipart = 0;
  while (ipart < npart) {
    const int pdgid = event_intermediate_[ipart].id();
    if (event_intermediate_[ipart].isFinal() &&
        !event_intermediate_[ipart].isParton() &&
        !pythia_parton_->particleData.isOctetHadron(pdgid)) {
      FourVector momentum = reorient(event_intermediate_[ipart], evecBasisAB_);
      log.debug("4-momentum from Pythia: ", momentum);
      bool found_ptype = append_intermediate_list(pdgid, momentum,
                                                  new_non_hadron_particles);
      if (!found_ptype) {
        log.warn("PDG ID ", pdgid,
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
  log.debug("Hard non-diff: partonic process gives ",
            event_intermediate_.size(), " partons.");
  // identify and fragment strings until there is no parton left.
  while (event_intermediate_.size() > 1) {
    /* identify string from a most forward or backward parton.
     * if there is no junction. */
    compose_string_parton(find_forward_string,
                          event_intermediate_, pythia_hadron_->event);
    // identify string from a junction if there is any.
    compose_string_junction(find_forward_string,
                            event_intermediate_, pythia_hadron_->event);

    pythia_hadron_->rndm.init(random::uniform_int(1,
                                  std::numeric_limits<int>::max()));
    // fragment the (identified) string into hadrons.
    hadronize_success = pythia_hadron_->next();
    log.debug("Pythia hadronized, success = ", hadronize_success);

    new_intermediate_particles.clear();
    if (hadronize_success) {
      for (int i = 0; i < event_hadron.size(); i++) {
        if (event_hadron[i].isFinal()) {
          int pythia_id = event_hadron[i].id();
          log.debug("PDG ID from Pythia:", pythia_id);
          /* K_short and K_long need to be converted to K0
           * since SMASH only knows K0 */
          convert_KaonLS(pythia_id);

          /* evecBasisAB_[0] is a unit 3-vector in the collision axis,
           * while evecBasisAB_[1] and evecBasisAB_[2] spans the transverse plane.
           * Given that PYTHIA assumes z-direction to be the collision axis,
           * pz from PYTHIA should be the momentum compoment in evecBasisAB_[0].
           * px and py are respectively the momentum components in two
           * transverse directions evecBasisAB_[1] and evecBasisAB_[2]. */
          FourVector momentum = reorient(event_hadron[i], evecBasisAB_);
          log.debug("4-momentum from Pythia: ", momentum);
          log.debug("appending the particle ", pythia_id,
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
            log.warn("PDG ID ", pythia_id,
                     " does not exist in ParticleType - start over.");
            hadronize_success = false;
          }
        }
      }
    }

    /* if hadronization is not successful,
     * reset the event records, return false and then start over. */
    if (!hadronize_success) {
      event_intermediate_.reset();
      pythia_hadron_->event.reset();
      new_intermediate_particles.clear();
      break;
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

void StringProcess::replace_constituent(Pythia8::Particle &particle,
                        std::array<int, 5> &excess_constituent) {
  const auto &log = logger<LogArea::Pythia>();

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
      const int l = random::uniform_int(0,
                        static_cast<int>(k_found.size()) - 1);
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
  log.debug("  parton id = ", particle.id(), " is converted to ", pdgid_new);

  // update the constituent mass and energy.
  Pythia8::Vec4 pquark = particle.p();
  double mass_new = pythia_hadron_->particleData.m0(pdgid_new);
  double e_new = std::sqrt(mass_new * mass_new + pquark.pAbs() * pquark.pAbs());
  // update the particle object.
  particle.id(pdgid_new);
  particle.e(e_new);
  particle.m(mass_new);
}

void StringProcess::restore_constituent(Pythia8::Event &event_intermediate,
                        std::array<std::array<int, 5>, 2> &excess_quark,
                        std::array<std::array<int, 5>, 2> &excess_antiq) {
  const auto &log = logger<LogArea::Pythia>();

  Pythia8::Vec4 pSum = event_intermediate[0].p();
  const double energy_init = pSum.e();
  log.debug("  initial total energy [GeV] : ", energy_init);

  bool recovered_quarks = false;
  while (!recovered_quarks) {
    /* Flavor conversion begins with the most forward and backward parton
     * respectively for incoming_particles_[0] and incoming_particles_[1]. */
    std::array<bool, 2> find_forward = {true, false};
    const std::array<int, 5> excess_null = {0, 0, 0, 0, 0};
    std::array<int, 5> excess_total = excess_null;

    for (int ih = 0; ih < 2; ih++) {  // loop over incoming hadrons
      int nfrag = event_intermediate.size();
      for (int np_end = 0; np_end < nfrag - 2; np_end++) {  // constituent loop
        int iforward = 1;
        /* select the np_end-th most forward or backward parton and
         * change its specie.
         * np_end = 0 corresponds to the most forward,
         * np_end = 1 corresponds to the second most forward and so on. */
        for (int ip = 2; ip < event_intermediate.size() - np_end; ip++) {
          const double y_quark_current = event_intermediate[ip].y();
          const double y_quark_forward = event_intermediate[iforward].y();
          if ((find_forward[ih] &&
              y_quark_current > y_quark_forward) ||
              (!find_forward[ih] &&
              y_quark_current < y_quark_forward)) {
            iforward = ip;
          }
        }
        pSum -= event_intermediate[iforward].p();

        if (event_intermediate[iforward].id() > 0) {  // quark and diquark
          replace_constituent(event_intermediate[iforward], excess_quark[ih]);
        } else {  // antiquark and anti-diquark
          replace_constituent(event_intermediate[iforward], excess_antiq[ih]);
        }

        const int pdgid = event_intermediate[iforward].id();
        Pythia8::Vec4 pquark = event_intermediate[iforward].p();
        const double mass = pythia_hadron_->particleData.m0(pdgid);

        const int status = event_intermediate[iforward].status();
        const int color = event_intermediate[iforward].col();
        const int anticolor = event_intermediate[iforward].acol();

        pSum += pquark;
        event_intermediate.append(pdgid, status,
                                  color, anticolor, pquark, mass);

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
  log.debug("  valence quark contents of hadons are recovered.");

  log.debug("  current total energy [GeV] : ", pSum.e());
  /* rescale momenta of all partons by a constant factor
   * to conserve the total energy. */
  while (true) {
    if (std::abs(pSum.e() - energy_init) < really_small * energy_init) {
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
    log.debug("  parton momenta are rescaled by factor of ", rescale_factor);
  }

  log.debug("  final total energy [GeV] : ", pSum.e());
  event_intermediate[0].p(pSum);
  event_intermediate[0].m(pSum.mCalc());
}

void StringProcess::compose_string_parton(bool find_forward_string,
                                          Pythia8::Event &event_intermediate,
                                          Pythia8::Event &event_hadronize) {
  const auto &log = logger<LogArea::Pythia>();

  // if there is a junction, process it first.
  if (event_intermediate.sizeJunction() > 0) {
    return;
  }

  Pythia8::Vec4 pSum = 0.;
  event_hadronize.reset();

  int iforward = 1;
  // select the most forward or backward parton.
  for (int i = 2; i < event_intermediate.size(); i++) {
    if ((find_forward_string &&
         event_intermediate[i].y() > event_intermediate[iforward].y()) ||
        (!find_forward_string &&
         event_intermediate[i].y() < event_intermediate[iforward].y())) {
      iforward = i;
    }
  }
  log.debug("Hard non-diff: iforward = ", iforward,
            "(", event_intermediate[iforward].id(), ")");

  pSum += event_intermediate[iforward].p();
  event_hadronize.append(event_intermediate[iforward]);

  int col_to_find = event_intermediate[iforward].acol();
  int acol_to_find = event_intermediate[iforward].col();
  event_intermediate.remove(iforward, iforward);
  log.debug("Hard non-diff: event_intermediate reduces in size to ",
            event_intermediate.size());

  // trace color and anti-color indices and find corresponding partons.
  while (col_to_find != 0 || acol_to_find != 0) {
    log.debug("  col_to_find = ", col_to_find,
              ", acol_to_find = ", acol_to_find);

    int ifound = -1;
    for (int i = 1; i < event_intermediate.size(); i++) {
      const int pdgid = event_intermediate[i].id();
      bool found_col = col_to_find != 0 &&
                       col_to_find == event_intermediate[i].col();
      bool found_acol = acol_to_find != 0 &&
                        acol_to_find == event_intermediate[i].acol();
      if (found_col) {
        log.debug("  col_to_find ", col_to_find,
                  " from i ", i, "(", pdgid, ") found");
      }
      if (found_acol) {
        log.debug("  acol_to_find ", acol_to_find,
                  " from i ", i, "(", pdgid, ") found");
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
      throw std::runtime_error("Hard string could not be identified.");
    } else {
      pSum += event_intermediate[ifound].p();
      // add a parton to the new event record.
      event_hadronize.append(event_intermediate[ifound]);
      // then remove from the original event record.
      event_intermediate.remove(ifound, ifound);
      log.debug("Hard non-diff: event_intermediate reduces in size to ",
                event_intermediate.size());
    }
  }

  event_hadronize[0].p(pSum);
  event_hadronize[0].m(pSum.mCalc());
}

void StringProcess::compose_string_junction(bool &find_forward_string,
                                            Pythia8::Event &event_intermediate,
                                            Pythia8::Event &event_hadronize) {
  const auto &log = logger<LogArea::Pythia>();

  if (event_intermediate.sizeJunction() == 0) {
    return;
  }

  event_hadronize.reset();

  const int kind = event_intermediate.kindJunction(0);
  bool sign_color = kind % 2 == 1;
  std::vector<int> col;
  col.clear();
  for (int j = 0; j < 3; j++) {
    col.push_back(event_intermediate.colJunction(0, j));
  }
  event_hadronize.appendJunction(kind, col[0], col[1], col[2]);
  event_intermediate.eraseJunction(0);
  log.debug("junction (", col[0], ", ", col[1], ", ", col[2],
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
      /* if there is any leg which is not associated with parton,
       * find another junction(s) connected. */
      log.debug("  still has leg(s) unfinished.");
      sign_color = !sign_color;
      std::vector<int> junction_to_move;
      junction_to_move.clear();
      for (int i = 0; i < event_intermediate.sizeJunction(); i++) {
        const int kind_new = event_intermediate.kindJunction(i);
        if (sign_color != (kind_new % 2 == 1)) {
          continue;
        }

        std::array<int, 3> col_new;
        for (int k = 0; k < 3; k++) {
          col_new[k] = event_intermediate.colJunction(i, k);
        }

        int n_legs_connected = 0;
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

        if (n_legs_connected > 0) {
          for (int k = 0; k < 3; k++) {
            if (col_new[k] != 0) {
              col.push_back(col_new[k]);
            }
          }
          log.debug("  junction ", i, " (",
                    event_intermediate.colJunction(i, 0), ", ",
                    event_intermediate.colJunction(i, 1), ", ",
                    event_intermediate.colJunction(i, 2),
                    ") with kind ", kind_new, " will be added.");
          junction_to_move.push_back(i);
        }
      }

      for (unsigned int i = 0; i < junction_to_move.size(); i++) {
        unsigned int imove = junction_to_move[i] - i;
        const int kind_add = event_intermediate.kindJunction(imove);
        std::array<int, 3> col_add;
        for (int k = 0; k < 3; k++) {
          col_add[k] = event_intermediate.colJunction(imove, k);
        }
        // add a junction to the new event record.
        event_hadronize.appendJunction(kind_add,
                                       col_add[0], col_add[1], col_add[2]);
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
  const auto &log = logger<LogArea::Pythia>();

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
          log.debug("  col[", j, "] = ", col[j],
                    " from i ", i, "(", pdgid, ") found");
          ifound = i;
          col[j] = event_intermediate[i].acol();
          break;
        } else if (!sign_color && col[j] == event_intermediate[i].acol()) {
          log.debug("  acol[", j, "] = ", col[j],
                    " from i ", i, "(", pdgid, ") found");
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
          throw std::runtime_error("Hard string could not be identified.");
        }
      } else {
        pSum += event_intermediate[ifound].p();
        // add a parton to the new event record.
        event_hadronize.append(event_intermediate[ifound]);
        // then remove from the original event record.
        event_intermediate.remove(ifound, ifound);
        log.debug("Hard non-diff: event_intermediate reduces in size to ",
                  event_intermediate.size());
        if (col[j] == 0) {
          found_leg = true;
        }
      }
    }
  }

  event_hadronize[0].p(pSum);
  event_hadronize[0].m(pSum.mCalc());
}

// baryon-antibaryon annihilation
bool StringProcess::next_BBbarAnn() {
  const auto &log = logger<LogArea::Pythia>();
  const std::array<FourVector, 2> ustrcom = {FourVector(1., 0., 0., 0.),
                                             FourVector(1., 0., 0., 0.)};

  NpartFinal_ = 0;
  NpartString_[0] = 0;
  NpartString_[1] = 0;
  final_state_.clear();

  log.debug("Annihilation occurs between ", PDGcodes_[0], "+", PDGcodes_[1],
            " at CM energy [GeV] ", sqrtsAB_);

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
    const int nfrag = fragment_string(
        remaining_quarks[i], remaining_antiquarks[i],
        mstr[i], evec, true, false, new_intermediate_particles);
    if (nfrag <= 0) {
      NpartString_[i] = 0;
      return false;
    }
    NpartString_[i] = append_final_state(new_intermediate_particles,
                                         ustrcom[i], evec);
  }
  NpartFinal_ = NpartString_[0] + NpartString_[1];
  return true;
}

void StringProcess::make_orthonormal_basis(
         ThreeVector &evec_polar, std::array<ThreeVector, 3> &evec_basis) {
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
    evec_basis[1].set_x1(cos(theta) * cos(phi));
    evec_basis[1].set_x2(cos(theta) * sin(phi));
    evec_basis[1].set_x3(-sin(theta));

    evec_basis[2].set_x1(-sin(phi));
    evec_basis[2].set_x2(cos(phi));
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
}

void StringProcess::compute_incoming_lightcone_momenta() {
  PPosA_ = (pcom_[0].x0() + evecBasisAB_[0] * pcom_[0].threevec()) / sqrt2_;
  PNegA_ = (pcom_[0].x0() - evecBasisAB_[0] * pcom_[0].threevec()) / sqrt2_;
  PPosB_ = (pcom_[1].x0() + evecBasisAB_[0] * pcom_[1].threevec()) / sqrt2_;
  PNegB_ = (pcom_[1].x0() - evecBasisAB_[0] * pcom_[1].threevec()) / sqrt2_;
}

void StringProcess::quarks_from_diquark(int diquark,
                                        int &q1, int &q2, int &deg_spin) {
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

void StringProcess::make_string_ends(const PdgCode &pdg, int &idq1, int &idq2) {
  std::array<int, 3> quarks = pdg.quark_content();

  if (pdg.is_meson()) {
    idq1 = quarks[1];
    idq2 = quarks[2];
    /* Some mesons with PDG id 11X are actually mixed state of uubar and ddbar.
     * have a random selection whether we have uubar or ddbar in this case. */
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
  // Fulfil the convention: idq1 should be quark or anti-diquark
  if (idq1 < 0) {
    std::swap(idq1, idq2);
  }
}

int StringProcess::fragment_string(int idq1, int idq2, double mString,
                                   ThreeVector &evecLong,
                                   bool flip_string_ends,
                                   bool separate_fragment_baryon,
                                   ParticleList &intermediate_particles) {
  const auto &log = logger<LogArea::Pythia>();
  pythia_hadron_->event.reset();
  intermediate_particles.clear();

  log.debug("initial quark content for fragment_string : ", idq1, ", ", idq2);
  std::array<int, 2> idqIn;
  idqIn[0] = idq1;
  idqIn[1] = idq2;

  int bstring = 0;
  std::array<double, 2> m_const;

  for (int i = 0; i < 2; i++) {
    // evaluate 3 times total baryon number of the string
    bstring += pythia_hadron_->particleData.baryonNumberType(idqIn[i]);

    m_const[i] = pythia_hadron_->particleData.m0(idqIn[i]);
  }
  log.debug("3 times baryon number of string : ", bstring);

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

  if (separate_fragment_baryon &&
      (std::abs(bstring) == 3) && (mString > m_const[0] + m_const[1] + 1.)) {
    const double ssbar_supp = pythia_hadron_->parm("StringFlav:probStoUD");
    std::array<ThreeVector, 3> evec_basis;
    make_orthonormal_basis(evecLong, evec_basis);

    double QTrx, QTry, QTrn;
    double ppos_string_new, pneg_string_new;
    double mTrn_string_new;
    std::array<double, 2> m_trans;

    const int niter_max = 10000;
    int iiter = 0;
    bool found_leading_baryon = false;
    while (!found_leading_baryon) {
      int idnew_qqbar = 0;
      if (random::uniform(0., 1. + ssbar_supp) < 1.) {
        if (random::uniform_int(0, 1) == 0) {
          idnew_qqbar = 1;
        } else {
          idnew_qqbar = 2;
        }
      } else {
        idnew_qqbar = 3;
      }

      int id_diquark = 0;
      if (bstring > 0) {
        id_diquark = std::abs(idq2);
        idqIn[1] = -idnew_qqbar;
      } else {
        id_diquark = std::abs(idq1);
        idqIn[0] = idnew_qqbar;
      }
      log.debug("quark constituents for leading baryon : ",
                idnew_qqbar, ", ", id_diquark);

      std::array<int, 5> frag_net_q;
      for (int iq = 0; iq < 5; iq++) {
        frag_net_q[iq] = (bstring > 0 ? 1 : -1) *
                         (pythia_hadron_->particleData.nQuarksInCode(
                          idnew_qqbar, iq + 1) +
                          pythia_hadron_->particleData.nQuarksInCode(
                          id_diquark, iq + 1));
      }
      const int frag_iso3 = frag_net_q[1] - frag_net_q[0];
      const int frag_strange = -frag_net_q[2];
      const int frag_charm = frag_net_q[3];
      const int frag_bottom = -frag_net_q[4];
      log.debug("  conserved charges of leading baryon : iso3 = ", frag_iso3,
                ", strangeness = ", frag_strange,
                ", charmness = ", frag_charm,
                ", bottomness = ", frag_bottom);

      std::vector<int> pdgid_possible;
      std::vector<double> weight_possible;
      std::vector<double> weight_summed;
      pdgid_possible.clear();
      weight_possible.clear();
      weight_summed.clear();
      for (auto &ptype : ParticleType::list_all()) {
        const int pdgid = (bstring > 0 ? 1 : -1) *
                          std::abs(ptype.pdgcode().get_decimal());
        if ((pythia_hadron_->particleData.isParticle(pdgid)) &&
            (bstring == 3 * ptype.pdgcode().baryon_number()) &&
            (frag_iso3 == ptype.pdgcode().isospin3()) &&
            (frag_strange == ptype.pdgcode().strangeness()) &&
            (frag_charm == ptype.pdgcode().charmness()) &&
            (frag_bottom == ptype.pdgcode().bottomness())) {
          const int spin_degeneracy = ptype.pdgcode().spin_degeneracy();
          const double mass_pole = ptype.mass();
          const double weight = static_cast<double>(spin_degeneracy) /
                                mass_pole;
          pdgid_possible.push_back(pdgid);
          weight_possible.push_back(weight);

          log.debug("  PDG id ", pdgid, " is possible with weight ", weight);
        }
      }
      const int n_possible = pdgid_possible.size();
      weight_summed.push_back(0.);
      for (int i = 0; i < n_possible; i++) {
        weight_summed.push_back(weight_summed[i] + weight_possible[i]);
      }

      int pdgid_frag = 0;
      const double uspc = random::uniform(0., weight_summed[n_possible]);
      for (int i = 0; i < n_possible; i++) {
        if ((uspc >= weight_summed[i]) && (uspc < weight_summed[i + 1])) {
          pdgid_frag = pdgid_possible[i];
          log.debug("PDG id ", pdgid_frag, " selected for leading baryon.");
          break;
        }
      }
      const double mass_frag = pythia_hadron_->particleData.mSel(pdgid_frag);

      QTrx = random::normal(0., sigma_qperp_ / sqrt2_);
      QTry = random::normal(0., sigma_qperp_ / sqrt2_);
      QTrn = std::sqrt(QTrx * QTrx + QTry * QTry);
      const double mTrn_frag = std::sqrt(QTrn * QTrn + mass_frag * mass_frag);

      double xfrac = 0.;
      const double xf_const_A = leading_frag_width_;
      const double xf_const_B = leading_frag_mean_;
      bool xfrac_accepted = false;
      while (!xfrac_accepted) {
        const double angle = random::uniform(0., 1.) *
                             (std::atan(0.5 * (1. - xf_const_B) / xf_const_A) +
                              std::atan(0.5 * xf_const_B / xf_const_A)) -
                             std::atan(0.5 * xf_const_B / xf_const_A);
        xfrac = xf_const_B + 2. * xf_const_A * std::tan(angle);

        const double xf_tmp = std::abs(xfrac - xf_const_B) / xf_const_A;
        const double xf_env = (1. + really_small) / (1. + xf_tmp * xf_tmp / 4.);
        const double xf_pdf = std::exp(-xf_tmp * xf_tmp / 2.);
        if (random::uniform(0., xf_env) < xf_pdf &&
            xfrac > 0. && xfrac < 1.) {
          xfrac_accepted = true;
        }
      }

      const double ppos_frag = xfrac * mString / sqrt2_;
      const double pneg_frag = 0.5 * mTrn_frag * mTrn_frag / ppos_frag;
      ppos_string_new = mString / sqrt2_ - ppos_frag;
      pneg_string_new = mString / sqrt2_ - pneg_frag;
      mTrn_string_new = std::sqrt(std::max(0.,
                        2. * ppos_string_new * pneg_string_new));

      for (int i = 0; i < 2; i++) {
        m_const[i] = pythia_hadron_->particleData.m0(idqIn[i]);
      }
      if (bstring > 0) {
        m_trans[0] = m_const[0];
        m_trans[1] = std::sqrt(m_const[1] * m_const[1] + QTrn * QTrn);
      } else {
        m_trans[0] = std::sqrt(m_const[0] * m_const[0] + QTrn * QTrn);
        m_trans[1] = m_const[1];
      }

      if (mTrn_string_new > m_trans[0] + m_trans[1]) {
        found_leading_baryon = true;

        FourVector mom_frag((ppos_frag + pneg_frag) / sqrt2_,
                            evec_basis[0] * (ppos_frag - pneg_frag) / sqrt2_ +
                            evec_basis[1] * QTrx + evec_basis[2] * QTry);
        log.debug("appending the leading baryon ", pdgid_frag,
                  " to the intermediate particle list.");
        bool found_ptype = append_intermediate_list(pdgid_frag, mom_frag,
                                                    intermediate_particles);
        if (!found_ptype) {
          log.error("PDG ID ", pdgid_frag, " should exist in ParticleType.");
          throw std::runtime_error("string fragmentation failed.");
        }
        number_of_fragments++;
        log.debug("proceed to the next step");
      }

      if (iiter == niter_max) {
        return 0;
      }
      iiter += 1;
    }

    std::array<double, 2> ppos_parton;
    std::array<double, 2> pneg_parton;

    const double pb_const = (mTrn_string_new * mTrn_string_new +
                            m_trans[0] * m_trans[0] - m_trans[1] * m_trans[1]) /
                           (4. * pneg_string_new);
    const double pc_const = 0.5 * m_trans[0] * m_trans[0] *
                           ppos_string_new / pneg_string_new;
    ppos_parton[0] = pb_const + (bstring > 0 ? -1. : 1.) *
                     std::sqrt(pb_const * pb_const - pc_const);
    ppos_parton[1] = ppos_string_new - ppos_parton[0];

    for (int i = 0; i < 2; i++) {
      pneg_parton[i] = 0.5 * m_trans[i] * m_trans[i] / ppos_parton[i];
    }

    const int status = 1;
    int color, anticolor;
    ThreeVector three_mom;
    Pythia8::Vec4 pquark;

    color = 1;
    anticolor = 0;
    if (bstring > 0) {
      three_mom = evec_basis[0] * (ppos_parton[0] - pneg_parton[0]) / sqrt2_;
    } else {
      three_mom = evec_basis[0] * (ppos_parton[0] - pneg_parton[0]) / sqrt2_ -
                  evec_basis[1] * QTrx - evec_basis[2] * QTry;
    }
    pquark = set_Vec4((ppos_parton[0] + pneg_parton[0]) / sqrt2_,
                      three_mom);
    pSum += pquark;
    pythia_hadron_->event.append(idqIn[0], status, color, anticolor,
                                 pquark, m_const[0]);

    color = 0;
    anticolor = 1;
    if (bstring > 0) {
      three_mom = evec_basis[0] * (ppos_parton[1] - pneg_parton[1]) / sqrt2_ -
                  evec_basis[1] * QTrx - evec_basis[2] * QTry;
    } else {
      three_mom = evec_basis[0] * (ppos_parton[1] - pneg_parton[1]) / sqrt2_;
    }
    pquark = set_Vec4((ppos_parton[1] + pneg_parton[1]) / sqrt2_,
                      three_mom);
    pSum += pquark;
    pythia_hadron_->event.append(idqIn[1], status, color, anticolor,
                                 pquark, m_const[1]);
  } else {
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
    pythia_hadron_->event.append(idqIn[0], status1, color1, anticolor1,
                                 pquark, m_const[0]);

    const int status2 = 1, color2 = 0, anticolor2 = 1;
    pquark = set_Vec4(E2, direction * pCMquark);
    pSum += pquark;
    pythia_hadron_->event.append(idqIn[1], status2, color2, anticolor2,
                                 pquark, m_const[1]);
  }

  log.debug("fragmenting a string with ", idqIn[0], ", ", idqIn[1]);
  // implement PYTHIA fragmentation
  pythia_hadron_->event[0].p(pSum);
  pythia_hadron_->event[0].m(pSum.mCalc());
  pythia_hadron_->rndm.init(random::uniform_int(1,
                                std::numeric_limits<int>::max()));
  const bool successful_hadronization = pythia_hadron_->next();
  if (successful_hadronization) {
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
      log.debug("appending the fragmented hadron ", pythia_id,
                " to the intermediate particle list.");
      bool found_ptype = append_intermediate_list(pythia_id, momentum,
                                                  intermediate_particles);
      if (!found_ptype) {
        log.warn("PDG ID ", pythia_id,
                 " does not exist in ParticleType - start over.");
        intermediate_particles.clear();
        return 0;
      }

      number_of_fragments++;
    }

    return number_of_fragments;
  } else {
    return 0;
  }
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
    if (pdg.charge() > 0) {
      pdg_mapped = PdgCode(pdg::pi_p);
    } else if (pdg.charge() < 0) {
      pdg_mapped = PdgCode(pdg::pi_m);
    } else {
      pdg_mapped = PdgCode(pdg::pi_z);
    }
  } else if (pdg.is_lepton()) {  // lepton
    pdg_mapped = pdg.charge() < 0 ? PdgCode(0x11) : PdgCode(-0x11);
  } else {
    log.error("PdgCode", pdg,
              " could not be mapped onto appropriate PDG id for PYTHIA.");
    throw std::runtime_error("StringProcess::pdg_map_for_pythia failed.");
  }

  return pdg_mapped.get_decimal();
}

}  // namespace smash
