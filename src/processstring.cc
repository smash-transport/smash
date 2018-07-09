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
                             double sigma_perp, double stringz_a,
                             double stringz_b, double string_sigma_T,
                             double factor_t_form, bool use_yoyo_model)
    : pmin_gluon_lightcone_(gluon_pmin),
      pow_fgluon_beta_(gluon_beta),
      pow_fquark_alpha_(quark_alpha),
      pow_fquark_beta_(quark_beta),
      sigma_qperp_(sigma_perp),
      kappa_tension_string_(string_tension),
      additional_xsec_supp_(0.7),
      time_formation_const_(time_formation),
      soft_t_form_(factor_t_form),
      time_collision_(0.),
      gamma_factor_com_(1.),
      use_yoyo_model_(use_yoyo_model) {
  // setup and initialize pythia for hard string process
  pythia_parton_ = make_unique<Pythia8::Pythia>(PYTHIA_XML_DIR, false);
  /* select only non-diffractive events
   * diffractive ones are implemented in a separate routine */
  pythia_parton_->readString("SoftQCD:nonDiffractive = on");
  pythia_parton_->readString("MultipartonInteractions:pTmin = 1.5");
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
int StringProcess::append_final_state(const FourVector &uString,
                                      const ThreeVector &evecLong) {
  struct fragment_type {
    FourVector momentum;
    double pparallel;  // longitudinal component of momentum
    double y;          // momentum rapidity
    double xtotfac;    // cross-section scaling factor
    PdgCode pdg;
    bool is_leading;  // is leading hadron or not
  };

  std::vector<fragment_type> fragments;
  double p_pos_tot = 0.0, p_neg_tot = 0.0;
  int bstring = 0;

  for (int ipyth = 0; ipyth < pythia_hadron_->event.size(); ipyth++) {
    if (!pythia_hadron_->event[ipyth].isFinal()) {
      continue;
    }
    int pythia_id = pythia_hadron_->event[ipyth].id();
    /* K_short and K_long need are converted to K0 since SMASH only knows K0 */
    convert_KaonLS(pythia_id);
    PdgCode pdg = PdgCode::from_decimal(pythia_id);
    if (!pdg.is_hadron()) {
      throw std::invalid_argument(
          "StringProcess::append_final_state warning :"
          " particle is not meson or baryon.");
    }
    FourVector mom(
        pythia_hadron_->event[ipyth].e(), pythia_hadron_->event[ipyth].px(),
        pythia_hadron_->event[ipyth].py(), pythia_hadron_->event[ipyth].pz());
    double pparallel = mom.threevec() * evecLong;
    double y = 0.5 * std::log((mom.x0() + pparallel) / (mom.x0() - pparallel));
    fragments.push_back({mom, pparallel, y, 0.0, pdg, false});
    // total lightcone momentum
    p_pos_tot += (mom.x0() + pparallel) / sqrt2_;
    p_neg_tot += (mom.x0() - pparallel) / sqrt2_;
    bstring += pythia_hadron_->particleData.baryonNumberType(pythia_id);
  }
  const int nfrag = fragments.size();
  assert(nfrag > 0);
  // Sort the particles in descending order of momentum rapidity
  std::sort(fragments.begin(), fragments.end(),
            [&](fragment_type i, fragment_type j) { return i.y > j.y; });

  std::vector<double> xvertex_pos, xvertex_neg;
  xvertex_pos.resize(nfrag + 1);
  xvertex_neg.resize(nfrag + 1);
  // x^{+} coordinates of the forward end
  xvertex_pos[0] = p_pos_tot / kappa_tension_string_;
  for (int i = 0; i < nfrag; i++) {
    // recursively compute x^{+} coordinates of q-qbar formation vertex
    xvertex_pos[i + 1] =
        xvertex_pos[i] - (fragments[i].momentum.x0() + fragments[i].pparallel) /
                             (kappa_tension_string_ * sqrt2_);
  }
  // x^{-} coordinates of the backward end
  xvertex_neg[nfrag] = p_neg_tot / kappa_tension_string_;
  for (int i = nfrag - 1; i >= 0; i--) {
    // recursively compute x^{-} coordinates of q-qbar formation vertex
    xvertex_neg[i] = xvertex_neg[i + 1] -
                     (fragments[i].momentum.x0() - fragments[i].pparallel) /
                         (kappa_tension_string_ * sqrt2_);
  }

  /* compute the cross section suppression factor for leading hadrons
   * based on the number of valence quarks. */
  if (bstring == 0) {  // mesonic string
    for (int i : {0, nfrag - 1}) {
      fragments[i].xtotfac = fragments[i].pdg.is_baryon() ? 1. / 3. : 0.5;
      fragments[i].is_leading = true;
    }
  } else if (bstring == 3 || bstring == -3) {  // baryonic string
    // The first baryon in forward direction
    int i = 0;
    while (i < nfrag && 3 * fragments[i].pdg.baryon_number() != bstring) {
      i++;
    }
    fragments[i].xtotfac = 2. / 3.;
    fragments[i].is_leading = true;
    // The most backward meson
    i = nfrag - 1;
    while (i >= 0 && 3 * fragments[i].pdg.baryon_number() == bstring) {
      i--;
    }
    fragments[i].xtotfac = fragments[i].pdg.is_baryon() ? 1. / 3. : 0.5;
    fragments[i].is_leading = true;
  } else {
    throw std::invalid_argument(
        "StringProcess::append_final_state"
        " encountered bstring != 0, 3, -3");
  }

  // Velocity three-vector to perform Lorentz boost.
  const ThreeVector vstring = uString.velocity();

  /* compute the formation times of hadrons
   * from the lightcone coordinates of q-qbar formation vertices. */
  for (int i = 0; i < nfrag; i++) {
    /* create new particle with specific PDG id
     * and assign momentum. */
    ParticleData new_particle(ParticleType::find(fragments[i].pdg));
    fragments[i].momentum = fragments[i].momentum.LorentzBoost(-vstring);
    new_particle.set_4momentum(fragments[i].momentum);

    /* additional suppression factor to take
     * the quantum coherence effect into account. */
    new_particle.set_cross_section_scaling_factor(
        fragments[i].is_leading ? additional_xsec_supp_ * fragments[i].xtotfac
                                : 0.);

    if (use_yoyo_model_) {
      // set the formation time and position in the rest frame of string
      double t_prod = (xvertex_pos[i] + xvertex_neg[i + 1]) / sqrt2_;
      FourVector fragment_position = FourVector(
          t_prod, evecLong * (xvertex_pos[i] - xvertex_neg[i + 1]) / sqrt2_);
      // boost into the center of mass frame
      fragment_position = fragment_position.LorentzBoost(-vstring);
      t_prod = fragment_position.x0();
      new_particle.set_formation_time(
          time_collision_ + soft_t_form_ * gamma_factor_com_ * t_prod);
    } else {
      ThreeVector v_calc =
          new_particle.momentum().LorentzBoost(-vcomAB_).velocity();
      double gamma_factor = 1.0 / std::sqrt(1 - (v_calc).sqr());
      new_particle.set_formation_time(time_formation_const_ * gamma_factor +
                                      time_collision_);
    }
    final_state_.push_back(new_particle);
  }

  return nfrag;
}

void StringProcess::init(const ParticleList &incoming, double tcoll,
                         double gamma) {
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

  make_orthonormal_basis();
  compute_incoming_lightcone_momenta();

  time_collision_ = tcoll;
  gamma_factor_com_ = gamma;
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
  int nfrag = fragment_string(idqX1, idqX2, massX, evec, true);
  if (nfrag < 1) {
    NpartString_[0] = 0;
    return false;
  }
  NpartString_[0] = append_final_state(ustrXcom, evec);

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
    const std::array<ThreeVector, 2> &evec_str, bool flip_string_ends) {
  const std::array<FourVector, 2> ustr_com = {pstr_com[0] / m_str[0],
                                              pstr_com[1] / m_str[1]};
  for (int i = 0; i < 2; i++) {
    // determine direction in which string i is stretched.
    ThreeVector evec = evec_str[i];
    // perform fragmentation and add particles to final_state.
    int nfrag = fragment_string(quarks[i][0], quarks[i][1], m_str[i], evec,
                                flip_string_ends);
    if (nfrag <= 0) {
      NpartString_[i] = 0;
      return false;
    }
    NpartString_[i] = append_final_state(ustr_com[i], evec);
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
  const bool success = make_final_state_2strings(quarks, pstr_com, m_str,
                                                 evec_str, flip_string_ends);
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
  const bool success = make_final_state_2strings(quarks, pstr_com, m_str,
                                                 evec_str, flip_string_ends);
  return success;
}

// hard non-diffractive
bool StringProcess::next_NDiffHard() {
  const auto &log = logger<LogArea::Pythia>();
  final_state_.clear();

  std::array<bool, 2> accepted_by_pythia;
  for (int i = 0; i < 2; i++) {
    int pdgid = PDGcodes_[i].get_decimal();
    accepted_by_pythia[i] = pdgid == 2212 || pdgid == -2212 || pdgid == 2112 ||
                            pdgid == -2112 || pdgid == 211 || pdgid == 111 ||
                            pdgid == -211;
  }
  if (!accepted_by_pythia[0] || !accepted_by_pythia[1]) {
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

  int previous_idA = pythia_parton_->mode("Beams:idA"),
      previous_idB = pythia_parton_->mode("Beams:idB");
  double previous_eCM = pythia_parton_->parm("Beams:eCM");

  bool same_initial_state = previous_idA == PDGcodes_[0].get_decimal() &&
                            previous_idB == PDGcodes_[1].get_decimal() &&
                            std::abs(previous_eCM - sqrtsAB_) < really_small;

  /* Perform PYTHIA initialization if it was not previously initialized
   * or the initial state changed. */
  if (!pythia_parton_initialized_ || !same_initial_state) {
    pythia_parton_->settings.mode("Beams:idA", PDGcodes_[0].get_decimal());
    pythia_parton_->settings.mode("Beams:idB", PDGcodes_[1].get_decimal());
    pythia_parton_->settings.parm("Beams:eCM", sqrtsAB_);

    pythia_parton_initialized_ = pythia_parton_->init();
    log.info("Pythia initialized with ", PDGcodes_[0], "+", PDGcodes_[1],
             " at CM energy [GeV] ", sqrtsAB_);
    if (!pythia_parton_initialized_) {
      throw std::runtime_error("Pythia failed to initialize.");
    }
  }
  /* Set the random seed of the Pythia random Number Generator.
   * Pythia's random is controlled by SMASH in every single collision.
   * In this way we ensure that the results are reproducible
   * for every event if one knows SMASH random seed. */
  const int maxint = std::numeric_limits<int>::max();
  pythia_parton_->rndm.init(random::uniform_int(1, maxint));

  // Short notation for Pythia event
  Pythia8::Event &event = pythia_parton_->event;
  log.debug("Pythia hard event created");
  bool final_state_success = false;
  while (!final_state_success) {
    final_state_success = pythia_parton_->next();
    log.debug("Pythia final state computed, success = ", final_state_success);
  }

  ParticleList new_intermediate_particles;
  ParticleList new_non_hadron_particles;
  for (int i = 0; i < event.size(); i++) {
    if (event[i].isFinal()) {
      int pythia_id = event[i].id();
      log.debug("PDG ID from Pythia:", pythia_id);
      /* K_short and K_long need to be converted to K0
       * since SMASH only knows K0 */
      convert_KaonLS(pythia_id);
      const std::string s = std::to_string(pythia_id);
      PdgCode pythia_code(s);
      ParticleData new_particle(ParticleType::find(pythia_code));
      /* evecBasisAB_[0] is a unit 3-vector in the collision axis,
       * while evecBasisAB_[1] and evecBasisAB_[2] spans the transverse plane.
       * Given that PYTHIA assumes z-direction to be the collision axis,
       * pz from PYTHIA should be the momentum compoment in evecBasisAB_[0].
       * px and py are respectively the momentum components in two
       * transverse directions evecBasisAB_[1] and evecBasisAB_[2]. */
      ThreeVector threeMomentum = evecBasisAB_[0] * event[i].pz() +
                                  evecBasisAB_[1] * event[i].px() +
                                  evecBasisAB_[2] * event[i].py();
      FourVector momentum = FourVector(event[i].e(), threeMomentum);
      new_particle.set_4momentum(momentum);
      log.debug("4-momentum from Pythia: ", momentum);
      if (event[i].isHadron()) {
        new_intermediate_particles.push_back(new_particle);
      } else {
        new_non_hadron_particles.push_back(new_particle);
      }
    }
  }

  /* Additional suppression factor to mimic coherence taken as 0.7
   * from UrQMD (CTParam(59) */
  const int baryon_string =
      PDGcodes_[random::uniform_int(0, 1)].baryon_number();
  assign_all_scaling_factors(baryon_string, new_intermediate_particles,
                             evecBasisAB_[0], additional_xsec_supp_);
  for (ParticleData data : new_intermediate_particles) {
    log.debug("Particle momenta after sorting: ", data.momentum());
    /* The hadrons are not immediately formed, currently a formation time of
     * 1 fm is universally applied and cross section is reduced to zero and
     * to a fraction corresponding to the valence quark content. Hadrons
     * containing a valence quark are determined by highest z-momentum. */
    log.debug("The formation time is: ", time_formation_const_, "fm/c.");
    ThreeVector v_calc =
        (data.momentum().LorentzBoost(-1.0 * vcomAB_)).velocity();
    /* Set formation time: actual time of collision + time to form the
     * particle */
    double gamma_factor = 1.0 / std::sqrt(1 - (v_calc).sqr());
    data.set_slow_formation_times(
        time_collision_,
        time_formation_const_ * gamma_factor + time_collision_);
    final_state_.push_back(data);
  }

  for (ParticleData data : new_non_hadron_particles) {
    data.set_cross_section_scaling_factor(1.);
    data.set_formation_time(time_collision_);
    final_state_.push_back(data);
  }

  return final_state_success;
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
    ThreeVector evec = pcom_[i].threevec() / pcom_[i].threevec().abs();
    const int nfrag = fragment_string(
        remaining_quarks[i], remaining_antiquarks[i], mstr[i], evec, true);
    if (nfrag <= 0) {
      NpartString_[i] = 0;
      return false;
    }
    NpartString_[i] = append_final_state(ustrcom[i], evec);
  }
  NpartFinal_ = NpartString_[0] + NpartString_[1];
  return true;
}

void StringProcess::make_orthonormal_basis() {
  const double pabscomAB = pCM(sqrtsAB_, massA_, massB_);
  if (std::abs(pcom_[0].x3()) < (1. - 1.0e-8) * pabscomAB) {
    double ex, ey, et;
    double theta, phi;

    /* evecBasisAB_[0] is set to be longitudinal direction
     * which is parallel to the collision axis. */
    evecBasisAB_[0] = pcom_[0].threevec() / pabscomAB;

    theta = std::acos(evecBasisAB_[0].x3());

    ex = evecBasisAB_[0].x1();
    ey = evecBasisAB_[0].x2();
    et = std::sqrt(ex * ex + ey * ey);
    if (ey > 0.) {
      phi = std::acos(ex / et);
    } else {
      phi = -std::acos(ex / et);
    }

    /* The transverse plane is spanned
     * by evecBasisAB_[1] and evecBasisAB_[2]. */
    evecBasisAB_[1].set_x1(cos(theta) * cos(phi));
    evecBasisAB_[1].set_x2(cos(theta) * sin(phi));
    evecBasisAB_[1].set_x3(-sin(theta));

    evecBasisAB_[2].set_x1(-sin(phi));
    evecBasisAB_[2].set_x2(cos(phi));
    evecBasisAB_[2].set_x3(0.);
  } else {
    // if pcom_[0].threevec() is very close to the z axis
    if (pcom_[0].x3() > 0.) {
      evecBasisAB_[1] = ThreeVector(1., 0., 0.);
      evecBasisAB_[2] = ThreeVector(0., 1., 0.);
      evecBasisAB_[0] = ThreeVector(0., 0., 1.);
    } else {
      evecBasisAB_[1] = ThreeVector(0., 1., 0.);
      evecBasisAB_[2] = ThreeVector(1., 0., 0.);
      evecBasisAB_[0] = ThreeVector(0., 0., -1.);
    }
  }
}

void StringProcess::compute_incoming_lightcone_momenta() {
  PPosA_ = (pcom_[0].x0() + evecBasisAB_[0] * pcom_[0].threevec()) / sqrt2_;
  PNegA_ = (pcom_[0].x0() - evecBasisAB_[0] * pcom_[0].threevec()) / sqrt2_;
  PPosB_ = (pcom_[1].x0() + evecBasisAB_[0] * pcom_[1].threevec()) / sqrt2_;
  PNegB_ = (pcom_[1].x0() - evecBasisAB_[0] * pcom_[1].threevec()) / sqrt2_;
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
                                   bool flip_string_ends) {
  pythia_hadron_->event.reset();
  // evaluate 3 times total baryon number of the string
  const int bstring = pythia_hadron_->particleData.baryonNumberType(idq1) +
                      pythia_hadron_->particleData.baryonNumberType(idq2);
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

  const double m1 = pythia_hadron_->particleData.m0(idq1);
  const double m2 = pythia_hadron_->particleData.m0(idq2);
  if (m1 + m2 > mString) {
    throw std::runtime_error("String fragmentation: m1 + m2 > mString");
  }

  // evaluate momenta of quarks
  const double pCMquark = pCM(mString, m1, m2);
  const double E1 = std::sqrt(m1 * m1 + pCMquark * pCMquark);
  const double E2 = std::sqrt(m2 * m2 + pCMquark * pCMquark);

  if (flip_string_ends && random::uniform_int(0, 1) == 0) {
    /* in the case where we flip the string ends,
     * we need to flip the longitudinal unit vector itself
     * since it is set to be direction of diquark (anti-quark)
     * or anti-diquark. */
    evecLong = -evecLong;
  }
  ThreeVector direction = sign_direction * evecLong;

  Pythia8::Vec4 pSum = 0.;

  // For status and (anti)color see \iref{Sjostrand:2007gs}.
  const int status1 = 1, color1 = 1, anticolor1 = 0;
  Pythia8::Vec4 pquark = set_Vec4(E1, -direction * pCMquark);
  pSum += pquark;
  pythia_hadron_->event.append(idq1, status1, color1, anticolor1, pquark, m1);

  const int status2 = 1, color2 = 0, anticolor2 = 1;
  pquark = set_Vec4(E2, direction * pCMquark);
  pSum += pquark;
  pythia_hadron_->event.append(idq2, status2, color2, anticolor2, pquark, m2);

  // implement PYTHIA fragmentation
  pythia_hadron_->event[0].p(pSum);
  pythia_hadron_->event[0].m(pSum.mCalc());
  const bool successful_hadronization = pythia_hadron_->forceHadronLevel();
  int number_of_fragments = 0;
  if (successful_hadronization) {
    for (int ipart = 0; ipart < pythia_hadron_->event.size(); ipart++) {
      if (pythia_hadron_->event[ipart].isFinal()) {
        number_of_fragments++;
      }
    }
  }

  return number_of_fragments;
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
                                               ThreeVector &evec_coll,
                                               double suppression_factor) {
  // Set each particle's cross section scaling factor to 0 first
  for (ParticleData &data : outgoing_particles) {
    data.set_cross_section_scaling_factor(0.0);
  }
  // sort outgoing particles according to z-velocity
  std::sort(outgoing_particles.begin(), outgoing_particles.end(),
            [&](ParticleData i, ParticleData j) {
              return i.momentum().velocity() * evec_coll <
                     j.momentum().velocity() * evec_coll;
            });
  int nq1, nq2;  // number of quarks at both ends of the string
  switch (baryon_string) {
    case 0:
      nq1 = 1;
      nq2 = -1;
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
  if (i.second - i.first > j.second - j.first) {
    assign_scaling_factor(nq1, outgoing_particles[i.first], suppression_factor);
    assign_scaling_factor(nq2, outgoing_particles[i.second],
                          suppression_factor);
  } else {
    assign_scaling_factor(nq2, outgoing_particles[j.first], suppression_factor);
    assign_scaling_factor(nq1, outgoing_particles[j.second],
                          suppression_factor);
  }
}

}  // namespace smash
