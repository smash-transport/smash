/*
 *
 *    Copyright (c) 2017-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <array>

#include "include/angles.h"
#include "include/kinematics.h"
#include "include/random.h"
#include "include/stringprocess.h"

namespace Smash {

StringProcess::StringProcess()
    : pmin_gluon_lightcone_(0.001),
      pow_fgluon_beta_(0.5),
      pow_fquark_alpha_(1.),
      pow_fquark_beta_(2.5),
      sigma_qperp_(0.5),
      kappa_tension_string_(1.0),
      time_collision_(0.),
      gamma_factor_com_(1.) {
  // setup and initialize pythia
  pythia_ = make_unique<Pythia8::Pythia>(PYTHIA_XML_DIR, false);
  /* select only inelastic events: */
  pythia_->readString("SoftQCD:inelastic = on");
  /* suppress unnecessary output */
  pythia_->readString("Print:quiet = on");
  /* No resonance decays, since the resonances will be handled by SMASH */
  pythia_->readString("HadronLevel:Decay = off");
  /* transverse momentum spread in string fragmentation */
  pythia_->readString("StringPT:sigma = 0.25");
  /* manually set the parton distribution function */
  pythia_->readString("PDF:pSet = 13");
  pythia_->readString("PDF:pSetB = 13");
  pythia_->readString("PDF:piSet = 1");
  pythia_->readString("PDF:piSetB = 1");
  pythia_->readString("Beams:idA = 2212");
  pythia_->readString("Beams:idB = 2212");
  pythia_->readString("Beams:eCM = 10.");
  /* set particle masses and widths in PYTHIA to be same with those in SMASH */
  for (auto &ptype : ParticleType::list_all()) {
    int pdgid = ptype.pdgcode().get_decimal();
    double mass_pole = ptype.mass();
    double width_pole = ptype.width_at_pole();
    /* check if the particle specie is in PYTHIA */
    if (pythia_->particleData.isParticle(pdgid)) {
      /* set mass and width in PYTHIA */
      pythia_->particleData.m0(pdgid, mass_pole);
      pythia_->particleData.mWidth(pdgid, width_pole);
    }
  }
  pythia_->init();
  pythia_sigmatot_.init(&pythia_->info, pythia_->settings,
                        &pythia_->particleData);

  sqrt2_ = std::sqrt(2.);

  for (int imu = 0; imu < 3; imu++) {
    evecBasisAB_[imu] = ThreeVector(0., 0., 0.);
  }

  final_state_.clear();
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

  for (int ipyth = 0; ipyth < pythia_->event.size(); ipyth++) {
    if (!pythia_->event[ipyth].isFinal()) {
      continue;
    }
    int pythia_id = pythia_->event[ipyth].id();
    /* K_short and K_long need are converted to K0 since SMASH only knows K0 */
    if (pythia_id == 310 || pythia_id == 130) {
      pythia_id = (Random::uniform_int(0, 1) == 0) ? 311 : -311;
    }
    PdgCode pdg = PdgCode::from_decimal(pythia_id);
    if (!pdg.is_hadron()) {
      throw std::invalid_argument(
          "StringProcess::append_final_state warning :"
          " particle is not meson or baryon.");
    }
    FourVector mom(pythia_->event[ipyth].e(), pythia_->event[ipyth].px(),
                   pythia_->event[ipyth].py(), pythia_->event[ipyth].pz());
    double pparallel = mom.threevec() * evecLong;
    double y = 0.5 * std::log((mom.x0() + pparallel) / (mom.x0() - pparallel));
    fragments.push_back({mom, pparallel, y, 0.0, pdg, false});
    // total lightcone momentum
    p_pos_tot += (mom.x0() + pparallel) / sqrt2_;
    p_neg_tot += (mom.x0() - pparallel) / sqrt2_;
    bstring += pythia_->particleData.baryonNumberType(pythia_id);
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
    xvertex_pos[i + 1] = xvertex_pos[i] -
                         (fragments[i].momentum.x0() + fragments[i].pparallel) /
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
    // set the formation time and position in the rest frame of string
    double t_prod = (xvertex_pos[i] + xvertex_neg[i + 1]) / sqrt2_;
    FourVector fragment_position = FourVector(
        t_prod,
        evecLong * (xvertex_pos[i] - xvertex_neg[i + 1]) / sqrt2_);
    // boost into the center of mass frame
    fragments[i].momentum = fragments[i].momentum.LorentzBoost(-vstring);
    fragment_position = fragment_position.LorentzBoost(-vstring);
    t_prod = fragment_position.x0();

    ParticleData new_particle(ParticleType::find(fragments[i].pdg));
    new_particle.set_4momentum(fragments[i].momentum);

    constexpr double suppression_factor = 0.7;
    new_particle.set_cross_section_scaling_factor(
        fragments[i].is_leading ? suppression_factor * fragments[i].xtotfac
                                : 0.);
    new_particle.set_formation_time(time_collision_ +
                                    gamma_factor_com_ * t_prod);
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

  sqrtsAB_ = (plab_[0] + plab_[1]).abs();
  /* Transverse momentum transferred to strings,
     parametrization to fit the experimental data */
  sigma_qperp_ = (sqrtsAB_ < 4.)
                     ? 0.5
                     : 0.5 + 0.2 * std::log(sqrtsAB_ / 4.) / std::log(5.);

  ucomAB_ = (plab_[0] + plab_[1]) / sqrtsAB_;
  vcomAB_ = ucomAB_.velocity();

  pcom_[0] = plab_[0].LorentzBoost(vcomAB_);
  pcom_[1] = plab_[1].LorentzBoost(vcomAB_);

  make_orthonormal_basis();
  compute_incoming_lightcone_momenta();

  time_collision_ = tcoll;
  gamma_factor_com_ = gamma;
}

/**
 * single diffractive
 * channel = 1 : A + B -> A + X
 * channel = 2 : A + B -> X + B
 */
bool StringProcess::next_SDiff(bool is_AB_to_AX) {
  NpartFinal_ = 0;
  NpartString_[0] = 0;
  NpartString_[1] = 0;
  final_state_.clear();

  double massH = is_AB_to_AX ? massA_ : massB_;
  double mstrMin = is_AB_to_AX ? massB_ : massA_;
  double mstrMax = sqrtsAB_ - massH;

  bool foundPabsX = false;
  int ntry = 0;
  int idqX1, idqX2;
  double QTrn, QTrx, QTry;
  double pabscomHX_sqr, massX;
  while (!foundPabsX && ntry < 100) {
    ntry++;
    // decompose hadron into quarks
    make_string_ends(is_AB_to_AX ? PDGcodes_[1] : PDGcodes_[0], idqX1, idqX2);
    // string mass must be larger than threshold set by PYTHIA.
    mstrMin = pythia_->particleData.m0(idqX1) +
              pythia_->particleData.m0(idqX2);
    // this threshold cannot be larger than maximum of allowed string mass.
    if (mstrMin > mstrMax) {
      continue;
    }
    // sample the transverse momentum transfer
    QTrx = Random::normal(0., sigma_qperp_ / sqrt2_);
    QTry = Random::normal(0., sigma_qperp_ / sqrt2_);
    QTrn = std::sqrt(QTrx * QTrx + QTry * QTry);
    // sample the string mass and evaluate the three-momenta of hadron and
    // string.
    massX = Random::power(-1.0, mstrMin, mstrMax);
    pabscomHX_sqr = pCM_sqr(sqrtsAB_, massH, massX);
    // magnitude of the three momentum must be larger than the transverse
    // momentum.
    foundPabsX = pabscomHX_sqr > QTrn * QTrn;
  }

  if (!foundPabsX) {
    return false;
  }
  double sign_direction = is_AB_to_AX ? 1. : -1.;
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
   * this is set to be same with the three-momentum of string
   * in the center of mass frame. */
  const ThreeVector threeMomentum = pstrXcom.threevec();
  const FourVector pnull = FourVector(threeMomentum.abs(), threeMomentum);
  const FourVector prs = pnull.LorentzBoost(ustrXcom.velocity());
  ThreeVector evec = prs.threevec() / prs.threevec().abs();
  // perform fragmentation and add particles to final_state.
  int nfrag = fragment_string(idqX1, idqX2, massX, evec, false);
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
  new_particle.set_formation_time(0.);
  final_state_.push_back(new_particle);

  NpartFinal_ = NpartString_[0] + NpartString_[1];
  return true;
}

bool StringProcess::make_final_state_2strings(
    const std::array<std::array<int, 2>, 2> &quarks,
    const std::array<FourVector, 2> &pstr_com,
    const std::array<double, 2> &m_str) {
  const std::array<FourVector, 2> ustr_com = {pstr_com[0] / m_str[0],
                                              pstr_com[1] / m_str[1]};
  for (int i = 0; i < 2; i++) {
    /* determine direction in which string i is stretched.
     * this is set to be same with the three-momentum of string
     * in the center of mass frame. */
    const ThreeVector mom = pstr_com[i].threevec();
    const FourVector pnull(mom.abs(), mom);
    const FourVector prs = pnull.LorentzBoost(ustr_com[i].velocity());
    ThreeVector evec = prs.threevec() / prs.threevec().abs();
    // perform fragmentation and add particles to final_state.
    int nfrag =
        fragment_string(quarks[i][0], quarks[i][1], m_str[i], evec, false);
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

/** double-diffractive : A + B -> X + X */
bool StringProcess::next_DDiff() {
  NpartFinal_ = 0;
  NpartString_[0] = 0;
  NpartString_[1] = 0;
  final_state_.clear();

  int ntry = 0;
  std::array<bool, 2> found_mass = {false, false};
  std::array<std::array<int, 2>, 2> quarks;
  std::array<FourVector, 2> pstr_com;
  std::array<double, 2> m_str;
  ThreeVector threeMomentum;
  while ((!found_mass[0] || !found_mass[1]) && (ntry < 100)) {
    ntry++;

    make_string_ends(PDGcodes_[0], quarks[0][0], quarks[0][1]);
    make_string_ends(PDGcodes_[1], quarks[1][0], quarks[1][1]);
    // sample the lightcone momentum fraction carried by gluons
    const double xmin_gluon_fraction = pmin_gluon_lightcone_ / sqrtsAB_;
    const double xfracA =
        Random::beta_a0(xmin_gluon_fraction, pow_fgluon_beta_ + 1.);
    const double xfracB =
        Random::beta_a0(xmin_gluon_fraction, pow_fgluon_beta_ + 1.);
    // sample the transverse momentum transfer
    const double QTrx = Random::normal(0., sigma_qperp_ / sqrt2_);
    const double QTry = Random::normal(0., sigma_qperp_ / sqrt2_);
    const double QTrn = std::sqrt(QTrx * QTrx + QTry * QTry);
    // evaluate the lightcone momentum transfer
    const double QPos = -QTrn * QTrn / (2. * xfracB * PNegB_);
    const double QNeg = QTrn * QTrn / (2. * xfracA * PPosA_);
    // compute four-momentum of string 1
    threeMomentum =
        evecBasisAB_[0] * (PPosA_ + QPos - PNegA_ - QNeg) / sqrt2_ +
        evecBasisAB_[1] * QTrx + evecBasisAB_[2] * QTry;
    pstr_com[0] = FourVector((PPosA_ + QPos + PNegA_ + QNeg) / sqrt2_,
                             threeMomentum);
    // compute four-momentum of string 2
    threeMomentum =
        evecBasisAB_[0] * (PPosB_ - QPos - PNegB_ + QNeg) / sqrt2_ -
        evecBasisAB_[1] * QTrx - evecBasisAB_[2] * QTry;
    pstr_com[1] = FourVector((PPosB_ - QPos + PNegB_ - QNeg) / sqrt2_,
                             threeMomentum);
    found_mass[0] = false;
    found_mass[1] = false;
    for (int i = 0; i < 2; i++) {
      m_str[i] = pstr_com[i].sqr();
      m_str[i] = (m_str[i] > 0.) ? std::sqrt(m_str[i]) : 0.;
      const double threshold = pythia_->particleData.m0(quarks[i][0]) +
                               pythia_->particleData.m0(quarks[i][1]);
      // string mass must be larger than threshold set by PYTHIA.
      if (m_str[i] > threshold) {
        found_mass[i] = true;
      }
    }
  }

  if (!found_mass[0] || !found_mass[1]) {
    return false;
  }
  const bool success = make_final_state_2strings(quarks, pstr_com, m_str);
  return success;
}

/** soft non-diffractive */
bool StringProcess::next_NDiffSoft() {
  NpartFinal_ = 0;
  NpartString_[0] = 0;
  NpartString_[1] = 0;
  final_state_.clear();

  int ntry = 0;
  std::array<bool, 2> found_mass = {false, false};
  std::array<std::array<int, 2>, 2> quarks;
  std::array<FourVector, 2> pstr_com;
  std::array<double, 2> m_str;
  while ((!found_mass[0] || !found_mass[1]) && (ntry < 100)) {
    ntry++;

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
      ss << "  StringProcess::next_NDiff : baryonA = " << bar_a <<
            ", baryonB = " << bar_b;
      throw std::runtime_error(ss.str());
    }
    // sample the lightcone momentum fraction carried by quarks
    const double xfracA = Random::beta(pow_fquark_alpha_, pow_fquark_beta_);
    const double xfracB = Random::beta(pow_fquark_alpha_, pow_fquark_beta_);
    // sample the transverse momentum transfer
    const double QTrx = Random::normal(0., sigma_qperp_ / sqrt2_);
    const double QTry = Random::normal(0., sigma_qperp_ / sqrt2_);
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
    pstr_com[0] = FourVector((PPosA_ + dPPos + PNegA_ + dPNeg) / sqrt2_,
                             threeMomentum);
    m_str[0] = pstr_com[0].sqr();
    // compute four-momentum of string 2
    threeMomentum =
        evecBasisAB_[0] * (PPosB_ - dPPos - PNegB_ + dPNeg) / sqrt2_ -
        evecBasisAB_[1] * QTrx - evecBasisAB_[2] * QTry;
    pstr_com[1] = FourVector((PPosB_ - dPPos + PNegB_ - dPNeg) / sqrt2_,
                             threeMomentum);

    found_mass[0] = false;
    found_mass[1] = false;
    for (int i = 0; i < 2; i++) {
      m_str[i] = pstr_com[i].sqr();
      m_str[i] = (m_str[i] > 0.) ? std::sqrt(m_str[i]) : 0.;
      const double threshold = pythia_->particleData.m0(quarks[i][0]) +
                               pythia_->particleData.m0(quarks[i][1]);
      // string mass must be larger than threshold set by PYTHIA.
      if (m_str[i] > threshold) {
        found_mass[i] = true;
      }
    }
  }

  if (!found_mass[0] || !found_mass[1]) {
    return false;
  }
  const bool success = make_final_state_2strings(quarks, pstr_com, m_str);
  return success;
}

/** baryon-antibaryon annihilation */
bool StringProcess::next_BBbarAnn() {
  const std::array<FourVector, 2> ustrcom = {FourVector(1.,0.,0.,0.),
                                             FourVector(1.,0.,0.,0.)};

  NpartFinal_ = 0;
  NpartString_[0] = 0;
  NpartString_[1] = 0;
  final_state_.clear();

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
    qcount_bar.push_back(baryon.net_quark_number(i));
    qcount_antibar.push_back(-antibaryon.net_quark_number(i));
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
      new_particle.set_formation_time(0.);
      final_state_.push_back(new_particle);
    }
    NpartFinal_ = NpartString_[0] + NpartString_[1];
    return true;
  }

  // Select qqbar pair to annihilate and remove it away
  auto discrete_distr = Random::discrete_dist<int>(n_combinations);
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

  const int max_ntry = 100;
  const std::array<double, 2> mstr = {0.5 * sqrtsAB_, 0.5 * sqrtsAB_};
  for (int ntry = 0; ntry < max_ntry; ntry++) {
    // Randomly select two quark-antiquark pairs
    if (Random::uniform_int(0, 1) == 0) {
      std::swap(remaining_quarks[0], remaining_quarks[1]);
    }
    if (Random::uniform_int(0, 1) == 0) {
      std::swap(remaining_antiquarks[0], remaining_antiquarks[1]);
    }
    // Make sure it satisfies kinematical threshold constraint
    bool kin_threshold_satisfied = true;
    for (int i = 0; i < 2; i++) {
      const double mstr_min = pythia_->particleData.m0(remaining_quarks[i]) +
                              pythia_->particleData.m0(remaining_antiquarks[i]);
      if (mstr_min > mstr[i]) {
        kin_threshold_satisfied = false;
      }
    }
    if (!kin_threshold_satisfied) {
      continue;
    }
    // Fragment two strings
    for (int i = 0; i < 2; i++) {
      ThreeVector evec = pcom_[i].threevec() / pcom_[i].threevec().abs();
      const int nfrag = fragment_string(
          remaining_quarks[i], remaining_antiquarks[i], mstr[i], evec, false);
      if (nfrag <= 0) {
        NpartString_[i] = 0;
        return false;
      }
      NpartString_[i] = append_final_state(ustrcom[i], evec);
    }
    NpartFinal_ = NpartString_[0] + NpartString_[1];
    return true;
  }

  return false;
}

void StringProcess::make_orthonormal_basis() {
  const double pabscomAB = pCM(sqrtsAB_, massA_, massB_);
  if (std::abs(pcom_[0].x3()) < (1. - 1.0e-8) * pabscomAB) {
    double ex, ey, et;
    double theta, phi;

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

    evecBasisAB_[1].set_x1(cos(theta) * cos(phi));
    evecBasisAB_[1].set_x2(cos(theta) * sin(phi));
    evecBasisAB_[1].set_x3(-sin(theta));

    evecBasisAB_[2].set_x1(-sin(phi));
    evecBasisAB_[2].set_x2(cos(phi));
    evecBasisAB_[2].set_x3(0.);
  } else {
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
  PPosA_ =
      (pcom_[0].x0() + evecBasisAB_[0] * pcom_[0].threevec()) / sqrt2_;
  PNegA_ =
      (pcom_[0].x0() - evecBasisAB_[0] * pcom_[0].threevec()) / sqrt2_;
  PPosB_ =
      (pcom_[1].x0() + evecBasisAB_[0] * pcom_[1].threevec()) / sqrt2_;
  PNegB_ =
      (pcom_[1].x0() - evecBasisAB_[0] * pcom_[1].threevec()) / sqrt2_;
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
  diquark += (q1 != q2 && Random::uniform_int(0, 3) == 0) ? 1 : 3;
  return (q1 < 0) ? -diquark : diquark;
}

void StringProcess::make_string_ends(const PdgCode &pdg, int &idq1, int &idq2) {
  std::array<int, 3> quarks = pdg.quark_content();

  if (pdg.is_meson()) {
    idq1 = quarks[1];
    idq2 = quarks[2];
    /* Some mesons with PDG id 11X are actually mixed state of uubar and ddbar.
     * have a random selection whether we have uubar or ddbar in this case. */
    if (idq1 == 1 && idq2 == -1 && Random::uniform_int(0, 1) == 0) {
      idq1 = 2;
      idq2 = -2;
    }
  } else {
    assert(pdg.is_baryon());
    // Get random quark to position 0
    std::swap(quarks[Random::uniform_int(0, 2)], quarks[0]);
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
                                   bool random_rotation) {
  pythia_->event.reset();
  // evaluate 3 times total baryon number of the string
  const int bstring = pythia_->particleData.baryonNumberType(idq1) +
                      pythia_->particleData.baryonNumberType(idq2);
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

  const double m1 = pythia_->particleData.m0(idq1);
  const double m2 = pythia_->particleData.m0(idq2);
  if (m1 + m2 > mString) {
    throw std::runtime_error("String fragmentation: m1 + m2 > mString");
  }

  // evaluate momenta of quarks
  const double pCMquark = pCM(mString, m1, m2);
  const double E1 = std::sqrt(m1 * m1 + pCMquark * pCMquark);
  const double E2 = std::sqrt(m2 * m2 + pCMquark * pCMquark);

  ThreeVector direction = sign_direction * evecLong;
  if (random_rotation) {
    Angles phitheta;
    phitheta.distribute_isotropically();
    direction = phitheta.threevec();
  } else if (Random::uniform_int(0, 1) == 0) {
    /* in the case where we flip the string ends,
     * we need to flip the longitudinal unit vector itself
     * since it is set to be direction of diquark (anti-quark) or anti-diquark.
     */
    evecLong = -evecLong;
    direction = sign_direction * evecLong;
  }

  // For status and (anti)color see \iref{Sjostrand:2007gs}.
  const int status1 = 1, color1 = 1, anticolor1 = 0;
  Pythia8::Vec4 pquark = set_Vec4(E1, -direction * pCMquark);
  pythia_->event.append(idq1, status1, color1, anticolor1, pquark, m1);

  const int status2 = 1, color2 = 0, anticolor2 = 1;
  pquark = set_Vec4(E2, direction * pCMquark);
  pythia_->event.append(idq2, status2, color2, anticolor2, pquark, m2);

  // implement PYTHIA fragmentation
  const bool successful_hadronization = pythia_->forceHadronLevel();
  int number_of_fragments = 0;
  if (successful_hadronization) {
    for (int ipart = 0; ipart < pythia_->event.size(); ipart++) {
      if (pythia_->event[ipart].isFinal()) {
        number_of_fragments++;
      }
    }
  }

  return number_of_fragments;
}

}  // namespace Smash
