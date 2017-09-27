#include "include/stringprocess.h"
#include "include/angles.h"
#include "include/kinematics.h"
#include "include/random.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace Smash {

// constructor
StringProcess::StringProcess() {
  for (int imu = 0; imu < 4; imu++) {
    evecBasisAB[imu] = ThreeVector(0., 0., 0.);
  }

  baryonA = baryonB = 0;
  //chargeA = chargeB = 0;
  PPosA = 0.;
  PNegA = 0.;
  PPosB = 0.;
  PNegB = 0.;
  massA = massB = 0.;
  sqrtsAB = 0.;
  pabscomAB = 0.;

  pmin_gluon_lightcone = 0.001;
  pow_fgluon_beta = 0.5;
  pow_fquark_alpha = 1.;
  pow_fquark_beta = 2.5;

  sigma_qperp = 0.5;
  kappa_tension_string = 1.;

  time_collision = 0.;
  gamma_factor_com = 1.;

  NpartFinal = 0;
  NpartString1 = 0;
  NpartString2 = 0;
  final_state.clear();
}

// destructor
StringProcess::~StringProcess() {
}

void StringProcess::set_pythia(Pythia8::Pythia *pythiaIn) {
  pythia = pythiaIn;
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
    bool is_leading;   // is leading hadron or not
  };

  std::vector<fragment_type> fragments;
  double p_pos_tot = 0.0, p_neg_tot = 0.0;
  int bstring = 0;

  for (int ipyth = 0; ipyth < pythia->event.size(); ipyth++) {
    if (!pythia->event[ipyth].isFinal()) {
      continue;
    }
    int pythia_id = pythia->event[ipyth].id();
    /* K_short and K_long need are converted to K0 since SMASH only knows K0 */
    if (pythia_id == 310 || pythia_id == 130) {
      pythia_id = (Random::uniform_int(0, 1) == 0) ? 311 : -311;
    }
    PdgCode pdg = PdgCode::from_decimal(pythia_id);
    if (!fragments[ipyth].pdg.is_hadron()) {
      throw std::invalid_argument("StringProcess::append_final_state warning :"
             " particle is not meson or baryon.");
    }
    FourVector mom(pythia->event[ipyth].e(),
                   pythia->event[ipyth].px(),
                   pythia->event[ipyth].py(),
                   pythia->event[ipyth].pz());
    double pparallel = mom.threevec() * evecLong;
    double y = 0.5 * std::log((mom.x0() + pparallel) / (mom.x0() - pparallel));
    fragments.push_back({mom, pparallel, y, 0.0, pdg, false});
    // total lightcone momentum
    p_pos_tot += (mom.x0() + pparallel) / std::sqrt(2.);
    p_neg_tot += (mom.x0() - pparallel) / std::sqrt(2.);
    bstring += pythia->particleData.baryonNumberType(pythia_id);
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
  xvertex_pos[0] = p_pos_tot / kappa_tension_string;
  for (int i = 0; i < nfrag; i++) {
    // recursively compute x^{+} coordinates of q-qbar formation vertex
    xvertex_pos[i + 1] = xvertex_pos[i] -
        (fragments[i].momentum.x0() + fragments[i].pparallel) /
        (kappa_tension_string * std::sqrt(2.));
  }
  // x^{-} coordinates of the backward end
  xvertex_neg[nfrag] = p_neg_tot / kappa_tension_string;
  for (int i = nfrag - 1; i >= 0; i--) {
    // recursively compute x^{-} coordinates of q-qbar formation vertex
    xvertex_neg[i] = xvertex_neg[i + 1] -
        (fragments[i].momentum.x0() - fragments[i].pparallel) /
        (kappa_tension_string * std::sqrt(2.));
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
    while (i < nfrag && fragments[i].pdg.is_meson()) { i++; }
    fragments[i].xtotfac = 2. /  3.;
    fragments[i].is_leading = true;
    // The most backward meson
    i = nfrag - 1;
    while (i >= 0 && fragments[i].pdg.is_baryon()) { i--; }
    fragments[i].xtotfac = 0.5;
    fragments[i].is_leading = true;
  } else {
    throw std::invalid_argument("StringProcess::append_final_state"
        " encountered bstring != 0, 3, -3");
  }

  // Velocity three-vector to perform Lorentz boost.
  const ThreeVector vstring = uString.velocity();

  /* compute the formation times of hadrons
   * from the lightcone coordinates of q-qbar formation vertices. */
  for (int i = 0; i < nfrag; i++) {
    // set the formation time and position in the rest frame of string
    double t_prod = (xvertex_pos[i] + xvertex_neg[i + 1]) / std::sqrt(2.);
    FourVector fragment_position = FourVector(t_prod, evecLong *
            (xvertex_pos[i] - xvertex_neg[i + 1]) / std::sqrt(2.));
    // boost into the lab frame
    fragments[i].momentum = fragments[i].momentum.LorentzBoost(  -vstring );
    fragment_position = fragment_position.LorentzBoost(-vstring );
    t_prod = fragment_position.x0();

    ParticleData new_particle(ParticleType::find(fragments[i].pdg));
    new_particle.set_4momentum(fragments[i].momentum);

    constexpr double suppression_factor = 0.7;
    new_particle.set_cross_section_scaling_factor(fragments[i].is_leading ?
        suppression_factor * fragments[i].xtotfac : 0.);
    new_particle.set_formation_time(time_collision +
                                    gamma_factor_com * t_prod);
    final_state.push_back(new_particle);
  }

  return nfrag;
}

bool StringProcess::init(const ParticleList &incomingList,
                         double tcollIn, double gammaFacIn){
  PDGcodeA = incomingList[0].pdgcode();
  PDGcodeB = incomingList[1].pdgcode();
  massA = incomingList[0].effective_mass();
  massB = incomingList[1].effective_mass();

  plabA = incomingList[0].momentum();
  plabB = incomingList[1].momentum();

  sqrtsAB = ( plabA + plabB ).abs();
  pabscomAB = pCM(sqrtsAB, massA, massB);

  make_incoming_com_momenta();
  make_orthonormal_basis();
  compute_incoming_lightcone_momenta();

  xmin_gluon_fraction = pmin_gluon_lightcone / sqrtsAB;

  // quantum numbers of hadron A
  baryonA = 3*PDGcodeA.baryon_number();
  //chargeA = 3*PDGcodeA.charge();
  // quantum numbers of hadron B
  baryonB = 3*PDGcodeB.baryon_number();
  //chargeB = 3*PDGcodeB.charge();

  time_collision = tcollIn;
  gamma_factor_com = gammaFacIn;

  return true;
}

/**
 * single diffractive
 * channel = 1 : A + B -> A + X
 * channel = 2 : A + B -> X + B
 */
bool StringProcess::next_SDiff(int channel) {
  bool ret;

  int ntry;
  bool foundPabsX, foundMassX;

  int pdgidH;
  double mstrMin = 0.;
  double mstrMax = 0.;
  double pabscomHX, massH, massX, rmass;
  double QTrn, QTrx, QTry;

  int nfrag;
  int idqX1, idqX2;

  FourVector pstrHcom;
  FourVector pstrHlab;
  FourVector pstrXcom;
  FourVector pstrXlab;
  ThreeVector threeMomentum;

  FourVector ustrHcom;
  FourVector ustrHlab;
  FourVector ustrXcom;
  FourVector ustrXlab;

  double pabs;
  FourVector pnull;
  FourVector prs;
  ThreeVector evec;

  NpartFinal = 0;
  NpartString1 = 0;
  NpartString2 = 0;
  final_state.clear();

  ntry = 0;
  foundPabsX = false;
  foundMassX = false;

  if( channel == 1 ) { // AB > AX
    mstrMin = massB;
    mstrMax = sqrtsAB - massA;
    pdgidH = PDGcodeA.get_decimal();
    massH = massA;
  } else if( channel == 2 ) { // AB > XB
    mstrMin = massA;
    mstrMax = sqrtsAB - massB;
    pdgidH = PDGcodeB.get_decimal();
    massH = massB;
  } else {
    throw std::runtime_error("invalid argument for StringProcess::next_SDiff");
  }

  while (((foundPabsX == false) || (foundMassX == false)) && (ntry < 100)) {
    ntry = ntry + 1;
    // decompose hadron into quarks
    if( channel == 1 ) { // AB > AX
      make_string_ends(PDGcodeB, idqX1, idqX2);
    } else if( channel == 2 ) { // AB > XB
      make_string_ends(PDGcodeA, idqX1, idqX2);
    }
    // sample the transverse momentum transfer
    QTrx = Random::normal(0., sigma_qperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigma_qperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);
    // sample the string mass and evaluate the three-momenta of hadron and string.
    rmass = std::log(mstrMax / mstrMin) * Random::uniform(0., 1.);
    massX = mstrMin * exp(rmass);
    pabscomHX = pCM( sqrtsAB, massH, massX );
    // magnitude of the three momentum must be larger than the transvers momentum.
    foundPabsX = pabscomHX > QTrn;
    // string mass must be larger than threshold set by PYTHIA.
    foundMassX = massX > (pythia->particleData.m0(idqX1) +
                          pythia->particleData.m0(idqX2));
  }

  ret = false;
  if ((foundPabsX == true) && (foundMassX == true)) {
    double sign_direction = 0.;
    if( channel == 1 ) { // AB > AX
      sign_direction = 1.;
    } else if( channel == 2 ) { // AB > XB
      sign_direction = -1.;
    }
    threeMomentum = sign_direction * (
                        evecBasisAB[3] * std::sqrt(pabscomHX*pabscomHX - QTrn*QTrn) +
                        evecBasisAB[1] * QTrx +
                        evecBasisAB[2] * QTry );
    pstrHcom = FourVector( std::sqrt(pabscomHX*pabscomHX + massH*massH), threeMomentum );
    threeMomentum = -sign_direction * (
                        evecBasisAB[3] * std::sqrt(pabscomHX*pabscomHX - QTrn*QTrn) +
                        evecBasisAB[1] * QTrx +
                        evecBasisAB[2] * QTry );
    pstrXcom = FourVector( std::sqrt(pabscomHX*pabscomHX + massX*massX), threeMomentum );

    pstrHlab = pstrHcom.LorentzBoost( -vcomAB );
    pstrXlab = pstrXcom.LorentzBoost( -vcomAB );

    ustrHcom = pstrHcom / massH;
    ustrXcom = pstrXcom / massX;
    ustrHlab = pstrHlab / massH;
    ustrXlab = pstrXlab / massX;
    /* determin direction in which the string is stretched.
     * this is set to be same with the three-momentum of string
     * in the center of mass frame. */
    threeMomentum = pstrXcom.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustrXcom.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;
    // perform fragmentation and add particles to final_state.
    nfrag = fragment_string(idqX1, idqX2, massX, evec, false);
    if (nfrag > 0) {
      NpartString1 = append_final_state(ustrXlab, evec);
    } else {
      nfrag = 0;
      NpartString1 = 0;
      ret = false;
    }

    NpartString2 = 1;
    const std::string s = std::to_string(pdgidH);
    PdgCode hadron_code(s);
    ParticleData new_particle(ParticleType::find(hadron_code));
    new_particle.set_4momentum(pstrHlab);
    new_particle.set_cross_section_scaling_factor(1.);
    new_particle.set_formation_time(0.);
    final_state.push_back(new_particle);

    if ((NpartString1 > 0) && (NpartString2 > 0) && (nfrag == NpartString1)) {
      NpartFinal = NpartString1 + NpartString2;
      ret = true;
    }
  }

  return ret;
}

/** double-diffractive : A + B -> X + X */
bool StringProcess::next_DDiff() {
  bool ret;

  int ntry;
  bool foundMass1, foundMass2;

  double xfracA, xfracB;
  double QPos, QNeg;
  double QTrn, QTrx, QTry;

  int nfrag1, nfrag2;
  int idq11, idq12;
  int idq21, idq22;
  double mstr1, mstr2;

  FourVector pstr1com;
  FourVector pstr1lab;
  FourVector pstr2com;
  FourVector pstr2lab;
  ThreeVector threeMomentum;

  FourVector ustr1com;
  FourVector ustr1lab;
  FourVector ustr2com;
  FourVector ustr2lab;

  double pabs;
  FourVector pnull;
  FourVector prs;
  ThreeVector evec;

  NpartFinal = 0;
  NpartString1 = 0;
  NpartString2 = 0;
  final_state.clear();

  ntry = 0;
  foundMass1 = false;
  foundMass2 = false;
  while (((foundMass1 == false) || (foundMass2 == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    make_string_ends(PDGcodeA, idq11, idq12);
    make_string_ends(PDGcodeB, idq21, idq22);
    // sample the lightcone momentum fraction carried by gluons
    xfracA = Random::beta_a0(xmin_gluon_fraction, pow_fgluon_beta + 1.);
    xfracB = Random::beta_a0(xmin_gluon_fraction, pow_fgluon_beta + 1.);
    // sample the transverse momentum transfer
    QTrx = Random::normal(0., sigma_qperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigma_qperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);
    // evaluate the lightcone momentum transfer
    QPos = -QTrn*QTrn / (2. * xfracB * PNegB);
    QNeg = QTrn*QTrn / (2. * xfracA * PPosA);
    // compute four-momentum of string 1
    threeMomentum = evecBasisAB[3] * (PPosA + QPos - PNegA - QNeg) / std::sqrt(2.) +
                        evecBasisAB[1] * QTrx + evecBasisAB[2] * QTry;
    pstr1com = FourVector( (PPosA + QPos + PNegA + QNeg) / std::sqrt(2.), threeMomentum );
    mstr1 = pstr1com.sqr();
    // compute four-momentum of string 2
    threeMomentum = evecBasisAB[3] * (PPosB - QPos - PNegB + QNeg) / std::sqrt(2.) -
                        evecBasisAB[1] * QTrx - evecBasisAB[2] * QTry;
    pstr2com = FourVector( (PPosB - QPos + PNegB - QNeg) / std::sqrt(2.), threeMomentum );
    mstr2 = pstr2com.sqr();
    // string mass must be larger than threshold set by PYTHIA.
    mstr1 = (mstr1 > 0.) ? std::sqrt(mstr1) : 0.;
    foundMass1 = mstr1 > (pythia->particleData.m0(idq11) +
                          pythia->particleData.m0(idq12));
    mstr2 = (mstr2 > 0.) ? std::sqrt(mstr2) : 0.;
    foundMass2 = mstr2 > (pythia->particleData.m0(idq21) +
                          pythia->particleData.m0(idq22));
  }

  ret = false;
  bool both_masses_above_pythia_threshold = foundMass1 && foundMass2;
  if ( both_masses_above_pythia_threshold ) {
    pstr1lab = pstr1com.LorentzBoost( -vcomAB );
    pstr2lab = pstr2com.LorentzBoost( -vcomAB );

    ustr1com = pstr1com / mstr1;
    ustr2com = pstr2com / mstr2;
    ustr1lab = pstr1lab / mstr1;
    ustr2lab = pstr2lab / mstr2;
    /* determin direction in which string 1 is stretched.
     * this is set to be same with the three-momentum of string
     * in the center of mass frame. */
    threeMomentum = pstr1com.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustr1com.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;
    // perform fragmentation and add particles to final_state.
    nfrag1 = fragment_string(idq11, idq12, mstr1, evec, false);
    if (nfrag1 > 0) {
      NpartString1 = append_final_state(ustr1lab, evec);
    } else {
      nfrag1 = 0;
      NpartString1 = 0;
      ret = false;
    }
    /* determin direction in which string 2 is stretched.
     * this is set to be same with the three-momentum of string
     * in the center of mass frame. */
    threeMomentum = pstr2com.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustr2com.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;
    // perform fragmentation and add particles to final_state.
    nfrag2 = fragment_string(idq21, idq22, mstr2, evec, false);
    if (nfrag2 > 0) {
      NpartString2 = append_final_state(ustr2lab, evec);
    } else {
      nfrag2 = 0;
      NpartString2 = 0;
      ret = false;
    }

    if ((NpartString1 > 0) && (NpartString2 > 0) && (nfrag1 == NpartString1) &&
        (nfrag2 == NpartString2)) {
      NpartFinal = NpartString1 + NpartString2;
      ret = true;
    }
  }

  return ret;
}

/** non-diffractive */
bool StringProcess::next_NDiff() {
  bool ret;

  int ntry;
  bool foundMass1, foundMass2;

  double xfracA, xfracB;
  double QPos, QNeg;
  double dPPos, dPNeg;
  double QTrn, QTrx, QTry;

  int nfrag1, nfrag2;
  int idqA1, idqA2;
  int idqB1, idqB2;
  int idq11, idq12;
  int idq21, idq22;
  double mstr1, mstr2;

  FourVector pstr1com;
  FourVector pstr1lab;
  FourVector pstr2com;
  FourVector pstr2lab;
  ThreeVector threeMomentum;

  FourVector ustr1com;
  FourVector ustr1lab;
  FourVector ustr2com;
  FourVector ustr2lab;

  double pabs;
  FourVector pnull;
  FourVector prs;
  ThreeVector evec;

  NpartFinal = 0;
  NpartString1 = 0;
  NpartString2 = 0;
  final_state.clear();

  ntry = 0;
  foundMass1 = false;
  foundMass2 = false;
  while (((foundMass1 == false) || (foundMass2 == false)) && (ntry < 100)) {
    ntry = ntry + 1;

    make_string_ends(PDGcodeA, idqA1, idqA2);
    make_string_ends(PDGcodeB, idqB1, idqB2);

    if ((baryonA == 3) && (baryonB == 3)) {  // baryon-baryon
      idq11 = idqB1;
      idq12 = idqA2;
      idq21 = idqA1;
      idq22 = idqB2;
    } else if ((baryonA == 3) && (baryonB == 0)) {  // baryon-meson
      idq11 = idqB1;
      idq12 = idqA2;
      idq21 = idqA1;
      idq22 = idqB2;
    } else if ((baryonA == 3) && (baryonB == -3)) {  // baryon-antibaryon
      idq11 = idqB1;
      idq12 = idqA2;
      idq21 = idqA1;
      idq22 = idqB2;
    } else if ((baryonA == 0) && (baryonB == 3)) {  // meson-baryon
      idq11 = idqB1;
      idq12 = idqA2;
      idq21 = idqA1;
      idq22 = idqB2;
    } else if ((baryonA == 0) && (baryonB == 0)) {  // meson-meson
      idq11 = idqB1;
      idq12 = idqA2;
      idq21 = idqA1;
      idq22 = idqB2;
    } else if ((baryonA == 0) && (baryonB == -3)) {  // meson-antibaryon
      idq11 = idqA1;
      idq12 = idqB2;
      idq21 = idqB1;
      idq22 = idqA2;
    } else if ((baryonA == -3) && (baryonB == 3)) {  // antibaryon-baryon
      idq11 = idqA1;
      idq12 = idqB2;
      idq21 = idqB1;
      idq22 = idqA2;
    } else if ((baryonA == -3) && (baryonB == 0)) {  // antibaryon-meson
      idq11 = idqA1;
      idq12 = idqB2;
      idq21 = idqB1;
      idq22 = idqA2;
    } else if ((baryonA == -3) && (baryonB == -3)) {  // antibaryon-antibaryon
      idq11 = idqA1;
      idq12 = idqB2;
      idq21 = idqB1;
      idq22 = idqA2;
    } else {
      fprintf(stderr,
              "  StringProcess::next_NDiff : incorrect baryon number of incoming "
              "hadrons.\n");
      fprintf(stderr, "  StringProcess::next_NDiff : baryonA = %d, baryonB = %d\n",
              baryonA, baryonB);
      exit(1);
    }
    // sample the lightcone momentum fraction carried by quarks
    xfracA = Random::beta(pow_fquark_alpha, pow_fquark_beta);
    xfracB = Random::beta(pow_fquark_alpha, pow_fquark_beta);
    // sample the transverse momentum transfer
    QTrx = Random::normal(0., sigma_qperp/std::sqrt(2.) );
    QTry = Random::normal(0., sigma_qperp/std::sqrt(2.) );
    QTrn = std::sqrt(QTrx*QTrx + QTry*QTry);
    // evaluate the lightcone momentum transfer
    QPos = -QTrn*QTrn / (2. * xfracB * PNegB);
    QNeg = QTrn*QTrn / (2. * xfracA * PPosA);
    dPPos = -xfracA * PPosA - QPos;
    dPNeg = xfracB * PNegB - QNeg;
    // compute four-momentum of string 1
    threeMomentum = evecBasisAB[3] * (PPosA + dPPos - PNegA - dPNeg) / std::sqrt(2.) +
                        evecBasisAB[1] * QTrx + evecBasisAB[2] * QTry;
    pstr1com = FourVector( (PPosA + dPPos + PNegA + dPNeg) / std::sqrt(2.), threeMomentum );
    mstr1 = pstr1com.sqr();
    // compute four-momentum of string 2
    threeMomentum = evecBasisAB[3] * (PPosB - dPPos - PNegB + dPNeg) / std::sqrt(2.) -
                        evecBasisAB[1] * QTrx - evecBasisAB[2] * QTry;
    pstr2com = FourVector( (PPosB - dPPos + PNegB - dPNeg) / std::sqrt(2.), threeMomentum );
    mstr2 = pstr2com.sqr();
    // string mass must be larger than threshold set by PYTHIA.
    mstr1 = (mstr1 > 0.) ? std::sqrt(mstr1) : 0.;
    foundMass1 = mstr1 > (pythia->particleData.m0(idq11) +
                          pythia->particleData.m0(idq12));
    mstr2 = (mstr2 > 0.) ? std::sqrt(mstr2) : 0.;
    foundMass2 = mstr2 > (pythia->particleData.m0(idq21) +
                          pythia->particleData.m0(idq22));
  }

  ret = false;
  bool both_masses_above_pythia_threshold = foundMass1 && foundMass2;
  if ( both_masses_above_pythia_threshold ) {
    pstr1lab = pstr1com.LorentzBoost( -vcomAB );
    pstr2lab = pstr2com.LorentzBoost( -vcomAB );

    ustr1com = pstr1com / mstr1;
    ustr2com = pstr2com / mstr2;
    ustr1lab = pstr1lab / mstr1;
    ustr2lab = pstr2lab / mstr2;
    /* determin direction in which string 1 is stretched.
     * this is set to be same with the three-momentum of string
     * in the center of mass frame. */
    threeMomentum = pstr1com.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustr1com.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;
    // perform fragmentation and add particles to final_state.
    nfrag1 = fragment_string(idq11, idq12, mstr1, evec, false);
    if (nfrag1 > 0) {
      NpartString1 = append_final_state(ustr1lab, evec);
    } else {
      nfrag1 = 0;
      NpartString1 = 0;
      ret = false;
    }
    /* determin direction in which string 2 is stretched.
     * this is set to be same with the three-momentum of string
     * in the center of mass frame. */
    threeMomentum = pstr2com.threevec();
    pnull = FourVector( threeMomentum.abs(), threeMomentum );
    prs = pnull.LorentzBoost( ustr2com.velocity() );
    pabs = prs.threevec().abs();
    evec = prs.threevec() / pabs;
    // perform fragmentation and add particles to final_state.
    nfrag2 = fragment_string(idq21, idq22, mstr2, evec, false);
    if (nfrag2 > 0) {
      NpartString2 = append_final_state(ustr2lab, evec);
    } else {
      nfrag2 = 0;
      NpartString2 = 0;
      ret = false;
    }

    if ((NpartString1 > 0) && (NpartString2 > 0) && (nfrag1 == NpartString1) &&
        (nfrag2 == NpartString2)) {
      NpartFinal = NpartString1 + NpartString2;
      ret = true;
    }
  }

  return ret;
}

/** baryon-antibaryon annihilation */
bool StringProcess::next_BBbarAnn(){
	bool ret;

	int ntry;
	bool isBBbarpair, isAnnihilating;

	int ic, jc;
	int ijc, ipr, npr;
	std::array<int,3> quark_content_A;
	std::array<int,3> quark_content_B;
	std::vector<int> indexAnn;

	int idq11, idq12, idq12prev;
	int idq21, idq22, idq22prev;
	int nfrag1, nfrag2;
	double mstr1, mstr2;
	double mstr1Min, mstr2Min;

	FourVector ustr1lab;
	FourVector ustr2lab;

	double pabs;
	ThreeVector evec;

	ustr1lab = ucomAB;
	ustr2lab = ucomAB;

	indexAnn.resize(0);

	NpartFinal = 0;
	NpartString1 = 0;
	NpartString2 = 0;
	final_state.clear();

	quark_content_A = PDGcodeA.quark_content();
	quark_content_B = PDGcodeB.quark_content();

	isBBbarpair = ( (baryonA == 3) && (baryonB == -3) ) || ( (baryonA == -3) && (baryonB == 3) );
	isAnnihilating = false;

	// if it is baryon-antibaryon pair
	if( isBBbarpair == true ){ // if it is
		mstr1 = 0.5*sqrtsAB;
		mstr2 = 0.5*sqrtsAB;

		for(ic = 0; ic < 3; ic++){
			for(jc = 0; jc < 3; jc++){
				if( quark_content_A[ic] == -quark_content_B[jc] ){
					ijc = ic*10 + jc;
					indexAnn.push_back( ijc );
				}
			}
		}
		npr = indexAnn.size();
		fprintf(stderr,"  StringProcess::next_BBarAnn : %d possible pairs for qqbar annihilation\n", npr);
		/* if it is a BBbar pair but there is no qqbar pair to annihilate,
		 * nothing happens */
		if( npr == 0 ){
			NpartString1 = 1;
			ParticleData new_particle1(ParticleType::find(PDGcodeA));
			new_particle1.set_4momentum(plabA);
			new_particle1.set_cross_section_scaling_factor(1.);
			new_particle1.set_formation_time(0.);
			final_state.push_back(new_particle1);

			NpartString2 = 1;
			ParticleData new_particle2(ParticleType::find(PDGcodeB));
			new_particle2.set_4momentum(plabB);
			new_particle2.set_cross_section_scaling_factor(1.);
			new_particle2.set_formation_time(0.);
			final_state.push_back(new_particle2);

			isAnnihilating = false;

			NpartFinal = NpartString1 + NpartString2;
			ret = true;
		}// endif no qqbar pair to annihilate

		ntry = 0;
		while( ( npr > 0 ) && ( isAnnihilating == false ) && ( ntry < 100 ) ){
			ntry = ntry + 1;

			// randomly choose a qqbar pair to annihilate
			ipr = Random::uniform_int(0, npr - 1);
			ijc = indexAnn.at(ipr);
			ic = ( ijc - (ijc%10) )/10;
			jc = ijc%10;
			fprintf(stderr,"  StringProcess::next_BBarAnn : ic = %d, jc = %d chosen\n", ic, jc);
			// make two qqbar pairs to excite strings
			idq11 = 0;
			idq12 = 0;
			idq21 = 0;
			idq22 = 0;
			if( (baryonA == 3) && (baryonB == -3) ){
				idq11 = quark_content_A[(ic + 1)%3];
				idq12 = quark_content_B[(jc + 1)%3];
				idq21 = quark_content_A[(ic + 2)%3];
				idq22 = quark_content_B[(jc + 2)%3];
			}
			else if( (baryonA == -3) && (baryonB == 3) ){
				idq11 = quark_content_B[(ic + 1)%3];
				idq12 = quark_content_A[(jc + 1)%3];
				idq21 = quark_content_B[(ic + 2)%3];
				idq22 = quark_content_A[(jc + 2)%3];
			}
			// randomly choose if we flip the antiquark contents
			if( Random::uniform_int(0, 1) == 0 ){
				idq12prev = idq12;
				idq22prev = idq22;
				idq12 = idq22prev;
				idq22 = idq12prev;
			}
			fprintf(stderr,"  StringProcess::next_BBarAnn : string 1 with %d, %d\n", idq11, idq12);
			fprintf(stderr,"  StringProcess::next_BBarAnn : string 2 with %d, %d\n", idq21, idq22);

			mstr1Min = pythia->particleData.m0(idq11) + pythia->particleData.m0(idq12);
			mstr2Min = pythia->particleData.m0(idq21) + pythia->particleData.m0(idq22);
			isAnnihilating = ( mstr1 > mstr1Min ) && ( mstr2 > mstr2Min );
		}
	}
	else{ // if it is not
		fprintf(stderr,"  StringProcess::next_BBarAnn failure : it is not BBbar pair.\n");
		isAnnihilating = false;
	}
	// endif baryon-antibaryon pair

	// implement collision in the case of annihilating BBbar pair
	if( isAnnihilating == true ){
		ret = false;

		// string 1
		pabs = pcomA.threevec().abs();
		evec = pcomA.threevec() / pabs;

		nfrag1 = fragment_string(idq11, idq12, mstr1, evec, false);
		if( nfrag1 > 0 ){
			NpartString1 = append_final_state(ustr1lab, evec);
		}
		else{
			nfrag1 = 0;
			NpartString1 = 0;
			ret = false;
		}

		// string 2
		pabs = pcomB.threevec().abs();
		evec = pcomB.threevec() / pabs;

		nfrag2 = fragment_string(idq21, idq22, mstr2, evec, false);
		if( nfrag2 > 0 ){
			NpartString2 = append_final_state(ustr2lab, evec);
		}
		else{
			nfrag2 = 0;
			NpartString2 = 0;
			ret = false;
		}

		if( ( NpartString1 > 0 ) && ( NpartString2 > 0 )
			&& ( nfrag1 == NpartString1 ) && ( nfrag2 == NpartString2 ) ){
			NpartFinal = NpartString1 + NpartString2;
			ret = true;
		}
	}

	return ret;
}

void StringProcess::make_incoming_com_momenta(){
  ucomAB = ( plabA + plabB )/sqrtsAB;
  vcomAB = ucomAB.velocity();

  pcomA = plabA.LorentzBoost(vcomAB);
  pcomB = plabB.LorentzBoost(vcomAB);
}

void StringProcess::make_orthonormal_basis(){
  if (std::abs(pcomA.x3()) < (1. - 1.0e-8) * pabscomAB) {
    double ex, ey, et;
    double theta, phi;

    evecBasisAB[3] = pcomA.threevec() / pabscomAB;

    theta = std::acos(evecBasisAB[3].x3());

    ex = evecBasisAB[3].x1();
    ey = evecBasisAB[3].x2();
    et = std::sqrt(ex*ex + ey*ey);
    if (ey > 0.) {
      phi = std::acos(ex / et);
    } else {
      phi = -std::acos(ex / et);
    }

    evecBasisAB[1].set_x1( cos(theta) * cos(phi) );
    evecBasisAB[1].set_x2( cos(theta) * sin(phi) );
    evecBasisAB[1].set_x3( -sin(theta) );

    evecBasisAB[2].set_x1( -sin(phi) );
    evecBasisAB[2].set_x2( cos(phi) );
    evecBasisAB[2].set_x3( 0. );
  } else {
    if (pcomA.x3() > 0.) {
      evecBasisAB[1] = ThreeVector(1., 0., 0.);
      evecBasisAB[2] = ThreeVector(0., 1., 0.);
      evecBasisAB[3] = ThreeVector(0., 0., 1.);
    } else {
      evecBasisAB[1] = ThreeVector(0., 1., 0.);
      evecBasisAB[2] = ThreeVector(1., 0., 0.);
      evecBasisAB[3] = ThreeVector(0., 0., -1.);
    }
  }
}

void StringProcess::compute_incoming_lightcone_momenta(){
  PPosA = ( pcomA.x0() + evecBasisAB[3] * pcomA.threevec() ) / std::sqrt(2.);
  PNegA = ( pcomA.x0() - evecBasisAB[3] * pcomA.threevec() ) / std::sqrt(2.);
  PPosB = ( pcomB.x0() + evecBasisAB[3] * pcomB.threevec() ) / std::sqrt(2.);
  PNegB = ( pcomB.x0() - evecBasisAB[3] * pcomB.threevec() ) / std::sqrt(2.);
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

void StringProcess::make_string_ends(const PdgCode &pdg,
                                     int &idq1, int &idq2) {
  std::array<int, 3> quarks = pdg.quark_content();

  if (pdg.is_meson()) {
    idq1 = quarks[1];
    idq2 = quarks[2];
    /* Some mesons with PDG id 11X are actually mixed state of uubar and ddbar.
     * have a random selection whether we have uubar or ddbar in this case. */
    if (idq1 == 1 && idq2 == -1 && Random::uniform_int(0, 1) == 0) {
      idq1 =  2;
      idq2 = -2;
    }
  } else {
    assert(pdg.is_baryon());
    // Get random quark to position 0
    std::swap(quarks[Random::uniform_int(0, 2)], quarks[0]);
    idq1 = quarks[0];
    idq2 = diquark_from_quarks(quarks[1], quarks[2]);
  }
}

int StringProcess::fragment_string(int idq1, int idq2, double mString,
                            ThreeVector &evecLong, bool random_rotation) {
  pythia->event.reset();
  // evaluate 3 times total baryon number of the string
  const int bstring = pythia->particleData.baryonNumberType(idq1) +
                      pythia->particleData.baryonNumberType(idq2);
  /* diquark (anti-quark) with PDG id idq2 is going in the direction of evecLong.
   * quark with PDG id idq1 is going in the direction opposite to evecLong. */
  double sign_direction = 1.;
  if (bstring == -3) {  // anti-baryonic string
    /* anti-diquark with PDG id idq1 is going in the direction of evecLong.
     * anti-quark with PDG id idq2 is going in the direction
     * opposite to evecLong. */
    sign_direction = -1;
  }

  const double m1 = pythia->particleData.m0(idq1);
  const double m2 = pythia->particleData.m0(idq2);
  if (m1 + m2 > mString) {
    throw std::runtime_error("String fragmentation: m1 + m2 > mString");
  }

  // evaluate momenta of quarks
  const double pCMquark = pCM(mString, m1, m2);
  const double E1 = std::sqrt(m1*m1 + pCMquark*pCMquark);
  const double E2 = std::sqrt(m2*m2 + pCMquark*pCMquark);

  ThreeVector direction;
  if (random_rotation) {
    Angles phitheta;
    phitheta.distribute_isotropically();
    direction = phitheta.threevec();
  } else if (Random::uniform_int(0, 1) == 0) {
    /* in the case where we flip the string ends,
     * we need to flip the longitudinal unit vector itself
     * since it is set to be direction of diquark (anti-quark) or anti-diquark. */
    evecLong = -evecLong;
    direction = sign_direction * evecLong;
  }

  // For status and (anti)color see \iref{Sjostrand:2007gs}.
  const int status1 = 1, color1 = 1, anticolor1 = 0;
  Pythia8::Vec4 pquark = set_Vec4(E1, -direction * pCMquark);
  pythia->event.append(idq1, status1, color1, anticolor1, pquark, m1);

  const int status2 = 1, color2 = 0, anticolor2 = 1;
  pquark = set_Vec4(E2, direction * pCMquark);
  pythia->event.append(idq2, status2, color2, anticolor2, pquark, m2);

  // implement PYTHIA fragmentation
  const bool successful_hadronization = pythia->forceHadronLevel();
  int number_of_fragments = 0;
  if (successful_hadronization) {
    for (int ipart = 0; ipart < pythia->event.size(); ipart++) {
      if (pythia->event[ipart].isFinal()) {
        number_of_fragments++;
      }
    }
  }

  return number_of_fragments;
}

}  // namespace Smash
