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
#include <cmath>
#include <limits>
#include <string>

#include "smash/configuration.h"
#include "smash/forwarddeclarations.h"
#include "smash/input_keys.h"
#include "smash/kinematics.h"
#include "smash/pow.h"
#include "smash/processbranch.h"
#include "smash/random.h"

namespace smash {

bool SmashFragHook::doChangeFragPar(
    [[maybe_unused]] Pythia8::StringFlav* flavSelPtr, Pythia8::StringZ* zSelPtr,
    [[maybe_unused]] Pythia8::StringPT* pTSelPtr, [[maybe_unused]] int idEnd,
    [[maybe_unused]] double m2Had, std::vector<int> iParton,
    const Pythia8::StringEnd* nowEnd) {
  return true;
}

bool SmashFragHook::doVetoFragmentation(const Pythia8::Particle had,
                                        const Pythia8::StringEnd* sEnd) {
  return false;
}

StringProcess::StringProcess(Configuration& config)
    : pmin_gluon_lightcone_(
          config.take(InputKeys::collTerm_stringParam_gluonPMin)),
      pow_fgluon_beta_(config.take(InputKeys::collTerm_stringParam_gluonBeta)),
      pow_fquark_alpha_(
          config.take(InputKeys::collTerm_stringParam_quarkAlpha)),
      pow_fquark_beta_(config.take(InputKeys::collTerm_stringParam_quarkBeta)),
      sigma_qperp_(config.take(InputKeys::collTerm_stringParam_sigmaPerp)),
      stringz_a_leading_(
          config.take(InputKeys::collTerm_stringParam_stringZALeading)),
      stringz_b_leading_(
          config.take(InputKeys::collTerm_stringParam_stringZBLeading)),
      stringz_a_produce_(config.take(InputKeys::collTerm_stringParam_stringZA)),
      stringz_b_produce_(config.take(InputKeys::collTerm_stringParam_stringZB)),
      strange_supp_(
          config.take(InputKeys::collTerm_stringParam_strangeSuppression)),
      diquark_supp_(
          config.take(InputKeys::collTerm_stringParam_diquarkSuppression)),
      popcorn_rate_(config.take(InputKeys::collTerm_stringParam_popcornRate)),
      damp_popcorn_(config.take(InputKeys::collTerm_stringParam_dampPopcorn)),
      string_sigma_T_(
          config.take(InputKeys::collTerm_stringParam_stringSigmaT)),
      kappa_tension_string_(
          config.take(InputKeys::collTerm_stringParam_stringTension)),
      time_formation_const_(
          config.take(InputKeys::collTerm_stringParam_formationTime)),
      soft_t_form_(config.take(InputKeys::collTerm_stringParam_formTimeFactor)),
      mass_dependent_formation_times_(config.take(
          InputKeys::collTerm_stringParam_mDependentFormationTimes)),
      prob_proton_to_d_uu_(
          config.take(InputKeys::collTerm_stringParam_probabilityPToDUU)),
      separate_fragment_baryon_(
          config.take(InputKeys::collTerm_stringParam_separateFragmentBaryon)),
      use_monash_tune_(
          config.take(InputKeys::collTerm_stringParam_useMonashTune, false)),
      additional_xsec_supp_(
          config.take(InputKeys::collTerm_stringParam_unformedXsecSuppression)),
      prob_sq_to_qq_(

          config.take(InputKeys::collTerm_stringParam_probSQtoQQ)),
      popcorn_spair_(

          config.take(InputKeys::collTerm_stringParam_popcornSpair)),
      prob_qq1_to_qq0_(

          config.take(InputKeys::collTerm_stringParam_probQQ1toQQ0)),

      popcorn_smeson_(
          config.take(InputKeys::collTerm_stringParam_popcornSmeson)) {
  // setup and initialize pythia for fragmentation
  pythia_hadron_ = std::make_unique<Pythia8::Pythia>(PYTHIA_XML_DIR, false);
  /* turn off all parton-level processes to implement only hadronization */
  pythia_hadron_->readString("ProcessLevel:all = off");
  common_setup_pythia(pythia_hadron_.get(), strange_supp_, diquark_supp_,
                      popcorn_rate_, stringz_a_produce_, stringz_b_produce_,
                      string_sigma_T_);

  if (separate_fragment_baryon_) {
    frag_hook = std::make_shared<SmashFragHook>(
        pythia_hadron_.get(),
        static_cast<int>(StringProcess::LeadingStatus::LEADING_PARTON),
        stringz_a_produce_, stringz_b_produce_, stringz_a_leading_,
        stringz_b_leading_, popcorn_rate_);
    pythia_hadron_->addUserHooksPtr(frag_hook);
  }

  if (config.take(InputKeys::collTerm_stringParam_fragmentationModel) ==
      PythiaFragmentationModel::Thermal) {
    pythia_hadron_->readString("Fragmentation:model = 1");

    pythia_hadron_->readString(
        "StringPT:temperature = " +
        std::to_string(
            config.take(InputKeys::collTerm_stringParam_thermalTemperature)));

    pythia_hadron_->readString(
        "StringPT:tempPreFactor = " +
        std::to_string(
            config.take(InputKeys::collTerm_stringParam_thermalTempPrefactor)));

    pythia_hadron_->readString(
        "StringFlav:BtoMratio = " +
        std::to_string(config.take(
            InputKeys::collTerm_stringParam_thermalBaryonToMesonRatio)));

    pythia_hadron_->readString(
        "StringFlav:StrangeSuppression = " +
        std::to_string(config.take(
            InputKeys::collTerm_stringParam_thermalStrangeSuppression)));

    pythia_hadron_->readString(
        "StringFlav:nQuark = " +
        std::to_string(
            config.take(InputKeys::collTerm_stringParam_thermalNQuark)));

    pythia_hadron_->readString(
        std::string("StringFlav:mesonNonetL1 = ") +
        (config.take(InputKeys::collTerm_stringParam_thermalMesonNonetL1)
             ? "on"
             : "off"));

    logg[LPythia].warn(
        "Using Pythia thermal fragmentation model.\n"
        "Please cite: https://inspirehep.net/literature/1495274\n"
        "\n"
        "This model uses thermal fragmentation settings that differ\n"
        "partially from the standard tuning employed in SMASH.\n"
        "\n"
        "While the model itself is implemented in PYTHIA, its behaviour\n"
        "within the SMASH framework has not yet been systematically\n"
        "tuned or validated.\n"
        "\n"
        "The thermal fragmentation model may significantly modify\n"
        "hadronization observables such as particle yields, spectra,\n"
        "and strangeness production.\n"
        "\n"
        "Please validate carefully before drawing quantitative\n"
        "physics conclusions.");
  }
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
      const_cast<Pythia8::Info&>(pythia_hadron_->info));
  pythia_sigmatot_.init();

  pythia_stringflav_.initInfoPtr(
      const_cast<Pythia8::Info&>(pythia_hadron_->info));
  pythia_stringflav_.init();

  event_intermediate_.init("intermediate partons",
                           &pythia_hadron_->particleData);

  for (int imu = 0; imu < 3; imu++) {
    evecBasisAB_[imu] = ThreeVector(0., 0., 0.);
  }

  final_state_.clear();
}
void StringProcess::common_setup_pythia(Pythia8::Pythia* pythia_in,
                                        double strange_supp,
                                        double diquark_supp,
                                        double popcorn_rate, double stringz_a,
                                        double stringz_b,
                                        double string_sigma_T) {
  // choose parametrization for mass-dependent width
  pythia_in->readString("ParticleData:modeBreitWigner = 4");
  /* choose minimum transverse momentum scale
   * involved in partonic interactions */

  // transverse momentum spread in string fragmentation
  // Global Lund fragmentation
  pythia_in->readString("StringZ:aLund = " + std::to_string(stringz_a));
  pythia_in->readString("StringZ:bLund = " + std::to_string(stringz_b));

  pythia_in->readString("BeamRemnants:dampPopcorn = " +
                        std::to_string(damp_popcorn_));
  pythia_in->readString("BeamRemnants:hardRemnantBaryon = on");
  pythia_in->readString("BeamRemnants:aRemnantBaryon = " +
                        std::to_string(stringz_a_leading_));
  pythia_in->readString("BeamRemnants:bRemnantBaryon = " +
                        std::to_string(stringz_b_leading_));
  pythia_in->readString("StringPT:sigma = " + std::to_string(string_sigma_T));
  // diquark suppression factor in string fragmentation
  pythia_in->readString("StringFlav:probQQtoQ = " +
                        std::to_string(diquark_supp));
  // strangeness suppression factor in string fragmentation
  pythia_in->readString("StringFlav:probStoUD = " +
                        std::to_string(strange_supp));
  pythia_in->readString("StringFlav:popcornRate = " +
                        std::to_string(popcorn_rate));
  pythia_in->readString("StringFlav:probSQtoQQ = " +
                        std::to_string(prob_sq_to_qq_));

  pythia_in->readString("StringFlav:popcornSpair = " +
                        std::to_string(popcorn_spair_));

  pythia_in->readString("StringFlav:popcornSmeson = " +
                        std::to_string(popcorn_smeson_));

  pythia_in->readString("StringFlav:probQQ1toQQ0 = " +
                        std::to_string(prob_qq1_to_qq0_));

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
  for (auto& ptype : ParticleType::list_all()) {
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

void StringProcess::form_intermediate_particles(
    ParticleList& intermediate_particles, const FourVector& uString,
    const ThreeVector& evecLong, double additional_xsec_supp,
    bool find_and_scale_leading) {
  const int nfrag = intermediate_particles.size();
  assert(nfrag > 0);

  int bstring = 0;
  for (const ParticleData& data : intermediate_particles) {
    bstring += data.pdgcode().baryon_number();
  }

  if (find_and_scale_leading) {
    assign_all_scaling_factors(bstring, intermediate_particles, evecLong,
                               additional_xsec_supp);
  }

  const ThreeVector vstring = uString.velocity();

  for (ParticleData& particle : intermediate_particles) {
    const FourVector p_string = particle.momentum();
    const ThreeVector velocity_string = p_string.velocity();
    const double gamma_string = 1.0 / particle.inverse_gamma();

    const FourVector p_com = p_string.lorentz_boost(-vstring);
    particle.set_4momentum(p_com);

    if (mass_dependent_formation_times_) {
      const double tau_prod =
          M_SQRT2 * particle.effective_mass() / kappa_tension_string_;

      const double t_prod_string = tau_prod * gamma_string;

      FourVector fragment_position(t_prod_string,
                                   t_prod_string * velocity_string);

      fragment_position = fragment_position.lorentz_boost(-vstring);
      fragment_position = fragment_position.lorentz_boost(-vcomAB_);

      particle.set_slow_formation_times(
          time_collision_,
          soft_t_form_ * fragment_position.x0() + time_collision_);
    } else {
      const ThreeVector v_lab = p_com.lorentz_boost(-vcomAB_).velocity();
      const double gamma_lab = 1.0 / std::sqrt(1.0 - v_lab.sqr());

      particle.set_slow_formation_times(
          time_collision_, time_formation_const_ * gamma_lab + time_collision_);
    }
  }
}
void StringProcess::init(const ParticleList& incoming, double tcoll) {
  PDGcodes_[0] = incoming[0].pdgcode();
  PDGcodes_[1] = incoming[1].pdgcode();
  massA_ = incoming[0].effective_mass();
  massB_ = incoming[1].effective_mass();

  plab_[0] = incoming[0].momentum();
  plab_[1] = incoming[1].momentum();

  sqrtsAB_ = (plab_[0] + plab_[1]).abs();
  ucomAB_ = (plab_[0] + plab_[1]) / sqrtsAB_;
  vcomAB_ = ucomAB_.velocity();

  const Pythia8::Vec4 pA_lab = make_pythia_4vec(plab_[0]);
  const Pythia8::Vec4 pB_lab = make_pythia_4vec(plab_[1]);

  to_cm_.reset();
  to_cm_.toCMframe(pA_lab, pB_lab);

  const Pythia8::Vec4 pA_cm = to_cm_ * pA_lab;
  const Pythia8::Vec4 pB_cm = to_cm_ * pB_lab;

  pcom_[0] = make_smash_4vec(pA_cm);
  pcom_[1] = make_smash_4vec(pB_cm);
  ThreeVector evec_polar(pA_cm.px(), pA_cm.py(), pA_cm.pz());
  evec_polar /= std::sqrt(evec_polar.sqr());

  make_orthonormal_basis(evec_polar, evecBasisAB_);
  compute_incoming_lightcone_momenta();

  time_collision_ = tcoll;
}

bool StringProcess::next(ProcessType type) {
  string_parton_events_.clear();
  final_state_.clear();

  bool parton_level_success = false;

  switch (type) {
    case ProcessType::StringSoftNonDiffractive:
      parton_level_success = next_NDiffSoft();
      break;
    case ProcessType::StringSoftSingleDiffractiveAX:
      parton_level_success = next_SDiff(true);
      break;
    case ProcessType::StringSoftSingleDiffractiveXB:
      parton_level_success = next_SDiff(false);
      break;
    case ProcessType::StringSoftDoubleDiffractive:
      parton_level_success = next_DDiff();
      break;
    case ProcessType::StringSoftAnnihilation:
      parton_level_success = next_BBbarAnn();
      break;
    case ProcessType::StringHardNonDiffractive:
    case ProcessType::StringHardDoubleDiffractive:
    case ProcessType::StringHardSingleDiffractiveAX:
    case ProcessType::StringHardSingleDiffractiveXB:
      parton_level_success = next_Hard(type);
      break;
    default:
      logg[LPythia].error("Unknown string process required.");
      return false;
  }

  if (!parton_level_success) {
    return false;
  }

  for (const Pythia8::Event& string_event : string_parton_events_) {
    auto hadrons = hadronize(string_event);

    if (!hadrons) {
      final_state_.clear();
      return false;
    }

    final_state_.insert(final_state_.end(), hadrons->begin(), hadrons->end());
  }

  return true;
}

std::optional<ParticleList> StringProcess::hadronize(
    const Pythia8::Event& string_evt) {
  ParticleList intermediate_particles;

  bool has_leading_parton = false;
  bool has_parton = false;

  const bool has_junction = string_evt.sizeJunction() > 0;

  for (const Pythia8::Particle& particle : string_evt) {
    if (is_leading_parton(particle)) {
      has_leading_parton = true;

      if (!particle.isGluon()) {
        has_parton = true;
      }
    }

    if (has_leading_parton && has_parton) {
      break;
    }
  }

  const bool has_leading_and_is_open = has_leading_parton && !has_junction;

  // usually 0 is "system", so need at least indices 1 and 2
  if (string_evt.size() < 3) {
    logg[LPythia].error("String event too small to hadronize.");
    return std::nullopt;
  }

  pythia_hadron_->event = string_evt;

  const Pythia8::Vec4 p_str = string_evt[0].p();

  Pythia8::RotBstMatrix to_string_rest;
  to_string_rest.bstback(p_str);

  pythia_hadron_->event.rotbst(to_string_rest);

  if (!pythia_hadron_->forceHadronLevel(false)) {
    logg[LPythia].error("Pythia fragmentation failed for one string.");
    return std::nullopt;
  }

  if (has_leading_and_is_open) {
    tag_leading_hadron(pythia_hadron_->event, pythia_hadron_->particleData);
  }

  const FourVector pString_smash(
      p_str.e(), ThreeVector(p_str.px(), p_str.py(), p_str.pz()));

  const double m2 = pString_smash.sqr();

  if (m2 <= 0.0) {
    logg[LPythia].error("String has non-positive invariant mass.");
    return std::nullopt;
  }

  const double m = std::sqrt(m2);
  const FourVector uString = pString_smash / m;

  // evecLong: unit vector along string momentum
  // fallback if |p| ~ 0
  ThreeVector evecLong(p_str.px(), p_str.py(), p_str.pz());

  const double pabs = evecLong.abs();

  if (pabs > 1e-12) {
    evecLong = evecLong / pabs;
  } else {
    evecLong = ThreeVector(0.0, 0.0, 1.0);
  }

  // Build intermediate SMASH particles from final hadrons
  intermediate_particles.reserve(pythia_hadron_->event.size());

  for (int i = 0; i < pythia_hadron_->event.size(); ++i) {
    const auto& particle = pythia_hadron_->event[i];

    if (!(particle.isFinal() && particle.isHadron())) {
      continue;
    }

    const FourVector mom_smash = make_smash_4vec(particle.p());

    int particle_id = particle.id();
    convert_KaonLS(particle_id);

    if (!append_intermediate_list(particle_id, mom_smash,
                                  intermediate_particles)) {
      logg[LPythia].error("Unknown hadron in SMASH during hadronization: PDG=",
                          particle.id());

      return std::nullopt;
    }

    if (is_leading_from_diquark(particle)) {
      intermediate_particles.back().set_cross_section_scaling_factor(
          2.0 / 3.0 * additional_xsec_supp_);

    } else if (is_leading_from_quark(particle)) {
      intermediate_particles.back().set_cross_section_scaling_factor(
          (pythia_hadron_->particleData.isBaryon(particle.id()) ? 1.0 / 3.0
                                                                : 0.5) *
          additional_xsec_supp_);

    } else {
      intermediate_particles.back().set_cross_section_scaling_factor(0.0);
    }
  }

  if (intermediate_particles.empty()) {
    logg[LPythia].error("Hadronization produced no final hadrons.");
    return std::nullopt;
  }

  const bool should_assign_scaling = has_junction && has_leading_parton;

  form_intermediate_particles(intermediate_particles, uString, evecLong,
                              additional_xsec_supp_, should_assign_scaling);

  return intermediate_particles;
}
bool StringProcess::append_string(const Pythia8::Vec4& p_str,
                                  const std::array<int, 2>& ends, int color_tag,
                                  bool from_projectile) {
  const int id1 = ends[0];
  const int id2 = ends[1];

  const auto& pd = pythia_hadron_->particleData;

  const double m_string = p_str.mCalc();

  if (m_string < estimate_string_threshold(id1, id2))
    return false;

  const double m1 = pd.m0(id1);
  const double m2 = pd.m0(id2);

  const double p_cm = pCM(m_string, m1, m2);
  if (p_cm < really_small)
    return false;

  const double E1 = std::sqrt(p_cm * p_cm + m1 * m1);
  const double E2 = std::sqrt(p_cm * p_cm + m2 * m2);

  const int ibeam = from_projectile ? 0 : 1;

  const ThreeVector mom = pcom_[ibeam].threevec();
  Pythia8::Vec4 beam_axis_cm(mom.x1(), mom.x2(), mom.x3(), mom.abs());

  Pythia8::RotBstMatrix to_string_rest;
  to_string_rest.bstback(p_str);

  Pythia8::Vec4 beam_axis_rest = to_string_rest * beam_axis_cm;

  // Boost beam direction into the string rest frame.
  to_string_rest.bstback(p_str);

  const double norm = std::sqrt(beam_axis_rest.px() * beam_axis_rest.px() +
                                beam_axis_rest.py() * beam_axis_rest.py() +
                                beam_axis_rest.pz() * beam_axis_rest.pz());

  if (norm < really_small)
    return false;

  const double nx = beam_axis_rest.px() / norm;
  const double ny = beam_axis_rest.py() / norm;
  const double nz = beam_axis_rest.pz() / norm;

  const Pythia8::Vec4 p1_rest(nx * p_cm, ny * p_cm, nz * p_cm, E1);
  const Pythia8::Vec4 p2_rest(-nx * p_cm, -ny * p_cm, -nz * p_cm, E2);
  Pythia8::RotBstMatrix boost;
  boost.bst(p_str);

  const Pythia8::Vec4 p1_cm = boost * p1_rest;
  const Pythia8::Vec4 p2_cm = boost * p2_rest;

  // Create event for this string.
  string_parton_events_.emplace_back();
  Pythia8::Event& evt = string_parton_events_.back();
  evt.init("Soft SMASH string", &pythia_hadron_->particleData);

  evt.append(90, -11, 0, 0, p1_cm + p2_cm, (p1_cm + p2_cm).mCalc());

  const int status_quark = static_cast<int>(LeadingStatus::LEADING_PARTON);
  const int status_diquark = static_cast<int>(LeadingStatus::LEADING_DIQUARK);

  const bool id1_is_quark = pythia_hadron_->particleData.isQuark(id1);
  const bool id2_is_quark = pythia_hadron_->particleData.isQuark(id2);

  const int i1 = evt.append(id1, id1_is_quark ? status_quark : status_diquark,
                            0, 0, p1_cm, m1);
  const int i2 = evt.append(id2, id2_is_quark ? status_quark : status_diquark,
                            0, 0, p2_cm, m2);

  set_color_by_type(evt[i1], color_tag);
  set_color_by_type(evt[i2], color_tag);
  evt.rotbst(to_cm_.inverse());
  return true;
}

bool StringProcess::next_SDiff(bool is_AB_to_AX) {
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

  // determine three momentum of the final state hadron
  double sign_direction = is_AB_to_AX ? 1. : -1.;

  // longitudinal component (along +z by convention)
  const double pL = std::sqrt(pabscomHX_sqr - QTrn * QTrn);

  // momentum of the outgoing hadron H in the CM frame
  const ThreeVector cm_momentum(sign_direction * QTrx, sign_direction * QTry,
                                sign_direction * pL);

  const FourVector pstrHcom(std::sqrt(pabscomHX_sqr + massH * massH),
                            cm_momentum);
  const FourVector pstrXcom(std::sqrt(pabscomHX_sqr + massX * massX),
                            -cm_momentum);

  if (!append_string(make_pythia_4vec(pstrXcom), {idqX2, idqX1}, 101,
                     !is_AB_to_AX))
    return false;
  // string_parton_events_.back().list();
  PdgCode hadron_code = is_AB_to_AX ? PDGcodes_[0] : PDGcodes_[1];
  ParticleData new_particle(ParticleType::find(hadron_code));

  auto had_mom = make_pythia_4vec(pstrHcom);
  had_mom.rotbst(to_cm_.inverse());

  new_particle.set_4momentum(make_smash_4vec(had_mom));
  new_particle.set_cross_section_scaling_factor(1.);
  new_particle.set_formation_time(time_collision_);
  final_state_.push_back(new_particle);

  return true;
}

// double-diffractive : A + B -> X + X
bool StringProcess::next_DDiff() {
  std::array<std::array<int, 2>, 2> quarks;

  // decompose hadron into quark (and diquark) contents
  make_string_ends(PDGcodes_[0], quarks[0][1], quarks[0][0],
                   prob_proton_to_d_uu_);
  make_string_ends(PDGcodes_[1], quarks[1][1], quarks[1][0],
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

  const double pz1 = ((PPosA_ + QPos) - (PNegA_ + QNeg)) * M_SQRT1_2;
  const double E1 = ((PPosA_ + QPos) + (PNegA_ + QNeg)) * M_SQRT1_2;
  Pythia8::Vec4 p_str1(QTrx, QTry, pz1, E1);

  const double pz2 = ((PPosB_ - QPos) - (PNegB_ - QNeg)) * M_SQRT1_2;
  const double E2 = ((PPosB_ - QPos) + (PNegB_ - QNeg)) * M_SQRT1_2;
  Pythia8::Vec4 p_str2(-QTrx, -QTry, pz2, E2);
  if (random::uniform_int(0, 1) == 0) {
    std::swap(quarks[0][0], quarks[0][1]);
  }
  if (random::uniform_int(0, 1) == 0) {
    std::swap(quarks[1][0], quarks[1][1]);
  }
  const bool added = append_string(p_str1, quarks[0], 101, true) &&
                     append_string(p_str2, quarks[1], 102, false);

  if (!added)
    return false;

  return true;
}

bool StringProcess::next_BBbarAnn() {
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
      ParticleData new_particle(ParticleType::find(PDGcodes_[i]));
      Pythia8::Vec4 pcom_pyth = make_pythia_4vec(pcom_[i]);
      pcom_pyth.rotbst(to_cm_.inverse());
      new_particle.set_4momentum(make_smash_4vec(pcom_pyth));
      new_particle.set_cross_section_scaling_factor(1.);
      new_particle.set_formation_time(time_collision_);
      final_state_.push_back(new_particle);
    }
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
    append_string(make_pythia_4vec(pcom_[i]),
                  {remaining_quarks[i], remaining_antiquarks[i]}, 101,
                  random::canonical() > 0.5);
  }

  return true;
}

// soft non-diffractive
bool StringProcess::next_NDiffSoft() {
  // ids of the two ends of each string
  std::array<int, 2> endsA;  // string built with p_strA
  std::array<int, 2> endsB;  // string built with p_strB

  // decompose hadron into quark (and diquark) contents
  int idqA1, idqA2, idqB1, idqB2;
  make_string_ends(PDGcodes_[0], idqA1, idqA2, prob_proton_to_d_uu_);
  make_string_ends(PDGcodes_[1], idqB1, idqB2, prob_proton_to_d_uu_);

  const int bar_a = PDGcodes_[0].baryon_number();
  const int bar_b = PDGcodes_[1].baryon_number();

  if (bar_a == 1 ||                  // baryon-*
      (bar_a == 0 && bar_b == 1) ||  // meson-baryon
      (bar_a == 0 && bar_b == 0)) {  // meson-meson
    endsA = {idqA2, idqB1};
    endsB = {idqB2, idqA1};
  } else if ((bar_a == 0 && bar_b == -1) ||  // meson-antibaryon
             (bar_a == -1)) {                // antibaryon-*
    endsA = {idqB2, idqA1};
    endsB = {idqA2, idqB1};
  } else {
    std::stringstream ss;
    ss << "StringProcess::next_NDiffSoft: baryonA = " << bar_a
       << ", baryonB = " << bar_b;
    throw std::runtime_error(ss.str());
  }

  // CM frame with A along +z and B along -z; PPos/PNeg are lightcone components
  // in this frame.

  // sample the lightcone momentum fraction carried by quarks
  const double xfracA = random::beta(pow_fquark_alpha_, pow_fquark_beta_);
  const double xfracB = random::beta(pow_fquark_alpha_, pow_fquark_beta_);

  // sample the transverse momentum transfer
  const double qx = random::normal(0., sigma_qperp_ * M_SQRT1_2);
  const double qy = random::normal(0., sigma_qperp_ * M_SQRT1_2);
  const double qT2 = qx * qx + qy * qy;

  // evaluate the lightcone momentum transfer
  const double QPos = -qT2 / (2. * xfracB * PNegB_);
  const double QNeg = qT2 / (2. * xfracA * PPosA_);
  const double dPPos = -xfracA * PPosA_ - QPos;
  const double dPNeg = xfracB * PNegB_ - QNeg;

  // string from hadron A side (with +qT)
  const double pzA = ((PPosA_ + dPPos) - (PNegA_ + dPNeg)) * M_SQRT1_2;
  const double EA = ((PPosA_ + dPPos) + (PNegA_ + dPNeg)) * M_SQRT1_2;
  Pythia8::Vec4 p_strA(qx, qy, pzA, EA);

  // string from hadron B side (with -qT)
  const double pzB = ((PPosB_ - dPPos) - (PNegB_ - dPNeg)) * M_SQRT1_2;
  const double EB = ((PPosB_ - dPPos) + (PNegB_ - dPNeg)) * M_SQRT1_2;
  Pythia8::Vec4 p_strB(-qx, -qy, pzB, EB);
  if (!append_string(p_strA, endsA, 101, true) ||
      !append_string(p_strB, endsB, 102, false))
    return false;

  return true;
}

std::vector<bool> StringProcess::compute_beam_valence_flags(
    Pythia8::Pythia& pythia) {
  const Pythia8::Event& event = pythia.event;

  auto find_final_copy = [&](int iPos) -> int {
    if (iPos <= 0 || iPos >= event.size())
      return -1;
    if (event[iPos].isFinal())
      return iPos;

    const int id = event[iPos].id();
    const auto ds = event[iPos].daughterListRecursive();
    for (int j : ds) {
      if (j > 0 && j < event.size() && event[j].isFinal() &&
          event[j].id() == id)
        return j;
    }
    return -1;
  };

  auto tag_from_beam = [&](const Pythia8::BeamParticle& beam,
                           std::vector<bool>& isValenceFinal) -> void {
    for (int i = 0; i < beam.size(); ++i) {
      const int j = find_final_copy(beam[i].iPos());
      if (j >= 0)
        isValenceFinal[j] = beam[i].isValence();
    }
  };

  std::vector<bool> isValenceFinal(event.size(), false);
  tag_from_beam(pythia.beamA, isValenceFinal);
  tag_from_beam(pythia.beamB, isValenceFinal);
  return isValenceFinal;
}

void StringProcess::tag_leading_hadron(Pythia8::Event& event,
                                       const Pythia8::ParticleData& pd) {
  // 1) Find the first two (di)quark endpoints in the event record
  int idx1 = -1, idx2 = -1;
  bool is_leading1 = false, is_leading2 = false;

  for (const auto& part : event) {
    if (part.isQuark() || part.isDiquark()) {
      if (idx1 == -1) {
        idx1 = part.index();
        is_leading1 = (part.statusAbs() ==
                       static_cast<int>(LeadingStatus::LEADING_PARTON));
      } else {
        idx2 = part.index();
        is_leading2 = (part.statusAbs() ==
                       static_cast<int>(LeadingStatus::LEADING_PARTON));
        break;
      }
    }
  }

  if (idx1 < 0 || idx2 < 0) {
    return;
  }

  // 2) Boost to the CM frame of the two endpoints
  Pythia8::RotBstMatrix toRest;
  toRest.toCMframe(event[idx1].p(), event[idx2].p());
  event.rotbst(toRest);

  // Helper: pick final hadron with pz closest to endpoint pz.
  // If endpoint is diquark -> require baryon.
  auto endpoint_required_flavours =
      [&](const Pythia8::Particle& end) -> std::pair<std::array<int, 2>, int> {
    std::array<int, 2> flavours{{0, 0}};

    if (end.isDiquark()) {
      int q1 = 0, q2 = 0, degeneracy = 0;
      StringProcess::quarks_from_diquark(end.id(), q1, q2, degeneracy);

      flavours[0] = std::abs(q1);
      flavours[1] = std::abs(q2);
      return {flavours, 2};
    }

    // single (anti)quark endpoint
    flavours[0] = std::abs(end.id());
    return {flavours, 1};
  };

  auto find_best_hadron = [&](int iEnd) -> int {
    const auto& end = event[iEnd];
    const double pzEnd = end.pz();
    const bool endIsDiquark = end.isDiquark();

    auto [req, nreq] = endpoint_required_flavours(end);

    int bestIdx = -1;
    double bestScore = std::numeric_limits<double>::infinity();

    for (int i = 0; i < event.size(); ++i) {
      const auto& p = event[i];

      if (!p.isFinal() || !p.isHadron())
        continue;

      if (endIsDiquark && !pd.isBaryon(p.id()))
        continue;

      // Flavour containment check (ignore sign)
      const PdgCode had = PdgCode::from_decimal(p.id());
      bool ok = true;
      for (int k = 0; k < nreq; ++k) {
        if (!had.contains_quark(req[k])) {
          ok = false;
          break;
        }
      }
      if (!ok)
        continue;

      const double score = std::abs(p.pz() - pzEnd);
      if (score < bestScore) {
        bestScore = score;
        bestIdx = i;
      }
    }
    return bestIdx;
  };

  // 3) Tag matching hadrons (still in CM frame, but indices/status are the
  // same)

  if (is_leading1) {
    const int h1 = find_best_hadron(idx1);
    if (h1 >= 0) {
      event[h1].status(leading_hadron_status_from_endpoint(event[idx1]));
    }
  }

  if (is_leading2) {
    const int h2 = find_best_hadron(idx2);
    if (h2 >= 0) {
      event[h2].status(leading_hadron_status_from_endpoint(event[idx2]));
    }
  }

  // 4) Boost back
  event.rotbst(toRest.inverse());
}

// hard non-diffractive
bool StringProcess::next_Hard(ProcessType type) {
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
    hard_map_[idAB]->readString("SoftQCD:singleDiffractiveXB = on");
    hard_map_[idAB]->readString("SoftQCD:singleDiffractiveAX = on");
    hard_map_[idAB]->readString("SoftQCD:doubleDiffractive = on");
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
  hard_map_[idAB]->setKinematics(sqrtsAB_);

  bool final_state_success = false;
  /* Hard-process codes follow the Pythia convention:
   * https://www.pythia.org/latest-manual/QCDSoftProcesses.html
   *
   * Note: Pythia effectively classifies processes by the last digit
   * of the code (e.g. 101 -> 1, 102 -> 2, ...).
   *
   * The terminology "soft" vs "hard" is somewhat misleading:
   * Pythia's SoftQCD corresponds to a semi-perturbative (semi-hard)
   * model valid down to low pT, while HardQCD represents purely
   * perturbative processes and is not applicable over the full pT range.
   */

  switch (type) {
    case ProcessType::StringHardNonDiffractive:
      final_state_success = hard_map_[idAB]->next(1);
      break;
    case ProcessType::StringHardDoubleDiffractive:
      final_state_success = hard_map_[idAB]->next(5);
      break;
    case ProcessType::StringHardSingleDiffractiveAX:
      final_state_success = hard_map_[idAB]->next(4);
      break;
    case ProcessType::StringHardSingleDiffractiveXB:
      final_state_success = hard_map_[idAB]->next(3);
      break;
    default:
      logg[LPythia].error("Unknown string process required.");
      final_state_success = false;
      break;
  }

  logg[LPythia].debug("Pythia final state computed, success = ",
                      final_state_success);
  if (!final_state_success) {
    return false;
  }
  auto valance_tags = compute_beam_valence_flags(*hard_map_[idAB]);
  for (auto& p : hard_map_[idAB]->event) {
    if (valance_tags[p.index()] && p.isFinal()) {
      p.statusCode(static_cast<int>(p.isDiquark()
                                        ? LeadingStatus::LEADING_DIQUARK
                                        : LeadingStatus::LEADING_PARTON));
    }
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
      const double mass = pquark.mCalc();
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
  // Boost after the constituents have been corrected!
  event_intermediate_.rotbst(to_cm_.inverse());
  int npart = event_intermediate_.size();
  int ipart = 0;
  while (ipart < npart) {
    const int pdgid = event_intermediate_[ipart].id();
    if (event_intermediate_[ipart].isFinal() &&
        !event_intermediate_[ipart].isParton() &&
        !hard_map_[idAB]->particleData.isOctetHadron(pdgid)) {
      logg[LPythia].debug("PDG ID from Pythia: ", pdgid);
      Pythia8::Vec4 momentum_pythia = event_intermediate_[ipart].p();
      FourVector momentum = make_smash_4vec(momentum_pythia);

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

  for (ParticleData& non_hadron : new_non_hadron_particles) {
    non_hadron.set_cross_section_scaling_factor(1.);
    non_hadron.set_formation_time(time_collision_);
    final_state_.push_back(non_hadron);
  }
  logg[LPythia].debug("Hard non-diff: partonic process gives ",
                      event_intermediate_.size(), " partons.");
  bool find_forward_string = true;
  while (event_intermediate_.size() > 1) {
    string_parton_events_.emplace_back();
    string_parton_events_.back().init("Hard string",
                                      &pythia_hadron_->particleData);
    if (event_intermediate_.sizeJunction() > 0) {
      // identify string from a junction if there is any.
      compose_string_junction(find_forward_string, event_intermediate_,
                              string_parton_events_.back());
    } else {
      /* identify string from a most forward or backward parton.
       * if there is no junction. */
      compose_string_parton(find_forward_string, event_intermediate_,
                            string_parton_events_.back());
    }
    if (!string_above_threshold(string_parton_events_.back())) {
      return false;
    }

    find_forward_string = !find_forward_string;
  }

  return true;
}

double StringProcess::estimate_string_threshold(int p_left, int p_right) {
  Pythia8::ParticleData& pd = pythia_hadron_->particleData;
  bool left_is_dq = pd.isDiquark(p_left);
  bool right_is_dq = pd.isDiquark(p_right);

  bool left_is_q_or_dq = pd.isQuark(p_left) || pd.isDiquark(p_left);
  bool right_is_q_or_dq = pd.isQuark(p_right) || pd.isDiquark(p_right);

  if (!(left_is_q_or_dq && right_is_q_or_dq)) {
    std::cout << p_left << ", " << p_right << std::endl;
    string_parton_events_.back().list();
    throw std::runtime_error(
        "Threshold can got input that is not quark or diquark.");
  }
  // Case: both ends are diquarks → estimate as two separate baryons
  if (left_is_dq && right_is_dq) {
    int h1 = pythia_stringflav_.combineToLightest(p_left, 2);    // u
    int h2 = pythia_stringflav_.combineToLightest(p_right, -2);  // ū

    double m1 = pd.m0(h1);
    double m2 = pd.m0(h2);

    return std::min(m1, m2);
  }

  int h = pythia_stringflav_.combineToLightest(p_left, p_right);
  PdgCode pdg = PdgCode::from_decimal(h);
  double mass = ParticleType::find(pdg).mass() + 0.150;
  return mass;
}

void StringProcess::set_color_by_type(Pythia8::Particle& p, int color) {
  const int id = p.id();
  const bool is_anti = (id < 0);

  const auto& pd = pythia_hadron_->particleData;
  const bool is_quark = pd.isQuark(id);
  const bool is_diquark = pd.isDiquark(id);

  if ((is_quark && !is_anti) || (is_diquark && is_anti)) {
    // quark or antidiquark: color flows out
    p.cols(color, 0);
  } else if ((is_quark && is_anti) || (is_diquark && !is_anti)) {
    // antiquark or diquark: color flows in
    p.cols(0, color);
  } else {
    throw std::runtime_error("Cannot set color: not a quark or diquark (id=" +
                             std::to_string(id) + ")");
  }
}
bool StringProcess::string_above_threshold(const Pythia8::Event& event) {
  const double string_mass = event[0].mCalc();
  const auto& particle_data = pythia_hadron_->particleData;

  std::vector<int> endpoints;
  endpoints.reserve(event.size());

  for (const Pythia8::Particle& p : event) {
    if (!p.isFinal() || p.isGluon() || !p.isParton()) {
      continue;
    }

    const int id = p.id();

    if (!particle_data.isQuark(id) && !particle_data.isDiquark(id)) {
      logg[LPythia].error("Colored non-quark/diquark in threshold check: id=",
                          id, ", col=", p.col(), ", acol=", p.acol());
      return false;
    }

    endpoints.push_back(id);
  }

  // Closed gluon string.
  if (endpoints.size() < 2) {
    return string_mass > 0.28;
  }

  // Ordinary open string.
  if (endpoints.size() == 2) {
    return string_mass > estimate_string_threshold(endpoints[0], endpoints[1]);
  }

  // Baryonic/antibaryonic endpoint: combine two same-sign quarks.
  if (endpoints.size() == 3) {
    for (std::size_t i = 0; i < endpoints.size(); ++i) {
      for (std::size_t j = i + 1; j < endpoints.size(); ++j) {
        const int id_i = endpoints[i];
        const int id_j = endpoints[j];

        if (!particle_data.isQuark(id_i) || !particle_data.isQuark(id_j)) {
          continue;
        }

        if (id_i * id_j < 0) {
          continue;
        }

        const int diquark_id = diquark_from_quarks(id_i, id_j);

        const std::size_t k = 3 - i - j;
        const int remaining_id = endpoints[k];

        return string_mass >
               estimate_string_threshold(diquark_id, remaining_id);
      }
    }

    logg[LPythia].error(
        "Cannot reduce 3 string endpoints to quark-diquark "
        "threshold: ",
        endpoints[0], ", ", endpoints[1], ", ", endpoints[2]);
    return false;
  }

  logg[LPythia].error(
      "Unexpected number of string endpoints in threshold check: ",
      endpoints.size());
  return false;
}
void StringProcess::find_excess_constituent(PdgCode& pdg_actual,
                                            PdgCode& pdg_mapped,
                                            std::array<int, 5>& excess_quark,
                                            std::array<int, 5>& excess_antiq) {
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
    Pythia8::Particle& particle, std::array<int, 5>& excess_constituent) {
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
    Pythia8::Event& event_intermediate, std::array<int, 5>& nquark_total,
    std::array<int, 5>& nantiq_total) {
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
    Pythia8::Event& event_intermediate, std::array<int, 5>& nquark_total,
    std::array<int, 5>& nantiq_total, bool sign_constituent,
    std::array<std::array<int, 5>, 2>& excess_constituent) {
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
    std::array<int, 5>& nquark_total,
    std::array<std::array<int, 5>, 2>& excess_quark,
    std::array<std::array<int, 5>, 2>& excess_antiq) {
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
    Pythia8::Event& event_intermediate,
    std::array<std::array<int, 5>, 2>& excess_quark,
    std::array<std::array<int, 5>, 2>& excess_antiq) {
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
        const double mass = pquark.mCalc();

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
                                          Pythia8::Event& event_intermediate,
                                          Pythia8::Event& event_hadronize) {
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

void StringProcess::compose_string_junction(bool& find_forward_string,
                                            Pythia8::Event& event_intermediate,
                                            Pythia8::Event& event_hadronize) {
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

void StringProcess::find_junction_leg(bool sign_color, std::vector<int>& col,
                                      Pythia8::Event& event_intermediate,
                                      Pythia8::Event& event_hadronize) {
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

void StringProcess::make_orthonormal_basis(
    ThreeVector& evec_polar, std::array<ThreeVector, 3>& evec_basis) {
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

void StringProcess::quarks_from_diquark(int diquark, int& q1, int& q2,
                                        int& deg_spin) {
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

void StringProcess::make_string_ends(const PdgCode& pdg, int& idq1, int& idq2,
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

void StringProcess::assign_scaling_factor(int nquark, ParticleData& data,
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
                                                ParticleList& list) {
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
                                               ParticleList& outgoing_particles,
                                               const ThreeVector& evecLong,
                                               double suppression_factor) {
  // Set each particle's cross section scaling factor to 0 first
  for (ParticleData& data : outgoing_particles) {
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

int StringProcess::pdg_map_for_pythia(PdgCode& pdg) {
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
