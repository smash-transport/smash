/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/particletype.h"

#include <assert.h>
#include <algorithm>
#include <map>
#include <vector>

#include "smash/constants.h"
#include "smash/cxx14compat.h"
#include "smash/decaymodes.h"
#include "smash/distributions.h"
#include "smash/formfactors.h"
#include "smash/inputfunctions.h"
#include "smash/integrate.h"
#include "smash/iomanipulators.h"
#include "smash/isoparticletype.h"
#include "smash/kinematics.h"
#include "smash/logging.h"
#include "smash/numerics.h"
#include "smash/particledata.h"
#include "smash/pdgcode.h"
#include "smash/pow.h"
#include "smash/processbranch.h"
#include "smash/stringfunctions.h"

namespace smash {

#ifdef SMASH_INLINE_LIST_ALL
const ParticleTypeList *all_particle_types = nullptr;
#else
namespace {
/// Global pointer to the Particle Type list.
const ParticleTypeList *all_particle_types = nullptr;
/// Global pointer to the Particle Type list of nucleons
ParticleTypePtrList nucleons_list;
/// Global pointer to the Particle Type list of anti-nucleons
ParticleTypePtrList anti_nucs_list;
/// Global pointer to the Particle Type list of deltas
ParticleTypePtrList deltas_list;
/// Global pointer to the Particle Type list of anti-deltas
ParticleTypePtrList anti_deltas_list;
/// Global pointer to the Particle Type list of baryon resonances
ParticleTypePtrList baryon_resonances_list;
/// Global pointer to the Particle Type list of light nuclei
ParticleTypePtrList light_nuclei_list;
}  // unnamed namespace

const ParticleTypeList &ParticleType::list_all() {
  assert(all_particle_types);
  return *all_particle_types;
}

ParticleTypePtr ParticleType::operator&() const {
  // Calculate the offset via pointer subtraction:
  const auto offset = this - std::addressof(list_all()[0]);
  // Since we're using uint16_t for storing the index better be safe than sorry:
  // The offset must fit into the data type. If this ever fails we got a lot
  // more particle types than initially expected and you have to increase the
  // ParticleTypePtr storage to uint32_t.
  assert(offset >= 0 && offset < 0xffff);
  // After the assertion above the down-cast to uint16_t is safe:
  return ParticleTypePtr(static_cast<uint16_t>(offset));
}
#endif

ParticleTypePtrList &ParticleType::list_nucleons() { return nucleons_list; }

ParticleTypePtrList &ParticleType::list_anti_nucleons() {
  return anti_nucs_list;
}

ParticleTypePtrList &ParticleType::list_Deltas() { return deltas_list; }

ParticleTypePtrList &ParticleType::list_anti_Deltas() {
  return anti_deltas_list;
}

ParticleTypePtrList &ParticleType::list_baryon_resonances() {
  return baryon_resonances_list;
}

ParticleTypePtrList &ParticleType::list_light_nuclei() {
  return light_nuclei_list;
}

const ParticleTypePtr ParticleType::try_find(PdgCode pdgcode) {
  const auto found = std::lower_bound(
      all_particle_types->begin(), all_particle_types->end(), pdgcode,
      [](const ParticleType &l, const PdgCode &r) { return l.pdgcode() < r; });
  if (found == all_particle_types->end() || found->pdgcode() != pdgcode) {
    return {};  // The default constructor creates an invalid pointer.
  }
  return &*found;
}

const ParticleType &ParticleType::find(PdgCode pdgcode) {
  const auto found = ParticleType::try_find(pdgcode);
  if (!found) {
    throw PdgNotFoundFailure("PDG code " + pdgcode.string() + " not found!");
  }
  return *found;
}

bool ParticleType::exists(PdgCode pdgcode) {
  const auto found = ParticleType::try_find(pdgcode);
  return found;
}

bool ParticleType::exists(const std::string &name) {
  const auto found =
      std::find_if(all_particle_types->begin(), all_particle_types->end(),
                   [&](const ParticleType &p) { return p.name() == name; });
  if (found == all_particle_types->end()) {
    return false;
  }
  return true;
}

ParticleType::ParticleType(std::string n, double m, double w, Parity p,
                           PdgCode id)
    : name_(n),
      mass_(m),
      width_(w),
      parity_(p),
      pdgcode_(id),
      min_mass_kinematic_(-1.),
      min_mass_spectral_(-1.),
      charge_(pdgcode_.charge()),
      isospin_(-1),
      I3_(pdgcode_.isospin3()) {}

/**
 * Construct an antiparticle name-string from the given name-string for the
 * particle and its PDG code.
 *
 * \param[in] name the name-string of the particle to convert
 * \param[in] code the pdgcode of the particle to convert
 * \return the name-string of the converted antiparticle
 */
static std::string antiname(const std::string &name, PdgCode code) {
  std::string basename, charge;

  if (name.find("⁺⁺") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁺⁺") + 1);
    charge = "⁻⁻";
  } else if (name.find("⁺") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁺") + 1);
    charge = "⁻";
  } else if (name.find("⁻⁻") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁻⁻") + 1);
    charge = "⁺⁺";
  } else if (name.find("⁻") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁻") + 1);
    charge = "⁺";
  } else if (name.find("⁰") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁰") + 1);
    charge = "⁰";
  } else {
    basename = name;
    charge = "";
  }

  // baryons & strange mesons: insert a bar
  if (code.baryon_number() != 0 || code.strangeness() != 0) {
    constexpr char bar[] = "\u0305";
    basename.insert(utf8::sequence_length(basename.begin()), bar);
  }

  return basename + charge;
}

/**
 * Construct a charge string, given the charge as integer.
 *
 * \param[in] charge charge of a particle
 * \return the corresponding string to write out this charge
 * \throw runtime_error if the charge is not an integer between -2 and 2
 */
static std::string chargestr(int charge) {
  switch (charge) {
    case 2:
      return "⁺⁺";
    case 1:
      return "⁺";
    case 0:
      return "⁰";
    case -1:
      return "⁻";
    case -2:
      return "⁻⁻";
    default:
      throw std::runtime_error("Invalid charge " + std::to_string(charge));
  }
}

void ParticleType::create_type_list(const std::string &input) {  // {{{
  const auto &log = logger<LogArea::ParticleType>();
  static ParticleTypeList type_list;
  type_list.clear();  // in case LoadFailure was thrown and caught and we should
                      // try again
  for (const Line &line : line_parser(input)) {
    std::istringstream lineinput(line.text);
    std::string name;
    double mass, width;
    std::string parity_string;
    std::vector<std::string> pdgcode_strings;
    // We expect at most 4 PDG codes per multiplet.
    pdgcode_strings.reserve(4);
    lineinput >> name >> mass >> width >> parity_string;
    Parity parity;
    bool fail = false;
    if (parity_string == "+") {
      parity = Parity::Pos;
    } else if (parity_string == "-") {
      parity = Parity::Neg;
    } else {
      fail = true;
    }
    if (lineinput.fail() || fail) {
      throw ParticleType::LoadFailure(build_error_string(
          "While loading the ParticleType data:\nFailed to convert the input "
          "string to the expected data types.",
          line));
    }
    // read additional PDG codes (if present)
    while (!lineinput.eof()) {
      pdgcode_strings.push_back("");
      lineinput >> pdgcode_strings.back();
      if (lineinput.fail()) {
        throw ParticleType::LoadFailure(build_error_string(
            "While loading the ParticleType data:\nFailed to convert the input "
            "string to the expected data types.",
            line));
      }
    }
    if (pdgcode_strings.size() < 1) {
      throw ParticleType::LoadFailure(build_error_string(
          "While loading the ParticleType data:\nFailed to convert the input "
          "string due to missing PDG code.",
          line));
    }
    std::vector<PdgCode> pdgcode;
    pdgcode.resize(pdgcode_strings.size());
    std::transform(pdgcode_strings.begin(), pdgcode_strings.end(),
                   pdgcode.begin(),
                   [](const std::string &s) { return PdgCode(s); });
    ensure_all_read(lineinput, line);

    /* Check if nucleon, kaon, and delta masses are
     * the same as hardcoded ones, if present */
    if (pdgcode[0].is_nucleon() && !almost_equal(mass, nucleon_mass)) {
      throw std::runtime_error(
          "Nucleon mass in input file different from 0.938");
    }
    if (pdgcode[0].is_kaon() && !almost_equal(mass, kaon_mass)) {
      throw std::runtime_error("Kaon mass in input file different from 0.494");
    }
    if (pdgcode[0].is_Delta() && !almost_equal(mass, delta_mass)) {
      throw std::runtime_error("Delta mass in input file different from 1.232");
    }
    if (pdgcode[0].is_deuteron() && !almost_equal(mass, deuteron_mass)) {
      throw std::runtime_error("d mass in input file different from 1.8756");
    }

    // add all states to type list
    for (size_t i = 0; i < pdgcode.size(); i++) {
      std::string full_name = name;
      if (pdgcode.size() > 1) {
        // for multiplets: add charge string to name
        full_name += chargestr(pdgcode[i].charge());
      }
      type_list.emplace_back(full_name, mass, width, parity, pdgcode[i]);
      log.debug() << "Setting     particle type: " << type_list.back();
      if (pdgcode[i].has_antiparticle()) {
        /* add corresponding antiparticle */
        PdgCode anti = pdgcode[i].get_antiparticle();
        // For bosons the parity does not change, for fermions it gets inverted.
        const auto anti_parity =
            (anti.spin() % 2 == 0) ? parity : invert_parity(parity);
        full_name = antiname(full_name, pdgcode[i]);
        type_list.emplace_back(full_name, mass, width, anti_parity, anti);
        log.debug() << "Setting antiparticle type: " << type_list.back();
      }
    }
  }
  type_list.shrink_to_fit();

  /* Sort the type list by PDG code. */
  std::sort(type_list.begin(), type_list.end());

  /* Look for duplicates. */
  PdgCode prev_pdg = 0;
  for (const auto &t : type_list) {
    if (t.pdgcode() == prev_pdg) {
      throw ParticleType::LoadFailure("Duplicate PdgCode in particles.txt: " +
                                      t.pdgcode().string());
    }
    prev_pdg = t.pdgcode();
  }

  if (all_particle_types != nullptr) {
    throw std::runtime_error("Error: Type list was already built!");
  }
  all_particle_types = &type_list;  // note that type_list is a function-local
                                    // static and thus will live on until after
                                    // main().

  // create all isospin multiplets
  for (const auto &t : type_list) {
    IsoParticleType::create_multiplet(t);
  }
  // link the multiplets to the types
  for (auto &t : type_list) {
    t.iso_multiplet_ = IsoParticleType::find(t);
  }

  // Create nucleons/anti-nucleons list
  if (IsoParticleType::exists("N")) {
    for (const auto state : IsoParticleType::find("N").get_states()) {
      nucleons_list.push_back(state);
      anti_nucs_list.push_back(state->get_antiparticle());
    }
  }

  // Create deltas list
  if (IsoParticleType::exists("Δ")) {
    for (const auto state : IsoParticleType::find("Δ").get_states()) {
      deltas_list.push_back(state);
      anti_deltas_list.push_back(state->get_antiparticle());
    }
  }

  // Create baryon resonances list
  for (const ParticleType &type_resonance : ParticleType::list_all()) {
    /* Only loop over baryon resonances. */
    if (type_resonance.is_stable() ||
        type_resonance.pdgcode().baryon_number() != 1) {
      continue;
    }
    baryon_resonances_list.push_back(&type_resonance);
    baryon_resonances_list.push_back(type_resonance.get_antiparticle());
  }

  for (const ParticleType &type : ParticleType::list_all()) {
    if (type.is_nucleus()) {
      light_nuclei_list.push_back(&type);
    }
  }
} /*}}}*/

double ParticleType::min_mass_kinematic() const {
  if (unlikely(min_mass_kinematic_ < 0.)) {
    /* If the particle is stable, min. mass is just the mass. */
    min_mass_kinematic_ = mass_;
    /* Otherwise, find the lowest mass value needed in any decay mode */
    if (!is_stable()) {
      for (const auto &mode : decay_modes().decay_mode_list()) {
        min_mass_kinematic_ = std::min(min_mass_kinematic_, mode->threshold());
      }
    }
  }
  return min_mass_kinematic_;
}

double ParticleType::min_mass_spectral() const {
  if (unlikely(min_mass_spectral_ < 0.)) {
    /* If the particle is stable or it has a non-zero spectral function value at
     * the minimum mass that is allowed by kinematics, min_mass_spectral is just
     * the min_mass_kinetic. */
    min_mass_spectral_ = min_mass_kinematic();
    /* Otherwise, find the lowest mass value where spectral function has a
     * non-zero value by bisection.*/
    if (!is_stable() &&
        this->spectral_function(min_mass_kinematic()) < really_small) {
      // find a right bound that has non-zero spectral function for bisection
      const double m_step = 0.01;
      double right_bound_bis;
      for (unsigned int i = 0;; i++) {
        right_bound_bis = min_mass_kinematic() + m_step * i;
        if (this->spectral_function(right_bound_bis) > really_small) {
          break;
        }
      }
      // bisection
      const double precision = 1E-6;
      double left_bound_bis = right_bound_bis - m_step;
      while (right_bound_bis - left_bound_bis > precision) {
        const double mid = (left_bound_bis + right_bound_bis) / 2.0;
        if (this->spectral_function(mid) > really_small) {
          right_bound_bis = mid;
        } else {
          left_bound_bis = mid;
        }
      }
      min_mass_spectral_ = right_bound_bis;
    }
  }
  return min_mass_spectral_;
}

int ParticleType::isospin() const {
  if (isospin_ < 0) {
    isospin_ = (pdgcode_.is_hadron() && iso_multiplet_)
                   ? iso_multiplet_->isospin()
                   : 0;
  }
  return isospin_;
}

double ParticleType::partial_width(const double m,
                                   const DecayBranch *mode) const {
  if (m < mode->threshold()) {
    return 0.;
  }
  double partial_width_at_pole = width_at_pole() * mode->weight();
  return mode->type().width(mass(), partial_width_at_pole, m);
}

const DecayModes &ParticleType::decay_modes() const {
  const auto offset = this - std::addressof(list_all()[0]);
  const auto &modes = (*DecayModes::all_decay_modes)[offset];
  assert(is_stable() || !modes.is_empty());
  return modes;
}

double ParticleType::total_width(const double m) const {
  double w = 0.;
  if (is_stable()) {
    return w;
  }
  /* Loop over decay modes and sum up all partial widths. */
  const auto &modes = decay_modes().decay_mode_list();
  for (unsigned int i = 0; i < modes.size(); i++) {
    w = w + partial_width(m, modes[i].get());
  }
  if (w < width_cutoff) {
    return 0.;
  }
  return w;
}

void ParticleType::check_consistency() {
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (!ptype.is_stable() && ptype.decay_modes().is_empty()) {
      throw std::runtime_error("Unstable particle " + ptype.name() +
                               " has no decay chanels!");
    }
  }
}

DecayBranchList ParticleType::get_partial_widths(const double m) const {
  const auto &decay_mode_list = decay_modes().decay_mode_list();
  if (decay_mode_list.size() == 0) {
    return {};
  }

  /* Loop over decay modes and calculate all partial widths. */
  DecayBranchList partial;
  partial.reserve(decay_mode_list.size());
  for (unsigned int i = 0; i < decay_mode_list.size(); i++) {
    const double w = partial_width(m, decay_mode_list[i].get());
    if (w > 0.) {
      partial.push_back(
          make_unique<DecayBranch>(decay_mode_list[i]->type(), w));
    }
  }
  return partial;
}

DecayBranchList ParticleType::get_partial_widths_hadronic(
    const double m) const {
  if (is_stable()) {
    return {};
  }
  /* Loop over decay modes and calculate all partial widths. */
  const auto &decay_mode_list = decay_modes().decay_mode_list();
  DecayBranchList partial;
  partial.reserve(decay_mode_list.size());
  for (unsigned int i = 0; i < decay_mode_list.size(); i++) {
    switch (decay_mode_list[i]->type().particle_number()) {
      case 2: {
        if (!(is_dilepton(
                decay_mode_list[i]->type().particle_types()[0]->pdgcode(),
                decay_mode_list[i]->type().particle_types()[1]->pdgcode()))) {
          const double w = partial_width(m, decay_mode_list[i].get());
          if (w > 0.) {
            partial.push_back(
                make_unique<DecayBranch>(decay_mode_list[i]->type(), w));
          }
        }
        break;
      }
      case 3: {
        if (!(has_lepton_pair(
                decay_mode_list[i]->type().particle_types()[0]->pdgcode(),
                decay_mode_list[i]->type().particle_types()[1]->pdgcode(),
                decay_mode_list[i]->type().particle_types()[2]->pdgcode()))) {
          const double w = partial_width(m, decay_mode_list[i].get());
          if (w > 0.) {
            partial.push_back(
                make_unique<DecayBranch>(decay_mode_list[i]->type(), w));
          }
        }
        break;
      }
      default:
        throw std::runtime_error("Problem in get_partial_widths_hadronic()");
    }
  }
  return partial;
}

DecayBranchList ParticleType::get_partial_widths_dilepton(
    const double m) const {
  const auto &decay_mode_list = decay_modes().decay_mode_list();
  if (decay_mode_list.size() == 0) {
    return {};
  }
  /* Loop over decay modes and calculate all partial widths. */
  DecayBranchList partial;
  partial.reserve(decay_mode_list.size());
  for (unsigned int i = 0; i < decay_mode_list.size(); i++) {
    switch (decay_mode_list[i]->type().particle_number()) {
      case 2: {
        if (is_dilepton(
                decay_mode_list[i]->type().particle_types()[0]->pdgcode(),
                decay_mode_list[i]->type().particle_types()[1]->pdgcode())) {
          const double w = partial_width(m, decay_mode_list[i].get());
          if (w > 0.) {
            partial.push_back(
                make_unique<DecayBranch>(decay_mode_list[i]->type(), w));
          }
        }
        break;
      }
      case 3: {
        if (has_lepton_pair(
                decay_mode_list[i]->type().particle_types()[0]->pdgcode(),
                decay_mode_list[i]->type().particle_types()[1]->pdgcode(),
                decay_mode_list[i]->type().particle_types()[2]->pdgcode())) {
          const double w = partial_width(m, decay_mode_list[i].get());
          if (w > 0.) {
            partial.push_back(
                make_unique<DecayBranch>(decay_mode_list[i]->type(), w));
          }
        }
        break;
      }
      default:
        throw std::runtime_error("Problem in get_partial_widths_dilepton()");
    }
  }
  return partial;
}

double ParticleType::get_partial_width(const double m, const ParticleType &t_a,
                                       const ParticleType &t_b) const {
  /* Get all decay modes. */
  const auto &decaymodes = decay_modes().decay_mode_list();

  /* Find the right one(s) and add up corresponding widths. */
  double w = 0.;
  for (const auto &mode : decaymodes) {
    double partial_width_at_pole = width_at_pole() * mode->weight();
    const ParticleTypePtrList l = {&t_a, &t_b};
    if (mode->type().has_particles(l)) {
      w += mode->type().width(mass(), partial_width_at_pole, m);
    }
  }
  return w;
}

double ParticleType::get_partial_in_width(const double m,
                                          const ParticleData &p_a,
                                          const ParticleData &p_b) const {
  /* Get all decay modes. */
  const auto &decaymodes = decay_modes().decay_mode_list();

  /* Find the right one(s) and add up corresponding widths. */
  double w = 0.;
  for (const auto &mode : decaymodes) {
    double partial_width_at_pole = width_at_pole() * mode->weight();
    const ParticleTypePtrList l = {&p_a.type(), &p_b.type()};
    if (mode->type().has_particles(l)) {
      w += mode->type().in_width(mass(), partial_width_at_pole, m,
                                 p_a.effective_mass(), p_b.effective_mass());
    }
  }
  return w;
}

double ParticleType::spectral_function(double m) const {
  if (norm_factor_ < 0.) {
    /* Initialize the normalization factor
     * by integrating over the unnormalized spectral function. */
    static /*thread_local (see #3075)*/ Integrator integrate;
    const double width = width_at_pole();
    const double m_pole = mass();
    // We transform the integral using m = m_min + width_pole * tan(x), to
    // make it definite and to avoid numerical issues.
    const double x_min = std::atan((min_mass_kinematic() - m_pole) / width);
    norm_factor_ = 1. / integrate(x_min, M_PI / 2., [&](double x) {
                     const double tanx = std::tan(x);
                     const double m_x = m_pole + width * tanx;
                     const double jacobian = width * (1.0 + tanx * tanx);
                     return spectral_function_no_norm(m_x) * jacobian;
                   });
  }
  return norm_factor_ * spectral_function_no_norm(m);
}

double ParticleType::spectral_function_no_norm(double m) const {
  /* The spectral function is a relativistic Breit-Wigner function
   * with mass-dependent width. Here: without normalization factor. */
  const double resonance_width = total_width(m);
  if (resonance_width < ParticleType::width_cutoff) {
    return 0.;
  }
  return breit_wigner(m, mass(), resonance_width);
}

double ParticleType::spectral_function_const_width(double m) const {
  /* The spectral function is a relativistic Breit-Wigner function.
   * This variant is using a constant width (evaluated at the pole mass). */
  const double resonance_width = width_at_pole();
  if (resonance_width < ParticleType::width_cutoff) {
    return 0.;
  }
  return breit_wigner(m, mass(), resonance_width);
}

double ParticleType::spectral_function_simple(double m) const {
  return breit_wigner_nonrel(m, mass(), width_at_pole());
}

/* Resonance mass sampling for 2-particle final state */
double ParticleType::sample_resonance_mass(const double mass_stable,
                                           const double cms_energy,
                                           int L) const {
  /* largest possible mass: Use 'nextafter' to make sure it is not above the
   * physical limit by numerical error. */
  const double max_mass = std::nextafter(cms_energy - mass_stable, 0.);
  // largest possible cm momentum (from smallest mass)
  const double pcm_max =
      pCM(cms_energy, mass_stable, this->min_mass_spectral());
  const double blw_max = pcm_max * blatt_weisskopf_sqr(pcm_max, L);
  /* The maximum of the spectral-function ratio 'usually' happens at the
   * largest mass. However, this is not always the case, therefore we need
   * and additional fudge factor (determined automatically). Additionally,
   * a heuristic knowledge is used that usually such mass exist that
   * spectral_function(m) > spectral_function_simple(m). */
  const double sf_ratio_max =
      std::max(1., this->spectral_function(max_mass) /
                       this->spectral_function_simple(max_mass));

  double mass_res, val;
  // outer loop: repeat if maximum is too small
  do {
    const double q_max = sf_ratio_max * this->max_factor1_;
    const double max = blw_max * q_max;  // maximum value for rejection sampling
    // inner loop: rejection sampling
    do {
      // sample mass from a simple Breit-Wigner (aka Cauchy) distribution
      mass_res = random::cauchy(this->mass(), this->width_at_pole() / 2.,
                                this->min_mass_spectral(), max_mass);
      // determine cm momentum for this case
      const double pcm = pCM(cms_energy, mass_stable, mass_res);
      const double blw = pcm * blatt_weisskopf_sqr(pcm, L);
      // determine ratio of full to simple spectral function
      const double q = this->spectral_function(mass_res) /
                       this->spectral_function_simple(mass_res);
      val = q * blw;
    } while (val < random::uniform(0., max));

    // check that we are using the proper maximum value
    if (val > max) {
      const auto &log = logger<LogArea::Resonances>();
      log.debug("maximum is being increased in sample_resonance_mass: ",
                this->max_factor1_, " ", val / max, " ", this->pdgcode(), " ",
                mass_stable, " ", cms_energy, " ", mass_res);
      this->max_factor1_ *= val / max;
    } else {
      break;  // maximum ok, exit loop
    }
  } while (true);

  return mass_res;
}

/* Resonance mass sampling for 2-particle final state with two resonances. */
std::pair<double, double> ParticleType::sample_resonance_masses(
    const ParticleType &t2, const double cms_energy, int L) const {
  const ParticleType &t1 = *this;
  /* Sample resonance mass from the distribution
   * used for calculating the cross section. */
  const double max_mass_1 =
      std::nextafter(cms_energy - t2.min_mass_spectral(), 0.);
  const double max_mass_2 =
      std::nextafter(cms_energy - t1.min_mass_spectral(), 0.);
  // largest possible cm momentum (from smallest mass)
  const double pcm_max =
      pCM(cms_energy, t1.min_mass_spectral(), t2.min_mass_spectral());
  const double blw_max = pcm_max * blatt_weisskopf_sqr(pcm_max, L);

  double mass_1, mass_2, val;
  // outer loop: repeat if maximum is too small
  do {
    // maximum value for rejection sampling (determined automatically)
    const double max = blw_max * t1.max_factor2_;
    // inner loop: rejection sampling
    do {
      // sample mass from a simple Breit-Wigner (aka Cauchy) distribution
      mass_1 = random::cauchy(t1.mass(), t1.width_at_pole() / 2.,
                              t1.min_mass_spectral(), max_mass_1);
      mass_2 = random::cauchy(t2.mass(), t2.width_at_pole() / 2.,
                              t2.min_mass_spectral(), max_mass_2);
      // determine cm momentum for this case
      const double pcm = pCM(cms_energy, mass_1, mass_2);
      const double blw = pcm * blatt_weisskopf_sqr(pcm, L);
      // determine ratios of full to simple spectral function
      const double q1 =
          t1.spectral_function(mass_1) / t1.spectral_function_simple(mass_1);
      const double q2 =
          t2.spectral_function(mass_2) / t2.spectral_function_simple(mass_2);
      val = q1 * q2 * blw;
    } while (val < random::uniform(0., max));

    if (val > max) {
      const auto &log = logger<LogArea::Resonances>();
      log.debug("maximum is being increased in sample_resonance_masses: ",
                t1.max_factor2_, " ", val / max, " ", t1.pdgcode(), " ",
                t2.pdgcode(), " ", cms_energy, " ", mass_1, " ", mass_2);
      t1.max_factor2_ *= val / max;
    } else {
      break;  // maximum ok, exit loop
    }
  } while (true);

  return {mass_1, mass_2};
}

void ParticleType::dump_width_and_spectral_function() const {
  if (is_stable()) {
    std::stringstream err;
    err << "Particle " << *this << " is stable, so it makes no"
        << " sense to print its spectral function, etc.";
    throw std::runtime_error(err.str());
  }

  double rightmost_pole = 0.0;
  const auto &decaymodes = decay_modes().decay_mode_list();
  for (const auto &mode : decaymodes) {
    double pole_mass_sum = 0.0;
    for (const ParticleTypePtr p : mode->type().particle_types()) {
      pole_mass_sum += p->mass();
    }
    if (pole_mass_sum > rightmost_pole) {
      rightmost_pole = pole_mass_sum;
    }
  }

  std::cout << "# mass m[GeV], width w(m) [GeV],"
            << " spectral function(m^2)*m [GeV^-1] of " << *this << std::endl;
  constexpr double m_step = 0.02;
  const double m_min = min_mass_spectral();
  // An empirical value used to stop the printout. Assumes that spectral
  // function decays at high mass, which is true for all known resonances.
  constexpr double spectral_function_threshold = 8.e-3;
  std::cout << std::fixed << std::setprecision(5);
  for (unsigned int i = 0;; i++) {
    const double m = m_min + m_step * i;
    const double w = total_width(m), sf = spectral_function(m);
    if (m > rightmost_pole * 2 && sf < spectral_function_threshold) {
      break;
    }
    std::cout << m << " " << w << " " << sf << std::endl;
  }
}

std::ostream &operator<<(std::ostream &out, const ParticleType &type) {
  const PdgCode &pdg = type.pdgcode();
  return out << type.name() << std::setfill(' ') << std::right
             << "[ mass:" << field<6> << type.mass()
             << ", width:" << field<6> << type.width_at_pole()
             << ", PDG:" << field<6> << pdg
             << ", charge:" << field<3> << pdg.charge()
             << ", spin:" << field<2> << pdg.spin() << "/2 ]";
}

}  // namespace smash
