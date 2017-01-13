/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/particletype.h"

#include <assert.h>
#include <algorithm>
#include <map>
#include <vector>

#include "include/cxx14compat.h"
#include "include/decaymodes.h"
#include "include/distributions.h"
#include "include/formfactors.h"
#include "include/inputfunctions.h"
#include "include/integrate.h"
#include "include/iomanipulators.h"
#include "include/isoparticletype.h"
#include "include/kinematics.h"
#include "include/logging.h"
#include "include/numerics.h"
#include "include/particledata.h"
#include "include/pdgcode.h"
#include "include/pow.h"
#include "include/processbranch.h"
#include "include/stringfunctions.h"

namespace Smash {

#ifdef SMASH_INLINE_LIST_ALL
const ParticleTypeList *all_particle_types = nullptr;
#else
namespace {
/// Global pointer to the Particle Type list.
const ParticleTypeList *all_particle_types = nullptr;
ParticleTypePtrList nucleons_list;
ParticleTypePtrList deltas_list;
ParticleTypePtrList baryon_resonances_list;
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

ParticleTypePtrList &ParticleType::list_nucleons() {
  return nucleons_list;
}

ParticleTypePtrList &ParticleType::list_Deltas() {
  return deltas_list;
}

ParticleTypePtrList &ParticleType::list_baryon_resonances() {
  return baryon_resonances_list;
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

ParticleType::ParticleType(std::string n, float m, float w, PdgCode id)
    : name_(n),
      mass_(m),
      width_(w),
      pdgcode_(id),
      minimum_mass_(-1.f),
      charge_(pdgcode_.charge()),
      isospin_(-1),
      I3_(pdgcode_.isospin3()) {}

/* Construct an antiparticle name-string from the given name-string for the
 * particle and its PDG code. */
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

  return basename+charge;
}

/* Construct a charge string, given the charge as integer. */
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
    float mass, width;
    std::array<PdgCode, 4> pdgcode;
    lineinput >> name >> mass >> width >> pdgcode[0];
    if (lineinput.fail()) {
      throw ParticleType::LoadFailure(build_error_string(
          "While loading the ParticleType data:\nFailed to convert the input "
          "string to the expected data types.", line));
    }
    // read additional PDG codes (if present)
    unsigned int n = 1;  // number of PDG codes found
    while (!lineinput.eof() && n < pdgcode.size()) {
      lineinput >> pdgcode[n++];
      if (lineinput.fail()) {
        throw ParticleType::LoadFailure(build_error_string(
            "While loading the ParticleType data:\nFailed to convert the input "
            "string to the expected data types.", line));
      }
    }
    ensure_all_read(lineinput, line);

    //Check if nucleon, kaon, and delta masses are the same as hardcoded ones, if present
    if (pdgcode[0].is_nucleon() && !almost_equal(mass, nucleon_mass)) {
      throw std::runtime_error("Nucleon mass in input file different from 0.938");
    }
    if (pdgcode[0].is_kaon() && !almost_equal(mass, kaon_mass)) {
      throw std::runtime_error("Kaon mass in input file different from 0.494"); 
    }
    if (pdgcode[0].is_Delta() && !almost_equal(mass, delta_mass)) {
      throw std::runtime_error("Delta mass in input file different from 1.232"); 
    }

    // add all states to type list
    for (unsigned int i = 0; i < n; i++) {
      std::string full_name = name;
      if (n > 1) {
        // for multiplets: add charge string to name
        full_name += chargestr(pdgcode[i].charge());
      }
      type_list.emplace_back(full_name, mass, width, pdgcode[i]);
      log.debug() << "Setting     particle type: " << type_list.back();
      if (pdgcode[i].has_antiparticle()) {
        /* add corresponding antiparticle */
        PdgCode anti = pdgcode[i].get_antiparticle();
        full_name = antiname(full_name, pdgcode[i]);
        type_list.emplace_back(full_name, mass, width, anti);
        log.debug() << "Setting antiparticle type: " << type_list.back();
      }
    }
  }
  type_list.shrink_to_fit();

  /* Sort the type list by PDG code. */
  std::sort(type_list.begin(), type_list.end());

  /* Look for duplicates. */
  PdgCode prev_pdg = 0;
  for (const auto& t : type_list) {
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

  // Create nucleons list
  if (IsoParticleType::exists("N")) {
    for (const auto state : IsoParticleType::find("N").get_states()) {
      nucleons_list.push_back(state);
    }
  }

  // Create deltas list
  if (IsoParticleType::exists("Δ")) {
    for (const auto state : IsoParticleType::find("Δ").get_states()) {
      deltas_list.push_back(state);
    }
  }

  // Create baryon resonances list
  for (const ParticleType &type_resonance : ParticleType::list_all()) {
    /* Only loop over baryon resonances. */
    if (type_resonance.is_stable()
        || type_resonance.pdgcode().baryon_number() != 1) {
      continue;
    }
    baryon_resonances_list.push_back(&type_resonance);
  }

}/*}}}*/


float ParticleType::minimum_mass() const {
  if (unlikely(minimum_mass_ < 0.f)) {
    /* If the particle is stable, min. mass is just the mass. */
    minimum_mass_ = mass_;
    /* Otherwise, find the lowest mass value needed in any decay mode */
    if (!is_stable()) {
      for (const auto &mode : decay_modes().decay_mode_list()) {
        minimum_mass_ = std::min(minimum_mass_, mode->threshold());
      }
    }
  }
  return minimum_mass_;
}

int ParticleType::isospin() const {
  if (isospin_ < 0) {
    isospin_ =  (pdgcode_.is_hadron() && iso_multiplet_) ?
                iso_multiplet_->isospin() : 0;
  }
  return isospin_;
}

float ParticleType::partial_width(const float m,
                                  const DecayBranch *mode) const {
  if (m < mode->threshold()) {
    return 0.;
  }
  float mode_weight = mode->weight();
  float wap = width_at_pole();
  float partial_width_at_pole = wap*mode_weight;
  return mode->type().width(mass(), partial_width_at_pole, m);
}

const DecayModes &ParticleType::decay_modes() const {
  const auto offset = this - std::addressof(list_all()[0]);
  const auto &modes = (*DecayModes::all_decay_modes)[offset];
  assert(is_stable() || !modes.is_empty());
  return modes;
}

float ParticleType::total_width(const float m) const {
  float w = 0.;
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

DecayBranchList ParticleType::get_partial_widths(const float m) const {
  const auto &decay_mode_list = decay_modes().decay_mode_list();
  if (decay_mode_list.size() == 0) {
    return {};
  }

  /* Loop over decay modes and calculate all partial widths. */
  DecayBranchList partial;
  partial.reserve(decay_mode_list.size());
  for (unsigned int i = 0; i < decay_mode_list.size(); i++) {
    const float w = partial_width(m, decay_mode_list[i].get());
    if (w > 0.) {
      partial.push_back(
          make_unique<DecayBranch>(decay_mode_list[i]->type(), w));
    }
  }
  return partial;
}

DecayBranchList ParticleType::get_partial_widths_hadronic(const float m) const {
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
          const float w = partial_width(m, decay_mode_list[i].get());
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
          const float w = partial_width(m, decay_mode_list[i].get());
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

DecayBranchList ParticleType::get_partial_widths_dilepton(const float m) const {
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
          const float w = partial_width(m, decay_mode_list[i].get());
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
          const float w = partial_width(m, decay_mode_list[i].get());
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

float ParticleType::get_partial_width(const float m,
                                      const ParticleType &t_a,
                                      const ParticleType &t_b) const {
  /* Get all decay modes. */
  const auto &decaymodes = decay_modes().decay_mode_list();

  /* Find the right one(s) and add up corresponding widths. */
  float w = 0.;
  for (const auto &mode : decaymodes) {
    float partial_width_at_pole = width_at_pole()*mode->weight();
    const ParticleTypePtrList l = {&t_a, &t_b};
    if (mode->type().has_particles(l)) {
      w += mode->type().width(mass(), partial_width_at_pole, m);
    }
  }
  return w;
}

float ParticleType::get_partial_in_width(const float m,
                                         const ParticleData &p_a,
                                         const ParticleData &p_b) const {
  /* Get all decay modes. */
  const auto &decaymodes = decay_modes().decay_mode_list();

  /* Find the right one(s) and add up corresponding widths. */
  float w = 0.;
  for (const auto &mode : decaymodes) {
    float partial_width_at_pole = width_at_pole()*mode->weight();
    const ParticleTypePtrList l = {&p_a.type(), &p_b.type()};
    if (mode->type().has_particles(l)) {
      w += mode->type().in_width(mass(), partial_width_at_pole, m,
                                 p_a.effective_mass(), p_b.effective_mass());
    }
  }
  return w;
}


float ParticleType::spectral_function(float m) const {
  if (norm_factor_ < 0.) {
    /* Initialize the normalization factor
     * by integrating over the unnormalized spectral function. */
    Integrator integrate;
    //^ This should be static, but for some reason then the integrals sometimes
    //  yield different results. See #4299.
    const auto width = width_at_pole();
    // We transform the integral using m = m_min + width_pole * tan(x), to
    // make it definite and to avoid numerical issues.
    norm_factor_ = 1./integrate(std::atan((minimum_mass() - mass())/width), M_PI/2.,
        [&](double x) {
          return spectral_function_no_norm(mass() + width*std::tan(x)) * width
                 * (1 + square(std::tan(x))) ;
    });
  }
  return norm_factor_ * spectral_function_no_norm(m);
}

float ParticleType::spectral_function_no_norm(float m) const {
  /* The spectral function is a relativistic Breit-Wigner function
   * with mass-dependent width. Here: without normalization factor. */
  const float resonance_width = total_width(m);
  if (resonance_width < ParticleType::width_cutoff) {
    return 0.;
  }
  return breit_wigner(m, mass(), resonance_width);
}

float ParticleType::spectral_function_const_width(float m) const {
  /* The spectral function is a relativistic Breit-Wigner function.
   * This variant is using a constant width (evaluated at the pole mass). */
  const float resonance_width = width_at_pole();
  if (resonance_width < ParticleType::width_cutoff) {
    return 0.;
  }
  return breit_wigner(m, mass(), resonance_width);
}

float ParticleType::spectral_function_simple(float m) const {
  return breit_wigner_nonrel(m, mass(), width_at_pole());
}


/* Resonance mass sampling for 2-particle final state */
float ParticleType::sample_resonance_mass(const float mass_stable,
                                          const float cms_energy, int L) const {
  /* largest possible mass: Use 'nextafter' to make sure it is not above the
   * physical limit by numerical error. */
  const float max_mass = std::nextafter(cms_energy - mass_stable, 0.f);
  // largest possible cm momentum (from smallest mass)
  const float pcm_max = pCM(cms_energy, mass_stable, this->minimum_mass());
  const float blw_max = pcm_max * blatt_weisskopf_sqr(pcm_max, L);

  float mass_res, val;
  // outer loop: repeat if maximum is too small
  do {
    /* The maximum of the spectral-function ratio 'usually' happens at the
     * largest mass. However, this is not always the case, therefore we need
     * and additional fudge factor (determined automatically). */
    const float q_max = this->spectral_function(max_mass)
                      / this->spectral_function_simple(max_mass)
                      * this->max_factor1_;
    const float max = blw_max * q_max;  // maximum value for rejection sampling
    // inner loop: rejection sampling
    do {
      // sample mass from a simple Breit-Wigner (aka Cauchy) distribution
      mass_res = Random::cauchy(this->mass(), this->width_at_pole()/2.f,
                                this->minimum_mass(), max_mass);
      // determine cm momentum for this case
      const float pcm = pCM(cms_energy, mass_stable, mass_res);
      const float blw = pcm * blatt_weisskopf_sqr(pcm, L);
      // determine ratio of full to simple spectral function
      const float q = this->spectral_function(mass_res)
                    / this->spectral_function_simple(mass_res);
      val = q * blw;
    } while (val < Random::uniform(0.f, max));

    // check that we are using the proper maximum value
    if (val > max) {
      const auto &log = logger<LogArea::Resonances>();
      log.debug("maximum is being increased in sample_resonance_mass: ",
                this->max_factor1_, " ", val/max, " ", this->pdgcode(),
                " ", mass_stable, " ", cms_energy, " ", mass_res);
      this->max_factor1_ *= val/max;
    } else {
      break;  // maximum ok, exit loop
    }
  } while (true);

  return mass_res;
}


/* Resonance mass sampling for 2-particle final state with two resonances. */
std::pair<float, float> ParticleType::sample_resonance_masses(
                  const ParticleType &t2, const float cms_energy, int L) const {
  const ParticleType &t1 = *this;
  /* Sample resonance mass from the distribution
   * used for calculating the cross section. */
  const float max_mass_1 = std::nextafter(cms_energy - t2.minimum_mass(), 0.f);
  const float max_mass_2 = std::nextafter(cms_energy - t1.minimum_mass(), 0.f);
  // largest possible cm momentum (from smallest mass)
  const float pcm_max = pCM(cms_energy, t1.minimum_mass(), t2.minimum_mass());
  const float blw_max = pcm_max * blatt_weisskopf_sqr(pcm_max, L);

  float mass_1, mass_2, val;
  // outer loop: repeat if maximum is too small
  do {
    // maximum value for rejection sampling (determined automatically)
    const float max = blw_max * t1.max_factor2_;
    // inner loop: rejection sampling
    do {
      // sample mass from a simple Breit-Wigner (aka Cauchy) distribution
      mass_1 = Random::cauchy(t1.mass(), t1.width_at_pole()/2.f,
                              t1.minimum_mass(), max_mass_1);
      mass_2 = Random::cauchy(t2.mass(), t2.width_at_pole()/2.f,
                              t2.minimum_mass(), max_mass_2);
      // determine cm momentum for this case
      const float pcm = pCM(cms_energy, mass_1, mass_2);
      const float blw = pcm * blatt_weisskopf_sqr(pcm, L);
      // determine ratios of full to simple spectral function
      const float q1 = t1.spectral_function(mass_1)
                    / t1.spectral_function_simple(mass_1);
      const float q2 = t2.spectral_function(mass_2)
                    / t2.spectral_function_simple(mass_2);
      val = q1 * q2 * blw;
    } while (val < Random::uniform(0.f, max));

    if (val > max) {
      const auto &log = logger<LogArea::Resonances>();
      log.debug("maximum is being increased in sample_resonance_masses: ",
                t1.max_factor2_, " ", val/max, " ", t1.pdgcode(), " ",
                t2.pdgcode(), " ", cms_energy, " ", mass_1, " ", mass_2);
      t1.max_factor2_ *= val/max;
    } else {
      break;  // maximum ok, exit loop
    }
  } while (true);

  return {mass_1, mass_2};
}


std::ostream &operator<<(std::ostream &out, const ParticleType &type) {
  const PdgCode &pdg = type.pdgcode();
  return out << type.name() << std::setfill(' ') << std::right
             << "[ mass:"   << field<6> << type.mass()
             << ", width:"  << field<6> << type.width_at_pole()
             << ", PDG:"    << field<6> << pdg
             << ", charge:" << field<3> << pdg.charge()
             << ", spin:"   << field<2> << pdg.spin() << "/2 ]";
}

}  // namespace Smash
