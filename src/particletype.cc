/*
 *
 *    Copyright (c) 2014
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

#include "include/decaymodes.h"
#include "include/inputfunctions.h"
#include "include/iomanipulators.h"
#include "include/logging.h"
#include "include/outputroutines.h"
#include "include/pdgcode.h"
#include "include/stringfunctions.h"
#include "include/width.h"

namespace Smash {

namespace {
/// Global pointer to the Particle Type list.
const ParticleTypeList *all_particle_types = nullptr;
}  // unnamed namespace

const ParticleTypeList &ParticleType::list_all() {
  assert(all_particle_types);
  return *all_particle_types;
}

const ParticleTypeList ParticleType::list_nucleons() {
  return {find(0x2212), find(0x2112)};
}

const ParticleTypeList ParticleType::list_baryon_resonances() {
  ParticleTypeList list;

  for (const ParticleType &type_resonance : ParticleType::list_all()) {
    /* Only loop over baryon resonances. */
    if (type_resonance.is_stable()
        || type_resonance.pdgcode().baryon_number() != 1) {
      continue;
    }
   list.push_back(type_resonance);
  }
  return list;
}

SMASH_CONST const ParticleType &ParticleType::find(PdgCode pdgcode) {
  const auto found = std::lower_bound(
      all_particle_types->begin(), all_particle_types->end(), pdgcode,
      [](const ParticleType &l, const PdgCode &r) { return l.pdgcode() < r; });
  if (found == all_particle_types->end()) {
    throw std::runtime_error("PDG code "+pdgcode.string()+ " not found!");
  }
  assert(found->pdgcode() == pdgcode);
  return *found;
}

SMASH_CONST bool ParticleType::exists(PdgCode pdgcode) {
  const auto found = std::lower_bound(
      all_particle_types->begin(), all_particle_types->end(), pdgcode,
      [](const ParticleType &l, const PdgCode &r) { return l.pdgcode() < r; });
  if (found != all_particle_types->end()) {
    return found->pdgcode() == pdgcode;
  }
  return false;
}

#ifdef NDEBUG
ParticleType::ParticleType(std::string, float m, float w, PdgCode id)
    :
#else
ParticleType::ParticleType(std::string n, float m, float w, PdgCode id)
    : name_(fill_right(n, 3)),
#endif
      mass_(m),
      width_(w),
      pdgcode_(id),
      isospin_(pdgcode_.isospin_total()),
      charge_(pdgcode_.charge()) {}

/* Construct an antiparticle name-string from the given name-string for the
 * particle and its PDG code. */
static std::string antiname(const std::string &name, PdgCode code) {
  std::string basename, charge;

  if (name.find("++") != std::string::npos) {
    basename = name.substr(0, name.length()-2);
    charge = "--";
  } else if (name.find("⁺⁺") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁺⁺") + 1);
    charge = "⁻⁻";
  } else if (name.find("+") != std::string::npos) {
    basename = name.substr(0, name.length()-1);
    charge = "-";
  } else if (name.find("⁺") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁺") + 1);
    charge = "⁻";
  } else if (name.find("0") != std::string::npos) {
    basename = name.substr(0, name.length()-1);
    charge = "0";
  } else if (name.find("⁰") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁰") + 1);
    charge = "⁰";
  } else if (name.find("-") != std::string::npos) {
    basename = name.substr(0, name.length()-1);
    charge = "+";
  } else if (name.find("⁻") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁻") + 1);
    charge = "⁺";
  } else if (name.find("--") != std::string::npos) {
    basename = name.substr(0, name.length()-2);
    charge = "++";
  } else if (name.find("⁻⁻") != std::string::npos) {
    basename = name.substr(0, name.length() - sizeof("⁻⁻") + 1);
    charge = "⁺⁺";
  } else {
    basename = name;
    charge = "";
  }

  constexpr char bar[] = "\u0305";
  if (code.baryon_number() != 0) {
    return basename+bar+charge;  // baryon
  } else if (code.charge() != 0) {
    return basename+charge;        // charged meson
  } else {
    return basename+bar+charge;  // neutral meson
  }
};

void ParticleType::create_type_list(const std::string &input) {  //{{{
  const auto &log = logger<LogArea::ParticleType>();
  static ParticleTypeList type_list;
  type_list.clear();  // in case LoadFailure was thrown and caught and we should
                      // try again
  for (const Line &line : line_parser(input)) {
    std::istringstream lineinput(line.text);
    std::string name;
    float mass, width;
    PdgCode pdgcode;
    lineinput >> name >> mass >> width >> pdgcode;
    if (lineinput.fail()) {
      throw Particles::LoadFailure(build_error_string(
          "While loading the ParticleType data:\nFailed to convert the input "
          "string to the expected data types.",
          line));
    }
    ensure_all_read(lineinput, line);

    type_list.emplace_back(name, mass, width, pdgcode);
    log.debug() << "Setting     particle type: " << type_list.back();
    if (pdgcode.has_antiparticle()) {
      /* add corresponding antiparticle */
      PdgCode anti = pdgcode.get_antiparticle();
      name = antiname(name, pdgcode);
      type_list.emplace_back(name, mass, width, anti);
      log.debug() << "Setting antiparticle type: " << type_list.back();
    }
  }
  type_list.shrink_to_fit();

  std::sort(type_list.begin(), type_list.end(),
            [](const ParticleType &l,
               const ParticleType &r) { return l.pdgcode() < r.pdgcode(); });

  assert(nullptr == all_particle_types);
  all_particle_types = &type_list;  // note that type_list is a function-local
                                    // static and thus will live on until after
                                    // main().
}/*}}}*/


float ParticleType::minimum_mass() const {
  float minmass = mass();

  /* If the particle happens to be stable, just return the mass. */
  if (is_stable()) {
    return minmass;
  }
  /* Otherwise, find the lowest mass value needed in any decay mode */
  for (const auto &mode : DecayModes::find(pdgcode()).decay_mode_list()) {
    minmass = std::min(minmass, mode.threshold());
  }
  return minmass;
}


float ParticleType::partial_width(const float m,
                                  const DecayBranch &mode) const {
  if (m < mode.threshold()) {
    return 0.;
  }
  float partial_width_at_pole = width_at_pole()*mode.weight();
  const ParticleType &t1 = ParticleType::find(mode.pdg_list()[0]);
  const ParticleType &t2 = ParticleType::find(mode.pdg_list()[1]);
  if (mode.pdg_list().size() == 2 && t1.is_stable() && t2.is_stable()) {
    /* mass-dependent width for 2-body decays with stable decay products */
    return width_Manley_stable(m, mass(), t1.mass(), t2.mass(),
                               mode.angular_momentum(),
                               partial_width_at_pole);
  } else {
    /* constant width for three-body decays
     * (and two-body decays with unstable products) */
    return partial_width_at_pole;
  }
}


float ParticleType::total_width(const float m) const {
  float w = 0.;
  if (is_stable()) {
    return w;
  }
  /* Loop over decay modes and sum up all partial widths. */
  for (const auto &mode : DecayModes::find(pdgcode()).decay_mode_list()) {
    w = w + partial_width(m, mode);
  }
  return w;
}


ProcessBranchList ParticleType::get_partial_widths(const float m) const {
  float w = 0.;
  if (is_stable()) {
    return {};
  }
  /* Loop over decay modes and calculate all partial widths. */
  const auto &decay_mode_list = DecayModes::find(pdgcode()).decay_mode_list();
  ProcessBranchList partial;
  partial.reserve(decay_mode_list.size());
  for (const auto &mode : decay_mode_list) {
    w = partial_width(m, mode);
    if (w > 0.) {
      partial.emplace_back(mode.pdg_list(), w);
    }
  }
  return std::move(partial);
}

std::ostream &operator<<(std::ostream &out, const ParticleType &type) {
  const PdgCode &pdg = type.pdgcode();
  return out << type.name() << std::setfill(' ') << std::right
             << "[mass:" << field<6> << type.mass()
             << ", width:" << field<6> << type.width_at_pole()
             << ", PDG:" << field<6> << pdg
             << ", Isospin:" << field<2> << pdg.isospin_total()
             << "/2, Charge:" << field<3> << pdg.charge()
             << ", Spin:" << field<2> << pdg.spin() << "/2]";
}

}  // namespace Smash
