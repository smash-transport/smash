/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/particletype.h"

#include "include/lineparser.h"
#include "include/outputroutines.h"
#include "include/pdgcode.h"
#include "include/width.h"
#include "include/decaymodes.h"

#include <algorithm>
#include <assert.h>
#include <map>
#include <vector>

namespace Smash {

namespace {
/// Global pointer to the Particle Type list.
const ParticleTypeList *all_particle_types = nullptr;
}  // unnamed namespace

const ParticleTypeList &ParticleType::list_all() {
  assert(all_particle_types);
  return *all_particle_types;
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

/* Construct an antiparticle name-string from the given name-string for the
 * particle and its PDG code. */
static std::string antiname(std::string name, PdgCode code) {
  std::string basename, charge;

  if (name.find("++") != std::string::npos) {
    basename = name.substr(0, name.length()-2);
    charge = "--";
  } else if (name.find("+") != std::string::npos) {
    basename = name.substr(0, name.length()-1);
    charge = "-";
  } else if (name.find("0") != std::string::npos) {
    basename = name.substr(0, name.length()-1);
    charge = "0";
  } else if (name.find("-") != std::string::npos) {
    basename = name.substr(0, name.length()-1);
    charge = "+";
  } else if (name.find("--") != std::string::npos) {
    basename = name.substr(0, name.length()-2);
    charge = "++";
  } else {
    basename = name;
    charge = "";
  }

  if (code.baryon_number() != 0) {
    return basename+"bar"+charge;  // baron
  } else if (code.charge() != 0) {
    return basename+charge;        // charged meson
  } else {
    return basename+charge+"bar";  // neutral meson
  }
};

void ParticleType::create_type_list(const std::string &input) {  //{{{
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
          "While loading the Particle data:\nFailed to convert the input "
          "string to the expected data types.",
          line));
    }
    ensure_all_read(lineinput, line);

    printd("Setting particle type %s mass %g width %g pdgcode %s\n",
           name.c_str(), mass, width, pdgcode.string().c_str());
    printd("Setting particle type %s isospin %i/2 charge %i spin %i/2\n",
           name.c_str(), pdgcode.isospin_total(), pdgcode.charge(),
                                                  pdgcode.spin());

    type_list.emplace_back(name, mass, width, pdgcode);
    if (pdgcode.has_antiparticle()) {
      /* add corresponding antiparticle */
      PdgCode anti = pdgcode.anti();
      name = antiname (name, pdgcode);
      type_list.emplace_back(name, mass, width, anti);
      printd("Setting antiparticle type %s mass %g width %g pdgcode %s\n",
             name.c_str(), mass, width, anti.string().c_str());
      printd("Setting antiparticle type %s isospin %i/2 charge %i spin %i/2\n",
             name.c_str(), anti.isospin_total(), anti.charge(), anti.spin());
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
                                  const DecayBranch &mode) const{
  if (m < mode.threshold()) {
    return 0.;
  }
  float partial_width_at_pole = width_at_pole()*mode.weight();
  const ParticleType &t1 = ParticleType::find(mode.pdg_list()[0]);
  const ParticleType &t2 = ParticleType::find(mode.pdg_list()[1]);
  if (mode.pdg_list().size()==2 && t1.is_stable() && t2.is_stable()) {
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
  ProcessBranchList partial;
  if (is_stable()) {
    return partial;
  }
  /* Loop over decay modes and calculate all partial widths. */
  for (const auto &mode : DecayModes::find(pdgcode()).decay_mode_list()) {
    w = partial_width(m, mode);
    if (w>0.) {
      partial.push_back(ProcessBranch(mode.pdg_list(),w));
    }
  }
  return partial;
}


}  // namespace Smash
