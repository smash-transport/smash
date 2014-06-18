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
  assert(found != all_particle_types->end());
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


float ParticleType::width_total(const float m) const {
  float w = 0.;
  const ProcessBranchList partial = width_partial(m);
  // loop over decay modes and sum up all partial widths
  for (const auto &mode : partial) {
    w = w + mode.weight();
  }
  return w;
}


ProcessBranchList ParticleType::width_partial(const float m) const {
  float w = 0.;
  ProcessBranchList partial;
  if (is_stable()) {
    return partial;
  }
  const DecayBranchList decaymodes
        = DecayModes::find(pdgcode()).decay_mode_list();
  // loop over decay modes and sum up all partial widths
  for (const auto &mode : decaymodes) {
    if (m < mode.threshold()) {
      continue; }
    float partial_width_at_pole = width_at_pole()*mode.weight();
    const ParticleType &t1 = ParticleType::find(mode.pdg_list()[0]);
    const ParticleType &t2 = ParticleType::find(mode.pdg_list()[1]);
    if (mode.pdg_list().size()==2 && t1.is_stable() && t2.is_stable()) {
      // mass-dependent width for 2-body decays
      w = width_Manley_stable(m, mass(), t1.mass(), t2.mass(),
                              mode.angular_momentum(), partial_width_at_pole);
    }
    else {
      // constant width for three-body decays
      w = partial_width_at_pole;
    }
    partial.push_back(ProcessBranch(mode.pdg_list(),w));
  }

  return partial;
}


}  // namespace Smash
