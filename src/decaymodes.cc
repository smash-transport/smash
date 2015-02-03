/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/decaymodes.h"

#include <assert.h>
#include <cstdio>
#include <map>

#include "include/constants.h"
#include "include/inputfunctions.h"
#include "include/logging.h"
#include "include/pdgcode.h"
#include "include/processbranch.h"
#include "include/stringfunctions.h"

namespace Smash {

const std::vector<DecayModes> *DecayModes::all_decay_modes = nullptr;

void DecayModes::add_mode(float ratio, int L,
                          std::vector<ParticleTypePtr> particle_types) {
  switch (particle_types.size()) {
  case 2:
    if (!particle_types[0]->is_hadron() || !particle_types[1]->is_hadron()) {
      logger<LogArea::DecayModes>().warn(
          "decay products A: ", *particle_types[0], " B: ", *particle_types[1]);
    }
    break;
  case 3:
    if (!particle_types[0]->is_hadron() || !particle_types[1]->is_hadron() ||
        !particle_types[2]->is_hadron()) {
      logger<LogArea::DecayModes>().warn(
          "decay products A: ", *particle_types[0], " B: ", *particle_types[1],
          " C: ", *particle_types[2]);
    }
    break;
  default:
    throw InvalidDecay(
        "DecayModes::add_mode was instructed to add a decay mode with " +
        std::to_string(particle_types.size()) +
        " particles. This is an invalid input.");
  }
  decay_modes_.emplace_back(L, std::move(particle_types), ratio);
}

void DecayModes::renormalize(float renormalization_constant) {
  const auto &log = logger<LogArea::DecayModes>();
  if (renormalization_constant < really_small) {
    log.warn("Extremely small renormalization constant: ",
             renormalization_constant,
             "\n=> Skipping the renormalization.");
  } else {
    log.info("Renormalizing decay modes with ", renormalization_constant);
    float new_sum = 0.0;
    for (auto &mode : decay_modes_) {
      mode.set_weight(mode.weight() / renormalization_constant);
      new_sum += mode.weight();
    }
    log.info("After renormalization sum of ratios is ", new_sum);
  }
}

namespace {
inline std::size_t find_offset(PdgCode pdg) {
  return std::addressof(ParticleType::find(pdg)) -
         std::addressof(ParticleType::list_all()[0]);
}
}  // unnamed namespace

void DecayModes::load_decaymodes(const std::string &input) {
  static std::vector<DecayModes> decaymodes;
  decaymodes.clear();  // in case an exception was thrown and should try again
  decaymodes.resize(ParticleType::list_all().size());

  PdgCode pdgcode = PdgCode::invalid();
  DecayModes decay_modes_to_add;
  float ratio_sum = 0.0;

  const auto end_of_decaymodes = [&]() {
    if (pdgcode == PdgCode::invalid()) {  // at the start of the file
      return;
    }
    if (decay_modes_to_add.is_empty()) {
      throw MissingDecays("No decay modes found for particle " +
                          pdgcode.string());
    }
    // TODO: why not just unconditionally call renormalize? (mkretz)
    /* Check if ratios add to 1 */
    if (fabs(ratio_sum - 1.0) > really_small) {
      /* They didn't; renormalize */
      logger<LogArea::DecayModes>().info("Particle ", pdgcode);
      decay_modes_to_add.renormalize(ratio_sum);
    }

    if (pdgcode.has_antiparticle()) {
      /* Construct and add the list of decay modes for the antiparticle.  */
      DecayModes decay_modes_anti;
      for (const auto &mode : decay_modes_to_add.decay_mode_list()) {
        std::vector<ParticleTypePtr> list = mode.particle_types();
        for (auto &type : list) {
          if (type->has_antiparticle()) {
            type = type->get_antiparticle();
          }
        }
        decay_modes_anti.add_mode(mode.weight(), mode.angular_momentum(),
                                  list);
      }
      decaymodes[find_offset(pdgcode.get_antiparticle())] =
          std::move(decay_modes_anti);
    }

    /* Add the list of decay modes for this particle type */
    decaymodes[find_offset(pdgcode)] = std::move(decay_modes_to_add);

    /* Clean up the list for the next particle type */
    decay_modes_to_add.clear();
    ratio_sum = 0.0;
  };

  for (const Line &line : line_parser(input)) {
    const auto trimmed = trim(line.text);
    assert(!trimmed.empty());  // trim(line.text) is never empty - else
                               // line_parser is broken
    // if (trimmed.find_first_not_of("-0123456789") ==
    if (trimmed.find_first_of(" \t") ==
        std::string::npos) {  // a single record on one line signifies a new
                              // decay mode section
      end_of_decaymodes();
      pdgcode = PdgCode(trim(line.text));
      if (!ParticleType::exists(pdgcode)) {
        throw ReferencedParticleNotFound(build_error_string(
            "Inconsistency: The particle with PDG id " +
                pdgcode.string() +
                " was not registered through particles.txt, but "
                "decaymodes.txt referenced it.",
            line));
      }
      assert(pdgcode != PdgCode::invalid());  // special value for start of file
    } else {
      std::istringstream lineinput(line.text);
      std::vector<ParticleTypePtr> decay_particles;
      decay_particles.reserve(4);
      float ratio;
      lineinput >> ratio;

      int L;
      lineinput >> L;

      PdgCode pdg;
      lineinput >> pdg;
      while (lineinput) {
        try {
          decay_particles.emplace_back(&ParticleType::find(pdg));
          lineinput >> pdg;
        }
        catch ( ParticleType::PdgNotFoundFailure ) {
          throw ReferencedParticleNotFound(build_error_string(
              "Inconsistency: The particle with PDG id " + pdg.string() +
                  " was not registered through particles.txt, but "
                  "decaymodes.txt referenced it.",
              line));
        }
      }
      if (pdg != PdgCode::invalid()) {
        decay_particles.emplace_back(&ParticleType::find(pdg));
      }
      if (lineinput.fail() && !lineinput.eof()) {
        throw LoadFailure(
            build_error_string("Parse error: expected a PdgCode ", line));
      }
      decay_particles.shrink_to_fit();
      decay_modes_to_add.add_mode(ratio, L, std::move(decay_particles));
      ratio_sum += ratio;
    }
  }
  end_of_decaymodes();
  assert(nullptr == all_decay_modes);
  all_decay_modes = &decaymodes;
}


}  // namespace Smash
