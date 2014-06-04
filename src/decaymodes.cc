/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/decaymodes.h"

#include "include/constants.h"
#include "include/pdgcode.h"
#include "include/processbranch.h"
#include "include/lineparser.h"

#include <cstdio>
#include <map>
#include <assert.h>

namespace Smash {

using DecayModesMap = std::map<PdgCode, DecayModes>;

namespace {
  /// a map between pdg and corresponding decay modes
  const DecayModesMap *all_decay_modes = nullptr;
}  // unnamed namespace

void DecayModes::add_mode(float ratio, int L, std::vector<PdgCode> pdg_list) {
  if (pdg_list.size() < 2) {
    throw InvalidDecay(
        "DecayModes::add_mode was instructed to add a decay mode with less "
        "than 2 particles. This is an invalid input.");
  }
  DecayBranch branch;
  branch.set_weight(ratio);
  branch.set_angular_momentum(L);
  branch.set_particles(std::move(pdg_list));
  decay_modes_.push_back(branch);
}

void DecayModes::renormalize(float renormalization_constant) {
  if (renormalization_constant < really_small) {
    printf("Warning: Extremely small renormalization constant: %g\n",
           renormalization_constant);
    printf("Skipping the renormalization.\n");
  } else {
    printf("Renormalizing decay modes with %g \n", renormalization_constant);
    float new_sum = 0.0;
    for (std::vector<DecayBranch>::iterator mode
           = decay_modes_.begin(); mode != decay_modes_.end(); ++mode) {
      mode->set_weight(mode->weight() / renormalization_constant);
      new_sum += mode->weight();
    }
    printf("After renormalization sum of ratios is %g. \n", new_sum);
  }
}

/* return the decay modes of specific type */
const DecayModes &DecayModes::find(PdgCode pdg) {
  return all_decay_modes->at(pdg);
}


void DecayModes::load_decaymodes(const std::string &input) {
  static DecayModesMap decaymodes;
  decaymodes.clear();  // in case an exception was thrown and we should try again
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
    // XXX: why not just unconditionally call renormalize? (mkretz)
    /* Check if ratios add to 1 */
    if (fabs(ratio_sum - 1.0) > really_small) {
      /* They didn't; renormalize */
      printf("Particle %s:\n", pdgcode.string().c_str());
      decay_modes_to_add.renormalize(ratio_sum);
    }
    /* Add the list of decay modes for this particle type */
    decaymodes.insert(std::make_pair(pdgcode, decay_modes_to_add));
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
      std::vector<PdgCode> decay_particles;
      decay_particles.reserve(4);
      float ratio;
      lineinput >> ratio;

      int L;
      lineinput >> L;

      PdgCode pdg;
      lineinput >> pdg;
      while (lineinput) {
        if (!ParticleType::exists(pdg)) {
          throw ReferencedParticleNotFound(build_error_string(
              "Inconsistency: The particle with PDG id " +
                  pdg.string() +
                  " was not registered through particles.txt, but "
                  "decaymodes.txt referenced it.",
              line));
        }
        decay_particles.push_back(pdg);
        lineinput >> pdg;
      }
      if (pdg != PdgCode::invalid()) {
        decay_particles.push_back(pdg);
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
