/*
 *
 *    Copyright (c) 2014-2015
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
#include "include/cxx14compat.h"
#include "include/inputfunctions.h"
#include "include/isoparticletype.h"
#include "include/logging.h"
#include "include/pdgcode.h"
#include "include/processbranch.h"
#include "include/resonances.h"
#include "include/stringfunctions.h"

namespace Smash {

std::vector<DecayModes> *DecayModes::all_decay_modes = nullptr;

std::vector<std::unique_ptr<DecayType>> *all_decay_types = nullptr;

void DecayModes::add_mode(float ratio, int L,
                          ParticleTypePtrList particle_types) {
  const auto &log = logger<LogArea::DecayModes>();
  assert(all_decay_types != nullptr);
  switch (particle_types.size()) {
  case 2:
    if (!particle_types[0]->is_hadron() || !particle_types[1]->is_hadron()) {
      log.warn("decay products A: ", *particle_types[0],
               " B: ", *particle_types[1]);
    }
    if (particle_types[0]->is_stable() && particle_types[1]->is_stable()) {
      all_decay_types->emplace_back(
          make_unique<TwoBodyDecayStable>(particle_types, L));
    } else if (particle_types[0]->is_stable() ||
               particle_types[1]->is_stable()) {
      all_decay_types->emplace_back(
          make_unique<TwoBodyDecaySemistable>(particle_types, L));
    } else {
      all_decay_types->emplace_back(
          make_unique<TwoBodyDecayUnstable>(particle_types, L));
    }
    break;
  case 3:
    if (!particle_types[0]->is_hadron() || !particle_types[1]->is_hadron() ||
        !particle_types[2]->is_hadron()) {
      log.warn("decay products A: ", *particle_types[0],
               " B: ", *particle_types[1], " C: ", *particle_types[2]);
    }
    all_decay_types->emplace_back(
        make_unique<ThreeBodyDecay>(particle_types, L));
    break;
  default:
    throw InvalidDecay(
        "DecayModes::add_mode was instructed to add a decay mode with " +
        std::to_string(particle_types.size()) +
        " particles. This is an invalid input.");
  }
  decay_modes_.push_back(
      make_unique<DecayBranch>(*all_decay_types->back(), ratio));
}

void DecayModes::renormalize(std::string name) {
  const auto &log = logger<LogArea::DecayModes>();
  float sum = 0.;
  for (auto &mode : decay_modes_) {
    sum += mode->weight();
  }
  if (std::abs(sum - 1.) < really_small) {
    log.debug("Particle ", name, ": Extremely small renormalization constant: ",
              sum, "\n=> Skipping the renormalization.");
  } else {
    log.warn("Particle ", name, ": Renormalizing decay modes with ", sum);
    float new_sum = 0.0;
    for (auto &mode : decay_modes_) {
      mode->set_weight(mode->weight() / sum);
      new_sum += mode->weight();
    }
    log.debug("After renormalization sum of ratios is ", new_sum);
  }
}

namespace {
inline std::size_t find_offset(PdgCode pdg) {
  return std::addressof(ParticleType::find(pdg)) -
         std::addressof(ParticleType::list_all()[0]);
}
}  // unnamed namespace

void DecayModes::load_decaymodes(const std::string &input) {
  const auto &log = logger<LogArea::DecayModes>();
  // create the DecayType vector first, then it outlives the DecayModes vector,
  // which references the DecayType objects.
  static std::vector<std::unique_ptr<DecayType>> decaytypes;
  decaytypes.clear();  // in case an exception was thrown and should try again
  // ten decay types per decay mode should be a good guess.
  decaytypes.reserve(10 * ParticleType::list_all().size());
  all_decay_types = &decaytypes;

  static std::vector<DecayModes> decaymodes;
  decaymodes.clear();  // in case an exception was thrown and should try again
  decaymodes.resize(ParticleType::list_all().size());
  all_decay_modes = &decaymodes;

  const IsoParticleType *isotype_mother = NULL;
  ParticleTypePtrList mother_states;
  std::vector<DecayModes> decay_modes_to_add;  // one for each mother state

  const auto end_of_decaymodes = [&]() {
    if (isotype_mother == NULL) {  // at the start of the file
      return;
    }
    // Loop over all states in the mother multiplet and add modes
    for (unsigned int m = 0; m < mother_states.size(); m++) {
      if (decay_modes_to_add[m].is_empty()) {
        throw MissingDecays("No decay modes found for particle " +
                            mother_states[m]->name());
      }
      decay_modes_to_add[m].renormalize(mother_states[m]->name());
      PdgCode pdgcode = mother_states[m]->pdgcode();
      /* Add the list of decay modes for this particle type */
      decaymodes[find_offset(pdgcode)] = std::move(decay_modes_to_add[m]);
    }
    if (isotype_mother->has_anti_multiplet()) {
      /* Construct the decay modes for the anti-multiplet.  */
      log.debug("generating decay modes for anti-multiplet: " +
                isotype_mother->name());
      for (unsigned int m = 0; m < mother_states.size(); m++) {
        PdgCode pdg = mother_states[m]->pdgcode();
        PdgCode pdg_anti = pdg.get_antiparticle();
        DecayModes &decay_modes_orig = decaymodes[find_offset(pdg)];
        DecayModes &decay_modes_anti = decaymodes[find_offset(pdg_anti)];
        for (const auto &mode : decay_modes_orig.decay_mode_list()) {
          ParticleTypePtrList list = mode->particle_types();
          for (auto &type : list) {
            if (type->has_antiparticle()) {
              type = type->get_antiparticle();
            }
          }
          decay_modes_anti.add_mode(mode->weight(), mode->angular_momentum(),
                                    list);
        }
      }
    }
  };

  for (const Line &line : line_parser(input)) {
    const auto trimmed = trim(line.text);
    assert(!trimmed.empty());  // trim(line.text) is never empty,
                               // else line_parser is broken
    if (trimmed.find_first_of(" \t") == std::string::npos) {
      // a single record on one line signifies a new decay mode section
      end_of_decaymodes();
      std::string name = trim(line.text);
      isotype_mother = &IsoParticleType::find(name);
      mother_states = isotype_mother->get_states();
      decay_modes_to_add.clear();
      decay_modes_to_add.resize(mother_states.size());
      log.debug("reading decay modes for " + name);
    } else {
      std::istringstream lineinput(line.text);
      std::vector<std::string> decay_particles;
      decay_particles.reserve(3);
      float ratio;
      lineinput >> ratio;

      int L;
      lineinput >> L;
      if (L < 0 || L > 3) {  // at some point we might need to support L up to 4
                             // (cf. BlattWeisskopf in decaytype.cc)
        throw LoadFailure("Invalid angular momentum '" + std::to_string(L) +
                          "' in decaymodes.txt:" + std::to_string(line.number) +
                          ": '" + line.text + "'");
      }

      std::string name;
      lineinput >> name;
      bool multi = true;  // does the decay channel refer to whole multiplets?
      while (lineinput) {
        decay_particles.emplace_back(name);
        multi &= IsoParticleType::exists(name);
        lineinput >> name;
      }
      if (multi) {
        // set up multiplet types
        assert(decay_particles.size() == 2);
        const IsoParticleType &isotype_daughter_1
                              = IsoParticleType::find(decay_particles[0]);
        const IsoParticleType &isotype_daughter_2
                              = IsoParticleType::find(decay_particles[1]);
        // loop through multiplets
        for (unsigned int m = 0; m < mother_states.size(); m++) {
          for (const auto daughter1 : isotype_daughter_1.get_states()) {
            for (const auto daughter2 : isotype_daughter_2.get_states()) {
              // calculate Clebsch-Gordan factor
              const double cg = isospin_clebsch_gordan(*daughter1, *daughter2,
                                                       *mother_states[m]);
              const double cg_sqr = cg*cg;
              if (cg_sqr > 0.) {
                // add mode
                log.debug("decay mode generated: " + mother_states[m]->name() +
                          " -> " + daughter1->name() + " " + daughter2->name() +
                          " (" + std::to_string(ratio*cg_sqr) + ")");
                decay_modes_to_add[m].add_mode(ratio*cg_sqr, L,
                                            {daughter1, daughter2});
              }
            }
          }
        }
      } else {
        assert(mother_states.size() == 1);
        ParticleTypePtrList types;
        for (unsigned int i = 0; i < decay_particles.size(); i++) {
          types.push_back(IsoParticleType::find_state(decay_particles[i]));
        }
        log.debug("decay mode found: " + isotype_mother->name() + " -> " +
                  std::to_string(decay_particles.size()));
        decay_modes_to_add[0].add_mode(ratio, L, types);
      }
    }
  }
  end_of_decaymodes();
}

}  // namespace Smash
