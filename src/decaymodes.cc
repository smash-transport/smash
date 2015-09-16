/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/decaymodes.h"

#include <vector>

#include "include/clebschgordan.h"
#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/inputfunctions.h"
#include "include/isoparticletype.h"
#include "include/logging.h"
#include "include/stringfunctions.h"

namespace Smash {

std::vector<DecayModes> *DecayModes::all_decay_modes = nullptr;

std::vector<std::unique_ptr<DecayType>> *all_decay_types = nullptr;

void DecayModes::add_mode(float ratio, int L,
                          ParticleTypePtrList particle_types) {
  assert(all_decay_types != nullptr);
  switch (particle_types.size()) {
    case 2:
      if (is_dilepton(particle_types[0]->pdgcode(),
                      particle_types[1]->pdgcode())) {
        all_decay_types->emplace_back(
            make_unique<TwoBodyDecayDilepton>(particle_types, L));
      } else if (particle_types[0]->is_stable() &&
                 particle_types[1]->is_stable()) {
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
      if (has_lepton_pair(particle_types[0]->pdgcode(),
                          particle_types[1]->pdgcode(),
                          particle_types[2]->pdgcode())) {
        all_decay_types->emplace_back(
        make_unique<ThreeBodyDecayDilepton>(particle_types, L));
      } else {
      all_decay_types->emplace_back(
          make_unique<ThreeBodyDecay>(particle_types, L));
      }
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
    if (std::abs(sum - 1.) < 0.01) {
      // Reasonably small correction: do not warn, only give a debug message
      log.debug("Particle ", name, ": Renormalizing decay modes with ", sum);
    } else {
      log.warn("Particle ", name, ": Renormalizing decay modes with ", sum);
    }
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

  const IsoParticleType *isotype_mother = nullptr;
  ParticleTypePtrList mother_states;
  std::vector<DecayModes> decay_modes_to_add;  // one for each mother state

  const auto end_of_decaymodes = [&]() {
    if (isotype_mother == nullptr) {  // at the start of the file
      return;
    }
    // Loop over all states in the mother multiplet and add modes
    for (size_t m = 0; m < mother_states.size(); m++) {
      if (decay_modes_to_add[m].is_empty() && !mother_states[m]->is_stable()) {
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
      for (const auto &state : mother_states) {
        PdgCode pdg = state->pdgcode();
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
        /* References to isospin multiplets: Automatically determine all valid
         * combinations and calculate Clebsch-Gordan factors */
        switch (decay_particles.size()) {
          case 2: {
            const IsoParticleType &isotype_daughter_1 =
                IsoParticleType::find(decay_particles[0]);
            const IsoParticleType &isotype_daughter_2 =
                IsoParticleType::find(decay_particles[1]);
            // loop through multiplets
            for (size_t m = 0; m < mother_states.size(); m++) {
              for (const auto &daughter1 : isotype_daughter_1.get_states()) {
                for (const auto &daughter2 : isotype_daughter_2.get_states()) {
                  // calculate Clebsch-Gordan factor
                  const double cg = isospin_clebsch_gordan(
                      *daughter1, *daughter2, *mother_states[m]);
                  const double cg_sqr = cg * cg;
                  if (cg_sqr > 0.) {
                    // add mode
                    log.debug("decay mode generated: " +
                              mother_states[m]->name() + " -> " +
                              daughter1->name() + " " + daughter2->name() +
                              " (" + std::to_string(ratio * cg_sqr) + ")");
                    decay_modes_to_add[m].add_mode(ratio * cg_sqr, L,
                                                   {daughter1, daughter2});
                  }
                }
              }
            }
            break;
          }
          case 3: {
            const IsoParticleType &isotype_daughter_1 =
                IsoParticleType::find(decay_particles[0]);
            const IsoParticleType &isotype_daughter_2 =
                IsoParticleType::find(decay_particles[1]);
            const IsoParticleType &isotype_daughter_3 =
                IsoParticleType::find(decay_particles[2]);
            // loop through multiplets
            for (size_t m = 0; m < mother_states.size(); m++) {
              for (const auto &daughter1 : isotype_daughter_1.get_states()) {
                for (const auto &daughter2 : isotype_daughter_2.get_states()) {
                  for (const auto &daughter3 :
                       isotype_daughter_3.get_states()) {
                    // calculate allowed I_12
                    const auto min_I_12 =
                        std::abs(daughter1->isospin() - daughter2->isospin());
                    const auto max_I_12 =
                        daughter1->isospin() + daughter2->isospin();
                    std::vector<int> possible_I_12(max_I_12 - min_I_12 + 1);
                    std::iota(possible_I_12.begin(), possible_I_12.end(),
                              min_I_12);
                    std::vector<int> allowed_I_12;
                    allowed_I_12.reserve(possible_I_12.size());
                    const int target_I = 0;
                    for (const auto I_12 : possible_I_12) {
                      const auto min_I = std::abs(I_12 - daughter3->isospin());
                      const auto max_I = I_12 + daughter3->isospin();
                      if (min_I <= target_I && target_I <= max_I) {
                        allowed_I_12.push_back(I_12);
                      }
                    }
                    if (allowed_I_12.size() != 1) {
                      throw std::runtime_error(
                          "The coupled 3-body isospin state is not uniquely "
                          "defined for " +
                          mother_states[m]->name() + " -> " +
                          daughter1->name() + " " + daughter2->name() + " " +
                          daughter3->name());
                    }
                    const auto I_12 = allowed_I_12[0];
                    const auto &mother = *mother_states[m];
                    const double cg = isospin_clebsch_gordan(
                        *daughter1, *daughter2, *daughter3, mother.isospin(),
                        mother.isospin3(), I_12);
                    const double cg_sqr = cg * cg;
                    if (cg_sqr > 0.) {
                      // add mode
                      log.debug("decay mode generated: " +
                                mother_states[m]->name() + " -> " +
                                daughter1->name() + " " + daughter2->name() +
                                " " + daughter3->name() + " (" +
                                std::to_string(ratio * cg_sqr) + ")");
                      decay_modes_to_add[m].add_mode(
                          ratio * cg_sqr, L, {daughter1, daughter2, daughter3});
                    }
                  }
                }
              }
            }
            break;
          }
          default:
            throw std::runtime_error(
                "References to isospin multiplets only "
                "allowed in two-body or three-body decays: " +
                line.text);
        }
      } else {
        /* References to specific states, not multiplets:
         * Loop over all mother states and check charge conservation. */
        ParticleTypePtrList types;
        int charge = 0;
        for (auto part : decay_particles) {
          types.push_back(IsoParticleType::find_state(part));
          charge += types.back()->charge();
        }
        bool no_decays = true;
        for (size_t m = 0; m < mother_states.size(); m++) {
          if (mother_states[m]->charge() == charge) {
            log.debug("decay mode found: " + mother_states[m]->name() + " -> " +
                      std::to_string(decay_particles.size()));
            decay_modes_to_add[m].add_mode(ratio, L, types);
            no_decays = false;
          }
        }
        if (no_decays) {
          throw InvalidDecay(isotype_mother->name() +
                             " decay mode violates charge conservation: \"" +
                             line.text + "\"");
        }
      }
    }
  }
  end_of_decaymodes();
}

}  // namespace Smash
