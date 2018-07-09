/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/decaymodes.h"

#include <vector>

#include "smash/clebschgordan.h"
#include "smash/constants.h"
#include "smash/cxx14compat.h"
#include "smash/inputfunctions.h"
#include "smash/isoparticletype.h"
#include "smash/logging.h"
#include "smash/stringfunctions.h"

namespace smash {

/// Global pointer to the decay modes list
std::vector<DecayModes> *DecayModes::all_decay_modes = nullptr;
/// Global pointer to the decay types list
std::vector<DecayTypePtr> *all_decay_types = nullptr;

void DecayModes::add_mode(ParticleTypePtr mother, double ratio, int L,
                          ParticleTypePtrList particle_types) {
  DecayType *type = get_decay_type(mother, particle_types, L);
  // Check if mode already exists: if yes, add weight.
  for (auto &mode : decay_modes_) {
    if (type == &mode->type()) {
      mode->set_weight(mode->weight() + ratio);
      return;
    }
  }
  // Add new mode.
  decay_modes_.push_back(make_unique<DecayBranch>(*type, ratio));
}

DecayType *DecayModes::get_decay_type(ParticleTypePtr mother,
                                      ParticleTypePtrList particle_types,
                                      int L) {
  assert(all_decay_types != nullptr);

  // check if the decay type already exisits
  for (const auto &type : *all_decay_types) {
    if (type->has_mother(mother) && type->has_particles(particle_types) &&
        type->angular_momentum() == L) {
      return type.get();
    }
  }

  // if the type does not exist yet, create a new one
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
            make_unique<ThreeBodyDecayDilepton>(mother, particle_types, L));
      } else {
        all_decay_types->emplace_back(
            make_unique<ThreeBodyDecay>(particle_types, L));
      }
      break;
    default:
      throw InvalidDecay(
          "DecayModes::get_decay_type was instructed to add a decay mode "
          "with " +
          std::to_string(particle_types.size()) +
          " particles. This is an invalid input.");
  }

  return all_decay_types->back().get();
}

void DecayModes::renormalize(std::string name) {
  const auto &log = logger<LogArea::DecayModes>();
  double sum = 0.;
  for (auto &mode : decay_modes_) {
    sum += mode->weight();
  }
  if (std::abs(sum - 1.) < really_small) {
    log.debug("Particle ", name,
              ": Extremely small renormalization constant: ", sum,
              "\n=> Skipping the renormalization.");
  } else {
    if (std::abs(sum - 1.) < 0.01) {
      // Reasonably small correction: do not warn, only give a debug message
      log.debug("Particle ", name, ": Renormalizing decay modes with ", sum);
    } else {
      log.warn("Particle ", name, ": Renormalizing decay modes with ", sum);
    }
    double new_sum = 0.0;
    for (auto &mode : decay_modes_) {
      mode->set_weight(mode->weight() / sum);
      new_sum += mode->weight();
    }
    log.debug("After renormalization sum of ratios is ", new_sum);
  }
}

/* unnamed namespace; its content is accessible only in this file.
 * This is identically the same as have a static function, but if one
 * wanted to specify a file-local type, this would be the way to go. */
namespace {
/**
 * Passes back the address offset of a particletype in the list of all particles
 *
 * \param[in] pdg the pdg code of the sought particle type
 * \return the address offset of this particle type
 */
inline std::size_t find_offset(PdgCode pdg) {
  return std::addressof(ParticleType::find(pdg)) -
         std::addressof(ParticleType::list_all()[0]);
}
}  // unnamed namespace

void DecayModes::load_decaymodes(const std::string &input) {
  const auto &log = logger<LogArea::DecayModes>();
  // create the DecayType vector first, then it outlives the DecayModes vector,
  // which references the DecayType objects.
  static std::vector<DecayTypePtr> decaytypes;
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
        const ParticleType &type_anti = ParticleType::find(pdg_anti);
        DecayModes &decay_modes_orig = decaymodes[find_offset(pdg)];
        DecayModes &decay_modes_anti = decaymodes[find_offset(pdg_anti)];
        for (const auto &mode : decay_modes_orig.decay_mode_list()) {
          ParticleTypePtrList list = mode->particle_types();
          for (auto &type : list) {
            if (type->has_antiparticle()) {
              type = type->get_antiparticle();
            }
          }
          decay_modes_anti.add_mode(&type_anti, mode->weight(),
                                    mode->angular_momentum(), list);
        }
      }
    }
  };

  // Track the line number for better error messages.
  // FIXME: At the moment this does not include comments and empty lines.
  uint64_t linenumber = 1;
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
      // check if any of the states have decay modes already
      for (size_t m = 0; m < mother_states.size(); m++) {
        PdgCode pdgcode = mother_states[m]->pdgcode();
        if (!decaymodes[find_offset(pdgcode)].is_empty()) {
          throw LoadFailure("Duplicate entry for " + name +
                            " in decaymodes.txt:" + std::to_string(linenumber));
        }
      }
    } else {
      std::istringstream lineinput(line.text);
      std::vector<std::string> decay_particles;
      decay_particles.reserve(3);
      double ratio;
      lineinput >> ratio;

      int L;
      lineinput >> L;
      if (L < 0 || L > 4) {  // at some point we might need to support L > 4 ?
                             // (cf. blatt_weisskopf_sqr in formfactors.h)
        throw LoadFailure("Invalid angular momentum '" + std::to_string(L) +
                          "' in decaymodes.txt:" + std::to_string(line.number) +
                          ": '" + line.text + "'");
      }

      std::string name;
      lineinput >> name;
      bool multi = true;
      while (lineinput) {
        decay_particles.emplace_back(name);
        const auto isotype = IsoParticleType::try_find(name);
        const bool is_multiplet = isotype;
        const bool is_state = ParticleType::exists(name);
        if (!is_multiplet && !is_state) {
          throw InvalidDecay(
              "Daughter " + name +
              " is neither an isospin multiplet nor a particle." + " (line " +
              std::to_string(linenumber) + ": \"" + trimmed + "\")");
        }
        const bool is_hadronic_multiplet =
            is_multiplet && isotype->get_states()[0]->is_hadron();
        multi &= is_hadronic_multiplet;
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
            bool forbidden_by_isospin = true;
            for (size_t m = 0; m < mother_states.size(); m++) {
              for (const auto &daughter1 : isotype_daughter_1.get_states()) {
                for (const auto &daughter2 : isotype_daughter_2.get_states()) {
                  // calculate Clebsch-Gordan factor
                  const double cg_sqr = isospin_clebsch_gordan_sqr_2to1(
                      *daughter1, *daughter2, *mother_states[m]);
                  if (cg_sqr > 0.) {
                    // add mode
                    log.debug(
                        "decay mode generated: " + mother_states[m]->name() +
                        " -> " + daughter1->name() + " " + daughter2->name() +
                        " (" + std::to_string(ratio * cg_sqr) + ")");
                    decay_modes_to_add[m].add_mode(mother_states[m],
                                                   ratio * cg_sqr, L,
                                                   {daughter1, daughter2});
                    forbidden_by_isospin = false;
                  }
                }
              }
            }
            if (forbidden_by_isospin) {
              std::stringstream s;
              s << ",\nwhere isospin mother: " << isotype_mother->isospin()
                << ", daughters: " << isotype_daughter_1.isospin() << " "
                << isotype_daughter_2.isospin();
              throw InvalidDecay(isotype_mother->name() +
                                 " decay mode is forbidden by isospin: \"" +
                                 line.text + "\"" + s.str());
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
                    const double cg_sqr = isospin_clebsch_gordan_sqr_3to1(
                        *daughter1, *daughter2, *daughter3, *mother_states[m]);
                    if (cg_sqr > 0.) {
                      // add mode
                      log.debug(
                          "decay mode generated: " + mother_states[m]->name() +
                          " -> " + daughter1->name() + " " + daughter2->name() +
                          " " + daughter3->name() + " (" +
                          std::to_string(ratio * cg_sqr) + ")");
                      decay_modes_to_add[m].add_mode(
                          mother_states[m], ratio * cg_sqr, L,
                          {daughter1, daughter2, daughter3});
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
                line.text + " (line " + std::to_string(linenumber) + ": \"" +
                trimmed + "\")");
        }
      } else {
        /* References to specific states, not multiplets:
         * Loop over all mother states and check charge conservation. */
        ParticleTypePtrList types;
        int charge = 0;
        for (auto part : decay_particles) {
          try {
            types.push_back(IsoParticleType::find_state(part));
          } catch (std::runtime_error &e) {
            throw std::runtime_error(std::string() + e.what() + " (line " +
                                     std::to_string(linenumber) + ": \"" +
                                     trimmed + "\")");
          }
          charge += types.back()->charge();
        }
        bool no_decays = true;
        for (size_t m = 0; m < mother_states.size(); m++) {
          if (mother_states[m]->charge() == charge) {
            log.debug("decay mode found: " + mother_states[m]->name() + " -> " +
                      std::to_string(decay_particles.size()));
            decay_modes_to_add[m].add_mode(mother_states[m], ratio, L, types);
            no_decays = false;
          }
        }
        if (no_decays) {
          throw InvalidDecay(isotype_mother->name() +
                             " decay mode violates charge conservation: \"" +
                             line.text + "\" (line " +
                             std::to_string(linenumber) + ": \"" + trimmed +
                             "\"");
        }
      }
    }
    linenumber++;
  }
  end_of_decaymodes();

  // Check whether the mother's pole mass is strictly larger than the minimal
  // masses of the daughters. This is required by the Manley-Saleski ansatz.
  const auto &particles = ParticleType::list_all();
  for (const auto &mother : particles) {
    if (mother.is_stable()) {
      continue;
    }
    const auto &decays = mother.decay_modes().decay_mode_list();
    for (const auto &decay : decays) {
      if (mother.mass() <= decay->threshold()) {
        std::stringstream s;
        s << mother.name() << " →  ";
        for (const auto p : decay->particle_types()) {
          s << p->name();
        }
        s << " with " << mother.mass() << " ≤ " << decay->threshold();
        throw InvalidDecay(
            "For all decays, the minimum mass of daughters"
            "must be smaller\nthan the mother's pole mass "
            "(Manley-Saleski Ansatz)\n"
            "Violated by the following decay: " +
            s.str());
      }
    }
  }
}

}  // namespace smash
