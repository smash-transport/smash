/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DECAYMODES_H_
#define SRC_INCLUDE_DECAYMODES_H_

#include <stdexcept>
#include <string>
#include <vector>

#include "processbranch.h"

namespace Smash {

/**
 * \ingroup data
 *
 * If you want to find a DecayModes object for a specific particle type use
 * ParticleType::decay_modes().
 */
class DecayModes {
 public:
  /* Add a decay mode */
  void add_mode(float ratio, int L, ParticleTypePtrList particle_types);
  void add_mode(DecayBranchPtr branch) {
    decay_modes_.push_back(std::move(branch));
  }

  /* Make sure branching ratios add up to 1. */
  void renormalize(std::string name);

  /* Check if empty */
  bool is_empty() const { return decay_modes_.empty(); }

  /* Pass out the decay modes */
  const DecayBranchList &decay_mode_list() const {
    return decay_modes_;
  }

  /**
   * Loads the DecayModes map as described in the \p input string.
   *
   * It does sanity checking - that the particles it talks about are in the
   * ParticleType map.
   */
  static void load_decaymodes(const std::string &input);

  /// \ingroup exception
  struct InvalidDecay : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };
  /// \ingroup exception
  struct LoadFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };
  /// \ingroup exception
  struct MissingDecays : public LoadFailure {
    using LoadFailure::LoadFailure;
  };
  /// \ingroup exception
  struct ReferencedParticleNotFound : public LoadFailure {
    using LoadFailure::LoadFailure;
  };

 private:
  /* Vector of decay modes.
   * Each mode consists of a vector of the pdg codes of decay products
   * and a ratio of this decay mode compared to all possible modes
   */
  DecayBranchList decay_modes_;

  /// allow ParticleType::decay_modes to access all_decay_modes
  friend const DecayModes &ParticleType::decay_modes() const;

  /**
   * A list of all DecayModes objects using the same indexing as
   * all_particle_types.
   */
  static std::vector<DecayModes> *all_decay_modes;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYMODES_H_
