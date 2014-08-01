/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DECAYMODES_H_
#define SRC_INCLUDE_DECAYMODES_H_

#include "processbranch.h"

#include <stdexcept>
#include <vector>

namespace Smash {

class DecayModes {
 public:
  /* Add a decay mode */
  void add_mode(float ratio, int L, std::vector<PdgCode> pdg_list);
  void add_mode(DecayBranch branch) { decay_modes_.push_back(branch); }

  /* Make sure ratios add to 1 */
  void renormalize(float renormalization_constant);

  /* Remove all modes */
  void clear() { decay_modes_.clear(); }

  /* Check if empty */
  bool is_empty() const { return decay_modes_.empty(); }

  /* Pass out the decay modes */
  const std::vector<DecayBranch> &decay_mode_list(void) const {
    return decay_modes_;
  }

  /// Return decay modes of this particle type
  static const DecayModes &find(PdgCode pdg);

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
  std::vector<DecayBranch> decay_modes_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYMODES_H_
