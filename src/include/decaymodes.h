/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DECAYMODES_H_
#define SRC_INCLUDE_DECAYMODES_H_

#include "forwarddeclarations.h"
#include "processbranch.h"

#include <stdexcept>
#include <vector>

namespace Smash {

class DecayModes {
 public:
  /* Add a decay mode */
  void add_mode(std::vector<PdgCode> pdg_list, float ratio);
  void add_mode(ProcessBranch branch) { decay_modes_.push_back(branch); }

  /* Make sure ratios add to 1 */
  void renormalize(float renormalization_constant);

  /* Remove all modes */
  void clear() { decay_modes_.clear(); }

  /* Check if empty */
  bool empty() const { return decay_modes_.empty(); }

  /* Pass out the decay modes */
  const std::vector<ProcessBranch> &decay_mode_list(void) const {
    return decay_modes_;
  }

  struct InvalidDecay : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

 private:
  /* Vector of decay modes.
   * Each mode consists of a vector of the pdg codes of decay products
   * and a ratio of this decay mode compared to all possible modes
   */
  std::vector<ProcessBranch> decay_modes_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYMODES_H_
