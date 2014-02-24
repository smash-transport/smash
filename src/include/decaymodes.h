/*
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_DECAYMODES_H_
#define SRC_INCLUDE_DECAYMODES_H_

#include <cstdio>
#include <vector>

#include "include/constants.h"
#include "include/processbranch.h"
#include <stdexcept>

class DecayModes {
 public:
  /* Add a decay mode */
  inline void add_mode(std::vector<int> particles, float ratio);
  inline void add_mode(ProcessBranch mode);
  /* Make sure ratios add to 1 */
  inline void renormalize(float renormalization_constant);
  /* Remove all modes */
  inline void clear(void);
  /* Check if empty */
  inline bool empty(void) const;
  /* Pass out the decay modes */
  inline std::vector<ProcessBranch> decay_mode_list(void) const;

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

/* Add a decay mode */
inline void DecayModes::add_mode(std::vector<int> particles, float ratio) {
  if (particles.size() < 2) {
    throw InvalidDecay(
        "DecayModes::add_mode was instructed to add a decay mode with less "
        "than 2 branches. This is an invalid input.");
  }
  ProcessBranch branch;
  branch.set_particles(std::move(particles));
  branch.set_weight(ratio);
  decay_modes_.push_back(branch);
}

inline void DecayModes::add_mode(ProcessBranch branch) {
  decay_modes_.push_back(branch);
}

/* Make sure ratios add to 1 */
inline void DecayModes::renormalize(float renormalization_constant) {
  if (renormalization_constant < really_small) {
    printf("Warning: Extremely small renormalization constant: %g\n",
           renormalization_constant);
    printf("Skipping the renormalization.\n");
  } else {
    printf("Renormalizing decay modes with %g \n", renormalization_constant);
    float new_sum = 0.0;
    for (std::vector<ProcessBranch>::iterator mode
           = decay_modes_.begin(); mode != decay_modes_.end(); ++mode) {
      mode->set_weight(mode->weight() / renormalization_constant);
      new_sum += mode->weight();
    }
    printf("After renormalization sum of ratios is %g. \n", new_sum);
  }
}

/* Remove all modes */
inline void DecayModes::clear(void) {
  decay_modes_.clear();
}

/* Check if empty */
inline bool DecayModes::empty(void) const {
  return decay_modes_.empty();
}

/* Pass out the decay modes */
inline std::vector<ProcessBranch> DecayModes::decay_mode_list(void) const {
  return decay_modes_;
}

#endif  // SRC_INCLUDE_DECAYMODES_H_
