/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_DECAYMODES_H_
#define SRC_INCLUDE_DECAYMODES_H_

#include <cstdio>
#include <utility>
#include <vector>

class DecayModes {
 public:
  /* Add a decay mode */
  inline void add_mode(std::vector<int> particles, float ratio);
  /* Make sure ratios add to 1 */
  inline void renormalize(const float renormalization_constant);
  /* Remove all modes */
  inline void clear(void);
  /* Pass out the decay modes */
  inline std::vector< std::pair<std::vector<int>, float> >
    decay_mode_list(void) const;
 private:
  /* Vector of decay modes.
   * Each mode consists of a vector of the pdg codes of decay products
   * and a ratio of this decay mode compared to all possible modes
   */
  std::vector< std::pair<std::vector<int>, float> > decay_modes_;
};

/* Add a decay mode */
inline void DecayModes::add_mode(std::vector<int> particles, float ratio) {
  std::pair<std::vector<int>, float> decay_mode
    = std::make_pair(particles, ratio);
  decay_modes_.push_back(decay_mode);
}

/* Make sure ratios add to 1 */
inline void DecayModes::renormalize(const float renormalization_constant) {
  printf("Renormalizing decay modes with %g \n", renormalization_constant);
  float new_sum = 0.0;
  for (std::vector< std::pair<std::vector<int>, float> >::iterator mode
         = decay_modes_.begin(); mode != decay_modes_.end(); ++mode) {
    mode->second = mode->second / renormalization_constant;
    new_sum += mode->second;
  }
  printf("After renormalization sum of ratios is %g. \n", new_sum);
}

/* Remove all modes */
inline void DecayModes::clear(void) {
  decay_modes_.clear();
}

/* Pass out the decay modes */
inline std::vector< std::pair<std::vector<int>, float> >
  DecayModes::decay_mode_list(void) const {
  return decay_modes_;
}

#endif  // SRC_INCLUDE_DECAYMODES_H_
