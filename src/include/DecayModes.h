/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_DECAYMODES_H_
#define SRC_INCLUDE_DECAYMODES_H_

#include <utility>
#include <vector>

class DecayModes {
 public:
  /* Add a decay mode */
  inline void add_mode(std::vector<int> particles, float ratio);
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
