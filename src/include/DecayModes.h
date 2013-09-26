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
  /* pass out the decay modes */
  inline std::vector< std::pair<std::vector<int>, float> >
    decay_modes(void) const;
 private:
  /* Vector of decay modes.
   * Each mode consists of a vector of the pdg codes of decay products
   * and a ratio of this decay mode with respect to total
   */
  std::vector< std::pair<std::vector<int>, float> > decay_modes_;
}

inline std::vector< std::pair<std::vector<int>, float> >
  decay_modes(void) const {
  return decay_modes_;
}

#endif  // SRC_INCLUDE_DECAYMODES_H_
