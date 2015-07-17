/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_
#define SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_

#include "decayactionsfinder.h"

namespace Smash {

/**
 * \ingroup action
 * A dilepton decay finder:
 * Loops through all particles and if they can decay into dileptons, it
 * treats the decays with the shining method. This means in every timestep
 * it stores every possible decay in an extra action and later in the
 * write_dilepton_action routine in experiment.cc these decays are written
 * in the DileptonOutput. Because too many dileptons are produced this way,
 * every decay is weighted with the so called "shining_weight".
 * see \iref{Schmidt:2008hm}, chapter 2D
 * Currently the finder only works with two-body-dilepton decays.
 * Dalitz dilepton decays are not supported.
 */

class DecayActionsFinderDilepton : public DecayActionsFinder {
 public:
  /** Initialize the finder */
  DecayActionsFinderDilepton() {}

  /** Check the whole particle list for possible dilepton decays
   * and return a list with the corrsponding Action objects. */
  ActionList find_actions_in_cell(
      const ParticleList &search_list,
      float dt) const override;

  /** All resonances decay at the end of the simulation. Take all the shining
   * weight for the remaining time.
   */
  ActionList find_final_actions(const Particles &search_list) const override;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_
