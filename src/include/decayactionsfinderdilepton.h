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

#include "outputinterface.h"

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
 * See \iref{Schmidt:2008hm}, chapter 2D.
 * The finder works with two body dilepton decays as well as with dalitz
 * dilepton decays.
 */

class DecayActionsFinderDilepton {
 public:
  /** Initialize the finder */
  DecayActionsFinderDilepton() {}

  /// Check the whole particles list and print out possible dilepton decays
  void shine(const Particles &search_list,
             OutputInterface* output,
             double dt) const;

  /** Shine dileptons from resonances at the end of the simulation. This is
   * special, because the shining time is now until resonance decays and not
   * some fixed dt interval.
   */
  void shine_final(const Particles &search_list,
                   OutputInterface* output,
                   bool only_res = false) const;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_
