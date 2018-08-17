/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_
#define SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_

#include "outputinterface.h"

namespace smash {

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
  /// Initialize the finder
  DecayActionsFinderDilepton() {}

  /**
   * Check the whole particles list and print out possible dilepton decays.
   *
   * \param[in] search_list List of all particles.
   * \param[in] output Pointer to the dilepton output.
   * \param[in] dt Length of timestep [fm]
   */
  void shine(const Particles& search_list, OutputInterface* output,
             double dt) const;

  /**
   * Shine dileptons from resonances at the end of the simulation.
   *
   * This is special, because the shining time is now until the resonance would
   * decay and not for some fixed dt interval.
   *
   * \param[in] search_list List of all particles.
   * \param[in] output Pointer to the dilepton output.
   * \param[in] only_res optional parameter that requests that only actions
   *                     regarding resonances are considered (disregarding
   *                     stable particles)
   */
  void shine_final(const Particles& search_list, OutputInterface* output,
                   bool only_res = false) const;
};

}  // namespace smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_
