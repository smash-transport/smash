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

#include <boost/filesystem.hpp>
#include <vector>

#include "decayactionsfinder.h"
#include "dileptonoutput.h"
#include "forwarddeclarations.h"

namespace Smash {

/**
 * \ingroup action
 * A simple dilepton decay finder:
 * Just loops through all particles and checks if they can decay into dileptons.
 */
class DecayActionsFinderDilepton : public DecayActionsFinder {
 public:
  /** Initialize the finder */
  DecayActionsFinderDilepton(bf::path output_path)
     :dil_out_{new DileptonOutput(output_path)} {}
  /** Check the whole particle list for decays
   * and return a list with the corrsponding Action objects. */
  ActionList find_possible_actions(
      const ParticleList &search_list,
      float dt) const override;

  /** Force all resonances to decay at the end of the simulation. */
  ActionList find_final_actions(const Particles &search_list) const override;

 private:
  DileptonOutput* dil_out_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_
