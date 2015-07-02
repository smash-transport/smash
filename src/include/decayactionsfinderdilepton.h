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
#include <string>
#include <memory>

#include "decayactionsfinder.h"
#include "oscaroutput.h"
#include "forwarddeclarations.h"
#include "cxx14compat.h"

namespace Smash {

/**
 * \ingroup action
 * A simple dilepton decay finder:
 * Just loops through all particles and checks if they can decay into dileptons.
 */
class DecayActionsFinderDilepton : public DecayActionsFinder {
 public:

  /** Initialize the finder */
  DecayActionsFinderDilepton(bf::path output_path, Configuration config, std::string name)
     {
       dil_out_=create_select_format<OscarInteractions>(std::move(output_path), std::move(config), "full_event_history")
       //dil_out_=make_unique<OscarOutput<OscarFormat2013, OscarInteractions>>(std::move(output_path), std::move(name));
     }
  /** Check the whole particle list for decays
   * and return a list with the corrsponding Action objects. */
  ActionList find_possible_actions(
      const ParticleList &search_list,
      float dt) const override;

  /** Force all resonances to decay at the end of the simulation. */
  ActionList find_final_actions(const Particles &search_list) const override;

 private:
  build_unique_ptr_<OutputInterface> dil_out_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONSFINDERDILEPTON_H_
