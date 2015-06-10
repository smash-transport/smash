/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */


#include "include/dileptondecayactionsfinder.h"


#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/decayaction.h"
#include "include/experimentparameters.h"
#include "include/fourvector.h"
#include "include/particles.h"
#include "include/random.h"

#include <typeinfo>
#include <iostream>
#include <fstream>

namespace Smash {

ActionList DileptonDecayActionsFinder::find_possible_actions(
      const ParticleList &search_list,
      const std::vector<const ParticleList *> &,  // neighbors irrelevant
      float dt) const {
  ActionList actions;

  std::FILE* dileptons;
  dileptons = std::fopen("dilepton_ouput.txt", "a");


  for (const auto &p : search_list) {
    if (p.type().is_stable()) {
      continue;
    }

    // work in progress //

    /* current (first) implementation: everything happens in finder
     * It slows down smash probably because of the extra output
     * file.
     * The plan is to use the standard output routine in next versions.
     */

    DecayBranchList all_modes = p.type().get_partial_widths(p.effective_mass());

    /* dileptons should be radiated at the end of the timestep,
        so they wont radiate if the resonance decays (so dt)*/
    auto act = make_unique<DecayAction>(p, dt);

    for (DecayBranchPtr & mode : all_modes) {
        if (is_dilepton(mode->type().particle_types()[0]->pdgcode(),
                        mode->type().particle_types()[1]->pdgcode())) {
          act->add_decay(std::move(mode));
          /* problem for more than one dilepton mode per particle
           * later generate_final_state would choose one
           * not conform with shinning method
           * doesn't matter for know only e+e- decay channel switched
           * on
           */
      }
    }



    if (act->total_width() > 0.0) {  // check if their are any (dil.) decays
      actions.emplace_back(std::move(act));
      actions[0]->generate_final_state();  // should only be one action
      // for a first version the finder writes its own dilepton output
      // current outputformat pdg p0 p1 p2 p3
      std::fprintf(dileptons, "# --------------\n");
      std::fprintf(dileptons, "%s %g %g %g %g\n",
               actions[0]->outgoing_particles()[0].pdgcode().string().c_str(),
               actions[0]->outgoing_particles()[0].momentum().x0(),
               actions[0]->outgoing_particles()[0].momentum().x1(),
               actions[0]->outgoing_particles()[0].momentum().x2(),
               actions[0]->outgoing_particles()[0].momentum().x3());
      std::fprintf(dileptons, "%s %g %g %g %g\n",
               actions[0]->outgoing_particles()[1].pdgcode().string().c_str(),
               actions[0]->outgoing_particles()[1].momentum().x0(),
               actions[0]->outgoing_particles()[1].momentum().x1(),
               actions[0]->outgoing_particles()[1].momentum().x2(),
               actions[0]->outgoing_particles()[1].momentum().x3());
    }
  }

  fclose(dileptons);

  /** in the current impl. nothing else should be done by the experiment,
      so I return an empty action list*/
  ActionList empty_actionlist;

  return std::move(empty_actionlist);
}

ActionList DileptonDecayActionsFinder::find_final_actions(
                  const Particles &) const {   // temp. rmvd search_list (warn.)
  // not done yet
  ActionList empty_actionlist;

  return std::move(empty_actionlist);
}

}  // namespace Smash
