/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/decayactionsfinderdilepton.h"

#include <boost/filesystem.hpp>

#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/decayactiondilepton.h"
#include "include/experimentparameters.h"
#include "include/fourvector.h"
#include "include/particles.h"
#include "include/random.h"

#include <typeinfo>
#include <iostream>
#include <fstream>

namespace Smash {

ActionList DecayActionsFinderDilepton::find_possible_actions(
      const ParticleList &search_list,
      float dt) const {

  for (const auto &p : search_list) {
    if (p.type().is_stable()) {
      continue;
    }

    // work in progress //

    /* current (first) implementation: everything happens in finder
     * It slows down smash probably because of the extra output
     * file.
     */

    float inv_gamma = p.inverse_gamma();

    DecayBranchList all_modes =
                  p.type().get_partial_widths_dilepton(p.effective_mass());

    for (DecayBranchPtr & mode : all_modes) {
      float partial_width = mode->weight();
      // SHINNING as described in \iref{Schmidt:2008hm}, chapter 2D
      float sh_weight = dt * partial_width * inv_gamma;
      auto act = make_unique<DecayActionDilepton>(p, 0.f, sh_weight);
      act->add_decay(std::move(mode));

      if (act->total_width() > 0.0) {  // check if their are any (dil.) decays

        act->generate_final_state();

        // output in oscar format
        dil_out_->at_interaction(act->incoming_particles(),
                                 act->outgoing_particles(),
                                 0.0,
                                 act->raw_weight_value(),
                                 act->get_type());
      }
    }
  }


  /** in the current impl. nothing else should be done by the experiment,
      so I return an empty action list*/
  ActionList empty_actionlist;

  return std::move(empty_actionlist);
}

ActionList DecayActionsFinderDilepton::find_final_actions(
                  const Particles &search_list) const {
  for (const auto &p : search_list) {
    if (p.type().is_stable()) {
      continue;      /* particle doesn't decay */
    }
    float inv_gamma = p.inverse_gamma();

    DecayBranchList dil_modes =
                  p.type().get_partial_widths_dilepton(p.effective_mass());

    // total decay width, also hadronic decays
    const float width_tot = total_weight<DecayBranch>(
                               p.type().get_partial_widths(p.effective_mass()));

    for (DecayBranchPtr & mode : dil_modes) {
      float partial_width = mode->weight();

      float sh_weight = partial_width * inv_gamma / width_tot;
      auto act = make_unique<DecayActionDilepton>(p, 0.f, sh_weight);
      act->add_decay(std::move(mode));

      if (act->total_width() > 0.0) {  // check if their are any (dil.) decays

        act->generate_final_state();

        // output in oscar format
        dil_out_->at_interaction(act->incoming_particles(),
                                 act->outgoing_particles(),
                                 1.0,
                                 act->raw_weight_value(),
                                 act->get_type());
      }
    }

  }

  ActionList empty_actionlist;

  return std::move(empty_actionlist);
}

}  // namespace Smash
