/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/decayactionsfinderdilepton.h"

#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/decayactiondilepton.h"
#include "include/decaymodes.h"
#include "include/particles.h"


namespace Smash {

ActionList DecayActionsFinderDilepton::find_actions_in_cell(
      const ParticleList &search_list,
      float dt) const {
  ActionList actions;

  for (const auto &p : search_list) {
    unsigned long n_all_modes = p.type().decay_modes().decay_mode_list().size();
    if (n_all_modes == 0) {
      continue;
    }

    const float inv_gamma = p.inverse_gamma();
    DecayBranchList dil_modes =
                  p.type().get_partial_widths_dilepton(p.effective_mass());

    // if particle can only decay into dileptons, use shining in only in
    // find_final_actions and ignore them here
    if (dil_modes.size() == n_all_modes) {
      continue;
    }

    for (DecayBranchPtr & mode : dil_modes) {
      const float partial_width = mode->weight();
      // SHINING as described in \iref{Schmidt:2008hm}, chapter 2D
      const float sh_weight = dt * partial_width * inv_gamma / hbarc;

      auto act = make_unique<DecayActionDilepton>(p, 0.f, sh_weight);
      act->add_decay(std::move(mode));
      actions.emplace_back(std::move(act));
    }
  }

  return std::move(actions);
}


ActionList DecayActionsFinderDilepton::find_final_actions(
                  const Particles &search_list) const {
  ActionList actions;

  for (const auto &p : search_list) {
    if (p.type().is_stable()) {
      continue;
    }

    DecayBranchList dil_modes =
                  p.type().get_partial_widths_dilepton(p.effective_mass());

    // total decay width, also hadronic decays
    const float width_tot = total_weight<DecayBranch>(
                               p.type().get_partial_widths(p.effective_mass()));

    for (DecayBranchPtr & mode : dil_modes) {
      const float partial_width = mode->weight();
      const float sh_weight = partial_width / width_tot;

      auto act = make_unique<DecayActionDilepton>(p, 0.f, sh_weight);
      act->add_decay(std::move(mode));
      actions.emplace_back(std::move(act));
    }
  }

  return std::move(actions);
}

}  // namespace Smash
