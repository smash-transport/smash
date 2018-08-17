/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/decayactionsfinderdilepton.h"

#include "smash/constants.h"
#include "smash/cxx14compat.h"
#include "smash/decayactiondilepton.h"
#include "smash/decaymodes.h"
#include "smash/particles.h"

namespace smash {

void DecayActionsFinderDilepton::shine(const Particles &search_list,
                                       OutputInterface *output,
                                       double dt) const {
  if (!output->is_dilepton_output()) {
    return;
  }
  for (const auto &p : search_list) {
    // effective mass of decaying particle
    const double m_eff = p.effective_mass();
    const auto n_all_modes = p.type().get_partial_widths(m_eff).size();
    if (n_all_modes == 0) {
      continue;
    }

    const double inv_gamma = p.inverse_gamma();
    DecayBranchList dil_modes = p.type().get_partial_widths_dilepton(m_eff);

    /* if particle can only decay into dileptons or is stable, use shining only
     * in find_final_actions and ignore them here, also unformed
     * resonances cannot decay */
    if (dil_modes.size() == n_all_modes || p.type().is_stable() ||
        (p.formation_time() > p.position().x0())) {
      continue;
    }

    for (DecayBranchPtr &mode : dil_modes) {
      // SHINING as described in \iref{Schmidt:2008hm}, chapter 2D
      const double shining_weight = dt * inv_gamma * mode->weight() / hbarc;

      if (shining_weight > 0.0) {  // decays that can happen
        DecayActionDilepton act(p, 0., shining_weight);
        act.add_decay(std::move(mode));
        act.generate_final_state();
        output->at_interaction(act, 0.0);
      }
    }
  }
}

void DecayActionsFinderDilepton::shine_final(const Particles &search_list,
                                             OutputInterface *output,
                                             bool only_res) const {
  if (!output->is_dilepton_output()) {
    return;
  }
  for (const auto &p : search_list) {
    const ParticleType &t = p.type();
    if (t.decay_modes().decay_mode_list().empty() ||
        (only_res && t.is_stable())) {
      continue;
    }

    // effective mass of decaying particle
    const double m_eff = p.effective_mass();
    DecayBranchList dil_modes = t.get_partial_widths_dilepton(m_eff);

    // total decay width, also hadronic decays
    const double width_tot =
        total_weight<DecayBranch>(t.get_partial_widths(m_eff));

    for (DecayBranchPtr &mode : dil_modes) {
      const double shining_weight = mode->weight() / width_tot;

      if (shining_weight > 0.0) {  // decays that can happen
        DecayActionDilepton act(p, 0., shining_weight);
        act.add_decay(std::move(mode));
        act.generate_final_state();
        output->at_interaction(act, 0.0);
      }
    }
  }
}

}  // namespace smash
