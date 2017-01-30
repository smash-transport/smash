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

void DecayActionsFinderDilepton::shine(
      const Particles &search_list,
      OutputInterface* output,
      float dt) const {

  for (const auto &p : search_list) {
    // effective mass of decaying particle
    const float m_eff = p.effective_mass();
    const auto n_all_modes = p.type().get_partial_widths(m_eff).size();
    if (n_all_modes == 0) {
      continue;
    }

    const float inv_gamma = p.inverse_gamma();
    DecayBranchList dil_modes = p.type().get_partial_widths_dilepton(m_eff);

    // if particle can only decay into dileptons or is stable, use shining only
    // in find_final_actions and ignore them here, also unformed
    // resonances cannot decay
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wfloat-equal"
    if (dil_modes.size() == n_all_modes || p.type().is_stable()
        || (p.formation_time() > p.position().x0() &&
            p.cross_section_scaling_factor() == 0.0)) {
      continue;
    }
    #pragma GCC diagnostic pop

    for (DecayBranchPtr & mode : dil_modes) {
      // SHINING as described in \iref{Schmidt:2008hm}, chapter 2D
      const float shining_weight = dt * inv_gamma * mode->weight() / hbarc;

      if (shining_weight > 0.0) {  // decays that can happen
        DecayActionDilepton act(p, 0.f, shining_weight);
        act.add_decay(std::move(mode));
        act.generate_final_state();
        output->at_interaction(act, 0.0);
      }
    }
  }
}


void DecayActionsFinderDilepton::shine_final(
                  const Particles &search_list,
                  OutputInterface* output,
                  bool only_res) const {

  for (const auto &p : search_list) {
    const ParticleType &t = p.type();
    if (t.decay_modes().decay_mode_list().empty() ||
        (only_res && t.is_stable())) {
      continue;
    }

    // effective mass of decaying particle
    const float m_eff = p.effective_mass();
    DecayBranchList dil_modes = t.get_partial_widths_dilepton(m_eff);

    // total decay width, also hadronic decays
    const float width_tot = total_weight<DecayBranch>(
                                                  t.get_partial_widths(m_eff));

    for (DecayBranchPtr & mode : dil_modes) {
      const float shining_weight = mode->weight() / width_tot;

      if (shining_weight > 0.0) {  // decays that can happen
        DecayActionDilepton act(p, 0.f, shining_weight);
        act.add_decay(std::move(mode));
        act.generate_final_state();
        output->at_interaction(act, 0.0);
      }
    }
  }
}

}  // namespace Smash
