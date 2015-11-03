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
    const auto n_all_modes = p.type().get_partial_widths(p.effective_mass()).size();
    if (n_all_modes == 0) {
      continue;
    }

    const float inv_gamma = p.inverse_gamma();
    DecayBranchList dil_modes =
                  p.type().get_partial_widths_dilepton(p.effective_mass());

    // if particle can only decay into dileptons or is stable, use shining only
    // in find_final_actions and ignore them here
    if (dil_modes.size() == n_all_modes || p.type().is_stable()) {
      continue;
    }

    for (DecayBranchPtr & mode : dil_modes) {
      float sh_weight = 0.0;
      float dilepton_mass = 0.f;

      switch (mode->particle_number()) {
        case 2: {
          const float partial_width = mode->weight();
          // SHINING as described in \iref{Schmidt:2008hm}, chapter 2D
          sh_weight = dt * inv_gamma * partial_width / hbarc;
          break;
        }
        case 3: {
          // find the non lepton particle position
          int non_lepton_position = -1;
          for (int i = 0; i < 3; ++i) {
            if (!(mode->particle_types()[i]->is_lepton())) {
              non_lepton_position = i;
              break;
            }
          }

          // mass of non-lepton final state particle
          const float m_nl =
                            mode->particle_types()[non_lepton_position]->mass();
          // mass of leptons in final state
          const float m_l =
                      mode->particle_types()[(non_lepton_position+1)%3]->mass();

          // randomly select a mass
          dilepton_mass = Random::uniform(2*m_l, p.effective_mass()-m_nl);

          const float diff_width = ThreeBodyDecayDilepton::diff_width(
                                              p.effective_mass(), dilepton_mass,
                                              m_nl, p.type().pdgcode());

          sh_weight = dt * inv_gamma * diff_width / hbarc;
          break;
          }
        default:
          throw std::runtime_error("Error in DecayActionFinderDilepton");
      }

      if (sh_weight > 0.0) {  // decays that can happen
        auto act = make_unique<DecayActionDilepton>(p, 0.f, sh_weight,
                                                                 dilepton_mass);
        act->add_decay(std::move(mode));
        actions.emplace_back(std::move(act));
      }
    }
  }

  return std::move(actions);
}


ActionList DecayActionsFinderDilepton::find_final_actions(
                  const Particles &search_list) const {
  ActionList actions;

  for (const auto &p : search_list) {
    if (p.type().decay_modes().decay_mode_list().size() == 0) {
      continue;
    }

    DecayBranchList dil_modes =
                  p.type().get_partial_widths_dilepton(p.effective_mass());

    // total decay width, also hadronic decays
    const float width_tot = total_weight<DecayBranch>(
                               p.type().get_partial_widths(p.effective_mass()));

    for (DecayBranchPtr & mode : dil_modes) {
      float sh_weight;
      float dilepton_mass = 0.f;

      switch (mode->particle_number()) {
        case 2: {
          const float partial_width = mode->weight();
          sh_weight = partial_width / width_tot;
          break;
        }
        case 3: {
          // find the non lepton particle position
          int non_lepton_position = -1;
          for (int i = 0; i < 3; ++i) {
            if (!(mode->particle_types()[i]->is_lepton())) {
              non_lepton_position = i;
              break;
            }
          }
          // mass of non-lepton final state particle
          const float m_nl =
                            mode->particle_types()[non_lepton_position]->mass();
          // mass of leptons in final state
          const float m_l =
                      mode->particle_types()[(non_lepton_position+1)%3]->mass();

          // randomly select a mass
          dilepton_mass = Random::uniform(2*m_l, p.effective_mass()-m_nl);

          const float diff_width = ThreeBodyDecayDilepton::diff_width(
                                              p.effective_mass(), dilepton_mass,
                                              m_nl, p.type().pdgcode());

          sh_weight = diff_width / width_tot;
          break;
        }
        default:
          throw std::runtime_error("Error in DecayActionFinderDilepton");
      }
      if (sh_weight > 0.0) {  // decays that can happen
        auto act = make_unique<DecayActionDilepton>(p, 0.f, sh_weight,
                                                                 dilepton_mass);
        act->add_decay(std::move(mode));
        actions.emplace_back(std::move(act));
      }
    }
  }

  return std::move(actions);
}

}  // namespace Smash
