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
#include "include/decaymodes.h"
#include "include/experimentparameters.h"
#include "include/logging.h"
#include "include/fourvector.h"
#include "include/particles.h"
#include "include/random.h"

#include <typeinfo>
#include <iostream>
#include <fstream>

namespace Smash {

ActionList DecayActionsFinderDilepton::find_actions_in_cell(
      const ParticleList &search_list,
      float dt) const {

  ActionList actions;

  for (const auto &p : search_list) {
    if (p.type().decay_modes().decay_mode_list().size() == 0) {
      continue;
    }

    // work in progress // #CleanUp

    /* current (first) implementation: everything happens in finder
     * It slows down smash probably because of the extra output
     * file.
     */

    float inv_gamma = p.inverse_gamma();

  // for three body decays partial wodth = on_shell width (dummy)

    DecayBranchList dil_modes =
                  p.type().get_partial_widths_dilepton(p.effective_mass());

    for (DecayBranchPtr & mode : dil_modes) {
      float sh_weight = 0.0;
      float dilepton_mass = 0.f;

      switch (mode->particle_number()) {
        case 2:
          {
          float partial_width = mode->weight();
          // SHINNING as described in \iref{Schmidt:2008hm}, chapter 2D
          sh_weight = dt * inv_gamma * partial_width;
          break;
          }
        case 3:
          {
          // find the non lepton particle position
          int non_lepton_position = -1;
          for (int i=0; i<3; ++i) {
            if (!(mode->particle_types()[i]->is_lepton())) {
              non_lepton_position = i;
              break;
            }
          }
          float m_nl = mode->particle_types()[non_lepton_position]->mass();  // mass of non-lepton final state particle
          float m_l = mode->particle_types()[(non_lepton_position+1)%3]->mass(); // mass of leptons in final state
          // randomly select a mass
          dilepton_mass = Random::uniform(2*m_l,p.effective_mass());
          float diff_width = mode->type().diff_width(p.effective_mass(), dilepton_mass, m_nl,
                                              p.type().pdgcode());  // #CleanUp
          sh_weight = dt * inv_gamma * diff_width;
          break;
          }
        default:
          throw std::runtime_error("Error in DecayActionFinderDilepton");
      }

      auto act = make_unique<DecayActionDilepton>(p, 0.f, sh_weight, dilepton_mass);
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
      continue;      /* particle doesn't decay */
    }

    float inv_gamma = p.inverse_gamma();

    DecayBranchList dil_modes =
                  p.type().get_partial_widths_dilepton(p.effective_mass());

    // total decay width, also hadronic decays
    const float width_tot = total_weight<DecayBranch>(
                               p.type().get_partial_widths(p.effective_mass()));

    for (DecayBranchPtr & mode : dil_modes) {
      float sh_weight = 0.0;
      float dilepton_mass = 0.f;

      switch (mode->particle_number()) {
        case 2:
          {
          float partial_width = mode->weight();
          sh_weight = partial_width * inv_gamma / width_tot;
          break;
          }
        case 3:
          {
          // find the non lepton particle position
          int non_lepton_position = -1;
          for (int i=0; i<3; ++i) {
            if (!(mode->particle_types()[i]->is_lepton())) {
              non_lepton_position = i;
              break;
            }
          }
          float m_nl = mode->particle_types()[non_lepton_position]->mass();  // mass of non-lepton final state particle
          float m_l = mode->particle_types()[(non_lepton_position+1)%3]->mass(); // mass of leptons in final state
          // randomly select a mass
          dilepton_mass = Random::uniform(2*m_l,p.effective_mass());
          float diff_width = mode->type().diff_width(p.effective_mass(), dilepton_mass, m_nl,
                                              p.type().pdgcode());  // #CleanUp
          sh_weight = diff_width * inv_gamma / width_tot;;
          break;
          }
        default:
          throw std::runtime_error("Error in DecayActionFinderDilepton");
      }

      auto act = make_unique<DecayActionDilepton>(p, 0.f, sh_weight, dilepton_mass);
      act->add_decay(std::move(mode));
      actions.emplace_back(std::move(act));
    }
  }

  return std::move(actions);
}

}  // namespace Smash
