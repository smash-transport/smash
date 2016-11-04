/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionhyperonpion.h"

#include "include/cxx14compat.h"
#include "include/parametrizations.h"
#include "include/pdgcode_constants.h"

namespace Smash {

void ScatterActionHyperonPion::format_debug_output(std::ostream &out) const {
  out << "Hyperon-Pion  ";
  ScatterAction::format_debug_output(out);
}

CollisionBranchList ScatterActionHyperonPion::two_to_two_cross_sections() {
  const ParticleType &type_particle_a = incoming_particles_[0].type();
  const ParticleType &type_particle_b = incoming_particles_[1].type();

  CollisionBranchList process_list = two_to_two_inel(type_particle_a,
                                                     type_particle_b);

  return process_list;
}

CollisionBranchList ScatterActionHyperonPion::two_to_two_inel(
                            const ParticleType &type_particle_a,
                            const ParticleType &type_particle_b) {
  CollisionBranchList process_list;

  const ParticleType &type_hyperon =
      type_particle_a.pdgcode().is_hyperon() ? type_particle_a : type_particle_b;
  const ParticleType &type_pion =
      type_particle_a.pdgcode().is_hyperon() ? type_particle_b : type_particle_a;

  const auto pdg_hyperon = type_hyperon.pdgcode().code();
  const auto pdg_pion = type_pion.pdgcode().code();

  const double s = mandelstam_s();
  const double sqrts = sqrt_s();

  switch (pack(pdg_hyperon, pdg_pion)) {
    case pack(pdg::Sigma_z, pdg::pi_m): {
      const auto& neutron = ParticleType::find(pdg::n);
      const auto& kaon = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, neutron, kaon)
                               * kminusn_piminussigma0(sqrts); },
                  sqrts, neutron, kaon);
      break;
    }
    case pack(pdg::Sigma_m, pdg::pi_z): {
      const auto& neutron = ParticleType::find(pdg::n);
      const auto& kaon = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, neutron, kaon)
                               * kminusn_pi0sigmaminus(sqrts); },
                  sqrts, neutron, kaon);
      break;
    }
    case pack(pdg::Lambda, pdg::pi_m): {
      const auto& neutron = ParticleType::find(pdg::n);
      const auto& kaon = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, neutron, kaon)
                               * kminusn_piminuslambda(sqrts); },
                  sqrts, neutron, kaon);
      break;
    }
    case pack(pdg::Sigma_z, pdg::pi_z): {
      const auto& proton = ParticleType::find(pdg::p);
      const auto& kaon = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, proton, kaon)
                               * kminusp_pi0sigma0(sqrts); },
                  sqrts, proton, kaon);
      break;
    }
    case pack(pdg::Sigma_m, pdg::pi_p): {
      const auto& proton = ParticleType::find(pdg::p);
      const auto& kaon = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, proton, kaon)
                               * kminusp_piplussigmaminus(sqrts); },
                  sqrts, proton, kaon);
      break;
    }
    case pack(pdg::Lambda, pdg::pi_z): {
      const auto& proton = ParticleType::find(pdg::p);
      const auto& kaon = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, proton, kaon)
                               * kminusp_pi0lambda(sqrts); },
                  sqrts, proton, kaon);
      break;
    }
    case pack(pdg::Sigma_p, pdg::pi_m): {
      const auto& proton = ParticleType::find(pdg::p);
      const auto& kaon = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, proton, kaon)
                               * kminusp_piminussigmaplus(sqrts); },
                  sqrts, proton, kaon);
      break;
    }
    default:
      break;
  }

  return process_list;
}

}  // namespace Smash
