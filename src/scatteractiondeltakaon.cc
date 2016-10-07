/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractiondeltakaon.h"

#include "include/cxx14compat.h"
#include "include/parametrizations.h"
#include "include/pdgcode_constants.h"

namespace Smash {

void ScatterActionDeltaKaon::format_debug_output(std::ostream &out) const {
  out << "Delta-Kaon  ";
  ScatterAction::format_debug_output(out);
}

CollisionBranchList ScatterActionDeltaKaon::two_to_two_cross_sections() {
  const ParticleType &type_particle_a = incoming_particles_[0].type();
  const ParticleType &type_particle_b = incoming_particles_[1].type();

  CollisionBranchList process_list = two_to_two_inel(type_particle_a,
                                                     type_particle_b);

  return process_list;
}

CollisionBranchList ScatterActionDeltaKaon::two_to_two_inel(
                            const ParticleType &type_particle_a,
                            const ParticleType &type_particle_b) {
  CollisionBranchList process_list;

  const ParticleType &type_delta =
      type_particle_a.pdgcode().is_Delta() ? type_particle_a : type_particle_b;
  const ParticleType &type_kaon =
      type_particle_a.pdgcode().is_Delta() ? type_particle_b : type_particle_a;

  const auto pdg_delta = type_delta.pdgcode().code();
  const auto pdg_kaon = type_kaon.pdgcode().code();

  const auto s = mandelstam_s();
  const auto sqrts = sqrt_s();

  // calculate cross section
  auto add_channel
    = [&](float xsection, const ParticleType &type_a, const ParticleType &type_b) {
      const float sqrt_s_min = type_a.minimum_mass() + type_b.minimum_mass();
      if ((xsection > really_small) && (sqrts > sqrt_s_min)) {
        process_list.push_back(make_unique<CollisionBranch>(
          type_a, type_b, xsection, ProcessType::TwoToTwo));
      }
  };

  switch (pack(pdg_delta, pdg_kaon)) {
    case pack(pdg::Delta_pp, pdg::K_z):
    case pack(pdg::Delta_p, pdg::K_p): {
      const auto& proton = ParticleType::find(pdg::p);
      const auto& kaon = ParticleType::find(pdg::K_p);
      const auto factor = detailed_balance_factor(s, proton, kaon, type_delta, type_kaon);
      add_channel(0.5 * factor * kplusp_inelastic(s), proton, kaon);
      break;
    }
    case pack(pdg::Delta_p, pdg::K_z):
    case pack(pdg::Delta_z, pdg::K_p): {
      const auto& neutron = ParticleType::find(pdg::n);
      const auto& kaon_p = ParticleType::find(pdg::K_p);
      const auto factor1 = detailed_balance_factor(s, neutron, kaon_p, type_delta, type_kaon);
      add_channel(0.5 * factor1 * kplusn_inelastic(s), neutron, kaon_p);

      const auto& proton = ParticleType::find(pdg::p);
      const auto& kaon_z = ParticleType::find(pdg::K_z);
      const auto factor2 = detailed_balance_factor(s, proton, kaon_z, type_delta, type_kaon);
      add_channel(0.5 * factor2 * kplusp_inelastic(s), proton, kaon_z);
      break;
    }
    case pack(pdg::Delta_z, pdg::K_z):
    case pack(pdg::Delta_m, pdg::K_p): {
      const auto& neutron = ParticleType::find(pdg::n);
      const auto& kaon = ParticleType::find(pdg::K_z);
      const auto factor = detailed_balance_factor(s, neutron, kaon, type_delta, type_kaon);
      add_channel(0.5 * factor * kplusn_inelastic(s), neutron, kaon);
      break;
    }
    default:
      break;
  }

  return process_list;
}

}  // namespace Smash
