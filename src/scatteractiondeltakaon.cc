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

  const double s = mandelstam_s();
  const double sqrts = sqrt_s();

  // The cross sections are determined from the backward reactions via
  // detailed balance. The same isospin factors as for the backward reaction
  // are used.
  switch (pack(pdg_delta, pdg_kaon)) {
    case pack(pdg::Delta_pp, pdg::K_z):
    case pack(pdg::Delta_p, pdg::K_p): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_RK(s,
                               type_p, type_K_p, type_delta, type_kaon)
                               * kplusn_ratios.get_ratio(type_p, type_K_p, type_kaon, type_delta)
                               * kplusp_inelastic(s); },
                  sqrts, type_p, type_K_p);
      break;
    }
    case pack(pdg::Delta_p, pdg::K_z):
    case pack(pdg::Delta_z, pdg::K_p): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_RK(s,
                               type_n, type_K_p, type_delta, type_kaon)
                               * kplusn_ratios.get_ratio(type_n, type_K_p, type_kaon, type_delta)
                               * kplusn_inelastic(s); },
                  sqrts, type_n, type_K_p);

      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_RK(s,
                               type_p, type_K_z, type_delta, type_kaon)
                               * kplusn_ratios.get_ratio(type_p, type_K_z, type_kaon, type_delta)
                               * kplusp_inelastic(s); },
                  sqrts, type_p, type_K_z);
      break;
    }
    case pack(pdg::Delta_z, pdg::K_z):
    case pack(pdg::Delta_m, pdg::K_p): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_RK(s,
                               type_n, type_K_z, type_delta, type_kaon)
                               * kplusn_ratios.get_ratio(type_n, type_K_z, type_kaon, type_delta)
                               * kplusn_inelastic(s); },
                  sqrts, type_n, type_K_z);
      break;
    }
    default:
      break;
  }

  return process_list;
}

}  // namespace Smash
