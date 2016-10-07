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
#include "include/kinematics.h"

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

static float detailed_balance_factor(float s, const ParticleType& particle_a, const ParticleType& particle_b, const ParticleType& particle_c, const ParticleType& particle_d) {
    float spin_factor = (particle_c.spin() + 1)*(particle_d.spin() + 1);
    spin_factor /= (particle_a.spin() + 1)*(particle_b.spin() + 1);
    float symmetry_factor = (1 + (particle_a == particle_b));
    symmetry_factor /= (1 + (particle_c == particle_d));
    const float momentum_factor = pCM_sqr_from_s(s, particle_c.mass(), particle_d.mass())
        / pCM_sqr_from_s(s, particle_a.mass(), particle_b.mass());
    return spin_factor * symmetry_factor * momentum_factor;
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

  // calculate cross section
  auto add_channel
    = [&](float xsection, const ParticleType &type_a, const ParticleType &type_b) {
      const float sqrt_s_min = type_a.minimum_mass() + type_b.minimum_mass();
      if ((xsection > really_small) && (sqrts > sqrt_s_min)) {
        process_list.push_back(make_unique<CollisionBranch>(
          type_a, type_b, xsection, ProcessType::TwoToTwo));
      }
  };

  switch (pack(pdg_hyperon, pdg_pion)) {
    case pack(pdg::Sigma_z, pdg::pi_m): {
      const auto& neutron = ParticleType::find(pdg::n);
      const auto& kaon = ParticleType::find(pdg::K_m);
      const auto factor = detailed_balance_factor(s, neutron, kaon, type_hyperon, type_pion);
      add_channel(factor * kminusn_piminussigma0(sqrts), neutron, kaon);
      break;
    }
    case pack(pdg::Sigma_m, pdg::pi_z): {
      const auto& neutron = ParticleType::find(pdg::n);
      const auto& kaon = ParticleType::find(pdg::K_m);
      const auto factor = detailed_balance_factor(s, neutron, kaon, type_hyperon, type_pion);
      add_channel(factor * kminusn_pi0sigmaminus(sqrts), neutron, kaon);
      break;
    }
    case pack(pdg::Lambda, pdg::pi_m): {
      const auto& neutron = ParticleType::find(pdg::n);
      const auto& kaon = ParticleType::find(pdg::K_m);
      const auto factor = detailed_balance_factor(s, neutron, kaon, type_hyperon, type_pion);
      add_channel(factor * kminusn_piminuslambda(sqrts), neutron, kaon);
      break;
    }
    case pack(pdg::Sigma_z, pdg::pi_z): {
      const auto& proton = ParticleType::find(pdg::p);
      const auto& kaon = ParticleType::find(pdg::K_m);
      const auto factor = detailed_balance_factor(s, proton, kaon, type_hyperon, type_pion);
      add_channel(factor * kminusp_pi0sigma0(sqrts), proton, kaon);
      break;
    }
    case pack(pdg::Sigma_m, pdg::pi_p): {
      const auto& proton = ParticleType::find(pdg::p);
      const auto& kaon = ParticleType::find(pdg::K_m);
      const auto factor = detailed_balance_factor(s, proton, kaon, type_hyperon, type_pion);
      add_channel(factor * kminusp_piplussigmaminus(sqrts), proton, kaon);
      break;
    }
    case pack(pdg::Lambda, pdg::pi_z): {
      const auto& proton = ParticleType::find(pdg::p);
      const auto& kaon = ParticleType::find(pdg::K_m);
      const auto factor = detailed_balance_factor(s, proton, kaon, type_hyperon, type_pion);
      add_channel(factor * kminusp_pi0lambda(sqrts), proton, kaon);
      break;
    }
    case pack(pdg::Sigma_p, pdg::pi_m): {
      const auto& proton = ParticleType::find(pdg::p);
      const auto& kaon = ParticleType::find(pdg::K_m);
      const auto factor = detailed_balance_factor(s, proton, kaon, type_hyperon, type_pion);
      add_channel(factor * kminusp_piminussigmaplus(sqrts), proton, kaon);
      break;
    }
    default:
      break;
  }

  return process_list;
}

}  // namespace Smash
