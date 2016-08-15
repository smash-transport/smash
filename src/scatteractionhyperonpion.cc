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

  const double sqrts = sqrt_s();

  // calculate cross section
  auto add_channel
    = [&](float xsection, const ParticleType &type_a, const ParticleType &type_b) {
      if (xsection > really_small) {
        process_list.push_back(make_unique<CollisionBranch>(
          type_a, type_b, xsection, ProcessType::TwoToTwo));
      }
  };

  switch (pack(pdg_hyperon, pdg_pion)) {
    case pack(pdg::Sigma_z, pdg::pi_m):
      add_channel(kminusn_piminussigma0(sqrts),
                  ParticleType::find(pdg::n), ParticleType::find(pdg::K_m));
      break;
    /*
    case pack(pdg::Sigma_m, pdg::pi_z):
      add_channel(kminusn_pi0sigmaminus(sqrts),
                  ParticleType::find(pdg::n), ParticleType::find(pdg::K_m));
      break;
    case pack(pdg::Lambda, pdg::pi_m):
      add_channel(kminusn_piminuslambda(sqrts),
                  ParticleType::find(pdg::n), ParticleType::find(pdg::K_m));
      break;
    case pack(pdg::Sigma_z, pdg::pi_z):
      add_channel(kminusp_pi0sigma0(sqrts),
                  ParticleType::find(pdg::p), ParticleType::find(pdg::K_m));
      break;
    case pack(pdg::Sigma_m, pdg::pi_p):
      add_channel(kminusp_piplussigmaminus(sqrts),
                  ParticleType::find(pdg::p), ParticleType::find(pdg::K_m));
      break;
    case pack(pdg::Lambda, pdg::pi_z):
      add_channel(kminusp_pi0lambda(sqrts),
                  ParticleType::find(pdg::p), ParticleType::find(pdg::K_m));
      break;
    case pack(pdg::Sigma_p, pdg::pi_m):
      add_channel(kminusp_piminussigmaplus(sqrts),
                  ParticleType::find(pdg::p), ParticleType::find(pdg::K_m));
      break;
    */
    default:
      break;
  }

  return process_list;
}

}  // namespace Smash
