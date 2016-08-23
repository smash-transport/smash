/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionnucleonkaon.h"

#include "include/cxx14compat.h"
#include "include/parametrizations.h"
#include "include/pdgcode_constants.h"

namespace Smash {

float ScatterActionNucleonKaon::elastic_parametrization() {
  const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();

  const PdgCode &nucleon = pdg_a.is_nucleon() ? pdg_a : pdg_b;
  const PdgCode &kaon = pdg_a.is_nucleon() ? pdg_b : pdg_a;
  assert(kaon != nucleon);

  const double s = mandelstam_s();

  float sig_el = 0.f;
  switch (nucleon.code()) {
      case pdg::p:
      switch (kaon.code()) {
        case pdg::K_p:
          sig_el = kplusp_elastic(s);
          break;
        case pdg::K_m:
          sig_el = kminusp_elastic(s);
          break;
        case pdg::K_z:
          sig_el = k0p_elastic(s);
          break;
        case pdg::Kbar_z:
          sig_el = kbar0p_elastic(s);
          break;
      }
      break;
      case pdg::n:
      switch (kaon.code()) {
        case pdg::K_p:
          sig_el = kplusn_elastic(s);
          break;
        case pdg::K_m:
          sig_el = kminusn_elastic(s);
          break;
        case pdg::K_z:
          sig_el = k0n_elastic(s);
          break;
        case pdg::Kbar_z:
          sig_el = kbar0n_elastic(s);
          break;
      }
      break;
    default:
      throw std::runtime_error("elastic cross section for antinucleon-kaon "
                               "not implemented");
  }

  if (sig_el > 0) {
    return sig_el;
  } else {
    std::stringstream ss;
    const auto name_a = incoming_particles_[0].type().name();
    const auto name_b = incoming_particles_[1].type().name();
    ss << "problem in CrossSections::elastic: a=" << name_a
       << " b=" << name_b << " j_a=" << pdg_a.spin() << " j_b="
       << pdg_b.spin() << " sigma=" << sig_el << " s=" << s;
    throw std::runtime_error(ss.str());
  }
}

void ScatterActionNucleonKaon::format_debug_output(std::ostream &out) const {
  out << "Nucleon-Kaon  ";
  ScatterAction::format_debug_output(out);
}

CollisionBranchList ScatterActionNucleonKaon::two_to_two_cross_sections() {
  const ParticleType &type_particle_a = incoming_particles_[0].type();
  const ParticleType &type_particle_b = incoming_particles_[1].type();

  CollisionBranchList process_list = two_to_two_inel(type_particle_a,
                                                     type_particle_b);

  return process_list;
}

CollisionBranchList ScatterActionNucleonKaon::two_to_two_inel(
                            const ParticleType &type_particle_a,
                            const ParticleType &type_particle_b) {
  CollisionBranchList process_list;

  const ParticleType &type_nucleon =
      type_particle_a.pdgcode().is_nucleon() ? type_particle_a : type_particle_b;
  const ParticleType &type_kaon =
      type_particle_a.pdgcode().is_nucleon() ? type_particle_b : type_particle_a;

  const auto pdg_nucleon = type_nucleon.pdgcode().code();
  const auto pdg_kaon = type_kaon.pdgcode().code();

  const double s = mandelstam_s();
  const double sqrts = sqrt_s();

  // calculate cross section
  auto add_channel
    = [&](float xsection, const ParticleType &type_a, const ParticleType &type_b) {
      if (xsection > really_small) {
        process_list.push_back(make_unique<CollisionBranch>(
          type_a, type_b, xsection, ProcessType::TwoToTwo));
      }
  };
  switch (pdg_kaon) {
    case pdg::K_m: {
      // All inelatic K- N channels here are strangeness exchange.
      switch (pdg_nucleon) {
        case pdg::p: {
          const ParticleType &type_pi0 = ParticleType::find(pdg::pi_z);
          add_channel(kminusp_piminussigmaplus(sqrts),
                      ParticleType::find(pdg::pi_m), ParticleType::find(pdg::Sigma_p));
          add_channel(kminusp_piplussigmaminus(sqrts),
                      ParticleType::find(pdg::pi_p), ParticleType::find(pdg::Sigma_m));
          add_channel(kminusp_pi0sigma0(sqrts),
                      type_pi0, ParticleType::find(pdg::Sigma_z));
          add_channel(kminusp_pi0lambda(sqrts),
                      type_pi0, ParticleType::find(pdg::Lambda));
          break;
        }
        case pdg::n: {
          const ParticleType &type_piminus = ParticleType::find(pdg::pi_m);
          add_channel(kminusn_piminussigma0(sqrts),
                      type_piminus, ParticleType::find(pdg::Sigma_z));
          add_channel(kminusn_pi0sigmaminus(sqrts),
                      ParticleType::find(pdg::pi_z), ParticleType::find(pdg::Sigma_m));
          add_channel(kminusn_piminuslambda(sqrts),
                      type_piminus, ParticleType::find(pdg::Lambda));
          break;
        }
      }
      break;
    }
    case pdg::K_p: {
      // All inelastic channels are K+ N -> K Delta -> K pi N, with identical
      // cross section.
      const auto sigma_kplusp = kplusp_inelastic(s);
      switch (pdg_nucleon) {
        case pdg::p: {
          add_channel(sigma_kplusp * 0.5, ParticleType::find(pdg::K_z),
                      ParticleType::find(pdg::Delta_pp));
          add_channel(sigma_kplusp * 0.5, ParticleType::find(pdg::K_p),
                      ParticleType::find(pdg::Delta_p));
          break;
        }
        case pdg::n: {
          add_channel(sigma_kplusp / 3, ParticleType::find(pdg::K_m),
                      ParticleType::find(pdg::Delta_pp));
          add_channel(sigma_kplusp / 3, ParticleType::find(pdg::K_z),
                      ParticleType::find(pdg::Delta_p));
          add_channel(sigma_kplusp / 3, ParticleType::find(pdg::K_p),
                      ParticleType::find(pdg::Delta_z));
          break;
        }
      }
      break;
    }
  }

  return process_list;
}

}  // namespace Smash
