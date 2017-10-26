/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionnucleonkaon.h"

#include "include/clebschgordan.h"
#include "include/cxx14compat.h"
#include "include/parametrizations.h"
#include "include/pdgcode_constants.h"

namespace Smash {

double ScatterActionNucleonKaon::elastic_parametrization() {
  const PdgCode& pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode& pdg_b = incoming_particles_[1].type().pdgcode();

  const PdgCode& nucleon = pdg_a.is_nucleon() ? pdg_a : pdg_b;
  const PdgCode& kaon = pdg_a.is_nucleon() ? pdg_b : pdg_a;
  assert(kaon != nucleon);

  const double s = mandelstam_s();

  double sig_el = 0.;
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
    case -pdg::p:
      switch (kaon.code()) {
        case pdg::K_p:
          sig_el = kminusp_elastic(s);
          break;
        case pdg::K_m:
          sig_el = kplusp_elastic(s);
          break;
        case pdg::K_z:
          sig_el = kbar0p_elastic(s);
          break;
        case pdg::Kbar_z:
          sig_el = k0p_elastic(s);
          break;
      }
      break;
    case -pdg::n:
      switch (kaon.code()) {
        case pdg::K_p:
          sig_el = kminusn_elastic(s);
          break;
        case pdg::K_m:
          sig_el = kplusn_elastic(s);
          break;
        case pdg::K_z:
          sig_el = kbar0n_elastic(s);
          break;
        case pdg::Kbar_z:
          sig_el = k0n_elastic(s);
          break;
      }
      break;
    default:
      throw std::runtime_error(
          "elastic cross section for antinucleon-kaon "
          "not implemented");
  }

  if (sig_el > 0) {
    return sig_el;
  } else {
    std::stringstream ss;
    const auto name_a = incoming_particles_[0].type().name();
    const auto name_b = incoming_particles_[1].type().name();
    ss << "problem in CrossSections::elastic: a=" << name_a << " b=" << name_b
       << " j_a=" << pdg_a.spin() << " j_b=" << pdg_b.spin()
       << " sigma=" << sig_el << " s=" << s;
    throw std::runtime_error(ss.str());
  }
}

void ScatterActionNucleonKaon::format_debug_output(std::ostream& out) const {
  out << "Nucleon-Kaon  ";
  ScatterAction::format_debug_output(out);
}

CollisionBranchList ScatterActionNucleonKaon::two_to_two_cross_sections() {
  const ParticleType& a = incoming_particles_[0].type();
  const ParticleType& b = incoming_particles_[1].type();
  const ParticleType& type_nucleon = a.pdgcode().is_nucleon() ? a : b;
  const ParticleType& type_kaon = a.pdgcode().is_nucleon() ? b : a;

  const auto pdg_nucleon = type_nucleon.pdgcode().code();
  const auto pdg_kaon = type_kaon.pdgcode().code();

  const double s = mandelstam_s();
  const double sqrts = sqrt_s();

  // Some variable declarations for frequently used quantities
  const auto sigma_kplusp = kplusp_inelastic(s);
  const auto sigma_kplusn = kplusn_inelastic(s);

  CollisionBranchList process_list;
  switch (pdg_kaon) {
    case pdg::K_m: {
      // All inelastic K- N channels here are strangeness exchange, plus one
      // charge exchange.
      switch (pdg_nucleon) {
        case pdg::p: {
          const auto& type_n = ParticleType::find(pdg::n);
          const auto& type_pi_z = ParticleType::find(pdg::pi_z);
          const auto& type_pi_m = ParticleType::find(pdg::pi_m);
          const auto& type_pi_p = ParticleType::find(pdg::pi_p);
          const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
          const auto& type_Sigma_p = ParticleType::find(pdg::Sigma_p);
          const auto& type_Sigma_m = ParticleType::find(pdg::Sigma_m);
          const auto& type_Sigma_z = ParticleType::find(pdg::Sigma_z);
          const auto& type_Lambda = ParticleType::find(pdg::Lambda);
          add_channel(process_list,
                      [&] { return kminusp_piminussigmaplus(sqrts); }, sqrts,
                      type_pi_m, type_Sigma_p);
          add_channel(process_list,
                      [&] { return kminusp_piplussigmaminus(sqrts); }, sqrts,
                      type_pi_p, type_Sigma_m);
          add_channel(process_list, [&] { return kminusp_pi0sigma0(sqrts); },
                      sqrts, type_pi_z, type_Sigma_z);
          add_channel(process_list, [&] { return kminusp_pi0lambda(sqrts); },
                      sqrts, type_pi_z, type_Lambda);
          add_channel(process_list, [&] { return kminusp_kbar0n(s); }, sqrts,
                      type_Kbar_z, type_n);
          break;
        }
        case pdg::n: {
          const auto& type_pi_z = ParticleType::find(pdg::pi_z);
          const auto& type_pi_m = ParticleType::find(pdg::pi_m);
          const auto& type_Sigma_m = ParticleType::find(pdg::Sigma_m);
          const auto& type_Sigma_z = ParticleType::find(pdg::Sigma_z);
          const auto& type_Lambda = ParticleType::find(pdg::Lambda);
          add_channel(process_list,
                      [&] { return kminusn_piminussigma0(sqrts); }, sqrts,
                      type_pi_m, type_Sigma_z);
          add_channel(process_list,
                      [&] { return kminusn_pi0sigmaminus(sqrts); }, sqrts,
                      type_pi_z, type_Sigma_m);
          add_channel(process_list,
                      [&] { return kminusn_piminuslambda(sqrts); }, sqrts,
                      type_pi_m, type_Lambda);
          break;
        }
        case -pdg::p: {
          const auto& type_K_m = ParticleType::find(pdg::K_m);
          const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
          const auto& type_Delta_pp_bar = ParticleType::find(-pdg::Delta_pp);
          const auto& type_Delta_p_bar = ParticleType::find(-pdg::Delta_p);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_Kbar_z,
                                                       type_Delta_pp_bar);
                      },
                      sqrts, type_Kbar_z, type_Delta_pp_bar);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_m,
                                                       type_Delta_p_bar);
                      },
                      sqrts, type_K_m, type_Delta_p_bar);
          break;
        }
        case -pdg::n: {
          const auto& type_K_m = ParticleType::find(pdg::K_m);
          const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
          const auto& type_Delta_p_bar = ParticleType::find(-pdg::Delta_p);
          const auto& type_Delta_z_bar = ParticleType::find(-pdg::Delta_z);
          const auto& type_p_bar = ParticleType::find(-pdg::p);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_Kbar_z,
                                                       type_Delta_p_bar);
                      },
                      sqrts, type_Kbar_z, type_Delta_p_bar);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_m,
                                                       type_Delta_z_bar);
                      },
                      sqrts, type_K_m, type_Delta_z_bar);
          add_channel(process_list,
                      [&] {
                        return kplusn_k0p(s);
                      },
                      sqrts, type_Kbar_z, type_p_bar);
          break;
        }
      }
      break;
    }
    case pdg::K_p: {
      // All inelastic channels are K+ N -> K Delta -> K pi N, with identical
      // cross section, weighted by the isospin factor.
      switch (pdg_nucleon) {
        case pdg::p: {
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          const auto& type_Delta_pp = ParticleType::find(pdg::Delta_pp);
          const auto& type_Delta_p = ParticleType::find(pdg::Delta_p);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_pp);
                      },
                      sqrts, type_K_z, type_Delta_pp);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_p);
                      },
                      sqrts, type_K_p, type_Delta_p);
          break;
        }
        case pdg::n: {
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          const auto& type_p = ParticleType::find(pdg::p);
          const auto& type_Delta_p = ParticleType::find(pdg::Delta_p);
          const auto& type_Delta_z = ParticleType::find(pdg::Delta_z);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_p);
                      },
                      sqrts, type_K_z, type_Delta_p);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_z);
                      },
                      sqrts, type_K_p, type_Delta_z);
          add_channel(process_list,
                      [&] {
                        return kplusn_k0p(s);
                      },
                      sqrts, type_K_z, type_p);
          break;
        }
        case -pdg::p: {
          const auto& type_n_bar = ParticleType::find(-pdg::n);
          const auto& type_pi_z = ParticleType::find(pdg::pi_z);
          const auto& type_pi_m = ParticleType::find(pdg::pi_m);
          const auto& type_pi_p = ParticleType::find(pdg::pi_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          const auto& type_Sigma_p_bar = ParticleType::find(-pdg::Sigma_p);
          const auto& type_Sigma_m_bar = ParticleType::find(-pdg::Sigma_m);
          const auto& type_Sigma_z_bar = ParticleType::find(-pdg::Sigma_z);
          const auto& type_Lambda_bar = ParticleType::find(-pdg::Lambda);
          add_channel(process_list,
                      [&] { return kminusp_piminussigmaplus(sqrts); }, sqrts,
                      type_pi_p, type_Sigma_p_bar);
          add_channel(process_list,
                      [&] { return kminusp_piplussigmaminus(sqrts); }, sqrts,
                      type_pi_m, type_Sigma_m_bar);
          add_channel(process_list, [&] { return kminusp_pi0sigma0(sqrts); },
                      sqrts, type_pi_z, type_Sigma_z_bar);
          add_channel(process_list, [&] { return kminusp_pi0lambda(sqrts); },
                      sqrts, type_pi_z, type_Lambda_bar);
          add_channel(process_list, [&] { return kminusp_kbar0n(s); }, sqrts,
                      type_K_z, type_n_bar);
          break;
        }
        case -pdg::n: {
          const auto& type_pi_z = ParticleType::find(pdg::pi_z);
          const auto& type_pi_p = ParticleType::find(pdg::pi_p);
          const auto& type_Sigma_m_bar = ParticleType::find(-pdg::Sigma_m);
          const auto& type_Sigma_z_bar = ParticleType::find(-pdg::Sigma_z);
          const auto& type_Lambda_bar = ParticleType::find(-pdg::Lambda);
          add_channel(process_list,
                      [&] { return kminusn_piminussigma0(sqrts); }, sqrts,
                      type_pi_p, type_Sigma_z_bar);
          add_channel(process_list,
                      [&] { return kminusn_pi0sigmaminus(sqrts); }, sqrts,
                      type_pi_z, type_Sigma_m_bar);
          add_channel(process_list,
                      [&] { return kminusn_piminuslambda(sqrts); }, sqrts,
                      type_pi_p, type_Lambda_bar);
          break;
        }
      }
      break;
    }
    case pdg::K_z: {
      // K+ and K0 have the same isospin projection, they are assumed to have
      // the same cross section here.

      switch (pdg_nucleon) {
        case pdg::p: {
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          const auto& type_n = ParticleType::find(pdg::n);
          const auto& type_Delta_p = ParticleType::find(pdg::Delta_p);
          const auto& type_Delta_z = ParticleType::find(pdg::Delta_z);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_p);
                      },
                      sqrts, type_K_z, type_Delta_p);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_z);
                      },
                      sqrts, type_K_p, type_Delta_z);
          add_channel(process_list,
                      [&] {
                        return kplusn_k0p(s) *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_n);
                      },
                      sqrts, type_K_p, type_n);
          break;
        }
        case pdg::n: {
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          const auto& type_Delta_z = ParticleType::find(pdg::Delta_z);
          const auto& type_Delta_m = ParticleType::find(pdg::Delta_m);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_z);
                      },
                      sqrts, type_K_z, type_Delta_z);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_m);
                      },
                      sqrts, type_K_p, type_Delta_m);
          break;
        }
        case -pdg::n: {
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_p_bar = ParticleType::find(-pdg::p);
          add_channel(process_list, [&] { return kminusp_kbar0n(s); }, sqrts,
                      type_K_p, type_p_bar);
          break;
        }
      }
      break;
    }
    case pdg::Kbar_z:
      switch (pdg_nucleon) {
        case pdg::n: {
          const auto& type_p = ParticleType::find(pdg::p);
          const auto& type_K_m = ParticleType::find(pdg::K_m);
          add_channel(process_list, [&] { return kminusp_kbar0n(s); }, sqrts,
                      type_K_m, type_p);
          break;
        }
        case -pdg::p: {
          const auto& type_K_m = ParticleType::find(pdg::K_m);
          const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
          const auto& type_Delta_p_bar = ParticleType::find(-pdg::Delta_p);
          const auto& type_Delta_z_bar = ParticleType::find(-pdg::Delta_z);
          const auto& type_n_bar = ParticleType::find(-pdg::n);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_Kbar_z,
                                                       type_Delta_p_bar);
                      },
                      sqrts, type_Kbar_z, type_Delta_p_bar);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_m,
                                                       type_Delta_z_bar);
                      },
                      sqrts, type_K_m, type_Delta_z_bar);
          add_channel(process_list,
                      [&] {
                        return kplusn_k0p(s) *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_m, type_n_bar);
                      },
                      sqrts, type_K_m, type_n_bar);
          break;
        }
        case -pdg::n: {
          const auto& type_K_m = ParticleType::find(pdg::K_m);
          const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
          const auto& type_Delta_z_bar = ParticleType::find(-pdg::Delta_z);
          const auto& type_Delta_m_bar = ParticleType::find(-pdg::Delta_m);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_Kbar_z,
                                                       type_Delta_z_bar);
                      },
                      sqrts, type_Kbar_z, type_Delta_z_bar);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_m,
                                                       type_Delta_m_bar);
                      },
                      sqrts, type_K_m, type_Delta_m_bar);
          break;
        }
      }
      break;
  }

  return process_list;
}

}  // namespace Smash
