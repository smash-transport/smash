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
  const ParticleType &a = incoming_particles_[0].type();
  const ParticleType &b = incoming_particles_[1].type();
  const ParticleType &type_hyperon = a.pdgcode().is_hyperon() ? a : b;
  const ParticleType &type_pion    = a.pdgcode().is_hyperon() ? b : a;

  const auto pdg_hyperon = type_hyperon.pdgcode().code();
  const auto pdg_pion = type_pion.pdgcode().code();

  const double s = mandelstam_s();
  const double sqrts = sqrt_s();

  CollisionBranchList process_list;

  switch (pack(pdg_hyperon, pdg_pion)) {
    case pack(pdg::Sigma_z, pdg::pi_m): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_n, type_K_m)
                               * kminusn_piminussigma0(sqrts); },
                  sqrts, type_n, type_K_m);
      break;
    }
    case pack(-pdg::Sigma_z, pdg::pi_p): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_n_bar, type_K_p)
                               * kminusn_piminussigma0(sqrts); },
                  sqrts, type_n_bar, type_K_p);
      break;
    }
    case pack(pdg::Sigma_m, pdg::pi_z): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_n, type_K_m)
                               * kminusn_pi0sigmaminus(sqrts); },
                  sqrts, type_n, type_K_m);
      break;
    }
    case pack(-pdg::Sigma_m, pdg::pi_z): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_n_bar, type_K_p)
                               * kminusn_pi0sigmaminus(sqrts); },
                  sqrts, type_n_bar, type_K_p);
      break;
    }
    case pack(pdg::Lambda, pdg::pi_m): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_n, type_K_m)
                               * kminusn_piminuslambda(sqrts); },
                  sqrts, type_n, type_K_m);
      break;
    }
    case pack(-pdg::Lambda, pdg::pi_p): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_n_bar, type_K_p)
                               * kminusn_piminuslambda(sqrts); },
                  sqrts, type_n_bar, type_K_p);
      break;
    }
    case pack(pdg::Sigma_z, pdg::pi_z): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_p, type_K_m)
                               * kminusp_pi0sigma0(sqrts); },
                  sqrts, type_p, type_K_m);
      break;
    }
    case pack(-pdg::Sigma_z, pdg::pi_z): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_p_bar, type_K_p)
                               * kminusp_pi0sigma0(sqrts); },
                  sqrts, type_p_bar, type_K_p);
      break;
    }
    case pack(pdg::Sigma_m, pdg::pi_p): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_p, type_K_m)
                               * kminusp_piplussigmaminus(sqrts); },
                  sqrts, type_p, type_K_m);
      break;
    }
    case pack(-pdg::Sigma_m, pdg::pi_m): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_p_bar, type_K_p)
                               * kminusp_piplussigmaminus(sqrts); },
                  sqrts, type_p_bar, type_K_p);
      break;
    }
    case pack(pdg::Lambda, pdg::pi_z): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_p, type_K_m)
                               * kminusp_pi0lambda(sqrts); },
                  sqrts, type_p, type_K_m);
      break;
    }
    case pack(-pdg::Lambda, pdg::pi_z): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_p_bar, type_K_p)
                               * kminusp_pi0lambda(sqrts); },
                  sqrts, type_p_bar, type_K_p);
      break;
    }
    case pack(pdg::Sigma_p, pdg::pi_m): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_p, type_K_m)
                               * kminusp_piminussigmaplus(sqrts); },
                  sqrts, type_p, type_K_m);
      break;
    }
    case pack(-pdg::Sigma_p, pdg::pi_p): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] { return detailed_balance_factor_stable(s,
                               type_hyperon, type_pion, type_p_bar, type_K_p)
                               * kminusp_piminussigmaplus(sqrts); },
                  sqrts, type_p_bar, type_K_p);
      break;
    }
    default:
      break;
  }

  return process_list;
}

}  // namespace Smash
