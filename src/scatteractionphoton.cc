/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <fstream>
#include <iostream>

#include "include/scatteractionphoton.h"

#include "include/angles.h"
#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/integrate.h"
#include "include/kinematics.h"
#include "include/particletype.h"
#include "include/pdgcode.h"
#include "include/photoncrosssections.h"
#include "include/pow.h"
#include "include/random.h"
#include "include/tabulation.h"

using std::atan;
using std::pow;
using std::sqrt;

namespace Smash {

std::unique_ptr<Tabulation> tabulation_pi_pi_rho0 = nullptr;
std::unique_ptr<Tabulation> tabulation_pi0_pi_rho = nullptr;

void ScatterActionPhoton::generate_final_state() {
  /* Decide for a particular final state. */
  const CollisionBranch *proc = choose_channel<CollisionBranch>(
      collision_channels_photons_, cross_section_photons_);
  process_type_ = proc->get_type();
  outgoing_particles_ = proc->particle_list();

  /* The production point of the new particles.  */
  FourVector middle_point = get_interaction_point();

  /* 2->2 inelastic scattering */
  /* Sample the particle momenta in CM system. */
  const std::pair<double, double> masses = sample_masses();
  const double m1 = incoming_particles_[0].effective_mass();
  const double m2 = incoming_particles_[1].effective_mass();
  const double m3 = masses.first;
  const double s = mandelstam_s();
  const double sqrts = sqrt_s();
  std::array<double, 2> mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
  const double t1 = mandelstam_t[1];
  const double t2 = mandelstam_t[0];
  const double pcm_in = cm_momentum();
  const double pcm_out = pCM(sqrts, m3, 0.0);

  assert(t1 < t2);
  const double stepsize = (t2 - t1) / 100.0;
  for (double t = t1; t < t2; t += stepsize) {
    double diff_xsection_max =
        std::max(diff_cross_section(t, m3, t2, t1), diff_xsection_max);
  }

  double t = Random::uniform(t1, t2);
  double diff_xsection_max = 0;
  int iteration_number = 0;
  // Bug: diff_xsection_max set to 0 but then used in comparison
  do {
    t = Random::uniform(t1, t2);
    iteration_number++;
  } while (diff_cross_section(t, m3, t2, t1) <
               Random::uniform(0., diff_xsection_max) &&
           iteration_number < 100);

  // TODO(schaefer): this should be moved to kinematics.h and tested
  double costheta =
      (t - pow_int(m2, 2) +
       0.5 * (s + pow_int(m2, 2) - pow_int(m1, 2)) * (s - pow_int(m3, 2)) / s) /
      (pcm_in * (s - pow_int(m3, 2)) / sqrts);

  Angles phitheta(Random::uniform(0.0, twopi), costheta);
  outgoing_particles_[0].set_4momentum(masses.first,
                                       phitheta.threevec() * pcm_out);
  outgoing_particles_[1].set_4momentum(masses.second,
                                       -phitheta.threevec() * pcm_out);

  /* Weighing of the fractional photons */
  if (number_of_fractional_photons_ > 1) {
    weight_ = diff_cross_section(t, m3, t2, t1) * (t2 - t1) /
              (number_of_fractional_photons_ * cross_section());
  } else {
    weight_ = proc->weight() / cross_section();
  }

  /* Set positions & boost to computational frame. */
  for (ParticleData &new_particle : outgoing_particles_) {
    new_particle.set_4position(middle_point);
    new_particle.boost_momentum(-beta_cm());
  }

  /* Inlcusion of form factors:
  Usual procedure would be the multplication of the photon cross section
  with the corresponding form factor. This form factor is however energy
  dependent, such that the energy of the generated photon in the computational
  frame is a necessary to determine FF. Yet this is not directly accessible in
  ScatterActionPhoton::photon_cross_section().
  The alternative solution is to multiply the weighting factor (proportional to
  cross section) by the form factor, which is equivalent to multiplying the
  cross section directly.

  The modification is as follows:
  weight_FF = weight_noFF * FF^4
  The actual value of the form factor is determined in
  ScatterActionPhoton::form_factor */

  double E_Photon_Comp = outgoing_particles_[1].momentum()[0];

  weight_ *= pow(form_factor(E_Photon_Comp), 4);

  // Photons are not really part of the normal processes, so we have to set a
  // constant arbitrary number.
  const auto id_process = ID_PROCESS_PHOTON;
  Action::check_conservation(id_process);
}

void ScatterActionPhoton::add_dummy_hadronic_channels(
    double reaction_cross_section) {
  CollisionBranchPtr dummy_process = make_unique<CollisionBranch>(
      incoming_particles_[0].type(), incoming_particles_[1].type(),
      reaction_cross_section, ProcessType::TwoToTwo);
  add_collision(std::move(dummy_process));
}

ScatterActionPhoton::ReactionType ScatterActionPhoton::is_photon_reaction(
    const ParticleList &in) {
  if (in.size() != 2) {
    return ReactionType::no_reaction;
  }

  PdgCode a = in[0].pdgcode();
  PdgCode b = in[1].pdgcode();

  // swap so that pion should be first and there are less cases to be listed

  if (!a.is_pion()) {
    std::swap(a, b);
  }

  switch (pack(a.code(), b.code())) {
    case (pack(pdg::pi_p, pdg::pi_z)):
    case (pack(pdg::pi_z, pdg::pi_p)):
    case (pack(pdg::pi_m, pdg::pi_z)):
    case (pack(pdg::pi_z, pdg::pi_m)):
      return ReactionType::pi0_pi;
    case (pack(pdg::pi_p, pdg::rho_z)):
    case (pack(pdg::pi_m, pdg::rho_z)):
      return ReactionType::pi_rho0;
    case (pack(pdg::pi_m, pdg::rho_p)):
    case (pack(pdg::pi_p, pdg::rho_m)):
      return ReactionType::pi_rho;
    case (pack(pdg::pi_z, pdg::rho_p)):
    case (pack(pdg::pi_z, pdg::rho_m)):
      return ReactionType::pi0_rho;
    /*case(pack(pdg::pi_p, pdg::eta)):
    case(pack(pdg::pi_m, pdg::eta)):
      return ReactionType::pi_eta;*/
    case (pack(pdg::pi_p, pdg::pi_m)):
    case (pack(pdg::pi_m, pdg::pi_p)):
      return ReactionType::pi_pi;
    case (pack(pdg::pi_z, pdg::rho_z)):
      return ReactionType::pi0_rho0;
    default:
      return ReactionType::no_reaction;
  }
}

CollisionBranchList ScatterActionPhoton::photon_cross_sections() {
  CollisionBranchList process_list;
  static ParticleTypePtr rho0_particle = &ParticleType::find(pdg::rho_z);
  static ParticleTypePtr rho_plus_particle = &ParticleType::find(pdg::rho_p);
  static ParticleTypePtr rho_minus_particle = &ParticleType::find(pdg::rho_m);
  static ParticleTypePtr pi0_particle = &ParticleType::find(pdg::pi_z);
  static ParticleTypePtr pi_plus_particle = &ParticleType::find(pdg::pi_p);
  static ParticleTypePtr pi_minus_particle = &ParticleType::find(pdg::pi_m);
  static ParticleTypePtr photon_particle = &ParticleType::find(pdg::photon);

  const double m_rho = rho0_particle->mass();
  const double m_pi = pi0_particle->mass();

  /*
    const double to_mb = 0.3894;
    const double Const = 0.059;
    const double g_POR = 11.93;
    const double ma1 = 1.26;
    const double ghat = 6.4483;
    const double eta1 = 2.3920;
    const double eta2 = 1.9430;
    const double delta = -0.6426;
    const double C4 = -0.14095;
    const double Gammaa1 = 0.4;
    const double Pi = M_PI;
    double m_omega = 0.783;
    double momega = m_omega;
    double mrho = m_rho;
    double mpion = m_pi;
  */

  PhotonCrossSection<ComputationMethod::Lookup> xs_object;

  //std::cout << xs_object.s_min << " " << xs_object.s_max << " " <<
  //    xs_object.t_min << " " << xs_object.t_max << std::endl;

  ParticleData part_a = incoming_particles_[0];
  ParticleData part_b = incoming_particles_[1];

  bool pion_found = true;

  if (!part_a.type().pdgcode().is_pion()) {
    if (part_b.type().pdgcode().is_pion()) {
      part_a = incoming_particles_[1];
      part_b = incoming_particles_[0];
    } else {
      pion_found = false;
    }
  }

  if (pion_found) {
    // do a check according to incoming_particles_ and calculate the
    // cross sections (xsection) for all possible reactions

    const double s = mandelstam_s();
    double sqrts = sqrt_s();
    const double &m1 = part_a.effective_mass();
    const double &m2 = part_b.effective_mass();
    double m3 = 0.0;  // will be fixed according to reaction outcome
    ParticleTypePtr part_out = photon_particle;
    ParticleTypePtr photon_out = photon_particle;

    reac = is_photon_reaction(Action::incoming_particles());

    if (sqrts <= m1 + m2) {
      reac = ReactionType::no_reaction;
    }

    if (reac != ReactionType::no_reaction) {
      std::array<double, 2> mandelstam_t;
      double t1, t2, xsection = 0.0;

      switch (reac) {
        case ReactionType::pi_pi:
          // there are three possible reaction channels
          // the first possible reaction produces eta (removed on this branch)

          // the second possible reaction (produces rho0)
          part_out = rho0_particle;
          m3 = part_out->mass();

          if (sqrts > m3) {
            mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
            t1 = mandelstam_t[1];
            t2 = mandelstam_t[0];

            xsection = xs_object.xs_pi_pi_rho0(s);

            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          }

          // the third possible reaction (produces photon)
          part_out = photon_particle;
          m3 = 0.0;

          mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
          t1 = mandelstam_t[1];
          t2 = mandelstam_t[0];

          // replace by really small constant
          xsection = 0.0000000000001 * to_mb;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;

        case ReactionType::pi0_pi:
          if (part_a.type().pdgcode() == pdg::pi_p ||
              part_b.type().pdgcode() == pdg::pi_p) {
            part_out = rho_plus_particle;
          } else {
            part_out = rho_minus_particle;
          }
          m3 = part_out->mass();

          if (sqrts > m3) {
            mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
            t1 = mandelstam_t[1];
            t2 = mandelstam_t[0];

            xsection = xs_object.xs_pi_pi0_rho(s);

            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          }

          // dummy: just for stable rho
          // could be also in block above?
          if (part_a.type().pdgcode() == pdg::pi_p ||
              part_b.type().pdgcode() == pdg::pi_p) {
            part_out = pi_plus_particle;
          } else {
            part_out = pi_minus_particle;
          }
          m3 = 0.0;

          mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
          t1 = mandelstam_t[1];
          t2 = mandelstam_t[0];

          xsection = 0.0000000000000001;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;

        case ReactionType::pi_rho0:
          if (part_a.type().pdgcode() == pdg::pi_p) {
            part_out = pi_plus_particle;
          } else {
            part_out = pi_minus_particle;
          }

          m3 = part_out->mass();
          if (sqrts > m3) {
            mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
            t1 = mandelstam_t[1];
            t2 = mandelstam_t[0];

            xsection = xs_object.xs_pi_rho0_pi(s);

            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          } else {
            xsection = 0.0000000000001 * to_mb;
            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          }
          break;

        case ReactionType::pi_rho:
          part_out = pi0_particle;

          m3 = part_out->mass();
          if (sqrts > m3) {
            mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
            t1 = mandelstam_t[1];
            t2 = mandelstam_t[0];

            // omega:
            xsection = xs_object.xs_pi_rho_pi0(s);

            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          } else {
            xsection = 0.0000000000001 * to_mb;
            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          }
          break;

        case ReactionType::pi0_rho:
          if (part_b.type().pdgcode() == pdg::rho_p) {
            part_out = pi_plus_particle;
          } else {
            part_out = pi_minus_particle;
          }
          m3 = part_out->mass();
          if (sqrts > m3) {
            mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
            t1 = mandelstam_t[1];
            t2 = mandelstam_t[0];

            // omega:
            xsection = xs_object.xs_pi0_rho_pi(s);

            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          } else {
            xsection = 0.0000000000001 * to_mb;
            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          }
          break;

        case ReactionType::pi0_rho0:
          part_out = pi0_particle;

          m3 = part_out->mass();
          if (sqrts > m3) {
            mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
            t1 = mandelstam_t[1];
            t2 = mandelstam_t[0];

            xsection = xs_object.xs_pi0_rho0_pi0(s);

            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          } else {
            xsection = 0.0000000000001 * to_mb;
            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          }

          break;

        case ReactionType::no_reaction:
          // never reached
          break;
      }
    }
  }
  return process_list;
}

double ScatterActionPhoton::diff_cross_section(double t, double m3, double t2,
                                               double t1) const {
  const double to_mb = 0.3894;
  static const float m_rho = ParticleType::find(pdg::rho_z).mass();
  static const float m_pi = ParticleType::find(pdg::pi_z).mass();
  double s = mandelstam_s();
  double diff_xsection = 0.0;

  /*
  const double Const = 0.059;
  const double g_POR = 11.93;
  const double ma1 = 1.26;
  const double ghat = 6.4483;
  const double eta1 = 2.3920;
  const double eta2 = 1.9430;
  const double delta = -0.6426;
  const double C4 = -0.14095;
  const double Gammaa1 = 0.4;
  const double Pi = M_PI;
  double m_omega = 0.783;
  double momega = m_omega;
  double mrho = m_rho;
  double mpion = m_pi;
*/
  PhotonCrossSection<ComputationMethod::Lookup> xs_object;

  switch (reac) {
    case ReactionType::pi_pi:
      if (outgoing_particles_[0].type().pdgcode().is_rho()) {
        diff_xsection = xs_object.xs_diff_pi_pi_rho0(s, t);

      } else if (outgoing_particles_[0].type().pdgcode() == pdg::photon) {
        diff_xsection = 0.0000000000001 / to_mb / (t2 - t1);
      }
      break;
    case ReactionType::pi0_pi:
      if (outgoing_particles_[0].type().pdgcode().is_rho()) {
        diff_xsection = xs_object.xs_diff_pi_pi0_rho(s, t);

      } else if (outgoing_particles_[0].type().pdgcode().is_pion()) {
        diff_xsection = 0.0000000000001 / to_mb / (t2 - t1);
      }
      break;

    case ReactionType::pi_rho0:

      diff_xsection = xs_object.xs_diff_pi_rho0_pi(s, t);

      break;

    case ReactionType::pi_rho:

      // omega:
      diff_xsection = xs_object.xs_diff_pi_rho_pi0(s, t);

      break;
    case ReactionType::pi0_rho:

      // omega:
      diff_xsection = xs_object.xs_diff_pi0_rho_pi(s, t);

      break;

    case ReactionType::pi0_rho0:
      diff_xsection = xs_object.xs_diff_pi0_rho0_pi0(s, t);

      break;
    case ReactionType::no_reaction:
      // never reached
      break;
  }
  return diff_xsection;
  // conversionfactor to_mb in functions. this makes it analog to total cross
  // sections.
  // return diff_xsection * to_mb;
}

double ScatterActionPhoton::form_factor(double E_photon) {
  double form_factor = 1.0;
  double t_ff = 0.0;
  double Lambda = 1.0;
  switch (reac) {
      /* The form factor is assumed to be a hadronic dipole form factor which
      takes the shape of: FF = (2*Lambda^2/(2*Lambda^2 - t))^2 with
      Lambda = 1.0 GeV. t depends on the lightest possible exchange particle in
      the different channels. This could either be a pion or an omega meson. For
      the computation the parametrizations given in REF! are used. */

    case ReactionType::pi_pi:
    case ReactionType::pi0_pi:
    case ReactionType::pi_rho0:
      // case ReactionType::pi_rho:
      // case ReactionType::pi0_rho:
      t_ff = 34.5096 * pow(E_photon, 0.737) - 67.557 * pow(E_photon, 0.7584) +
             32.858 * pow(E_photon, 0.7806);
      break;
    // lightest exchange particle: omega
    case ReactionType::pi_rho:   // for omega
    case ReactionType::pi0_rho:  // for omega
    case ReactionType::pi0_rho0:
      t_ff = -61.595 * pow(E_photon, 0.9979) + 28.592 * pow(E_photon, 1.1579) +
             37.738 * pow(E_photon, 0.9317) - 5.282 * pow(E_photon, 1.3686);
      break;

    case ReactionType::no_reaction:
      // never reached
      break;
  }
  form_factor = pow(2.0 * pow(Lambda, 2) / (2.0 * pow(Lambda, 2) - t_ff), 2);
  return form_factor;
}

}  // namespace Smash
