/*
 *
 *    Copyright (c) 2016-2017
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
#include "include/pow.h"
#include "include/random.h"
#include "include/tabulation.h"

using std::sqrt;
using std::pow;
using std::atan;

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
        std::max(diff_cross_section(t, m3), diff_xsection_max);
  }

  double t = Random::uniform(t1, t2);
  double diff_xsection_max = 0;
  int iteration_number = 0;
  do {
    t = Random::uniform(t1, t2);
    iteration_number++;
  } while (diff_cross_section(t, m3) < Random::uniform(0., diff_xsection_max) &&
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
    weight_ = diff_cross_section(t, m3) * (t2 - t1) /
              (number_of_fractional_photons_ * cross_section());
  } else {
    weight_ = proc->weight() / cross_section();
  }
  /* Set positions & boost to computational frame. */
  for (ParticleData &new_particle : outgoing_particles_) {
    new_particle.set_4position(middle_point);
    new_particle.boost_momentum(-beta_cm());
  }
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
    case (pack(pdg::pi_p, pdg::eta)):
    case (pack(pdg::pi_m, pdg::eta)):
      return ReactionType::pi_eta;
    case (pack(pdg::pi_p, pdg::pi_m)):
    case (pack(pdg::pi_m, pdg::pi_p)):
      return ReactionType::pi_pi;
    default:
      return ReactionType::no_reaction;
  }
}

CollisionBranchList ScatterActionPhoton::photon_cross_sections() {
  CollisionBranchList process_list;
  ParticleTypePtr rho0_particle = &ParticleType::find(pdg::rho_z);
  ParticleTypePtr rho_plus_particle = &ParticleType::find(pdg::rho_p);
  ParticleTypePtr rho_minus_particle = &ParticleType::find(pdg::rho_m);
  ParticleTypePtr eta_particle = &ParticleType::find(pdg::eta);
  ParticleTypePtr pi0_particle = &ParticleType::find(pdg::pi_z);
  ParticleTypePtr pi_plus_particle = &ParticleType::find(pdg::pi_p);
  ParticleTypePtr pi_minus_particle = &ParticleType::find(pdg::pi_m);
  ParticleTypePtr photon_particle = &ParticleType::find(pdg::photon);
  const double m_rho = rho0_particle->mass();
  const double m_rho_2 = pow_int(m_rho, 2);
  const double m_pi = pi0_particle->mass();
  const double m_pi_2 = pow_int(m_pi, 2);
  const double m_eta = eta_particle->mass();
  const double m_eta_2 = pow_int(m_eta, 2);
  const double gamma_rho_tot = rho0_particle->width_at_pole();
  const double g_rho_2 =
      24 * twopi * gamma_rho_tot * pow_int(m_rho, 2) /
      pow(pow_int(m_rho, 2) - 4 * pow_int(m_pi, 2), 3.0 / 2.0);
  const double to_mb = 0.3894;

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
    const double p_cm_2 = cm_momentum_squared();
    ParticleTypePtr part_out = photon_particle;
    ParticleTypePtr photon_out = photon_particle;

    reac = is_photon_reaction(Action::incoming_particles());

    if (sqrts <= m1 + m2) {
      reac = ReactionType::no_reaction;
    }

    if (reac != ReactionType::no_reaction) {
      std::array<double, 2> mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
      double t1 = mandelstam_t[1];
      double t2 = mandelstam_t[0];

      double u1 = pow_int(m1, 2) + pow_int(m2, 2) + pow_int(m3, 2) - s - t1;
      double u2 = pow_int(m1, 2) + pow_int(m2, 2) + pow_int(m3, 2) - s - t2;

      double e, I0, I1;
      double xsection = 0.0;

      Integrator1dMonte integrate;

      switch (reac) {
        case ReactionType::pi_pi:  // there are three possible reaction channels
          // the first possible reaction has part_out = photon_particle with
          // m3 = 0, which is the default declared above

          xsection = twopi * pow_int(alpha, 2) / (s * p_cm_2);
          t1 += -m_pi_2;  // is t+
          t2 += -m_pi_2;
          u1 = -s - t1;
          u2 = -s - t2;
          e = t2 - t1 +
              2 * m_pi_2 * ((1 - 2 * m_pi_2 / s) * std::log(t2 / t1) +
                            m_pi_2 * (t2 - t1) / (t1 * t2));
          e += 2 * m_pi_2 * ((1 - 2 * m_pi_2 / s) * std::log(u1 / u2) +
                             m_pi_2 * (u1 - u2) / (u1 * u2));
          xsection = xsection * e * to_mb;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));

          // now the second possible reaction (produces eta)
          part_out = eta_particle;
          m3 = part_out->mass();

          if (sqrts > m3) {
            mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
            t1 = mandelstam_t[1];
            t2 = mandelstam_t[0];
            u1 = pow_int(m1, 2) + pow_int(m2, 2) + pow_int(m3, 2) - s - t1;
            u2 = pow_int(m1, 2) + pow_int(m2, 2) + pow_int(m3, 2) - s - t2;
            xsection = M_PI * alpha * 4.7 * pow(m_rho, 4) /
                       (pow_int(s - m_rho_2, 2) +
                        pow_int(gamma_rho_tot, 2) * m_rho_2) /
                       (16 * m_eta_2 * pow_int(m_rho, 4) * s * p_cm_2);
            xsection = xsection *
                       ((2 * m_pi_2 + m_eta_2 - s) * s / 2 *
                            (pow_int(t2, 2) - pow_int(t1, 2)) -
                        s / 3 * (pow_int(t2, 3) - pow_int(t1, 3)) -
                        (t2 - t1) * m_pi_2 *
                            (pow_int(m_eta, 4) + s * (m_pi_2 - m_eta_2))) *
                       to_mb;
            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          }

          // and the third possible reaction (produces rho0)

          part_out = rho0_particle;
          m3 = part_out->mass();

          if (gamma_rho_tot > really_small) {
            if (tabulation_pi_pi_rho0 == nullptr) {
              tabulation_pi_pi_rho0 = make_unique<Tabulation>(
                  2. * m_pi, 15. - 2. * m_pi, num_tab_pts_, [&](double sqrts1) {
                    return integrate(2. * m_pi, sqrts1, [&](double M) {
                      return pi_pi_rho0(M, pow_int(sqrts1, 2)) *
                             part_out->spectral_function(M);
                    });
                  });
            }
            xsection = tabulation_pi_pi_rho0->get_value_linear(sqrts);
          } else {
            xsection = pi_pi_rho0(m3, pow_int(sqrts, 2));
          }
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

          if (gamma_rho_tot > really_small) {
            if (tabulation_pi0_pi_rho == nullptr) {
              tabulation_pi0_pi_rho = make_unique<Tabulation>(
                  2. * m_pi, 15. - 2. * m_pi, num_tab_pts_, [&](double sqrts1) {
                    return integrate(2. * m_pi, sqrts1, [&](double M) {
                      return pi_pi0_rho(M, pow_int(sqrts1, 2)) *
                             part_out->spectral_function(M);
                    });
                  });
            }
            xsection = tabulation_pi0_pi_rho->get_value_linear(sqrts);
          } else {
            xsection = pi_pi0_rho(m3, pow_int(sqrts, 2));
          }
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

          xsection = alpha * g_rho_2 / (12 * s * p_cm_2);
          t1 += -m_pi_2;
          t2 += -m_pi_2;
          xsection = xsection * (2 * (t2 - t1) -
                                 s * (pow_int(m2, 2) - 4 * m_pi_2) /
                                     pow_int(s - m_pi_2, 2) * (t2 - t1) -
                                 (pow_int(m2, 2) - 4 * m_pi_2) *
                                     ((s - pow_int(m2, 2) + m_pi_2) /
                                          (s - m_pi_2) * std::log(t2 / t1) +
                                      m_pi_2 * (t2 - t1) / (t1 * t2))) *
                     to_mb;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;

        case ReactionType::pi_rho:
          part_out = pi0_particle;
          m3 = part_out->mass();

          xsection = -alpha * g_rho_2 / (48 * s * p_cm_2);
          t1 += -m_pi_2;
          t2 += -m_pi_2;
          u1 += -m_rho_2;
          u2 += -m_rho_2;
          e = 4 * (pow_int(m2, 2) - 4 * m_pi_2) *
              (pow_int(m2, 2) * (t2 - t1) / (u1 * u2) +
               m_pi_2 * (t2 - t1) / (t2 * t1) + std::log(u1 / u2 * t2 / t1) -
               pow_int(m2, 2) / (s - m_pi_2) * std::log(t2 / t1 * u1 / u2));
          e += (s - m_pi_2) * (3.0 + (s - m_pi_2) / pow_int(m2, 2)) *
               std::log(u1 / u2);
          e += (t2 - t1) *
               (s / pow_int(m2, 2) - 0.5 - pow_int(s - m_pi_2, 2) / (u1 * u2));
          xsection = xsection * e * to_mb;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;

        case ReactionType::pi0_rho:
          if (part_b.type().pdgcode() == pdg::rho_p) {
            part_out = pi_plus_particle;
          } else {
            part_out = pi_minus_particle;
          }
          m3 = part_out->mass();

          xsection = alpha * g_rho_2 / (48 * s * p_cm_2);
          u1 += -pow_int(m2, 2);  // is u+
          u2 += -pow_int(m2, 2);
          e = (t2 - t1) *
              (4.5 - s / pow_int(m2, 2) -
               4 * s * (pow_int(m2, 2) - 4 * m_pi_2) / pow_int(s - m_pi_2, 2) +
               (pow_int(s - m_pi_2, 2) -
                4 * pow_int(m2, 2) * (pow_int(m2, 2) - 4 * m_pi_2)) /
                   (u1 * u2));
          e += std::log(u1 / u2) *
               (5 * (s - m_pi_2) - pow_int(s - m_pi_2, 2) / pow_int(m2, 2) -
                4 * (pow_int(m2, 2) - 4 * m_pi_2) *
                    (s - m_pi_2 + pow_int(m2, 2)) / (s - m_pi_2));
          xsection = xsection * e * to_mb;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;

        case ReactionType::pi_eta:
          if (part_a.type().pdgcode() == pdg::pi_p) {
            part_out = pi_plus_particle;
          } else {
            part_out = pi_minus_particle;
          }
          m3 = part_out->mass();

          xsection = M_PI * alpha * 4.7 / (16 * m_eta_2 * s * p_cm_2);
          I0 = 1 / (m_rho * gamma_rho_tot) *
               (atan((u1 - m_rho_2) / (m_rho * gamma_rho_tot)) -
                atan((u2 - m_rho_2) / (m_rho * gamma_rho_tot)));
          I1 = std::log(
              (pow_int(u2 - m_rho_2, 2) + m_rho_2 * pow_int(gamma_rho_tot, 2)) /
              (pow_int(u1 - m_rho_2, 2) + m_rho_2 * pow_int(gamma_rho_tot, 2)));
          e = -m_pi_2 *
              ((t2 + u2) * (s - m_pi_2) + pow_int(2 * m_pi_2 - s, 2)) * I0;
          e += ((s - m_pi_2) * (m_pi_2 + t2 + u2) -
                2 * m_pi_2 * (s - 2 * m_pi_2)) *
               ((t2 + u2 - m_rho_2) * I0 + 0.5 * I1);
          e += -s * (t2 - t1 + (t2 + u2 - m_rho_2) * I1 +
                     pow_int(t2 + u2 - m_rho_2, 2) * I0 -
                     m_rho_2 * pow_int(gamma_rho_tot, 2) * I0);
          xsection = xsection * e * to_mb;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;

        case ReactionType::no_reaction:
          // never reached
          break;
      }
    }
  }
  return process_list;
}

double ScatterActionPhoton::pi_pi_rho0(const double M, const double s) const {
  const double to_mb = 0.3894;
  const double m_pi = ParticleType::find(pdg::pi_z).mass();
  const double m_pi_2 = pow_int(m_pi, 2);
  const double m_rho = ParticleType::find(pdg::rho_z).mass();
  const double gamma_rho_tot = ParticleType::find(pdg::rho_z).width_at_pole();
  const double g_rho_2 =
      24 * twopi * gamma_rho_tot * pow_int(m_rho, 2) /
      pow(pow_int(m_rho, 2) - 4 * pow_int(m_pi, 2), 3.0 / 2.0);
  const double DM = pow_int(M, 2) - 4 * pow_int(m_pi, 2);
  const double sqrts = sqrt(s);
  const double p_cm_2 = 0.25 * s - m_pi_2;
  if (sqrts <= M) {
    return 0;
  }
  std::array<double, 2> mandelstam_t = get_t_range(sqrts, m_pi, m_pi, M, 0.);
  double t1 = mandelstam_t[1];
  double t2 = mandelstam_t[0];
  double u1 = 2 * m_pi_2 + pow_int(M, 2) - s - t1;
  double u2 = 2 * m_pi_2 + pow_int(M, 2) - s - t2;
  double xsection = alpha * g_rho_2 / (4 * s * p_cm_2);

  t1 += -m_pi_2;
  t2 += -m_pi_2;
  if (std::abs(t1) < really_small) {
    t1 = -really_small;
  }
  if (t2 / t1 <= 0) {
    return 0;
  }
  u1 += -m_pi_2;
  u2 += -m_pi_2;
  if (std::abs(u2) < really_small) {
    return 0;
  }
  if (u1 / u2 <= 0) {
    return 0;
  }
  xsection = xsection *
             (2 * (t2 - t1) -
              DM * ((s - 2 * m_pi_2) / (s - pow_int(M, 2)) * std::log(t2 / t1) +
                    m_pi_2 * (t2 - t1) / (t2 * t1) +
                    (s - 2 * m_pi_2) / (s - pow_int(M, 2)) * std::log(u1 / u2) +
                    m_pi_2 * (u1 - u2) / (u1 * u2))) *
             to_mb;
  if (xsection > 0) {
    return xsection;
  } else {
    return really_small;
  }
}

double ScatterActionPhoton::pi_pi0_rho(const double M, const double s) const {
  const double to_mb = 0.3894;
  const double m_pi = ParticleType::find(pdg::pi_z).mass();
  const double m_pi_2 = pow_int(m_pi, 2);
  const double m_rho = ParticleType::find(pdg::rho_z).mass();
  const double gamma_rho_tot = ParticleType::find(pdg::rho_z).width_at_pole();
  const double g_rho_2 =
      24 * twopi * gamma_rho_tot * pow_int(m_rho, 2) /
      pow(pow_int(m_rho, 2) - 4 * pow_int(m_pi, 2), 3.0 / 2.0);
  const double DM = pow_int(M, 2) - 4 * pow_int(m_pi, 2);
  const double sqrts = sqrt(s);
  const double p_cm_2 = 0.25 * s - m_pi_2;
  if (sqrts <= M) {
    return 0;
  }
  std::array<double, 2> mandelstam_t = get_t_range(sqrts, m_pi, m_pi, M, 0.);
  double t1 = mandelstam_t[1];
  double t2 = mandelstam_t[0];
  double xsection = -alpha * g_rho_2 / (16 * s * p_cm_2);
  double e = 1.0 / 3.0 * (s - 2 * pow_int(M, 2)) / pow_int(M, 2) /
             pow_int(s - pow_int(M, 2), 2) * (pow_int(t2, 3) - pow_int(t1, 3));
  e += 0.5 * (s - 6 * pow_int(M, 2)) / pow_int(M, 2) / (s - pow_int(M, 2)) *
       (pow_int(t2, 2) - pow_int(t1, 2));
  t1 += -m_pi_2;
  t2 += -m_pi_2;
  if (std::abs(t1) < really_small) {
    t1 = -really_small;
  }
  if (t2 / t1 <= 0) {
    return 0;
  }
  e += (4 * s * DM / pow_int(s - pow_int(M, 2), 2) + m_pi_2 / pow_int(M, 2) -
        4.5) *
       (t2 - t1);
  e += 4 * s * DM / (s - pow_int(M, 2)) * std::log(t2 / t1);
  e += 4 * m_pi_2 * DM * ((t2 - t1) / (t2 * t1));
  xsection = xsection * e * to_mb;
  if (xsection > 0) {
    return xsection;
  } else {
    return really_small;
  }
}

double ScatterActionPhoton::diff_cross_section(double t, double m3) const {
  const double to_mb = 0.3894;
  const double m_rho = ParticleType::find(pdg::rho_z).mass();
  const double m_rho_2 = pow_int(m_rho, 2);
  const double m_pi = ParticleType::find(pdg::pi_z).mass();
  const double m_pi_2 = pow_int(m_pi, 2);
  const double m_eta = ParticleType::find(pdg::eta).mass();
  const double m_eta_2 = pow_int(m_eta, 2);
  const double gamma_rho_tot = ParticleType::find(pdg::rho_z).width_at_pole();
  const double g_rho_2 =
      24 * twopi * gamma_rho_tot * pow_int(m_rho, 2) /
      pow(pow_int(m_rho, 2) - 4 * pow_int(m_pi, 2), 3.0 / 2.0);
  double s = mandelstam_s();
  const double p_cm_2 = cm_momentum_squared();
  const double m1 = incoming_particles_[0].effective_mass();
  const double m2 = incoming_particles_[1].effective_mass();
  const double m3_2 = pow_int(m3, 2);
  const double DM = pow_int(m3, 2) - 4 * pow_int(m_pi, 2);
  double u = pow_int(m1, 2) + pow_int(m2, 2) + pow_int(m3, 2) - s - t;
  double diff_xsection = 0.0;
  double e = 0.0;
  switch (reac) {
    case ReactionType::pi_pi:
      if (outgoing_particles_[0].type().pdgcode().is_rho()) {
        diff_xsection = alpha * g_rho_2 / (4 * s * p_cm_2);
        diff_xsection =
            diff_xsection *
            (2 -
             DM / (t - m_pi_2) * ((s - 2 * m_pi_2) / (s - pow_int(m3, 2)) +
                                  m_pi_2 / (t - m_pi_2)) -
             DM / (u - m_pi_2) * ((s - 2 * m_pi_2) / (s - pow_int(m3, 2)) +
                                  m_pi_2 / (u - m_pi_2)));
      } else if (outgoing_particles_[0].type().pdgcode() == pdg::eta) {
        diff_xsection =
            twopi * alpha * 4.7 * pow_int(m_rho, 4) /
            (pow_int(s - m_rho_2, 2) + pow_int(gamma_rho_tot, 2) * m_rho_2) /
            (32 * m_eta_2 * pow_int(m_rho, 4) * s * p_cm_2);
        diff_xsection = diff_xsection * (s * (u - m_pi_2) * (t - m_pi_2) -
                                         m_pi_2 * pow_int(s - m_eta_2, 2));
      } else if (outgoing_particles_[0].type().pdgcode() == pdg::photon) {
        diff_xsection = twopi * pow_int(alpha, 2) / (s * p_cm_2);
        u += -m_pi_2;
        t += -m_pi_2;
        diff_xsection =
            diff_xsection * (1 + 2 * m_pi_2 * (1 / t + 1 / u) +
                             2 * pow_int(m_pi, 4) * pow_int(1 / t + 1 / u, 2));
      }
      break;
    case ReactionType::pi0_pi:
      diff_xsection = -alpha * g_rho_2 / (16 * s * p_cm_2);
      e = (s - 2 * m3_2) * pow_int(t - m_pi_2, 2) / m3_2 /
              pow_int(s - m3_2, 2) +
          (s - 2 * m3_2) * pow_int(u - m_pi_2, 2) / m3_2 / pow_int(s - m3_2, 2);
      e += (s - 6 * m3_2) * (t - m_pi_2) / m3_2 / (s - m3_2) +
           (s - 6 * m3_2) * (u - m_pi_2) / m3_2 / (s - m3_2);
      e += 2 * 4 * s * DM / pow_int(s - m3_2, 2);
      e += 4 * DM / (t - m_pi_2) * (s / (s - m3_2) + m_pi_2 / (t - m_pi_2)) +
           4 * DM / (u - m_pi_2) * (s / (s - m3_2) + m_pi_2 / (u - m_pi_2));
      e += 2 * (m_pi_2 / m3_2 - 4.5);
      diff_xsection = diff_xsection * e;
      break;
    case ReactionType::pi_rho0:
      diff_xsection = alpha * g_rho_2 / (12 * s * p_cm_2);
      diff_xsection = diff_xsection *
                      (2 -
                       s * (pow_int(m2, 2) - 4 * pow_int(m_pi, 2)) /
                           pow_int(s - m_pi_2, 2) -
                       (pow_int(m2, 2) - 4 * pow_int(m_pi, 2)) / (t - m_pi_2) *
                           ((s - pow_int(m2, 2) + m_pi_2) / (s - m_pi_2) +
                            m_pi_2 / (t - m_pi_2)));
      break;
    case ReactionType::pi_rho:
      diff_xsection = -alpha * g_rho_2 / (48 * s * p_cm_2);
      e = 4 * (pow_int(m2, 2) - 4 * pow_int(m_pi, 2)) *
          (t / pow_int(t - m_pi_2, 2) + u / pow_int(u - m2, 2) -
           m2 / (s - m_pi_2) * (1 / (t - m_pi_2) + 1 / (u - m2)));
      e += (3.0 + (s - m_pi_2) / m2) * (s - m_pi_2) / (u - m2) - 0.5 + s / m2 -
           pow_int((s - m_pi_2) / (u - m2), 2);
      diff_xsection = diff_xsection * e;
      break;
    case ReactionType::pi0_rho:
      diff_xsection = alpha * g_rho_2 / (48 * s * p_cm_2);
      e = 4.5 - s / m2 -
          4 * s * (pow_int(m2, 2) - 4 * pow_int(m_pi, 2)) /
              pow_int(s - m_pi_2, 2) +
          (pow_int(s - m_pi_2, 2) -
           4 * m2 * (pow_int(m2, 2) - 4 * pow_int(m_pi, 2))) /
              pow_int(u - m2, 2);
      e += 1 / (u - m2) * (5 * (s - m_pi_2) - pow_int(s - m_pi_2, 2) / m2 -
                           4 * (pow_int(m2, 2) - 4 * pow_int(m_pi, 2)) *
                               (s - m_pi_2 + m2) / (s - m_pi_2));
      diff_xsection = diff_xsection * e;
      break;
    case ReactionType::pi_eta:
      diff_xsection =
          twopi * alpha * 4.7 *
          (pow_int(m_rho, 4) /
           (pow_int(u - m_rho_2, 2) + pow_int(gamma_rho_tot, 2) * m_rho_2)) /
          (32 * m_eta_2 * pow_int(m_rho, 4) * s * p_cm_2);
      diff_xsection = diff_xsection * (u * (s - m_pi_2) * (t - m_pi_2) -
                                       m_pi_2 * pow_int(u - m_eta_2, 2));
      break;
    case ReactionType::no_reaction:
      // never reached
      break;
  }
  return diff_xsection * to_mb;
}

}  // namespace Smash
