/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionphoton.h"

#include "include/angles.h"
#include "include/cxx14compat.h"
#include "include/integrate.h"
#include "include/kinematics.h"
#include "include/particletype.h"
#include "include/pdgcode.h"
#include "include/random.h"
#include "include/tabulation.h"

#include <fstream>
#include <iostream>

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

  const double pcm = cm_momentum();
  double t = Random::uniform(t1, t2);
  double costheta =
      (t - pow(m1, 2) +
       0.5 * (s + pow(m1, 2) - pow(m2, 2)) * (s - pow(m3, 2)) / s) /
      (pcm * (s - pow(m3, 2)) / sqrts);
  if (costheta > 1)
    costheta = 1;
  if (costheta < -1)
    costheta = -1;
  Angles phitheta(Random::uniform(0.0, twopi), costheta);
  outgoing_particles_[0].set_4momentum(masses.first, phitheta.threevec() * pcm);
  outgoing_particles_[1].set_4momentum(masses.second,
                                       -phitheta.threevec() * pcm);

  /* Weighing of the fractional photons */
  weight_ = diff_cross_section(t) * (t2 - t1) / number_of_fractional_photons /
            total_cross_section_;

  /* Set positions & boost to computational frame. */
  for (ParticleData &new_particle : outgoing_particles_) {
    new_particle.set_4position(middle_point);
    new_particle.boost_momentum(-beta_cm());
  }
}

CollisionBranchList ScatterActionPhoton::two_to_two_cross_sections() {
  CollisionBranchList process_list;
  const float m_rho = ParticleType::find(0x113).mass();
  const float m_rho_2 = pow(m_rho, 2);
  const float m_pi = ParticleType::find(0x111).mass();
  const float m_pi_2 = pow(m_pi, 2);
  const float m_eta = ParticleType::find(0x221).mass();
  const float m_eta_2 = pow(m_eta, 2);
  const float gamma_rho_tot = ParticleType::find(0x113).width_at_pole();
  const float g_rho_2 = 24 * twopi * gamma_rho_tot * pow(m_rho, 2) /
                        pow(pow(m_rho, 2) - 4 * pow(m_pi, 2), 3.0 / 2.0);
  const float to_mb = 0.3894;

  ParticleData &part_a = incoming_particles_[0];
  ParticleData &part_b = incoming_particles_[1];

  bool no_pion = false;

  if (!part_a.type().pdgcode().is_pion()) {
    if (part_b.type().pdgcode().is_pion()) {
      ParticleData dummy = part_a;
      part_a = part_b;
      part_b = dummy;
    } else
      no_pion = true;
  }

  if (no_pion) {
  } else {  // do a check according to incoming_particles_ and calculate the
            // cross sections (xsection) for all possible reactions

    const double s = mandelstam_s();
    double sqrts = sqrt_s();
    const double &m1 = part_a.effective_mass();
    const double &m2 = part_b.effective_mass();
    double m3 = 0.0;  // will be fixed according to reaction outcome
    const double p_cm_2 = cm_momentum_squared();
    ParticleTypePtr part_out = &ParticleType::find(0x022);
    ParticleTypePtr photon_out = &ParticleType::find(0x022);

    reac = no_reaction;
    if (part_a.type().charge() == 0) {
      if (part_b.type().charge() != 0) {
        if (part_b.type().pdgcode().is_pion()) {
          reac = pi0_pi;
          part_out = &ParticleType::find(0x213);
        }
        if (part_b.type().pdgcode().is_rho()) {
          reac = pi0_rho;
          part_out = &ParticleType::find(0x211);
        }
      }
    } else {  // so part_a.type().charge()!=0
      if (part_b.type().charge() == 0) {
        if (part_b.type().pdgcode().is_pion()) {
          reac = pi0_pi;
          part_out = &ParticleType::find(0x213);
        }
        if (part_b.type().pdgcode().is_rho()) {
          reac = piplus_rho0;
          part_out = &ParticleType::find(0x211);
        }
        if (part_b.type().pdgcode() == 0x221) {  // corresponds to eta meson
          reac = piplus_eta;
          part_out = &ParticleType::find(0x211);
        }
      } else if (part_b.type().charge() == -part_a.type().charge()) {
        if (part_b.type().pdgcode().is_pion()) {
          reac = pi_pi;  // actually three reactions (afh), start with h
        }
        if (part_b.type().pdgcode().is_rho()) {
          reac = pi_rho;
          part_out = &ParticleType::find(0x111);
        }
      }
    }

    m3 = part_out->mass();
    if (sqrts <= m1 + m2)
      reac = no_reaction;
    if (reac != no_reaction) {
      std::array<double, 2> mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
      double t1 = mandelstam_t[1];
      double t2 = mandelstam_t[0];

      double u1 = pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s - t1;
      double u2 = pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s - t2;

      double e, I0, I1;
      float xsection = 0.0;

      Integrator1dMonte integrate;

      switch (reac) {
        case pi_pi:
          xsection = twopi * pow(alpha, 2) / (s * p_cm_2);
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
          part_out = &ParticleType::find(0x221);
          m3 = part_out->mass();
          if (sqrts > m3) {
            mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
            t1 = mandelstam_t[1];
            t2 = mandelstam_t[0];
            u1 = pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s - t1;
            u2 = pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s - t2;
            xsection = M_PI * alpha * 4.7 * pow(m_rho, 4) /
                       (pow(s - m_rho_2, 2) + pow(gamma_rho_tot, 2) * m_rho_2) /
                       (16 * m_eta_2 * pow(m_rho, 4) * s * p_cm_2);
            xsection = xsection *
                       ((2 * m_pi_2 + m_eta_2 - s) * s / 2 *
                            (pow(t2, 2) - pow(t1, 2)) -
                        s / 3 * (pow(t2, 3) - pow(t1, 3)) -
                        (t2 - t1) * m_pi_2 *
                            (pow(m_eta, 4) + s * (m_pi_2 - m_eta_2))) *
                       to_mb;
            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          }

          // and the third possible reaction (produces rho0)

          part_out = &ParticleType::find(0x113);
          m3 = part_out->mass();
          if (tabulation_pi_pi_rho0 == nullptr) {
            tabulation_pi_pi_rho0 = make_unique<Tabulation>(
                2.0f * m_pi, 15.0f - 2.0f * m_pi, num_tab_pts,
                [&](float sqrts1) {
                  return integrate(2.0f * m_pi, sqrts1, [&](float M) {
                    return pi_pi_rho0(M, pow(sqrts1, 2)) *
                           part_out->spectral_function(M);
                  });
                });
          }
          xsection = tabulation_pi_pi_rho0->get_value_linear(sqrts);

          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;
        case pi0_pi:
          if (tabulation_pi0_pi_rho == nullptr) {
            tabulation_pi0_pi_rho = make_unique<Tabulation>(
                2.0f * m_pi, 15.0f - 2.0f * m_pi, num_tab_pts,
                [&](float sqrts1) {
                  return integrate(2.0f * m_pi, sqrts1, [&](float M) {
                    return pi_pi0_rho(M, pow(sqrts1, 2)) *
                           part_out->spectral_function(M);
                  });
                });
          }
          xsection = tabulation_pi0_pi_rho->get_value_linear(sqrts);
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;
        case piplus_rho0:
          xsection = alpha * g_rho_2 / (12 * s * p_cm_2);
          t1 += -m_pi_2;
          t2 += -m_pi_2;
          xsection = xsection * (2 * (t2 - t1) -
                                 s * (pow(m2, 2) - 4 * m_pi_2) /
                                     pow(s - m_pi_2, 2) * (t2 - t1) -
                                 (pow(m2, 2) - 4 * m_pi_2) *
                                     ((s - pow(m2, 2) + m_pi_2) / (s - m_pi_2) *
                                          std::log(t2 / t1) +
                                      m_pi_2 * (t2 - t1) / (t1 * t2))) *
                     to_mb;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;
        case pi_rho:
          xsection = -alpha * g_rho_2 / (48 * s * p_cm_2);
          t1 += -m_pi_2;
          t2 += -m_pi_2;
          u1 += -m_rho_2;
          u2 += -m_rho_2;
          e = 4 * (pow(m2, 2) - 4 * m_pi_2) *
              (pow(m2, 2) * (t2 - t1) / (u1 * u2) +
               m_pi_2 * (t2 - t1) / (t2 * t1) + std::log(u1 / u2 * t2 / t1) -
               pow(m2, 2) / (s - m_pi_2) * std::log(t2 / t1 * u1 / u2));
          e += (s - m_pi_2) * (3.0 + (s - m_pi_2) / pow(m2, 2)) *
               std::log(u1 / u2);
          e += (t2 - t1) *
               (s / pow(m2, 2) - 0.5 - pow(s - m_pi_2, 2) / (u1 * u2));
          xsection = xsection * e * to_mb;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;
        case pi0_rho:
          xsection = alpha * g_rho_2 / (48 * s * p_cm_2);
          u1 += -pow(m2, 2);  // is u+
          u2 += -pow(m2, 2);
          e = (t2 - t1) *
              (4.5 - s / pow(m2, 2) -
               4 * s * (pow(m2, 2) - 4 * m_pi_2) / pow(s - m_pi_2, 2) +
               (pow(s - m_pi_2, 2) -
                4 * pow(m2, 2) * (pow(m2, 2) - 4 * m_pi_2)) /
                   (u1 * u2));
          e += std::log(u1 / u2) *
               (5 * (s - m_pi_2) - pow(s - m_pi_2, 2) / pow(m2, 2) -
                4 * (pow(m2, 2) - 4 * m_pi_2) * (s - m_pi_2 + pow(m2, 2)) /
                    (s - m_pi_2));
          xsection = xsection * e * to_mb;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;
        case piplus_eta:
          xsection = M_PI * alpha * 4.7 / (16 * m_eta_2 * s * p_cm_2);
          I0 = 1 / (m_rho * gamma_rho_tot) *
               (atan((u1 - m_rho_2) / (m_rho * gamma_rho_tot)) -
                atan((u2 - m_rho_2) / (m_rho * gamma_rho_tot)));
          I1 = std::log(
              (pow(u2 - m_rho_2, 2) + m_rho_2 * pow(gamma_rho_tot, 2)) /
              (pow(u1 - m_rho_2, 2) + m_rho_2 * pow(gamma_rho_tot, 2)));
          e = -m_pi_2 * ((t2 + u2) * (s - m_pi_2) + pow(2 * m_pi_2 - s, 2)) *
              I0;
          e += ((s - m_pi_2) * (m_pi_2 + t2 + u2) -
                2 * m_pi_2 * (s - 2 * m_pi_2)) *
               ((t2 + u2 - m_rho_2) * I0 + 0.5 * I1);
          e += -s * (t2 - t1 + (t2 + u2 - m_rho_2) * I1 +
                     pow(t2 + u2 - m_rho_2, 2) * I0 -
                     m_rho_2 * pow(gamma_rho_tot, 2) * I0);
          xsection = xsection * e * to_mb;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;
        case no_reaction:
          // never reached
          break;
      }
    }

    // add to extra CollisionBranch only for photon producing reactions!
    add_processes<CollisionBranch>(std::move(process_list),
                                   collision_channels_photons_,
                                   cross_section_photons_);
  }

  return process_list;
}

float ScatterActionPhoton::pi_pi_rho0(const float M, const float s) const {
  const float to_mb = 0.3894;
  const float m_pi = ParticleType::find(0x111).mass();
  const float m_pi_2 = pow(m_pi, 2);
  const float m_rho = ParticleType::find(0x113).mass();
  const float gamma_rho_tot = ParticleType::find(0x113).width_at_pole();
  const float g_rho_2 = 24 * twopi * gamma_rho_tot * pow(m_rho, 2) /
                        pow(pow(m_rho, 2) - 4 * pow(m_pi, 2), 3.0 / 2.0);
  const float DM = pow(M, 2) - 4 * pow(m_pi, 2);
  const float sqrts = sqrt(s);
  const float p_cm_2 = 0.25 * s - m_pi_2;
  if (sqrts <= M) {
    return 0;
  }
  std::array<float, 2> mandelstam_t = get_t_range(sqrts, m_pi, m_pi, M, 0.0f);
  float t1 = mandelstam_t[1];
  float t2 = mandelstam_t[0];
  float u1 = 2 * m_pi_2 + pow(M, 2) - s - t1;
  float u2 = 2 * m_pi_2 + pow(M, 2) - s - t2;
  float xsection = alpha * g_rho_2 / (4 * s * p_cm_2);

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
              DM * ((s - 2 * m_pi_2) / (s - pow(M, 2)) * std::log(t2 / t1) +
                    m_pi_2 * (t2 - t1) / (t2 * t1) +
                    (s - 2 * m_pi_2) / (s - pow(M, 2)) * std::log(u1 / u2) +
                    m_pi_2 * (u1 - u2) / (u1 * u2))) *
             to_mb;
  if (xsection > 0) {
    return xsection;
  } else {
    return really_small;
  }
}

float ScatterActionPhoton::pi_pi0_rho(const float M, const float s) const {
  const float to_mb = 0.3894;
  const float m_pi = ParticleType::find(0x111).mass();
  const float m_pi_2 = pow(m_pi, 2);
  const float m_rho = ParticleType::find(0x113).mass();
  const float gamma_rho_tot = ParticleType::find(0x113).width_at_pole();
  const float g_rho_2 = 24 * twopi * gamma_rho_tot * pow(m_rho, 2) /
                        pow(pow(m_rho, 2) - 4 * pow(m_pi, 2), 3.0 / 2.0);
  const float DM = pow(M, 2) - 4 * pow(m_pi, 2);
  const float sqrts = sqrt(s);
  const float p_cm_2 = 0.25 * s - m_pi_2;
  if (sqrts <= M) {
    return 0;
  }
  std::array<float, 2> mandelstam_t = get_t_range(sqrts, m_pi, m_pi, M, 0.0f);
  float t1 = mandelstam_t[1];
  float t2 = mandelstam_t[0];
  float xsection = -alpha * g_rho_2 / (16 * s * p_cm_2);
  float e = 1.0 / 3.0 * (s - 2 * pow(M, 2)) / pow(M, 2) /
            pow(s - pow(M, 2), 2) * (pow(t2, 3) - pow(t1, 3));
  e += 0.5 * (s - 6 * pow(M, 2)) / pow(M, 2) / (s - pow(M, 2)) *
       (pow(t2, 2) - pow(t1, 2));
  t1 += -m_pi_2;
  t2 += -m_pi_2;
  if (std::abs(t1) < really_small) {
    t1 = -really_small;
  }
  if (t2 / t1 <= 0) {
    return 0;
  }
  e += (4 * s * DM / pow(s - pow(M, 2), 2) + m_pi_2 / pow(M, 2) - 4.5) *
       (t2 - t1);
  e += 4 * s * DM / (s - pow(M, 2)) * std::log(t2 / t1);
  e += 4 * m_pi_2 * DM * ((t2 - t1) / (t2 * t1));
  xsection = xsection * e * to_mb;
  if (xsection > 0) {
    return xsection;
  } else {
    return really_small;
  }
}

float ScatterActionPhoton::diff_cross_section(float t) const {
  const float m_rho = ParticleType::find(0x113).mass();
  const float m_rho_2 = pow(m_rho, 2);
  const float m_pi = ParticleType::find(0x111).mass();
  const float m_pi_2 = pow(m_pi, 2);
  const float m_eta = ParticleType::find(0x221).mass();
  const float m_eta_2 = pow(m_eta, 2);
  const float gamma_rho_tot = ParticleType::find(0x113).width_at_pole();
  const float g_rho_2 = 24 * twopi * gamma_rho_tot * pow(m_rho, 2) /
                        pow(pow(m_rho, 2) - 4 * pow(m_pi, 2), 3.0 / 2.0);
  const float s = mandelstam_s();
  const float p_cm_2 = cm_momentum_squared();
  const float m1 = incoming_particles_[0].effective_mass();
  const float m2 = incoming_particles_[1].effective_mass();
  const float m3 = outgoing_particles_[0].effective_mass();
  const float m3_2 = pow(m3, 2);
  const float DM = pow(m3, 2) - 4 * pow(m_pi, 2);
  float u = pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s - t;
  float diff_xsection = 0.0;
  float e = 0.0;
  switch (reac) {
    case pi_pi:
      if (outgoing_particles_[0].type().pdgcode().is_rho()) {
        diff_xsection = alpha * g_rho_2 / (4 * s * p_cm_2);
        diff_xsection =
            diff_xsection *
            (2 -
             DM / (t - m_pi_2) *
                 ((s - 2 * m_pi_2) / (s - pow(m3, 2)) + m_pi_2 / (t - m_pi_2)) -
             DM / (u - m_pi_2) *
                 ((s - 2 * m_pi_2) / (s - pow(m3, 2)) + m_pi_2 / (u - m_pi_2)));
      }
      if (outgoing_particles_[0].type().pdgcode() == 0x221) {
        diff_xsection =
            twopi * alpha * 4.7 * pow(m_rho, 4) /
            (pow(s - m_rho_2, 2) + pow(gamma_rho_tot, 2) * m_rho_2) /
            (32 * m_eta_2 * pow(m_rho, 4) * s * p_cm_2);
        diff_xsection = diff_xsection * (s * (u - m_pi_2) * (t - m_pi_2) -
                                         m_pi_2 * pow(s - m_eta_2, 2));
      }
      if (outgoing_particles_[0].type().pdgcode() == 0x022) {
        diff_xsection = twopi * pow(alpha, 2) / (s * p_cm_2);
        u += -m_pi_2;
        t += -m_pi_2;
        diff_xsection =
            diff_xsection * (1 + 2 * m_pi_2 * (1 / t + 1 / u) +
                             2 * pow(m_pi, 4) * pow(1 / t + 1 / u, 2));
      }
      break;
    case pi0_pi:
      diff_xsection = -alpha * g_rho_2 / (16 * s * p_cm_2);
      e = (s - 2 * m3_2) * pow(t - m_pi_2, 2) / m3_2 / pow(s - m3_2, 2) +
          (s - 2 * m3_2) * pow(u - m_pi_2, 2) / m3_2 / pow(s - m3_2, 2);
      e += (s - 6 * m3_2) * (t - m_pi_2) / m3_2 / (s - m3_2) +
           (s - 6 * m3_2) * (u - m_pi_2) / m3_2 / (s - m3_2);
      e += 2 * 4 * s * DM / pow(s - m3_2, 2);
      e += 4 * DM / (t - m_pi_2) * (s / (s - m3_2) + m_pi_2 / (t - m_pi_2)) +
           4 * DM / (u - m_pi_2) * (s / (s - m3_2) + m_pi_2 / (u - m_pi_2));
      e += 2 * (m_pi_2 / m3_2 - 4.5);
      diff_xsection = diff_xsection * e;
      break;
    case piplus_rho0:
      diff_xsection = alpha * g_rho_2 / (12 * s * p_cm_2);
      diff_xsection =
          diff_xsection *
          (2 - s * (pow(m2, 2) - 4 * pow(m_pi, 2)) / pow(s - m_pi_2, 2) -
           (pow(m2, 2) - 4 * pow(m_pi, 2)) / (t - m_pi_2) *
               ((s - pow(m2, 2) + m_pi_2) / (s - m_pi_2) +
                m_pi_2 / (t - m_pi_2)));
      break;
    case pi_rho:
      diff_xsection = -alpha * g_rho_2 / (48 * s * p_cm_2);
      e = 4 * (pow(m2, 2) - 4 * pow(m_pi, 2)) *
          (t / pow(t - m_pi_2, 2) + u / pow(u - m2, 2) -
           m2 / (s - m_pi_2) * (1 / (t - m_pi_2) + 1 / (u - m2)));
      e += (3.0 + (s - m_pi_2) / m2) * (s - m_pi_2) / (u - m2) - 0.5 + s / m2 -
           pow((s - m_pi_2) / (u - m2), 2);
      diff_xsection = diff_xsection * e;
      break;
    case pi0_rho:
      diff_xsection = alpha * g_rho_2 / (48 * s * p_cm_2);
      e = 4.5 - s / m2 -
          4 * s * (pow(m2, 2) - 4 * pow(m_pi, 2)) / pow(s - m_pi_2, 2) +
          (pow(s - m_pi_2, 2) - 4 * m2 * (pow(m2, 2) - 4 * pow(m_pi, 2))) /
              pow(u - m2, 2);
      e += 1 / (u - m2) * (5 * (s - m_pi_2) - pow(s - m_pi_2, 2) / m2 -
                           4 * (pow(m2, 2) - 4 * pow(m_pi, 2)) *
                               (s - m_pi_2 + m2) / (s - m_pi_2));
      diff_xsection = diff_xsection * e;
      break;
    case piplus_eta:
      diff_xsection = twopi * alpha * 4.7 *
                      (pow(m_rho, 4) / (pow(u - m_rho_2, 2) +
                                        pow(gamma_rho_tot, 2) * m_rho_2)) /
                      (32 * m_eta_2 * pow(m_rho, 4) * s * p_cm_2);
      diff_xsection = diff_xsection * (u * (s - m_pi_2) * (t - m_pi_2) -
                                       m_pi_2 * pow(u - m_eta_2, 2));
      break;
    case no_reaction:
      // never reached
      break;
  }
  return diff_xsection;
}

}  // namespace Smash
