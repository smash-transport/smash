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
  const double stepsize = (t2-t1)/100.0;
  for (double t = t1; t < t2; t += stepsize) {
    float diff_xsection_max = std::max(diff_cross_section(t, m3,t2,t1),
                                       diff_xsection_max);
  }

  float t = Random::uniform(t1, t2);
  float diff_xsection_max = 0;
  int iteration_number = 0;
  do {
    t = Random::uniform(t1, t2);
    iteration_number++;
  } while (diff_cross_section(t, m3,t2,t1) < Random::uniform(0.f, diff_xsection_max)
           && iteration_number < 100);

  // todo: this should move to kinematics.h and tested
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
    weight_ = diff_cross_section(t, m3,t2,t1) * (t2 - t1)
          / (number_of_fractional_photons_ * cross_section());
  } else {
    weight_ = proc->weight() / cross_section();
  }
  /* Set positions & boost to computational frame. */
  for (ParticleData &new_particle : outgoing_particles_) {
    new_particle.set_4position(middle_point);
    new_particle.boost_momentum(-beta_cm());
  }
}

void ScatterActionPhoton::add_dummy_hadronic_channels(
                            float reaction_cross_section) {
  CollisionBranchPtr dummy_process = make_unique<CollisionBranch>(
    incoming_particles_[0].type(),
    incoming_particles_[1].type(),
    reaction_cross_section,
    ProcessType::TwoToTwo);
  add_collision(std::move(dummy_process));
}

ScatterActionPhoton::ReactionType
  ScatterActionPhoton::is_photon_reaction(const ParticleList &in) {
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
    case(pack(pdg::pi_p, pdg::pi_z)):
    case(pack(pdg::pi_z, pdg::pi_p)):
    case(pack(pdg::pi_m, pdg::pi_z)):
    case(pack(pdg::pi_z, pdg::pi_m)):
      return ReactionType::pi0_pi;
    case(pack(pdg::pi_p, pdg::rho_z)):
    case(pack(pdg::pi_m, pdg::rho_z)):
      return ReactionType::pi_rho0;
    case(pack(pdg::pi_m, pdg::rho_p)):
    case(pack(pdg::pi_p, pdg::rho_m)):
      return ReactionType::pi_rho;
    case(pack(pdg::pi_z, pdg::rho_p)):
    case(pack(pdg::pi_z, pdg::rho_m)):
      return ReactionType::pi0_rho;
    /*case(pack(pdg::pi_p, pdg::eta)):
    case(pack(pdg::pi_m, pdg::eta)):
      return ReactionType::pi_eta;*/
    case(pack(pdg::pi_p, pdg::pi_m)):
    case(pack(pdg::pi_m, pdg::pi_p)):
      return ReactionType::pi_pi;
    case(pack(pdg::pi_z, pdg::rho_z)):
      return ReactionType::pi0_rho0;
    default:
      return ReactionType::no_reaction;
  }
}

CollisionBranchList ScatterActionPhoton::photon_cross_sections() {
  CollisionBranchList process_list;
  ParticleTypePtr rho0_particle = &ParticleType::find(pdg::rho_z);
  ParticleTypePtr rho_plus_particle = &ParticleType::find(pdg::rho_p);
  ParticleTypePtr rho_minus_particle = &ParticleType::find(pdg::rho_m);
  ParticleTypePtr pi0_particle = &ParticleType::find(pdg::pi_z);
  ParticleTypePtr pi_plus_particle = &ParticleType::find(pdg::pi_p);
  ParticleTypePtr pi_minus_particle = &ParticleType::find(pdg::pi_m);
  ParticleTypePtr photon_particle = &ParticleType::find(pdg::photon);
  const float m_rho = rho0_particle->mass();
  const float m_pi = pi0_particle->mass();

  const float to_mb = 0.3894;
  const float Const = 0.059;
  const float g_POR = 11.93;
  const float ma1 = 1.26;
  const float ghat = 6.4483;
  const float eta1 = 2.3920;
  const float eta2 = 1.9430;
  const float delta = -0.6426;
  const float C4 = -0.14095;
  const float Gammaa1 = 0.4;
  const float Pi = M_PI;
  float m_omega = 0.783;
  float momega = m_omega;
  float mrho = m_rho;
  float mpion = m_pi;

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
      std::array<double, 2> mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
      double t1;
      double t2;
      float xsection = 0.0;

      switch (reac) {
         case ReactionType::pi_pi:
        // there are three possible reaction channels
        // the first possible reaction produces eta
        /*  part_out = eta_particle;
          m3 = part_out->mass();

          if (sqrts > m3) {
            mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
            t1 = mandelstam_t[1];
            t2 = mandelstam_t[0];

            xsection = to_be_determined * to_mb;
            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          }*/

          // the second possible reaction (produces rho0)
          part_out = rho0_particle;
          m3 = part_out->mass();

          if (sqrts > m3) {
            mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
            t1 = mandelstam_t[1];
            t2 = mandelstam_t[0];

            xsection = to_mb*((pow(Const,2)*pow(ghat,4)*((0.5*pow(-2. + delta,2)*pow(mpion,2)*(0. - 1.*t2))/pow(mrho,2) +
             0.25*(-2. + delta)*(eta1 - 1.*eta2)*(-0.5*eta2*pow(ma1,2) + 1.*eta1*pow(mpion,2) - 1.*eta2*pow(mpion,2) + 0.5*eta1*s -
                0.5*eta2*s)*(0. - 1.*t2) - (2.*(pow(mpion,2)*
                   (1.5 + 0.125*pow(delta,2) - 2.*C4*pow(mrho,2) + delta*(-1. + 1.*C4*pow(mrho,2))) +
                  (-0.5 - 0.375*pow(delta,2) - 2.*C4*pow(mrho,2) + delta*(1. + 1.*C4*pow(mrho,2)))*s)*(0. - 1.*t2))/pow(mrho,2) -
             0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(1.*pow(mpion,2) - 0.5*s) +
                eta2*(-0.5*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 1.*s))*(0. - 1.*t2) -
             (0.25*(-2. + 1.*delta)*(-8.*C4*pow(mrho,4) + pow(mpion,2)*(2. + 1.*delta - 8.*C4*pow(mrho,2)) + (-2. - 3.*delta)*s +
                  pow(mrho,2)*(2. + 1.*delta + 16.*C4*s))*(0. - 1.*t2))/pow(mrho,2) -
             0.09375*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-2.*pow(ma1,4) - 4.*pow(mpion,4) + 0.6666666666666666*pow(mrho,4) +
                   pow(ma1,2)*(5.333333333333333*pow(mpion,2) - 2.6666666666666665*s) + 2.6666666666666665*pow(mpion,2)*s +
                   1.3333333333333333*pow(mrho,2)*s - 1.3333333333333333*pow(s,2)) +
                pow(eta1,2)*(1.*pow(ma1,4) + 2.*pow(mpion,4) + 0.3333333333333333*pow(mrho,4) +
                   pow(mpion,2)*(2.*pow(mrho,2) - 1.3333333333333333*s) - 1.3333333333333333*pow(mrho,2)*s + 0.6666666666666666*pow(s,2) +
                   pow(ma1,2)*(-2.6666666666666665*pow(mpion,2) - 1.3333333333333333*pow(mrho,2) + 1.3333333333333333*s)) +
                pow(eta2,2)*(1.*pow(ma1,4) + 2.*pow(mpion,4) + 0.3333333333333333*pow(mrho,4) +
                   pow(mpion,2)*(-2.*pow(mrho,2) - 1.3333333333333333*s) - 0.6666666666666666*pow(mrho,2)*s + 0.6666666666666666*pow(s,2) +
                   pow(ma1,2)*(-2.6666666666666665*pow(mpion,2) + 1.3333333333333333*pow(mrho,2) + 1.3333333333333333*s)))*(0. - 1.*t2) -
             0.0625*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(pow(Gammaa1,2)*pow(ma1,2) - 1.*pow(ma1,4) - 2.*pow(mpion,4) -
                   2.*pow(mpion,2)*pow(mrho,2) + 2.*pow(mrho,4) + pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) -
                   3.*pow(mrho,2)*s + pow(s,2)) + pow(eta1,2)*
                 (pow(Gammaa1,2)*pow(ma1,2) - 1.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) +
                   pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) + pow(mrho,2)*s + pow(s,2)) +
                eta1*eta2*(-2.*pow(Gammaa1,2)*pow(ma1,2) + 2.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) +
                   2.*pow(mrho,4) - 2.*pow(s,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)))*(0. - 1.*t2) +
             (1.*(pow(eta1,2)*(0.5*pow(mrho,4) - 1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                     pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) - 0.5*pow(mrho,2)*s + 0.25*delta*pow(mrho,2)*s -
                     0.25*delta*pow(s,2) + 1.*C4*pow(mrho,2)*pow(s,2)) +
                  pow(eta2,2)*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                     pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) - 0.5*pow(mrho,2)*s + 0.75*delta*pow(mrho,2)*s -
                     2.*C4*pow(mrho,4)*s - 0.25*delta*pow(s,2) + 1.*C4*pow(mrho,2)*pow(s,2)) +
                  eta1*eta2*(pow(ma1,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) + pow(mpion,2)*(-2.*pow(mrho,2) + 4.*C4*pow(mrho,4)) +
                     s*(2.*C4*pow(mrho,4) + 0.5*delta*s + pow(mrho,2)*(1. - 1.*delta - 2.*C4*s))))*(0. - 1.*t2))/pow(mrho,2) +
             (0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
                     pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) + pow(ma1,2)*pow(mpion,2)*
                      (8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
                     pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
                  pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
                     pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                        2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s +
                        pow(mpion,4)*(6.*pow(mrho,2) + 2.*s))) +
                  pow(eta1,2)*(1.*pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                        2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
                        pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
                     pow(mpion,2)*(1.*pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                        pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/(0. + 1.*pow(ma1,2) - 1.*t2) +
             (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(0. + 1.*pow(mpion,2) - 1.*t2) +
             0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(1.*pow(Gammaa1,2)*pow(ma1,2) - 3.*pow(ma1,4) - 2.*pow(mpion,4) +
                   2.*pow(mpion,2)*pow(mrho,2) - 1.*pow(mrho,4) + pow(ma1,2)*(4.*pow(mpion,2) - 4.*pow(mrho,2) - 4.*s) +
                   2.*pow(mrho,2)*s - 2.*pow(s,2)) + pow(eta1,2)*
                 (1.*pow(Gammaa1,2)*pow(ma1,2) - 3.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) - 1.*pow(mrho,4) +
                   pow(ma1,2)*(4.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) + 4.*pow(mrho,2)*s - 2.*pow(s,2)) +
                eta1*eta2*(-2.*pow(Gammaa1,2)*pow(ma1,2) + 6.*pow(ma1,4) + 4.*pow(mpion,4) - 2.*pow(mrho,4) - 4.*pow(mrho,2)*s +
                   4.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) + 8.*s)))*(0. + 1.*pow(mrho,2) - 1.*s - 1.*t2) -
             0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) - 1.*pow(mrho,2) - 1.*s) + eta1*(pow(ma1,2) - 1.*pow(mrho,2) + s))*
              pow(0. + 1.*pow(mrho,2) - 1.*s - 1.*t2,2) -
             0.010416666666666666*pow(eta1 - 1.*eta2,4)*pow(0. + 1.*pow(mrho,2) - 1.*s - 1.*t2,3) +
             (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/
              (0. + 1.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s - 1.*t2) -
             (1.*(1.*eta1 - 1.*eta2)*(eta1*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6) + pow(ma1,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                     pow(mpion,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) + 0.25*delta*pow(s,2) - 1.*C4*pow(mrho,2)*pow(s,2)) +
                  eta2*(-0.75*pow(mrho,4) + 0.125*delta*pow(mrho,4) + 1.*C4*pow(mrho,6) +
                     pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) + pow(ma1,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) +
                     0.25*pow(mrho,2)*s + 0.375*delta*pow(mrho,2)*s - 2.*C4*pow(mrho,4)*s - 0.25*delta*pow(s,2) +
                     1.*C4*pow(mrho,2)*pow(s,2)))*(0. - 1.*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s - 1.*t2))/pow(mrho,2)\
              + 0.5*pow(1.*eta1 - 1.*eta2,2)*(-0.5 + 1.*C4*pow(mrho,2))*
              pow(0. - 1.*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s - 1.*t2,2) +
             (0.25*(32.*pow(C4,2)*pow(mrho,8) + 2.*pow(delta,2)*pow(s,2) + 8.*C4*pow(mrho,6)*(-6. + delta - 8.*C4*s) +
                  2.*delta*pow(mrho,2)*s*(-6. + delta - 8.*C4*s) +
                  pow(mrho,4)*(12. - 1.*pow(delta,2) + 8.*C4*(6. + delta)*s + 32.*pow(C4,2)*pow(s,2)))*(0. + t2))/pow(mrho,4) -
             (0.125*(-2. + 1.*delta)*(2. + 1.*delta - 8.*C4*pow(mrho,2))*(0. + 1.*pow(t2,2)))/pow(mrho,2) -
             (1.*(0.5 - 0.125*pow(delta,2) - 2.*C4*pow(mrho,2) + 1.*C4*delta*pow(mrho,2))*(0. + 1.*pow(t2,2)))/pow(mrho,2) -
             0.5*(eta1*eta2*(1. - 2.*C4*pow(mrho,2)) + pow(eta1,2)*(-0.5 + 1.*C4*pow(mrho,2)) + pow(eta2,2)*(-0.5 + 1.*C4*pow(mrho,2)))*
              (0. + 1.*pow(t2,2)) + 0.0625*pow(1.*eta1 - 1.*eta2,4)*(1.*pow(mpion,2) + 0.5*pow(mrho,2) - 0.5*s)*(0. + 1.*pow(t2,2)) +
             0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*(0. + 1.*pow(t2,2)) +
             0.010416666666666666*pow(eta1 - 1.*eta2,4)*(0. + 1.*pow(t2,3)) -
             0.020833333333333332*pow(1.*eta1 - 1.*eta2,4)*(0. + 1.*pow(t2,3)) -
             (0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
                     pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) + pow(ma1,2)*pow(mpion,2)*
                      (8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
                     pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
                  pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
                     pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                        2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s +
                        pow(mpion,4)*(6.*pow(mrho,2) + 2.*s))) +
                  pow(eta1,2)*(1.*pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                        2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
                        pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
                     pow(mpion,2)*(1.*pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                        pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/(1.*pow(ma1,2) - 1.*t1) -
             (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*t1) -
             0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(1.*pow(Gammaa1,2)*pow(ma1,2) - 3.*pow(ma1,4) - 2.*pow(mpion,4) +
                   2.*pow(mpion,2)*pow(mrho,2) - 1.*pow(mrho,4) + pow(ma1,2)*(4.*pow(mpion,2) - 4.*pow(mrho,2) - 4.*s) +
                   2.*pow(mrho,2)*s - 2.*pow(s,2)) + pow(eta1,2)*
                 (1.*pow(Gammaa1,2)*pow(ma1,2) - 3.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) - 1.*pow(mrho,4) +
                   pow(ma1,2)*(4.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) + 4.*pow(mrho,2)*s - 2.*pow(s,2)) +
                eta1*eta2*(-2.*pow(Gammaa1,2)*pow(ma1,2) + 6.*pow(ma1,4) + 4.*pow(mpion,4) - 2.*pow(mrho,4) - 4.*pow(mrho,2)*s +
                   4.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) + 8.*s)))*(1.*pow(mrho,2) - 1.*s - 1.*t1) +
             0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) - 1.*pow(mrho,2) - 1.*s) + eta1*(pow(ma1,2) - 1.*pow(mrho,2) + s))*
              pow(1.*pow(mrho,2) - 1.*s - 1.*t1,2) + 0.010416666666666666*pow(eta1 - 1.*eta2,4)*pow(1.*pow(mrho,2) - 1.*s - 1.*t1,3) -
             (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/
              (1.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s - 1.*t1) + (0.5*pow(-2. + delta,2)*pow(mpion,2)*t1)/pow(mrho,2) +
             0.25*(-2. + delta)*(eta1 - 1.*eta2)*(-0.5*eta2*pow(ma1,2) + 1.*eta1*pow(mpion,2) - 1.*eta2*pow(mpion,2) + 0.5*eta1*s -
                0.5*eta2*s)*t1 - (0.25*(pow(mpion,2)*(12. + 1.*pow(delta,2) - 16.*C4*pow(mrho,2) + delta*(-8. + 8.*C4*pow(mrho,2))) +
                  (-4. - 3.*pow(delta,2) - 16.*C4*pow(mrho,2) + delta*(8. + 8.*C4*pow(mrho,2)))*s)*t1)/pow(mrho,2) -
             0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(1.*pow(mpion,2) - 0.5*s) +
                eta2*(-0.5*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 1.*s))*t1 -
             (0.25*(-2. + 1.*delta)*(-8.*C4*pow(mrho,4) + pow(mpion,2)*(2. + 1.*delta - 8.*C4*pow(mrho,2)) + (-2. - 3.*delta)*s +
                  pow(mrho,2)*(2. + 1.*delta + 16.*C4*s))*t1)/pow(mrho,2) -
             (0.25*(32.*pow(C4,2)*pow(mrho,8) + 2.*pow(delta,2)*pow(s,2) + 8.*C4*pow(mrho,6)*(-6. + delta - 8.*C4*s) +
                  2.*delta*pow(mrho,2)*s*(-6. + delta - 8.*C4*s) +
                  pow(mrho,4)*(12. - 1.*pow(delta,2) + 8.*C4*(6. + delta)*s + 32.*pow(C4,2)*pow(s,2)))*t1)/pow(mrho,4) -
             0.09375*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-2.*pow(ma1,4) - 4.*pow(mpion,4) + 0.6666666666666666*pow(mrho,4) +
                   pow(ma1,2)*(5.333333333333333*pow(mpion,2) - 2.6666666666666665*s) + 2.6666666666666665*pow(mpion,2)*s +
                   1.3333333333333333*pow(mrho,2)*s - 1.3333333333333333*pow(s,2)) +
                pow(eta1,2)*(1.*pow(ma1,4) + 2.*pow(mpion,4) + 0.3333333333333333*pow(mrho,4) +
                   pow(mpion,2)*(2.*pow(mrho,2) - 1.3333333333333333*s) - 1.3333333333333333*pow(mrho,2)*s + 0.6666666666666666*pow(s,2) +
                   pow(ma1,2)*(-2.6666666666666665*pow(mpion,2) - 1.3333333333333333*pow(mrho,2) + 1.3333333333333333*s)) +
                pow(eta2,2)*(1.*pow(ma1,4) + 2.*pow(mpion,4) + 0.3333333333333333*pow(mrho,4) +
                   pow(mpion,2)*(-2.*pow(mrho,2) - 1.3333333333333333*s) - 0.6666666666666666*pow(mrho,2)*s + 0.6666666666666666*pow(s,2) +
                   pow(ma1,2)*(-2.6666666666666665*pow(mpion,2) + 1.3333333333333333*pow(mrho,2) + 1.3333333333333333*s)))*t1 -
             0.0625*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(pow(Gammaa1,2)*pow(ma1,2) - 1.*pow(ma1,4) - 2.*pow(mpion,4) -
                   2.*pow(mpion,2)*pow(mrho,2) + 2.*pow(mrho,4) + pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) -
                   3.*pow(mrho,2)*s + pow(s,2)) + pow(eta1,2)*
                 (pow(Gammaa1,2)*pow(ma1,2) - 1.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) +
                   pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) + pow(mrho,2)*s + pow(s,2)) +
                eta1*eta2*(-2.*pow(Gammaa1,2)*pow(ma1,2) + 2.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) +
                   2.*pow(mrho,4) - 2.*pow(s,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)))*t1 +
             (1.*(pow(eta1,2)*(0.5*pow(mrho,4) - 1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                     pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) - 0.5*pow(mrho,2)*s + 0.25*delta*pow(mrho,2)*s -
                     0.25*delta*pow(s,2) + 1.*C4*pow(mrho,2)*pow(s,2)) +
                  pow(eta2,2)*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                     pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) - 0.5*pow(mrho,2)*s + 0.75*delta*pow(mrho,2)*s -
                     2.*C4*pow(mrho,4)*s - 0.25*delta*pow(s,2) + 1.*C4*pow(mrho,2)*pow(s,2)) +
                  eta1*eta2*(pow(ma1,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) + pow(mpion,2)*(-2.*pow(mrho,2) + 4.*C4*pow(mrho,4)) +
                     s*(2.*C4*pow(mrho,4) + 0.5*delta*s + pow(mrho,2)*(1. - 1.*delta - 2.*C4*s))))*t1)/pow(mrho,2) +
             (0.125*(-2. + 1.*delta)*(2. + 1.*delta - 8.*C4*pow(mrho,2))*pow(t1,2))/pow(mrho,2) +
             (1.*(0.5 - 0.125*pow(delta,2) - 2.*C4*pow(mrho,2) + 1.*C4*delta*pow(mrho,2))*pow(t1,2))/pow(mrho,2) +
             0.5*(eta1*eta2*(1. - 2.*C4*pow(mrho,2)) + pow(eta1,2)*(-0.5 + 1.*C4*pow(mrho,2)) + pow(eta2,2)*(-0.5 + 1.*C4*pow(mrho,2)))*
              pow(t1,2) - 0.0625*pow(1.*eta1 - 1.*eta2,4)*(1.*pow(mpion,2) + 0.5*pow(mrho,2) - 0.5*s)*pow(t1,2) -
             0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(t1,2) -
             0.010416666666666666*pow(eta1 - 1.*eta2,4)*pow(t1,3) + 0.020833333333333332*pow(1.*eta1 - 1.*eta2,4)*pow(t1,3) +
             (2.*(1.*eta1 - 1.*eta2)*(eta2*(0.375*pow(mrho,4) - 0.0625*delta*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                     pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(mpion,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                     0.125*pow(mrho,2)*s - 0.1875*delta*pow(mrho,2)*s + 1.*C4*pow(mrho,4)*s + 0.125*delta*pow(s,2) -
                     0.5*C4*pow(mrho,2)*pow(s,2)) + eta1*(0.25*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                     pow(mpion,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                     0.125*delta*pow(s,2) + 0.5*C4*pow(mrho,2)*pow(s,2)))*
                (1.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s + 1.*t1))/pow(mrho,2) -
             0.5*pow(1.*eta1 - 1.*eta2,2)*(-0.5 + 1.*C4*pow(mrho,2))*
              pow(1.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s + 1.*t1,2) +
             (2.*(1.*eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*(0.375*pow(mrho,4) - 0.0625*delta*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                     pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(mpion,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                     0.125*pow(mrho,2)*s - 0.1875*delta*pow(mrho,2)*s + 1.*C4*pow(mrho,4)*s + 0.125*delta*pow(s,2) -
                     0.5*C4*pow(mrho,2)*pow(s,2)) + eta1*(0.25*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                     pow(mpion,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                     0.125*delta*pow(s,2) + 0.5*C4*pow(mrho,2)*pow(s,2)))*
                atan((0. + pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + t2)/(Gammaa1*ma1)))/pow(mrho,2) +
             (0.25*(-2. + delta)*(eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*
                   (-1.*pow(ma1,6) + pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(ma1,2) + 0.5*pow(mrho,2) - 1.*s) +
                     pow(ma1,4)*(2.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) + pow(mpion,4)*(-1.5*pow(mrho,2) + 1.*s) +
                     pow(ma1,2)*pow(mpion,2)*(-1.*pow(mpion,2) - 1.*pow(mrho,2) + 2.*s)) +
                  eta1*(pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) + 0.5*s) + pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
                     pow(ma1,2)*(-2.*pow(mpion,4) - 1.*pow(mpion,2)*s) +
                     pow(mpion,2)*(1.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) + 1.*pow(mrho,2)*s)))*
                atan((0. + pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + t2)/(Gammaa1*ma1)))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) -
             (0.25*(-2. + delta)*(eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*
                   (-1.*pow(ma1,6) - 2.*pow(mpion,4)*pow(mrho,2) - 1.*pow(mpion,2)*pow(mrho,4) +
                     pow(ma1,4)*(2.*pow(mpion,2) + 2.*pow(mrho,2) - 1.5*s) + 2.5*pow(mpion,4)*s + 3.*pow(mpion,2)*pow(mrho,2)*s +
                     0.5*pow(mrho,4)*s - 2.*pow(mpion,2)*pow(s,2) - 1.*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) +
                     pow(Gammaa1,2)*(-1.*pow(ma1,4) + 0.5*pow(ma1,2)*s) +
                     pow(ma1,2)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) + 1.*pow(mrho,2)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 1.*s))) +
                  eta1*(1.*pow(mpion,6) + 4.*pow(mpion,4)*pow(mrho,2) + 1.*pow(mpion,2)*pow(mrho,4) +
                     pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) + pow(ma1,4)*(1.*pow(mpion,2) - 0.5*s) - 4.5*pow(mpion,4)*s -
                     4.*pow(mpion,2)*pow(mrho,2)*s - 0.5*pow(mrho,4)*s + 3.*pow(mpion,2)*pow(s,2) + 1.*pow(mrho,2)*pow(s,2) -
                     0.5*pow(s,3) + pow(ma1,2)*(-2.*pow(mpion,4) + (1.*pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 3.*s))))*
                atan((0. + pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + t2)/(Gammaa1*ma1)))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
                2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
             (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8) + pow(mpion,8) -
                     2.*pow(mpion,6)*pow(mrho,2) + pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                        2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s +
                        pow(mpion,4)*(6.*pow(mrho,2) + 2.*s)) +
                     pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) - 6.*pow(mpion,4) - 1.*pow(mrho,4) +
                        pow(ma1,2)*(12.*pow(mpion,2) - 6.*pow(mrho,2) - 6.*s) + 2.*pow(mrho,2)*s - 2.*pow(s,2) +
                        pow(mpion,2)*(6.*pow(mrho,2) + 4.*s))) +
                  eta1*eta2*(-2.*pow(Gammaa1,4)*pow(ma1,4) - 2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
                     pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) + pow(ma1,2)*pow(mpion,2)*
                      (8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
                     pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
                     pow(Gammaa1,2)*pow(ma1,2)*(12.*pow(ma1,4) + 12.*pow(mpion,4) - 2.*pow(mrho,4) - 8.*pow(mpion,2)*s -
                        4.*pow(mrho,2)*s + 4.*pow(s,2) + pow(ma1,2)*(-24.*pow(mpion,2) + 12.*s))) +
                  pow(eta1,2)*(pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                        2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
                        pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
                     pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) - 6.*pow(mpion,4) - 1.*pow(mrho,4) +
                        pow(ma1,2)*(12.*pow(mpion,2) + 6.*pow(mrho,2) - 6.*s) + 4.*pow(mrho,2)*s - 2.*pow(s,2) +
                        pow(mpion,2)*(-6.*pow(mrho,2) + 4.*s)) +
                     pow(mpion,2)*(pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                        pow(mpion,2)*(pow(mrho,4) - 2.*pow(mrho,2)*s))))*
                atan((0. + pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + t2)/(Gammaa1*ma1)))/(Gammaa1*ma1) -
             (0.0625*pow(eta1 - 1.*eta2,2)*Gammaa1*ma1*(eta1*eta2*
                   (-2.*pow(Gammaa1,4)*pow(ma1,4) + 14.*pow(ma1,8) + 14.*pow(mpion,8) + 28.*pow(mpion,6)*pow(mrho,2) +
                     20.*pow(mpion,4)*pow(mrho,4) + 10.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) - 16.*pow(mpion,6)*s -
                     16.*pow(mpion,4)*pow(mrho,2)*s - 12.*pow(mpion,2)*pow(mrho,4)*s - 4.*pow(mrho,6)*s - 4.*pow(mpion,4)*pow(s,2) -
                     6.*pow(mpion,2)*pow(mrho,2)*pow(s,2) + 8.*pow(mpion,2)*pow(s,3) + 4.*pow(mrho,2)*pow(s,3) - 2.*pow(s,4) +
                     pow(ma1,6)*(-56.*pow(mpion,2) - 28.*pow(mrho,2) + 28.*s) +
                     pow(ma1,4)*(84.*pow(mpion,4) + 24.*pow(mrho,4) + pow(mpion,2)*(84.*pow(mrho,2) - 72.*s) - 36.*pow(mrho,2)*s +
                        12.*pow(s,2)) + pow(Gammaa1,2)*pow(ma1,2)*
                      (-4.*pow(ma1,4) - 4.*pow(mpion,4) + pow(ma1,2)*(8.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) +
                        (4.*pow(mrho,2) - 4.*s)*s + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)) +
                     pow(ma1,2)*(-56.*pow(mpion,6) - 10.*pow(mrho,6) + 18.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3) +
                        pow(mpion,4)*(-84.*pow(mrho,2) + 60.*s) + pow(mpion,2)*(-48.*pow(mrho,4) + 60.*pow(mrho,2)*s - 12.*pow(s,2)))) +
                  pow(eta1,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) - 7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
                     7.*pow(mpion,4)*pow(mrho,4) - 2.*pow(mpion,2)*pow(mrho,6) +
                     pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) + 8.*pow(mpion,6)*s + 11.*pow(mpion,4)*pow(mrho,2)*s +
                     6.*pow(mpion,2)*pow(mrho,4)*s + 1.*pow(mrho,6)*s + 2.*pow(mpion,4)*pow(s,2) - 1.*pow(mrho,4)*pow(s,2) -
                     4.*pow(mpion,2)*pow(s,3) - 1.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
                     pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,4) + 2.*pow(mpion,4) + 1.*pow(mrho,4) +
                        pow(mpion,2)*(2.*pow(mrho,2) - 4.*s) - 1.*pow(mrho,2)*s + 2.*pow(s,2) +
                        pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
                     pow(ma1,4)*(-42.*pow(mpion,4) - 9.*pow(mrho,4) + 21.*pow(mrho,2)*s - 6.*pow(s,2) +
                        pow(mpion,2)*(-42.*pow(mrho,2) + 36.*s)) +
                     pow(ma1,2)*(28.*pow(mpion,6) + 2.*pow(mrho,6) + pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) - 9.*pow(mrho,4)*s +
                        6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(18.*pow(mrho,4) - 36.*pow(mrho,2)*s + 6.*pow(s,2)))) +
                  pow(eta2,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) - 7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
                     1.*pow(mpion,4)*pow(mrho,4) + 6.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) +
                     pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) + 8.*pow(mpion,6)*s - 1.*pow(mpion,4)*pow(mrho,2)*s -
                     16.*pow(mpion,2)*pow(mrho,4)*s - 7.*pow(mrho,6)*s + 2.*pow(mpion,4)*pow(s,2) +
                     14.*pow(mpion,2)*pow(mrho,2)*pow(s,2) + 9.*pow(mrho,4)*pow(s,2) - 4.*pow(mpion,2)*pow(s,3) -
                     5.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
                     pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,4) + 2.*pow(mpion,4) + 3.*pow(mrho,4) +
                        pow(mpion,2)*(2.*pow(mrho,2) - 4.*s) - 5.*pow(mrho,2)*s + 2.*pow(s,2) +
                        pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
                     pow(ma1,4)*(-42.*pow(mpion,4) - 3.*pow(mrho,4) + 9.*pow(mrho,2)*s - 6.*pow(s,2) +
                        pow(mpion,2)*(-42.*pow(mrho,2) + 36.*s)) +
                     pow(ma1,2)*(28.*pow(mpion,6) - 4.*pow(mrho,6) + pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) + 9.*pow(mrho,4)*s -
                        6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s + 6.*pow(s,2)))))*
                atan((0. + pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + t2)/(Gammaa1*ma1)))/
              (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
                4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) -
             (2.*(1.*eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*(0.375*pow(mrho,4) - 0.0625*delta*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                     pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(mpion,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                     0.125*pow(mrho,2)*s - 0.1875*delta*pow(mrho,2)*s + 1.*C4*pow(mrho,4)*s + 0.125*delta*pow(s,2) -
                     0.5*C4*pow(mrho,2)*pow(s,2)) + eta1*(0.25*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                     pow(mpion,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                     0.125*delta*pow(s,2) + 0.5*C4*pow(mrho,2)*pow(s,2)))*
                atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + t1)/(Gammaa1*ma1)))/pow(mrho,2) -
             (0.25*(-2. + delta)*(eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*
                   (-1.*pow(ma1,6) + pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(ma1,2) + 0.5*pow(mrho,2) - 1.*s) +
                     pow(ma1,4)*(2.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) + pow(mpion,4)*(-1.5*pow(mrho,2) + 1.*s) +
                     pow(ma1,2)*pow(mpion,2)*(-1.*pow(mpion,2) - 1.*pow(mrho,2) + 2.*s)) +
                  eta1*(pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) + 0.5*s) + pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
                     pow(ma1,2)*(-2.*pow(mpion,4) - 1.*pow(mpion,2)*s) +
                     pow(mpion,2)*(1.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) + 1.*pow(mrho,2)*s)))*
                atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + t1)/(Gammaa1*ma1)))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) +
             (0.25*(-2. + delta)*(eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*
                   (-1.*pow(ma1,6) - 2.*pow(mpion,4)*pow(mrho,2) - 1.*pow(mpion,2)*pow(mrho,4) +
                     pow(ma1,4)*(2.*pow(mpion,2) + 2.*pow(mrho,2) - 1.5*s) + 2.5*pow(mpion,4)*s + 3.*pow(mpion,2)*pow(mrho,2)*s +
                     0.5*pow(mrho,4)*s - 2.*pow(mpion,2)*pow(s,2) - 1.*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) +
                     pow(Gammaa1,2)*(-1.*pow(ma1,4) + 0.5*pow(ma1,2)*s) +
                     pow(ma1,2)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) + 1.*pow(mrho,2)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 1.*s))) +
                  eta1*(1.*pow(mpion,6) + 4.*pow(mpion,4)*pow(mrho,2) + 1.*pow(mpion,2)*pow(mrho,4) +
                     pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) + pow(ma1,4)*(1.*pow(mpion,2) - 0.5*s) - 4.5*pow(mpion,4)*s -
                     4.*pow(mpion,2)*pow(mrho,2)*s - 0.5*pow(mrho,4)*s + 3.*pow(mpion,2)*pow(s,2) + 1.*pow(mrho,2)*pow(s,2) -
                     0.5*pow(s,3) + pow(ma1,2)*(-2.*pow(mpion,4) + (1.*pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 3.*s))))*
                atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + t1)/(Gammaa1*ma1)))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
                2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) -
             (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8) + pow(mpion,8) -
                     2.*pow(mpion,6)*pow(mrho,2) + pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                        2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s +
                        pow(mpion,4)*(6.*pow(mrho,2) + 2.*s)) +
                     pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) - 6.*pow(mpion,4) - 1.*pow(mrho,4) +
                        pow(ma1,2)*(12.*pow(mpion,2) - 6.*pow(mrho,2) - 6.*s) + 2.*pow(mrho,2)*s - 2.*pow(s,2) +
                        pow(mpion,2)*(6.*pow(mrho,2) + 4.*s))) +
                  eta1*eta2*(-2.*pow(Gammaa1,4)*pow(ma1,4) - 2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
                     pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) + pow(ma1,2)*pow(mpion,2)*
                      (8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
                     pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
                     pow(Gammaa1,2)*pow(ma1,2)*(12.*pow(ma1,4) + 12.*pow(mpion,4) - 2.*pow(mrho,4) - 8.*pow(mpion,2)*s -
                        4.*pow(mrho,2)*s + 4.*pow(s,2) + pow(ma1,2)*(-24.*pow(mpion,2) + 12.*s))) +
                  pow(eta1,2)*(pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                        2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
                        pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
                     pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) - 6.*pow(mpion,4) - 1.*pow(mrho,4) +
                        pow(ma1,2)*(12.*pow(mpion,2) + 6.*pow(mrho,2) - 6.*s) + 4.*pow(mrho,2)*s - 2.*pow(s,2) +
                        pow(mpion,2)*(-6.*pow(mrho,2) + 4.*s)) +
                     pow(mpion,2)*(pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                        pow(mpion,2)*(pow(mrho,4) - 2.*pow(mrho,2)*s))))*
                atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + t1)/(Gammaa1*ma1)))/(Gammaa1*ma1) +
             (0.0625*pow(eta1 - 1.*eta2,2)*Gammaa1*ma1*(eta1*eta2*
                   (-2.*pow(Gammaa1,4)*pow(ma1,4) + 14.*pow(ma1,8) + 14.*pow(mpion,8) + 28.*pow(mpion,6)*pow(mrho,2) +
                     20.*pow(mpion,4)*pow(mrho,4) + 10.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) - 16.*pow(mpion,6)*s -
                     16.*pow(mpion,4)*pow(mrho,2)*s - 12.*pow(mpion,2)*pow(mrho,4)*s - 4.*pow(mrho,6)*s - 4.*pow(mpion,4)*pow(s,2) -
                     6.*pow(mpion,2)*pow(mrho,2)*pow(s,2) + 8.*pow(mpion,2)*pow(s,3) + 4.*pow(mrho,2)*pow(s,3) - 2.*pow(s,4) +
                     pow(ma1,6)*(-56.*pow(mpion,2) - 28.*pow(mrho,2) + 28.*s) +
                     pow(ma1,4)*(84.*pow(mpion,4) + 24.*pow(mrho,4) + pow(mpion,2)*(84.*pow(mrho,2) - 72.*s) - 36.*pow(mrho,2)*s +
                        12.*pow(s,2)) + pow(Gammaa1,2)*pow(ma1,2)*
                      (-4.*pow(ma1,4) - 4.*pow(mpion,4) + pow(ma1,2)*(8.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) +
                        (4.*pow(mrho,2) - 4.*s)*s + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)) +
                     pow(ma1,2)*(-56.*pow(mpion,6) - 10.*pow(mrho,6) + 18.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3) +
                        pow(mpion,4)*(-84.*pow(mrho,2) + 60.*s) + pow(mpion,2)*(-48.*pow(mrho,4) + 60.*pow(mrho,2)*s - 12.*pow(s,2)))) +
                  pow(eta1,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) - 7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
                     7.*pow(mpion,4)*pow(mrho,4) - 2.*pow(mpion,2)*pow(mrho,6) +
                     pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) + 8.*pow(mpion,6)*s + 11.*pow(mpion,4)*pow(mrho,2)*s +
                     6.*pow(mpion,2)*pow(mrho,4)*s + 1.*pow(mrho,6)*s + 2.*pow(mpion,4)*pow(s,2) - 1.*pow(mrho,4)*pow(s,2) -
                     4.*pow(mpion,2)*pow(s,3) - 1.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
                     pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,4) + 2.*pow(mpion,4) + 1.*pow(mrho,4) +
                        pow(mpion,2)*(2.*pow(mrho,2) - 4.*s) - 1.*pow(mrho,2)*s + 2.*pow(s,2) +
                        pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
                     pow(ma1,4)*(-42.*pow(mpion,4) - 9.*pow(mrho,4) + 21.*pow(mrho,2)*s - 6.*pow(s,2) +
                        pow(mpion,2)*(-42.*pow(mrho,2) + 36.*s)) +
                     pow(ma1,2)*(28.*pow(mpion,6) + 2.*pow(mrho,6) + pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) - 9.*pow(mrho,4)*s +
                        6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(18.*pow(mrho,4) - 36.*pow(mrho,2)*s + 6.*pow(s,2)))) +
                  pow(eta2,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) - 7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
                     1.*pow(mpion,4)*pow(mrho,4) + 6.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) +
                     pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) + 8.*pow(mpion,6)*s - 1.*pow(mpion,4)*pow(mrho,2)*s -
                     16.*pow(mpion,2)*pow(mrho,4)*s - 7.*pow(mrho,6)*s + 2.*pow(mpion,4)*pow(s,2) +
                     14.*pow(mpion,2)*pow(mrho,2)*pow(s,2) + 9.*pow(mrho,4)*pow(s,2) - 4.*pow(mpion,2)*pow(s,3) -
                     5.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
                     pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,4) + 2.*pow(mpion,4) + 3.*pow(mrho,4) +
                        pow(mpion,2)*(2.*pow(mrho,2) - 4.*s) - 5.*pow(mrho,2)*s + 2.*pow(s,2) +
                        pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
                     pow(ma1,4)*(-42.*pow(mpion,4) - 3.*pow(mrho,4) + 9.*pow(mrho,2)*s - 6.*pow(s,2) +
                        pow(mpion,2)*(-42.*pow(mrho,2) + 36.*s)) +
                     pow(ma1,2)*(28.*pow(mpion,6) - 4.*pow(mrho,6) + pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) + 9.*pow(mrho,4)*s -
                        6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s + 6.*pow(s,2)))))*
                atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + t1)/(Gammaa1*ma1)))/
              (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
                4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
             0.0625*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
                   pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) - 2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
                   pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
                pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) + pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s +
                   pow(mpion,4)*(-3.*pow(mrho,2) + s) + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
                   pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2)))\
                 + pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
                   pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
                   pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2)))
                )*log(fabs(-pow(ma1,2) + t2)) - (0.03125*pow(eta1 - 1.*eta2,2)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s)*
                (eta1*eta2*(-2.*pow(ma1,8) + pow(ma1,6)*(8.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) +
                     pow(ma1,4)*(-12.*pow(mpion,4) - 4.*pow(mrho,4) + 4.*pow(mrho,2)*s + pow(mpion,2)*(-12.*pow(mrho,2) + 8.*s)) +
                     pow(mpion,2)*(-2.*pow(mpion,6) - 4.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 4.*pow(mrho,4)*s -
                        2.*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(-8.*pow(mrho,4) + 8.*pow(mrho,2)*s)) +
                     pow(ma1,2)*(8.*pow(mpion,6) + 2.*pow(mrho,6) + pow(mpion,4)*(12.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,4)*s -
                        2.*pow(mrho,2)*pow(s,2) + 2.*pow(s,3) + pow(mpion,2)*(8.*pow(mrho,4) - 4.*pow(mrho,2)*s - 4.*pow(s,2)))) +
                  pow(eta2,2)*(pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                     pow(mpion,4)*(pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 1.*pow(mrho,2)*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) + pow(mrho,2)*s) +
                     pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mrho,6) - 5.*pow(mrho,4)*s + 4.*pow(mrho,2)*pow(s,2) - 1.*pow(s,3) +
                        pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) + pow(mpion,2)*(2.*pow(mrho,4) - 4.*pow(mrho,2)*s + 2.*pow(s,2)))) +
                  pow(eta1,2)*(pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 3.*pow(mrho,2)*s) +
                     pow(mpion,2)*(pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) + pow(mrho,2)*s*(-2.*pow(mrho,2) + 2.*s) +
                        pow(mpion,2)*(3.*pow(mrho,4) - 5.*pow(mrho,2)*s)) +
                     pow(ma1,2)*(-4.*pow(mpion,6) + pow(mrho,4)*s - 1.*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) +
                        pow(mpion,2)*(-2.*pow(mrho,4) + 4.*pow(mrho,2)*s + 2.*pow(s,2)))))*log(fabs(-pow(ma1,2) + t2)))/
              (0.25*pow(Gammaa1,2)*pow(ma1,2) + 1.*pow(ma1,4) + 1.*pow(mpion,4) + 1.*pow(mpion,2)*pow(mrho,2) + 0.25*pow(mrho,4) -
                1.*pow(mpion,2)*s - 0.5*pow(mrho,2)*s + 0.25*pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s)) -
             (1.*(eta1*eta2*(pow(ma1,8)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                     pow(ma1,6)*(-1.*pow(mrho,4) + 2.*C4*pow(mrho,6) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*C4*pow(mrho,4)) +
                        0.5*delta*pow(s,2) + pow(mrho,2)*s*(2. - 1.*delta - 2.*C4*s)) +
                     pow(ma1,2)*(pow(mpion,6)*(-4.*pow(mrho,2) + 8.*C4*pow(mrho,4)) +
                        pow(mrho,4)*s*((0.5 - 0.25*delta)*pow(mrho,2) + (-0.5 + 0.25*delta)*s) +
                        pow(mpion,4)*(10.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,4)*(-3. - 1.*delta - 8.*C4*s) +
                           pow(mrho,2)*s*(2. + 1.*delta - 2.*C4*s)) +
                        pow(mpion,2)*(2.*C4*pow(mrho,8) - 0.5*delta*pow(s,3) + pow(mrho,6)*(1. - 1.*delta - 2.*C4*s) +
                           pow(mrho,4)*s*(-1. + 1.*delta - 2.*C4*s) + pow(mrho,2)*pow(s,2)*(1. + 2.*C4*s))) +
                     pow(mpion,2)*pow(mrho,2)*(pow(mpion,6)*(1. - 2.*C4*pow(mrho,2)) +
                        pow(mrho,4)*((-0.5 + 0.25*delta)*pow(mrho,2) + (0.5 - 0.25*delta)*s) +
                        pow(mpion,2)*(-2.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,2)*s*(-1.5 - 0.25*delta - 2.*C4*s) +
                           pow(mrho,4)*(1. + 4.*C4*s)) + pow(mpion,4)*
                         (-4.*C4*pow(mrho,4) - 1.*delta*s + pow(mrho,2)*(1. + 0.5*delta + 4.*C4*s))) +
                     pow(ma1,4)*(pow(mpion,4)*(6.*pow(mrho,2) - 12.*C4*pow(mrho,4)) +
                        pow(mpion,2)*(-8.*C4*pow(mrho,6) - 1.*delta*pow(s,2) + pow(mrho,4)*(3. + 0.5*delta + 4.*C4*s) +
                           pow(mrho,2)*s*(-4. + 1.*delta + 4.*C4*s)) +
                        s*(-2.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,2)*s*(1. - 1.5*delta - 2.*C4*s) +
                           pow(mrho,4)*(-1.5 + 1.25*delta + 4.*C4*s)))) +
                  pow(eta1,2)*(pow(ma1,8)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                     pow(ma1,6)*(-2.*C4*pow(mrho,6) + pow(mpion,2)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) - 0.25*delta*pow(s,2) +
                        pow(mrho,4)*(1. + 1.*C4*s) + pow(mrho,2)*s*(-1. + 0.25*delta + 1.*C4*s)) +
                     pow(ma1,4)*(1.*C4*pow(mrho,8) + pow(mpion,4)*(-3.*pow(mrho,2) + 6.*C4*pow(mrho,4)) - 0.25*delta*pow(s,3) +
                        pow(mrho,6)*(-0.5 - 1.*C4*s) + pow(mrho,4)*s*(1.5 - 0.5*delta - 1.*C4*s) +
                        pow(mrho,2)*pow(s,2)*(-0.5 + 0.5*delta + 1.*C4*s) +
                        pow(mpion,2)*(7.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,4)*(-3. - 0.25*delta - 5.*C4*s) +
                           pow(mrho,2)*s*(2. + 0.25*delta - 2.*C4*s))) +
                     pow(mpion,2)*pow(mrho,2)*(pow(mpion,6)*(-0.5 + 1.*C4*pow(mrho,2)) +
                        pow(mrho,4)*((0.5 - 0.25*delta)*pow(mrho,2) + (-0.5 + 0.25*delta)*s) +
                        pow(mpion,4)*(3.*C4*pow(mrho,4) + 0.75*delta*s + pow(mrho,2)*(-1. - 0.25*delta - 3.*C4*s)) +
                        pow(mpion,2)*(2.*C4*pow(mrho,6) - 0.5*delta*pow(s,2) + pow(mrho,4)*(-1. - 4.*C4*s) +
                           pow(mrho,2)*s*(1.5 + 0.25*delta + 2.*C4*s))) +
                     pow(ma1,2)*(pow(mpion,6)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                        pow(mrho,4)*s*((-0.5 + 0.25*delta)*pow(mrho,2) + (0.5 - 0.25*delta)*s) +
                        pow(mpion,2)*(-3.*C4*pow(mrho,8) + 0.25*delta*pow(s,3) + pow(mrho,4)*s*(-1. - 0.75*delta - 1.*C4*s) +
                           pow(mrho,2)*pow(s,2)*(-0.5 + 0.5*delta - 1.*C4*s) + pow(mrho,6)*(0.5 + 0.5*delta + 5.*C4*s)) +
                        pow(mpion,4)*(-8.*C4*pow(mrho,6) - 0.25*delta*pow(s,2) + pow(mrho,2)*s*(-1. - 1.25*delta + 1.*C4*s) +
                           pow(mrho,4)*(3. + 0.5*delta + 7.*C4*s)))) +
                  pow(eta2,2)*(pow(ma1,8)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                     pow(mpion,6)*pow(mrho,2)*(-0.25*delta*pow(mrho,2) + 1.*C4*pow(mrho,4) + pow(mpion,2)*(-0.5 + 1.*C4*pow(mrho,2)) +
                        0.25*delta*s - 1.*C4*pow(mrho,2)*s) +
                     pow(ma1,6)*(pow(mpion,2)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                        s*(-1.*C4*pow(mrho,4) - 0.25*delta*s + pow(mrho,2)*(-1. + 0.75*delta + 1.*C4*s))) +
                     pow(ma1,4)*(-1.*C4*pow(mrho,8) + pow(mpion,4)*(-3.*pow(mrho,2) + 6.*C4*pow(mrho,4)) - 0.25*delta*pow(s,3) +
                        pow(mrho,4)*s*(-0.75*delta - 3.*C4*s) + pow(mrho,2)*pow(s,2)*(-0.5 + 1.*delta + 1.*C4*s) +
                        pow(mrho,6)*(0.5 + 3.*C4*s) + pow(mpion,2)*
                         (delta*(-0.25*pow(mrho,4) - 1.25*pow(mrho,2)*s + 0.5*pow(s,2)) +
                           pow(mrho,2)*(1.*C4*pow(mrho,4) + 2.*s + 1.*C4*pow(mrho,2)*s - 2.*C4*pow(s,2)))) +
                     pow(ma1,2)*pow(mpion,2)*(1.*C4*pow(mrho,8) + pow(mpion,4)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                        0.25*delta*pow(s,3) + pow(mrho,6)*(-1.5 + 0.5*delta - 3.*C4*s) + pow(mrho,2)*pow(s,2)*(-0.5 - 0.5*delta - 1.*C4*s) +
                        pow(mrho,4)*s*(2. - 0.25*delta + 3.*C4*s) +
                        pow(mpion,2)*(delta*(0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
                           pow(mrho,2)*(-2.*C4*pow(mrho,4) - 1.*s + 1.*C4*pow(mrho,2)*s + 1.*C4*pow(s,2))))))*log(fabs(-pow(ma1,2) + t2)))/
              ((pow(ma1,2) - 1.*pow(mpion,2))*pow(mrho,2)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s)) +
             (0.5*pow(mpion,2)*(pow(eta2,2)*pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s +
                     (-2. + 1.*delta)*pow(s,2)) + eta1*eta2*(pow(mpion,2)*
                      ((4. - 2.*delta)*pow(mrho,4) + (-8. + 4.*delta)*pow(mrho,2)*s + (4. - 2.*delta)*pow(s,2)) +
                     pow(mrho,2)*((-1. + 0.5*delta)*pow(mrho,4) + (2. - 1.*delta)*pow(mrho,2)*s + (-1. + 0.5*delta)*pow(s,2))) +
                  pow(eta1,2)*(pow(mrho,2)*((1. - 0.5*delta)*pow(mrho,4) + (-2. + 1.*delta)*pow(mrho,2)*s + (1. - 0.5*delta)*pow(s,2)) +
                     pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))))*
                log(fabs(-pow(mpion,2) + t2)))/((-1.*pow(ma1,2) + 1.*pow(mpion,2))*(pow(mrho,2) - 1.*s)) -
             (0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*
                (eta1*(pow(mpion,4)*(-1.*pow(mrho,2) + 1.*s) +
                     pow(ma1,2)*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + (-0.5*pow(mrho,2) + 0.5*s)*s) +
                     pow(mpion,2)*(-1.*pow(mrho,4) + 2.5*pow(mrho,2)*s - 1.5*pow(s,2)) +
                     s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2))) +
                  eta2*(0.5*pow(mrho,6) + pow(mpion,4)*(1.*pow(mrho,2) - 1.*s) - 1.5*pow(mrho,4)*s + 1.5*pow(mrho,2)*pow(s,2) -
                     0.5*pow(s,3) + pow(mpion,2)*(1.5*pow(mrho,4) - 3.*pow(mrho,2)*s + 1.5*pow(s,2)) +
                     pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*pow(mrho,2)*s - 0.5*pow(s,2) + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*
                log(fabs(-pow(mpion,2) + t2)))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
                2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) -
             (0.5*pow(mpion,2)*(eta1*eta2*((1. - 0.5*delta)*pow(mrho,6) + (-4. + 2.*delta)*pow(mrho,4)*s +
                     (5. - 2.5*delta)*pow(mrho,2)*pow(s,2) + (-2. + 1.*delta)*pow(s,3) +
                     pow(mpion,2)*((4. - 2.*delta)*pow(mrho,4) + (-8. + 4.*delta)*pow(mrho,2)*s + (4. - 2.*delta)*pow(s,2))) +
                  pow(eta2,2)*((-1. + 0.5*delta)*pow(mrho,6) + (3. - 1.5*delta)*pow(mrho,4)*s + (-3. + 1.5*delta)*pow(mrho,2)*pow(s,2) +
                     (1. - 0.5*delta)*pow(s,3) + pow(mpion,2)*
                      ((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))) +
                  pow(eta1,2)*(s*((1. - 0.5*delta)*pow(mrho,4) + (-2. + 1.*delta)*pow(mrho,2)*s + (1. - 0.5*delta)*pow(s,2)) +
                     pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))))*
                log(fabs(-pow(mpion,2) - pow(mrho,2) + s + t2)))/
              ((pow(mrho,2) - 1.*s)*(-1.*pow(ma1,2) + pow(mpion,2) + pow(mrho,2) - 1.*s)) +
             0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + t2)) +
             (0.25000000000000006*pow(2. - 1.*delta,2)*(7.999999999999998*pow(mpion,4) - 5.999999999999998*pow(mpion,2)*s + 1.*pow(s,2))*
                log(fabs(-pow(mpion,2) - pow(mrho,2) + s + t2)))/(pow(mrho,2) - 1.*s) -
             (1.*(-2. + 1.*delta)*((0.5 - 0.25*delta)*pow(mrho,2)*s +
                  pow(mpion,2)*(4.*C4*pow(mrho,4) + 1.*delta*s + pow(mrho,2)*(-2. - 4.*C4*s)))*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + t2)))/pow(mrho,2)
               - 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + t2)) -
             (2.*(1.*pow(2. - 1.*delta,2)*pow(mpion,4)*pow(mrho,2) + 0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,4)*s +
                  pow(mpion,2)*(C4*(4. - 2.*delta)*pow(mrho,6) + (-1. + 0.5*delta)*delta*pow(s,2) +
                     pow(mrho,2)*s*(-1. + 3.*delta - 1.25*pow(delta,2) + 4.*C4*s - 2.*C4*delta*s) +
                     pow(mrho,4)*(-2. + 1.*delta - 8.*C4*s + 4.*C4*delta*s)))*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + t2)))/
              (pow(mrho,4) - 1.*pow(mrho,2)*s) + (0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*
                (eta2*pow(mpion,2)*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + pow(ma1,2)*(-1.*pow(mrho,2) + 1.*s)) +
                  eta1*(pow(ma1,2)*(-0.5*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s) +
                     pow(mpion,2)*(0.5*pow(mrho,4) - 0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*
                log(fabs(-pow(mpion,2) - pow(mrho,2) + s + t2)))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) +
             (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8) + 0.5*pow(mpion,8) +
                     0.5*pow(ma1,4)*pow(mpion,2)*pow(mrho,2) - 0.5*pow(mpion,6)*pow(mrho,2) +
                     pow(Gammaa1,2)*pow(ma1,2)*(pow(mpion,2)*(1.*pow(mpion,2) + 1.5*pow(mrho,2) - 2.*s) +
                        pow(ma1,2)*(-1.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s)) +
                     pow(ma1,6)*(1.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) +
                     pow(ma1,2)*pow(mpion,4)*(-1.*pow(mpion,2) - 0.5*pow(mrho,2) + 1.*s)) +
                  eta1*(pow(ma1,6)*(1.*pow(mpion,2) + 0.5*s) +
                     pow(ma1,2)*(3.*pow(mpion,6) + 1.*pow(mpion,2)*pow(mrho,4) - 0.5*pow(mpion,4)*s) +
                     pow(ma1,4)*(-3.*pow(mpion,4) + pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
                     pow(mpion,4)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 0.5*s) + 0.5*pow(mrho,2)*s) +
                     pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) + pow(ma1,2)*(1.*pow(mpion,2) + 0.5*s) - 0.5*pow(mrho,2)*s +
                        pow(mpion,2)*(-1.*pow(mrho,2) + 1.5*s))))*
                log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*t2 - 2.*pow(mrho,2)*t2 + 2.*s*t2 + pow(t2,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*t2))))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) -
             (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(-1.*pow(mpion,8) + 1.*pow(mpion,4)*pow(mrho,4) +
                     pow(ma1,6)*(1.*pow(mpion,2) - 0.5*s) - 0.5*pow(mpion,6)*s - 4.*pow(mpion,4)*pow(mrho,2)*s -
                     1.5*pow(mpion,2)*pow(mrho,4)*s + 3.5*pow(mpion,4)*pow(s,2) + 4.*pow(mpion,2)*pow(mrho,2)*pow(s,2) +
                     0.5*pow(mrho,4)*pow(s,2) - 2.5*pow(mpion,2)*pow(s,3) - 1.*pow(mrho,2)*pow(s,3) + 0.5*pow(s,4) +
                     pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) + pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) - 0.5*pow(mpion,2)*s +
                        0.5*pow(s,2)) + pow(ma1,4)*(-3.*pow(mpion,4) + (1.*pow(mrho,2) - 0.5*s)*s +
                        pow(mpion,2)*(-2.*pow(mrho,2) + 2.5*s)) +
                     pow(ma1,2)*(3.*pow(mpion,6) + pow(mpion,4)*(2.*pow(mrho,2) - 1.5*s) - 0.5*pow(mrho,4)*s + 0.5*pow(s,3) +
                        pow(mpion,2)*(1.*pow(mrho,4) - 1.*pow(mrho,2)*s - 1.*pow(s,2)))) +
                  eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8) + pow(ma1,6)*(1.*pow(mpion,2) + 1.*pow(mrho,2) - 0.5*s) +
                     pow(ma1,4)*(-0.5*pow(mrho,4) + (-0.5*pow(mpion,2) + 0.5*s)*s) +
                     pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) -
                        1.*pow(mrho,2)*s + 0.5*pow(s,2) + pow(ma1,2)*(-1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.5*s)) +
                     pow(mpion,2)*(0.5*pow(mpion,6) + 0.5*pow(mpion,4)*s +
                        pow(mpion,2)*(-0.5*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.5*pow(s,2)) +
                        s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2))) +
                     pow(ma1,2)*(-1.*pow(mpion,6) + pow(mpion,4)*(-1.*pow(mrho,2) + 0.5*s) +
                        pow(mpion,2)*(-1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
                        s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2)))))*
                log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*t2 - 2.*pow(mrho,2)*t2 + 2.*s*t2 + pow(t2,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*t2))))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
                2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) -
             (0.0625*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(-1.*pow(ma1,10) + pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
                     pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
                     pow(ma1,4)*(10.*pow(mpion,6) - 2.5*pow(mrho,6) + pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) + 6.*pow(mrho,4)*s -
                        4.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3)) +
                     pow(ma1,6)*(-10.*pow(mpion,4) + (1.*pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-10.*pow(mrho,2) + 8.*s)) +
                     pow(mpion,4)*(1.*pow(mpion,6) + 0.5*pow(mrho,6) + pow(mpion,4)*(2.5*pow(mrho,2) - 0.5*s) - 1.*pow(mrho,4)*s +
                        0.5*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(2.*pow(mrho,4) - 2.*pow(mrho,2)*s)) +
                     pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) - 4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
                        1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
                        pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
                        pow(mpion,2)*(-3.*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.*pow(s,2)) +
                        pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4) + pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
                           3.*pow(s,2))) + pow(ma1,2)*(-5.*pow(mpion,8) + 1.*pow(mrho,8) - 3.5*pow(mrho,6)*s +
                        4.5*pow(mrho,4)*pow(s,2) - 2.5*pow(mrho,2)*pow(s,3) + 0.5*pow(s,4) + pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
                        pow(mpion,4)*(-2.*pow(mrho,4) + 1.*pow(mrho,2)*s + 1.*pow(s,2)) +
                        pow(mpion,2)*(3.*pow(mrho,6) - 8.*pow(mrho,4)*s + 7.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)))) +
                  pow(eta1,2)*(-1.*pow(ma1,10) + pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
                     pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
                     pow(ma1,6)*(-10.*pow(mpion,4) - 2.*pow(mrho,4) + 5.*pow(mrho,2)*s - 1.*pow(s,2) +
                        pow(mpion,2)*(-10.*pow(mrho,2) + 8.*s)) +
                     pow(ma1,4)*(10.*pow(mpion,6) + 0.5*pow(mrho,6) + pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) - 3.*pow(mrho,4)*s +
                        1.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s)) +
                     pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) - 4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
                        1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
                        pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
                        pow(mpion,2)*(-3.*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.*pow(s,2)) +
                        pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4) + pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
                           3.*pow(s,2))) + pow(mpion,2)*(1.*pow(mpion,8) + pow(mpion,6)*(2.5*pow(mrho,2) - 0.5*s) +
                        pow(mpion,4)*(4.*pow(mrho,4) - 6.*pow(mrho,2)*s) +
                        pow(mrho,2)*s*(-1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
                        pow(mpion,2)*(1.5*pow(mrho,6) - 6.*pow(mrho,4)*s + 4.5*pow(mrho,2)*pow(s,2))) +
                     pow(ma1,2)*(-5.*pow(mpion,8) + pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
                        pow(mpion,4)*(-8.*pow(mrho,4) + 13.*pow(mrho,2)*s + 1.*pow(s,2)) +
                        pow(mpion,2)*(-1.*pow(mrho,6) + 6.*pow(mrho,4)*s - 3.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)) +
                        s*(0.5*pow(mrho,6) - 0.5*pow(mrho,4)*s - 0.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3)))) +
                  eta1*eta2*(2.*pow(ma1,10) + pow(Gammaa1,4)*pow(ma1,4)*(-2.*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s) +
                     pow(ma1,8)*(-10.*pow(mpion,2) - 5.*pow(mrho,2) + 5.*s) +
                     pow(ma1,6)*(20.*pow(mpion,4) + 6.*pow(mrho,4) + pow(mpion,2)*(20.*pow(mrho,2) - 16.*s) - 8.*pow(mrho,2)*s +
                        2.*pow(s,2)) + pow(ma1,4)*(-20.*pow(mpion,6) - 4.*pow(mrho,6) + 6.*pow(mrho,4)*s - 2.*pow(s,3) +
                        pow(mpion,4)*(-30.*pow(mrho,2) + 18.*s) + pow(mpion,2)*(-18.*pow(mrho,4) + 18.*pow(mrho,2)*s)) +
                     pow(mpion,2)*(-2.*pow(mpion,8) - 1.*pow(mrho,8) + 3.*pow(mrho,6)*s - 3.*pow(mrho,4)*pow(s,2) +
                        1.*pow(mrho,2)*pow(s,3) + pow(mpion,6)*(-5.*pow(mrho,2) + 1.*s) +
                        pow(mpion,4)*(-10.*pow(mrho,4) + 10.*pow(mrho,2)*s) +
                        pow(mpion,2)*(-6.*pow(mrho,6) + 12.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2))) +
                     pow(ma1,2)*(10.*pow(mpion,8) + 1.*pow(mrho,8) + pow(mpion,6)*(20.*pow(mrho,2) - 8.*s) - 2.*pow(mrho,6)*s +
                        2.*pow(mrho,2)*pow(s,3) - 1.*pow(s,4) + pow(mpion,4)*(22.*pow(mrho,4) - 20.*pow(mrho,2)*s - 2.*pow(s,2)) +
                        pow(mpion,2)*(8.*pow(mrho,6) - 12.*pow(mrho,4)*s + 4.*pow(s,3))) +
                     pow(Gammaa1,2)*pow(ma1,2)*(-8.*pow(ma1,6) + 8.*pow(mpion,6) + 1.*pow(mrho,6) +
                        pow(mpion,4)*(12.*pow(mrho,2) - 12.*s) + pow(ma1,4)*(24.*pow(mpion,2) + 12.*pow(mrho,2) - 12.*s) -
                        3.*pow(mrho,4)*s + 3.*pow(mrho,2)*pow(s,2) - 1.*pow(s,3) +
                        pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s + 6.*pow(s,2)) +
                        pow(ma1,2)*(-24.*pow(mpion,4) - 6.*pow(mrho,4) + 12.*pow(mrho,2)*s - 6.*pow(s,2) +
                           pow(mpion,2)*(-24.*pow(mrho,2) + 24.*s)))))*
                log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*t2 - 2.*pow(mrho,2)*t2 + 2.*s*t2 + pow(t2,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*t2))))/
              (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
                4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
             0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(4.*pow(ma1,6) +
                   pow(Gammaa1,2)*pow(ma1,2)*(-4.*pow(ma1,2) + 4.*pow(mpion,2) - 2.*s) + pow(ma1,4)*(-12.*pow(mpion,2) + 6.*s) +
                   pow(mpion,2)*(-4.*pow(mpion,4) + 4.*pow(mrho,4) + 2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s) +
                   pow(ma1,2)*(12.*pow(mpion,4) - 2.*pow(mrho,4) - 8.*pow(mpion,2)*s - 4.*pow(mrho,2)*s + 4.*pow(s,2))) +
                pow(eta1,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) + 3.*pow(mpion,4)*pow(mrho,2) +
                   pow(ma1,4)*(6.*pow(mpion,2) + 3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s - 1.*pow(mpion,2)*pow(mrho,2)*s -
                   1.*pow(mrho,4)*s + pow(mrho,2)*pow(s,2) +
                   pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s) +
                   pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) + 4.*pow(mrho,2)*s - 2.*pow(s,2) +
                      pow(mpion,2)*(-6.*pow(mrho,2) + 4.*s))) +
                pow(eta2,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) - 3.*pow(mpion,4)*pow(mrho,2) +
                   pow(ma1,4)*(6.*pow(mpion,2) - 3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s + pow(mpion,2)*pow(mrho,2)*s +
                   pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) - 2.*pow(mpion,2) + pow(mrho,2) + s) +
                   pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 2.*pow(s,2) +
                      pow(mpion,2)*(6.*pow(mrho,2) + 4.*s))))*
              log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + 1.*pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*t2 - 2.*pow(mrho,2)*t2 + 2.*s*t2 + pow(t2,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*t2))) -
             (0.5*(1.*eta1 - 1.*eta2)*(eta2*(pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(0.5 - 1.*C4*pow(mrho,2)) +
                     pow(ma1,4)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                     pow(mpion,2)*pow(mrho,2)*(pow(mpion,2)*(-0.5 + 1.*C4*pow(mrho,2)) + (0.25 - 0.125*delta)*(pow(mrho,2) + s)) +
                     pow(ma1,2)*(1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) - 0.25*delta*pow(s,2) +
                        pow(mrho,4)*(-0.75 + 0.125*delta - 2.*C4*s) + pow(mrho,2)*s*(0.25 + 0.375*delta + 1.*C4*s))) +
                  eta1*(pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(-0.5 + 1.*C4*pow(mrho,2)) +
                     pow(ma1,4)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
                     pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6) + pow(mpion,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) +
                        0.25*delta*pow(s,2) - 1.*C4*pow(mrho,2)*pow(s,2)) +
                     pow(mrho,2)*(pow(mpion,4)*(0.5 - 1.*C4*pow(mrho,2)) + s*((-0.25 + 0.125*delta)*pow(mrho,2) + (0.25 - 0.125*delta)*s) +
                        pow(mpion,2)*(-2.*C4*pow(mrho,4) + (-0.5 - 0.25*delta)*s + pow(mrho,2)*(1. + 2.*C4*s)))))*
                log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + 1.*pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*t2 - 2.*pow(mrho,2)*t2 + 2.*s*t2 + pow(t2,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*t2))))/pow(mrho,2) -
             0.0625*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
                   pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) - 2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
                   pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
                pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) + pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s +
                   pow(mpion,4)*(-3.*pow(mrho,2) + s) + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
                   pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2)))\
                 + pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
                   pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
                   pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2)))
                )*log(fabs(-pow(ma1,2) + t1)) + (0.03125*pow(eta1 - 1.*eta2,2)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s)*
                (eta1*eta2*(-2.*pow(ma1,8) + pow(ma1,6)*(8.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) +
                     pow(ma1,4)*(-12.*pow(mpion,4) - 4.*pow(mrho,4) + 4.*pow(mrho,2)*s + pow(mpion,2)*(-12.*pow(mrho,2) + 8.*s)) +
                     pow(mpion,2)*(-2.*pow(mpion,6) - 4.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 4.*pow(mrho,4)*s -
                        2.*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(-8.*pow(mrho,4) + 8.*pow(mrho,2)*s)) +
                     pow(ma1,2)*(8.*pow(mpion,6) + 2.*pow(mrho,6) + pow(mpion,4)*(12.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,4)*s -
                        2.*pow(mrho,2)*pow(s,2) + 2.*pow(s,3) + pow(mpion,2)*(8.*pow(mrho,4) - 4.*pow(mrho,2)*s - 4.*pow(s,2)))) +
                  pow(eta2,2)*(pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                     pow(mpion,4)*(pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 1.*pow(mrho,2)*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) + pow(mrho,2)*s) +
                     pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mrho,6) - 5.*pow(mrho,4)*s + 4.*pow(mrho,2)*pow(s,2) - 1.*pow(s,3) +
                        pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) + pow(mpion,2)*(2.*pow(mrho,4) - 4.*pow(mrho,2)*s + 2.*pow(s,2)))) +
                  pow(eta1,2)*(pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                     pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 3.*pow(mrho,2)*s) +
                     pow(mpion,2)*(pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) + pow(mrho,2)*s*(-2.*pow(mrho,2) + 2.*s) +
                        pow(mpion,2)*(3.*pow(mrho,4) - 5.*pow(mrho,2)*s)) +
                     pow(ma1,2)*(-4.*pow(mpion,6) + pow(mrho,4)*s - 1.*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) +
                        pow(mpion,2)*(-2.*pow(mrho,4) + 4.*pow(mrho,2)*s + 2.*pow(s,2)))))*log(fabs(-pow(ma1,2) + t1)))/
              (0.25*pow(Gammaa1,2)*pow(ma1,2) + 1.*pow(ma1,4) + 1.*pow(mpion,4) + 1.*pow(mpion,2)*pow(mrho,2) + 0.25*pow(mrho,4) -
                1.*pow(mpion,2)*s - 0.5*pow(mrho,2)*s + 0.25*pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s)) +
             (1.*(eta1*eta2*(pow(ma1,8)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                     pow(ma1,6)*(-1.*pow(mrho,4) + 2.*C4*pow(mrho,6) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*C4*pow(mrho,4)) +
                        0.5*delta*pow(s,2) + pow(mrho,2)*s*(2. - 1.*delta - 2.*C4*s)) +
                     pow(ma1,2)*(pow(mpion,6)*(-4.*pow(mrho,2) + 8.*C4*pow(mrho,4)) +
                        pow(mrho,4)*s*((0.5 - 0.25*delta)*pow(mrho,2) + (-0.5 + 0.25*delta)*s) +
                        pow(mpion,4)*(10.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,4)*(-3. - 1.*delta - 8.*C4*s) +
                           pow(mrho,2)*s*(2. + 1.*delta - 2.*C4*s)) +
                        pow(mpion,2)*(2.*C4*pow(mrho,8) - 0.5*delta*pow(s,3) + pow(mrho,6)*(1. - 1.*delta - 2.*C4*s) +
                           pow(mrho,4)*s*(-1. + 1.*delta - 2.*C4*s) + pow(mrho,2)*pow(s,2)*(1. + 2.*C4*s))) +
                     pow(mpion,2)*pow(mrho,2)*(pow(mpion,6)*(1. - 2.*C4*pow(mrho,2)) +
                        pow(mrho,4)*((-0.5 + 0.25*delta)*pow(mrho,2) + (0.5 - 0.25*delta)*s) +
                        pow(mpion,2)*(-2.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,2)*s*(-1.5 - 0.25*delta - 2.*C4*s) +
                           pow(mrho,4)*(1. + 4.*C4*s)) + pow(mpion,4)*
                         (-4.*C4*pow(mrho,4) - 1.*delta*s + pow(mrho,2)*(1. + 0.5*delta + 4.*C4*s))) +
                     pow(ma1,4)*(pow(mpion,4)*(6.*pow(mrho,2) - 12.*C4*pow(mrho,4)) +
                        pow(mpion,2)*(-8.*C4*pow(mrho,6) - 1.*delta*pow(s,2) + pow(mrho,4)*(3. + 0.5*delta + 4.*C4*s) +
                           pow(mrho,2)*s*(-4. + 1.*delta + 4.*C4*s)) +
                        s*(-2.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,2)*s*(1. - 1.5*delta - 2.*C4*s) +
                           pow(mrho,4)*(-1.5 + 1.25*delta + 4.*C4*s)))) +
                  pow(eta1,2)*(pow(ma1,8)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                     pow(ma1,6)*(-2.*C4*pow(mrho,6) + pow(mpion,2)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) - 0.25*delta*pow(s,2) +
                        pow(mrho,4)*(1. + 1.*C4*s) + pow(mrho,2)*s*(-1. + 0.25*delta + 1.*C4*s)) +
                     pow(ma1,4)*(1.*C4*pow(mrho,8) + pow(mpion,4)*(-3.*pow(mrho,2) + 6.*C4*pow(mrho,4)) - 0.25*delta*pow(s,3) +
                        pow(mrho,6)*(-0.5 - 1.*C4*s) + pow(mrho,4)*s*(1.5 - 0.5*delta - 1.*C4*s) +
                        pow(mrho,2)*pow(s,2)*(-0.5 + 0.5*delta + 1.*C4*s) +
                        pow(mpion,2)*(7.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,4)*(-3. - 0.25*delta - 5.*C4*s) +
                           pow(mrho,2)*s*(2. + 0.25*delta - 2.*C4*s))) +
                     pow(mpion,2)*pow(mrho,2)*(pow(mpion,6)*(-0.5 + 1.*C4*pow(mrho,2)) +
                        pow(mrho,4)*((0.5 - 0.25*delta)*pow(mrho,2) + (-0.5 + 0.25*delta)*s) +
                        pow(mpion,4)*(3.*C4*pow(mrho,4) + 0.75*delta*s + pow(mrho,2)*(-1. - 0.25*delta - 3.*C4*s)) +
                        pow(mpion,2)*(2.*C4*pow(mrho,6) - 0.5*delta*pow(s,2) + pow(mrho,4)*(-1. - 4.*C4*s) +
                           pow(mrho,2)*s*(1.5 + 0.25*delta + 2.*C4*s))) +
                     pow(ma1,2)*(pow(mpion,6)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                        pow(mrho,4)*s*((-0.5 + 0.25*delta)*pow(mrho,2) + (0.5 - 0.25*delta)*s) +
                        pow(mpion,2)*(-3.*C4*pow(mrho,8) + 0.25*delta*pow(s,3) + pow(mrho,4)*s*(-1. - 0.75*delta - 1.*C4*s) +
                           pow(mrho,2)*pow(s,2)*(-0.5 + 0.5*delta - 1.*C4*s) + pow(mrho,6)*(0.5 + 0.5*delta + 5.*C4*s)) +
                        pow(mpion,4)*(-8.*C4*pow(mrho,6) - 0.25*delta*pow(s,2) + pow(mrho,2)*s*(-1. - 1.25*delta + 1.*C4*s) +
                           pow(mrho,4)*(3. + 0.5*delta + 7.*C4*s)))) +
                  pow(eta2,2)*(pow(ma1,8)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                     pow(mpion,6)*pow(mrho,2)*(-0.25*delta*pow(mrho,2) + 1.*C4*pow(mrho,4) + pow(mpion,2)*(-0.5 + 1.*C4*pow(mrho,2)) +
                        0.25*delta*s - 1.*C4*pow(mrho,2)*s) +
                     pow(ma1,6)*(pow(mpion,2)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                        s*(-1.*C4*pow(mrho,4) - 0.25*delta*s + pow(mrho,2)*(-1. + 0.75*delta + 1.*C4*s))) +
                     pow(ma1,4)*(-1.*C4*pow(mrho,8) + pow(mpion,4)*(-3.*pow(mrho,2) + 6.*C4*pow(mrho,4)) - 0.25*delta*pow(s,3) +
                        pow(mrho,4)*s*(-0.75*delta - 3.*C4*s) + pow(mrho,2)*pow(s,2)*(-0.5 + 1.*delta + 1.*C4*s) +
                        pow(mrho,6)*(0.5 + 3.*C4*s) + pow(mpion,2)*
                         (delta*(-0.25*pow(mrho,4) - 1.25*pow(mrho,2)*s + 0.5*pow(s,2)) +
                           pow(mrho,2)*(1.*C4*pow(mrho,4) + 2.*s + 1.*C4*pow(mrho,2)*s - 2.*C4*pow(s,2)))) +
                     pow(ma1,2)*pow(mpion,2)*(1.*C4*pow(mrho,8) + pow(mpion,4)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                        0.25*delta*pow(s,3) + pow(mrho,6)*(-1.5 + 0.5*delta - 3.*C4*s) + pow(mrho,2)*pow(s,2)*(-0.5 - 0.5*delta - 1.*C4*s) +
                        pow(mrho,4)*s*(2. - 0.25*delta + 3.*C4*s) +
                        pow(mpion,2)*(delta*(0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
                           pow(mrho,2)*(-2.*C4*pow(mrho,4) - 1.*s + 1.*C4*pow(mrho,2)*s + 1.*C4*pow(s,2))))))*log(fabs(-pow(ma1,2) + t1)))/
              ((pow(ma1,2) - 1.*pow(mpion,2))*pow(mrho,2)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s)) -
             (0.5*pow(mpion,2)*(pow(eta2,2)*pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s +
                     (-2. + 1.*delta)*pow(s,2)) + eta1*eta2*(pow(mpion,2)*
                      ((4. - 2.*delta)*pow(mrho,4) + (-8. + 4.*delta)*pow(mrho,2)*s + (4. - 2.*delta)*pow(s,2)) +
                     pow(mrho,2)*((-1. + 0.5*delta)*pow(mrho,4) + (2. - 1.*delta)*pow(mrho,2)*s + (-1. + 0.5*delta)*pow(s,2))) +
                  pow(eta1,2)*(pow(mrho,2)*((1. - 0.5*delta)*pow(mrho,4) + (-2. + 1.*delta)*pow(mrho,2)*s + (1. - 0.5*delta)*pow(s,2)) +
                     pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))))*
                log(fabs(-pow(mpion,2) + t1)))/((-1.*pow(ma1,2) + 1.*pow(mpion,2))*(pow(mrho,2) - 1.*s)) +
             (0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*
                (eta1*(pow(mpion,4)*(-1.*pow(mrho,2) + 1.*s) +
                     pow(ma1,2)*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + (-0.5*pow(mrho,2) + 0.5*s)*s) +
                     pow(mpion,2)*(-1.*pow(mrho,4) + 2.5*pow(mrho,2)*s - 1.5*pow(s,2)) +
                     s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2))) +
                  eta2*(0.5*pow(mrho,6) + pow(mpion,4)*(1.*pow(mrho,2) - 1.*s) - 1.5*pow(mrho,4)*s + 1.5*pow(mrho,2)*pow(s,2) -
                     0.5*pow(s,3) + pow(mpion,2)*(1.5*pow(mrho,4) - 3.*pow(mrho,2)*s + 1.5*pow(s,2)) +
                     pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*pow(mrho,2)*s - 0.5*pow(s,2) + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*
                log(fabs(-pow(mpion,2) + t1)))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
                2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
             (0.5*pow(mpion,2)*(eta1*eta2*((1. - 0.5*delta)*pow(mrho,6) + (-4. + 2.*delta)*pow(mrho,4)*s +
                     (5. - 2.5*delta)*pow(mrho,2)*pow(s,2) + (-2. + 1.*delta)*pow(s,3) +
                     pow(mpion,2)*((4. - 2.*delta)*pow(mrho,4) + (-8. + 4.*delta)*pow(mrho,2)*s + (4. - 2.*delta)*pow(s,2))) +
                  pow(eta2,2)*((-1. + 0.5*delta)*pow(mrho,6) + (3. - 1.5*delta)*pow(mrho,4)*s + (-3. + 1.5*delta)*pow(mrho,2)*pow(s,2) +
                     (1. - 0.5*delta)*pow(s,3) + pow(mpion,2)*
                      ((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))) +
                  pow(eta1,2)*(s*((1. - 0.5*delta)*pow(mrho,4) + (-2. + 1.*delta)*pow(mrho,2)*s + (1. - 0.5*delta)*pow(s,2)) +
                     pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))))*
                log(fabs(-pow(mpion,2) - pow(mrho,2) + s + t1)))/
              ((pow(mrho,2) - 1.*s)*(-1.*pow(ma1,2) + pow(mpion,2) + pow(mrho,2) - 1.*s)) -
             0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) + t1)) -
             (0.25000000000000006*pow(2. - 1.*delta,2)*(7.999999999999998*pow(mpion,4) - 5.999999999999998*pow(mpion,2)*s + 1.*pow(s,2))*
                log(fabs(-pow(mpion,2) + t1)))/(pow(mrho,2) - 1.*s) +
             (1.*(-2. + 1.*delta)*((0.5 - 0.25*delta)*pow(mrho,2)*s +
                  pow(mpion,2)*(4.*C4*pow(mrho,4) + 1.*delta*s + pow(mrho,2)*(-2. - 4.*C4*s)))*log(fabs(-pow(mpion,2) + t1)))/pow(mrho,2)
               + 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + t1)) +
             (2.*(1.*pow(2. - 1.*delta,2)*pow(mpion,4)*pow(mrho,2) + 0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,4)*s +
                  pow(mpion,2)*(C4*(4. - 2.*delta)*pow(mrho,6) + (-1. + 0.5*delta)*delta*pow(s,2) +
                     pow(mrho,2)*s*(-1. + 3.*delta - 1.25*pow(delta,2) + 4.*C4*s - 2.*C4*delta*s) +
                     pow(mrho,4)*(-2. + 1.*delta - 8.*C4*s + 4.*C4*delta*s)))*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + t1)))/
              (pow(mrho,4) - 1.*pow(mrho,2)*s) - (0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*
                (eta2*pow(mpion,2)*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + pow(ma1,2)*(-1.*pow(mrho,2) + 1.*s)) +
                  eta1*(pow(ma1,2)*(-0.5*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s) +
                     pow(mpion,2)*(0.5*pow(mrho,4) - 0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*
                log(fabs(-pow(mpion,2) - pow(mrho,2) + s + t1)))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) -
             (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8) + 0.5*pow(mpion,8) +
                     0.5*pow(ma1,4)*pow(mpion,2)*pow(mrho,2) - 0.5*pow(mpion,6)*pow(mrho,2) +
                     pow(Gammaa1,2)*pow(ma1,2)*(pow(mpion,2)*(1.*pow(mpion,2) + 1.5*pow(mrho,2) - 2.*s) +
                        pow(ma1,2)*(-1.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s)) +
                     pow(ma1,6)*(1.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) +
                     pow(ma1,2)*pow(mpion,4)*(-1.*pow(mpion,2) - 0.5*pow(mrho,2) + 1.*s)) +
                  eta1*(pow(ma1,6)*(1.*pow(mpion,2) + 0.5*s) +
                     pow(ma1,2)*(3.*pow(mpion,6) + 1.*pow(mpion,2)*pow(mrho,4) - 0.5*pow(mpion,4)*s) +
                     pow(ma1,4)*(-3.*pow(mpion,4) + pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
                     pow(mpion,4)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 0.5*s) + 0.5*pow(mrho,2)*s) +
                     pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) + pow(ma1,2)*(1.*pow(mpion,2) + 0.5*s) - 0.5*pow(mrho,2)*s +
                        pow(mpion,2)*(-1.*pow(mrho,2) + 1.5*s))))*
                log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*t1 - 2.*pow(mrho,2)*t1 + 2.*s*t1 + pow(t1,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*t1))))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) +
             (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(-1.*pow(mpion,8) + 1.*pow(mpion,4)*pow(mrho,4) +
                     pow(ma1,6)*(1.*pow(mpion,2) - 0.5*s) - 0.5*pow(mpion,6)*s - 4.*pow(mpion,4)*pow(mrho,2)*s -
                     1.5*pow(mpion,2)*pow(mrho,4)*s + 3.5*pow(mpion,4)*pow(s,2) + 4.*pow(mpion,2)*pow(mrho,2)*pow(s,2) +
                     0.5*pow(mrho,4)*pow(s,2) - 2.5*pow(mpion,2)*pow(s,3) - 1.*pow(mrho,2)*pow(s,3) + 0.5*pow(s,4) +
                     pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) + pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) - 0.5*pow(mpion,2)*s +
                        0.5*pow(s,2)) + pow(ma1,4)*(-3.*pow(mpion,4) + (1.*pow(mrho,2) - 0.5*s)*s +
                        pow(mpion,2)*(-2.*pow(mrho,2) + 2.5*s)) +
                     pow(ma1,2)*(3.*pow(mpion,6) + pow(mpion,4)*(2.*pow(mrho,2) - 1.5*s) - 0.5*pow(mrho,4)*s + 0.5*pow(s,3) +
                        pow(mpion,2)*(1.*pow(mrho,4) - 1.*pow(mrho,2)*s - 1.*pow(s,2)))) +
                  eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8) + pow(ma1,6)*(1.*pow(mpion,2) + 1.*pow(mrho,2) - 0.5*s) +
                     pow(ma1,4)*(-0.5*pow(mrho,4) + (-0.5*pow(mpion,2) + 0.5*s)*s) +
                     pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) -
                        1.*pow(mrho,2)*s + 0.5*pow(s,2) + pow(ma1,2)*(-1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.5*s)) +
                     pow(mpion,2)*(0.5*pow(mpion,6) + 0.5*pow(mpion,4)*s +
                        pow(mpion,2)*(-0.5*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.5*pow(s,2)) +
                        s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2))) +
                     pow(ma1,2)*(-1.*pow(mpion,6) + pow(mpion,4)*(-1.*pow(mrho,2) + 0.5*s) +
                        pow(mpion,2)*(-1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
                        s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2)))))*
                log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*t1 - 2.*pow(mrho,2)*t1 + 2.*s*t1 + pow(t1,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*t1))))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
                2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
             (0.0625*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(-1.*pow(ma1,10) + pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
                     pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
                     pow(ma1,4)*(10.*pow(mpion,6) - 2.5*pow(mrho,6) + pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) + 6.*pow(mrho,4)*s -
                        4.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3)) +
                     pow(ma1,6)*(-10.*pow(mpion,4) + (1.*pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-10.*pow(mrho,2) + 8.*s)) +
                     pow(mpion,4)*(1.*pow(mpion,6) + 0.5*pow(mrho,6) + pow(mpion,4)*(2.5*pow(mrho,2) - 0.5*s) - 1.*pow(mrho,4)*s +
                        0.5*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(2.*pow(mrho,4) - 2.*pow(mrho,2)*s)) +
                     pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) - 4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
                        1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
                        pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
                        pow(mpion,2)*(-3.*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.*pow(s,2)) +
                        pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4) + pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
                           3.*pow(s,2))) + pow(ma1,2)*(-5.*pow(mpion,8) + 1.*pow(mrho,8) - 3.5*pow(mrho,6)*s +
                        4.5*pow(mrho,4)*pow(s,2) - 2.5*pow(mrho,2)*pow(s,3) + 0.5*pow(s,4) + pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
                        pow(mpion,4)*(-2.*pow(mrho,4) + 1.*pow(mrho,2)*s + 1.*pow(s,2)) +
                        pow(mpion,2)*(3.*pow(mrho,6) - 8.*pow(mrho,4)*s + 7.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)))) +
                  pow(eta1,2)*(-1.*pow(ma1,10) + pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
                     pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
                     pow(ma1,6)*(-10.*pow(mpion,4) - 2.*pow(mrho,4) + 5.*pow(mrho,2)*s - 1.*pow(s,2) +
                        pow(mpion,2)*(-10.*pow(mrho,2) + 8.*s)) +
                     pow(ma1,4)*(10.*pow(mpion,6) + 0.5*pow(mrho,6) + pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) - 3.*pow(mrho,4)*s +
                        1.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s)) +
                     pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) - 4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
                        1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
                        pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
                        pow(mpion,2)*(-3.*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.*pow(s,2)) +
                        pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4) + pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
                           3.*pow(s,2))) + pow(mpion,2)*(1.*pow(mpion,8) + pow(mpion,6)*(2.5*pow(mrho,2) - 0.5*s) +
                        pow(mpion,4)*(4.*pow(mrho,4) - 6.*pow(mrho,2)*s) +
                        pow(mrho,2)*s*(-1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
                        pow(mpion,2)*(1.5*pow(mrho,6) - 6.*pow(mrho,4)*s + 4.5*pow(mrho,2)*pow(s,2))) +
                     pow(ma1,2)*(-5.*pow(mpion,8) + pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
                        pow(mpion,4)*(-8.*pow(mrho,4) + 13.*pow(mrho,2)*s + 1.*pow(s,2)) +
                        pow(mpion,2)*(-1.*pow(mrho,6) + 6.*pow(mrho,4)*s - 3.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)) +
                        s*(0.5*pow(mrho,6) - 0.5*pow(mrho,4)*s - 0.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3)))) +
                  eta1*eta2*(2.*pow(ma1,10) + pow(Gammaa1,4)*pow(ma1,4)*(-2.*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s) +
                     pow(ma1,8)*(-10.*pow(mpion,2) - 5.*pow(mrho,2) + 5.*s) +
                     pow(ma1,6)*(20.*pow(mpion,4) + 6.*pow(mrho,4) + pow(mpion,2)*(20.*pow(mrho,2) - 16.*s) - 8.*pow(mrho,2)*s +
                        2.*pow(s,2)) + pow(ma1,4)*(-20.*pow(mpion,6) - 4.*pow(mrho,6) + 6.*pow(mrho,4)*s - 2.*pow(s,3) +
                        pow(mpion,4)*(-30.*pow(mrho,2) + 18.*s) + pow(mpion,2)*(-18.*pow(mrho,4) + 18.*pow(mrho,2)*s)) +
                     pow(mpion,2)*(-2.*pow(mpion,8) - 1.*pow(mrho,8) + 3.*pow(mrho,6)*s - 3.*pow(mrho,4)*pow(s,2) +
                        1.*pow(mrho,2)*pow(s,3) + pow(mpion,6)*(-5.*pow(mrho,2) + 1.*s) +
                        pow(mpion,4)*(-10.*pow(mrho,4) + 10.*pow(mrho,2)*s) +
                        pow(mpion,2)*(-6.*pow(mrho,6) + 12.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2))) +
                     pow(ma1,2)*(10.*pow(mpion,8) + 1.*pow(mrho,8) + pow(mpion,6)*(20.*pow(mrho,2) - 8.*s) - 2.*pow(mrho,6)*s +
                        2.*pow(mrho,2)*pow(s,3) - 1.*pow(s,4) + pow(mpion,4)*(22.*pow(mrho,4) - 20.*pow(mrho,2)*s - 2.*pow(s,2)) +
                        pow(mpion,2)*(8.*pow(mrho,6) - 12.*pow(mrho,4)*s + 4.*pow(s,3))) +
                     pow(Gammaa1,2)*pow(ma1,2)*(-8.*pow(ma1,6) + 8.*pow(mpion,6) + 1.*pow(mrho,6) +
                        pow(mpion,4)*(12.*pow(mrho,2) - 12.*s) + pow(ma1,4)*(24.*pow(mpion,2) + 12.*pow(mrho,2) - 12.*s) -
                        3.*pow(mrho,4)*s + 3.*pow(mrho,2)*pow(s,2) - 1.*pow(s,3) +
                        pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s + 6.*pow(s,2)) +
                        pow(ma1,2)*(-24.*pow(mpion,4) - 6.*pow(mrho,4) + 12.*pow(mrho,2)*s - 6.*pow(s,2) +
                           pow(mpion,2)*(-24.*pow(mrho,2) + 24.*s)))))*
                log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*t1 - 2.*pow(mrho,2)*t1 + 2.*s*t1 + pow(t1,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*t1))))/
              (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
                4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) -
             0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(4.*pow(ma1,6) +
                   pow(Gammaa1,2)*pow(ma1,2)*(-4.*pow(ma1,2) + 4.*pow(mpion,2) - 2.*s) + pow(ma1,4)*(-12.*pow(mpion,2) + 6.*s) +
                   pow(mpion,2)*(-4.*pow(mpion,4) + 4.*pow(mrho,4) + 2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s) +
                   pow(ma1,2)*(12.*pow(mpion,4) - 2.*pow(mrho,4) - 8.*pow(mpion,2)*s - 4.*pow(mrho,2)*s + 4.*pow(s,2))) +
                pow(eta1,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) + 3.*pow(mpion,4)*pow(mrho,2) +
                   pow(ma1,4)*(6.*pow(mpion,2) + 3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s - 1.*pow(mpion,2)*pow(mrho,2)*s -
                   1.*pow(mrho,4)*s + pow(mrho,2)*pow(s,2) +
                   pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s) +
                   pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) + 4.*pow(mrho,2)*s - 2.*pow(s,2) +
                      pow(mpion,2)*(-6.*pow(mrho,2) + 4.*s))) +
                pow(eta2,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) - 3.*pow(mpion,4)*pow(mrho,2) +
                   pow(ma1,4)*(6.*pow(mpion,2) - 3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s + pow(mpion,2)*pow(mrho,2)*s +
                   pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) - 2.*pow(mpion,2) + pow(mrho,2) + s) +
                   pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 2.*pow(s,2) +
                      pow(mpion,2)*(6.*pow(mrho,2) + 4.*s))))*
              log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + 1.*pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*t1 - 2.*pow(mrho,2)*t1 + 2.*s*t1 + pow(t1,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*t1))) +
             (0.5*(1.*eta1 - 1.*eta2)*(eta2*(pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(0.5 - 1.*C4*pow(mrho,2)) +
                     pow(ma1,4)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                     pow(mpion,2)*pow(mrho,2)*(pow(mpion,2)*(-0.5 + 1.*C4*pow(mrho,2)) + (0.25 - 0.125*delta)*(pow(mrho,2) + s)) +
                     pow(ma1,2)*(1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) - 0.25*delta*pow(s,2) +
                        pow(mrho,4)*(-0.75 + 0.125*delta - 2.*C4*s) + pow(mrho,2)*s*(0.25 + 0.375*delta + 1.*C4*s))) +
                  eta1*(pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(-0.5 + 1.*C4*pow(mrho,2)) +
                     pow(ma1,4)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
                     pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6) + pow(mpion,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) +
                        0.25*delta*pow(s,2) - 1.*C4*pow(mrho,2)*pow(s,2)) +
                     pow(mrho,2)*(pow(mpion,4)*(0.5 - 1.*C4*pow(mrho,2)) + s*((-0.25 + 0.125*delta)*pow(mrho,2) + (0.25 - 0.125*delta)*s) +
                        pow(mpion,2)*(-2.*C4*pow(mrho,4) + (-0.5 - 0.25*delta)*s + pow(mrho,2)*(1. + 2.*C4*s)))))*
                log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + 1.*pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*t1 - 2.*pow(mrho,2)*t1 + 2.*s*t1 + pow(t1,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*t1))))/pow(mrho,2))))/
                  (16.*M_PI*s*(-4*pow(mpion,2) + s));
            process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          }

          // the third possible reaction (produces photon)
          part_out = photon_particle;
          m3 = 0.0;

          mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
          t1 = mandelstam_t[1];
          t2 = mandelstam_t[0];

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

           xsection = to_mb*(-(pow(Const,2)*pow(ghat,4)*((0.03125*pow(eta1 - 1.*eta2,2)*
        				   (eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
        						pow(ma1,2)*pow(mpion,2)*(8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
        						pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
        					 pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
        						pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
        						pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
        						   2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(6.*pow(mrho,2) + 2.*s)))
        					   + pow(eta1,2)*(1.*pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
        						pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
        						   2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
        						   pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
        						pow(mpion,2)*(1.*pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
        						   pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/(1.*pow(ma1,2) - 1.*t2) +
        				(1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*t2) -
        				(0.25*pow(-2. + delta,2)*pow(mpion,2)*t2)/pow(mrho,2) -
        				0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,2) + pow(mrho,2) - 2.*s) + eta1*(2.*pow(mpion,2) + s))*t2 +
        				(0.5*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*(4.*pow(mpion,4)*pow(mrho,2) + 1.*pow(mrho,6) - 3.5*pow(mrho,4)*s + 0.5*pow(s,3) +
        					 pow(mpion,2)*(10.*pow(mrho,4) - 2.*pow(s,2)))*t2)/(pow(mrho,6)*pow(pow(mrho,2) - 1.*s,2)) -
        				(0.25*(eta1 - 1.*eta2)*(1.*pow(mrho,2) - 0.5*delta*s)*
        				   (eta2*(-2.*pow(ma1,4) - 6.*pow(mpion,4) + 1.5*pow(mrho,4) + pow(ma1,2)*(6.*pow(mpion,2) - 2.*pow(mrho,2) - 2.*s) -
        						1.*pow(mrho,2)*s - 0.5*pow(s,2) + pow(mpion,2)*(2.*pow(mrho,2) + 2.*s)) +
        					 eta1*(2.*pow(ma1,4) + 6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
        						1.*pow(s,2) + pow(ma1,2)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s)))*t2)/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
        				0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-6.*pow(ma1,4) - 12.*pow(mpion,4) + 2.*pow(mrho,4) +
        					  pow(ma1,2)*(16.*pow(mpion,2) - 8.*s) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
        				   pow(eta1,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
        					  2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
        				   pow(eta2,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
        					  2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) + 4.*s)))*t2 -
        				(1.*(pow(mpion,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (1.5 - 2.*delta + 0.625*pow(delta,2))*pow(mrho,2)*s +
        						(0.25 - 0.125*delta)*delta*pow(s,2) + pow(mrho,4)*(2.5 - 2.25*delta + 0.5*pow(delta,2) + 2.*C4*s - 1.*C4*delta*s)) +
        					 pow(mrho,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (0.75 - 0.375*delta)*delta*pow(s,2) +
        						pow(mrho,4)*(0.5 - 0.25*delta + 6.*C4*s - 3.*C4*delta*s) +
        						pow(mrho,2)*s*(-0.5 - 0.5*delta + 0.375*pow(delta,2) - 4.*C4*s + 2.*C4*delta*s)))*t2)/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
        				(0.25*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*(eta1*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(ma1,2)*(1. - 2.*C4*pow(mrho,2)) +
        						   pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) - 2.*C4*pow(s,2)) +
        						eta2*(-1.5*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) +
        						   pow(ma1,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s - 4.*C4*pow(mrho,2)*s + 2.*C4*pow(s,2))) +
        					 delta*(eta2*(-1.*pow(ma1,4) - 3.*pow(mpion,4) + 1.*pow(mrho,4) + pow(ma1,2)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
        						   0.25*pow(mrho,2)*s - 0.75*pow(s,2) + pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)) +
        						eta1*(1.*pow(ma1,4) + 3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s +
        						   1.*pow(s,2) + pow(ma1,2)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s))))*t2)/pow(mrho,2) +
        				(0.5*(pow(delta,2)*(1.*pow(mpion,4)*pow(mrho,2) + 0.25*pow(mrho,6) - 0.75*pow(mrho,4)*s + 0.125*pow(mrho,2)*pow(s,2) +
        						0.25*pow(s,3) + pow(mpion,2)*(2.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2))) +
        					 pow(mrho,6)*(1.5 + C4*(-6.*pow(mrho,2) + 6.*s) + pow(C4,2)*(4.*pow(mrho,4) - 8.*pow(mrho,2)*s + 4.*pow(s,2))) +
        					 delta*pow(mrho,2)*(4.*C4*pow(mrho,6) - 0.5*pow(s,2) + pow(mrho,4)*(-1.5 - 3.*C4*s) + pow(mrho,2)*s*(0.5 - 1.*C4*s) +
        						pow(mpion,2)*(6.*C4*pow(mrho,4) + 0.5*s + pow(mrho,2)*(-2.5 - 2.*C4*s))))*t2)/pow(mrho,6) -
        				(3.*(1.*pow(mrho,2) - 0.5*delta*s)*(delta*(0.666667*pow(mpion,4)*pow(mrho,2) + 0.166667*pow(mrho,6) - 0.541667*pow(mrho,4)*s -
        						0.0833333*pow(mrho,2)*pow(s,2) + 0.125*pow(s,3) +
        						pow(mpion,2)*(1.66667*pow(mrho,4) + 0.0833333*pow(mrho,2)*s - 0.416667*pow(s,2))) +
        					 pow(mrho,2)*(1.*C4*pow(mrho,6) - 0.0833333*pow(s,2) + pow(mrho,4)*(-0.416667 - 1.33333*C4*s) +
        						pow(mrho,2)*s*(0.5 + 0.333333*C4*s) + pow(mpion,2)*(2.*C4*pow(mrho,4) + 0.166667*s + pow(mrho,2)*(-0.833333 - 0.666667*C4*s))
        						))*t2)/(pow(mrho,8) - 1.*pow(mrho,6)*s) - 1.*C4*pow(t2,2) - 1.*C4*delta*pow(t2,2) +
        				0.0625*(-2. + delta)*(eta1 - 1.*eta2)*eta2*pow(t2,2) - (0.5*pow(delta,2)*pow(mpion,2)*pow(t2,2))/pow(mrho,4) +
        				(0.25*pow(t2,2))/pow(mrho,2) + (0.5*delta*pow(t2,2))/pow(mrho,2) - (0.25*pow(delta,2)*pow(t2,2))/pow(mrho,2) -
        				(0.25*delta*s*pow(t2,2))/pow(mrho,4) + (0.25*pow(delta,2)*s*pow(t2,2))/pow(mrho,4) +
        				(0.5*C4*delta*s*pow(t2,2))/pow(mrho,2) + (0.0625*pow(delta,2)*pow(s,2)*pow(t2,2))/pow(mrho,6) -
        				(1.*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s)*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*pow(t2,2))/
        				 (pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) + (0.375*(eta1 - 1.*eta2)*
        				   (eta1*(-0.6666666666666666*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s) +
        					 eta2*(0.6666666666666666*pow(ma1,2) - 2.*pow(mpion,2) + 0.6666666666666666*pow(mrho,2) + 0.6666666666666666*s))*
        				   (1.*pow(mrho,2) - 0.5*delta*s)*pow(t2,2))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
        				0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
        				   eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(t2,2) +
        				(3.*(1.*pow(mrho,2) - 0.5*delta*s)*(1.*C4*pow(mrho,6) + 0.0833335*pow(mrho,2)*s + pow(mrho,4)*(-0.416667 - 0.333334*C4*s) +
        					 delta*(0.666665*pow(mpion,2)*pow(mrho,2) + 0.333334*pow(mrho,4) - 0.291667*pow(mrho,2)*s - 0.0416667*pow(s,2)))*pow(t2,2))
        				  /(pow(mrho,8) - 1.*pow(mrho,6)*s) + (0.125*(1.*eta1 - 1.*eta2)*
        				   (pow(mrho,2)*(eta1*(1. - 2.*C4*pow(mrho,2)) + eta2*(-1. + 2.*C4*pow(mrho,2))) +
        					 delta*(eta2*(-1.*pow(ma1,2) + 3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
        						eta1*(1.*pow(ma1,2) - 3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s)))*pow(t2,2))/pow(mrho,2) +
        				0.0104167*pow(eta1 - 1.*eta2,4)*pow(t2,3) + (0.166667*pow(delta,2)*pow(t2,3))/pow(mrho,4) +
        				(0.0833333*delta*pow(1.*eta1 - 1.*eta2,2)*pow(t2,3))/pow(mrho,2) +
        				(0.666667*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*pow(t2,3))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) -
        				(0.166667*pow(1.*eta1 - 1.*eta2,2)*(1.*pow(mrho,2) - 0.5*delta*s)*pow(t2,3))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
        				(0.333334*delta*(-2.*pow(mrho,2) + 1.*delta*s)*pow(t2,3))/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
        				(0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
        						pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) + pow(ma1,2)*pow(mpion,2)*
        						 (8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
        						pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
        					 pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
        						pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
        						pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
        						   2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(6.*pow(mrho,2) + 2.*s)))
        					   + pow(eta1,2)*(1.*pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
        						pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
        						   2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
        						   pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
        						pow(mpion,2)*(1.*pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
        						   pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/(1.*pow(ma1,2) - 1.*t1) -
        				(1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*t1) +
        				(0.25*pow(-2. + delta,2)*pow(mpion,2)*t1)/pow(mrho,2) +
        				0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,2) + pow(mrho,2) - 2.*s) + eta1*(2.*pow(mpion,2) + s))*t1 -
        				(0.5*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*(4.*pow(mpion,4)*pow(mrho,2) + 1.*pow(mrho,6) - 3.5*pow(mrho,4)*s + 0.5*pow(s,3) +
        					 pow(mpion,2)*(10.*pow(mrho,4) - 2.*pow(s,2)))*t1)/(pow(mrho,6)*pow(pow(mrho,2) - 1.*s,2)) +
        				(0.25*(eta1 - 1.*eta2)*(1.*pow(mrho,2) - 0.5*delta*s)*
        				   (eta2*(-2.*pow(ma1,4) - 6.*pow(mpion,4) + 1.5*pow(mrho,4) + pow(ma1,2)*(6.*pow(mpion,2) - 2.*pow(mrho,2) - 2.*s) -
        						1.*pow(mrho,2)*s - 0.5*pow(s,2) + pow(mpion,2)*(2.*pow(mrho,2) + 2.*s)) +
        					 eta1*(2.*pow(ma1,4) + 6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
        						1.*pow(s,2) + pow(ma1,2)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s)))*t1)/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
        				0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-6.*pow(ma1,4) - 12.*pow(mpion,4) + 2.*pow(mrho,4) +
        					  pow(ma1,2)*(16.*pow(mpion,2) - 8.*s) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
        				   pow(eta1,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
        					  2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
        				   pow(eta2,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
        					  2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) + 4.*s)))*t1 +
        				(1.*(pow(mpion,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (1.5 - 2.*delta + 0.625*pow(delta,2))*pow(mrho,2)*s +
        						(0.25 - 0.125*delta)*delta*pow(s,2) + pow(mrho,4)*(2.5 - 2.25*delta + 0.5*pow(delta,2) + 2.*C4*s - 1.*C4*delta*s)) +
        					 pow(mrho,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (0.75 - 0.375*delta)*delta*pow(s,2) +
        						pow(mrho,4)*(0.5 - 0.25*delta + 6.*C4*s - 3.*C4*delta*s) +
        						pow(mrho,2)*s*(-0.5 - 0.5*delta + 0.375*pow(delta,2) - 4.*C4*s + 2.*C4*delta*s)))*t1)/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
        				(0.25*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*(eta1*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(ma1,2)*(1. - 2.*C4*pow(mrho,2)) +
        						   pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) - 2.*C4*pow(s,2)) +
        						eta2*(-1.5*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) +
        						   pow(ma1,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s - 4.*C4*pow(mrho,2)*s + 2.*C4*pow(s,2))) +
        					 delta*(eta2*(-1.*pow(ma1,4) - 3.*pow(mpion,4) + 1.*pow(mrho,4) + pow(ma1,2)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
        						   0.25*pow(mrho,2)*s - 0.75*pow(s,2) + pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)) +
        						eta1*(1.*pow(ma1,4) + 3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s +
        						   1.*pow(s,2) + pow(ma1,2)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s))))*t1)/pow(mrho,2) -
        				(0.5*(pow(delta,2)*(1.*pow(mpion,4)*pow(mrho,2) + 0.25*pow(mrho,6) - 0.75*pow(mrho,4)*s + 0.125*pow(mrho,2)*pow(s,2) +
        						0.25*pow(s,3) + pow(mpion,2)*(2.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2))) +
        					 pow(mrho,6)*(1.5 + C4*(-6.*pow(mrho,2) + 6.*s) + pow(C4,2)*(4.*pow(mrho,4) - 8.*pow(mrho,2)*s + 4.*pow(s,2))) +
        					 delta*pow(mrho,2)*(4.*C4*pow(mrho,6) - 0.5*pow(s,2) + pow(mrho,4)*(-1.5 - 3.*C4*s) + pow(mrho,2)*s*(0.5 - 1.*C4*s) +
        						pow(mpion,2)*(6.*C4*pow(mrho,4) + 0.5*s + pow(mrho,2)*(-2.5 - 2.*C4*s))))*t1)/pow(mrho,6) +
        				(3.*(1.*pow(mrho,2) - 0.5*delta*s)*(delta*(0.666667*pow(mpion,4)*pow(mrho,2) + 0.166667*pow(mrho,6) - 0.541667*pow(mrho,4)*s -
        						0.0833333*pow(mrho,2)*pow(s,2) + 0.125*pow(s,3) +
        						pow(mpion,2)*(1.66667*pow(mrho,4) + 0.0833333*pow(mrho,2)*s - 0.416667*pow(s,2))) +
        					 pow(mrho,2)*(1.*C4*pow(mrho,6) - 0.0833333*pow(s,2) + pow(mrho,4)*(-0.416667 - 1.33333*C4*s) +
        						pow(mrho,2)*s*(0.5 + 0.333333*C4*s) + pow(mpion,2)*(2.*C4*pow(mrho,4) + 0.166667*s + pow(mrho,2)*(-0.833333 - 0.666667*C4*s))
        						))*t1)/(pow(mrho,8) - 1.*pow(mrho,6)*s) + 1.*C4*pow(t1,2) + 1.*C4*delta*pow(t1,2) -
        				0.0625*(-2. + delta)*(eta1 - 1.*eta2)*eta2*pow(t1,2) + (0.5*pow(delta,2)*pow(mpion,2)*pow(t1,2))/pow(mrho,4) -
        				(0.25*pow(t1,2))/pow(mrho,2) - (0.5*delta*pow(t1,2))/pow(mrho,2) + (0.25*pow(delta,2)*pow(t1,2))/pow(mrho,2) +
        				(0.25*delta*s*pow(t1,2))/pow(mrho,4) - (0.25*pow(delta,2)*s*pow(t1,2))/pow(mrho,4) -
        				(0.5*C4*delta*s*pow(t1,2))/pow(mrho,2) - (0.0625*pow(delta,2)*pow(s,2)*pow(t1,2))/pow(mrho,6) +
        				(1.*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s)*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*pow(t1,2))/
        				 (pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) - (0.375*(eta1 - 1.*eta2)*
        				   (eta1*(-0.6666666666666666*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s) +
        					 eta2*(0.6666666666666666*pow(ma1,2) - 2.*pow(mpion,2) + 0.6666666666666666*pow(mrho,2) + 0.6666666666666666*s))*
        				   (1.*pow(mrho,2) - 0.5*delta*s)*pow(t1,2))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
        				0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
        				   eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(t1,2) -
        				(3.*(1.*pow(mrho,2) - 0.5*delta*s)*(1.*C4*pow(mrho,6) + 0.0833335*pow(mrho,2)*s + pow(mrho,4)*(-0.416667 - 0.333334*C4*s) +
        					 delta*(0.666665*pow(mpion,2)*pow(mrho,2) + 0.333334*pow(mrho,4) - 0.291667*pow(mrho,2)*s - 0.0416667*pow(s,2)))*pow(t1,2))
        				  /(pow(mrho,8) - 1.*pow(mrho,6)*s) - (0.125*(1.*eta1 - 1.*eta2)*
        				   (pow(mrho,2)*(eta1*(1. - 2.*C4*pow(mrho,2)) + eta2*(-1. + 2.*C4*pow(mrho,2))) +
        					 delta*(eta2*(-1.*pow(ma1,2) + 3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
        						eta1*(1.*pow(ma1,2) - 3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s)))*pow(t1,2))/pow(mrho,2) -
        				0.0104167*pow(eta1 - 1.*eta2,4)*pow(t1,3) - (0.166667*pow(delta,2)*pow(t1,3))/pow(mrho,4) -
        				(0.0833333*delta*pow(1.*eta1 - 1.*eta2,2)*pow(t1,3))/pow(mrho,2) -
        				(0.666667*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*pow(t1,3))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) +
        				(0.166667*pow(1.*eta1 - 1.*eta2,2)*(1.*pow(mrho,2) - 0.5*delta*s)*pow(t1,3))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
        				(0.333334*delta*(-2.*pow(mrho,2) + 1.*delta*s)*pow(t1,3))/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
        				0.0625*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
        					  pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) - 2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
        					  pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
        				   pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) + pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s +
        					  pow(mpion,4)*(-3.*pow(mrho,2) + s) + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
        					  pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2))) +
        				   pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
        					  pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
        					  pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2))))*
        				 log(fabs(-pow(ma1,2) + t2)) - (0.25*(-2. + delta)*(eta1 - 1.*eta2)*
        				   (eta2*(-0.5*pow(ma1,6) - 0.5*pow(mpion,6) + 0.5*pow(mpion,4)*pow(mrho,2) +
        						pow(ma1,4)*(0.5*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) +
        						pow(ma1,2)*pow(mpion,2)*(0.5*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s)) +
        					 eta1*(pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
        						pow(mpion,2)*(1.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
        						pow(ma1,2)*(-2.*pow(mpion,4) - 0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*log(fabs(-pow(ma1,2) + t2)))/
        				 (pow(ma1,2) - 1.*pow(mpion,2)) - (0.125*(eta1 - 1.*eta2)*(1.*pow(mrho,2) - 0.5*delta*s)*
        				   (eta1*(4.*pow(ma1,6) - 4.*pow(mpion,6) + pow(mrho,4)*s + 4.*pow(mpion,2)*pow(s,2) - 1.*pow(s,3) +
        						pow(mpion,4)*(-10.*pow(mrho,2) + 2.*s) + pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
        						pow(ma1,2)*(12.*pow(mpion,4) + 2.*pow(mrho,4) + pow(mpion,2)*(16.*pow(mrho,2) - 8.*s) - 8.*pow(mrho,2)*s + 2.*pow(s,2))
        						) + eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 4.*pow(mrho,2) - 4.*s) +
        						pow(mpion,2)*(4.*pow(mpion,4) - 1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
        						pow(ma1,2)*(-12.*pow(mpion,4) + 3.*pow(mrho,4) - 2.*pow(mrho,2)*s - 1.*pow(s,2) + pow(mpion,2)*(4.*pow(mrho,2) + 4.*s))
        						))*log(fabs(-pow(ma1,2) + t2)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
        				(0.25*(1.*eta1 - 1.*eta2)*(delta*(eta1*(1.*pow(ma1,6) - 1.*pow(mpion,6) + pow(mpion,4)*(-2.5*pow(mrho,2) + 0.5*s) +
        						   pow(mpion,2)*s*(-0.5*pow(mrho,2) + 1.*s) + pow(ma1,4)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s) +
        						   s*(0.5*pow(mrho,4) - 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
        						   pow(ma1,2)*(3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s +
        							  1.*pow(s,2))) + eta2*(-1.*pow(ma1,6) + pow(ma1,4)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
        						   pow(mpion,2)*(1.*pow(mpion,4) - 0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
        						   pow(ma1,2)*(-3.*pow(mpion,4) + 1.*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2) +
        							  pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)))) +
        					 pow(mrho,2)*(eta2*(pow(ma1,4)*(-1. + 2.*C4*pow(mrho,2)) +
        						   pow(mpion,2)*(0.5*pow(mrho,2) + pow(mpion,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s) +
        						   pow(ma1,2)*(2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) + pow(mrho,2)*(-1.5 - 4.*C4*s) +
        							  s*(0.5 + 2.*C4*s))) + eta1*(pow(ma1,4)*(1. - 2.*C4*pow(mrho,2)) + pow(mpion,4)*(1. - 2.*C4*pow(mrho,2)) +
        						   (-0.5*pow(mrho,2) + 0.5*s)*s + pow(ma1,2)*
        							(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) - 2.*C4*pow(s,2)) +
        						   pow(mpion,2)*(-4.*C4*pow(mrho,4) - 1.*s + pow(mrho,2)*(2. + 4.*C4*s)))))*log(fabs(-pow(ma1,2) + t2)))/pow(mrho,2) +
        				0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) + t2)) +
        				(0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*(eta2*pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s) +
        					 eta1*(-0.5*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s))*log(fabs(-pow(mpion,2) + t2)))/
        				 (-1.*pow(ma1,2) + 1.*pow(mpion,2)) - (2.*(-0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,4)*s +
        					 pow(mpion,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (0.5 - 0.25*delta)*delta*pow(s,2) +
        						pow(mrho,4)*(1. - 0.5*delta + 4.*C4*s - 2.*C4*delta*s) +
        						pow(mrho,2)*s*(1. - 2.*delta + 0.75*pow(delta,2) - 2.*C4*s + 1.*C4*delta*s)))*log(fabs(-pow(mpion,2) + t2)))/
        				 (pow(mrho,4) - 1.*pow(mrho,2)*s) - 0.0625*pow(eta1 - 1.*eta2,2)*
        				 (eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
        					  pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) - 2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
        					  pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
        				   pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) + pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s +
        					  pow(mpion,4)*(-3.*pow(mrho,2) + s) + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
        					  pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2))) +
        				   pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
        					  pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
        					  pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2))))*
        				 log(fabs(-pow(ma1,2) + t1)) + (0.25*(-2. + delta)*(eta1 - 1.*eta2)*
        				   (eta2*(-0.5*pow(ma1,6) - 0.5*pow(mpion,6) + 0.5*pow(mpion,4)*pow(mrho,2) +
        						pow(ma1,4)*(0.5*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) +
        						pow(ma1,2)*pow(mpion,2)*(0.5*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s)) +
        					 eta1*(pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
        						pow(mpion,2)*(1.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
        						pow(ma1,2)*(-2.*pow(mpion,4) - 0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*log(fabs(-pow(ma1,2) + t1)))/
        				 (pow(ma1,2) - 1.*pow(mpion,2)) + (0.125*(eta1 - 1.*eta2)*(1.*pow(mrho,2) - 0.5*delta*s)*
        				   (eta1*(4.*pow(ma1,6) - 4.*pow(mpion,6) + pow(mrho,4)*s + 4.*pow(mpion,2)*pow(s,2) - 1.*pow(s,3) +
        						pow(mpion,4)*(-10.*pow(mrho,2) + 2.*s) + pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
        						pow(ma1,2)*(12.*pow(mpion,4) + 2.*pow(mrho,4) + pow(mpion,2)*(16.*pow(mrho,2) - 8.*s) - 8.*pow(mrho,2)*s + 2.*pow(s,2))
        						) + eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 4.*pow(mrho,2) - 4.*s) +
        						pow(mpion,2)*(4.*pow(mpion,4) - 1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
        						pow(ma1,2)*(-12.*pow(mpion,4) + 3.*pow(mrho,4) - 2.*pow(mrho,2)*s - 1.*pow(s,2) + pow(mpion,2)*(4.*pow(mrho,2) + 4.*s))
        						))*log(fabs(-pow(ma1,2) + t1)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
        				(0.25*(1.*eta1 - 1.*eta2)*(delta*(eta1*(1.*pow(ma1,6) - 1.*pow(mpion,6) + pow(mpion,4)*(-2.5*pow(mrho,2) + 0.5*s) +
        						   pow(mpion,2)*s*(-0.5*pow(mrho,2) + 1.*s) + pow(ma1,4)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s) +
        						   s*(0.5*pow(mrho,4) - 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
        						   pow(ma1,2)*(3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s +
        							  1.*pow(s,2))) + eta2*(-1.*pow(ma1,6) + pow(ma1,4)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
        						   pow(mpion,2)*(1.*pow(mpion,4) - 0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
        						   pow(ma1,2)*(-3.*pow(mpion,4) + 1.*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2) +
        							  pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)))) +
        					 pow(mrho,2)*(eta2*(pow(ma1,4)*(-1. + 2.*C4*pow(mrho,2)) +
        						   pow(mpion,2)*(0.5*pow(mrho,2) + pow(mpion,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s) +
        						   pow(ma1,2)*(2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) + pow(mrho,2)*(-1.5 - 4.*C4*s) +
        							  s*(0.5 + 2.*C4*s))) + eta1*(pow(ma1,4)*(1. - 2.*C4*pow(mrho,2)) + pow(mpion,4)*(1. - 2.*C4*pow(mrho,2)) +
        						   (-0.5*pow(mrho,2) + 0.5*s)*s + pow(ma1,2)*
        							(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) - 2.*C4*pow(s,2)) +
        						   pow(mpion,2)*(-4.*C4*pow(mrho,4) - 1.*s + pow(mrho,2)*(2. + 4.*C4*s)))))*log(fabs(-pow(ma1,2) + t1)))/pow(mrho,2) -
        				0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) + t1)) -
        				(0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*(eta2*pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s) +
        					 eta1*(-0.5*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s))*log(fabs(-pow(mpion,2) + t1)))/
        				 (-1.*pow(ma1,2) + 1.*pow(mpion,2)) + (2.*(-0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,4)*s +
        					 pow(mpion,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (0.5 - 0.25*delta)*delta*pow(s,2) +
        						pow(mrho,4)*(1. - 0.5*delta + 4.*C4*s - 2.*C4*delta*s) +
        						pow(mrho,2)*s*(1. - 2.*delta + 0.75*pow(delta,2) - 2.*C4*s + 1.*C4*delta*s)))*log(fabs(-pow(mpion,2) + t1)))/
        				 (pow(mrho,4) - 1.*pow(mrho,2)*s)))/(16.*Pi*(4*pow(mpion,2) - s)*s));

             process_list.push_back(make_unique<CollisionBranch>(
                *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          }

          //dummy: just for stable rho
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

          xsection = 0.0000000000000001 * to_mb;
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

          mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
          t1 = mandelstam_t[1];
          t2 = mandelstam_t[0];

           xsection = to_mb*1/3.0*(pow(Const,2)*pow(ghat,4)*((pow(eta1 - eta2,2)*(-2*eta1*eta2*
                      (pow(ma1,8) + pow(m_pi,8) - pow(m_pi,4)*pow(mrho,4) - 2*pow(ma1,2)*pow(m_pi,2)*(pow(m_pi,2) - pow(mrho,2))*(pow(mrho,2) + s) +
                        pow(ma1,6)*(-4*pow(m_pi,2) + 2*s) + pow(ma1,4)*
                         (4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2))) +
                     pow(eta2,2)*(pow(ma1,8) + pow(m_pi,4)*pow(pow(m_pi,2) - pow(mrho,2),2) + 2*pow(ma1,6)*(-2*pow(m_pi,2) + pow(mrho,2) + s) +
                        2*pow(ma1,2)*pow(m_pi,2)*(-pow(mrho,4) + pow(m_pi,2)*(2*pow(mrho,2) - s) + pow(mrho,2)*s) +
                        pow(ma1,4)*(4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s))) +
                     pow(eta1,2)*(pow(ma1,8) + pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) - 2*pow(ma1,6)*(2*pow(m_pi,2) + pow(mrho,2) - s) -
                        2*pow(m_pi,2)*pow(mrho,4)*s + pow(m_pi,4)*(3*pow(mrho,4) + 2*pow(mrho,2)*s) +
                        pow(ma1,4)*(4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) - 4*pow(mrho,2)*s + 2*pow(s,2)) -
                        2*pow(ma1,2)*(pow(mrho,2)*s*(-pow(mrho,2) + s) + pow(m_pi,4)*(3*pow(mrho,2) + s) + pow(m_pi,2)*(2*pow(mrho,4) - 3*pow(mrho,2)*s)))))
                  /((pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*(pow(ma1,2) - t2)) +
                (8*pow(-2 + delta,2)*pow(m_pi,2)*(4*pow(m_pi,2) - pow(mrho,2)))/
                 ((pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*(pow(m_pi,2) - t2)) -
                (8*pow(-2 + delta,2)*pow(m_pi,2)*t2)/(pow(mrho,2)*pow(pow(m_pi,2) - s,2)) -
                (8*pow(-2 + delta,2)*pow(m_pi,2)*t2)/(pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (8*(-2 + delta)*(-8*C4*pow(mrho,4) + pow(m_pi,2)*(2 + delta - 8*C4*pow(mrho,2)) - (2 + 3*delta)*s + pow(mrho,2)*(-2 + 3*delta + 16*C4*s))*t2)/
                 (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (4*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(eta2*(pow(m_pi,2) + pow(mrho,2) - 2*s)*(pow(m_pi,2) + s) +
                     eta1*(-2*pow(m_pi,4) + pow(mrho,4) - 3*pow(mrho,2)*s + 2*pow(s,2) + pow(m_pi,2)*(pow(mrho,2) + s)))*t2)/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (pow(eta1 - eta2,2)*(pow(eta1,2)*(3*pow(ma1,4) + 4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) -
                        4*pow(ma1,2)*(2*pow(m_pi,2) + pow(mrho,2) - s) - 4*pow(mrho,2)*s + 2*pow(s,2)) +
                     pow(eta2,2)*(3*pow(ma1,4) + 4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s) +
                        4*pow(ma1,2)*(-2*pow(m_pi,2) + pow(mrho,2) + s)) -
                     2*eta1*eta2*(3*pow(ma1,4) + 4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2) +
                        pow(ma1,2)*(-8*pow(m_pi,2) + 4*s)))*t2)/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) +
                (8*(pow(delta,2)*(8*pow(m_pi,4) + 3*pow(mrho,4) + 4*pow(m_pi,2)*(3*pow(mrho,2) - 2*s) - 6*pow(mrho,2)*s + 2*pow(s,2)) +
                     4*pow(mrho,4)*(3 + 12*C4*(2*pow(m_pi,2) - s) + 8*pow(C4,2)*pow(-2*pow(m_pi,2) + s,2)) -
                     4*delta*pow(mrho,2)*(16*C4*pow(m_pi,4) + 2*pow(m_pi,2)*(3 + 6*C4*pow(mrho,2) - 8*C4*s) + pow(mrho,2)*(3 - 6*C4*s) + s*(-3 + 4*C4*s)))*t2)/
                 (pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (8*(-2 + delta)*(pow(m_pi,4)*(-2 + 3*delta - 8*C4*pow(mrho,2)) + (pow(mrho,2) - s)*((-2 + 3*delta)*s + pow(mrho,2)*(-2 + delta - 8*C4*s)) +
                     4*pow(m_pi,2)*(2*C4*pow(mrho,4) + delta*s - pow(mrho,2)*(-1 + delta + 4*C4*s)))*t2)/
                 (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (4*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(eta2*(pow(m_pi,2) + s)*(pow(m_pi,4) - pow(m_pi,2)*(pow(mrho,2) - 2*s) + (pow(mrho,2) - s)*s) +
                     eta1*(-4*pow(m_pi,6) + pow(pow(mrho,2) - s,2)*s + pow(m_pi,4)*(3*pow(mrho,2) + s) -
                        pow(m_pi,2)*(pow(mrho,4) - pow(mrho,2)*s + 2*pow(s,2))))*t2)/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,2) - s)*
                   (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (pow(eta1 - eta2,2)*(-2*eta1*eta2*(pow(m_pi,8) - pow(mrho,4)*pow(s,2) + pow(s,4) -
                        pow(m_pi,4)*(pow(mrho,4) + 2*pow(mrho,2)*s - 4*pow(s,2)) + 2*pow(m_pi,2)*s*(pow(mrho,4) + pow(mrho,2)*s - 2*pow(s,2))) +
                     pow(eta2,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + pow(s,2)*pow(pow(mrho,2) + s,2) + pow(m_pi,4)*pow(pow(mrho,2) + 2*s,2) -
                        2*pow(m_pi,2)*s*(pow(mrho,4) + 2*pow(mrho,2)*s + 2*pow(s,2))) +
                     pow(eta1,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) - 4*pow(m_pi,2)*pow(pow(mrho,2) - s,2)*s + pow(pow(mrho,2) - s,2)*pow(s,2) +
                        pow(m_pi,4)*(3*pow(mrho,4) - 6*pow(mrho,2)*s + 4*pow(s,2))))*t2)/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (2*pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(pow(eta1,2)*(pow(ma1,4)*s + pow(m_pi,4)*(-3*pow(mrho,2) + 2*s) +
                        s*(2*pow(mrho,4) - 3*pow(mrho,2)*s + pow(s,2)) - 2*pow(m_pi,2)*(pow(mrho,4) - 4*pow(mrho,2)*s + 2*pow(s,2)) +
                        pow(ma1,2)*(2*pow(m_pi,2)*(pow(mrho,2) - 2*s) + 3*s*(-pow(mrho,2) + s))) -
                     2*eta1*eta2*(pow(ma1,4)*s + s*(2*pow(m_pi,4) + 4*pow(m_pi,2)*(pow(mrho,2) - s) + s*(-2*pow(mrho,2) + s)) +
                        pow(ma1,2)*(pow(m_pi,2)*(pow(mrho,2) - 4*s) + s*(-2*pow(mrho,2) + 3*s))) +
                     pow(eta2,2)*(-4*pow(m_pi,2)*s*(pow(ma1,2) + pow(mrho,2) + s) + pow(m_pi,4)*(pow(mrho,2) + 2*s) +
                        s*(pow(ma1,4) + s*(pow(mrho,2) + s) + pow(ma1,2)*(pow(mrho,2) + 3*s))))*t2)/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (4*(eta1 - eta2)*(pow(ma1,2) - s)*(eta1*(-4*pow(m_pi,4)*(6*C4*pow(mrho,4) + 2*delta*s + pow(mrho,2)*(1 - 2*delta - 8*C4*s)) +
                        2*pow(m_pi,2)*(4*delta*pow(s,2) + pow(mrho,2)*s*(6 - 7*delta - 16*C4*s) + 2*pow(mrho,4)*(-2 + delta + 8*C4*s)) -
                        (pow(mrho,2) - s)*s*(-2*delta*s + pow(mrho,2)*(-6 + 3*delta + 8*C4*s))) +
                     eta2*(delta*(2*pow(m_pi,4)*(pow(mrho,2) + 4*s) + pow(m_pi,2)*(2*pow(mrho,4) + pow(mrho,2)*s - 8*pow(s,2)) +
                           s*(-2*pow(mrho,4) - pow(mrho,2)*s + 2*pow(s,2))) -
                        2*pow(mrho,2)*(4*C4*pow(m_pi,4)*(pow(mrho,2) + 4*s) + pow(m_pi,2)*(s*(5 - 16*C4*s) + pow(mrho,2)*(2 - 8*C4*s)) +
                           s*(s*(-3 + 4*C4*s) + pow(mrho,2)*(-2 + 4*C4*s)))))*t2)/
                 (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                   (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (8*(eta1 - eta2)*(delta*(eta1*(4*pow(m_pi,6) + pow(m_pi,4)*(7*pow(mrho,2) - 8*s) + pow(ma1,4)*(pow(m_pi,2) - s) -
                           pow(ma1,2)*(2*pow(m_pi,2) + pow(mrho,2) - 2*s)*(2*pow(m_pi,2) - s) + pow(m_pi,2)*s*(-8*pow(mrho,2) + 5*s) +
                           s*(pow(mrho,4) + pow(mrho,2)*s - pow(s,2))) +
                        eta2*(-4*pow(m_pi,6) - pow(m_pi,4)*(pow(mrho,2) - 8*s) + pow(ma1,4)*(-pow(m_pi,2) + s) +
                           pow(m_pi,2)*(2*pow(mrho,4) - 5*pow(s,2)) + s*(-2*pow(mrho,4) + pow(mrho,2)*s + pow(s,2)) +
                           pow(ma1,2)*(4*pow(m_pi,4) - 6*pow(m_pi,2)*s + s*(pow(mrho,2) + 2*s)))) -
                     2*pow(mrho,2)*(eta1*(8*C4*pow(m_pi,6) + pow(m_pi,4)*(3 + 8*C4*(pow(mrho,2) - 2*s)) + 2*C4*pow(ma1,4)*(pow(m_pi,2) - s) +
                           2*pow(m_pi,2)*s*(-1 - 6*C4*pow(mrho,2) + 5*C4*s) -
                           pow(ma1,2)*(8*C4*pow(m_pi,4) + pow(m_pi,2)*(1 + 2*C4*(pow(mrho,2) - 6*s)) + 2*C4*s*(-pow(mrho,2) + 2*s)) +
                           s*(-(s*(1 + 2*C4*s)) + pow(mrho,2)*(1 + 4*C4*s))) +
                        eta2*(2*C4*pow(ma1,4)*(-pow(m_pi,2) + s) - (pow(m_pi,2) - s)*
                            (8*C4*pow(m_pi,4) - 2*pow(mrho,2) + s + 2*C4*pow(s,2) + pow(m_pi,2)*(3 - 4*C4*(pow(mrho,2) + 2*s))) +
                           pow(ma1,2)*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(mrho,2) + 2*s) + pow(m_pi,2)*(1 - 2*C4*(pow(mrho,2) + 6*s))))))*t2)/
                 (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (8*(-2 + delta)*(delta - 4*C4*pow(mrho,2))*pow(t2,2))/
                 (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (16*(-2 + delta)*(delta - 4*C4*pow(mrho,2))*s*pow(t2,2))/
                 (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (pow(eta1 - eta2,2)*(pow(eta1,2)*(pow(mrho,2) - s) + 2*eta1*eta2*s - pow(eta2,2)*s)*
                   (pow(m_pi,4) - pow(m_pi,2)*(pow(mrho,2) - 2*s) + (pow(mrho,2) - s)*s)*pow(t2,2))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (2*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(-(eta1*(pow(m_pi,2) + 2*pow(mrho,2) - 3*s)) - eta2*(pow(m_pi,2) + s))*pow(t2,2))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (2*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(pow(m_pi,2) + s)*(-2*eta2*s + eta1*(pow(m_pi,2) - pow(mrho,2) + s))*pow(t2,2))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,2) - s)*
                   (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (pow(eta1 - eta2,3)*(eta1*(pow(ma1,2) - 2*pow(m_pi,2) - pow(mrho,2) + s) - eta2*(pow(ma1,2) - 2*pow(m_pi,2) + pow(mrho,2) + s))*
                   pow(t2,2))/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) -
                (8*(delta - 4*C4*pow(mrho,2))*(delta*(4*pow(m_pi,2) + 3*pow(mrho,2) - 2*s) - 2*pow(mrho,2)*(3 + 8*C4*pow(m_pi,2) - 4*C4*s))*pow(t2,2))/
                 (pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(pow(eta2,2)*s*(pow(ma1,2) - 4*pow(m_pi,2) + pow(mrho,2) + 3*s) +
                     pow(eta1,2)*(2*pow(m_pi,2)*(pow(mrho,2) - 2*s) + s*(pow(ma1,2) - 3*pow(mrho,2) + 3*s)) -
                     2*eta1*eta2*(pow(m_pi,2)*(pow(mrho,2) - 4*s) + s*(pow(ma1,2) - 2*pow(mrho,2) + 3*s)))*pow(t2,2))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (2*(eta1 - eta2)*(pow(ma1,2) - s)*(eta1*(4*delta*pow(s,2) - 2*pow(mrho,2)*s*(-2 + 3*delta + 8*C4*s) + pow(mrho,4)*(-2 + delta + 16*C4*s) -
                        2*pow(m_pi,2)*(8*C4*pow(mrho,4) + 4*delta*s + pow(mrho,2)*(2 - 3*delta - 16*C4*s))) +
                     eta2*(pow(m_pi,2)*(8*delta*s + pow(mrho,2)*(-2 + delta - 32*C4*s)) + s*(-4*delta*s + pow(mrho,2)*(-2 + delta + 16*C4*s))))*pow(t2,2))/
                 (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                   (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (4*(eta1 - eta2)*(delta*(eta1*(pow(ma1,2)*(pow(m_pi,2) - s) - (2*pow(m_pi,2) + pow(mrho,2) - 2*s)*(2*pow(m_pi,2) - s)) +
                        eta2*(4*pow(m_pi,4) - 6*pow(m_pi,2)*s + pow(ma1,2)*(-pow(m_pi,2) + s) + s*(pow(mrho,2) + 2*s))) +
                     2*pow(mrho,2)*(eta1*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(ma1,2) - pow(mrho,2) + 2*s) +
                           pow(m_pi,2)*(1 - 2*C4*(pow(ma1,2) - pow(mrho,2) + 6*s))) -
                        eta2*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(ma1,2) + pow(mrho,2) + 2*s) - pow(m_pi,2)*(-1 + 2*C4*(pow(ma1,2) + pow(mrho,2) + 6*s)))))*
                   pow(t2,2))/(pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (pow(eta1 - eta2,4)*pow(t2,3))/(3.*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (8*pow(eta1 - eta2,2)*(delta - 4*C4*pow(mrho,2))*pow(t2,3))/
                 (3.*pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (16*pow(delta - 4*C4*pow(mrho,2),2)*pow(t2,3))/
                 (3.*pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (4*(-2 + delta)*eta1*(eta1 - eta2)*(pow(ma1,2) - s)*pow(t2,3))/
                 (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (2*pow(eta1 - eta2,4)*(pow(ma1,2) - s)*s*pow(t2,3))/
                 (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (2*pow(eta1 - eta2,2)*s*(-2*eta1*eta2*s + pow(eta2,2)*s + pow(eta1,2)*(-pow(mrho,2) + s))*pow(t2,3))/
                 (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (4*(eta1 - eta2)*(pow(ma1,2) - s)*(2*eta2*(delta - 4*C4*pow(mrho,2))*s + eta1*(-2*delta*s + pow(mrho,2)*(-2 + delta + 8*C4*s)))*pow(t2,3))/
                 (3.*pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                   (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (pow(eta1 - eta2,2)*(-2*eta1*eta2*(pow(ma1,8) + pow(m_pi,8) - pow(m_pi,4)*pow(mrho,4) -
                        2*pow(ma1,2)*pow(m_pi,2)*(pow(m_pi,2) - pow(mrho,2))*(pow(mrho,2) + s) + pow(ma1,6)*(-4*pow(m_pi,2) + 2*s) +
                        pow(ma1,4)*(4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2))) +
                     pow(eta2,2)*(pow(ma1,8) + pow(m_pi,4)*pow(pow(m_pi,2) - pow(mrho,2),2) + 2*pow(ma1,6)*(-2*pow(m_pi,2) + pow(mrho,2) + s) +
                        2*pow(ma1,2)*pow(m_pi,2)*(-pow(mrho,4) + pow(m_pi,2)*(2*pow(mrho,2) - s) + pow(mrho,2)*s) +
                        pow(ma1,4)*(4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s))) +
                     pow(eta1,2)*(pow(ma1,8) + pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) - 2*pow(ma1,6)*(2*pow(m_pi,2) + pow(mrho,2) - s) -
                        2*pow(m_pi,2)*pow(mrho,4)*s + pow(m_pi,4)*(3*pow(mrho,4) + 2*pow(mrho,2)*s) +
                        pow(ma1,4)*(4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) - 4*pow(mrho,2)*s + 2*pow(s,2)) -
                        2*pow(ma1,2)*(pow(mrho,2)*s*(-pow(mrho,2) + s) + pow(m_pi,4)*(3*pow(mrho,2) + s) + pow(m_pi,2)*(2*pow(mrho,4) - 3*pow(mrho,2)*s)))))
                  /((pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*(pow(ma1,2) - t1)) -
                (8*pow(-2 + delta,2)*pow(m_pi,2)*(4*pow(m_pi,2) - pow(mrho,2)))/
                 ((pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*(pow(m_pi,2) - t1)) +
                (8*pow(-2 + delta,2)*pow(m_pi,2)*t1)/(pow(mrho,2)*pow(pow(m_pi,2) - s,2)) +
                (8*pow(-2 + delta,2)*pow(m_pi,2)*t1)/(pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (8*(-2 + delta)*(-8*C4*pow(mrho,4) + pow(m_pi,2)*(2 + delta - 8*C4*pow(mrho,2)) - (2 + 3*delta)*s + pow(mrho,2)*(-2 + 3*delta + 16*C4*s))*t1)/
                 (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (4*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(eta2*(pow(m_pi,2) + pow(mrho,2) - 2*s)*(pow(m_pi,2) + s) +
                     eta1*(-2*pow(m_pi,4) + pow(mrho,4) - 3*pow(mrho,2)*s + 2*pow(s,2) + pow(m_pi,2)*(pow(mrho,2) + s)))*t1)/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (pow(eta1 - eta2,2)*(pow(eta1,2)*(3*pow(ma1,4) + 4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) -
                        4*pow(ma1,2)*(2*pow(m_pi,2) + pow(mrho,2) - s) - 4*pow(mrho,2)*s + 2*pow(s,2)) +
                     pow(eta2,2)*(3*pow(ma1,4) + 4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s) +
                        4*pow(ma1,2)*(-2*pow(m_pi,2) + pow(mrho,2) + s)) -
                     2*eta1*eta2*(3*pow(ma1,4) + 4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2) +
                        pow(ma1,2)*(-8*pow(m_pi,2) + 4*s)))*t1)/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) -
                (8*(pow(delta,2)*(8*pow(m_pi,4) + 3*pow(mrho,4) + 4*pow(m_pi,2)*(3*pow(mrho,2) - 2*s) - 6*pow(mrho,2)*s + 2*pow(s,2)) +
                     4*pow(mrho,4)*(3 + 12*C4*(2*pow(m_pi,2) - s) + 8*pow(C4,2)*pow(-2*pow(m_pi,2) + s,2)) -
                     4*delta*pow(mrho,2)*(16*C4*pow(m_pi,4) + 2*pow(m_pi,2)*(3 + 6*C4*pow(mrho,2) - 8*C4*s) + pow(mrho,2)*(3 - 6*C4*s) + s*(-3 + 4*C4*s)))*t1)/
                 (pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (8*(-2 + delta)*(pow(m_pi,4)*(-2 + 3*delta - 8*C4*pow(mrho,2)) + (pow(mrho,2) - s)*((-2 + 3*delta)*s + pow(mrho,2)*(-2 + delta - 8*C4*s)) +
                     4*pow(m_pi,2)*(2*C4*pow(mrho,4) + delta*s - pow(mrho,2)*(-1 + delta + 4*C4*s)))*t1)/
                 (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (4*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(eta2*(pow(m_pi,2) + s)*(pow(m_pi,4) - pow(m_pi,2)*(pow(mrho,2) - 2*s) + (pow(mrho,2) - s)*s) +
                     eta1*(-4*pow(m_pi,6) + pow(pow(mrho,2) - s,2)*s + pow(m_pi,4)*(3*pow(mrho,2) + s) -
                        pow(m_pi,2)*(pow(mrho,4) - pow(mrho,2)*s + 2*pow(s,2))))*t1)/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,2) - s)*
                   (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (pow(eta1 - eta2,2)*(-2*eta1*eta2*(pow(m_pi,8) - pow(mrho,4)*pow(s,2) + pow(s,4) -
                        pow(m_pi,4)*(pow(mrho,4) + 2*pow(mrho,2)*s - 4*pow(s,2)) + 2*pow(m_pi,2)*s*(pow(mrho,4) + pow(mrho,2)*s - 2*pow(s,2))) +
                     pow(eta2,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + pow(s,2)*pow(pow(mrho,2) + s,2) + pow(m_pi,4)*pow(pow(mrho,2) + 2*s,2) -
                        2*pow(m_pi,2)*s*(pow(mrho,4) + 2*pow(mrho,2)*s + 2*pow(s,2))) +
                     pow(eta1,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) - 4*pow(m_pi,2)*pow(pow(mrho,2) - s,2)*s + pow(pow(mrho,2) - s,2)*pow(s,2) +
                        pow(m_pi,4)*(3*pow(mrho,4) - 6*pow(mrho,2)*s + 4*pow(s,2))))*t1)/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (2*pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(pow(eta1,2)*(pow(ma1,4)*s + pow(m_pi,4)*(-3*pow(mrho,2) + 2*s) +
                        s*(2*pow(mrho,4) - 3*pow(mrho,2)*s + pow(s,2)) - 2*pow(m_pi,2)*(pow(mrho,4) - 4*pow(mrho,2)*s + 2*pow(s,2)) +
                        pow(ma1,2)*(2*pow(m_pi,2)*(pow(mrho,2) - 2*s) + 3*s*(-pow(mrho,2) + s))) -
                     2*eta1*eta2*(pow(ma1,4)*s + s*(2*pow(m_pi,4) + 4*pow(m_pi,2)*(pow(mrho,2) - s) + s*(-2*pow(mrho,2) + s)) +
                        pow(ma1,2)*(pow(m_pi,2)*(pow(mrho,2) - 4*s) + s*(-2*pow(mrho,2) + 3*s))) +
                     pow(eta2,2)*(-4*pow(m_pi,2)*s*(pow(ma1,2) + pow(mrho,2) + s) + pow(m_pi,4)*(pow(mrho,2) + 2*s) +
                        s*(pow(ma1,4) + s*(pow(mrho,2) + s) + pow(ma1,2)*(pow(mrho,2) + 3*s))))*t1)/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (4*(eta1 - eta2)*(pow(ma1,2) - s)*(eta1*(4*pow(m_pi,4)*(6*C4*pow(mrho,4) + 2*delta*s + pow(mrho,2)*(1 - 2*delta - 8*C4*s)) -
                        2*pow(m_pi,2)*(4*delta*pow(s,2) + pow(mrho,2)*s*(6 - 7*delta - 16*C4*s) + 2*pow(mrho,4)*(-2 + delta + 8*C4*s)) +
                        (pow(mrho,2) - s)*s*(-2*delta*s + pow(mrho,2)*(-6 + 3*delta + 8*C4*s))) +
                     eta2*(-(delta*(2*pow(m_pi,4)*(pow(mrho,2) + 4*s) + pow(m_pi,2)*(2*pow(mrho,4) + pow(mrho,2)*s - 8*pow(s,2)) +
                             s*(-2*pow(mrho,4) - pow(mrho,2)*s + 2*pow(s,2)))) +
                        2*pow(mrho,2)*(4*C4*pow(m_pi,4)*(pow(mrho,2) + 4*s) + pow(m_pi,2)*(s*(5 - 16*C4*s) + pow(mrho,2)*(2 - 8*C4*s)) +
                           s*(s*(-3 + 4*C4*s) + pow(mrho,2)*(-2 + 4*C4*s)))))*t1)/
                 (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                   (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (8*(eta1 - eta2)*(delta*(eta1*(4*pow(m_pi,6) + pow(m_pi,4)*(7*pow(mrho,2) - 8*s) + pow(ma1,4)*(pow(m_pi,2) - s) -
                           pow(ma1,2)*(2*pow(m_pi,2) + pow(mrho,2) - 2*s)*(2*pow(m_pi,2) - s) + pow(m_pi,2)*s*(-8*pow(mrho,2) + 5*s) +
                           s*(pow(mrho,4) + pow(mrho,2)*s - pow(s,2))) +
                        eta2*(-4*pow(m_pi,6) - pow(m_pi,4)*(pow(mrho,2) - 8*s) + pow(ma1,4)*(-pow(m_pi,2) + s) +
                           pow(m_pi,2)*(2*pow(mrho,4) - 5*pow(s,2)) + s*(-2*pow(mrho,4) + pow(mrho,2)*s + pow(s,2)) +
                           pow(ma1,2)*(4*pow(m_pi,4) - 6*pow(m_pi,2)*s + s*(pow(mrho,2) + 2*s)))) -
                     2*pow(mrho,2)*(eta1*(8*C4*pow(m_pi,6) + pow(m_pi,4)*(3 + 8*C4*(pow(mrho,2) - 2*s)) + 2*C4*pow(ma1,4)*(pow(m_pi,2) - s) +
                           2*pow(m_pi,2)*s*(-1 - 6*C4*pow(mrho,2) + 5*C4*s) -
                           pow(ma1,2)*(8*C4*pow(m_pi,4) + pow(m_pi,2)*(1 + 2*C4*(pow(mrho,2) - 6*s)) + 2*C4*s*(-pow(mrho,2) + 2*s)) +
                           s*(-(s*(1 + 2*C4*s)) + pow(mrho,2)*(1 + 4*C4*s))) +
                        eta2*(2*C4*pow(ma1,4)*(-pow(m_pi,2) + s) - (pow(m_pi,2) - s)*
                            (8*C4*pow(m_pi,4) - 2*pow(mrho,2) + s + 2*C4*pow(s,2) + pow(m_pi,2)*(3 - 4*C4*(pow(mrho,2) + 2*s))) +
                           pow(ma1,2)*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(mrho,2) + 2*s) + pow(m_pi,2)*(1 - 2*C4*(pow(mrho,2) + 6*s))))))*t1)/
                 (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (8*(-2 + delta)*(delta - 4*C4*pow(mrho,2))*pow(t1,2))/
                 (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (16*(-2 + delta)*(delta - 4*C4*pow(mrho,2))*s*pow(t1,2))/
                 (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (pow(eta1 - eta2,2)*(pow(eta1,2)*(pow(mrho,2) - s) + 2*eta1*eta2*s - pow(eta2,2)*s)*
                   (pow(m_pi,4) - pow(m_pi,2)*(pow(mrho,2) - 2*s) + (pow(mrho,2) - s)*s)*pow(t1,2))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (2*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(-(eta1*(pow(m_pi,2) + 2*pow(mrho,2) - 3*s)) - eta2*(pow(m_pi,2) + s))*pow(t1,2))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (2*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(pow(m_pi,2) + s)*(-2*eta2*s + eta1*(pow(m_pi,2) - pow(mrho,2) + s))*pow(t1,2))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,2) - s)*
                   (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (pow(eta1 - eta2,3)*(-(eta1*(pow(ma1,2) - 2*pow(m_pi,2) - pow(mrho,2) + s)) + eta2*(pow(ma1,2) - 2*pow(m_pi,2) + pow(mrho,2) + s))*
                   pow(t1,2))/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) +
                (8*(delta - 4*C4*pow(mrho,2))*(delta*(4*pow(m_pi,2) + 3*pow(mrho,2) - 2*s) - 2*pow(mrho,2)*(3 + 8*C4*pow(m_pi,2) - 4*C4*s))*pow(t1,2))/
                 (pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(pow(eta2,2)*s*(pow(ma1,2) - 4*pow(m_pi,2) + pow(mrho,2) + 3*s) +
                     pow(eta1,2)*(2*pow(m_pi,2)*(pow(mrho,2) - 2*s) + s*(pow(ma1,2) - 3*pow(mrho,2) + 3*s)) -
                     2*eta1*eta2*(pow(m_pi,2)*(pow(mrho,2) - 4*s) + s*(pow(ma1,2) - 2*pow(mrho,2) + 3*s)))*pow(t1,2))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (2*(eta1 - eta2)*(pow(ma1,2) - s)*(eta1*(4*delta*pow(s,2) - 2*pow(mrho,2)*s*(-2 + 3*delta + 8*C4*s) + pow(mrho,4)*(-2 + delta + 16*C4*s) -
                        2*pow(m_pi,2)*(8*C4*pow(mrho,4) + 4*delta*s + pow(mrho,2)*(2 - 3*delta - 16*C4*s))) +
                     eta2*(pow(m_pi,2)*(8*delta*s + pow(mrho,2)*(-2 + delta - 32*C4*s)) + s*(-4*delta*s + pow(mrho,2)*(-2 + delta + 16*C4*s))))*pow(t1,2))/
                 (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                   (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (4*(eta1 - eta2)*(delta*(eta1*(pow(ma1,2)*(pow(m_pi,2) - s) - (2*pow(m_pi,2) + pow(mrho,2) - 2*s)*(2*pow(m_pi,2) - s)) +
                        eta2*(4*pow(m_pi,4) - 6*pow(m_pi,2)*s + pow(ma1,2)*(-pow(m_pi,2) + s) + s*(pow(mrho,2) + 2*s))) +
                     2*pow(mrho,2)*(eta1*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(ma1,2) - pow(mrho,2) + 2*s) +
                           pow(m_pi,2)*(1 - 2*C4*(pow(ma1,2) - pow(mrho,2) + 6*s))) -
                        eta2*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(ma1,2) + pow(mrho,2) + 2*s) - pow(m_pi,2)*(-1 + 2*C4*(pow(ma1,2) + pow(mrho,2) + 6*s)))))*
                   pow(t1,2))/(pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (pow(eta1 - eta2,4)*pow(t1,3))/(3.*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (8*pow(eta1 - eta2,2)*(delta - 4*C4*pow(mrho,2))*pow(t1,3))/
                 (3.*pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (16*pow(delta - 4*C4*pow(mrho,2),2)*pow(t1,3))/
                 (3.*pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (4*(-2 + delta)*eta1*(eta1 - eta2)*(pow(ma1,2) - s)*pow(t1,3))/
                 (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (2*pow(eta1 - eta2,4)*(pow(ma1,2) - s)*s*pow(t1,3))/
                 (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (2*pow(eta1 - eta2,2)*s*(-2*eta1*eta2*s + pow(eta2,2)*s + pow(eta1,2)*(-pow(mrho,2) + s))*pow(t1,3))/
                 (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (4*(eta1 - eta2)*(pow(ma1,2) - s)*(2*eta2*(delta - 4*C4*pow(mrho,2))*s + eta1*(-2*delta*s + pow(mrho,2)*(-2 + delta + 8*C4*s)))*pow(t1,3))/
                 (3.*pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                   (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (2*pow(eta1 - eta2,2)*(pow(eta1,2)*(2*pow(ma1,6) - 3*pow(ma1,4)*(2*pow(m_pi,2) + pow(mrho,2) - s) + pow(mrho,2)*(pow(mrho,2) - s)*s -
                        pow(m_pi,4)*(3*pow(mrho,2) + s) + pow(m_pi,2)*(-2*pow(mrho,4) + 3*pow(mrho,2)*s) +
                        pow(ma1,2)*(4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) - 4*pow(mrho,2)*s + 2*pow(s,2))) -
                     2*eta1*eta2*(2*pow(ma1,6) - pow(m_pi,2)*(pow(m_pi,2) - pow(mrho,2))*(pow(mrho,2) + s) + pow(ma1,4)*(-6*pow(m_pi,2) + 3*s) +
                        pow(ma1,2)*(4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2))) +
                     pow(eta2,2)*(2*pow(ma1,6) + 3*pow(ma1,4)*(-2*pow(m_pi,2) + pow(mrho,2) + s) +
                        pow(m_pi,2)*(-pow(mrho,4) + pow(m_pi,2)*(2*pow(mrho,2) - s) + pow(mrho,2)*s) +
                        pow(ma1,2)*(4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s))))*log(fabs(-pow(ma1,2) + t2)))
                  /(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) -
                (2*pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(-2*eta1*eta2*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + 2*pow(ma1,2)*pow(m_pi,4)*s +
                        pow(m_pi,2)*(pow(ma1,4)*(pow(mrho,2) - 4*s) + 4*pow(ma1,2)*(pow(mrho,2) - s)*s + pow(mrho,2)*pow(s,2)) +
                        pow(ma1,2)*s*(pow(ma1,4) + s*(-2*pow(mrho,2) + s) + pow(ma1,2)*(-2*pow(mrho,2) + 3*s))) +
                     pow(eta2,2)*(pow(m_pi,8) - 4*pow(ma1,2)*pow(m_pi,2)*s*(pow(ma1,2) + pow(mrho,2) + s) +
                        pow(m_pi,4)*(pow(mrho,2)*s + pow(ma1,2)*(pow(mrho,2) + 2*s)) +
                        pow(ma1,2)*s*(pow(ma1,4) + s*(pow(mrho,2) + s) + pow(ma1,2)*(pow(mrho,2) + 3*s))) +
                     pow(eta1,2)*(pow(m_pi,8) + pow(ma1,2)*s*(pow(ma1,4) + 2*pow(mrho,4) - 3*pow(ma1,2)*(pow(mrho,2) - s) - 3*pow(mrho,2)*s + pow(s,2)) +
                        pow(m_pi,4)*(2*pow(mrho,4) - 3*pow(mrho,2)*s + pow(ma1,2)*(-3*pow(mrho,2) + 2*s)) +
                        2*pow(m_pi,2)*(pow(ma1,4)*(pow(mrho,2) - 2*s) + pow(mrho,2)*s*(-pow(mrho,2) + s) -
                           pow(ma1,2)*(pow(mrho,4) - 4*pow(mrho,2)*s + 2*pow(s,2)))))*log(fabs(-pow(ma1,2) + t2)))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (8*(eta1 - eta2)*(delta*(eta2*(pow(m_pi,6)*pow(mrho,2)*(2*pow(m_pi,2) - s) + pow(ma1,8)*(-pow(m_pi,2) + s) +
                           pow(ma1,6)*(5*pow(m_pi,4) - 7*pow(m_pi,2)*s + s*(pow(mrho,2) + 2*s)) +
                           pow(ma1,4)*(-8*pow(m_pi,6) - pow(m_pi,4)*(pow(mrho,2) - 14*s) + pow(m_pi,2)*(2*pow(mrho,4) - pow(mrho,2)*s - 7*pow(s,2)) +
                              s*(-2*pow(mrho,4) + pow(mrho,2)*s + pow(s,2))) +
                           pow(ma1,2)*pow(m_pi,2)*(4*pow(m_pi,6) + pow(m_pi,4)*(pow(mrho,2) - 8*s) + s*(2*pow(mrho,4) + pow(mrho,2)*s - pow(s,2)) +
                              pow(m_pi,2)*(-2*pow(mrho,4) - 3*pow(mrho,2)*s + 5*pow(s,2)))) +
                        eta1*(pow(ma1,8)*(pow(m_pi,2) - s) + pow(ma1,6)*(-5*pow(m_pi,4) + (pow(mrho,2) - 2*s)*s + pow(m_pi,2)*(-2*pow(mrho,2) + 7*s)) +
                           pow(m_pi,2)*pow(mrho,2)*(2*pow(m_pi,6) + pow(m_pi,4)*(4*pow(mrho,2) - 5*s) + pow(mrho,4)*s -
                              pow(m_pi,2)*(pow(mrho,4) + 3*pow(mrho,2)*s - 2*pow(s,2))) +
                           pow(ma1,4)*(8*pow(m_pi,6) + pow(m_pi,4)*(9*pow(mrho,2) - 14*s) + pow(m_pi,2)*s*(-9*pow(mrho,2) + 7*s) +
                              s*(pow(mrho,4) + pow(mrho,2)*s - pow(s,2))) +
                           pow(ma1,2)*(-4*pow(m_pi,8) + pow(mrho,4)*s*(-pow(mrho,2) + s) + pow(m_pi,6)*(-11*pow(mrho,2) + 8*s) +
                              pow(m_pi,4)*(-3*pow(mrho,4) + 17*pow(mrho,2)*s - 5*pow(s,2)) + pow(m_pi,2)*(pow(mrho,6) - 5*pow(mrho,2)*pow(s,2) + pow(s,3))
                              ))) - 2*pow(mrho,2)*(eta2*(pow(m_pi,8)*(1 + 2*C4*pow(mrho,2)) - 2*C4*pow(m_pi,6)*pow(mrho,2)*s +
                           2*C4*pow(ma1,8)*(-pow(m_pi,2) + s) + pow(ma1,4)*
                            (-16*C4*pow(m_pi,6) + pow(m_pi,4)*(-4 + 6*C4*pow(mrho,2) + 28*C4*s) +
                              2*pow(m_pi,2)*(pow(mrho,2) + s - 3*C4*pow(mrho,2)*s - 7*C4*pow(s,2)) + s*(-2*pow(mrho,2) + s + 2*C4*pow(s,2))) +
                           pow(ma1,6)*(10*C4*pow(m_pi,4) + 2*C4*s*(pow(mrho,2) + 2*s) + pow(m_pi,2)*(1 - 2*C4*(pow(mrho,2) + 7*s))) +
                           pow(ma1,2)*pow(m_pi,2)*(8*C4*pow(m_pi,6) - 2*pow(m_pi,4)*(-2 + 3*C4*pow(mrho,2) + 8*C4*s) +
                              s*(2*pow(mrho,2) + s - 2*C4*pow(s,2)) + 2*pow(m_pi,2)*(pow(mrho,2)*(-1 + 3*C4*s) + s*(-3 + 5*C4*s)))) +
                        eta1*(pow(m_pi,8)*(-1 + 6*C4*pow(mrho,2)) + 2*C4*pow(ma1,8)*(pow(m_pi,2) - s) + pow(m_pi,2)*pow(mrho,4)*s +
                           2*pow(m_pi,6)*pow(mrho,2)*(2 - 5*C4*s) - pow(ma1,6)*
                            (10*C4*pow(m_pi,4) + pow(m_pi,2)*(1 + 2*C4*(pow(mrho,2) - 7*s)) + 2*C4*s*(-pow(mrho,2) + 2*s)) -
                           pow(m_pi,4)*pow(mrho,2)*(pow(mrho,2) + s*(3 - 4*C4*s)) +
                           pow(ma1,4)*(16*C4*pow(m_pi,6) + 2*pow(m_pi,4)*(2 + 5*C4*pow(mrho,2) - 14*C4*s) + 2*pow(m_pi,2)*s*(-1 - 7*C4*pow(mrho,2) + 7*C4*s) +
                              s*(-(s*(1 + 2*C4*s)) + pow(mrho,2)*(1 + 4*C4*s))) -
                           pow(ma1,2)*(8*C4*pow(m_pi,8) + pow(mrho,2)*(pow(mrho,2) - s)*s + 2*pow(m_pi,6)*(2 + 7*C4*pow(mrho,2) - 8*C4*s) +
                              pow(m_pi,2)*(-pow(mrho,4) + pow(s,2) + 8*C4*pow(mrho,2)*pow(s,2) - 2*C4*pow(s,3)) +
                              pow(m_pi,4)*(pow(mrho,2)*(3 - 22*C4*s) + 2*s*(-3 + 5*C4*s))))))*log(fabs(-pow(ma1,2) + t2)))/
                 ((pow(ma1,2) - pow(m_pi,2))*pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (16*pow(-2 + delta,2)*pow(m_pi,2)*log(fabs(-pow(m_pi,2) + t2)))/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) -
                (8*pow(-2 + delta,2)*(3*pow(m_pi,4) - 4*pow(m_pi,2)*(pow(mrho,2) - s) + pow(pow(mrho,2) - s,2))*log(fabs(-pow(m_pi,2) + t2)))/
                 ((pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (8*(-2 + delta)*(eta1 - eta2)*pow(m_pi,2)*(2*eta1*pow(m_pi,2) - 2*eta2*pow(m_pi,2) - eta1*pow(mrho,2))*(pow(m_pi,2) - s)*
                   log(fabs(-pow(m_pi,2) + t2)))/((pow(ma1,2) - pow(m_pi,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (8*(-2 + delta)*(eta1 - eta2)*pow(m_pi,2)*(pow(ma1,2) - s)*(pow(m_pi,2) - s)*
                   (-(eta2*(pow(m_pi,2) + s)) + eta1*(pow(m_pi,2) - pow(mrho,2) + s))*log(fabs(-pow(m_pi,2) + t2)))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (8*(-2 + delta)*(-(delta*(4*pow(m_pi,2) - pow(mrho,2))*(pow(m_pi,2) + pow(mrho,2) - s)) +
                     2*pow(mrho,2)*(8*C4*pow(m_pi,4) - pow(mrho,2) + s + pow(m_pi,2)*(3 - 8*C4*s)))*log(fabs(-pow(m_pi,2) + t2)))/
                 (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (2*pow(eta1 - eta2,2)*(pow(eta1,2)*(2*pow(ma1,6) - 3*pow(ma1,4)*(2*pow(m_pi,2) + pow(mrho,2) - s) + pow(mrho,2)*(pow(mrho,2) - s)*s -
                        pow(m_pi,4)*(3*pow(mrho,2) + s) + pow(m_pi,2)*(-2*pow(mrho,4) + 3*pow(mrho,2)*s) +
                        pow(ma1,2)*(4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) - 4*pow(mrho,2)*s + 2*pow(s,2))) -
                     2*eta1*eta2*(2*pow(ma1,6) - pow(m_pi,2)*(pow(m_pi,2) - pow(mrho,2))*(pow(mrho,2) + s) + pow(ma1,4)*(-6*pow(m_pi,2) + 3*s) +
                        pow(ma1,2)*(4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2))) +
                     pow(eta2,2)*(2*pow(ma1,6) + 3*pow(ma1,4)*(-2*pow(m_pi,2) + pow(mrho,2) + s) +
                        pow(m_pi,2)*(-pow(mrho,4) + pow(m_pi,2)*(2*pow(mrho,2) - s) + pow(mrho,2)*s) +
                        pow(ma1,2)*(4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s))))*log(fabs(-pow(ma1,2) + t1)))
                  /(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) +
                (2*pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(-2*eta1*eta2*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + 2*pow(ma1,2)*pow(m_pi,4)*s +
                        pow(m_pi,2)*(pow(ma1,4)*(pow(mrho,2) - 4*s) + 4*pow(ma1,2)*(pow(mrho,2) - s)*s + pow(mrho,2)*pow(s,2)) +
                        pow(ma1,2)*s*(pow(ma1,4) + s*(-2*pow(mrho,2) + s) + pow(ma1,2)*(-2*pow(mrho,2) + 3*s))) +
                     pow(eta2,2)*(pow(m_pi,8) - 4*pow(ma1,2)*pow(m_pi,2)*s*(pow(ma1,2) + pow(mrho,2) + s) +
                        pow(m_pi,4)*(pow(mrho,2)*s + pow(ma1,2)*(pow(mrho,2) + 2*s)) +
                        pow(ma1,2)*s*(pow(ma1,4) + s*(pow(mrho,2) + s) + pow(ma1,2)*(pow(mrho,2) + 3*s))) +
                     pow(eta1,2)*(pow(m_pi,8) + pow(ma1,2)*s*(pow(ma1,4) + 2*pow(mrho,4) - 3*pow(ma1,2)*(pow(mrho,2) - s) - 3*pow(mrho,2)*s + pow(s,2)) +
                        pow(m_pi,4)*(2*pow(mrho,4) - 3*pow(mrho,2)*s + pow(ma1,2)*(-3*pow(mrho,2) + 2*s)) +
                        2*pow(m_pi,2)*(pow(ma1,4)*(pow(mrho,2) - 2*s) + pow(mrho,2)*s*(-pow(mrho,2) + s) -
                           pow(ma1,2)*(pow(mrho,4) - 4*pow(mrho,2)*s + 2*pow(s,2)))))*log(fabs(-pow(ma1,2) + t1)))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (8*(eta1 - eta2)*(delta*(eta2*(pow(m_pi,6)*pow(mrho,2)*(2*pow(m_pi,2) - s) + pow(ma1,8)*(-pow(m_pi,2) + s) +
                           pow(ma1,6)*(5*pow(m_pi,4) - 7*pow(m_pi,2)*s + s*(pow(mrho,2) + 2*s)) +
                           pow(ma1,4)*(-8*pow(m_pi,6) - pow(m_pi,4)*(pow(mrho,2) - 14*s) + pow(m_pi,2)*(2*pow(mrho,4) - pow(mrho,2)*s - 7*pow(s,2)) +
                              s*(-2*pow(mrho,4) + pow(mrho,2)*s + pow(s,2))) +
                           pow(ma1,2)*pow(m_pi,2)*(4*pow(m_pi,6) + pow(m_pi,4)*(pow(mrho,2) - 8*s) + s*(2*pow(mrho,4) + pow(mrho,2)*s - pow(s,2)) +
                              pow(m_pi,2)*(-2*pow(mrho,4) - 3*pow(mrho,2)*s + 5*pow(s,2)))) +
                        eta1*(pow(ma1,8)*(pow(m_pi,2) - s) + pow(ma1,6)*(-5*pow(m_pi,4) + (pow(mrho,2) - 2*s)*s + pow(m_pi,2)*(-2*pow(mrho,2) + 7*s)) +
                           pow(m_pi,2)*pow(mrho,2)*(2*pow(m_pi,6) + pow(m_pi,4)*(4*pow(mrho,2) - 5*s) + pow(mrho,4)*s -
                              pow(m_pi,2)*(pow(mrho,4) + 3*pow(mrho,2)*s - 2*pow(s,2))) +
                           pow(ma1,4)*(8*pow(m_pi,6) + pow(m_pi,4)*(9*pow(mrho,2) - 14*s) + pow(m_pi,2)*s*(-9*pow(mrho,2) + 7*s) +
                              s*(pow(mrho,4) + pow(mrho,2)*s - pow(s,2))) +
                           pow(ma1,2)*(-4*pow(m_pi,8) + pow(mrho,4)*s*(-pow(mrho,2) + s) + pow(m_pi,6)*(-11*pow(mrho,2) + 8*s) +
                              pow(m_pi,4)*(-3*pow(mrho,4) + 17*pow(mrho,2)*s - 5*pow(s,2)) + pow(m_pi,2)*(pow(mrho,6) - 5*pow(mrho,2)*pow(s,2) + pow(s,3))
                              ))) - 2*pow(mrho,2)*(eta2*(pow(m_pi,8)*(1 + 2*C4*pow(mrho,2)) - 2*C4*pow(m_pi,6)*pow(mrho,2)*s +
                           2*C4*pow(ma1,8)*(-pow(m_pi,2) + s) + pow(ma1,4)*
                            (-16*C4*pow(m_pi,6) + pow(m_pi,4)*(-4 + 6*C4*pow(mrho,2) + 28*C4*s) +
                              2*pow(m_pi,2)*(pow(mrho,2) + s - 3*C4*pow(mrho,2)*s - 7*C4*pow(s,2)) + s*(-2*pow(mrho,2) + s + 2*C4*pow(s,2))) +
                           pow(ma1,6)*(10*C4*pow(m_pi,4) + 2*C4*s*(pow(mrho,2) + 2*s) + pow(m_pi,2)*(1 - 2*C4*(pow(mrho,2) + 7*s))) +
                           pow(ma1,2)*pow(m_pi,2)*(8*C4*pow(m_pi,6) - 2*pow(m_pi,4)*(-2 + 3*C4*pow(mrho,2) + 8*C4*s) +
                              s*(2*pow(mrho,2) + s - 2*C4*pow(s,2)) + 2*pow(m_pi,2)*(pow(mrho,2)*(-1 + 3*C4*s) + s*(-3 + 5*C4*s)))) +
                        eta1*(pow(m_pi,8)*(-1 + 6*C4*pow(mrho,2)) + 2*C4*pow(ma1,8)*(pow(m_pi,2) - s) + pow(m_pi,2)*pow(mrho,4)*s +
                           2*pow(m_pi,6)*pow(mrho,2)*(2 - 5*C4*s) - pow(ma1,6)*
                            (10*C4*pow(m_pi,4) + pow(m_pi,2)*(1 + 2*C4*(pow(mrho,2) - 7*s)) + 2*C4*s*(-pow(mrho,2) + 2*s)) -
                           pow(m_pi,4)*pow(mrho,2)*(pow(mrho,2) + s*(3 - 4*C4*s)) +
                           pow(ma1,4)*(16*C4*pow(m_pi,6) + 2*pow(m_pi,4)*(2 + 5*C4*pow(mrho,2) - 14*C4*s) + 2*pow(m_pi,2)*s*(-1 - 7*C4*pow(mrho,2) + 7*C4*s) +
                              s*(-(s*(1 + 2*C4*s)) + pow(mrho,2)*(1 + 4*C4*s))) -
                           pow(ma1,2)*(8*C4*pow(m_pi,8) + pow(mrho,2)*(pow(mrho,2) - s)*s + 2*pow(m_pi,6)*(2 + 7*C4*pow(mrho,2) - 8*C4*s) +
                              pow(m_pi,2)*(-pow(mrho,4) + pow(s,2) + 8*C4*pow(mrho,2)*pow(s,2) - 2*C4*pow(s,3)) +
                              pow(m_pi,4)*(pow(mrho,2)*(3 - 22*C4*s) + 2*s*(-3 + 5*C4*s))))))*log(fabs(-pow(ma1,2) + t1)))/
                 ((pow(ma1,2) - pow(m_pi,2))*pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (16*pow(-2 + delta,2)*pow(m_pi,2)*log(fabs(-pow(m_pi,2) + t1)))/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) +
                (8*pow(-2 + delta,2)*(3*pow(m_pi,4) - 4*pow(m_pi,2)*(pow(mrho,2) - s) + pow(pow(mrho,2) - s,2))*log(fabs(-pow(m_pi,2) + t1)))/
                 ((pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (8*(-2 + delta)*(eta1 - eta2)*pow(m_pi,2)*(2*eta1*pow(m_pi,2) - 2*eta2*pow(m_pi,2) - eta1*pow(mrho,2))*(pow(m_pi,2) - s)*
                   log(fabs(-pow(m_pi,2) + t1)))/((pow(ma1,2) - pow(m_pi,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                (8*(-2 + delta)*(eta1 - eta2)*pow(m_pi,2)*(pow(ma1,2) - s)*(pow(m_pi,2) - s)*
                   (-(eta2*(pow(m_pi,2) + s)) + eta1*(pow(m_pi,2) - pow(mrho,2) + s))*log(fabs(-pow(m_pi,2) + t1)))/
                 ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                (8*(-2 + delta)*(delta*(4*pow(m_pi,2) - pow(mrho,2))*(pow(m_pi,2) + pow(mrho,2) - s) -
                     2*pow(mrho,2)*(8*C4*pow(m_pi,4) - pow(mrho,2) + s + pow(m_pi,2)*(3 - 8*C4*s)))*log(fabs(-pow(m_pi,2) + t1)))/
                 (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)))))/(512.*Pi);

          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;

        case ReactionType::pi_rho:
          part_out = pi0_particle;
          m3 = part_out->mass();

          mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
          t1 = mandelstam_t[1];
          t2 = mandelstam_t[0];

          xsection = to_mb*(1/3.0*(pow(Const,2)*pow(ghat,4)*((0.03125*pow(eta1 - 1.*eta2,2)*
                    (pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
                         pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                         pow(ma1,2)*pow(mpion,2)*(-2.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) + 2.*pow(mrho,2)*s) +
                         pow(ma1,4)*(4.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2)))\
                       + eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
                         pow(ma1,2)*pow(mpion,2)*(-4.*pow(mrho,4) - 4.*pow(mrho,2)*s + pow(mpion,2)*(4.*pow(mrho,2) + 4.*s)) +
                         pow(ma1,4)*(-8.*pow(mpion,4) + 2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 4.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)))
                        + pow(eta1,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) - 2.*pow(mpion,2)*pow(mrho,4)*s +
                         pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) + pow(mpion,4)*(3.*pow(mrho,4) + 2.*pow(mrho,2)*s) +
                         pow(ma1,4)*(4.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2)) +
                         pow(ma1,2)*(pow(mpion,4)*(-6.*pow(mrho,2) - 2.*s) + pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s +
                            pow(mpion,2)*(-4.*pow(mrho,4) + 6.*pow(mrho,2)*s)))))/(1.*pow(ma1,2) - 1.*t2) +
                 (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*t2) +
                 (pow(2. - 1.*delta,2)*(0.25*pow(mpion,4) - 0.12500000000000003*pow(mrho,4) +
                      pow(mpion,2)*(0.5000000000000001*pow(mrho,2) - 0.5000000000000001*s) + 0.25*pow(s,2)))/(1.*pow(mpion,2) - 0.5*s - 0.5*t2) -
                 0.5*t2 + 1.*delta*t2 - 0.375*pow(delta,2)*t2 - 2.*C4*pow(mpion,2)*t2 + 1.*C4*delta*pow(mpion,2)*t2 +
                 (0.5*delta*pow(mpion,4)*t2)/pow(mrho,4) - (0.25*pow(delta,2)*pow(mpion,4)*t2)/pow(mrho,4) +
                 (1.5*pow(mpion,2)*t2)/pow(mrho,2) - (0.25*pow(-2. + delta,2)*pow(mpion,2)*t2)/pow(mrho,2) -
                 (1.*delta*pow(mpion,2)*t2)/pow(mrho,2) + (0.125*pow(delta,2)*pow(mpion,2)*t2)/pow(mrho,2) - 2.*C4*pow(mrho,2)*t2 +
                 1.*C4*delta*pow(mrho,2)*t2 + 4.*C4*s*t2 - 2.*C4*delta*s*t2 - (0.25*delta*pow(mpion,2)*s*t2)/pow(mrho,4) +
                 (0.125*pow(delta,2)*pow(mpion,2)*s*t2)/pow(mrho,4) - (0.75*delta*s*t2)/pow(mrho,2) +
                 (0.375*pow(delta,2)*s*t2)/pow(mrho,2) + 0.125*(-2. + delta)*(eta1 - 1.*eta2)*
                  (eta2*(-1.*pow(ma1,2) + pow(mrho,2) - 2.*s) + eta1*(pow(ma1,2) - 1.*pow(mpion,2) - 2.*pow(mrho,2) + s))*t2 +
                 0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(3.*pow(ma1,4) + 4.*pow(mpion,4) + pow(mrho,4) +
                       pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2) +
                       pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
                    pow(eta2,2)*(3.*pow(ma1,4) + 4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                       2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) + 4.*s)) +
                    eta1*eta2*(-6.*pow(ma1,4) - 8.*pow(mpion,4) + 2.*pow(mrho,4) + pow(ma1,2)*(16.*pow(mpion,2) - 8.*s) + 4.*pow(mrho,2)*s -
                       4.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)))*t2 -
                 (0.5000000625000001*(pow(mrho,4)*(1.*pow(mpion,2) - 1.499996062500492*pow(mrho,2) + 0.4999994375000703*s) +
                      delta*pow(mrho,2)*(-2.0000022499997185*pow(mpion,4) - 2.0000022499997185*pow(mpion,2)*pow(mrho,2) -
                         0.4999994375000703*pow(mrho,4) + 2.0000022499997185*pow(mrho,2)*s + 0.4999994375000703*pow(s,2)) +
                      pow(delta,2)*(1.*pow(mpion,6) + 1.1249998593750175*pow(mrho,6) +
                         pow(mpion,4)*(0.4999994375000703*pow(mrho,2) - 0.4999994375000703*s) - 1.6249947968756502*pow(mrho,4)*s -
                         0.6250004218749473*pow(mrho,2)*pow(s,2) + 0.12499985937501758*pow(s,3) +
                         pow(mpion,2)*(-0.24999971875003515*pow(mrho,4) + 0.4999994375000703*pow(mrho,2)*s - 0.24999971875003515*pow(s,2))))*t2)/
                  pow(mrho,6) - (0.0624999375*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*
                       (eta1*(2.0000030000030002*pow(ma1,2) + 2.0000030000030002*pow(mrho,2)) +
                         eta2*(-2.0000030000030002*pow(ma1,2) + 2.0000030000030002*pow(mpion,2) - 8.000018000018*pow(mrho,2) + 6.0000060000060005*s)) +
                      delta*(eta1*(1.*pow(ma1,4) - 2.0000030000030002*pow(mpion,4) - 3.0000030000030002*pow(mrho,4) +
                            pow(mpion,2)*(3.999993999994*pow(mrho,2) - 3.999993999994*s) - 2.0000030000030002*pow(mrho,2)*s +
                            5.000015000015001*pow(s,2) + pow(ma1,2)*(-2.0000030000030002*pow(mpion,2) + 1.*s)) +
                         eta2*(-1.*pow(ma1,4) - 3.999993999994*pow(mpion,4) + 3.999993999994*pow(mrho,4) +
                            pow(mpion,2)*(-1.*pow(mrho,2) - 2.0000030000030002*s) + 1.*pow(mrho,2)*s - 1.*pow(s,2) +
                            pow(ma1,2)*(3.0000030000030002*pow(mpion,2) - 3.0000030000030002*pow(mrho,2) + 2.0000030000030002*s))))*t2)/pow(mrho,2)\
                  - (0.5*(pow(mrho,6)*(-1.5 + C4*(-12.*pow(mpion,2) + 6.*s) + pow(C4,2)*(-16.*pow(mpion,4) + 16.*pow(mpion,2)*s - 4.*pow(s,2))) +
                      pow(delta,2)*(1.*pow(mpion,6) - 2.*pow(mpion,4)*pow(mrho,2) + 0.125*pow(mrho,6) + 0.25*pow(mrho,4)*s -
                         0.875*pow(mrho,2)*pow(s,2) + 0.25*pow(s,3) + pow(mpion,2)*(-2.5*pow(mrho,4) + 2.25*pow(mrho,2)*s - 0.75*pow(s,2))) +
                      delta*pow(mrho,2)*(pow(mpion,4)*(1. + 8.*C4*pow(mrho,2)) + 0.5*pow(s,2) + pow(mrho,4)*(1.5 - 5.*C4*s) +
                         pow(mrho,2)*s*(-0.5 + 1.*C4*s) + pow(mpion,2)*(6.*C4*pow(mrho,4) - 1.5*s + pow(mrho,2)*(3. - 6.*C4*s))))*t2)/pow(mrho,6)\
                  - (0.5*(-4.*C4*pow(mrho,8) - 0.5*pow(mrho,4)*s + pow(mrho,6)*(2. + 2.*C4*s) +
                      pow(delta,2)*(-2.*pow(mpion,6) - 2.*pow(mrho,6) + 0.5*pow(mpion,4)*s + 2.125*pow(mrho,4)*s + 1.25*pow(mrho,2)*pow(s,2) -
                         0.375*pow(s,3) + pow(mpion,2)*(1.5*pow(mrho,4) - 1.5*pow(mrho,2)*s + 1.*pow(s,2))) +
                      delta*pow(mrho,2)*(2.*pow(mpion,4) + 2.*C4*pow(mrho,6) - 1.*pow(s,2) + pow(mrho,2)*s*(-3. + 1.*C4*s) +
                         pow(mrho,4)*(1. + 1.*C4*s) + pow(mpion,2)*(1.*s + pow(mrho,2)*(1. - 2.*C4*s))))*t2)/pow(mrho,6) +
                 (0.0625*(1.*eta1 - 1.*eta2)*(delta*(eta1*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(18.*pow(mrho,2) - 12.*s) -
                            8.*pow(mrho,2)*s + 7.*pow(s,2) + pow(ma1,2)*(-10.*pow(mpion,2) - 4.*pow(mrho,2) + 5.*s)) +
                         eta2*(-3.*pow(ma1,4) - 12.*pow(mpion,4) + 2.*pow(mrho,4) + pow(ma1,2)*(11.*pow(mpion,2) - 3.*pow(mrho,2) - 2.*s) +
                            5.*pow(mrho,2)*s - 3.*pow(s,2) + pow(mpion,2)*(-1.*pow(mrho,2) + 6.*s))) +
                      pow(mrho,2)*(eta1*(-8.*C4*pow(ma1,4) - 32.*C4*pow(mpion,4) - 6.*pow(mrho,2) +
                            pow(ma1,2)*(6. + C4*(32.*pow(mpion,2) + 8.*pow(mrho,2) - 16.*s)) + 4.*s + 16.*C4*pow(mrho,2)*s - 8.*C4*pow(s,2) +
                            pow(mpion,2)*(-12. - 32.*C4*pow(mrho,2) + 32.*C4*s)) +
                         eta2*(8.*C4*pow(ma1,4) + 32.*C4*pow(mpion,4) - 4.*pow(mrho,2) - 2.*s + 8.*C4*pow(s,2) +
                            pow(mpion,2)*(10. - 16.*C4*pow(mrho,2) - 32.*C4*s) + pow(ma1,2)*(-6. + C4*(-32.*pow(mpion,2) + 8.*pow(mrho,2) + 16.*s))))
                      )*t2)/pow(mrho,2) + 0.0625*(-2. + delta)*pow(eta1 - 1.*eta2,2)*pow(t2,2) +
                 ((0.25*pow(mrho,2) + 1.*C4*pow(mrho,4) + pow(delta,2)*(0.1875*pow(mpion,2) + 0.25*pow(mrho,2) - 0.0625*s) +
                      delta*(-0.375*pow(mpion,2) - 0.625*pow(mrho,2) - 0.5*C4*pow(mrho,4) + 0.125*s))*pow(t2,2))/pow(mrho,4) +
                 0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                    eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(t2,2) -
                 (0.375*(0.3333333333333333*pow(mrho,4) - 1.3333333333333333*C4*pow(mrho,6) +
                      delta*pow(mrho,2)*(1.3333333333333333*pow(mrho,2) - 0.6666666666666666*C4*pow(mrho,4) +
                         pow(mpion,2)*(-0.6666666666666666 + 1.3333333333333333*C4*pow(mrho,2)) - 0.6666666666666666*s) +
                      pow(delta,2)*(1.*pow(mpion,4) + 0.25*pow(mrho,4) + pow(mpion,2)*(-0.3333333333333333*pow(mrho,2) + 0.6666666666666666*s) -
                         0.5833333333333334*pow(s,2)))*pow(t2,2))/pow(mrho,6) +
                 ((-0.125*pow(mrho,4) + delta*(0.5*pow(mrho,4) - 0.25*pow(mrho,2)*s) +
                      pow(delta,2)*(0.125*pow(mpion,4) - 0.09375*pow(mrho,4) + pow(mpion,2)*(-0.125*pow(mrho,2) + 0.125*s) +
                         0.0625*pow(mrho,2)*s - 0.09375*pow(s,2)))*pow(t2,2))/pow(mrho,6) -
                 (0.03125*(1.*eta1 - 1.*eta2)*((2.*eta1 - 2.*eta2)*pow(mrho,2) +
                      delta*(eta1*(1.*pow(ma1,2) - 2.*pow(mpion,2) + 1.*s) + eta2*(-1.*pow(ma1,2) + 3.*pow(mpion,2) - 3.*pow(mrho,2) + 2.*s)))*
                    pow(t2,2))/pow(mrho,2) + (0.25*(C4*pow(mrho,6)*(-6. - 16.*C4*pow(mpion,2) + 8.*C4*s) +
                      pow(delta,2)*(1.*pow(mpion,4) - 0.25*pow(mrho,4) + pow(mpion,2)*(-1.75*pow(mrho,2) + 0.5*s) + 0.5*pow(mrho,2)*s -
                         0.5*pow(s,2)) + delta*pow(mrho,2)*(1.*C4*pow(mrho,4) + pow(mpion,2)*(0.5 + 10.*C4*pow(mrho,2)) - 0.5*s +
                         pow(mrho,2)*(2.5 - 4.*C4*s)))*pow(t2,2))/pow(mrho,6) +
                 (0.09375*(1.*eta1 - 1.*eta2)*(delta*(eta2*(-1.*pow(ma1,2) + 3.6666666666666665*pow(mpion,2) - 1.*pow(mrho,2) - 0.6666666666666666*s) +
                         eta1*(1.*pow(ma1,2) - 3.3333333333333335*pow(mpion,2) - 1.3333333333333333*pow(mrho,2) + 1.6666666666666667*s)) +
                      pow(mrho,2)*(eta1*(2. + C4*(-2.6666666666666665*pow(ma1,2) + 10.666666666666666*pow(mpion,2) + 2.6666666666666665*pow(mrho,2) -
                               5.333333333333333*s)) + eta2*(-2. + C4*
                             (2.6666666666666665*pow(ma1,2) - 10.666666666666666*pow(mpion,2) + 2.6666666666666665*pow(mrho,2) + 5.333333333333333*s))))*
                    pow(t2,2))/pow(mrho,2) + 0.666667*pow(C4,2)*pow(t2,3) + 0.0104167*pow(eta1 - 1.*eta2,4)*pow(t2,3) -
                 (0.0416667*pow(delta,2)*pow(mpion,2)*pow(t2,3))/pow(mrho,6) - (0.0208333*pow(delta,2)*pow(t2,3))/pow(mrho,4) -
                 (0.0416667*delta*(-2. + 1.*delta)*pow(t2,3))/pow(mrho,4) - (0.333333*C4*delta*pow(t2,3))/pow(mrho,2) -
                 (0.0208333*delta*pow(1.*eta1 - 1.*eta2,2)*pow(t2,3))/pow(mrho,2) -
                 (0.166667*pow(1.*eta1 - 1.*eta2,2)*(-0.375*delta + 1.*C4*pow(mrho,2))*pow(t2,3))/pow(mrho,2) +
                 (0.104167*delta*(-0.8*pow(mrho,2) + delta*(0.4*pow(mpion,2) + 1.*pow(mrho,2) - 0.6*s))*pow(t2,3))/pow(mrho,6) +
                 (0.0625*pow(delta,2)*s*pow(t2,3))/pow(mrho,6) -
                 (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
                         1.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                         pow(ma1,2)*pow(mpion,2)*(-2.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) + 2.*pow(mrho,2)*s) +
                         pow(ma1,4)*(4.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2)))\
                       + eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
                         pow(ma1,2)*pow(mpion,2)*(-4.*pow(mrho,4) - 4.*pow(mrho,2)*s + pow(mpion,2)*(4.*pow(mrho,2) + 4.*s)) +
                         pow(ma1,4)*(-8.*pow(mpion,4) + 2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 4.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)))
                        + pow(eta1,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) - 2.*pow(mpion,2)*pow(mrho,4)*s +
                         pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) + pow(mpion,4)*(3.*pow(mrho,4) + 2.*pow(mrho,2)*s) +
                         pow(ma1,4)*(4.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2)) +
                         pow(ma1,2)*(pow(mpion,4)*(-6.*pow(mrho,2) - 2.*s) + pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s +
                            pow(mpion,2)*(-4.*pow(mrho,4) + 6.*pow(mrho,2)*s)))))/(1.*pow(ma1,2) - 1.*t1) -
                 (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*t1) +
                 (pow(2. - 1.*delta,2)*(-0.25*pow(mpion,4) + 0.12500000000000003*pow(mrho,4) +
                      pow(mpion,2)*(-0.5000000000000001*pow(mrho,2) + 0.5000000000000001*s) - 0.25*pow(s,2)))/(1.*pow(mpion,2) - 0.5*s - 0.5*t1) +
                 0.5*t1 - 1.*delta*t1 + 0.375*pow(delta,2)*t1 + 2.*C4*pow(mpion,2)*t1 - 1.*C4*delta*pow(mpion,2)*t1 -
                 (0.5*delta*pow(mpion,4)*t1)/pow(mrho,4) + (0.25*pow(delta,2)*pow(mpion,4)*t1)/pow(mrho,4) -
                 (1.5*pow(mpion,2)*t1)/pow(mrho,2) + (0.25*pow(-2. + delta,2)*pow(mpion,2)*t1)/pow(mrho,2) +
                 (1.*delta*pow(mpion,2)*t1)/pow(mrho,2) - (0.125*pow(delta,2)*pow(mpion,2)*t1)/pow(mrho,2) + 2.*C4*pow(mrho,2)*t1 -
                 1.*C4*delta*pow(mrho,2)*t1 - 4.*C4*s*t1 + 2.*C4*delta*s*t1 + (0.25*delta*pow(mpion,2)*s*t1)/pow(mrho,4) -
                 (0.125*pow(delta,2)*pow(mpion,2)*s*t1)/pow(mrho,4) + (0.75*delta*s*t1)/pow(mrho,2) -
                 (0.375*pow(delta,2)*s*t1)/pow(mrho,2) - 0.125*(-2. + delta)*(eta1 - 1.*eta2)*
                  (eta2*(-1.*pow(ma1,2) + pow(mrho,2) - 2.*s) + eta1*(pow(ma1,2) - 1.*pow(mpion,2) - 2.*pow(mrho,2) + s))*t1 -
                 0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(3.*pow(ma1,4) + 4.*pow(mpion,4) + pow(mrho,4) +
                       pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2) +
                       pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
                    pow(eta2,2)*(3.*pow(ma1,4) + 4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                       2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) + 4.*s)) +
                    eta1*eta2*(-6.*pow(ma1,4) - 8.*pow(mpion,4) + 2.*pow(mrho,4) + pow(ma1,2)*(16.*pow(mpion,2) - 8.*s) + 4.*pow(mrho,2)*s -
                       4.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)))*t1 +
                 (0.5000000625000001*(pow(mrho,4)*(1.*pow(mpion,2) - 1.499996062500492*pow(mrho,2) + 0.4999994375000703*s) +
                      delta*pow(mrho,2)*(-2.0000022499997185*pow(mpion,4) - 2.0000022499997185*pow(mpion,2)*pow(mrho,2) -
                         0.4999994375000703*pow(mrho,4) + 2.0000022499997185*pow(mrho,2)*s + 0.4999994375000703*pow(s,2)) +
                      pow(delta,2)*(1.*pow(mpion,6) + 1.1249998593750175*pow(mrho,6) +
                         pow(mpion,4)*(0.4999994375000703*pow(mrho,2) - 0.4999994375000703*s) - 1.6249947968756502*pow(mrho,4)*s -
                         0.6250004218749473*pow(mrho,2)*pow(s,2) + 0.12499985937501758*pow(s,3) +
                         pow(mpion,2)*(-0.24999971875003515*pow(mrho,4) + 0.4999994375000703*pow(mrho,2)*s - 0.24999971875003515*pow(s,2))))*t1)/
                  pow(mrho,6) + (0.0624999375*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*
                       (eta1*(2.0000030000030002*pow(ma1,2) + 2.0000030000030002*pow(mrho,2)) +
                         eta2*(-2.0000030000030002*pow(ma1,2) + 2.0000030000030002*pow(mpion,2) - 8.000018000018*pow(mrho,2) + 6.0000060000060005*s)) +
                      delta*(eta1*(1.*pow(ma1,4) - 2.0000030000030002*pow(mpion,4) - 3.0000030000030002*pow(mrho,4) +
                            pow(mpion,2)*(3.999993999994*pow(mrho,2) - 3.999993999994*s) - 2.0000030000030002*pow(mrho,2)*s +
                            5.000015000015001*pow(s,2) + pow(ma1,2)*(-2.0000030000030002*pow(mpion,2) + 1.*s)) +
                         eta2*(-1.*pow(ma1,4) - 3.999993999994*pow(mpion,4) + 3.999993999994*pow(mrho,4) +
                            pow(mpion,2)*(-1.*pow(mrho,2) - 2.0000030000030002*s) + 1.*pow(mrho,2)*s - 1.*pow(s,2) +
                            pow(ma1,2)*(3.0000030000030002*pow(mpion,2) - 3.0000030000030002*pow(mrho,2) + 2.0000030000030002*s))))*t1)/pow(mrho,2)\
                  + (0.5*(pow(mrho,6)*(-1.5 + C4*(-12.*pow(mpion,2) + 6.*s) + pow(C4,2)*(-16.*pow(mpion,4) + 16.*pow(mpion,2)*s - 4.*pow(s,2))) +
                      pow(delta,2)*(1.*pow(mpion,6) - 2.*pow(mpion,4)*pow(mrho,2) + 0.125*pow(mrho,6) + 0.25*pow(mrho,4)*s -
                         0.875*pow(mrho,2)*pow(s,2) + 0.25*pow(s,3) + pow(mpion,2)*(-2.5*pow(mrho,4) + 2.25*pow(mrho,2)*s - 0.75*pow(s,2))) +
                      delta*pow(mrho,2)*(pow(mpion,4)*(1. + 8.*C4*pow(mrho,2)) + 0.5*pow(s,2) + pow(mrho,4)*(1.5 - 5.*C4*s) +
                         pow(mrho,2)*s*(-0.5 + 1.*C4*s) + pow(mpion,2)*(6.*C4*pow(mrho,4) - 1.5*s + pow(mrho,2)*(3. - 6.*C4*s))))*t1)/pow(mrho,6)\
                  + (0.5*(-4.*C4*pow(mrho,8) - 0.5*pow(mrho,4)*s + pow(mrho,6)*(2. + 2.*C4*s) +
                      pow(delta,2)*(-2.*pow(mpion,6) - 2.*pow(mrho,6) + 0.5*pow(mpion,4)*s + 2.125*pow(mrho,4)*s + 1.25*pow(mrho,2)*pow(s,2) -
                         0.375*pow(s,3) + pow(mpion,2)*(1.5*pow(mrho,4) - 1.5*pow(mrho,2)*s + 1.*pow(s,2))) +
                      delta*pow(mrho,2)*(2.*pow(mpion,4) + 2.*C4*pow(mrho,6) - 1.*pow(s,2) + pow(mrho,2)*s*(-3. + 1.*C4*s) +
                         pow(mrho,4)*(1. + 1.*C4*s) + pow(mpion,2)*(1.*s + pow(mrho,2)*(1. - 2.*C4*s))))*t1)/pow(mrho,6) -
                 (0.0625*(1.*eta1 - 1.*eta2)*(delta*(eta1*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(18.*pow(mrho,2) - 12.*s) -
                            8.*pow(mrho,2)*s + 7.*pow(s,2) + pow(ma1,2)*(-10.*pow(mpion,2) - 4.*pow(mrho,2) + 5.*s)) +
                         eta2*(-3.*pow(ma1,4) - 12.*pow(mpion,4) + 2.*pow(mrho,4) + pow(ma1,2)*(11.*pow(mpion,2) - 3.*pow(mrho,2) - 2.*s) +
                            5.*pow(mrho,2)*s - 3.*pow(s,2) + pow(mpion,2)*(-1.*pow(mrho,2) + 6.*s))) +
                      pow(mrho,2)*(eta1*(-8.*C4*pow(ma1,4) - 32.*C4*pow(mpion,4) - 6.*pow(mrho,2) +
                            pow(ma1,2)*(6. + C4*(32.*pow(mpion,2) + 8.*pow(mrho,2) - 16.*s)) + 4.*s + 16.*C4*pow(mrho,2)*s - 8.*C4*pow(s,2) +
                            pow(mpion,2)*(-12. - 32.*C4*pow(mrho,2) + 32.*C4*s)) +
                         eta2*(8.*C4*pow(ma1,4) + 32.*C4*pow(mpion,4) - 4.*pow(mrho,2) - 2.*s + 8.*C4*pow(s,2) +
                            pow(mpion,2)*(10. - 16.*C4*pow(mrho,2) - 32.*C4*s) + pow(ma1,2)*(-6. + C4*(-32.*pow(mpion,2) + 8.*pow(mrho,2) + 16.*s))))
                      )*t1)/pow(mrho,2) - 0.0625*(-2. + delta)*pow(eta1 - 1.*eta2,2)*pow(t1,2) +
                 ((-0.25*pow(mrho,2) - 1.*C4*pow(mrho,4) + delta*(0.375*pow(mpion,2) + 0.625*pow(mrho,2) + 0.5*C4*pow(mrho,4) - 0.125*s) +
                      pow(delta,2)*(-0.1875*pow(mpion,2) - 0.25*pow(mrho,2) + 0.0625*s))*pow(t1,2))/pow(mrho,4) -
                 0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                    eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(t1,2) +
                 (0.375*(0.3333333333333333*pow(mrho,4) - 1.3333333333333333*C4*pow(mrho,6) +
                      delta*pow(mrho,2)*(1.3333333333333333*pow(mrho,2) - 0.6666666666666666*C4*pow(mrho,4) +
                         pow(mpion,2)*(-0.6666666666666666 + 1.3333333333333333*C4*pow(mrho,2)) - 0.6666666666666666*s) +
                      pow(delta,2)*(1.*pow(mpion,4) + 0.25*pow(mrho,4) + pow(mpion,2)*(-0.3333333333333333*pow(mrho,2) + 0.6666666666666666*s) -
                         0.5833333333333334*pow(s,2)))*pow(t1,2))/pow(mrho,6) +
                 ((0.125*pow(mrho,4) + delta*(-0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s) +
                      pow(delta,2)*(-0.125*pow(mpion,4) + 0.09375*pow(mrho,4) + pow(mpion,2)*(0.125*pow(mrho,2) - 0.125*s) -
                         0.0625*pow(mrho,2)*s + 0.09375*pow(s,2)))*pow(t1,2))/pow(mrho,6) +
                 (0.03125*(1.*eta1 - 1.*eta2)*((2.*eta1 - 2.*eta2)*pow(mrho,2) +
                      delta*(eta1*(1.*pow(ma1,2) - 2.*pow(mpion,2) + 1.*s) + eta2*(-1.*pow(ma1,2) + 3.*pow(mpion,2) - 3.*pow(mrho,2) + 2.*s)))*
                    pow(t1,2))/pow(mrho,2) - (0.25*(C4*pow(mrho,6)*(-6. - 16.*C4*pow(mpion,2) + 8.*C4*s) +
                      pow(delta,2)*(1.*pow(mpion,4) - 0.25*pow(mrho,4) + pow(mpion,2)*(-1.75*pow(mrho,2) + 0.5*s) + 0.5*pow(mrho,2)*s -
                         0.5*pow(s,2)) + delta*pow(mrho,2)*(1.*C4*pow(mrho,4) + pow(mpion,2)*(0.5 + 10.*C4*pow(mrho,2)) - 0.5*s +
                         pow(mrho,2)*(2.5 - 4.*C4*s)))*pow(t1,2))/pow(mrho,6) -
                 (0.09375*(1.*eta1 - 1.*eta2)*(delta*(eta2*(-1.*pow(ma1,2) + 3.6666666666666665*pow(mpion,2) - 1.*pow(mrho,2) - 0.6666666666666666*s) +
                         eta1*(1.*pow(ma1,2) - 3.3333333333333335*pow(mpion,2) - 1.3333333333333333*pow(mrho,2) + 1.6666666666666667*s)) +
                      pow(mrho,2)*(eta1*(2. + C4*(-2.6666666666666665*pow(ma1,2) + 10.666666666666666*pow(mpion,2) + 2.6666666666666665*pow(mrho,2) -
                               5.333333333333333*s)) + eta2*(-2. + C4*
                             (2.6666666666666665*pow(ma1,2) - 10.666666666666666*pow(mpion,2) + 2.6666666666666665*pow(mrho,2) + 5.333333333333333*s))))*
                    pow(t1,2))/pow(mrho,2) - 0.666667*pow(C4,2)*pow(t1,3) - 0.0104167*pow(eta1 - 1.*eta2,4)*pow(t1,3) +
                 (0.0416667*pow(delta,2)*pow(mpion,2)*pow(t1,3))/pow(mrho,6) + (0.0208333*pow(delta,2)*pow(t1,3))/pow(mrho,4) +
                 (0.0416667*delta*(-2. + 1.*delta)*pow(t1,3))/pow(mrho,4) + (0.333333*C4*delta*pow(t1,3))/pow(mrho,2) +
                 (0.0208333*delta*pow(1.*eta1 - 1.*eta2,2)*pow(t1,3))/pow(mrho,2) +
                 (0.166667*pow(1.*eta1 - 1.*eta2,2)*(-0.375*delta + 1.*C4*pow(mrho,2))*pow(t1,3))/pow(mrho,2) -
                 (0.104167*delta*(-0.8*pow(mrho,2) + delta*(0.4*pow(mpion,2) + 1.*pow(mrho,2) - 0.6*s))*pow(t1,3))/pow(mrho,6) -
                 (0.0625*pow(delta,2)*s*pow(t1,3))/pow(mrho,6) -
                 (0.0625*(1.*eta1 - 1.*eta2)*(2.*pow(mrho,2) + delta*(1.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s))*
                    (eta2*(-1.*pow(ma1,6) + pow(mpion,2)*(4.*pow(mpion,2) - 1.*s)*s + pow(ma1,4)*(3.*pow(mpion,2) - 4.*pow(mrho,2) + 2.*s) +
                         pow(ma1,2)*(-4.*pow(mpion,4) - 2.*pow(mpion,2)*s + (4.*pow(mrho,2) - 1.*s)*s)) +
                      eta1*(1.*pow(ma1,6) + 8.*pow(mpion,6) + pow(mpion,4)*(-4.*pow(mrho,2) - 6.*s) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s)*s +
                         pow(ma1,4)*(-2.*pow(mpion,2) + 1.*pow(mrho,2) + 1.*s) + s*(2.*pow(mrho,4) - 3.*pow(mrho,2)*s + 1.*pow(s,2)) +
                         pow(ma1,2)*(-2.*pow(mpion,4) - 2.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 5.*pow(s,2))))
                     *log(fabs(-pow(ma1,2) + t2)))/(pow(mrho,2)*(pow(ma1,2) - 2.*pow(mpion,2) + s)) +
                 (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,6) + pow(mpion,6) - 1.*pow(mpion,4)*pow(mrho,2) +
                         pow(ma1,4)*(pow(mpion,2) + pow(mrho,2) - 2.*s) + pow(ma1,2)*(3.*pow(mpion,4) - 2.*pow(mpion,2)*s)) +
                      eta1*(pow(ma1,6) + pow(ma1,4)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + s) +
                         pow(mpion,2)*(-4.*pow(mpion,4) - 1.*pow(mrho,4) - 1.*pow(mrho,2)*s + pow(mpion,2)*(3.*pow(mrho,2) + s)) +
                         pow(ma1,2)*(pow(mpion,4) + pow(mrho,4) - 1.*pow(mrho,2)*s + pow(mpion,2)*(pow(mrho,2) + 2.*s))))*
                    log(fabs(-pow(ma1,2) + t2)))/(pow(ma1,2) - 1.*pow(mpion,2)) +
                 0.0625*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(2.*pow(ma1,6) + pow(mpion,4)*(-3.*pow(mrho,2) - 1.*s) +
                       pow(mrho,2)*(pow(mrho,2) - 1.*s)*s + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
                       pow(mpion,2)*(-2.*pow(mrho,4) + 3.*pow(mrho,2)*s) +
                       pow(ma1,2)*(4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2))) +
                    pow(eta2,2)*(2.*pow(ma1,6) + pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
                       pow(mpion,2)*(-1.*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.*s) + pow(mrho,2)*s) +
                       pow(ma1,2)*(4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2))) +
                    eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
                       pow(mpion,2)*(-2.*pow(mrho,4) - 2.*pow(mrho,2)*s + pow(mpion,2)*(2.*pow(mrho,2) + 2.*s)) +
                       pow(ma1,2)*(-8.*pow(mpion,4) + 2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 4.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s))))*
                  log(fabs(-pow(ma1,2) + t2)) + (0.0625*(1.*eta1 - 1.*eta2)*
                    (delta*(eta1*(3.*pow(ma1,6) + 8.*pow(mpion,6) + pow(mpion,4)*(-12.*pow(mrho,2) - 6.*s) +
                            pow(ma1,4)*(-10.*pow(mpion,2) - 4.*pow(mrho,2) + 5.*s) +
                            pow(mpion,2)*(-4.*pow(mrho,4) + 10.*pow(mrho,2)*s - 2.*pow(s,2)) +
                            s*(3.*pow(mrho,4) - 4.*pow(mrho,2)*s + pow(s,2)) +
                            pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(18.*pow(mrho,2) - 12.*s) - 8.*pow(mrho,2)*s + 7.*pow(s,2)))
                           + eta2*(-3.*pow(ma1,6) + pow(ma1,4)*(11.*pow(mpion,2) - 3.*pow(mrho,2) - 2.*s) +
                            pow(mpion,2)*(-2.*pow(mrho,4) + pow(mrho,2)*s - 1.*pow(s,2) + pow(mpion,2)*(-2.*pow(mrho,2) + 4.*s)) +
                            pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 5.*pow(mrho,2)*s - 3.*pow(s,2) +
                               pow(mpion,2)*(-1.*pow(mrho,2) + 6.*s)))) +
                      pow(mrho,2)*(eta2*(8.*C4*pow(ma1,6) + pow(mpion,2)*((4. + 8.*C4*pow(mpion,2))*pow(mrho,2) - 2.*s) +
                            pow(ma1,4)*(-6. + C4*(-32.*pow(mpion,2) + 8.*pow(mrho,2) + 16.*s)) +
                            pow(ma1,2)*(32.*C4*pow(mpion,4) - 4.*pow(mrho,2) + pow(mpion,2)*(10. - 16.*C4*pow(mrho,2) - 32.*C4*s) +
                               s*(-2. + 8.*C4*s))) + eta1*(-8.*C4*pow(ma1,6) + pow(mpion,4)*(4. + 24.*C4*pow(mrho,2)) +
                            pow(ma1,4)*(6. + C4*(32.*pow(mpion,2) + 8.*pow(mrho,2) - 16.*s)) + s*(-2.*pow(mrho,2) + 2.*s) +
                            pow(mpion,2)*(-4.*s + pow(mrho,2)*(8. - 16.*C4*s)) +
                            pow(ma1,2)*(-32.*C4*pow(mpion,4) + s*(4. - 8.*C4*s) + pow(mrho,2)*(-6. + 16.*C4*s) +
                               pow(mpion,2)*(-12. - 32.*C4*pow(mrho,2) + 32.*C4*s)))))*log(fabs(-pow(ma1,2) + t2)))/pow(mrho,2) +
                 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) + t2)) +
                 (0.5*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(mpion,6) + 1.*pow(mpion,4)*s) +
                      eta1*(1.*pow(mpion,6) + pow(mpion,4)*(-0.5*pow(mrho,2) - 1.*s) + 0.5*pow(mpion,2)*pow(mrho,2)*s))*
                    log(fabs(-pow(mpion,2) + t2)))/(pow(ma1,2) - 1.*pow(mpion,2)) -
                 (0.5*(pow(delta,2)*(1.*pow(mpion,6) + 0.5*pow(mrho,6) + pow(mpion,4)*(-1.*pow(mrho,2) - 2.*s) - 0.5*pow(mrho,4)*s +
                         pow(mpion,2)*(-1.5*pow(mrho,4) + 1.*pow(mrho,2)*s + 1.*pow(s,2))) +
                      delta*(-2.*pow(mrho,6) + pow(mpion,6)*(-2. - 4.*C4*pow(mrho,2)) + 2.*pow(mrho,4)*s +
                         pow(mpion,2)*(6.*pow(mrho,4) - 2.*pow(s,2) + pow(mrho,2)*s*(-4. - 4.*C4*s)) +
                         pow(mpion,4)*(4.*s + pow(mrho,2)*(4. + 8.*C4*s))) +
                      pow(mrho,2)*(8.*C4*pow(mpion,6) + 2.*pow(mrho,4) - 2.*pow(mrho,2)*s + pow(mpion,4)*(-4. - 16.*C4*s) +
                         pow(mpion,2)*(-6.*pow(mrho,2) + s*(4. + 8.*C4*s))))*log(fabs(-pow(mpion,2) + t2)))/(pow(mrho,2)*(1.*pow(mpion,2) - 1.*s)) -
                 (0.25*(1.*eta1 - 1.*eta2)*(eta2*((2. - 1.*delta)*pow(mpion,6) + pow(mpion,2)*s*((-12. + 6.*delta)*pow(mrho,2) + (6. - 3.*delta)*s) +
                         pow(s,2)*((4. - 2.*delta)*pow(mrho,2) + (-2. + 1.*delta)*s) + pow(mpion,4)*((8. - 4.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s))
                        + eta1*((-2. + 1.*delta)*pow(mpion,6) + (-2. + 1.*delta)*pow(mrho,4)*s + (2. - 1.*delta)*pow(s,3) +
                         pow(mpion,4)*((-4. + 2.*delta)*pow(mrho,2) + (6. - 3.*delta)*s) +
                         pow(mpion,2)*((2. - 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-6. + 3.*delta)*pow(s,2))))*
                    log(fabs(-2.*pow(mpion,2) + s + t2)))/(pow(ma1,2) - 2.*pow(mpion,2) + s) -
                 (1.*pow(2. - 1.*delta,2)*(1.*pow(mpion,4) - 0.25*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s))*
                    log(fabs(-2.*pow(mpion,2) + 1.*s + 1.*t2)))/(1.*pow(mpion,2) - 1.*s) -
                 (0.9999974999999999*(pow(mrho,2)*(-2.00000750001875*pow(mpion,2) - 1.*pow(mrho,2) + 2.00000750001875*s) +
                      pow(delta,2)*(1.*pow(mpion,4) - 0.7500018750046875*pow(mrho,4) +
                         pow(mpion,2)*(1.500003750009375*pow(mrho,2) - 2.00000750001875*s) + 0.5000015000037501*pow(mrho,2)*s + 1.*pow(s,2)) +
                      delta*(-2.00000750001875*pow(mpion,4) + 2.00000750001875*pow(mrho,4) - 2.00000750001875*pow(mrho,2)*s -
                         2.00000750001875*pow(s,2) + pow(mpion,2)*(-2.00000750001875*pow(mrho,2) + 4.00000750001875*s)))*
                    log(fabs(-2.*pow(mpion,2) + 1.*s + 1.*t2)))/pow(mrho,2) +
                 (1.*(pow(mrho,2)*(pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) + (2. - 4.*C4*pow(mrho,2))*s) +
                      pow(delta,2)*(1.*pow(mpion,4) - 0.5*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 2.*s) + 1.*pow(s,2)) +
                      delta*(-2.*pow(mpion,4) - 1.*pow(mrho,2)*s - 2.*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) - 2.*C4*pow(mrho,4) + 4.*s) +
                         pow(mrho,4)*(1. + 2.*C4*s)))*log(fabs(-2.*pow(mpion,2) + 1.*s + 1.*t2)))/pow(mrho,2) +
                 (0.0625*(1.*eta1 - 1.*eta2)*(2.*pow(mrho,2) + delta*(1.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s))*
                    (eta2*(-1.*pow(ma1,6) + pow(mpion,2)*(4.*pow(mpion,2) - 1.*s)*s + pow(ma1,4)*(3.*pow(mpion,2) - 4.*pow(mrho,2) + 2.*s) +
                         pow(ma1,2)*(-4.*pow(mpion,4) - 2.*pow(mpion,2)*s + (4.*pow(mrho,2) - 1.*s)*s)) +
                      eta1*(1.*pow(ma1,6) + 8.*pow(mpion,6) + pow(mpion,4)*(-4.*pow(mrho,2) - 6.*s) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s)*s +
                         pow(ma1,4)*(-2.*pow(mpion,2) + 1.*pow(mrho,2) + 1.*s) + s*(2.*pow(mrho,4) - 3.*pow(mrho,2)*s + 1.*pow(s,2)) +
                         pow(ma1,2)*(-2.*pow(mpion,4) - 2.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 5.*pow(s,2))))
                     *log(fabs(-pow(ma1,2) + t1)))/(pow(mrho,2)*(pow(ma1,2) - 2.*pow(mpion,2) + s)) -
                 (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,6) + pow(mpion,6) - 1.*pow(mpion,4)*pow(mrho,2) +
                         pow(ma1,4)*(pow(mpion,2) + pow(mrho,2) - 2.*s) + pow(ma1,2)*(3.*pow(mpion,4) - 2.*pow(mpion,2)*s)) +
                      eta1*(pow(ma1,6) + pow(ma1,4)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + s) +
                         pow(mpion,2)*(-4.*pow(mpion,4) - 1.*pow(mrho,4) - 1.*pow(mrho,2)*s + pow(mpion,2)*(3.*pow(mrho,2) + s)) +
                         pow(ma1,2)*(pow(mpion,4) + pow(mrho,4) - 1.*pow(mrho,2)*s + pow(mpion,2)*(pow(mrho,2) + 2.*s))))*
                    log(fabs(-pow(ma1,2) + t1)))/(pow(ma1,2) - 1.*pow(mpion,2)) -
                 0.0625*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(2.*pow(ma1,6) + pow(mpion,4)*(-3.*pow(mrho,2) - 1.*s) +
                       pow(mrho,2)*(pow(mrho,2) - 1.*s)*s + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
                       pow(mpion,2)*(-2.*pow(mrho,4) + 3.*pow(mrho,2)*s) +
                       pow(ma1,2)*(4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2))) +
                    pow(eta2,2)*(2.*pow(ma1,6) + pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
                       pow(mpion,2)*(-1.*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.*s) + pow(mrho,2)*s) +
                       pow(ma1,2)*(4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2))) +
                    eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
                       pow(mpion,2)*(-2.*pow(mrho,4) - 2.*pow(mrho,2)*s + pow(mpion,2)*(2.*pow(mrho,2) + 2.*s)) +
                       pow(ma1,2)*(-8.*pow(mpion,4) + 2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 4.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s))))*
                  log(fabs(-pow(ma1,2) + t1)) - (0.0625*(1.*eta1 - 1.*eta2)*
                    (delta*(eta1*(3.*pow(ma1,6) + 8.*pow(mpion,6) + pow(mpion,4)*(-12.*pow(mrho,2) - 6.*s) +
                            pow(ma1,4)*(-10.*pow(mpion,2) - 4.*pow(mrho,2) + 5.*s) +
                            pow(mpion,2)*(-4.*pow(mrho,4) + 10.*pow(mrho,2)*s - 2.*pow(s,2)) +
                            s*(3.*pow(mrho,4) - 4.*pow(mrho,2)*s + pow(s,2)) +
                            pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(18.*pow(mrho,2) - 12.*s) - 8.*pow(mrho,2)*s + 7.*pow(s,2)))
                           + eta2*(-3.*pow(ma1,6) + pow(ma1,4)*(11.*pow(mpion,2) - 3.*pow(mrho,2) - 2.*s) +
                            pow(mpion,2)*(-2.*pow(mrho,4) + pow(mrho,2)*s - 1.*pow(s,2) + pow(mpion,2)*(-2.*pow(mrho,2) + 4.*s)) +
                            pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 5.*pow(mrho,2)*s - 3.*pow(s,2) +
                               pow(mpion,2)*(-1.*pow(mrho,2) + 6.*s)))) +
                      pow(mrho,2)*(eta2*(8.*C4*pow(ma1,6) + pow(mpion,2)*((4. + 8.*C4*pow(mpion,2))*pow(mrho,2) - 2.*s) +
                            pow(ma1,4)*(-6. + C4*(-32.*pow(mpion,2) + 8.*pow(mrho,2) + 16.*s)) +
                            pow(ma1,2)*(32.*C4*pow(mpion,4) - 4.*pow(mrho,2) + pow(mpion,2)*(10. - 16.*C4*pow(mrho,2) - 32.*C4*s) +
                               s*(-2. + 8.*C4*s))) + eta1*(-8.*C4*pow(ma1,6) + pow(mpion,4)*(4. + 24.*C4*pow(mrho,2)) +
                            pow(ma1,4)*(6. + C4*(32.*pow(mpion,2) + 8.*pow(mrho,2) - 16.*s)) + s*(-2.*pow(mrho,2) + 2.*s) +
                            pow(mpion,2)*(-4.*s + pow(mrho,2)*(8. - 16.*C4*s)) +
                            pow(ma1,2)*(-32.*C4*pow(mpion,4) + s*(4. - 8.*C4*s) + pow(mrho,2)*(-6. + 16.*C4*s) +
                               pow(mpion,2)*(-12. - 32.*C4*pow(mrho,2) + 32.*C4*s)))))*log(fabs(-pow(ma1,2) + t1)))/pow(mrho,2) -
                 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) + t1)) -
                 (0.5*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(mpion,6) + 1.*pow(mpion,4)*s) +
                      eta1*(1.*pow(mpion,6) + pow(mpion,4)*(-0.5*pow(mrho,2) - 1.*s) + 0.5*pow(mpion,2)*pow(mrho,2)*s))*
                    log(fabs(-pow(mpion,2) + t1)))/(pow(ma1,2) - 1.*pow(mpion,2)) +
                 (0.5*(pow(delta,2)*(1.*pow(mpion,6) + 0.5*pow(mrho,6) + pow(mpion,4)*(-1.*pow(mrho,2) - 2.*s) - 0.5*pow(mrho,4)*s +
                         pow(mpion,2)*(-1.5*pow(mrho,4) + 1.*pow(mrho,2)*s + 1.*pow(s,2))) +
                      delta*(-2.*pow(mrho,6) + pow(mpion,6)*(-2. - 4.*C4*pow(mrho,2)) + 2.*pow(mrho,4)*s +
                         pow(mpion,2)*(6.*pow(mrho,4) - 2.*pow(s,2) + pow(mrho,2)*s*(-4. - 4.*C4*s)) +
                         pow(mpion,4)*(4.*s + pow(mrho,2)*(4. + 8.*C4*s))) +
                      pow(mrho,2)*(8.*C4*pow(mpion,6) + 2.*pow(mrho,4) - 2.*pow(mrho,2)*s + pow(mpion,4)*(-4. - 16.*C4*s) +
                         pow(mpion,2)*(-6.*pow(mrho,2) + s*(4. + 8.*C4*s))))*log(fabs(-pow(mpion,2) + t1)))/(pow(mrho,2)*(1.*pow(mpion,2) - 1.*s)) +
                 (0.25*(1.*eta1 - 1.*eta2)*(eta2*((2. - 1.*delta)*pow(mpion,6) + pow(mpion,2)*s*((-12. + 6.*delta)*pow(mrho,2) + (6. - 3.*delta)*s) +
                         pow(s,2)*((4. - 2.*delta)*pow(mrho,2) + (-2. + 1.*delta)*s) + pow(mpion,4)*((8. - 4.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s))
                        + eta1*((-2. + 1.*delta)*pow(mpion,6) + (-2. + 1.*delta)*pow(mrho,4)*s + (2. - 1.*delta)*pow(s,3) +
                         pow(mpion,4)*((-4. + 2.*delta)*pow(mrho,2) + (6. - 3.*delta)*s) +
                         pow(mpion,2)*((2. - 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-6. + 3.*delta)*pow(s,2))))*
                    log(fabs(-2.*pow(mpion,2) + s + t1)))/(pow(ma1,2) - 2.*pow(mpion,2) + s) +
                 (1.*pow(2. - 1.*delta,2)*(1.*pow(mpion,4) - 0.25*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s))*
                    log(fabs(-2.*pow(mpion,2) + s + t1)))/(1.*pow(mpion,2) - 1.*s) +
                 (0.9999974999999999*(pow(mrho,2)*(-2.00000750001875*pow(mpion,2) - 1.*pow(mrho,2) + 2.00000750001875*s) +
                      pow(delta,2)*(1.*pow(mpion,4) - 0.7500018750046875*pow(mrho,4) +
                         pow(mpion,2)*(1.500003750009375*pow(mrho,2) - 2.00000750001875*s) + 0.5000015000037501*pow(mrho,2)*s + 1.*pow(s,2)) +
                      delta*(-2.00000750001875*pow(mpion,4) + 2.00000750001875*pow(mrho,4) - 2.00000750001875*pow(mrho,2)*s -
                         2.00000750001875*pow(s,2) + pow(mpion,2)*(-2.00000750001875*pow(mrho,2) + 4.00000750001875*s)))*
                    log(fabs(-2.*pow(mpion,2) + s + t1)))/pow(mrho,2) -
                 (1.*(pow(mrho,2)*(pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) + (2. - 4.*C4*pow(mrho,2))*s) +
                      pow(delta,2)*(1.*pow(mpion,4) - 0.5*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 2.*s) + 1.*pow(s,2)) +
                      delta*(-2.*pow(mpion,4) - 1.*pow(mrho,2)*s - 2.*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) - 2.*C4*pow(mrho,4) + 4.*s) +
                         pow(mrho,4)*(1. + 2.*C4*s)))*log(fabs(-2.*pow(mpion,2) + s + t1)))/pow(mrho,2)))/(16.*Pi*s*(-4*pow(mpion,2) + s)) );

          /* omega: xsection = to_mb*1/3.0*(0.0024867959858108648*pow(Const,2)*pow(g_POR,4)*(pow(mpion,8)*(t2 - t1) + pow(mpion,6)*pow(mrho,2)*(-2.*t2 + 2.*t1) +
               pow(mpion,4)*(pow(mrho,4)*(t2 - t1) + s*(4.*s*t2 - pow(t2,2) - 4.*s*t1 + pow(t1,2))) +
               pow(s,2)*(pow(s,2)*t2 + s*pow(t2,2) + 0.6666666666666666*pow(t2,3) + pow(mrho,4)*(t2 - t1) -
                  pow(s,2)*t1 - s*pow(t1,2) - 0.6666666666666666*pow(t1,3) +
                  pow(mrho,2)*(-2.*s*t2 - pow(t2,2) + 2.*s*t1 + pow(t1,2))) +
               pow(mpion,2)*s*(pow(mrho,4)*(-2.*t2 + 2.*t1) + pow(mrho,2)*(4.*s*t2 + pow(t2,2) - 4.*s*t1 - pow(t1,2)) +
                  s*(-4.*s*t2 - 2.*pow(t2,2) + 4.*s*t1 + 2.*pow(t1,2)))))/
           (pow(pow(momega,2) - s,2)*(pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-2.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s + pow(s,2)));*/

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

          mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
          t1 = mandelstam_t[1];
          t2 = mandelstam_t[0];

          xsection = to_mb*(1/3.0*(pow(Const,2)*pow(ghat,4)*((pow(2. - 1.*delta,2)*(0.25*pow(mpion,4) - 0.12500000000000003*pow(mrho,4) +
                      pow(mpion,2)*(0.5000000000000001*pow(mrho,2) - 0.5000000000000001*s) + 0.25*pow(s,2)))/(1.*pow(mpion,2) - 0.5*s - 0.5*t2) -
                 (0.25*pow(-2. + delta,2)*pow(mpion,2)*(pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-2.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s +
                      1.*pow(s,2))*t2)/(pow(mrho,2)*pow(pow(mpion,2) - 1.*s,2)) -
                 (0.5000000625000001*(pow(mrho,4)*(1.*pow(mpion,2) - 1.499996062500492*pow(mrho,2) + 0.4999994375000703*s) +
                      delta*pow(mrho,2)*(-2.0000022499997185*pow(mpion,4) - 2.0000022499997185*pow(mpion,2)*pow(mrho,2) -
                         0.4999994375000703*pow(mrho,4) + 2.0000022499997185*pow(mrho,2)*s + 0.4999994375000703*pow(s,2)) +
                      pow(delta,2)*(1.*pow(mpion,6) + 1.1249998593750175*pow(mrho,6) +
                         pow(mpion,4)*(0.4999994375000703*pow(mrho,2) - 0.4999994375000703*s) - 1.6249947968756502*pow(mrho,4)*s -
                         0.6250004218749473*pow(mrho,2)*pow(s,2) + 0.12499985937501758*pow(s,3) +
                         pow(mpion,2)*(-0.24999971875003515*pow(mrho,4) + 0.4999994375000703*pow(mrho,2)*s - 0.24999971875003515*pow(s,2))))*t2)/
                  pow(mrho,6) - (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
                    (eta2*(pow(mpion,6) + pow(mpion,2)*pow(s,2) + (pow(mrho,2) - 1.*s)*pow(s,2) + pow(mpion,4)*(-1.*pow(mrho,2) + 3.*s)) +
                      eta1*(-4.*pow(mpion,6) + pow(mpion,4)*(3.*pow(mrho,2) + s) +
                         pow(mpion,2)*(-1.*pow(mrho,4) + pow(mrho,2)*s - 2.*pow(s,2)) + s*(pow(mrho,4) - 2.*pow(mrho,2)*s + pow(s,2))))*t2)/
                  ((-1.*pow(mpion,2) + s)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
                 (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
                         pow(mpion,2)*s*(-4.*pow(mrho,4) + 8.*pow(mrho,2)*s - 4.*pow(s,2)) +
                         pow(s,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s + 1.*pow(s,2)) +
                         pow(mpion,4)*(3.*pow(mrho,4) - 6.*pow(mrho,2)*s + 4.*pow(s,2))) +
                      pow(eta2,2)*(1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
                         pow(mpion,2)*s*(-2.*pow(mrho,4) - 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
                         pow(s,2)*(1.*pow(mrho,4) + 2.*pow(mrho,2)*s + 1.*pow(s,2)) +
                         pow(mpion,4)*(1.*pow(mrho,4) + 4.*pow(mrho,2)*s + 4.*pow(s,2))) +
                      eta1*eta2*(-2.*pow(mpion,8) + 2.*pow(mrho,4)*pow(s,2) - 2.*pow(s,4) +
                         pow(mpion,4)*(2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 8.*pow(s,2)) +
                         pow(mpion,2)*s*(-4.*pow(mrho,4) - 4.*pow(mrho,2)*s + 8.*pow(s,2))))*t2)/
                  (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) +
                 (1.*(-2.*C4*pow(mrho,8) - 0.25*pow(mrho,4)*s + pow(mrho,6)*(1. + 1.*C4*s) +
                      pow(delta,2)*(1.*pow(mpion,6) + 1.*pow(mrho,6) - 0.75*pow(mpion,4)*s - 1.1875*pow(mrho,4)*s - 0.375*pow(mrho,2)*pow(s,2) +
                         0.0625*pow(s,3) + pow(mpion,2)*(-0.75*pow(mrho,4) + 0.25*pow(mrho,2)*s)) +
                      delta*pow(mrho,2)*(-1.*pow(mpion,4) + pow(mpion,2)*(0.5*s + pow(mrho,2)*(-0.5 - 1.*C4*s)) +
                         pow(mrho,2)*(1.*C4*pow(mrho,4) + pow(mrho,2)*(-1.5 + 0.5*C4*s) + s*(1. + 0.5*C4*s))))*t2)/pow(mrho,6) -
                 (0.5*(pow(mrho,6)*(-1.5 + C4*(-12.*pow(mpion,2) + 6.*s) + pow(C4,2)*(-16.*pow(mpion,4) + 16.*pow(mpion,2)*s - 4.*pow(s,2))) +
                      pow(delta,2)*(1.*pow(mpion,6) + 0.125*pow(mrho,6) + pow(mpion,4)*(-2.*pow(mrho,2) - 1.*s) + 0.25*pow(mrho,4)*s -
                         0.625*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(-2.5*pow(mrho,4) + 1.75*pow(mrho,2)*s + 0.25*pow(s,2))) +
                      delta*pow(mrho,2)*(pow(mpion,4)*(1. + 8.*C4*pow(mrho,2)) +
                         pow(mpion,2)*(6.*C4*pow(mrho,4) - 0.5*s + pow(mrho,2)*(3. - 10.*C4*s)) +
                         pow(mrho,2)*(pow(mrho,2)*(1.5 - 1.*C4*s) + s*(-2.5 + 3.*C4*s))))*t2)/pow(mrho,6) -
                 (0.25*(pow(delta,2)*(1.*pow(mpion,6) - 1.*pow(mrho,6) + pow(mpion,4)*(-2.5*pow(mrho,2) - 2.5*s) - 1.5*pow(mrho,4)*s +
                         2.*pow(mrho,2)*pow(s,2) - 0.5*pow(s,3) + pow(mpion,2)*(3.5*pow(mrho,4) - 1.5*pow(mrho,2)*s + 2.*pow(s,2))) +
                      pow(mrho,2)*(pow(mpion,4)*(-6. - 8.*C4*pow(mrho,2)) + 2.*pow(s,2) + pow(mrho,4)*(-4. - 8.*C4*s) +
                         pow(mrho,2)*s*(-2. + 8.*C4*s) + pow(mpion,2)*(8.*C4*pow(mrho,4) + 4.*s + pow(mrho,2)*(10. - 16.*C4*s))) +
                      delta*(-2.*pow(mpion,6) - 5.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
                         pow(mpion,4)*(8.*pow(mrho,2) + 4.*C4*pow(mrho,4) + 5.*s) + pow(mrho,4)*s*(4. - 4.*C4*s) + pow(mrho,6)*(4. + 4.*C4*s) +
                         pow(mpion,2)*(-4.*C4*pow(mrho,6) + 1.*pow(mrho,2)*s - 4.*pow(s,2) + pow(mrho,4)*(-12. + 8.*C4*s))))*t2)/
                  (pow(mrho,4)*(pow(mpion,2) - 1.*s)) + (0.0625*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
                    (pow(mrho,2)*(eta2*(4.*pow(mpion,4) - 6.*pow(mpion,2)*s + s*(8.*pow(mrho,2) + 6.*s)) +
                         eta1*(-12.*pow(mpion,4) + 4.*pow(mrho,4) + 2.*pow(mrho,2)*s - 6.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s))) +
                      delta*(eta1*(8.*pow(mpion,6) - 2.*pow(mrho,6) + pow(mpion,4)*(2.*pow(mrho,2) - 2.*s) - 3.*pow(mrho,4)*s +
                            4.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(2.*pow(mrho,4) - 2.*pow(s,2))) +
                         eta2*(pow(mpion,4)*(-2.*pow(mrho,2) - 4.*s) + pow(mpion,2)*s*(3.*pow(mrho,2) + 3.*s) +
                            s*(-4.*pow(mrho,4) - 7.*pow(mrho,2)*s - 1.*pow(s,2)))))*t2)/
                  (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
                 (0.1875*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(delta*
                       (eta1*(2.6666666666666665*pow(mpion,6) + pow(mpion,4)*(-4.*pow(mrho,2) + 2.*s) +
                            pow(mpion,2)*(-1.3333333333333333*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.3333333333333335*pow(s,2)) +
                            s*(0.3333333333333333*pow(mrho,4) - 1.3333333333333333*pow(mrho,2)*s + 1.*pow(s,2))) +
                         eta2*(pow(mpion,4)*(-0.6666666666666666*pow(mrho,2) - 4.*s) +
                            s*(0.6666666666666666*pow(mrho,4) - 1.*pow(mrho,2)*s - 1.*pow(s,2)) +
                            pow(mpion,2)*(-0.6666666666666666*pow(mrho,4) - 0.3333333333333333*pow(mrho,2)*s + 3.6666666666666665*pow(s,2)))) +
                      pow(mrho,2)*(eta2*(C4*pow(mpion,4)*(2.6666666666666665*pow(mrho,2) + 10.666666666666666*s) +
                            pow(mpion,2)*(s*(3.3333333333333335 - 10.666666666666666*C4*s) + pow(mrho,2)*(1.3333333333333333 - 5.333333333333333*C4*s)) +
                            s*(s*(-2. + 2.6666666666666665*C4*s) + pow(mrho,2)*(-1.3333333333333333 + 2.6666666666666665*C4*s))) +
                         eta1*(pow(mpion,4)*(1.3333333333333333 + 8.*C4*pow(mrho,2) - 10.666666666666666*C4*s) +
                            s*(s*(2. - 2.6666666666666665*C4*s) + pow(mrho,2)*(-2. + 2.6666666666666665*C4*s)) +
                            pow(mpion,2)*(pow(mrho,2)*(2.6666666666666665 - 10.666666666666666*C4*s) + s*(-4. + 10.666666666666666*C4*s)))))*t2)/
                  (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
                 (0.0625*(-2. + delta)*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(pow(mpion,2) + s)*(-2.*eta2*s + eta1*(pow(mpion,2) - 1.*pow(mrho,2) + s))*
                    pow(t2,2))/((-1.*pow(mpion,2) + s)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
                 (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(pow(mrho,2) - 1.*s) + 2.*eta1*eta2*s - 1.*pow(eta2,2)*s)*
                    (pow(mpion,4) + (pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 2.*s))*pow(t2,2))/
                  (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) +
                 ((-0.125*pow(mrho,4) + delta*(0.5*pow(mrho,4) - 0.25*pow(mrho,2)*s) +
                      pow(delta,2)*(0.125*pow(mpion,4) - 0.09375*pow(mrho,4) + pow(mpion,2)*(-0.125*pow(mrho,2) + 0.125*s) +
                         0.0625*pow(mrho,2)*s - 0.09375*pow(s,2)))*pow(t2,2))/pow(mrho,6) -
                 (0.125*(-1.*pow(mrho,4) + 4.*C4*pow(mrho,6) + delta*pow(mrho,2)*
                       (2.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) - 2.*s) +
                      pow(delta,2)*(1.*pow(mpion,4) + 0.25*pow(mrho,4) - 1.25*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 2.*s)))*pow(t2,2))/
                  pow(mrho,6) + (0.0625*(-32.*C4*pow(mrho,4)*s + pow(delta,2)*
                       (1.*pow(mpion,4) + pow(mpion,2)*(-1.*pow(mrho,2) - 2.*s) + s*(-3.*pow(mrho,2) + 1.*s)) +
                      delta*(-2.*pow(mpion,4) + (6.*pow(mrho,2) + 16.*C4*pow(mrho,4) - 2.*s)*s + pow(mpion,2)*(2.*pow(mrho,2) + 4.*s)))*pow(t2,2)
                    )/(pow(mrho,4)*(pow(mpion,2) - 1.*s)) - (0.5625*
                    (C4*pow(mrho,6)*(2.6666666666666665 + 7.111111111111112*C4*pow(mpion,2) - 3.555555555555556*C4*s) +
                      pow(delta,2)*(0.11111111111111112*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 0.22222222222222224*s) -
                         0.22222222222222224*pow(mrho,2)*s + 0.11111111111111112*pow(s,2)) +
                      delta*pow(mrho,2)*(-2.2222222222222223*C4*pow(mrho,4) + pow(mpion,2)*(-0.6666666666666666 - 2.6666666666666665*C4*pow(mrho,2)) +
                         0.22222222222222224*s + pow(mrho,2)*(-0.22222222222222224 + 1.777777777777778*C4*s)))*pow(t2,2))/pow(mrho,6) +
                 (0.03125*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(pow(mrho,2)*
                       (-2.*eta2*pow(mpion,2) - 5.999999999999999*eta1*pow(mrho,2) + 8.*eta1*s - 2.*eta2*s) +
                      delta*(eta1*(-5.999999999999999*pow(mpion,4) + 5.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) -
                            5.999999999999999*pow(mrho,2)*s + 1.*pow(s,2)) +
                         eta2*(4.*pow(mpion,4) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*s) + s*(5.*pow(mrho,2) + 2.*s))))*pow(t2,2))/
                  (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
                 (0.15625*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(delta*
                       (eta1*(-1.2*pow(mpion,4) + 0.6*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 2.4*s) - 1.6*pow(mrho,2)*s + 1.*pow(s,2)) +
                         eta2*(0.8*pow(mpion,4) + (1.*pow(mrho,2) - 0.4*s)*s + pow(mpion,2)*(0.2*pow(mrho,2) + 1.2*s))) +
                      pow(mrho,2)*(eta2*(pow(mpion,2)*(-0.4 - 6.4*C4*s) + s*(-0.4 + 3.2*C4*s)) +
                         eta1*(s*(0.8 - 3.2*C4*s) + pow(mrho,2)*(-0.4 + 3.2*C4*s) + pow(mpion,2)*(-0.8 - 3.2*C4*pow(mrho,2) + 6.4*C4*s))))*pow(t2,2))
                   /(pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) + 0.666667*pow(C4,2)*pow(t2,3) +
                 (0.125*pow(delta,2)*pow(mpion,2)*pow(t2,3))/pow(mrho,6) - (0.0833333*delta*pow(t2,3))/pow(mrho,4) +
                 (0.145833*pow(delta,2)*pow(t2,3))/pow(mrho,4) - (0.166667*C4*delta*pow(t2,3))/pow(mrho,2) -
                 (0.208333*delta*(-0.8*pow(mrho,2) + 0.8*C4*pow(mrho,4) + delta*(0.8*pow(mpion,2) + 1.*pow(mrho,2) - 0.7*s))*pow(t2,3))/
                  pow(mrho,6) + (0.104167*delta*(-0.8*pow(mrho,2) + delta*(0.4*pow(mpion,2) + 1.*pow(mrho,2) - 0.6*s))*pow(t2,3))/pow(mrho,6) -
                 (0.0833333*pow(delta,2)*s*pow(t2,3))/pow(mrho,6) +
                 (0.0208333*pow(eta1 - 1.*eta2,2)*s*(-2.*eta1*eta2*s + pow(eta2,2)*s + pow(eta1,2)*(-1.*pow(mrho,2) + s))*pow(t2,3))/
                  (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) +
                 (0.104166625*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(0.4*eta1*pow(mrho,2) +
                      delta*(-0.2*eta2*pow(mpion,2) - 0.2*eta2*s + eta1*(-0.4*pow(mpion,2) - 0.8*pow(mrho,2) + 1.*s)))*pow(t2,3))/
                  (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
                 (0.14583362500000002*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
                    (delta*(-0.14285714285714285*eta2*pow(mpion,2) - 0.42857142857142855*eta2*s +
                         eta1*(-0.2857142857142857*pow(mpion,2) - 0.5714285714285714*pow(mrho,2) + 1.*s)) +
                      pow(mrho,2)*(1.1428571428571428*C4*eta2*s + eta1*(0.2857142857142857 - 1.1428571428571428*C4*s)))*pow(t2,3))/
                  (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
                 (pow(2. - 1.*delta,2)*(-0.25*pow(mpion,4) + 0.12500000000000003*pow(mrho,4) +
                      pow(mpion,2)*(-0.5000000000000001*pow(mrho,2) + 0.5000000000000001*s) - 0.25*pow(s,2)))/(1.*pow(mpion,2) - 0.5*s - 0.5*t1) +
                 (0.25*pow(-2. + delta,2)*pow(mpion,2)*(pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-2.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s +
                      1.*pow(s,2))*t1)/(pow(mrho,2)*pow(pow(mpion,2) - 1.*s,2)) +
                 (0.5000000625000001*(pow(mrho,4)*(1.*pow(mpion,2) - 1.499996062500492*pow(mrho,2) + 0.4999994375000703*s) +
                      delta*pow(mrho,2)*(-2.0000022499997185*pow(mpion,4) - 2.0000022499997185*pow(mpion,2)*pow(mrho,2) -
                         0.4999994375000703*pow(mrho,4) + 2.0000022499997185*pow(mrho,2)*s + 0.4999994375000703*pow(s,2)) +
                      pow(delta,2)*(1.*pow(mpion,6) + 1.1249998593750175*pow(mrho,6) +
                         pow(mpion,4)*(0.4999994375000703*pow(mrho,2) - 0.4999994375000703*s) - 1.6249947968756502*pow(mrho,4)*s -
                         0.6250004218749473*pow(mrho,2)*pow(s,2) + 0.12499985937501758*pow(s,3) +
                         pow(mpion,2)*(-0.24999971875003515*pow(mrho,4) + 0.4999994375000703*pow(mrho,2)*s - 0.24999971875003515*pow(s,2))))*t1)/
                  pow(mrho,6) + (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
                    (eta2*(pow(mpion,6) + pow(mpion,2)*pow(s,2) + (pow(mrho,2) - 1.*s)*pow(s,2) + pow(mpion,4)*(-1.*pow(mrho,2) + 3.*s)) +
                      eta1*(-4.*pow(mpion,6) + pow(mpion,4)*(3.*pow(mrho,2) + s) +
                         pow(mpion,2)*(-1.*pow(mrho,4) + pow(mrho,2)*s - 2.*pow(s,2)) + s*(pow(mrho,4) - 2.*pow(mrho,2)*s + pow(s,2))))*t1)/
                  ((-1.*pow(mpion,2) + s)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
                 (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
                         pow(mpion,2)*s*(-4.*pow(mrho,4) + 8.*pow(mrho,2)*s - 4.*pow(s,2)) +
                         pow(s,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s + 1.*pow(s,2)) +
                         pow(mpion,4)*(3.*pow(mrho,4) - 6.*pow(mrho,2)*s + 4.*pow(s,2))) +
                      pow(eta2,2)*(1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
                         pow(mpion,2)*s*(-2.*pow(mrho,4) - 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
                         pow(s,2)*(1.*pow(mrho,4) + 2.*pow(mrho,2)*s + 1.*pow(s,2)) +
                         pow(mpion,4)*(1.*pow(mrho,4) + 4.*pow(mrho,2)*s + 4.*pow(s,2))) +
                      eta1*eta2*(-2.*pow(mpion,8) + 2.*pow(mrho,4)*pow(s,2) - 2.*pow(s,4) +
                         pow(mpion,4)*(2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 8.*pow(s,2)) +
                         pow(mpion,2)*s*(-4.*pow(mrho,4) - 4.*pow(mrho,2)*s + 8.*pow(s,2))))*t1)/
                  (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) -
                 (1.*(-2.*C4*pow(mrho,8) - 0.25*pow(mrho,4)*s + pow(mrho,6)*(1. + 1.*C4*s) +
                      pow(delta,2)*(1.*pow(mpion,6) + 1.*pow(mrho,6) - 0.75*pow(mpion,4)*s - 1.1875*pow(mrho,4)*s - 0.375*pow(mrho,2)*pow(s,2) +
                         0.0625*pow(s,3) + pow(mpion,2)*(-0.75*pow(mrho,4) + 0.25*pow(mrho,2)*s)) +
                      delta*pow(mrho,2)*(-1.*pow(mpion,4) + pow(mpion,2)*(0.5*s + pow(mrho,2)*(-0.5 - 1.*C4*s)) +
                         pow(mrho,2)*(1.*C4*pow(mrho,4) + pow(mrho,2)*(-1.5 + 0.5*C4*s) + s*(1. + 0.5*C4*s))))*t1)/pow(mrho,6) +
                 (0.5*(pow(mrho,6)*(-1.5 + C4*(-12.*pow(mpion,2) + 6.*s) + pow(C4,2)*(-16.*pow(mpion,4) + 16.*pow(mpion,2)*s - 4.*pow(s,2))) +
                      pow(delta,2)*(1.*pow(mpion,6) + 0.125*pow(mrho,6) + pow(mpion,4)*(-2.*pow(mrho,2) - 1.*s) + 0.25*pow(mrho,4)*s -
                         0.625*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(-2.5*pow(mrho,4) + 1.75*pow(mrho,2)*s + 0.25*pow(s,2))) +
                      delta*pow(mrho,2)*(pow(mpion,4)*(1. + 8.*C4*pow(mrho,2)) +
                         pow(mpion,2)*(6.*C4*pow(mrho,4) - 0.5*s + pow(mrho,2)*(3. - 10.*C4*s)) +
                         pow(mrho,2)*(pow(mrho,2)*(1.5 - 1.*C4*s) + s*(-2.5 + 3.*C4*s))))*t1)/pow(mrho,6) +
                 (0.25*(pow(delta,2)*(1.*pow(mpion,6) - 1.*pow(mrho,6) + pow(mpion,4)*(-2.5*pow(mrho,2) - 2.5*s) - 1.5*pow(mrho,4)*s +
                         2.*pow(mrho,2)*pow(s,2) - 0.5*pow(s,3) + pow(mpion,2)*(3.5*pow(mrho,4) - 1.5*pow(mrho,2)*s + 2.*pow(s,2))) +
                      pow(mrho,2)*(pow(mpion,4)*(-6. - 8.*C4*pow(mrho,2)) + 2.*pow(s,2) + pow(mrho,4)*(-4. - 8.*C4*s) +
                         pow(mrho,2)*s*(-2. + 8.*C4*s) + pow(mpion,2)*(8.*C4*pow(mrho,4) + 4.*s + pow(mrho,2)*(10. - 16.*C4*s))) +
                      delta*(-2.*pow(mpion,6) - 5.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
                         pow(mpion,4)*(8.*pow(mrho,2) + 4.*C4*pow(mrho,4) + 5.*s) + pow(mrho,4)*s*(4. - 4.*C4*s) + pow(mrho,6)*(4. + 4.*C4*s) +
                         pow(mpion,2)*(-4.*C4*pow(mrho,6) + 1.*pow(mrho,2)*s - 4.*pow(s,2) + pow(mrho,4)*(-12. + 8.*C4*s))))*t1)/
                  (pow(mrho,4)*(pow(mpion,2) - 1.*s)) - (0.0625*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
                    (pow(mrho,2)*(eta2*(4.*pow(mpion,4) - 6.*pow(mpion,2)*s + s*(8.*pow(mrho,2) + 6.*s)) +
                         eta1*(-12.*pow(mpion,4) + 4.*pow(mrho,4) + 2.*pow(mrho,2)*s - 6.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s))) +
                      delta*(eta1*(8.*pow(mpion,6) - 2.*pow(mrho,6) + pow(mpion,4)*(2.*pow(mrho,2) - 2.*s) - 3.*pow(mrho,4)*s +
                            4.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(2.*pow(mrho,4) - 2.*pow(s,2))) +
                         eta2*(pow(mpion,4)*(-2.*pow(mrho,2) - 4.*s) + pow(mpion,2)*s*(3.*pow(mrho,2) + 3.*s) +
                            s*(-4.*pow(mrho,4) - 7.*pow(mrho,2)*s - 1.*pow(s,2)))))*t1)/
                  (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
                 (0.1875*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(delta*
                       (eta1*(2.6666666666666665*pow(mpion,6) + pow(mpion,4)*(-4.*pow(mrho,2) + 2.*s) +
                            pow(mpion,2)*(-1.3333333333333333*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.3333333333333335*pow(s,2)) +
                            s*(0.3333333333333333*pow(mrho,4) - 1.3333333333333333*pow(mrho,2)*s + 1.*pow(s,2))) +
                         eta2*(pow(mpion,4)*(-0.6666666666666666*pow(mrho,2) - 4.*s) +
                            s*(0.6666666666666666*pow(mrho,4) - 1.*pow(mrho,2)*s - 1.*pow(s,2)) +
                            pow(mpion,2)*(-0.6666666666666666*pow(mrho,4) - 0.3333333333333333*pow(mrho,2)*s + 3.6666666666666665*pow(s,2)))) +
                      pow(mrho,2)*(eta2*(C4*pow(mpion,4)*(2.6666666666666665*pow(mrho,2) + 10.666666666666666*s) +
                            pow(mpion,2)*(s*(3.3333333333333335 - 10.666666666666666*C4*s) + pow(mrho,2)*(1.3333333333333333 - 5.333333333333333*C4*s)) +
                            s*(s*(-2. + 2.6666666666666665*C4*s) + pow(mrho,2)*(-1.3333333333333333 + 2.6666666666666665*C4*s))) +
                         eta1*(pow(mpion,4)*(1.3333333333333333 + 8.*C4*pow(mrho,2) - 10.666666666666666*C4*s) +
                            s*(s*(2. - 2.6666666666666665*C4*s) + pow(mrho,2)*(-2. + 2.6666666666666665*C4*s)) +
                            pow(mpion,2)*(pow(mrho,2)*(2.6666666666666665 - 10.666666666666666*C4*s) + s*(-4. + 10.666666666666666*C4*s)))))*t1)/
                  (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
                 (0.0625*(-2. + delta)*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(pow(mpion,2) + s)*(-2.*eta2*s + eta1*(pow(mpion,2) - 1.*pow(mrho,2) + s))*
                    pow(t1,2))/((-1.*pow(mpion,2) + s)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
                 (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(pow(mrho,2) - 1.*s) + 2.*eta1*eta2*s - 1.*pow(eta2,2)*s)*
                    (pow(mpion,4) + (pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 2.*s))*pow(t1,2))/
                  (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) +
                 ((0.125*pow(mrho,4) + delta*(-0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s) +
                      pow(delta,2)*(-0.125*pow(mpion,4) + 0.09375*pow(mrho,4) + pow(mpion,2)*(0.125*pow(mrho,2) - 0.125*s) -
                         0.0625*pow(mrho,2)*s + 0.09375*pow(s,2)))*pow(t1,2))/pow(mrho,6) +
                 (0.125*(-1.*pow(mrho,4) + 4.*C4*pow(mrho,6) + delta*pow(mrho,2)*
                       (2.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) - 2.*s) +
                      pow(delta,2)*(1.*pow(mpion,4) + 0.25*pow(mrho,4) - 1.25*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 2.*s)))*pow(t1,2))/
                  pow(mrho,6) - (0.0625*(-32.*C4*pow(mrho,4)*s + pow(delta,2)*
                       (1.*pow(mpion,4) + pow(mpion,2)*(-1.*pow(mrho,2) - 2.*s) + s*(-3.*pow(mrho,2) + 1.*s)) +
                      delta*(-2.*pow(mpion,4) + (6.*pow(mrho,2) + 16.*C4*pow(mrho,4) - 2.*s)*s + pow(mpion,2)*(2.*pow(mrho,2) + 4.*s)))*pow(t1,2)
                    )/(pow(mrho,4)*(pow(mpion,2) - 1.*s)) + (0.5625*
                    (C4*pow(mrho,6)*(2.6666666666666665 + 7.111111111111112*C4*pow(mpion,2) - 3.555555555555556*C4*s) +
                      pow(delta,2)*(0.11111111111111112*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 0.22222222222222224*s) -
                         0.22222222222222224*pow(mrho,2)*s + 0.11111111111111112*pow(s,2)) +
                      delta*pow(mrho,2)*(-2.2222222222222223*C4*pow(mrho,4) + pow(mpion,2)*(-0.6666666666666666 - 2.6666666666666665*C4*pow(mrho,2)) +
                         0.22222222222222224*s + pow(mrho,2)*(-0.22222222222222224 + 1.777777777777778*C4*s)))*pow(t1,2))/pow(mrho,6) -
                 (0.03125*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(pow(mrho,2)*
                       (-2.*eta2*pow(mpion,2) - 5.999999999999999*eta1*pow(mrho,2) + 8.*eta1*s - 2.*eta2*s) +
                      delta*(eta1*(-5.999999999999999*pow(mpion,4) + 5.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) -
                            5.999999999999999*pow(mrho,2)*s + 1.*pow(s,2)) +
                         eta2*(4.*pow(mpion,4) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*s) + s*(5.*pow(mrho,2) + 2.*s))))*pow(t1,2))/
                  (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
                 (0.15625*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(delta*
                       (eta1*(-1.2*pow(mpion,4) + 0.6*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 2.4*s) - 1.6*pow(mrho,2)*s + 1.*pow(s,2)) +
                         eta2*(0.8*pow(mpion,4) + (1.*pow(mrho,2) - 0.4*s)*s + pow(mpion,2)*(0.2*pow(mrho,2) + 1.2*s))) +
                      pow(mrho,2)*(eta2*(pow(mpion,2)*(-0.4 - 6.4*C4*s) + s*(-0.4 + 3.2*C4*s)) +
                         eta1*(s*(0.8 - 3.2*C4*s) + pow(mrho,2)*(-0.4 + 3.2*C4*s) + pow(mpion,2)*(-0.8 - 3.2*C4*pow(mrho,2) + 6.4*C4*s))))*pow(t1,2))
                   /(pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) - 0.666667*pow(C4,2)*pow(t1,3) -
                 (0.125*pow(delta,2)*pow(mpion,2)*pow(t1,3))/pow(mrho,6) + (0.0833333*delta*pow(t1,3))/pow(mrho,4) -
                 (0.145833*pow(delta,2)*pow(t1,3))/pow(mrho,4) + (0.166667*C4*delta*pow(t1,3))/pow(mrho,2) +
                 (0.208333*delta*(-0.8*pow(mrho,2) + 0.8*C4*pow(mrho,4) + delta*(0.8*pow(mpion,2) + 1.*pow(mrho,2) - 0.7*s))*pow(t1,3))/
                  pow(mrho,6) - (0.104167*delta*(-0.8*pow(mrho,2) + delta*(0.4*pow(mpion,2) + 1.*pow(mrho,2) - 0.6*s))*pow(t1,3))/pow(mrho,6) +
                 (0.0833333*pow(delta,2)*s*pow(t1,3))/pow(mrho,6) -
                 (0.0208333*pow(eta1 - 1.*eta2,2)*s*(-2.*eta1*eta2*s + pow(eta2,2)*s + pow(eta1,2)*(-1.*pow(mrho,2) + s))*pow(t1,3))/
                  (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) -
                 (0.104166625*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(0.4*eta1*pow(mrho,2) +
                      delta*(-0.2*eta2*pow(mpion,2) - 0.2*eta2*s + eta1*(-0.4*pow(mpion,2) - 0.8*pow(mrho,2) + 1.*s)))*pow(t1,3))/
                  (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
                 (0.14583362500000002*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
                    (delta*(-0.14285714285714285*eta2*pow(mpion,2) - 0.42857142857142855*eta2*s +
                         eta1*(-0.2857142857142857*pow(mpion,2) - 0.5714285714285714*pow(mrho,2) + 1.*s)) +
                      pow(mrho,2)*(1.1428571428571428*C4*eta2*s + eta1*(0.2857142857142857 - 1.1428571428571428*C4*s)))*pow(t1,3))/
                  (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
                 (0.25*(pow(mrho,2)*(16.*pow(mpion,4) + 4.*pow(mrho,4) + pow(mpion,2)*(-16.*pow(mrho,2) - 16.*s)) +
                      delta*pow(mrho,2)*(-16.*pow(mpion,4) - 4.*pow(mrho,4) + pow(mpion,2)*(16.*pow(mrho,2) + 16.*s)) +
                      pow(delta,2)*(1.77636e-15*pow(mpion,6) + 1.*pow(mrho,6) + pow(mpion,4)*(4.*pow(mrho,2) - 1.77636e-15*s) +
                         pow(mpion,2)*(-4.*pow(mrho,4) - 4.*pow(mrho,2)*s)))*log(fabs(-pow(mpion,2) + 0.5*s + 0.5*t2)))/
                  (pow(mrho,2)*(pow(mpion,2) - 1.*s)) - (0.25*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
                    (eta2*((-2. + 1.*delta)*pow(mpion,6) + (6. - 3.*delta)*pow(mpion,4)*s +
                         pow(s,2)*((4. - 2.*delta)*pow(mrho,2) + (2. - 1.*delta)*s) +
                         pow(mpion,2)*s*((-4. + 2.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s)) +
                      eta1*((2. - 1.*delta)*pow(mpion,6) + (2. - 1.*delta)*pow(mrho,4)*s + (-2. + 1.*delta)*pow(s,3) +
                         pow(mpion,4)*((4. - 2.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s) +
                         pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (-4. + 2.*delta)*pow(mrho,2)*s + (6. - 3.*delta)*pow(s,2))))*
                    log(fabs(-2.*pow(mpion,2) + s + t2)))/(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) -
                 (0.9999974999999999*(pow(mrho,2)*(-2.00000750001875*pow(mpion,2) - 1.*pow(mrho,2) + 2.00000750001875*s) +
                      pow(delta,2)*(1.*pow(mpion,4) - 0.7500018750046875*pow(mrho,4) +
                         pow(mpion,2)*(1.500003750009375*pow(mrho,2) - 2.00000750001875*s) + 0.5000015000037501*pow(mrho,2)*s + 1.*pow(s,2)) +
                      delta*(-2.00000750001875*pow(mpion,4) + 2.00000750001875*pow(mrho,4) - 2.00000750001875*pow(mrho,2)*s -
                         2.00000750001875*pow(s,2) + pow(mpion,2)*(-2.00000750001875*pow(mrho,2) + 4.00000750001875*s)))*
                    log(fabs(-2.*pow(mpion,2) + s + t2)))/pow(mrho,2) +
                 (1.*(pow(mrho,2)*(pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) + (-2. + 4.*C4*pow(mrho,2))*s) +
                      pow(delta,2)*(1.*pow(mpion,4) - 0.5*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 2.*s) + 1.*pow(s,2)) +
                      delta*(-2.*pow(mpion,4) + 1.*pow(mrho,2)*s - 2.*pow(s,2) + pow(mpion,2)*(-5.*pow(mrho,2) + 2.*C4*pow(mrho,4) + 4.*s) +
                         pow(mrho,4)*(1. - 2.*C4*s)))*log(fabs(-2.*pow(mpion,2) + s + t2)))/pow(mrho,2) +
                 (0.25*(pow(mrho,2)*(16.*pow(mpion,4) + 4.*pow(mrho,4) + pow(mpion,2)*(-16.*pow(mrho,2) - 16.*s)) +
                      delta*pow(mrho,2)*(-16.*pow(mpion,4) - 4.*pow(mrho,4) + pow(mpion,2)*(16.*pow(mrho,2) + 16.*s)) +
                      pow(delta,2)*(1.77636e-15*pow(mpion,6) + 1.*pow(mrho,6) + pow(mpion,4)*(4.*pow(mrho,2) - 1.77636e-15*s) +
                         pow(mpion,2)*(-4.*pow(mrho,4) - 4.*pow(mrho,2)*s)))*log(fabs(-pow(mpion,2) + 0.5*s + 0.5*t1)))/
                  (pow(mrho,2)*(pow(mpion,2) - 1.*s)) + (0.25*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
                    (eta2*((-2. + 1.*delta)*pow(mpion,6) + (6. - 3.*delta)*pow(mpion,4)*s +
                         pow(s,2)*((4. - 2.*delta)*pow(mrho,2) + (2. - 1.*delta)*s) +
                         pow(mpion,2)*s*((-4. + 2.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s)) +
                      eta1*((2. - 1.*delta)*pow(mpion,6) + (2. - 1.*delta)*pow(mrho,4)*s + (-2. + 1.*delta)*pow(s,3) +
                         pow(mpion,4)*((4. - 2.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s) +
                         pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (-4. + 2.*delta)*pow(mrho,2)*s + (6. - 3.*delta)*pow(s,2))))*
                    log(fabs(-2.*pow(mpion,2) + s + t1)))/(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) +
                 (0.9999974999999999*(pow(mrho,2)*(-2.00000750001875*pow(mpion,2) - 1.*pow(mrho,2) + 2.00000750001875*s) +
                      pow(delta,2)*(1.*pow(mpion,4) - 0.7500018750046875*pow(mrho,4) +
                         pow(mpion,2)*(1.500003750009375*pow(mrho,2) - 2.00000750001875*s) + 0.5000015000037501*pow(mrho,2)*s + 1.*pow(s,2)) +
                      delta*(-2.00000750001875*pow(mpion,4) + 2.00000750001875*pow(mrho,4) - 2.00000750001875*pow(mrho,2)*s -
                         2.00000750001875*pow(s,2) + pow(mpion,2)*(-2.00000750001875*pow(mrho,2) + 4.00000750001875*s)))*
                    log(fabs(-2.*pow(mpion,2) + s + t1)))/pow(mrho,2) -
                 (1.*(pow(mrho,2)*(pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) + (-2. + 4.*C4*pow(mrho,2))*s) +
                      pow(delta,2)*(1.*pow(mpion,4) - 0.5*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 2.*s) + 1.*pow(s,2)) +
                      delta*(-2.*pow(mpion,4) + 1.*pow(mrho,2)*s - 2.*pow(s,2) + pow(mpion,2)*(-5.*pow(mrho,2) + 2.*C4*pow(mrho,4) + 4.*s) +
                         pow(mrho,4)*(1. - 2.*C4*s)))*log(fabs(-2.*pow(mpion,2) + s + t1)))/pow(mrho,2)))/(16.*Pi*s*(-4*pow(mpion,2) + s)));

          /* omega: xsection = to_mb*1/3.0*(0.0024868*pow(Const,2)*pow(g_POR,4)*((pow(momega,8) + pow(mpion,4)*pow(pow(mpion,2) - pow(mrho,2),2) -
              2*pow(momega,6)*(2*pow(mpion,2) + pow(mrho,2) - s) -
              2*pow(momega,2)*pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
              pow(momega,4)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))/
            (pow(momega,2) - t2) + 3*pow(momega,4)*t2 - 8*pow(momega,2)*pow(mpion,2)*t2 + 4*pow(mpion,4)*t2 -
           4*pow(momega,2)*pow(mrho,2)*t2 + 4*pow(mpion,2)*pow(mrho,2)*t2 + pow(mrho,4)*t2 + 4*pow(momega,2)*s*t2 -
           4*pow(mpion,2)*s*t2 - 2*pow(mrho,2)*s*t2 + 2*pow(s,2)*t2 + pow(momega,2)*pow(t2,2) - 2*pow(mpion,2)*pow(t2,2) -
           pow(mrho,2)*pow(t2,2) + s*pow(t2,2) + pow(t2,3)/3. -
           (pow(momega,8) + pow(mpion,4)*pow(pow(mpion,2) - pow(mrho,2),2) - 2*pow(momega,6)*(2*pow(mpion,2) + pow(mrho,2) - s) -
              2*pow(momega,2)*pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
              pow(momega,4)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))/
            (pow(momega,2) - t1) - 3*pow(momega,4)*t1 + 8*pow(momega,2)*pow(mpion,2)*t1 - 4*pow(mpion,4)*t1 +
           4*pow(momega,2)*pow(mrho,2)*t1 - 4*pow(mpion,2)*pow(mrho,2)*t1 - pow(mrho,4)*t1 - 4*pow(momega,2)*s*t1 +
           4*pow(mpion,2)*s*t1 + 2*pow(mrho,2)*s*t1 - 2*pow(s,2)*t1 - pow(momega,2)*pow(t1,2) + 2*pow(mpion,2)*pow(t1,2) +
           pow(mrho,2)*pow(t1,2) - s*pow(t1,2) - pow(t1,3)/3. +
           2*(2*pow(momega,6) - 3*pow(momega,4)*(2*pow(mpion,2) + pow(mrho,2) - s) -
              pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
              pow(momega,2)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))*
            log(fabs(-pow(momega,2) + t2)) - 2*(2*pow(momega,6) - 3*pow(momega,4)*(2*pow(mpion,2) + pow(mrho,2) - s) -
              pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
              pow(momega,2)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))*
            log(fabs(-pow(momega,2) + t1))))/(pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s)); */
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break;

        /*case ReactionType::pi_eta:
          if (part_a.type().pdgcode() == pdg::pi_p) {
            part_out = pi_plus_particle;
          } else {
            part_out = pi_minus_particle;
          }
          m3 = part_out->mass();

          mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
          t1 = mandelstam_t[1];
          t2 = mandelstam_t[0];

          xsection = to_be_determined * to_mb;
          process_list.push_back(make_unique<CollisionBranch>(
              *part_out, *photon_out, xsection, ProcessType::TwoToTwo));
          break; */

        case ReactionType::pi0_rho0:
          part_out = pi0_particle;
          m3 = part_out->mass();

          mandelstam_t = get_t_range(sqrts, m1, m2, m3, 0.0);
          t1 = mandelstam_t[1];
          t2 = mandelstam_t[0];

          xsection = to_mb*1/3.0*(pow(Const,2)*pow(g_POR,4)*((pow(pow(m_omega,2) - s,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(m_rho,2) + pow(m_pi,4)*(pow(m_rho,4) + 4*pow(m_omega,4) - 2*pow(m_omega,2)*s) +
                  pow(m_omega,4)*(pow(m_rho,4) + pow(m_omega,4) + 2*pow(m_omega,2)*s + 2*pow(s,2) - 2*pow(m_rho,2)*(pow(m_omega,2) + s)) -
                  2*pow(m_pi,2)*pow(m_omega,2)*(pow(m_rho,4) + 2*pow(m_omega,2)*(pow(m_omega,2) + s) - pow(m_rho,2)*(2*pow(m_omega,2) + s))))/(pow(m_omega,2) - t2) +
             (pow(m_pi,8) - 2*pow(m_pi,6)*pow(m_rho,2) + 3*pow(m_omega,8) - 4*pow(m_omega,6)*s - 7*pow(m_omega,4)*pow(s,2) + 4*pow(m_omega,2)*pow(s,3) + 5*pow(s,4) +
                pow(m_rho,4)*(pow(m_omega,4) - 2*pow(m_omega,2)*s + 2*pow(s,2)) + pow(m_rho,2)*(-4*pow(m_omega,6) + 8*pow(m_omega,4)*s - 6*pow(s,3)) -
                2*pow(m_pi,2)*(4*pow(m_omega,6) - 2*pow(m_rho,2)*pow(pow(m_omega,2) - 2*s,2) + pow(m_rho,4)*s - 10*pow(m_omega,4)*s + 8*pow(s,3)) +
                pow(m_pi,4)*(pow(m_rho,4) + 2*pow(m_rho,2)*(pow(m_omega,2) - s) + 4*(pow(m_omega,4) - 3*pow(m_omega,2)*s + 3*pow(s,2))))*t2 -
             2*pow(m_pi,2)*pow(m_omega,4)*pow(t2,2) - pow(m_rho,2)*pow(m_omega,4)*pow(t2,2) + pow(m_omega,6)*pow(t2,2) - pow(m_pi,4)*s*pow(t2,2) +
             pow(m_pi,2)*pow(m_rho,2)*s*pow(t2,2) + 8*pow(m_pi,2)*pow(m_omega,2)*s*pow(t2,2) + 3*pow(m_rho,2)*pow(m_omega,2)*s*pow(t2,2) -
             2*pow(m_omega,4)*s*pow(t2,2) - 8*pow(m_pi,2)*pow(s,2)*pow(t2,2) - 3*pow(m_rho,2)*pow(s,2)*pow(t2,2) - 3*pow(m_omega,2)*pow(s,2)*pow(t2,2) +
             5*pow(s,3)*pow(t2,2) + ((pow(m_omega,4) - 4*pow(m_omega,2)*s + 5*pow(s,2))*pow(t2,3))/3. -
             (pow(pow(m_omega,2) - s,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(m_rho,2) + pow(m_pi,4)*(pow(m_rho,4) + 4*pow(m_omega,4) - 2*pow(m_omega,2)*s) +
                  pow(m_omega,4)*(pow(m_rho,4) + pow(m_omega,4) + 2*pow(m_omega,2)*s + 2*pow(s,2) - 2*pow(m_rho,2)*(pow(m_omega,2) + s)) -
                  2*pow(m_pi,2)*pow(m_omega,2)*(pow(m_rho,4) + 2*pow(m_omega,2)*(pow(m_omega,2) + s) - pow(m_rho,2)*(2*pow(m_omega,2) + s))))/(pow(m_omega,2) - t1) -
             (pow(m_pi,8) - 2*pow(m_pi,6)*pow(m_rho,2) + 3*pow(m_omega,8) - 4*pow(m_omega,6)*s - 7*pow(m_omega,4)*pow(s,2) + 4*pow(m_omega,2)*pow(s,3) + 5*pow(s,4) +
                pow(m_rho,4)*(pow(m_omega,4) - 2*pow(m_omega,2)*s + 2*pow(s,2)) + pow(m_rho,2)*(-4*pow(m_omega,6) + 8*pow(m_omega,4)*s - 6*pow(s,3)) -
                2*pow(m_pi,2)*(4*pow(m_omega,6) - 2*pow(m_rho,2)*pow(pow(m_omega,2) - 2*s,2) + pow(m_rho,4)*s - 10*pow(m_omega,4)*s + 8*pow(s,3)) +
                pow(m_pi,4)*(pow(m_rho,4) + 2*pow(m_rho,2)*(pow(m_omega,2) - s) + 4*(pow(m_omega,4) - 3*pow(m_omega,2)*s + 3*pow(s,2))))*t1 +
             2*pow(m_pi,2)*pow(m_omega,4)*pow(t1,2) + pow(m_rho,2)*pow(m_omega,4)*pow(t1,2) - pow(m_omega,6)*pow(t1,2) + pow(m_pi,4)*s*pow(t1,2) -
             pow(m_pi,2)*pow(m_rho,2)*s*pow(t1,2) - 8*pow(m_pi,2)*pow(m_omega,2)*s*pow(t1,2) - 3*pow(m_rho,2)*pow(m_omega,2)*s*pow(t1,2) +
             2*pow(m_omega,4)*s*pow(t1,2) + 8*pow(m_pi,2)*pow(s,2)*pow(t1,2) + 3*pow(m_rho,2)*pow(s,2)*pow(t1,2) + 3*pow(m_omega,2)*pow(s,2)*pow(t1,2) -
             5*pow(s,3)*pow(t1,2) - ((pow(m_omega,4) - 4*pow(m_omega,2)*s + 5*pow(s,2))*pow(t1,3))/3. +
             2*(pow(m_omega,2) - s)*(-pow(m_pi,8) + pow(m_pi,4)*(4*pow(m_omega,4) - 7*pow(m_omega,2)*s + pow(s,2) + pow(m_rho,2)*(pow(m_omega,2) + s)) +
                pow(m_pi,2)*(-6*pow(m_omega,6) + 6*pow(m_omega,4)*s + 8*pow(m_omega,2)*pow(s,2) + pow(m_rho,4)*(-pow(m_omega,2) + s) +
                   pow(m_rho,2)*(4*pow(m_omega,4) - 7*pow(m_omega,2)*s - pow(s,2))) +
                pow(m_omega,2)*(2*pow(m_omega,6) + pow(m_rho,4)*(pow(m_omega,2) - s) - 4*pow(m_omega,2)*pow(s,2) - 3*pow(s,3) +
                   pow(m_rho,2)*(-3*pow(m_omega,4) + 2*pow(m_omega,2)*s + 3*pow(s,2))))*log((-pow(m_omega,2) + t2)/(-pow(m_omega,2) + t1))))/
                   (128.0*M_PI*pow(pow(m_omega,2) - s,2)*(pow(m_pi,4) + pow(pow(m_rho,2) - s,2) - 2*pow(m_pi,2)*(pow(m_rho,2) + s)));

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

float ScatterActionPhoton::diff_cross_section(float t, float m3, float t2, float t1) const {
  const float to_mb = 0.3894;
  const float m_rho = ParticleType::find(pdg::rho_z).mass();
  const float m_pi = ParticleType::find(pdg::pi_z).mass();
  float s = mandelstam_s();
  float diff_xsection = 0.0;

  const float Const = 0.059;
  const float g_POR = 11.93;
  const float ma1 = 1.26;
  const float ghat = 6.4483;
  const float eta1 = 2.3920;
  const float eta2 = 1.9430;
  const float delta = -0.6426;
  const float C4 = -0.14095;
  const float Gammaa1 = 0.4;
  const float Pi = M_PI;
  float m_omega = 0.783;
  float momega = m_omega;
  float mrho = m_rho;
  float mpion = m_pi;

  switch (reac) {
    case ReactionType::pi_pi:
      if (outgoing_particles_[0].type().pdgcode().is_rho()) {
        diff_xsection = (pow(Const,2)*pow(ghat,4)*((0.25*(32*pow(C4,2)*pow(mrho,8) + 2*pow(delta,2)*pow(s,2) + 8*C4*pow(mrho,6)*(-6 + delta - 8*C4*s) +
                    2*delta*pow(mrho,2)*s*(-6 + delta - 8*C4*s) +
                    pow(mrho,4)*(12 - pow(delta,2) + 8*C4*(6 + delta)*s + 32*pow(C4,2)*pow(s,2))))/pow(mrho,4) -
               (0.25*pow(-2 + delta,2)*pow(mpion,2)*(pow(mpion,4) + pow(pow(mrho,2) - t,2) - 2*pow(mpion,2)*(pow(mrho,2) + t)))/
                (pow(mrho,2)*pow(pow(mpion,2) - t,2)) - (0.25*pow(-2 + delta,2)*pow(mpion,2)*
                  (pow(mpion,4) + pow(s + t,2) - 2*pow(mpion,2)*(2*pow(mrho,2) + s + t)))/
                (pow(mrho,2)*pow(pow(mpion,2) + pow(mrho,2) - s - t,2)) +
               (0.125*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - 2*pow(mpion,2) - pow(mrho,2) + s + t)*
                  (eta1*(2*pow(mpion,2) - s) + eta2*(-3*pow(mpion,2) - pow(mrho,2) + s + t))*
                  (pow(mpion,4) + t*(-pow(mrho,2) + 2*s + t) - pow(mpion,2)*(pow(mrho,2) + 2*t)))/
                ((-pow(mpion,2) + t)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - 2*pow(mpion,2) - pow(mrho,2) + s + t,2))) +
               (0.25*(-2. + delta)*(pow(mpion,4)*(2. + delta - 8.*C4*pow(mrho,2)) + 8.*C4*pow(mrho,4)*t +
                    t*((2. + 3.*delta)*s + (2. + delta)*t) + pow(mrho,2)*(s*(2. - 1.*delta - 16.*C4*t) + t*(-2. - 1.*delta - 8.*C4*t)) +
                    pow(mpion,2)*(8.*C4*pow(mrho,4) + (-2. + delta)*s + (-4. - 2.*delta)*t + pow(mrho,2)*(-6. + delta + 16.*C4*t))))/
                (pow(mrho,2)*(pow(mpion,2) - 1.*t)) - (0.125*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - 2*pow(mpion,2) - pow(mrho,2) + s + t)*
                  (-(eta2*(3*pow(mpion,2) + pow(mrho,2) - s - t)*
                       (pow(mpion,4) + (pow(mrho,2) - s - t)*(s - t) - pow(mpion,2)*(pow(mrho,2) - 2*s + 2*t))) +
                    eta1*(2*pow(mpion,6) + pow(mpion,4)*(-2*pow(mrho,2) + 5*s - 4*t) + s*(s + t)*(-pow(mrho,2) + s + t) +
                       pow(mpion,2)*(2*pow(mrho,4) + pow(mrho,2)*(s - 2*t) - 2*(2*s - t)*(s + t)))))/
                ((-pow(mpion,2) - pow(mrho,2) + s + t)*(pow(Gammaa1,2)*pow(ma1,2) +
                    pow(pow(ma1,2) - 2*pow(mpion,2) - pow(mrho,2) + s + t,2))) +
               (0.03125*pow(eta1 - eta2,2)*(-2*eta1*eta2*(pow(mpion,8) - 4*pow(mpion,6)*t +
                       pow(t,2)*(-pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) + 2*s*t + pow(t,2)) -
                       2*pow(mpion,2)*t*(-2*pow(mrho,4) + pow(mrho,2)*s + 2*t*(s + t)) + pow(mpion,4)*(-pow(mrho,4) + 2*t*(s + 3*t))) +
                    pow(eta2,2)*(pow(mpion,8) - 2*pow(mpion,6)*(pow(mrho,2) + 2*t) +
                       pow(t,2)*(pow(mrho,4) + 2*pow(s,2) + 2*s*t + pow(t,2) + 2*pow(mrho,2)*(-s + t)) -
                       2*pow(mpion,2)*t*(2*t*(s + t) + pow(mrho,2)*(s + 3*t)) + pow(mpion,4)*(pow(mrho,4) + 6*pow(mrho,2)*t + 2*t*(s + 3*t)))
                      + pow(eta1,2)*(pow(mpion,8) + 2*pow(mpion,6)*(pow(mrho,2) - 2*t) -
                       2*pow(mpion,2)*(pow(mrho,2) - s - t)*(pow(mrho,4) + pow(mrho,2)*t - 2*pow(t,2)) +
                       t*(-pow(mrho,2) + t)*(2*pow(s,2) + 2*s*t + pow(t,2) - pow(mrho,2)*(2*s + t)) +
                       pow(mpion,4)*(pow(mrho,4) - 2*pow(mrho,2)*(s + 3*t) + 2*t*(s + 3*t)))))/pow(pow(ma1,2) - t,2) +
               (2*((0.125*pow(-2 + delta,2)*(2*pow(mpion,2) - s)*
                       (pow(mpion,4) + pow(mrho,2)*(s - t) + t*(s + t) - pow(mpion,2)*(3*pow(mrho,2) + s + 2*t)))/
                     ((pow(mpion,2) - t)*(pow(mpion,2) + pow(mrho,2) - s - t)) -
                    (0.125*(-2. + delta)*(pow(mpion,4)*(2. + delta - 8.*C4*pow(mrho,2)) - 2.*delta*pow(s,2) + 2.*s*t - 1.*delta*s*t +
                         2.*pow(t,2) + delta*pow(t,2) + C4*pow(mrho,4)*(-8.*s + 8.*t) +
                         pow(mrho,2)*((2. + delta)*s + 8.*C4*pow(s,2) + t*(-2. - 1.*delta - 8.*C4*t)) +
                         pow(mpion,2)*(8.*C4*pow(mrho,4) - 2.*s + 5.*delta*s - 4.*t - 2.*delta*t +
                            pow(mrho,2)*(-6. + delta - 16.*C4*s + 16.*C4*t))))/(pow(mpion,2) + pow(mrho,2) - 1.*s - 1.*t)))/pow(mrho,2) +
               (0.03125*pow(eta1 - eta2,2)*(-2*eta1*eta2*(pow(mpion,8) + 4*pow(mpion,6)*(pow(mrho,2) - t) +
                       pow(-pow(mrho,2) + s + t,2)*(pow(s,2) + pow(t,2) - 2*pow(mrho,2)*(s + t)) +
                       pow(mpion,4)*(9*pow(mrho,4) + 4*pow(s,2) + 2*s*t + 6*pow(t,2) - 2*pow(mrho,2)*(7*s + 6*t)) +
                       2*pow(mpion,2)*(pow(mrho,2) - s - t)*(2*pow(mrho,4) - pow(mrho,2)*(5*s + 4*t) + 2*(pow(s,2) + pow(t,2)))) +
                    pow(eta2,2)*(pow(mpion,8) + pow(mpion,6)*(6*pow(mrho,2) - 4*t) +
                       pow(-pow(mrho,2) + s + t,2)*(4*pow(mrho,4) + pow(s,2) + pow(t,2) - 4*pow(mrho,2)*(s + t)) +
                       pow(mpion,4)*(17*pow(mrho,4) + 4*pow(s,2) + 2*s*t + 6*pow(t,2) - 2*pow(mrho,2)*(10*s + 9*t)) +
                       2*pow(mpion,2)*(pow(mrho,2) - s - t)*(7*pow(mrho,4) - pow(mrho,2)*(8*s + 7*t) + 2*(pow(s,2) + pow(t,2)))) +
                    pow(eta1,2)*(pow(mpion,8) + 2*pow(mpion,6)*(pow(mrho,2) - 2*t) +
                       (s + t)*(-pow(mrho,2) + s + t)*(pow(s,2) + pow(t,2) - pow(mrho,2)*(s + t)) +
                       pow(mpion,4)*(5*pow(mrho,4) + 4*pow(s,2) + 2*s*t + 6*pow(t,2) - 2*pow(mrho,2)*(5*s + 3*t)) -
                       2*pow(mpion,2)*(2*pow(mrho,4)*(s + t) + 2*(s + t)*(pow(s,2) + pow(t,2)) -
                          pow(mrho,2)*(4*pow(s,2) + 5*s*t + 3*pow(t,2))))))/
                (pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - 2*pow(mpion,2) - pow(mrho,2) + s + t,2)) +
               (0.0625*pow(eta1 - eta2,2)*(pow(ma1,2) - 2*pow(mpion,2) - pow(mrho,2) + s + t)*
                  (-(pow(eta2,2)*(pow(mpion,8) + 2*pow(mpion,6)*(pow(mrho,2) - 2*t) +
                         2*pow(mpion,2)*t*(pow(pow(mrho,2) - s,2) + (3*pow(mrho,2) - 2*s)*t - 2*pow(t,2)) +
                         (pow(mrho,2) - s - t)*t*(2*pow(mrho,4) + pow(s,2) - s*t - pow(t,2) + pow(mrho,2)*(-3*s + t)) +
                         pow(mpion,4)*(pow(mrho,4) + 2*t*(s + 3*t) - pow(mrho,2)*(s + 6*t)))) -
                    pow(eta1,2)*(pow(mpion,8) + 2*pow(mpion,6)*(pow(mrho,2) - 2*t) +
                       (pow(mrho,2) - s - t)*t*(pow(s,2) - s*t - pow(t,2) + pow(mrho,2)*(s + t)) +
                       pow(mpion,4)*(3*pow(mrho,4) + 2*t*(s + 3*t) - pow(mrho,2)*(5*s + 6*t)) +
                       2*pow(mpion,2)*(-(pow(mrho,4)*(s + t)) + t*(pow(s,2) - 2*s*t - 2*pow(t,2)) +
                          pow(mrho,2)*(pow(s,2) + 2*s*t + 3*pow(t,2)))) +
                    2*eta1*eta2*(pow(mpion,8) + 2*pow(mpion,6)*(pow(mrho,2) - 2*t) -
                       (pow(mrho,2) - s - t)*t*(pow(mrho,4) - pow(s,2) - pow(mrho,2)*t + t*(s + t)) +
                       2*pow(mpion,4)*(2*pow(mrho,4) + t*(s + 3*t) - pow(mrho,2)*(2*s + 3*t)) +
                       pow(mpion,2)*(pow(mrho,6) - 2*pow(mrho,4)*(s + 2*t) + 2*t*(pow(s,2) - 2*s*t - 2*pow(t,2)) +
                          pow(mrho,2)*(pow(s,2) + 2*s*t + 6*pow(t,2))))))/
                ((-pow(ma1,2) + t)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - 2*pow(mpion,2) - pow(mrho,2) + s + t,2))) +
               2*((-0.0625*(-2 + delta)*(eta1 - eta2)*(eta1*(-2*pow(mpion,2) + s) + eta2*(pow(mpion,2) + t))*
                     (-pow(mpion,4) + pow(s,2) - pow(t,2) + pow(mrho,2)*(-s + t) + pow(mpion,2)*(pow(mrho,2) - 2*s + 2*t)))/
                   ((pow(mpion,2) + pow(mrho,2) - s - t)*(-pow(ma1,2) + t)) +
                  (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-0.5*pow(mpion,6) + pow(mpion,4)*(0.5*pow(mrho,2) + 0.5*t) +
                          pow(mpion,2)*(1.*pow(mrho,2) - 1.*s + 0.5*t)*t + (0.5*pow(mrho,2) - 1.*s - 0.5*t)*pow(t,2)) +
                       eta1*(1.*pow(mpion,6) + pow(mpion,4)*(-1.*pow(mrho,2) + 0.5*s - 2.*t) + s*(-0.5*pow(mrho,2) + 0.5*t)*t +
                          pow(mpion,2)*(1.*pow(mrho,4) + pow(mrho,2)*(-0.5*s - 1.*t) + t*(1.*s + 1.*t)))))/
                   ((pow(ma1,2) - 1.*t)*(-1.*pow(mpion,2) + t)) +
                  (0.0625*(eta1 - eta2)*(eta2*(8*C4*pow(mrho,6)*t - 2*delta*pow(s,2)*t +
                          pow(mrho,2)*(-4*pow(mpion,4) + (s*(2 + 3*delta + 8*C4*s) - 4*t)*t + pow(mpion,2)*(-((-2 + delta)*s) + 8*t)) +
                          pow(mrho,4)*(8*C4*pow(mpion,4) - pow(mpion,2)*(-2 + delta + 16*C4*t) + t*(-6 + delta + 8*C4*(-2*s + t)))) +
                       eta1*(2*delta*pow(s,2)*t + 8*C4*pow(mrho,6)*(-2*pow(mpion,2) + t) -
                          pow(mrho,2)*(-4*pow(mpion,4) - 4*pow(t,2) + 2*pow(mpion,2)*((2 + delta)*s + 4*t) +
                             pow(s,2)*(-2 + delta + 8*C4*t)) + pow(mrho,4)*
                           (-8*C4*pow(mpion,4) + (-2 + delta)*s - 4*t*(1 + 2*C4*t) + 8*pow(mpion,2)*(1 + 2*C4*(s + t))))))/
                   (pow(mrho,2)*(-pow(ma1,2) + t))) - (0.125*(eta1 - eta2)*(pow(ma1,2) - 2*pow(mpion,2) - pow(mrho,2) + s + t)*
                  (eta1*(pow(mpion,4)*(4*pow(mrho,2) - 8*C4*pow(mrho,4)) + 8*C4*pow(mrho,6)*(s + t) - 2*delta*pow(s,2)*(s + t) +
                       pow(mrho,2)*((6 + delta)*pow(s,2) + 8*s*t + 4*pow(t,2) + 8*C4*pow(s,2)*(s + t)) -
                       pow(mrho,4)*(-((-6 + delta)*s) + 4*t + 8*C4*(2*pow(s,2) + 2*s*t + pow(t,2))) +
                       2*pow(mpion,2)*(-8*C4*pow(mrho,6) + 2*delta*pow(s,2) - pow(mrho,2)*(s*(6 + delta + 8*C4*s) + 4*t) +
                          4*pow(mrho,4)*(1 + 2*C4*(2*s + t)))) +
                    eta2*(pow(mpion,4)*(-4*pow(mrho,2) + 8*C4*pow(mrho,4)) -
                       (-pow(mrho,2) + s + t)*(16*C4*pow(mrho,6) - 2*delta*pow(s,2) + pow(mrho,2)*(s*(6 + 3*delta + 8*C4*s) + 4*t) +
                          pow(mrho,4)*(-10 + delta - 8*C4*(3*s + t))) +
                       pow(mpion,2)*(32*C4*pow(mrho,6) - 4*delta*pow(s,2) + pow(mrho,2)*(s*(14 + 5*delta + 16*C4*s) + 8*t) +
                          pow(mrho,4)*(delta - 2*(9 + 8*C4*(3*s + t)))))))/
                (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - 2*pow(mpion,2) - pow(mrho,2) + s + t,2)))))/
           (16.*M_PI*s*(-4*pow(mpion,2) + s));
    //  } else if (outgoing_particles_[0].type().pdgcode() == pdg::eta) {
    //    diff_xsection = to_be_determined;
      } else if (outgoing_particles_[0].type().pdgcode() == pdg::photon) {
        diff_xsection = 0.0000000000001/to_mb/(t2-t1);
      }
      break;
    case ReactionType::pi0_pi:
      if (outgoing_particles_[0].type().pdgcode().is_rho()) {
        diff_xsection = (pow(Const,2)*pow(ghat,4)*((-0.25*pow(-2 + delta,2)*pow(mpion,2)*
      				  (pow(mpion,4) + pow(pow(mrho,2) - t,2) - 2*pow(mpion,2)*(pow(mrho,2) + t)))/(pow(mrho,2)*pow(pow(mpion,2) - t,2)) +
      			   (0.0625*pow(-2*pow(mrho,2) + delta*s,2)*(8*pow(mpion,4)*pow(mrho,2) + 2*pow(mrho,6) + pow(s,3) + 8*pow(mrho,2)*t*(s + t) -
      					pow(mrho,4)*(7*s + 8*t) + 4*pow(mpion,2)*(5*pow(mrho,4) - pow(s,2) - 4*pow(mrho,2)*t)))/
      				(pow(mrho,6)*pow(pow(mrho,2) - s,2)) - (0.0625*(eta1 - eta2)*(2*pow(mrho,2) - delta*s)*
      				  (-(eta2*(2*pow(mpion,2) + pow(mrho,2) - s - 2*t)*
      					   (2*pow(mpion,4) + pow(mpion,2)*(-pow(mrho,2) + s - 4*t) + t*(3*pow(mrho,2) + s + 2*t))) +
      					eta1*(4*pow(mpion,6) - pow(mrho,4)*s + pow(s,3) + 2*pow(mpion,4)*(5*pow(mrho,2) - s - 6*t) -
      					   2*(pow(mrho,4) - 4*pow(mrho,2)*s + pow(s,2))*t + 6*(pow(mrho,2) - s)*pow(t,2) - 4*pow(t,3) -
      					   4*pow(mpion,2)*(4*pow(mrho,2)*t + (s - 3*t)*(s + t)))))/(pow(mrho,2)*(pow(mrho,2) - s)*(pow(ma1,2) - t)) -
      			   (0.125*(-2 + delta)*(eta1 - eta2)*(-(eta2*(pow(mpion,2) + t)*
      					   (pow(mpion,4) + t*(-pow(mrho,2) + 2*s + t) - pow(mpion,2)*(pow(mrho,2) + 2*t))) +
      					eta1*(2*pow(mpion,6) + pow(mpion,4)*(-2*pow(mrho,2) + s - 4*t) + s*t*(-pow(mrho,2) + t) +
      					   pow(mpion,2)*(2*pow(mrho,4) + 2*t*(s + t) - pow(mrho,2)*(s + 2*t)))))/((-pow(ma1,2) + t)*(-pow(mpion,2) + t)) +
      			   (0.03125*pow(eta1 - eta2,2)*(-2*eta1*eta2*(pow(mpion,8) - 4*pow(mpion,6)*t +
      					   pow(t,2)*(-pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) + 2*s*t + pow(t,2)) -
      					   2*pow(mpion,2)*t*(-2*pow(mrho,4) + pow(mrho,2)*s + 2*t*(s + t)) + pow(mpion,4)*(-pow(mrho,4) + 2*t*(s + 3*t))) +
      					pow(eta2,2)*(pow(mpion,8) - 2*pow(mpion,6)*(pow(mrho,2) + 2*t) +
      					   pow(t,2)*(pow(mrho,4) + 2*pow(s,2) + 2*s*t + pow(t,2) + 2*pow(mrho,2)*(-s + t)) -
      					   2*pow(mpion,2)*t*(2*t*(s + t) + pow(mrho,2)*(s + 3*t)) + pow(mpion,4)*(pow(mrho,4) + 6*pow(mrho,2)*t + 2*t*(s + 3*t))) +
      					pow(eta1,2)*(pow(mpion,8) + 2*pow(mpion,6)*(pow(mrho,2) - 2*t) -
      					   2*pow(mpion,2)*(pow(mrho,2) - s - t)*(pow(mrho,4) + pow(mrho,2)*t - 2*pow(t,2)) +
      					   t*(-pow(mrho,2) + t)*(2*pow(s,2) + 2*s*t + pow(t,2) - pow(mrho,2)*(2*s + t)) +
      					   pow(mpion,4)*(pow(mrho,4) - 2*pow(mrho,2)*(s + 3*t) + 2*t*(s + 3*t)))))/pow(pow(ma1,2) - t,2) -
      			   (3.*(1.*pow(mrho,2) - 0.5*delta*s)*(delta*(0.666667*pow(mpion,4)*pow(mrho,2) + 0.166667*pow(mrho,6) +
      					   pow(mpion,2)*(1.66667*pow(mrho,4) - 0.416667*pow(s,2) + pow(mrho,2)*(0.0833333*s - 1.33333*t)) +
      					   pow(mrho,4)*(-0.541667*s - 0.666667*t) + pow(s,2)*(0.125*s + 0.0833333*t) +
      					   pow(mrho,2)*(-0.0833333*pow(s,2) + 0.583333*s*t + 0.666667*pow(t,2))) +
      					pow(mrho,2)*(1.*C4*pow(mrho,6) + pow(mpion,2)*(2.*C4*pow(mrho,4) + 0.166667*s + pow(mrho,2)*(-0.833333 - 0.666667*C4*s)) +
      					   s*(-0.0833333*s - 0.166667*t) + pow(mrho,4)*(-0.416667 - 1.33333*C4*s - 2.*C4*t) +
      					   pow(mrho,2)*(0.833333*t + s*(0.5 + 0.333333*C4*s + 0.666667*C4*t)))))/(pow(mrho,8) - 1.*pow(mrho,6)*s) +
      			   (pow(mrho,6)*(0.75 + C4*(2.*C4*pow(mrho,4) + pow(mrho,2)*(-3. - 4.*C4*s) + s*(3. + 2.*C4*s))) +
      				  pow(delta,2)*(0.5*pow(mpion,4)*pow(mrho,2) + 0.125*pow(mrho,6) +
      					 pow(mpion,2)*(1.25*pow(mrho,4) - 0.375*pow(s,2) + pow(mrho,2)*(0.125*s - 1.*t)) + pow(mrho,4)*(-0.375*s - 0.5*t) +
      					 pow(s,2)*(0.125*s + 0.125*t) + pow(mrho,2)*(0.0625*pow(s,2) + 0.375*s*t + 0.5*pow(t,2))) +
      				  delta*pow(mrho,2)*(2.*C4*pow(mrho,6) + pow(mpion,2)*(3.*C4*pow(mrho,4) + 0.25*s + pow(mrho,2)*(-1.25 - 1.*C4*s)) +
      					 s*(-0.25*s - 0.25*t) + pow(mrho,4)*(-0.75 - 1.5*C4*s - 3.*C4*t) + pow(mrho,2)*(1.25*t + s*(0.25 - 0.5*C4*s + 1.*C4*t))))/
      				pow(mrho,6) + (2*((-0.0625*(-2. + delta)*(-2.*pow(mrho,2) + delta*s)*
      					   (pow(mpion,4)*(4.*pow(mrho,2) + 4.*s) + pow(mpion,2)*(pow(mrho,2)*(-7.*s - 4.*t) + s*(-1.*s - 4.*t)) +
      						 s*(pow(mrho,4) + pow(mrho,2)*(s - 1.*t) + s*t)))/((pow(mrho,2) - 1.*s)*(pow(mpion,2) - 1.*t)) +
      					(0.0625*(-2 + delta)*(pow(mpion,4)*((-2 + 4*delta)*pow(mrho,2) + 8*C4*pow(mrho,4) + 5*delta*s) - 8*C4*pow(mrho,6)*t +
      						 delta*s*t*(s + t) + pow(mrho,2)*(delta*s*(s - 3*t) - 2*t*(s + t)) + 2*pow(mrho,4)*((-1 + delta)*s + t + 4*C4*t*(2*s + t)) -
      						 pow(mpion,2)*(8*C4*pow(mrho,6) + delta*s*(s + 6*t) + 2*pow(mrho,4)*(-3 + 8*C4*t) +
      							pow(mrho,2)*((-2 + 9*delta)*s + 4*(-1 + delta)*t))))/(-pow(mpion,2) + t)))/pow(mrho,4) -
      			   (0.0625*(eta1 - eta2)*(eta2*(-4*delta*pow(mpion,6) + 4*pow(mpion,4)*(pow(mrho,2) - 2*C4*pow(mrho,4) + 3*delta*t) +
      					   pow(mpion,2)*(delta*(s - 6*t)*(s + 2*t) - (2 + delta)*pow(mrho,2)*(s + 4*t) + 2*pow(mrho,4)*(-1 + delta + 8*C4*t)) +
      					   t*(-8*C4*pow(mrho,6) + pow(mrho,4)*(6 - 4*delta + 16*C4*s - 8*C4*t) +
      						  pow(mrho,2)*(-(s*(2 + delta + 8*C4*s)) + 4*(1 + delta)*t) + delta*(3*pow(s,2) + 4*s*t + 4*pow(t,2)))) +
      					eta1*(4*delta*pow(mpion,6) - 8*C4*pow(mrho,6)*t + delta*(pow(s,3) - 4*pow(s,2)*t - 6*s*pow(t,2) - 4*pow(t,3)) -
      					   2*pow(mpion,4)*((2 - 5*delta)*pow(mrho,2) - 4*C4*pow(mrho,4) + delta*(s + 6*t)) +
      					   2*pow(mrho,4)*(s - delta*s + t*(2 - delta + 4*C4*t)) +
      					   pow(mrho,2)*(8*delta*s*t + 2*(-2 + 3*delta)*pow(t,2) + pow(s,2)*(-2 + delta + 8*C4*t)) -
      					   2*pow(mpion,2)*(-8*C4*pow(mrho,6) + 2*delta*(s - 3*t)*(s + t) - pow(mrho,2)*((2 + delta)*s + (4 - 8*delta)*t) +
      						  pow(mrho,4)*(4 + 8*C4*(s + t))))))/(pow(mrho,2)*(-pow(ma1,2) + t))))/(16.*Pi*s*(-4*pow(mpion,2) + s));

      } else if (outgoing_particles_[0].type().pdgcode().is_pion()) {
        diff_xsection = 0.0000000000001/to_mb/(t2-t1);
      }
      break;
    case ReactionType::pi_rho0:

      diff_xsection = 1/3.0*(pow(Const,2)*pow(ghat,4)*((-8*pow(-2 + delta,2)*pow(m_pi,2))/(pow(mrho,2)*pow(pow(m_pi,2) - s,2)) -
       (8*pow(-2 + delta,2)*pow(m_pi,2)*(pow(m_pi,4) + pow(pow(mrho,2) - t,2) - 2*pow(m_pi,2)*(pow(mrho,2) + t)))/
        (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*pow(pow(m_pi,2) - t,2)) +
       (4*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(-(eta2*(pow(m_pi,2) + s)) + eta1*(-pow(mrho,2) + s + t))*
          (-pow(m_pi,4) + pow(m_pi,2)*(pow(mrho,2) - 2*t) + t*(-pow(mrho,2) + 2*s + t)))/
        ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*
          (pow(m_pi,2) - t)) - (8*(-2 + delta)*(pow(m_pi,4)*(2 - 3*delta + 8*C4*pow(mrho,2)) + pow(mrho,4)*(-2 + delta + 8*C4*t) +
            t*((2 + 3*delta)*s + 2*delta*t) + pow(m_pi,2)*(-8*C4*pow(mrho,4) + (-2 + delta)*s - (2 + 3*delta)*t + 4*pow(mrho,2)*(1 + 4*C4*t)) -
            pow(mrho,2)*(t*(-2 + 3*delta + 8*C4*t) + s*(-2 + delta + 16*C4*t))))/
        (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*(pow(m_pi,2) - t)) +
       (4*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(eta2*(pow(m_pi,2) + s)*
             (pow(m_pi,4) - pow(m_pi,2)*(pow(mrho,2) - 2*s) + s*(pow(mrho,2) - s - 2*t)) +
            eta1*(-4*pow(m_pi,6) + s*(-pow(mrho,2) + s)*(-pow(mrho,2) + s + t) + pow(m_pi,4)*(3*pow(mrho,2) + s + t) -
               pow(m_pi,2)*(pow(mrho,4) + 2*s*(s - t) + pow(mrho,2)*(-s + t)))))/
        ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,2) - s)*
          (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
       (pow(eta1 - eta2,2)*(pow(eta2,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + pow(m_pi,4)*(pow(pow(mrho,2) + 2*s,2) - 2*s*t) +
               pow(s,2)*(pow(pow(mrho,2) + s,2) + 2*(-pow(mrho,2) + s)*t + 2*pow(t,2)) -
               2*pow(m_pi,2)*s*(pow(mrho,4) + pow(mrho,2)*(2*s - t) + 2*s*(s + t))) +
            2*eta1*eta2*(-pow(m_pi,8) + pow(m_pi,4)*(pow(mrho,4) + 2*pow(mrho,2)*s + 2*s*(-2*s + t)) -
               2*pow(m_pi,2)*s*(pow(mrho,4) + pow(mrho,2)*(s + t) - 2*s*(s + t)) + pow(s,2)*(pow(mrho,4) - pow(s,2) + 2*pow(mrho,2)*t - 2*t*(s + t)))\
             + pow(eta1,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + pow(m_pi,4)*(3*pow(mrho,4) + 2*s*(2*s - t) + 2*pow(mrho,2)*(-3*s + t)) -
               2*pow(m_pi,2)*(pow(mrho,2) - s)*(-2*s*(s + t) + pow(mrho,2)*(2*s + t)) +
               s*(-pow(mrho,2) + s)*(pow(s,2) + 2*s*t + 2*pow(t,2) - pow(mrho,2)*(s + 2*t)))))/
        ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
       (pow(eta1 - eta2,2)*(-2*eta1*eta2*(pow(m_pi,8) - pow(m_pi,4)*(pow(mrho,4) + 2*(pow(mrho,2) + s)*t - 4*pow(t,2)) +
               pow(t,2)*(-pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) + 2*s*t + pow(t,2)) +
               2*pow(m_pi,2)*t*(pow(mrho,4) + pow(mrho,2)*(s + t) - 2*t*(s + t))) +
            pow(eta2,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + pow(m_pi,4)*(pow(mrho,4) + 4*pow(mrho,2)*t - 2*(s - 2*t)*t) +
               pow(t,2)*(pow(mrho,4) + 2*pow(s,2) + 2*s*t + pow(t,2) + 2*pow(mrho,2)*(-s + t)) -
               2*pow(m_pi,2)*t*(pow(mrho,4) - pow(mrho,2)*(s - 2*t) + 2*t*(s + t))) +
            pow(eta1,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + pow(m_pi,4)*(3*pow(mrho,4) + 2*pow(mrho,2)*(s - 3*t) - 2*(s - 2*t)*t) +
               t*(-pow(mrho,2) + t)*(2*pow(s,2) + 2*s*t + pow(t,2) - pow(mrho,2)*(2*s + t)) -
               2*pow(m_pi,2)*(-pow(mrho,2) + t)*(2*t*(s + t) - pow(mrho,2)*(s + 2*t)))))/
        ((pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*pow(pow(ma1,2) - t,2)) +
       (8*(-2 + delta)*((-2 + delta)*pow(mrho,6) + pow(m_pi,6)*(-2 + 3*delta - 8*C4*pow(mrho,2)) + s*t*((-2 + 3*delta)*s + 4*delta*t) +
            pow(m_pi,4)*(8*C4*pow(mrho,4) + 4*delta*s + 2*t - 3*delta*t - pow(mrho,2)*(2 + delta + 16*C4*s - 8*C4*t)) +
            pow(mrho,4)*(-((-2 + delta)*t) + s*(4 - 2*delta + 8*C4*t)) + pow(mrho,2)*s*(s*(-2 + delta - 8*C4*t) - 2*t*(delta + 8*C4*t)) +
            pow(m_pi,2)*(s*((2 - 3*delta)*s - 8*delta*t) - pow(mrho,4)*(-6 + 3*delta + 8*C4*(s + t)) +
               pow(mrho,2)*(8*C4*pow(s,2) + 4*(-1 + delta)*t + s*(-8 + 6*delta + 32*C4*t)))))/
        (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*(pow(m_pi,2) - t)) +
       (2*pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(pow(eta1,2)*(pow(m_pi,8) + pow(m_pi,4)*(2*pow(mrho,4) + 2*s*t - 3*pow(mrho,2)*(s + t)) +
               s*t*(2*pow(mrho,4) + pow(s,2) + 3*s*t + pow(t,2) - 3*pow(mrho,2)*(s + t)) -
               2*pow(m_pi,2)*(pow(mrho,2) - s - t)*(-2*s*t + pow(mrho,2)*(s + t))) +
            pow(eta2,2)*(pow(m_pi,8) - 4*pow(m_pi,2)*s*t*(pow(mrho,2) + s + t) + pow(m_pi,4)*(2*s*t + pow(mrho,2)*(s + t)) +
               s*t*(pow(s,2) + 3*s*t + pow(t,2) + pow(mrho,2)*(s + t))) +
            2*eta1*eta2*(-pow(m_pi,8) + 2*pow(m_pi,6)*pow(mrho,2) - 2*pow(m_pi,4)*s*t - s*t*(pow(s,2) + 3*s*t + pow(t,2) - 2*pow(mrho,2)*(s + t)) -
               pow(m_pi,2)*(-4*s*t*(s + t) + pow(mrho,2)*(pow(s,2) + 4*s*t + pow(t,2))))))/
        ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*
          (pow(ma1,2) - t)) + (8*(pow(delta,2)*(8*pow(m_pi,4) + 3*pow(mrho,4) - 6*pow(mrho,2)*(s + t) + 2*pow(s + t,2) +
               4*pow(m_pi,2)*(3*pow(mrho,2) - 2*(s + t))) - 4*delta*pow(mrho,2)*
             (16*C4*pow(m_pi,4) + pow(mrho,2)*(3 - 6*C4*(s + t)) + (s + t)*(-3 + 4*C4*(s + t)) + 2*pow(m_pi,2)*(3 + C4*(6*pow(mrho,2) - 8*(s + t)))) +
            4*pow(mrho,4)*(3 + 4*C4*(2*pow(m_pi,2) - s - t)*(3 + C4*(4*pow(m_pi,2) - 2*(s + t))))))/
        (pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
       (4*(eta1 - eta2)*(-pow(ma1,2) + s)*(eta2*(-2*pow(m_pi,4)*(delta - 4*C4*pow(mrho,2))*(pow(mrho,2) + 4*s) +
               pow(m_pi,2)*(-2*pow(mrho,4)*(-2 + delta + 8*C4*s) + 8*delta*s*(s + t) - pow(mrho,2)*((-10 + delta)*s - (-2 + delta)*t + 32*C4*s*(s + t))) +
               s*(2*pow(mrho,4)*(-2 + delta + 4*C4*s) - 2*delta*pow(s + t,2) + pow(mrho,2)*((-6 + delta)*s + (-2 + delta)*t + 8*C4*pow(s + t,2)))) +
            eta1*(4*pow(m_pi,4)*(6*C4*pow(mrho,4) + 2*delta*s + pow(mrho,2)*(1 - 2*delta - 8*C4*s)) + 2*delta*s*pow(s + t,2) -
               pow(mrho,2)*((-6 + 5*delta)*pow(s,2) + 2*(-2 + 3*delta)*s*t + (-2 + delta)*pow(t,2) + 8*C4*s*pow(s + t,2)) +
               pow(mrho,4)*((-2 + delta)*(3*s + t) + 8*C4*s*(s + 2*t)) -
               2*pow(m_pi,2)*(4*delta*s*(s + t) - pow(mrho,2)*(-6*s + 7*delta*s - 2*t + 3*delta*t + 16*C4*s*(s + t)) +
                  2*pow(mrho,4)*(-2 + delta + 4*C4*(2*s + t))))))/
        (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
          (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
       (4*(eta1 - eta2)*(((-2 + delta)*(pow(m_pi,4) - pow(m_pi,2)*(pow(mrho,2) - 2*s) + s*(pow(mrho,2) - s - 2*t))*
               (eta1*(pow(mrho,2) - s - t) + eta2*(pow(m_pi,2) + t)))/((pow(m_pi,2) - s)*(pow(ma1,2) - t)) +
            ((-2 + delta)*(eta2*(pow(m_pi,2) + t)*(pow(m_pi,4) - pow(m_pi,2)*(pow(mrho,2) - 2*t) + (pow(mrho,2) - 2*s - t)*t) +
                 eta1*(-4*pow(m_pi,6) + (pow(mrho,2) - t)*(pow(mrho,2) - s - t)*t + pow(m_pi,4)*(3*pow(mrho,2) + s + t) -
                    pow(m_pi,2)*(pow(mrho,4) + pow(mrho,2)*(s - t) + 2*t*(-s + t)))))/((-pow(ma1,2) + t)*(-pow(m_pi,2) + t)) +
            (eta2*(-2*pow(m_pi,4)*(delta - 4*C4*pow(mrho,2))*(pow(mrho,2) + 4*t) +
                  pow(m_pi,2)*(8*delta*t*(s + t) - 2*pow(mrho,4)*(-2 + delta + 8*C4*t) -
                     pow(mrho,2)*(-((-2 + delta)*s) + (-10 + delta)*t + 32*C4*t*(s + t))) +
                  t*(-2*delta*pow(s + t,2) + 2*pow(mrho,4)*(-2 + delta + 4*C4*t) + pow(mrho,2)*((-2 + delta)*s + (-6 + delta)*t + 8*C4*pow(s + t,2)))) +
               eta1*(2*delta*t*pow(s + t,2) - pow(mrho,2)*((-2 + delta)*pow(s,2) + 2*(-2 + 3*delta)*s*t + (-6 + 5*delta)*pow(t,2) + 8*C4*t*pow(s + t,2)) +
                  pow(mrho,4)*(8*C4*t*(2*s + t) + (-2 + delta)*(s + 3*t)) +
                  4*pow(m_pi,4)*(6*C4*pow(mrho,4) + 2*delta*t + pow(mrho,2)*(1 - 2*delta - 8*C4*t)) -
                  2*pow(m_pi,2)*(4*delta*t*(s + t) - pow(mrho,2)*(-2*s + 3*delta*s - 6*t + 7*delta*t + 16*C4*t*(s + t)) +
                     2*pow(mrho,4)*(-2 + delta + 4*C4*(s + 2*t)))))/(pow(mrho,2)*(-pow(ma1,2) + t))))/
        (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))))/(512.*Pi);

      break;
    case ReactionType::pi_rho:
      diff_xsection = (pow(Const,2)*pow(ghat,4)*((-0.25*pow(-2 + delta,2)*pow(mpion,2)*
                    (pow(mpion,4) + pow(pow(mrho,2) - t,2) - 2*pow(mpion,2)*(pow(mrho,2) + t)))/(pow(mrho,2)*pow(pow(mpion,2) - t,2)) -
                 (0.0625*(eta1 - eta2)*(2*pow(mrho,2) + delta*(-2*pow(mpion,2) - pow(mrho,2) + s + t))*
                    (eta1*(8*pow(mpion,6) + pow(s,3) + 2*pow(mrho,4)*(s - t) + 5*pow(s,2)*t + s*pow(t,2) + pow(t,3) +
                         2*pow(mpion,2)*(2*pow(mrho,2) - s - t)*(s + t) - pow(mrho,2)*(3*s - t)*(s + t) - 2*pow(mpion,4)*(2*pow(mrho,2) + 3*s + t))\
                       + eta2*(s - t)*(4*pow(mpion,4) + t*(4*pow(mrho,2) - s + t) - pow(mpion,2)*(s + 3*t))))/
                  (pow(mrho,2)*(-pow(ma1,2) + t)*(-2*pow(mpion,2) + s + t)) -
                 (0.0625*pow(-2.*pow(mrho,2) + delta*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s - 1.*t),2)*
                    (8.*pow(mpion,6) + 4.*pow(mrho,6) + pow(s,3) + pow(mrho,4)*(-4.*s - 4.*t) + pow(mpion,4)*(-4.*pow(mrho,2) - 4.*s - 4.*t) +
                      3.*pow(s,2)*t + 3.*s*pow(t,2) + pow(t,3) + pow(mrho,2)*(-3.*pow(s,2) + 2.*s*t - 3.*pow(t,2)) +
                      pow(mpion,2)*(-8.*pow(mrho,4) - 2.*pow(s,2) - 4.*s*t - 2.*pow(t,2) + pow(mrho,2)*(4.*s + 4.*t))))/
                  (pow(mrho,6)*pow(2.*pow(mpion,2) - 1.*s - 1.*t,2)) +
                 (0.125*(-2 + delta)*(eta1 - eta2)*(eta2*(pow(mpion,2) + t)*
                       (pow(mpion,4) - pow(mpion,2)*(pow(mrho,2) - 2*t) + (pow(mrho,2) - 2*s - t)*t) +
                      eta1*(-4*pow(mpion,6) + (pow(mrho,2) - t)*(pow(mrho,2) - s - t)*t + pow(mpion,4)*(3*pow(mrho,2) + s + t) -
                         pow(mpion,2)*(pow(mrho,4) + pow(mrho,2)*(s - t) + 2*t*(-s + t)))))/((-pow(ma1,2) + t)*(-pow(mpion,2) + t)) +
                 (0.03125*pow(eta1 - eta2,2)*(-2*eta1*eta2*(pow(mpion,8) - pow(mpion,4)*(pow(mrho,4) + 2*(pow(mrho,2) + s)*t - 4*pow(t,2)) +
                         pow(t,2)*(-pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) + 2*s*t + pow(t,2)) +
                         2*pow(mpion,2)*t*(pow(mrho,4) + pow(mrho,2)*(s + t) - 2*t*(s + t))) +
                      pow(eta2,2)*(pow(mpion,8) - 2*pow(mpion,6)*pow(mrho,2) + pow(mpion,4)*(pow(mrho,4) + 4*pow(mrho,2)*t - 2*(s - 2*t)*t) +
                         pow(t,2)*(pow(mrho,4) + 2*pow(s,2) + 2*s*t + pow(t,2) + 2*pow(mrho,2)*(-s + t)) -
                         2*pow(mpion,2)*t*(pow(mrho,4) - pow(mrho,2)*(s - 2*t) + 2*t*(s + t))) +
                      pow(eta1,2)*(pow(mpion,8) - 2*pow(mpion,6)*pow(mrho,2) +
                         pow(mpion,4)*(3*pow(mrho,4) + 2*pow(mrho,2)*(s - 3*t) - 2*(s - 2*t)*t) +
                         t*(-pow(mrho,2) + t)*(2*pow(s,2) + 2*s*t + pow(t,2) - pow(mrho,2)*(2*s + t)) -
                         2*pow(mpion,2)*(-pow(mrho,2) + t)*(2*t*(s + t) - pow(mrho,2)*(s + 2*t)))))/pow(pow(ma1,2) - t,2) -
                 (0.5*(-2.*pow(mrho,2) + delta*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s - 1.*t))*
                    (delta*(-1.*pow(mpion,6) - 0.5*pow(mrho,6) - 0.1875*pow(s,3) +
                         pow(mpion,2)*(1.*pow(mrho,4) + pow(mrho,2)*(-0.625*s - 0.375*t) + s*(0.5*s + 0.5*t)) + pow(mrho,4)*(0.5*s + 0.5*t) +
                         pow(mpion,4)*(0.5*pow(mrho,2) + 0.25*s + 0.75*t) - 0.4375*pow(s,2)*t - 0.3125*s*pow(t,2) - 0.0625*pow(t,3) +
                         pow(mrho,2)*(0.4375*pow(s,2) - 0.25*s*t + 0.3125*pow(t,2))) +
                      pow(mrho,2)*(-0.125*pow(s,2) + C4*pow(mrho,4)*(1.*s - 1.*t) + 0.125*pow(t,2) +
                         pow(mpion,2)*((0.25 - 1.*C4*pow(mrho,2))*s + (-0.25 + 1.*C4*pow(mrho,2))*t) +
                         pow(mrho,2)*(-0.5*s + 0.5*C4*pow(s,2) + t*(0.5 - 0.5*C4*t)))))/(pow(mrho,6)*(1.*pow(mpion,2) - 0.5*s - 0.5*t)) +
                 (pow(delta,2)*(-0.5*pow(mpion,6) - 0.0625*pow(mrho,6) + pow(mrho,4)*(-0.125*s - 0.125*t) +
                       pow(mpion,4)*(1.*pow(mrho,2) + 0.5*t) + s*(-0.125*pow(s,2) - 0.25*s*t - 0.125*pow(t,2)) +
                       pow(mpion,2)*(1.25*pow(mrho,4) + 0.375*pow(s,2) + pow(mrho,2)*(-1.125*s - 0.875*t) + 0.25*s*t - 0.125*pow(t,2)) +
                       pow(mrho,2)*(0.4375*pow(s,2) + 0.25*s*t + 0.3125*pow(t,2))) +
                    pow(mrho,6)*(0.75 + C4*(8.*C4*pow(mpion,4) + 2.*C4*pow(s,2) + pow(mpion,2)*(6. - 8.*C4*s - 8.*C4*t) + t*(-3. + 2.*C4*t) +
                          s*(-3. + 4.*C4*t))) + delta*pow(mrho,2)*(pow(mpion,4)*(-0.5 - 4.*C4*pow(mrho,2)) + s*(-0.25*s - 0.25*t) +
                       pow(mrho,4)*(-0.75 + 2.5*C4*s + 0.5*C4*t) + pow(mrho,2)*(-0.5*C4*pow(s,2) + s*(0.25 - 2.*C4*t) + t*(1.25 - 1.5*C4*t)) +
                       pow(mpion,2)*(-3.*C4*pow(mrho,4) + 0.75*s + 0.25*t + pow(mrho,2)*(-1.5 + 3.*C4*s + 5.*C4*t))))/pow(mrho,6) +
                 (2*((0.0625*(-2. + delta)*(-2.*pow(mrho,2) + delta*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s - 1.*t))*
                         (2.*pow(mpion,6) + 1.*pow(mrho,6) + pow(mpion,4)*(-3.*pow(mrho,2) - 2.*t) + pow(mrho,4)*(-1.5*s - 1.5*t) +
                           pow(mrho,2)*s*(0.5*s + 0.5*t) + pow(mpion,2)*
                            (-1.*pow(mrho,4) - 0.5*pow(s,2) + pow(mrho,2)*(2.5*s - 0.5*t) - 1.*s*t - 0.5*pow(t,2)) +
                           t*(0.5*pow(s,2) + 1.*s*t + 0.5*pow(t,2))))/((pow(mpion,2) - 1.*t)*(1.*pow(mpion,2) - 0.5*s - 0.5*t)) +
                      (0.0625*(-2 + delta)*(6*delta*pow(mpion,6) + delta*s*t*(s + t) + pow(mrho,6)*(-2 + 3*delta + 8*C4*t) -
                           pow(mpion,4)*((-2 + 9*delta)*pow(mrho,2) - 8*C4*pow(mrho,4) + delta*(s + 9*t)) -
                           2*pow(mrho,4)*(t*(-1 + 3*delta + 4*C4*t) + s*(-1 + 2*delta + 8*C4*t)) -
                           pow(mpion,2)*(8*C4*pow(mrho,6) + 2*pow(mrho,4)*(-2 + delta - 8*C4*t) + pow(mrho,2)*((2 - 7*delta)*s + (2 + 5*delta)*t) +
                              delta*(pow(s,2) - 3*pow(t,2))) + pow(mrho,2)*(2*s*t + delta*(pow(s,2) + 3*s*t + 3*pow(t,2)))))/(-pow(mpion,2) + t)))/
                  pow(mrho,4) + (0.0625*(eta1 - eta2)*(-(eta2*(-2*pow(mpion,4)*
                            (4*C4*pow(mrho,2)*(pow(mrho,2) + 4*t) - delta*(pow(mrho,2) - 2*s + 6*t)) +
                           pow(mpion,2)*(2*pow(mrho,4)*(-2 + delta + 8*C4*t) + delta*(pow(s,2) - 6*s*t - 11*pow(t,2)) +
                              pow(mrho,2)*(-((-2 + delta)*s) + (-10 + delta)*t + 32*C4*t*(s + t))) +
                           t*(-2*pow(mrho,4)*(-2 + delta + 4*C4*t) + delta*(3*pow(s,2) + 2*s*t + 3*pow(t,2)) +
                              pow(mrho,2)*((2 - 5*delta)*s + 3*(2 + delta)*t - 8*C4*pow(s + t,2))))) +
                      eta1*(8*delta*pow(mpion,6) + delta*(pow(s,3) + 7*pow(s,2)*t + 5*s*pow(t,2) + 3*pow(t,3)) -
                         2*pow(mrho,2)*((-1 + 2*delta)*pow(s,2) + 2*(-1 + 2*delta)*s*t + (-3 + 2*delta)*pow(t,2) + 4*C4*t*pow(s + t,2)) +
                         pow(mpion,4)*(24*C4*pow(mrho,4) + 6*delta*(-s + t) - 4*pow(mrho,2)*(-1 + 3*delta + 8*C4*t)) +
                         pow(mrho,4)*(t*(-6 + delta + 8*C4*t) + s*(-2 + 3*delta + 16*C4*t)) -
                         2*pow(mpion,2)*(delta*(s + t)*(s + 5*t) - pow(mrho,2)*(-2*s + 5*delta*s - 6*t + 9*delta*t + 16*C4*t*(s + t)) +
                            2*pow(mrho,4)*(-2 + delta + 4*C4*(s + 2*t))))))/(pow(mrho,2)*(-pow(ma1,2) + t))))/(16.*Pi*s*(-4*pow(mpion,2) + s));

      /* omega: diff_xsection = (0.0024867959858108648*pow(Const,2)*pow(g_POR,4)*(pow(mpion,8) - 2*pow(mpion,6)*pow(mrho,2) +
        pow(mpion,4)*(pow(mrho,4) + 4*pow(s,2) - 2*s*t) +
        pow(s,2)*(pow(mrho,4) + pow(s,2) + 2*s*t + 2*pow(t,2) - 2*pow(mrho,2)*(s + t)) -
        2*pow(mpion,2)*s*(pow(mrho,4) + 2*s*(s + t) - pow(mrho,2)*(2*s + t))))/(pow(pow(momega,2) - s,2)*(pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s))); */
      break;
    case ReactionType::pi0_rho:
      diff_xsection = (pow(Const,2)*pow(ghat,4)*((-0.25*pow(-2 + delta,2)*pow(mpion,2)*
                (pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s)))/(pow(mrho,2)*pow(pow(mpion,2) - s,2)) -
             (0.0625*(eta1 - eta2)*(-pow(ma1,2) + s)*(2*pow(mrho,2) + delta*(-2*pow(mpion,2) - pow(mrho,2) + s + t))*
                (-(eta2*(s - t)*(4*pow(mpion,4) + s*(4*pow(mrho,2) + s - t) - pow(mpion,2)*(3*s + t))) +
                  eta1*(8*pow(mpion,6) + pow(s,3) + pow(s,2)*t + 5*s*pow(t,2) + pow(t,3) + 2*pow(mrho,4)*(-s + t) +
                     pow(mrho,2)*(s - 3*t)*(s + t) + 2*pow(mpion,2)*(2*pow(mrho,2) - s - t)*(s + t) - 2*pow(mpion,4)*(2*pow(mrho,2) + s + 3*t))))
               /(pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(-2*pow(mpion,2) + s + t)) -
             (0.0625*pow(-2.*pow(mrho,2) + delta*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s - 1.*t),2)*
                (8.*pow(mpion,6) + 4.*pow(mrho,6) + pow(s,3) + pow(mrho,4)*(-4.*s - 4.*t) + pow(mpion,4)*(-4.*pow(mrho,2) - 4.*s - 4.*t) +
                  3.*pow(s,2)*t + 3.*s*pow(t,2) + pow(t,3) + pow(mrho,2)*(-3.*pow(s,2) + 2.*s*t - 3.*pow(t,2)) +
                  pow(mpion,2)*(-8.*pow(mrho,4) - 2.*pow(s,2) - 4.*s*t - 2.*pow(t,2) + pow(mrho,2)*(4.*s + 4.*t))))/
              (pow(mrho,6)*pow(2.*pow(mpion,2) - 1.*s - 1.*t,2)) +
             (0.125*(-2 + delta)*(eta1 - eta2)*(-pow(ma1,2) + s)*
                (-(eta2*(pow(mpion,2) + s)*(-pow(mpion,4) + pow(mpion,2)*(pow(mrho,2) - 2*s) + s*(-pow(mrho,2) + s + 2*t))) +
                  eta1*(-4*pow(mpion,6) + s*(-pow(mrho,2) + s)*(-pow(mrho,2) + s + t) + pow(mpion,4)*(3*pow(mrho,2) + s + t) -
                     pow(mpion,2)*(pow(mrho,4) + 2*s*(s - t) + pow(mrho,2)*(-s + t)))))/
              ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(-pow(mpion,2) + s)) +
             (0.03125*pow(eta1 - eta2,2)*(pow(eta2,2)*(pow(mpion,8) - 2*pow(mpion,6)*pow(mrho,2) +
                     pow(mpion,4)*(pow(pow(mrho,2) + 2*s,2) - 2*s*t) +
                     pow(s,2)*(pow(pow(mrho,2) + s,2) + 2*(-pow(mrho,2) + s)*t + 2*pow(t,2)) -
                     2*pow(mpion,2)*s*(pow(mrho,4) + pow(mrho,2)*(2*s - t) + 2*s*(s + t))) -
                  2*eta1*eta2*(pow(mpion,8) - pow(mpion,4)*(pow(mrho,4) + 2*pow(mrho,2)*s + 2*s*(-2*s + t)) +
                     2*pow(mpion,2)*s*(pow(mrho,4) + pow(mrho,2)*(s + t) - 2*s*(s + t)) +
                     pow(s,2)*(-pow(mrho,4) + pow(s,2) - 2*pow(mrho,2)*t + 2*t*(s + t))) +
                  pow(eta1,2)*(pow(mpion,8) - 2*pow(mpion,6)*pow(mrho,2) +
                     pow(mpion,4)*(3*pow(mrho,4) + 2*s*(2*s - t) + 2*pow(mrho,2)*(-3*s + t)) -
                     2*pow(mpion,2)*(-pow(mrho,2) + s)*(2*s*(s + t) - pow(mrho,2)*(2*s + t)) +
                     s*(-pow(mrho,2) + s)*(pow(s,2) + 2*s*t + 2*pow(t,2) - pow(mrho,2)*(s + 2*t)))))/
              (pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2)) +
             (0.5*(-2.*pow(mrho,2) + delta*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s - 1.*t))*
                (delta*(1.*pow(mpion,6) + 0.5*pow(mrho,6) + 0.0625*pow(s,3) + pow(mrho,4)*(-0.5*s - 0.5*t) +
                     pow(mpion,4)*(-0.5*pow(mrho,2) - 0.75*s - 0.25*t) + 0.3125*pow(s,2)*t + 0.4375*s*pow(t,2) + 0.1875*pow(t,3) +
                     pow(mpion,2)*(-1.*pow(mrho,4) + pow(mrho,2)*(0.375*s + 0.625*t) + (-0.5*s - 0.5*t)*t) +
                     pow(mrho,2)*(-0.3125*pow(s,2) + 0.25*s*t - 0.4375*pow(t,2))) +
                  pow(mrho,2)*(-0.125*pow(s,2) + C4*pow(mrho,4)*(1.*s - 1.*t) + 0.125*pow(t,2) +
                     pow(mpion,2)*((0.25 - 1.*C4*pow(mrho,2))*s + (-0.25 + 1.*C4*pow(mrho,2))*t) +
                     pow(mrho,2)*(-0.5*s + 0.5*C4*pow(s,2) + t*(0.5 - 0.5*C4*t)))))/(pow(mrho,6)*(1.*pow(mpion,2) - 0.5*s - 0.5*t)) +
             (pow(delta,2)*(-0.5*pow(mpion,6) - 0.0625*pow(mrho,6) + pow(mpion,4)*(1.*pow(mrho,2) + 0.5*s) +
                   pow(mrho,4)*(-0.125*s - 0.125*t) + t*(-0.125*pow(s,2) - 0.25*s*t - 0.125*pow(t,2)) +
                   pow(mpion,2)*(1.25*pow(mrho,4) - 0.125*pow(s,2) + pow(mrho,2)*(-0.875*s - 1.125*t) + 0.25*s*t + 0.375*pow(t,2)) +
                   pow(mrho,2)*(0.3125*pow(s,2) + 0.25*s*t + 0.4375*pow(t,2))) +
                delta*pow(mrho,2)*(pow(mpion,4)*(-0.5 - 4.*C4*pow(mrho,2)) + (-0.25*s - 0.25*t)*t + pow(mrho,4)*(-0.75 + 0.5*C4*s + 2.5*C4*t) +
                   pow(mrho,2)*(-1.5*C4*pow(s,2) + s*(1.25 - 2.*C4*t) + t*(0.25 - 0.5*C4*t)) +
                   pow(mpion,2)*(-3.*C4*pow(mrho,4) + 0.25*s + 0.75*t + pow(mrho,2)*(-1.5 + 5.*C4*s + 3.*C4*t))) +
                pow(mrho,6)*(0.75 + C4*(8.*C4*pow(mpion,4) + 2.*C4*pow(s,2) + pow(mpion,2)*(6. - 8.*C4*s - 8.*C4*t) + t*(-3. + 2.*C4*t) +
                      s*(-3. + 4.*C4*t))))/pow(mrho,6) + (0.0625*(eta1 - eta2)*(-pow(ma1,2) + s)*
                (-(eta2*(2*pow(mpion,4)*(-4*C4*pow(mrho,2)*(pow(mrho,2) + 4*s) + delta*(pow(mrho,2) + 6*s - 2*t)) +
                       pow(mpion,2)*(2*pow(mrho,4)*(-2 + delta + 8*C4*s) + delta*(-11*pow(s,2) - 6*s*t + pow(t,2)) +
                          pow(mrho,2)*((-10 + delta)*s - (-2 + delta)*t + 32*C4*s*(s + t))) +
                       s*(-2*pow(mrho,4)*(-2 + delta + 4*C4*s) + delta*(3*pow(s,2) + 2*s*t + 3*pow(t,2)) +
                          pow(mrho,2)*(3*(2 + delta)*s + (2 - 5*delta)*t - 8*C4*pow(s + t,2))))) +
                  eta1*(8*delta*pow(mpion,6) + 2*pow(mpion,4)*(12*C4*pow(mrho,4) - 2*pow(mrho,2)*(-1 + 3*delta + 8*C4*s) + 3*delta*(s - t)) +
                     delta*(3*pow(s,3) + 5*pow(s,2)*t + 7*s*pow(t,2) + pow(t,3)) -
                     2*pow(mrho,2)*((-3 + 2*delta)*pow(s,2) + 2*(-1 + 2*delta)*s*t + (-1 + 2*delta)*pow(t,2) + 4*C4*s*pow(s + t,2)) +
                     pow(mrho,4)*((-6 + delta)*s + (-2 + 3*delta)*t + 8*C4*s*(s + 2*t)) -
                     2*pow(mpion,2)*(delta*(s + t)*(5*s + t) - pow(mrho,2)*(-6*s + 9*delta*s - 2*t + 5*delta*t + 16*C4*s*(s + t)) +
                        2*pow(mrho,4)*(-2 + delta + 4*C4*(2*s + t))))))/(pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))) +
             (2*((0.0625*(-2. + delta)*(-2.*pow(mrho,2) + delta*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s - 1.*t))*
                     (2.*pow(mpion,6) + 1.*pow(mrho,6) + pow(mpion,4)*(-3.*pow(mrho,2) - 2.*s) + pow(mrho,4)*(-1.5*s - 1.5*t) +
                       pow(mrho,2)*(0.5*s + 0.5*t)*t + s*(0.5*pow(s,2) + 1.*s*t + 0.5*pow(t,2)) +
                       pow(mpion,2)*(-1.*pow(mrho,4) - 0.5*pow(s,2) - 1.*s*t - 0.5*pow(t,2) + pow(mrho,2)*(-0.5*s + 2.5*t))))/
                   ((pow(mpion,2) - 1.*s)*(1.*pow(mpion,2) - 0.5*s - 0.5*t)) +
                  (0.0625*(-2 + delta)*(delta*(6*pow(mpion,6) - pow(mpion,4)*(9*(pow(mrho,2) + s) + t) -
                          pow(mpion,2)*(2*pow(mrho,4) - 3*pow(s,2) + pow(mrho,2)*(5*s - 7*t) + pow(t,2)) +
                          (pow(mrho,2) - s - t)*(3*pow(mrho,4) - s*t - pow(mrho,2)*(3*s + t))) +
                       2*pow(mrho,2)*(pow(mpion,4)*(1 + 4*C4*pow(mrho,2)) + pow(mrho,4)*(-1 + 4*C4*s) + s*t -
                          pow(mpion,2)*(4*C4*pow(mrho,4) + s - 2*pow(mrho,2)*(1 + 4*C4*s) + t) + pow(mrho,2)*(t + s*(1 - 4*C4*(s + 2*t))))))/
                   (-pow(mpion,2) + s)))/pow(mrho,4)))/(16.*Pi*s*(-4*pow(mpion,2) + s));

      /* omega: diff_xsection = (0.0024867959858108648*pow(Const,2)*pow(g_POR,4)*(pow(mpion,8) - 2*pow(mpion,6)*pow(mrho,2) +
       pow(mpion,4)*(pow(mrho,4) - 2*(s - 2*t)*t) + pow(t,2)*(pow(mrho,4) + 2*pow(s,2) + 2*s*t + pow(t,2) - 2*pow(mrho,2)*(s + t)) -
       2*pow(mpion,2)*t*(pow(mrho,4) + 2*t*(s + t) - pow(mrho,2)*(s + 2*t))))/
       ((pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s))*pow(pow(momega,2) - t,2)); */
      break; 
    /*case ReactionType::pi_eta:
      diff_xsection = to_be_determined;
      break;*/
    case ReactionType::pi0_rho0:

      diff_xsection = 1/3.0*(pow(Const,2)*pow(g_POR,4)*(pow(m_omega,4)*pow(s,4) + 4*pow(m_omega,4)*pow(s,3)*t - 4*pow(m_omega,2)*pow(s,4)*t + 10*pow(m_omega,4)*pow(s,2)*pow(t,2) -
             16*pow(m_omega,2)*pow(s,3)*pow(t,2) + 5*pow(s,4)*pow(t,2) + 4*pow(m_omega,4)*s*pow(t,3) - 16*pow(m_omega,2)*pow(s,2)*pow(t,3) +
             10*pow(s,3)*pow(t,3) + pow(m_omega,4)*pow(t,4) - 4*pow(m_omega,2)*s*pow(t,4) + 5*pow(s,2)*pow(t,4) + pow(m_pi,8)*pow(-2*pow(m_omega,2) + s + t,2) -
             2*pow(m_pi,6)*pow(m_rho,2)*(2*pow(m_omega,4) + pow(s,2) + pow(t,2) - 2*pow(m_omega,2)*(s + t)) +
             pow(m_rho,4)*(2*pow(s,2)*pow(t,2) - 2*pow(m_omega,2)*s*t*(s + t) + pow(m_omega,4)*(pow(s,2) + pow(t,2))) -
             2*pow(m_rho,2)*(3*pow(s,2)*pow(t,2)*(s + t) - 3*pow(m_omega,2)*s*t*pow(s + t,2) +
                pow(m_omega,4)*(pow(s,3) + 2*pow(s,2)*t + 2*s*pow(t,2) + pow(t,3))) +
             pow(m_pi,4)*(-2*pow(m_rho,2)*(pow(m_omega,2) - s)*(pow(m_omega,2) - t)*(s + t) - 8*pow(m_omega,2)*s*t*(s + t) + 4*pow(m_omega,4)*(pow(s,2) + pow(t,2)) -
                2*s*t*(pow(s,2) - 6*s*t + pow(t,2)) + pow(m_rho,4)*(2*pow(m_omega,4) + pow(s,2) + pow(t,2) - 2*pow(m_omega,2)*(s + t))) -
             2*pow(m_pi,2)*(2*(s + t)*pow(-2*s*t + pow(m_omega,2)*(s + t),2) + pow(m_rho,4)*(-4*pow(m_omega,2)*s*t + pow(m_omega,4)*(s + t) + s*t*(s + t)) -
                pow(m_rho,2)*(-10*pow(m_omega,2)*s*t*(s + t) + 2*pow(m_omega,4)*(pow(s,2) + 3*s*t + pow(t,2)) + s*t*(pow(s,2) + 8*s*t + pow(t,2))))))/
         (128.*Pi*pow(pow(m_omega,2) - s,2)*(pow(pow(m_pi,2) - pow(m_rho,2),2) - 2*(pow(m_pi,2) + pow(m_rho,2))*s + pow(s,2))*pow(pow(m_omega,2) - t,2));
      break;
    case ReactionType::no_reaction:
      // never reached
      break;
  }
  return diff_xsection*to_mb;
}

}  // namespace Smash
