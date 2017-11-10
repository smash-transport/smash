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
#include "include/photoncrosssections.h"
#include "include/pow.h"
#include "include/random.h"
#include "include/tabulation.h"

using std::atan;
using std::pow;
using std::sqrt;

namespace Smash {

ScatterActionPhoton::ReactionType ScatterActionPhoton::photon_reaction_type(
    const ParticleList &in) {
  // ToDo: is this check here necessary?
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
      return ReactionType::pi_z_pi_p_rho_p;

    case (pack(pdg::pi_m, pdg::pi_z)):
    case (pack(pdg::pi_z, pdg::pi_m)):
      return ReactionType::pi_z_pi_m_rho_m;

    case (pack(pdg::pi_p, pdg::rho_z)):
      return ReactionType::pi_p_rho_z_pi_p;

    case (pack(pdg::pi_m, pdg::rho_z)):
      return ReactionType::pi_m_rho_z_pi_m;

    case (pack(pdg::pi_m, pdg::rho_p)):
      return ReactionType::pi_m_rho_p_pi_z;

    case (pack(pdg::pi_p, pdg::rho_m)):
      return ReactionType::pi_p_rho_m_pi_z;

    case (pack(pdg::pi_z, pdg::rho_p)):
      return ReactionType::pi_z_rho_p_pi_p;

    case (pack(pdg::pi_z, pdg::rho_m)):
      return ReactionType::pi_z_rho_m_pi_m;

    case (pack(pdg::pi_p, pdg::pi_m)):
    case (pack(pdg::pi_m, pdg::pi_p)):
      return ReactionType::pi_p_pi_m_rho_z;

    case (pack(pdg::pi_z, pdg::rho_z)):
      return ReactionType::pi_z_rho_z_pi_z;

    default:
      return ReactionType::no_reaction;
  }
}

ParticleTypePtr ScatterActionPhoton::outgoing_hadron(const ParticleList &in) {
  // ToDo: not sure how to handle the returning of pointers. i have to think
  // about this
  const static ParticleTypePtr rho_z_particle_ptr =
      &ParticleType::find(pdg::rho_z);
  const static ParticleTypePtr rho_p_particle_ptr =
      &ParticleType::find(pdg::rho_p);
  const static ParticleTypePtr rho_m_particle_ptr =
      &ParticleType::find(pdg::rho_m);
  const static ParticleTypePtr pi_z_particle_ptr =
      &ParticleType::find(pdg::pi_z);
  const static ParticleTypePtr pi_p_particle_ptr =
      &ParticleType::find(pdg::pi_p);
  const static ParticleTypePtr pi_m_particle_ptr =
      &ParticleType::find(pdg::pi_m);
  // const static ParticleTypePtr photon_particle =
  // &ParticleType::find(pdg::photon);

  auto reac = photon_reaction_type(in);

  switch (reac) {
    case ReactionType::pi_z_pi_p_rho_p:
      return rho_p_particle_ptr;
      break;
    case ReactionType::pi_z_pi_m_rho_m:
      return rho_m_particle_ptr;
      break;
    case ReactionType::pi_p_pi_m_rho_z:
      return rho_z_particle_ptr;
      break;

    case ReactionType::pi_p_rho_z_pi_p:
    case ReactionType::pi_z_rho_p_pi_p:
      return pi_p_particle_ptr;

    case ReactionType::pi_m_rho_z_pi_m:
    case ReactionType::pi_z_rho_m_pi_m:
      return pi_m_particle_ptr;

    case ReactionType::pi_m_rho_p_pi_z:
    case ReactionType::pi_p_rho_m_pi_z:
    case ReactionType::pi_z_rho_z_pi_z:
      return pi_z_particle_ptr;
      break;
    default:
      // ToDo: Think of something better, now returns particletypeptr with
      // ivnalid index
      ParticleTypePtr p;
      return p;
  }
}

bool ScatterActionPhoton::is_kinematically_possible(const double s,
                                                    const ParticleList &in) {
  auto hadron = outgoing_hadron(in);
  if (hadron && hadron->mass() > s) {
    return true;
  } else {
    return false;
  }
}

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
  double diff_xsection_max = 0.0;
  const double stepsize = (t2 - t1) / 100.0;
  for (double t = t1; t < t2; t += stepsize) {
    diff_xsection_max =
        std::max(diff_cross_section(t, t2, t1), diff_xsection_max);
  }

  double t = 0.0;
  int iteration_number = 0;
  // Bug: diff_xsection_max set to 0 but then used in comparison
  do {
    t = Random::uniform(t1, t2);
    iteration_number++;
  } while (diff_cross_section(t, t2, t1) <
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
    weight_ = diff_cross_section(t, t2, t1) * (t2 - t1) /
              (number_of_fractional_photons_ * cross_section());
  } else {
    weight_ = proc->weight() / cross_section();
  }

  /* Set positions & boost to computational frame. */
  for (ParticleData &new_particle : outgoing_particles_) {
    new_particle.set_4position(middle_point);
    new_particle.boost_momentum(-beta_cm());
  }

  /* Inlcusion of form factors (FF):
  The usual procedure would be the multplication of the photon cross section
  by the corresponding form factor. This form factor is however energy
  dependent, such that the energy of the generated photon in the computational
  frame is necessary to determine FF. Yet this is not directly accessible in
  ScatterActionPhoton::photon_cross_section().
  The alternative solution is to multiply the weighting factor (proportional to
  cross section) by FF. This is equivalent to multiplying the cross section
  directly.

  The modification is as follows:
  weight_FF = weight_noFF * FF^4
  The actual value of the form factor is determined in
  ScatterActionPhoton::form_factor(E_photon) */

  double E_Photon = outgoing_particles_[1].momentum()[0];

  weight_ *= pow(form_factor(E_Photon), 4);

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

/*
ScatterActionPhoton::ReactionType ScatterActionPhoton::photon_reaction_type(
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
    case (pack(pdg::pi_p, pdg::pi_z)
    case (pack(pdg::pi_z, pdg::pi_p)):
      return ReactionType::pi_z_pi_p_rho_p;

    case (pack(pdg::pi_m, pdg::pi_z)):
    case (pack(pdg::pi_z, pdg::pi_m)):
      return ReactionType::pi_z_pi_m_rho_m;

    case (pack(pdg::pi_p, pdg::rho_z)):
      return ReactionType::pi_p_rho_z_pi_p;

    case (pack(pdg::pi_m, pdg::rho_z)):
      return ReactionType::pi_m_rho_z_pi_m;

    case (pack(pdg::pi_m, pdg::rho_p)):
      return ReactionType::pi_m_rho_p_pi_z;

    case (pack(pdg::pi_p, pdg::rho_m)):
      return ReactionType::pi_p_rho_m_pi_z;

    case (pack(pdg::pi_z, pdg::rho_p)):
      return ReactionType::pi_z_rho_p_pi_p;

    case (pack(pdg::pi_z, pdg::rho_m)):
      return ReactionType::pi_z_rho_m_pi_m;

      case(pack(pdg::pi_p, pdg::pi_m)):
    case(pack(pdg::pi_m, pdg::pi_p)):
      return ReactionType::pi_p_pi_m_rho_z;

    case(pack(pdg::pi_z, pdg::rho_z)):
      return ReactionType::pi_z_rho_z_pi_z;

      default:
      return ReactionType::no_reaction;
  }
}
*/

CollisionBranchList ScatterActionPhoton::photon_cross_sections() {
  CollisionBranchList process_list;
  PhotonCrossSection<ComputationMethod::Lookup> xs_object;
  
  reac = photon_reaction_type(Action::incoming_particles());

  auto hadron_out = outgoing_hadron(incoming_particles_);
  static ParticleTypePtr photon_particle = &ParticleType::find(pdg::photon);
  


  ParticleData part_a = incoming_particles_[0];
  ParticleData part_b = incoming_particles_[1];
  const double &m1 = part_a.effective_mass();
  const double &m2 = part_b.effective_mass();
  const double &m3 = hadron_out->mass();

  const double s = mandelstam_s();
  const double sqrts = sqrt_s();

  ParticleTypePtr photon_out = photon_particle;
  double xsection = 0.0;
  switch (reac) {
    case ReactionType::pi_p_pi_m_rho_z:
        xsection = xs_object.xs_pi_pi_rho0(s);
        process_list.push_back(make_unique<CollisionBranch>(
            *hadron_out, *photon_out, xsection, ProcessType::TwoToTwo));
      break;

    case ReactionType::pi_z_pi_m_rho_m:
    case ReactionType::pi_z_pi_p_rho_p:
        xsection = xs_object.xs_pi_pi0_rho(s);
        process_list.push_back(make_unique<CollisionBranch>(
            *hadron_out, *photon_out, xsection, ProcessType::TwoToTwo));
        break;

    case ReactionType::pi_m_rho_z_pi_m:
    case ReactionType::pi_p_rho_z_pi_p:
        xsection = xs_object.xs_pi_rho0_pi(s);
        process_list.push_back(make_unique<CollisionBranch>(
            *hadron_out, *photon_out, xsection, ProcessType::TwoToTwo));
      break;

    case ReactionType::pi_m_rho_p_pi_z:
    case ReactionType::pi_p_rho_m_pi_z:

        xsection = xs_object.xs_pi_rho_pi0(s);
        process_list.push_back(make_unique<CollisionBranch>(
            *hadron_out, *photon_out, xsection, ProcessType::TwoToTwo));
      break;

    case ReactionType::pi_z_rho_m_pi_m:
    case ReactionType::pi_z_rho_p_pi_p:
        xsection = xs_object.xs_pi0_rho_pi(s);
        process_list.push_back(make_unique<CollisionBranch>(
            *hadron_out, *photon_out, xsection, ProcessType::TwoToTwo));
      break;

    case ReactionType::pi_z_rho_z_pi_z:
        xsection = xs_object.xs_pi0_rho0_pi0(s);
        process_list.push_back(make_unique<CollisionBranch>(
            *hadron_out, *photon_out, xsection, ProcessType::TwoToTwo));
      break;
  
    case ReactionType::no_reaction:
      // never reached
      break;
  }

  return process_list;
}

double ScatterActionPhoton::diff_cross_section(double t, double t2,
                                               double t1) const {
  const double to_mb = 0.3894;
  static const float m_rho = ParticleType::find(pdg::rho_z).mass();
  static const float m_pi = ParticleType::find(pdg::pi_z).mass();
  double s = mandelstam_s();
  double diff_xsection = 0.0;

  PhotonCrossSection<ComputationMethod::Analytic> xs_object;

  switch (reac) {
    case ReactionType::pi_p_pi_m_rho_z:
      if (outgoing_particles_[0].type().pdgcode().is_rho()) {
        diff_xsection = xs_object.xs_diff_pi_pi_rho0(s, t);

      } else if (outgoing_particles_[0].type().pdgcode() == pdg::photon) {
        diff_xsection = 0.0000000000001 / to_mb / (t2 - t1);
      }
      break;

    case ReactionType::pi_z_pi_m_rho_m:
    case ReactionType::pi_z_pi_p_rho_p:

      if (outgoing_particles_[0].type().pdgcode().is_rho()) {
        diff_xsection = xs_object.xs_diff_pi_pi0_rho(s, t);

      } else if (outgoing_particles_[0].type().pdgcode().is_pion()) {
        diff_xsection = 0.0000000000001 / to_mb / (t2 - t1);
      }
      break;

    case ReactionType::pi_m_rho_z_pi_m:
    case ReactionType::pi_p_rho_z_pi_p:
      diff_xsection = xs_object.xs_diff_pi_rho0_pi(s, t);

      break;

    case ReactionType::pi_m_rho_p_pi_z:
    case ReactionType::pi_p_rho_m_pi_z:
      // omega:
      diff_xsection = xs_object.xs_diff_pi_rho_pi0(s, t);

      break;

    case ReactionType::pi_z_rho_m_pi_m:
    case ReactionType::pi_z_rho_p_pi_p:
      // omega:
      diff_xsection = xs_object.xs_diff_pi0_rho_pi(s, t);

      break;
    case ReactionType::pi_z_rho_z_pi_z:

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
      takes the shape: FF = (2*Lambda^2/(2*Lambda^2 - t))^2 with
      Lambda = 1.0 GeV. t depends on the lightest possible exchange particle in
      the different channels. This could either be a pion or an omega meson. For
      the computation the parametrizations given in \ref! are used. */
      // TODO (schaefer): Include reference for FF!

    case ReactionType::pi_p_pi_m_rho_z:
    case ReactionType::pi_z_pi_p_rho_p:
    case ReactionType::pi_z_pi_m_rho_m:
    case ReactionType::pi_p_rho_z_pi_p:
    case ReactionType::pi_m_rho_z_pi_m:
      // case ReactionType::pi_rho:
      // case ReactionType::pi0_rho:
      t_ff = 34.5096 * pow(E_photon, 0.737) - 67.557 * pow(E_photon, 0.7584) +
             32.858 * pow(E_photon, 0.7806);
      break;
    // lightest exchange particle: omega
    case ReactionType::pi_m_rho_p_pi_z:
    case ReactionType::pi_p_rho_m_pi_z:  // for omeaga

    case ReactionType::pi_z_rho_m_pi_m:
    case ReactionType::pi_z_rho_p_pi_p:  // for omega

    case ReactionType::pi_z_rho_z_pi_z:
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
