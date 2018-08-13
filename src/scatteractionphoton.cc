/*
 *
 *    Copyright (c) 2016-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/scatteractionphoton.h"

#include <algorithm>

#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/crosssectionsphoton.h"
#include "smash/cxx14compat.h"
#include "smash/forwarddeclarations.h"
#include "smash/kinematics.h"
#include "smash/outputinterface.h"
#include "smash/particletype.h"
#include "smash/pdgcode.h"
#include "smash/pow.h"
#include "smash/random.h"
#include "smash/tabulation.h"

namespace smash {

ScatterActionPhoton::ScatterActionPhoton(const ParticleList &in,
                    const double time,
                    const int n_frac_photons,
                    const double hadronic_cross_section_input)
    : ScatterAction(in[0], in[1], time),
      reac_(photon_reaction_type(in)),
      number_of_fractional_photons_(n_frac_photons),
      hadron_out_t_(outgoing_hadron_type(in)),
      hadron_out_mass_(sample_out_hadron_mass(hadron_out_t_)),
      hadronic_cross_section_(hadronic_cross_section_input) {}

ScatterActionPhoton::ReactionType ScatterActionPhoton::photon_reaction_type(
    const ParticleList &in) {
  if (in.size() != 2) {
    return ReactionType::no_reaction;
  }

  PdgCode a = in[0].pdgcode();
  PdgCode b = in[1].pdgcode();

  // swap so that pion is first and there are less cases to be listed
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

void ScatterActionPhoton::perform_photons(const OutputsList &outputs) {
  for (int i = 0; i < number_of_fractional_photons_; i++) {
    generate_final_state();
    for (const auto &output : outputs) {
      if (output->is_photon_output()) {
        // we do not care about the local density
        output->at_interaction(*this, 0.0);
      }
    }
  }
}

ParticleTypePtr ScatterActionPhoton::outgoing_hadron_type(
    const ReactionType reaction) {
  static const ParticleTypePtr rho_z_particle_ptr =
      &ParticleType::find(pdg::rho_z);
  static const ParticleTypePtr rho_p_particle_ptr =
      &ParticleType::find(pdg::rho_p);
  static const ParticleTypePtr rho_m_particle_ptr =
      &ParticleType::find(pdg::rho_m);
  static const ParticleTypePtr pi_z_particle_ptr =
      &ParticleType::find(pdg::pi_z);
  static const ParticleTypePtr pi_p_particle_ptr =
      &ParticleType::find(pdg::pi_p);
  static const ParticleTypePtr pi_m_particle_ptr =
      &ParticleType::find(pdg::pi_m);

  switch (reaction) {
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
      // default constructor constructs p with invalid index
      ParticleTypePtr p{};
      return p;
  }
}

ParticleTypePtr ScatterActionPhoton::outgoing_hadron_type(
    const ParticleList &in) {
  auto reac = photon_reaction_type(in);
  return outgoing_hadron_type(reac);
}

bool ScatterActionPhoton::is_kinematically_possible(const double s_sqrt,
                                                    const ParticleList &in) {
  auto reac = photon_reaction_type(in);
  auto hadron = outgoing_hadron_type(in);

  if (reac == ReactionType::no_reaction)
    return false;

  if (default_mediator_ == MediatorType::PION &&
      reac == ReactionType::pi_z_rho_z_pi_z) {
    return false;
  }

  // C15 has only s-channel. Make sure that CM-energy is high
  // enough to produce mediating omega meson
  if ((reac == ReactionType::pi_m_rho_p_pi_z ||
       reac == ReactionType::pi_p_rho_m_pi_z) &&
      default_mediator_ == MediatorType::OMEGA) {
    if (s_sqrt < omega_mass) {
      return false;
    }
  }

  // for all other processes: if cm-energy is not high enough to produce final
  // state particle reject the collision.
  if (hadron->is_stable() && s_sqrt < hadron->mass()) {
    return false;
  } else {
    return true;
  }
}

void ScatterActionPhoton::generate_final_state() {
  // we have only one reaction per incoming particle pair
  if (collision_processes_photons_.size() != 1) {
    const auto &log = logger<LogArea::ScatterAction>();
    log.fatal() << "Problem in ScatterActionPhoton::generate_final_state().\n";
    throw std::runtime_error("");
  }
  auto *proc = collision_processes_photons_[0].get();

  outgoing_particles_ = proc->particle_list();

  FourVector middle_point = get_interaction_point();

  // 2->2 inelastic scattering
  // Sample the particle momenta in CM system
  const double m1 = incoming_particles_[0].effective_mass();
  const double m2 = incoming_particles_[1].effective_mass();

  const double &m_out = hadron_out_mass_;

  const double s = mandelstam_s();
  const double sqrts = sqrt_s();
  std::array<double, 2> mandelstam_t = get_t_range(sqrts, m1, m2, m_out, 0.0);
  const double t1 = mandelstam_t[1];
  const double t2 = mandelstam_t[0];
  const double pcm_in = cm_momentum();
  const double pcm_out = pCM(sqrts, m_out, 0.0);

  const double t = random::uniform(t1, t2);

  double costheta = (t - pow_int(m2, 2) +
                     0.5 * (s + pow_int(m2, 2) - pow_int(m1, 2)) *
                         (s - pow_int(m_out, 2)) / s) /
                    (pcm_in * (s - pow_int(m_out, 2)) / sqrts);

  Angles phitheta(random::uniform(0.0, twopi), costheta);
  outgoing_particles_[0].set_4momentum(hadron_out_mass_,
                                       phitheta.threevec() * pcm_out);
  outgoing_particles_[1].set_4momentum(0.0, -phitheta.threevec() * pcm_out);

  // Set positions & boost to computational frame.
  for (ParticleData &new_particle : outgoing_particles_) {
    new_particle.set_4position(middle_point);
    new_particle.boost_momentum(-beta_cm());
  }

  double E_Photon = outgoing_particles_[1].momentum()[0];

  // if rho in final state take already sampled mass (same as m_out). If rho is
  // incoming take the mass of the incoming particle
  const double m_rho = rho_mass();

  // compute the differential cross section with form factor included
  const double diff_xs = diff_cross_section_w_ff(t, m_rho, E_Photon);

  // Weighing of the fractional photons
  if (number_of_fractional_photons_ > 1) {
    weight_ = diff_xs * (t2 - t1) /
              (number_of_fractional_photons_ * hadronic_cross_section());
  } else {
    weight_ = proc->weight() / hadronic_cross_section();
  }
  // Photons are not really part of the normal processes, so we have to set a
  // constant arbitrary number.
  const auto id_process = ID_PROCESS_PHOTON;
  Action::check_conservation(id_process);
}

void ScatterActionPhoton::add_dummy_hadronic_process(
    double reaction_cross_section) {
  CollisionBranchPtr dummy_process = make_unique<CollisionBranch>(
      incoming_particles_[0].type(), incoming_particles_[1].type(),
      reaction_cross_section, ProcessType::TwoToTwo);
  add_collision(std::move(dummy_process));
}

double ScatterActionPhoton::sample_out_hadron_mass(
    const ParticleTypePtr out_t) {
  double mass = out_t->mass();
  const double cms_energy = sqrt_s();
  if (cms_energy <= out_t->min_mass_kinematic()) {
    throw InvalidResonanceFormation(
        "Problem in ScatterActionPhoton::sample_hadron_mass");
  }

  if (!out_t->is_stable()) {
    mass = out_t->sample_resonance_mass(0, cms_energy);
  }

  return mass;
}

double ScatterActionPhoton::rho_mass() const {
  assert(reac_ != ReactionType::no_reaction);
  switch (reac_) {
    // rho in final state. use already sampled mass
    case ReactionType::pi_p_pi_m_rho_z:
    case ReactionType::pi_z_pi_m_rho_m:
    case ReactionType::pi_z_pi_p_rho_p:
      return hadron_out_mass_;
    // rho in initial state, use its mass
    case ReactionType::pi_m_rho_p_pi_z:
    case ReactionType::pi_p_rho_m_pi_z:
    case ReactionType::pi_p_rho_z_pi_p:
    case ReactionType::pi_m_rho_z_pi_m:
    case ReactionType::pi_z_rho_m_pi_m:
    case ReactionType::pi_z_rho_p_pi_p:
    case ReactionType::pi_z_rho_z_pi_z:
      return (incoming_particles_[0].is_rho())
                 ? incoming_particles_[0].effective_mass()
                 : incoming_particles_[1].effective_mass();
    case ReactionType::no_reaction:
    default:
      throw std::runtime_error(
          "Invalid ReactionType in ScatterActionPhoton::rho_mass()");
  }
}

CollisionBranchList ScatterActionPhoton::photon_cross_sections(
    MediatorType mediator) {
  CollisionBranchList process_list;
  CrosssectionsPhoton<ComputationMethod::Analytic> xs_object;

  static ParticleTypePtr photon_particle = &ParticleType::find(pdg::photon);

  const double s = mandelstam_s();
  // the mass of the mediating particle depends on the channel. For an incoming
  // rho it is the mass of the incoming particle, for an outgoing rho it is the
  // sampled mass
  const double m_rho = rho_mass();
  double xsection = 0.0;

  switch (reac_) {
    case ReactionType::pi_p_pi_m_rho_z:
      xsection = xs_object.xs_pi_pi_rho0(s, m_rho);
      break;

    case ReactionType::pi_z_pi_m_rho_m:
    case ReactionType::pi_z_pi_p_rho_p:
      xsection = xs_object.xs_pi_pi0_rho(s, m_rho);
      break;

    case ReactionType::pi_m_rho_z_pi_m:
    case ReactionType::pi_p_rho_z_pi_p:
      xsection = xs_object.xs_pi_rho0_pi(s, m_rho);
      break;

    case ReactionType::pi_m_rho_p_pi_z:
    case ReactionType::pi_p_rho_m_pi_z:
      if (mediator == MediatorType::SUM) {
        xsection = xs_object.xs_pi_rho_pi0(s, m_rho);
        break;
      } else if (mediator == MediatorType::PION) {
        xsection = xs_object.xs_pi_rho_pi0_rho_mediated(s, m_rho);
        break;
      } else if (mediator == MediatorType::OMEGA) {
        xsection = xs_object.xs_pi_rho_pi0_omega_mediated(s, m_rho);
        break;
      } else {
        throw std::runtime_error("");
      }
    case ReactionType::pi_z_rho_m_pi_m:
    case ReactionType::pi_z_rho_p_pi_p:
      if (mediator == MediatorType::SUM) {
        xsection = xs_object.xs_pi0_rho_pi(s, m_rho);
        break;
      } else if (mediator == MediatorType::PION) {
        xsection = xs_object.xs_pi0_rho_pi_rho_mediated(s, m_rho);
        break;
      } else if (mediator == MediatorType::OMEGA) {
        xsection = xs_object.xs_pi0_rho_pi_omega_mediated(s, m_rho);
        break;
      } else {
        throw std::runtime_error("");
      }

    case ReactionType::pi_z_rho_z_pi_z:
      xsection = xs_object.xs_pi0_rho0_pi0(s, m_rho);
      break;

    case ReactionType::no_reaction:
      // never reached
      break;
  }

  // Due to numerical reasons it can happen that the calculated cross sections
  // are negative (approximately -1e-15) if sqrt(s) is close to the threshold
  // energy. In those cases the cross section is manually set to 0.1 mb, which
  // is a reasonable value for the processes we are looking at (C14,C15,C16).

  if (xsection <= 0) {
    xsection = 0.1;
    const auto &log = logger<LogArea::ScatterAction>();
    log.warn("Calculated negative cross section.\nParticles ",
             incoming_particles_, " mass rho particle: ", m_rho,
             ", sqrt_s: ", std::sqrt(s));
  }
  process_list.push_back(make_unique<CollisionBranch>(
      *hadron_out_t_, *photon_particle, xsection, ProcessType::TwoToTwo));
  return process_list;
}

double ScatterActionPhoton::diff_cross_section(const double t,
                                               const double m_rho,
                                               MediatorType mediator) const {
  const double s = mandelstam_s();
  double diff_xsection = 0.0;

  CrosssectionsPhoton<ComputationMethod::Analytic> xs_object;

  switch (reac_) {
    case ReactionType::pi_p_pi_m_rho_z:
      diff_xsection = xs_object.xs_diff_pi_pi_rho0(s, t, m_rho);
      break;

    case ReactionType::pi_z_pi_m_rho_m:
    case ReactionType::pi_z_pi_p_rho_p:
      diff_xsection = xs_object.xs_diff_pi_pi0_rho(s, t, m_rho);
      break;

    case ReactionType::pi_m_rho_z_pi_m:
    case ReactionType::pi_p_rho_z_pi_p:
      diff_xsection = xs_object.xs_diff_pi_rho0_pi(s, t, m_rho);
      break;

    case ReactionType::pi_m_rho_p_pi_z:
    case ReactionType::pi_p_rho_m_pi_z:
      if (mediator == MediatorType::SUM) {
        diff_xsection =
            xs_object.xs_diff_pi_rho_pi0_rho_mediated(s, t, m_rho) +
            xs_object.xs_diff_pi_rho_pi0_omega_mediated(s, t, m_rho);
      } else if (mediator == MediatorType::OMEGA) {
        diff_xsection =
            xs_object.xs_diff_pi_rho_pi0_omega_mediated(s, t, m_rho);
      } else if (mediator == MediatorType::PION) {
        diff_xsection = xs_object.xs_diff_pi_rho_pi0_rho_mediated(s, t, m_rho);
      }
      break;

    case ReactionType::pi_z_rho_m_pi_m:
    case ReactionType::pi_z_rho_p_pi_p:
      if (mediator == MediatorType::SUM) {
        diff_xsection =
            xs_object.xs_diff_pi0_rho_pi_rho_mediated(s, t, m_rho) +
            xs_object.xs_diff_pi0_rho_pi_omega_mediated(s, t, m_rho);
      } else if (mediator == MediatorType::OMEGA) {
        diff_xsection =
            xs_object.xs_diff_pi0_rho_pi_omega_mediated(s, t, m_rho);
      } else if (mediator == MediatorType::PION) {
        diff_xsection = xs_object.xs_diff_pi0_rho_pi_rho_mediated(s, t, m_rho);
      }
      break;

    case ReactionType::pi_z_rho_z_pi_z:
      diff_xsection = xs_object.xs_diff_pi0_rho0_pi0(s, t, m_rho);
      break;
    case ReactionType::no_reaction:
      // never reached
      break;
  }
  return diff_xsection;
}

double ScatterActionPhoton::diff_cross_section_w_ff(const double t,
                                                    const double m_rho,
                                                    const double E_photon) {
  /**
   * The form factor is assumed to be a hadronic dipole form factor which
   * takes the shape: FF = (2*Lambda^2/(2*Lambda^2 - t))^2 with
   * Lambda = 1.0 GeV. t depends on the lightest possible exchange particle in
   * the different channels. This could either be a pion or an omega meson. For
   * the computation the parametrizations given in (\iref{Turbide:2006})
   * are used.
   */

  /* C12, C13, C15, C16 need special treatment. These processes have identical
     incoming and outgoing particles, but diffrent mediating particles and
     hence different form factors. If both channels are added up
     (MediatorType::SUM), each contribution is corrected by the corresponding
     form factor.
  */
  switch (reac_) {
    case ReactionType::pi_m_rho_p_pi_z:
    case ReactionType::pi_p_rho_m_pi_z:
    case ReactionType::pi_z_rho_m_pi_m:
    case ReactionType::pi_z_rho_p_pi_p: {
      if (default_mediator_ == MediatorType::SUM) {
        std::pair<double, double> FF = form_factor_single(E_photon);
        std::pair<double, double> diff_xs = diff_cross_section_single(t, m_rho);
        const double xs_ff = pow_int(FF.first, 4) * diff_xs.first +
                             pow_int(FF.second, 4) * diff_xs.second;
        return xs_ff;
      } else if (default_mediator_ == MediatorType::PION) {
        const double FF = form_factor_pion(E_photon);
        const double diff_xs = diff_cross_section(t, m_rho);
        return pow_int(FF, 4) * diff_xs;
      } else if (default_mediator_ == MediatorType::OMEGA) {
        const double FF = form_factor_omega(E_photon);
        const double diff_xs = diff_cross_section(t, m_rho);
        return pow_int(FF, 4) * diff_xs;
      }
      break;
    }
    case ReactionType::pi_z_pi_p_rho_p:
    case ReactionType::pi_z_pi_m_rho_m:
    case ReactionType::pi_p_rho_z_pi_p:
    case ReactionType::pi_m_rho_z_pi_m:
    case ReactionType::pi_p_pi_m_rho_z: {
      const double FF = form_factor_pion(E_photon);
      const double xs = diff_cross_section(t, m_rho);
      const double xs_ff = pow_int(FF, 4) * xs;
      return xs_ff;
      break;
    }

    case ReactionType::pi_z_rho_z_pi_z: {
      const double FF = form_factor_omega(E_photon);
      const double xs = diff_cross_section(t, m_rho);
      const double xs_ff = pow_int(FF, 4) * xs;
      return xs_ff;
    }

    case ReactionType::no_reaction:
    default:
      throw std::runtime_error("");
      return 0;
  }
}

double ScatterActionPhoton::form_factor_pion(const double E_photon) const {
  const double Lambda = 1.0;
  const double Lambda2 = Lambda * Lambda;

  const double t_ff = 34.5096 * pow(E_photon, 0.737) -
                      67.557 * pow(E_photon, 0.7584) +
                      32.858 * pow(E_photon, 0.7806);
  const double ff = 2 * Lambda2 / (2 * Lambda2 - t_ff);

  return ff * ff;
}

double ScatterActionPhoton::form_factor_omega(const double E_photon) const {
  const double Lambda = 1.0;
  const double Lambda2 = Lambda * Lambda;

  const double t_ff =
      -61.595 * pow(E_photon, 0.9979) + 28.592 * pow(E_photon, 1.1579) +
      37.738 * pow(E_photon, 0.9317) - 5.282 * pow(E_photon, 1.3686);
  const double ff = 2 * Lambda2 / (2 * Lambda2 - t_ff);

  return ff * ff;
}

std::pair<double, double> ScatterActionPhoton::form_factor_single(
    const double E_photon) {
  return std::pair<double, double>(form_factor_pion(E_photon),
                                   form_factor_omega(E_photon));
}

std::pair<double, double> ScatterActionPhoton::diff_cross_section_single(
    const double t, const double m_rho) {
  const double diff_xs_rho = diff_cross_section(t, m_rho, MediatorType::PION);
  const double diff_xs_omega =
      diff_cross_section(t, m_rho, MediatorType::OMEGA);

  return std::pair<double, double>(diff_xs_rho, diff_xs_omega);
}

}  // namespace smash
