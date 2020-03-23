/*
 *
 *    Copyright (c) 2016-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/bremsstrahlungaction.h"

#include "smash/outputinterface.h"
#include "smash/random.h"

#include "smash/interpolation2D.h"

namespace smash {
static constexpr int LScatterAction = LogArea::ScatterAction::id;

BremsstrahlungAction::BremsstrahlungAction(
    const ParticleList &in, const double time, const int n_frac_photons,
    const double hadronic_cross_section_input)
    : ScatterAction(in[0], in[1], time),
      reac_(bremsstrahlung_reaction_type(in)),
      number_of_fractional_photons_(n_frac_photons),
      hadronic_cross_section_(hadronic_cross_section_input) {}

BremsstrahlungAction::ReactionType
BremsstrahlungAction::bremsstrahlung_reaction_type(const ParticleList &in) {
  if (in.size() != 2) {
    return ReactionType::no_reaction;
  }

  PdgCode a = in[0].pdgcode();
  PdgCode b = in[1].pdgcode();

  switch (pack(a.code(), b.code())) {
    case (pack(pdg::pi_z, pdg::pi_m)):
    case (pack(pdg::pi_m, pdg::pi_z)):
      return ReactionType::pi_z_pi_m;

    case (pack(pdg::pi_z, pdg::pi_p)):
    case (pack(pdg::pi_p, pdg::pi_z)):
      return ReactionType::pi_z_pi_p;

    case (pack(pdg::pi_m, pdg::pi_p)):
    case (pack(pdg::pi_p, pdg::pi_m)):
      return ReactionType::pi_p_pi_m;

    case (pack(pdg::pi_m, pdg::pi_m)):
      return ReactionType::pi_m_pi_m;

    case (pack(pdg::pi_p, pdg::pi_p)):
      return ReactionType::pi_p_pi_p;

    case (pack(pdg::pi_z, pdg::pi_z)):
      return ReactionType::pi_z_pi_z;

    default:
      return ReactionType::no_reaction;
  }
}

void BremsstrahlungAction::perform_bremsstrahlung(const OutputsList &outputs) {
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

void BremsstrahlungAction::generate_final_state() {
  // we have only one reaction per incoming particle pair
  if (collision_processes_bremsstrahlung_.size() != 1) {
    logg[LScatterAction].fatal()
        << "Problem in BremsstrahlungAction::generate_final_state().\nThe "
           "brocess branch has "
        << collision_processes_bremsstrahlung_.size()
        << " entries. It should however have 1.";
    throw std::runtime_error("");
  }

  auto *proc = collision_processes_bremsstrahlung_[0].get();

  outgoing_particles_ = proc->particle_list();
  process_type_ = proc->get_type();
  FourVector interaction_point = get_interaction_point();

  // Sample k and theta:
  double k_max =
      (sqrt_s() * sqrt_s() - 2 * outgoing_particles_[0].type().mass() * 2 *
                                 outgoing_particles_[1].type().mass()) /
      (2 * sqrt_s());
  k_ = random::uniform(0.001, k_max);
  theta_ = random::uniform(0.0, M_PI);

  // Sample the phase space anisotropically in the local rest frame
  sample_3body_phasespace();

  // ToDo: Update implementation and saveguards if sqrts, k, theta are out of
  // range
  if (number_of_fractional_photons_ > 1) {
    double diff_xs_theta = 0.0;
    double diff_xs_k = 0.0;
    double energy = sqrt_s();
    if (reac_ == ReactionType::pi_m_pi_m || reac_ == ReactionType::pi_p_pi_p) {
      diff_xs_k = (*pipi_same_charge_interpolation_diff_sigma_k)(k_, energy);
      diff_xs_theta =
          (*pipi_same_charge_interpolation_diff_sigma_theta)(theta_, energy);
      if (diff_xs_k < 0.0) {
        diff_xs_k = really_small;
      }

      if (diff_xs_theta < 0.0) {
        diff_xs_theta = really_small;
      }
    }

    // Assign weighting factor
    const double W_theta = diff_xs_theta * (M_PI - 0.0);
    const double W_k = diff_xs_k * (k_max - 0.0);
    weight_ = sqrt(W_theta * W_k) /
              (number_of_fractional_photons_ * hadronic_cross_section());
  } else {
    weight_ = proc->weight() / hadronic_cross_section();
  }

  // Set position and formation time and boost back to computational frame
  for (auto &new_particle : outgoing_particles_) {
    // assuming decaying particles are always fully formed
    new_particle.set_formation_time(time_of_execution_);
    new_particle.set_4position(interaction_point);
    new_particle.boost_momentum(
        -total_momentum_of_outgoing_particles().velocity());
  }

  // Photons are not really part of the normal processes, so we have to set a
  // constant arbitrary number.
  const auto id_process = ID_PROCESS_PHOTON;
  Action::check_conservation(id_process);
}

void BremsstrahlungAction::sample_3body_phasespace() {
  assert(outgoing_part.size() == 3);
  const double m_a = outgoing_particles_[0].type().mass(),
               m_b = outgoing_particles_[1].type().mass(),
               m_c = outgoing_particles_[2].type().mass();

  const double sqrts = sqrt_s();
  const double E_ab = sqrts - m_c - k_;  // Ekin of the pion pair in cm frame
  const double pcm =
      pCM(sqrts, E_ab, m_c);  // cm momentum of (pion pair - photon)
  const double pcm_pions = pCM(E_ab, m_a, m_b);  // cm momentum within pion pair

  // Photon angle: Phi random, theta from theta_ sampled above
  const Angles phitheta_photon(random::uniform(0.0, twopi), std::cos(theta_));
  outgoing_particles_[2].set_4momentum(m_c, pcm * phitheta_photon.threevec());
  const ThreeVector beta_cm =
      pcm * phitheta_photon.threevec() / std::sqrt(pcm * pcm + E_ab * E_ab);

  // Sample pion pair isotropically
  Angles phitheta;
  phitheta.distribute_isotropically();
  outgoing_particles_[0].set_4momentum(m_a, pcm_pions * phitheta.threevec());
  outgoing_particles_[1].set_4momentum(m_b, -pcm_pions * phitheta.threevec());
  outgoing_particles_[0].boost_momentum(beta_cm);
  outgoing_particles_[1].boost_momentum(beta_cm);
}

void BremsstrahlungAction::add_dummy_hadronic_process(
    double reaction_cross_section) {
  CollisionBranchPtr dummy_process = make_unique<CollisionBranch>(
      incoming_particles_[0].type(), incoming_particles_[1].type(),
      reaction_cross_section, ProcessType::Bremsstrahlung);
  add_collision(std::move(dummy_process));
}

CollisionBranchList BremsstrahlungAction::brems_cross_sections() {
  CollisionBranchList process_list;
  // ParticleList final_state_particles;
  static const ParticleTypePtr photon_particle =
      &ParticleType::find(pdg::photon);
  static const ParticleTypePtr pi_z_particle = &ParticleType::find(pdg::pi_z);
  static const ParticleTypePtr pi_p_particle = &ParticleType::find(pdg::pi_p);
  static const ParticleTypePtr pi_m_particle = &ParticleType::find(pdg::pi_m);

  // Create interpolation object, if not yet existent
  if (pipi_opp_charge_interpolation == nullptr ||
      pipi_same_charge_interpolation == nullptr ||
      pi0pi_interpolation == nullptr) {
    create_interpolations();
  }

  // Find cross section corresponding to given sqrt(s)
  double sqrts = sqrt_s();
  double xsection;

  if (reac_ == ReactionType::pi_p_pi_m) {
    // Here the final state is determined by the the final state provided by the
    // picked process

    // In the case of two oppositely charged pions as incoming particles,
    // there are two potential final states: pi+ + pi- and pi0 + pi0
    double xsection_pipi = (*pipi_opp_charge_interpolation)(sqrts);
    double xsection_pi0pi0 = (*pip_pim_pi0_pi0_interpolation)(sqrts);

    if (xsection_pipi <= 0.0) {
      xsection_pipi = really_small;
    }

    if (xsection_pi0pi0 <= 0.0) {
      xsection_pi0pi0 = really_small;
    }
    // Necessary only to decide for a final state with pi+ and pi- as incoming
    // particles.
    CollisionBranchList process_list_pipi;

    // Add both processes to the process_list
    process_list_pipi.push_back(make_unique<CollisionBranch>(
        incoming_particles_[0].type(), incoming_particles_[1].type(),
        *photon_particle, xsection_pipi, ProcessType::Bremsstrahlung));
    process_list_pipi.push_back(make_unique<CollisionBranch>(
        *pi_z_particle, *pi_z_particle, *photon_particle, xsection_pi0pi0,
        ProcessType::Bremsstrahlung));

    // Decide for one of the possible final states
    double total_cross_section = xsection_pipi + xsection_pi0pi0;
    const CollisionBranch *proc =
        choose_channel<CollisionBranch>(process_list_pipi, total_cross_section);

    xsection = proc->weight();

    if (xsection <= 0.0) {
      xsection = really_small;
    }

    process_list.push_back(make_unique<CollisionBranch>(
        proc->particle_list()[0].type(), proc->particle_list()[1].type(),
        *photon_particle, xsection, ProcessType::Bremsstrahlung));

  } else if (reac_ == ReactionType::pi_m_pi_m ||
             reac_ == ReactionType::pi_p_pi_p ||
             reac_ == ReactionType::pi_z_pi_m ||
             reac_ == ReactionType::pi_z_pi_p) {
    // Here the final state hadrons are identical to the initial state hadrons
    if (reac_ == ReactionType::pi_m_pi_m || reac_ == ReactionType::pi_p_pi_p) {
      xsection = (*pipi_same_charge_interpolation)(sqrts);
    } else {
      // Pi0 in initial state
      xsection = (*pi0pi_interpolation)(sqrts);
    }

    if (xsection <= 0.0) {
      xsection = really_small;
    }

    process_list.push_back(make_unique<CollisionBranch>(
        incoming_particles_[0].type(), incoming_particles_[1].type(),
        *photon_particle, xsection, ProcessType::Bremsstrahlung));

  } else if (reac_ == ReactionType::pi_z_pi_z) {
    // Here we have a hard-coded final state that differs from the initial
    // state, namely: pi0 + pi0 -> pi+ + pi- + gamma
    xsection = (*pi0_pi0_pip_pim_interpolation)(sqrts);

    if (xsection <= 0.0) {
      xsection = really_small;
    }
    process_list.push_back(make_unique<CollisionBranch>(
        *pi_p_particle, *pi_m_particle, *photon_particle, xsection,
        ProcessType::Bremsstrahlung));
  } else {
    throw std::runtime_error("Unknown ReactionType in BremsstrahlungAction.");
  }

  return process_list;
}

}  // namespace smash
