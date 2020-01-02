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

namespace smash {
static constexpr int LScatterAction = LogArea::ScatterAction::id;

BremsstrahlungAction::BremsstrahlungAction(
    const ParticleList &in, const double time, const int n_frac_photons,
    const double hadronic_cross_section_input)
    : ScatterAction(in[0], in[1], time),
      reac_(bremsstrahlung_reaction_type(in)),
      number_of_fractional_photons_(n_frac_photons),
      hadronic_cross_section_(hadronic_cross_section_input) {}

BremsstrahlungAction::ReactionType BremsstrahlungAction::bremsstrahlung_reaction_type(
    const ParticleList &in) {
    if (in.size() != 2) {
      return ReactionType::no_reaction;
    }

    PdgCode a = in[0].pdgcode();
    PdgCode b = in[1].pdgcode();

    switch (pack(a.code(), b.code())) {
      case (pack(pdg::pi_z, pdg::pi_z)):
        return ReactionType::pi_z_pi_z;

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

      default:
        return ReactionType::no_reaction;
      }
    };

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
};

void BremsstrahlungAction::generate_final_state() {
  // we have only one reaction per incoming particle pair
  if (collision_processes_bremsstrahlung_.size() != 1) {
    logg[LScatterAction].fatal()
        << "Problem in BremsstrahlungAction::generate_final_state().\n";
    throw std::runtime_error("");
  }

  auto *proc = collision_processes_bremsstrahlung_[0].get();

  outgoing_particles_ = proc->particle_list();
  process_type_ = proc->get_type();
  FourVector interaction_point = get_interaction_point();

  // This samples the phase space isotropically in the local rest frame
  sample_3body_phasespace();

  // Set position and formation time and boost back to computational frame
  for (auto &new_particle : outgoing_particles_) {
    // assuming decaying particles are always fully formed
    new_particle.set_formation_time(time_of_execution_);
    new_particle.set_4position(interaction_point);
    new_particle.boost_momentum(-total_momentum_of_outgoing_particles().velocity());
  }

  // Weighing of the fractional photons
  if (number_of_fractional_photons_ > 1) {
    throw std::runtime_error("Fractional photons currently not implemented for bremsstrahlung.");
  } else {
    weight_ = proc->weight() / hadronic_cross_section();
  }

  // Photons are not really part of the normal processes, so we have to set a
  // constant arbitrary number.
  const auto id_process = ID_PROCESS_PHOTON;
  Action::check_conservation(id_process);
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
  static ParticleTypePtr photon_particle = &ParticleType::find(pdg::photon);

  double xsection = 10.0;
  process_list.push_back(make_unique<CollisionBranch>(
      incoming_particles_[0].type(), incoming_particles_[1].type(), *photon_particle, xsection, ProcessType::Bremsstrahlung));
  return process_list;

};

} // namespace smash
