/*
 *
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/scatteractionthree.h"


namespace smash {

ScatterActionThree::ScatterActionThree(const ParticleData &in_part_a,
                             const ParticleData &in_part_b, const ParticleData &in_part_c, const ParticleData &out_part, double time)
    : Action({in_part_a, in_part_b, in_part_c}, {out_part}, time, ProcessType::MultiParticle) {}


void ScatterActionThree::generate_final_state() {

  // TODO Implementation //

  // logg[LScatterActionThree].debug("Incoming particles: ", incoming_particles_);
  //
  // // /* Decide for a particular final state. */
  // // const CollisionBranch *proc = choose_channel<CollisionBranch>(
  // //     collision_channels_, total_cross_section_);
  // // process_type_ = proc->get_type();
  // // outgoing_particles_ = proc->particle_list();
  // // partial_cross_section_ = proc->weight();
  //
  // logg[LScatterAction].debug("Chosen channel: ", process_type_,
  //                            outgoing_particles_);
  //
  // /* The production point of the new particles.  */
  // FourVector middle_point = get_interaction_point();
  //
  // switch (process_type_) {
  //   case ProcessType::Elastic:
  //     /* 2->2 elastic scattering */
  //     elastic_scattering();
  //     break;
  //   case ProcessType::TwoToOne:
  //     /* resonance formation */
  //     resonance_formation();
  //     break;
  //   case ProcessType::TwoToTwo:
  //     /* 2->2 inelastic scattering */
  //     /* Sample the particle momenta in CM system. */
  //     inelastic_scattering();
  //     break;
  //   case ProcessType::StringSoftSingleDiffractiveAX:
  //   case ProcessType::StringSoftSingleDiffractiveXB:
  //   case ProcessType::StringSoftDoubleDiffractive:
  //   case ProcessType::StringSoftAnnihilation:
  //   case ProcessType::StringSoftNonDiffractive:
  //   case ProcessType::StringHard:
  //     string_excitation();
  //     break;
  //   default:
  //     throw InvalidScatterAction(
  //         "ScatterActionThree::generate_final_state: Invalid process type " +
  //         std::to_string(static_cast<int>(process_type_)) + " was requested. " +
  //         "(PDGcode1=" + incoming_particles_[0].pdgcode().string() +
  //         ", PDGcode2=" + incoming_particles_[1].pdgcode().string() + ")");
  // }
  //
  // for (ParticleData &new_particle : outgoing_particles_) {
  //   // Boost to the computational frame
  //   new_particle.boost_momentum(
  //       -total_momentum_of_outgoing_particles().velocity());
  //   /* Set positions of the outgoing particles */
  //   if (proc->get_type() != ProcessType::Elastic) {
  //     new_particle.set_4position(middle_point);
  //   }
  // }
}

double ScatterActionThree::get_total_weight() const {
  // No cross section for mulitparticle scatterings
  // TODO Maybe put probability here?
  return 0.0;
}

double ScatterActionThree::get_partial_weight() const {
  // No cross section for mulitparticle scatterings
  // TODO Maybe put probability here?
  return 0.0;
}

void ScatterActionThree::format_debug_output(std::ostream &out) const {
  out << "MultiParticleScatter of " << incoming_particles_;
  out << " to " << outgoing_particles_;
}

}  // namespace smash
