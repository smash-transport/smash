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
                                       const ParticleData &in_part_b,
                                       const ParticleData &in_part_c,
                                       const ParticleData &out_part,
                                       double time)
    : Action({in_part_a, in_part_b, in_part_c}, {out_part}, time,
             ProcessType::MultiParticle) {}

void ScatterActionThree::generate_final_state() {
  /* The production point of the new particles.  */
  FourVector middle_point = get_interaction_point();

  // TODO Check that incoming and outgoing are correct (3-to-1) + maybe check
  // type

  // if (outgoing_particles_.size() != 1) {
  //   std::string s =
  //       "resonance_formation: "
  //       "Incorrect number of particles in final state: ";
  //   s += std::to_string(outgoing_particles_.size()) + " (";
  //   s += incoming_particles_[0].pdgcode().string() + " + ";
  //   s += incoming_particles_[1].pdgcode().string() + ")";
  //   throw InvalidResonanceFormation(s);
  // }

  // Set the momentum of the formed resonance in its rest frame.
  outgoing_particles_[0].set_4momentum(
      total_momentum_of_outgoing_particles().abs(), 0., 0., 0.);

  // TODO Allow formation of particles (?)

  for (ParticleData &new_particle : outgoing_particles_) {
    // Boost to the computational frame
    new_particle.boost_momentum(
        -total_momentum_of_outgoing_particles().velocity());
    /* Set positions of the outgoing particles */
    if (get_type() != ProcessType::Elastic) {
      new_particle.set_4position(middle_point);
    }
  }
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
  // TODO Distinguish rejected from accepted reaction (Possible with set final
  // state?)
  out << "MultiParticleScatter of " << incoming_particles_;
  out << " to " << outgoing_particles_;
}

}  // namespace smash
