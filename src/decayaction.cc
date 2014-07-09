/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/action.h"

#include "include/angles.h"
#include "include/outputroutines.h"
#include "include/pdgcode.h"
#include "include/resonances.h"
#include "include/decaymodes.h"

namespace Smash {


DecayAction::DecayAction(const ParticleData &in_part, float time_of_execution)
    : Action ({in_part}, time_of_execution) {}

DecayAction::DecayAction(const ParticleData &p) : Action ({p}, 0.) {
  add_processes(p.type().get_partial_widths(p.mass()));
}

void DecayAction::one_to_two(const ParticleData &incoming0) {
  ParticleData &outgoing0 = outgoing_particles_[0];
  ParticleData &outgoing1 = outgoing_particles_[1];
  const ParticleType &outgoing0_type = outgoing0.type();
  const ParticleType &outgoing1_type = outgoing1.type();

  // TODO(weil) This may be checked already when reading in the possible
  // decay channels.
  if (!outgoing0.is_hadron() || !outgoing1.is_hadron()) {
    printf("Warning: decay products A: %s B: %s\n",
           outgoing0.pdgcode().string().c_str(),
           outgoing1.pdgcode().string().c_str());
  }

  double mass_a = outgoing0_type.mass();
  double mass_b = outgoing1_type.mass();
  const double total_energy = incoming0.momentum().x0();

  /* If one of the particles is resonance, sample its mass */
  /* XXX: Other particle assumed stable! */
  if (!outgoing0_type.is_stable()) {
    mass_a = sample_resonance_mass (outgoing0_type, outgoing1_type,
                                    total_energy);
  } else if (!outgoing1_type.is_stable()) {
    mass_b = sample_resonance_mass (outgoing1_type, outgoing0_type,
                                    total_energy);
  }

  /* Sample the momenta */
  sample_cms_momenta(&outgoing0, &outgoing1, total_energy, mass_a, mass_b);

  /* Both decay products begin from the same point */
  outgoing0.set_position(incoming0.position());
  outgoing1.set_position(incoming0.position());
}


void DecayAction::one_to_three(const ParticleData &incoming0) {
  ParticleData &outgoing0 = outgoing_particles_[0];
  ParticleData &outgoing1 = outgoing_particles_[1];
  ParticleData &outgoing2 = outgoing_particles_[2];
  const ParticleType &outgoing0_type = outgoing0.type();
  const ParticleType &outgoing1_type = outgoing1.type();
  const ParticleType &outgoing2_type = outgoing2.type();

  // TODO(weil) This may be checked already when reading in the possible
  // decay channels.
  if (!outgoing0.is_hadron() || !outgoing1.is_hadron() ||
      !outgoing2.is_hadron()) {
    printf("Warning: decay products A: %s B: %s C: %s\n",
           outgoing0.pdgcode().string().c_str(),
           outgoing1.pdgcode().string().c_str(),
           outgoing2.pdgcode().string().c_str());
  }
  printd("Note: Doing 1->3 decay!\n");

  const FourVector momentum_resonance = incoming0.momentum();
  const double mass_a = outgoing0_type.mass();
  const double mass_b = outgoing1_type.mass();
  const double mass_c = outgoing2_type.mass();
  const double mass_resonance =
      sqrt(momentum_resonance.Dot(momentum_resonance));

  /* mandelstam-s limits for pairs ab and bc */
  const double s_ab_max = (mass_resonance - mass_c) * (mass_resonance - mass_c);
  const double s_ab_min = (mass_a + mass_b) * (mass_a + mass_b);
  const double s_bc_max = (mass_resonance - mass_a) * (mass_resonance - mass_a);
  const double s_bc_min = (mass_b + mass_c) * (mass_b + mass_c);

  printd("s_ab limits: %g %g \n", s_ab_min, s_ab_max);
  printd("s_bc limits: %g %g \n", s_bc_min, s_bc_max);

  /* randomly pick values for s_ab and s_bc
   * until the pair is within the Dalitz plot */
  double dalitz_bc_max = 0.0, dalitz_bc_min = 1.0;
  double s_ab = 0.0, s_bc = 0.5;
  while (s_bc > dalitz_bc_max || s_bc < dalitz_bc_min) {
    s_ab = Random::uniform(s_ab_min, s_ab_max);
    s_bc = Random::uniform(s_bc_min, s_bc_max);
    const double e_b_rest =
        (s_ab - mass_a * mass_a + mass_b * mass_b) / (2 * sqrt(s_ab));
    const double e_c_rest =
        (mass_resonance * mass_resonance - s_ab - mass_c * mass_c) /
        (2 * sqrt(s_ab));
    dalitz_bc_max = (e_b_rest + e_c_rest) * (e_b_rest + e_c_rest) -
                    (sqrt(e_b_rest * e_b_rest - mass_b * mass_b) -
                     sqrt(e_c_rest * e_c_rest - mass_c * mass_c)) *
                        (sqrt(e_b_rest * e_b_rest - mass_b * mass_b) -
                         sqrt(e_c_rest * e_c_rest - mass_c * mass_c));
    dalitz_bc_min = (e_b_rest + e_c_rest) * (e_b_rest + e_c_rest) -
                    (sqrt(e_b_rest * e_b_rest - mass_b * mass_b) +
                     sqrt(e_c_rest * e_c_rest - mass_c * mass_c)) *
                        (sqrt(e_b_rest * e_b_rest - mass_b * mass_b) +
                         sqrt(e_c_rest * e_c_rest - mass_c * mass_c));
  }

  printd("s_ab: %g s_bc: %g min: %g max: %g\n", s_ab, s_bc, dalitz_bc_min,
         dalitz_bc_max);

  /* Compute energy and momentum magnitude */
  const double energy_a =
      (mass_resonance * mass_resonance + mass_a * mass_a - s_bc) /
      (2 * mass_resonance);
  const double energy_c =
      (mass_resonance * mass_resonance + mass_c * mass_c - s_ab) /
      (2 * mass_resonance);
  const double energy_b =
      (s_ab + s_bc - mass_a * mass_a - mass_c * mass_c) / (2 * mass_resonance);
  const double momentum_a = sqrt(energy_a * energy_a - mass_a * mass_a);
  const double momentum_c = sqrt(energy_c * energy_c - mass_c * mass_c);
  const double momentum_b = sqrt(energy_b * energy_b - mass_b * mass_b);

  const double total_energy = momentum_resonance.x0();
  if (fabs(energy_a + energy_b + energy_c - total_energy) > really_small)
    printf("1->3 warning: Ea + Eb + Ec: %g Total E: %g",
           energy_a + energy_b + energy_c, total_energy);
  printd("Calculating the angles...\n");

  /* momentum_a direction is random */
  Angles phitheta;
  phitheta.distribute_isotropically();
  /* This is the angle of the plane of the three decay particles */
  outgoing0.set_momentum(mass_a, phitheta.threevec() * momentum_a);

  /* Angle between a and b */
  double theta_ab = acos(
      (energy_a * energy_b - 0.5 * (s_ab - mass_a * mass_a - mass_b * mass_b)) /
      (momentum_a * momentum_b));
  printd("theta_ab: %g Ea: %g Eb: %g sab: %g pa: %g pb: %g\n", theta_ab,
         energy_a, energy_b, s_ab, momentum_a, momentum_b);
  bool phi_has_changed = phitheta.add_to_theta(theta_ab);
  outgoing1.set_momentum(mass_b, phitheta.threevec() * momentum_b);

  /* Angle between b and c */
  double theta_bc = acos(
      (energy_b * energy_c - 0.5 * (s_bc - mass_b * mass_b - mass_c * mass_c)) /
      (momentum_b * momentum_c));
  printd("theta_bc: %g Eb: %g Ec: %g sbc: %g pb: %g pc: %g\n", theta_bc,
         energy_b, energy_c, s_bc, momentum_b, momentum_c);
  // pass information on whether phi has changed during the last adding
  // on to add_to_theta:
  phitheta.add_to_theta(theta_bc, phi_has_changed);
  outgoing2.set_momentum(mass_c, phitheta.threevec() * momentum_c);

  /* Momentum check */
  FourVector ptot = outgoing0.momentum() + outgoing1.momentum() +
                    outgoing2.momentum();

  if (fabs(ptot.x0() - total_energy) > really_small) {
    printf("1->3 energy not conserved! Before: %g After: %g\n", total_energy,
           ptot.x0());
  }
  if (fabs(ptot.x1()) > really_small || fabs(ptot.x2()) > really_small ||
      fabs(ptot.x3()) > really_small) {
    printf("1->3 momentum check failed. Total momentum: %g %g %g\n", ptot.x1(),
           ptot.x2(), ptot.x3());
  }

  /* All decay products begin from the same point */
  outgoing0.set_position(incoming0.position());
  outgoing1.set_position(incoming0.position());
  outgoing2.set_position(incoming0.position());

  printd("p0: %g %g %g \n", outgoing0.momentum().x0(),
         outgoing1.momentum().x0(), outgoing2.momentum().x0());
  printd("p1: %g %g %g \n", outgoing0.momentum().x1(),
         outgoing1.momentum().x1(), outgoing2.momentum().x1());
  printd("p2: %g %g %g \n", outgoing0.momentum().x2(),
         outgoing1.momentum().x2(), outgoing2.momentum().x2());
  printd("p3: %g %g %g \n", outgoing0.momentum().x3(),
         outgoing1.momentum().x3(), outgoing2.momentum().x3());
}


void DecayAction::perform(Particles *particles, size_t &id_process) {
  /* Check if particle still exists. */
  if (!is_valid(*particles)) {
    printf("DecayAction::perform: ID %i not found!\n", incoming_particles_[0].id());
    return;
  }

  /* local copy of the initial state */
  ParticleData particle0 = incoming_particles_[0];

  printd("Process: Resonance decay. ");
  printd_momenta("Resonance momenta before decay", particle0);

  /* boost to rest frame */
  ThreeVector velocity_CM = particle0.velocity();
  particle0.boost(velocity_CM);

  printd_momenta("Boosted resonance momenta before decay", particle0);

  /*
   * Execute a decay process for the selected particle.
   *
   * Randomly select one of the decay modes of the particle
   * according to their relative weights. Then decay the particle
   * by calling function one_to_two or one_to_three.
   */
  outgoing_particles_ = choose_channel();

  if (outgoing_particles_.size() == 2) {
    one_to_two(particle0);
  } else if (outgoing_particles_.size() == 3) {
    one_to_three(particle0);
  } else {
    throw InvalidDecay(
        "DecayAction::perform: Only 1->2 or 1->3 processes are supported. "
        "Decay from 1->" +
        std::to_string(outgoing_particles_.size()) + " was requested.");
  }

  for (auto &p : outgoing_particles_) {
    printd_momenta("particle momenta in lrf", p);
    p.boost(-velocity_CM);
    printd_momenta("particle momenta in comp", p);
    // unset collision time for both particles + keep id + unset partner
    p.set_collision_past(id_process);
  }

  id_process++;

  check_conservation(id_process);

  /* Remove decayed particle */
  particles->remove(particle0.id());
  printd("ID %i has decayed and removed from the list.\n", particle0.id());

  for (auto &p : outgoing_particles_) {
    p.set_id(particles->add_data(p));
  }
  printd("Particle map has now %zu elements. \n", particles->size());
}

}
