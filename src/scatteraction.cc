/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/action.h"

#include "include/constants.h"
#include "include/outputroutines.h"
#include "include/pdgcode.h"
#include "include/random.h"
#include "include/resonances.h"
#include "include/angles.h"

namespace Smash {

ScatterAction::ScatterAction(const ParticleData &in_part1,
                             const ParticleData &in_part2,
                             float time_of_execution)
    : Action({in_part1, in_part2}, time_of_execution) {}


void ScatterAction::perform(Particles *particles, size_t &id_process) {
  /* Relevant particle IDs for the collision. */
  int id_a = incoming_particles_[0].id();
  int id_b = incoming_particles_[1].id();

  printd("Process %zu particle %s<->%s colliding %d<->%d time %g\n",
         id_process, incoming_particles_[0].type().name().c_str(),
         incoming_particles_[1].type().name().c_str(), id_a, id_b,
         incoming_particles_[0].position().x0());
  printd_momenta("particle 1 momenta before", incoming_particles_[0]);
  printd_momenta("particle 2 momenta before", incoming_particles_[1]);

  /* Decide for a particular final state. */
  outgoing_particles_ = choose_channel();

  if (is_elastic()) {
    /* 2->2 elastic scattering */
    printd("Process: Elastic collision.\n");

    momenta_exchange();

    /* unset collision time for both particles + keep id + unset partner */
    outgoing_particles_[0].set_collision_past(id_process);
    outgoing_particles_[1].set_collision_past(id_process);

    particles->data(id_a) = outgoing_particles_[0];
    particles->data(id_b) = outgoing_particles_[1];
  } else {
    /* resonance formation */
    printd("Process: Resonance formation. ");

    /* The starting point of resonance is between the two initial particles:
     * x_middle = (x_a + x_b) / 2   */
    FourVector middle_point = (incoming_particles_[0].position()
                               + incoming_particles_[1].position()) / 2.;

    /* processes computed in the center of momenta */
    resonance_formation();

    /* Set positions & boost to computational frame. */
    for (ParticleData &new_particle : outgoing_particles_) {
      new_particle.set_position(middle_point);

      new_particle.set_momentum(
          new_particle.momentum().LorentzBoost(-beta_cm()));

      /* unset collision time for particles + keep id + unset partner */
      new_particle.set_collision_past(id_process);

      printd("Resonance %s with ID %i \n",
             new_particle.type().name().c_str(), new_particle.id());
      printd_momenta("momentum in comp frame", new_particle);
      printd_position("position in comp frame", new_particle);

      new_particle.set_id(particles->add_data(new_particle));
    }

    /* Remove the initial particles */
    particles->remove(id_a);
    particles->remove(id_b);

    printd("Particle map has now %zu elements. \n", particles->size());
  }

  check_conservation(id_process);

  id_process++;
}


ThreeVector ScatterAction::beta_cm() const {
  FourVector mom = incoming_particles_[0].momentum() +
                   incoming_particles_[1].momentum();
  return mom.threevec() / mom.x0();
}


double ScatterAction::sqrt_s() const {
  FourVector mom = incoming_particles_[0].momentum() +
                   incoming_particles_[1].momentum();
  return mom.abs();
}


double ScatterAction::particle_distance() const {
  // local copy of particles (since we need to boost them)
  ParticleData p1 = incoming_particles_[0];
  ParticleData p2 = incoming_particles_[1];
  /* Boost particles to center-of-momentum frame. */
  ThreeVector velocity = beta_cm();
  p1.boost(velocity);
  p2.boost(velocity);
  ThreeVector pos_diff = p1.position().threevec() - p2.position().threevec();
  printd("Particle %d<->%d position difference: %g %g %g %g [fm]\n",
         p1.id(), p2.id(), pos_diff.x1(), pos_diff.x2(), pos_diff.x3());
  ThreeVector mom_diff = p1.momentum().threevec() - p2.momentum().threevec();
  printd("Particle %d<->%d momentum difference: %g %g %g %g [fm]\n",
         p1.id(), p2.id(), mom_diff.x1(), mom_diff.x2(), mom_diff.x3());
  /* Zero momentum leads to infite distance. */
  if (fabs(mom_diff.sqr()) < really_small)
    return  pos_diff.sqr();

  /* UrQMD squared distance criteria:
   * arXiv:nucl-th/9803035 (3.27): in center of momemtum frame
   * position of particle a: x_a
   * position of particle b: x_b
   * velocity of particle a: v_a
   * velocity of particle b: v_b
   * d^2_{coll} = (x_a - x_b)^2 - ((x_a - x_a) . (v_a - v_b))^2 / (v_a - v_b)^2
   */
  return pos_diff.sqr()
         - (pos_diff * mom_diff) * (pos_diff * mom_diff) / mom_diff.sqr();
}


bool ScatterAction::is_elastic() const {
  return outgoing_particles_.size() == 2 &&
         outgoing_particles_[0].pdgcode() == incoming_particles_[0].pdgcode() &&
         outgoing_particles_[1].pdgcode() == incoming_particles_[1].pdgcode();
}


void ScatterAction::momenta_exchange() {
  outgoing_particles_[0] = incoming_particles_[0];
  outgoing_particles_[1] = incoming_particles_[1];

  ParticleData *p1 = &outgoing_particles_[0];
  ParticleData *p2 = &outgoing_particles_[1];

  /* Boost to CM frame. */
  ThreeVector velocity_CM = beta_cm();
  p1->boost(velocity_CM);
  p2->boost(velocity_CM);

  /* debug output */
  printd_momenta("center of momenta 1", *p1);
  printd_momenta("center of momenta 2", *p2);

  /* We are in the center of momentum,
     hence this is equal for both particles. */
  const double momentum_radial = p1->momentum().abs3();

  /* Particle exchange momenta and scatter to random direction.
   * XXX: Angles should be sampled from differential cross section
   * of this process. */
  Angles phitheta;
  phitheta.distribute_isotropically();
  printd("Random momentum: %g %g %g %g \n", momentum_radial, phitheta.phi(),
        phitheta.costheta(), phitheta.sintheta());

  /* Only direction of 3-momentum, not magnitude, changes in CM frame.
   * Thus particle energies remain the same (Lorentz boost will change them for
   * computational frame, however). */
  p1->set_3momentum(phitheta.threevec() * momentum_radial);
  p2->set_3momentum(-phitheta.threevec() * momentum_radial);

  /* debug output */
  printd_momenta("exchanged momenta 1", *p1);
  printd_momenta("exchanged momenta 2", *p2);

  /* Boost back. */
  p1->boost(-velocity_CM);
  p2->boost(-velocity_CM);
}


void ScatterAction::resonance_formation() {
  const double cms_energy = sqrt_s();

  switch (outgoing_particles_.size()) {
  case 1:
    /* 1 particle in final state: Center-of-momentum frame of initial
     * particles is the rest frame of the resonance.
     */
    outgoing_particles_[0].set_momentum(FourVector(cms_energy, 0., 0., 0.));

    printd("Momentum of the new particle: %g %g %g %g \n",
           outgoing_particles_[0].momentum().x0(),
           outgoing_particles_[0].momentum().x1(),
           outgoing_particles_[0].momentum().x2(),
           outgoing_particles_[0].momentum().x3());
    break;
  case 2:
    /* 2 particles in final state: Sample the particle momenta. */
    sample_cms_momenta(cms_energy);
    break;
  default:
    std::string s = "resonance_formation: "
                    "Incorrect number of particles in final state: ";
    s += std::to_string(outgoing_particles_.size()) + " (";
    s += incoming_particles_[0].pdgcode().string() + " + ";
    s += incoming_particles_[1].pdgcode().string() + ")";
    throw InvalidResonanceFormation(s);
  }
}

}  // namespace Smash
