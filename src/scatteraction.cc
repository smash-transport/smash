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

ScatterAction::ScatterAction(const ParticleList &in_part,
                             float time_of_execution)
    : Action(in_part, time_of_execution) {}


void ScatterAction::perform (Particles *particles, size_t &id_process)
{
  /* Relevant particle IDs for the collision. */
  int id_a = incoming_particles_[0].id();
  int id_b = incoming_particles_[1].id();

  /* Check if particles still exist. */
  if (!is_valid(*particles)) {
    printd("ScatterAction::perform: ID %i or %i not found!\n", id_a, id_b);
    return;
  }

  printd("Process %zu particle %s<->%s colliding %d<->%d time %g\n",
         id_process, incoming_particles_[0].type().name().c_str(),
         incoming_particles_[1].type().name().c_str(), id_a, id_b, incoming_particles_[0].position().x0());
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

    ParticleData data_a = incoming_particles_[0];
    ParticleData data_b = incoming_particles_[1];

    /* The starting point of resonance is between the two initial particles:
     * x_middle = (x_a + x_b) / 2   */
    FourVector middle_point = (data_a.position() + data_b.position()) / 2.;

    /* processes computed in the center of momenta */
    ThreeVector velocity_CM = boost_CM(&data_a, &data_b);
    resonance_formation(data_a, data_b);

    /* Set positions & boost to computational frame. */
    for (ParticleData &new_particle : outgoing_particles_) {
      new_particle.set_position(middle_point);
 
      new_particle.set_momentum(
          new_particle.momentum().LorentzBoost(-velocity_CM));

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


bool ScatterAction::is_elastic() const {
  return outgoing_particles_.size()==2 &&
         outgoing_particles_[0].pdgcode() == incoming_particles_[0].pdgcode() &&
         outgoing_particles_[1].pdgcode() == incoming_particles_[1].pdgcode();
}


void ScatterAction::momenta_exchange() {

  outgoing_particles_[0] = incoming_particles_[0];
  outgoing_particles_[1] = incoming_particles_[1];

  ParticleData *p1 = &outgoing_particles_[0];
  ParticleData *p2 = &outgoing_particles_[1];

  ThreeVector velocity_CM = boost_CM(p1, p2);

  /* debug output */
  printd_momenta("center of momenta 1", *p1);
  printd_momenta("center of momenta 2", *p2);

  /* We are in the center of momentum, hence this is equal for both particles. */
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
  p1->set_3momentum(  phitheta.threevec() * momentum_radial);
  p2->set_3momentum(- phitheta.threevec() * momentum_radial);

  /* debug output */
  printd_momenta("exchanged momenta 1", *p1);
  printd_momenta("exchanged momenta 2", *p2);

  boost_back_CM(p1, p2, velocity_CM);
}


void ScatterAction::resonance_formation(const ParticleData &particle0,
                                        const ParticleData &particle1) {
  const double cms_energy =
      particle0.momentum().x0() + particle1.momentum().x0();

  if (outgoing_particles_.size() == 1) {
    ParticleData &resonance = outgoing_particles_.at(0);
    /* Center-of-momentum frame of initial particles
     * is the rest frame of the resonance
     *
     * We use fourvector to set 4-momentum, as setting it
     * with doubles requires that particle is on
     * mass shell, which is not generally true for resonances
     */
    FourVector resonance_momentum(cms_energy, 0.0, 0.0, 0.0);
    resonance.set_momentum(resonance_momentum);

    printd("Momentum of the new particle: %g %g %g %g \n",
      resonance.momentum().x0(),
      resonance.momentum().x1(),
      resonance.momentum().x2(),
      resonance.momentum().x3());

  } else if (outgoing_particles_.size() == 2) {
    /* 2 particles in final state. */

    /* Sample the particle momenta */
    sample_cms_momenta(cms_energy);

  } else {
    std::string s = "resonance_formation: Incorrect number of particles in final state: ";
    s += std::to_string(outgoing_particles_.size()) + " (";
    s += particle0.pdgcode().string() + " + " + particle1.pdgcode().string() + ")";
    throw InvalidResonanceFormation(s);
  }
}

}
