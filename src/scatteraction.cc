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

namespace Smash {

ScatterAction::ScatterAction(const std::vector<int> &in_part,
                             float time_of_execution)
    : Action(in_part, time_of_execution) {}


void ScatterAction::choose_channel () {
  if (total_weight_ > really_small) {
    double random_interaction = Random::canonical();
    float interaction_probability = 0.0;
    std::vector<ProcessBranch>::const_iterator proc = subprocesses_.begin();
    while (outgoing_particles_.size() == 0 && proc != subprocesses_.end()) {
      if (proc->pdg_list().size() > 1
          || proc->pdg_list().at(0) != PdgCode::invalid()) {
        interaction_probability += proc->weight() / total_weight_;
        if (random_interaction < interaction_probability) {
          outgoing_particles_ = proc->particle_list();
        }
      }
      ++proc;
      printd("ScatterAction::choose_channel: particle %d <-> %d time: %g\n",
             incoming_particles_[0], incoming_particles_[1],
             time_of_execution_);
    }
  }
}

void ScatterAction::perform (Particles *particles, size_t &id_process)
{
  ThreeVector velocity_CM;

  /* Relevant particle IDs for the collision. */
  int id_a = incoming_particles_[0];
  int id_b = incoming_particles_[1];

  /* Check if particles still exist. */
  if (!is_valid(*particles)) {
    printd("ScatterAction::perform: ID %i or %i not found!\n", id_a, id_b);
    return;
  }

  ParticleData data_a = particles->data(id_a);
  ParticleData data_b = particles->data(id_b);

  printd("Process %zu particle %s<->%s colliding %d<->%d time %g\n",
         id_process, data_a.type().name().c_str(),
         data_a.type().name().c_str(), id_a, id_b, data_a.position().x0());
  printd_momenta("particle 1 momenta before", data_a);
  printd_momenta("particle 2 momenta before", data_b);

  if (is_elastic(particles)) {
    /* 2->2 elastic scattering */
    printd("Process: Elastic collision.\n");

    /* processes computed in the center of momenta */
    velocity_CM = boost_CM(&data_a, &data_b);
    momenta_exchange(&data_a, &data_b);
    boost_back_CM(&data_a, &data_b, velocity_CM);

    printd_momenta("particle 1 momenta after", data_a);
    printd_momenta("particle 2 momenta after", data_b);

    /* unset collision time for both particles + keep id + unset partner */
    data_a.set_collision_past(id_process);
    data_b.set_collision_past(id_process);

    particles->data(id_a) = data_a;
    particles->data(id_b) = data_b;
    outgoing_particles_[0] = data_a;
    outgoing_particles_[1] = data_b;
  } else {
    /* resonance formation */
    printd("Process: Resonance formation. ");

    /* processes computed in the center of momenta */
    velocity_CM = boost_CM(&data_a, &data_b);
    resonance_formation(data_a, data_b);
    boost_back_CM(&data_a, &data_b, velocity_CM);  // TODO(mkretz) why? can't
                                                    // we just boost a copy of
                                                    // the ParticleData objects?

    /* Boost the new particle to computational frame */
    for (ParticleData &new_particle : outgoing_particles_) {
      new_particle.set_momentum(
          new_particle.momentum().LorentzBoost(-velocity_CM));

      /* The starting point of resonance is between
       * the two initial particles
       * x_middle = (x_a + x_b) / 2
       */
      FourVector middle_point = (data_a.position() + data_b.position()) / 2.;
      new_particle.set_position(middle_point);
      /* unset collision time for particles + keep id + unset partner */
      new_particle.set_collision_past(id_process);

      printd("Resonance %s with ID %zu \n",
             new_particle.type(*particles).name().c_str(), new_particle.id());
      printd_momenta("momentum in comp frame", new_particle);
      printd_position("position in comp frame", new_particle);

      new_particle.set_id(particles->add_data(new_particle));
    }

    check_conservation(*particles, id_process);

    /* Remove the initial particles */
    particles->remove(data_a.id());
    particles->remove(data_b.id());

    printd("Particle map has now %zu elements. \n", particles->size());
  }

  id_process++;
}


bool ScatterAction::is_elastic(Particles *particles) const {
  return outgoing_particles_.size()==2 &&
         outgoing_particles_[0].pdgcode() == particles->data(incoming_particles_[0]).pdgcode() &&
         outgoing_particles_[1].pdgcode() == particles->data(incoming_particles_[1]).pdgcode();
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

    /* Initialize position */
    resonance.set_position(FourVector(1., 0., 0., 0.));
  } else if (outgoing_particles_.size() == 2) {
    /* 2 particles in final state. Need another particle template */
    /* XXX: For now, it is assumed that the other particle is stable! */
    ParticleData *resonance = &outgoing_particles_.at(0);
    ParticleData *stable_product;
    if (!resonance->type().is_stable()) {
      stable_product = &outgoing_particles_.at(1);
    } else {
      stable_product = resonance;
      resonance = &outgoing_particles_.at(1);
    }
    float mass_stable = stable_product->type().mass();
    /* Sample resonance mass */
    double mass_resonance = sample_resonance_mass(resonance->pdgcode(),
                                                  stable_product->pdgcode(),
                                                  cms_energy);

    /* Sample the particle momenta */
    sample_cms_momenta(resonance, stable_product, cms_energy, mass_resonance,
                       mass_stable);

    /* Initialize positions */
    resonance->set_position(FourVector(1., 0., 0., 0.));
    stable_product->set_position(FourVector(1., 0., 0., 0.));
  } else {
    throw InvalidResonanceFormation(
        "resonance_formation: Incorrect number of particles in final state. 1 "
        "or 2 are expected. Called with " +
        std::to_string(outgoing_particles_.size()));
  }
}

}
