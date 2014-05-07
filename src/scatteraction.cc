/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/action.h"

#include "include/outputroutines.h"
#include "include/pdgcode.h"
#include "include/random.h"
#include "include/resonances.h"

namespace Smash {

ScatterAction::ScatterAction(const std::vector<int> &in_part,
                             float time_of_execution)
    : Action(in_part, time_of_execution) {}


void ScatterAction::choose_channel () {
  interaction_type_ = 0;
  if (total_weight_ > really_small) {
    double random_interaction = Random::canonical();
    float interaction_probability = 0.0;
    std::vector<ProcessBranch>::const_iterator proc = subprocesses_.begin();
    while (interaction_type_ == 0 && proc != subprocesses_.end()) {
      if (proc->pdg_list().size() > 1
          || proc->pdg_list().at(0) != PdgCode::invalid()) {
        interaction_probability += proc->weight() / total_weight_;
        if (random_interaction < interaction_probability) {
          interaction_type_ = proc->type();
          outgoing_particles_ = proc->particle_list();
        }
      }
      ++proc;
      printd("ScatterAction::choose_channel: collision type %d particle %d <-> %d time: %g\n",
             interaction_type_, incoming_particles_[0], incoming_particles_[1],
             time_of_execution_);
    }
  }
}

void ScatterAction::perform (Particles *particles, size_t &id_process)
{
  FourVector velocity_CM, neg_velocity_CM;

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

  FourVector initial_momentum(data_a.momentum() + data_b.momentum());
  FourVector final_momentum;

  printd("Process %zu type %i particle %s<->%s colliding %d<->%d time %g\n",
         id_process, interaction_type_, data_a.type(*particles).name().c_str(),
         data_a.type(*particles).name().c_str(), id_a, id_b,
         data_a.position().x0());
  printd_momenta("particle 1 momenta before", data_a);
  printd_momenta("particle 2 momenta before", data_b);

  /* 2->2 elastic scattering */
  switch (interaction_type_) {
  case 0: {
    printd("Process: Elastic collision.\n");

    /* processes computed in the center of momenta */
    boost_CM(&data_a, &data_b, &velocity_CM);
    momenta_exchange(&data_a, &data_b);
    boost_back_CM(&data_a, &data_b, &velocity_CM);

    printd_momenta("particle 1 momenta after", data_a);
    printd_momenta("particle 2 momenta after", data_b);

    final_momentum = data_a.momentum() + data_b.momentum();

    /* unset collision time for both particles + keep id + unset partner */
    data_a.set_collision_past(id_process);
    data_b.set_collision_past(id_process);

    *particles->data_pointer(id_a) = data_a;
    *particles->data_pointer(id_b) = data_b;
    outgoing_particles_[0] = data_a;
    outgoing_particles_[1] = data_b;
  } break;

  case 1: {
    /* resonance formation */
    printd("Process: Resonance formation. ");

    /* processes computed in the center of momenta */
    boost_CM(&data_a, &data_b, &velocity_CM);
    resonance_formation(particles, data_a, data_b);
    boost_back_CM(&data_a, &data_b, &velocity_CM);  // TODO(mkretz) why? can't
                                                    // we just boost a copy of
                                                    // the ParticleData objects?

    /* Boost the new particle to computational frame */
    neg_velocity_CM.set_FourVector(1.0, -velocity_CM.x1(),
      -velocity_CM.x2(), -velocity_CM.x3());

    for (ParticleData &new_particle : outgoing_particles_) {
      new_particle.set_momentum(
          new_particle.momentum().LorentzBoost(neg_velocity_CM));
      final_momentum += new_particle.momentum();

      /* The starting point of resonance is between
       * the two initial particles
       * x_middle = x_a + (x_b - x_a) / 2
       */
      FourVector middle_point =
          data_a.position() + (data_b.position() - data_a.position()) / 2.0;
      new_particle.set_position(middle_point);
      /* unset collision time for particles + keep id + unset partner */
      new_particle.set_collision_past(id_process);

      printd("Resonance %s with ID %zu \n",
             new_particle.type(*particles).name().c_str(), new_particle.id());
      printd_momenta("momentum in comp frame", new_particle);
      printd_position("position in comp frame", new_particle);

      particles->add_data(new_particle);
    }

    /* Remove the initial particles */
    particles->remove(data_a.id());
    particles->remove(data_b.id());

    printd("Particle map has now %zu elements. \n", particles->size());
  } break;

  default:
    printf("Warning: ID %i (%s) has unspecified process type %i.\n",
           id_a, particles->type(id_a).name().c_str(), interaction_type_);
  } /* end switch (interaction_type_) */

  id_process++;

  FourVector momentum_difference;
  momentum_difference += initial_momentum;
  momentum_difference -= final_momentum;
  if (fabs(momentum_difference.x0()) > really_small) {
    printf("Process %zu type %i\n", id_process, interaction_type_);
    printf("Warning: Interaction type %i E conservation violation %g\n",
           interaction_type_, momentum_difference.x0());
  }
  if (fabs(momentum_difference.x1()) > really_small)
    printf("Warning: Interaction type %i px conservation violation %g\n",
           interaction_type_, momentum_difference.x1());
  if (fabs(momentum_difference.x2()) > really_small)
    printf("Warning: Interaction type %i py conservation violation %g\n",
           interaction_type_, momentum_difference.x2());
  if (fabs(momentum_difference.x3()) > really_small)
    printf("Warning: Interaction type %i pz conservation violation %g\n",
           interaction_type_, momentum_difference.x3());
}

void ScatterAction::resonance_formation(Particles *particles, const ParticleData &particle0,
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
    resonance.set_position(1.0, 0.0, 0.0, 0.0);
  } else if (outgoing_particles_.size() == 2) {
    /* 2 particles in final state. Need another particle template */
    /* XXX: For now, it is assumed that the other particle is stable! */
    ParticleData *resonance = &outgoing_particles_.at(0);
    ParticleData *stable_product;
    if (resonance->type(*particles).width() >
        0) {  // TODO: can we change this to an is_stable or is_unstable method,
              // please?
      stable_product = &outgoing_particles_.at(1);
    } else {
      stable_product = resonance;
      resonance = &outgoing_particles_.at(1);
    }
    float mass_stable = stable_product->type(*particles).mass();
    /* Sample resonance mass */
    double mass_resonance = sample_resonance_mass(
        particles, resonance->pdgcode(), stable_product->pdgcode(), cms_energy);

    /* Sample the particle momenta */
    sample_cms_momenta(resonance, stable_product, cms_energy, mass_resonance,
                       mass_stable);

    /* Initialize positions */
    resonance->set_position(1.0, 0.0, 0.0, 0.0);
    stable_product->set_position(1.0, 0.0, 0.0, 0.0);
  } else {
    throw InvalidResonanceFormation(
        "resonance_formation: Incorrect number of particles in final state. 1 "
        "or 2 are expected. Called with " +
        std::to_string(outgoing_particles_.size()));
  }
}

}
