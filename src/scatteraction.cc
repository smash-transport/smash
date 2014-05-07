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
  size_t new_particles, id_new;

  /* Relevant particle IDs for the collision. */
  int id_a = incoming_particles_[0];
  int id_b = incoming_particles_[1];

  /* Check if particles still exist. */
  if (!is_valid(*particles)) {
    printd("ScatterAction::perform: ID %i or %i not found!\n", id_a, id_b);
    return;
  }

  FourVector initial_momentum(particles->data(id_a).momentum()
    + particles->data(id_b).momentum());
  FourVector final_momentum;

  printd("Process %zu type %i particle %s<->%s colliding %d<->%d time %g\n",
         id_process, interaction_type_, particles->type(id_a).name().c_str(),
         particles->type(id_a).name().c_str(), id_a, id_b,
         particles->data(id_a).position().x0());
  printd_momenta("particle 1 momenta before", particles->data(id_a));
  printd_momenta("particle 2 momenta before", particles->data(id_b));

  /* 2->2 elastic scattering */
  switch (interaction_type_) {
  case 0:
    printd("Process: Elastic collision.\n");

    /* processes computed in the center of momenta */
    boost_CM(particles->data_pointer(id_a), particles->data_pointer(id_b),
             &velocity_CM);
    momenta_exchange(particles->data_pointer(id_a),
                     particles->data_pointer(id_b));
    boost_back_CM(particles->data_pointer(id_a),
                  particles->data_pointer(id_b), &velocity_CM);

    printd_momenta("particle 1 momenta after", particles->data(id_a));
    printd_momenta("particle 2 momenta after", particles->data(id_b));

    final_momentum = particles->data(id_a).momentum()
    + particles->data(id_b).momentum();

    /* unset collision time for both particles + keep id + unset partner */
    particles->data_pointer(id_a)->set_collision_past(id_process);
    particles->data_pointer(id_b)->set_collision_past(id_process);
    break;

  case 1:
    /* resonance formation */
    printd("Process: Resonance formation. ");
    new_particles = outgoing_particles_.size();
    /* processes computed in the center of momenta */
    boost_CM(particles->data_pointer(id_a), particles->data_pointer(id_b),
             &velocity_CM);

    id_new = resonance_formation(particles, id_a, id_b, outgoing_particles_);

    boost_back_CM(particles->data_pointer(id_a),
                  particles->data_pointer(id_b), &velocity_CM);

    /* Boost the new particle to computational frame */
    neg_velocity_CM.set_FourVector(1.0, -velocity_CM.x1(),
      -velocity_CM.x2(), -velocity_CM.x3());

    for (size_t id_value = id_new; id_value < id_new + new_particles;
         id_value++) {
      particles->data_pointer(id_value)->set_momentum(
        particles->data(id_value).momentum().LorentzBoost(neg_velocity_CM));
      final_momentum += particles->data(id_value).momentum();

      /* The starting point of resonance is between
       * the two initial particles
       * x_middle = x_a + (x_b - x_a) / 2
       */
      FourVector middle_point = particles->data(id_a).position()
        + (particles->data(id_b).position()
            - particles->data(id_a).position())
        / 2.0;
      particles->data_pointer(id_value)->set_position(middle_point);
      /* unset collision time for particles + keep id + unset partner */
      particles->data_pointer(id_value)->set_collision_past(id_process);

      printd("Resonance %s with ID %zu \n",
             particles->type(id_new).name().c_str(), id_new);
      printd_momenta("momentum in comp frame", particles->data(id_new));
      printd_position("position in comp frame", particles->data(id_new));
    }

    /* Remove the initial particles */
    particles->remove(id_a);
    particles->remove(id_b);

    printd("Particle map has now %zu elements. \n", particles->size());
    break;

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

int ScatterAction::resonance_formation(Particles *particles, int particle_id,
                                       int other_id,
                                       const ParticleList &produced_particles) {
  if (produced_particles.empty()) {
    printf("resonance_formation:\n");
    printf("Warning: No final state particles found!\n");
    printf("Resonance formation canceled. Returning -1.\n");
    return -1;
  }

  const double cms_energy = particles->data(particle_id).momentum().x0()
    + particles->data(other_id).momentum().x0();

  int id_first_new = -1;
  if (produced_particles.size() == 1) {
    ParticleData resonance = produced_particles.at(0);
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
    id_first_new = particles->add_data(resonance);
  } else if (produced_particles.size() == 2) {
    /* 2 particles in final state. Need another particle template */
    /* XXX: For now, it is assumed that the other particle is stable! */
    ParticleData stable_product;
    ParticleData resonance;
    if (produced_particles.at(0).type(*particles).width() > 0) {
      resonance = produced_particles.at(0);
      stable_product = produced_particles.at(1);
    } else {
      stable_product = produced_particles.at(0);
      resonance = produced_particles.at(1);
    }
    float mass_stable
      = particles->particle_type(stable_product.pdgcode()).mass();
    /* Sample resonance mass */
    double mass_resonance = sample_resonance_mass(particles,
      resonance.pdgcode(), stable_product.pdgcode(), cms_energy);

    /* Sample the particle momenta */
    sample_cms_momenta(&resonance, &stable_product, cms_energy, mass_resonance,
                       mass_stable);

    /* Initialize positions */
    resonance.set_position(1.0, 0.0, 0.0, 0.0);
    stable_product.set_position(1.0, 0.0, 0.0, 0.0);
    id_first_new = particles->add_data(resonance);
    particles->add_data(stable_product);
  } else {
    printf("resonance_formation:\n");
    printf("Warning: %zu particles in final state!\n",
           produced_particles.size());
    printf("Resonance formation canceled. Returning -1.\n");
    return -1;
  }
  /* Return the id of the first new particle */
  return id_first_new;
}

}
