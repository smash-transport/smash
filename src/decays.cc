/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */

#include <cstdio>
#include <cmath>
#include <list>
#include <map>
#include <vector>

#include "include/decays.h"

#include "include/Parameters.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/constants.h"
#include "include/FourVector.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/particles.h"

/* does_decay - does a resonance decay on this timestep? */
bool does_decay(ParticleData *particle, ParticleType *particle_type,
  std::list<int> *decay_list, const Parameters &parameters) {
  /* local rest frame velocity */
  FourVector velocity_lrf;
  velocity_lrf.set_x0(1.0);
  velocity_lrf.set_x1(particle->momentum().x1()
                      / particle->momentum().x0());
  velocity_lrf.set_x2(particle->momentum().x2()
                      / particle->momentum().x0());
  velocity_lrf.set_x3(particle->momentum().x3()
                      / particle->momentum().x0());

  /* The clock goes slower in the rest frame of the resonance */
  double inverse_gamma = sqrt(velocity_lrf.Dot(velocity_lrf));
  double resonance_frame_timestep = parameters.eps() * inverse_gamma;

  /* Exponential decay. Average lifetime t_avr = 1 / width
   * t / t_avr = width * t (remember GeV-fm conversion)
   * P(decay at Delta_t) = width * Delta_t
   * P(alive after n steps) = (1 - width * Delta_t)^n
   * = (1 - width * Delta_t)^(t / Delta_t)
   * -> exp(-width * t) when Delta_t -> 0
   */
  if (drand48() < resonance_frame_timestep * particle_type->width() / hbarc) {
    /* Time is up! Set the particle to decay at this timestep */
    particle->set_collision(2, 0.0, -1);
    decay_list->push_back(particle->id());
    return true;
  }
  return false;
}



/* colliding_particle - particle interaction */
size_t decay_particles(std::map<int, ParticleData> *particle,
  std::vector<ParticleType> *type, std::map<int, int> *map_type,
  std::list<int> *decay_list, size_t id_process, size_t *largest_id) {
  FourVector velocity_CM;

  for (std::list<int>::iterator id = decay_list->begin();
    id != decay_list->end(); ++id) {
    /* relevant particle id's for the collision */
    int id_a = *id;
    int interaction_type = (*particle)[id_a].process_type();

    if (interaction_type != 2)
      printf("Decays warning: ID %i (%s) has process type %i.\n",
           id_a, (*type)[(*map_type)[id_a]].name().c_str(), interaction_type);

    FourVector initial_momentum, final_momentum;
    initial_momentum += (*particle)[id_a].momentum();

    /* 1->2 resonance decay */
    printd("Process: Resonance decay. ");
    write_oscar((*particle)[id_a], (*type)[(*map_type)[id_a]], 1, 2);

    printd("Resonance momenta before decay: %g %g %g %g\n",
         (*particle)[id_a].momentum().x0(), (*particle)[id_a].momentum().x1(),
         (*particle)[id_a].momentum().x2(), (*particle)[id_a].momentum().x3());

    /* boost to rest frame */
    velocity_CM.set_x0(1.0);
    velocity_CM.set_x1((*particle)[id_a].momentum().x1()
                       / (*particle)[id_a].momentum().x0());
    velocity_CM.set_x2((*particle)[id_a].momentum().x2()
                       / (*particle)[id_a].momentum().x0());
    velocity_CM.set_x3((*particle)[id_a].momentum().x3()
                       / (*particle)[id_a].momentum().x0());
    (*particle)[id_a].set_momentum(
        (*particle)[id_a].momentum().LorentzBoost(velocity_CM));
    (*particle)[id_a].set_position(
        (*particle)[id_a].position().LorentzBoost(velocity_CM));

    printd("Boosted resonance momenta before decay: %g %g %g %g\n",
        (*particle)[id_a].momentum().x0(), (*particle)[id_a].momentum().x1(),
        (*particle)[id_a].momentum().x2(), (*particle)[id_a].momentum().x3());

    size_t id_new_a = resonance_decay(particle, type, map_type, &id_a,
                                        largest_id);
    size_t id_new_b = id_new_a + 1;

    printd("particle 1 momenta in lrf: %g %g %g %g\n",
           (*particle)[id_new_a].momentum().x0(),
           (*particle)[id_new_a].momentum().x1(),
           (*particle)[id_new_a].momentum().x2(),
           (*particle)[id_new_a].momentum().x3());
    printd("particle 2 momenta in lrf: %g %g %g %g\n",
           (*particle)[id_new_b].momentum().x0(),
           (*particle)[id_new_b].momentum().x1(),
           (*particle)[id_new_b].momentum().x2(),
           (*particle)[id_new_b].momentum().x3());

    boost_back_CM(&(*particle)[id_new_a], &(*particle)[id_new_b],
                  &velocity_CM);

    write_oscar((*particle)[id_new_a], (*type)[(*map_type)[id_new_a]]);
    write_oscar((*particle)[id_new_b], (*type)[(*map_type)[id_new_b]]);

    printd("particle 1 momenta in comp: %g %g %g %g\n",
           (*particle)[id_new_a].momentum().x0(),
           (*particle)[id_new_a].momentum().x1(),
           (*particle)[id_new_a].momentum().x2(),
           (*particle)[id_new_a].momentum().x3());
    printd("particle 2 momenta in comp: %g %g %g %g\n",
           (*particle)[id_new_b].momentum().x0(),
           (*particle)[id_new_b].momentum().x1(),
           (*particle)[id_new_b].momentum().x2(),
           (*particle)[id_new_b].momentum().x3());

    final_momentum += (*particle)[id_new_a].momentum();
    final_momentum += (*particle)[id_new_b].momentum();

    particle->erase(id_a);
    printd("ID %i has decayed and removed from the list.\n", id_a);

    /* unset collision time for both particles + keep id + unset partner */
    (*particle)[id_new_a].set_collision_past(id_process);
    (*particle)[id_new_b].set_collision_past(id_process);
    printd("Particle map has now %zu elements. \n", particle->size());

    id_process++;

    double numerical_tolerance = 1.0e-7;
    FourVector momentum_difference;
    momentum_difference += initial_momentum;
    momentum_difference -= final_momentum;
    if (fabs(momentum_difference.x0()) > numerical_tolerance) {
      printf("Process %lu type %i particle %s decay to %zu and %zu time %g\n",
        id_process, interaction_type, (*type)[(*map_type)[id_a]].name().c_str(),
             id_new_a, id_new_b, (*particle)[id_a].position().x0());
      printf("Warning: Interaction type %i E conservation violation %g\n",
             interaction_type, momentum_difference.x0());
    }
    if (fabs(momentum_difference.x1()) > numerical_tolerance)
      printf("Warning: Interaction type %i px conservation violation %g\n",
             interaction_type, momentum_difference.x1());
    if (fabs(momentum_difference.x2()) > numerical_tolerance)
      printf("Warning: Interaction type %i py conservation violation %g\n",
             interaction_type, momentum_difference.x2());
    if (fabs(momentum_difference.x3()) > numerical_tolerance)
      printf("Warning: Interaction type %i pz conservation violation %g\n",
             interaction_type, momentum_difference.x3());
  } /* end for std::list<int>::iterator id = decay_list->begin() */
  /* empty the decay table */
  decay_list->clear();
  printd("Decay list done.\n");

  /* return how many processes we have handled so far*/
  return id_process;
}
