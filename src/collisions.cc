/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/collisions.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <list>
#include <map>
#include <utility>
#include <vector>

#include "include/CollisionData.h"
#include "include/FourVector.h"
#include "include/Laboratory.h"
#include "include/Particles.h"
#include "include/ParticleData.h"
#include "include/constants.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/resonances.h"

/* collision_criteria_geometry - check by geometrical method if a collision
 *                               happens between particles
 */
void collision_criteria_geometry(Particles *particles,
  CrossSections *cross_sections,
  std::list<int> *collision_list, const Laboratory &parameters, int id_a,
  int id_b, size_t *rejection_conflict) {
  /* just collided with this particle */
  if (particles->data(id_a).id_process() >= 0
      && particles->data(id_a).id_process()
      == particles->data(id_b).id_process()) {
    printd("Skipping collided particle %d <-> %d at time %g due process %d\n",
        id_a, id_b, particles->data(id_a).position().x0(),
        particles->data(id_a).id_process());
    return;
  }

  /* Compute kinematic quantities needed for cross section calculations  */
  cross_sections->compute_kinematics(particles, id_a, id_b);
  /* Resonance production cross section */
  std::vector< std::pair<std::vector<int>, double> > resonance_xsections
    = resonance_cross_section(particles->data(id_a), particles->data(id_b),
      particles->type(id_a), particles->type(id_b), particles);

  /* Total cross section is elastic + resonance production  */
  /* (Ignore annihilation and total for now) */
  const double total_cross_section
    = cross_sections->elastic(particles, id_a, id_b)
      + resonance_xsections.at(0).second;

  {
    /* distance criteria according to cross_section */
    const double distance_squared = particle_distance(
               particles->data_pointer(id_a), particles->data_pointer(id_b));
    if (distance_squared >= total_cross_section * fm2_mb * M_1_PI)
      return;
    printd("distance squared particle %d <-> %d: %g \n", id_a, id_b,
           distance_squared);
  }

  /* check according timestep: positive and smaller */
  const double time_collision = collision_time(particles->data(id_a),
    particles->data(id_b));
  if (time_collision < 0.0 || time_collision >= parameters.eps())
    return;

  /* check for minimal collision time both particles */
  if ((particles->data(id_a).collision_time() > 0.0
        && time_collision > particles->data(id_a).collision_time())
      || (particles->data(id_b).collision_time() > 0
        && time_collision > particles->data(id_b).collision_time())) {
    printd("%g Not minimal particle %d <-> %d\n",
        particles->data(id_a).position().x0(), id_a, id_b);
    return;
  }

  /* handle minimal collision time of both particles */
  if (unlikely(particles->data(id_a).collision_time() > 0.0)) {
    int id_not = particles->data(id_a).id_partner();
    printd("Not colliding particle %d <-> %d\n", id_a, id_not);
    /* unset collision partner to zero time and unexisting id */
    if (particles->count(id_not) > 0)
      particles->data_pointer(id_not)->set_collision(-1, 0.0, -1);
    /* remove any of those partners from the list */
    if (id_a < id_not) {
      printd("Removing particle %d from collision list\n", id_a);
      collision_list->remove(id_a);
    } else {
      printd("Removing particle %d from collision list\n", id_not);
      collision_list->remove(id_not);
    }
    /* collect statistics of multiple possible collision partner */
    (*rejection_conflict)++;
  }
  if (unlikely(particles->data(id_b).collision_time() > 0.0)) {
    int id_not = particles->data(id_b).id_partner();
    printd("Not colliding particle %d <-> %d\n", id_b, id_not);
    /* unset collision partner to zero time and unexisting id */
    if (particles->count(id_not) > 0)
      particles->data_pointer(id_not)->set_collision(-1, 0.0, -1);
    /* remove any of those partners from the list */
    if (id_b < id_not) {
      printd("Removing particle %d from collision list\n", id_b);
      collision_list->remove(id_b);
    } else {
      printd("Removing particle %d from collision list\n", id_not);
      collision_list->remove(id_not);
    }
    /* collect statistics of multiple possible collision partner */
    (*rejection_conflict)++;
  }

  /* If resonance formation probability is high enough, do that,
   * otherwise do elastic collision
   */
  int interaction_type = 0;
  std::vector<int> final_particles;
  if ((resonance_xsections.at(0)).second > really_small) {
    double random_interaction = drand48();
    double interaction_probability = 0.0;
    std::vector< std::pair<std::vector<int>, double> >::iterator resonances
      = resonance_xsections.begin();
    while (interaction_type == 0 && resonances != resonance_xsections.end()) {
      if ((resonances->first).size() > 1 || (resonances->first).at(0) != 0) {
        interaction_probability += resonances->second / total_cross_section;
        if (random_interaction < interaction_probability) {
          interaction_type = 1;
          final_particles = resonances->first;
        }
      }
      ++resonances;
    }
  }

  /* setup collision partners */
  printd("collision type %d particle %d <-> %d time: %g\n", interaction_type,
     id_a, id_b, time_collision);
  particles->data_pointer(id_a)->set_collision(interaction_type,
                                   time_collision, id_b, final_particles);
  particles->data_pointer(id_b)->set_collision(interaction_type,
                                   time_collision, id_a, final_particles);
  printd("collision type %d particle %d <-> %d time: %g\n",
         particles->data(id_a).process_type(),
         particles->data(id_a).id_partner(),
         particles->data(id_b).id_partner(),
         particles->data(id_a).collision_time());
  /* add to collision list */
  collision_list->push_back(id_a);
  resonance_xsections.clear();
}

/* colliding_particle - particle interaction */
size_t collide_particles(Particles *particles, std::list<int> *collision_list,
  size_t id_process, int *resonance_formations) {
  FourVector velocity_CM;

  /* XXX: print debug output of collision list */
  /* collide: 2 <-> 2 soft momenta exchange */
  for (std::list<int>::iterator id = collision_list->begin();
    id != collision_list->end(); ++id) {
    /* relevant particle id's for the collision */
    int id_a = *id;
    int id_b = particles->data(id_a).id_partner();
    int interaction_type = particles->data(id_a).process_type();
    FourVector initial_momentum(particles->data(id_a).momentum()
      + particles->data(id_b).momentum());
    FourVector final_momentum;

    printd("Process %zu type %i particle %s<->%s colliding %d<->%d time %g\n",
      id_process, interaction_type, particles->type(id_a).name().c_str(),
           particles->type(id_a).name().c_str(), id_a, id_b,
           particles->data(id_a).position().x0());
    printd_momenta("particle 1 momenta before", particles->data(id_a));
    printd_momenta("particle 2 momenta before", particles->data(id_b));

    /* 2->2 elastic scattering */
    if (interaction_type == 0) {
      printd("Process: Elastic collision.\n");
      write_oscar(particles->data(id_a), particles->type(id_a), 2, 2);
      write_oscar(particles->data(id_b), particles->type(id_b));

      /* processes computed in the center of momenta */
      boost_CM(particles->data_pointer(id_a), particles->data_pointer(id_b),
               &velocity_CM);
      momenta_exchange(particles->data_pointer(id_a),
                       particles->data_pointer(id_b));
      boost_back_CM(particles->data_pointer(id_a),
                    particles->data_pointer(id_b), &velocity_CM);

      write_oscar(particles->data(id_a), particles->type(id_a));
      write_oscar(particles->data(id_b), particles->type(id_b));

      printd_momenta("particle 1 momenta after", particles->data(id_a));
      printd_momenta("particle 2 momenta after", particles->data(id_b));

      final_momentum = particles->data(id_a).momentum()
      + particles->data(id_b).momentum();

      /* unset collision time for both particles + keep id + unset partner */
      particles->data_pointer(id_a)->set_collision_past(id_process);
      particles->data_pointer(id_b)->set_collision_past(id_process);

    } else if (interaction_type == 1) {
      /* resonance formation */
      (*resonance_formations)++;
      printd("Process: Resonance formation. ");
      size_t new_particles = (particles->data(id_a).final_state()).size();
      write_oscar(particles->data(id_a), particles->type(id_a), 2,
                  new_particles);
      write_oscar(particles->data(id_b), particles->type(id_b));
      /* processes computed in the center of momenta */
      boost_CM(particles->data_pointer(id_a), particles->data_pointer(id_b),
               &velocity_CM);

      size_t id_new = resonance_formation(particles, id_a, id_b,
        particles->data(id_a).final_state());

      boost_back_CM(particles->data_pointer(id_a),
                    particles->data_pointer(id_b), &velocity_CM);

      /* Boost the new particle to computational frame */
      FourVector neg_velocity_CM;
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
        write_oscar(particles->data(id_value), particles->type(id_value));
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
    } else {
      printf("Warning: ID %i (%s) has unspecified process type %i.\n",
             id_a, particles->type(id_a).name().c_str(), interaction_type);
    } /* end if (interaction_type == 0) */
    id_process++;

    FourVector momentum_difference;
    momentum_difference += initial_momentum;
    momentum_difference -= final_momentum;
    if (fabs(momentum_difference.x0()) > really_small) {
      printf("Process %zu type %i particle %s<->%s colliding %d<->%d time %g\n",
        id_process, interaction_type, particles->type(id_a).name().c_str(),
             particles->type(id_b).name().c_str(), id_a, id_b,
             particles->data(id_a).position().x0());
      printf("Warning: Interaction type %i E conservation violation %g\n",
             interaction_type, momentum_difference.x0());
    }
    if (fabs(momentum_difference.x1()) > really_small)
      printf("Warning: Interaction type %i px conservation violation %g\n",
             interaction_type, momentum_difference.x1());
    if (fabs(momentum_difference.x2()) > really_small)
      printf("Warning: Interaction type %i py conservation violation %g\n",
             interaction_type, momentum_difference.x2());
    if (fabs(momentum_difference.x3()) > really_small)
      printf("Warning: Interaction type %i pz conservation violation %g\n",
             interaction_type, momentum_difference.x3());
  } /* end for std::list<int>::iterator id = collision_list->begin() */
  /* empty the collision table */
  collision_list->clear();
  printd("Collision list done.\n");

  /* return how many processes we have handled so far*/
  return id_process;
}
