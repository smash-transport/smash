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
#include <cmath>

#include "include/Parameters.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/constants.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/particles.h"

/* does_decay - does a resonance decay on this timestep? */
bool does_decay(std::vector<ParticleData> *particle,
  std::vector<ParticleType> *particle_type, std::map<int, int> *map_type,
  std::list<int> *collision_list, const Parameters &parameters, int id_res) {
  /* local rest frame velocity */
  FourVector velocity_lrf;
  velocity_lrf.set_x0(1.0);
  velocity_lrf.set_x1((*particle)[id_res].momentum().x1()
                      / (*particle)[id_res].momentum().x0());
  velocity_lrf.set_x2((*particle)[id_res].momentum().x2()
                      / (*particle)[id_res].momentum().x0());
  velocity_lrf.set_x3((*particle)[id_res].momentum().x3()
                      / (*particle)[id_res].momentum().x0());

  /* The clock goes slower in the rest frame of the resonance */
  double inverse_gamma = velocity_lrf.Dot(velocity_lrf);
  double resonance_frame_timestep = parameters.eps() * inverse_gamma;

  /* Exponential decay. Average lifetime t_avr = 1 / width
   * t / t_avr = width * t (remember GeV-fm conversion)
   * P(decay at Delta_t) = width * Delta_t
   * P(alive after n steps) = (1 - width * Delta_t)^n
   * = (1 - width * Delta_t)^(t / Delta_t)
   * -> exp(-width * t) when Delta_t -> 0
   */
  if (drand48() < resonance_frame_timestep
      * (*particle_type)[(*map_type)[id_res]].width() / hbarc) {
    /* Time is up! Set the particle to decay at this timestep */
    (*particle)[id_res].set_collision(2, 0.0, -1);
    collision_list->push_back(id_res);
    return true;
  } else {
    return false;
  }
}

/* collision_criteria_geometry - check by geometrical method if a collision
 *                               happens between particles
 */
void collision_criteria_geometry(std::vector<ParticleData> *particle,
  std::vector<ParticleType> *particle_type, std::map<int, int> *map_type,
  std::list<int> *collision_list, const Parameters &parameters, int id_a,
  int id_b) {
  /* just collided with this particle */
  if ((*particle)[id_a].id_process() >= 0
      && (*particle)[id_a].id_process() == (*particle)[id_b].id_process()) {
    printd("Skipping collided particle %d <-> %d at time %g due process %d\n",
        id_a, id_b, (*particle)[id_a].position().x0(),
        (*particle)[id_a].id_process());
    return;
  }

  /* Resonance production cross section */
  const double resonance_xsection = resonance_cross_section(
  &(*particle)[id_a], &(*particle)[id_b], &(*particle_type)[(*map_type)[id_a]],
  &(*particle_type)[(*map_type)[id_b]], particle_type);

  /* Total cross section is elastic + resonance production  */
  const double total_cross_section = parameters.cross_section()
    + resonance_xsection;

  {
    /* distance criteria according to cross_section */
    const double distance_squared = particle_distance(&(*particle)[id_a],
                                     &(*particle)[id_b]);
    if (distance_squared >= total_cross_section * fm2_mb * M_1_PI)
      return;
    printd("distance squared particle %d <-> %d: %g \n", id_a, id_b,
           distance_squared);
  }

  /* check according timestep: positive and smaller */
  const double time_collision = collision_time(&(*particle)[id_a],
    &(*particle)[id_b]);
  if (time_collision < 0.0 || time_collision >= parameters.eps())
    return;

  /* check for minimal collision time both particles */
  if (((*particle)[id_a].collision_time() > 0.0
        && time_collision > (*particle)[id_a].collision_time())
      || ((*particle)[id_b].collision_time() > 0
        && time_collision > (*particle)[id_b].collision_time())) {
    printd("%g Not minimal particle %d <-> %d\n",
        (*particle)[id_a].position().x0(), id_a, id_b);
    return;
  }

  /* handle minimal collision time of both particles */
  /* XXX: keep track of multiple possible collision partners */
  if (unlikely((*particle)[id_a].collision_time() > 0.0)) {
    int id_not = (*particle)[id_a].id_partner();
    printd("Not colliding particle %d <-> %d\n", id_a, id_not);
    /* unset collision partner to zero time and unexisting id */
    (*particle)[id_not].set_collision(-1, 0.0, -1);
    /* remove any of those partners from the list */
    if (id_a < id_not) {
      printd("Removing particle %d from collision list\n", id_a);
      collision_list->remove(id_a);
    } else {
      printd("Removing particle %d from collision list\n", id_not);
      collision_list->remove(id_not);
    }
  }
  if (unlikely((*particle)[id_b].collision_time() > 0.0)) {
    int id_not = (*particle)[id_b].id_partner();
    printd("Not colliding particle %d <-> %d\n", id_b, id_not);
    /* unset collision partner to zero time and unexisting id */
    (*particle)[id_not].set_collision(-1, 0.0, -1);
    /* remove any of those partners from the list */
    if (id_b < id_not) {
      printd("Removing particle %d from collision list\n", id_b);
      collision_list->remove(id_b);
    } else {
      printd("Removing particle %d from collision list\n", id_not);
      collision_list->remove(id_not);
    }
  }

  /* If resonance formation probability is high enough, do that,
   * otherwise do elastic collision
   */
  int interaction_type = 0;
  if (drand48() < resonance_xsection / total_cross_section)
    interaction_type = 1;

  /* setup collision partners */
  printd("collision time particle %d <-> %d: %g \n", id_a, id_b,
    time_collision);
  (*particle)[id_a].set_collision(interaction_type, time_collision, id_b);
  (*particle)[id_b].set_collision(interaction_type, time_collision, id_a);
  /* add to collision list */
  collision_list->push_back(id_a);
}

/* colliding_particle - particle interaction */
size_t collide_particles(std::vector<ParticleData> *particle,
  std::vector<ParticleType> *type, std::map<int, int> *map_type,
  std::list<int> *collision_list, size_t id_process) {
  FourVector velocity_CM;

  /* XXX: print debug output of collision list */
  /* collide: 2 <-> 2 soft momenta exchange */
  for (std::list<int>::iterator id = collision_list->begin();
    id != collision_list->end(); ++id) {
    /* relevant particle id's for the collision */
    int id_a = *id;
    int id_b = 0;
    int interaction_type = (*particle)[id_a].process_type();
    FourVector initial_momentum, final_momentum;
    initial_momentum += (*particle)[id_a].momentum();
    if (interaction_type < 2) {
      id_b = (*particle)[id_a].id_partner();
      initial_momentum += (*particle)[id_b].momentum();
      printd("Process %lu type %i particle %s<->%s colliding %d<->%d time %g\n",
        id_process, interaction_type, (*type)[(*map_type)[id_a]].name().c_str(),
             (*type)[(*map_type)[id_b]].name().c_str(), id_a, id_b,
             (*particle)[id_a].position().x0());
      write_oscar((*particle)[id_a], (*particle)[id_b],
                  (*type)[(*map_type)[id_a]], (*type)[(*map_type)[id_b]], 1);
      printd("particle 1 momenta before: %g %g %g %g\n",
          (*particle)[id_a].momentum().x0(), (*particle)[id_a].momentum().x1(),
          (*particle)[id_a].momentum().x2(), (*particle)[id_a].momentum().x3());
      printd("particle 2 momenta before: %g %g %g %g\n",
          (*particle)[id_b].momentum().x0(), (*particle)[id_b].momentum().x1(),
          (*particle)[id_b].momentum().x2(), (*particle)[id_b].momentum().x3());

      /* processes computed in the center of momenta */
      boost_CM(&(*particle)[id_a], &(*particle)[id_b], &velocity_CM);
    }

    if (interaction_type == 0) {
      /* 2->2 elastic scattering*/
      printd("Process: Elastic collision.\n");
      momenta_exchange(&(*particle)[id_a], &(*particle)[id_b],
       (*type)[(*map_type)[id_a]].mass(), (*type)[(*map_type)[id_b]].mass());

      boost_back_CM(&(*particle)[id_a], &(*particle)[id_b],
       &velocity_CM);
      write_oscar((*particle)[id_a], (*particle)[id_b],
       (*type)[(*map_type)[id_a]], (*type)[(*map_type)[id_b]], -1);
      printd("particle 1 momenta after: %g %g %g %g\n",
       (*particle)[id_a].momentum().x0(), (*particle)[id_a].momentum().x1(),
       (*particle)[id_a].momentum().x2(), (*particle)[id_a].momentum().x3());
      printd("particle 2 momenta after: %g %g %g %g\n",
       (*particle)[id_b].momentum().x0(), (*particle)[id_b].momentum().x1(),
       (*particle)[id_b].momentum().x2(), (*particle)[id_b].momentum().x3());

      final_momentum += (*particle)[id_a].momentum();
      final_momentum += (*particle)[id_b].momentum();

      /* unset collision time for both particles + keep id + unset partner */
      (*particle)[id_b].set_collision_past(id_process);

    } else if (interaction_type == 1) {
      /* 2->1 resonance formation */
      printd("Process: Resonance formation.\n");
      size_t id_new = resonance_formation(particle, type, map_type,
                                          &id_a, &id_b);
      /* Boost the new particle to computational frame */
      FourVector neg_velocity_CM;
      neg_velocity_CM.set_FourVector(1.0, -velocity_CM.x1(), -velocity_CM.x2(),
                                     -velocity_CM.x3());
      (*particle)[id_new].set_momentum(
          (*particle)[id_new].momentum().LorentzBoost(neg_velocity_CM));

      final_momentum += (*particle)[id_new].momentum();

      boost_back_CM(&(*particle)[id_a], &(*particle)[id_b],
       &velocity_CM);

      /* unset collision time for particles + keep id + unset partner */
      (*particle)[id_b].set_collision_past(id_process);

      /* The starting point of resonance is between the two initial particles */
      /* x_middle = x_a + (x_b - x_a) / 2 */
      FourVector middle_point;
      middle_point += (*particle)[id_b].position();
      middle_point -= (*particle)[id_a].position();
      middle_point /= 2.0;
      middle_point += (*particle)[id_a].position();
      (*particle)[id_new].set_position(middle_point);

      /* XXX: We need oscar output for 2->1 process:
       * write_oscar((*particle)[id_new],
       *  (*type)[(*map_type)[id_a]], (*type)[(*map_type)[id_b]], -1);
       */

      printd("Resonance %s with ID %lu \n",
       (*type)[(*map_type)[id_new]].name().c_str(), id_new);

      printd("has momentum in comp frame: %g %g %g %g\n",
      (*particle)[id_new].momentum().x0(), (*particle)[id_new].momentum().x1(),
      (*particle)[id_new].momentum().x2(), (*particle)[id_new].momentum().x3());

      printd("and position in comp frame: %g %g %g %g\n",
      (*particle)[id_new].position().x0(), (*particle)[id_new].position().x1(),
      (*particle)[id_new].position().x2(), (*particle)[id_new].position().x3());
    } else if (interaction_type == 2) {
      /* 1->2 resonance decay */
      printd("Process: Resonance decay.\n");
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

      size_t id_new_a = resonance_decay(particle, type, map_type, &id_a);
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
    } else {
       printf("Warning: Unspecified process type, nothing done.\n");
    }

    /* unset collision time for particles + keep id + unset partner */
    (*particle)[id_a].set_collision_past(id_process);
    id_process++;

    double numerical_tolerance = 1.0e-7;
    FourVector momentum_difference;
    momentum_difference += initial_momentum;
    momentum_difference -= final_momentum;
    if (fabs(momentum_difference.x0()) > numerical_tolerance)
      printf("Warning: Energy conservation violated by %g\n",
             momentum_difference.x0());
    if (fabs(momentum_difference.x1()) > numerical_tolerance)
      printf("Warning: x-momentum conservation violated by %g\n",
             momentum_difference.x1());
    if (fabs(momentum_difference.x2()) > numerical_tolerance)
      printf("Warning: y-momentum conservation violated by %g\n",
             momentum_difference.x2());
    if (fabs(momentum_difference.x3()) > numerical_tolerance)
      printf("Warning: z-momentum conservation violated by %g\n",
             momentum_difference.x3());
  }
  /* empty the collision table */
  collision_list->clear();
  printd("Collision list done.\n");

  /* return how many processes we have handled so far*/
  return id_process;
}
