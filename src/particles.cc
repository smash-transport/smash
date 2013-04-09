/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/particles.h"

#include <cstdio>
#include <list>

#include "include/FourVector.h"
#include "include/ParticleData.h"
#include "include/constants.h"
#include "include/box.h"
#include "include/outputroutines.h"

/* boost_COM - boost to center of momentum */
static void boost_COM(ParticleData *particle1, ParticleData *particle2,
  FourVector *velocity) {
  FourVector momentum1(particle1->momentum()), momentum2(particle2->momentum());
  FourVector position1(particle1->x()), position2(particle2->x());
  double cms_energy = momentum1.x0() + momentum2.x0();

  // CMS 4-velocity
  velocity->set_x0(1.0);
  velocity->set_x1((momentum1.x1() + momentum2.x1()) / cms_energy);
  velocity->set_x2((momentum1.x2() + momentum2.x2()) / cms_energy);
  velocity->set_x3((momentum1.x3() + momentum2.x3()) / cms_energy);

  // Boost the momenta into CMS frame
  momentum1 = momentum1.LorentzBoost(momentum1, *velocity);
  momentum2 = momentum2.LorentzBoost(momentum2, *velocity);

  // Boost the positions into CMS frame
  position1 = position1.LorentzBoost(position1, *velocity);
  position2 = position2.LorentzBoost(position2, *velocity);

  particle1->set_momentum(momentum1);
  particle1->set_position(position1);
  particle2->set_momentum(momentum2);
  particle2->set_position(position2);
}

/* boost_from_COM - boost back from center of momentum */
static void boost_from_COM(ParticleData *particle1, ParticleData *particle2,
  FourVector *velocity_orig) {
  FourVector momentum1(particle1->momentum()), momentum2(particle2->momentum());
  FourVector position1(particle1->x()), position2(particle2->x());
  FourVector velocity = *velocity_orig;

  /* To boost back set 1 + velocity */
  velocity *= -1;
  velocity.set_x0(1);

  /* Boost the momenta back to lab frame */
  momentum1 = momentum1.LorentzBoost(momentum1, velocity);
  momentum2 = momentum2.LorentzBoost(momentum2, velocity);

  /* Boost the positions back to lab frame */
  position1 = position1.LorentzBoost(position1, velocity);
  position2 = position2.LorentzBoost(position2, velocity);

  particle1->set_momentum(momentum1);
  particle1->set_position(position1);
  particle2->set_momentum(momentum2);
  particle2->set_position(position2);
}

/* particle_distance - measure distance between two particles */
static double particle_distance(ParticleData *particle_orig1,
  ParticleData *particle_orig2) {
  ParticleData particle1 = *particle_orig1, particle2 = *particle_orig2;
  FourVector velocity_com, position_diff, momentum_diff;
  double distance_squared;

  /* XXX: allow Kodama + UrQMD distance criteria */
  /* arXiv:nucl-th/9803035 (3.27): in center of momemtum frame
   * d^2_{coll} = (x1 - x2)^2 - ((x1 - x2) . (v1 - v2))^2 / (v1 - v2)^2
   */
  boost_COM(&particle1, &particle2, &velocity_com);
  position_diff = particle1.x() - particle2.x();
  momentum_diff = particle1.momentum() - particle2.momentum();
  distance_squared = - position_diff.DotThree(position_diff)
    + position_diff.DotThree(momentum_diff)
      * position_diff.DotThree(momentum_diff)
      / momentum_diff.DotThree(momentum_diff);
  return distance_squared;
}

/* time_collision - measure collision time of two particles */
static double collision_time(ParticleData *particle1,
  ParticleData *particle2) {
  FourVector position_diff, velocity_diff;
  double time;

  /* XXX: allow both Kodama + UrQMD distance criteria */
  /* arXiv:1203.4418 (5.15): t_{coll} = - (x1 - x2) . (v1 - v2) / (v1 - v2)^2 */
  position_diff = particle1->x() - particle2->x();
  velocity_diff = particle1->momentum() / particle1->momentum().x0()
    - particle2->momentum() / particle2->momentum().x0();
  time = - position_diff.DotThree(velocity_diff)
           / velocity_diff.DotThree(velocity_diff);
  return time;
}

/* momenta_exchange - soft scattering */
static void momenta_exchange(ParticleData *particle1, ParticleData *particle2) {
  FourVector momentum_copy;

  momentum_copy = particle1->momentum();
  particle1->set_momentum(particle2->momentum());
  particle2->set_momentum(momentum_copy);
}

/* check_collision - check if a collision can happen betwenn particles */
void check_collision(ParticleData *particle,
  std::list<ParticleData> *collision_list, box box, int id, int number) {
  double distance_squared, time_collision;

  /* check which particles interact:
   * This processes all particles above the certain id.
   */
  for (int i = id + 1; i < number; i++) {
    /* XXX: only check particles within nearest neighbour cells - size
     * according to cross_section */
    distance_squared = particle_distance(&particle[id], &particle[i]);

    /* particles are far apart */
    if (distance_squared >= box.cross_section() * fm2_mb / M_PI)
      continue;

    /* check according timestep: positive and smaller */
    time_collision = collision_time(&particle[id], &particle[i]);
    if (time_collision < 0 || time_collision >= box.eps())
      continue;

    /* check for minimal collision time */
    if (particle[id].collision_time() > 0
          && time_collision > particle[id].collision_time()) {
      printd("%g Not minimal particle %d <-> %d\n", particle[id].x().x0(), id,
        i);
      continue;
    }

    /* just collided with this particle */
    if (particle[id].collision_time() == 0
        && i == particle[id].collision_id()) {
      printd("%g Skipping particle %d <-> %d\n", particle[id].x().x0(), id,
        i);
      continue;
    }

    /* handle minimal collision time */
    if (unlikely(particle[id].collision_time() > 0)) {
      int not_id = particle[id].collision_id();
      printd("Not colliding particle %d <-> %d\n", id, not_id);
      /* unset collision partner */
      particle[not_id].set_collision(0, 0);
      /* remove any of those partners from the list */
      collision_list->remove(particle[id]);
      collision_list->remove(particle[not_id]);
    }

    /* setup collision partners */
    printd("distance particle %d <-> %d: %g \n", id, i, distance_squared);
    printd("t_coll particle %d <-> %d: %g \n", id, i, time_collision);
    particle[id].set_collision(time_collision, i);
    particle[i].set_collision(time_collision, id);
    /* add to collision list */
    collision_list->push_back(particle[id]);
  }
}

/* colliding_particle - particle interaction */
void collide_particles(ParticleData *particle,
  std::list<ParticleData> *collision_list) {
  FourVector velocity_com;

  /* collide: 2 <-> 2 soft momenta exchange */
  for (std::list<ParticleData>::iterator p = collision_list->begin();
    p != collision_list->end(); p++) {
    printd("particle id colliding %d with %d\n", p->id(), p->collision_id());
    write_oscar(particle[p->id()], particle[p->collision_id()], 1);

    /* exchange in center of momenta */
    boost_COM(&particle[p->id()], &particle[p->collision_id()], &velocity_com);
    momenta_exchange(&particle[p->id()], &particle[p->collision_id()]);
    boost_from_COM(&particle[p->id()], &particle[p->collision_id()],
      &velocity_com);
    write_oscar(particle[p->id()], particle[p->collision_id()], -1);

    /* unset collision time for both particles */
    particle[p->collision_id()].set_collision_time(0);
    particle[p->id()].set_collision_time(0);
  }
  /* empty the collision table */
  collision_list->clear();
}
