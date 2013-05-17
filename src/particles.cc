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
#include "include/ParticleType.h"
#include "include/constants.h"
#include "include/box.h"
#include "include/outputroutines.h"

/* boost_COM - boost to center of momentum */
static void boost_COM(ParticleData *particle1, ParticleData *particle2,
  FourVector *velocity) {
  FourVector momentum1(particle1->momentum()), momentum2(particle2->momentum());
  FourVector position1(particle1->position()), position2(particle2->position());
  double cms_energy = momentum1.x0() + momentum2.x0();

  // CMS 4-velocity
  velocity->set_x0(1.0);
  velocity->set_x1((momentum1.x1() + momentum2.x1()) / cms_energy);
  velocity->set_x2((momentum1.x2() + momentum2.x2()) / cms_energy);
  velocity->set_x3((momentum1.x3() + momentum2.x3()) / cms_energy);

  // Boost the momenta into CMS frame
  momentum1 = momentum1.LorentzBoost(*velocity);
  momentum2 = momentum2.LorentzBoost(*velocity);

  // Boost the positions into CMS frame
  position1 = position1.LorentzBoost(*velocity);
  position2 = position2.LorentzBoost(*velocity);

  particle1->set_momentum(momentum1);
  particle1->set_position(position1);
  particle2->set_momentum(momentum2);
  particle2->set_position(position2);
}

/* boost_from_COM - boost back from center of momentum */
static void boost_from_COM(ParticleData *particle1, ParticleData *particle2,
  FourVector *velocity_orig) {
  FourVector momentum1(particle1->momentum()), momentum2(particle2->momentum());
  FourVector position1(particle1->position()), position2(particle2->position());
  FourVector velocity = *velocity_orig;

  /* To boost back set 1 + velocity */
  velocity *= -1;
  velocity.set_x0(1.0);

  /* Boost the momenta back to lab frame */
  momentum1 = momentum1.LorentzBoost(velocity);
  momentum2 = momentum2.LorentzBoost(velocity);

  /* Boost the positions back to lab frame */
  position1 = position1.LorentzBoost(velocity);
  position2 = position2.LorentzBoost(velocity);

  particle1->set_momentum(momentum1);
  particle1->set_position(position1);
  particle2->set_momentum(momentum2);
  particle2->set_position(position2);
}

/* particle_distance - measure distance between two particles */
static double particle_distance(ParticleData *particle_orig1,
  ParticleData *particle_orig2) {
  ParticleData particle1 = *particle_orig1, particle2 = *particle_orig2;
  FourVector velocity_com;
  double distance_squared;

  /* boost particles in center of momenta frame */
  boost_COM(&particle1, &particle2, &velocity_com);
  FourVector position_diff = particle1.position() - particle2.position();
  printd("Particle %d<->%d position diff: %g %g %g %g [fm]\n",
    particle1.id(), particle2.id(), position_diff.x0(), position_diff.x1(),
    position_diff.x2(), position_diff.x3());

  FourVector momentum_diff = particle1.momentum() - particle2.momentum();
  /* zero momentum leads to infite distance */
  if (momentum_diff.x1() == 0 || momentum_diff.x2() == 0
      || momentum_diff.x3() == 0)
    return  - position_diff.DotThree(position_diff);
  /* UrQMD squared distance criteria:
   * arXiv:nucl-th/9803035 (3.27): in center of momemtum frame
   * position of particle a: x_a
   * position of particle b: x_b
   * velocity of particle a: v_a
   * velocity of particle b: v_b
   * d^2_{coll} = (x_a - x_b)^2 - ((x_a - x_a) . (v_a - v_b))^2 / (v_a - v_b)^2
   */
  distance_squared = - position_diff.DotThree(position_diff)
    + position_diff.DotThree(momentum_diff)
      * position_diff.DotThree(momentum_diff)
      / momentum_diff.DotThree(momentum_diff);
  return distance_squared;
}

/* time_collision - measure collision time of two particles */
static double collision_time(ParticleData *particle1,
  ParticleData *particle2) {
  /* UrQMD collision time
   * arXiv:1203.4418 (5.15): in computational frame
   * position of particle a: x_a
   * position of particle b: x_b
   * momentum of particle a: p_a
   * momentum of particle b: p_b
   * t_{coll} = - (x_a - x_b) . (p_a - p_b) / (p_a - p_b)^2
   */
  FourVector position_diff = particle1->position() - particle2->position();
  FourVector velocity_diff = particle1->momentum() / particle1->momentum().x0()
    - particle2->momentum() / particle2->momentum().x0();
  return - position_diff.DotThree(velocity_diff)
           / velocity_diff.DotThree(velocity_diff);
}

/* momenta_exchange - soft scattering */
static void momenta_exchange(ParticleData *particle1, ParticleData *particle2,
  const float &particle1_mass, const float &particle2_mass) {
  double phi, theta;
  /* center of momentum hence equal for both particles */
  double momentum_radial = sqrt(particle1->momentum().x0()
    * particle1->momentum().x0() - particle1_mass * particle1_mass);
  printd("center of momenta 1: %g %g %g %g \n", particle1->momentum().x0(),
    particle1->momentum().x1(), particle1->momentum().x2(),
    particle1->momentum().x3());
  printd("center of momenta 2: %g %g %g %g \n", particle2->momentum().x0(),
    particle2->momentum().x1(), particle2->momentum().x2(),
    particle2->momentum().x3());

  /* particle exchange momenta and scatter to random direction */
  phi =  2 * M_PI * drand48();
  theta = M_PI * drand48();
  particle1->set_momentum(particle1_mass,
     momentum_radial * cos(phi) * sin(theta),
     momentum_radial * sin(phi) * sin(theta), momentum_radial * cos(theta));
  particle2->set_momentum(particle2_mass,
     -momentum_radial * cos(phi) * sin(theta),
     -momentum_radial * sin(phi) * sin(theta), -momentum_radial * cos(theta));

  printd("exchanged momenta 1: %g %g %g %g \n", particle1->momentum().x0(),
    particle1->momentum().x1(), particle1->momentum().x2(),
    particle1->momentum().x3());
  printd("exchanged momenta 2: %g %g %g %g \n", particle2->momentum().x0(),
    particle2->momentum().x1(), particle2->momentum().x2(),
    particle2->momentum().x3());
}

/* check_collision_criteria - check if a collision happens between particles */
void check_collision_criteria(std::vector<ParticleData> *particle,
  std::list<int> *collision_list, box box, int id, int id_other) {
  double distance_squared, time_collision;

  /* distance criteria according to cross_section */
  distance_squared = particle_distance(&(*particle)[id],
    &(*particle)[id_other]);
  if (distance_squared >= box.cross_section() * fm2_mb * M_1_PI)
    return;

  /* check according timestep: positive and smaller */
  time_collision = collision_time(&(*particle)[id], &(*particle)[id_other]);
  if (time_collision < 0 || time_collision >= box.eps())
    return;

  /* check for minimal collision time */
  if ((*particle)[id].collision_time() > 0
        && time_collision > (*particle)[id].collision_time()) {
    printd("%g Not minimal particle %d <-> %d\n",
        (*particle)[id].position().x0(), id, id_other);
    return;
  }

  /* just collided with this particle */
  if ((*particle)[id].collision_time() == 0
      && id_other == (*particle)[id].collision_id()) {
    printd("%g Skipping particle %d <-> %d\n",
        (*particle)[id].position().x0(), id, id_other);
    return;
  }

  /* handle minimal collision time */
  if (unlikely((*particle)[id].collision_time() > 0)) {
    int not_id = (*particle)[id].collision_id();
    printd("Not colliding particle %d <-> %d\n", id, not_id);
    /* unset collision partner to zero time and unexisting id */
    (*particle)[not_id].set_collision(0.0, -1);
    /* remove any of those partners from the list */
    collision_list->remove(id);
    collision_list->remove(not_id);
    /* XXX: keep track of multiple possible collision partners */
  }

  /* setup collision partners */
  printd("distance particle %d <-> %d: %g \n", id, id_other, distance_squared);
  printd("t_coll particle %d <-> %d: %g \n", id, id_other, time_collision);
  (*particle)[id].set_collision(time_collision, id_other);
  (*particle)[id_other].set_collision(time_collision, id);
  /* add to collision list */
  collision_list->push_back(id);
}

/* colliding_particle - particle interaction */
void collide_particles(std::vector<ParticleData> *particle, ParticleType *type,
  std::map<int, int> *map_type, std::list<int> *collision_list) {
  FourVector velocity_com;

  /* collide: 2 <-> 2 soft momenta exchange */
  for (std::list<int>::iterator id = collision_list->begin();
    id != collision_list->end(); ++id) {
    int id_other = (*particle)[*id].collision_id();
    printd("particle types %s<->%s colliding %d<->%d %g\n",
      type[(*map_type)[*id]].name().c_str(),
      type[(*map_type)[id_other]].name().c_str(), *id, id_other,
      (*particle)[*id].position().x0());
    write_oscar((*particle)[*id], (*particle)[id_other], type[(*map_type)[*id]],
      type[(*map_type)[id_other]], 1);
    printd("particle 1 momenta before: %g %g %g %g\n",
      (*particle)[*id].momentum().x0(), (*particle)[*id].momentum().x1(),
      (*particle)[*id].momentum().x2(), (*particle)[*id].momentum().x3());
    printd("particle 2 momenta before: %g %g %g %g\n",
      (*particle)[id_other].momentum().x0(),
      (*particle)[id_other].momentum().x1(),
      (*particle)[id_other].momentum().x2(),
      (*particle)[id_other].momentum().x3());

    /* exchange in center of momenta */
    boost_COM(&(*particle)[*id], &(*particle)[id_other], &velocity_com);
    momenta_exchange(&(*particle)[*id], &(*particle)[id_other],
      type[(*map_type)[*id]].mass(), type[(*map_type)[id_other]].mass());
    boost_from_COM(&(*particle)[*id], &(*particle)[id_other],
      &velocity_com);
    write_oscar((*particle)[*id], (*particle)[id_other], type[(*map_type)[*id]],
      type[(*map_type)[id_other]], -1);
    printd("particle 1 momenta after: %g %g %g %g\n",
      (*particle)[*id].momentum().x0(), (*particle)[*id].momentum().x1(),
      (*particle)[*id].momentum().x2(), (*particle)[*id].momentum().x3());
    printd("particle 2 momenta after: %g %g %g %g\n",
      (*particle)[id_other].momentum().x0(),
      (*particle)[id_other].momentum().x1(),
      (*particle)[id_other].momentum().x2(),
      (*particle)[id_other].momentum().x3());

    /* unset collision time for both particles + keep id */
    (*particle)[*id].set_collision_time(0.0);
    (*particle)[id_other].set_collision_time(0.0);
  }
  /* empty the collision table */
  collision_list->clear();
}
