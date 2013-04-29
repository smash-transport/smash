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
#include <vector>

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
  ParticleData *particle_orig2, box box) {
  ParticleData particle1 = *particle_orig1, particle2 = *particle_orig2;
  FourVector velocity_com;
  double distance_squared;

  /* boost particles and box in center of momenta frame */
  FourVector box_com(particle1.x().x0(), box.a(), box.a(), box.a());
  boost_COM(&particle1, &particle2, &velocity_com);
  box_com = box_com.LorentzBoost(box_com, velocity_com);
  FourVector position_diff = particle1.x() - particle2.x();
  printd("Particle %d<->%d position diff: %g %g %g %g [fm]\n",
    particle1.id(), particle2.id(), position_diff.x0(), position_diff.x1(),
    position_diff.x2(), position_diff.x3());

  /* check for a periodic boundary cross */
  if (position_diff.x1() > box_com.x1() / 2)
    position_diff.set_x1(box_com.x1() - position_diff.x1());
  if (position_diff.x1() < -box_com.x1() / 2)
    position_diff.set_x1(box_com.x1() + position_diff.x1());
  if (position_diff.x2() > box_com.x2() / 2)
    position_diff.set_x2(box_com.x2() - position_diff.x2());
  if (position_diff.x2() < -box_com.x2() / 2)
    position_diff.set_x2(box_com.x2() + position_diff.x2());
  if (position_diff.x3() > box_com.x3() / 2)
    position_diff.set_x3(box_com.x3() - position_diff.x3());
  if (position_diff.x3() < -box_com.x3() / 2)
    position_diff.set_x3(box_com.x3() + position_diff.x3());
  printd("Particle %d<->%d boundary diff: %g %g %g %g [fm]\n",
    particle1.id(), particle2.id(), position_diff.x0(), position_diff.x1(),
    position_diff.x2(), position_diff.x3());

  FourVector momentum_diff = particle1.momentum() - particle2.momentum();
  /* zero momentum leads to infite distance */
  if (momentum_diff.x1() == 0 || momentum_diff.x2() == 0
      || momentum_diff.x3() == 0)
    return  - position_diff.DotThree(position_diff);
  /* UrQMD distance criteria:
   * arXiv:nucl-th/9803035 (3.27): in center of momemtum frame
   * d^2_{coll} = (x1 - x2)^2 - ((x1 - x2) . (v1 - v2))^2 / (v1 - v2)^2
   */
  distance_squared = - position_diff.DotThree(position_diff)
    + position_diff.DotThree(momentum_diff)
      * position_diff.DotThree(momentum_diff)
      / momentum_diff.DotThree(momentum_diff);
  return distance_squared;
}

/* time_collision - measure collision time of two particles */
static double collision_time(ParticleData *particle1,
  ParticleData *particle2, box box) {
  double time;

  /* UrQMD distance criteria
   * arXiv:1203.4418 (5.15): t_{coll} = - (x1 - x2) . (v1 - v2) / (v1 - v2)^2
   */
  FourVector position_diff = particle1->x() - particle2->x();
  /* check for a periodic boundary cross */
  if (position_diff.x1() > box.a() / 2)
    position_diff.set_x1(box.a() - position_diff.x1());
  if (position_diff.x2() > box.a() / 2)
    position_diff.set_x2(box.a() - position_diff.x2());
  if (position_diff.x3() > box.a() / 2)
    position_diff.set_x3(box.a() - position_diff.x3());
  FourVector velocity_diff = particle1->momentum() / particle1->momentum().x0()
    - particle2->momentum() / particle2->momentum().x0();
  time = - position_diff.DotThree(velocity_diff)
           / velocity_diff.DotThree(velocity_diff);
  return time;
}

/* momenta_exchange - soft scattering */
static void momenta_exchange(ParticleData *particle1, ParticleData *particle2) {
  FourVector momentum_copy = particle1->momentum();
  particle1->set_momentum(particle2->momentum());
  particle2->set_momentum(momentum_copy);
}

/* check_collision_criteria - check if a collision happens between particles */
static void check_collision_criteria(ParticleData *particle,
  std::list<int> *collision_list, box box, int id, int id_other) {
  double distance_squared, time_collision;

  /* criteria according to cross_section */
  distance_squared = particle_distance(&particle[id], &particle[id_other],
    box);

  /* particles are far apart */
  if (distance_squared >= box.cross_section() * fm2_mb / M_PI)
    return;

  /* check according timestep: positive and smaller */
  time_collision = collision_time(&particle[id], &particle[id_other], box);
  if (time_collision < 0 || time_collision >= box.eps())
    return;

  /* check for minimal collision time */
  if (particle[id].collision_time() > 0
        && time_collision > particle[id].collision_time()) {
    printd("%g Not minimal particle %d <-> %d\n", particle[id].x().x0(), id,
        id_other);
    return;
  }

  /* just collided with this particle */
  if (particle[id].collision_time() == 0
      && id_other == particle[id].collision_id()) {
    printd("%g Skipping particle %d <-> %d\n", particle[id].x().x0(), id,
        id_other);
    return;
  }

  /* handle minimal collision time */
  if (unlikely(particle[id].collision_time() > 0)) {
    int not_id = particle[id].collision_id();
    printd("Not colliding particle %d <-> %d\n", id, not_id);
    /* unset collision partner to zero time and unexisting id */
    particle[not_id].set_collision(0, -1);
    /* remove any of those partners from the list */
    collision_list->remove(id);
    collision_list->remove(not_id);
    /* XXX: keep track of multiple possible collision partners */
  }

  /* setup collision partners */
  printd("distance particle %d <-> %d: %g \n", id, id_other, distance_squared);
  printd("t_coll particle %d <-> %d: %g \n", id, id_other, time_collision);
  particle[id].set_collision(time_collision, id_other);
  particle[id_other].set_collision(time_collision, id);
  /* add to collision list */
  collision_list->push_back(id);
}

/* check_collision - check if a collision can happen betwenn particles */
void check_collision(ParticleData *particle,
  std::list<int> *collision_list, box box, int number) {
  std::vector<std::vector<std::vector<std::vector<int> > > > grid;
  int N;
  int x, y, z;

  /* calculate approximate grid size according to double interaction length */
  N = round(box.a() / sqrt(box.cross_section() * fm2_mb / M_PI) / 2);

  /* For small boxes no point in splitting up in grids */
  if (unlikely(N < 4 || number < 10)) {
    for (int id = 0; id < number - 1; id++)
      for (int id_other = id + 1; id_other < number; id_other++)
        check_collision_criteria(particle, collision_list, box, id, id_other);
    return;
  }

  /* allocate grid */
  grid.resize(N);
  for (int i = 0; i < N; i++) {
    grid[i].resize(N);
    for (int j = 0; j < N; j++)
      grid[i][j].resize(N);
  }

  /* populate grid */
  for (int id = 0; id < number; id++) {
    /* XXX: function - map particle position to grid number */
    z = round(particle[id].x().x1() / box.a() * (N - 1));
    x = round(particle[id].x().x2() / box.a() * (N - 1));
    y = round(particle[id].x().x3() / box.a() * (N - 1));
    printd_position(particle[id]);
    printd("grid cell %i: %i %i %i\n", N, z, x, y);
    grid[z][x][y].push_back(id);
  }

  /* semi optimised nearest neighbour search:
   * http://en.wikipedia.org/wiki/Cell_lists
   */
  for (int id = 0; id < number - 1; id++) {
    /* XXX: function - map particle position to grid number */
    z = round(particle[id].x().x1() / box.a() * (N - 1));
    x = round(particle[id].x().x2() / box.a() * (N - 1));
    y = round(particle[id].x().x3() / box.a() * (N - 1));
    /* check all neighbour grids */
    int sz, sy, sx;
    for (int cz = 0; cz < 3; cz++) {
      if (cz == 0)
        sz = sm(z, N);
      else if (cz == 2)
        sz = sp(z, N);
      else
        sz = z;
      for (int cy = 0; cy < 3; cy++) {
        if (cy == 0)
          sy = sm(y, N);
        else if (cy == 2)
          sy = sp(y, N);
        else
          sy = y;
        for (int cx = 0; cx < 3; cx++) {
          if (cx == 0)
            sx = sm(x, N);
          else if (cx == 2)
            sx = sp(x, N);
          else
            sx = x;
          /* empty grid cell */
          if (grid[sz][sx][sy].empty())
            continue;
          /* grid cell particle list */
          for (std::vector<int>::iterator id_other = grid[sz][sx][sy].begin();
               id_other != grid[sz][sx][sy].end(); ++id_other) {
	    /* only check against particles above current id
	     * to avoid double counting
	     */
            if (*id_other <= id)
              continue;

            printd("grid cell particle %i <-> %i\n", id, *id_other);
            check_collision_criteria(particle, collision_list, box, id,
              *id_other);
          }
        }
      }
    }
  }
}

/* colliding_particle - particle interaction */
void collide_particles(ParticleData *particle,
  std::list<int> *collision_list) {
  FourVector velocity_com;

  /* collide: 2 <-> 2 soft momenta exchange */
  for (std::list<int>::iterator id = collision_list->begin();
    id != collision_list->end(); ++id) {
    int id_other = particle[*id].collision_id();
    printd("particle id colliding %d with %d\n", *id, id_other);
    write_oscar(particle[*id], particle[id_other], 1);

    /* exchange in center of momenta */
    boost_COM(&particle[*id], &particle[id_other], &velocity_com);
    momenta_exchange(&particle[*id], &particle[id_other]);
    boost_from_COM(&particle[*id], &particle[id_other],
      &velocity_com);
    write_oscar(particle[*id], particle[id_other], -1);

    /* unset collision time for both particles + keep id */
    particle[*id].set_collision_time(0);
    particle[id_other].set_collision_time(0);
  }
  /* empty the collision table */
  collision_list->clear();
}
