/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionsfinder.h"

#include "include/action.h"
#include "include/constants.h"
#include "include/crosssections.h"
#include "include/experimentparameters.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/particles.h"
#include "include/resonances.h"

namespace Smash {


double ScatterActionsFinder::collision_time(const ParticleData &p1,
                                            const ParticleData &p2) {
  /* UrQMD collision time
   * arXiv:1203.4418 (5.15): in computational frame
   * position of particle a: x_a
   * position of particle b: x_b
   * momentum of particle a: p_a
   * momentum of particle b: p_b
   * t_{coll} = - (x_a - x_b) . (p_a - p_b) / (p_a - p_b)^2
   */
  ThreeVector pos_diff = p1.position().threevec() - p2.position().threevec();
  printd("Particle %d<->%d position difference: %g %g %g %g [fm]\n",
    p1.id(), p2.id(), pos_diff.x1(), pos_diff.x2(), pos_diff.x3());
  ThreeVector velo_diff = p1.velocity() - p2.velocity();
  printd("Particle %d<->%d velocity difference: %g %g %g %g [fm]\n",
    p1.id(), p2.id(), velo_diff.x1(), velo_diff.x2(), velo_diff.x3());
  /* Zero momentum leads to infite distance, particles are not approaching. */
  if (fabs(velo_diff.sqr()) < really_small) {
    return -1.0;
  } else {
    return -pos_diff * velo_diff/velo_diff.sqr();
  }
}


ActionPtr
ScatterActionsFinder::check_collision(const int id_a, const int id_b, Particles *particles,
                                      const ExperimentParameters &parameters,
                                      CrossSections *cross_sections) const {

  ScatterAction* act = nullptr;

  const ParticleData data_a = particles->data(id_a);
  const ParticleData data_b = particles->data(id_b);

  /* just collided with this particle */
  if (data_a.id_process() >= 0 && data_a.id_process() == data_b.id_process()) {
    printd("Skipping collided particle %d <-> %d at time %g due process %d\n",
           id_a, id_b, data_a.position().x0(), data_a.id_process());
    return nullptr;
  }

  /* check according timestep: positive and smaller */
  const double time_until_collision = collision_time(data_a, data_b);
  if (time_until_collision < 0.0 ||
      time_until_collision >= parameters.timestep_duration()) {
    return nullptr;
  }

  /* check for minimal collision time both particles */
  if ((data_a.collision_time() > 0.0 &&
       time_until_collision > data_a.collision_time()) ||
      (data_b.collision_time() > 0.0 &&
       time_until_collision > data_b.collision_time())) {
    printd("%g Not minimal particle %d <-> %d\n", data_a.position().x0(), id_a,
           id_b);
    return nullptr;
  }

  act = new ScatterAction(data_a, data_b, time_until_collision);

  /* Resonance production cross section */
  ProcessBranchList resonance_xsections = resonance_cross_section(data_a,
                                                                  data_b);
  act->add_processes(resonance_xsections);

  /* Add elastic process.  */
  act->add_process(ProcessBranch(data_a.pdgcode(), data_b.pdgcode(),
                                 cross_sections->elastic(data_a, data_b)));

  {
    /* distance criteria according to cross_section */
    const double distance_squared = act->particle_distance();
    if (distance_squared >= act->weight() * fm2_mb * M_1_PI) {
        delete act;
        return nullptr;
      }
    printd("distance squared particle %d <-> %d: %g \n", id_a, id_b,
           distance_squared);
  }

  /* Set up collision partners. */
  particles->data(id_a).set_collision_time(time_until_collision);
  particles->data(id_b).set_collision_time(time_until_collision);

  return ActionPtr(act);
}

std::vector<ActionPtr> ScatterActionsFinder::find_possible_actions(
    Particles *particles, const ExperimentParameters &parameters,
    CrossSections *cross_sections) const {
  std::vector<ActionPtr> actions;
  double neighborhood_radius_squared = parameters.cross_section * fm2_mb * M_1_PI * 4;

  for (const auto &p1 : particles->data()) {
    for (const auto &p2 : particles->data()) {

      int id_a = p1.id(), id_b = p2.id();

      /* Check for same particle and double counting. */
      if (id_a >= id_b) continue;

      /* Skip particles that are double interaction radius length away
       * (3-product gives negative values
       * with the chosen sign convention for the metric). */
      FourVector distance = p1.position() - p2.position();
      if (distance.sqr3() > neighborhood_radius_squared)
        continue;

      /* Check if collision is possible. */
      ActionPtr act = check_collision (id_a, id_b, particles, parameters, cross_sections);

      /* Add to collision list. */
      if (act != nullptr) {
        actions.push_back (std::move(act));
      }
    }
  }
  return std::move(actions);
}



#if 0
GridScatterFinder::GridScatterFinder(float length) : length_(length) {
}


void
GridScatterFinder::find_possible_actions (std::vector<ActionPtr> &actions,
        Particles *particles, const ExperimentParameters &parameters,
        CrossSections *cross_sections) const {

  std::vector<std::vector<std::vector<std::vector<int> > > > grid;
  int N;
  int x, y, z;

  /* For small boxes no point in splitting up in grids */
  /* calculate approximate grid size according to double interaction length */
  N = round(length_ / sqrt(parameters.cross_section * fm2_mb * M_1_PI) * 0.5);
  if (unlikely(N < 4 || particles->size() < 10)) {
    /* XXX: apply periodic boundary condition */
    ScatterActionsFinder::find_possible_actions (actions, particles, parameters,
                                                 cross_sections);
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
  for (const ParticleData &data : particles->data()) {
    /* XXX: function - map particle position to grid number */
    x = round(data.position().x1() / length_ * (N - 1));
    y = round(data.position().x2() / length_ * (N - 1));
    z = round(data.position().x3() / length_ * (N - 1));
    printd_position(data);
    printd("grid cell particle %i: %i %i %i of %i\n", data.id(), x, y, z, N);
    if (unlikely(x >= N || y >= N || z >= N))
      printf("W: Particle outside the box: %g %g %g \n",
             data.position().x1(), data.position().x2(),
             data.position().x3());
    grid[x][y][z].push_back(data.id());
  }
  /* semi optimised nearest neighbour search:
   * http://en.wikipedia.org/wiki/Cell_lists
   */
  FourVector shift;
  for (const ParticleData &data : particles->data()) {
    /* XXX: function - map particle position to grid number */
    x = round(data.position().x1() / length_ * (N - 1));
    y = round(data.position().x2() / length_ * (N - 1));
    z = round(data.position().x3() / length_ * (N - 1));
    if (unlikely(x >= N || y >= N || z >= N))
      printf("grid cell particle %i: %i %i %i of %i\n", data.id(), x, y, z, N);
    /* check all neighbour grids */
    for (int cx = -1; cx < 2; cx++) {
      int sx = cx + x;
      /* apply periodic boundary condition for particle positions */
      if (sx < 0) {
        sx = N - 1;
        shift.set_x1(-length_);
      } else if (sx > N - 1) {
        sx = 0;
        shift.set_x1(length_);
      } else {
        shift.set_x1(0);
      }
      for (int cy = -1; cy < 2; cy++) {
        int sy = cy + y;
        if (sy < 0) {
          sy = N - 1;
          shift.set_x2(-length_);
        } else if (sy > N - 1) {
          sy = 0;
          shift.set_x2(length_);
        } else {
          shift.set_x2(0);
        }
        for (int cz = -1; cz < 2; cz++) {
          int sz = cz + z;
          if (sz < 0) {
            sz = N - 1;
            shift.set_x3(-length_);
          } else if (sz > N - 1) {
            sz = 0;
            shift.set_x3(length_);
          } else {
            shift.set_x3(0);
          }
          /* empty grid cell */
          if (grid[sx][sy][sz].empty()) {
            continue;
          }
          /* grid cell particle list */
          for (auto id_other = grid[sx][sy][sz].begin();
               id_other != grid[sx][sy][sz].end(); ++id_other) {
            /* only check against particles above current id
             * to avoid double counting
             */
            if (*id_other <= data.id()) {
              continue;
            }
            printd("grid cell particle %i <-> %i\n", data.id(), *id_other);
            ActionPtr act;
            if (shift == 0)
              act = check_collision (data.id(), *id_other, particles, parameters, cross_sections);
            else {
              /* apply eventual boundary before and restore after */
              particles->data(*id_other)
                  .set_position(particles->data(*id_other).position() + shift);
              act = check_collision (data.id(), *id_other, particles, parameters, cross_sections);
              particles->data(*id_other)
                  .set_position(particles->data(*id_other).position() - shift);
            }
            /* Add to collision list. */
            if (act != nullptr) {
              actions.push_back (std::move(act));
            }
          } /* grid particles loop */
        }   /* grid sz */
      }     /* grid sy */
    }       /* grid sx */
  }         /* outer particle loop */
}
#endif

}  // namespace Smash
