/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/scatteractionsfinder.h"

#include <algorithm>    // std::sort

#include "include/resonances.h"
#include "include/macros.h"
#include "include/outputroutines.h"

namespace Smash {

ActionPtr
ScatterActionsFinder::check_collision(int id_a, int id_b, Particles *particles,
				      const ExperimentParameters &parameters,
				      CrossSections *cross_sections)
      const
{
  ActionPtr act = nullptr;
  std::vector<int> in_part;

  /* just collided with this particle */
  if (particles->data(id_a).id_process() >= 0
      && particles->data(id_a).id_process()
      == particles->data(id_b).id_process()) {
    printd("Skipping collided particle %d <-> %d at time %g due process %d\n",
	id_a, id_b, particles->data(id_a).position().x0(),
	particles->data(id_a).id_process());
    return nullptr;
  }

  /* check according timestep: positive and smaller */
  const double time_collision = collision_time(particles->data(id_a),
    particles->data(id_b));
  if (time_collision < 0.0 || time_collision >= parameters.eps)
    return nullptr;

  /* check for minimal collision time both particles */
  if ((particles->data(id_a).collision_time() > 0.0
	&& time_collision > particles->data(id_a).collision_time())
      || (particles->data(id_b).collision_time() > 0.0
	&& time_collision > particles->data(id_b).collision_time())) {
    printd("%g Not minimal particle %d <-> %d\n",
	particles->data(id_a).position().x0(), id_a, id_b);
    return nullptr;
  }

  in_part.push_back(id_a);
  in_part.push_back(id_b);
  act = ActionPtr(new ScatterAction(in_part, time_collision));

  /* Compute kinematic quantities needed for cross section calculations  */
  cross_sections->compute_kinematics(particles, id_a, id_b);

  /* Resonance production cross section */
  std::vector<ProcessBranch> resonance_xsections
    = resonance_cross_section(particles->data(id_a), particles->data(id_b),
      particles->type(id_a), particles->type(id_b), particles);
  act->add_processes(resonance_xsections);

  /* Add elastic process.  */
  act->add_process(ProcessBranch(particles->data(id_a).pdgcode(),
				  particles->data(id_b).pdgcode(),
				  cross_sections->elastic(particles,
							  id_a, id_b), 0));

  {
    /* distance criteria according to cross_section */
    const double distance_squared = particle_distance(
	      particles->data_pointer(id_a), particles->data_pointer(id_b));
    if (distance_squared >= act->weight() * fm2_mb * M_1_PI)
      {
	act.reset(nullptr);
	return nullptr;
      }
    printd("distance squared particle %d <-> %d: %g \n", id_a, id_b,
	  distance_squared);
  }

  /* TODO: handle minimal collision time of both particles */
//       if (unlikely(particles->data(id_a).collision_time() > 0.0)) {
// 	int id_not = particles->data(id_a).id_partner();
// 	printd("Not colliding particle %d <-> %d\n", id_a, id_not);
// 	/* unset collision partner to zero time and unexisting id */
// 	if (particles->count(id_not) > 0)
// 	  particles->data_pointer(id_not)->set_collision(-1, 0.0, -1);
// 	/* remove any of those partners from the list */
// 	if (id_a < id_not) {
// 	  printd("Removing particle %d from collision list\n", id_a);
// 	  collision_list->remove(id_a);
// 	} else {
// 	  printd("Removing particle %d from collision list\n", id_not);
// 	  collision_list->remove(id_not);
// 	}
// 	/* collect statistics of multiple possible collision partner */
// 	(*rejection_conflict)++;
//       }
//       if (unlikely(particles->data(id_b).collision_time() > 0.0)) {
// 	int id_not = particles->data(id_b).id_partner();
// 	printd("Not colliding particle %d <-> %d\n", id_b, id_not);
//     /* unset collision partner to zero time and unexisting id */
// 	if (particles->count(id_not) > 0)
// 	  particles->data_pointer(id_not)->set_collision(-1, 0.0, -1);
// 	/* remove any of those partners from the list */
// 	if (id_b < id_not) {
// 	  printd("Removing particle %d from collision list\n", id_b);
// 	  collision_list->remove(id_b);
// 	} else {
// 	  printd("Removing particle %d from collision list\n", id_not);
// 	  collision_list->remove(id_not);
// 	}
// 	/* collect statistics of multiple possible collision partner */
// 	(*rejection_conflict)++;
//       }

  /* Decide for a particular final state. */
  act->decide();

  /* Set up collision partners. */
  printd("collision type %d particle %d <-> %d time: %g\n",
	 act->process_type(), id_a, id_b, time_collision);
  particles->data_pointer(id_a)->set_collision_time(time_collision);
  particles->data_pointer(id_a)->set_collision_time(time_collision);
  printd("collision type %d particle %d <-> %d time: %g\n",
	particles->data(id_a).process_type(),
	particles->data(id_a).id_partner(),
	particles->data(id_b).id_partner(),
	particles->data(id_a).collision_time());

  return act;
}


std::vector<ActionPtr>
ScatterActionsFinder::find_possible_actions (Particles *particles,
					const ExperimentParameters &parameters,
					CrossSections *cross_sections) const
{
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
      if (-distance.DotThree() > neighborhood_radius_squared)
        continue;

      /* Check if collision is possible. */
      ActionPtr act = check_collision (id_a, id_b, particles, parameters, cross_sections);

      /* Add to collision list. */
      if (act != nullptr)
	actions.push_back (std::move(act));
    }
  }

  /* Sort action list by collision time. */
  std::sort (actions.begin(), actions.end());

  return actions;
}



GridScatterFinder::GridScatterFinder(float length) : length_(length)
{
}


std::vector<ActionPtr>
GridScatterFinder::find_possible_actions (Particles *particles,
					   const ExperimentParameters &parameters,
					   CrossSections *cross_sections) const
{
  std::vector<std::vector<std::vector<std::vector<int> > > > grid;
  int N;
  int x, y, z;
  std::vector<ActionPtr> actions;

  /* For small boxes no point in splitting up in grids */
  /* calculate approximate grid size according to double interaction length */
  N = round(length_ / sqrt(parameters.cross_section * fm2_mb * M_1_PI) * 0.5);
  if (unlikely(N < 4 || particles->size() < 10))
    /* XXX: apply periodic boundary condition */
    return ScatterActionsFinder::find_possible_actions (particles, parameters, cross_sections);

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
              particles->data_pointer(*id_other)
                  ->set_position(particles->data(*id_other).position() + shift);
	      act = check_collision (data.id(), *id_other, particles, parameters, cross_sections);
              particles->data_pointer(*id_other)
                  ->set_position(particles->data(*id_other).position() - shift);
            }
            /* Add to collision list. */
	    if (act != nullptr)
	      actions.push_back (std::move(act));
          } /* grid particles loop */
        }   /* grid sz */
      }     /* grid sy */
    }       /* grid sx */
  }         /* outer particle loop */

  /* Sort action list by collision time. */
  std::sort (actions.begin(), actions.end());

  return actions;
}


}  // namespace Smash
