/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */

#include "include/Box.h"
#include "include/CrossSections.h"
#include "include/Particles.h"
#include "include/constants.h"
#include "include/collisions.h"
#include "include/decays.h"
#include "include/initial-conditions.h"
#include "include/input-decaymodes.h"
#include "include/input-particles.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/param-reader.h"
#include "include/propagation.h"

/* check_collision_geometry - check if a collision happens between particles */
static void check_collision_geometry(Particles *particles,
  CrossSections *cross_sections,
  std::list<int> *collision_list, Box const &box,
  size_t *rejection_conflict) {
  std::vector<std::vector<std::vector<std::vector<int> > > > grid;
  int N;
  int x, y, z;

  /* For small boxes no point in splitting up in grids */
  /* calculate approximate grid size according to double interaction length */
  N = round(box.length()
            / sqrt(box.cross_section() * fm2_mb * M_1_PI) * 0.5);
  if (unlikely(N < 4 || particles->size() < 10)) {
    FourVector distance;
    double radial_interaction = sqrt(box.cross_section() * fm2_mb
                                     * M_1_PI) * 2;
    for (std::map<int, ParticleData>::iterator i = particles->begin();
         i != particles->end(); ++i) {
      for (std::map<int, ParticleData>::iterator j = particles->begin();
           j != particles->end(); ++j) {
        /* exclude check on same particle and double counting */
        if (i->first >= j->first)
          continue;

        /* XXX: apply periodic boundary condition */
        distance = i->second.position() - j->second.position();
        /* skip particles that are double interaction radius length away */
        if (distance > radial_interaction)
           continue;
        collision_criteria_geometry(particles, cross_sections, collision_list,
         box, i->first, j->first, rejection_conflict);
      }
    }
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
  for (std::map<int, ParticleData>::iterator i = particles->begin();
         i != particles->end(); ++i) {
    /* XXX: function - map particle position to grid number */
    x = round(i->second.position().x1() / box.length() * (N - 1));
    y = round(i->second.position().x2() / box.length() * (N - 1));
    z = round(i->second.position().x3() / box.length() * (N - 1));
    printd_position(i->second);
    printd("grid cell particle %i: %i %i %i of %i\n", i->first, x, y, z, N);
    if (unlikely(x >= N || y >= N || z >= N))
      printf("W: Particle outside the box: %g %g %g \n",
             i->second.position().x1(), i->second.position().x2(),
             i->second.position().x3());
    grid[x][y][z].push_back(i->first);
  }
  /* semi optimised nearest neighbour search:
   * http://en.wikipedia.org/wiki/Cell_lists
   */
  FourVector shift;
  for (std::map<int, ParticleData>::iterator i = particles->begin();
       i != particles->end(); ++i) {
    /* XXX: function - map particle position to grid number */
    x = round(i->second.position().x1() / box.length() * (N - 1));
    y = round(i->second.position().x2() / box.length() * (N - 1));
    z = round(i->second.position().x3() / box.length() * (N - 1));
    if (unlikely(x >= N || y >= N || z >= N))
      printf("grid cell particle %i: %i %i %i of %i\n", i->first, x, y, z, N);
    /* check all neighbour grids */
    for (int cx = -1; cx < 2; cx++) {
      int sx = cx + x;
      /* apply periodic boundary condition for particle positions */
      if (sx < 0) {
        sx = N - 1;
        shift.set_x1(-box.length());
      } else if (sx > N - 1) {
        sx = 0;
        shift.set_x1(box.length());
      } else {
        shift.set_x1(0);
      }
      for (int cy = -1; cy <  2; cy++) {
        int sy = cy + y;
        if (sy < 0) {
          sy = N - 1;
          shift.set_x2(-box.length());
        } else if (sy > N - 1) {
          sy = 0;
          shift.set_x2(box.length());
        } else {
          shift.set_x2(0);
        }
        for (int cz = -1; cz < 2; cz++) {
          int sz = cz + z;
          if (sz < 0) {
            sz = N - 1;
            shift.set_x3(-box.length());
          } else if (sz > N - 1) {
            sz = 0;
            shift.set_x3(box.length());
          } else {
            shift.set_x3(0);
          }
          /* empty grid cell */
          if (grid[sx][sy][sz].empty())
            continue;
          /* grid cell particle list */
          for (std::vector<int>::iterator id_other
               = grid[sx][sy][sz].begin(); id_other != grid[sx][sy][sz].end();
               ++id_other) {
            /* only check against particles above current id
             * to avoid double counting
             */
            if (*id_other <= i->first)
              continue;

            printd("grid cell particle %i <-> %i\n", i->first, *id_other);
            if (shift == 0) {
              collision_criteria_geometry(particles, cross_sections,
                collision_list, box, i->first, *id_other,
                rejection_conflict);
            } else {
              /* apply eventual boundary before and restore after */
              particles->data_pointer(*id_other)->set_position(
                particles->data(*id_other).position() + shift);
              collision_criteria_geometry(particles, cross_sections,
                collision_list, box, i->first, *id_other,
                rejection_conflict);
              particles->data_pointer(*id_other)->set_position(
                particles->data(*id_other).position() - shift);
            }
          } /* grid particles loop */
        } /* grid sz */
      } /* grid sy */
    } /* grid sx */
  } /* outer particle loop */
}

/* evolve - the core of the box, stepping forward in time */
int Box::evolve(Particles *particles, CrossSections *cross_sections) {
  std::list<int> collision_list, decay_list;
  size_t interactions_total = 0, previous_interactions_total = 0,
    interactions_this_interval = 0;
  size_t rejection_conflict = 0;
  int resonances = 0, decays = 0;

  /* fixup positions on startup, particles need to be *inside* the box */
  for (std::map<int, ParticleData>::iterator i = particles->begin();
       i != particles->end(); ++i) {
    bool boundary_hit = false;
    i->second.set_position(boundary_condition(i->second.position(), *this,
                           &boundary_hit));
  }

  /* startup values */
  print_measurements(*particles, interactions_total,
                     interactions_this_interval, *this);

  for (int step = 0; step < this->steps(); step++) {
    /* Check resonances for decays */
    check_decays(particles, &decay_list, *this);

    /* Do the decays */
    if (!decay_list.empty()) {
      decays += decay_list.size();
      interactions_total = decay_particles(particles,
        &decay_list, interactions_total);
    }

    /* fill collision table by cells */
    check_collision_geometry(particles, cross_sections,
      &collision_list, *this, &rejection_conflict);

    /* particle interactions */
    if (!collision_list.empty()) {
      printd_list(collision_list);
      interactions_total = collide_particles(particles, &collision_list,
        interactions_total, &resonances);
    }

    /* propagate all particles */
    propagate_particles(particles, *this);

    /* physics output during the run */
    if (step > 0 && (step + 1) % this->output_interval() == 0) {
      interactions_this_interval = interactions_total
        - previous_interactions_total;

      previous_interactions_total = interactions_total;

      print_measurements(*particles, interactions_total,
                         interactions_this_interval, *this);
      printd("Resonances: %i Decays: %i\n", resonances, decays);
      printd("Ignored collisions %zu\n", rejection_conflict);
      /* save evolution data */
      write_measurements(*particles, interactions_total,
        interactions_this_interval, resonances, decays, rejection_conflict);
      write_vtk(*particles);
    }
  }

  /* Guard against evolution */
  if (likely(this->steps() > 0)) {
    /* if there are not particles no interactions happened */
    if (likely(!particles->empty()))
      print_tail(*this, interactions_total * 2
                 / particles->time() / particles->size());
    else
      print_tail(*this, 0);
    printf("Total ignored collisions: %zu\n", rejection_conflict);
  }
  return 0;
}
