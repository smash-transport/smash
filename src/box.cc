/*
 *
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */

#include <getopt.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <list>
#include <map>
#include <vector>

#include "include/Box.h"
#include "include/CrossSections.h"
#include "include/FourVector.h"
#include "include/Parameters.h"
#include "include/Particles.h"
#include "include/ParticleData.h"
#include "include/collisions.h"
#include "include/constants.h"
#include "include/decays.h"
#include "include/input-decaymodes.h"
#include "include/input-particles.h"
#include "include/initial-conditions.h"
#include "include/Laboratory.h"
#include "include/macros.h"
#include "include/param-reader.h"
#include "include/outputroutines.h"
#include "include/propagation.h"

/* build dependent variables */
#include "include/Config.h"

char *progname;

static void usage(int rc) {
  printf("\nUsage: %s [option]\n\n", progname);
  printf("Calculate transport box\n"
         "  -e, --eps            time step\n"
         "  -h, --help           usage information\n"
         "  -l, --length         length of the box in fermi\n"
         "  -O, --output-interval          step interval between measurements\n"
         "  -R, --random         random number seed\n"
         "  -s, --sigma          cross section in mbarn\n"
         "  -S, --steps          number of steps\n"
         "  -T, --temperature    initial temperature\n"
         "  -V, --version\n\n");
  exit(rc);
}

/* boundary_condition - enforce specific type of boundaries
 *
 * This assumes that the particle is at most one box length
 * away from the boundary to shift it in.
 */
FourVector boundary_condition(FourVector position, const Box &box,
                              bool *boundary_hit) {
  /* Check positivity and box size */
  if (position.x1() > 0 && position.x2() > 0 && position.x3() > 0
      && position.x1() < box.length() && position.x2() < box.length()
      && position.x3() < box.length())
    goto out;

  *boundary_hit = true;

  /* Enforce periodic boundary condition */
  if (position.x1() < 0)
    position.set_x1(position.x1() + box.length());

  if (position.x2() < 0)
    position.set_x2(position.x2() + box.length());

  if (position.x3() < 0)
    position.set_x3(position.x3() + box.length());

  if (position.x1() > box.length())
    position.set_x1(position.x1() - box.length());

  if (position.x2() > box.length())
    position.set_x2(position.x2() - box.length());

  if (position.x3() > box.length())
    position.set_x3(position.x3() - box.length());

 out:
    return position;
}

/* check_collision_geometry - check if a collision happens between particles */
static void check_collision_geometry(Particles *particles,
  CrossSections *cross_sections,
  std::list<int> *collision_list, Laboratory const &parameters,
  Box const &box, size_t *rejection_conflict) {
  std::vector<std::vector<std::vector<std::vector<int> > > > grid;
  int N;
  int x, y, z;

  /* For small boxes no point in splitting up in grids */
  /* calculate approximate grid size according to double interaction length */
  N = round(box.length()
            / sqrt(parameters.cross_section() * fm2_mb * M_1_PI) * 0.5);
  if (unlikely(N < 4 || particles->size() < 10)) {
    FourVector distance;
    double radial_interaction = sqrt(parameters.cross_section() * fm2_mb
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
         parameters, i->first, j->first, rejection_conflict);
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
                collision_list, parameters, i->first, *id_other,
                rejection_conflict);
            } else {
              /* apply eventual boundary before and restore after */
              particles->data_pointer(*id_other)->set_position(
                particles->data(*id_other).position() + shift);
              collision_criteria_geometry(particles, cross_sections,
                collision_list, parameters, i->first, *id_other,
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

/* Evolve - the core of the box, stepping forward in time */
static int Evolve(Particles *particles, CrossSections *cross_sections,
                  const Laboratory &parameters, const Box &box) {
  std::list<int> collision_list, decay_list;
  size_t interactions_total = 0, previous_interactions_total = 0,
    interactions_this_interval = 0;
  size_t rejection_conflict = 0;
  int resonances = 0, decays = 0;

  /* fixup positions on startup, particles need to be *inside* the box */
  for (std::map<int, ParticleData>::iterator i = particles->begin();
       i != particles->end(); ++i) {
    bool boundary_hit = false;
    i->second.set_position(boundary_condition(i->second.position(), box,
                           &boundary_hit));
  }

  /* startup values */
  print_measurements(*particles, interactions_total,
                     interactions_this_interval, box);

  for (int steps = 0; steps < parameters.steps(); steps++) {
    /* Check resonances for decays */
    check_decays(particles, &decay_list, parameters);

    /* Do the decays */
    if (!decay_list.empty()) {
      decays += decay_list.size();
      interactions_total = decay_particles(particles,
        &decay_list, interactions_total);
    }

    /* fill collision table by cells */
    check_collision_geometry(particles, cross_sections,
      &collision_list, parameters, box, &rejection_conflict);

    /* particle interactions */
    if (!collision_list.empty())
      interactions_total = collide_particles(particles, &collision_list,
        interactions_total, &resonances);

    /* propagate all particles */
    propagate_particles(particles, parameters, box);

    /* physics output during the run */
    if (steps > 0 && (steps + 1) % parameters.output_interval() == 0) {
      interactions_this_interval = interactions_total
        - previous_interactions_total;

      previous_interactions_total = interactions_total;

      print_measurements(*particles, interactions_total,
                         interactions_this_interval, box);
      printd("Resonances: %i Decays: %i\n", resonances, decays);
      printd("Ignored collisions %zu\n", rejection_conflict);
      /* save evolution data */
      write_measurements(*particles, interactions_total,
        interactions_this_interval, resonances, decays, rejection_conflict);
      write_vtk(*particles);
    }
  }

  /* Guard against evolution */
  if (likely(parameters.steps() > 0)) {
    /* if there are not particles no interactions happened */
    if (likely(!particles->empty()))
      print_tail(box, interactions_total * 2
                 / particles->time() / particles->size());
    else
      print_tail(box, 0);
    printf("Total ignored collisions: %zu\n", rejection_conflict);
  }
  return 0;
}

int main(int argc, char *argv[]) {
  char *p, *path;
  int opt, rc;
  Particles *particles = new Particles;
  Box *cube = new Box;
  Laboratory *parameters = new Laboratory;

  struct option longopts[] = {
    { "eps",        required_argument,      0, 'e' },
    { "help",       no_argument,            0, 'h' },
    { "length",     required_argument,      0, 'l' },
    { "output-interval", required_argument,      0, 'O' },
    { "random",     required_argument,      0, 'R' },
    { "sigma",      required_argument,      0, 's' },
    { "steps",      required_argument,      0, 'S' },
    { "temperature", required_argument,     0, 'T' },
    { "version",    no_argument,            0, 'V' },
    { NULL,         0, 0, 0 }
  };

  /* strip any path to progname */
  progname = argv[0];
  if ((p = strrchr(progname, '/')) != NULL)
    progname = p + 1;
  printf("%s (%d)\n", progname, VERSION_MAJOR);

  /* Read config file overrides box constructor defaults */
  std::vector<Parameters> configuration;
  int len = 3;
  path = reinterpret_cast<char *>(malloc(len));
  /* XXX: make path configurable */
  snprintf(path, len, "./");
  process_params(path, &configuration);
  assign_params(&configuration, cube);
  assign_params(&configuration, parameters);

  /* parse the command line options, they override all previous */
  while ((opt = getopt_long(argc, argv, "e:hl:O:R:s:S:T:V", longopts,
    NULL)) != -1) {
    switch (opt) {
    case 'e':
      parameters->set_eps(fabs(atof(optarg)));
      break;
    case 'h':
      usage(EXIT_SUCCESS);
      break;
    case 'l':
      cube->set_length(fabs(atof(optarg)));
      break;
    case 'O':
      {
      /* guard output_interval to be positive and greater null */
      int output_interval = abs(atoi(optarg));
      if (output_interval > 0)
        parameters->set_output_interval(output_interval);
      }
      break;
    case 'R':
      /* negative seed is for time */
      if (atol(optarg) > 0)
        parameters->set_seed(atol(optarg));
      else
        parameters->set_seed(time(NULL));
      break;
    case 's':
      parameters->set_cross_section(fabs(atof(optarg)));
      break;
    case 'S':
      parameters->set_steps(abs(atoi(optarg)));
      break;
    case 'T':
      cube->set_temperature(fabs(atof(optarg)));
      break;
    case 'V':
      exit(EXIT_SUCCESS);
    default:
      usage(EXIT_FAILURE);
    }
  }

  /* Output IC values */
  print_startup(*parameters);
  print_startup(*cube);
  mkdir_data();
  write_oscar_header();

  /* Initialize box */
  input_particles(particles, path);
  initial_conditions(particles, parameters, cube);
  input_decaymodes(particles, path);
  CrossSections *cross_sections = new CrossSections;
  cross_sections->add_elastic_parameter(parameters->cross_section());

  write_measurements_header(*particles);
  print_header();
  write_particles(*particles);

  /* Compute stuff */
  rc = Evolve(particles, cross_sections, *parameters, *cube);

  /* tear down */
  delete particles;
  delete cross_sections;
  delete cube;
  delete parameters;
  free(path);
  return rc;
}
