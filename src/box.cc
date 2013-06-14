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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <list>
#include <map>
#include <vector>

#include "include/Box.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/collisions.h"
#include "include/constants.h"
#include "include/decays.h"
#include "include/input-particles.h"
#include "include/initial-conditions.h"
#include "include/macros.h"
#include "include/param-reader.h"
#include "include/Parameters.h"
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
         "  -s, --steps          number of steps\n"
         "  -T, --temperature    initial temperature\n"
         "  -u, --update         output measurements each nth steps\n"
         "  -V, --version\n\n");
  exit(rc);
}

/* boundary_condition - enforce specific type of boundaries */
FourVector boundary_condition(FourVector position, const Box &box) {
  /* Check positivity and box size */
  if (position.x1() > 0 && position.x2() > 0 && position.x3() > 0
      && position.x1() < box.length() && position.x2() < box.length()
      && position.x3() < box.length())
    goto out;

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
static void check_collision_geometry(std::map<int, ParticleData> *particle,
  std::vector<ParticleType> *particle_type, std::map<int, int> *map_type,
  std::list<int> *collision_list, Parameters const &parameters,
  Box const &box) {
  std::vector<std::vector<std::vector<std::vector<int> > > > grid;
  int N;
  int x, y, z;

  /* For small boxes no point in splitting up in grids */
  /* calculate approximate grid size according to double interaction length */
  N = round(box.length()
            / sqrt(parameters.cross_section() * fm2_mb * M_1_PI) * 0.5);
  if (unlikely(N < 4 || particle->size() < 10)) {
    FourVector distance;
    double radial_interaction = sqrt(parameters.cross_section() * fm2_mb
                                     * M_1_PI) * 2;
    for (std::map<int, ParticleData>::iterator i = particle->begin();
         i != particle->end(); ++i) {
      /* The particle has formed a resonance or has decayed
       * and is not active anymore
       */
      if (i->second.process_type() > 0)
        printf("Attention: i %i has process type %i \n", i->first,
               i->second.process_type());

      /* Check resonances for decays first */
      if ((*particle_type)[(*map_type)[i->first]].width() > 0.0) {
        if (does_decay(&(i->second),
                       &(*particle_type)[(*map_type)[i->first]],
                       collision_list, parameters))
          continue;
      }

      for (std::map<int, ParticleData>::iterator j; j != particle->end();
         ++j) {
        /* exclude check on same particle */
        if (i->first != j->first)
          continue;
        /* The other particle has formed a resonance or has decayed
         * and is not active anymore
         */
        if (j->second.process_type() > 0)
          printf("Attention: j %i has process type %i \n", j->first,
               j->second.process_type());

        /* Check resonances for decays here too */
        if ((*particle_type)[(*map_type)[j->first]].width() > 0.0) {
          if (does_decay(&(j->second), &(*particle_type)[(*map_type)[j->first]],
                         collision_list, parameters))
              continue;
        }

        /* XXX: apply periodic boundary condition */
        distance = i->second.position() - j->second.position();
        /* skip particles that are double interaction radius length away */
        if (distance > radial_interaction)
           continue;
        collision_criteria_geometry(particle, particle_type, map_type,
                              collision_list, parameters, i->first, j->first);
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
  for (std::map<int, ParticleData>::iterator i = particle->begin();
         i != particle->end(); ++i) {
    /* XXX: function - map particle position to grid number */
    z = round(i->second.position().x1() / box.length() * (N - 1));
    x = round(i->second.position().x2() / box.length() * (N - 1));
    y = round(i->second.position().x3() / box.length() * (N - 1));
    printd_position(i->second);
    printd("grid cell %i: %i %i %i of %i\n", i->first, z, x, y, N);
    if (unlikely(z >= N || x >= N || y >= N)) {
      printf("Particle position: %g %g %g \n", i->second.position().x1(),
             i->second.position().x2(), i->second.position().x3());
      double coord_time = i->second.position().x0();
      double coord[3];
      coord[0] = i->second.position().x1();
      coord[1] = i->second.position().x2();
      coord[2] = i->second.position().x3();
      for (int coordi = 0; coordi < 3; coordi++) {
        if (coord[coordi] > box.length())
          coord[coordi] -= box.length();
      }
      i->second.set_position(coord_time, coord[0], coord[1], coord[2]);
      z = round(i->second.position().x1() / box.length() * (N - 1));
      x = round(i->second.position().x2() / box.length() * (N - 1));
      y = round(i->second.position().x3() / box.length() * (N - 1));
    }
    grid[z][x][y].push_back(i->first);
  }
  /* semi optimised nearest neighbour search:
   * http://en.wikipedia.org/wiki/Cell_lists
   */
  FourVector shift;
  for (std::map<int, ParticleData>::iterator i = particle->begin();
       i != particle->end(); ++i) {
    /* The particle has formed a resonance or has decayed
     * and is not active anymore
     */
    if (i->second.process_type() > 0)
      printf("Attention: i %i has process type %i \n", i->first,
             i->second.process_type());

    /* Check resonances for decay */
    if ((*particle_type)[(*map_type)[i->first]].width() > 0.0) {
        if (does_decay(&(i->second), &(*particle_type)[(*map_type)[i->first]],
                       collision_list, parameters))
        continue;
    }

    /* XXX: function - map particle position to grid number */
    z = round(i->second.position().x1() / box.length() * (N - 1));
    x = round(i->second.position().x2() / box.length() * (N - 1));
    y = round(i->second.position().x3() / box.length() * (N - 1));
    if (unlikely(z >= N || x >= N || y >= N))
      printf("grid cell %i: %i %i %i of %i\n", i->first, z, x, y, N);
    /* check all neighbour grids */
    for (int cz = -1; cz < 2; cz++) {
      int sz = cz + z;
      /* apply periodic boundary condition for particle positions */
      if (sz < 0) {
        sz = N - 1;
        shift.set_x1(-box.length());
      } else if (sz > N - 1) {
        sz = 0;
        shift.set_x1(box.length());
      } else {
        shift.set_x1(0);
      }
      for (int cx = -1; cx <  2; cx++) {
        int sx = cx + x;
        if (sx < 0) {
          sx = N - 1;
          shift.set_x2(-box.length());
        } else if (sx > N - 1) {
          sx = 0;
          shift.set_x2(box.length());
        } else {
          shift.set_x2(0);
        }
        for (int cy = -1; cy < 2; cy++) {
          int sy = cy + y;
          if (sy < 0) {
            sy = N - 1;
            shift.set_x3(-box.length());
          } else if (sy > N - 1) {
            sy = 0;
            shift.set_x3(box.length());
          } else {
            shift.set_x3(0);
          }
          /* empty grid cell */
          if (grid[sz][sx][sy].empty())
            continue;
          /* grid cell particle list */
          for (std::vector<int>::iterator id_other
               = grid[sz][sx][sy].begin(); id_other != grid[sz][sx][sy].end();
               ++id_other) {
            /* only check against particles above current id
             * to avoid double counting
             */
            if (*id_other <= i->first)
              continue;

            /* The other particle has formed a resonance or has decayed
             * and is not active anymore
             */
            if ((*particle)[*id_other].process_type() > 0)
              printf("Attention: j %i has process type %i \n", *id_other,
               (*particle)[*id_other].process_type());

            /* Check resonances for decay */
            if ((*particle_type)[(*map_type)[*id_other]].width() > 0.0) {
              if (does_decay(&((*particle)[*id_other]),
                             &(*particle_type)[(*map_type)[*id_other]],
                             collision_list, parameters))
                continue;
            }

            printd("grid cell particle %i <-> %i\n", i->first, *id_other);
            if (shift == 0) {
              collision_criteria_geometry(particle, particle_type, map_type,
                             collision_list, parameters, i->first, *id_other);
            } else {
              /* apply eventual boundary before and restore after */
              (*particle)[*id_other].set_position(
                (*particle)[*id_other].position() + shift);
              collision_criteria_geometry(particle, particle_type, map_type,
                              collision_list, parameters, i->first, *id_other);
              (*particle)[*id_other].set_position(
                (*particle)[*id_other].position() - shift);
            }
          } /* grid particles loop */
        } /* grid sy */
      } /* grid sx */
    } /* grid sz */
  } /* outer particle loop */
}

/* Evolve - the core of the box, stepping forward in time */
static int Evolve(std::map<int, ParticleData> *particles,
  std::vector<ParticleType> *particle_type, std::map<int, int> *map_type,
  const Parameters &parameters, const Box &box, size_t *largest_id) {
  std::list<int> collision_list;
  size_t scatterings_total = 0, previous_scatterings_total = 0,
    scatterings_this_interval = 0;

  /* startup values */
  print_measurements(*particles, scatterings_total,
                     scatterings_this_interval, box);

  for (int steps = 0; steps < box.steps(); steps++) {
    /* fill collision table by cells */
    check_collision_geometry(particles, particle_type, map_type,
                                &collision_list, parameters, box);

    /* particle interactions */
    if (!collision_list.empty())
      scatterings_total = collide_particles(particles, particle_type,
        map_type, &collision_list, scatterings_total, largest_id);

    /* propagate all particles */
    propagate_particles(particles, parameters, box);

    /* physics output during the run */
    if (steps > 0 && (steps + 1) % parameters.output_interval() == 0) {
      scatterings_this_interval = scatterings_total
        - previous_scatterings_total;

      previous_scatterings_total = scatterings_total;

      print_measurements(*particles, scatterings_total,
                         scatterings_this_interval, box);
      /* save evolution data */
      write_particles(*particles);
      write_vtk(*particles);
    }
  }

  if (likely(box.steps() > 0))
    print_tail(box, scatterings_total * 2
     / ((*particles)[0].position().x0() - 1.0) / particles->size());
  return 0;
}

int main(int argc, char *argv[]) {
  char *p, *path;
  int opt, rc;
  std::map<int, ParticleData> particles;
  std::vector<ParticleType> particle_types;
  Box *cube = new Box;
  Parameters *parameters = new Parameters;
  std::map<int, int> map_type;
  size_t largest_id = 0;

  struct option longopts[] = {
    { "eps",        required_argument,      0, 'e' },
    { "help",       no_argument,            0, 'h' },
    { "length",     required_argument,      0, 'l' },
    { "output-interval", required_argument,      0, 'O' },
    { "random",     required_argument,      0, 'r' },
    { "steps",      required_argument,      0, 's' },
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
  int len = 3;
  path = reinterpret_cast<char *>(malloc(len));
  /* XXX: make path configurable */
  snprintf(path, len, "./");
  process_params(cube, parameters, path);

  /* parse the command line options, they override all previous */
  while ((opt = getopt_long(argc, argv, "e:hl:O:r:s:T:V", longopts,
    NULL)) != -1) {
    switch (opt) {
    case 'e':
      parameters->set_eps(atof(optarg));
      break;
    case 'h':
      usage(EXIT_SUCCESS);
      break;
    case 'l':
      cube->set_length(atof(optarg));
      break;
    case 'O':
      parameters->set_output_interval(abs(atoi(optarg)));
      break;
    case 'r':
      /* negative seed is for time */
      if (atol(optarg) > 0)
        parameters->set_seed(atol(optarg));
      else
        parameters->set_seed(time(NULL));
      break;
    case 's':
      cube->set_steps(abs(atoi(optarg)));
      break;
    case 'T':
      cube->set_temperature(atof(optarg));
      break;
    case 'V':
      exit(EXIT_SUCCESS);
    default:
      usage(EXIT_FAILURE);
    }
  }

  /* Output IC values */
  print_startup(*cube, *parameters);
  mkdir_data();
  write_oscar_header();

  /* Initialize box */
  input_particles(&particle_types, path);
  initial_conditions(&particles, &particle_types, &map_type, parameters, cube,
                     &largest_id);
  write_particles(particles);

  /* Compute stuff */
  rc = Evolve(&particles, &particle_types, &map_type, *parameters, *cube,
              &largest_id);

  /* tear down */
  particles.clear();
  particle_types.clear();
  delete cube;
  delete parameters;
  free(path);
  return rc;
}
