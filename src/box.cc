/*
 *
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/box.h"

#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <list>
#include <map>
#include <vector>

#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/collisions.h"
#include "include/constants.h"
#include "include/initial-conditions.h"
#include "include/param-reader.h"
#include "include/outputroutines.h"

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
static FourVector boundary_condition(FourVector position, const box &box) {
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
static void check_collision_geometry(std::vector<ParticleData> *particle,
  std::list<int> *collision_list, box box) {
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > grid;
  int N;
  int x, y, z;

  /* For small boxes no point in splitting up in grids */
  N = box.grid_number();
  if (unlikely(N < 4 || particle->size() < 10)) {
    FourVector distance;
    double radial_interaction = sqrt(box.cross_section() * fm2_mb * M_1_PI) * 2;
    for (size_t id = 0; id < particle->size() - 1; id++)
      for (size_t id_other = id + 1; id_other < particle->size(); id_other++) {
        /* XXX: apply periodic boundary condition */
        distance = (*particle)[id].position()
          - (*particle)[id_other].position();
        /* skip particles that are double interaction radius length away */
        if (distance > radial_interaction)
            continue;
        collision_criteria_geometry(particle, collision_list, box, id,
          id_other);
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
  for (size_t id = 0; id < particle->size(); id++) {
    /* XXX: function - map particle position to grid number */
    z = round((*particle)[id].position().x1() / box.length() * (N - 1));
    x = round((*particle)[id].position().x2() / box.length() * (N - 1));
    y = round((*particle)[id].position().x3() / box.length() * (N - 1));
    printd_position((*particle)[id]);
    printd("grid cell %lu: %i %i %i of %i\n", id, z, x, y, N);
    grid[z][x][y].push_back(id);
  }

  /* semi optimised nearest neighbour search:
   * http://en.wikipedia.org/wiki/Cell_lists
   */
  FourVector shift;
  for (size_t id = 0; id < particle->size() - 1; id++) {
    /* XXX: function - map particle position to grid number */
    z = round((*particle)[id].position().x1() / box.length() * (N - 1));
    x = round((*particle)[id].position().x2() / box.length() * (N - 1));
    y = round((*particle)[id].position().x3() / box.length() * (N - 1));
    printd("grid cell %lu: %i %i %i of %i\n", id, z, x, y, N);
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
          for (std::vector<unsigned int>::iterator id_other
               = grid[sz][sx][sy].begin(); id_other != grid[sz][sx][sy].end();
               ++id_other) {
	    /* only check against particles above current id
	     * to avoid double counting
	     */
            if (*id_other <= id)
              continue;

            printd("grid cell particle %lu <-> %i\n", id, *id_other);
            if (shift == 0) {
              collision_criteria_geometry(particle, collision_list, box, id,
                *id_other);
            } else {
              /* apply eventual boundary before and restore after */
              (*particle)[*id_other].set_position(
                (*particle)[*id_other].position() + shift);
              collision_criteria_geometry(particle, collision_list, box, id,
                *id_other);
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
static int Evolve(std::vector<ParticleData> *particles,
  ParticleType *particle_type, std::map<int, int> *map_type, const box &box) {
  FourVector distance, position;
  std::list<int> collision_list;
  size_t scatterings_total = 0;

  /* startup values */
  print_measurements(*particles, scatterings_total, box);

  for (int steps = 0; steps < box.steps(); steps++) {
    /* fill collision table by cells */
    check_collision_geometry(particles, &collision_list, box);

    /* particle interactions */
    if (!collision_list.empty()) 
      scatterings_total += collide_particles(particles, particle_type,
        map_type, &collision_list, scatterings_total);

    /* propagate all particles */
    for (size_t i = 0; i < particles->size(); i++) {
      distance.set_FourVector(1.0, (*particles)[i].velocity_x(),
        (*particles)[i].velocity_y(), (*particles)[i].velocity_z());
      distance *= box.eps();
      printd("Particle %d motion: %g %g %g %g\n", (*particles)[i].id(),
         distance.x0(), distance.x1(), distance.x2(), distance.x3());

      /* treat the box boundaries */
      position = (*particles)[i].position();
      position += distance;
      position = boundary_condition(position, box);
      (*particles)[i].set_position(position);
      printd_position((*particles)[i]);
    }

    /* physics output during the run */
    if (steps > 0 && (steps + 1) % box.update() == 0) {
      print_measurements(*particles, scatterings_total, box);
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
  std::vector<ParticleData> particles;
  ParticleType *particle_types = NULL;
  box *cube = new box;
  std::map<int, int> map_type;

  struct option longopts[] = {
    { "eps",        required_argument,      0, 'e' },
    { "help",       no_argument,            0, 'h' },
    { "length",     required_argument,      0, 'l' },
    { "random",     required_argument,      0, 'r' },
    { "steps",      required_argument,      0, 's' },
    { "temperature", required_argument,     0, 'T' },
    { "update",     required_argument,      0, 'u' },
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
  snprintf(path, len, "./");
  process_params(cube, path);

  /* parse the command line options, they override all previous */
  while ((opt = getopt_long(argc, argv, "e:hl:r:s:T:u:V", longopts,
    NULL)) != -1) {
    switch (opt) {
    case 'e':
      cube->set_eps(atof(optarg));
      break;
    case 'h':
      usage(EXIT_SUCCESS);
      break;
    case 'l':
      cube->set_length(atof(optarg));
      break;
    case 'r':
      /* negative seed is for time */
      if (atol(optarg) > 0)
        cube->set_seed(atol(optarg));
      else
        cube->set_seed(time(NULL));
      break;
    case 's':
      cube->set_steps(abs(atoi(optarg)));
      break;
    case 'T':
      cube->set_temperature(atof(optarg));
      break;
    case 'u':
      cube->set_update(abs(atoi(optarg)));
      break;
    case 'V':
      exit(EXIT_SUCCESS);
    default:
      usage(EXIT_FAILURE);
    }
  }

  /* Output IC values */
  print_startup(*cube);
  mkdir_data();
  write_oscar_header();

  /* Initialize box */
  particle_types = initial_particles(particle_types);
  initial_conditions(&particles, particle_types, &map_type, cube);
  write_particles(particles);

  /* Compute stuff */
  rc = Evolve(&particles, particle_types, &map_type, *cube);

  /* tear down */
  particles.clear();
  delete [] particle_types;
  delete cube;
  free(path);
  return rc;
}
