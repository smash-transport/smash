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
#include "include/constants.h"
#include "include/initial-conditions.h"
#include "include/param-reader.h"
#include "include/particles.h"
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
      && position.x1() < box.a() && position.x2() < box.a()
      && position.x3() < box.a())
    goto out;

  /* XXX: add hard wall conditions too */
  /* Enforce periodic boundary condition */
  if (position.x1() < 0)
    position.set_x1(position.x1() + box.a());

  if (position.x2() < 0)
    position.set_x2(position.x2() + box.a());

  if (position.x3() < 0)
    position.set_x3(position.x3() + box.a());

  if (position.x1() > box.a())
    position.set_x1(position.x1() - box.a());

  if (position.x2() > box.a())
    position.set_x2(position.x2() - box.a());

  if (position.x3() > box.a())
    position.set_x3(position.x3() - box.a());

 out:
    return position;
}

/* check_collision - check if a collision can happen betwenn particles */
static void check_collision(ParticleData *particle,
  std::list<int> *collision_list, box box, int number) {
  std::vector<std::vector<std::vector<std::vector<int> > > > grid;
  int N;
  int x, y, z;

  /* For small boxes no point in splitting up in grids */
  N = box.grid_number();
  if (unlikely(N < 4 || number < 10)) {
    FourVector distance;
    double radial_interaction = sqrt(box.cross_section() * fm2_mb * M_1_PI) * 2;
    for (int id = 0; id < number - 1; id++)
      for (int id_other = id + 1; id_other < number; id_other++) {
        /* XXX: apply periodic boundary condition */
        distance = particle[id].x() - particle[id_other].x();
        /* skip particles that are double interaction radius length away */
        if (distance > radial_interaction)
            continue;
        check_collision_criteria(particle, collision_list, box, id, id_other);
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
  for (int id = 0; id < number; id++) {
    /* XXX: function - map particle position to grid number */
    z = round(particle[id].x().x1() / box.a() * (N - 1));
    x = round(particle[id].x().x2() / box.a() * (N - 1));
    y = round(particle[id].x().x3() / box.a() * (N - 1));
    printd_position(particle[id]);
    printd("grid cell %i: %i %i %i of %i\n", id, z, x, y, N);
    grid[z][x][y].push_back(id);
  }

  /* semi optimised nearest neighbour search:
   * http://en.wikipedia.org/wiki/Cell_lists
   */
  FourVector shift;
  for (int id = 0; id < number - 1; id++) {
    /* XXX: function - map particle position to grid number */
    z = round(particle[id].x().x1() / box.a() * (N - 1));
    x = round(particle[id].x().x2() / box.a() * (N - 1));
    y = round(particle[id].x().x3() / box.a() * (N - 1));
    printd("grid cell %i: %i %i %i of %i\n", id, z, x, y, N);
    /* check all neighbour grids */
    for (int cz = -1; cz < 2; cz++) {
      int sz = cz + z;
      /* apply periodic boundary condition for particle positions */
      if (sz < 0) {
        sz = N - 1;
        shift.set_x1(-box.a());
      } else if (sz > N - 1) {
        sz = 0;
        shift.set_x1(box.a());
      } else {
        shift.set_x1(0);
      }
      for (int cx = -1; cx <  2; cx++) {
        int sx = cx + x;
        if (sx < 0) {
          sx = N - 1;
          shift.set_x2(-box.a());
        } else if (sx > N - 1) {
          sx = 0;
          shift.set_x2(box.a());
        } else {
          shift.set_x2(0);
        }
        for (int cy = -1; cy < 2; cy++) {
          int sy = cy + y;
          if (sy < 0) {
            sy = N - 1;
            shift.set_x3(-box.a());
          } else if (sy > N - 1) {
            sy = 0;
            shift.set_x3(box.a());
          } else {
            shift.set_x3(0);
          }
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
            if (shift == 0) {
              check_collision_criteria(particle, collision_list, box, id,
                *id_other);
            } else {
              /* apply eventual boundary before and restore after */
              particle[*id_other].set_position(particle[*id_other].x() + shift);
              check_collision_criteria(particle, collision_list, box, id,
                *id_other);
              particle[*id_other].set_position(particle[*id_other].x() - shift);
            }
          } /* grid particles loop */
        } /* grid sy */
      } /* grid sx */
    } /* grid sz */
  } /* outer particle loop */
}

/* Evolve - the core of the box, stepping forward in time */
static int Evolve(ParticleData *particles, ParticleType *particle_type,
  std::map<int, int> *map_type, int &number, const box &box) {
  FourVector distance, position;
  std::list<int> collision_list;
  size_t scatterings_total = 0;

  /* startup values */
  print_measurements(particles, number, scatterings_total, box);

  for (int steps = 0; steps < box.steps(); steps++) {
    /* fill collision table by cells */
    check_collision(particles, &collision_list, box, number);

    /* particle interactions */
    if (!collision_list.empty()) {
      scatterings_total += collision_list.size();
      collide_particles(particles, particle_type, map_type, &collision_list);
    }

    /* propagate all particles */
    for (int i = 0; i < number; i++) {
      distance.set_FourVector(1.0, particles[i].velocity_z(),
        particles[i].velocity_x(), particles[i].velocity_y());
      distance *= box.eps();
      printd("Particle %d motion: %g %g %g %g\n", particles[i].id(),
         distance.x0(), distance.x1(), distance.x2(), distance.x3());

      /* treat the box boundaries */
      position = particles[i].x();
      position += distance;
      position = boundary_condition(position, box);
      particles[i].set_position(position);
      printd_position(particles[i]);
    }

    /* physics output during the run */
    if (steps > 0 && (steps + 1) % box.update() == 0) {
      print_measurements(particles, number, scatterings_total, box);
      /* save evolution data */
      write_particles(particles, number);
      write_vtk(particles, number);
    }
  }

  if (likely(box.steps() > 0))
    print_tail(box,
      scatterings_total * 2 / (particles[0].x().x0() - 1.0) / number);
  return 0;
}

int main(int argc, char *argv[]) {
  char *p, *path;
  int opt, rc, number = 0;
  ParticleData *particles = NULL;
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
      cube->set_a(atof(optarg));
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
  particles = initial_conditions(particles, particle_types, &map_type, number,
    cube);
  write_particles(particles, number);

  /* Compute stuff */
  rc = Evolve(particles, particle_types, &map_type, number, *cube);

  /* tear down */
  delete [] particles;
  delete [] particle_types;
  delete cube;
  free(path);
  return rc;
}
