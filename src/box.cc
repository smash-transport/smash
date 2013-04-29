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

#include "include/ParticleData.h"
#include "include/ParticleType.h"
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

/* Evolve - the core of the box, stepping forward in time */
static int Evolve(ParticleData *particles, ParticleType *particle_type,
  int &number, const box &box) {
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
      collide_particles(particles, &collision_list);
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
  particles = initial_conditions(particles, particle_types, number, cube);
  write_particles(particles, number);

  /* Compute stuff */
  rc = Evolve(particles, particle_types, number, *cube);

  /* tear down */
  delete [] particles;
  delete [] particle_types;
  delete cube;
  free(path);
  return rc;
}
