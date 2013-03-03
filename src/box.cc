/*
 *
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/box.h"

#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/initial-conditions.h"
#include "include/param-reader.h"
#include "include/outputroutines.h"

/* build dependent variables */
#include "include/Config.h"

char *progname;

/* Be chatty on demand */
bool verbose = 0;

/* Default random seed */
unsigned int seedp = 1;

static void usage(int rc) {
  printf("\nUsage: %s [option]\n\n", progname);
  printf("Calculate transport box\n"
         "  -e, --eps            time step\n"
         "  -h, --help           usage information\n"
         "  -l, --length         length of the box in fermi\n"
         "  -s, --steps          number of steps\n"
         "  -T, --temperature    initial temperature\n"
         "  -v, --verbose        show debug info\n"
         "  -V, --version\n\n");
  exit(rc);
}


/* boundary_condition - enforce specific type of boundaries */
static FourVector boundary_condition(FourVector position, FourVector distance,
       box box) {
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
static int Evolve(ParticleData *particles, int number, box box) {
  FourVector distance, position;

  for (int steps = 0; steps < box.steps(); steps++)
    for (int i = 0; i < number; i++) {
       /* XXX: interactions */
       /* calculate particles motion */
       distance.set_FourVector(1.0, particles[i].velocity_x(),
         particles[i].velocity_y(), particles[i].velocity_z());
       distance *= box.eps();
       printd("Particle %d motion: %g %g %g %g\n", particles[i].id(),
         distance.x0(), distance.x1(), distance.x2(), distance.x3());

       /* treat the box boundaries */
       position = particles[i].x();
       position += distance;
       position = boundary_condition(position, distance, box);
       particles[i].set_position(position);
       printd_position(particles[i]);

       /* save evolution data */
       if (i > 0 && i % box.update() == 0)
         write_particles(particles, number);
    }
  return 0;
}

int main(int argc, char *argv[]) {
  char *p, *path;
  int opt, rc, number = 0;
  ParticleData *particles = NULL;
  box box;

  struct option longopts[] = {
    { "eps",        no_argument,            0, 'e' },
    { "help",       no_argument,            0, 'h' },
    { "length",     no_argument,            0, 'l' },
    { "steps",      no_argument,            0, 's' },
    { "temperature", no_argument,           0, 'T' },
    { "verbose",    no_argument,            0, 'v' },
    { "version",    no_argument,            0, 'V' },
    { NULL,         0, 0, 0 }
  };

  /* strip any path to progname */
  progname = argv[0];
  if ((p = strrchr(progname, '/')) != NULL)
    progname = p + 1;

  /* set default box configuration */
  box = init_box(box);

  /* Read config file overrides default */
  int len = 3;
  path = reinterpret_cast<char *>(malloc(len));
  snprintf(path, len, "./");
  process_params(box, path);

  /* parse the command line options, they override all previous */
  while ((opt = getopt_long(argc, argv, "e:hl:s:T:Vv", longopts, NULL)) != -1) {
    switch (opt) {
    case 'e':
      box.set_eps(atof(optarg));
      break;
    case 'h':
      usage(EXIT_SUCCESS);
      break;
    case 'l':
      box.set_a(atof(optarg));
      break;
    case 's':
      box.set_steps(abs(atoi(optarg)));
      break;
    case 'T':
      box.set_temperature(atof(optarg));
      break;
    case 'V':
      printf("%s (%d)\n", progname, VERSION_MAJOR);
      exit(EXIT_SUCCESS);
    default:
      usage(EXIT_FAILURE);
    }
  }

  /* Output IC values */
  print_startup(box);
  mkdir_data();

  /* Initialize box */
  particles = initial_conditions(particles, &number, box);
  write_particles(particles, number);

  /* Compute stuff */
  rc = Evolve(particles, number, box);

  delete [] particles;
  free(path);
  return rc;
}
