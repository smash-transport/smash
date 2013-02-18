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
         "  -T, --temperature    initial temperature\n"
         "  -v, --verbose        show debug info\n"
         "  -V, --version\n\n");
  exit(rc);
}

static int Evolve(ParticleData *particles, int number, box box) {
  FourVector distance;

  for (int steps = 0; steps < box.steps(); steps++)
    for (int i = 0; i < number; i++) {
       /* XXX: interactions */
       /* calculate particles motion */
       distance.set_FourVector(1.0, particles[i].velocity_x(),
         particles[i].velocity_y(), particles[i].velocity_z());
       distance *= box.eps();
       printd("Particle %d motion: %g %g %g %g\n", particles[i].id(),
         distance.x0(), distance.x1(), distance.x2(), distance.x3());

       /* XXX : treat properly boundaries */
       particles[i].add_position(distance);
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
    { "eps",    no_argument,                0, 'e' },
    { "help",       no_argument,            0, 'h' },
    { "length",    no_argument,             0, 'l' },
    { "temperature", no_argument,           0, 'T' },
    { "verbose",    no_argument,            0, 'v' },
    { "version",    no_argument,            0, 'V' },
    { NULL,         0, 0, 0 }
  };

  /* strip any path to progname */
  progname = argv[0];
  if ((p = strrchr(progname, '/')) != NULL)
    progname = p + 1;

  /* parse the command line options */
  while ((opt = getopt_long(argc, argv, "e:hl:T:Vv", longopts, NULL)) != -1) {
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

  /* Read config file */
  int len = 3;
  path = reinterpret_cast<char *>(malloc(len));
  snprintf(path, len, "./");
  process_params(box, path);

  /* Output IC values */
  box = init_box(box);
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
