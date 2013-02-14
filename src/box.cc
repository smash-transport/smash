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

/* Side length of the cube in fm */
float A = 5.0;

/* Time steps */
float EPS = 0.1;

/* Total number of steps */
int STEPS = 10000;

/* Temperature of the Boltzmann distribution for thermal initialization */
float temperature = 0.1;

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

static int Evolve(ParticleData *particles, int number) {
  FourVector distance;	

  for (int steps = 0; steps < STEPS; steps++)
    for (int i = 0; i < number; i++) {
       /* XXX: interactions */
       /* calculate particles motion */
       distance.set_FourVector(1.0, particles[i].velocity_x(),	    
         particles[i].velocity_y(), particles[i].velocity_z());
       distance *= EPS;
       printd("Particle %d motion: %g %g %g %g\n", particles[i].id(),
         distance.x0(), distance.x1(), distance.x2(), distance.x3());

       /* XXX : treat properly boundaries */
       particles[i].add_position(distance);
       printd("Particle %d position: %g %g %g %g\n", particles[i].id(),
                particles[i].x().x0(), particles[i].x().x1(),
                particles[i].x().x2(), particles[i].x().x3());
    }
  return 0;
}

int main(int argc, char *argv[]) {
  char *p, *path;
  int opt, rc, number = 0;
  ParticleData *particles = NULL;

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
      EPS = atof(optarg);
      break;
    case 'h':
      usage(EXIT_SUCCESS);
      break;
    case 'l':
      A = atof(optarg);
      break;
    case 'T':
      temperature = atof(optarg);
      break;
    case 'V':
      printf("%s (%d.%d)\n", progname, VERSION_MAJOR, VERSION_MINOR);
      exit(EXIT_SUCCESS);
    default:
      usage(EXIT_FAILURE);
    }
  }

  /* Read config file */
  int len = 3;
  path = reinterpret_cast<char *>(malloc(len));
  snprintf(path, len, "./");
  process_params(path);

  /* Output IC values */
  print_startup();

  /* Initialize box */
  particles = initial_conditions(particles, number);
  print_particles(particles, number);

  /* Compute stuff */
  rc = Evolve(particles, number);

  delete [] particles;
  free(path);
  return rc;
}
