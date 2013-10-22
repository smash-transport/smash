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
#include "include/Sphere.h"
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
         "  -m, --modus          modus of laboratory\n"
         "  -O, --output-interval          step interval between measurements\n"
         "  -R, --random         random number seed\n"
         "  -s, --sigma          cross section in mbarn\n"
         "  -S, --steps          number of steps\n"
         "  -V, --version\n\n");
  exit(rc);
}

/* main - do command line parsing and hence decides modus */
int main(int argc, char *argv[]) {
  char *p, *path;
  int opt, rc = 0;

  struct option longopts[] = {
    { "eps",        required_argument,      0, 'e' },
    { "help",       no_argument,            0, 'h' },
    { "modus",      required_argument,      0, 'm' },
    { "output-interval", required_argument,      0, 'O' },
    { "random",     required_argument,      0, 'R' },
    { "sigma",      required_argument,      0, 's' },
    { "steps",      required_argument,      0, 'S' },
    { "version",    no_argument,            0, 'V' },
    { NULL,         0, 0, 0 }
  };

  /* strip any path to progname */
  progname = argv[0];
  if ((p = strrchr(progname, '/')) != NULL)
    progname = p + 1;
  printf("%s (%d)\n", progname, VERSION_MAJOR);

  /* XXX: make path configurable */
  size_t len = strlen("./") + 1;
  path = reinterpret_cast<char *>(malloc(len));
  snprintf(path, len, "./");

  /* Read Laboratory config file parameters */
  Laboratory *lab = new Laboratory();
  process_config_laboratory(lab, path);

  /* parse the command line options, they override all previous */
  while ((opt = getopt_long(argc, argv, "e:hm:O:R:s:S:V", longopts,
    NULL)) != -1) {
    switch (opt) {
    case 'e':
      lab->set_eps(fabs(atof(optarg)));
      break;
    case 'h':
      usage(EXIT_SUCCESS);
      break;
    case 'm':
      lab->set_modus(fabs(atoi(optarg)));
      break;
    case 'O':
      {
      /* guard output_interval to be positive and greater null */
      int output_interval = abs(atoi(optarg));
      if (output_interval > 0)
        lab->set_output_interval(output_interval);
      }
      break;
    case 'R':
      /* negative seed is for time */
      if (atol(optarg) > 0)
        lab->set_seed(atol(optarg));
      else
        lab->set_seed(time(NULL));
      break;
    case 's':
      lab->set_cross_section(fabs(atof(optarg)));
      break;
    case 'S':
      lab->set_steps(abs(atoi(optarg)));
      break;
    case 'V':
      exit(EXIT_SUCCESS);
    default:
      usage(EXIT_FAILURE);
    }
  }

  /* Output IC values */
  print_startup(*lab);
  mkdir_data();
  write_oscar_header();

  /* initialize random seed */
  srand48(lab->seed());

  /* reducing cross section according to number of test particle */
  if (lab->testparticles() > 1) {
    printf("IC test particle: %i\n", lab->testparticles());
    lab->set_cross_section(lab->cross_section() / lab->testparticles());
    printf("Elastic cross section: %g [mb]\n", lab->cross_section());
  }

  /* modus operandi */
  if (lab->modus() == 1) {
    Box *box = new Box(*lab);
    lab = box;
  } else if (lab->modus() == 2) {
    Sphere *ball = new Sphere(*lab);
    lab = ball;
  }

  /* read specific config files */
  lab->process_config(path);

  /* initiate particles */
  Particles *particles = new Particles;
  input_particles(particles, path);
  input_decaymodes(particles, path);
  CrossSections *cross_sections = new CrossSections;
  cross_sections->add_elastic_parameter(lab->cross_section());
  lab->initial_conditions(particles);

  /* record IC startup */
  print_startup(*lab);
  write_measurements_header(*particles);
  print_header();
  write_particles(*particles);

  /* the time evolution of the relevant subsystem */
  rc = lab->evolve(particles, cross_sections);

  /* tear down */
  free(path);
  delete particles;
  delete cross_sections;
  return rc;
}
