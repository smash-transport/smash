/*
 *
 *    Copyright (c) 2012
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
#include "include/param-reader.h"

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

static void usage(int rc) {
  printf("\nUsage: %s [option]\n\n", progname);
  printf("Calculate transport box\n"
         "  -h, --help           usage information\n"
         "  -v, --verbose        show debug info\n"
         "  -V, --version\n\n");
  exit(rc);
}

static int Evolve(void) {
  /* do something */
  for (int i = 0; i < STEPS; i++)
    continue;
  return 0;
}

int main(int argc, char *argv[]) {
  char *p, *path;
  int opt, rc;

  struct option longopts[] = {
    { "help",       no_argument,            0, 'h' },
    { "verbose",    no_argument,            0, 'v' },
    { "version",    no_argument,            0, 'V' },
    { NULL,         0, 0, 0 }
  };

  /* strip any path to progname */
  progname = argv[0];
  if ((p = strrchr(progname, '/')) != NULL)
    progname = p + 1;

  /* parse the command line options */
  while ((opt = getopt_long(argc, argv, "hVv", longopts, NULL)) != -1) {
    switch (opt) {
    case 'h':
      usage(EXIT_SUCCESS);
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

  /* Compute stuff */
  rc = Evolve();

  free(path);
  return rc;
}
