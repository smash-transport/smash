/*
 *
 *    Copyright (c) 2012-2014
 *      SMASH Team
 * 
 *    GNU General Public License (GPLv3 or later)
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


#include "include/Experiment.h"
#include "include/Parameters.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/param-reader.h"

/* build dependent variables */
#include "include/Config.h"

char *progname;

static void usage(int rc) {
  printf("\nUsage: %s [option]\n\n", progname);
  printf("Calculate transport box\n"
         "  -h, --help           usage information\n"
         "  -m, --modus          modus of laboratory\n"
         "  -s, --steps          number of steps\n"
         "  -v, --version\n\n");
  exit(rc);
}


/* main - do command line parsing and hence decides modus */
int main(int argc, char *argv[]) {
  char *p, *path;
  int opt, rc = 0;
  int steps = 0, nevents = 0;
  struct option longopts[] = {
    { "help",       no_argument,            0, 'h' },
    { "modus",      required_argument,      0, 'm' },
    { "steps",      required_argument,      0, 's' },
    { "version",    no_argument,            0, 'v' },
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

  /* read in config file */
  std::list<Parameters> configuration;
  process_config(&configuration, path);
  bool match = false;
  std::string modus_chooser;
  std::list<Parameters>::iterator i = configuration.begin();
  while (i != configuration.end()) {
    char *key = i->key();
    char *value = i->value();
    printd("Looking for match %s %s\n", key, value);
    /* integer values */
    if (strcmp(key, "MODUS") == 0) {
      modus_chooser = value;
      match = true;
    }
    if (strcmp(key, "NEVENTS") == 0) {
      nevents = abs(atoi(value));
      match = true;
    }
    /* remove processed entry */
    if (match) {
      printd("Erasing %s %s\n", key, value);
      i = configuration.erase(i);
      match = false;
    } else {
      ++i;
    }
  }
  /* check for overriding command line arguments */
  while ((opt = getopt_long(argc, argv, "hm:s:v", longopts,
                              NULL)) != -1) {
    switch (opt) {
      case 'h':
        usage(EXIT_SUCCESS);
        break;
      case 'm':
        modus_chooser = optarg;
        printf("Modus set: %s\n", modus_chooser.c_str());
        break;
//    case 'r':
/* negative seed is for time */
//      if (atol(optarg) > 0)
//         lab->set_seed(atol(optarg));
//      else
//         lab->set_seed(time(NULL));
//      break;
      case 's':
        steps = abs(atoi(optarg));
        break;
      case 'v':
        exit(EXIT_SUCCESS);
      default:
        usage(EXIT_FAILURE);
    }
  }
  printf("Modus for this calculation: %s\n", modus_chooser.c_str());
  auto experiment = Experiment::create(modus_chooser);
  experiment->configure(configuration);
  mkdir_data();
  /* overriden arguments */
  if (steps > 0)
    experiment->commandline_arg(steps);

  for (int j = 1; j < nevents; j++) {
    write_oscar_header();
    experiment->initialize(path);
    /* the time evolution of the relevant subsystem */
    experiment->run_time_evolution();
  }

  /* tear down */
  free(path);
  experiment->end();
  return rc;
}
