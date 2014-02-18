/*
 *
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <string>

#include <boost/filesystem.hpp>

#include "include/configuration.h"
#include "include/experiment.h"
#include "include/parameters.h"
#include "include/macros.h"
#include "include/outputroutines.h"

/* build dependent variables */
#include "include/config.h"

char *progname;

static void usage(int rc) {
  printf("\nUsage: %s [option]\n\n", progname);
  printf("Calculate transport box\n"
         "  -h, --help              usage information\n"
         "  -i, --inputfile <file>  path to input configuration file\n"
         "  -c, --config <YAML>     specify config value overrides\n"
         "  -m, --modus <modus>     shortcut for -c 'General: { MODUS: <modus> }'\n"
         "  -s, --steps <steps>     shortcut for -c 'General: { STEPS: <modus> }'\n"
         "  -v, --version\n\n");
  exit(rc);
}


/* main - do command line parsing and hence decides modus */
int main(int argc, char *argv[]) {
  char *p, *path;
  int opt, rc = 0;
  struct option longopts[] = {
    { "config",     required_argument,      0, 'c' },
    { "help",       no_argument,            0, 'h' },
    { "inputfile",  required_argument,      0, 'i' },
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

  try {
    /* read in config file */
    Configuration configuration(path);

    /* check for overriding command line arguments */
    while ((opt = getopt_long(argc, argv, "c:hi:m:s:v", longopts,
                              nullptr)) != -1) {
      switch (opt) {
        case 'c':
          configuration.merge_yaml(optarg);
          break;
        case 'i': {
          const boost::filesystem::path file(optarg);
          configuration = Configuration(file.parent_path(), file.filename());
        } break;
        case 'h':
          usage(EXIT_SUCCESS);
          break;
        case 'm':
          configuration["General"]["MODUS"] = std::string(optarg);
          break;
        case 's':
          configuration["General"]["STEPS"] = abs(atoi(optarg));
          break;
        case 'v':
          exit(EXIT_SUCCESS);
        default:
          usage(EXIT_FAILURE);
      }
    }

    // TODO(mkretz): this needs to go into a file in the output
    printf("Configuration:\n%s\n", configuration.unused_values_report().c_str());

    mkdir_data();
    auto experiment = ExperimentBase::create(configuration);
    const std::string report = configuration.unused_values_report();
    if (!report.empty()) {
      printf("The following configuration values were not used:\n%s\n", report.c_str());
    }
    experiment->run(path);
  }
  catch(std::exception &e) {
    printf("SMASH failed with the following error:\n%s\n", e.what());
    rc = EXIT_FAILURE;
  }

  /* tear down */
  free(path);
  return rc;
}
