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

#include "include/Parameters.h"
#include "include/macros.h"
#include "include/param-reader.h"
#include "include/outputroutines.h"
#include "include/Experiment.h"
/* build dependent variables */
#include "include/Config.h"

char *progname;

static void usage(int rc) {
  printf("\nUsage: %s [option]\n\n", progname);
  printf("Calculate transport box\n"
         "  -h, --help           usage information\n"
         "  -m, --modus          modus of laboratory\n"
         "  -V, --version\n\n");
  exit(rc);
}


/* main - do command line parsing and hence decides modus */
int main(int argc, char *argv[]) {
  char *p, *path;
  int opt, rc = 0;
  char modus[20];
    
  struct option longopts[] = {
    { "help",       no_argument,            0, 'h' },
    { "modus",      required_argument,      0, 'm' },
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

// read in config file
    
    std::list<Parameters> configuration;
    process_config(&configuration, path);

    bool match = false;
    std::list<Parameters>::iterator i = configuration.begin();
    while (i != configuration.end()) {
        char *key = i->key();
        char *value = i->value();
        printd("Looking for match %s %s\n", key, value);
        
        /* integer values */
        if (strcmp(key, "MODUS") == 0) {
            strncpy(modus,value, sizeof(&modus));
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

    
    while ((opt = getopt_long(argc, argv, "hm:V", longopts,
                              NULL)) != -1) {
        switch (opt) {
            case 'h':
                usage(EXIT_SUCCESS);
                break;
            case 'm':
                strncpy(modus,optarg, sizeof(modus));
                printf("Modus read in: %s \n", modus);
                break;
                //    case 'R':
                /* negative seed is for time */
                //      if (atol(optarg) > 0)
                //         lab->set_seed(atol(optarg));
                //      else
                //         lab->set_seed(time(NULL));
                //      break;
                //    case 'S':
                //            lab->set_steps(abs(atoi(optarg)));
                //      break;
            case 'V':
                exit(EXIT_SUCCESS);
            default:
                usage(EXIT_FAILURE);
        }
    }

    printf("Modus for this calculation: %s \n", modus);
    
    auto experiment = Experiment::create(modus);
      
    experiment->config(configuration);

    mkdir_data();
    write_oscar_header();
    
    experiment->initialize(path);

  /* the time evolution of the relevant subsystem */
//  rc = lab->evolve(particles, cross_sections);

  /* tear down */
    free(path);
    experiment->end();
 return rc;
      
 
}
