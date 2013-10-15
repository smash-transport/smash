/*
 *
 *    Copyright (c) 2012
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <list>

#include "include/param-reader.h"
#include "include/Box.h"
#include "include/Laboratory.h"
#include "include/Parameters.h"
#include "include/outputroutines.h"

/* XXX: hardcoded length cap */
#define FILELEN 256

/* space separation between items */
const char *sep = " \t\n";

/* process_params - read in all params into vetor */
void process_params(char *path, std::list<Parameters> *configuration) {
  char *line = NULL, *saveptr = NULL, params[FILELEN];
  size_t len = 0;
  ssize_t read;
  char *key, *value;
  FILE *fp;

  /* Looking for parameters in config file
   * If none exists, we'll use default values.
   */
  snprintf(params, strlen(path) + 12, "%s/params.txt", path);
  fp = fopen(params, "r");
  if (!fp) {
    fprintf(stderr, "W: No params.txt at %s path.\n", path);
    return;
  }

  printf("Processing %s/params.txt.\n", path);

  while ((read = getline(&line, &len, fp)) != -1) {
    printd("Retrieved params.txt line of length %li :\n", read);
    printd("%s", line);
    /* Skip comments and blank lines */
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\t' || line[0] == '/')
      continue;

    key = strtok_r(line, sep, &saveptr);
    if (key == NULL)
      continue;
    value = strtok_r(NULL, sep, &saveptr);
    if (value == NULL)
      continue;

    Parameters parameter(key, value);
    configuration->push_back(parameter);
    printd("%s %s\n", configuration->back().key(),
           configuration->back().value());
  }
  free(line);
  fclose(fp);
  printd("Read all params.txt.\n");
}

void assign_params(std::list<Parameters> *configuration,
  Laboratory *parameters) {
  bool match = false;
  for (std::list<Parameters>::iterator i = configuration->begin();
       i != configuration->end(); ++i) {
    char *key = i->key();
    char *value = i->value();
    printd("%s %s\n", key, value);

    /* integer values */
    if (strcmp(key, "STEPS") == 0) {
      parameters->set_steps(abs(atoi(value)));
      match = true;
    }
    if (strcmp(key, "RANDOMSEED") == 0) {
      /* negative seed means random startup value */
      if (atol(value) > 0)
        parameters->set_seed(atol(value));
      else
        parameters->set_seed(time(NULL));
      match = true;
    }
    if (strcmp(key, "UPDATE") == 0) {
      parameters->set_output_interval(abs(atoi(value)));
      match = true;
    }
    if (strcmp(key, "TESTPARTICLES") == 0) {
      parameters->set_testparticles(abs(atoi(value)));
      match = true;
    }
    if (strcmp(key, "INITIAL_CONDITION") == 0) {
      parameters->set_initial_condition(abs(atoi(value)));
      match = true;
    }


    /* double or float values */
    if (strcmp(key, "EPS") == 0) {
      parameters->set_eps(fabs(atof(value)));
      match = true;
    }
    if (strcmp(key, "SIGMA") == 0) {
      parameters->set_cross_section(fabs(atof(value)));
      match = true;
    }

    /* remove processed entry */
    if (match) {
      i = configuration->erase(i);
      match = false;
    }
  }
}

void assign_params(std::list<Parameters> *configuration, Box *box) {
  bool match = false;
  for (std::list<Parameters>::iterator i = configuration->begin();
       i != configuration->end(); ++i) {
    char *key = i->key();
    char *value = i->value();
    printd("%s %s\n", key, value);

    /* double or float values */
    if (strcmp(key, "LENGTH") == 0) {
      box->set_length(fabs(atof(value)));
      match = true;
    }
    if (strcmp(key, "TEMPERATURE") == 0) {
      box->set_temperature(fabs(atof(value)));
      match = true;
    }
    /* remove processed entry */
    if (match) {
      i = configuration->erase(i);
      match = false;
    }
  }
}
