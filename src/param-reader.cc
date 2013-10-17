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
#include "include/Sphere.h"
#include "include/outputroutines.h"

/* space separation between items */
const char *sep = " \t\n";

/* process_params - read in all params into list of key value paris */
static void process_params(char *file_path,
  std::list<Parameters> *configuration) {
  char *line = NULL, *saveptr = NULL;
  size_t len = 0;
  ssize_t read;
  char *key, *value;
  FILE *fp;

  /* Looking for parameters in config file
   * If none exists, we'll use default values.
   */
  fp = fopen(file_path, "r");
  if (!fp) {
    fprintf(stderr, "W: No config at %s path.\n", file_path);
    return;
  }

  printf("Processing config %s.\n", file_path);

  while ((read = getline(&line, &len, fp)) != -1) {
    printd("Retrieved %s line of length %li :\n", file_path, read);
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
  printd("Read all %s.\n", file_path);
}

static void assign_params(std::list<Parameters> *configuration,
  Laboratory *parameters) {
  bool match = false;
  std::list<Parameters>::iterator i = configuration->begin();
  while (i != configuration->end()) {
    char *key = i->key();
    char *value = i->value();
    printd("Looking for match %s %s\n", key, value);

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
    if (strcmp(key, "MODUS") == 0) {
      parameters->set_modus(abs(atoi(value)));
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
      printd("Erasing %s %s\n", key, value);
      i = configuration->erase(i);
      match = false;
    } else {
      ++i;
    }
  }
}

static void assign_params(std::list<Parameters> *configuration, Box *box) {
  bool match = false;
  std::list<Parameters>::iterator i = configuration->begin();
  while (i != configuration->end()) {
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
    /* int values */
    if (strcmp(key, "INITIAL_CONDITION") == 0) {
      box->set_initial_condition(abs(atoi(value)));
      match = true;
    }
    /* remove processed entry */
    if (match) {
      i = configuration->erase(i);
      match = false;
    } else {
      ++i;
    }
  }
}

static void assign_params(std::list<Parameters> *configuration, Sphere *ball) {
  bool match = false;
  std::list<Parameters>::iterator i = configuration->begin();
  while (i != configuration->end()) {
    char *key = i->key();
    char *value = i->value();
    printd("%s %s\n", key, value);

    /* double or float values */
    if (strcmp(key, "RADIUS") == 0) {
      ball->set_radius(fabs(atof(value)));
      match = true;
    }

    /* remove processed entry */
    if (match) {
      i = configuration->erase(i);
      match = false;
    } else {
      ++i;
    }
  }
}

/* process_laboratory_config - configuration handling */
void process_config_laboratory(Laboratory *parameters, char *path) {
  std::list<Parameters> configuration;
  size_t len = strlen("./config_laboratory.txt") + strlen(path) + 1;
  char *config_path = reinterpret_cast<char *>(malloc(len));
  snprintf(config_path, len, "%s/config_laboratory.txt", path);
  process_params(config_path, &configuration);
  assign_params(&configuration, parameters);
  warn_wrong_params(&configuration);
  free(config_path);
}

/* process_box_config - configuration handling */
void process_config_box(Box *cube, char *path) {
  std::list<Parameters> configuration;
  size_t len = strlen("./config_box.txt") + strlen(path) + 1;
  char *config_path = reinterpret_cast<char *>(malloc(len));
  snprintf(config_path, len, "%s/config_box.txt", path);
  process_params(config_path, &configuration);
  assign_params(&configuration, cube);
  warn_wrong_params(&configuration);
  free(config_path);
}

/* process_sphere_config - configuration handling */
void process_config_sphere(Sphere *ball, char *path) {
  std::list<Parameters> configuration;
  size_t len = strlen("./config_sphere.txt") + strlen(path) + 1;
  char *config_path = reinterpret_cast<char *>(malloc(len));
  snprintf(config_path, len, "%s/config_sphere.txt", path);
  process_params(config_path, &configuration);
  assign_params(&configuration, ball);
  warn_wrong_params(&configuration);
  free(config_path);
}
