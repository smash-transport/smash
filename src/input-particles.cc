/*
 *
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "include/input-particles.h"
#include "include/initial-conditions.h"
#include "include/ParticleType.h"
#include "include/outputroutines.h"

/* XXX: hardcoded length cap */
#define FILELEN 256

/* input_particles - read in particle types */
void input_particles(std::vector<ParticleType> *type, char *path) {
  char *line = NULL, *saveptr = NULL, *characters, input_particles[FILELEN];
  size_t len = 0, type_number = 0;
  ssize_t read;
  FILE *fp;

  /* Looking for parameters in config file
   * If none exists, we'll use default values.
   */
  snprintf(input_particles, strlen(path) + 15, "%s/particles.txt", path);
  fp = fopen(input_particles, "r");
  if (!fp) {
    fprintf(stderr, "W: No particles.txt at %s path.\n", path);
    /* use just pions in that case */
    initial_particles(type);
    return;
  }

  printf("Processing %s/particles.txt.\n", path);

  while ((read = getline(&line, &len, fp)) != -1) {
    printd("Retrieved particles.txt line of length %zu :\n", read);
    printd("%s", line);
    /* Skip comments and blank lines */
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\t' || line[0] == '/')
      continue;

    char *particle_name = strtok_r(line, sep, &saveptr);
    if (particle_name == NULL)
      continue;
    characters = strtok_r(NULL, sep, &saveptr);
    if (characters == NULL)
      continue;
    float mass = atof(characters);
    characters = strtok_r(NULL, sep, &saveptr);
    if (characters == NULL)
      continue;
    float width = atof(characters);
    characters = strtok_r(NULL, sep, &saveptr);
    if (characters == NULL)
      continue;
    int pdgcode = atoi(characters);
    characters = strtok_r(NULL, sep, &saveptr);
    if (characters == NULL)
      continue;
    int isospin = atoi(characters);
    characters = strtok_r(NULL, sep, &saveptr);
    int charge;
    if (characters != NULL)
      charge = atoi(characters);
    else
      charge = 0;

    /* Have a new particle type */
    (*type).resize(type_number + 1);
    printf("Setting particle type %s mass %g width %g"
           " pdgcode %i isospin %i charge %i \n", particle_name,
           mass, width, pdgcode, isospin, charge);
    std::string name(particle_name);
    (*type)[type_number].set(name, mass, width, pdgcode, isospin, charge);
    type_number++;
  }
  free(line);
  fclose(fp);
  return;
}
