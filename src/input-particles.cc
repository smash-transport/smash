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

#include "include/input-particles.h"
#include "include/initial-conditions.h"
#include "include/Particles.h"
#include "include/outputroutines.h"

/* XXX: hardcoded length cap */
#define FILELEN 256

/* input_particles - read in particle types */
void input_particles(Particles *particles, char *path) {
  char *line = NULL, *saveptr = NULL, *characters, input_particles[FILELEN];
  size_t len = 0;
  ssize_t read;
  FILE *fp = NULL;

  /* Looking for parameters in config file
   * If none exists, we'll use default values.
   */
  snprintf(input_particles, strlen(path) + 15, "%s/particles.txt", path);
  fp = fopen(input_particles, "r");
  if (!fp) {
    fprintf(stderr, "E: No particles.txt at %s path.\n", path);
    exit(EXIT_FAILURE);
  }

  printf("Processing %s/particles.txt.\n", path);

  while ((read = getline(&line, &len, fp)) != -1) {
    printd("Retrieved particles.txt line of length %li:\n", read);
    /* Skip comments and blank lines */
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\t'
        || line[0] == '/') {
      printd("Skipping line: %s", line);
      continue;
    }
    printd("line: %s", line);

    char *particle_name = strtok_r(line, sep, &saveptr);
    if (particle_name == NULL)
      continue;
    printd("particle name: %s\n", particle_name);
    characters = strtok_r(NULL, sep, &saveptr);
    if (characters == NULL)
      continue;
    float mass = atof(characters);
    printd("mass: %f\n", mass);
    characters = strtok_r(NULL, sep, &saveptr);
    if (characters == NULL)
      continue;
    float width = atof(characters);
    printd("width: %f\n", width);
    characters = strtok_r(NULL, sep, &saveptr);
    if (characters == NULL)
      continue;
    int pdgcode = atoi(characters);
    printd("pdgcode: %d\n", pdgcode);
    characters = strtok_r(NULL, sep, &saveptr);
    if (characters == NULL)
      continue;
    int isospin = atoi(characters);
    printd("isospin: %d\n", isospin);
    characters = strtok_r(NULL, sep, &saveptr);
    if (characters == NULL)
      continue;
    int charge = atoi(characters);
    printd("charge: %d\n", charge);
    characters = strtok_r(NULL, sep, &saveptr);
    if (characters == NULL)
      continue;
    int spin = atoi(characters);
    printd("spin: %d\n", spin);

    /* Have a new particle type */
    printf("Setting particle type %s mass %g width %g pdgcode %i\n",
           particle_name, mass, width, pdgcode);
    printf("Setting particle type %s isospin %i charge %i spin %i\n",
           particle_name, isospin, charge, spin);
    std::string name(particle_name);
    ParticleType new_type(name, mass, width, pdgcode, isospin, charge, spin);
    particles->add_type(new_type, pdgcode);
  }
  free(line);
  fclose(fp);
  printd("Finished reading particles.txt\n");

  if (particles->types_empty()) {
    fprintf(stderr, "E: Empty particles.txt at %s path.\n", path);
    fprintf(stderr, "E: No particle types.\n");
    exit(EXIT_FAILURE);
  }

  return;
}
