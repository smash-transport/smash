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

#include "include/input-decaymodes.h"
#include "include/Particles.h"
#include "include/outputroutines.h"

/* XXX: hardcoded length cap */
#define FILELEN 256

/* input_decaymodes - read in particle decay modes */
void input_decaymodes(Particles *particles, char *path) {
  char input_decaymodes[FILELEN];
  FILE *file;

  /* Looking for decay mode list
   * If none exists, we'll use default values.
   */
  snprintf(input_decaymodes, strlen(path) + 15, "%s/decaymodes.txt", path);
  file = fopen(input_decaymodes, "r");
  if (!file) {
    fprintf(stderr, "W: No decaymodes.txt at %s path.\n", path);
    return;
  }

  printf("Processing %s/decaymodes.txt.\n", path);

  char *line = NULL, *line_position = NULL, *list_position = NULL;
  char *characters, *pdgs;
  size_t line_size = 0;
  ssize_t characters_read;
  const char mode_separator = ";";
  std::vector<int> decay_particles;
  std::pair<vector<int>, float> decay_mode;
  std::vector< std::pair<vector<int>, float> > decay_modes;
  while ((characters_read = getline(&line, &line_size, file)) != -1) {
    printd("Retrieved particles.txt line of length %li:\n", characters_read);
    /* Skip comments and blank lines */
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\t'
        || line[0] == '/') {
      printd("Skipping line: %s", line);
      continue;
    }
    printd("line: %s", line);

    characters = strtok_r(line, " ", &line_position);
    if (characters == NULL)
      continue;
    int pdgcode = atoi(characters);
    printd("pdgcode: %d\n", pdgcode);
    char *particle_name = strtok_r(NULL, " ", &line_position);
    if (particle_name == NULL)
      continue;
    printd("particle name: %s\n", particle_name);
    characters = strtok_r(NULL, &mode_separator, &line_position);
    while (characters != NULL) {
      pdgs = strtok_r(characters, " ", &list_position);
      while (pdgs != NULL) {
        int decay_particle = atoi(pdgs);
        printd("decay particle: %f ", decay_particle);
        decay_particles.push_back(decay_particle);
        pdgs = strtok_r(NULL, " ", &list_position);
      }
      characters = strtok_r(NULL, &mode_separator, &line_position);
      float ratio = atof(characters);
      printd("ratio: %f\n", ratio);
      decay_mode = std::make_pair(decay_particles, ratio);
      decay_modes.push_back(decay_mode);
    }
    particles->add_decay_modes(pdgcode, decay_modes);
  }
  free(line);
  fclose(file);
  printd("Finished reading decaymodes.txt\n");
  return;
}
