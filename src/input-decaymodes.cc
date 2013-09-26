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
#include <utility>
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
  snprintf(input_decaymodes, strlen(path) + 16, "%s/decaymodes.txt", path);
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
  const char mode_separator = ';';
  std::vector<int> decay_particles;
  std::pair<std::vector<int>, float> decay_mode;
  std::vector< std::pair<std::vector<int>, float> > decay_modes;
  while ((characters_read = getline(&line, &line_size, file)) != -1) {
    printf("Retrieved decaymodes.txt line of length %li:\n", characters_read);
    /* Skip comments and blank lines */
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\t'
        || line[0] == '/') {
      printf("Skipping line: %s", line);
      continue;
    }
    printf("line: %s", line);

    characters = strtok_r(line, " ", &line_position);
    if (characters == NULL)
      continue;
    int pdgcode = atoi(characters);
    printf("pdgcode: %d\n", pdgcode);
    char *particle_name = strtok_r(NULL, " ", &line_position);
    if (particle_name == NULL)
      continue;
    printf("particle name: %s\n", particle_name);
    characters = strtok_r(NULL, &mode_separator, &line_position);
    printf("characters: %s \n", characters);
    while (characters != NULL) {
      if (strlen(characters) > 2) {
        pdgs = strtok_r(characters, " ", &list_position);
        printf("pdgs: %s \n", pdgs);
        while (pdgs != NULL) {
          int decay_particle = atoi(pdgs);
          printf("decay particle: %i \n", decay_particle);
          decay_particles.push_back(decay_particle);
          pdgs = strtok_r(NULL, " ", &list_position);
          printf("pdgs: %s \n", pdgs);
        }
        characters = strtok_r(NULL, " ", &line_position);
        printf("characters: %s \n", characters);
        float ratio = atof(characters);
        printf("ratio: %f\n", ratio);
        decay_mode = std::make_pair(decay_particles, ratio);
        decay_modes.push_back(decay_mode);
      }
      characters = strtok_r(NULL, &mode_separator, &line_position);
      printf("characters: %s \n", characters);
    }
    //particles->add_decay_modes(pdgcode, decay_modes);
  }
  free(line);
  fclose(file);
  printf("Finished reading decaymodes.txt\n");
  return;
}
