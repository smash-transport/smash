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

#include "include/DecayModes.h"
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
  DecayModes decay_modes;
  while ((characters_read = getline(&line, &line_size, file)) != -1) {
    printf("Retrieved decaymodes.txt line of length %li:\n", characters_read);
    /* Skip comments and blank lines */
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\t'
        || line[0] == '/') {
      printd("Skipping line: %s", line);
      continue;
    }
    printd("line: %s", line);

     /* First item in a line should be the PDG code */
    characters = strtok_r(line, " ", &line_position);
    if (characters == NULL)
      continue;
    int pdgcode = atoi(characters);
    printd("pdgcode: %d\n", pdgcode);
    /* Second item should be the particle name */
    char *particle_name = strtok_r(NULL, " ", &line_position);
    if (particle_name == NULL)
      continue;
    printd("particle name: %s\n", particle_name);
    /* Rest of the line should contain the decay modes.
     * the list of decay particles (PDG codes) is marked with special separator
     * at the beginning and end of the list
     */
    characters = strtok_r(NULL, &mode_separator, &line_position);
    printd("characters: %s \n", characters);
    while (characters != NULL) {
      if (strlen(characters) > 2) {
         /* Read in the decay particle PDG codes for this mode */
        pdgs = strtok_r(characters, " ", &list_position);
        printd("pdgs: %s \n", pdgs);
        while (pdgs != NULL) {
          int decay_particle = atoi(pdgs);
          printd("decay particle: %i \n", decay_particle);
          decay_particles.push_back(decay_particle);
          pdgs = strtok_r(NULL, " ", &list_position);
          printd("pdgs: %s \n", pdgs);
        }
        /* Finished reading the particles; this should be
         * followed by the ratio of this decay mode
         * among all the possible modes
         */
        characters = strtok_r(NULL, " ", &line_position);
        printd("characters: %s \n", characters);
        float ratio = atof(characters);
        printd("ratio: %f\n", ratio);
        /* Add mode to the list of possible decays for this particle */
        decay_modes.add_mode(decay_particles, ratio);
        /* Clean the particle list for the next mode */
        decay_particles.clear();
      }
      characters = strtok_r(NULL, &mode_separator, &line_position);
      printd("characters: %s \n", characters);
    }
    /* Add the list of decay modes for this particle type */
    particles->add_decaymodes(decay_modes, pdgcode);
    /* Clean up the list for the next particle type */
    decay_modes.clear();
  }
  free(line);
  fclose(file);
  printf("Finished reading decaymodes.txt\n");
  return;
}
