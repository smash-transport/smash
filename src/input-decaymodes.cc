/*
 *
 *    Copyright (c) 2013
 *      SMASH Team
 * 
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "include/constants.h"
#include "include/DecayModes.h"
#include "include/input-decaymodes.h"
#include "include/Particles.h"
#include "include/ParticleType.h"
#include "include/outputroutines.h"

/* XXX: hardcoded length cap */
#define FILELEN 256

/* input_decaymodes - read in particle decay modes */
void input_decaymodes(Particles *particles, char *path) {
  char input_decaymodes[FILELEN];
  FILE *file = NULL;

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

  char *line = NULL, *line_position = NULL;
  char *characters, *pdgs;
  size_t line_size = 0;
  ssize_t characters_read;
  std::vector<int> decay_particles;
  DecayModes decay_modes;
  int modes = 0, pdgcode = 0;
  float ratio_sum = 0.0;
  bool unknown_particle = true;
  while ((characters_read = getline(&line, &line_size, file)) != -1) {
    printd("Retrieved decaymodes.txt line of length %li:\n", characters_read);
    /* Skip comments and blank lines */
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\t'
        || line[0] == '/' || characters_read == 0) {
      printd("Skipping line: %s", line);
      continue;
    }
    printd("line: %s", line);
    /* Check if we are reading in the decay modes */
    if (modes == 0) {
      /* Not reading modes; the first item in a line should be the PDG code */
      characters = strtok_r(line, " ", &line_position);
      if (characters == NULL)
        continue;
      pdgcode = atoi(characters);
      printd("pdgcode: %d\n", pdgcode);
      /* Check if SMASH knows this particle */
      unknown_particle = true;
      std::map<int, ParticleType>::const_iterator type_i
          = particles->types_cbegin();
      while (type_i != particles->types_cend()) {
        if (type_i->second.pdgcode() == pdgcode) {
            unknown_particle = false;
            break;
        }
        ++type_i;
      }
      /* Second item should be the number of decay modes */
      characters = strtok_r(NULL, " ", &line_position);
      if (characters == NULL)
        continue;
      modes = atoi(characters);
      /* Initialize the sum of the ratios */
      ratio_sum = 0.0;
      printd("Number of modes: %i\n", modes);
      if (unknown_particle) {
        printf("Warning: Could not find particle %i in the list!\n",
               pdgcode);
        printf("Can't add decay modes for particle %i.\n", pdgcode);
        printf("The following mode information will be ignored.\n");
        continue;
      }
    } else {
      if (unknown_particle) {
        modes--;
        continue;
      }
      /* Reading the list of the decay modes */
      characters = strtok_r(line, " ", &line_position);
      /* First number should be the ratio for this mode */
      float ratio = atof(characters);
      printd("Ratio: %g \n", ratio);
      /* Read in the decay particle PDG codes for this mode */
      pdgs = strtok_r(NULL, " ", &line_position);
      printd("pdgs: %s \n", pdgs);
      while (pdgs != NULL) {
        if (strlen(pdgs) > 2) {
          int decay_particle = atoi(pdgs);
          printd("decay particle: %i \n", decay_particle);
          decay_particles.push_back(decay_particle);
        }
        pdgs = strtok_r(NULL, " ", &line_position);
        printd("pdgs: %s \n", pdgs);
      }
      /* Check that the particles in the decay mode are known */
      bool unknown_decay = false;
      for (std::vector<int>::const_iterator decay_pdg
           = decay_particles.begin(); decay_pdg != decay_particles.end();
           ++decay_pdg) {
        bool not_found = true;
        std::map<int, ParticleType>::const_iterator type_i
          = particles->types_cbegin();
        while (type_i != particles->types_cend()) {
          if (type_i->second.pdgcode() == *decay_pdg) {
            not_found = false;
            break;
          }
          ++type_i;
        }
        if (not_found) {
          unknown_decay = true;
          printf("Warning: Could not find particle %i in the list!\n",
                 *decay_pdg);
          printf("Can't add this decay mode for particle %i.\n", pdgcode);
        }
      }
      if (!unknown_decay) {
        /* Add mode to the list of possible decays for this particle */
        decay_modes.add_mode(decay_particles, ratio);
        ratio_sum += ratio;
      }
      /* Clean the particle list for the next mode */
      decay_particles.clear();
      /* Number of modes remaining decreases by one */
      modes--;
      /* Check if this was the last mode */
      if (modes == 0) {
        if (decay_modes.empty()) {
          printf("Warning: No decay modes found for particle %i.\n",
                 pdgcode);
        } else {
          /* Check if ratios add to 1 */
          if (fabs(ratio_sum - 1.0) > really_small) {
            /* They didn't; renormalize */
            printf("Particle %i:\n", pdgcode);
            decay_modes.renormalize(ratio_sum);
          }
          /* Add the list of decay modes for this particle type */
          particles->add_decaymodes(decay_modes, pdgcode);
          /* Clean up the list for the next particle type */
          decay_modes.clear();
        }
      }
    }
  }
  free(line);
  fclose(file);
  printf("Finished reading decaymodes.txt\n");
  return;
}
