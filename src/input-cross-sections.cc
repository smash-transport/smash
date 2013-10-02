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

#include "include/CrossSections.h"
#include "include/constants.h"
#include "include/input-cross-sections.h"
#include "include/outputroutines.h"

/* XXX: hardcoded length cap */
#define FILELEN 256

/* input_cross_sections - read cross section parametrizations */
void input_cross_sections(CrossSections *cross_sections, char *path) {
  char input_cross_sections[FILELEN];
  FILE *file;

  /* Looking for cross section parameter list
   * If none exists, we'll use default values.
   */
  snprintf(input_cross_sections, strlen(path) + 16,
           "%s/cross-sections.txt", path);
  file = fopen(input_cross_sections, "r");
  if (!file) {
    fprintf(stderr, "W: No cross-sections.txt at %s path.\n", path);
    return;
  }

  printf("Processing %s/cross-sections.txt.\n", path);

  char *line = NULL, *line_position = NULL;
  char *characters;
  size_t line_size = 0;
  ssize_t characters_read;
  std::vector<float> parameters_values;
  bool pp_elastic = true, pp_total = true, pn_elastic = true,
    pn_total = true, ppbar_elastic = true, ppbar_annihilation = true,
    ppbar_total = true;
  while ((characters_read = getline(&line, &line_size, file)) != -1) {
    printf("Retrieved cross-sections.txt line of length %li:\n",
           characters_read);
    /* Skip comments and blank lines */
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\t'
        || line[0] == '/') {
      printd("Skipping line: %s", line);
      continue;
    }
    printd("line: %s", line);
    /* First item in a line should be the minimum p_lab */
    characters = strtok_r(line, " ", &line_position);
    if (characters == NULL)
      continue;
    float p_lab_min = atof(characters);
    parameters_values.push_back(p_lab_min);
    printd("p_lab: %d\n", p_lab_min);
    /* Rest of the line should contain the parameter values */
    characters = strtok_r(NULL, " ", &line_position);
    printd("characters: %s \n", characters);
    while (characters != NULL) {
      float parameter_value = atof(characters);
      printd("parameter value: %d \n", parameter_value);
      parameters_values.push_back(parameter_value);
      characters = strtok_r(NULL, " ", &line_position);
      printd("characters: %s \n", characters);
    }
    /* Add the parameter values to cross section */
    if (pp_elastic) {
      cross_sections->add_pp_elastic(parameters_values);
      if (p_lab_min < really_small)
        pp_elastic = false;
    } else if (pp_total) {
      cross_sections->add_pp_total(parameters_values);
      if (p_lab_min < really_small)
        pp_total = false;
    } else if (pn_elastic) {
      cross_sections->add_pn_elastic(parameters_values);
      if (p_lab_min < really_small)
        pn_elastic = false;
    } else if (pn_total) {
      cross_sections->add_pn_total(parameters_values);
      if (p_lab_min < really_small)
        pn_total = false;
    } else if (ppbar_elastic) {
      cross_sections->add_ppbar_elastic(parameters_values);
      if (p_lab_min < really_small)
        ppbar_elastic = false;
    } else if (ppbar_annihilation) {
      cross_sections->add_ppbar_annihilation(parameters_values);
      if (p_lab_min < really_small)
        ppbar_annihilation = false;
    } else if (ppbar_total) {
      cross_sections->add_ppbar_total(parameters_values);
      if (p_lab_min < really_small)
        ppbar_total = false;
    } else {
      printf("Warning: No place for additional parameters!\n");
    }
    /* Clean the parameter values */
    parameters_values.clear();
  }
  free(line);
  fclose(file);
  printf("Finished reading cross-sections.txt\n");
  return;
}
