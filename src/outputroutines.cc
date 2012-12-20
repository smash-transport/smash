/*
 *
 *    Copyright (c) 2012
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/outputroutines.h"

#include <stdio.h>

#include "include/box.h"

void print_startup(void) {
  printf("Size of the box: %g x %g x %g [fm]\n", A, A, A);
  printf("Initial temperature: %g [GeV]\n", temperature);
  printf("Using temporal stepsize: %g [GeV]\n", EPS);
  printf("Maximum number of steps: %i \n", STEPS);
  printf("Random number seed: %i \n", seedp);
}
