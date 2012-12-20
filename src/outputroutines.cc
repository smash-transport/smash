/*
 *
 *    Copyright (c) 2012
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/outputroutines.h"
#include "include/box.h"

#include <stdio.h>

void print_startup(void) {
  printf("Size of the box: %g [fm]\n", A);
  printf("Initial temperature: %g [GeV]\n", temperature);
  printf("Using temporal stepsize: %g [GeV]\n", EPS);
}
