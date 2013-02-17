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
#include <sys/stat.h>
#include <sys/types.h>

#include "include/box.h"
#include "include/ParticleData.h"

/* print_startup - console output on startup */
void print_startup(void) {
  printf("Size of the box: %g x %g x %g [fm]\n", A, A, A);
  printf("Initial temperature: %g [GeV]\n", temperature);
  printf("Using temporal stepsize: %g [GeV]\n", EPS);
  printf("Maximum number of steps: %i \n", STEPS);
  printf("Random number seed: %u \n", seedp);
}

/* mkdir_data - directory for data files */
void mkdir_data(void) {
  int ret;

  ret = mkdir("data", 0751);
  if (ret == 0) {
    printf("dir 'data' successfully created.\n");
    return;
  }
  fprintf(stderr, "mkdir 'data' failed.\n");
}

/* printd_momenta - print debug data of the specific particle */
void printd_momenta(ParticleData particle) {
  printd("Particle %d momenta: %g %g %g %g [GeV]\n", particle.id(),
      particle.momentum().x0(), particle.momentum().x1(),
      particle.momentum().x2(), particle.momentum().x3());
}

/* printd_position - print debug data of the specific particle */
void printd_position(ParticleData particle) {
  printd("Particle %d position: %g %g %g %g [fm]\n", particle.id(),
      particle.x().x0(), particle.x().x1(), particle.x().x2(),
      particle.x().x3());
}

/* write_particles - writes out data of the specific particles */
void write_particles(ParticleData *particles, const int number) {
  FILE *fp;

  fp = fopen("data/momenta.dat", "a");
  for (int i = 0; i < number; i++) {
     fprintf(fp, "%g %g %g %g\n", particles[i].momentum().x0(),
       particles[i].momentum().x1(),  particles[i].momentum().x2(),
       particles[i].momentum().x3());
  }
  fclose(fp);
  fp = fopen("data/position.dat", "a");
  for (int i = 0; i < number; i++) {
     fprintf(fp, "%g %g %g %g\n", particles[i].x().x0(),
       particles[i].x().x1(),  particles[i].x().x2(),
       particles[i].x().x3());
  }
  fclose(fp);
}
