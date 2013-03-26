/*
 *
 *    Copyright (c) 2012-2013
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
#include "include/FourVector.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"

/* print_line - output a visible seperator */
static void print_line(void) {
  int field_width = 13;

  for (int i = 0; i < 5 * field_width; i++)
    printf("-");
  printf("\n");
}

/* print_startup - console output on startup */
void print_startup(const box &box) {
  float A = box.a();
  printf("Size of the box: %g x %g x %g [fm]\n", A, A, A);
  printf("Initial temperature: %g [GeV]\n", box.temperature());
  printf("Elastic cross section: %g [mb]\n", box.cross_section());
  printf("Using temporal stepsize: %g [GeV]\n", box.eps());
  printf("Maximum number of steps: %i \n", box.steps());
  printf("Random number seed: %u \n", seedp);
}

/* print_startup - console output on startup */
void print_header(void) {
  print_line();
  printf("        Time        <Etot>      <ptot>  <scattering_rate>\n");
  print_line();
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

/* print_measurements - console output during simulation */
void print_measurements(const ParticleData *particles, const int &number,
  const size_t &scatterings_total) {
  FourVector momentum_total(0, 0, 0, 0);
  /* use the time from the first particle - startup time */
  double time = particles[0].x().x0() - 1.0;

  for (int i = 0; i < number; i++)
    momentum_total += particles[i].momentum();
  if (likely(time > 0))
    printf("%13g%13g%13g%13g\n", time, momentum_total.x0(),
          -1 * momentum_total.DotThree(),
          scatterings_total * 2 / (number * time));
  else
    printf("%13g%13g%13g%13g\n", time, momentum_total.x0(), 
          -1 * momentum_total.DotThree(), 0.0);
}

/* printd_momenta - print debug data of the specific particle */
void printd_momenta(const ParticleData &particle) {
  printd("Particle %d momenta: %g %g %g %g [GeV]\n", particle.id(),
      particle.momentum().x0(), particle.momentum().x1(),
      particle.momentum().x2(), particle.momentum().x3());
}

/* printd_position - print debug data of the specific particle */
void printd_position(const ParticleData &particle) {
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

/* write_oscar_header - OSCAR header format */
void write_oscar_header(void) {
  FILE *fp;

  fp = fopen("data/collision.dat", "w");
  fprintf(fp, "# OSC1999A\n");
  fprintf(fp, "# Elastic scattering history\n");
  fprintf(fp, "# Pionbox \n");
  fprintf(fp, "# \n");
  fclose(fp);
}

/* write_oscar - OSCAR file */
void write_oscar(const ParticleData &particle1,
  const ParticleData &particle2, int flag) {
  /* XXX: generalise particle types */
  ParticleType pi("pi", 0.13957);
  FourVector momentum, position;
  FILE *fp;

  fp = fopen("data/collision.dat", "a");
  /* flag adds OSCAR line prefix with the meanings:
   * particle creation -> 0 0
   * particle 2<-> collision -> 2 2
   */
  if (flag == 0)
    fprintf(fp, "%i %i \n", 0, 0);
  else if (flag == 1)
    fprintf(fp, "%i %i \n", 2, 2);
  /* particle_index, particle_pdgcode, ?, momenta, mass position */
  momentum = particle1.momentum(), position = particle1.x();
  fprintf(fp, "%i %i %i %g %g %g %g %g %g %g %g %g \n",
              particle1.id(), 0, 0,
              momentum.x1(), momentum.x2(), momentum.x3(), momentum.x0(),
              pi.mass(), position.x1(), position.x2(), position.x3(),
              position.x0());
  momentum = particle2.momentum(), position = particle2.x();
  fprintf(fp, "%i %i %i %g %g %g %g %g %g %g %g %g \n",
              particle2.id(), 0, 0,
              momentum.x1(), momentum.x2(), momentum.x3(), momentum.x0(),
              pi.mass(), position.x1(), position.x2(), position.x3(),
              position.x0());
  fclose(fp);
}
