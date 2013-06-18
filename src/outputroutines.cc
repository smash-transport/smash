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
#include <map>

#include "include/Box.h"
#include "include/FourVector.h"
#include "include/macros.h"
#include "include/Parameters.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"

/* print_line - output a visible seperator */
static void print_line(void) {
  int field_width = 80;

  for (int i = 0; i < field_width; i++)
    printf("-");
  printf("\n");
}

/* print_startup - console output on startup */
void print_startup(const Box &box, const Parameters &parameters) {
  printf("Size of the box: %g x %g x %g [fm]\n", box.length(), box.length(),
    box.length());
  printf("Initial temperature: %g [GeV]\n", box.temperature());
  printf("Elastic cross section: %g [mb]\n", parameters.cross_section());
  printf("Using temporal stepsize: %g [GeV]\n", parameters.eps());
  printf("Maximum number of steps: %i \n", box.steps());
  printf("Random number seed: %li \n", parameters.seed());
}

/* print_header - title for each row */
void print_header(void) {
  print_line();
  printf(" Time    <Etot>    <Ediff>        <ptot>     <scattrate>"
         "  <scatt>    <timing>   <particles>\n");
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

/* measure_timediff - time the simulation used */
double measure_timediff(const Box &box) {
  timespec now;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
  return (now.tv_sec + now.tv_nsec / 10.0E9
    - box.time_start().tv_sec -   box.time_start().tv_nsec / 10.0E9);
}

/* print_measurements - console output during simulation */
void print_measurements(const std::map<int, ParticleData> &particles,
                        const size_t &scatterings_total,
                        const size_t &scatterings_this_interval,
                        const Box &box) {
  FourVector momentum_total(0, 0, 0, 0);
  /* calculate elapsed time */
  double elapsed = measure_timediff(box);
  double time = 0.0;
  int current_particle_total = 0;

  for (std::map<int, ParticleData>::const_iterator i = particles.begin();
       i != particles.end(); ++i) {
    momentum_total += i->second.momentum();
    /* use the time from the last active particle - startup time */
    time = i->second.position().x0() - 1.0;
    current_particle_total++;
  }
  if (likely(time > 0))
    printf("%5g%10g%13g%13g%13g%10zu%13g%10i\n", time, momentum_total.x0(),
           box.energy_initial() - momentum_total.x0(),
           sqrt(-1 * momentum_total.DotThree()),
           scatterings_total * 2 / (particles.size() * time),
           scatterings_this_interval, elapsed, current_particle_total);
  else
    printf("%5g%10g%13g%13g%13g%10i%13g%10i\n", time, momentum_total.x0(),
          box.energy_initial() - momentum_total.x0(),
           sqrt(-1 * momentum_total.DotThree()), 0.0, 0, elapsed,
           current_particle_total);
}

/* print_tail - output at the end of the simulation */
void print_tail(const Box &box, const double &scattering_rate) {
  double time = measure_timediff(box);
  print_line();
  /* print finishing time in human readable way:
   * time < 10 min => seconds
   * 10 min < time < 3 h => minutes
   * time > 3h => hours
   */
  if (time < 600)
    printf("Time real: %g [s]\n", time);
  else if (time < 10800)
    printf("Time real: %g [min]\n", time / 60);
  else
    printf("Time real: %g [h]\n", time / 3600);
  printf("Final scattering rate: %g [fm-1]\n", scattering_rate);
}

/* printd_momenta - print debug data of the specific particle */
void printd_momenta(const ParticleData &particle __attribute__((unused))) {
  printd("Particle %d momenta: %g %g %g %g [GeV]\n", particle.id(),
      particle.momentum().x0(), particle.momentum().x1(),
      particle.momentum().x2(), particle.momentum().x3());
}

/* printd_position - print debug data of the specific particle */
void printd_position(const ParticleData &particle __attribute__((unused))) {
  printd("Particle %d position: %g %g %g %g [fm]\n", particle.id(),
      particle.position().x0(), particle.position().x1(),
      particle.position().x2(), particle.position().x3());
}

/* write_particles - writes out data of the specific particles */
void write_particles(const std::map<int, ParticleData> &particles) {
  FILE *fp;
  char filename[256];

  snprintf(filename, sizeof(filename), "data/momenta_%.5f.dat",
           particles.begin()->second.position().x0() - 1.0);
  fp = fopen(filename, "w");
  for (std::map<int, ParticleData>::const_iterator i = particles.begin();
       i != particles.end(); ++i) {
     fprintf(fp, "%g %g %g %g\n", i->second.momentum().x0(),
             i->second.momentum().x1(), i->second.momentum().x2(),
             i->second.momentum().x3());
  }
  fclose(fp);
  snprintf(filename, sizeof(filename), "data/position_%.5f.dat",
           particles.begin()->second.position().x0() - 1.0);
  fp = fopen(filename, "w");
  for (std::map<int, ParticleData>::const_iterator i = particles.begin();
       i != particles.end(); ++i) {
     fprintf(fp, "%g %g %g %g\n", i->second.position().x0(),
             i->second.position().x1(), i->second.position().x2(),
             i->second.position().x3());
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
  const ParticleData &particle2, const ParticleType &type1,
  const ParticleType &type2, int flag) {
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
  momentum = particle1.momentum(), position = particle1.position();
  fprintf(fp, "%i %i %i %g %g %g %g %g %g %g %g %g \n",
              particle1.id(), 0, 0,
              momentum.x1(), momentum.x2(), momentum.x3(), momentum.x0(),
              type1.mass(), position.x1(), position.x2(), position.x3(),
              position.x0());
  momentum = particle2.momentum(), position = particle2.position();
  fprintf(fp, "%i %i %i %g %g %g %g %g %g %g %g %g \n",
              particle2.id(), 0, 0,
              momentum.x1(), momentum.x2(), momentum.x3(), momentum.x0(),
              type2.mass(), position.x1(), position.x2(), position.x3(),
              position.x0());
  fclose(fp);
}

/* write_vtk - VTK file describing particle position */
void write_vtk(const std::map<int, ParticleData> &particles) {
  FILE *fp;
  char filename[256];

  snprintf(filename, sizeof(filename), "data/pos_%.5f.vtk",
           particles.begin()->second.position().x0() - 1.0);
  fp = fopen(filename, "w");
  /* Legacy VTK file format */
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "Generated from molecular-offset data\n");
  fprintf(fp, "ASCII\n");
  /* Unstructured data sets are composed of points, lines, polygons, .. */
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp, "POINTS %lu double\n", particles.size());
  for (std::map<int, ParticleData>::const_iterator i = particles.begin();
       i != particles.end(); ++i)
    fprintf(fp, "%g %g %g\n", i->second.position().x1(),
            i->second.position().x2(), i->second.position().x3());
  fprintf(fp, "POINT_DATA %lu\n", particles.size());
  fprintf(fp, "SCALARS momenta_x double 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for (std::map<int, ParticleData>::const_iterator i = particles.begin();
       i != particles.end(); ++i)
    fprintf(fp, "%g\n", i->second.momentum().x1());
  fclose(fp);
}
