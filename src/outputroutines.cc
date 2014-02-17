/*
 *
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/outputroutines.h"

#include <sys/stat.h>
#include <cstdio>
#include <ctime>
#include <list>
#include <map>
#include <string>
#include <utility>

#include "include/fourvector.h"
#include "include/parameters.h"
#include "include/particles.h"
#include "include/particledata.h"
#include "include/particletype.h"
#include "include/macros.h"
#include "include/time.h"

/* print_line - output a visible seperator */
static void print_line(void) {
  int field_width = 80;

  for (int i = 0; i < field_width; i++)
    printf("-");
  printf("\n");
}

/* warn_wrong_params - warn about unassigned parameters */
void warn_wrong_params(std::list<Parameters> *configuration) {
  if (configuration->empty())
    return;

  for (auto i = configuration->begin(); i != configuration->end(); ++i) {
    char *key = i->key();
    char *value = i->value();
    printf("W: Unused parameter %s with %s\n", key, value);
  }
  /* Remove as we warned about them and can't do something */
  configuration->clear();
}


/* print_header - title for each row */
void print_header(void) {
  print_line();
  printf(" Time    <Ediff>       <ptot>    <scattrate>"
         "  <scatt>   <particles>  <timing>\n");
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
double measure_timediff(const timespec time_start) {
  timespec now;
  clock_gettime(&now);
  return (now.tv_sec + now.tv_nsec / 10.0E9
    - time_start.tv_sec -   time_start.tv_nsec / 10.0E9);
}

/* print_measurements - console output during simulation */
void print_measurements(const Particles &particles,
                        const size_t &scatterings_total,
                        const size_t &scatterings_this_interval,
                        float energy_ini, timespec time_start) {
  FourVector momentum_total(0, 0, 0, 0);
  /* calculate elapsed time */
  double elapsed = measure_timediff(time_start);
  double time = 0.0;

  for (auto i = particles.cbegin(); i != particles.cend(); ++i) {
    momentum_total += i->second.momentum();
    /* use the time from the last active particle - startup time */
    time = i->second.position().x0() - 1.0;
  }

  if (likely(time > 0))
    printf("%5g%13g%13g%13g%10zu%10zu%13g\n", time,
           energy_ini - momentum_total.x0(),
           sqrt(-1 * momentum_total.DotThree()),
           scatterings_total * 2 / (particles.size() * time),
           scatterings_this_interval, particles.size(), elapsed);
  else
    printf("%5g%13g%13g%13g%10i%10zu%13g\n", time,
           energy_ini - momentum_total.x0(),
           sqrt(-1 * momentum_total.DotThree()), 0.0, 0, particles.size(),
           elapsed);
}

/* print_tail - output at the end of the simulation */
void print_tail(const timespec time_start, const double &scattering_rate) {
  double time = measure_timediff(time_start);
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

/* printd_momenta - print debug data of the specific particle with message */
void printd_momenta(const char *message __attribute__((unused)),
  const ParticleData &particle __attribute__((unused))) {
  printd("%s: %g %g %g %g [GeV]\n", message,
      particle.momentum().x0(), particle.momentum().x1(),
      particle.momentum().x2(), particle.momentum().x3());
}

/* printd_momenta - print debug data of the specific particle */
void printd_momenta(const ParticleData &particle __attribute__((unused))) {
  printd("Particle %d momenta: %g %g %g %g [GeV]\n", particle.id(),
      particle.momentum().x0(), particle.momentum().x1(),
      particle.momentum().x2(), particle.momentum().x3());
}

/* printd_position - print debug data of the specific particle with message */
void printd_position(const char *message __attribute__((unused)),
  const ParticleData &particle __attribute__((unused))) {
  printd("%s: %g %g %g %g [fm]\n", message,
      particle.position().x0(), particle.position().x1(),
      particle.position().x2(), particle.position().x3());
}

/* printd_position - print debug data of the specific particle */
void printd_position(const ParticleData &particle __attribute__((unused))) {
  printd("Particle %d position: %g %g %g %g [fm]\n", particle.id(),
      particle.position().x0(), particle.position().x1(),
      particle.position().x2(), particle.position().x3());
}

void printd_list(const std::list<int> &collision_list) {
  printd("Collision list contains:");
  for (std::list<int>::const_iterator id = collision_list.cbegin();
       id != collision_list.cend(); ++id)
    printd(" particle %d", *id);
  printd("\n");
}

/**
 *  write_particles - writes out data of the specific particles
 */
void write_particles(const Particles &particles) {
  FILE *fp;
  char filename[256];

  snprintf(filename, sizeof(filename), "data/momenta_%.5f.dat",
           particles.time());
  fp = fopen(filename, "w");
  for (auto i = particles.cbegin(); i != particles.cend(); ++i) {
    fprintf(fp, "%g %g %g %g %i %i\n", i->second.momentum().x0(),
            i->second.momentum().x1(), i->second.momentum().x2(),
            i->second.momentum().x3(), i->second.id(), i->second.pdgcode());
  }
  fclose(fp);
  snprintf(filename, sizeof(filename), "data/position_%.5f.dat",
           particles.time());
  fp = fopen(filename, "w");
  for (auto i = particles.cbegin(); i != particles.cend(); ++i) {
    fprintf(fp, "%g %g %g %g %i %i\n", i->second.position().x0(),
            i->second.position().x1(), i->second.position().x2(),
            i->second.position().x3(), i->second.id(), i->second.pdgcode());
  }
  fclose(fp);
}

/* write_oscar_header - OSCAR header format */
void write_oscar_header(void) {
  FILE *fp;

  fp = fopen("data/collision.dat", "w");
  fprintf(fp, "# OSC1999A\n");
  fprintf(fp, "# Interaction history\n");
  fprintf(fp, "# smash \n");
  fprintf(fp, "# \n");
  fclose(fp);
}

/**
 *  write_oscar_event_block
 *  - writes the initial and final particle information of an event
 */
void write_oscar_event_block(Particles *particles,
                             size_t initial, size_t final, int event_id) {
  FILE *fp;
  fp = fopen("data/collision.dat", "a");
  /* OSCAR line prefix : initial particles; final particles; event id
   * First block of an event: initial = 0, final = number of particles
   * Vice versa for the last block
   */
  fprintf(fp, "%zu %zu %i\n", initial, final, event_id);
  for (auto i = particles->cbegin(); i != particles->cend(); ++i) {
    fprintf(fp, "%i %i %i %g %g %g %g %g %g %g %g %g \n",
            i->first, i->second.pdgcode(), 0,
            i->second.momentum().x1(), i->second.momentum().x2(),
            i->second.momentum().x3(), i->second.momentum().x0(),
            particles->type(i->first).mass(),
            i->second.position().x1(), i->second.position().x2(),
            i->second.position().x3(), i->second.position().x0() - 1.0);
  }
  fclose(fp);
}

/* write_oscar - OSCAR file */
/* Use this for the first particle in a process */
void write_oscar(const ParticleData &particle_data,
                 const ParticleType &particle_type,
                 const int initial, const int final) {
  FILE *fp;
  fp = fopen("data/collision.dat", "a");
  /* OSCAR line prefix : initial final
   * particle creation: 0 1
   * particle 2<->2 collision: 2 2
   * resonance formation: 2 1
   * resonance decay: 1 2
   * etc.
   */
  if (initial > 0 || final > 0)
    fprintf(fp, "%i %i \n", initial, final);

  /* particle_index, particle_pdgcode, ?, momenta, mass position */
  FourVector momentum = particle_data.momentum(),
    position = particle_data.position();
  fprintf(fp, "%i %i %i %g %g %g %g %g %g %g %g %g \n",
          particle_data.id(), particle_type.pdgcode(), 0,
          momentum.x1(), momentum.x2(), momentum.x3(), momentum.x0(),
          particle_type.mass(), position.x1(), position.x2(), position.x3(),
          position.x0() - 1.0);

  fclose(fp);
}

/* write_oscar - use this for the other particles in same process */
void write_oscar(const ParticleData &particle_data,
                 const ParticleType &particle_type) {
  FILE *fp;
  fp = fopen("data/collision.dat", "a");

  /* particle_index, particle_pdgcode, ?, momenta, mass position */
  FourVector momentum = particle_data.momentum(),
             position = particle_data.position();
  fprintf(fp, "%i %i %i %g %g %g %g %g %g %g %g %g \n", particle_data.id(),
          particle_type.pdgcode(), 0, momentum.x1(), momentum.x2(),
          momentum.x3(), momentum.x0(), particle_type.mass(), position.x1(),
          position.x2(), position.x3(), position.x0() - 1.0);

  fclose(fp);
}

/* write_vtk - VTK file describing particle position */
void write_vtk(const Particles &particles) {
  FILE *fp = NULL;
  char filename[256];
  snprintf(filename, sizeof(filename), "data/pos_0.%05i.vtk",
           static_cast<int>((particles.cbegin()->second.position().x0() - 1.0) *
                            10));
  fp = fopen(filename, "w");
  /* Legacy VTK file format */
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "Generated from molecular-offset data\n");
  fprintf(fp, "ASCII\n");
  /* Unstructured data sets are composed of points, lines, polygons, .. */
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp, "POINTS %zu double\n", particles.size());
  for (auto i = particles.cbegin(); i != particles.cend(); ++i)
    fprintf(fp, "%g %g %g\n", i->second.position().x1(),
            i->second.position().x2(), i->second.position().x3());
  fprintf(fp, "CELLS %zu %zu\n", particles.size(), particles.size() * 2);
  for (size_t point_index = 0; point_index < particles.size(); point_index++)
    fprintf(fp, "1 %zu\n", point_index);
  fprintf(fp, "CELL_TYPES %zu\n", particles.size());
  for (size_t point_index = 0; point_index < particles.size(); point_index++)
    fprintf(fp, "1\n");
  fprintf(fp, "POINT_DATA %zu\n", particles.size());
  fprintf(fp, "SCALARS pdg_codes int 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for (auto i = particles.cbegin(); i != particles.cend(); ++i)
    fprintf(fp, "%i\n", i->second.pdgcode());
  fprintf(fp, "VECTORS momentum double\n");
  for (auto i = particles.cbegin(); i != particles.cend(); ++i)
    fprintf(fp, "%g %g %g\n", i->second.momentum().x1(),
            i->second.momentum().x2(), i->second.momentum().x3());
  fclose(fp);
}
