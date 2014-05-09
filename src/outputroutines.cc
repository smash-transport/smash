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
#include <cmath>
#include <cstdio>
#include <ctime>
#include <list>
#include <map>
#include <string>
#include <utility>

#include "include/chrono.h"
#include "include/fourvector.h"
#include "include/particles.h"
#include "include/particledata.h"
#include "include/particletype.h"
#include "include/macros.h"

namespace Smash {

/* print_line - output a visible seperator */
static void print_line(void) {
  int field_width = 80;

  for (int i = 0; i < field_width; i++)
    printf("-");
  printf("\n");
}

/* print_header - title for each row */
void print_header(void) {
  print_line();
  printf(" Time    <Ediff>       <ptot>    <scattrate>"
         "  <scatt>   <particles>  <timing>\n");
  print_line();
}


/* print_measurements - console output during simulation */
void print_measurements(const Particles &particles,
                        const size_t &scatterings_total,
                        const size_t &scatterings_this_interval,
                        float energy_ini,
                SystemTimePoint time_start) {
  FourVector momentum_total(0, 0, 0, 0);
  /* calculate elapsed time */
  SystemTimeSpan elapsed_seconds = SystemClock::now() - time_start;

  for (const ParticleData &data : particles.data()) {
    momentum_total += data.momentum();
  }
  double time = particles.time();

  if (likely(time > 0))
    printf("%5g%13g%13g%13g%10zu%10zu%13g\n", time,
           energy_ini - momentum_total.x0(),
           sqrt(-1 * momentum_total.DotThree()),
           scatterings_total * 2 / (particles.size() * time),
           scatterings_this_interval, particles.size(), elapsed_seconds.count());
  else
    printf("%5g%13g%13g%13g%10i%10zu%13g\n", time,
           energy_ini - momentum_total.x0(),
           sqrt(-1 * momentum_total.DotThree()), 0.0, 0, particles.size(),
           elapsed_seconds.count());
}

/* print_tail - output at the end of the simulation */
void print_tail(const
                SystemTimePoint time_start,
                const double &scattering_rate) {
  SystemTimeSpan time = SystemClock::now() - time_start;
  print_line();
  /* print finishing time in human readable way:
   * time < 10 min => seconds
   * 10 min < time < 3 h => minutes
   * time > 3h => hours
   */
  if (time.count() < 600)
    printf("Time real: %g [s]\n", time.count());
  else if (time.count() < 10800)
    printf("Time real: %g [min]\n", time.count() / 60);
  else
    printf("Time real: %g [h]\n", time.count() / 3600);
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
  for (const ParticleData &data : particles->data()) {
    fprintf(fp, "%i %s %i %g %g %g %g %g %g %g %g %g \n",
            data.id(), data.pdgcode().string().c_str(), 0,
            data.momentum().x1(), data.momentum().x2(),
            data.momentum().x3(), data.momentum().x0(),
            sqrt(data.momentum().Dot(data.momentum())),
            data.position().x1(), data.position().x2(),
            data.position().x3(), data.position().x0());
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
  float mass = sqrt(momentum.Dot(momentum));
  fprintf(fp, "%i %s %i %g %g %g %g %g %g %g %g %g \n", particle_data.id(),
          particle_type.pdgcode().string().c_str(), 0, momentum.x1(),
          momentum.x2(), momentum.x3(), momentum.x0(), mass,
          position.x1(), position.x2(), position.x3(),
          position.x0());

  fclose(fp);
}

}  // namespace Smash
