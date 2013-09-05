/*
 *
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 *
 */
#include "include/initial-conditions.h"

#include <stdio.h>
#include <stdint.h>
#include <gsl/gsl_sf_bessel.h>

#include "include/Box.h"
#include "include/constants.h"
#include "include/distributions.h"
#include "include/macros.h"
#include "include/Parameters.h"
#include "include/Particles.h"
#include "include/ParticleData.h"
#include "include/ParticleType.h"
#include "include/outputroutines.h"

/* initial_conditions - sets particle type */
void initial_particles(Particles *particles) {
  /* XXX: use nosql table for particle type values */
  ParticleType piplus("pi+", 0.13957f, -1.0, 211, 1, 1, 0);
  particles->add_type(piplus, 211);
  ParticleType piminus("pi-", 0.13957f, -1.0, -211, 1, -1, 0);
  particles->add_type(piminus, -211);
  ParticleType pi0("pi0", 0.134977f, -1.0, 111, 1, 0, 0);
  particles->add_type(pi0, 111);
}

/* initial_conditions - sets particle data for @particles */
void initial_conditions(Particles *particles, Parameters *parameters,
  Box *box) {
  double phi, cos_theta, sin_theta, momentum_radial, number_density_total = 0;
  FourVector momentum_total(0, 0, 0, 0);
  size_t number_total = 0, number = 0;

  /* initialize random seed */
  srand48(parameters->seed());

  if (particles->types().empty()) {
    fprintf(stderr, "E: No particle types\n");
    exit(EXIT_FAILURE);
  }

  /* Let's check how many non-resonances we have */
  unsigned int non_resonances = 0;
  printd("IC has %lu particle types\n", particles->types().size());
  for (std::map<int, ParticleType>::const_iterator
       i = particles->types().begin(); i != particles->types().end(); ++i) {
    if (i->second.width() < 0.0)
      non_resonances++;
  }

  /* loop over all the particle types */
  for (std::map<int, ParticleType>::const_iterator
       i = particles->types().begin(); i != particles->types().end(); ++i) {
    /* Particles with width > 0 (resonances) do not exist in the beginning */
    if (i->second.width() > 0.0)
      continue;

    /* Number of non-resonances left */
    non_resonances--;

    printd("%s mass: %g [GeV]\n", i->second.name().c_str(), i->second.mass());
    /*
     * The particle number depends on distribution function
     * (assumes Bose-Einstein):
     * Volume m^2 T BesselK[2, m/T] / (2\pi^2)
     */
    double number_density = i->second.mass()
      * i->second.mass() * box->temperature()
      * gsl_sf_bessel_Knu(2, i->second.mass() / box->temperature())
      * 0.5 * M_1_PI * M_1_PI / hbarc / hbarc / hbarc;

    /* particle number depending on IC geometry either sphere or box */
    if (unlikely(parameters->initial_condition() == 3)) {
      /* cast while reflecting probability of extra particle */
      number = 4.0 / 3.0 * M_PI * box->length() * box->length() * box->length()
        * number_density * parameters->testparticles();
      if (4.0 / 3.0 * M_PI * box->length() * box->length() * box->length()
        * number_density - number > drand48())
        number++;
    } else {
      /* cast while reflecting probability of extra particle */
      number = box->length() * box->length() * box->length() * number_density
        * parameters->testparticles();
      if (box->length() * box->length() * box->length() * number_density
        - number > drand48())
        number++;
    }
    printf("IC number density %.6g [fm^-3]\n", number_density);
    printf("IC %zu number of %s\n", number, i->second.name().c_str());
    number_density_total += number_density;

    /* Set random IC:
     * particles at random position in the box with thermal momentum
     */
    /* allocate the particles */
    ParticleData particle_new;
    for (size_t id = number_total; id < number_total + number; id++) {
      double x, y, z, time_start;
      /* ID uniqueness check */
      if (unlikely(particles->count(id) > 0))
        continue;

      /* back to back pair creation with random momenta direction */
      if (unlikely(id == number + number_total - 1 && !(id % 2)
          && non_resonances == 0)) {
        /* poor last guy just sits around */
        particle_new.set_momentum(i->second.mass(), 0, 0, 0);
      } else if (!(id % 2)) {
        if (parameters->initial_condition() != 2) {
          /* thermal momentum according Maxwell-Boltzmann distribution */
          momentum_radial = sample_momenta(*box, i->second);
        } else {
          /* IC == 2 initial thermal momentum is the average 3T */
          momentum_radial = 3.0 * box->temperature();
        }
        /* phi in the range from [0, 2 * pi) */
        phi = 2.0 * M_PI * drand48();
        /* cos(theta) in the range from [-1.0, 1.0) */
        cos_theta = -1.0 + 2.0 * drand48();
        sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        printd("Particle %zu radial momenta %g phi %g cos_theta %g\n", id,
          momentum_radial, phi, cos_theta);
        particle_new.set_momentum(i->second.mass(),
          momentum_radial * cos(phi) * sin_theta,
          momentum_radial * sin(phi) * sin_theta,
          momentum_radial * cos_theta);
      } else {
        particle_new.set_momentum(i->second.mass(),
          - particles->data(id - 1).momentum().x1(),
          - particles->data(id - 1).momentum().x2(),
          - particles->data(id - 1).momentum().x3());
      }
      momentum_total += particle_new.momentum();

      time_start = 1.0;
      /* ramdom position in a quadratic box */
      if (parameters->initial_condition() != 3) {
        x = drand48() * box->length();
        y = drand48() * box->length();
        z = drand48() * box->length();
      /* ramdom position in a sphere
       * box length here has the meaning of the sphere radius
       */
      } else {
        x = -box->length() + 2.0 * drand48() * box->length();
        y = -box->length() + 2.0 * drand48() * box->length();
        z = -box->length() + 2.0 * drand48() * box->length();
        /* sampling points inside of the sphere, rejected if outside */
        while (sqrt(x * x + y * y + z * z) > box->length()) {
          x = -box->length() + 2.0 * drand48() * box->length();
          y = -box->length() + 2.0 * drand48() * box->length();
          z = -box->length() + 2.0 * drand48() * box->length();
        }
      }
      particle_new.set_position(time_start, x, y, z);

      /* no collision yet hence zero time and set id */
      particle_new.set_collision(-1, 0, -1);
      particle_new.set_id(id);

      /* create new particle id is enhanced in the class itself */
      particles->add_data(particle_new);

      /* IC: debug checks */
      printd_momenta(particles->data(id));
      printd_position(particles->data(id));
    }
    number_total += number;
  }
  printf("IC total number density %.6g [fm^-3]\n", number_density_total);
  printf("IC contains %zu particles\n", number_total);
  /* loop over all particles */
  number = number_total;
  /* reducing cross section according to number of test particle */
  if (parameters->testparticles() > 1) {
    printf("IC test particle: %i\n", parameters->testparticles());
    parameters->set_cross_section(parameters->cross_section()
                                  / parameters->testparticles());
    printf("Elastic cross section: %g [mb]\n", parameters->cross_section());
  }

  /* Display on startup if pseudo grid is used */
  int const grid_number = round(box->length()
              / sqrt(parameters->cross_section() * fm2_mb * M_1_PI) * 0.5);
  /* pseudo grid not used for 3^3 or extremely small particle numbers */
  if (grid_number >= 4 && number > 10)
    printf("Simulation with pseudo grid: %d^3\n", grid_number);
  else
    printf("W: Not using pseudo grid: %d^3\n", grid_number);

  /* allows to check energy conservation */
  printf("IC total energy: %g [GeV]\n", momentum_total.x0());
  box->set_energy_initial(momentum_total.x0());
  box->set_number_density_inital(number_density_total);
}
