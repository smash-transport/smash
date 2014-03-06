/*
 *    Copyright (c) 2013-14
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <utility>
#include <vector>

#include "include/algorithms.h"
#include "include/angles.h"
#include "include/boxmodus.h"
#include "include/collisions.h"
#include "include/configuration.h"
#include "include/constants.h"
#include "include/crosssections.h"
#include "include/distributions.h"
#include "include/experimentparameters.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/parameters.h"
#include "include/particles.h"
#include "include/random.h"

namespace Smash {

BoxModus::BoxModus(Configuration modus_config, const ExperimentParameters &)
    : initial_condition_(modus_config.take({"Box", "INITIAL_CONDITION"})),
      length_           (modus_config.take({"Box", "LENGTH"})),
      temperature_      (modus_config.take({"Box", "TEMPERATURE"})) {
}

/* print_startup - console output on startup of box specific parameters */
void BoxModus::print_startup() {
  printf("Size of the box: %g x %g x %g fm\n", length_, length_, length_);
  printf("Initial temperature: %g GeV\n", temperature_);
  printf("IC type %d\n", initial_condition_);
}

/* initial_conditions - sets particle data for @particles */
void BoxModus::initial_conditions(Particles *particles,
                                  const ExperimentParameters &parameters) {
  double momentum_radial, number_density_total = 0;
  Angles phitheta;
  FourVector momentum_total(0, 0, 0, 0);
  size_t number_total = 0;
  /* loop over all the particle types */
  for (auto i = particles->types_cbegin(); i != particles->types_cend(); ++i) {
    /* Particles with width > 0 (resonances) do not exist in the beginning */
    if (i->second.width() > 0.0) {
      continue;
    }
    printd("%s mass: %g [GeV]\n", i->second.name().c_str(), i->second.mass());
    /* bose einstein distribution function */
    double number_density =
        number_density_bose(i->second.mass(), this->temperature_);
    // calculate the expected of particles in the box
    double real_number = this->length_ * this->length_ * this->length_ *
                         number_density * parameters.testparticles;
    size_t int_number = static_cast<size_t>(real_number);
    // decide if we have an extra particle: the probability for that is
    // equal to the fractional part of the number.
    if (real_number - int_number > rng.canonical())
      int_number++;
    printf("IC number density %.6g [fm^-3]\n", number_density);
    printf("IC %zu number of %s\n", int_number, i->second.name().c_str());
    number_density_total += number_density;
    /* create bunch of particles */
    printf("IC creating %zu particles\n", int_number);
    particles->create(int_number, i->second.pdgcode());
    number_total += int_number;
  }
  printf("IC total number density %.6g [fm^-3]\n", number_density_total);
  printf("IC contains %zu particles\n", number_total);
  rng.set_uniform(0.0, this->length_);
  /* Set paricles IC: */
  for (auto i = particles->begin(); i != particles->end(); ++i) {
    double x, y, z, time_begin;
    /* back to back pair creation with random momenta direction */
    if (unlikely(i->first == particles->id_max() && !(i->first % 2))) {
      /* poor last guy just sits around */
      i->second.set_momentum(particles->type(i->first).mass(), 0, 0, 0);
    } else if (!(i->first % 2)) {
      if (this->initial_condition_ != 2) {
        /* thermal momentum according Maxwell-Boltzmann distribution */
        momentum_radial = sample_momenta(this->temperature_,
                                         particles->type(i->first).mass());
      } else {
        /* IC == 2 initial thermal momentum is the average 3T */
        momentum_radial = 3.0 * this->temperature_;
      }
      phitheta.distribute_isotropically();
      printd("Particle %d radial momenta %g phi %g cos_theta %g\n", i->first,
             momentum_radial, phitheta.phi(), phitheta.costheta());
      i->second.set_momentum(
          particles->type(i->first).mass(), momentum_radial * phitheta.x(),
          momentum_radial * phitheta.y(), momentum_radial * phitheta.z());
    } else {
      i->second.set_momentum(particles->type(i->first).mass(),
                             -particles->data(i->first - 1).momentum().x1(),
                             -particles->data(i->first - 1).momentum().x2(),
                             -particles->data(i->first - 1).momentum().x3());
    }
    momentum_total += i->second.momentum();
    time_begin = 1.0;
    /* ramdom position in a quadratic box */
    x = rng.uniform();
    y = rng.uniform();
    z = rng.uniform();
    i->second.set_position(time_begin, x, y, z);
    /* IC: debug checks */
    printd_momenta(i->second);
    printd_position(i->second);
  }
  /* Display on startup if pseudo grid is used */
  size_t number = number_total;
  int const grid_number = round(
      this->length_ / sqrt(parameters.cross_section * fm2_mb * M_1_PI) * 0.5);
  /* pseudo grid not used for 3^3 or extremely small particle numbers */
  if (grid_number >= 4 && number > 10)
    printf("Simulation with pseudo grid: %d^3\n", grid_number);
  else
    printf("W: Not using pseudo grid: %d^3\n", grid_number);
  /* allows to check energy conservation */
  printf("IC total energy: %g [GeV]\n", momentum_total.x0());
  number_density_initial_ = number_density_total;
}

/* evolve - the core of the box, stepping forward in time */
int BoxModus::sanity_check(Particles *particles) {
  /* fixup positions on startup, particles need to be *inside* the box */
  for (auto i = particles->begin(); i != particles->end(); ++i) {
    FourVector p = i->second.position();
    enforce_periodic_boundaries(p.begin() + 1, p.end(), length_);
    i->second.set_position(p);
  }
  return 0;
}

/* check_collision_geometry - check if a collision happens between particles */
void BoxModus::check_collision_geometry(
    Particles *particles, CrossSections *cross_sections,
    std::list<int> *collision_list, size_t *rejection_conflict,
    const ExperimentParameters &parameters) {
  std::vector<std::vector<std::vector<std::vector<int> > > > grid;
  int N;
  int x, y, z;

  /* For small boxes no point in splitting up in grids */
  /* calculate approximate grid size according to double interaction length */
  N = round(length_ / sqrt(parameters.cross_section * fm2_mb * M_1_PI) * 0.5);
  if (unlikely(N < 4 || particles->size() < 10)) {
    /* XXX: apply periodic boundary condition */
    ModusDefault::check_collision_geometry(particles, cross_sections,
                                           collision_list, rejection_conflict,
                                           parameters);
    return;
  }
  /* allocate grid */
  grid.resize(N);
  for (int i = 0; i < N; i++) {
    grid[i].resize(N);
    for (int j = 0; j < N; j++)
      grid[i][j].resize(N);
  }
  /* populate grid */
  for (auto i = particles->begin(); i != particles->end(); ++i) {
    /* XXX: function - map particle position to grid number */
    x = round(i->second.position().x1() / length_ * (N - 1));
    y = round(i->second.position().x2() / length_ * (N - 1));
    z = round(i->second.position().x3() / length_ * (N - 1));
    printd_position(i->second);
    printd("grid cell particle %i: %i %i %i of %i\n", i->first, x, y, z, N);
    if (unlikely(x >= N || y >= N || z >= N))
      printf("W: Particle outside the box: %g %g %g \n",
             i->second.position().x1(), i->second.position().x2(),
             i->second.position().x3());
    grid[x][y][z].push_back(i->first);
  }
  /* semi optimised nearest neighbour search:
   * http://en.wikipedia.org/wiki/Cell_lists
   */
  FourVector shift;
  for (auto i = particles->begin(); i != particles->end(); ++i) {
    /* XXX: function - map particle position to grid number */
    x = round(i->second.position().x1() / length_ * (N - 1));
    y = round(i->second.position().x2() / length_ * (N - 1));
    z = round(i->second.position().x3() / length_ * (N - 1));
    if (unlikely(x >= N || y >= N || z >= N))
      printf("grid cell particle %i: %i %i %i of %i\n", i->first, x, y, z, N);
    /* check all neighbour grids */
    for (int cx = -1; cx < 2; cx++) {
      int sx = cx + x;
      /* apply periodic boundary condition for particle positions */
      if (sx < 0) {
        sx = N - 1;
        shift.set_x1(-length_);
      } else if (sx > N - 1) {
        sx = 0;
        shift.set_x1(length_);
      } else {
        shift.set_x1(0);
      }
      for (int cy = -1; cy < 2; cy++) {
        int sy = cy + y;
        if (sy < 0) {
          sy = N - 1;
          shift.set_x2(-length_);
        } else if (sy > N - 1) {
          sy = 0;
          shift.set_x2(length_);
        } else {
          shift.set_x2(0);
        }
        for (int cz = -1; cz < 2; cz++) {
          int sz = cz + z;
          if (sz < 0) {
            sz = N - 1;
            shift.set_x3(-length_);
          } else if (sz > N - 1) {
            sz = 0;
            shift.set_x3(length_);
          } else {
            shift.set_x3(0);
          }
          /* empty grid cell */
          if (grid[sx][sy][sz].empty()) {
            continue;
          }
          /* grid cell particle list */
          for (auto id_other = grid[sx][sy][sz].begin();
               id_other != grid[sx][sy][sz].end(); ++id_other) {
            /* only check against particles above current id
             * to avoid double counting
             */
            if (*id_other <= i->first) {
              continue;
            }
            printd("grid cell particle %i <-> %i\n", i->first, *id_other);
            if (shift == 0) {
              collision_criteria_geometry(
                  particles, cross_sections, collision_list, parameters.eps,
                  i->first, *id_other, rejection_conflict);
            } else {
              /* apply eventual boundary before and restore after */
              particles->data_pointer(*id_other)
                  ->set_position(particles->data(*id_other).position() + shift);
              collision_criteria_geometry(
                  particles, cross_sections, collision_list, parameters.eps,
                  i->first, *id_other, rejection_conflict);
              particles->data_pointer(*id_other)
                  ->set_position(particles->data(*id_other).position() - shift);
            }
          } /* grid particles loop */
        }   /* grid sz */
      }     /* grid sy */
    }       /* grid sx */
  }         /* outer particle loop */
}

/* propagate all particles */
void BoxModus::propagate(Particles *particles,
                         const ExperimentParameters &parameters) {
  FourVector distance, position;
  for (auto i = particles->begin(); i != particles->end(); ++i) {
    /* propagation for this time step */
    distance.set_FourVector(parameters.eps,
                            i->second.velocity_x() * parameters.eps,
                            i->second.velocity_y() * parameters.eps,
                            i->second.velocity_z() * parameters.eps);
    printd("Particle %d motion: %g %g %g %g\n", i->first, distance.x0(),
           distance.x1(), distance.x2(), distance.x3());
    /* treat the box boundaries */
    position = i->second.position();
    position += distance;
    bool wall_hit = enforce_periodic_boundaries(position.begin() + 1,
                                                position.end(), length_);
    if (wall_hit) {
      write_oscar(i->second, particles->type(i->first), 1, 1);
    }
    i->second.set_position(position);
    if (wall_hit) {
      write_oscar(i->second, particles->type(i->first));
    }
    printd_position(i->second);
  }
}

}  // namespace Smash
