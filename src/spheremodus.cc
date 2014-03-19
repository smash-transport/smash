/*
 *
 *    Copyright (c) 2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <cmath>

#include "include/box.h"
#include "include/collisions.h"
#include "include/constants.h"
#include "include/crosssections.h"
#include "include/decays.h"
#include "include/fourvector.h"
#include "include/initial-conditions.h"
#include "include/input-decaymodes.h"
#include "include/input-particles.h"
#include "include/laboratory.h"
#include "include/macros.h"
#include "include/outputroutines.h"
#include "include/particledata.h"
#include "include/particles.h"
#include "include/propagation.h"
#include "include/sphere.h"

/* build dependent variables */
#include "include/config.h"

namespace Smash {

SphereModus::SphereModus(Configuration modus_config)
    : radius_(modus_config.take({"Sphere", "RADIUS"})),
      timer_start_(set_timer_start()) {
}

/* print_startup - console output on startup of sphere specific parameters */
/* void print_startup(const SphereModus &ball) {
   /* printf("Volume of the sphere: 4 * pi * %g^2 [fm]\n", ball.radius);
   /* }


/* initial_conditions - sets particle data for @particles */
void SphereModus::initial_conditions(Particles *particles) {
  size_t number_total = 0;
  double time_start = 1.0;
  FourVector momentum_total(0, 0, 0, 0);
  /* loop over all the particle types creating each particles */
  for (auto i = particles->types_cbegin(); i != particles->types_cend(); ++i) {
    /* Particles with width > 0 (resonances) do not exist in the beginning */
    if (data.width() > 0.0) continue;
    printd("%s mass: %g [GeV]\n", data.name().c_str(), data.mass());
    /* bose einstein distribution function with temperature 0.3 GeV */
    double number_density = number_density_bose(data.mass(), 0.3);
    printf("IC number density %.6g [fm^-3]\n", number_density);
    /* cast while reflecting probability of extra particle */
    size_t number = 4.0 / 3.0 * M_PI * radius_ * radius_ * radius_ *
                    number_density * testparticles;
    if (4.0 / 3.0 * M_PI * radius_ * radius_ * radius_ * number_density -
            number >
        drand48())
      number++;
    /* create bunch of particles */
    printf("IC creating %zu particles\n", number);
    particles->create(number, data.pdgcode());
    number_total += number;
  }
  printf("IC contains %zu particles\n", number_total);
  /* now set position and momentum of the particles */
  double momentum_radial;
  Angles phitheta = Angles();
  for (ParticleData &data : particles->data()) {
    if (unlikely(data.id() == particles->id_max() && !(data.id() % 2))) {
      /* poor last guy just sits around */
      data.set_momentum(particles->particle_type(data.pdgcode()).mass(), 0, 0, 0);
    } else if (!(data.id() % 2)) {
      /* thermal momentum according Maxwell-Boltzmann distribution */
      momentum_radial = sample_momenta(0.3, particles->particle_type(data.pdgcode()).mass());
      phitheta = Angles().distribute_isotropically();
      printd("Particle %d radial momenta %g phi %g cos_theta %g\n", data.id(),
             momentum_radial, phitheta.phi(), phitheta.costheta());
      data.set_momentum(
          particles->particle_type(data.pdgcode()).mass(), momentum_radial * phitheta.x(),
          momentum_radial * phitheta.y(), momentum_radial * phitheta.z());
    } else {
      data.set_momentum(particles->particle_type(data.pdgcode()).mass(),
                             -particles->data(data.id() - 1).momentum().x1(),
                             -particles->data(data.id() - 1).momentum().x2(),
                             -particles->data(data.id() - 1).momentum().x3());
    }
    momentum_total += data.momentum();
    double x, y, z;
    /* ramdom position in a sphere
     * box length here has the meaning of the sphere radius
     */
    x = -radius_ + 2.0 * drand48() * radius_;
    y = -radius_ + 2.0 * drand48() * radius_;
    z = -radius_ + 2.0 * drand48() * radius_;
    /* sampling points inside of the sphere, rejected if outside */
    while (sqrt(x * x + y * y + z * z) > radius_) {
      x = -radius_ + 2.0 * drand48() * radius_;
      y = -radius_ + 2.0 * drand48() * radius_;
      z = -radius_ + 2.0 * drand48() * radius_;
    }
    data.set_position(time_start, x, y, z);
    /* IC: debug checks */
    printd_momenta(data);
    printd_position(data);
  }
  printf("IC total energy: %g [GeV]\n", momentum_total.x0());
}

/* boundary_condition - enforce specific type of boundaries */
// FourVector boundary_condition(FourVector position,
//  const SphereModus &sphere __attribute__((unused)), bool *boundary_hit) {
  /* no boundary */
//  *boundary_hit = false;
//  return position;
// }

/* check_collision_geometry - check if a collision happens between particles */
void SphereModus::check_collision_geometry(Particles *particles,
  CrossSections *cross_sections,
  std::list<int> *collision_list, ModusDefault const &parameters,
  BoxModus const &box, size_t *rejection_conflict) {
  std::vector<std::vector<std::vector<std::vector<int> > > > grid;
  int N, x, y, z;
  /* the maximal radial propagation for light particle */
  int a = box.length + particles->time();
  /* For small boxes no point in splitting up in grids */
  /* calculate approximate grid size according to double interaction length */
  N = round(2.0 * a / sqrt(parameters.cross_section() * fm2_mb * M_1_PI) * 0.5);
  /* for small boxes not possible to split upo */
  if (unlikely(N < 4 || particles->size() < 10)) {
    FourVector distance;
    double radial_interaction = sqrt(parameters.cross_section() * fm2_mb
                                     * M_1_PI) * 2;
  for (const ParticleData &data : particles->data()) {
    for (const ParticleData &data2 : particles->data()) {
        /* exclude check on same particle and double counting */
        if (data.id() >= data2.id())
          continue;
        distance = data.position() - data2.position();
        /* skip particles that are double interaction radius length away */
        if (distance > radial_interaction)
           continue;
        collision_criteria_geometry(particles, cross_sections, collision_list,
          parameters, data.id(), data2.id(), rejection_conflict);
      }
    }
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
  for (const ParticleData &data : particles->data()) {
    /* XXX: function - map particle position to grid number */
    x = round((a + data.position().x1()) / (N - 1));
    y = round((a + data.position().x2()) / (N - 1));
    z = round((a + data.position().x3()) / (N - 1));
    printd_position(data);
    printd("grid cell particle %i: %i %i %i of %i\n", data.id(), x, y, z, N);
    grid[x][y][z].push_back(data.id());
  }
  /* semi optimised nearest neighbour search:
   * http://en.wikipedia.org/wiki/Cell_lists
   */
  FourVector shift;
  for (const ParticleData &data : particles->data()) {
    /* XXX: function - map particle position to grid number */
    x = round((a + data.position().x1()) / (N - 1));
    y = round((a + data.position().x2()) / (N - 1));
    z = round((a + data.position().x3()) / (N - 1));
    if (unlikely(x >= N || y >= N || z >= N))
      printf("grid cell particle %i: %i %i %i of %i\n", data.id(), x, y, z, N);
    /* check all neighbour grids */
    for (int cx = -1; cx < 2; cx++) {
      int sx = cx + x;
      if (sx < 0 || sx >= N)
        continue;
      for (int cy = -1; cy <  2; cy++) {
        int sy = cy + y;
        if (sy < 0 || sy >= N)
          continue;
        for (int cz = -1; cz < 2; cz++) {
          int sz = cz + z;
          if (sz < 0 || sz >= N)
            continue;
          /* empty grid cell */
          if (grid[sx][sy][sz].empty())
            continue;
          /* grid cell particle list */
          for (auto id_b = grid[sx][sy][sz].begin();
               id_b != grid[sx][sy][sz].end(); ++id_b) {
            /* only check against particles above current id
             * to avoid double counting
             */
            if (*id_b <= data.id())
              continue;
            printd("grid cell particle %i <-> %i\n", data.id(), *id_b);
            collision_criteria_geometry(particles, cross_sections,
              collision_list, parameters, data.id(), *id_b, rejection_conflict);
          } /* grid particles loop */
        } /* grid sy */
      } /* grid sx */
    } /* grid sz */
  } /* outer particle loop */
}


/* Evolve - the core of the box, stepping forward in time */
int SphereModus::Evolve(Particles *particles, CrossSections *cross_sections,
                        int *resonances, int *decays) {
  std::list<int> collision_list, decay_list;
  size_t interactions_total = 0, previous_interactions_total = 0,
    interactions_this_interval = 0;
  size_t rejection_conflict = 0;
  /* startup values */
  print_measurements(*particles, interactions_total,
                     interactions_this_interval, ball);
  for (int steps = 0; steps < parameters.steps(); steps++) {
    /* Check resonances for decays */
    check_decays(particles, &decay_list, parameters);
    /* Do the decays */
    if (!decay_list.empty()) {
      (*decays) += decay_list.size();
      interactions_total = decay_particles(particles, &decay_list,
        interactions_total);
    }
    /* fill collision table by cells */
    check_collision_geometry(particles, cross_sections, &collision_list,
      parameters, lab, &rejection_conflict);
    /* particle interactions */
    if (!collision_list.empty())
      interactions_total = collide_particles(particles, &collision_list,
        interactions_total, resonances);
    /* propagate all particles */
    propagate_particles(particles, parameters, lab);
    /* physics output during the run */
    if (steps > 0 && (steps + 1) % parameters.output_interval() == 0) {
      interactions_this_interval = interactions_total
        - previous_interactions_total;
      previous_interactions_total = interactions_total;
      print_measurements(*particles, interactions_total,
                         interactions_this_interval, lab);
      printd("Resonances: %i Decays: %i\n", *resonances, *decays);
      printd("Ignored collisions %zu\n", rejection_conflict);
      /* save evolution data */
      write_particles(*particles);
      write_vtk(*particles);
    }
  }
  /* Guard against evolution */
  if (likely(parameters.steps > 0)) {
    /* if there are not particles no interactions happened */
    if (likely(!particles->empty()))
      print_tail(lab, interactions_total * 2
                 / particles->time() / particles->size());
    else
      print_tail(ball, 0);
    printf("Total ignored collisions: %zu\n", rejection_conflict);
  }
  return 0;
}

/* start up a sphere and run it */
// int Sphere::evolve(const Laboratory &lab, char *path) {
//  /* Read sphere config file parameters */
//  Sphere *ball = new Sphere(lab);
//  process_config_sphere(ball, path);
  /* Initialize box */
//  print_startup(*ball);
//  Particles *particles = new Particles;
//  input_particles(particles, path);
//  initial_conditions(particles, ball);
// input_decaymodes(particles, path);
//  CrossSections *cross_sections = new CrossSections;
//  cross_sections->add_elastic_parameter(lab.cross_section());
  /* Compute stuff */
//  int rc = Evolve(particles, cross_sections, lab);
  /* record IC startup */
//  write_measurements_header(*particles);
//  print_header();
//  write_particles(*particles);
  /* tear down */
//  delete particles;
// delete cross_sections;
// delete ball;
//  return rc;
// }


}  // namespace Smash
