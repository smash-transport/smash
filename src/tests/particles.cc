/*
 *    Copyright (c) 2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *    GNU General Public License (GPLv3)
 */

#include "tests/unittest.h"

#include <cstdio>

#include "include/particles.h"
#include "include/constants.h"
#include "include/particledata.h"
#include "include/outputroutines.h"
#include <algorithm>

TEST(everything) {
  /* checks for geometric distance criteria */
  ParticleData particle_a, particle_b;
  particle_a.set_id(0);
  particle_b.set_id(1);

  /* 2 particles with null momenta */
  particle_a.set_momentum(0.1, 0.0, 0.0, 0.0);
  particle_b.set_momentum(0.1, 0.0, 0.0, 0.0);
  particle_a.set_position(1.0, 1.0, 1.0, 1.0);
  particle_b.set_position(2.0, 2.0, 2.0, 2.0);

  /* check return of particle distance of null momenta particles */
  double distance_squared = particle_distance(&particle_a, &particle_b);
  VERIFY(!(distance_squared < 0.0));
  /* XXX: does this test NaN?? */
  VERIFY(!(distance_squared > 1000.0));

  /* check collision_time for parallel momenta => impossible collision */
  particle_a.set_momentum(0.1, 0.3, -0.1, 0.2);
  particle_b.set_momentum(0.1, 0.3, -0.1, 0.2);
  double time = collision_time(particle_a, particle_b);
  VERIFY(!(time >= 0.0));

  /* reset momenta for possible collision and compare to Particle class */
  particle_a.set_momentum(0.1, 10.0, 9.0, 8.0);
  particle_b.set_momentum(0.1, -10.0, -90.0, -80.0);
  distance_squared = particle_distance(&particle_a, &particle_b);
  VERIFY(!(distance_squared < 0.0));

  /* now check the Particles class itself */
  Particles particles;
  /* check addition of particle types */
  ParticleType piplus("pi+", 0.13957, -1.0, 211, 1, 1, 0);
  particles.add_type(piplus, 211);
  VERIFY(!(particles.types_size() != 1));
  size_t type_size = 0;
  for (std::map<int, ParticleType>::const_iterator
       i = particles.types_cbegin(); i != particles.types_cend(); ++i) {
    printd("pdg %d mass: %g [GeV]\n", i->first, i->second.mass());
    type_size++;
  }
  VERIFY(!(type_size != 1));
  type_size = 0;
  ParticleType piminus("pi-", 0.13957, -1.0, -211, 1, -1, 0);
  particles.add_type(piminus, -211);
  for (std::map<int, ParticleType>::const_iterator
       i = particles.types_cbegin(); i != particles.types_cend(); ++i) {
    printd("pdg %d mass: %g [GeV]\n", i->first, i->second.mass());
    type_size++;
  }
  VERIFY(!(type_size != 2));
  particles.add_type(piminus, -211);
  VERIFY(!(particles.types_size() != 2));

  /* check addition of particles */
  particles.add_data(particle_a);
  VERIFY(!(particles.size() != 1));
  particles.add_data(particle_b);
  VERIFY(!(particles.size() != 2));
  type_size = 0;
  for (std::map<int, ParticleData>::const_iterator
       i = particles.cbegin(); i != particles.cend(); ++i) {
    printd("id %d: pdg %d\n", i->first, i->second.pdgcode());
    /* check that id and and own id are the same */
    VERIFY(!(i->first != i->second.id()));
    type_size++;
  }
  VERIFY(!(type_size != 2));
  VERIFY(!(particles.empty()));
  VERIFY(!(particles.count(1) != 1));

  /* check usage particle data */
  double distance_squared_2 = particle_distance(particles.data_pointer(0),
    particles.data_pointer(1));
  printd("%g versus %g\n", distance_squared, distance_squared_2);
  VERIFY(!(distance_squared_2 < 0.0));
  VERIFY(!(distance_squared_2 - distance_squared > really_small));
}

TEST_CATCH(load_from_incorrect_string, Particles::LoadFailure) {
  Particles p;
  std::istringstream input("Hallo Welt! (wave)");
  p.load(input);
}

TEST(load_one_particle_no_extra_whitespace) {
  Particles p;
  std::istringstream input1("pi0 0.1350 -1.0 111 2 0 0");
  p.load(input1);
  COMPARE(p.types_size(), 1);
  int count = 0;
  std::for_each(p.types_cbegin(), p.types_cend(),
                [&count](std::pair<int, ParticleType> i) {
    ++count;
    const ParticleType &type = i.second;
    COMPARE(type.mass(), 0.135f);
    COMPARE(type.width(), -1.f);
    COMPARE(type.pdgcode(), 111);
    COMPARE(type.isospin(), 2);
    COMPARE(type.charge(), 0);
    COMPARE(type.spin(), 0);
  });
  COMPARE(count, 1);
}

TEST(load_one_particle_with_whitespace) {
  Particles p;
  std::istringstream input("\t\n\t  pi0  0.1350 \t -1.0 111\t2 0 0 \n ");
  p.load(input);
  COMPARE(p.types_size(), 1);
  int count = 0;
  std::for_each(p.types_cbegin(), p.types_cend(),
                [&count](std::pair<int, ParticleType> i) {
    ++count;
    const ParticleType &type = i.second;
    COMPARE(type.mass(), 0.135f);
    COMPARE(type.width(), -1.f);
    COMPARE(type.pdgcode(), 111);
    COMPARE(type.isospin(), 2);
    COMPARE(type.charge(), 0);
    COMPARE(type.spin(), 0);
  });
  COMPARE(count, 1);
}

TEST_CATCH(load_one_particle_with_incorrect_newline, Particles::LoadFailure) {
  Particles p;
  std::istringstream input("pi0 0.1350\n-1.0 111 2 0 0");
  p.load(input);
}

TEST(load_only_comments) {
  Particles p;
  std::istringstream input(
      "# Hello\n"
      "  # Will you ignore me? #### sldfkjsdf\n"
      "\t\t  \t # yes?");
  p.load(input);
  COMPARE(p.types_size(), 0);
}

TEST(load_one_particle_with_comment) {
  Particles p;
  std::istringstream input("pi0 0.1350  -1.0 111 2 0 0 # This is pi0. Swell.");
  p.load(input);
  COMPARE(p.types_size(), 1);
  int count = 0;
  std::for_each(p.types_cbegin(), p.types_cend(),
                [&count](std::pair<int, ParticleType> i) {
    ++count;
    const ParticleType &type = i.second;
    COMPARE(type.mass(), 0.135f);
    COMPARE(type.width(), -1.f);
    COMPARE(type.pdgcode(), 111);
    COMPARE(type.isospin(), 2);
    COMPARE(type.charge(), 0);
    COMPARE(type.spin(), 0);
  });
  COMPARE(count, 1);
}

TEST(load_many_particles) {
  Particles p;
  std::istringstream input(
      "# NAME MASS[GEV] WIDTH[GEV] PDG ISOSPIN CHARGE SPIN\n"
      "pi0 0.1350 -1.0 111 2 0 0\n"
      "pi+ 0.1396 -1.0 211 2 1 0\n"
      "pi- 0.1396 -1.0 -211 2 -1 0\n"
      "rho0 0.7755 0.149 113 2 0 2\n"
      "rho+ 0.7755 0.149 213 2 1 2\n"
      "rho- 0.7755 0.149 -213 2 -1 2\n"
      "eta 0.5479 1.0e-6 221 0 0 0\n"
      "omega 0.7827 0.0085 223 0 0 2\n"
      "p 0.9383 -1.0 2212 1 1 1\n"
      "pbar 0.9383 -1.0 -2212 1 -1 1\n"
      "n 0.9396 -1.0 2112 1 0 1\n"
      "nbar 0.9396 -1.0 -2112 1 0 1\n"
      "Delta++ 1.232 0.117 2224 3 2 3\n"
      "Delta+ 1.232 0.117 2214 3 1 3\n"
      "Delta0 1.232 0.117 2114 3 0 3\n"
      "Delta- 1.232 0.117 1114 3 -1 3\n"
      "Deltabar++ 1.232 0.117 -2224 3 -2 3\n"
      "Deltabar+ 1.232 0.117 -2214 3 -1 3\n"
      "Deltabar0 1.232 0.117 -2114 3 0 3\n"
      "Deltabar- 1.232 0.117 -1114 3 1 3\n");
  p.load(input);
  COMPARE(p.types_size(), 20);
  ParticleType type = p.particle_type(-1114);
  COMPARE(type.mass(), 1.232f);
  COMPARE(type.width(), .117f);
  COMPARE(type.pdgcode(), -1114);
  COMPARE(type.isospin(), 3);
  COMPARE(type.charge(), 1);
  COMPARE(type.spin(), 3);

  type = p.particle_type(2112);
  COMPARE(type.mass(), .9396f);
  COMPARE(type.width(), -1.f);
  COMPARE(type.pdgcode(), 2112);
  COMPARE(type.isospin(), 1);
  COMPARE(type.charge(), 0);
  COMPARE(type.spin(), 1);
}
