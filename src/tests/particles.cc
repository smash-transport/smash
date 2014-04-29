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

using namespace Smash;

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
  const std::string pis(
      "pi+ 0.13957 -1.0 211 1 1 0\n"
      "pi- 0.13957 -1.0 -211 1 -1 0\n");
  Particles particles(pis, {});

  /* check addition of particles */
  particles.add_data(particle_a);
  VERIFY(!(particles.size() != 1));
  particles.add_data(particle_b);
  VERIFY(!(particles.size() != 2));
  int type_size = 0;
  for (const ParticleData &data : particles.data()) {
    printd("id %d: pdg %d\n", data.id(), data.pdgcode());
    type_size++;
  }
  VERIFY(!(type_size != 2));
  VERIFY(!(particles.empty()));
  VERIFY(particles.has_data(1));

  /* check usage particle data */
  double distance_squared_2 = particle_distance(particles.data_pointer(0),
    particles.data_pointer(1));
  printd("%g versus %g\n", distance_squared, distance_squared_2);
  VERIFY(!(distance_squared_2 < 0.0));
  VERIFY(!(distance_squared_2 - distance_squared > really_small));
}

TEST_CATCH(load_from_incorrect_string, Particles::LoadFailure) {
  const std::string parts("Hallo Welt! (wave)");
  Particles p(parts, {});
}

TEST(load_one_particle_no_extra_whitespace) {
  const std::string parts("pi0 0.1350 -1.0 111 2 0 0");
  Particles p(parts, {});
  COMPARE(p.types_size(), 1u);
  int count = 0;
  for (const auto &type : p.types()) {
    ++count;
    COMPARE(type.mass(), 0.135f);
    COMPARE(type.width(), -1.f);
    COMPARE(type.pdgcode(), 111);
    COMPARE(type.isospin(), 2);
    COMPARE(type.charge(), 0);
    COMPARE(type.spin(), 0);
  }
  COMPARE(count, 1);
}

TEST(load_one_particle_with_whitespace) {
  const std::string parts("\t\n\t  pi0  0.1350 \t -1.0 111\t2 0 0 \n ");
  Particles p(parts, {});
  COMPARE(p.types_size(), 1u);
  int count = 0;
  for (const auto &type : p.types()) {
    ++count;
    COMPARE(type.mass(), 0.135f);
    COMPARE(type.width(), -1.f);
    COMPARE(type.pdgcode(), 111);
    COMPARE(type.isospin(), 2);
    COMPARE(type.charge(), 0);
    COMPARE(type.spin(), 0);
  }
  COMPARE(count, 1);
}

TEST_CATCH(load_one_particle_with_incorrect_newline, Particles::LoadFailure) {
  const std::string parts("pi0 0.1350\n-1.0 111 2 0 0");
  Particles p(parts, {});
}

TEST(load_only_comments) {
  const std::string parts(
      "# Hello\n"
      "  # Will you ignore me? #### sldfkjsdf\n"
      "\t\t  \t # yes?");
  Particles p(parts, {});
  COMPARE(p.types_size(), 0u);
}

TEST(load_one_particle_with_comment) {
  const std::string parts("pi0 0.1350  -1.0 111 2 0 0 # This is pi0. Swell.");
  Particles p(parts, {});
  COMPARE(p.types_size(), 1u);
  int count = 0;
  for (const auto &type : p.types()) {
    ++count;
    COMPARE(type.mass(), 0.135f);
    COMPARE(type.width(), -1.f);
    COMPARE(type.pdgcode(), 111);
    COMPARE(type.isospin(), 2);
    COMPARE(type.charge(), 0);
    COMPARE(type.spin(), 0);
  }
  COMPARE(count, 1);
}

namespace particles_txt {
#include "particles.txt.h"
}  // namespace particles_txt
namespace decaymodes_txt {
#include "decaymodes.txt.h"
}  // namespace decaymodes_txt

TEST(load_many_particles) {
  Particles p(particles_txt::data, {});
  COMPARE(p.types_size(), 20u);
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

TEST_CATCH(load_decaymodes_missing_pdg, Particles::ReferencedParticleNotFound) {
  const std::string decays_input(
      "113 \n"
      );
  Particles p({}, decays_input);
}

TEST_CATCH(load_decaymodes_no_decays, Particles::MissingDecays) {
  const std::string decays_input(
      "113 # rho0\n"
      );
  Particles p(particles_txt::data, decays_input);
}

TEST_CATCH(load_decaymodes_incorrect_start, Particles::ParseError) {
  const std::string decays_input(
      "113. # rho0\n"
      );
  Particles p(particles_txt::data, decays_input);
}

TEST(load_decaymodes_two_channels) {
  const std::string decays_input(
      " 113\t# rho0\n"
      "\n"
      " 1.0\t211 -211\t# pi+ pi- \n"
      " \n"
      "\n"
      "223	# omega\n"
      "0.33 111 113	# pi0 rho0\n"
      "\n"
      "0.33 211 -213	# pi+ rho-\n"
      "0.33 -211 213	# pi- rho+\n"
      );
  Particles p(particles_txt::data, decays_input);

  {
    const auto &rho0 = p.decay_modes(113);
    VERIFY(!rho0.empty());
    const auto &modelist = rho0.decay_mode_list();
    COMPARE(modelist.size(), 1u);
    COMPARE(modelist[0].weight(), 1.);
    COMPARE(modelist[0].particle_list().size(), 2u);
    COMPARE(modelist[0].particle_list()[0], 211);
    COMPARE(modelist[0].particle_list()[1], -211);
  }
  {
    const auto &omega = p.decay_modes(223);
    VERIFY(!omega.empty());
    const auto &modelist = omega.decay_mode_list();
    COMPARE(modelist.size(), 3u);
    FUZZY_COMPARE(float(modelist[0].weight()), 1.f/3.f);
    FUZZY_COMPARE(float(modelist[1].weight()), 1.f/3.f);
    FUZZY_COMPARE(float(modelist[2].weight()), 1.f/3.f);
    COMPARE(modelist[0].particle_list().size(), 2u);
    COMPARE(modelist[0].particle_list()[0], 111);
    COMPARE(modelist[0].particle_list()[1], 113);
    COMPARE(modelist[1].particle_list().size(), 2u);
    COMPARE(modelist[1].particle_list()[0], 211);
    COMPARE(modelist[1].particle_list()[1], -213);
    COMPARE(modelist[2].particle_list().size(), 2u);
    COMPARE(modelist[2].particle_list()[0], -211);
    COMPARE(modelist[2].particle_list()[1], 213);
  }
}

template <typename T>
void check_particle_data_iteration(T *p) {
  std::size_t count = 0;
  for (auto &data : p->data()) {
    const int id = data.id();
    const ParticleData &data2 = p->data(id);
    COMPARE(&data, &data2);
    ++count;
  }
  COMPARE(count, p->size());
}

template <typename T>
void check_particle_type_iteration(T *p) {
  std::size_t count = 0;
  for (const auto &type : p->types()) {
    const int pdg = type.pdgcode();
    const ParticleType &type2 = p->particle_type(pdg);
    COMPARE(&type, &type2);
    ++count;
  }
  COMPARE(count, p->types_size());
}

TEST(iterate_particle_data) {
  Particles p(particles_txt::data, decaymodes_txt::data);
  const Particles *p2 = &p;
  check_particle_type_iteration(&p);
  check_particle_type_iteration(p2);

  check_particle_data_iteration(&p);
  check_particle_data_iteration(p2);
  p.create(211);
  check_particle_data_iteration(&p);
  check_particle_data_iteration(p2);
  p.create(-211);
  check_particle_data_iteration(&p);
  check_particle_data_iteration(p2);
}

TEST(erase_particle) {
  Particles p(particles_txt::data, decaymodes_txt::data);
  p.create(211);
  p.create(-211);
  p.create(111);
  COMPARE(p.size(), 3u);
  VERIFY(p.has_data(0));
  VERIFY(p.has_data(1));
  VERIFY(p.has_data(2));
  VERIFY(!p.has_data(3));
  COMPARE(p.data(1).pdgcode(), -211);

  p.remove(0);
  COMPARE(p.size(), 2u);
  VERIFY(!p.has_data(0));
  VERIFY(p.has_data(1));
  VERIFY(p.has_data(2));
  VERIFY(!p.has_data(3));
  COMPARE(p.data(1).pdgcode(), -211);

  p.remove(2);
  COMPARE(p.size(), 1u);
  VERIFY(!p.has_data(0));
  VERIFY(p.has_data(1));
  VERIFY(!p.has_data(2));
  VERIFY(!p.has_data(3));
  COMPARE(p.data(1).pdgcode(), -211);
}
