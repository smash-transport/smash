/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"

#include <cstdio>

#include "../include/particles.h"
#include "../include/constants.h"
#include "../include/particledata.h"
#include "../include/pdgcode.h"
#include "../include/logging.h"
#include "../include/macros.h"
#include <algorithm>

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "smashon 0.123 1.2 661\n"
      "pi0 0.1350 -1.0 111\n"
      "pi+ 0.1396 -1.0 211\n"
      "rho0 0.7755 0.149 113\n"
      "rho+ 0.7755 0.149 213\n"
      "eta 0.5479 1.0e-6 221\n"
      "omega 0.7827 0.0085 223\n"
      "p 0.9383 -1.0 2212\n"
      "n 0.9396 -1.0 2112\n"
      "Delta++ 1.232 0.117 2224\n"
      "Delta+ 1.232 0.117 2214\n"
      "Delta0 1.232 0.117 2114\n"
      "Delta- 1.232 0.117 1114\n");
}

static ParticleData create_smashon_particle(int id = -1) {
  return ParticleData{ParticleType::find(0x661), id};
}

TEST(everything) {
  const einhard::Logger<> log(einhard::ALL);
  ParticleData particle_a = create_smashon_particle(0),
               particle_b = create_smashon_particle(1);

  /* 2 particles with null momenta */
  particle_a.set_4momentum(0.1, 0.0, 0.0, 0.0);
  particle_b.set_4momentum(0.1, 0.0, 0.0, 0.0);
  particle_a.set_4position(FourVector(1., 1., 1., 1.));
  particle_b.set_4position(FourVector(2., 2., 2., 2.));

  /* check the Particles class itself */
  Particles particles;

  /* check addition of particles */
  particles.add_data(particle_a);
  VERIFY(particles.size() == 1);
  particles.add_data(particle_b);
  VERIFY(particles.size() == 2);
  int data_size = 0;
  for (const ParticleData &data : particles.data()) {
    log.debug("id ", data.id(), ": pdg ", data.pdgcode());
    SMASH_UNUSED(data);
    data_size++;
  }
  VERIFY(data_size == 2);
  VERIFY(!particles.empty());
  VERIFY(particles.has_data(0));
  VERIFY(particles.has_data(1));

  /* check usage of particles.data */
  VERIFY(particles.data(0).id()==particle_a.id());
  VERIFY(particles.data(1).id()==particle_b.id());
  VERIFY(particles.data(0).pdgcode()==particle_a.pdgcode());
  VERIFY(particles.data(1).pdgcode()==particle_b.pdgcode());
  VERIFY(particles.data(0).momentum()==particle_a.momentum());
  VERIFY(particles.data(1).momentum()==particle_b.momentum());
  VERIFY(particles.data(0).position()==particle_a.position());
  VERIFY(particles.data(1).position()==particle_b.position());
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

static void check_particle_type_iteration() {
  std::size_t count = 0;
  for (const auto &type : ParticleType::list_all()) {
    const PdgCode pdg = type.pdgcode();
    const ParticleType &type2 = ParticleType::find(pdg);
    COMPARE(&type, &type2);
    ++count;
  }
  COMPARE(count, ParticleType::list_all().size());
}

TEST(iterate_particle_data) {
  Particles p;
  const Particles *p2 = &p;
  check_particle_type_iteration();
  check_particle_type_iteration();

  check_particle_data_iteration(&p);
  check_particle_data_iteration(p2);
  p.create(0x211);
  check_particle_data_iteration(&p);
  check_particle_data_iteration(p2);
  p.create(-0x211);
  check_particle_data_iteration(&p);
  check_particle_data_iteration(p2);
}

TEST(erase_particle) {
  Particles p;
  p.create(0x211);
  p.create(-0x211);
  p.create(0x111);
  COMPARE(p.size(), 3u);
  VERIFY(p.has_data(0));
  VERIFY(p.has_data(1));
  VERIFY(p.has_data(2));
  VERIFY(!p.has_data(3));
  COMPARE(p.data(1).pdgcode().dump(), 0x80000211u);

  p.remove(0);
  COMPARE(p.size(), 2u);
  VERIFY(!p.has_data(0));
  VERIFY(p.has_data(1));
  VERIFY(p.has_data(2));
  VERIFY(!p.has_data(3));
  COMPARE(p.data(1).pdgcode().dump(), 0x80000211u);

  p.remove(2);
  COMPARE(p.size(), 1u);
  VERIFY(!p.has_data(0));
  VERIFY(p.has_data(1));
  VERIFY(!p.has_data(2));
  VERIFY(!p.has_data(3));
  COMPARE(p.data(1).pdgcode().dump(), 0x80000211u);
}
