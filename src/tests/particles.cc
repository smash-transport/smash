/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "setup.h"

#include "../include/particles.h"
#include "../include/particledata.h"
#include "../include/pdgcode.h"

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

TEST(empty_particles) {
  Particles p;
  VERIFY(p.is_empty());
  COMPARE(p.size(), 0u);

  ParticleData pd = Test::smashon();
  VERIFY(!p.is_valid(pd));
}

TEST(create) {
  Particles p;
  p.create(0x661);

  VERIFY(!p.is_empty());
  COMPARE(p.size(), 1u);

  ParticleData pd = Test::smashon();
  VERIFY(!p.is_valid(pd));

  pd = p.front();
  VERIFY(p.is_valid(pd));
  COMPARE(pd.id(), 0);
}

TEST(replace) {
  Particles p;
  p.create(0x661);
  COMPARE(p.size(), 1u);
  ParticleList to_remove = {p.front()};
  ParticleList to_add = {Test::smashon(-1)};
  VERIFY(p.is_valid(to_remove.front()));
  VERIFY(!p.is_valid(to_add.front()));

  COMPARE(to_remove.front().id(), 0);
  COMPARE(to_add.front().id(), -1);
  p.replace(to_remove, to_add);
  COMPARE(p.size(), 1u);
  COMPARE(to_add.front().id(), -1);
  COMPARE(p.front().id(), 1);
  VERIFY(!p.is_valid(to_add.front()));
  VERIFY(!p.is_valid(to_remove.front()));
  VERIFY(p.is_valid(p.front()));

  to_remove = {p.front()};
  to_add = {Test::smashon(), Test::smashon(), Test::smashon()};
  p.replace(to_remove, to_add);
  COMPARE(p.size(), 3u);
  COMPARE(p.front().id(), 2);
  COMPARE(p.back().id(), 4);
}

TEST(insert) {
  Particles p;
  p.insert(Test::smashon());
  COMPARE(p.size(), 1u);
  VERIFY(p.is_valid(p.front()));
  p.insert(Test::smashon());
  COMPARE(p.size(), 2u);
  VERIFY(p.is_valid(p.front()));
  VERIFY(p.is_valid(p.back()));
  p.create(0x661);
  COMPARE(p.size(), 3u);
  VERIFY(p.is_valid(p.front()));
  VERIFY(p.is_valid(p.back()));
  p.insert(Test::smashon());
  COMPARE(p.size(), 4u);
  VERIFY(p.is_valid(p.front()));
  VERIFY(p.is_valid(p.back()));

  ParticleData smashon = Test::smashon();
  smashon.set_id_process(1);
  p.insert(smashon);
  COMPARE(p.back().id_process(), 1);
  smashon.set_id_process(2);
  p.insert(smashon);
  COMPARE(p.back().id_process(), 2);
}

TEST(insert_2) {
  /* 2 particles with null momenta */
  const auto particle_a = Test::smashon(Test::Position{1., 1., 1., 1.},
                                        Test::Momentum{0.1, 0.0, 0.0, 0.0}, 0);
  const auto particle_b = Test::smashon(Test::Position{2., 2., 2., 2.},
                                        Test::Momentum{0.1, 0.0, 0.0, 0.0}, 1);

  Particles particles;

  particles.insert(particle_a);
  COMPARE(particles.size(), 1u);
  particles.insert(particle_b);
  COMPARE(particles.size(), 2u);
  VERIFY(particles.is_valid(particles.front())) << particles.front();
  VERIFY(particles.is_valid(particles.back())) << particles.back();

  COMPARE(particles.front().id(), particle_a.id());
  COMPARE(particles.back() .id(), particle_b.id());
  COMPARE(particles.front().pdgcode(), particle_a.pdgcode());
  COMPARE(particles.back() .pdgcode(), particle_b.pdgcode());
  COMPARE(particles.front().momentum(), particle_a.momentum());
  COMPARE(particles.back() .momentum(), particle_b.momentum());
  COMPARE(particles.front().position(), particle_a.position());
  COMPARE(particles.back() .position(), particle_b.position());
}

TEST(create_multiple) {
  Particles p;
  p.create(4, 0x661);
  COMPARE(p.size(), 4u);
  COMPARE(p.front().id(), 0);
  p.remove(p.front());
  COMPARE(p.size(), 3u);
  COMPARE(p.front().id(), 1);
  p.remove(p.front());
  COMPARE(p.size(), 2u);
  COMPARE(p.front().id(), 2);
}

template <typename T>
void check_particle_data_iteration(T *p, std::size_t expected_size) {
  std::size_t count = 0;
  for (auto &data : *p) {
    VERIFY(p->is_valid(data)) << count;
    ++count;
  }
  COMPARE(count, p->size());
  COMPARE(count, expected_size);
}

TEST(iterate_particle_data) {
  Particles p;
  const Particles *p2 = &p;
  check_particle_data_iteration(&p, 0);
  check_particle_data_iteration(p2, 0);
  p.create(0x211);
  check_particle_data_iteration(&p, 1);
  check_particle_data_iteration(p2, 1);
  p.create(-0x211);
  check_particle_data_iteration(&p, 2);
  check_particle_data_iteration(p2, 2);
  p.create(0x211);
  check_particle_data_iteration(&p, 3);
  check_particle_data_iteration(p2, 3);
  p.replace({p.front(), p.back()}, {Test::smashon()});
  check_particle_data_iteration(&p, 2);
  check_particle_data_iteration(p2, 2);
  p.create(-0x211);
  p.replace({p.front(), *(++p.begin())}, {Test::smashon()});
  check_particle_data_iteration(&p, 2);
  check_particle_data_iteration(p2, 2);
}

TEST(erase_particle) {
  Particles p;
  p.create(0x211);
  p.create(-0x211);
  p.create(0x111);
  COMPARE(p.size(), 3u);
  for (auto &&x : p) {
    VERIFY(p.is_valid(x));
  }

  auto copy = p.front();
  VERIFY(p.is_valid(copy));
  p.remove(copy);
  COMPARE(p.size(), 2u);
  VERIFY(!p.is_valid(copy));
  for (auto &&x : p) {
    VERIFY(p.is_valid(x));
  }

  auto copy2 = p.back();
  p.remove(copy2);
  COMPARE(p.size(), 1u);
  VERIFY(!p.is_valid(copy));
  VERIFY(!p.is_valid(copy2));
  for (auto &&x : p) {
    VERIFY(p.is_valid(x));
  }
}

TEST(reset) {
  Particles p;
  std::size_t count = 1000;
  for (auto i = count; i; --i) {
    p.insert(Test::smashon());
  }
  COMPARE(p.size(), count);
  p.remove(p.front());
  COMPARE(p.size(), count - 1);
  VERIFY(!p.is_empty());
  p.reset();
  COMPARE(p.size(), 0u);
  VERIFY(p.is_empty());
  for (auto i = count; i; --i) {
    p.insert(Test::smashon());
  }
  COMPARE(p.size(), count);
  int n = 0;
  for (auto &&pd : p) {
    COMPARE(pd.id(), n);
    ++n;
  }
}

TEST(copy_to_vector) {
  Particles p;
  p.create(100, 0x661);
  auto copy = p.copy_to_vector();
  COMPARE(copy.size(), 100u);
  for (auto &&x : copy) {
    VERIFY(p.is_valid(x));
  }

  p.remove(copy[5]);
  p.remove(copy[40]);
  p.remove(copy[41]);
  copy = p.copy_to_vector();
  COMPARE(copy.size(), 97u);
  for (auto &&x : copy) {
    VERIFY(p.is_valid(x));
  }
}

TEST(exceed_capacity) {
  Particles p;
  p.create(50, 0x661);
  COMPARE(p.size(), 50u);
  p.create(150, 0x661);
  COMPARE(p.size(), 200u);
  p.create(450, 0x661);
  COMPARE(p.size(), 650u);
  p.create(1350, 0x661);
  COMPARE(p.size(), 2000u);
  p.create(4050, 0x661);
  COMPARE(p.size(), 6050u);
  p.create(12150, 0x661);
  COMPARE(p.size(), 18200u);
  int n = 0;
  for (auto &&x : p) {
    COMPARE(x.id(), n);
    ++n;
  }
}

TEST(update) {
  Particles p;
  auto pd =
      Test::smashon(Test::Momentum{1, 1, 1, 1}, Test::Position{1, 1, 1, 1});
  pd.set_id_process(1);
  p.insert(pd);
  p.insert(pd);
  p.insert(pd);
  COMPARE(p.size(), 3u);
  COMPARE(p.front().momentum(), FourVector(1, 1, 1, 1));
  COMPARE(p.front().position(), FourVector(1, 1, 1, 1));
  COMPARE(p.front().id_process(), 1);
  pd.set_id_process(2);
  pd.set_4momentum({2, 2, 2, 2});
  pd.set_4position({3, 3, 3, 3});
  p.update(p.front(), pd);
  COMPARE(p.size(), 3u);
  COMPARE(p.front().momentum(), FourVector(2, 2, 2, 2));
  COMPARE(p.front().position(), FourVector(3, 3, 3, 3));
  COMPARE(p.front().id_process(), 2);
}
