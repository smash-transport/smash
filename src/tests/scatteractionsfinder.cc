/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"

#include <cstdio>

#include "../include/particledata.h"
#include "../include/pdgcode.h"
#include "../include/action.h"
#include "../include/scatteractionsfinder.h"

using namespace Smash;

// ------
// everything copied from the test action.cc
// ------

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "smashon 0.123 1.2 661\n");
}

static ParticleData create_smashon_particle(int id = -1) {
  return ParticleData{ParticleType::find(0x661), id};
}

TEST(two_particles) {
  ParticleData particle_a = create_smashon_particle(0),
               particle_b = create_smashon_particle(1);
  particle_a.set_4position(FourVector(0., 1., 0., 0.));
  particle_b.set_4position(FourVector(0., 0., 1., 0.));

  // set momenta
  particle_a.set_4momentum(0.1, 0., 0.5, 0.2);
  particle_b.set_4momentum(0.1, 0.5, 0., 0.2);

  // put in list
  ParticleList l;
  l.emplace_back(particle_a);
  l.emplace_back(particle_b);

  // prepare scatteractionsfinder
  float elastic_parameter = 1;
  int testparticles = 1;
  // delta t in fermi
  float dt;
  ScatterActionsFinder finder(elastic_parameter, testparticles);

  // no neighbors
  std::vector<const ParticleList*> neighbors_list;
  std::vector<ActionPtr> actions;

  dt = 0.9f;
  actions = finder.find_possible_actions(
          l, neighbors_list, dt);
  COMPARE(actions.size(), 0) << "there shouldn't be a collision yet";

  dt = 1.1f;
  actions = finder.find_possible_actions(
          l, neighbors_list, dt);
  COMPARE(actions.size(), 1) << "there should already be a collision";
}
