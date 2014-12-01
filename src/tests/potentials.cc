/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <map>
#include "unittest.h"

#include "../include/configuration.h"
#include "../include/potentials.h"
#include "../include/experiment.h"
#include "../include/modusdefault.h"
#include "../include/nucleus.h"

#include <boost/filesystem.hpp>

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "proton 0.938 0.0 2212\n"
      "neutron 0.938 0.0 2112\n"
      "pion    0.138 0.0 211\n");
}

static ParticleData create_proton(int id = -1) {
  return ParticleData{ParticleType::find(0x2212), id};
}

static ParticleData create_pion(int id = -1) {
  return ParticleData{ParticleType::find(0x211), id};
}

// check that analytical and numerical results for gradient of potential coincide
TEST(potential_gradient) {
  // create two protons
  ParticleData part1 = create_proton();
  ParticleData part2 = create_proton();
  ParticleData part3 = create_pion();
  double mass = 0.938;
  // set momenta (just randomly):
  part1.set_4momentum(mass, ThreeVector(1.0, 4.0, 3.0));
  part2.set_4momentum(mass, ThreeVector(4.0, 2.0, 5.0));
  part3.set_4momentum(0.138, ThreeVector(2.0, -1.0, 1.0));
  // set coordinates (not too far from each other)
  part1.set_4position(FourVector(0.0, -0.5, 0.0, 0.0));
  part2.set_4position(FourVector(0.0, 0.5, 0.0, 0.0));
  part3.set_4position(FourVector(0.0, 0.0, 0.4, 0.0));
  // make particle list out of them
  ParticleList P;
  P.push_back(part1);
  P.push_back(part2);
  P.push_back(part3);

  ThreeVector r,dr;
  Configuration conf(TEST_CONFIG_PATH);
  conf["Potentials"]["Skyrme"]["Skyrme_A"] = -209.2;
  conf["Potentials"]["Skyrme"]["Skyrme_B"] = 156.4;
  conf["Potentials"]["Skyrme"]["Skyrme_Tau"] = 1.35;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 1, 1.0};
  Potentials* pot = new Potentials(conf["Potentials"], param);

  ThreeVector num_grad, analit_grad;
  r = ThreeVector(0.2, 0.0, 0.0);

  // analytical gradient
  analit_grad = pot->potential_gradient(r, P);
  // numerical gradient
  const double U = pot->potential(r, P);
  dr = ThreeVector(1.e-4, 0.0, 0.0);
  num_grad.set_x1((pot->potential(r + dr, P) - U)/dr.x1());
  dr = ThreeVector(0.0, 1.e-4, 0.0);
  num_grad.set_x2((pot->potential(r + dr, P) - U)/dr.x2());
  dr = ThreeVector(0.0, 0.0, 1.e-4);
  num_grad.set_x3((pot->potential(r + dr, P) - U)/dr.x3());
  // compare them with: accuracy should not be worse than |dr|
  std::cout << num_grad << analit_grad << std::endl;
  COMPARE_ABSOLUTE_ERROR(num_grad.x1(), analit_grad.x1(), 1.e-4);
  COMPARE_ABSOLUTE_ERROR(num_grad.x2(), analit_grad.x2(), 1.e-4);
  COMPARE_ABSOLUTE_ERROR(num_grad.x3(), analit_grad.x3(), 1.e-4);
}
