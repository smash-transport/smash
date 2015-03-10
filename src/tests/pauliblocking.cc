/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/boxmodus.h"
#include "../include/configuration.h"
#include "../include/experiment.h"
#include "../include/pauliblocking.h"
#include "../include/potentials.h"

#include <boost/filesystem.hpp>

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "mock_De 0.1 0.0 2114\n");
}

TEST(phase_space_density) {
  Configuration conf(TEST_CONFIG_PATH);
  conf["General"]["Randomseed"] = -1;
  conf["Modi"]["Box"]["Initial_Condition"] = 1;
  conf["Modi"]["Box"]["Length"] = 5.0;
  conf["Modi"]["Box"]["Temperature"] = 1.0;
  conf["Modi"]["Box"]["Start_Time"] = 0.0;
  conf.take({"Modi", "Box", "Init_Multiplicities"});
  conf["Modi"]["Box"]["Init_Multiplicities"]["2114"] = 10000;
  conf["Collision_Term"]["Pauli_Blocking"]["Spatial_Averaging_Radius"] = 1.86;
  conf["Collision_Term"]["Pauli_Blocking"]["Momentum_Averaging_Radius"] = 2.0;
  conf["Collision_Term"]["Pauli_Blocking"]["Gaussian_Cutoff"] = 2.2;
  // Number of test-particles
  const size_t ntest = 100;
  ExperimentParameters param{{0.f, 1.f}, 1.f, ntest, 1.0};
  BoxModus *b = new BoxModus(conf["Modi"], param);
  Particles Pbox;
  b->initial_conditions(&Pbox, param);
  PauliBlocker *pb = new PauliBlocker(conf, param);
  ThreeVector r(0.0, 0.0, 0.0);
  ThreeVector p;
  PdgCode pdg = 0x2114;
  float f;
  for (int i = 1; i < 100; i++) {
    p = ThreeVector(0.0, 0.0, 5.0/100*i);
    f = pb->phasespace_dens(r, p, Pbox, pdg);
    std::cout << 5.0/100*i << "  " << f << std::endl;
  }

}
