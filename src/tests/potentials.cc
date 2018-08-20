/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include <fstream>
#include <map>

#include "../include/smash/collidermodus.h"
#include "../include/smash/configuration.h"
#include "../include/smash/cxx14compat.h"
#include "../include/smash/experiment.h"
#include "../include/smash/modusdefault.h"
#include "../include/smash/nucleus.h"
#include "../include/smash/potentials.h"
#include "../include/smash/propagation.h"
#include "../include/smash/spheremodus.h"

#include <boost/filesystem.hpp>

using namespace smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "N+ 0.938 0.0 + 2212\n"
      "N0 0.938 0.0 + 2112\n"
      "π+ 0.138 0.0 - 211\n");
}
/*
static ParticleData create_proton(int id = -1) {
  return ParticleData{ParticleType::find(0x2212), id};
}

static ParticleData create_pion(int id = -1) {
  return ParticleData{ParticleType::find(0x211), id};
}
*/

// Create nuclear potential profile in XY plane
TEST(nucleus_potential_profile) {
  // Create a nucleus
  Configuration conf = Test::configuration();
  // All interactions off
  conf["Collision_Term"]["Decays"] = "False";
  conf["Collision_Term"]["Collisions"] = "False";
  conf["Collision_Term"]["Sigma"] = 0.0;
  // Fixed target: Copper
  conf["Modi"]["Collider"]["Calculation_Frame"] = "fixed target";
  conf["Modi"]["Collider"]["E_Kin"] = 1.23;
  conf["Modi"]["Collider"]["Projectile"]["Particles"]["211"] = 1;
  conf["Modi"]["Collider"]["Target"]["Particles"]["2212"] = 29;
  conf["Modi"]["Collider"]["Target"]["Particles"]["2112"] = 34;
  conf["Modi"]["Collider"]["Target"]["Automatic"] = "True";

  ExperimentParameters param = smash::Test::default_parameters();
  ColliderModus c(conf["Modi"], param);
  Particles P;
  c.initial_conditions(&P, param);
  ParticleList plist;

  // Create potentials
  conf["Potentials"]["Skyrme"]["Skyrme_A"] = -209.2;
  conf["Potentials"]["Skyrme"]["Skyrme_B"] = 156.4;
  conf["Potentials"]["Skyrme"]["Skyrme_Tau"] = 1.35;
  std::unique_ptr<Potentials> pot =
      make_unique<Potentials>(conf["Potentials"], param);

  // Write potential XY map in a vtk output
  ThreeVector r;
  const int nx = 50, ny = 50;
  const double dx = 0.2, dy = 0.2;
  double pot_value;
  const ParticleType &proton = ParticleType::find(0x2212);

  std::ofstream a_file;
  const double timestep = param.labclock.timestep_duration();
  for (auto it = 0; it < 20; it++) {
    {
      a_file.open(("Nucleus_U_xy.vtk." + std::to_string(it)).c_str(),
                  std::ios::out);
      plist = P.copy_to_vector();
      a_file << "# vtk DataFile Version 2.0\n"
             << "potential\n"
             << "ASCII\n"
             << "DATASET STRUCTURED_POINTS\n"
             << "DIMENSIONS " << 2 * nx + 1 << " " << 2 * ny + 1 << " 1\n"
             << "SPACING 1 1 1\n"
             << "ORIGIN " << -nx << " " << -ny << " 0\n"
             << "POINT_DATA " << (2 * nx + 1) * (2 * ny + 1) << "\n"
             << "SCALARS potential double 1\n"
             << "LOOKUP_TABLE default\n";

      a_file << std::setprecision(8);
      a_file << std::fixed;
      for (auto iy = -ny; iy <= ny; iy++) {
        for (auto ix = -nx; ix <= nx; ix++) {
          r = ThreeVector(ix * dx, iy * dy, 8.0);
          pot_value = pot->potential(r, plist, proton);
          a_file << pot_value << " ";
        }
        a_file << "\n";
      }
    }

    for (auto i = 0; i < 50; i++) {
      const double time_to = 5.0 * it + i * timestep;
      const double dt = propagate_straight_line(&P, time_to, {});
      update_momenta(&P, dt, *pot, nullptr, nullptr);
    }
  }
}
