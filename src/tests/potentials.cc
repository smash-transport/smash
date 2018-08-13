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
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "N+ 0.938 0.0 2212\n"
      "N0 0.938 0.0 2112\n"
      "π+  0.138 0.0 211\n");
}

static ParticleData create_proton(int id = -1) {
  return ParticleData{ParticleType::find(0x2212), id};
}

static ParticleData create_pion(int id = -1) {
  return ParticleData{ParticleType::find(0x211), id};
}

// check that analytical and numerical results for gradient of potential
// coincide
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

  ThreeVector r, dr;
  Configuration conf = Test::configuration();
  conf["Potentials"]["Skyrme"]["Skyrme_A"] = -209.2;
  conf["Potentials"]["Skyrme"]["Skyrme_B"] = 156.4;
  conf["Potentials"]["Skyrme"]["Skyrme_Tau"] = 1.35;
  ExperimentParameters param = smash::Test::default_parameters();
  std::unique_ptr<Potentials> pot =
      make_unique<Potentials>(conf["Potentials"], param);

  ThreeVector num_grad, analit_grad;
  r = ThreeVector(0.2, 0.0, 0.0);

  // analytical gradient
  const ParticleType &proton = ParticleType::find(0x2212);
  auto tmp = pot->potential_gradient(r, P);
  std::cout << "Grads:" << tmp.first << " " << tmp.second << std::endl;
  analit_grad = tmp.first + tmp.second;
  // numerical gradient
  const double U = pot->potential(r, P, proton);
  dr = ThreeVector(1.e-4, 0.0, 0.0);
  num_grad.set_x1((pot->potential(r + dr, P, proton) - U) / dr.x1());
  dr = ThreeVector(0.0, 1.e-4, 0.0);
  num_grad.set_x2((pot->potential(r + dr, P, proton) - U) / dr.x2());
  dr = ThreeVector(0.0, 0.0, 1.e-4);
  num_grad.set_x3((pot->potential(r + dr, P, proton) - U) / dr.x3());
  // compare them with: accuracy should not be worse than |dr|
  std::cout << num_grad << analit_grad << std::endl;
  COMPARE_ABSOLUTE_ERROR(num_grad.x1(), analit_grad.x1(), 1.e-4);
  COMPARE_ABSOLUTE_ERROR(num_grad.x2(), analit_grad.x2(), 1.e-4);
  COMPARE_ABSOLUTE_ERROR(num_grad.x3(), analit_grad.x3(), 1.e-4);
}

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

TEST(propagation_in_test_potential) {
  /* A dummy potential is created: U(x) = U_0/(1 + exp(x/d))
     A particle is propagated through this potential and
     it's momentum and energy are checked against analytically expected
     from conservation laws.
   */

  // Create a dummy potential
  class Dummy_Pot : public Potentials {
   public:
    Dummy_Pot(Configuration conf, const ExperimentParameters &param,
              const double U0, const double d)
        : Potentials(conf, param), U0_(U0), d_(d) {}
    double potential(const ThreeVector &r, const ParticleList & /*plist*/,
                     const ParticleType & /*acts_on*/) const override {
      return U0_ / (1.0 + std::exp(r.x1() / d_));
    }

    std::pair<ThreeVector, ThreeVector> potential_gradient(
        const ThreeVector &r, const ParticleList &) const override {
      const double tmp = std::exp(r.x1() / d_);
      return std::make_pair(
          ThreeVector(-U0_ / d_ * tmp / ((1.0 + tmp) * (1.0 + tmp)), 0.0, 0.0),
          ThreeVector(0.0, 0.0, 0.0));
    }

    bool use_skyrme() const override { return true; }
    bool use_symmetry() const override { return true; }

   private:
    const double U0_, d_;
  };

  // Create spheremodus with arbitrary parameters
  // Do not initialize particles: just artificially put one particle to list
  const double p_mass = 0.938;
  Configuration conf = Test::configuration();
  ExperimentParameters param = smash::Test::default_parameters();

  // Create dummy outputs and our test potential
  const double U0 = 0.5;
  const double d = 4.0;
  std::unique_ptr<Dummy_Pot> pot =
      make_unique<Dummy_Pot>(conf["Potentials"], param, U0, d);

  // Create one particle
  ParticleData part = create_proton();
  part.set_4momentum(p_mass, ThreeVector(2.0, -1.0, 1.0));
  part.set_4position(FourVector(0.0, -20 * d, 0.0, 0.0));
  Particles P;
  P.insert(part);
  COMPARE(P.back().id(), 0);

  // Propagate, until particle is at x>>d, where d is parameter of potential
  const double timestep = param.labclock.timestep_duration();
  double time_to = 0.0;
  while (P.front().position().x1() < 20 * d) {
    time_to += timestep;
    const double dt = propagate_straight_line(&P, time_to, {});
    update_momenta(&P, dt, *pot, nullptr, nullptr);
  }
  // Calculate 4-momentum, expected from conservation laws
  const FourVector pm = part.momentum();
  FourVector expected_p = FourVector(
      pm.x0() + U0, std::sqrt(pm.x1() * pm.x1() + 2 * pm.x0() * U0 + U0 * U0),
      pm.x2(), pm.x3());

  COMPARE_ABSOLUTE_ERROR(expected_p.x0(), P.front().momentum().x0(), 1.e-4)
      << "Expected energy " << expected_p.x0() << ", obtained "
      << P.front().momentum().x0();
  COMPARE_ABSOLUTE_ERROR(expected_p.x1(), P.front().momentum().x1(), 1.e-4)
      << "Expected px " << expected_p.x1() << ", obtained "
      << P.front().momentum().x1();
  // y and z components did not have to change at all, so check is precise
  COMPARE(expected_p.x2(), P.front().momentum().x2());
  COMPARE(expected_p.x3(), P.front().momentum().x3());
}
