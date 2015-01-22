/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <map>
#include <fstream>
#include "unittest.h"

#include "../include/collidermodus.h"
#include "../include/configuration.h"
#include "../include/experiment.h"
#include "../include/modusdefault.h"
#include "../include/nucleus.h"
#include "../include/potentials.h"
#include "../include/spheremodus.h"

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

// Create nuclear potential profile in XY plane
TEST(nucleus_potential_profile) {
  // Create a nucleus
  std::map<PdgCode, int> nuc_list = {{0x2212, 79}, {0x2112, 118}};
  const int Ntest = 1;
  const double sigma = 1.0;
  const double dt = 0.1;

  // Create a nucleus
  Configuration conf(TEST_CONFIG_PATH);
  // All interactions off
  conf["Collision_Term"]["Decays"] = "False";
  conf["Collision_Term"]["Collisions"] = "False";
  conf["Collision_Term"]["Sigma"] = 0.0;
  // Fixed target: Copper
  conf["Modi"]["Collider"]["Calculation_Frame"] = 3;
  conf["Modi"]["Collider"]["E_Kin"] = 1.23;
  conf.take({"Modi", "Collider", "Sqrtsnn"});
  conf.take({"Modi", "Collider", "Projectile"});
  conf.take({"Modi", "Collider", "Target"});
  conf["Modi"]["Collider"]["Projectile"]["Particles"]["211"] = 1;
  conf["Modi"]["Collider"]["Target"]["Particles"]["2212"] = 29;
  conf["Modi"]["Collider"]["Target"]["Particles"]["2112"] = 34;
  conf["Modi"]["Collider"]["Target"]["Automatic"] = "True";

  ExperimentParameters param{{0.f, dt}, 1.f, Ntest, sigma};
  ColliderModus c(conf["Modi"], param);
  Particles P;
  c.initial_conditions(&P, param);
  OutputsList out;
  ParticleList plist;

  // Create potentials
  conf["Potentials"]["Skyrme"]["Skyrme_A"] = -209.2;
  conf["Potentials"]["Skyrme"]["Skyrme_B"] = 156.4;
  conf["Potentials"]["Skyrme"]["Skyrme_Tau"] = 1.35;
  Potentials* pot = new Potentials(conf["Potentials"], param);

  // Write potential XY map in a vtk output
  ThreeVector r;
  const int nx = 50, ny = 50;
  const double dx = 0.2, dy = 0.2;
  double pot_value;

  std::ofstream a_file;
  for (auto it = 0; it < 20; it++) {
    a_file.open(("Nucleus_U_xy.vtk." + std::to_string(it)).c_str(),
                                                     std::ios::out);
    plist = ParticleList(P.data().begin(), P.data().end());
    a_file << "# vtk DataFile Version 2.0\n" <<
              "potential\n" <<
              "ASCII\n" <<
              "DATASET STRUCTURED_POINTS\n" <<
              "DIMENSIONS " << 2*nx+1 << " " << 2*ny+1 << " 1\n" <<
              "SPACING 1 1 1\n" <<
              "ORIGIN " << -nx << " " << -ny << " 0\n" <<
              "POINT_DATA " << (2*nx+1)*(2*ny+1) << "\n" <<
              "SCALARS potential float 1\n" <<
              "LOOKUP_TABLE default\n";

    a_file << std::setprecision(8);
    a_file << std::fixed;
    for (auto iy = -ny; iy <= ny; iy++) {
      for (auto ix = -nx; ix <= nx; ix++) {
        r = ThreeVector(ix*dx, iy*dy, 8.0);
        pot_value = pot->potential(r, plist);
        a_file << pot_value << " ";
      }
      a_file << "\n";
    }
    a_file.close();
    for (auto i = 0; i < 50; i++) {
      c.propagate(&P, param, out, pot);
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
  class Dummy_Pot: public Potentials{
   public:
    Dummy_Pot(Configuration conf, const ExperimentParameters &param,
              const double U0, const double d) :
       Potentials(conf, param),
       U0_(U0), d_(d) {}
    double potential(const ThreeVector &r,
                    const ParticleList &/*plist*/) const {
      return U0_/(1.0 + std::exp(r.x1()/d_));
    }

    ThreeVector potential_gradient(const ThreeVector &r,
                        const ParticleList &/*plist*/) const {
      const double tmp = std::exp(r.x1()/d_);
      return ThreeVector(- U0_/d_ * tmp / ((1.0 + tmp)*(1.0 + tmp)), 0.0, 0.0);
    }
   private:
    const double U0_, d_;
  };

  // Create spheremodus with arbitrary parameters
  // Do not initialize particles: just artificially put one particle to list
  const int Ntest = 1;
  const double sigma = 1.0;
  const double dt = 0.1;
  const double p_mass = 0.938;
  Configuration conf(TEST_CONFIG_PATH);
  ExperimentParameters param{{0.f, dt}, 1.f, Ntest, sigma};
  SphereModus c(conf["Modi"], param);

  // Create dummy outputs and our test potential
  OutputsList out;
  const double U0 = 0.5;
  const double d = 4.0;
  Dummy_Pot* pot = new Dummy_Pot(conf["Potentials"], param, U0, d);

  // Create one particle
  ParticleData part = create_proton();
  part.set_4momentum(p_mass, ThreeVector(2.0, -1.0, 1.0));
  part.set_4position(FourVector(0.0, -20*d, 0.0, 0.0));
  Particles P;
  COMPARE(P.add_data(part), 0);

  // Propagate, until particle is at x>>d, where d is parameter of potential
  while (P.data(0).position().x1() < 20*d) {
    c.propagate(&P, param, out, pot);
  }
  // Calculate 4-momentum, expected from conservation laws
  const FourVector pm = part.momentum();
  FourVector expected_p = FourVector(
         pm.x0() + U0,
         std::sqrt(pm.x1()*pm.x1() + 2*pm.x0()*U0 + U0*U0),
         pm.x2(),
         pm.x3());

  COMPARE_ABSOLUTE_ERROR(expected_p.x0(), P.data(0).momentum().x0(), 1.e-4)<<
    "Expected energy " << expected_p.x0() <<
    ", obtained " << P.data(0).momentum().x0();
  COMPARE_ABSOLUTE_ERROR(expected_p.x1(), P.data(0).momentum().x1(), 1.e-4)<<
    "Expected px " << expected_p.x1() <<
    ", obtained " << P.data(0).momentum().x1();
  // y and z components did not have to change at all, so check is precise
  COMPARE(expected_p.x2(), P.data(0).momentum().x2());
  COMPARE(expected_p.x3(), P.data(0).momentum().x3());

}
