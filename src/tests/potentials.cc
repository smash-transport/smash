/*
 *
 *    Copyright (c) 2014-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/potentials.h"

#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <vector>

#include "setup.h"
#include "smash/algorithms.h"
#include "smash/collidermodus.h"
#include "smash/configuration.h"
#include "smash/constants.h"
#include "smash/density.h"
#include "smash/experiment.h"
#include "smash/modusdefault.h"
#include "smash/nucleus.h"
#include "smash/propagation.h"
#include "smash/quantumsampling.h"
#include "smash/spheremodus.h"

using namespace smash;

TEST(set_random_seed) {
  std::random_device rd;
  int64_t seed = rd();
  random::set_seed(seed);
}

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "π 0.138 7.7e-9 - 111 211\n"
      "N 0.938 0 + 2112 2212\n"
      "Λ 1.116 0 + 3122\n"
      "Σ 1.189 0 + 3112 3212 3222\n"
      "Ξ 1.318 0 + 3312 3322\n"
      "Ω⁻ 1.672 0 + 3334\n");
}

static void compare_pair(const std::pair<double, int>& a,
                         const std::pair<double, int>& b) {
  COMPARE(a.first, b.first);
  COMPARE(a.second, b.second);
}

TEST(force_scale) {
  const auto& p = ParticleType::find(0x2212);
  const auto& antip = ParticleType::find(-0x2212);
  const auto& n = ParticleType::find(0x2112);
  const auto& lambda = ParticleType::find(0x3122);
  const auto& antilambda = ParticleType::find(-0x3122);
  const auto& sigma = ParticleType::find(0x3112);
  const auto& xi = ParticleType::find(0x3312);
  const auto& omega = ParticleType::find(0x3334);
  compare_pair(Potentials::force_scale(p), std::make_pair(1., 1));
  compare_pair(Potentials::force_scale(antip), std::make_pair(-1., -1));
  compare_pair(Potentials::force_scale(n), std::make_pair(1., 1));
  compare_pair(Potentials::force_scale(lambda), std::make_pair(2. / 3., 1));
  compare_pair(Potentials::force_scale(antilambda),
               std::make_pair(-2. / 3., -1));
  compare_pair(Potentials::force_scale(sigma), std::make_pair(2. / 3., 1));
  compare_pair(Potentials::force_scale(xi), std::make_pair(1. / 3., 1));
  compare_pair(Potentials::force_scale(omega), std::make_pair(0., 1));
}

static ParticleData create_proton(int id = -1) {
  return ParticleData{ParticleType::find(0x2212), id};
}

// Create nuclear potential profile in XY plane
TEST(nucleus_potential_profile) {
  // Create a nucleus
  Configuration conf{R"(
    Modi:
      Collider:
        Calculation_Frame: "fixed target"
        E_Kin: 1.23
        Projectile:
          Particles:
            211: 1
        Target:
          Particles:
            2212: 29
            2112: 34
    Potentials:
      Skyrme:
        Skyrme_A: -209.2
        Skyrme_B: 156.4
        Skyrme_Tau: 1.35
  )"};
  conf.validate();
  ExperimentParameters param = smash::Test::default_parameters();
  ColliderModus c(conf.extract_complete_sub_configuration(InputSections::modi),
                  param);
  std::vector<Particles> P(1);
  c.initial_conditions(&(P[0]), param);
  ParticleList plist;
  Potentials pot = Potentials(
      conf.extract_complete_sub_configuration(InputSections::potentials),
      param);

  // Write potential XY map in a vtk output
  ThreeVector r;
  const int nx = 50, ny = 50;
  const double dx = 0.2, dy = 0.2;
  double pot_value;
  const ParticleType& proton = ParticleType::find(0x2212);

  std::ofstream a_file;
  const double timestep = param.labclock->timestep_duration();
  for (auto it = 0; it < 20; it++) {
    {
      a_file.open(("Nucleus_U_xy.vtk." + std::to_string(it)).c_str(),
                  std::ios::out);
      plist = P[0].copy_to_vector();
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
          pot_value = pot.potential(r, plist, proton);
          a_file << pot_value << " ";
        }
        a_file << "\n";
      }
    }

    for (auto i = 0; i < 50; i++) {
      const double time_to = 5.0 * it + i * timestep;
      const double dt = propagate_straight_line(&(P[0]), time_to, {});
      update_momenta(P, dt, pot, nullptr, nullptr, nullptr, nullptr);
    }
  }
}

TEST(propagation_in_test_potential) {
  /* Two dummy potentials are created:
   * One has only the time component: U(x) = U_0/(1 + exp(x/d))
   * A particle is propagated through this stationary potential and
   * its momentum and energy are checked against analytically expected
   * from conservation laws.
   * The other gives rise to a constant magnetic field along z-axis.
   * A particle is expected to do a circular motion with a constant
   * speed in the magnetic field. Its final velocity after one period is
   * compared with the initial one.*/

  // Create a dummy potential
  class Dummy_Pot : public Potentials {
   public:
    Dummy_Pot(const ExperimentParameters& param, const double U0,
              const double d, const double B0)
        : Potentials(Configuration{""}, param), U0_(U0), d_(d), B0_(B0) {}

    std::tuple<ThreeVector, ThreeVector, ThreeVector, ThreeVector> all_forces(
        const ThreeVector& r, const ParticleList&) const override {
      const double tmp = std::exp(r.x1() / d_);
      return std::make_tuple(
          ThreeVector(U0_ / d_ * tmp / ((1.0 + tmp) * (1.0 + tmp)), 0.0, 0.0),
          ThreeVector(0., 0., B0_), ThreeVector(), ThreeVector());
    }

    bool use_skyrme() const override { return true; }
    bool use_symmetry() const override { return true; }

   private:
    const double U0_, d_, B0_;
  };

  // Create spheremodus with arbitrary parameters
  // Do not initialize particles: just artificially put one particle to list
  const double p_mass = 0.938;
  ExperimentParameters param = smash::Test::default_parameters();

  /* Create two dummy test potentials: one has only the electrical
   * force, the other has only magnetic force. */
  const double U0 = 0.5;
  const double d = 4.0;
  const double B0 = 0.5;
  std::unique_ptr<Dummy_Pot> pot1 =
      std::make_unique<Dummy_Pot>(param, U0, d, 0.);
  std::unique_ptr<Dummy_Pot> pot2 =
      std::make_unique<Dummy_Pot>(param, 0., d, B0);

  /* Create two particles: one flies in the pure electrical field, while
   * the other flies in the pure magnetical field. */
  ParticleData part1 = create_proton();
  part1.set_4momentum(p_mass, 2.0, -1.0, 1.0);
  part1.set_4position(FourVector(0.0, -20 * d, 0.0, 0.0));
  std::vector<Particles> P1(1);
  P1[0].insert(part1);
  COMPARE(P1[0].back().id(), 0);

  // This particle is expected to do circular motion with a constant speed.
  ParticleData part2 = create_proton();
  // The speed of the particle is 0.6
  part2.set_4momentum(4., 3., 0., 0.);
  part2.set_4position(FourVector(0.0, 0.0, 2.0, 0.0));
  std::vector<Particles> P2(1);
  P2[0].insert(part2);
  COMPARE(P2[0].back().id(), 0);

  /* Propagate the first particle, until particle1 is at x>>d,
   * where d is parameter of potential */
  const double timestep = 0.01;
  double time_to = 0.0;
  while (P1[0].front().position().x1() < 20 * d) {
    time_to += timestep;
    const double dt1 = propagate_straight_line(&(P1[0]), time_to, {});
    update_momenta(P1, dt1, *pot1, nullptr, nullptr, nullptr, nullptr);
  }

  // Propagate the second particle for one period.
  const double period = twopi * 5. / B0;
  time_to = 0.0;
  while (P2[0].front().position().x0() < period) {
    time_to += timestep;
    const double dt2 = propagate_straight_line(&(P2[0]), time_to, {});
    update_momenta(P2, dt2, *pot2, nullptr, nullptr, nullptr, nullptr);
  }

  // Calculate 4-momentum, expected from conservation laws
  const FourVector pm = part1.momentum();
  FourVector expected_p = FourVector(
      pm.x0() + U0, std::sqrt(pm.x1() * pm.x1() + 2 * pm.x0() * U0 + U0 * U0),
      pm.x2(), pm.x3());

  COMPARE_ABSOLUTE_ERROR(expected_p.x0(), P1[0].front().momentum().x0(), 1.e-6)
      << "Expected energy " << expected_p.x0() << ", obtained "
      << P1[0].front().momentum().x0();
  COMPARE_ABSOLUTE_ERROR(expected_p.x1(), P1[0].front().momentum().x1(), 1.e-6)
      << "Expected px " << expected_p.x1() << ", obtained "
      << P1[0].front().momentum().x1();
  // y and z components did not have to change at all, so check is precise
  COMPARE(expected_p.x2(), P1[0].front().momentum().x2());
  COMPARE(expected_p.x3(), P1[0].front().momentum().x3());

  /* compare the final velocity of the second particle, component by component,
   * with the initial values,
   * the test is passed if the errors are within 0.003, which is 0.5 percent of
   * the initial speed. */

  COMPARE_ABSOLUTE_ERROR(P2[0].front().momentum().velocity().x1(), 0.6, 0.003)
      << "Expected x-component of velocity " << 0.6 << ", obtained "
      << P2[0].front().momentum().velocity().x1();
  COMPARE_ABSOLUTE_ERROR(P2[0].front().momentum().velocity().x2(), 0.0, 0.003)
      << "Expected y-component of velocity " << 0 << ", obtained "
      << P2[0].front().momentum().velocity().x2();
  COMPARE_ABSOLUTE_ERROR(P2[0].front().momentum().velocity().x3(), 0.0, 0.003)
      << "Expected z-component of velocity " << 0 << ", obtained "
      << P2[0].front().momentum().velocity().x3();
}

/*
 * The idea is to compute potentials from the same set of particles,
 * but in one case they are testparticles in one ensemble, while in the
 * other case they are particles in many ensembles. The density from both
 * calculations should be the same, so the potentials should be identical too.
 *
 */
TEST(ensembles_vs_testparticles) {
  auto random_value = random::make_uniform_distribution(-2.0, +2.0);
  int Ntest = 17;
  int N = 39;
  ParticleList plist;

  for (int id = 0; id < N * Ntest; id++) {
    ParticleData p{ParticleType::find(0x2212), id};
    p.set_4position(
        {random_value(), random_value(), random_value(), random_value()});
    p.set_4momentum(smash::nucleon_mass,
                    {random_value(), random_value(), random_value()});
    plist.push_back(p);
  }

  const char* conf_pot{R"(
    Potentials:
      Skyrme:
          Skyrme_A: -209.2
          Skyrme_B: 156.4
          Skyrme_Tau: 1.35
      Symmetry:
          S_Pot: 18.0
  )"};
  Configuration conf1{conf_pot}, conf2{conf_pot};
  ExperimentParameters param1 = smash::Test::default_parameters(),
                       param2 = smash::Test::default_parameters();
  param1.testparticles = Ntest;
  param1.n_ensembles = 1;
  param2.testparticles = 1;
  param2.n_ensembles = Ntest;
  Potentials pot1(std::move(conf1), param1), pot2(std::move(conf2), param2);

  const ThreeVector r = ThreeVector(0., 0., 0.);
  const auto forces1 = pot1.all_forces(r, plist),
             forces2 = pot2.all_forces(r, plist);
  ThreeVector a, b;
  a = std::get<0>(forces1);
  b = std::get<0>(forces2);
  VERIFY(a == b) << a << " " << b;
  a = std::get<1>(forces1);
  b = std::get<1>(forces2);
  VERIFY(a == b) << a << " " << b;
  a = std::get<2>(forces1);
  b = std::get<2>(forces2);
  VERIFY(a == b) << a << " " << b;
  a = std::get<3>(forces1);
  b = std::get<3>(forces2);
  VERIFY(a == b) << a << " " << b;
}

/*
 * Compare the calculation of the energy gradient in the calculation frame using
 * single_particle_energy_gradient to the gradient of the Skyrme potential
 * calculated using the chain rule. They should coincide when no boosts are
 * involved, so when the local rest-frame and the calculation frame coincide.
 */
TEST(energy_gradient_vs_pot_gradient) {
  double length = 5.0;  // consider a 5 fm box
  // create some protons but make sure they are at rest
  std::vector<Particles> ensembles(1);
  Particles particles;
  auto uniform_length = random::make_uniform_distribution(0.0, length);
  auto uniform_zero_to_one = random::make_uniform_distribution(0.0, 1.0);
  const int nparticles = 500;
  int isampled = 0;
  for (int itry = 0; (isampled < nparticles) && (itry < 1000 * nparticles);
       itry++) {
    const ThreeVector position_sample = {uniform_length(), uniform_length(),
                                         uniform_length()};
    // sample particles according to P ~ z/length to create a finite gradient
    const double weight = position_sample.x3() / length;
    if (weight < uniform_zero_to_one()) {
      continue;
    }
    ParticleData p(ParticleType::find(0x2212), itry);
    p.set_4momentum(0.938, 0.0, 0.0, 0.0);
    p.set_3position(position_sample);
    ensembles[0].insert(p);
    isampled++;
  }
  if (isampled < nparticles) {
    throw std::runtime_error("Failed to sample enough particles, sampled " +
                             std::to_string(isampled) + " instead of " +
                             std::to_string(nparticles));
  }
  // set up lattice
  std::array<double, 3> lengths = {length, length, length};
  std::array<int, 3> ncells = {
      11, 11, 51};  // lattice needs to be fine because we compare different
                    // types of derivatives in the following
  std::array<double, 3> origin = {0.0, 0.0, 0.0};
  std::unique_ptr<RectangularLattice<DensityOnLattice>> lat =
      std::make_unique<RectangularLattice<DensityOnLattice>>(
          lengths, ncells, origin, false, LatticeUpdate::EveryTimestep);
  // update lattice to calculate densities
  ExperimentParameters exp_par = Test::default_parameters();
  exp_par.testparticles = 50;
  DensityParameters denspar(exp_par);
  update_lattice_accumulating_ensembles(lat.get(), LatticeUpdate::EveryTimestep,
                                        DensityType::Baryon, denspar, ensembles,
                                        true);
  DensityOnLattice jB = (*lat)[lat->index1d(
      ncells[0] / 2, ncells[1] / 2,
      ncells[2] / 2)];  // density in the center (int division intended)
  // set up a potentials object without momentum dependence
  const char* conf_pot{R"(
    Potentials:
      Skyrme:
          Skyrme_A: -209.2
          Skyrme_B: 156.4
          Skyrme_Tau: 1.35
      Momentum_Dependence:
          C: 0.0
          Lambda: 2.13
  )"};
  Configuration conf{conf_pot};
  Potentials pot(std::move(conf), denspar);

  // calculate needed densities
  double baryon_density = jB.rho();
  ThreeVector baryon_grad_j0 = jB.grad_j0();
  ThreeVector baryon_dvecj_dt = jB.dvecj_dt();
  ThreeVector baryon_curl_vecj = jB.curl_vecj();
  ThreeVector force_chain_rule =
      pot.skyrme_force(baryon_density, baryon_grad_j0, baryon_dvecj_dt,
                       baryon_curl_vecj)
          .first;
  // calculate Force via the gradient of the single particle energy
  ParticleList dummy_plist;
  ThreeVector energy_grad = pot.single_particle_energy_gradient(
      lat.get(), {2.5, 2.5, 2.5}, {0.0, 0.0, 0.0}, 0.938, dummy_plist);
  /* compare the two with relatively large tolerance as the gradients are
   * calculated in a different way, energy_grad uses finite difference but the
   * density gradient is done with Gaussian derivatives. They would converge for
   * a very fine lattice. */
  COMPARE_ABSOLUTE_ERROR((-energy_grad - force_chain_rule).abs(), 0.0, 0.001);
  COMPARE_RELATIVE_ERROR(-energy_grad[2], force_chain_rule[2], 0.001);

  // do the same calculation without using the lattice for the energy gradient
  ParticleList real_plist = ensembles[0].copy_to_vector();
  energy_grad = pot.single_particle_energy_gradient(
      nullptr, {2.5, 2.5, 2.5}, {0.0, 0.0, 0.0}, 0.938, real_plist);
  COMPARE_ABSOLUTE_ERROR((-energy_grad - force_chain_rule).abs(), 0.0, 0.001);
  COMPARE_RELATIVE_ERROR(-energy_grad[2], force_chain_rule[2], 0.001);
};

// create experiment parameters for tests with VDF
static ExperimentParameters default_parameters_vdf(
    int testparticles = 1, double dt = 0.1,
    double triangular_smearing_range = 2.0) {
  return ExperimentParameters{
      std::make_unique<UniformClock>(0., dt, 300.0),  // labclock
      std::make_unique<UniformClock>(0., 1., 300.0),  // outputclock
      1,                                              // ensembles
      testparticles,                                  // testparticles
      DerivativesMode::FiniteDifference,              // derivatives mode
      // both the rest frame and the direct derivatives need to be on for the
      // test of forces calculated using chain rule and direct derivatives
      RestFrameDensityDerivativesMode::On,  // rest frame derivatives mode
      FieldDerivativesMode::Direct,         // field derivatives mode
      SmearingMode::Triangular,             // smearing modeing mode
      1.0,                                  // Gaussian smearing width
      4.0,                                  // Gaussian smearing cut-off
      0.333333,                             // discrete smearing weight
      triangular_smearing_range,            // triangular smearing range
      CollisionCriterion::Geometric,
      false,  // two_to_one
      false, Test::no_multiparticle_reactions(),
      false,  // strings switch
      1.0, NNbarTreatment::NoAnnihilation,
      0.,           // low energy sigma_NN cut-off
      false,        // potential_affect_threshold
      -1.0,         // box_length
      200.0,        // max. cross section
      2.5,          // fixed min. cell length
      1.0,          // cross section scaling
      false,        // in thermodynamics outputs spectators are included
      false,        // do non-strong decays
      true,         // can decay initial particles
      std::nullopt  // use monash tune, not known
  };
}

/*
 * Testing the values of the vdf forces for two ways of calculating gradients of
 * the vector field: with chain rule derivatives and with field derivatives. The
 * test initializes matter in a cubic space of size L, with density distributed
 * uniformly in the x- and y-directions, and distributed linearly in the
 * z-direction. This results in a constant density gradient in the z-direction
 * (up to corrections due to fluctuations). For simplicity, we initialize the
 * matter at zero temperature (sampling the momenta from the Fermi sphere). The
 * initialization resembles an initialization of a Box Modus generalized to a
 * case where there is a constant gradient of density, however, there are no
 * boundary conditions applied. For this reason the calculation of density on
 * the edge of the space will be faulty, and it is the safest to check the
 * values of the forces in the middle of the cubic space.
 */
TEST(vdf_chain_rule_derivatives_vs_vdf_direct_derivatives) {
  // a large Ntest is necessary for high precision
  const int Ntest = 1000;
  // initialize the experiment parameters and potentials
  Configuration conf{R"(
    Potentials:
      VDF:
        Sat_rhoB: 0.168
        Powers: [2.0, 2.35]
        Coeffs: [-209.2, 156.5]
  )"};
  ExperimentParameters param = default_parameters_vdf(Ntest);
  Potentials pot(std::move(conf), param);

  // the side length of the cubic space
  const double length = 10.0;

  // the gradient in density rises linearly from some minimal to some maximal
  // value; here given in units of the saturation density
  const double rho_min = 0.25;
  const double rho_max = 0.75;
  const double rho_0 = pot.saturation_density();
  // the average density in the cubic space, in fm^{-3}
  const double rho_avg = rho_0 * (rho_max + rho_min) / 2.0;

  // calculate how many particles of one type are needed in a cubic space of
  // given length to reproduce this density profile
  const int N = numeric_cast<int>(rho_avg * length * length * length);

  // generate numbers uniformly on the intervals (0,1) and (0, length)
  auto uniform_one = random::make_uniform_distribution(0.0, 1.0);
  auto uniform_length = random::make_uniform_distribution(0.0, length);

  ParticleList plist;
  std::vector<Particles> P(1);
  // orientation of the momentum vector
  Angles phitheta;
  for (int id = 0; id < N * Ntest; id++) {
    ParticleData p{ParticleType::find(0x2212), id};

    // the density distribution is uniform in the x- and y-directions, but in
    // the z-direction we want the density to start at n_min for z=0 and rise
    // linearly to n_max for z=length; we use rejection sampling to sample the
    // z-position
    double sampled_z = 0.0;
    bool success = false;
    while (!success) {
      sampled_z = uniform_length();
      const double distribution_at_sampled_z =
          rho_min + (rho_max - rho_min) * sampled_z / length;
      const double sampled_ratio = distribution_at_sampled_z / rho_max;
      const double accept_or_reject = uniform_one();
      if (sampled_ratio > accept_or_reject) {
        success = true;
      }
    }
    const ThreeVector pos{uniform_length(), uniform_length(), sampled_z};
    p.set_4position({0, pos});

    // particle's momentum needs to be sampled from the Fermi sphere, whose
    // radius pF varies with the z-position of the particle; the density needs
    // to be expressed in GeV to yield a momentum in GeV
    double density_at_z = rho_min + (rho_max - rho_min) * sampled_z / length;
    density_at_z = rho_0 * density_at_z * hbarc * hbarc * hbarc;
    // spin degeneracy: we treat the protons that we initialized (2212) as
    // standing in for nucleons (protons and neutrons); this is so that the
    // behavior of the initialized system can be easily compared to that of
    // uniform nuclear matter at the given density
    const double spin_degeneracy = 4.0;
    const double pF = cbrt(6.0 * M_PI * M_PI * density_at_z / spin_degeneracy);

    double sample = uniform_one();
    double radial_momentum = cbrt(sample) * pF;
    double mass = p.type().mass();
    phitheta.distribute_isotropically();
    p.set_4momentum(mass, phitheta.threevec() * radial_momentum);

    plist.push_back(p);
    P[0].insert(p);
  }

  // parameters of the lattices; 1 lattice point per 1 fm^3
  const std::array<double, 3> l = {10, 10, 10};
  const std::array<int, 3> n = {10, 10, 10};
  const std::array<double, 3> origin = {0.0, 0.0, 0.0};
  const bool periodic = false;
  const DensityParameters par(param);

  // lattices for calculation of the VDF force using the chain rule derivatives
  std::unique_ptr<RectangularLattice<FourVector>> old_jmu_aux =
      std::make_unique<RectangularLattice<FourVector>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);
  std::unique_ptr<RectangularLattice<FourVector>> new_jmu_aux =
      std::make_unique<RectangularLattice<FourVector>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);
  std::unique_ptr<RectangularLattice<std::array<FourVector, 4>>>
      four_gradient_aux =
          std::make_unique<RectangularLattice<std::array<FourVector, 4>>>(
              l, n, origin, periodic, LatticeUpdate::EveryTimestep);

  std::unique_ptr<DensityLattice> jmu_B_lattice =
      std::make_unique<DensityLattice>(l, n, origin, periodic,
                                       LatticeUpdate::EveryTimestep);

  std::unique_ptr<RectangularLattice<FourVector>> UB_lattice_density =
      std::make_unique<RectangularLattice<FourVector>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);
  std::unique_ptr<RectangularLattice<std::pair<ThreeVector, ThreeVector>>>
      FB_lattice_density = std::make_unique<
          RectangularLattice<std::pair<ThreeVector, ThreeVector>>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);

  // lattices for calculation of the VDF force using the direct field
  // derivatives
  std::unique_ptr<RectangularLattice<FourVector>> old_fields_aux =
      std::make_unique<RectangularLattice<FourVector>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);
  std::unique_ptr<RectangularLattice<FourVector>> new_fields_aux =
      std::make_unique<RectangularLattice<FourVector>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);
  std::unique_ptr<RectangularLattice<std::array<FourVector, 4>>>
      fields_four_gradient_aux =
          std::make_unique<RectangularLattice<std::array<FourVector, 4>>>(
              l, n, origin, periodic, LatticeUpdate::EveryTimestep);
  std::unique_ptr<FieldsLattice> fields_lattice =
      std::make_unique<FieldsLattice>(l, n, origin, periodic,
                                      LatticeUpdate::EveryTimestep);

  std::unique_ptr<RectangularLattice<FourVector>> UB_lattice_field =
      std::make_unique<RectangularLattice<FourVector>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);
  std::unique_ptr<RectangularLattice<std::pair<ThreeVector, ThreeVector>>>
      FB_lattice_field = std::make_unique<
          RectangularLattice<std::pair<ThreeVector, ThreeVector>>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);

  update_lattice(jmu_B_lattice.get(), old_jmu_aux.get(), new_jmu_aux.get(),
                 four_gradient_aux.get(), LatticeUpdate::EveryTimestep,
                 DensityType::Baryon, par, P,
                 param.labclock->timestep_duration(), true);

  update_fields_lattice(fields_lattice.get(), old_fields_aux.get(),
                        new_fields_aux.get(), fields_four_gradient_aux.get(),
                        jmu_B_lattice.get(), LatticeUpdate::EveryTimestep, pot,
                        param.labclock->timestep_duration());

  const size_t lattice_size = UB_lattice_density->size();

  for (auto& node : *jmu_B_lattice) {
    node.overwrite_drho_dt_to_zero();
    node.overwrite_djmu_dt_to_zero();
  }

  for (auto& node : *fields_lattice) {
    node.overwrite_dAmu_dt_to_zero();
  }

  // we will calculate the average of the relative difference between the two
  // ways of calculating the forces over a number of nodes (samples)
  double relative_diff = 0.0;
  int samples = 0;

  for (size_t i = 0; i < lattice_size; i++) {
    auto jB = (*jmu_B_lattice)[i];
    auto Amu = (*fields_lattice)[i];
    ThreeVector cell_center = (*FB_lattice_density).cell_center(i);

    // get force using chain rule derivatives
    (*FB_lattice_density)[i] =
        pot.vdf_force(jB.rho(), jB.drho_dxnu().x0(), jB.drho_dxnu().threevec(),
                      jB.grad_rho_cross_vecj(), jB.jmu_net().x0(), jB.grad_j0(),
                      jB.jmu_net().threevec(), jB.dvecj_dt(), jB.curl_vecj());

    // get force using direct derivatives
    (*FB_lattice_field)[i] =
        pot.vdf_force(Amu.grad_A0(), Amu.dvecA_dt(), Amu.curl_vecA());

    // calculate average relative difference error over all cells in the
    // relative middle of the cubic space
    if ((cell_center.x1() > 3.0) && (cell_center.x1() < 7.0)) {
      if ((cell_center.x2() > 3.0) && (cell_center.x1() < 7.0)) {
        if ((cell_center.x3() > 3.0) && (cell_center.x3() < 7.0)) {
          double aux = std::abs((*FB_lattice_density)[i].first.x3() -
                                (*FB_lattice_field)[i].first.x3());
          aux = aux / (std::abs((*FB_lattice_density)[i].first.x3() +
                                (*FB_lattice_field)[i].first.x3()) /
                       2.0);

          relative_diff += aux;
          samples++;
        }
      }
    }
  }

  relative_diff = relative_diff / samples;

  VERIFY(relative_diff < 0.025)
      << "The relative difference between forces calculated with chain rule "
         "derivatives and with direct derivatives is too large. \nThe "
         "expectation is for the two forces to be comparable.";
}

/*
 * Comparing the computation of Skyrme and VDF potentials in the case where
 * latice is not used; both should give the same force for the VDF potential
 * defined to reproduce the Skyrme potential. Details of the initialization
 * mirror those used in the vdf_chain_rule_derivatives_vs_vdf_direct_derivatives
 * test.
 */
TEST(skyrme_vs_vdf_wo_lattice) {
  // a large Ntest is necessary for high precision
  const int Ntest = 1000;
  // initialize the experiment parameters and potentials
  Configuration conf_pot1{R"(
    Potentials:
      Skyrme:
          Skyrme_A: -209.2
          Skyrme_B: 156.4
          Skyrme_Tau: 1.35
  )"};

  Configuration conf_pot2{R"(
    Potentials:
      VDF:
        Sat_rhoB: 0.168
        Powers: [2.0, 2.35]
        Coeffs: [-209.2, 156.5]
  )"};
  ExperimentParameters param = default_parameters_vdf(Ntest);
  Potentials pot1(std::move(conf_pot1), param),
      pot2(std::move(conf_pot2), param);

  // the side length of the cubic space
  const double length = 10.0;

  // the gradient in density rises linearly from some minimal to some maximal
  // value; here given in units of the saturation density
  const double rho_min = 0.25;
  const double rho_max = 0.75;
  const double rho_0 = pot2.saturation_density();
  // the average density in the cubic space, in fm^{-3}
  const double rho_avg = rho_0 * (rho_max + rho_min) / 2.0;

  // calculate how many particles of one type are needed in a cubic space of
  // given length to reproduce this density profile
  const int N = numeric_cast<int>(rho_avg * length * length * length);

  // generate numbers uniformly on the intervals (0,1) and (0, length)
  auto uniform_one = random::make_uniform_distribution(0.0, 1.0);
  auto uniform_length = random::make_uniform_distribution(0.0, length);

  ParticleList plist;
  std::vector<Particles> P(1);
  // orientation of the momentum vector
  Angles phitheta;
  for (int id = 0; id < N * Ntest; id++) {
    ParticleData p{ParticleType::find(0x2212), id};

    // the density distribution is uniform in the x- and y-directions, but in
    // the z-direction we want the density to start at n_min for z=0 and rise
    // linearly to n_max for z=length; we use rejection sampling to sample the
    // z-position
    double sampled_z = 0.0;
    bool success = false;
    while (!success) {
      sampled_z = uniform_length();
      const double distribution_at_sampled_z =
          rho_min + (rho_max - rho_min) * sampled_z / length;
      const double sampled_ratio = distribution_at_sampled_z / rho_max;
      const double accept_or_reject = uniform_one();
      if (sampled_ratio > accept_or_reject) {
        success = true;
      }
    }
    const ThreeVector pos{uniform_length(), uniform_length(), sampled_z};
    p.set_4position({0, pos});

    // particle's momentum needs to be sampled from the Fermi sphere, whose
    // radius pF varies with the z-position of the particle; the density needs
    // to be expressed in GeV to yield a momentum in GeV
    double density_at_z = rho_min + (rho_max - rho_min) * sampled_z / length;
    density_at_z = rho_0 * density_at_z * hbarc * hbarc * hbarc;
    // spin degeneracy: we treat the protons that we initialized (2212) as
    // standing in for nucleons (protons and neutrons); this is so that the
    // behavior of the initialized system can be easily compared to that of
    // uniform nuclear matter at the given density
    const double spin_degeneracy = 4.0;
    const double pF = cbrt(6.0 * M_PI * M_PI * density_at_z / spin_degeneracy);

    double sample = uniform_one();
    double radial_momentum = cbrt(sample) * pF;
    double mass = p.type().mass();
    phitheta.distribute_isotropically();
    p.set_4momentum(mass, phitheta.threevec() * radial_momentum);

    plist.push_back(p);
    P[0].insert(p);
  }

  // calculate the average of the relative difference between the forces
  // calculated from the two potentials over a number of nodes (samples)
  double relative_diff = 0;
  int samples = 0;
  for (int x = 4; x < 7; x++) {
    for (int y = 4; y < 7; y++) {
      for (int z = 4; z < 7; z++) {
        const ThreeVector r = ThreeVector(x, y, z);
        const auto forces1 = pot1.all_forces(r, plist),
                   forces2 = pot2.all_forces(r, plist);
        double aux =
            std::abs(std::get<0>(forces1).x3() - std::get<0>(forces2).x3());
        aux =
            aux /
            ((std::abs(std::get<0>(forces1).x3() + std::get<0>(forces2).x3())) /
             2.0);

        relative_diff += aux;
        samples++;
      }
    }
  }

  relative_diff = relative_diff / samples;

  VERIFY(relative_diff < 0.025)
      << "The relative difference between forces calculated with the Skyrme "
         "and with the VDF potential is too large. \nThe expectation is for "
         "the Skyrme and VDF forces to be comparable for this initialization.";
}

/*
 * Creates a histogram of densities in cells of a (cubic) box, and checks
 * whether the expected relation between histogram densities was found:
 * returns true if a) there are more cells with density test_left than cells
 * with density mean_density; b) there is at least one cell with density
 * test_right. Such relations are expected when spinodal decomposition occurs.
 */
static bool density_hist(Particles* P, int box_length, int cell_length,
                         int test_p, double mean_density, double test_left,
                         double test_right, double step,
                         double saturation_density) {
  // we split the box into cubic cells; for each of the cells, the local density
  // (scaled by the saturation density) is computed
  VERIFY(box_length % cell_length == 0)
      << "We expect the cell length of the density histogram to be a divisor "
         "of the box length. Instead we got "
      << box_length << " " << cell_length;
  static const int num_cell = int(box_length * box_length * box_length /
                                  (cell_length * cell_length * cell_length));
  // number of cells per a side of the box
  const double factor = box_length / cell_length;

  std::vector<double> cells(num_cell, 0.0);

  int count = 0;
  // go through the particle list
  for (Particles::iterator it = P->begin(); it != P->end(); ++it) {
    count = count + 1;
    // find the cell in which a given particle resides
    double x1 = it->position().x1();
    double x2 = it->position().x2();
    double x3 = it->position().x3();
    const int index = static_cast<int>(
        std::floor(x1 / cell_length) + std::floor(x2 / cell_length) * factor +
        std::floor(x3 / cell_length) * factor * factor);
    // increase the local density of this cell
    cells[index] = cells[index] + 1 / (cell_length * cell_length * cell_length *
                                       saturation_density * test_p);
  }
  double min_bound = DBL_MAX;
  double max_bound = 0;

  // construct a histogram over densities, counting how many cells have density
  // that falls within density range for a given histogram bin
  for (double c : cells) {
    min_bound = std::min(min_bound, c);
    max_bound = std::max(max_bound, c);
  }

  unsigned int hist_bin =
      numeric_cast<int>(std::ceil((max_bound - min_bound) / step));
  min_bound = min_bound / step;
  max_bound = max_bound / step;
  std::vector<unsigned int> histogram(hist_bin + 1, 0);

  for (double c : cells) {
    double position = std::max(0.0, (c / step - min_bound));
    position = std::min(max_bound - min_bound, position);
    histogram[int(std::round(position))]++;
  }
  // get the number of cells which have the two test densities and the mean
  // density
  double position_a = std::max(0.0, (test_left / step - min_bound));
  position_a = std::min(max_bound - min_bound, position_a);
  double position_b = std::max(0.0, (mean_density / step - min_bound));
  position_b = std::min(max_bound - min_bound, position_b);
  double position_c = std::max(0.0, (test_right / step - min_bound));
  position_c = std::min(max_bound - min_bound, position_c);
  bool peak_test = false;
  // we expect a change of the initial distribution in the box for spinodal
  // decomposition: in a box in which the average density is mean_density, many
  // cells are expected to have a lower density test_left, and we ecpect to have
  // at least one cell with a higher density test_right
  if (histogram[int(round(position_b))] <= histogram[int(round(position_a))] and
      histogram[int(round(position_c))] > 0) {
    peak_test = true;
  } else {
    std::cout << " Density histogram " << std::endl;
    for (unsigned int i = 0; i < hist_bin; i++) {
      std::cout << i * step - min_bound << "\t " << histogram[i] << std::endl;
    }
  }
  return peak_test;
}

/*
 * Test for the occurence of a spinodal decomposition in the spinodal region
 * of the liquid/gas phase transition at zero temperature in a cubic box with
 * periodic boundary conditions. There are two parts to the test: The first
 * checks, using a histogram of densities in the box, whether the expected
 * separation of densities has occurred. The second checks whether an
 * expected change in the value of the mean-field energy has taken place.
 */
TEST(spinodal_dilute) {
  // initialise the system with a mean density of 0.25 saturation density;
  // distribute the particles using a Fermi distribution;
  // we use isospin symmetric matter with 21 protons and 21 neutrons which
  // gives 0.25*0.168 density in the 10 fm box
  const int Ntest = 400, N_p = 21, box_length = 10.0, cell_length = 2;
  std::vector<Particles> P(1);

  // The baryon density includes both protons and neutrons, of which we
  // initialise N_p each.
  double density = 2.0 * N_p / (box_length * box_length * box_length);
  // degeneracy of isospin symmetric nuclear matter
  const double nucl_deg = 4.0;
  // calculate Fermi momentum
  const double pF =
      cbrt(6.0 * M_PI * M_PI * density * hbarc * hbarc * hbarc / nucl_deg);

  auto uniform_length =
      random::make_uniform_distribution(0.0, static_cast<double>(box_length));
  double momentum_radial = 0.0, start_time_ = 0, mass = smash::nucleon_mass;
  Angles phitheta;
  double time_to = start_time_;

  auto uniform_one = random::make_uniform_distribution(0.0, 1.0);
  // sample protons
  for (int id = 0; id < N_p * Ntest; id++) {
    ParticleData p{ParticleType::find(0x2212), id};
    momentum_radial = (cbrt(uniform_one()) * pF);
    phitheta.distribute_isotropically();
    p.set_4momentum(mass, phitheta.threevec() * momentum_radial);
    ThreeVector pos{uniform_length(), uniform_length(), uniform_length()};
    p.set_4position(FourVector(start_time_, pos));
    /// Initialize formation time
    p.set_formation_time(start_time_);
    P[0].insert(p);
  }

  // sample neutrons
  for (int id = 0; id < N_p * Ntest; id++) {
    ParticleData p{ParticleType::find(0x2112), id};
    momentum_radial = (cbrt(uniform_one()) * pF);
    phitheta.distribute_isotropically();
    p.set_4momentum(mass, phitheta.threevec() * momentum_radial);
    ThreeVector pos{uniform_length(), uniform_length(), uniform_length()};
    p.set_4position(FourVector(start_time_, pos));
    /// Initialize formation time
    p.set_formation_time(start_time_);
    P[0].insert(p);
  }

  // initialize potential
  Configuration conf{R"(
    Potentials:
      VDF:
        Sat_rhoB: 0.168
        Powers: [2.0, 2.35]
        Coeffs: [-209.2, 156.5]
  )"};
  ExperimentParameters param = default_parameters_vdf(Ntest, 0.1, 2.0);

  Potentials pot(std::move(conf), param);
  const std::array<double, 3> l = {10, 10, 10};
  const std::array<int, 3> n = {10, 10, 10};
  const std::array<double, 3> origin = {0, 0, 0};
  const bool periodic = true;
  const DensityParameters par(param);

  std::unique_ptr<RectangularLattice<FourVector>> old_jmu_auxiliary_df =
      std::make_unique<RectangularLattice<FourVector>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);
  std::unique_ptr<RectangularLattice<FourVector>> new_jmu_auxiliary_df =
      std::make_unique<RectangularLattice<FourVector>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);
  std::unique_ptr<RectangularLattice<std::array<FourVector, 4>>>
      four_gradient_auxiliary_df =
          std::make_unique<RectangularLattice<std::array<FourVector, 4>>>(
              l, n, origin, periodic, LatticeUpdate::EveryTimestep);

  std::unique_ptr<DensityLattice> jmu_B_lat_df =
      std::make_unique<DensityLattice>(l, n, origin, periodic,
                                       LatticeUpdate::EveryTimestep);
  std::unique_ptr<RectangularLattice<FourVector>> UB_lat_df =
      std::make_unique<RectangularLattice<FourVector>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);
  std::unique_ptr<RectangularLattice<std::pair<ThreeVector, ThreeVector>>>
      FB_lat_df = std::make_unique<
          RectangularLattice<std::pair<ThreeVector, ThreeVector>>>(
          l, n, origin, periodic, LatticeUpdate::EveryTimestep);

  // the mean field energy at the beginning and end of the evolution
  double E_init, E_final;
  // time evolution of the system
  for (int i = 0; i < 2000; i++) {
    time_to = time_to + 0.1;

    const double dt = propagate_straight_line(&P[0], time_to, {});
    for (ParticleData& data : P[0]) {
      FourVector position = data.position();
      bool wall_hit = enforce_periodic_boundaries(position.begin() + 1,
                                                  position.end(), box_length);
      if (wall_hit) {
        data.set_4position(position);
      }
    }

    update_lattice(jmu_B_lat_df.get(), old_jmu_auxiliary_df.get(),
                   new_jmu_auxiliary_df.get(), four_gradient_auxiliary_df.get(),
                   LatticeUpdate::EveryTimestep, DensityType::Baryon, par, P,
                   dt, true);
    if (i == 0) {
      for (auto& node : *jmu_B_lat_df) {
        node.overwrite_drho_dt_to_zero();
        node.overwrite_djmu_dt_to_zero();
      }
      // initial mean field energy
      E_init = calculate_mean_field_energy(pot, *jmu_B_lat_df, nullptr, param);
    }

    const size_t UBlattice_size_df = UB_lat_df->size();
    for (size_t j = 0; j < UBlattice_size_df; j++) {
      auto jB_df = (*jmu_B_lat_df)[j];
      (*UB_lat_df)[j] = pot.vdf_pot(jB_df.rho(), jB_df.jmu_net());
      (*FB_lat_df)[j] = pot.vdf_force(
          jB_df.rho(), jB_df.drho_dxnu().x0(), jB_df.drho_dxnu().threevec(),
          jB_df.grad_rho_cross_vecj(), jB_df.jmu_net().x0(), jB_df.grad_j0(),
          jB_df.jmu_net().threevec(), jB_df.dvecj_dt(), jB_df.curl_vecj());
    }
    update_momenta(P, dt, pot, FB_lat_df.get(), nullptr, nullptr, nullptr);
  }
  // final mean field energy
  E_final = calculate_mean_field_energy(pot, *jmu_B_lat_df, nullptr, param);
  /*
   * We compare the density histogram at 0.1*saturation density,0.5*saturation
   * density and 0.25*saturation density, which is the mean density used for
   * initialisation. Due to spinodal decomposition, the system should have moved
   * out of equilibrium and we expect to observe cells with around
   * 0.5*saturation density which were not present before as well as more cells
   * at 0.1*saturation density than at 0.25*saturation density
   */
  bool final_state = density_hist(&P[0], box_length, cell_length, Ntest,
                                  density / pot.saturation_density(), 0.1, 0.5,
                                  0.1, pot.saturation_density());

  VERIFY(final_state) << "The expected peak structure was not found. Spinodal "
                         "decomposition has probably not happened.";
  // for spinodal decomposition, we expect the absolute value of the mean-field
  // energy to increase substantially
  VERIFY(E_final / E_init > 1.3)
      << "Mean-field energy did not increase enough. This means that spinodal "
         "decomposition has probably not happened";
}
