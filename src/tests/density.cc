/*
 *
 *    Copyright (c) 2014-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/density.h"

#include <filesystem>
#include <map>

#include "setup.h"
#include "smash/boxmodus.h"
#include "smash/configuration.h"
#include "smash/experiment.h"
#include "smash/modusdefault.h"
#include "smash/nucleus.h"
#include "smash/thermodynamicoutput.h"

using namespace smash;

static const std::filesystem::path testoutputpath =
    std::filesystem::absolute(SMASH_TEST_OUTPUT_PATH);

TEST(directory_is_created) {
  std::filesystem::create_directories(testoutputpath);
  VERIFY(std::filesystem::exists(testoutputpath));
}

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "N+ 0.938 0.0 + 2212\n"
      "N0 0.938 0.0 + 2112\n"
      "π⁺ 0.138 0.0 -  211\n"
      "π⁰ 0.138 0.0 -  111\n");
}

static ParticleData create_proton(int id = -1) {
  return ParticleData{ParticleType::find(0x2212), id};
}

static ParticleData create_antiproton(int id = -1) {
  return ParticleData{ParticleType::find(-0x2212), id};
}

TEST(density_type) {
  // pions
  const ParticleType &pi_zero = ParticleType::find(0x111);
  const ParticleType &pi_plus = ParticleType::find(0x211);
  const ParticleType &pi_minus = ParticleType::find(-0x211);
  // baryons
  const ParticleType &proton = ParticleType::find(0x2212);

  // verify that pions are recognized as pions
  COMPARE(density_factor(pi_zero, DensityType::Pion), 1.);
  COMPARE(density_factor(pi_plus, DensityType::Pion), 1.);
  COMPARE(density_factor(pi_minus, DensityType::Pion), 1.);

  // verify that pions are not recognized as baryons
  COMPARE(density_factor(pi_zero, DensityType::Baryon), 0.);

  // verify that protons are recognized as baryons
  COMPARE(density_factor(proton, DensityType::Baryon), 1.);

  // verify that protons are not recognized as pions
  COMPARE(density_factor(proton, DensityType::Pion), 0.);

  // verify that all are recognized as particles
  COMPARE(density_factor(proton, DensityType::Hadron), 1.);
  COMPARE(density_factor(pi_zero, DensityType::Hadron), 1.);
  COMPARE(density_factor(pi_plus, DensityType::Hadron), 1.);
  COMPARE(density_factor(pi_minus, DensityType::Hadron), 1.);
}

// create one particle moving along x axis and check density in comp. frame
// check if density in the comp. frame gets contracted as expected
TEST(density_value) {
  ParticleData part_x = create_proton();
  part_x.set_4momentum(FourVector(1.0, 0.95, 0.0, 0.0));
  part_x.set_4position(FourVector(0.0, 0.0, 0.0, 0.0));
  ParticleList P;
  P.push_back(part_x);
  ThreeVector r;
  const ExperimentParameters exp_par = smash::Test::default_parameters();
  const DensityParameters par(exp_par);
  double rho;
  DensityType bar_dens = DensityType::Baryon;

  /* Got numbers executing mathematica code:
     rhoYZ = SetPrecision[Exp[-1.0/2.0]/(2*Pi)^(3/2), 16]
     rhoX = SetPrecision[Exp[-1.0/2.0/(1.0 - 0.95*0.95)]/(2*Pi)^(3/2), 16]
   */

  const double f =
      1.0 / smearing_factor_rcut_correction(exp_par.gauss_cutoff_in_sigma);
  r = ThreeVector(1.0, 0.0, 0.0);
  rho = std::get<0>(current_eckart(r, P, par, bar_dens, false, true));
  COMPARE_RELATIVE_ERROR(rho, 0.0003763388107782538 * f, 1.e-5);

  r = ThreeVector(0.0, 1.0, 0.0);
  rho = std::get<0>(current_eckart(r, P, par, bar_dens, false, true));
  COMPARE_RELATIVE_ERROR(rho, 0.03851083689074894 * f, 1.e-5);

  r = ThreeVector(0.0, 0.0, 1.0);
  rho = std::get<0>(current_eckart(r, P, par, bar_dens, false, true));
  COMPARE_RELATIVE_ERROR(rho, 0.03851083689074894 * f, 1.e-5);
}

TEST(density_eckart_special_cases) {
  /* This one checks one especially nasty case:
  Eckart rest frame of baryon density for proton and antiproton
  flying with the same speed in the opposite directions. */
  ParticleData pr = create_proton();
  ParticleData apr = create_antiproton();
  pr.set_4position(FourVector(0.0, 0.0, 0.0, 0.0));
  apr.set_4position(FourVector(0.0, 0.0, 0.0, 0.0));
  pr.set_4momentum(FourVector(1.0, 0.95, 0.0, 0.0));
  apr.set_4momentum(FourVector(1.0, -0.95, 0.0, 0.0));
  ParticleList P;
  P.push_back(pr);
  P.push_back(apr);
  ThreeVector r(1.0, 0.0, 0.0);
  const ExperimentParameters exp_par = smash::Test::default_parameters();
  const DensityParameters par(exp_par);
  const double f =
      1.0 / smearing_factor_rcut_correction(exp_par.gauss_cutoff_in_sigma);
  double rho =
      std::get<0>(current_eckart(r, P, par, DensityType::Baryon, false, true));
  COMPARE_ABSOLUTE_ERROR(rho, 0.0, 1.e-15) << rho;

  /* Now check for negative baryon density from antiproton */
  P.erase(P.begin());  // Remove proton
  rho =
      std::get<0>(current_eckart(r, P, par, DensityType::Baryon, false, true));
  COMPARE_RELATIVE_ERROR(rho, -0.0003763388107782538 * f, 1.e-5) << rho;
}

TEST(smearing_factor_normalization) {
  // Create density lattice with small lattice spacing
  const std::array<double, 3> l = {10., 10., 10.};
  const std::array<int, 3> n = {50, 60, 70};
  const std::array<double, 3> origin = {0., 0., 0.};
  bool periodicity = true;
  auto lat = std::make_unique<DensityLattice>(l, n, origin, periodicity,
                                              LatticeUpdate::EveryTimestep);
  // Create box with 1 proton
  const int N = 1;
  const double L = 10.;
  Configuration conf{R"(
    Modi:
      Box:
        Initial_Condition: "thermal momenta"
        Temperature: 0.2
        Start_Time: 0.0
  )"};
  // Note that it is not possible here to use set_value passing in the database
  // InputKeys::modi_box_initialMultiplicities key as its value is of type
  // std::map<PdgCode, int> and this cannot be set by YAML library automatically
  conf.merge_yaml(InputKeys::modi_box_initialMultiplicities.as_yaml(
      "{2212: " + std::to_string(N) + "}"));
  conf.set_value(InputKeys::modi_box_length, L);
  ExperimentParameters par = smash::Test::default_parameters();
  par.box_length = L;
  const DensityParameters dens_par = DensityParameters(par);
  std::unique_ptr<BoxModus> b =
      std::make_unique<BoxModus>(std::move(conf), par);
  std::vector<Particles> ensembles(1);
  b->initial_conditions(&ensembles[0], par);
  // Fill lattice from particles
  update_lattice_accumulating_ensembles(lat.get(), LatticeUpdate::EveryTimestep,
                                        DensityType::Baryon, dens_par,
                                        ensembles, false);
  // Compute integral rho(r) d^r. Should be equal to N.
  double int_rho_r_d3r = 0.0;
  for (auto &node : *lat) {
    int_rho_r_d3r += node.rho();
  }
  int_rho_r_d3r *= L * L * L / n[0] / n[1] / n[2] / N;
  COMPARE_RELATIVE_ERROR(int_rho_r_d3r, 1.0, 3.e-6);
}

TEST(smearing_factor_rcut_correction) {
  FUZZY_COMPARE(smearing_factor_rcut_correction(3.0), 0.97070911346511177);
  FUZZY_COMPARE(smearing_factor_rcut_correction(4.0), 0.99886601571021467);
}

// check that analytical and numerical results for gradient of density coincide
TEST(density_gradient) {
  // create two protons
  ParticleData part1 = create_proton();
  ParticleData part2 = create_proton();
  double mass = 0.938;
  // set momenta (just randomly):
  part1.set_4momentum(mass, ThreeVector());
  part2.set_4momentum(mass, ThreeVector());
  // set coordinates (not too far from each other)
  part1.set_4position(FourVector(0.0, 0.0, 0.0, 0.0));
  part2.set_4position(FourVector(0.0, 0.5, 0.0, 0.0));
  // make particle list out of them
  ParticleList P;
  P.push_back(part1);
  P.push_back(part2);

  const ExperimentParameters par = smash::Test::default_parameters();
  ThreeVector r, dr;
  FourVector jmu;
  DensityType dtype = DensityType::Baryon;

  ThreeVector num_grad, analit_grad;
  r = ThreeVector(0.5, 0.3, 0.0);
  auto rho_and_grad = current_eckart(r, P, par, dtype, true, true);
  double rho_r = std::get<0>(rho_and_grad);

  // analytical gradient
  analit_grad = std::get<2>(rho_and_grad);
  // numerical gradient
  dr = ThreeVector(1.e-4, 0.0, 0.0);
  double rho_rdr =
      std::get<0>(current_eckart(r + dr, P, par, dtype, false, true));
  num_grad.set_x1((rho_rdr - rho_r) / dr.x1());
  dr = ThreeVector(0.0, 1.e-4, 0.0);
  rho_rdr = std::get<0>(current_eckart(r + dr, P, par, dtype, false, true));
  num_grad.set_x2((rho_rdr - rho_r) / dr.x2());
  dr = ThreeVector(0.0, 0.0, 1.e-4);
  rho_rdr = std::get<0>(current_eckart(r + dr, P, par, dtype, false, true));
  num_grad.set_x3((rho_rdr - rho_r) / dr.x3());
  // compare them with: accuracy should not be worse than |dr|
  COMPARE_ABSOLUTE_ERROR(num_grad.x1(), analit_grad.x1(), 1.e-4);
  COMPARE_ABSOLUTE_ERROR(num_grad.x2(), analit_grad.x2(), 1.e-4);
  COMPARE_ABSOLUTE_ERROR(num_grad.x3(), analit_grad.x3(), 1.e-4);
}

TEST(density_gradient_in_linear_box) {
  // set parameters for the test
  ExperimentParameters par = smash::Test::default_parameters();
  par.testparticles = 10000;
  par.gaussian_sigma = 0.5;
  par.gauss_cutoff_in_sigma = 8.0;
  DensityType dtype = DensityType::Baryon;
  /* Prepare a box of protons whose density is linearly increasing
   * along z-axis. The box is 2*2*2 fm^3 large. There're 10000000
   * particles, including test particles, inside.*/
  ParticleList P;
  for (int i = 0; i < 100 * par.testparticles; i++) {
    const double x = random::uniform(-2., 2.);
    const double y = random::uniform(-2., 2.);
    const double z = random::power(1., 0., 4.) - 2.;
    ParticleData part = create_proton();
    const double mass = 0.938;
    part.set_4momentum(mass, ThreeVector());
    part.set_4position(FourVector(0., x, y, z));
    P.push_back(part);
  }
  // set the position where the gradient is measured.
  const double x0 = random::uniform(-0.4, 0.4);
  const double y0 = random::uniform(-0.4, 0.4);
  const double z0 = random::uniform(-0.4, 0.4);
  const ThreeVector r0 = ThreeVector(x0, y0, z0);
  // calculate the gradients
  const auto rho_and_grad = current_eckart(r0, P, par, dtype, true, true);
  const auto drho_dr = std::get<2>(rho_and_grad);
  /* Theoretically, the gradient should be (0., 0., DU/Dz) at any
   * point inside the box, where DU/Dz can be caluclated by taking
   * the difference between the potentials evaluated at (0., 0., -0.4)
   * and (0., 0., 0.4). Compare the z-component of the gradient with the
   * theoretical value. */
  const double theo_drho_dz =
      (std::get<0>(current_eckart(ThreeVector(0., 0., 0.4), P, par, dtype,
                                  false, true)) -
       std::get<0>(current_eckart(ThreeVector(0., 0., -0.4), P, par, dtype,
                                  false, true))) /
      0.8;
  COMPARE_RELATIVE_ERROR(drho_dr.x3(), theo_drho_dz, 0.02);
  /* Meanwhile, the gradient should be along z-axis, which means its
   * transverse component should be much smaller than its longitudinal
   * component*/
  const double drho_T_over_z =
      sqrt(drho_dr.x1() * drho_dr.x1() + drho_dr.x2() * drho_dr.x2()) /
      drho_dr.x3();
  COMPARE_ABSOLUTE_ERROR(drho_T_over_z, 0., 0.02);
}

TEST(current_curl_in_rotating_box) {
  // set parameters for the test
  ExperimentParameters par = smash::Test::default_parameters();
  par.testparticles = 10000;
  par.gaussian_sigma = 0.65;
  par.gauss_cutoff_in_sigma = 8.0;
  DensityType dtype = DensityType::Baryon;
  /* Prepare a box of uniformly distributed protons which are rotating about
   * the z-axis with a constant angular velocity. The box is 2*2*2 fm^3 large.
   * There're 10000000 particles, including test particles, inside.*/
  ParticleList P;
  const double omega = 0.2;
  for (int i = 0; i < 100 * par.testparticles; i++) {
    const double x = random::uniform(-2., 2.);
    const double y = random::uniform(-2., 2.);
    const double z = random::uniform(-2., 2.);
    const double r2 = x * x + y * y + z * z;
    const double gamma = 1. / sqrt(1. - omega * omega * r2);
    ParticleData part = create_proton();
    const double mass = 0.9;
    part.set_4momentum(mass, -mass * omega * gamma * y,
                       mass * omega * gamma * x, 0.);
    part.set_4position(FourVector(0., x, y, z));
    P.push_back(part);
  }
  // set the location where the curl is measured.
  const double x0 = random::uniform(-0.4, 0.4);
  const double y0 = random::uniform(-0.4, 0.4);
  const double z0 = random::uniform(-0.4, 0.4);
  const ThreeVector r0 = ThreeVector(x0, y0, z0);
  // calculate the curl
  const auto rho_and_grad = current_eckart(r0, P, par, dtype, true, true);
  const auto rot_j = std::get<3>(rho_and_grad);
  /* Theoretically, the curl should be (0., 0., 2. * density * omega) at any
   * point inside the box, Compare the z-component of the curl with the
   * theoretical value. */
  const double theo_rot_j_z = 2. * 1.5 * omega;
  COMPARE_RELATIVE_ERROR(rot_j.x3(), theo_rot_j_z, 0.05);
  /* Meanwhile, the curl should be along z-axis, which means its
   * transverse component should be much smaller than its longitudinal
   * component*/
  const double rot_j_T_over_z =
      sqrt(rot_j.x1() * rot_j.x1() + rot_j.x2() * rot_j.x2()) / rot_j.x3();
  COMPARE_ABSOLUTE_ERROR(rot_j_T_over_z, 0., 0.01);
}

TEST(baryon_current_j_B_smearing_false) {
  // create a proton
  ParticleData particle = create_proton();
  double mass = 0.938;
  double px = 1.0;
  double E0 = sqrt(px * px + mass * mass);
  // set momenta:
  particle.set_4momentum(mass, px, 0.0, 0.0);
  particle.set_4position(FourVector(0.0, 0.0, 0.0, 0.0));
  // make particle list out of them
  ParticleList P;
  P.push_back(particle);
  // set default parameters
  const ExperimentParameters par = smash::Test::default_parameters();
  // define the position for the calculation. The result should not depend on r0
  const ThreeVector r0 = ThreeVector(0.0, 0.0, 0.0);
  // calculate j_B
  DensityType dtype = DensityType::Baryon;
  // calculate the density
  bool comp_gradient = false;
  bool smearing = false;
  const auto j_mu_B =
      std::get<1>(current_eckart(r0, P, par, dtype, comp_gradient, smearing));
  /* In the case of no smearing, the four current j_B is calculated as j_B^k =
   * \sum_i B_i * p^k_i / p^0_i with the sum running over all particles. Since
   * in the above system there is only one proton, the result should simply be
   * j^\mu = (1, p^1 / p^0, 0, 0). Note that in the case of no smearing the
   * function returns just the current and not current density. As a result the
   * volume plays no role here.
   */
  COMPARE_ABSOLUTE_ERROR(j_mu_B[0], 1.0, 1e-12);
  COMPARE_ABSOLUTE_ERROR(j_mu_B[1], px / E0, 1e-12);
  COMPARE_ABSOLUTE_ERROR(j_mu_B[2], 0.0, 1e-12);
  COMPARE_ABSOLUTE_ERROR(j_mu_B[3], 0.0, 1e-12);
}

/*
   This test does not compare anything. It only prints density map versus
   time to vtk files, so that one can open it with paraview and make sure
   that Lorentz contraction works properly for density.
 */
/*
TEST(density_eckart_frame) {

  // create a particle list with particles moving in different directions.

  // proton that moves along the x axis:
  ParticleData proton_x = create_proton();
  // proton that moves along the y axis:
  ParticleData proton_y = create_proton();
  // proton that moves along the z axis:
  ParticleData proton_z = create_proton();
  // proton that moves in neither of xyz planes, v = v0:
  ParticleData proton = create_proton();
  // proton that moves in neither of xyz planes, v = - v0:
  ParticleData antiproton = create_antiproton();

  double mass = 0.938;
  // set momenta:
  proton_x.set_4momentum(mass, ThreeVector(2.0, 0.0, 0.0));
  proton_y.set_4momentum(mass, ThreeVector(0.0, -2.0, 0.0));
  proton_z.set_4momentum(mass, ThreeVector(0.0, 0.0, 2.0));
  proton.set_4momentum(mass, ThreeVector(0.4, 0.5, 0.7));
  antiproton.set_4momentum(mass, ThreeVector(-0.4, -0.5, -0.7));

  // set positions:
  proton_x.set_4position(FourVector(0.0, -5.0, 0.0, 0.0));
  proton_y.set_4position(FourVector(0.0, 0.0, 5.0, 0.0));
  proton_z.set_4position(FourVector(0.0, 0.0, 0.0, -5.0));
  proton.set_4position(FourVector(0.0,-4.0,-5.0,-7.0));
  antiproton.set_4position(FourVector(0.0, 4.0, 5.0, 7.0));

  // add particles (and make sure the particles get the correct ID):
  COMPARE(P.add_data(proton_x), 0);
  COMPARE(P.add_data(proton_y), 1);
  COMPARE(P.add_data(proton_z), 2);
  COMPARE(P.add_data(proton), 3);
  COMPARE(P.add_data(antiproton), 4);


  ModusDefault m;
  Particles Pdef;
  create_particle_list(Pdef);
  OutputsList out;
  // clock, output interval, cross-section, testparticles, gauss. sigma
  ExperimentParameters param = smash::Test::default_parameters();
  double dx = 0.3;
  double dy = 0.3;
  double dz = 0.3;
  int nx = 20;
  int ny = 20;
  int nz = 20;
  double sigma = 0.8;
  ThreeVector r;
  DensityType bar_dens = baryon;
  ParticleList plist;

  for (auto it = 0; it < 30; it++) {
    plist = ParticleList(Pdef.data().begin(), Pdef.data().end());
    vtk_density_map( ("bdens.vtk." + std::to_string(it)).c_str(),
                     plist, sigma, bar_dens, 1,
                     nx, ny, nz, dx, dy, dz);
    m.propagate(&Pdef, param, out);
  }
}
*/

/*
  Generates nuclei and prints their density profiles to vtk files
*/
/*
TEST(nucleus_density) {
  std::string configfilename = "densconfig.yaml";
  std::filesystem::ofstream(testoutputpath / configfilename) << "x: 0\ny: 0\nz:
0\n"; VERIFY(std::filesystem::exists(testoutputpath / configfilename));

  // Lead nuclei with 1000 test-particles
  std::map<PdgCode, int> lead_list = {{0x2212, 79}, {0x2112, 118}};
  int Ntest = 1000;
  Nucleus lead;
  lead.fill_from_list(lead_list, Ntest);
  lead.set_parameters_automatic();
  lead.arrange_nucleons();

  // make particle list out of generated nucleus
  Particles p;
  lead.copy_particles(&p);
  ParticleList plist = p.copy_to_vector();

  // write density profile to file, time-consuming!
  DensityType dens_type = DensityType::Baryon;
  double sigma = 0.5; // fm
//  vtk_density_map("lead_density.vtk", plist, sigma, dens_type, Ntest,
//                     20, 20, 20, 0.5, 0.5, 0.5);
  ThreeVector lstart(-10.0, 0.0, 0.0);
  ThreeVector lend(10.0, 0.0, 0.0);
  const int npoints = 100;

  Configuration&& conf{testoutputpath, configfilename};
  ExperimentParameters par = smash::Test::default_parameters(Ntest);
  std::unique_ptr<ThermodynamicOutput> out =
std::make_unique<ThermodynamicOutput>(testoutputpath, std::move(conf));
  out->density_along_line("lead_densityX.dat", plist, par, dens_type,
                          lstart, lend, npoints);

  lstart = ThreeVector(0.0, -10.0, 0.0);
  lend = ThreeVector(0.0, 10.0, 0.0);
  out->density_along_line("lead_densityY.dat", plist, par, dens_type,
                          lstart, lend, npoints);

  lstart = ThreeVector(0.0, 0.0, -10.0);
  lend = ThreeVector(0.0, 0.0, 10.0);
  out->density_along_line("lead_densityZ.dat", plist, par, dens_type,
                          lstart, lend, npoints);
}*/

/*TEST(box_density) {
ParticleType::create_type_list(
    "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
    "proton 0.938 0.0 + 2212\n");
const int Ntest = 1000;
const double L = 10.0;
Configuration conf(TEST_CONFIG_PATH);
conf["Modus"] = "Box";
conf.take({"Modi", "Box", "Init_Multiplicities"});
conf["Modi"]["Box"]["Init_Multiplicities"]["2212"] = 1000;
conf["Modi"]["Box"]["Length"] = L;
const ExperimentParameters par = smash::Test::default_parameters(Ntest);
std::unique_ptr<BoxModus> b = std::make_unique<BoxModus>(conf["Modi"], par);
Particles P;
b->initial_conditions(&P, par);
ParticleList plist = P.copy_to_vector();
conf["Output"]["Density"]["x"] = 0.0;
conf["Output"]["Density"]["y"] = 0.0;
conf["Output"]["Density"]["z"] = 0.0;
std::unique_ptr<ThermodynamicOutput> out =
std::make_unique<ThermodynamicOutput>(testoutputpath,
                                         conf["Output"]["Density"]);
const ThreeVector lstart = ThreeVector(0.0, 0.0, 0.0);
const ThreeVector lend = ThreeVector(L, L, L);
const int npoints = 100;
ExperimentParameters par = smash::Test::default_parameters(Ntest);
out->density_along_line("box_density.dat", plist, par, DensityType::Baryon,
                        lstart, lend, npoints);
}*/
