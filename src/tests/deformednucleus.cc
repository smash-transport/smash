/*
 *
 *    Copyright (c) 2014-2015,2017-2020,2022,2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/deformednucleus.h"

#include <map>
#include <vector>

#include "setup.h"
#include "smash/constants.h"
#include "smash/fourvector.h"
#include "smash/nucleus.h"
#include "smash/particledata.h"
#include "smash/pdgcode.h"
#include "smash/pow.h"

namespace particles_txt {
#include <particles.txt.h>
}

using namespace smash;

std::map<PdgCode, int> small_list = {{0x2212, 1}};

TEST(init_particle_types) {
  ParticleType::create_type_list(particles_txt::data);
}

TEST(rotate_phi) {
  DeformedNucleus dnucleus(small_list, 1);
  // Plan is to rotate the (0, 1, 0, 1) vector by phi=pi/2.
  // Rotation by pi/2 means (0, 1, 0, 1) -> (0, 0, 1, 1)
  Configuration config{
      InputKeys::modi_collider_projectile_orientation_phi
          .as_yaml("1.57079632679489661923")  // as string to retain all digits
          .c_str()};
  dnucleus.set_orientation_from_config(config);
  FourVector expectation = FourVector(0., 0., 1., 1.);
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    i->set_4position(FourVector(0., 1., 0., 1.));
  }
  dnucleus.rotate();
  FourVector actual;
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    actual = i->position();
  }
  COMPARE_ABSOLUTE_ERROR(actual.x0(), expectation.x0(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x1(), expectation.x1(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x2(), expectation.x2(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x3(), expectation.x3(), 1e-14);
}

TEST(rotate_theta) {
  DeformedNucleus dnucleus(small_list, 1);
  // Plan is to rotate the (0, 0, 0, -1) vector by theta=pi/2
  // Rotation by pi/2 means (0, 0, 0, -1) -> (0, 0, 1, 0)
  Configuration config{
      InputKeys::modi_collider_projectile_orientation_theta
          .as_yaml("1.57079632679489661923")  // as string to retain all digits
          .c_str()};
  dnucleus.set_orientation_from_config(config);
  FourVector expectation = FourVector(0., 0., 1., 0.);
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    i->set_4position(FourVector(0., 0., 0., -1.));
  }
  dnucleus.rotate();
  FourVector actual;
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    actual = i->position();
  }
  COMPARE_ABSOLUTE_ERROR(actual.x0(), expectation.x0(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x1(), expectation.x1(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x2(), expectation.x2(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x3(), expectation.x3(), 1e-14);
}

TEST(rotate_both) {
  DeformedNucleus dnucleus(small_list, 1);
  // Plan is to rotate the (0, 1, 1, 0) vector by phi=pi
  // and then by theta=pi around the rotated x-axis
  // Result: (0, 1, 1, 0) -> (0, -1, -1, 0) -> (0, -1, 1, 0)
  Configuration config{R"(
    Modi:
      Collider:
        Projectile:
          Orientation:
            Theta: 3.14159265358979323846
            Phi: 3.14159265358979323846
    )"};
  dnucleus.set_orientation_from_config(config);
  FourVector expectation = FourVector(0., -1., 1., 0.);
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    i->set_4position(FourVector(0., 1., 1., 0.));
  }
  dnucleus.rotate();
  FourVector actual;
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    actual = i->position();
  }
  COMPARE_ABSOLUTE_ERROR(actual.x0(), expectation.x0(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x1(), expectation.x1(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x2(), expectation.x2(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x3(), expectation.x3(), 1e-14);
}

TEST(rotate_all_three) {
  DeformedNucleus dnucleus(small_list, 1);
  // Plan is to rotate the (0, 1, 1, 0) vector by phi=pi around the z-axis, then
  // by theta=pi around the rotated x-axis and finally by psi=pi around the
  // rotated z-axis. Result: (0, 1, 1, 0) -> (0, -1, -1, 0) -> (0, -1, 1, 0) ->
  // (0, 1,-1,0)
  Configuration config{R"(
    Modi:
      Collider:
        Target:
          Orientation:
            Theta: 3.14159265358979323846
            Phi: 3.14159265358979323846
            Psi: 3.14159265358979323846
    )"};
  dnucleus.set_orientation_from_config(config);
  FourVector expectation = FourVector(0., 1., -1., 0.);
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    i->set_4position(FourVector(0., 1., 1., 0.));
  }
  dnucleus.rotate();
  FourVector actual;
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    actual = i->position();
  }
  COMPARE_ABSOLUTE_ERROR(actual.x0(), expectation.x0(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x1(), expectation.x1(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x2(), expectation.x2(), 1e-14);
  COMPARE_ABSOLUTE_ERROR(actual.x3(), expectation.x3(), 1e-14);
}

// Tests if the function for spherical harmonics
// returns correct values for given l and cosx
// simple values for cosx for demonstration purposes
TEST(ylm) {
  COMPARE_ABSOLUTE_ERROR(y_l_m(2, 0, 1., 0.), std::sqrt(5. / M_PI) / 2., 1e-7);
  COMPARE_ABSOLUTE_ERROR(y_l_m(2, 2, 1., 0.), 0.0, 1e-7);
  COMPARE_ABSOLUTE_ERROR(y_l_m(4, 0, 1., 0.), 3. / (2. * std::sqrt(M_PI)),
                         1e-7);
}

TEST(deformation_parameters_from_config) {
  // creates config for arbitrary nucleus (Gold in this case)
  Configuration conf{R"(
    Modi:
      Collider:
        Target:
          Particles:
            2112: 118
            2212: 79
          Saturation_Density: 0.1968
          Diffusiveness: 0.8
          Radius: 2.0
          Deformed:
            Beta_2: 1
            Beta_4: 2
  )"};
  // verifies if the beta values have been transcribed correctly
  DeformedNucleus dnucleus(conf, 1, false);
  VERIFY(dnucleus.get_beta2() == 1);
  VERIFY(dnucleus.get_beta4() == 2);
}

TEST(set_deformation_parameters_automatic) {
  auto create_conf = [](int n1, int n2) {
    std::string tmp{R"(
    Modi:
      Collider:
        Projectile:
          Saturation_Density: 0.1968
          Diffusiveness: 1.0
          Radius: 1.0
          Particles: )"};
    std::string particles{"{2112: " + std::to_string(n1) +
                          ", 2212: " + std::to_string(n2) + "}"};
    return Configuration{(tmp + particles).c_str()};
  };
  // config for uranium nucleus
  Configuration conf1 = create_conf(146, 92);
  // verifies that the values were automatically set
  DeformedNucleus dnucleus1(conf1, 1, true);
  VERIFY(dnucleus1.get_beta2() == 0.28);
  VERIFY(dnucleus1.get_beta4() == 0.093);

  // config for copper nucleus
  Configuration conf2 = create_conf(34, 29);
  // verifies that the values were automatically set
  DeformedNucleus dnucleus2(conf2, 1, true);
  VERIFY(dnucleus2.get_beta2() == 0.162);
  VERIFY(dnucleus2.get_beta4() == -0.006);

  // config for Zirconium nucleus
  Configuration conf3 = create_conf(56, 40);
  // verifies that the values were automatically set
  DeformedNucleus dnucleus3(conf3, 1, true);
  VERIFY(dnucleus3.get_beta2() == 0.0);
  VERIFY(dnucleus3.get_beta4() == 0.0);

  // config for Ruthenium nucleus
  Configuration conf4 = create_conf(52, 44);
  // verifies that the values were automatically set
  DeformedNucleus dnucleus4(conf4, 1, true);
  VERIFY(dnucleus4.get_beta2() == 0.158);
  VERIFY(dnucleus4.get_beta4() == 0.0);
}

TEST(nucleon_density) {
  // config with values for an easy analytic deformed-woods-saxon value
  // Uranium core with default values
  Configuration conf1{R"(
    Modi:
      Collider:
        Projectile:
          Particles:
            2112: 146
            2212: 92
          Saturation_Density: 0.166
          Diffusiveness: 0.556
          Radius: 6.86
  )"};
  // verifies that deformed Woods-Saxon is indeed 0 for some arbitrary values
  DeformedNucleus dnucleus1(conf1, 1, false);
  COMPARE_ABSOLUTE_ERROR(dnucleus1.nucleon_density(.0892, .1802, 0.),
                         0.16599914, 1e-7);

  // config with values for an easy analytic deformed Woods-Saxon value
  // Lead core with default values
  Configuration conf2{R"(
    Modi:
      Collider:
        Projectile:
          Particles:
            2112: 126
            2212: 82
          Saturation_Density: 0.161
          Diffusiveness: 0.54
          Radius: 6.67
  )"};
  // verifies that deformed Woods-Saxon is indeed 0.5
  DeformedNucleus dnucleus2(conf2, 1, false);
  COMPARE_ABSOLUTE_ERROR(dnucleus2.nucleon_density(.0892, .1802, 0.0),
                         0.16099917, 1e-7);
}

TEST(nucleon_density_norm) {
  const std::map<PdgCode, int> copper = {{pdg::p, 29}, {pdg::n, 63 - 29}},
                               zirconium = {{pdg::p, 40}, {pdg::n, 96 - 40}},
                               ruthenium = {{pdg::p, 44}, {pdg::n, 96 - 44}},
                               gold = {{pdg::p, 79}, {pdg::n, 197 - 79}},
                               lead = {{pdg::p, 82}, {pdg::n, 208 - 82}},
                               uranium = {{pdg::p, 92}, {pdg::n, 238 - 92}};
  std::vector<DeformedNucleus> deformed_nuclei{{copper, 1},    {zirconium, 1},
                                               {ruthenium, 1}, {gold, 1},
                                               {lead, 1},      {uranium, 1}};
  std::vector<double> allowed_errors = {1.0, 1.0, 1.0, 1.0, 6.0, 2.0};

  Integrator2d integrate;
  for (const DeformedNucleus &nucl : deformed_nuclei) {
    // Transform integral from (0, oo) to (0, 1) via r = (1 - t) / t.
    const auto result = integrate(0, 1, -1, 1, [&](double t, double cosx) {
      const double r = (1 - t) / t;
      return twopi * r * r * nucl.nucleon_density(r, cosx, 0.0) / (t * t);
    });
    const size_t Z = nucl.number_of_protons(), A = nucl.number_of_particles();
    logg[0].debug() << "Z: " << Z << "  A: " << A << '\n'
                    << result.value() << " ± " << result.error() << '\n';
    size_t index = &nucl - &deformed_nuclei[0];
    COMPARE_ABSOLUTE_ERROR(result.value(), static_cast<double>(A),
                           allowed_errors[index]);
  }
}
