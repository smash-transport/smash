/*
 *
 *    Copyright (c) 2014-2015,2017-2020
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
  dnucleus.set_azimuthal_angle(M_PI / 2);
  FourVector expectation = FourVector(0., 1., 1., 0.);
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    i->set_4position(FourVector(0., 1., 0., 1.));
  }
  dnucleus.rotate();
  FourVector actual;
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    actual = i->position();
  }
  COMPARE_ABSOLUTE_ERROR(actual.x0(), expectation.x0(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x1(), expectation.x1(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x2(), expectation.x2(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x3(), expectation.x3(), 1e-7);
}

TEST(rotate_theta) {
  DeformedNucleus dnucleus(small_list, 1);
  // Plan is to rotate the (0, 0, 0, -1) vector by theta=pi/2
  // Rotation by pi/2 means (0, 0, 0, -1) -> (0, 0, 1, 0)
  dnucleus.set_polar_angle(M_PI / 2);
  FourVector expectation = FourVector(0., 0., 1., 0.);
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    i->set_4position(FourVector(0., 0., 0., -1.));
  }
  dnucleus.rotate();
  FourVector actual;
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    actual = i->position();
  }
  COMPARE_ABSOLUTE_ERROR(actual.x0(), expectation.x0(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x1(), expectation.x1(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x2(), expectation.x2(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x3(), expectation.x3(), 1e-7);
}

TEST(rotate_both) {
  DeformedNucleus dnucleus(small_list, 1);
  // Plan is to rotate the (0, 1, 1, 0) vector by phi=pi
  // and then by theta=pi around the rotated x-axis
  // Result: (0, 1, 1, 0) -> (0, -1, -1, 0) -> (0, -1, 1, 0)
  dnucleus.set_azimuthal_angle(M_PI);
  dnucleus.set_polar_angle(M_PI);
  FourVector expectation = FourVector(0., -1., 1., 0.);
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    i->set_4position(FourVector(0., 1., 1., 0.));
  }
  dnucleus.rotate();
  FourVector actual;
  for (auto i = dnucleus.begin(); i != dnucleus.end(); i++) {
    actual = i->position();
  }
  COMPARE_ABSOLUTE_ERROR(actual.x0(), expectation.x0(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x1(), expectation.x1(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x2(), expectation.x2(), 1e-7);
  COMPARE_ABSOLUTE_ERROR(actual.x3(), expectation.x3(), 1e-7);
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
  Configuration conf = Test::configuration();
  conf.set_value({"Modi", "Collider", "Projectile", "Particles", "2112"}, 118);
  conf.set_value({"Modi", "Collider", "Projectile", "Particles", "2212"}, 79);
  conf.set_value({"Modi", "Collider", "Projectile", "Saturation_Density"},
                 0.1968);
  conf.set_value({"Modi", "Collider", "Projectile", "Diffusiveness"}, 0.8);
  conf.set_value({"Modi", "Collider", "Projectile", "Radius"}, 2.0);
  // inserts beta2_ and beta4_ values
  conf.set_value({"Modi", "Collider", "Projectile", "Deformed", "Beta_2"}, 1);
  conf.set_value({"Modi", "Collider", "Projectile", "Deformed", "Beta_4"}, 2);

  // verifies if the beta values have been transcribed correctly
  Configuration mod_conf = conf["Modi"];
  Configuration col_conf = mod_conf["Collider"];
  Configuration proj_conf = col_conf["Projectile"];
  DeformedNucleus dnucleus(proj_conf, 1, 0);
  VERIFY(dnucleus.get_beta2() == 1);
  VERIFY(dnucleus.get_beta4() == 2);
}

TEST(set_deformation_parameters_automatic) {
  // config for uranium nucleus
  Configuration conf1 = Test::configuration();
  conf1.set_value({"Modi", "Collider", "Projectile", "Particles", "2112"}, 146);
  conf1.set_value({"Modi", "Collider", "Projectile", "Particles", "2212"}, 92);
  conf1.set_value({"Modi", "Collider", "Projectile", "Saturation_Density"},
                  0.1968);
  conf1.set_value({"Modi", "Collider", "Projectile", "Diffusiveness"}, 1.0);
  conf1.set_value({"Modi", "Collider", "Projectile", "Radius"}, 1.0);
  conf1.set_value({"Modi", "Collider", "Projectile", "Deformed", "Automatic"},
                  "True");

  // verifies that the values were automatically set
  Configuration mod_conf1 = conf1["Modi"];
  Configuration col_conf1 = mod_conf1["Collider"];
  Configuration proj_conf1 = col_conf1["Projectile"];
  DeformedNucleus dnucleus1(proj_conf1, 1, 1);
  VERIFY(dnucleus1.get_beta2() == 0.28);
  VERIFY(dnucleus1.get_beta4() == 0.093);

  // config for copper nucleus
  Configuration conf2 = Test::configuration();
  conf2.set_value({"Modi", "Collider", "Projectile", "Particles", "2112"}, 34);
  conf2.set_value({"Modi", "Collider", "Projectile", "Particles", "2212"}, 29);
  conf2.set_value({"Modi", "Collider", "Projectile", "Saturation_Density"},
                  0.1968);
  conf2.set_value({"Modi", "Collider", "Projectile", "Diffusiveness"}, 1.0);
  conf2.set_value({"Modi", "Collider", "Projectile", "Radius"}, 1.0);
  conf2.set_value({"Modi", "Collider", "Projectile", "Deformed", "Automatic"},
                  "True");

  // verifies that the values were automatically set
  Configuration mod_conf2 = conf2["Modi"];
  Configuration col_conf2 = mod_conf2["Collider"];
  Configuration proj_conf2 = col_conf2["Projectile"];
  DeformedNucleus dnucleus2(proj_conf2, 1, 1);
  VERIFY(dnucleus2.get_beta2() == 0.162);
  VERIFY(dnucleus2.get_beta4() == -0.006);

  // config for Zirconium nucleus
  Configuration conf3 = Test::configuration();
  conf3.set_value({"Modi", "Collider", "Projectile", "Particles", "2112"}, 56);
  conf3.set_value({"Modi", "Collider", "Projectile", "Particles", "2212"}, 40);
  conf3.set_value({"Modi", "Collider", "Projectile", "Saturation_Density"},
                  0.1968);
  conf3.set_value({"Modi", "Collider", "Projectile", "Diffusiveness"}, 1.0);
  conf3.set_value({"Modi", "Collider", "Projectile", "Radius"}, 1.0);
  conf3.set_value({"Modi", "Collider", "Projectile", "Deformed", "Automatic"},
                  "True");

  // verifies that the values were automatically set
  Configuration mod_conf3 = conf3["Modi"];
  Configuration col_conf3 = mod_conf3["Collider"];
  Configuration proj_conf3 = col_conf3["Projectile"];
  DeformedNucleus dnucleus3(proj_conf3, 1, 1);
  VERIFY(dnucleus3.get_beta2() == 0.0);
  VERIFY(dnucleus3.get_beta4() == 0.0);

  // config for Ruthenium nucleus
  Configuration conf4 = Test::configuration();
  conf4.set_value({"Modi", "Collider", "Projectile", "Particles", "2112"}, 52);
  conf4.set_value({"Modi", "Collider", "Projectile", "Particles", "2212"}, 44);
  conf4.set_value({"Modi", "Collider", "Projectile", "Saturation_Density"},
                  0.1968);
  conf4.set_value({"Modi", "Collider", "Projectile", "Diffusiveness"}, 1.0);
  conf4.set_value({"Modi", "Collider", "Projectile", "Radius"}, 1.0);
  conf4.set_value({"Modi", "Collider", "Projectile", "Deformed", "Automatic"},
                  "True");

  // verifies that the values were automatically set
  Configuration mod_conf4 = conf4["Modi"];
  Configuration col_conf4 = mod_conf4["Collider"];
  Configuration proj_conf4 = col_conf4["Projectile"];
  DeformedNucleus dnucleus4(proj_conf4, 1, 1);
  VERIFY(dnucleus4.get_beta2() == 0.158);
  VERIFY(dnucleus4.get_beta4() == 0.0);
}

TEST(nucleon_density) {
  // config with values for an easy analytic deformed-woods-saxon value
  // Uranium core with default values
  Configuration conf1 = Test::configuration();
  conf1.set_value({"Modi", "Collider", "Projectile", "Particles", "2112"}, 146);
  conf1.set_value({"Modi", "Collider", "Projectile", "Particles", "2212"}, 92);
  conf1.set_value({"Modi", "Collider", "Projectile", "Saturation_Density"},
                  0.166);
  conf1.set_value({"Modi", "Collider", "Projectile", "Diffusiveness"}, 0.556);
  conf1.set_value({"Modi", "Collider", "Projectile", "Radius"}, 6.86);

  // verifies that deformed Woods-Saxon is indeed 0 for some arbitrary values
  Configuration mod_conf1 = conf1["Modi"];
  Configuration col_conf1 = mod_conf1["Collider"];
  Configuration proj_conf1 = col_conf1["Projectile"];
  DeformedNucleus dnucleus1(proj_conf1, 1, 0);
  COMPARE_ABSOLUTE_ERROR(dnucleus1.nucleon_density(.0892, .1802, 0.),
                         0.16599914, 1e-7);

  // config with values for an easy analytic deformed Woods-Saxon value
  // Lead core with default values
  Configuration conf2 = Test::configuration();
  conf2.set_value({"Modi", "Collider", "Projectile", "Particles", "2112"}, 126);
  conf2.set_value({"Modi", "Collider", "Projectile", "Particles", "2212"}, 82);
  conf2.set_value({"Modi", "Collider", "Projectile", "Saturation_Density"},
                  0.161);
  conf2.set_value({"Modi", "Collider", "Projectile", "Diffusiveness"}, 0.54);
  conf2.set_value({"Modi", "Collider", "Projectile", "Radius"}, 6.67);

  // verifies that deformed Woods-Saxon is indeed 0.5
  Configuration mod_conf2 = conf2["Modi"];
  Configuration col_conf2 = mod_conf2["Collider"];
  Configuration proj_conf2 = col_conf2["Projectile"];
  DeformedNucleus dnucleus2(proj_conf2, 1, 0);
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
      return twopi * square(r) * nucl.nucleon_density(r, cosx, 0.0) / square(t);
    });
    const size_t Z = nucl.number_of_protons(), A = nucl.number_of_particles();
    std::cout << "Z: " << Z << "  A: " << A << std::endl;
    std::cout << result.value() << " ± " << result.error() << std::endl;
    size_t index = &nucl - &deformed_nuclei[0];
    COMPARE_ABSOLUTE_ERROR(result.value(), static_cast<double>(A),
                           allowed_errors[index]);
  }
}
