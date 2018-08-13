/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include <map>

#include "../include/smash/nucleus.h"
#include "../include/smash/particles.h"
#include "../include/smash/pdgcode.h"
#include "../include/smash/threevector.h"

namespace particles_txt {
#include <particles.txt.h>
}

using namespace smash;

std::map<PdgCode, int> list = {{0x2212, 82}, {0x2112, 126}};

TEST(init_particle_types) {
  ParticleType::create_type_list(particles_txt::data);
}

TEST(initialize_realparticles) {
  Nucleus lead(list, 1);  // fill with 208 nucleons
  COMPARE(lead.number_of_particles(), 208u);
  COMPARE(lead.size(), 208u);
}

TEST(initialize_testparticles) {
  constexpr int N_TEST = 10;
  Nucleus lead(list, N_TEST);
  COMPARE(lead.number_of_particles(), 208u);
  COMPARE(lead.size(), 208u * N_TEST);
}

TEST(initialize_testparticles_multiple) {
  constexpr int N_TEST = 10;
  Nucleus lead(list, N_TEST);
  lead.fill_from_list(list, N_TEST);
  COMPARE(lead.number_of_particles(), 416u);
  COMPARE(lead.size(), 416u * N_TEST);
}

TEST_CATCH(initialize_testparticles_wrong, Nucleus::TestparticleConfusion) {
  constexpr int N_TEST = 10;
  Nucleus lead(list, 1);
  lead.fill_from_list(list, N_TEST);
  COMPARE(lead.number_of_particles(), 208u);
  // this should throw an error: list size is (N_TEST+1)*208 = 2288, not
  // divisible by N_TEST (unless someone set N_TEST to 1).
  size_t size = lead.size();
  std::printf("size: %g\n", 0.0 * size);
}

TEST(nuclear_radius) {
  Nucleus lead(list, 1);
  FUZZY_COMPARE(lead.default_nuclear_radius(), 1.2 * std::pow(208., 1. / 3.));
}

// check that center is at (0/0/0):
TEST(center) {
  constexpr int N_TEST = 20000;
  Nucleus lead(list, N_TEST);
  lead.set_nuclear_radius(lead.default_nuclear_radius());
  lead.arrange_nucleons();
  FourVector middle = lead.center();
  /** \f$\sqrt(\frac{\int_0^\infinity \frac{dr
   *  r^4}{\exp\left(\frac{x-R}{d}\right)+1}}{\int_0^\infinity \frac{dr
   *  r^2}{\exp\left(\frac{x-R}{d}\right)+1}\f$ with \f$R = 7.11 =
   *  1.2\sqrt[3]{208}\f$ and \f$d = 0.545\f$ is 5.86817. That is the
   *  standard deviation of our distribution. The distribution of
   *  centers is a gaussian with a width of \f$\sigma/\sqrt{N}\f$
   *  centered at 0.  I want the result to be within
   *  \f$3\frac{\sigma}{\sqrt{N}}\f$.
   *
   *  This code has been used at Wolfram Alpha to find the numerical
   *  values:
   *  \code
   *  sqrt(int( x**4/(exp((x-7.11)/.545)+1),x=0..infinity)/int(
   *  x**2/(exp((x-7.2)/.001)+1),x=0..infinity))
   *  \code
   **/
  double threesigma = 3 * 5.86817 / std::sqrt(N_TEST);
  VERIFY(std::abs(middle.x1()) < threesigma)
      << " x=" << middle.x1() << " vs. 3σ=" << threesigma
      << " (chance 1 in 370)";
  VERIFY(std::abs(middle.x2()) < threesigma)
      << " x=" << middle.x2() << " vs. 3σ=" << threesigma
      << " (chance 1 in 370)";
  VERIFY(std::abs(middle.x3()) < threesigma)
      << " x=" << middle.x3() << " vs. 3σ=" << threesigma
      << " (chance 1 in 370)";
}

TEST(center_hard_sphere) {
  constexpr int N_TEST = 20000;
  Nucleus lead(list, N_TEST);
  lead.set_diffusiveness(0.0);
  lead.set_nuclear_radius(lead.default_nuclear_radius());
  lead.arrange_nucleons();
  FourVector middle = lead.center();
  /**
   * Here, we can actually calculate the exact value for the width:
   * \f$\sigma = R\sqrt{\frac{3}{5}}\f$.
   **/
  double threesigma =
      3 * lead.default_nuclear_radius() * std::sqrt(0.6) / std::sqrt(N_TEST);
  VERIFY(std::abs(middle.x1()) < threesigma)
      << " x=" << middle.x1() << " vs. 3σ=" << threesigma
      << " (chance 1 in 370)";
  VERIFY(std::abs(middle.x2()) < threesigma)
      << " x=" << middle.x2() << " vs. 3σ=" << threesigma
      << " (chance 1 in 370)";
  VERIFY(std::abs(middle.x3()) < threesigma)
      << " x=" << middle.x3() << " vs. 3σ=" << threesigma
      << " (chance 1 in 370)";
}

// shift tests: here, in z direction the shift always depends on the
// maximum z value of all particles. Using many test particles and a
// hard sphere, I try to make this less random.
TEST(shift_zero) {
  constexpr int N_TEST = 1;
  Nucleus lead(list, N_TEST);
  lead.set_nuclear_radius(lead.default_nuclear_radius());
  lead.arrange_nucleons();
  FourVector precenter = lead.center();
  // shift with zero displacement: shouldn't change x and y, but note
  // that the z-parameter is the distance between the outer edges of the
  // nucleus!
  lead.shift(0, 0.0, 0.0);
  FourVector postcenter = lead.center();
  UnitTest::setFuzzyness<double>(30);
  FUZZY_COMPARE(postcenter.x1(), precenter.x1());
  FUZZY_COMPARE(postcenter.x2(), precenter.x2());
  FUZZY_COMPARE(postcenter.x3(), precenter.x3());
}

TEST(shift_x) {
  constexpr int N_TEST = 1;
  Nucleus lead(list, N_TEST);
  lead.set_nuclear_radius(lead.default_nuclear_radius());
  lead.arrange_nucleons();
  FourVector precenter = lead.center();
  // shift only in x.
  lead.shift(0, 4.0, 0.0);
  FourVector postcenter = lead.center();
  UnitTest::setFuzzyness<double>(30);
  FUZZY_COMPARE(postcenter.x1(), precenter.x1() + 4.0);
  FUZZY_COMPARE(postcenter.x2(), precenter.x2());
  FUZZY_COMPARE(postcenter.x3(), precenter.x3());
}

TEST(shift_z) {
  constexpr int N_TEST = 1;
  Nucleus lead(list, N_TEST);
  lead.set_nuclear_radius(lead.default_nuclear_radius());
  lead.arrange_nucleons();
  FourVector precenter = lead.center();
  lead.shift(4.0, 0.0, 0.0);
  FourVector postcenter = lead.center();
  UnitTest::setFuzzyness<double>(30);
  FUZZY_COMPARE(postcenter.x1(), precenter.x1());
  FUZZY_COMPARE(postcenter.x2(), precenter.x2());
  FUZZY_COMPARE(postcenter.x3(), precenter.x3() + 4.0);
}

// test the woods-saxon distribution at various discrete points:
TEST(woods_saxon) {
  // this is where we store the distribution.
  std::map<int, int> histogram{};
  // binning width for the distribution:
  constexpr double dx = 0.01;
  // the nucleus. Fill it from list with 1 testparticle.
  Nucleus projectile(list, 1);
  double R = projectile.default_nuclear_radius();
  projectile.set_nuclear_radius(R);
  // this is the number of times we access the distribution.
  constexpr int N_TEST = 10000000;
  // fill the histogram
  for (int i = 0; i < N_TEST; i++) {
    ThreeVector pos = projectile.distribute_nucleon();
    int bin = pos.abs() / dx;
    ++histogram[bin];
  }
  // We'll compare to relative values (I don't know what the integral
  // is)
  double value_at_radius = histogram.at(R / dx);
  double expected_at_radius = projectile.woods_saxon(R);
  // we'll probe at these values:
  double probes[9] = {1.0,    5.0,     7.2,     8.0,    8.5,
                      .5 * R, 1.1 * R, 1.2 * R, 1.3 * R};
  // now do probe these values:
  for (int i = 0; i < 9; ++i) {
    // value we have simulated:
    double value = histogram.at(probes[i] / dx) / value_at_radius;
    // value we have expected:
    double expec = projectile.woods_saxon(probes[i]) / expected_at_radius;
    // standard error we expect the histogram to have is 1/sqrt(N); we
    // give 3 sigma "space".
    double margin = 3. / std::sqrt(value);
    VERIFY(std::abs(value - expec) < margin)
        << " x = " << probes[i] << ": simulated: " << value
        << " vs. calculated: " << expec << " (allowed distance: " << margin
        << ")";
  }
}

TEST(Fermi_motion) {
  std::map<PdgCode, int> myfunnylist = {{0x2212, 22},  // protons
                                        {0x2112, 35},  // neutrons
                                        {0x111, 5},    // pions
                                        {0x3122, 1},   // Lambda
                                        {0x13, 1}};    // muon
  Nucleus myfunnynucleus(myfunnylist, 1);
  COMPARE(myfunnynucleus.size(), 22u + 35u + 5u + 1u + 1u);
  // Set some arbitrary radius and diffusiveness
  myfunnynucleus.set_nuclear_radius(3.);
  myfunnynucleus.set_diffusiveness(0.5);
  // Arrange coordinate space, because Fermi momenta depend on positions
  myfunnynucleus.arrange_nucleons();
  myfunnynucleus.generate_fermi_momenta();

  /* Check that: (i) only protons and neutrons get Fermi momenta
   *             (ii) total momentum is 0
   *             (iii) particles are on the mass shell
   */
  Particles particles;
  myfunnynucleus.copy_particles(&particles);
  ThreeVector ptot = ThreeVector();
  for (const auto &p : particles) {
    const ThreeVector mom3 = p.momentum().threevec();
    ptot += mom3;
    if (p.pdgcode() != 0x2212 && p.pdgcode() != 0x2112) {
      COMPARE(mom3.x1(), 0.0);
      COMPARE(mom3.x2(), 0.0);
      COMPARE(mom3.x3(), 0.0);
    }
    UnitTest::setFuzzyness<double>(2);
    FUZZY_COMPARE(p.momentum().sqr(), p.pole_mass() * p.pole_mass());
  }
  COMPARE_ABSOLUTE_ERROR(ptot.x1(), 0.0, 1.0e-15) << ptot.x1();
  COMPARE_ABSOLUTE_ERROR(ptot.x2(), 0.0, 1.0e-15) << ptot.x2();
  COMPARE_ABSOLUTE_ERROR(ptot.x3(), 0.0, 1.0e-15) << ptot.x3();
}
