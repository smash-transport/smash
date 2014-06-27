/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "unittest.h"

#include <map>
#include "../include/nucleus.h"
#include "../include/pdgcode.h"
#include "../include/lineparser.h"

#include <boost/filesystem/fstream.hpp>

using namespace Smash;

std::map<PdgCode, int> list = {{0x2212, 82}, {0x2112, 126}};

TEST(init_particle_types) {
  ParticleType::create_type_list(read_all(bf::ifstream{"particles.txt"}));
}

TEST(initialize_realparticles) {
  Nucleus lead;
  // fill with 208 nucleons:
  lead.fill_from_list(list, 1);
  COMPARE(lead.number_of_particles(), 208u);
  COMPARE(lead.size(), 208u);
}

TEST(initialize_testparticles) {
  Nucleus lead;
  constexpr int N_TEST = 10;
  lead.fill_from_list(list, N_TEST);
  COMPARE(lead.number_of_particles(), 208u);
  COMPARE(lead.size(), 208u * N_TEST);
}

TEST(initialize_testparticles_multiple) {
  Nucleus lead;
  constexpr int N_TEST = 10;
  lead.fill_from_list(list, N_TEST);
  lead.fill_from_list(list, N_TEST);
  COMPARE(lead.number_of_particles(), 416u);
  COMPARE(lead.size(), 416u * N_TEST);
}

TEST_CATCH(initialize_testparticles_wrong, Nucleus::TestparticleConfusion) {
  Nucleus lead;
  constexpr int N_TEST = 10;
  lead.fill_from_list(list, 1);
  lead.fill_from_list(list, N_TEST);
  COMPARE(lead.number_of_particles(), 208u);
  // this should throw an error: list size is (N_TEST+1)*208 = 2288, not
  // divisible by N_TEST (unless someone set N_TEST to 1).
  size_t size = lead.size();
  printf("size: %g\n", 0.0*size);
}

TEST(nuclear_radius) {
  Nucleus lead;
  lead.fill_from_list(list, 1);
  FUZZY_COMPARE(lead.nuclear_radius(), static_cast<float>(1.2f*pow(208,1./3.)));
}

// check that center is at (0/0/0):
TEST(center) {
  Nucleus lead;
  constexpr int N_TEST = 20000;
  lead.fill_from_list(list, N_TEST);
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
  double threesigma = 3*5.86817 / sqrt(N_TEST);
  VERIFY(std::abs(middle.x1()) < threesigma) << " x=" << middle.x1() << " vs. 3σ=" << threesigma << " (chance 1 in 370)";
  VERIFY(std::abs(middle.x2()) < threesigma) << " x=" << middle.x2() << " vs. 3σ=" << threesigma << " (chance 1 in 370)";
  VERIFY(std::abs(middle.x3()) < threesigma) << " x=" << middle.x3() << " vs. 3σ=" << threesigma << " (chance 1 in 370)";
}

TEST(center_hard_sphere) {
  Nucleus lead;
  constexpr int N_TEST = 20000;
  lead.fill_from_list(list, N_TEST);
  lead.set_diffusiveness(0.0);
  lead.arrange_nucleons();
  FourVector middle = lead.center();
  /**
   * Here, we can actually calculate the exact value for the width:
   * \f$\sigma = R\sqrt{\frac{3}{5}}\f$.
   **/
  double threesigma = 3* lead.nuclear_radius() * sqrt(0.6) / sqrt(N_TEST);
  VERIFY(std::abs(middle.x1()) < threesigma) << " x=" << middle.x1() << " vs. 3σ=" << threesigma << " (chance 1 in 370)";
  VERIFY(std::abs(middle.x2()) < threesigma) << " x=" << middle.x2() << " vs. 3σ=" << threesigma << " (chance 1 in 370)";
  VERIFY(std::abs(middle.x3()) < threesigma) << " x=" << middle.x3() << " vs. 3σ=" << threesigma << " (chance 1 in 370)";
}

// shift tests: here, in z direction the shift always depends on the
// maximum z value of all particles. Using many test particles and a
// hard sphere, I try to make this less random.
TEST(shift_zero) {
  Nucleus lead;
  constexpr int N_TEST = 1000;
  lead.fill_from_list(list, N_TEST);
  lead.set_diffusiveness(0.0);
  lead.arrange_nucleons();
  FourVector precenter = lead.center();
  // shift with zero displacement: shouldn't change x and y, but note
  // that the z-parameter is the distance between the outer edges of the
  // nucleus!
  lead.shift(true, 0, 0.0, 0.0);
  FourVector postcenter = lead.center();
  UnitTest::setFuzzyness<double>(10);
  FUZZY_COMPARE(postcenter.x1(), precenter.x1());
  FUZZY_COMPARE(postcenter.x2(), precenter.x2());
  COMPARE_RELATIVE_ERROR(postcenter.x3(), precenter.x3()-lead.nuclear_radius(),
                         0.2);
}

TEST(shift_x) {
  Nucleus lead;
  constexpr int N_TEST = 1000;
  lead.fill_from_list(list, N_TEST);
  lead.set_diffusiveness(0.0);
  lead.arrange_nucleons();
  FourVector precenter = lead.center();
  // shift only in x.
  lead.shift(true, 0, 4.0, 0.0);
  FourVector postcenter = lead.center();
  UnitTest::setFuzzyness<double>(150);
  FUZZY_COMPARE(postcenter.x1(), precenter.x1()+4.0);
  FUZZY_COMPARE(postcenter.x2(), precenter.x2());
  COMPARE_RELATIVE_ERROR(postcenter.x3(), precenter.x3()-lead.nuclear_radius(),
                        0.2);
}

TEST(shift_z) {
  Nucleus lead;
  constexpr int N_TEST = 1000;
  lead.fill_from_list(list, N_TEST);
  lead.set_diffusiveness(0.0);
  lead.arrange_nucleons();
  FourVector precenter = lead.center();
  // shift in z. Here, we cannot exactly predict what happens, because
  // the shift depends on the outermost particle. That's why we use many
  // test particles and a hard sphere.
  lead.shift(true, 4.0, 0.0, 0.0);
  FourVector postcenter = lead.center();
  UnitTest::setFuzzyness<double>(10);
  FUZZY_COMPARE(postcenter.x1(), precenter.x1());
  FUZZY_COMPARE(postcenter.x2(), precenter.x2());
  COMPARE_RELATIVE_ERROR(postcenter.x3(),
      precenter.x3() - lead.nuclear_radius() + 4.0, 0.2);
}

// test the woods-saxon distribution at various discrete points:
TEST(woods_saxon) {
  // this is where we store the distribution.
  std::map<int, int> histogram {};
  // binning width for the distribution:
  constexpr float dx = 0.01;
  // the nucleus. Fill it from list with 1 testparticle.
  Nucleus projectile;
  projectile.fill_from_list(list, 1);
  float R = projectile.nuclear_radius();
  // this is the number of times we access the distribution.
  constexpr int N_TEST = 10000000;
  // fill the histogram
  for (int i = 0; i < N_TEST; i++) {
    float radius = projectile.distribution_nucleons();
    int bin = radius/dx;
    ++histogram[bin];
  }
  // We'll compare to relative values (I don't know what the integral
  // is)
  float value_at_radius = histogram.at(R/dx);
  float expected_at_radius =
                           projectile.woods_saxon(R);
  // we'll probe at these values:
  float probes[9] = { 1.0, 5.0, 7.2, 8.0, 8.5, .5f*R, 1.1f*R, 1.2f*R, 1.3f*R };
  // now do probe these values:
  for (int i = 0; i < 9; ++i) {
    // value we have simulated:
    float value = histogram.at(probes[i]/dx)/value_at_radius;
    // value we have expected:
    float expec = projectile.woods_saxon(probes[i])/expected_at_radius;
    // standard error we expect the histogram to have is 1/sqrt(N); we
    // give 3 sigma "space".
    float margin = 3.0/sqrt(value);
    VERIFY(std::abs(value - expec) < margin) << " x = " << probes[i]
            << ": simulated: " << value
            << " vs. calculated: " << expec
            << " (allowed distance: " << margin << ")";
  }
}
