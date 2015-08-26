/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "histogram.h"

#include "../include/random.h"

using namespace Smash;

/**
 * Compare a random distribution to an analytical function.
 *
 * The distribution 'get_chi' is sampled 'n_test' times and put into a histogram
 * with bin size 'dx'. Then it is compared to the analytical function
 * 'get_analyt', multiplied by the normalization function 'norm'.
 */
template <typename Normalization, typename Chi, typename Analytical>
void test_distribution(int n_test, double dx, Normalization norm, Chi get_chi,
                       Analytical get_analyt) {
  Histogram1d hist(dx);
  // sample distribution and populate histogram
  hist.populate(n_test, get_chi);
  // test with the analytical function
  hist.test(norm, get_analyt);
}


constexpr int N_TEST = 1E7;

TEST(exponential) {
  constexpr double dx = 0.1;
  test_distribution(N_TEST, dx,
                    [](double) { return (1 - exp(-dx)); },
                    []() { return Random::exponential(1.0); },
                    [](double x) { return exp(-x); });
}

TEST(x_exponential) {
  constexpr double dx = 0.05;
  test_distribution(
      N_TEST, dx,
      [](double x) {
        const double normalize_1 = (1 - exp(-dx));
        const double normalize_2 = (1 + dx) * normalize_1 - dx;
        return x * normalize_1 + normalize_2;
      },
      []() { return Random::exponential(1.0) + Random::exponential(1.0); },
      [](double x) { return exp(-x); });
}

TEST(xsquared_exponential) {
  constexpr double dx = 0.1;
  test_distribution(
      N_TEST, dx,
      [](double x) {
        const double normalize_1 = (1 - exp(-dx)) / 2.0;
        const double normalize_2 = 2.0 * (dx + 1.0);
        const double normalize_3 =
            (1.0 - exp(-dx) * (dx * (dx / 2.0 + 1.0) + 1.0));
        return (normalize_1 * (x * x + x * normalize_2) - x * dx + normalize_3);
      },
      []() {
        return Random::exponential(1.0) + Random::exponential(1.0) +
               Random::exponential(1.0);
      },
      [](double x) { return exp(-x); });
}

TEST(canonical) {
  constexpr double dx = 0.0001;
  test_distribution(N_TEST, dx,
                    [](double) { return dx; },
                    []() { return Random::canonical(); },
                    [](double) { return 1.0; });
}

TEST(uniform) {
  constexpr double dx = 0.01;
  auto random_4_6 = Random::make_uniform_distribution(-4.0, 6.0);
  test_distribution(N_TEST, dx,
                    [](double) { return dx / 10.0; },
                    [&]() { return random_4_6(); },
                    [](double) { return 1.0; });
}
