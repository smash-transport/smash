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
#include "../include/random.h"

using namespace Smash;

// print the distribution (for plotting)
constexpr bool PRINT = false;

constexpr int N_TEST = 10000000;
constexpr int sigmabins = 4;
// We will check that at least the following ratios of numbers of samples lie
// within the corresponding sigma environment.
constexpr double allowed[sigmabins] = {.682 * .99, .954 * .99, .997 * .99, 1.0};

template <typename Normalization, typename Chi, typename Analytical>
void test_distribution(double dx, Normalization norm, Chi get_chi,
                       Analytical get_analyt) {
  std::map<int, int> hist{};
  for (int i = 0; i < N_TEST; i++) {
    const double chi = get_chi();
    ++hist[chi / dx];
  }
  int diffbad[sigmabins] = {0};
  int total = 0;
  for (auto b : hist) {
    // the bin
    const double x = dx * b.first;
    // this comes from the integral of the distribution inside the bins.
    const double normalize = norm(x);
    // the result that I got:
    const int result = b.second;
    // statistical error is sqrt(N)
    const double stat_err = sqrt(result);
    // the analytical value (stuff that is independent of x is only
    // calculated once before the loop)
    const double analyt = get_analyt(x) * normalize;
    // if we want to print the distribution, print!
    if (PRINT) {
      printf("%7.3f %7d %7.3f %7.3f\n", x, result, stat_err, analyt);
    }
    // normalize the deviation
    int diffbin = std::abs(result - analyt) / stat_err;
    diffbin = (diffbin >= sigmabins) ? sigmabins - 1 : diffbin;
    // this will be checked: how many points have been correct up to
    // one, two, three, ... sigmas?
    ++diffbad[diffbin];
    ++total;
  }
  // the integral (in each loop, this holds how many particles have a
  // normalized deviation less than this).
  int totalbad = 0;
  for (int unit = 0; unit < sigmabins - 1; ++unit) {
    totalbad += diffbad[unit];
    const double fraction = (totalbad + 0.0) / (total + 0.0);
    VERIFY(fraction > allowed[unit])
        << " too few entries have less than " << unit + 1
        << " sigma deviation (" << totalbad << "/" << total << "=" << fraction
        << ", required minimal fraction: " << allowed[unit];
  }
}

TEST(exponential) {
  constexpr double dx = 0.1;
  test_distribution(dx,
                    [](double) { return (1 - exp(-dx)) * N_TEST; },
                    []() { return Random::exponential(1.0); },
                    [](double x) { return exp(-x); });
}

TEST(x_exponential) {
  constexpr double dx = 0.05;
  test_distribution(
      dx,
      [](double x) {
        const double normalize_1 = (1 - exp(-dx)) * N_TEST;
        const double normalize_2 = (1 + dx) * normalize_1 - dx * N_TEST;
        return x * normalize_1 + normalize_2;
      },
      []() { return Random::exponential(1.0) + Random::exponential(1.0); },
      [](double x) { return exp(-x); });
}

TEST(xsquared_exponential) {
  constexpr double dx = 0.1;
  test_distribution(
      dx,
      [](double x) {
        const double normalize_1 = (1 - exp(-dx)) * N_TEST / 2.0;
        const double normalize_2 = 2.0 * (dx + 1.0);
        const double normalize_3 =
            N_TEST * (1.0 - exp(-dx) * (dx * (dx / 2.0 + 1.0) + 1.0));
        return (normalize_1 * (x * x + x * normalize_2) - x * dx * N_TEST +
                normalize_3);
      },
      []() {
        return Random::exponential(1.0) + Random::exponential(1.0) +
               Random::exponential(1.0);
      },
      [](double x) { return exp(-x); });
}

TEST(canonical) {
  constexpr double dx = 0.0001;
  test_distribution(dx,
                    [](double) { return N_TEST * dx; },
                    []() { return Random::canonical(); },
                    [](double) { return 1.0; });
}

TEST(uniform) {
  constexpr double dx = 0.01;
  auto random_4_6 = Random::make_uniform_distribution(-4.0, 6.0);
  test_distribution(dx,
                    [](double) { return N_TEST * dx / 10.0; },
                    [&]() { return random_4_6(); },
                    [](double) { return 1.0; });
}
