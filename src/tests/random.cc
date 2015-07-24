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

// we might want to print the distribution (for plotting) instead of
// testing it!
constexpr bool PRINT = true;

constexpr int N_TEST = 10000000;
constexpr int sigmabins = 4;
constexpr double allowed[sigmabins] = {.682*.99, .954*.99, .997*.99, 1.0};

// all tests work similarly, only the first is extensively commented.
TEST(exponential) {
  std::map<int, int> hist {};
  constexpr float dx = 0.10;
  for (int i = 0; i < N_TEST; i++) {
    double chi = Random::exponential(1.0);
    ++hist[chi/dx];
  }
  int diffbad[sigmabins] = {0};
  int total = 0;
  // this comes from the integral of the exponential inside the bins.
  double normalize = (1 - exp(-dx)) * N_TEST;
  for (auto b : hist) {
    // the result that I got:
    int result = b.second;
    // statistical error is sqrt(N)
    int stat_err = sqrt(b.second);
    // the analytical value (stuff that is independent of x is only
    // calculated once before the loop)
    int analyt = exp(- dx * b.first) * normalize;
    // if we want to print the distribution, print!
    if (PRINT) {
      printf("%7.3f %d %d %d\n", dx * b.first, result, stat_err, analyt);
    }
    // normalize the deviation
    int diffbin = std::abs(result - analyt)/stat_err;
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
    double fraction = (totalbad + 0.0) / (total + 0.0);
    VERIFY(fraction > allowed[unit]) << " too few entries have less than "
              << unit+1 << " sigma deviation ("
              << totalbad << "/" << total << "=" << fraction
              << ", required minimal fraction: " << allowed[unit];
  }
}

TEST(x_exponential) {
  std::map<int, int> hist {};
  constexpr float dx = 0.10;
  for (int i = 0; i < N_TEST; i++) {
    double chi = Random::exponential(1.0);
    chi += Random::exponential(1.0);
    ++hist[chi/dx];
  }
  int diffbad[sigmabins] = {0};
  int total = 0;
  double normalize_1 = (1 - exp(-dx)) * N_TEST;
  double normalize_2 = (1 + dx) * normalize_1 - dx * N_TEST;

  for (auto b : hist) {
    int result = b.second;
    int stat_err = sqrt(b.second);
    double x = dx * b.first;
    double analyt = exp(-x) * (x * normalize_1 + normalize_2);
    if (PRINT) {
      printf("%7.3f %d %d %13.6e\n", dx * b.first, result, stat_err, analyt);
    }
    int diffbin = std::abs(result - analyt)/stat_err;
    diffbin = (diffbin >= sigmabins) ? sigmabins - 1 : diffbin;
    ++diffbad[diffbin];
    ++total;
  }
  int totalbad = 0;
  for (int unit = 0; unit < sigmabins - 1; ++unit) {
    totalbad += diffbad[unit];
    double fraction = (totalbad + 0.0) / (total + 0.0);
    VERIFY(fraction > allowed[unit]) << " too few entries have less than "
              << unit+1 << " sigma deviation ("
              << totalbad << "/" << total << "=" << fraction
              << ", required minimal fraction: " << allowed[unit];
  }
}

TEST(xsquared_exponential) {
  std::map<int, int> hist {};
  constexpr float dx = 0.10;
  for (int i = 0; i < N_TEST; i++) {
    double chi  = Random::exponential(1.0);
           chi += Random::exponential(1.0);
           chi += Random::exponential(1.0);
    ++hist[chi/dx];
  }
  int diffbad[sigmabins] = {0};
  int total = 0;
  // this comes from the integral of the exponential inside the bins.
  double normalize_1 = (1 - exp(-dx)) * N_TEST / 2.0;
  double normalize_2 = 2.0 * (dx + 1.0);
  double normalize_3 = N_TEST * (1.0 - exp(-dx) * ( dx * (dx / 2.0 + 1.0) + 1.0));

  for (auto b : hist) {
    int result = b.second;
    int stat_err = sqrt(b.second);
    double x = dx * b.first;
    double analyt = exp(-x) * (normalize_1 * (x * x + x * normalize_2)
                              - x * dx * N_TEST + normalize_3);
    if (PRINT) {
      printf("%7.3f %d %d %13.6e\n", dx * b.first, result, stat_err, analyt);
    }
    int diffbin = std::abs(result - analyt)/stat_err;
    diffbin = (diffbin >= sigmabins) ? sigmabins - 1 : diffbin;
    ++diffbad[diffbin];
    ++total;
  }
  int totalbad = 0;
  for (int unit = 0; unit < sigmabins - 1; ++unit) {
    totalbad += diffbad[unit];
    double fraction = (totalbad + 0.0) / (total + 0.0);
    VERIFY(fraction > allowed[unit]) << " too few entries have less than "
              << unit+1 << " sigma deviation ("
              << totalbad << "/" << total << "=" << fraction
              << ", required minimal fraction: " << allowed[unit];
  }
}

TEST(canonical) {
  std::map<int, int> hist {};
  constexpr float dx = 0.001;
  for (int i = 0; i < N_TEST; i++) {
    double chi = Random::canonical();
    ++hist[chi/dx];
  }
  int diffbad[sigmabins] = {0};
  int total = 0;
  // here, all bins are equally likely.
  constexpr double analyt = N_TEST * dx;
  for (auto b : hist) {
    int result = b.second;
    int stat_err = sqrt(b.second);
    if (PRINT) {
      printf("%7.3f %d %d %13.6e\n", dx * b.first, result, stat_err, analyt);
    }
    int diffbin = std::abs(result - analyt)/stat_err;
    diffbin = (diffbin >= sigmabins) ? sigmabins - 1 : diffbin;
    ++diffbad[diffbin];
    ++total;
  }
  int totalbad = 0;
  for (int unit = 0; unit < sigmabins - 1; ++unit) {
    totalbad += diffbad[unit];
    double fraction = (totalbad + 0.0) / (total + 0.0);
    VERIFY(fraction > allowed[unit]) << " too few entries have less than "
              << unit+1 << " sigma deviation ("
              << totalbad << "/" << total << "=" << fraction
              << ", required minimal fraction: " << allowed[unit];
  }
}

TEST(uniform) {
  std::map<int, int> hist {};
  constexpr float dx = 0.01;
  auto random_4_6 = Random::make_uniform_distribution(-4.0,6.0);
  for (int i = 0; i < N_TEST; i++) {
    double chi = random_4_6();
    ++hist[floor(chi/dx)];
  }
  int diffbad[sigmabins] = {0};
  int total = 0;
  // here, all bins are equally likely.
  constexpr double analyt = N_TEST * dx / 10.0;
  for (auto b : hist) {
    int result = b.second;
    int stat_err = sqrt(b.second);
    if (PRINT) {
      printf("%7.3f %d %d %13.6e\n", dx * b.first, result, stat_err, analyt);
    }
    int diffbin = std::abs(result - analyt)/stat_err;
    diffbin = (diffbin >= sigmabins) ? sigmabins - 1 : diffbin;
    ++diffbad[diffbin];
    ++total;
  }
  int totalbad = 0;
  for (int unit = 0; unit < sigmabins - 1; ++unit) {
    totalbad += diffbad[unit];
    double fraction = (totalbad + 0.0) / (total + 0.0);
    VERIFY(fraction > allowed[unit]) << " too few entries have less than "
              << unit+1 << " sigma deviation ("
              << totalbad << "/" << total << "=" << fraction
              << ", required minimal fraction: " << allowed[unit];
  }
}

TEST(cos_like) {
  std::map<int, int> hist {};
  constexpr float dx = 0.01;
  auto cos_like = Random::make_uniform_distribution(-1.0, +1.0);
  for (int i = 0; i < N_TEST; i++) {
    double chi = cos_like();
    ++hist[floor(chi/dx)];
  }
  int diffbad[sigmabins] = {0};
  int total = 0;
  // here, all bins are equally likely.
  constexpr double analyt = N_TEST * dx / 2.0;
  for (auto b : hist) {
    int result = b.second;
    int stat_err = sqrt(b.second);
    if (PRINT) {
      printf("%7.3f %d %d %13.6e\n", dx * b.first, result, stat_err, analyt);
    }
    int diffbin = std::abs(result - analyt)/stat_err;
    diffbin = (diffbin >= sigmabins) ? sigmabins - 1 : diffbin;
    ++diffbad[diffbin];
    ++total;
  }
  int totalbad = 0;
  for (int unit = 0; unit < sigmabins - 1; ++unit) {
    totalbad += diffbad[unit];
    double fraction = (totalbad + 0.0) / (total + 0.0);
    VERIFY(fraction > allowed[unit]) << " too few entries have less than "
              << unit+1 << " sigma deviation ("
              << totalbad << "/" << total << "=" << fraction
              << ", required minimal fraction: " << allowed[unit];
  }
}

TEST(phi_like) {
  std::map<int, int> hist {};
  constexpr float dx = 2.0 * M_PI / 10000.0;
  auto phi_like = Random::make_uniform_distribution(0.0, 2 * M_PI);
  for (int i = 0; i < N_TEST; i++) {
    double chi = phi_like();
    ++hist[floor(chi/dx)];
  }
  int diffbad[sigmabins] = {0};
  int total = 0;
  // here, all bins are equally likely.
  constexpr double analyt = N_TEST * dx / (2.0 * M_PI);
  for (auto b : hist) {
    int result = b.second;
    int stat_err = sqrt(b.second);
    if (PRINT) {
      printf("%7.3f %d %d %13.6e\n", dx * b.first, result, stat_err, analyt);
    }
    int diffbin = std::abs(result - analyt)/stat_err;
    diffbin = (diffbin >= sigmabins) ? sigmabins - 1 : diffbin;
    ++diffbad[diffbin];
    ++total;
  }
  int totalbad = 0;
  for (int unit = 0; unit < sigmabins - 1; ++unit) {
    totalbad += diffbad[unit];
    double fraction = (totalbad + 0.0) / (total + 0.0);
    VERIFY(fraction > allowed[unit]) << " too few entries have less than "
              << unit+1 << " sigma deviation ("
              << totalbad << "/" << total << "=" << fraction
              << ", required minimal fraction: " << allowed[unit];
  }
}
