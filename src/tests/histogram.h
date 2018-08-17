/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <map>
#include <string>

#include "../include/smash/file.h"
#include "../include/smash/integrate.h"

namespace smash {

/** A class for representing a one-dimensional histogram. */
class Histogram1d {
 public:
  /// construct empty histogram with a given bin size
  Histogram1d(double d) : dx_(d), data_({}) {}
  /// get the bin size
  double dx() const { return dx_; }

  /// get the total number of entries in the histogram
  int num_entries() const { return entries_; }
  /// get the minimum value observed in the histogram
  double get_min() const { return min_ * dx_; }
  /// get the maximum value observed in the histogram
  double get_max() const { return max_ * dx_; }

  /// add an entry at position x
  void add(double x) {
    int n = std::floor(x / dx_);
    ++data_[n];
    ++entries_;
    if (n < min_)
      min_ = n;
    if (n >= max_)
      max_ = n + 1;
  }

  /// populate with 'n_test' entries from distribution 'chi'
  template <typename Chi>
  void populate(int n_test, Chi chi) {
    for (int i = 0; i < n_test; i++) {
      add(chi());
    }
  }

  // print the histogram to a file
  void print_to_file(const std::string& fname) const;

  // compare a histogram of a sampled distribution to an analytical function
  template <typename Analytical>
  void test(Analytical analyt, std::string dbg_file = "") const;

 private:
  double dx_;                          // bin size
  std::map<int, int> data_;            // data storage
  int entries_ = 0;                    // number of entries
  int min_ = INT_MAX, max_ = INT_MIN;  // minimum & maxiumum values
};

/** Print the histogram to a file. */
void Histogram1d::print_to_file(const std::string& fname) const {
  FilePtr file = fopen(fname, "w");
  for (auto b : data_) {
    const double m = dx() * (b.first + 0.5);  // mass bin
    const int result = b.second;              // number of counts
    const double stat_err = sqrt(result);     // statistical error is sqrt(N)
    fprintf(file.get(), "%7.3f %7d %7.3f\n", m, result, stat_err);
  }
}

/**
 * Compare a histogram of a sampled distribution to an analytical function.
 *
 * The histogram 'hist' (populated with sampled distribution values) is
 * compared to the analytical function 'analyt' (whose normalization is
 * being determined automatically).
 */
template <typename Analytical>
void Histogram1d::test(Analytical analyt, std::string dbg_file) const {
  // We will check that at least the following ratios of numbers of samples lie
  // within the corresponding sigma environment.
  constexpr int sigmabins = 4;
  constexpr double tolerance = 0.90;
  constexpr double allowed[sigmabins] = {.682 * tolerance, .954 * tolerance,
                                         .997 * tolerance, 1.00};

  int diffbad[sigmabins] = {0};
  int total = 0;

  /* Determine the normalization factor for the analytical function,
   * by integrating it from the minimum to the maximum value of the histogram.
   */
  Integrator integrate;
  const double min = get_min();
  const double max = get_max();
  const double itg = integrate(min, max, [&](double x) { return analyt(x); });
  /* The normalization factor also accounts for the bin width
   * and the total number of emtries in the histogram. */
  const double norm = dx() * num_entries() / itg;
  printf("norm: %f %f %f %f\n", min, max, itg, norm);

  double chi_sqr = 0.;
  {
    FilePtr file = nullptr;
    if (dbg_file != "") {
      file = fopen(dbg_file, "w");
    }

    for (const auto b : data_) {
      // center of the bin
      const double x = (b.first + 0.5) * dx();
      // the number N of counts in that bin
      const int counts = b.second;
      // statistical error is sqrt(N)
      const double stat_err = sqrt(counts);
      // the analytical value times normalization times number of samples
      const double ana = analyt(x) * norm;
      // if we want to print the distribution, print!
      if (file) {
        fprintf(file.get(), "%7.3f %7d %7.3f %7.3f\n", x, counts, stat_err,
                ana);
      }
      const double diff = (counts - ana) / stat_err;  // normalized deviation
      chi_sqr += diff * diff;                         // chi-squared
      // absolute value of deviation in units of sigma
      int abs_diff = std::abs(diff);
      abs_diff = (abs_diff >= sigmabins) ? sigmabins - 1 : abs_diff;
      // this will be checked: how many points have been correct up to
      // one, two, three, ... sigmas?
      ++diffbad[abs_diff];
      ++total;
    }
  }

  // divide chi^2 by d.o.f.
  chi_sqr = chi_sqr / data_.size();
  printf("chi_sqr per d.o.f: %f\n", chi_sqr);
  VERIFY(chi_sqr < 1.4) << "Error: chi_squared is too large!";

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

}  // namespace smash
