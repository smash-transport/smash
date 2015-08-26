/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <cstdio>
#include <map>
#include <string>

namespace Smash {

/** A class for representing a one-dimensional histogram. */
class Histogram1d {
 public:
  // construct empty histogram with a given bin size
  Histogram1d(double d) : dx_(d), data_({}) {}
  double dx () const { return dx_; }
  int num_entries() const { return entries_; }
  // add an entry at position x
  void add(double x) {
    int n = x / dx_;
    ++data_[n];
    ++entries_;
  }
  // populate with 'n_test' entries from distribution 'chi'
  template <typename Chi>
  void populate(int n_test, Chi chi) {
    for (int i = 0; i < n_test; i++) {
      add(chi());
    }
  }
  // print the histogram to a file
  void print_to_file(std::string fname) const;
  // compare a histogram of a sampled distribution to an analytical function
  template <typename Normalization, typename Analytical>
  void test(Normalization norm, Analytical get_analyt) const;
 private:
  double dx_;                // bin size
  std::map<int, int> data_;  // data storage
  int entries_ = 0;          // count number of entries;
};


/* Print the histogram to a file. */
void Histogram1d::print_to_file(std::string fname) const {
  std::FILE *file = std::fopen(fname.c_str(), "w");
  for (auto b : data_) {
    const double m = dx() * b.first;         // mass bin
    const int result = b.second;             // number of counts
    const double stat_err = sqrt(result);    // statistical error is sqrt(N)
    fprintf(file, "%7.3f %7d %7.3f\n", m, result, stat_err);
  }
  std::fclose(file);
}


/**
 * Compare a histogram of a sampled distribution to an analytical function.
 *
 * The histogram 'hist' (populated with sampled distribution values) is
 * compared to the analytical function 'get_analyt', multiplied by the
 * normalization function 'norm'.
 */
template <typename Normalization, typename Analytical>
void Histogram1d::test(Normalization norm, Analytical get_analyt) const {
  // print the distribution (for plotting)
  constexpr bool PRINT = false;

  // We will check that at least the following ratios of numbers of samples lie
  // within the corresponding sigma environment.
  constexpr int sigmabins = 4;
  constexpr double allowed[sigmabins] = {.682 * .99, .954 * .99, .997 * .99, 1.0};

  int diffbad[sigmabins] = {0};
  int total = 0;

  for (const auto b : data_) {
    // the bin
    const double x = dx() * b.first;
    // the result that I got:
    const int result = b.second;
    // statistical error is sqrt(N)
    const double stat_err = sqrt(result);
    // the analytical value times normalization times number of samples
    const double analyt = get_analyt(x) * norm(x) * num_entries();
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

}  // namespace Smash
