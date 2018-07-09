/*
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_AVERAGE_H_
#define SRC_INCLUDE_AVERAGE_H_

#include <cassert>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace smash {

/** Calculate an average value incrementally.
 *
 * \tparam T Type of the values (should be floating point).
 */
template <typename T>
class Average {
 public:
  /// Create a new object to calculate an average.
  Average() : avg_(0), n_(0) {}

  /// Add a value \p x to the set of numbers defining the average.
  void add(T x) {
    avg_ = (avg_ * n_ + x) / (n_ + 1);
    n_++;
  }

  /// \return Current average.
  T average() const { return avg_; }

  /// \return Number of values used to calculate the average.
  uint64_t number_of_values() const { return n_; }

  /// Reset the average to 0.
  void clear() {
    avg_ = 0;
    n_ = 0;
  }

 private:
  /// Average.
  T avg_;
  /// Sample size.
  uint64_t n_;
};

/** Remove duplicates from data (x, y) by averaging y.
 *
 * Assumes (x, y) is sorted.
 *
 * \tparam T Type of the values (should be floating point).
 * \param x x-values.
 * \param y y-values.
 * \return New x and y values as a pair of vectors.
 */
template <typename T>
std::pair<std::vector<T>, std::vector<T>> dedup_avg(const std::vector<T>& x,
                                                    const std::vector<T>& y) {
  if (x.size() != y.size()) {
    std::stringstream ss;
    ss << "x and y have to be of same size: " << x.size() << " != " << y.size();
    throw std::runtime_error(ss.str());
  }
  if (x.size() < 1) {
    throw std::runtime_error("x cannot be empty.");
  }
  std::vector<T> new_x;
  new_x.reserve(x.size());
  std::vector<T> new_y;
  new_y.reserve(y.size());
  Average<T> avg;
  T prev_x = x[0];
  for (size_t i = 0; i < x.size(); i++) {
    if (x[i] == prev_x) {
      avg.add(y[i]);
    } else {
      assert(i != 0);
      new_x.push_back(x[i - 1]);
      new_y.push_back(avg.average());
      avg.clear();
      avg.add(y[i]);
      prev_x = x[i];
    }
  }
  new_x.push_back(x.back());
  new_y.push_back(avg.average());
  return std::make_pair(std::move(new_x), std::move(new_y));
}

}  // namespace smash

#endif  // SRC_INCLUDE_AVERAGE_H_
