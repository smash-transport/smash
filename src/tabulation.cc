/*
 *    Copyright (c) 2015-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/tabulation.h"

namespace smash {

Tabulation::Tabulation(double x_min, double range, size_t num,
                       std::function<double(double)> f)
    : x_min_(x_min), x_max_(x_min + range), inv_dx_(num / range) {
  if (num < 2) {
    throw std::runtime_error("Tabulation needs at least two values");
  }
  values_.resize(num + 1);
  const double dx = range / num;
  for (size_t i = 0; i <= num; i++) {
    values_[i] = f(x_min_ + i * dx);
  }
}

double Tabulation::get_value_step(double x) const {
  if (x < x_min_) {
    return 0.;
  }
  // this rounds correctly because double -> int conversion truncates
  const unsigned int n = (x - x_min_) * inv_dx_ + 0.5;
  if (n >= values_.size()) {
    return values_.back();
  } else {
    return values_[n];
  }
}

double Tabulation::get_value_linear(double x, Extrapolation extrapol) const {
  if (x < x_min_) {
    return 0.;
  }
  if (extrapol == Extrapolation::Zero && x > x_max_) {
    return 0.0;
  }
  if (extrapol == Extrapolation::Const && x > x_max_) {
    return values_.back();
  }
  const double index_double = (x - x_min_) * inv_dx_;
  // here n is the lower index
  const size_t n =
      std::min(static_cast<size_t>(index_double), values_.size() - 2);
  const double r = index_double - n;
  return values_[n] + (values_[n + 1] - values_[n]) * r;
}

/**
 * Write binary representation to stream.
 *
 * \param stream Output stream.
 * \param x Value to be written.
 */
static void swrite(std::ofstream& stream, double x) {
  stream.write(reinterpret_cast<const char*>(&x), sizeof(x));
}

/**
 * Read binary representation of a double.
 *
 * \param Input stream.
 * \return Read value.
 */
static double sread_double(std::ifstream& stream) {
  double x;
  stream.read(reinterpret_cast<char*>(&x), sizeof(x));
  return x;
}

/**
 * Write binary representation to stream.
 *
 * \param stream Output stream.
 * \param x Value to be written.
 */
static void swrite(std::ofstream& stream, size_t x) {
  // We want to support 32-bit and 64-bit platforms, so we store a 64-bit
  // integer on all platforms.
  const auto const_size_x = static_cast<uint64_t>(x);
  stream.write(reinterpret_cast<const char*>(&const_size_x),
               sizeof(const_size_x));
}

/**
 * Read binary representation of a size_t.
 *
 * \param Input stream.
 * \return Read value.
 */
static size_t sread_size(std::ifstream& stream) {
  uint64_t x;
  stream.read(reinterpret_cast<char*>(&x), sizeof(x));
  if (x > std::numeric_limits<size_t>::max()) {
    throw std::runtime_error("trying to read vector larger than supported");
  }
  return x;
}

/**
 * Write binary representation to stream.
 *
 * \param stream Output stream.
 * \param x Value to be written.
 */
static void swrite(std::ofstream& stream, const std::vector<double> x) {
  swrite(stream, x.size());
  if (x.size() > 0) {
    stream.write(reinterpret_cast<const char*>(x.data()),
                 sizeof(x[0]) * x.size());
  }
}

/**
 * Read binary representation of a vector of doubles.
 *
 * \param Input stream.
 * \return Read value.
 */
static std::vector<double> sread_vector(std::ifstream& stream) {
  const size_t n = sread_size(stream);
  std::vector<double> x;
  x.resize(n);
  stream.read(reinterpret_cast<char*>(x.data()), sizeof(double) * n);
  return x;
}

/**
 * Write binary representation to stream.
 *
 * \param stream Output stream.
 * \param x Value to be written.
 */
static void swrite(std::ofstream& stream, sha256::Hash x) {
  // The size is always the same, so there is no need to write it.
  stream.write(reinterpret_cast<const char*>(x.data()),
               sizeof(x[0]) * x.size());
}

/**
 * Read binary representation of a SHA256 hash.
 *
 * \param Input stream.
 * \return Read value.
 */
static sha256::Hash sread_hash(std::ifstream& stream) {
  sha256::Hash x;
  stream.read(reinterpret_cast<char*>(x.data()), x.size());
  return x;
}

void Tabulation::write(std::ofstream& stream, sha256::Hash hash) const {
  swrite(stream, hash);
  swrite(stream, x_min_);
  swrite(stream, x_max_);
  swrite(stream, inv_dx_);
  swrite(stream, values_);
}

Tabulation Tabulation::from_file(std::ifstream& stream, sha256::Hash hash) {
  sha256::Hash hash_from_stream = sread_hash(stream);
  Tabulation t;
  if (hash != hash_from_stream) {
    return t;
  }
  t.x_min_ = sread_double(stream);
  t.x_max_ = sread_double(stream);
  t.inv_dx_ = sread_double(stream);
  t.values_ = sread_vector(stream);
  return t;
}

}  // namespace smash
