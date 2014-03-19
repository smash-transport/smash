/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_RANDOM_H_
#define SRC_INCLUDE_RANDOM_H_

#include<random>

namespace Smash {
namespace Random {

using Engine = std::ranlux48;
extern /*thread_local*/ Engine engine;

template <typename T> class uniform_dist {
 public:
  uniform_dist(T min, T max)
    : distribution(min, max) {
  }
  T operator ()() {
    return distribution(engine);
  }
  std::uniform_real_distribution<T> distribution;
};

template <typename T> void set_seed(T &&seed) {
  engine.seed(std::forward<T>(seed));
}
template <typename T> T uniform(T min, T max) {
  return std::uniform_real_distribution<T>(min, max)(engine);
}
template <typename T = double> T canonical() {
  static uniform_dist<T> canonical(0.0, 1.0);
  return canonical();
}
template <typename T>
uniform_dist<T> make_uniform_distribution(T min, T max) {
  return uniform_dist<T>(min, max);
}
template <typename T = double> T exponential() {
  return std::exponential_distribution<T>(1)(engine);
}

}  // namespace Random
}  // namespace Smash

#endif  // SRC_INCLUDE_RANDOM_H_
