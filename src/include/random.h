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

using RandomEngine = std::ranlux48;
extern /*thread_local*/ RandomEngine random_engine;

template <typename T> class uniform_dist {
 public:
  uniform_dist(T min, T max)
    : distribution(min, max) {
  }
  T operator ()() {
    return distribution(random_engine);
  }
  std::uniform_real_distribution<T> distribution;
};

template <typename T> void set_random_seed(T seed) {
  random_engine.seed(seed);
}
template <typename T> T random_uniform(T min, T max) {
  return std::uniform_real_distribution<T>(min, max)(random_engine);
}
template <typename T>
uniform_dist<T> make_uniform_distribution(T min, T max) {
  return uniform_dist<T>(min, max);
}
template <typename T> T random_exponential() {
  return std::exponential_distribution<T>(1)(random_engine);
}

}  // namespace Smash

#endif  // SRC_INCLUDE_RANDOM_H_
