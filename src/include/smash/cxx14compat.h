/*
 *
 *    Copyright (c) 2014-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_CXX14COMPAT_H_
#define SRC_INCLUDE_SMASH_CXX14COMPAT_H_

#include <memory>
#include <utility>

namespace smash {

/**
 * Definition for make_unique
 * Is in C++14's standard library; necessary for older compilers
 *
 * \see http://en.cppreference.com/w/cpp/memory/unique_ptr/make_unique
 */
#if __cplusplus < 201402L
template <typename T, typename... Args>
inline std::unique_ptr<T> make_unique(Args &&... args) {
  return std::unique_ptr<T>{new T{std::forward<Args>(args)...}};
}
#else
  using std::make_unique;
#endif
}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CXX14COMPAT_H_
