/*
 *
 *    Copyright (C) 2009-2014 Matthias Kretz <kretz@kde.org>
 *    Copyright (c) 2014-2015,2017-2018,2020-2021,2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_TSC_H_
#define SRC_INCLUDE_SMASH_TSC_H_

#ifndef USE_NANOBENCHMARKING_CODE
#error "Nanobenchmark code tried to be used without enabling it via CMake."
#else

#include <cstdint>
#include <iosfwd>

#ifdef _MSC_VER
#include <intrin.h>
#pragma intrinsic(__rdtsc)
#endif

namespace smash {

/**
 * A low-overhead timer for benchmarking small regions of code.
 *
 * This is usually unused in SMASH, because it should only be temporarily
 * used when benchmarking.
 */
class TimeStampCounter {
 public:
  /// Start the counter.
  void start();

  /// Stop the counter.
  void stop();

  /// \return Number of counted CPU cycles.
  uint64_t cycles() const;

 private:
  /// Union that stores cycles, \todo Why this data type?
  union Data {
    /// Either one 64-bit integer
    uint64_t a;
    /// Or two 32-bit integers
    unsigned int b[2];
  };
  /// Stores start of benchmarking
  Data m_start;
  /// Stores end of benchmarking
  Data m_end;
};

inline void TimeStampCounter::start() {
#ifdef VC_IMPL_MIC
  asm volatile("xor %%eax,%%eax\n\tcpuid\n\trdtsc"
               : "=a"(m_start.b[0]), "=d"(m_start.b[1])::"ebx", "ecx");
#elif defined _MSC_VER
  unsigned int tmp;
  m_start.a = __rdtscp(&tmp);
#else
  asm volatile("rdtscp" : "=a"(m_start.b[0]), "=d"(m_start.b[1])::"ecx");
#endif
}

inline void TimeStampCounter::stop() {
#ifdef VC_IMPL_MIC
  asm volatile("xor %%eax,%%eax\n\tcpuid\n\trdtsc"
               : "=a"(m_end.b[0]), "=d"(m_end.b[1])::"ebx", "ecx");
#elif defined _MSC_VER
  unsigned int tmp;
  m_end.a = __rdtscp(&tmp);
#else
  asm volatile("rdtscp" : "=a"(m_end.b[0]), "=d"(m_end.b[1])::"ecx");
#endif
}

inline uint64_t TimeStampCounter::cycles() const { return m_end.a - m_start.a; }

std::ostream &operator<<(std::ostream &out, const TimeStampCounter &tsc);

}  // namespace smash

#endif  // USE_NANOBENCHMARKING_CODE
#endif  // SRC_INCLUDE_SMASH_TSC_H_
