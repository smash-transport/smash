/*
    Copyright (C) 2009-2014 Matthias Kretz <kretz@kde.org>

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) version 3.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public License
    along with this library; see the file COPYING.LIB.  If not, write to
    the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
    Boston, MA 02110-1301, USA.

*/

#ifndef SRC_INCLUDE_TSC_H_
#define SRC_INCLUDE_TSC_H_

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

#endif  // SRC_INCLUDE_TSC_H_
