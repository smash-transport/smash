/*
 *    Copyright (c) 2013-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/clebschgordan.h"

#include <numeric>
#include <unordered_map>

#include "gsl/gsl_sf_coupling.h"

#include "smash/constants.h"
#include "smash/logging.h"

namespace {

/**
 * Auxiliary struct to be used as key in the look up table of Clebsch-Gordan
 * coefficients. It basically contains the input to retrieve one coefficient.
 */
struct ThreeSpins {
  int j1;  /// First spin
  int j2;  /// Second spin
  int j3;  /// Third spin
  int m1;  /// First isospin
  int m2;  /// Second isospin
  int m3;  /// Third isospin

  /**
   * Comparison operator between two set of spin information. This is needed in
   * order to use this object in a \c std::unordered_map container.
   *
   * @param other The object to be compared to
   * @return \c true If all 6 spins value are identical
   * @return \c false otherwise
   */
  bool operator==(const ThreeSpins &other) const {
    return std::tie(j1, j2, j3, m1, m2, m3) ==
           std::tie(other.j1, other.j2, other.j3, other.m1, other.m2, other.m3);
  }
};

/**
 * This is one of the possible ways to prepare a hashing mechanism to use a
 * custom object in a \c std::unordered_map container. It has been preferred
 * here to use a new \c struct instead of injecting a specialization into the
 * \c std namespace, because we are in an anonymous namespace and working in
 * one specific cpp file only. Since the hash is not trivial, we also preferred
 * this approach to using a lambda function to declare the hashing function.
 */
struct ThreeSpinHash {
  /**
   * The overload of the \c operator() is the only needed ingredient to make
   * this class ready to be used as hashing algorithm.
   *
   * Since there is not (yet) a C++ standard way of combining hashes, to
   * implement functions like this one, it is necessary to choose a hashing
   * algorithm. The algorithm should minimize the hash collision (different
   * input giving the same output), but at the same time be as fast as possible.
   * One could use boost library approach, which defines
   * \code {.cpp}
   * template <typename T>
   * inline void hash_combine(std::size_t &seed, const T &val) {
   *   seed ^= std::hash<T>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
   * }
   * \endcode
   * as function to combine hashes. Then this should be used here on each \c in
   * member to produce the final hash. Although it has been tested to work,
   * there is a more physics driven approach. In \cite Rasch2004 an algorithm
   * to efficiently store Clebsch-Gordan coefficient is discussed. In it a way
   * to map the three spins information to a single integer is proposed and this
   * can be used here as hashing function, basically offering the guarantee that
   * no hash collision will occur.
   *
   * @param in The spin information (meant to be used as input to be hashed)
   * @return \c std::size_t The calculated hash value
   */
  std::size_t operator()(const ThreeSpins &in) const noexcept {
    const int S = -in.j1 + in.j2 + in.j3;
    const int L = +in.j1 - in.j2 + in.j3;
    const int X = +in.j1 - in.m1;
    const int B = +in.j2 - in.m2;
    const int T = +in.j3 + in.m3;
    auto hash = L * (24 + L * (50 + L * (35 + L * (10 + L)))) / 120 +
                X * (6 + X * (11 + X * (6 + X))) / 24 +
                T * (2 + T * (3 + T)) / 6 + B * (B + 1) / 2 + S + 1;
    return hash;
  }
};

}  // namespace

namespace smash {
static constexpr int LResonances = LogArea::Resonances::id;

static double clebsch_gordan_calculation(const int j_a, const int j_b,
                                         const int j_c, const int m_a,
                                         const int m_b, const int m_c) {
  const double wigner_3j = gsl_sf_coupling_3j(j_a, j_b, j_c, m_a, m_b, -m_c);
  if (std::abs(wigner_3j) < really_small) {
    return 0.;
  }
  assert((j_a - j_b + m_c) % 2 == 0);
  const int j = (j_a - j_b + m_c) / 2;
  double result = std::sqrt(j_c + 1) * wigner_3j;
  result *= (j % 2 == 0) * 2 - 1;  // == (-1)**j

  logg[LResonances].debug("CG: ", result, " I1: ", j_a, " I2: ", j_b,
                          " IR: ", j_c, " iz1: ", m_a, " iz2: ", m_b,
                          " izR: ", m_c);

  return result;
}

double clebsch_gordan(const int j_a, const int j_b, const int j_c,
                      const int m_a, const int m_b, const int m_c) {
  /* Use a standard unordered map here as look up container. Unordered map is an
   associative container where search, insertion, and removal of elements have
   average constant-time complexity (worst case linear in the size of the
   container). */
  static std::unordered_map<ThreeSpins, double, ThreeSpinHash>
      coefficients_look_up_table{};
  ThreeSpins spin_information = {j_a, j_b, j_c, m_a, m_b, m_c};
  if (auto search = coefficients_look_up_table.find(spin_information);
      search != coefficients_look_up_table.end()) {
    return search->second;
  } else {
    double result = clebsch_gordan_calculation(j_a, j_b, j_c, m_a, m_b, m_c);
    coefficients_look_up_table[spin_information] = result;
    return result;
  }
}

/**
 * Calculate isospin Clebsch-Gordan coefficient for two particles p_a and p_b
 * coupling to a total isospin \see clebsch_gordan for details (I_tot, I_z).
 * \param[in] p_a Information of particle type for first particle
 * \param[in] p_b Information of particle type for second particle
 * \param[out] I_tot Total isospin of the reaction
 * \param[out] I_z Total isospin 3 component of the reaction
 */
static double isospin_clebsch_gordan_2to1(const ParticleType &p_a,
                                          const ParticleType &p_b,
                                          const int I_tot, const int I_z) {
  return clebsch_gordan(p_a.isospin(), p_b.isospin(), I_tot, p_a.isospin3(),
                        p_b.isospin3(), I_z);
}

double isospin_clebsch_gordan_sqr_3to1(const ParticleType &p_a,
                                       const ParticleType &p_b,
                                       const ParticleType &p_c,
                                       const ParticleType &Res) {
  // Calculate allowed isospin range for 3->1 reaction I_ab
  const auto min_I_ab = std::abs(p_a.isospin() - p_b.isospin());
  const auto max_I_ab = p_a.isospin() + p_b.isospin();
  std::vector<int> possible_I_ab(max_I_ab - min_I_ab + 1);
  std::iota(possible_I_ab.begin(), possible_I_ab.end(), min_I_ab);
  std::vector<int> allowed_I_ab;
  allowed_I_ab.reserve(possible_I_ab.size());
  for (const auto Iab : possible_I_ab) {
    const auto min_I = std::abs(Iab - p_c.isospin());
    const auto max_I = Iab + p_c.isospin();
    if (min_I <= Res.isospin() && Res.isospin() <= max_I) {
      allowed_I_ab.push_back(Iab);
    }
  }
  if (allowed_I_ab.size() != 1) {
    throw std::runtime_error(
        "The coupled 3-body isospin state is not uniquely defined for " +
        Res.name() + " -> " + p_a.name() + " " + p_b.name() + " " + p_c.name());
  }
  const auto I_ab = allowed_I_ab[0];

  const int I_abz = p_a.isospin3() + p_b.isospin3();
  const double cg = clebsch_gordan(I_ab, p_c.isospin(), Res.isospin(), I_abz,
                                   p_c.isospin3(), Res.isospin3()) *
                    clebsch_gordan(p_a.isospin(), p_b.isospin(), I_ab,
                                   p_a.isospin3(), p_b.isospin3(), I_abz);
  return cg * cg;
}

double isospin_clebsch_gordan_sqr_2to2(const ParticleType &p_a,
                                       const ParticleType &p_b,
                                       const ParticleType &p_c,
                                       const ParticleType &p_d, const int I) {
  const int I_z = p_a.isospin3() + p_b.isospin3();

  /* Loop over total isospin in allowed range. */
  double isospin_factor = 0.;
  for (const int I_tot : I_tot_range(p_a, p_b, p_c, p_d)) {
    if (I < 0 || I_tot == I) {
      const double cg_in = isospin_clebsch_gordan_2to1(p_a, p_b, I_tot, I_z);
      const double cg_out = isospin_clebsch_gordan_2to1(p_c, p_d, I_tot, I_z);
      isospin_factor = isospin_factor + cg_in * cg_in * cg_out * cg_out;
    }
  }
  return isospin_factor;
}

}  // namespace smash
