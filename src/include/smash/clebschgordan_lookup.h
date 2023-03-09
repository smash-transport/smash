/*
 *    Copyright (c) 2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_SMASH_CLEBSCHGORDAN_LOOKUP_H_
#define SRC_INCLUDE_SMASH_CLEBSCHGORDAN_LOOKUP_H_

#include <cassert>
#include <iostream>
#include <unordered_map>

#include "smash/iomanipulators.h"

namespace smash {

/**
 * Class to store and retrieve/calculate Clebsch-Gordan coefficients.
 */
class ClebschGordan {
 public:
  /**
   * Check in the Clebsch-Gordan lookup table if the requested coefficient is
   * available. If so, return it, otherwise calculate the requested one, store
   * it in the lookup table and return it.
   *
   * \see calculate_coefficient for a description of function arguments and
   * return value.
   */
  static double coefficient(const int j_a, const int j_b, const int j_c,
                            const int m_a, const int m_b, const int m_c);

  /**
   * Auxiliary struct to be used as key in the look up table of Clebsch-Gordan
   * coefficients. It basically contains the input to retrieve one coefficient.
   * Note that this is public since it is useful to be used from client code,
   * e.g. in tests.
   */
  struct ThreeSpins {
    int j1;  ///< First isospin
    int j2;  ///< Second isospin
    int j3;  ///< Third isospin
    int m1;  ///< z component of first isospin
    int m2;  ///< z component of second isospin
    int m3;  ///< z component of third isospin

   private:
    /**
     * A utility function to avoid duplication in comparison operator(s).
     * Note that in order to use \c auto deduced returned type, this member has
     * to be defined before using it.
     *
     * @return A tuple of constant references to members.
     */
    auto tied() const { return std::tie(j1, j2, j3, m1, m2, m3); }

   public:
    /**
     * Comparison operator between two set of spin information. This is needed
     * in order to use this object in a \c std::unordered_map container.
     *
     * @param other The object to be compared to
     * @return \c true If all 6 spins value are identical
     * @return \c false otherwise
     */
    bool operator==(const ThreeSpins &other) const {
      return tied() == other.tied();
    }
  };

 private:
  /**
   * Calculate Clebsch-Gordan coefficient
   * \f$(-1)^{j_a - j_b + m_c} \sqrt{(2 j_c + 1)} \cdot [Wigner 3J symbol] \f$
   * \param[in] j_a isospin of first particle
   * \param[in] j_b isospin of second particle
   * \param[in] j_c isospin of resonance
   * \param[in] m_a z-component of isospin of first particle
   * \param[in] m_b z-component of isospin of second particle
   * \param[in] m_c z-component of isospin of resonance
   * \return Clebsch-Gordan coefficient for coupling of particles a, b and c
   *
   * Note that the calculation assumes that the isospin values (j/m) have been
   * multiplied by two (in order to be integer).
   */
  static double calculate_coefficient(const int j_a, const int j_b,
                                      const int j_c, const int m_a,
                                      const int m_b, const int m_c);

  /**
   * This is one of the possible ways to prepare a hashing mechanism to use a
   * custom object in a \c std::unordered_map container. It has been preferred
   * here to use a new \c struct instead of injecting a specialization into the
   * \c std namespace, because we are here in the \c smash namespace and the
   * object is going to be possibly better localized. Since the hash is not
   * trivial, we also preferred this approach to using a lambda function to
   * declare the hashing function.
   */
  struct ThreeSpinHash {
    /**
     * The overload of the \c operator() is the only needed ingredient to make
     * this class ready to be used as hashing algorithm.
     *
     * Since there is not (yet) a C++ standard way of combining hashes, to
     * implement functions like this one, it is necessary to choose a hashing
     * algorithm. The algorithm should minimize the hash collision (different
     * input giving the same output), but at the same time be as fast as
     * possible. One could use boost library approach, which defines
     * \code{.cpp}
     * template <typename T>
     * inline void hash_combine(std::size_t &seed, const T &val)
     * {
     *   seed ^= std::hash<T>{}(val) + 0x9e3779b9 +
     *           (seed << 6) + (seed >> 2);
     * }
     * \endcode
     * as function to combine hashes. Then this should be used here on each
     * \c in member to produce the final hash. Although it has been tested to
     * work, there is a simpler approach.
     *
     * -# It can be assumed (asserted in the code), that all \c ThreeSpins
     *    members are numbers smaller than 16 (in absolute value). Hence 5 bits
     *    are enough to represent them (numbers between -16 and 15). Using
     *    bit-shift operations, it is then possible to build a bit-sequence that
     *    in its last 6*5=30 bits contains the 6 integer information,
     *    appropriately "converted" to 5-bit sequences that are then
     *    concatenated.
     * -# To get an \c int as a 5 bits integer into a \c std::size_t variable,
     *    one way is to cast to it first and then get rid of all but 5 rightmost
     *    bits. If e.g. the result variable is 64 bits large, this means to
     *    bit-shift to the left by 64-5=59 positions and then again to the right
     *    to bring the 5 bits back to the least significant positions (this is
     *    only an example, since the size of \c std::size_t is architecture
     *    dependent).
     * -# Finally, once having the 6 5-bits integer as \c std::size_t variables,
     *    these can be combined shifting them to the left in a way that the 5
     *    relevant bits do not "overlap" (i.e. occupy different positions) and
     *    then summing them. In the implementation we use the minimum needed
     *    shift, i.e. having N=6 numbers we shift the first number by (N-1)*5,
     *    the second by (N-2)*5 and so on till the last one that does not get
     *    shifted. Since both the cast to "5-bits" numbers and the preparation
     *    of the sum involve bit-shifts, these are combined together.
     *
     * \note
     * A couple of remarks a worth:
     * - Using \c std::bitset would probably make the code easier to read, but
     *   it has been benchmarked to be some % slower at low energies.
     * - Yet another possibility would be to use mathematics/physics to come up
     *   with a map between the 6 integers and a super-index. This has been
     *   devised in \cite Rasch2004, where an algorithm to efficiently store
     *   Clebsch-Gordan coefficient is discussed. In it a way to map the three
     *   spins information to a single integer is proposed and this can be used
     *   here as hashing function, basically offering the guarantee that no hash
     *   collision will occur. However, the simpler hand-made hash discussed
     *   above and implemented in the following turned out to be more efficient
     *   with GNU compiler (equivalent with LLVM).
     *
     * @param in The spin information (meant to be used as input to be hashed)
     * @return \c std::size_t The calculated hash value
     */
    std::size_t operator()(const ThreeSpins &in) const noexcept {
      assert(in.j1 >= 0 && in.j1 < 16);
      assert(in.j2 >= 0 && in.j2 < 16);
      assert(in.j3 >= 0 && in.j3 < 16);
      assert(std::abs(in.m1) < 16);
      assert(std::abs(in.m2) < 16);
      assert(std::abs(in.m3) < 16);
      /*
       * Although strictly speaking, this would only generate a very poor hash,
       * we prefer making compilation fail if size_t has less than 32 bits.
       * This is after all very unlikely and we'll deal with it only if needed.
       */
      static_assert(sizeof(std::size_t) >= 4);
      // This is the amount to shift to obtain "5-bits numbers"
      constexpr auto bitshift = sizeof(size_t) * 8 - 5;
      // The different shift to the right make the 5-bit occupy different bits
      return (static_cast<std::size_t>(in.j1) << bitshift >> (bitshift - 25)) +
             (static_cast<std::size_t>(in.j2) << bitshift >> (bitshift - 20)) +
             (static_cast<std::size_t>(in.j3) << bitshift >> (bitshift - 15)) +
             (static_cast<std::size_t>(in.m1) << bitshift >> (bitshift - 10)) +
             (static_cast<std::size_t>(in.m2) << bitshift >> (bitshift - 5)) +
             (static_cast<std::size_t>(in.m3) << bitshift >> bitshift);
    }
  };

  /**
   * Tabulation of Clebsch-Gordan coefficients. The C++ code to produce this
   * member declaration can be found in the "tabulate" unit test of this file.
   */
  inline static std::unordered_map<ThreeSpins, double, ThreeSpinHash>
      lookup_table = {
          {{0, 0, 0, +0, +0, +0}, 1.00000000000000000},
          {{0, 1, 1, +0, -1, -1}, 1.00000000000000022},
          {{0, 1, 1, +0, +1, +1}, 1.00000000000000022},
          {{0, 2, 2, +0, -2, -2}, 0.99999999999999989},
          {{0, 2, 2, +0, +0, +0}, 0.99999999999999989},
          {{0, 2, 2, +0, +2, +2}, 0.99999999999999989},
          {{0, 3, 3, +0, -3, -3}, 1.00000000000000000},
          {{0, 3, 3, +0, -1, -1}, 0.99999999999999989},
          {{0, 3, 3, +0, +1, +1}, 0.99999999999999989},
          {{0, 3, 3, +0, +3, +3}, 1.00000000000000000},
          {{1, 0, 1, -1, +0, -1}, 1.00000000000000022},
          {{1, 0, 1, +1, +0, +1}, 1.00000000000000022},
          {{1, 1, 0, -1, +1, +0}, -0.70710678118654757},
          {{1, 1, 0, +1, -1, +0}, 0.70710678118654757},
          {{1, 1, 2, -1, -1, -2}, 0.99999999999999989},
          {{1, 1, 2, -1, +1, +0}, 0.70710678118654746},
          {{1, 1, 2, +1, -1, +0}, 0.70710678118654746},
          {{1, 1, 2, +1, +1, +2}, 0.99999999999999989},
          {{1, 2, 1, -1, +0, -1}, -0.57735026918962584},
          {{1, 2, 1, -1, +2, +1}, -0.81649658092772615},
          {{1, 2, 1, +1, -2, -1}, 0.81649658092772615},
          {{1, 2, 1, +1, +0, +1}, 0.57735026918962584},
          {{1, 2, 3, -1, -2, -3}, 1.00000000000000000},
          {{1, 2, 3, -1, +0, -1}, 0.81649658092772615},
          {{1, 2, 3, -1, +2, +1}, 0.57735026918962584},
          {{1, 2, 3, +1, -2, -1}, 0.57735026918962584},
          {{1, 2, 3, +1, +0, +1}, 0.81649658092772615},
          {{1, 2, 3, +1, +2, +3}, 1.00000000000000000},
          {{1, 3, 2, -1, -1, -2}, -0.49999999999999983},
          {{1, 3, 2, -1, +1, +0}, -0.70710678118654724},
          {{1, 3, 2, -1, +3, +2}, -0.86602540378443837},
          {{1, 3, 2, +1, -3, -2}, 0.86602540378443837},
          {{1, 3, 2, +1, -1, +0}, 0.70710678118654724},
          {{1, 3, 2, +1, +1, +2}, 0.49999999999999983},
          {{1, 3, 4, -1, -3, -4}, 1.00000000000000022},
          {{1, 3, 4, -1, -1, -2}, 0.86602540378443871},
          {{1, 3, 4, -1, +1, +0}, 0.70710678118654746},
          {{1, 3, 4, -1, +3, +2}, 0.49999999999999994},
          {{1, 3, 4, +1, -3, -2}, 0.49999999999999994},
          {{1, 3, 4, +1, -1, +0}, 0.70710678118654746},
          {{1, 3, 4, +1, +1, +2}, 0.86602540378443871},
          {{1, 3, 4, +1, +3, +4}, 1.00000000000000022},
          {{2, 0, 2, -2, +0, -2}, 0.99999999999999989},
          {{2, 0, 2, +0, +0, +0}, 0.99999999999999989},
          {{2, 0, 2, +2, +0, +2}, 0.99999999999999989},
          {{2, 1, 1, -2, +1, -1}, -0.81649658092772615},
          {{2, 1, 1, +0, -1, -1}, 0.57735026918962584},
          {{2, 1, 1, +0, +1, +1}, -0.57735026918962584},
          {{2, 1, 1, +2, -1, +1}, 0.81649658092772615},
          {{2, 1, 3, -2, -1, -3}, 1.00000000000000000},
          {{2, 1, 3, -2, +1, -1}, 0.57735026918962584},
          {{2, 1, 3, +0, -1, -1}, 0.81649658092772615},
          {{2, 1, 3, +0, +1, +1}, 0.81649658092772615},
          {{2, 1, 3, +2, -1, +1}, 0.57735026918962584},
          {{2, 1, 3, +2, +1, +3}, 1.00000000000000000},
          {{2, 2, 0, -2, +2, +0}, 0.57735026918962584},
          {{2, 2, 0, +0, +0, +0}, -0.57735026918962573},
          {{2, 2, 0, +2, -2, +0}, 0.57735026918962584},
          {{2, 2, 2, -2, +0, -2}, -0.70710678118654735},
          {{2, 2, 2, -2, +2, +0}, -0.70710678118654735},
          {{2, 2, 2, +0, -2, -2}, 0.70710678118654735},
          {{2, 2, 2, +0, +2, +2}, -0.70710678118654735},
          {{2, 2, 2, +2, -2, +0}, 0.70710678118654735},
          {{2, 2, 2, +2, +0, +2}, 0.70710678118654735},
          {{2, 2, 4, -2, -2, -4}, 1.00000000000000022},
          {{2, 2, 4, -2, +0, -2}, 0.70710678118654746},
          {{2, 2, 4, -2, +2, +0}, 0.40824829046386313},
          {{2, 2, 4, +0, -2, -2}, 0.70710678118654746},
          {{2, 2, 4, +0, +0, +0}, 0.81649658092772615},
          {{2, 2, 4, +0, +2, +2}, 0.70710678118654746},
          {{2, 2, 4, +2, -2, +0}, 0.40824829046386313},
          {{2, 2, 4, +2, +0, +2}, 0.70710678118654746},
          {{2, 2, 4, +2, +2, +4}, 1.00000000000000022},
          {{2, 3, 1, -2, +1, -1}, 0.40824829046386302},
          {{2, 3, 1, -2, +3, +1}, 0.70710678118654746},
          {{2, 3, 1, +0, -1, -1}, -0.57735026918962573},
          {{2, 3, 1, +0, +1, +1}, -0.57735026918962573},
          {{2, 3, 1, +2, -3, -1}, 0.70710678118654746},
          {{2, 3, 1, +2, -1, +1}, 0.40824829046386302},
          {{2, 3, 3, -2, -1, -3}, -0.63245553203367610},
          {{2, 3, 3, -2, +1, -1}, -0.73029674334022165},
          {{2, 3, 3, -2, +3, +1}, -0.63245553203367610},
          {{2, 3, 3, +0, -3, -3}, 0.77459666924148352},
          {{2, 3, 3, +0, -1, -1}, 0.25819888974716126},
          {{2, 3, 3, +0, +1, +1}, -0.25819888974716126},
          {{2, 3, 3, +0, +3, +3}, -0.77459666924148352},
          {{2, 3, 3, +2, -3, -1}, 0.63245553203367610},
          {{2, 3, 3, +2, -1, +1}, 0.73029674334022165},
          {{2, 3, 3, +2, +1, +3}, 0.63245553203367610},
          {{2, 3, 5, -2, -3, -5}, 0.99999999999999989},
          {{2, 3, 5, -2, -1, -3}, 0.77459666924148318},
          {{2, 3, 5, -2, +1, -1}, 0.54772255750516596},
          {{2, 3, 5, -2, +3, +1}, 0.31622776601683794},
          {{2, 3, 5, +0, -3, -3}, 0.63245553203367599},
          {{2, 3, 5, +0, -1, -1}, 0.77459666924148318},
          {{2, 3, 5, +0, +1, +1}, 0.77459666924148318},
          {{2, 3, 5, +0, +3, +3}, 0.63245553203367599},
          {{2, 3, 5, +2, -3, -1}, 0.31622776601683794},
          {{2, 3, 5, +2, -1, +1}, 0.54772255750516596},
          {{2, 3, 5, +2, +1, +3}, 0.77459666924148318},
          {{2, 3, 5, +2, +3, +5}, 0.99999999999999989},
          {{3, 0, 3, -3, +0, -3}, 1.00000000000000000},
          {{3, 0, 3, -1, +0, -1}, 0.99999999999999989},
          {{3, 0, 3, +1, +0, +1}, 0.99999999999999989},
          {{3, 0, 3, +3, +0, +3}, 1.00000000000000000},
          {{3, 1, 2, -3, +1, -2}, -0.86602540378443837},
          {{3, 1, 2, -1, -1, -2}, 0.49999999999999983},
          {{3, 1, 2, -1, +1, +0}, -0.70710678118654724},
          {{3, 1, 2, +1, -1, +0}, 0.70710678118654724},
          {{3, 1, 2, +1, +1, +2}, -0.49999999999999983},
          {{3, 1, 2, +3, -1, +2}, 0.86602540378443837},
          {{3, 1, 4, -3, -1, -4}, 1.00000000000000022},
          {{3, 1, 4, -3, +1, -2}, 0.49999999999999994},
          {{3, 1, 4, -1, -1, -2}, 0.86602540378443871},
          {{3, 1, 4, -1, +1, +0}, 0.70710678118654746},
          {{3, 1, 4, +1, -1, +0}, 0.70710678118654746},
          {{3, 1, 4, +1, +1, +2}, 0.86602540378443871},
          {{3, 1, 4, +3, -1, +2}, 0.49999999999999994},
          {{3, 1, 4, +3, +1, +4}, 1.00000000000000022},
          {{3, 2, 1, -3, +2, -1}, 0.70710678118654746},
          {{3, 2, 1, -1, +0, -1}, -0.57735026918962573},
          {{3, 2, 1, -1, +2, +1}, 0.40824829046386302},
          {{3, 2, 1, +1, -2, -1}, 0.40824829046386302},
          {{3, 2, 1, +1, +0, +1}, -0.57735026918962573},
          {{3, 2, 1, +3, -2, +1}, 0.70710678118654746},
          {{3, 2, 3, -3, +0, -3}, -0.77459666924148352},
          {{3, 2, 3, -3, +2, -1}, -0.63245553203367610},
          {{3, 2, 3, -1, -2, -3}, 0.63245553203367610},
          {{3, 2, 3, -1, +0, -1}, -0.25819888974716126},
          {{3, 2, 3, -1, +2, +1}, -0.73029674334022165},
          {{3, 2, 3, +1, -2, -1}, 0.73029674334022165},
          {{3, 2, 3, +1, +0, +1}, 0.25819888974716126},
          {{3, 2, 3, +1, +2, +3}, -0.63245553203367610},
          {{3, 2, 3, +3, -2, +1}, 0.63245553203367610},
          {{3, 2, 3, +3, +0, +3}, 0.77459666924148352},
          {{3, 2, 5, -3, -2, -5}, 0.99999999999999989},
          {{3, 2, 5, -3, +0, -3}, 0.63245553203367599},
          {{3, 2, 5, -3, +2, -1}, 0.31622776601683794},
          {{3, 2, 5, -1, -2, -3}, 0.77459666924148318},
          {{3, 2, 5, -1, +0, -1}, 0.77459666924148318},
          {{3, 2, 5, -1, +2, +1}, 0.54772255750516596},
          {{3, 2, 5, +1, -2, -1}, 0.54772255750516596},
          {{3, 2, 5, +1, +0, +1}, 0.77459666924148318},
          {{3, 2, 5, +1, +2, +3}, 0.77459666924148318},
          {{3, 2, 5, +3, -2, +1}, 0.31622776601683794},
          {{3, 2, 5, +3, +0, +3}, 0.63245553203367599},
          {{3, 2, 5, +3, +2, +5}, 0.99999999999999989},
          {{3, 3, 0, -3, +3, +0}, -0.49999999999999994},
          {{3, 3, 0, -1, +1, +0}, 0.49999999999999994},
          {{3, 3, 0, +1, -1, +0}, -0.49999999999999994},
          {{3, 3, 0, +3, -3, +0}, 0.49999999999999994},
          {{3, 3, 2, -3, +1, -2}, 0.54772255750516596},
          {{3, 3, 2, -3, +3, +0}, 0.67082039324993692},
          {{3, 3, 2, -1, -1, -2}, -0.63245553203367599},
          {{3, 3, 2, -1, +1, +0}, -0.22360679774997907},
          {{3, 3, 2, -1, +3, +2}, 0.54772255750516596},
          {{3, 3, 2, +1, -3, -2}, 0.54772255750516596},
          {{3, 3, 2, +1, -1, +0}, -0.22360679774997907},
          {{3, 3, 2, +1, +1, +2}, -0.63245553203367599},
          {{3, 3, 2, +3, -3, +0}, 0.67082039324993692},
          {{3, 3, 2, +3, -1, +2}, 0.54772255750516596},
          {{3, 3, 4, -3, -1, -4}, -0.70710678118654746},
          {{3, 3, 4, -3, +1, -2}, -0.70710678118654746},
          {{3, 3, 4, -3, +3, +0}, -0.49999999999999994},
          {{3, 3, 4, -1, -3, -4}, 0.70710678118654746},
          {{3, 3, 4, -1, +1, +0}, -0.49999999999999994},
          {{3, 3, 4, -1, +3, +2}, -0.70710678118654746},
          {{3, 3, 4, +1, -3, -2}, 0.70710678118654746},
          {{3, 3, 4, +1, -1, +0}, 0.49999999999999994},
          {{3, 3, 4, +1, +3, +4}, -0.70710678118654746},
          {{3, 3, 4, +3, -3, +0}, 0.49999999999999994},
          {{3, 3, 4, +3, -1, +2}, 0.70710678118654746},
          {{3, 3, 4, +3, +1, +4}, 0.70710678118654746},
          {{3, 3, 6, -3, -3, -6}, 1.00000000000000022},
          {{3, 3, 6, -3, -1, -4}, 0.70710678118654746},
          {{3, 3, 6, -3, +1, -2}, 0.44721359549995793},
          {{3, 3, 6, -3, +3, +0}, 0.22360679774997894},
          {{3, 3, 6, -1, -3, -4}, 0.70710678118654746},
          {{3, 3, 6, -1, -1, -2}, 0.77459666924148352},
          {{3, 3, 6, -1, +1, +0}, 0.67082039324993670},
          {{3, 3, 6, -1, +3, +2}, 0.44721359549995793},
          {{3, 3, 6, +1, -3, -2}, 0.44721359549995793},
          {{3, 3, 6, +1, -1, +0}, 0.67082039324993670},
          {{3, 3, 6, +1, +1, +2}, 0.77459666924148352},
          {{3, 3, 6, +1, +3, +4}, 0.70710678118654746},
          {{3, 3, 6, +3, -3, +0}, 0.22360679774997894},
          {{3, 3, 6, +3, -1, +2}, 0.44721359549995793},
          {{3, 3, 6, +3, +1, +4}, 0.70710678118654746},
          {{3, 3, 6, +3, +3, +6}, 1.00000000000000022},
  };
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CLEBSCHGORDAN_LOOKUP_H_
