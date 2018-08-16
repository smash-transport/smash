/*
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_CLEBSCHGORDAN_H_
#define SRC_INCLUDE_CLEBSCHGORDAN_H_

#include <algorithm>

#include "particletype.h"

namespace smash {

/**
 * Calculate Clebsch-Gordan coefficient
 * \f$(-1)^{j_a - j_b + m_c} \sqrt{(2 j_c + 1)} \cdot [Wigner 3J symbol] \f$
 * \param[in] j_a spin of first particle
 * \param[in] j_b spin of second particle
 * \param[in] j_c spin of resonance
 * \param[in] m_a isospin of first particle
 * \param[in] m_b isospin of second particle
 * \param[in] m_c isospin of resonance
 * \return Clebsch-Gordan coefficient for coupling of particles a, b and c
 *
 * Note that the calculation assumes that the spin/isospin values (j/m)
 * have been multiplied by two (in order to be integer).
 */
double clebsch_gordan(const int j_a, const int j_b, const int j_c,
                      const int m_a, const int m_b, const int m_c);

/**
 * Calculate the squared isospin Clebsch-Gordan coefficient for two particles
 * p_a and p_b coupling to a resonance Res.
 * \param[in] p_a Information on spin/isospin of particle a
 * \param[in] p_b Information on spin/isospin of particle b
 * \param[in] Res Information on spin/isospin of resonance
 * \return Clebsch-Gordan squared for 2->1 reaction
 */
inline double isospin_clebsch_gordan_sqr_2to1(const ParticleType &p_a,
                                              const ParticleType &p_b,
                                              const ParticleType &Res) {
  const double cg =
      clebsch_gordan(p_a.isospin(), p_b.isospin(), Res.isospin(),
                     p_a.isospin3(), p_b.isospin3(), Res.isospin3());
  return cg * cg;
}

/**
 * Calculate the squared isospin Clebsch-Gordan coefficient for three particles
 * p_a, p_b and p_c coupling to a resonance Res.
 * \param[in] p_a Information on spin/isospin of particle a
 * \param[in] p_b Information on spin/isospin of particle b
 * \param[in] p_c Information on spin/isospin of particle c
 * \param[in] Res Information on spin/isospin of resonance
 * \return Clebsch-Gordan squared for 3->1 reaction
 */
double isospin_clebsch_gordan_sqr_3to1(const ParticleType &p_a,
                                       const ParticleType &p_b,
                                       const ParticleType &p_c,
                                       const ParticleType &Res);

/**
 * Calculate the squared isospin Clebsch-Gordan coefficient for a
 * 2-to-2 reaction A + B -> C + D. If a total isospin value I is given
 * (doubled in order to be integer), then only contributions with that total
 * isospin will be counted.
 * \param[in] p_a Information on spin/isospin of particle a
 * \param[in] p_b Information on spin/isospin of particle b
 * \param[in] p_c Information on spin/isospin of particle c
 * \param[in] p_d Information on spin/isospin of particle d
 * \param[in] I total isospin of the reaction
 * \return Clebsch-Gordan squared for 2->2 reaction
 */
double isospin_clebsch_gordan_sqr_2to2(const ParticleType &p_a,
                                       const ParticleType &p_b,
                                       const ParticleType &p_c,
                                       const ParticleType &p_d,
                                       const int I = -1);

/// Range of total isospin for reaction of particle a with particle b.
class I_tot_range {
 private:
  /// Value of minimum total isospin
  int I_min_;
  /// Value of maximum total isospin
  int I_max_;

 public:
  /**
   * Get the allowed range of total isospin for a collision a + b.
   * \param a Particle a.
   * \param a Particle b.
   * \return Maximum and minimum of allowed values.
   */
  I_tot_range(const ParticleType &p_a, const ParticleType &p_b) {
    // Compute total isospin range with given particles.
    const int I_z_abs = std::abs(p_a.isospin3() + p_b.isospin3());
    I_max_ = p_a.isospin() + p_b.isospin();
    I_min_ = std::max(std::abs(p_a.isospin() - p_b.isospin()), I_z_abs);
  }

  /**
   * Get the allowed range of total isospin for a collision a + b <-> c + d.
   * \param a Particle a.
   * \param a Particle b.
   * \param a Particle c.
   * \param a Particle d.
   * \return Maximum and minimum of allowed values or empty range, if reaction
   * is forbidden due to isospin.
   */
  I_tot_range(const ParticleType &p_a, const ParticleType &p_b,
              const ParticleType &p_c, const ParticleType &p_d) {
    // Compute total isospin range with given initial and final particles.
    const int I_z = p_a.isospin3() + p_b.isospin3();
    if (I_z != p_c.isospin3() + p_d.isospin3()) {
      /* This reaction is forbidden by isospin conservation.
       * Set impossible values to make sure an empty range is returned. */
      I_min_ = 1;
      I_max_ = 0;
      return;
    }
    I_max_ =
        std::min(p_a.isospin() + p_b.isospin(), p_c.isospin() + p_d.isospin());
    I_min_ = std::max(std::abs(p_a.isospin() - p_b.isospin()),
                      std::abs(p_c.isospin() - p_d.isospin()));
    I_min_ = std::max(I_min_, std::abs(I_z));
  }

  /// Iterator class for determination of total isospin.
  class iterator : public std::iterator<std::forward_iterator_tag, int> {
   private:
    int c_;
    I_tot_range &parent_;

   public:
    /**
     * Construct an iterator.
     * \param start Initial value.
     * \param parent Parent class giving the total isospin range.
     * \return The constructed iterator.
     */
    iterator(int start, I_tot_range &parent) : c_(start), parent_(parent) {}
    /// \return Current element of the iterator.
    int operator*() { return c_; }
    /// \return Next element of the iterator.
    const iterator *operator++() {
      c_ -= 2;
      return this;
    }
    /// \return Next element of the iterator.
    iterator operator++(int) {
      c_ -= 2;
      return iterator(c_ + 2, parent_);
    }
    /**
     * \param other Other iterator.
     * \return Whether both iterators are equal.
     */
    bool operator==(const iterator &other) { return c_ == other.c_; }
    /**
     * \param other Other iterator.
     * \return Whether both iterators are not equal.
     */
    bool operator!=(const iterator &other) { return c_ != other.c_; }
  };

  /// \return Beginning of iterator.
  iterator begin() { return iterator(I_max_, *this); }
  /// \return End of iterator.
  iterator end() {
    if (I_min_ > I_max_) {
      return begin();
    }
    return iterator(I_min_ - 2, *this);
  }
};

}  // namespace smash

#endif  // SRC_INCLUDE_CLEBSCHGORDAN_H_
