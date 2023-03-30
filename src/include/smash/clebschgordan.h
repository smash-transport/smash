/*
 *    Copyright (c) 2013-2018,2020,2022-2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_SMASH_CLEBSCHGORDAN_H_
#define SRC_INCLUDE_SMASH_CLEBSCHGORDAN_H_

#include <algorithm>

#include "clebschgordan_lookup.h"
#include "particletype.h"

namespace smash {

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
  const double cg = ClebschGordan::coefficient(p_a.isospin(), p_b.isospin(),
                                               Res.isospin(), p_a.isospin3(),
                                               p_b.isospin3(), Res.isospin3());
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
   * \param p_a Particle a.
   * \param p_b Particle b.
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
   * \param p_a Particle a.
   * \param p_b Particle b.
   * \param p_c Particle c.
   * \param p_d Particle d.
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

  /**
   * Iterator class for determination of total isospin.
   *
   * See documentation of smash::Particles::GenericIterator class for more
   * information about exposed types.
   */
  class iterator {
   private:
    /// Element of the iterator
    int c_;
    /// Parent class giving the total isospin range.
    I_tot_range &parent_;

   public:
    /// <b>Required by STL:</b> expose iterator_category
    using iterator_category = std::forward_iterator_tag;
    /// <b>Required by STL:</b> expose value_type
    using value_type = int;
    /// <b>Required by STL:</b> expose difference_type
    using difference_type = int;
    /// <b>Required by STL:</b> expose pointer
    using pointer = int *;
    /// <b>Required by STL:</b> expose reference
    using reference = int &;
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

#endif  // SRC_INCLUDE_SMASH_CLEBSCHGORDAN_H_
