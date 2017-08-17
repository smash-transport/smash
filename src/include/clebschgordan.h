/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_CLEBSCHGORDAN_H_
#define SRC_INCLUDE_CLEBSCHGORDAN_H_

#include <algorithm>

#include "particletype.h"

namespace Smash {

/**
 * Calculate Clebsch-Gordan coefficient
 *
 * \f$(-1)^{j_a - j_b + m_c} \sqrt(2 j_c + 1) \cdot [Wigner 3J symbol] \f$
 * Note that the calculation assumes that the spin/isospin values (j/m)
 * have been multiplied by two (in order to be integer).
 */
double clebsch_gordan(const int j_a, const int j_b, const int j_c,
                      const int m_a, const int m_b, const int m_c);

/**
 * Calculate the squared isospin Clebsch-Gordan coefficient for two particles
 * p_a and p_b coupling to a resonance Res.
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
 */
double isospin_clebsch_gordan_sqr_2to2(const ParticleType &t_a,
                                       const ParticleType &t_b,
                                       const ParticleType &t_c,
                                       const ParticleType &t_d,
                                       const int I = -1);

class I_tot_range {
 private:
  int I_min_;
  int I_max_;

 public:
  /** Get the allowed range of total isospin for a collision A + B. Returns
   * maximum and minimum of allowed values. */
  I_tot_range(const ParticleType &t_a, const ParticleType &t_b) {
    // Compute total isospin range with given particles.
    const int I_z_abs = std::abs(t_a.isospin3() + t_b.isospin3());
    I_max_ = t_a.isospin() + t_b.isospin();
    I_min_ = std::max(std::abs(t_a.isospin() - t_b.isospin()), I_z_abs);
  }

  /** Get the allowed range of total isospin for a collision A + B <-> C + D.
   * Returns maximum and minimum of allowed values. */
  I_tot_range(const ParticleType &t_a, const ParticleType &t_b,
              const ParticleType &t_c, const ParticleType &t_d) {
    // Compute total isospin range with given initial and final particles.
    const int I_z = t_a.isospin3() + t_b.isospin3();
    if (I_z != t_c.isospin3() + t_d.isospin3()) {
      // This reaction is forbidden by isospin conservation.
      // Set impossible values to make sure an empty range is returned.
      I_min_ = 1;
      I_max_ = 0;
      return;
    }
    I_max_ =
        std::min(t_a.isospin() + t_b.isospin(), t_c.isospin() + t_d.isospin());
    I_min_ = std::max(std::abs(t_a.isospin() - t_b.isospin()),
                      std::abs(t_c.isospin() - t_d.isospin()));
    I_min_ = std::max(I_min_, std::abs(I_z));
  }

  class iterator : public std::iterator<std::forward_iterator_tag, int> {
   private:
    int c_;
    I_tot_range &parent_;

   public:
    iterator(int start, I_tot_range &parent) : c_(start), parent_(parent) {}
    int operator*() { return c_; }
    const iterator *operator++() {
      c_ -= 2;
      return this;
    }
    iterator operator++(int) {
      c_ -= 2;
      return iterator(c_ + 2, parent_);
    }
    bool operator==(const iterator &other) { return c_ == other.c_; }
    bool operator!=(const iterator &other) { return c_ != other.c_; }
  };

  iterator begin() { return iterator(I_max_, *this); }
  iterator end() {
    if (I_min_ > I_max_) {
      return begin();
    }
    return iterator(I_min_ - 2, *this);
  }
};

}  // namespace Smash

#endif  // SRC_INCLUDE_CLEBSCHGORDAN_H_
