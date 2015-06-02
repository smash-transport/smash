/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/lattice.h"
#include "include/logging.h"

namespace Smash {

template <typename T>
RectangularLattice<T>::RectangularLattice(std::array<float, 3> l,
                                       std::array<size_t, 3> n,
                                       std::array<float, 3> origin,
                                       bool per)
  : l_(l),
    n_(n),
    csize_{l[0]/n[0], l[1]/n[1], l[2]/n[2]},
    origin_(origin),
    periodic_(per) {
  lattice_.resize(n_[0]*n_[1]*n_[2]);
  const auto &log = logger<LogArea::Lattice>();
  log.info("Rectangular lattice created: sizes[fm] = (",
           l_[0], ",", l_[1], ",", l_[2], "), dims = (",
           n_[0], ",", n_[1], ",", n_[2], "), origin = (",
           origin_[0], ",", origin_[1], ",", origin_[2],
           "), periodic: ", periodic_);
}

template <typename T>
void RectangularLattice<T>::iterate_sublattice(
                 const std::array<int, 3> &lower_bounds,
                 const std::array<int, 3> &upper_bounds,
                 const std::function<void(T&, int, int, int)> &func) {
  const auto &log = logger<LogArea::Lattice>();
  log.debug("Iterating sublattice with lower bound index (",
            lower_bounds[0], ",", lower_bounds[1], ",", lower_bounds[2],
            "), upper bound index (",
            upper_bounds[0], ",", upper_bounds[1], ",", upper_bounds[2], ")");

  if (periodic_) {
    for (int iz = lower_bounds[2]; iz < upper_bounds[2]; iz++) {
      const int z_offset = positive_modulo(iz, n_[2]) * n_[1];
      for (int iy = lower_bounds[1]; iy < upper_bounds[1]; iy++) {
        const int y_offset = n_[0] * (positive_modulo(iy, n_[1]) + z_offset);
        for (int ix = lower_bounds[0]; ix < upper_bounds[0]; ix++) {
          const int index = positive_modulo(ix, n_[0]) + y_offset;
          func(lattice_[index], ix, iy, iz);
        }
      }
    }
  } else {
    for (int iz = lower_bounds[2]; iz < upper_bounds[2]; iz++) {
      const int z_offset = iz * n_[1];
      for (int iy = lower_bounds[1]; iy < upper_bounds[1]; iy++) {
        const int y_offset = n_[0] * (iy + z_offset);
        for (int ix = lower_bounds[0]; ix < upper_bounds[0]; ix++) {
          func(lattice_[ix + y_offset], ix, iy, iz);
        }
      }
    }
  }
}

template <typename T>
void RectangularLattice<T>::iterate_in_radius(
        const ThreeVector& point, const double r_cut,
        const std::function<void(T&, int, int, int)> &func) {
  std::array<int, 3> l_bounds, u_bounds;

  // Array holds value at the cell center: r_center = r_0 + (i+0.5)cell_size,
  // where i is index in any direction. Therefore we want cells with condition
  // (r-r_cut)*csize - 0.5 < i < (r+r_cut)*csize - 0.5, r = r_center - r_0
  for (int i = 0; i < 3; i++) {
    l_bounds[i] = int_ceil((point[i] - origin_[i] - r_cut)/csize_[i] - 0.5);
    u_bounds[i] = int_ceil((point[i] - origin_[i] + r_cut)/csize_[i] - 0.5);
  }

  if (!periodic_) {
    for (int i = 0; i < 3; i++) {
      if (l_bounds[i] < 0) {
        l_bounds[i] = 0;
      }
      if (u_bounds[i] > n_[i]) {
        u_bounds[i] = n_[i];
      }
      if (l_bounds[i] > n_[i] || u_bounds[i] < 0) {
        return;
      }
    }
  }
  iterate_sublattice(l_bounds, u_bounds, func);
}

}  // namespace Smash
