/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_LATTICE_H_
#define SRC_INCLUDE_LATTICE_H_

#include <cstring>
#include <array>
#include <functional>
#include <utility>
#include <vector>

#include "forwarddeclarations.h"
#include "fourvector.h"
#include "logging.h"

namespace Smash {

/** Enumerator option for lattice updates.
 *  Lattice update is a costly operation and should be performed only if
 *  necessary. Possible needs are: output - then it is enough to update
 *  lattice just before output, need for physics - update every timestep
 *  is unavoidable. Other needs may occur - that's why enum, not bool.
 */
enum class LatticeUpdate {
  AtOutput = 0,
  EveryTimestep = 1
};

/**
  * A container class to hold all the arrays on the lattice and access them.
  */
template <typename T>
class RectangularLattice {
 public:
  /**
    * Creates rectangular lattice of sizes (lx,ly,lz) fm with
    * nx, ny, nz cells in x,y,z directions respectively. Cell
    * i,j,k comprises volume ((i,i+1)lx/nx; (j,j+1)ly/ny; (k,k+1)lz/nz),
    * i: 0,nx-1; j: 0, ny-1; k: 0, nz-1;
    */
  RectangularLattice(const std::array<float, 3> &l,
                     const std::array<int, 3> &n,
                     const std::array<float, 3> &orig, bool per,
                     const LatticeUpdate upd)
  : lattice_sizes_(l),
    n_cells_(n),
    cell_sizes_{l[0]/n[0], l[1]/n[1], l[2]/n[2]},
    origin_(orig),
    periodic_(per),
    when_update_(upd) {
    lattice_.resize(n_cells_[0]*
                    n_cells_[1]*
                    n_cells_[2]);
    const auto &log = logger<LogArea::Lattice>();
    log.debug("Rectangular lattice created: sizes[fm] = (",
             lattice_sizes_[0], ",", lattice_sizes_[1], ",",
             lattice_sizes_[2], "), dims = (",
             n_cells_[0], ",", n_cells_[1], ",", n_cells_[2], "), origin = (",
             origin_[0], ",", origin_[1], ",", origin_[2],
             "), periodic: ", periodic_);
  }

  /// Sets all values on lattice to zeros
  void reset() {
    std::fill(lattice_.begin(), lattice_.end(), T());
  }

  /// Checks if 3D index is out of lattice bounds
  inline bool out_of_bounds(int ix, int iy, int iz) {
    return !periodic_ &&
           (ix < 0 || ix >= n_cells_[0] ||
            iy < 0 || iy >= n_cells_[1] ||
            iz < 0 || iz >= n_cells_[2]);
  }

  /// Returns coordinate of cell center given its index
  inline ThreeVector cell_center(int ix, int iy, int iz) {
    return ThreeVector(origin_[0] + cell_sizes_[0] * (ix + 0.5f),
                       origin_[1] + cell_sizes_[1] * (iy + 0.5f),
                       origin_[2] + cell_sizes_[2] * (iz + 0.5f));
  }

  /// Returns lengths of the lattice in x,y,z directions
  const std::array<float, 3>& lattice_sizes() const { return lattice_sizes_; }

  /// Returns number of cells in x,y,z directions
  const std::array<int, 3>& dimensions() const { return n_cells_; }

  /// Returns lengths of one cell in x,y,z directions
  const std::array<float, 3>& cell_sizes() const { return cell_sizes_; }

  /// Returns lattice origin: left, down, near corner coordinates
  const std::array<float, 3>& origin() const { return origin_; }

  /// Returns if lattice is periodic or not
  bool periodic() const { return periodic_; }

  /// Returns the enum, which tells at which time lattice wants to be updated
  LatticeUpdate when_update() const { return when_update_; }

  /// Iterators and accessors
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  iterator begin() { return lattice_.begin(); }
  const_iterator begin() const { return lattice_.begin(); }
  iterator end() { return lattice_.end(); }
  const_iterator end() const { return lattice_.end(); }
  T &operator[](std::size_t i) { return lattice_[i]; }
  const T& operator[](std::size_t i) const { return lattice_[i]; }
  std::size_t size() const { return lattice_.size(); }

  T& node(int ix, int iy, int iz) {
    return periodic_ ?
           lattice_[positive_modulo(ix, n_cells_[0]) + n_cells_[0] *
                      (positive_modulo(iy, n_cells_[1]) + n_cells_[1] *
                         positive_modulo(iz, n_cells_[2]))] :
           lattice_[ix + n_cells_[0] * (iy + n_cells_[1] * iz)];
  }

  /**
   * A sub-lattice iterator, which iterates in a 3D-structured manner.
   * Gives index of the node it goes through: ix, iy, iz.
   */
  template <typename F>
  void iterate_sublattice(const std::array<int, 3> &lower_bounds,
                          const std::array<int, 3> &upper_bounds,
                          F &&func) {
    const auto &log = logger<LogArea::Lattice>();
    log.debug("Iterating sublattice with lower bound index (",
              lower_bounds[0], ",", lower_bounds[1], ",", lower_bounds[2],
              "), upper bound index (",
              upper_bounds[0], ",", upper_bounds[1], ",", upper_bounds[2], ")");

    if (periodic_) {
      for (int iz = lower_bounds[2]; iz < upper_bounds[2]; iz++) {
        const int z_offset = positive_modulo(iz, n_cells_[2]) * n_cells_[1];
        for (int iy = lower_bounds[1]; iy < upper_bounds[1]; iy++) {
          const int y_offset = n_cells_[0] *
                    (positive_modulo(iy, n_cells_[1]) + z_offset);
          for (int ix = lower_bounds[0]; ix < upper_bounds[0]; ix++) {
            const int index = positive_modulo(ix, n_cells_[0]) + y_offset;
            func(lattice_[index], ix, iy, iz);
          }
        }
      }
    } else {
      for (int iz = lower_bounds[2]; iz < upper_bounds[2]; iz++) {
        const int z_offset = iz * n_cells_[1];
        for (int iy = lower_bounds[1]; iy < upper_bounds[1]; iy++) {
          const int y_offset = n_cells_[0] * (iy + z_offset);
          for (int ix = lower_bounds[0]; ix < upper_bounds[0]; ix++) {
            func(lattice_[ix + y_offset], ix, iy, iz);
          }
        }
      }
    }

  };

  /**
   * Iterates only nodes, whose cell centers lie not further than r_cut in
   * x,y,z directions from the given point. Useful for adding quantities
   * from one particle to the lattice.
   */
  template <typename F>
  void iterate_in_radius(const ThreeVector& point, const double r_cut,
                         F &&func) {
    std::array<int, 3> l_bounds, u_bounds;

    // Array holds value at the cell center: r_center = r_0 + (i+0.5)cell_size,
    // where i is index in any direction. Therefore we want cells with condition
    // (r-r_cut)*csize - 0.5 < i < (r+r_cut)*csize - 0.5, r = r_center - r_0
    for (int i = 0; i < 3; i++) {
      l_bounds[i] = std::ceil( (point[i] - origin_[i] - r_cut) / cell_sizes_[i]
                               - 0.5f);
      u_bounds[i] = std::ceil( (point[i] - origin_[i] + r_cut) / cell_sizes_[i]
                               - 0.5f);
    }

    if (!periodic_) {
      for (int i = 0; i < 3; i++) {
        if (l_bounds[i] < 0) {
          l_bounds[i] = 0;
        }
        if (u_bounds[i] > n_cells_[i]) {
          u_bounds[i] = n_cells_[i];
        }
        if (l_bounds[i] > n_cells_[i] || u_bounds[i] < 0) {
          return;
        }
      }
    }
    iterate_sublattice(l_bounds, u_bounds, std::forward<F>(func));
  };

 protected:
  /// The lattice itself, array containing physical quantities
  std::vector<T> lattice_;
  /// Lattice sizes in x,y,z directions
  const std::array<float, 3> lattice_sizes_;
  /// Number of cells in x,y,z directions
  const std::array<int, 3> n_cells_;
  /// Cell sizes in x,y,z directions
  const std::array<float, 3> cell_sizes_;
  /// Coordinates of the left down nearer corner
  const std::array<float, 3> origin_;
  /// Periodicity
  const bool periodic_;
  /// when lattice should be recalculated
  const LatticeUpdate when_update_;

 private:
  /// Returns division modulo, which is always between 0 and n-1
  // i%n is not suitable, because it returns results from -(n-1) to n-1
  inline int positive_modulo(int i, int n) {
    /* (i % n + n) % n would be correct, but slow.
       Instead I rely on the fact that i should never go too far
       in negative region and replace i%n + n by i + 256 * n = i + (n << 8) */
    return (i + (n << 8)) % n;
  }
};

}  // namespace Smash

#endif  // SRC_INCLUDE_LATTICE_H_
