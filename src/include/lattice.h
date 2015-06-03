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
                     const std::array<float, 3> &origin, bool per,
                     const LatticeUpdate upd)
  : l_(l),
    n_(n),
    csize_{l[0]/n[0], l[1]/n[1], l[2]/n[2]},
    origin_(origin),
    periodic_(per),
    when_update_(upd) {
    lattice_.resize(n_[0]*n_[1]*n_[2]);
    const auto &log = logger<LogArea::Lattice>();
    log.info("Rectangular lattice created: sizes[fm] = (",
             l_[0], ",", l_[1], ",", l_[2], "), dims = (",
             n_[0], ",", n_[1], ",", n_[2], "), origin = (",
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
           (ix < 0 || ix >= n_[0] ||
            iy < 0 || iy >= n_[1] ||
            iz < 0 || iz >= n_[2]);
  }

  /// Returns coordinate of cell center given its index
  inline ThreeVector cell_center(int ix, int iy, int iz) {
    return ThreeVector(origin_[0] + csize_[0]*(ix + 0.5),
                       origin_[1] + csize_[1]*(iy + 0.5),
                       origin_[2] + csize_[2]*(iz + 0.5));
  }

  /// Getter for lattice sizes
  const std::array<float, 3>& lattice_sizes() const { return l_; }

  /// Getter for lattice dimensions
  const std::array<int, 3>& dimensions() const { return n_; }

  /// Getter for cell sizes
  const std::array<float, 3>& cell_sizes() const { return csize_; }

  /// Getter for periodicity
  bool periodic() const { return periodic_; }

  /// Getter for update case
  LatticeUpdate when_update() const { return when_update_; }

  /// Iterators and accessors
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  iterator begin() { return lattice_.begin(); }
  const_iterator begin() const { return lattice_.begin(); }
  iterator end() { return lattice_.end(); }
  const_iterator end() const { return lattice_.end(); }
  T &operator[](std::size_t i) { return lattice_[i]; }
  T operator[](std::size_t i) const { return lattice_[i]; }
  std::size_t size() const { return lattice_.size(); }

  T& node(int ix, int iy, int iz) {
    return periodic_ ?
           lattice_[positive_modulo(ix, n_[0])   + n_[0] *
                    (positive_modulo(iy, n_[1])  + n_[1] *
                     positive_modulo(iz, n_[2]))] :
           lattice_[ix + n_[0] * (iy + n_[1] * iz)];
  }

  /**
   * A sub-lattice iterator, which iterates in a 3D-structured manner.
   * Gives index of the node it goes through: ix, iy, iz.
   */
  void iterate_sublattice(const std::array<int, 3> &lower_bounds,
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

  };

  /**
   * Iterates only nodes, whose cell centers lie not further than r_cut in
   * x,y,z directions from the given point. Useful for adding quantities
   * from one particle to the lattice.
   */
  void iterate_in_radius(const ThreeVector& point, const double r_cut,
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

  };

 protected:
  /// The lattice itself, array containing physical quantities
  std::vector<T> lattice_;
  /// Lattice sizes in x,y,z directions
  const std::array<float, 3> l_;
  /// Number of cells in x,y,z directions
  const std::array<int, 3> n_;
  /// Cell sizes in x,y,z directions
  const std::array<float, 3> csize_;
  /// Coordinates of the left down nearer corner
  const std::array<float, 3> origin_;
  /// Periodicity
  const bool periodic_;
  /// when lattice should be recalculated
  const LatticeUpdate when_update_;

 private:
  /// Returns ceiling as integer
  inline int int_ceil(double x) {
    return static_cast<int>(std::ceil(x));
  }
  /// Returns division modulo, which is always between 0 and n-1
  // i%n is not suitable, because it returns results from -(n-1) to n-1
  inline int positive_modulo(int i, int n) {
    return (i % n + n) % n;
  }
};

}  // namespace Smash

#endif  // SRC_INCLUDE_LATTICE_H_
