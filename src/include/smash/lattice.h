/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_LATTICE_H_
#define SRC_INCLUDE_LATTICE_H_

#include <array>
#include <cstring>
#include <functional>
#include <utility>
#include <vector>

#include "forwarddeclarations.h"
#include "fourvector.h"
#include "logging.h"
#include "numerics.h"

namespace smash {

/**
 * Enumerator option for lattice updates.
 *
 * Updating the lattice is a costly operation and should be performed only if
 * necessary. Possible needs are:
 *
 * - output: then it is enough to update lattice just before output,
 * - physics: update every time step is unavoidable.
 */
enum class LatticeUpdate {
  AtOutput = 0,
  EveryTimestep = 1,
  EveryFixedInterval = 2,
};

/**
 * A container class to hold all the arrays on the lattice and access them.
 * \tparam T The type of the contained values.
 */
template <typename T>
class RectangularLattice {
 public:
  /**
   * Rectangular lattice constructor.
   *
   * \param[in] l 3-dimensional array (lx,ly,lz) indicates the size of
   *            the lattice in x, y, z directions respectively [fm].
   * \param[in] n 3-dimensional array (nx,ny,nz) indicates the number of
   *            cells of the lattice in x, y, z directions respectively.
   *            Each cell in the lattice is labeled by three integers i, j, k
   *            where \f$i\in[0, nx-1]\f$, \f$j\in[0, ny-1]\f$,
   *            \f$k\in[0, nz-1]\f$. The sizes of each cell are given by
   *            lx/nx, ly/ny, lz/nz in x,y,z directions respectively.
   * \param[in] orig A 3-dimensional array indicating the coordinates of the
   *            origin [fm].
   * \param[in] per Boolean indicating whether a periodic boundary condition
   *            is applied.
   * \param[in] upd Enum indicating how frequently the lattice is updated.
   */
  RectangularLattice(const std::array<double, 3>& l,
                     const std::array<int, 3>& n,
                     const std::array<double, 3>& orig, bool per,
                     const LatticeUpdate upd)
      : lattice_sizes_(l),
        n_cells_(n),
        cell_sizes_{l[0] / n[0], l[1] / n[1], l[2] / n[2]},
        origin_(orig),
        periodic_(per),
        when_update_(upd) {
    lattice_.resize(n_cells_[0] * n_cells_[1] * n_cells_[2]);
    const auto& log = logger<LogArea::Lattice>();
    log.debug("Rectangular lattice created: sizes[fm] = (", lattice_sizes_[0],
              ",", lattice_sizes_[1], ",", lattice_sizes_[2], "), dims = (",
              n_cells_[0], ",", n_cells_[1], ",", n_cells_[2], "), origin = (",
              origin_[0], ",", origin_[1], ",", origin_[2],
              "), periodic: ", periodic_);
    if (n_cells_[0] < 1 || n_cells_[1] < 1 || n_cells_[2] < 1 ||
        lattice_sizes_[0] < 0.0 || lattice_sizes_[1] < 0.0 ||
        lattice_sizes_[2] < 0.0) {
      throw std::invalid_argument(
          "Lattice sizes should be positive, "
          "lattice dimensions should be > 0.");
    }
  }

  /// Sets all values on lattice to zeros.
  void reset() { std::fill(lattice_.begin(), lattice_.end(), T()); }

  /**
   * Checks if 3D index is out of lattice bounds.
   *
   * \param[in] ix The index of the cell in x direction.
   * \param[in] iy The index of the cell in y direction.
   * \param[in] iz The index of the cell in z direction.
   * \return Whether the cell is out of the lattice.
   */
  inline bool out_of_bounds(int ix, int iy, int iz) const {
    // clang-format off
    return !periodic_ &&
        (ix < 0 || ix >= n_cells_[0] ||
         iy < 0 || iy >= n_cells_[1] ||
         iz < 0 || iz >= n_cells_[2]);
    // clang-format on
  }

  /**
   * Find the coordinates of a given cell.
   *
   * \param[in] ix The index of the cell in x direction.
   * \param[in] iy The index of the cell in y direction.
   * \param[in] iz The index of the cell in z direction.
   * \return Coordinates of the center of the given cell [fm].
   */
  inline ThreeVector cell_center(int ix, int iy, int iz) const {
    return ThreeVector(origin_[0] + cell_sizes_[0] * (ix + 0.5),
                       origin_[1] + cell_sizes_[1] * (iy + 0.5),
                       origin_[2] + cell_sizes_[2] * (iz + 0.5));
  }

  /**
   * Find the coordinate of cell center given the 1d index of the cell.
   *
   * \param[in] index 1-dimensional index of the given cell. It can be
   *            related to  a 3-dimensional one by
   *            index = ix + nx (iy + iz * ny).
   * \return Coordinates of the center of the given cell [fm].
   */
  inline ThreeVector cell_center(int index) const {
    const int ix = index % n_cells_[0];
    index = index / n_cells_[0];
    const int iy = index % n_cells_[1];
    const int iz = index / n_cells_[1];
    return cell_center(ix, iy, iz);
  }

  /// \return Lengths of the lattice in x, y, z directions.
  const std::array<double, 3>& lattice_sizes() const { return lattice_sizes_; }

  /// \return Number of cells in x, y, z directions.
  const std::array<int, 3>& dimensions() const { return n_cells_; }

  /// \return Lengths of one cell in x, y, z directions.
  const std::array<double, 3>& cell_sizes() const { return cell_sizes_; }

  /// \return Lattice origin: left, down, near corner coordinates.
  const std::array<double, 3>& origin() const { return origin_; }

  /// \return If lattice is periodic or not.
  bool periodic() const { return periodic_; }

  /// \return The enum, which tells at which time lattice needs to be updated.
  LatticeUpdate when_update() const { return when_update_; }

  /// Iterator of lattice.
  using iterator = typename std::vector<T>::iterator;
  /// Const interator of lattice.
  using const_iterator = typename std::vector<T>::const_iterator;
  /// \return First element of lattice.
  iterator begin() { return lattice_.begin(); }
  /// \return First element of lattice (const).
  const_iterator begin() const { return lattice_.begin(); }
  /// \return Last element of lattice.
  iterator end() { return lattice_.end(); }
  /// \return Last element of lattice (const).
  const_iterator end() const { return lattice_.end(); }
  /// \return ith element of lattice.
  T& operator[](std::size_t i) { return lattice_[i]; }
  /// \return ith element of lattice (const).
  const T& operator[](std::size_t i) const { return lattice_[i]; }
  /// \return Size of lattice.
  std::size_t size() const { return lattice_.size(); }

  /**
   * Take the value of a cell given its 3-D indices.
   *
   * \param[in] ix The index of the cell in x direction.
   * \param[in] iy The index of the cell in y direction.
   * \param[in] iz The index of the cell in z direction.
   * \return Physical quantity evaluated at the cell center.
   */
  T& node(int ix, int iy, int iz) {
    return periodic_
               ? lattice_[positive_modulo(ix, n_cells_[0]) +
                          n_cells_[0] *
                              (positive_modulo(iy, n_cells_[1]) +
                               n_cells_[1] * positive_modulo(iz, n_cells_[2]))]
               : lattice_[ix + n_cells_[0] * (iy + n_cells_[1] * iz)];
  }

  /**
   * Interpolates lattice quantity to coordinate r. Result is stored
   * in the value variable. Returns true if coordinate r is on the
   * lattice, false if out of the lattice. In the latter case, the
   * value is set to the default value (usually 0).
   *
   * \param[in] r Position where the physical quantity would be evaluated.
   * \param[out] value Physical quantity evaluated at the nearest cell
   *             to the given position.
   * \return Boolean indicates whether the position r is located inside
   *         the lattice.
   *
   * \todo (oliiny): maybe 1-order interpolation instead of 0-order?
   */
  bool value_at(const ThreeVector& r, T& value) {
    const int ix = std::floor((r.x1() - origin_[0]) / cell_sizes_[0]);
    const int iy = std::floor((r.x2() - origin_[1]) / cell_sizes_[1]);
    const int iz = std::floor((r.x3() - origin_[2]) / cell_sizes_[2]);
    if (out_of_bounds(ix, iy, iz)) {
      value = T();
      return false;
    } else {
      value = node(ix, iy, iz);
      return true;
    }
  }

  /**
   * A sub-lattice iterator, which iterates in a 3D-structured manner and
   * calls a function on every cell.
   *
   * \tparam F Type of the function. Arguments are the current node and the 3
   * integer indices of the cell.
   * \param[in] lower_bounds Starting numbers for iterating ix, iy, iz.
   * \param[in] upper_bounds Ending numbers for iterating ix, iy, iz.
   * \param[in] func Function acting on the cells (such as taking value).
   */
  template <typename F>
  void iterate_sublattice(const std::array<int, 3>& lower_bounds,
                          const std::array<int, 3>& upper_bounds, F&& func) {
    const auto& log = logger<LogArea::Lattice>();
    log.debug("Iterating sublattice with lower bound index (", lower_bounds[0],
              ",", lower_bounds[1], ",", lower_bounds[2],
              "), upper bound index (", upper_bounds[0], ",", upper_bounds[1],
              ",", upper_bounds[2], ")");

    if (periodic_) {
      for (int iz = lower_bounds[2]; iz < upper_bounds[2]; iz++) {
        const int z_offset = positive_modulo(iz, n_cells_[2]) * n_cells_[1];
        for (int iy = lower_bounds[1]; iy < upper_bounds[1]; iy++) {
          const int y_offset =
              n_cells_[0] * (positive_modulo(iy, n_cells_[1]) + z_offset);
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
  }

  /**
   * Iterates only nodes, whose cell centers lie not further than r_cut in x, y,
   * z directions from the given point and applies a function to each node.
   * Useful for adding quantities from one particle to the lattice.
   *
   * \tparam F Type of the function. Arguments are the current node and the 3
   * integer indices of the cell.
   * \param[in] point Position, usually the position of particle [fm].
   * \param[in] r_cut Maximum distance from the cell center to the
   *            given position. [fm]
   * \param[in] func Function acting on the cells (such as taking value).
   */
  template <typename F>
  void iterate_in_radius(const ThreeVector& point, const double r_cut,
                         F&& func) {
    std::array<int, 3> l_bounds, u_bounds;

    /* Array holds value at the cell center: r_center = r_0 + (i+0.5)cell_size,
     * where i is index in any direction. Therefore we want cells with condition
     * (r-r_cut)*csize - 0.5 < i < (r+r_cut)*csize - 0.5, r = r_center - r_0 */
    for (int i = 0; i < 3; i++) {
      l_bounds[i] =
          std::ceil((point[i] - origin_[i] - r_cut) / cell_sizes_[i] - 0.5);
      u_bounds[i] =
          std::ceil((point[i] - origin_[i] + r_cut) / cell_sizes_[i] - 0.5);
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
  }

  /**
   * Checks if lattices of possibly different types have identical structure.
   *
   * \tparam L Type of the other lattice.
   * \param[in] lat The other lattice being compared with the current one
   * \return Whether the two lattices have the same sizes, cell numbers,
   *         origins, and boundary conditions.
   */
  template <typename L>
  bool identical_to_lattice(const L* lat) const {
    return n_cells_[0] == lat->dimensions()[0] &&
           n_cells_[1] == lat->dimensions()[1] &&
           n_cells_[2] == lat->dimensions()[2] &&
           std::abs(lattice_sizes_[0] - lat->lattice_sizes()[0]) <
               really_small &&
           std::abs(lattice_sizes_[1] - lat->lattice_sizes()[1]) <
               really_small &&
           std::abs(lattice_sizes_[2] - lat->lattice_sizes()[2]) <
               really_small &&
           std::abs(origin_[0] - lat->origin()[0]) < really_small &&
           std::abs(origin_[1] - lat->origin()[1]) < really_small &&
           std::abs(origin_[2] - lat->origin()[2]) < really_small &&
           periodic_ == lat->periodic();
  }

 protected:
  /// The lattice itself, array containing physical quantities.
  std::vector<T> lattice_;
  /// Lattice sizes in x, y, z directions.
  const std::array<double, 3> lattice_sizes_;
  /// Number of cells in x,y,z directions.
  const std::array<int, 3> n_cells_;
  /// Cell sizes in x, y, z directions.
  const std::array<double, 3> cell_sizes_;
  /// Coordinates of the left down nearer corner.
  const std::array<double, 3> origin_;
  /// Whether the lattice is periodic.
  const bool periodic_;
  /// When the lattice should be recalculated.
  const LatticeUpdate when_update_;

 private:
  /**
   * Returns division modulo, which is always between 0 and n-1
   * i%n is not suitable, because it returns results from -(n-1) to n-1
   *
   * \param[in] i Dividend.
   * \param[in] n Divisor.
   * \return Positive remainder.
   */
  inline int positive_modulo(int i, int n) const {
    /* (i % n + n) % n would be correct, but slow.
     * Instead I rely on the fact that i should never go too far
     * in negative region and replace i%n + n by i + 256 * n = i + (n << 8) */
    // FIXME: This should use asserts, also checking for under- or overflows.
    return (i + (n << 8)) % n;
  }
};

}  // namespace smash

#endif  // SRC_INCLUDE_LATTICE_H_
