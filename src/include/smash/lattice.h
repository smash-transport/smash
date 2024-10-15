/*
 *
 *    Copyright (c) 2015-2021,2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_LATTICE_H_
#define SRC_INCLUDE_SMASH_LATTICE_H_

#include <array>
#include <cstring>
#include <functional>
#include <utility>
#include <vector>

#include "forwarddeclarations.h"
#include "fourvector.h"
#include "logging.h"
#include "numeric_cast.h"
#include "numerics.h"

namespace smash {
static constexpr int LLattice = LogArea::Lattice::id;

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
        cell_volume_{cell_sizes_[0] * cell_sizes_[1] * cell_sizes_[2]},
        origin_(orig),
        periodic_(per),
        when_update_(upd) {
    lattice_.resize(n_cells_[0] * n_cells_[1] * n_cells_[2]);
    logg[LLattice].debug(
        "Rectangular lattice created: sizes[fm] = (", lattice_sizes_[0], ",",
        lattice_sizes_[1], ",", lattice_sizes_[2], "), dims = (", n_cells_[0],
        ",", n_cells_[1], ",", n_cells_[2], "), origin = (", origin_[0], ",",
        origin_[1], ",", origin_[2], "), periodic: ", periodic_);
    if (n_cells_[0] < 1 || n_cells_[1] < 1 || n_cells_[2] < 1 ||
        lattice_sizes_[0] < 0.0 || lattice_sizes_[1] < 0.0 ||
        lattice_sizes_[2] < 0.0) {
      throw std::invalid_argument(
          "Lattice sizes should be positive, "
          "number of lattice cells should be > 0.");
    }
  }

  /// Copy-constructor
  RectangularLattice(RectangularLattice<T> const& rl)
      : lattice_(rl.lattice_),
        lattice_sizes_(rl.lattice_sizes_),
        n_cells_(rl.n_cells_),
        cell_sizes_(rl.cell_sizes_),
        cell_volume_(rl.cell_volume_),
        origin_(rl.origin_),
        periodic_(rl.periodic_),
        when_update_(rl.when_update_) {}

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
  const std::array<int, 3>& n_cells() const { return n_cells_; }

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
   * Overwrite with a template value T at a given node
   */
  void assign_value(int lattice_index, T value) {
    lattice_[lattice_index] = value;
  }

  /**
   * Return the index of a given cell.
   *
   * \param[in] ix index of a cell in the x-direction
   * \param[in] iy index of a cell in the y-direction
   * \param[in] iz index of a cell in the z-direction
   * \return index of cell
   */
  int index1d(int ix, int iy, int iz) {
    assert(ix < n_cells_[0]);
    assert(iy < n_cells_[1]);
    assert(iz < n_cells_[2]);
    return (ix + n_cells_[0] * (iy + iz * n_cells_[1]));
  }
  /**
   * Given the indices of a cell in the x, y, and z directions, return index of
   * the nearest cell in the -x direction ("left").
   *
   * \param[in] ix index of a cell in the x-direction
   * \param[in] iy index of a cell in the y-direction
   * \param[in] iz index of a cell in the z-direction
   * \return index of the nearest cell in the "left" direction (ix-1) with
   * respect to the current cell
   */
  int index_left(int ix, int iy, int iz) {
    if (unlikely(ix == 0)) {
      return periodic_ ? index1d(n_cells_[0] - 1, iy, iz) : index1d(ix, iy, iz);
    } else {
      return index1d(ix - 1, iy, iz);
    }
  }
  /**
   * Given the indices of a cell in the x, y, and z directions, return index of
   * the nearest cell in the +x direction ("right").
   *
   * \param[in] ix index of a cell in the x-direction
   * \param[in] iy index of a cell in the y-direction
   * \param[in] iz index of a cell in the z-direction
   * \return index of the nearest cell in the "right" direction (ix+1) with
   * respect to the current cell
   */
  int index_right(int ix, int iy, int iz) {
    if (unlikely(ix == (n_cells_[0] - 1))) {
      return periodic_ ? index1d(0, iy, iz) : index1d(ix, iy, iz);
    } else {
      return index1d(ix + 1, iy, iz);
    }
  }
  /**
   * Given the indices of a cell in the x, y, and z directions, return index of
   * the nearest cell in the -y direction ("down").
   *
   * \param[in] ix index of a cell in the x-direction
   * \param[in] iy index of a cell in the y-direction
   * \param[in] iz index of a cell in the z-direction
   * \return index of the nearest cell in the "down" direction (iy-1) with
   * respect to the current cell
   */
  int index_down(int ix, int iy, int iz) {
    if (unlikely(iy == 0)) {
      return periodic_ ? index1d(ix, n_cells_[1] - 1, iz) : index1d(ix, iy, iz);
    } else {
      return index1d(ix, iy - 1, iz);
    }
  }
  /**
   * Given the indices of a cell in the x, y, and z directions, return index of
   * the nearest cell in the +y direction ("up").
   *
   * \param[in] ix index of a cell in the x-direction
   * \param[in] iy index of a cell in the y-direction
   * \param[in] iz index of a cell in the z-direction
   * \return index of the nearest cell in the "up" direction (iy+1) with
   * respect to the current cell
   */
  int index_up(int ix, int iy, int iz) {
    if (unlikely(iy == (n_cells_[1] - 1))) {
      return periodic_ ? index1d(ix, 0, iz) : index1d(ix, iy, iz);
    } else {
      return index1d(ix, iy + 1, iz);
    }
  }
  /**
   * Given the indices of a cell in the x, y, and z directions, return index of
   * the nearest cell in the -z direction ("backward").
   *
   * \param[in] ix index of a cell in the x-direction
   * \param[in] iy index of a cell in the y-direction
   * \param[in] iz index of a cell in the z-direction
   * \return index of the nearest cell in the "backward" direction (iz-1) with
   * respect to the current cell
   */
  int index_backward(int ix, int iy, int iz) {
    if (unlikely(iz == 0)) {
      return periodic_ ? index1d(ix, iy, n_cells_[2] - 1) : index1d(ix, iy, iz);
    } else {
      return index1d(ix, iy, iz - 1);
    }
  }
  /**
   * Given the indices of a cell in the x, y, and z directions, return index of
   * the nearest cell in the +z direction ("forward").
   *
   * \param[in] ix index of a cell in the x-direction
   * \param[in] iy index of a cell in the y-direction
   * \param[in] iz index of a cell in the z-direction
   * \return index of the nearest cell in the "forward" direction (iz+1) with
   * respect to the current cell
   */
  int index_forward(int ix, int iy, int iz) {
    if (unlikely(iz == (n_cells_[2] - 1))) {
      return periodic_ ? index1d(ix, iy, 0) : index1d(ix, iy, iz);
    } else {
      return index1d(ix, iy, iz + 1);
    }
  }

  /**
   * Compute a gradient on a lattice of doubles via the finite difference method
   *
   * return a lattice of ThreeVectors which are gradients of the values on the
   * original lattice
   */
  void compute_gradient_lattice(
      RectangularLattice<ThreeVector>& grad_lat) const {
    if (n_cells_[0] < 2 || n_cells_[1] < 2 || n_cells_[2] < 2) {
      // Gradient calculation is impossible
      throw std::runtime_error(
          "Lattice is too small for gradient calculation"
          " (should be at least 2x2x2)");
    }
    if (!identical_to_lattice(&grad_lat)) {
      // Lattice for gradient should have identical origin/dims/periodicity
      throw std::invalid_argument(
          "Lattice for gradient should have the"
          " same origin/dims/periodicity as the original one.");
    }
    const double inv_2dx = 0.5 / cell_sizes_[0];
    const double inv_2dy = 0.5 / cell_sizes_[1];
    const double inv_2dz = 0.5 / cell_sizes_[2];
    const int dix = 1;
    const int diy = n_cells_[0];
    const int diz = n_cells_[0] * n_cells_[1];
    const int d = diz * n_cells_[2];

    for (int iz = 0; iz < n_cells_[2]; iz++) {
      const int z_offset = diz * iz;
      for (int iy = 0; iy < n_cells_[1]; iy++) {
        const int y_offset = diy * iy + z_offset;
        for (int ix = 0; ix < n_cells_[0]; ix++) {
          const int index = ix + y_offset;
          if (unlikely(ix == 0)) {
            (grad_lat)[index].set_x1(
                periodic_
                    ? (lattice_[index + dix] - lattice_[index + diy - dix]) *
                          inv_2dx
                    : (lattice_[index + dix] - lattice_[index]) * 2.0 *
                          inv_2dx);
          } else if (unlikely(ix == n_cells_[0] - 1)) {
            (grad_lat)[index].set_x1(
                periodic_
                    ? (lattice_[index - diy + dix] - lattice_[index - dix]) *
                          inv_2dx
                    : (lattice_[index] - lattice_[index - dix]) * 2.0 *
                          inv_2dx);
          } else {
            (grad_lat)[index].set_x1(
                (lattice_[index + dix] - lattice_[index - dix]) * inv_2dx);
          }

          if (unlikely(iy == 0)) {
            (grad_lat)[index].set_x2(
                periodic_
                    ? (lattice_[index + diy] - lattice_[index + diz - diy]) *
                          inv_2dy
                    : (lattice_[index + diy] - lattice_[index]) * 2.0 *
                          inv_2dy);
          } else if (unlikely(iy == n_cells_[1] - 1)) {
            (grad_lat)[index].set_x2(
                periodic_
                    ? (lattice_[index - diz + diy] - lattice_[index - diy]) *
                          inv_2dy
                    : (lattice_[index] - lattice_[index - diy]) * 2.0 *
                          inv_2dy);
          } else {
            (grad_lat)[index].set_x2(
                (lattice_[index + diy] - lattice_[index - diy]) * inv_2dy);
          }

          if (unlikely(iz == 0)) {
            (grad_lat)[index].set_x3(
                periodic_
                    ? (lattice_[index + diz] - lattice_[index + d - diz]) *
                          inv_2dz
                    : (lattice_[index + diz] - lattice_[index]) * 2.0 *
                          inv_2dz);
          } else if (unlikely(iz == n_cells_[2] - 1)) {
            (grad_lat)[index].set_x3(
                periodic_
                    ? (lattice_[index - d + diz] - lattice_[index - diz]) *
                          inv_2dz
                    : (lattice_[index] - lattice_[index - diz]) * 2.0 *
                          inv_2dz);
          } else {
            (grad_lat)[index].set_x3(
                (lattice_[index + diz] - lattice_[index - diz]) * inv_2dz);
          }
        }
      }
    }
  }

  /**
   * Compute a fourgradient on a lattice of FourVectors jmu via the finite
   * difference method.
   *
   * \param[in] old_lat the lattice of FourVectors jmu at a previous time step
   * \param[in] time_step the used time step, needed for the time derivative
   * \param[out] grad_lat a lattice of 4-arrays of 4-vectors with the following
   * structure: [djmu_dt, djmu_dx, djmu_dy, djmu_dz]
   */
  void compute_four_gradient_lattice(
      RectangularLattice<FourVector>& old_lat, double time_step,
      RectangularLattice<std::array<FourVector, 4>>& grad_lat) const {
    if (n_cells_[0] < 2 || n_cells_[1] < 2 || n_cells_[2] < 2) {
      // Gradient calculation is impossible
      throw std::runtime_error(
          "Lattice is too small for gradient calculation"
          " (should be at least 2x2x2)");
    }
    if (!identical_to_lattice(&grad_lat)) {
      // Lattice for gradient should have identical origin/dims/periodicity
      throw std::invalid_argument(
          "Lattice for gradient should have the"
          " same origin/dims/periodicity as the original one.");
    }
    const double inv_2dx = 0.5 / cell_sizes_[0];
    const double inv_2dy = 0.5 / cell_sizes_[1];
    const double inv_2dz = 0.5 / cell_sizes_[2];
    const int dix = 1;
    const int diy = n_cells_[0];
    const int diz = n_cells_[0] * n_cells_[1];
    const int d = diz * n_cells_[2];

    for (int iz = 0; iz < n_cells_[2]; iz++) {
      const int z_offset = diz * iz;
      for (int iy = 0; iy < n_cells_[1]; iy++) {
        const int y_offset = diy * iy + z_offset;
        for (int ix = 0; ix < n_cells_[0]; ix++) {
          const int index = ix + y_offset;

          // auxiliary vectors used for the calculation of gradients
          FourVector grad_t_jmu(0.0, 0.0, 0.0, 0.0);
          FourVector grad_x_jmu(0.0, 0.0, 0.0, 0.0);
          FourVector grad_y_jmu(0.0, 0.0, 0.0, 0.0);
          FourVector grad_z_jmu(0.0, 0.0, 0.0, 0.0);
          // t direction
          grad_t_jmu = (lattice_[index] - (old_lat)[index]) * (1.0 / time_step);
          // x direction
          if (unlikely(ix == 0)) {
            grad_x_jmu =
                periodic_
                    ? (lattice_[index + dix] - lattice_[index + diy - dix]) *
                          inv_2dx
                    : (lattice_[index + dix] - lattice_[index]) * 2.0 * inv_2dx;
          } else if (unlikely(ix == n_cells_[0] - 1)) {
            grad_x_jmu =
                periodic_
                    ? (lattice_[index - diy + dix] - lattice_[index - dix]) *
                          inv_2dx
                    : (lattice_[index] - lattice_[index - dix]) * 2.0 * inv_2dx;
          } else {
            grad_x_jmu =
                (lattice_[index + dix] - lattice_[index - dix]) * inv_2dx;
          }
          // y direction
          if (unlikely(iy == 0)) {
            grad_y_jmu =
                periodic_
                    ? (lattice_[index + diy] - lattice_[index + diz - diy]) *
                          inv_2dy
                    : (lattice_[index + diy] - lattice_[index]) * 2.0 * inv_2dy;
          } else if (unlikely(iy == n_cells_[1] - 1)) {
            grad_y_jmu =
                periodic_
                    ? (lattice_[index - diz + diy] - lattice_[index - diy]) *
                          inv_2dy
                    : (lattice_[index] - lattice_[index - diy]) * 2.0 * inv_2dy;
          } else {
            grad_y_jmu =
                (lattice_[index + diy] - lattice_[index - diy]) * inv_2dy;
          }
          // z direction
          if (unlikely(iz == 0)) {
            grad_z_jmu =
                periodic_
                    ? (lattice_[index + diz] - lattice_[index + d - diz]) *
                          inv_2dz
                    : (lattice_[index + diz] - lattice_[index]) * 2.0 * inv_2dz;
          } else if (unlikely(iz == n_cells_[2] - 1)) {
            grad_z_jmu =
                periodic_
                    ? (lattice_[index - d + diz] - lattice_[index - diz]) *
                          inv_2dz
                    : (lattice_[index] - lattice_[index - diz]) * 2.0 * inv_2dz;
          } else {
            grad_z_jmu =
                (lattice_[index + diz] - lattice_[index - diz]) * inv_2dz;
          }
          // fill
          (grad_lat)[index][0] = grad_t_jmu;
          (grad_lat)[index][1] = grad_x_jmu;
          (grad_lat)[index][2] = grad_y_jmu;
          (grad_lat)[index][3] = grad_z_jmu;
        }
      }
    }
  }

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
    const int ix =
        numeric_cast<int>(std::floor((r.x1() - origin_[0]) / cell_sizes_[0]));
    const int iy =
        numeric_cast<int>(std::floor((r.x2() - origin_[1]) / cell_sizes_[1]));
    const int iz =
        numeric_cast<int>(std::floor((r.x3() - origin_[2]) / cell_sizes_[2]));
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
    logg[LLattice].debug(
        "Iterating sublattice with lower bound index (", lower_bounds[0], ",",
        lower_bounds[1], ",", lower_bounds[2], "), upper bound index (",
        upper_bounds[0], ",", upper_bounds[1], ",", upper_bounds[2], ")");

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
   * Iterates only nodes whose cell centers lie not further than r_cut in x, y,
   * z directions from the given point, that is iterates within a cube of side
   * length 2*r_cut, and applies a function to each node.
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
  void iterate_in_cube(const ThreeVector& point, const double r_cut, F&& func) {
    std::array<int, 3> l_bounds, u_bounds;

    /* Array holds value at the cell center: r_center = r_0 + (i+0.5)cell_size,
     * where i is index in any direction. Therefore we want cells with condition
     * (r-r_cut)/csize - 0.5 < i < (r+r_cut)/csize - 0.5, r = r_center - r_0 */
    for (int i = 0; i < 3; i++) {
      l_bounds[i] = numeric_cast<int>(
          std::ceil((point[i] - origin_[i] - r_cut) / cell_sizes_[i] - 0.5));
      u_bounds[i] = numeric_cast<int>(
          std::ceil((point[i] - origin_[i] + r_cut) / cell_sizes_[i] - 0.5));
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
   * Calculate a volume integral with given integrand
   *
   * \tparam F return type of the integrand
   * \param[out] integral variable to store rsult in.
   *             Input should be 0.
   * \param[in] integrand Function to be integrated
   * \param[in] rcut size of the integration volume. In total the intgration
   *                 volume will be a cube with edge length 2*rcut
   * \param[in] point center of the integration volume
   **/
  template <typename F>
  void integrate_volume(F& integral,
                        F (*integrand)(ThreeVector, T&, ThreeVector),
                        const double rcut, const ThreeVector& point) {
    iterate_in_rectangle(
        point, {rcut, rcut, rcut},
        [&point, &integral, &integrand, this](T value, int ix, int iy, int iz) {
          ThreeVector pos = this->cell_center(ix, iy, iz);
          integral += integrand(pos, value, point) * this->cell_volume_;
        });
  }
  /**
   * Iterates only nodes whose cell centers lie not further than d_x in x-,
   * d_y in y-, and d_z in z-direction from the given point, that is iterates
   * within a rectangle of side lengths (2*dx, 2*dy, 2*dz), and applies a
   * function to each node.
   * Useful for adding quantities from one particle to the lattice.
   *
   * \tparam F Type of the function. Arguments are the current node and the 3
   * integer indices of the cell.
   * \param[in] point Position, usually the position of particle [fm].
   * \param[in] rectangle Maximum distances in the x-, y-, and z-directions
   * from the cell center to the given position. [fm]
   * \param[in] func Function acting on the cells (such as taking value).
   */
  template <typename F>
  void iterate_in_rectangle(const ThreeVector& point,
                            const std::array<double, 3>& rectangle, F&& func) {
    std::array<int, 3> l_bounds, u_bounds;

    /* Array holds value at the cell center: r_center = r_0 + (i+0.5)cell_size,
     * where i is index in any direction. Therefore we want cells with condition
     * (r[i]-rectangle[i])/csize - 0.5 < i < (r[i]+rectangle[i])/csize - 0.5,
     * r[i] = r_center[i] - r_0[i]
     */
    for (int i = 0; i < 3; i++) {
      l_bounds[i] = numeric_cast<int>(std::ceil(
          (point[i] - origin_[i] - rectangle[i]) / cell_sizes_[i] - 0.5));
      u_bounds[i] = numeric_cast<int>(std::ceil(
          (point[i] - origin_[i] + rectangle[i]) / cell_sizes_[i] - 0.5));
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
   * Iterates only over nodes corresponding to the center cell (the cell
   * containing the given point) and its nearest neighbors in the -x, +x, -y,
   * +y, -z, +z directions, and applies a function to each node.
   * Useful for adding quantities from one particle to the lattice.
   *
   * \tparam F Type of the function. Arguments are the current node and the 3
   * integer indices of the cell.
   * \param[in] point Position, usually the position of particle [fm].
   * \param[in] func Function acting on the cells (such as taking value).
   */
  template <typename F>
  void iterate_nearest_neighbors(const ThreeVector& point, F&& func) {
    // get the 3D indices of the cell containing the given point
    const int ix = numeric_cast<int>(
        std::floor((point.x1() - origin_[0]) / cell_sizes_[0]));
    const int iy = numeric_cast<int>(
        std::floor((point.x2() - origin_[1]) / cell_sizes_[1]));
    const int iz = numeric_cast<int>(
        std::floor((point.x3() - origin_[2]) / cell_sizes_[2]));

    logg[LLattice].debug(
        "Iterating over nearest neighbors of the cell at ix = ", ix,
        ", iy = ", iy, ", iz = ", iz);

    // determine the 1D index of the center cell
    const int i = index1d(ix, iy, iz);
    func(lattice_[i], i, i);

    // determine the indeces of nearby cells, perform function on them
    int ileft = index_left(ix, iy, iz);
    func(lattice_[ileft], ileft, i);

    int iright = index_right(ix, iy, iz);
    func(lattice_[iright], iright, i);

    int idown = index_down(ix, iy, iz);
    func(lattice_[idown], idown, i);

    int iup = index_up(ix, iy, iz);
    func(lattice_[iup], iup, i);

    int ibackward = index_backward(ix, iy, iz);
    func(lattice_[ibackward], ibackward, i);

    int iforward = index_forward(ix, iy, iz);
    func(lattice_[iforward], iforward, i);
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
    return n_cells_[0] == lat->n_cells()[0] &&
           n_cells_[1] == lat->n_cells()[1] &&
           n_cells_[2] == lat->n_cells()[2] &&
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
  /// Volume of a cell.
  const double cell_volume_;
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

#endif  // SRC_INCLUDE_SMASH_LATTICE_H_
