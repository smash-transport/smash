/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_GRID_H_
#define SRC_INCLUDE_GRID_H_

#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "forwarddeclarations.h"
#include "particles.h"

namespace Smash {

/**
 * Identifies the mode of the Grid.
 */
enum class GridOptions : char {
  /// without ghost cells
  Normal = 0,
  /// with ghost cells for periodic boundaries
  PeriodicBoundaries = 1
};

/**
 * Indentifies the strategy of determining the cell size.
 */
enum class CellSizeStrategy : char {

  /// Look for optimal cell size.
  Optimal,

  /// Make cells as large as possible.
  ///
  /// This means a single cell for normal boundary conditions and 8 cells
  /// for periodic boundary conditions.
  Largest
};

/**
 * Base class for Grid to host common functions that do not depend on the
 * GridOptions parameter.
 */
class GridBase {
 public:
  typedef int size_type;

  /**
   * The minimum cell length for the given testparticles (defaults to 1).
   */
  static float min_cell_length(int testparticles = 1) {
    // 2.5 fm corresponds to maximal cross-section of 200 mb = 20 fm^2
    // sqrt(20 fm^2/N_{test}/pi) is approximately 2.5/sqrt(N_{test})
    return 2.5f / std::sqrt(static_cast<float>(testparticles));
  }

 protected:
  /**
   * Returns the minimum x,y,z coordinates and the largest dx,dy,dz distances of
   * the particles in \p particles.
   */
  static std::pair<std::array<float, 3>, std::array<float, 3>>
  find_min_and_length(const Particles &particles);
};

/**
 * Abstracts a list of cells that partition the particles in the experiment into
 * regions of space that can interact / cannot interact.
 *
 * This class is used to construct a helper data structure to reduce the
 * combinatorics of finding particle pairs that could interact (scatter). It
 * takes a list of ParticleData objects and sorts them in such a way that it is
 * easy to look only at lists of particles that have a chance of interacting.
 *
 * \tparam Options This policy parameter determines whether ghost cells are
 * created to support periodic boundaries, or not.
 */
template <GridOptions Options = GridOptions::Normal>
class Grid : public GridBase {
 public:
  /**
   * Constructs a grid from the given particle list \p particles. It
   * automatically determines the necessary size for the grid from the positions
   * of the particles.
   *
   * \param particles The particles to place onto the grid.
   * \param testparticles Number of testparticles used in this event
   * \param strategy The strategy for determining the cell size
   */
  Grid(const Particles &particles, int testparticles,
       CellSizeStrategy strategy = CellSizeStrategy::Optimal)
      : Grid{find_min_and_length(particles), std::move(particles),
             testparticles, strategy} {}

  /**
   * Constructs a grid with the given minimum grid coordinates and grid length.
   * If you need periodic boundaries you have to use this constructor to set the
   * correct length to use for wrapping particles around the borders.
   *
   * \param min_and_length A pair consisting of the three min coordinates and
   * the three lengths.
   * \param particles The particles to place onto the grid.
   * \param testparticles Number of testparticles used in this event
   * \param strategy The strategy for determining the cell size
   */
  Grid(const std::pair<std::array<float, 3>, std::array<float, 3>> &
           min_and_length,
       const Particles &particles, int testparticles,
       CellSizeStrategy strategy = CellSizeStrategy::Optimal)
      : min_position_(min_and_length.first), length_(min_and_length.second) {
    /**
     * This normally equals 1/max_interaction_length, but if the number of cells
     * is reduced (because of low density) then this value is smaller.
     */
    std::array<float, 3> index_factor;

    switch (strategy) {
      case CellSizeStrategy::Optimal:
        std::tie(index_factor, number_of_cells_) =
            determine_cell_sizes(particles.size(), length_, testparticles);
        break;
      case CellSizeStrategy::Largest:
        // set number of cells
        if (Options == GridOptions::PeriodicBoundaries) {
          number_of_cells_ = {2, 2, 2};
        } else {
          number_of_cells_ = {1, 1, 1};
        }

        // set index factor
        for (std::size_t i = 0; i < index_factor.size(); ++i) {
          index_factor[i] = number_of_cells_[i] / length_[i];
          while (index_factor[i] * length_[i] >= number_of_cells_[i]) {
            index_factor[i] = std::nextafter(index_factor[i], 0.f);
          }
          assert(index_factor[i] * length_[i] < number_of_cells_[i]);
        }
        break;
      default:
        throw std::domain_error("Unhandled CellSizeStrategy!");
    }

    build_cells(index_factor, particles);
  }

  /**
   * Calculates the factor that, if multiplied with a x/y/z
   * coordinate, yields the 3-dim cell index and the required number of cells
   * (without ghost cells).
   *
   * \return A tuple of two 3-dim values:
   * \li The first tuple entry stores the conversion factors to turn a
   * normalized coordinate (particle position minus minimum position) into a
   * cell index.
   * \li The second tuple entry stores the dimensions of the grid (i.e. the
   * number of cells in each spatial direction).
   *
   * \param particle_count The number of particles to be placed in the grid.
   * \param length         Three lengths that identify the total dimensions of
   *                       the grid.
   * \param testparticles  The number of testparticles is used to scale the cell
   *                       size down, since the interaction length is also
   *                       reduced.
   */
  static std::tuple<std::array<float, 3>, std::array<int, 3>>
      determine_cell_sizes(size_type particle_count,
                           const std::array<float, 3> &length,
                           const int testparticles);

  /**
   * Iterates over all cells in the grid and calls the callback arguments with a
   * search cell and 0 to 13 neighbor cells.
   *
   * The neighbor cells are constructed like this:
   * - one cell at x+1
   * - three cells (x-1,x,x+1) at y+1
   * - nine cells (x-1, y-1)...(x+1, y+1) at z+1
   *
   * \param search_cell_callback A callable called for/with every non-empty cell
   *                             in the grid.
   * \param neighbor_cell_callback A callable called for/with every non-empty
   *                               cell and adjacent cell combination. For a
   *                               periodic grid, the first argument will be
   *                               adjusted to wrap around the grid.
   */
  void iterate_cells(
      const std::function<void(const ParticleList &)> &search_cell_callback,
      const std::function<void(const ParticleList &, const ParticleList &)> &
          neighbor_cell_callback) const;

 private:
  /**
   * Allocates the cells and fills them with ParticleData objects.
   *
   * This is different for the Normal and PeriodicBoundaries cases.
   */
  void build_cells(const std::array<float, 3> &index_factor,
                   const Particles &particles);

  /**
   * Returns the one-dimensional cell-index from the 3-dim index \p x, \p y, \p
   * z.
   */
  size_type make_index(size_type x, size_type y, size_type z) const;

  /**
   * Returns the one-dimensional cell-index from the 3-dim index \p idx.
   * This is a convenience overload for the above function.
   */
  size_type make_index(std::array<size_type, 3> idx) const {
    return make_index(idx[0], idx[1], idx[2]);
  }

  /// The lower bound of the cell coordinates.
  const std::array<float, 3> min_position_;

  /// The 3 lengths of the complete grid. Used for periodic boundary wrapping.
  const std::array<float, 3> length_;

  /// The number of cells in x, y, and z direction.
  std::array<int, 3> number_of_cells_;

  /// The cell storage.
  std::vector<ParticleList> cells_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_GRID_H_
