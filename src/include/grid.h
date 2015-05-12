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
#include <utility>
#include <vector>

#include "forwarddeclarations.h"

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
 * Base class for Grid to host common functions that do not depend on the
 * GridOptions parameter.
 */
class GridBase {
 public:
  typedef int size_type;

  /**
   * The minimum cell length for the given testparticles (defaults to 1).
   */
  static constexpr float min_cell_length(int testparticles = 1) {
    // 2.5 fm corresponds to maximal cross-section of 200 mb = 20 fm^2
    // sqrt(20 fm^2/N_{test}/pi) is approximately 2.5/sqrt(N_{test})
    return 2.5f / std::sqrt(static_cast<float>(testparticles));
  }

 protected:
  /**
   * Returns the minimum x,y,z coordinates and the largest dx,dy,dz distances of
   * the particles in \p all_particles.
   */
  static std::pair<std::array<float, 3>, std::array<float, 3>>
      find_min_and_length(const ParticleList &all_particles);
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
   * Constructs a grid from the given particle list \p all_particles. It
   * automatically determines the necessary size for the grid from the positions
   * of the particles.
   *
   * \param all_particles The particles to place onto the grid.
   * \param testparticles Number of testparticles used in this event
   */
  Grid(ParticleList &&all_particles, const int testparticles)
      : Grid{find_min_and_length(all_particles), std::move(all_particles),
             testparticles} {}

  /**
   * Constructs a grid with the given minimum grid coordinates and grid length.
   * If you need periodic boundaries you have to use this constructor to set the
   * correct length to use for wrapping particles around the borders.
   *
   * \param min_and_length A pair consisting of the three min coordinates and
   * the three lengths.
   * \param all_particles The particles to place onto the grid.
   * \param testparticles Number of testparticles used in this event
   */
  Grid(const std::pair<std::array<float, 3>, std::array<float, 3>> &
           min_and_length,
       ParticleList &&all_particles, const int testparticles)
      : min_position_(min_and_length.first) {
    const auto &length = min_and_length.second;

    std::tie(index_factor_, number_of_cells_) =
        determine_cell_sizes(all_particles.size(), length, testparticles);

    build_cells(std::move(all_particles), length);
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
   * Iterates over all cells in the grid and calls \p call_finder with a search
   * cell and 0 to 13 neighbor cells.
   *
   * The neighbor cells are constructed like this:
   * - one cell at x+1
   * - three cells (x-1,x,x+1) at y+1
   * - nine cells (x-1, y-1)...(x+1, y+1) at z+1
   *
   * \param call_finder A lambda (or other functor) that is called as often as
   * there are search cells. Experiment uses it to call the Action finders from
   * it.
   */
  void iterate_cells(const std::function<
      void(const ParticleList &, const std::vector<const ParticleList *> &)> &
                         call_finder) const;

 private:
  /**
   * Allocates the cells and fills them with ParticleData objects.
   *
   * This is different for the Normal and PeriodicBoundaries cases.
   */
  void build_cells(ParticleList &&all_particles,
                   const std::array<float, 3> &length);

  /**
   * Returns the one-dimensional cell-index from the 3-dim index \p x, \p y, \p
   * z.
   */
  size_type make_index(size_type x, size_type y, size_type z) const;

  /**
   * Returns the one-dimensional cell-index from the position vector inside the
   * grid.
   *
   * In Normal mode this simply calculates the distance to min_position_ and
   * multiplies it with index_factor_ to determine the 3 x,y,z indexes to pass
   * to the make_index overload above.
   *
   * In PeriodicBoundaries mode the x and y indexes are incremented by one to
   * adjust for the ghost cells.
   */
  size_type make_index(const ThreeVector &position) const;

  /// The lower bound of the cell coordinates.
  const std::array<float, 3> min_position_;

  /**
   * This normally equals 1/max_interaction_length, but if the number of cells
   * is reduced (because of low density) then this value is smaller.
   */
  std::array<float, 3> index_factor_;

  /// The number of cells in x, y, and z direction.
  std::array<int, 3> number_of_cells_;

  /// The cell storage.
  std::vector<ParticleList> cells_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_GRID_H_
