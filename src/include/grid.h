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

#include <functional>
#include <vector>
#include <array>

#include "forwarddeclarations.h"

namespace Smash {

enum class GridOptions : char {
  Normal = 0,
  PeriodicBoundaries = 1
};

class GridBase {
 public:
  typedef int size_type;

 protected:
  static constexpr std::array<float, 3> max_interaction_length = {
      {2.5f, 2.5f, 2.5f}};

  static std::pair<std::array<float, 3>, std::array<float, 3>>
      find_min_and_length(const ParticleList &all_particles);

  static std::tuple<std::array<float, 3>, std::array<int, 3>>
      determine_cell_sizes(size_type particle_count,
                           const std::array<float, 3> &length);

};

template <GridOptions Options = GridOptions::Normal>
class Grid : public GridBase {
 public:
  /**
   * Constructs a grid from the given particle list \p all_particles. It
   * automatically determines the necessary size for the grid from the positions
   * of the particles.
   */
  Grid(ParticleList &&all_particles)
      : Grid{find_min_and_length(all_particles), std::move(all_particles)} {}

  /**
   * Constructs a grid with the given minimum grid coordinates and grid length.
   * If you need periodic boundaries you have to use this constructor to set the
   * correct length to use for wrapping particles around the borders.
   */
  Grid(const std::pair<std::array<float, 3>, std::array<float, 3>> &
           min_and_length,
       ParticleList &&all_particles)
      : min_position_(min_and_length.first) {
    const auto &length = min_and_length.second;

    std::tie(index_factor_, number_of_cells_) =
        determine_cell_sizes(all_particles.size(), length);

    build_cells(std::move(all_particles), length);
  }

  /**
   * Iterates over all cells in the grid and calls \p call_finder with a search
   * cell and 0 to 13 neighbor cells.
   */
  void iterate_cells(const std::function<
      void(const ParticleList &, const std::vector<const ParticleList *> &)> &
                         call_finder) const;

 private:
  /**
   * Allocates the cells and fills them with ParticleData objects.
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
   * grid
   */
  size_type make_index(const ThreeVector &position) const;

  /**
   * Returns whether the cell at the given 3-dim index \p x, \p y, \p z is a
   * ghost cell.
   */
  bool is_ghost_cell(size_type x, size_type y, size_type z) const;

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
