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

#include <algorithm>
#include <vector>
#include <array>
#include <assert.h>

#include "forwarddeclarations.h"
#include "fourvector.h"
#include "logging.h"
#include "particledata.h"
#include "threevector.h"

namespace std {
template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  auto column = out.tellp();
  out << "vector{";
  for (const auto &x : v) {
    if (out.tellp() - column >= 100) {
      out << '\n';
      column = out.tellp();
    }
    out << x << ' ';
  }
  return out << '}';
}

template <typename T, std::size_t N>
std::ostream &operator<<(std::ostream &out, const std::array<T, N> &a) {
  out << "array{";
  for (const auto &x : a) {
    out << x << ' ';
  }
  return out << '}';
}
}  // namespace std

namespace Smash {

enum class GridOptions : char {
  Normal = 0,
  PeriodicBoundaries = 1
};

template <GridOptions Options = GridOptions::Normal>
class Grid {
  static constexpr std::array<float, 3> max_interaction_length = {
      {2.5f, 2.5f, 2.5f}};

 public:
  template <typename T>
  Grid(const T &all_particles){
    const auto &log = logger<LogArea::Grid>();

    const auto particle_count = all_particles.size();
    assert(particle_count > 0);

    // intialize min and max position arrays with the position of the first
    // particle in the list
    const auto first = all_particles.begin()->position().threevec();
    min_position_ = {{static_cast<float>(first[0]),
                      static_cast<float>(first[1]),
                      static_cast<float>(first[2])}};
    auto max_position = min_position_;
    for (const auto &p : all_particles) {
      const auto pos = p.position().threevec();
      min_position_[0] = std::min(min_position_[0], static_cast<float>(pos[0]));
      min_position_[1] = std::min(min_position_[1], static_cast<float>(pos[1]));
      min_position_[2] = std::min(min_position_[2], static_cast<float>(pos[2]));
      max_position [0] = std::max(max_position [0], static_cast<float>(pos[0]));
      max_position [1] = std::max(max_position [1], static_cast<float>(pos[1]));
      max_position [2] = std::max(max_position [2], static_cast<float>(pos[2]));
    }

    // The number of cells is determined by the min and max coordinates where
    // particles are positioned and the maximal interaction length (which equals
    // the length of a cell).
    // But don't let the number of cells exceed the actual number of particles.
    // That would be overkill. Let max_cells³ ≤ particle_count (conversion to
    // int truncates).
    const int max_cells = std::cbrt(particle_count);
    for (std::size_t i = 0; i < number_of_cells_.size(); ++i) {
      index_factor_[i] = 1.f / max_interaction_length[i];
      number_of_cells_[i] =
          std::ceil((max_position[i] - min_position_[i]) * index_factor_[i]);
      if (number_of_cells_[i] > max_cells) {
        number_of_cells_[i] = max_cells;
        index_factor_[i] = (max_cells - 0.1f)  // -0.1 for safety margin
                           / (max_position[i] - min_position_[i]);
      }
    }

    log.debug("min: ", min_position_, "\nmax: ", max_position, "\ncells: ",
              number_of_cells_, "\ninteraction length: ",
              max_interaction_length, "\nindex_factor: ", index_factor_);

    // After the grid parameters are determined, we can start placing the
    // particles in cells.
    cells_.resize(number_of_cells_[0] * number_of_cells_[1] *
                  number_of_cells_[2]);

    for (const auto &p : all_particles) {
      const auto idx = make_index(p.position().threevec());
      assert(idx < cells_.size());
      cells_[idx].push_back(p);
    }
    log.debug(cells_);
  }

  template <typename F>
  void iterate_cells(F &&closure) {
    std::vector<const ParticleList *> neighbors;
    neighbors.reserve(13);
    auto &&
    build_neighbors = [&](const std::initializer_list<std::size_t> &indexes)
                          -> const std::vector<const ParticleList *> &
    {
      neighbors.clear();
      for (const auto i : indexes) {
        neighbors.emplace_back(&cells_[i]);
      }
      return neighbors;
    };
    if (Options == GridOptions::PeriodicBoundaries) {
    } else {
      for (std::size_t z = 0; z < number_of_cells_[2] - 1; ++z) {
        closure(cells_[make_index(0, 0, z)], build_neighbors({
                make_index(1, 0, z + 0),
                make_index(0, 1, z + 0),
                make_index(1, 1, z + 0),
                make_index(0, 0, z + 1),
                make_index(1, 0, z + 1),
                make_index(0, 1, z + 1),
                make_index(1, 1, z + 1)
              }));
        for (std::size_t x = 1; x < number_of_cells_[0] - 1; ++x) {
          closure(cells_[make_index(x, 0, z)], build_neighbors({
                  make_index(x + 1, 0, z + 0),
                  make_index(x - 1, 1, z + 0),
                  make_index(x + 0, 1, z + 0),
                  make_index(x + 1, 1, z + 0),
                  make_index(x - 1, 0, z + 1),
                  make_index(x + 0, 0, z + 1),
                  make_index(x + 1, 0, z + 1),
                  make_index(x - 1, 1, z + 1),
                  make_index(x + 0, 1, z + 1),
                  make_index(x + 1, 1, z + 1)
                }));
        }
        {
          const std::size_t x = number_of_cells_[0] - 1;
          closure(cells_[make_index(x, 0, z)], build_neighbors({
                  make_index(x - 1, 1, z + 0),
                  make_index(x + 0, 1, z + 0),
                  make_index(x - 1, 0, z + 1),
                  make_index(x + 0, 0, z + 1),
                  make_index(x - 1, 1, z + 1),
                  make_index(x + 0, 1, z + 1)
                }));
        }
        for (std::size_t y = 1; y < number_of_cells_[1] - 1; ++y) {
          closure(cells_[make_index(0, y, z)], build_neighbors({
                  make_index(1, y + 0, z + 0),
                  make_index(0, y + 1, z + 0),
                  make_index(1, y + 1, z + 0),
                  make_index(0, y - 1, z + 1),
                  make_index(1, y - 1, z + 1),
                  make_index(0, y + 0, z + 1),
                  make_index(1, y + 0, z + 1),
                  make_index(0, y + 1, z + 1),
                  make_index(1, y + 1, z + 1)
                }));
          for (std::size_t x = 1; x < number_of_cells_[0] - 1; ++x) {
            neighbors.clear();
            neighbors.emplace_back(&cells_[make_index(x + 1, y + 0, z + 0)]);
            neighbors.emplace_back(&cells_[make_index(x - 1, y + 1, z + 0)]);
            neighbors.emplace_back(&cells_[make_index(x + 0, y + 1, z + 0)]);
            neighbors.emplace_back(&cells_[make_index(x + 1, y + 1, z + 0)]);
            neighbors.emplace_back(&cells_[make_index(x - 1, y - 1, z + 1)]);
            neighbors.emplace_back(&cells_[make_index(x + 0, y - 1, z + 1)]);
            neighbors.emplace_back(&cells_[make_index(x + 1, y - 1, z + 1)]);
            neighbors.emplace_back(&cells_[make_index(x - 1, y + 0, z + 1)]);
            neighbors.emplace_back(&cells_[make_index(x + 0, y + 0, z + 1)]);
            neighbors.emplace_back(&cells_[make_index(x + 1, y + 0, z + 1)]);
            neighbors.emplace_back(&cells_[make_index(x - 1, y + 1, z + 1)]);
            neighbors.emplace_back(&cells_[make_index(x + 0, y + 1, z + 1)]);
            neighbors.emplace_back(&cells_[make_index(x + 1, y + 1, z + 1)]);
            closure(cells_[make_index(x, y, z)], neighbors);
          }
          const std::size_t x = number_of_cells_[0] - 1;
          neighbors.clear();
          neighbors.emplace_back(&cells_[make_index(x - 1, y + 1, z + 0)]);
          neighbors.emplace_back(&cells_[make_index(x + 0, y + 1, z + 0)]);
          neighbors.emplace_back(&cells_[make_index(x - 1, y - 1, z + 1)]);
          neighbors.emplace_back(&cells_[make_index(x + 0, y - 1, z + 1)]);
          neighbors.emplace_back(&cells_[make_index(x - 1, y + 0, z + 1)]);
          neighbors.emplace_back(&cells_[make_index(x + 0, y + 0, z + 1)]);
          neighbors.emplace_back(&cells_[make_index(x - 1, y + 1, z + 1)]);
          neighbors.emplace_back(&cells_[make_index(x + 0, y + 1, z + 1)]);
          closure(cells_[make_index(x, y, z)], neighbors);
        }
        const std::size_t y = number_of_cells_[1] - 1;
        neighbors.clear();
        neighbors.emplace_back(&cells_[make_index(1, y + 0, z + 0)]);
        neighbors.emplace_back(&cells_[make_index(0, y - 1, z + 1)]);
        neighbors.emplace_back(&cells_[make_index(1, y - 1, z + 1)]);
        neighbors.emplace_back(&cells_[make_index(0, y + 0, z + 1)]);
        neighbors.emplace_back(&cells_[make_index(1, y + 0, z + 1)]);
        closure(cells_[make_index(0, y, z)], neighbors);
        for (std::size_t x = 1; x < number_of_cells_[0] - 1; ++x) {
          neighbors.clear();
          neighbors.emplace_back(&cells_[make_index(x + 1, y + 0, z + 0)]);
          neighbors.emplace_back(&cells_[make_index(x - 1, y - 1, z + 1)]);
          neighbors.emplace_back(&cells_[make_index(x + 0, y - 1, z + 1)]);
          neighbors.emplace_back(&cells_[make_index(x + 1, y - 1, z + 1)]);
          neighbors.emplace_back(&cells_[make_index(x - 1, y + 0, z + 1)]);
          neighbors.emplace_back(&cells_[make_index(x + 0, y + 0, z + 1)]);
          neighbors.emplace_back(&cells_[make_index(x + 1, y + 0, z + 1)]);
          closure(cells_[make_index(x, y, z)], neighbors);
        }
        const std::size_t x = number_of_cells_[0] - 1;
        neighbors.clear();
        neighbors.emplace_back(&cells_[make_index(x - 1, y - 1, z + 1)]);
        neighbors.emplace_back(&cells_[make_index(x + 0, y - 1, z + 1)]);
        neighbors.emplace_back(&cells_[make_index(x - 1, y + 0, z + 1)]);
        neighbors.emplace_back(&cells_[make_index(x + 0, y + 0, z + 1)]);
        closure(cells_[make_index(x, y, z)], neighbors);
      }
      const std::size_t z = number_of_cells_[2] - 1;
      neighbors.clear();
      neighbors.emplace_back(&cells_[make_index(1, 0, z + 0)]);
      neighbors.emplace_back(&cells_[make_index(0, 1, z + 0)]);
      neighbors.emplace_back(&cells_[make_index(1, 1, z + 0)]);
      closure(cells_[make_index(0, 0, z)], neighbors);
      for (std::size_t x = 1; x < number_of_cells_[0] - 1; ++x) {
        neighbors.clear();
        neighbors.emplace_back(&cells_[make_index(x + 1, 0, z + 0)]);
        neighbors.emplace_back(&cells_[make_index(x - 1, 1, z + 0)]);
        neighbors.emplace_back(&cells_[make_index(x + 0, 1, z + 0)]);
        neighbors.emplace_back(&cells_[make_index(x + 1, 1, z + 0)]);
        closure(cells_[make_index(x, 0, z)], neighbors);
      }
      {
        const std::size_t x = number_of_cells_[0] - 1;
        neighbors.clear();
        neighbors.emplace_back(&cells_[make_index(x - 1, 1, z + 0)]);
        neighbors.emplace_back(&cells_[make_index(x + 0, 1, z + 0)]);
        closure(cells_[make_index(x, 0, z)], neighbors);
      }
      for (std::size_t y = 1; y < number_of_cells_[1] - 1; ++y) {
        neighbors.clear();
        neighbors.emplace_back(&cells_[make_index(1, y + 0, z + 0)]);
        neighbors.emplace_back(&cells_[make_index(0, y + 1, z + 0)]);
        neighbors.emplace_back(&cells_[make_index(1, y + 1, z + 0)]);
        closure(cells_[make_index(0, y, z)], neighbors);
        for (std::size_t x = 1; x < number_of_cells_[0] - 1; ++x) {
          neighbors.clear();
          neighbors.emplace_back(&cells_[make_index(x + 1, y + 0, z + 0)]);
          neighbors.emplace_back(&cells_[make_index(x - 1, y + 1, z + 0)]);
          neighbors.emplace_back(&cells_[make_index(x + 0, y + 1, z + 0)]);
          neighbors.emplace_back(&cells_[make_index(x + 1, y + 1, z + 0)]);
          closure(cells_[make_index(x, y, z)], neighbors);
        }
        const std::size_t x = number_of_cells_[0] - 1;
        neighbors.clear();
        neighbors.emplace_back(&cells_[make_index(x - 1, y + 1, z + 0)]);
        neighbors.emplace_back(&cells_[make_index(x + 0, y + 1, z + 0)]);
        closure(cells_[make_index(x, y, z)], neighbors);
      }
      const std::size_t y = number_of_cells_[1] - 1;
      neighbors.clear();
      neighbors.emplace_back(&cells_[make_index(1, y + 0, z + 0)]);
      closure(cells_[make_index(0, y, z)], neighbors);
      for (std::size_t x = 1; x < number_of_cells_[0] - 1; ++x) {
        neighbors.clear();
        neighbors.emplace_back(&cells_[make_index(x + 1, y + 0, z + 0)]);
        closure(cells_[make_index(x, y, z)], neighbors);
      }
      const std::size_t x = number_of_cells_[0] - 1;
      neighbors.clear();
      closure(cells_[make_index(x, y, z)], neighbors);
    }
  }

 private:
  std::size_t make_index(std::size_t x, std::size_t y, std::size_t z) {
    return (z * number_of_cells_[1] + y) * number_of_cells_[0] + x;
  }

  std::size_t make_index(const ThreeVector &position) {
    return make_index(
        std::floor((static_cast<float>(position[0]) - min_position_[0]) *
                   index_factor_[0]),
        std::floor((static_cast<float>(position[1]) - min_position_[1]) *
                   index_factor_[1]),
        std::floor((static_cast<float>(position[2]) - min_position_[2]) *
                   index_factor_[2]));
  }

  std::array<float, 3> min_position_;

  /**
   * This normally equals 1/max_interaction_length, but if the number of cells
   * is reduced (because of low density) then this value is smaller.
   */
  std::array<float, 3> index_factor_;

  std::array<int, 3> number_of_cells_;
  std::vector<ParticleList> cells_;
};

template <GridOptions Options>
constexpr std::array<float, 3> Grid<Options>::max_interaction_length;

}  // namespace Smash

#endif  // SRC_INCLUDE_GRID_H_
