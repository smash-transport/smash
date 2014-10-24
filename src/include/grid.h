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
#include <stdexcept>
#include <vector>
#include <array>
#include <assert.h>

#include "algorithms.h"
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

template <typename T>
std::ostream &operator<<(std::ostream &out, const std::initializer_list<T> &v) {
  auto column = out.tellp();
  out << "initializer_list{";
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
  typedef std::size_t size_type;

  template <typename T>
  Grid(T &&all_particles) {
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
      number_of_cells_[i] = std::max(
          1, static_cast<int>(std::ceil((max_position[i] - min_position_[i]) *
                                        index_factor_[i])));
      if (number_of_cells_[i] > max_cells) {
        number_of_cells_[i] = max_cells;
        index_factor_[i] = (max_cells - 0.1f)  // -0.1 for safety margin
                           / (max_position[i] - min_position_[i]);
      }
    }

    if (all_of(number_of_cells_, [](size_type n) { return n <= 2; })) {
      // dilute limit:
      // the grid would have <= 2x2x2 cells, meaning every particle has to be
      // compared with every other particle anyway. Then we can just as well
      // fall back to not using the grid at all
      log.debug("There would only be ", number_of_cells_,
                " cells. Therefore the Grid falls back to a single cell / "
                "particle list.");
      number_of_cells_ = {1, 1, 1};
      cells_.reserve(1);
      cells_.emplace_back(std::forward<T>(all_particles));
    } else {
      // construct a normal grid
      log.debug("min: ", min_position_, "\nmax: ", max_position, "\ncells: ",
          number_of_cells_, "\ninteraction length: ",
          max_interaction_length, "\nindex_factor: ", index_factor_);

      // After the grid parameters are determined, we can start placing the
      // particles in cells.
      cells_.resize(number_of_cells_[0] * number_of_cells_[1] *
          number_of_cells_[2]);

      for (const auto &p : all_particles) {
        const auto idx = make_index(p.position().threevec());
        if (idx >= cells_.size()) {
          log.fatal(source_location,
              "\nan out-of-bounds access would be necessary for the particle ", p,
              "\nfor a grid with the following parameters:\nmin: ", min_position_,
              "\nmax: ", max_position, "\ncells: ", number_of_cells_,
              "\ninteraction length: ", max_interaction_length,
              "\nindex_factor: ", index_factor_, "\ncells_.size: ", cells_.size(),
              "\nrequested index: ", idx);
          throw std::runtime_error("out-of-bounds grid access on construction");
        }
        cells_[idx].push_back(p);
      }
    }

    //log.debug(cells_);
  }

  template <typename F>
  void iterate_cells(F &&call_finder) {
    const auto &log = logger<LogArea::Grid>();
    std::vector<const ParticleList *> neighbors;
    neighbors.reserve(13);

    auto &&call_closure = [&](size_type cell_index,
                              const std::initializer_list<int> &xoffsets,
                              const std::initializer_list<int> &yoffsets,
                              const std::initializer_list<int> &zoffsets) {
      neighbors.clear();
      log.debug("call_closure(", cell_index, ", ", xoffsets, ", ", yoffsets,
                ", ", zoffsets, ")");
      for (auto dz : zoffsets) {
        const auto cell_index_dz =
            cell_index + dz * number_of_cells_[1] * number_of_cells_[0];
        for (auto dy : yoffsets) {
          const auto cell_index_dzdy = cell_index_dz + dy * number_of_cells_[0];
          for (auto dx : xoffsets) {
            const auto cell_index_dzdydx = cell_index_dzdy + dx;
            if (cell_index_dzdydx > cell_index) {
              neighbors.emplace_back(&cells_[cell_index_dzdydx]);
            }
          }
        }
      }
      log.debug("iterate_cells calls closure with search_list: ", cells_[cell_index],
                " and neighbors_list: ", neighbors);
      call_finder(cells_[cell_index], neighbors);
    };

    auto &&build_neighbors_with_zy = [&](
        size_type y, size_type z, const std::initializer_list<int> &yoffsets,
        const std::initializer_list<int> &zoffsets) {
      if (number_of_cells_[0] > 1) {
        call_closure(make_index(0, y, z), {0, 1}, yoffsets, zoffsets);
        for (size_type x = 1; x < number_of_cells_[0] - 1; ++x) {
          call_closure(make_index(x, y, z), {-1, 0, 1}, yoffsets, zoffsets);
        }
        call_closure(make_index(number_of_cells_[0] - 1, y, z), {-1, 0},
                     yoffsets, zoffsets);
      } else {
        call_closure(make_index(0, y, z), {0}, yoffsets, zoffsets);
      }
    };

    auto &&build_neighbors_with_z = [&](
        size_type z, const std::initializer_list<int> &zoffsets) {
      if (number_of_cells_[1] > 1) {
        build_neighbors_with_zy(0, z, {0, 1}, zoffsets);
        for (size_type y = 1; y < number_of_cells_[1] - 1; ++y) {
          build_neighbors_with_zy(y, z, {-1, 0, 1}, zoffsets);
        }
        build_neighbors_with_zy(number_of_cells_[1] - 1, z, {-1, 0}, zoffsets);
      } else {
        build_neighbors_with_zy(0, z, {0}, zoffsets);
      }
    };

    if (Options == GridOptions::PeriodicBoundaries) {
    } else {
      for (size_type z = 0; z < number_of_cells_[2] - 1; ++z) {
        build_neighbors_with_z(z, {0, 1});
      }
      build_neighbors_with_z(number_of_cells_[2] - 1, {0});
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
