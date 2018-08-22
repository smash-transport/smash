/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/grid.h"

#include <stdexcept>

#include "smash/algorithms.h"
#include "smash/fourvector.h"
#include "smash/logging.h"
#include "smash/particledata.h"
#include "smash/threevector.h"

namespace std {
template <typename T>
static std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  auto column = out.tellp();
  out << "{ ";
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
static std::ostream &operator<<(std::ostream &out,
                                const std::initializer_list<T> &v) {
  auto column = out.tellp();
  out << "{ ";
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
static std::ostream &operator<<(std::ostream &out, const std::array<T, N> &a) {
  auto column = out.tellp();
  out << "{ ";
  for (const auto &x : a) {
    if (out.tellp() - column >= 100) {
      out << '\n';
      column = out.tellp();
    }
    out << x << ' ';
  }
  return out << '}';
}
}  // namespace std

namespace smash {

////////////////////////////////////////////////////////////////////////////////
// GridBase

std::pair<std::array<double, 3>, std::array<double, 3>>
GridBase::find_min_and_length(const Particles &particles) {
  std::pair<std::array<double, 3>, std::array<double, 3>> r;
  auto &min_position = r.first;
  auto &length = r.second;

  // intialize min and max position arrays with the position of the first
  // particle in the list
  const auto &first_position = particles.front().position();
  min_position = {{first_position[1], first_position[2], first_position[3]}};
  auto max_position = min_position;
  for (const auto &p : particles) {
    const auto &pos = p.position();
    min_position[0] = std::min(min_position[0], pos[1]);
    min_position[1] = std::min(min_position[1], pos[2]);
    min_position[2] = std::min(min_position[2], pos[3]);
    max_position[0] = std::max(max_position[0], pos[1]);
    max_position[1] = std::max(max_position[1], pos[2]);
    max_position[2] = std::max(max_position[2], pos[3]);
  }
  length[0] = max_position[0] - min_position[0];
  length[1] = max_position[1] - min_position[1];
  length[2] = max_position[2] - min_position[2];
  return r;
}

////////////////////////////////////////////////////////////////////////////////
// Grid

template <GridOptions O>
Grid<O>::Grid(const std::pair<std::array<double, 3>, std::array<double, 3>>
                  &min_and_length,
              const Particles &particles, double max_interaction_length,
              CellSizeStrategy strategy)
    : length_(min_and_length.second) {
  const auto min_position = min_and_length.first;
  const SizeType particle_count = particles.size();

  // very simple setup for non-periodic boundaries and largest cellsize strategy
  if (O == GridOptions::Normal && strategy == CellSizeStrategy::Largest) {
    number_of_cells_ = {1, 1, 1};
    cells_.clear();
    cells_.reserve(1);
    cells_.emplace_back(particles.copy_to_vector());
    return;
  }

  // The number of cells is determined by the min and max coordinates where
  // particles are positioned and the maximal interaction length (which equals
  // the length of a cell).
  // But don't let the number of cells exceed the actual number of particles.
  // That would be overkill. Let max_cells³ ≤ particle_count (conversion to
  // int truncates).
  // Consider that particle placement into cells uses half-open intervals. Thus
  // a cell includes particles in [0, a[. The next cell [a, 2a[. And so on. This
  // is important for calculating the number of cells. If length * index_factor
  // (equivalent to length / max_interaction_length) is integral, then
  // length * index_factor + 1 determines the number of required cells. That's
  // because the last cell will then store particles in the interval
  // [length, length + max_interaction_length[. The code below achieves this
  // effect by rounding down (floor) and adding 1 afterwards.
  const int max_cells =
      (O == GridOptions::Normal)
          ? std::cbrt(particle_count)
          : std::max(2, static_cast<int>(std::cbrt(particle_count)));

  // This normally equals 1/max_interaction_length, but if the number of cells
  // is reduced (because of low density) then this value is smaller.
  std::array<double, 3> index_factor = {1. / max_interaction_length,
                                        1. / max_interaction_length,
                                        1. / max_interaction_length};
  for (std::size_t i = 0; i < number_of_cells_.size(); ++i) {
    number_of_cells_[i] =
        (strategy == CellSizeStrategy::Largest)
            ? 2
            : static_cast<int>(std::floor(length_[i] * index_factor[i])) +
                  // The last cell in each direction can be smaller than
                  // max_interaction_length. In that case periodic boundaries
                  // will not work correctly. Thus, we need to reduce the number
                  // of cells in that direction by one and make the last cell
                  // larger. This basically merges a smaller boundary cell into
                  // a full cell inside the grid. There's a ~0% chance that the
                  // given boundaries create an integral number of cells with
                  // length of max_interaction_length. Therefore, just make the
                  // default number of cells one less than for non-periodic
                  // boundaries.
                  (O == GridOptions::Normal ? 1 : 0);

    // std::nextafter implements a safety margin so that no valid position
    // inside the grid can reference an out-of-bounds cell
    if (number_of_cells_[i] > max_cells) {
      number_of_cells_[i] = max_cells;
      index_factor[i] = number_of_cells_[i] / length_[i];
      while (index_factor[i] * length_[i] >= number_of_cells_[i]) {
        index_factor[i] = std::nextafter(index_factor[i], 0.);
      }
      assert(index_factor[i] * length_[i] < number_of_cells_[i]);
    } else if (O == GridOptions::PeriodicBoundaries) {
      if (number_of_cells_[i] == 1) {
        number_of_cells_[i] = 2;
      }
      index_factor[i] = number_of_cells_[i] / length_[i];
      while (index_factor[i] * length_[i] >= number_of_cells_[i]) {
        index_factor[i] = std::nextafter(index_factor[i], 0.);
      }
      assert(index_factor[i] * length_[i] < number_of_cells_[i]);
    }
  }

  const auto &log = logger<LogArea::Grid>();
  if (O == GridOptions::Normal &&
      all_of(number_of_cells_, [](SizeType n) { return n <= 2; })) {
    // dilute limit:
    // the grid would have <= 2x2x2 cells, meaning every particle has to be
    // compared with every other particle anyway. Then we can just as well
    // fall back to not using the grid at all
    // For a grid with periodic boundaries the situation is different and we
    // never want to have a grid smaller than 2x2x2.
    log.debug("There would only be ", number_of_cells_,
              " cells. Therefore the Grid falls back to a single cell / "
              "particle list.");
    number_of_cells_ = {1, 1, 1};
    cells_.resize(1);
    cells_.front().reserve(particles.size());
    std::copy_if(particles.begin(), particles.end(),
                 std::back_inserter(cells_.front()), [](const ParticleData &p) {
                   return p.xsec_scaling_factor() > 0.0;
                 });  // filter out the particles that can not interact
  } else {
    // construct a normal grid
    log.debug("min: ", min_position, "\nlength: ", length_,
              "\ncells: ", number_of_cells_, "\nindex_factor: ", index_factor);

    // After the grid parameters are determined, we can start placing the
    // particles in cells.
    cells_.resize(number_of_cells_[0] * number_of_cells_[1] *
                  number_of_cells_[2]);

    // Returns the one-dimensional cell-index from the position vector inside
    // the grid.
    // This simply calculates the distance to min_position and multiplies it
    // with index_factor to determine the 3 x,y,z indexes to pass to make_index.
    auto &&cell_index_for = [&](const ParticleData &p) {
      return make_index(
          std::floor((p.position()[1] - min_position[0]) * index_factor[0]),
          std::floor((p.position()[2] - min_position[1]) * index_factor[1]),
          std::floor((p.position()[3] - min_position[2]) * index_factor[2]));
    };

    for (const auto &p : particles) {
      if (p.xsec_scaling_factor() > 0.0) {
        const auto idx = cell_index_for(p);
#ifndef NDEBUG
        if (idx >= SizeType(cells_.size())) {
          log.fatal(
              source_location,
              "\nan out-of-bounds access would be necessary for the "
              "particle ",
              p, "\nfor a grid with the following parameters:\nmin: ",
              min_position, "\nlength: ", length_,
              "\ncells: ", number_of_cells_, "\nindex_factor: ", index_factor,
              "\ncells_.size: ", cells_.size(), "\nrequested index: ", idx);
          throw std::runtime_error("out-of-bounds grid access on construction");
        }
#endif
        cells_[idx].push_back(p);
      }
    }
  }

  log.debug(cells_);
}

template <GridOptions Options>
inline typename Grid<Options>::SizeType Grid<Options>::make_index(
    SizeType x, SizeType y, SizeType z) const {
  return (z * number_of_cells_[1] + y) * number_of_cells_[0] + x;
}

template <>
/// Specialization of iterate_cells
void Grid<GridOptions::Normal>::iterate_cells(
    const std::function<void(const ParticleList &)> &search_cell_callback,
    const std::function<void(const ParticleList &, const ParticleList &)>
        &neighbor_cell_callback) const {
  std::array<SizeType, 3> search_index;
  SizeType &x = search_index[0];
  SizeType &y = search_index[1];
  SizeType &z = search_index[2];
  SizeType search_cell_index = 0;
  for (z = 0; z < number_of_cells_[2]; ++z) {
    for (y = 0; y < number_of_cells_[1]; ++y) {
      for (x = 0; x < number_of_cells_[0]; ++x, ++search_cell_index) {
        assert(search_cell_index == make_index(search_index));
        assert(search_cell_index >= 0);
        assert(search_cell_index < SizeType(cells_.size()));
        const ParticleList &search = cells_[search_cell_index];
        search_cell_callback(search);

        const auto dz_list = z == number_of_cells_[2] - 1
                                 ? std::initializer_list<SizeType>{0}
                                 : std::initializer_list<SizeType>{0, 1};
        const auto dy_list =
            number_of_cells_[1] == 1
                ? std::initializer_list<SizeType>{0}
                : y == 0 ? std::initializer_list<SizeType>{0, 1}
                         : y == number_of_cells_[1] - 1
                               ? std::initializer_list<SizeType>{-1, 0}
                               : std::initializer_list<SizeType>{-1, 0, 1};
        const auto dx_list =
            number_of_cells_[0] == 1
                ? std::initializer_list<SizeType>{0}
                : x == 0 ? std::initializer_list<SizeType>{0, 1}
                         : x == number_of_cells_[0] - 1
                               ? std::initializer_list<SizeType>{-1, 0}
                               : std::initializer_list<SizeType>{-1, 0, 1};
        for (SizeType dz : dz_list) {
          for (SizeType dy : dy_list) {
            for (SizeType dx : dx_list) {
              const auto di = make_index(dx, dy, dz);
              if (di > 0) {
                neighbor_cell_callback(search, cells_[search_cell_index + di]);
              }
            }
          }
        }
      }
    }
  }
}

/**
 * The options determining what to do if a particle flies out of the grids
 * PlusLength:  Used if a periodic boundary condition is applied and a
 *              particle passes through the lower bound of the grid.
 * No:          Used if the boundary condition is not periodic.
 * MinusLength: Used if a periodic boundary condition is applied and a
 *              particle passes through the upper bound of the grid.
 */
enum class NeedsToWrap { PlusLength, No, MinusLength };

/// A strust containing the informations needed to search the neighboring cell
struct NeighborLookup {
  /// Index of the cell
  typename Grid<GridOptions::PeriodicBoundaries>::SizeType index = 0;
  /// Option to determine the neighbors of the cells on the boundary
  NeedsToWrap wrap = NeedsToWrap::No;
};

template <>
/// Specialization of iterate_cells
void Grid<GridOptions::PeriodicBoundaries>::iterate_cells(
    const std::function<void(const ParticleList &)> &search_cell_callback,
    const std::function<void(const ParticleList &, const ParticleList &)>
        &neighbor_cell_callback) const {
  const auto &log = logger<LogArea::Grid>();

  std::array<SizeType, 3> search_index;
  SizeType &x = search_index[0];
  SizeType &y = search_index[1];
  SizeType &z = search_index[2];
  SizeType search_cell_index = 0;

  // defaults:
  std::array<NeighborLookup, 2> dz_list;
  std::array<NeighborLookup, 3> dy_list;
  std::array<NeighborLookup, 3> dx_list;

  assert(number_of_cells_[2] >= 2);
  assert(number_of_cells_[1] >= 2);
  assert(number_of_cells_[0] >= 2);

  for (z = 0; z < number_of_cells_[2]; ++z) {
    dz_list[0].index = z;
    dz_list[1].index = z + 1;
    if (dz_list[1].index == number_of_cells_[2]) {
      dz_list[1].index = 0;
      // last z in the loop, so no need to reset wrap again
      dz_list[1].wrap = NeedsToWrap::MinusLength;
    }
    for (y = 0; y < number_of_cells_[1]; ++y) {
      dy_list[0].index = y;
      dy_list[1].index = y - 1;
      dy_list[2].index = y + 1;
      dy_list[2].wrap = NeedsToWrap::No;
      if (y == 0) {
        dy_list[1] = dy_list[2];
        dy_list[2].index = number_of_cells_[1] - 1;
        dy_list[2].wrap = NeedsToWrap::PlusLength;
      } else if (dy_list[2].index == number_of_cells_[1]) {
        dy_list[2].index = 0;
        dy_list[2].wrap = NeedsToWrap::MinusLength;
      }
      for (x = 0; x < number_of_cells_[0]; ++x, ++search_cell_index) {
        dx_list[0].index = x;
        dx_list[1].index = x - 1;
        dx_list[2].index = x + 1;
        dx_list[2].wrap = NeedsToWrap::No;
        if (x == 0) {
          dx_list[1] = dx_list[2];
          dx_list[2].index = number_of_cells_[0] - 1;
          dx_list[2].wrap = NeedsToWrap::PlusLength;
        } else if (dx_list[2].index == number_of_cells_[0]) {
          dx_list[2].index = 0;
          dx_list[2].wrap = NeedsToWrap::MinusLength;
        }

        assert(search_cell_index == make_index(search_index));
        assert(search_cell_index >= 0);
        assert(search_cell_index < SizeType(cells_.size()));
        ParticleList search = cells_[search_cell_index];
        search_cell_callback(search);

        auto virtual_search_index = search_index;
        ThreeVector wrap_vector = {};  // no change
        auto current_wrap_vector = wrap_vector;

        for (const auto &dz : dz_list) {
          if (dz.wrap == NeedsToWrap::MinusLength) {
            // last dz in the loop, so no need to undo the wrap
            wrap_vector[2] = -length_[2];
            virtual_search_index[2] = -1;
          }
          for (const auto &dy : dy_list) {
            // only the last dy in dy_list can wrap
            if (dy.wrap == NeedsToWrap::MinusLength) {
              wrap_vector[1] = -length_[1];
              virtual_search_index[1] = -1;
            } else if (dy.wrap == NeedsToWrap::PlusLength) {
              wrap_vector[1] = length_[1];
              virtual_search_index[1] = number_of_cells_[1];
            }
            for (const auto &dx : dx_list) {
              // only the last dx in dx_list can wrap
              if (dx.wrap == NeedsToWrap::MinusLength) {
                wrap_vector[0] = -length_[0];
                virtual_search_index[0] = -1;
              } else if (dx.wrap == NeedsToWrap::PlusLength) {
                wrap_vector[0] = length_[0];
                virtual_search_index[0] = number_of_cells_[0];
              }
              assert(dx.index >= 0);
              assert(dx.index < number_of_cells_[0]);
              assert(dy.index >= 0);
              assert(dy.index < number_of_cells_[1]);
              assert(dz.index >= 0);
              assert(dz.index < number_of_cells_[2]);
              const auto neighbor_cell_index =
                  make_index(dx.index, dy.index, dz.index);
              assert(neighbor_cell_index >= 0);
              assert(neighbor_cell_index < SizeType(cells_.size()));
              if (neighbor_cell_index <= make_index(virtual_search_index)) {
                continue;
              }

              if (wrap_vector != current_wrap_vector) {
                log.debug("translating search cell by ",
                          wrap_vector - current_wrap_vector);
                for_each(search, [&](ParticleData &p) {
                  p = p.translated(wrap_vector - current_wrap_vector);
                });
                current_wrap_vector = wrap_vector;
              }
              neighbor_cell_callback(search, cells_[neighbor_cell_index]);
            }
            virtual_search_index[0] = search_index[0];
            wrap_vector[0] = 0;
          }
          virtual_search_index[1] = search_index[1];
          wrap_vector[1] = 0;
        }
      }
    }
  }
}

template Grid<GridOptions::Normal>::Grid(
    const std::pair<std::array<double, 3>, std::array<double, 3>>
        &min_and_length,
    const Particles &particles, double max_interaction_length,
    CellSizeStrategy strategy);
template Grid<GridOptions::PeriodicBoundaries>::Grid(
    const std::pair<std::array<double, 3>, std::array<double, 3>>
        &min_and_length,
    const Particles &particles, double max_interaction_length,
    CellSizeStrategy strategy);
}  // namespace smash
