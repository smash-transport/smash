/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/grid.h"

#include <stdexcept>

#include "include/algorithms.h"
#include "include/fourvector.h"
#include "include/logging.h"
#include "include/particledata.h"
#include "include/threevector.h"

namespace std {
template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
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
std::ostream &operator<<(std::ostream &out, const std::initializer_list<T> &v) {
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
std::ostream &operator<<(std::ostream &out, const std::array<T, N> &a) {
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

namespace Smash {

////////////////////////////////////////////////////////////////////////////////
// GridBase

std::pair<std::array<float, 3>, std::array<float, 3>>
GridBase::find_min_and_length(const Particles &particles) {
  std::pair<std::array<float, 3>, std::array<float, 3>> r;
  auto &min_position = r.first;
  auto &length = r.second;

  // intialize min and max position arrays with the position of the first
  // particle in the list
  const auto &first_position = particles.front().position();
  min_position = {{static_cast<float>(first_position[1]),
                   static_cast<float>(first_position[2]),
                   static_cast<float>(first_position[3])}};
  auto max_position = min_position;
  for (const auto &p : particles) {
    const auto &pos = p.position();
    min_position[0] = std::min(min_position[0], static_cast<float>(pos[1]));
    min_position[1] = std::min(min_position[1], static_cast<float>(pos[2]));
    min_position[2] = std::min(min_position[2], static_cast<float>(pos[3]));
    max_position[0] = std::max(max_position[0], static_cast<float>(pos[1]));
    max_position[1] = std::max(max_position[1], static_cast<float>(pos[2]));
    max_position[2] = std::max(max_position[2], static_cast<float>(pos[3]));
  }
  length[0] = max_position[0] - min_position[0];
  length[1] = max_position[1] - min_position[1];
  length[2] = max_position[2] - min_position[2];
  return r;
}

template <GridOptions O>
std::tuple<std::array<float, 3>, std::array<int, 3>>
Grid<O>::determine_cell_sizes(size_type particle_count,
                              const std::array<float, 3> &length,
                              const int testparticles) {
  std::tuple<std::array<float, 3>, std::array<int, 3>> r;
  auto &index_factor = std::get<0>(r);
  auto &number_of_cells = std::get<1>(r);

  float max_interaction_length = min_cell_length(testparticles);

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
  const int max_cells = O == GridOptions::Normal
                            ? std::cbrt(particle_count)
                            : std::max(2, int(std::cbrt(particle_count)));
  for (std::size_t i = 0; i < number_of_cells.size(); ++i) {
    index_factor[i] = 1.f / max_interaction_length;
    number_of_cells[i] =
        static_cast<int>(std::floor(length[i] * index_factor[i])) +
        // The last cell in each direction can be smaller than
        // max_interaction_length. In that case periodic boundaries will not
        // work correctly. Thus, we need to reduce the number of cells in that
        // direction by one and make the last cell larger. This basically merges
        // a smaller boundary cell into a full cell inside the grid.
        // There's a ~0% chance that the given boundaries create an integral
        // number of cells with length of max_interaction_length. Therefore,
        // just make the default number of cells one less than for non-periodic
        // boundaries.
        (O == GridOptions::Normal ? 1 : 0);

    // std::nextafter implements a safety margin so that no valid position
    // inside the grid can reference an out-of-bounds cell
    if (number_of_cells[i] > max_cells) {
      number_of_cells[i] = max_cells;
      index_factor[i] = number_of_cells[i] / length[i];
      while (index_factor[i] * length[i] >= number_of_cells[i]) {
        index_factor[i] = std::nextafter(index_factor[i], 0.f);
      }
      assert(index_factor[i] * length[i] < number_of_cells[i]);
    } else if (O == GridOptions::PeriodicBoundaries) {
      if (number_of_cells[i] == 1) {
        number_of_cells[i] = 2;
      }
      index_factor[i] = number_of_cells[i] / length[i];
      while (index_factor[i] * length[i] >= number_of_cells[i]) {
        index_factor[i] = std::nextafter(index_factor[i], 0.f);
      }
      assert(index_factor[i] * length[i] < number_of_cells[i]);
    }
  }
  return r;
}

////////////////////////////////////////////////////////////////////////////////
// Grid<Options>

template <GridOptions Options>
inline typename Grid<Options>::size_type Grid<Options>::make_index(
    size_type x, size_type y, size_type z) const {
  return (z * number_of_cells_[1] + y) * number_of_cells_[0] + x;
}

template <GridOptions Options>
inline typename Grid<Options>::size_type Grid<Options>::make_index(
    const FourVector &position) const {
  return make_index(
      std::floor((static_cast<float>(position[1]) - min_position_[0]) *
                 index_factor_[0]),
      std::floor((static_cast<float>(position[2]) - min_position_[1]) *
                 index_factor_[1]),
      std::floor((static_cast<float>(position[3]) - min_position_[2]) *
                 index_factor_[2]));
}

template <GridOptions O>
void Grid<O>::build_cells(const Particles &particles) {
  const auto &log = logger<LogArea::Grid>();
  if (O == GridOptions::Normal &&
      all_of(number_of_cells_, [](size_type n) { return n <= 2; })) {
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
    index_factor_ = {0, 0, 0};
    cells_.resize(1);
    for (const auto &p : particles) {
      if (p.cross_section_scaling_factor() > 0.0) {
        const auto idx = make_index(p.position());
#ifndef NDEBUG
        if (idx >= size_type(cells_.size())) {
          log.fatal(source_location,
                    "\nan out-of-bounds access would be necessary for the "
                    "particle ",
                    p, "\nfor a grid with the following parameters:\nmin: ",
                    min_position_, "\nlength: ", length_, "\ncells: ",
                    number_of_cells_, "\nindex_factor: ", index_factor_,
                    "\ncells_.size: ", cells_.size(), "\nrequested index: ",
                    idx);
          throw std::runtime_error("out-of-bounds grid access on construction");
        }
#endif
        cells_[idx].push_back(p);
      } else {
        continue;
      }
    }
  } else {
    // construct a normal grid
    log.debug("min: ", min_position_, "\nlength: ", length_, "\ncells: ",
              number_of_cells_, "\nindex_factor: ", index_factor_);

    // After the grid parameters are determined, we can start placing the
    // particles in cells.
    cells_.resize(number_of_cells_[0] * number_of_cells_[1] *
                  number_of_cells_[2]);

    for (const auto &p : particles) {
      if (p.cross_section_scaling_factor() > 0.0) {
        const auto idx = make_index(p.position());
#ifndef NDEBUG
        if (idx >= size_type(cells_.size())) {
          log.fatal(source_location,
                    "\nan out-of-bounds access would be necessary for the "
                    "particle ",
                    p, "\nfor a grid with the following parameters:\nmin: ",
                    min_position_, "\nlength: ", length_, "\ncells: ",
                    number_of_cells_, "\nindex_factor: ", index_factor_,
                    "\ncells_.size: ", cells_.size(), "\nrequested index: ",
                    idx);
          throw std::runtime_error("out-of-bounds grid access on construction");
        }
#endif
        cells_[idx].push_back(p);
      } else {
        continue;
      }
    }
  }

  log.debug(cells_);
}

template <>
void Grid<GridOptions::Normal>::iterate_cells(
    const std::function<void(const ParticleList &)> &search_cell_callback,
    const std::function<void(const ParticleList &, const ParticleList &)> &
        neighbor_cell_callback) const {
  std::array<size_type, 3> search_index;
  size_type &x = search_index[0];
  size_type &y = search_index[1];
  size_type &z = search_index[2];
  size_type search_cell_index = 0;
  for (z = 0; z < number_of_cells_[2]; ++z) {
    for (y = 0; y < number_of_cells_[1]; ++y) {
      for (x = 0; x < number_of_cells_[0]; ++x, ++search_cell_index) {
        assert(search_cell_index == make_index(search_index));
        assert(search_cell_index >= 0);
        assert(search_cell_index < size_type(cells_.size()));
        const ParticleList &search = cells_[search_cell_index];
        search_cell_callback(search);

        const auto dz_list = z == number_of_cells_[2] - 1
                                 ? std::initializer_list<size_type>{0}
                                 : std::initializer_list<size_type>{0, 1};
        const auto dy_list =
            number_of_cells_[1] == 1
                ? std::initializer_list<size_type>{0}
                : y == 0 ? std::initializer_list<size_type>{0, 1}
                         : y == number_of_cells_[1] - 1
                               ? std::initializer_list<size_type>{-1, 0}
                               : std::initializer_list<size_type>{-1, 0, 1};
        const auto dx_list =
            number_of_cells_[0] == 1
                ? std::initializer_list<size_type>{0}
                : x == 0 ? std::initializer_list<size_type>{0, 1}
                         : x == number_of_cells_[0] - 1
                               ? std::initializer_list<size_type>{-1, 0}
                               : std::initializer_list<size_type>{-1, 0, 1};
        for (size_type dz : dz_list) {
          for (size_type dy : dy_list) {
            for (size_type dx : dx_list) {
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

enum class NeedsToWrap { PlusLength, No, MinusLength };
struct NeighborLookup {
  typename Grid<GridOptions::PeriodicBoundaries>::size_type index = 0;
  NeedsToWrap wrap = NeedsToWrap::No;
};

template <>
void Grid<GridOptions::PeriodicBoundaries>::iterate_cells(
    const std::function<void(const ParticleList &)> &search_cell_callback,
    const std::function<void(const ParticleList &, const ParticleList &)> &
        neighbor_cell_callback) const {
  const auto &log = logger<LogArea::Grid>();

  std::array<size_type, 3> search_index;
  size_type &x = search_index[0];
  size_type &y = search_index[1];
  size_type &z = search_index[2];
  size_type search_cell_index = 0;

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
        assert(search_cell_index < size_type(cells_.size()));
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
              assert(neighbor_cell_index < size_type(cells_.size()));
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

template std::tuple<std::array<float, 3>, std::array<int, 3>>
Grid<GridOptions::Normal>::determine_cell_sizes(size_type,
                                                const std::array<float, 3> &,
                                                const int);
template std::tuple<std::array<float, 3>, std::array<int, 3>>
Grid<GridOptions::PeriodicBoundaries>::determine_cell_sizes(
    size_type, const std::array<float, 3> &, const int);

template void Grid<GridOptions::Normal>::build_cells(const Particles &);
template void Grid<GridOptions::PeriodicBoundaries>::build_cells(
    const Particles &);
}  // namespace Smash
