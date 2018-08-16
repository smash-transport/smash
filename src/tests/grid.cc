/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include "../include/smash/grid.h"
#include "../include/smash/logging.h"

#include <set>
#include <unordered_set>

using namespace smash;

namespace std {
// helpers for printing sets/pairs on test failure
template <typename T, typename U>
ostream &operator<<(ostream &s, const pair<T, U> &p) {
  return s << '<' << p.first << '|' << p.second << '>';
}
template <typename T>
ostream &operator<<(ostream &s, const set<T> &set) {
  s << '{';
  for (auto &&x : set) {
    s << x << ' ';
  }
  return s << '}';
}
template <typename T>
static inline ostream &operator<<(ostream &s, const unordered_set<T> &data) {
  s << '{';
  for (auto &&x : data) {
    s << x << ' ';
  }
  return s << '}';
}
}  // namespace std

static inline double minimal_cell_length(int testparticles) {
  // In SMASH itself the minimal cell length is calculated by the function
  // ScatterActionsFinder::min_cell_length(). It uses the maximum cross section
  // of 200mb and the time step size to calculate the cell length. This test was
  // written before the maximum cross section was introduced and does not
  // accommodate for time step sizes. It assumes a minimal cell length based on
  // a maximal interaction length of 2.5fm.
  //
  // To make the test work with ScatterActionsFinder::min_cell_length() the
  // placements of the particles would need to be adapted to the different
  // minimal cell size.
  return 2.5 / std::sqrt(static_cast<double>(testparticles));
}

TEST(init) {
  set_default_loglevel(einhard::WARN);
  // create_all_loggers("Grid: DEBUG");
  Test::create_smashon_particletypes();
}

// Grid tests:
// - Test that the number of cells is determined correctly. I.e. what is
// expected for a specific volume of particles.
// - Test that the contents of each cell are correct.
// - Don't use random values! For a given set of positions an expected grid
// configuration can be written out. For random positions the expected grid must
// be determined programatically which can lead to errors in the the test code
// itself masking actual errors in the grid code.

TEST(grid_construction) {
  using NeighborsSet = std::set<std::pair<int, int>>;
  struct Parameter {
    // input:
    ParticleList particles;

    // expected results:
    std::size_t cellcount[3];  // per direction
    NeighborsSet neighbors;
    std::vector<std::unordered_set<int>> ids;
  };
  auto &&make_particle = [](double x, double y, double z) {
    return Test::smashon(Test::Position{0., x, y, z});
  };
  for (const int testparticles : {1, 5, 20, 100}) {
    const double min_cell_length = minimal_cell_length(testparticles);
    for (const Parameter &param : std::vector<Parameter>{
             Parameter{
                 {make_particle(0., 0., 0.), make_particle(1.9, 1.9, 1.9)},
                 {1, 1, 1},
                 NeighborsSet{},
                 {{0, 1}}},
             {{make_particle(0, 0, 0), make_particle(0, 0, 1),
               make_particle(0, 0, 2), make_particle(0, 1, 0),
               make_particle(0, 1, 1), make_particle(0, 1, 2),
               make_particle(0, 2, 0), make_particle(0, 2, 1),
               make_particle(0, 2, 2), make_particle(1, 0, 0),
               make_particle(1, 0, 1), make_particle(1, 0, 2),
               make_particle(1, 1, 0), make_particle(1, 1, 1),
               make_particle(1, 1, 2), make_particle(1, 2, 0),
               make_particle(1, 2, 1), make_particle(1, 2, 2),
               make_particle(2, 0, 0), make_particle(2, 0, 1),
               make_particle(2, 0, 2), make_particle(2, 1, 0),
               make_particle(2, 1, 1), make_particle(2, 1, 2),
               make_particle(2, 2, 0), make_particle(2, 2, 1),
               make_particle(2, 2, 2)},
              {3, 3, 3},
              {{0, 1},   {0, 3},   {0, 4},   {0, 9},   {0, 10},  {0, 12},
               {0, 13},  {1, 2},   {1, 3},   {1, 4},   {1, 5},   {1, 9},
               {1, 10},  {1, 11},  {1, 12},  {1, 13},  {1, 14},  {2, 4},
               {2, 5},   {2, 10},  {2, 11},  {2, 13},  {2, 14},  {3, 4},
               {3, 6},   {3, 7},   {3, 9},   {3, 10},  {3, 12},  {3, 13},
               {3, 15},  {3, 16},  {4, 5},   {4, 6},   {4, 7},   {4, 8},
               {4, 9},   {4, 10},  {4, 11},  {4, 12},  {4, 13},  {4, 14},
               {4, 15},  {4, 16},  {4, 17},  {5, 7},   {5, 8},   {5, 10},
               {5, 11},  {5, 13},  {5, 14},  {5, 16},  {5, 17},  {6, 7},
               {6, 12},  {6, 13},  {6, 15},  {6, 16},  {7, 8},   {7, 12},
               {7, 13},  {7, 14},  {7, 15},  {7, 16},  {7, 17},  {8, 13},
               {8, 14},  {8, 16},  {8, 17},  {9, 10},  {9, 12},  {9, 13},
               {9, 18},  {9, 19},  {9, 21},  {9, 22},  {10, 11}, {10, 12},
               {10, 13}, {10, 14}, {10, 18}, {10, 19}, {10, 20}, {10, 21},
               {10, 22}, {10, 23}, {11, 13}, {11, 14}, {11, 19}, {11, 20},
               {11, 22}, {11, 23}, {12, 13}, {12, 15}, {12, 16}, {12, 18},
               {12, 19}, {12, 21}, {12, 22}, {12, 24}, {12, 25}, {13, 14},
               {13, 15}, {13, 16}, {13, 17}, {13, 18}, {13, 19}, {13, 20},
               {13, 21}, {13, 22}, {13, 23}, {13, 24}, {13, 25}, {13, 26},
               {14, 16}, {14, 17}, {14, 19}, {14, 20}, {14, 22}, {14, 23},
               {14, 25}, {14, 26}, {15, 16}, {15, 21}, {15, 22}, {15, 24},
               {15, 25}, {16, 17}, {16, 21}, {16, 22}, {16, 23}, {16, 24},
               {16, 25}, {16, 26}, {17, 22}, {17, 23}, {17, 25}, {17, 26},
               {18, 19}, {18, 21}, {18, 22}, {19, 20}, {19, 21}, {19, 22},
               {19, 23}, {20, 22}, {20, 23}, {21, 22}, {21, 24}, {21, 25},
               {22, 23}, {22, 24}, {22, 25}, {22, 26}, {23, 25}, {23, 26},
               {24, 25}, {25, 26}},
              {{0}, {9},  {18}, {3}, {12}, {21}, {6}, {15}, {24},
               {1}, {10}, {19}, {4}, {13}, {22}, {7}, {16}, {25},
               {2}, {11}, {20}, {5}, {14}, {23}, {8}, {17}, {26}}},
         }) {
      Particles list;
      for (auto p : param.particles) {
        p.set_4position(min_cell_length * p.position());
        list.insert(p);
      }
      Grid<GridOptions::Normal> grid(list, min_cell_length);
      auto idsIt = param.ids.begin();
      auto neighbors = param.neighbors;
      grid.iterate_cells(
          [&](const ParticleList &search) {
            auto ids = *idsIt++;
            for (const auto &p : search) {
              COMPARE(ids.erase(p.id()), 1u)
                  << "p.id() = " << p.id() << ", ids = " << ids;
            }
            COMPARE(ids.size(), 0u);
          },
          [&](const ParticleList &search, const ParticleList &n) {
            for (const auto &p : search) {
              for (const auto &p2 : n) {
                COMPARE(neighbors.erase({std::min(p.id(), p2.id()),
                                         std::max(p.id(), p2.id())}),
                        1u)
                    << "<id|id>: <" << std::min(p.id(), p2.id()) << '|'
                    << std::max(p.id(), p2.id()) << '>';
              }
            }
          });
      COMPARE(neighbors.size(), 0u) << neighbors;
    }
  }
}

template <typename Container, typename T>
typename Container::const_iterator find(const Container &c, const T &value) {
  return std::find(std::begin(c), std::end(c), value);
}

namespace std {
template <typename T, typename U>
std::ostream &operator<<(std::ostream &s,
                         const std::vector<std::pair<T, U>> &l) {
  for (auto &&p : l) {
    s << "\n<" << p.first << ", " << p.second << '>';
  }
  return s;
}
}  // namespace std

TEST(periodic_grid) {
  using Test::Momentum;
  using Test::Position;
  for (const int testparticles : {1, 5}) {
    for (const int nparticles : {1, 5, 20, 75, 124, 125}) {
      const double min_cell_length = minimal_cell_length(testparticles);
      constexpr double length = 10;
      Particles list;
      auto random_value = random::make_uniform_distribution(0., 9.99);
      for (auto n = nparticles; n; --n) {
        list.insert(Test::smashon(
            Position{0., random_value(), random_value(), random_value()},
            Momentum{Test::smashon_mass,
                     {random_value(), random_value(), random_value()}},
            n));
      }
      Grid<GridOptions::PeriodicBoundaries> grid(
          make_pair(std::array<double, 3>{0, 0, 0},
                    std::array<double, 3>{length, length, length}),
          list, min_cell_length);

      // stores the neighbor pairs found via the grid:
      std::vector<std::pair<ParticleData, ParticleData>> neighbor_pairs;

      grid.iterate_cells(
          [&](const ParticleList &search) {
            for (const ParticleData &p : search) {
              {
                const auto it = find(list, p);
                VERIFY(it != list.cend());
                COMPARE(it->id(), p.id());
                COMPARE(it->position(), p.position());
              }

              for (const ParticleData &q : search) {
                if (p.id() <= q.id()) {
                  continue;
                }
                const auto sqrDistance =
                    (p.position().threevec() - q.position().threevec()).sqr();
                if (sqrDistance <= min_cell_length * min_cell_length) {
                  const auto pair = p.id() < q.id() ? std::make_pair(p, q)
                                                    : std::make_pair(q, p);
                  const auto it = find(neighbor_pairs, pair);
                  COMPARE(it, neighbor_pairs.end())
                      << "\np: " << p << "\nq: " << q << '\n'
                      << detailed(search);
                  neighbor_pairs.emplace_back(std::move(pair));
                }
              }
            }
          },
          [&](const ParticleList &search, const ParticleList &neighbors) {
            // for each particle in neighbors, find the same particle in list
            for (const ParticleData &p : neighbors) {
              const auto it = find(list, p);
              VERIFY(it != list.cend());
              COMPARE(it->id(), p.id());
              COMPARE(it->position(), p.position());
            }
            auto &&compareDiff = [length](double d) {
              if (d < -0.1 * length) {
                FUZZY_COMPARE(d, -length);
              } else if (d > 0.1 * length) {
                FUZZY_COMPARE(d, length);
              } else {
                COMPARE_ABSOLUTE_ERROR(
                    d, 0., length * std::numeric_limits<double>::epsilon());
              }
            };
            // for each particle in search, find the same particle in list
            for (const ParticleData &p : search) {
              const auto it = find(list, p);
              VERIFY(it != list.cend());
              COMPARE(it->id(), p.id());
              if (it->position() != p.position()) {
                // then the cell was wrapped around
                const auto diff = it->position() - p.position();
                COMPARE(diff[0], 0.);
                compareDiff(diff[1]);
                compareDiff(diff[2]);
                compareDiff(diff[3]);
                VERIFY(diff != FourVector(0, 0, 0, 0));
              }
            }

            // for each particle in search, search through the complete list of
            // neighbors to find those closer than 2.5fm
            for (const ParticleData &p : search) {
              for (const ParticleData &q : neighbors) {
                VERIFY(!(p == q)) << "\np: " << p << "\nq: " << q << '\n'
                                  << search << '\n'
                                  << neighbors;
                const auto sqrDistance =
                    (p.position().threevec() - q.position().threevec()).sqr();
                if (sqrDistance <= min_cell_length * min_cell_length) {
                  auto pair = p.id() < q.id() ? std::make_pair(p, q)
                                              : std::make_pair(q, p);
                  const auto it = find(neighbor_pairs, pair);
                  COMPARE(it, neighbor_pairs.end())
                      << "\np: " << p << "\nq: " << q << '\n'
                      << neighbor_pairs;
                  neighbor_pairs.emplace_back(std::move(pair));
                }
              }
            }
          });

      // Now search through the original list to verify the grid search found
      // everything.
      auto &&wrap = [](ParticleData p, int i) {
        auto pos = p.position();
        if (i > 0) {
          pos[i] += length;
        } else {
          pos[-i] -= length;
        }
        p.set_4position(pos);
        return p;
      };
      for (ParticleData p0 : list) {
        ParticleList p_periodic;
        p_periodic.push_back(p0);
        if (p0.position()[3] < min_cell_length) {
          p_periodic.push_back(wrap(p0, 3));
        }
        if (p0.position()[3] > length - min_cell_length) {
          p_periodic.push_back(wrap(p0, -3));
        }
        if (p0.position()[2] < min_cell_length) {
          p_periodic.push_back(wrap(p0, 2));
          if (p0.position()[3] < min_cell_length) {
            p_periodic.push_back(wrap(p0, 3));
          }
          if (p0.position()[3] > length - min_cell_length) {
            p_periodic.push_back(wrap(p0, -3));
          }
        }
        if (p0.position()[2] > length - min_cell_length) {
          p_periodic.push_back(wrap(p0, -2));
          if (p0.position()[3] < min_cell_length) {
            p_periodic.push_back(wrap(p0, 3));
          }
          if (p0.position()[3] > length - min_cell_length) {
            p_periodic.push_back(wrap(p0, -3));
          }
        }
        if (p0.position()[1] < min_cell_length) {
          p_periodic.push_back(wrap(p0, 1));
          if (p0.position()[3] < min_cell_length) {
            p_periodic.push_back(wrap(p0, 3));
          }
          if (p0.position()[3] > length - min_cell_length) {
            p_periodic.push_back(wrap(p0, -3));
          }
          if (p0.position()[2] < min_cell_length) {
            p_periodic.push_back(wrap(p0, 2));
            if (p0.position()[3] < min_cell_length) {
              p_periodic.push_back(wrap(p0, 3));
            }
            if (p0.position()[3] > length - min_cell_length) {
              p_periodic.push_back(wrap(p0, -3));
            }
          }
          if (p0.position()[2] > length - min_cell_length) {
            p_periodic.push_back(wrap(p0, -2));
            if (p0.position()[3] < min_cell_length) {
              p_periodic.push_back(wrap(p0, 3));
            }
            if (p0.position()[3] > length - min_cell_length) {
              p_periodic.push_back(wrap(p0, -3));
            }
          }
        }
        if (p0.position()[1] > length - min_cell_length) {
          p_periodic.push_back(wrap(p0, -1));
          if (p0.position()[3] < min_cell_length) {
            p_periodic.push_back(wrap(p0, 3));
          }
          if (p0.position()[3] > length - min_cell_length) {
            p_periodic.push_back(wrap(p0, -3));
          }
          if (p0.position()[2] < min_cell_length) {
            p_periodic.push_back(wrap(p0, 2));
            if (p0.position()[3] < min_cell_length) {
              p_periodic.push_back(wrap(p0, 3));
            }
            if (p0.position()[3] > length - min_cell_length) {
              p_periodic.push_back(wrap(p0, -3));
            }
          }
          if (p0.position()[2] > length - min_cell_length) {
            p_periodic.push_back(wrap(p0, -2));
            if (p0.position()[3] < min_cell_length) {
              p_periodic.push_back(wrap(p0, 3));
            }
            if (p0.position()[3] > length - min_cell_length) {
              p_periodic.push_back(wrap(p0, -3));
            }
          }
        }

        sort(neighbor_pairs.begin(), neighbor_pairs.end());
        for (const ParticleData &p : p_periodic) {
          for (const ParticleData &q : list) {
            if (p == q) {
              continue;
            }
            const auto sqrDistance =
                (p.position().threevec() - q.position().threevec()).sqr();
            if (sqrDistance <= min_cell_length * min_cell_length) {
              // (p,q) must be in neighbor_pairs
              auto pair =
                  p.id() > q.id() ? std::make_pair(q, p) : std::make_pair(p, q);
              const auto it = find(neighbor_pairs, pair);
              VERIFY(it != neighbor_pairs.end())
                  << "\ntestparticles: " << testparticles
                  << "\nnparticles: " << nparticles << "\np: " << p
                  << "\nq: " << q << "\n"
                  << neighbor_pairs;
            }
          }
        }
      }
    }
  }
}

TEST(max_positions_periodic_grid) {
  constexpr int testparticles = 1;
  const double min_cell_length = minimal_cell_length(testparticles);
  using Test::Position;
  Particles list;
  list.insert(Test::smashon(Position{0, 0, 0, 0}));
  list.insert(Test::smashon(Position{0, 0, 0, 0}));
  list.insert(Test::smashon(Position{0, 0, 0, 0}));
  list.insert(Test::smashon(Position{0, 0, 0, 0}));
  list.insert(Test::smashon(Position{0, 0, 0, 0}));
  list.insert(Test::smashon(Position{0, 0, 0, 0}));
  list.insert(Test::smashon(Position{0, 0, 0, 0}));
  list.insert(Test::smashon(Position{0, 0, 0, 0}));
  list.insert(Test::smashon(Position{
      0, 2 * min_cell_length, 2 * min_cell_length, 2 * min_cell_length}));
  // This grid construction is fragile because there are particles at 0 and 2 *
  // min_cell_length. A Normal grid simply would try to create a 3x3x3 grid and
  // be
  // fine. The PeriodicBoundaries grid cannot do so as it must fit the cells to
  // the total length. Thus it would create a 2x2x2 grid and the last particle
  // might result in an out-of-bounds cell index. This constructor call ensures
  // that no assertion/exception in the construction code is hit.
  Grid<GridOptions::PeriodicBoundaries> grid(list, min_cell_length);
}

TEST(max_positions_normal_grid) {
  constexpr int testparticles = 1;
  const double min_cell_length = minimal_cell_length(testparticles);
  using Test::Position;
  Particles list;
  list.insert(Test::smashon(Position{0, 0, 0, -6.2470569610595703125}));
  list.insert(Test::smashon(Position{
      0, 2.5 * min_cell_length, 2.5 * min_cell_length, 8.0611705780029296875}));
  for (int i = 5 * 5 * 5; i; --i) {
    list.insert(Test::smashon(Position(0, 0, 0, 0)));
  }
  // This grid construction uses fragile numbers in the z min/max coordinates,
  // which lead to an index_factor_ that even after one std::nextafter call
  // still generates an out-of-bounds cell index.
  Grid<GridOptions::Normal> grid2(list, testparticles);
}
