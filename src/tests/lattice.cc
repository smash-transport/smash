/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "../include/smash/cxx14compat.h"
#include "../include/smash/fourvector.h"
#include "../include/smash/lattice.h"

using namespace smash;

static std::unique_ptr<RectangularLattice<FourVector>> create_lattice(bool p) {
  const std::array<double, 3> l = {10., 6., 2.};
  const std::array<int, 3> n = {4, 8, 3};
  const std::array<double, 3> origin = {0., 0., 0.};
  return make_unique<RectangularLattice<FourVector>>(
      l, n, origin, p, LatticeUpdate::EveryTimestep);
}

TEST(getters) {
  auto lattice = create_lattice(true);
  COMPARE(lattice->lattice_sizes()[0], 10.);
  COMPARE(lattice->lattice_sizes()[1], 6.);
  COMPARE(lattice->lattice_sizes()[2], 2.);
  COMPARE(lattice->dimensions()[0], 4);
  COMPARE(lattice->dimensions()[1], 8);
  COMPARE(lattice->dimensions()[2], 3);
  COMPARE(lattice->cell_sizes()[0], 2.5);
  COMPARE(lattice->cell_sizes()[1], 0.75);
  COMPARE(lattice->cell_sizes()[2], 2. / 3.);
  VERIFY(lattice->periodic());
}

TEST(reset) {
  auto lattice = create_lattice(true);
  const double dummy = 3.1415926;
  // Fill the lattice with some rubbish, call reset and see if everything is 0
  for (auto &node : *lattice) {
    node = FourVector(dummy, 0.0, 0.0, 0.0);
  }
  lattice->reset();
  for (const auto &node : *lattice) {
    COMPARE(node, FourVector());
  }
}

TEST(out_of_bounds) {
  auto lattice1 = create_lattice(true);
  // For periodic lattice nothing is out of bounds
  VERIFY(!lattice1->out_of_bounds(0, 0, 0));
  VERIFY(!lattice1->out_of_bounds(1, 3, 2));
  VERIFY(!lattice1->out_of_bounds(5, 6, 7));
  VERIFY(!lattice1->out_of_bounds(155, 267, 731));
  // Now remove periodicity
  auto lattice2 = create_lattice(false);
  VERIFY(!lattice2->out_of_bounds(0, 0, 0));
  VERIFY(!lattice2->out_of_bounds(1, 3, 2));
  VERIFY(!lattice2->out_of_bounds(3, 7, 2));
  VERIFY(lattice2->out_of_bounds(1, 2, 3));
  VERIFY(lattice2->out_of_bounds(2, -1, 2));
  VERIFY(lattice2->out_of_bounds(1, 8, 2));
  VERIFY(lattice2->out_of_bounds(5, 5, 0));
  VERIFY(lattice2->out_of_bounds(-2, 5, 0));
  VERIFY(lattice2->out_of_bounds(1, 1, -1));
  VERIFY(lattice2->out_of_bounds(999, 666, 999999));
}

TEST(cell_center) {
  auto lattice = create_lattice(true);
  COMPARE(lattice->cell_center(0, 0, 0).x1(), 1.25);
  COMPARE(lattice->cell_center(0, 0, 0).x2(), 0.375);
  FUZZY_COMPARE(lattice->cell_center(0, 0, 0).x3(), 1.0 / 3.0);
  COMPARE(lattice->cell_center(1, 1, 1).x1(), 3.75);
  COMPARE(lattice->cell_center(1, 1, 1).x2(), 0.375 * 3);
  FUZZY_COMPARE(lattice->cell_center(1, 1, 1).x3(), 1.0);
  // Cell center from 1d index
  auto dims = lattice->dimensions();
  const ThreeVector r1 = lattice->cell_center(1, 3, 2);
  const ThreeVector r2 = lattice->cell_center(1 + dims[0] * (3 + dims[1] * 2));
  VERIFY(r1 == r2);
}

TEST(iterators) {
  auto lattice = create_lattice(false);
  // 1) Check that lattice size is as expected
  COMPARE(lattice->size(), 4u * 8u * 3u);
  // Check node and index accessors
  lattice->node(1, 3, 2) = FourVector(1., 2., 3., 4.);
  COMPARE((*lattice)[8 * 4 * 2 + 4 * 3 + 1], FourVector(1., 2., 3., 4.));

  // 2) Iterate lattice around a point, but take radius so big, that all
  //    the lattice is iterated. Check if really all the lattice was iterated
  ThreeVector r(0.0, 0.0, 0.0);
  double r_cut = 1.0e3;  // This is huge compared to lattice sizes

  int count_nodes = 0;
  lattice->iterate_in_radius(
      r, r_cut, [&](FourVector &, int, int, int) { count_nodes++; });
  COMPARE(count_nodes, 4 * 8 * 3);

  // 3) Iterate around a point, but take r_cut small enough that only one cell
  // is iterated.
  r = ThreeVector(1.26, 0.377, 0.333);  // close to (0,0,0) cell center
  r_cut = 0.2;
  count_nodes = 0;
  lattice->iterate_in_radius(
      r, r_cut, [&](FourVector &, int, int, int) { count_nodes++; });
  COMPARE(count_nodes, 1);

  // 4) Make lattice periodic and repeat 2). Now no node is out of lattice,
  //  because lattice is periodically continued. So expected number of nodes is
  //  (2*rcut)^3 / cell volume
  lattice = create_lattice(true);
  r = ThreeVector(0.0, 0.0, 0.0);
  r_cut = 1.0e2;
  count_nodes = 0;
  lattice->iterate_in_radius(
      r, r_cut, [&](FourVector &, int, int, int) { count_nodes++; });
  int nodes_expect =
      8 * static_cast<int>(std::ceil(r_cut / lattice->cell_sizes()[0] - 0.5) *
                           std::ceil(r_cut / lattice->cell_sizes()[1] - 0.5) *
                           std::ceil(r_cut / lattice->cell_sizes()[2] - 0.5));
  COMPARE(count_nodes, nodes_expect) << count_nodes << " " << nodes_expect;
}

TEST(iterate_in_radius) {
  // 1) Lattice is not periodic
  auto lattice = create_lattice(false);
  // Set all FourVectors to zeros
  lattice->reset();
  FourVector mark = FourVector(1.0, 2.0, 3.0, 4.0);
  ThreeVector r0 = ThreeVector(2.0, 2.0, 1.0), r;
  double r_cut = 1.0;
  lattice->iterate_in_radius(
      r0, r_cut, [&](FourVector &node, int, int, int) { node = mark; });
  /* Iterate all the lattice and check that marked nodes are within r_cut cube,
     while not marked nodes are out of r_cut cube */
  lattice->iterate_sublattice(
      {0, 0, 0}, lattice->dimensions(),
      [&](FourVector &node, int ix, int iy, int iz) {
        r = lattice->cell_center(ix, iy, iz);
        if (std::abs(r[0] - r0[0]) <= r_cut &&
            std::abs(r[1] - r0[1]) <= r_cut &&
            std::abs(r[2] - r0[2]) <= r_cut) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else {
          COMPARE(node, FourVector()) << ix << " " << iy << " " << iz;
        }
      });

  // 2) Lattice is periodic: here d(x1, x2) = |x2 - x1 - int((x2-x1)/l)*l|
  lattice = create_lattice(true);
  lattice->reset();
  r_cut = 2.0;
  lattice->iterate_in_radius(
      r0, r_cut, [&](FourVector &node, int, int, int) { node = mark; });
  double l;
  std::array<double, 3> d;
  lattice->iterate_sublattice(
      {0, 0, 0}, lattice->dimensions(),
      [&](FourVector &node, int ix, int iy, int iz) {
        r = lattice->cell_center(ix, iy, iz);
        for (int i = 0; i < 3; i++) {
          d[i] = r[i] - r0[i];
          l = lattice->lattice_sizes()[i];
          d[i] -= static_cast<int>(d[i] / l) * l;
          d[i] = std::abs(d[i]);
        }
        if (d[0] <= r_cut && d[1] <= r_cut && d[2] <= r_cut) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else {
          COMPARE(node, FourVector()) << ix << " " << iy << " " << iz;
        }
      });
}
