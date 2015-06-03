/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "unittest.h"

#include "../include/cxx14compat.h"
#include "../include/lattice.h"
#include "../include/fourvector.h"

using namespace Smash;

static std::unique_ptr<RectangularLattice<FourVector>> create_lattice(bool p) {
  const std::array<float,3> l = {10.0f, 6.0f, 2.0f};
  const std::array<int,3> n = {4, 8, 3};
  const std::array<float, 3> origin = {0.0f, 0.0f, 0.0f};
  return make_unique<RectangularLattice<FourVector> >(l, n, origin, p,
                                                  LatticeUpdate::EveryTimestep);
}

TEST(getters) {
  auto lattice = create_lattice(true);
  VERIFY(lattice->lattice_sizes()[0] == 10.0f);
  VERIFY(lattice->lattice_sizes()[1] ==  6.0f);
  VERIFY(lattice->lattice_sizes()[2] ==  2.0f);
  VERIFY(lattice->dimensions()[0] == 4);
  VERIFY(lattice->dimensions()[1] == 8);
  VERIFY(lattice->dimensions()[2] == 3);
  VERIFY(lattice->cell_sizes()[0] == 2.5f);
  VERIFY(lattice->cell_sizes()[1] == 0.75f);
  VERIFY(lattice->cell_sizes()[2] == 2.0f/3.0f);
  VERIFY(lattice->periodic());
}

TEST(reset) {
  auto lattice = create_lattice(true);
  const double dummy = 3.1415926;
  // Fill the lattice with some rubbish, call reset and see if everything is 0
  for (auto &node: *lattice) {
    node = FourVector(dummy, 0.0, 0.0, 0.0);
  }
  lattice->reset();
  for (const auto &node: *lattice) {
    VERIFY(node == FourVector());
  }
}

TEST(out_of_bounds) {
  auto lattice1 = create_lattice(true);
  // For periodic lattice nothing is out of bounds
  VERIFY(!lattice1->out_of_bounds(0,0,0));
  VERIFY(!lattice1->out_of_bounds(1,3,2));
  VERIFY(!lattice1->out_of_bounds(5,6,7));
  VERIFY(!lattice1->out_of_bounds(155,267,731));
  // Now remove periodicity
  auto lattice2 = create_lattice(false);
  VERIFY(!lattice2->out_of_bounds(0,0,0));
  VERIFY(!lattice2->out_of_bounds(1,3,2));
  VERIFY(!lattice2->out_of_bounds(3,7,2));
  VERIFY(lattice2->out_of_bounds(1,2,3));
  VERIFY(lattice2->out_of_bounds(2,-1,2));
  VERIFY(lattice2->out_of_bounds(1,8,2));
  VERIFY(lattice2->out_of_bounds(5,5,0));
  VERIFY(lattice2->out_of_bounds(-2,5,0));
  VERIFY(lattice2->out_of_bounds(1,1,-1));
  VERIFY(lattice2->out_of_bounds(999,666,999999));
}

TEST(cell_center) {
  auto lattice = create_lattice(true);
  VERIFY(lattice->cell_center(0,0,0).x1() == 1.25);
  VERIFY(lattice->cell_center(0,0,0).x2() == 0.375);
  // cell center is calculated from float, though ThreeVector contains doubles
  // so accuracy is not better then the float one
  FUZZY_COMPARE(float(lattice->cell_center(0,0,0).x3()), 1.0f/3.0f);
  VERIFY(lattice->cell_center(1,1,1).x1() == 3.75);
  VERIFY(lattice->cell_center(1,1,1).x2() == 0.375*3);
  FUZZY_COMPARE(float(lattice->cell_center(1,1,1).x3()), 1.0f);
}

TEST(iterators) {
  auto lattice = create_lattice(false);
  // 1) Check that lattice size is as expected
  VERIFY(lattice->size() == 4*8*3);
  // Check node and index accessors
  lattice->node(1,3,2) = FourVector(1., 2., 3., 4.);
  VERIFY((*lattice)[8*4*2 + 4*3 + 1] == FourVector(1., 2., 3., 4.));

  // 2) Iterate lattice around a point, but take radius so big, that all
  //    the lattice is iterated. Check if really all the lattice was iterated
  ThreeVector r(0.0, 0.0, 0.0);
  double r_cut = 1.0e3; // This is huge compared to lattice sizes

  int count_nodes = 0;
  lattice->iterate_in_radius(r, r_cut, [&](FourVector&, int, int, int) {
    count_nodes++;
  });
  VERIFY(count_nodes == 4*8*3);

  // 3) Iterate around a point, but take r_cut small enough that only one cell
  // is iterated.
  r = ThreeVector(1.26, 0.377, 0.333);  // close to (0,0,0) cell center
  r_cut = 0.2;
  count_nodes = 0;
  lattice->iterate_in_radius(r, r_cut, [&](FourVector &, int, int, int) {
    count_nodes++;
  });
  VERIFY(count_nodes == 1);

  // 4) Make lattice periodic and repeat 2). Now no node is out of lattice,
  //  because lattice is periodically continued. So expected number of nodes is
  //  (2*rcut)^3 / cell volume
  lattice = create_lattice(true);
  r = ThreeVector(0.0, 0.0, 0.0);
  r_cut = 1.0e2;
  count_nodes = 0;
  lattice->iterate_in_radius(r, r_cut, [&](FourVector &, int, int, int) {
    count_nodes++;
  });
  int nodes_expect = 8 * static_cast<int>(
                     std::ceil(r_cut/lattice->cell_sizes()[0] - 0.5) *
                     std::ceil(r_cut/lattice->cell_sizes()[1] - 0.5) *
                     std::ceil(r_cut/lattice->cell_sizes()[2] - 0.5));
  VERIFY(count_nodes == nodes_expect) << count_nodes << " " << nodes_expect;
}
