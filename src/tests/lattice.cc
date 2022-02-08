/*
 *
 *    Copyright (c) 2015-2021
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <vir/test.h>  // This include has to be first

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
  COMPARE(lattice->n_cells()[0], 4);
  COMPARE(lattice->n_cells()[1], 8);
  COMPARE(lattice->n_cells()[2], 3);
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
  auto dims = lattice->n_cells();
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
  lattice->iterate_in_cube(r, r_cut,
                           [&](FourVector &, int, int, int) { count_nodes++; });
  COMPARE(count_nodes, 4 * 8 * 3);

  // 3) Iterate around a point, but take r_cut small enough that only one cell
  // is iterated.
  r = ThreeVector(1.26, 0.377, 0.333);  // close to (0,0,0) cell center
  r_cut = 0.2;
  count_nodes = 0;
  lattice->iterate_in_cube(r, r_cut,
                           [&](FourVector &, int, int, int) { count_nodes++; });
  COMPARE(count_nodes, 1);

  // 4) Make lattice periodic and repeat 2). Now no node is out of lattice,
  //  because lattice is periodically continued. So expected number of nodes is
  //  (2*rcut)^3 / cell volume
  lattice = create_lattice(true);
  r = ThreeVector(0.0, 0.0, 0.0);
  r_cut = 1.0e2;
  count_nodes = 0;
  lattice->iterate_in_cube(r, r_cut,
                           [&](FourVector &, int, int, int) { count_nodes++; });
  int nodes_expect =
      8 * static_cast<int>(std::ceil(r_cut / lattice->cell_sizes()[0] - 0.5) *
                           std::ceil(r_cut / lattice->cell_sizes()[1] - 0.5) *
                           std::ceil(r_cut / lattice->cell_sizes()[2] - 0.5));
  COMPARE(count_nodes, nodes_expect) << count_nodes << " " << nodes_expect;
}

/*
 * Test that index_left(), index_right(), etc. return correct values; case of a
 * non-periodic lattice.
 */
TEST(neighbor_indeces) {
  const std::array<double, 3> lattice_l = {6.0, 8.0, 7.0};
  const std::array<int, 3> cell_n = {10, 4, 10};
  const std::array<double, 3> origin = {0.0, 0.0, 0.0};
  bool periodicity = false;
  RectangularLattice<double> l(lattice_l, cell_n, origin, periodicity,
                               LatticeUpdate::EveryTimestep);

  // shifts by one lattice spacing in x,y,z directions
  ThreeVector dx(l.cell_sizes()[0], 0.0, 0.0);  // = 6/10 = 0.6
  ThreeVector dy(0.0, l.cell_sizes()[0], 0.0);  // = 8/4 = 1.0
  ThreeVector dz(0.0, 0.0, l.cell_sizes()[0]);  // = 7/10 = 0.7

  // 1) position vector such that shifts according to dx, dy, dz will be within
  // the bounds of the lattice
  ThreeVector r1(4.0, 5.0, 6.0);
  // get indeces for r1
  const int ix1 = std::floor((r1.x1() - (l.origin())[0]) / (l.cell_sizes())[0]);
  const int iy1 = std::floor((r1.x2() - (l.origin())[1]) / (l.cell_sizes())[1]);
  const int iz1 = std::floor((r1.x3() - (l.origin())[2]) / (l.cell_sizes())[2]);
  // get shifted coordinates in the case of a non-periodic lattice
  double ix1_l =
      ((r1.x1() - l.cell_sizes()[0]) < l.origin()[0]) ? ix1 : (ix1 - 1);
  double ix1_r =
      ((r1.x1() + l.cell_sizes()[0]) >= (l.origin()[0] + l.lattice_sizes()[0]))
          ? ix1
          : (ix1 + 1);
  double iy1_d =
      ((r1.x2() - l.cell_sizes()[1]) < l.origin()[1]) ? iy1 : (iy1 - 1);
  double iy1_u =
      ((r1.x2() + l.cell_sizes()[1]) >= (l.origin()[1] + l.lattice_sizes()[1]))
          ? iy1
          : (iy1 + 1);
  double iz1_n =
      ((r1.x3() - l.cell_sizes()[2]) < l.origin()[2]) ? iz1 : (iz1 - 1);
  double iz1_f =
      ((r1.x3() + l.cell_sizes()[2]) >= (l.origin()[2] + l.lattice_sizes()[2]))
          ? iz1
          : (iz1 + 1);
  // compare
  COMPARE(l.cell_center(ix1_l, iy1, iz1),
          l.cell_center(l.index_left(ix1, iy1, iz1)));
  COMPARE(l.cell_center(ix1_r, iy1, iz1),
          l.cell_center(l.index_right(ix1, iy1, iz1)));
  COMPARE(l.cell_center(ix1, iy1_d, iz1),
          l.cell_center(l.index_down(ix1, iy1, iz1)));
  COMPARE(l.cell_center(ix1, iy1_u, iz1),
          l.cell_center(l.index_up(ix1, iy1, iz1)));
  COMPARE(l.cell_center(ix1, iy1, iz1_n),
          l.cell_center(l.index_backward(ix1, iy1, iz1)));
  COMPARE(l.cell_center(ix1, iy1, iz1_f),
          l.cell_center(l.index_forward(ix1, iy1, iz1)));

  // 2) position vector such that some of the shifts according to dx, dy, dz
  // will be beyond the bounds of the lattice; when that happens, the index
  // remains unchanged (for a non-periodic lattice)
  ThreeVector r2(0.25, 6.7, 6.5);
  // get indeces for r2
  const int ix2 = std::floor((r2.x1() - (l.origin())[0]) / (l.cell_sizes())[0]);
  const int iy2 = std::floor((r2.x2() - (l.origin())[1]) / (l.cell_sizes())[1]);
  const int iz2 = std::floor((r2.x3() - (l.origin())[2]) / (l.cell_sizes())[2]);
  // get shifted coordinates in the case of a non-periodic lattice
  double ix_l =
      ((r2.x1() - l.cell_sizes()[0]) < l.origin()[0]) ? ix2 : (ix2 - 1);
  double ix_r =
      ((r2.x1() + l.cell_sizes()[0]) >= (l.origin()[0] + l.lattice_sizes()[0]))
          ? ix2
          : (ix2 + 1);
  double iy_d =
      ((r2.x2() - l.cell_sizes()[1]) < l.origin()[1]) ? iy2 : (iy2 - 1);
  double iy_u =
      ((r2.x2() + l.cell_sizes()[1]) >= (l.origin()[1] + l.lattice_sizes()[1]))
          ? iy2
          : (iy2 + 1);
  double iz_n =
      ((r2.x3() - l.cell_sizes()[2]) < l.origin()[2]) ? iz2 : (iz2 - 1);
  double iz_f =
      ((r2.x3() + l.cell_sizes()[2]) >= (l.origin()[2] + l.lattice_sizes()[2]))
          ? iz2
          : (iz2 + 1);
  // compare
  COMPARE(l.cell_center(ix_l, iy2, iz2),
          l.cell_center(l.index_left(ix2, iy2, iz2)));
  COMPARE(l.cell_center(ix_r, iy2, iz2),
          l.cell_center(l.index_right(ix2, iy2, iz2)));
  COMPARE(l.cell_center(ix2, iy_d, iz2),
          l.cell_center(l.index_down(ix2, iy2, iz2)));
  COMPARE(l.cell_center(ix2, iy_u, iz2),
          l.cell_center(l.index_up(ix2, iy2, iz2)));
  COMPARE(l.cell_center(ix2, iy2, iz_n),
          l.cell_center(l.index_backward(ix2, iy2, iz2)));
  COMPARE(l.cell_center(ix2, iy2, iz_f),
          l.cell_center(l.index_forward(ix2, iy2, iz2)));
}

/*
 * Test that index_left(), index_right(), etc. return correct values; case of a
 * periodic lattice.
 */
TEST(neighbor_indeces_periodic) {
  const std::array<double, 3> lattice_l = {6.0, 8.0, 7.0};
  const std::array<int, 3> cell_n = {10, 4, 10};
  const std::array<double, 3> origin = {0.0, 0.0, 0.0};
  bool periodicity = true;
  RectangularLattice<double> l(lattice_l, cell_n, origin, periodicity,
                               LatticeUpdate::EveryTimestep);

  // shifts in x,y,z directions
  double dx = (l.cell_sizes())[0];  // 6/10 = 0.6
  double dy = (l.cell_sizes())[1];  // 8/4 = 2.0
  double dz = (l.cell_sizes())[2];  // 7/10 = 0.7

  // position vector such that some of the shifts according to dx, dy, dz will
  // be beyond the bounds of the lattice; when that happens, the index remains
  // unchanged
  ThreeVector r(0.25, 6.7, 6.5);
  // get indeces for r
  const int ix = std::floor((r.x1() - (l.origin())[0]) / (l.cell_sizes())[0]);
  const int iy = std::floor((r.x2() - (l.origin())[1]) / (l.cell_sizes())[1]);
  const int iz = std::floor((r.x3() - (l.origin())[2]) / (l.cell_sizes())[2]);

  // get shifted coordinates
  double ix_l =
      ((r.x1() - dx) < l.origin()[0]) ? (ix + l.n_cells()[0] - 1) : (ix - 1);
  double ix_r = ((r.x1() + dx) >= l.lattice_sizes()[0])
                    ? (ix - l.n_cells()[0] + 1)
                    : (ix + 1);
  double iy_d =
      ((r.x2() - dy) < l.origin()[1]) ? (iy + l.n_cells()[1] - 1) : (iy - 1);
  double iy_u = ((r.x2() + dy) >= l.lattice_sizes()[1])
                    ? (iy - l.n_cells()[1] + 1)
                    : (iy + 1);
  double iz_n =
      ((r.x3() - dz) < l.origin()[2]) ? (iz + l.n_cells()[2] - 1) : (iz - 1);
  double iz_f = ((r.x3() + dz) >= l.lattice_sizes()[2])
                    ? (iz - l.n_cells()[2] + 1)
                    : (iz + 1);

  COMPARE(l.cell_center(ix_l, iy, iz), l.cell_center(l.index_left(ix, iy, iz)));
  COMPARE(l.cell_center(ix_r, iy, iz),
          l.cell_center(l.index_right(ix, iy, iz)));
  COMPARE(l.cell_center(ix, iy_d, iz), l.cell_center(l.index_down(ix, iy, iz)));
  COMPARE(l.cell_center(ix, iy_u, iz), l.cell_center(l.index_up(ix, iy, iz)));
  COMPARE(l.cell_center(ix, iy, iz_n),
          l.cell_center(l.index_backward(ix, iy, iz)));
  COMPARE(l.cell_center(ix, iy, iz_f),
          l.cell_center(l.index_forward(ix, iy, iz)));
}

/*
 * Create lattice and fill it with 1/r function.
 * Check if gradient is -\vec{r}/r^3, pay attention to grid edges.
 */
TEST(gradient) {
  // Set a lattice with random parameters, but a relatively fine one.
  const std::array<double, 3> l = {9.0, 7.0, 13.0};
  const std::array<int, 3> n = {50, 60, 70};
  const std::array<double, 3> origin = {-5.2, -4.3, -6.7};
  bool periodicity = false;
  RectangularLattice<double> lat(l, n, origin, periodicity,
                                 LatticeUpdate::EveryTimestep);
  ThreeVector r;
  double d;
  // Fill lattice with 1/r function.
  lat.iterate_sublattice({0, 0, 0}, lat.n_cells(),
                         [&](double &node, int ix, int iy, int iz) {
                           r = lat.cell_center(ix, iy, iz);
                           d = r.abs();
                           node = (d > 0.0) ? (1.0 / d) : 0.0;
                         });
  ThreeVector expected_grad;
  RectangularLattice<ThreeVector> grad_lat(l, n, origin, periodicity,
                                           LatticeUpdate::EveryTimestep);
  lat.compute_gradient_lattice(grad_lat);
  /* Error of the derivative calculation is proportional to the
     lattice spacing squared and to the third derivative of the function.
     Third derivative of 1/r is getting larger for small r.
     That is why I only test r larger than some r0. For them the expected
     precision is better, and test with better precision has better
     capability to catch potential bugs. Furthermore, precision at the
     edges of the non-periodic grid is proportional to lattice spacing,
     not to lattice spacing squared.
  */
  grad_lat.iterate_sublattice(
      {0, 0, 0}, grad_lat.n_cells(),
      [&](ThreeVector &node, int ix, int iy, int iz) {
        r = grad_lat.cell_center(ix, iy, iz);
        d = r.abs();
        if (d > 2.0) {
          expected_grad = -r / (d * d * d);
          COMPARE_RELATIVE_ERROR(node.x1(), expected_grad.x1(), 6.e-2)
              << "node: (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << d;
          COMPARE_RELATIVE_ERROR(node.x2(), expected_grad.x2(), 6.e-2)
              << "node: (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << d;
          COMPARE_RELATIVE_ERROR(node.x3(), expected_grad.x3(), 6.e-2)
              << "node: (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << d;
        }
      });
}

/*
 * Create periodic lattice and fill it with
 * cos(2 pi x/lx) cos(2 pi y/ly) cos(2 pi z/lz) function.
 * Function is periodic and limited, error of computed gradient should be
 * proportional to lattice spacing squared everywhere, including lattice edges,
 * so computational errors should be smaller than in the "gradient" test.
 * This test checks the gradient on periodic lattice. The only difference
 * with non-periodic is treatment of derivatives on the edges.
 */
TEST(gradient_periodic) {
  const std::array<double, 3> l = {9.0, 7.0, 13.0};
  const std::array<int, 3> n = {50, 60, 70};
  const std::array<double, 3> origin = {-5.2, -4.3, -6.7};
  bool periodicity = true;
  RectangularLattice<double> lat(l, n, origin, periodicity,
                                 LatticeUpdate::EveryTimestep);
  ThreeVector r;
  // Fill lattice with (2 pi x/lx) cos(2 pi y/ly) cos(2 pi z/lz) function.
  lat.iterate_sublattice({0, 0, 0}, lat.n_cells(),
                         [&](double &node, int ix, int iy, int iz) {
                           r = lat.cell_center(ix, iy, iz);
                           node = std::cos(2 * M_PI * r.x1() / l[0]) *
                                  std::cos(2 * M_PI * r.x2() / l[1]) *
                                  std::cos(2 * M_PI * r.x3() / l[2]);
                         });
  ThreeVector expected_grad;
  RectangularLattice<ThreeVector> grad_lat(l, n, origin, periodicity,
                                           LatticeUpdate::EveryTimestep);
  lat.compute_gradient_lattice(grad_lat);
  grad_lat.iterate_sublattice(
      {0, 0, 0}, grad_lat.n_cells(),
      [&](ThreeVector &node, int ix, int iy, int iz) {
        r = grad_lat.cell_center(ix, iy, iz);
        expected_grad.set_x1(-2 * M_PI / l[0] *
                             std::sin(2 * M_PI * r.x1() / l[0]) *
                             std::cos(2 * M_PI * r.x2() / l[1]) *
                             std::cos(2 * M_PI * r.x3() / l[2]));
        expected_grad.set_x2(-2 * M_PI / l[1] *
                             std::cos(2 * M_PI * r.x1() / l[0]) *
                             std::sin(2 * M_PI * r.x2() / l[1]) *
                             std::cos(2 * M_PI * r.x3() / l[2]));
        expected_grad.set_x3(-2 * M_PI / l[2] *
                             std::cos(2 * M_PI * r.x1() / l[0]) *
                             std::cos(2 * M_PI * r.x2() / l[1]) *
                             std::sin(2 * M_PI * r.x3() / l[2]));

        COMPARE_RELATIVE_ERROR(node.x1(), expected_grad.x1(), 3.e-3)
            << "node: (" << ix << ", " << iy << ", " << iz << ")";
        COMPARE_RELATIVE_ERROR(node.x2(), expected_grad.x2(), 3.e-3)
            << "node: (" << ix << ", " << iy << ", " << iz << ")";
        COMPARE_RELATIVE_ERROR(node.x3(), expected_grad.x3(), 3.e-3)
            << "node: (" << ix << ", " << iy << ", " << iz << ")";
      });
}

/*
 * Test gradient for 2x2x2 lattice. The test is that it doesn't segfault.
 */
TEST(gradient_2x2x2lattice) {
  const std::array<double, 3> l = {9.0, 7.0, 13.0};
  const std::array<int, 3> n = {2, 2, 2};
  const std::array<double, 3> origin = {-5.2, -4.3, -6.7};
  bool periodicity = false;
  RectangularLattice<double> lat(l, n, origin, periodicity,
                                 LatticeUpdate::EveryTimestep);
  ThreeVector r;
  // Fill lattice with (2 pi x/lx) cos(2 pi y/ly) cos(2 pi z/lz) function.
  lat.iterate_sublattice({0, 0, 0}, lat.n_cells(),
                         [&](double &node, int ix, int iy, int iz) {
                           r = lat.cell_center(ix, iy, iz);
                           node = std::cos(2 * M_PI * r.x1() / l[0]) *
                                  std::cos(2 * M_PI * r.x2() / l[1]) *
                                  std::cos(2 * M_PI * r.x3() / l[2]);
                         });
  RectangularLattice<ThreeVector> grad_lat(l, n, origin, periodicity,
                                           LatticeUpdate::EveryTimestep);
  lat.compute_gradient_lattice(grad_lat);
}

// If one of the dimensions is 1, then 3D gradient calculation is impossible
TEST_CATCH(gradient_impossible_lattice, std::runtime_error) {
  const std::array<double, 3> l = {9.0, 7.0, 13.0};
  const std::array<int, 3> n = {5, 42, 1};
  const std::array<double, 3> origin = {-5.2, -4.3, -6.7};
  bool periodicity = false;
  RectangularLattice<double> lat(l, n, origin, periodicity,
                                 LatticeUpdate::EveryTimestep);
  RectangularLattice<ThreeVector> grad_lat(l, n, origin, periodicity,
                                           LatticeUpdate::EveryTimestep);
  lat.compute_gradient_lattice(grad_lat);
}

/*
 * Based on the compute_gradient test. Create two lattices of FourVectors: "old"
 * and "new" one. Fill the "old" lattice with zeros, and fill while for the
 * "new" lattice fill each components of its FourVectors with 1/r function.
 * Check if time derivative is 1/r. Check if gradient is -\vec{r}/r^3, pay
 * attention to grid edges. Case of a non-periodic lattice.
 */
TEST(compute_four_gradient_lattice) {
  // Set a lattice with random parameters, but a relatively fine one.
  const std::array<double, 3> l = {9.0, 7.0, 13.0};
  const std::array<int, 3> n = {50, 60, 70};
  const std::array<double, 3> origin = {-5.2, -4.3, -6.7};
  bool periodicity = false;
  RectangularLattice<FourVector> old_lat(l, n, origin, periodicity,
                                         LatticeUpdate::EveryTimestep);
  RectangularLattice<FourVector> new_lat(l, n, origin, periodicity,
                                         LatticeUpdate::EveryTimestep);

  // Fill the lattice old_lat with 0s
  old_lat.reset();

  // Fill the lattice new_lat with 1/r function.
  ThreeVector position;
  double r;
  new_lat.iterate_sublattice({0, 0, 0}, new_lat.n_cells(),
                             [&](FourVector &node, int ix, int iy, int iz) {
                               position = new_lat.cell_center(ix, iy, iz);
                               r = position.abs();
                               // get inverse of r unless r is zero:
                               double inv_r = (r > 0.0) ? (1.0 / r) : 0.0;

                               // filling each element of the four-vector at the
                               // node with the same function
                               node = FourVector(inv_r, inv_r, inv_r, inv_r);
                             });

  // time step for the time derivative
  const double time_step = 0.1;
  // create the lattice storing the four-gradient
  RectangularLattice<std::array<FourVector, 4>> grad_lat(
      l, n, origin, periodicity, LatticeUpdate::EveryTimestep);
  // compute the four-gradient
  new_lat.compute_four_gradient_lattice(old_lat, time_step, grad_lat);

  double expected_dt = 0.0;
  ThreeVector expected_grad(0.0, 0.0, 0.0);
  /* Error of the derivative calculation is proportional to the
     lattice spacing squared and to the third derivative of the function.
     Third derivative of 1/r is getting larger for small r.
     That is why I only test r larger than some r0. For them the expected
     precision is better, and test with better precision has better
     capability to catch potential bugs. Furthermore, precision at the
     edges of the non-periodic grid is proportional to lattice spacing,
     not to lattice spacing squared.
  */
  grad_lat.iterate_sublattice(
      {0, 0, 0}, grad_lat.n_cells(),
      [&](std::array<FourVector, 4> &node, int ix, int iy, int iz) {
        position = grad_lat.cell_center(ix, iy, iz);
        r = position.abs();
        // get inverse of r unless r is zero:
        double inv_r = (r > 0.0) ? (1.0 / r) : 0.0;

        // the time derivative can be realiably compared for any r:
        expected_dt = inv_r / time_step;
        // dA0 / dt
        COMPARE_RELATIVE_ERROR(node[0].x0(), expected_dt, 6.e-2)
            << "dA0 / dt \t node: (" << ix << ", " << iy << ", " << iz
            << "), |r| = " << r;
        // dA1 / dt
        COMPARE_RELATIVE_ERROR(node[0].x1(), expected_dt, 6.e-2)
            << "dA1 / dt \t node: (" << ix << ", " << iy << ", " << iz
            << "), |r| = " << r;
        // dA2 / dt
        COMPARE_RELATIVE_ERROR(node[0].x2(), expected_dt, 6.e-2)
            << "dA2 / dt \t node: (" << ix << ", " << iy << ", " << iz
            << "), |r| = " << r;
        // dA3 / dt
        COMPARE_RELATIVE_ERROR(node[0].x3(), expected_dt, 6.e-2)
            << "dA3 / dt \t node: (" << ix << ", " << iy << ", " << iz
            << "), |r| = " << r;

        // the spatial derivatives we only compare for r > 2.0:
        if (r > 2.0) {
          expected_grad = -position / (r * r * r);
          // grad x:
          COMPARE_RELATIVE_ERROR(node[1].x0(), expected_grad.x1(), 6.e-2)
              << "dA0 / dx \t node[1].x0(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
          COMPARE_RELATIVE_ERROR(node[1].x1(), expected_grad.x1(), 6.e-2)
              << "dA1 / dx \t node[1].x1(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
          COMPARE_RELATIVE_ERROR(node[1].x2(), expected_grad.x1(), 6.e-2)
              << "dA2 / dx \t node[1].x2(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
          COMPARE_RELATIVE_ERROR(node[1].x3(), expected_grad.x1(), 6.e-2)
              << "dA3 / dx \t node[1].x3(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
          // grad y:
          COMPARE_RELATIVE_ERROR(node[2].x0(), expected_grad.x2(), 6.e-2)
              << "dA0 / dy \t node[2].x0(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
          COMPARE_RELATIVE_ERROR(node[2].x1(), expected_grad.x2(), 6.e-2)
              << "dA1 / dx \t node[2].x1(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
          COMPARE_RELATIVE_ERROR(node[2].x2(), expected_grad.x2(), 6.e-2)
              << "dA2 / dx \t node[2].x2(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
          COMPARE_RELATIVE_ERROR(node[2].x3(), expected_grad.x2(), 6.e-2)
              << "dA3 / dx \t node[2].x3(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
          // grad z:
          COMPARE_RELATIVE_ERROR(node[3].x0(), expected_grad.x3(), 6.e-2)
              << "dA0 / dz \t node[3].x0(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
          COMPARE_RELATIVE_ERROR(node[3].x1(), expected_grad.x3(), 6.e-2)
              << "dA1 / dz \t node[3].x1(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
          COMPARE_RELATIVE_ERROR(node[3].x2(), expected_grad.x3(), 6.e-2)
              << "dA2 / dz \t node[3].x2(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
          COMPARE_RELATIVE_ERROR(node[3].x3(), expected_grad.x3(), 6.e-2)
              << "dA3 / dz \t node[3].x3(): (" << ix << ", " << iy << ", " << iz
              << "), |r| = " << r;
        }
      });
}

/*
 * Based on the compute_gradient_periodic test. Create two lattices of
 * FourVectors: "old" and "new" one. Fill the "old" lattice with zeros, and
 * fill while for the "new" lattice fill each components of its FourVectors with
 * 1/r function. Check if time derivative is 1/r. Check if gradient is
 * -\vec{r}/r^3, pay attention to grid edges. Case of a periodic lattice.
 */
TEST(compute_four_gradient_lattice_periodic) {
  const std::array<double, 3> l = {9.0, 7.0, 13.0};
  const std::array<int, 3> n = {50, 60, 70};
  const std::array<double, 3> origin = {-5.2, -4.3, -6.7};
  bool periodicity = true;
  RectangularLattice<FourVector> old_lat(l, n, origin, periodicity,
                                         LatticeUpdate::EveryTimestep);
  RectangularLattice<FourVector> new_lat(l, n, origin, periodicity,
                                         LatticeUpdate::EveryTimestep);

  // Fill the lattice old_lat with 0s
  old_lat.reset();

  // Fill the lattice new_lat with a function periodic in x, y, and z
  ThreeVector pos;
  new_lat.iterate_sublattice({0, 0, 0}, new_lat.n_cells(),
                             [&](FourVector &node, int ix, int iy, int iz) {
                               pos = new_lat.cell_center(ix, iy, iz);
                               double func =
                                   std::cos(2 * M_PI * pos.x1() / l[0]) *
                                   std::cos(2 * M_PI * pos.x2() / l[1]) *
                                   std::cos(2 * M_PI * pos.x3() / l[2]);

                               // filling each element of the four-vector at the
                               // node with the same function
                               node = FourVector(func, func, func, func);
                             });

  // time step for the time derivative
  const double time_step = 0.1;
  // create the lattice storing the four-gradient
  RectangularLattice<std::array<FourVector, 4>> grad_lat(
      l, n, origin, periodicity, LatticeUpdate::EveryTimestep);
  // compute the four-gradient
  new_lat.compute_four_gradient_lattice(old_lat, time_step, grad_lat);

  double expected_dt = 0.0;
  ThreeVector expected_grad(0.0, 0.0, 0.0);
  /* Error of the derivative calculation is proportional to the
     lattice spacing squared and to the third derivative of the function.
     Third derivative of 1/r is getting larger for small r.
     That is why I only test r larger than some r0. For them the expected
     precision is better, and test with better precision has better
     capability to catch potential bugs. Furthermore, precision at the
     edges of the non-periodic grid is proportional to lattice spacing,
     not to lattice spacing squared.
  */
  grad_lat.iterate_sublattice(
      {0, 0, 0}, grad_lat.n_cells(),
      [&](std::array<FourVector, 4> &node, int ix, int iy, int iz) {
        pos = grad_lat.cell_center(ix, iy, iz);
        double func = std::cos(2 * M_PI * pos.x1() / l[0]) *
                      std::cos(2 * M_PI * pos.x2() / l[1]) *
                      std::cos(2 * M_PI * pos.x3() / l[2]);

        // the time derivatives:
        expected_dt = func / time_step;
        // dA0 / dt
        COMPARE_RELATIVE_ERROR(node[0].x0(), expected_dt, 3.e-3)
            << "node: (" << ix << ", " << iy << ", " << iz << ")";
        // dA1 / dt
        COMPARE_RELATIVE_ERROR(node[0].x1(), expected_dt, 3.e-3)
            << "node: (" << ix << ", " << iy << ", " << iz << ")";
        // dA2 / dt
        COMPARE_RELATIVE_ERROR(node[0].x2(), expected_dt, 3.e-3)
            << "node: (" << ix << ", " << iy << ", " << iz << ")";
        // dA3 / dt
        COMPARE_RELATIVE_ERROR(node[0].x3(), expected_dt, 3.e-3)
            << "node: (" << ix << ", " << iy << ", " << iz << ")";

        // the spatial derivatives:
        expected_grad.set_x1(-2 * M_PI / l[0] *
                             std::sin(2 * M_PI * pos.x1() / l[0]) *
                             std::cos(2 * M_PI * pos.x2() / l[1]) *
                             std::cos(2 * M_PI * pos.x3() / l[2]));
        expected_grad.set_x2(-2 * M_PI / l[1] *
                             std::cos(2 * M_PI * pos.x1() / l[0]) *
                             std::sin(2 * M_PI * pos.x2() / l[1]) *
                             std::cos(2 * M_PI * pos.x3() / l[2]));
        expected_grad.set_x3(-2 * M_PI / l[2] *
                             std::cos(2 * M_PI * pos.x1() / l[0]) *
                             std::cos(2 * M_PI * pos.x2() / l[1]) *
                             std::sin(2 * M_PI * pos.x3() / l[2]));

        // grad x:
        COMPARE_RELATIVE_ERROR(node[1].x0(), expected_grad.x1(), 3.e-3)
            << "node[1].x0(): (" << ix << ", " << iy << ", " << iz << ")";
        COMPARE_RELATIVE_ERROR(node[1].x1(), expected_grad.x1(), 3.e-3)
            << "node[1].x1(): (" << ix << ", " << iy << ", " << iz << ")";
        COMPARE_RELATIVE_ERROR(node[1].x2(), expected_grad.x1(), 3.e-3)
            << "node[1].x2(): (" << ix << ", " << iy << ", " << iz << ")";
        COMPARE_RELATIVE_ERROR(node[1].x3(), expected_grad.x1(), 3.e-3)
            << "node[1].x3(): (" << ix << ", " << iy << ", " << iz << ")";
        // grad y:
        COMPARE_RELATIVE_ERROR(node[2].x0(), expected_grad.x2(), 3.e-3)
            << "node[2].x0(): (" << ix << ", " << iy << ", " << iz << ")";
        COMPARE_RELATIVE_ERROR(node[2].x1(), expected_grad.x2(), 3.e-3)
            << "node[2].x1(): (" << ix << ", " << iy << ", " << iz << ")";
        COMPARE_RELATIVE_ERROR(node[2].x2(), expected_grad.x2(), 3.e-3)
            << "node[2].x2(): (" << ix << ", " << iy << ", " << iz << ")";
        COMPARE_RELATIVE_ERROR(node[2].x3(), expected_grad.x2(), 3.e-3)
            << "node[2].x3(): (" << ix << ", " << iy << ", " << iz << ")";
        // grad z:
        COMPARE_RELATIVE_ERROR(node[3].x0(), expected_grad.x3(), 3.e-3)
            << "node[3].x0(): (" << ix << ", " << iy << ", " << iz << ")";
        COMPARE_RELATIVE_ERROR(node[3].x1(), expected_grad.x3(), 3.e-3)
            << "node[3].x1(): (" << ix << ", " << iy << ", " << iz << ")";
        COMPARE_RELATIVE_ERROR(node[3].x2(), expected_grad.x3(), 3.e-3)
            << "node[3].x2(): (" << ix << ", " << iy << ", " << iz << ")";
        COMPARE_RELATIVE_ERROR(node[3].x3(), expected_grad.x3(), 3.e-3)
            << "node[3].x3(): (" << ix << ", " << iy << ", " << iz << ")";
      });
}

/*
 * Test fourgradient for 2x2x2 lattice. The test is that it doesn't segfault.
 */
TEST(fourgradient_2x2x2lattice) {
  const std::array<double, 3> l = {9.0, 7.0, 13.0};
  const std::array<int, 3> n = {2, 2, 2};
  const std::array<double, 3> origin = {-5.2, -4.3, -6.7};
  bool periodicity = false;
  RectangularLattice<FourVector> old_lat(l, n, origin, periodicity,
                                         LatticeUpdate::EveryTimestep);
  RectangularLattice<FourVector> new_lat(l, n, origin, periodicity,
                                         LatticeUpdate::EveryTimestep);

  // Fill the lattice old_lat with 0s
  old_lat.reset();

  // Fill the lattice new_lat with 1/r function.
  ThreeVector position;
  double r;
  new_lat.iterate_sublattice({0, 0, 0}, new_lat.n_cells(),
                             [&](FourVector &node, int ix, int iy, int iz) {
                               position = new_lat.cell_center(ix, iy, iz);
                               r = position.abs();
                               // get inverse of r unless r is zero:
                               double inv_r = (r > 0.0) ? (1.0 / r) : 0.0;

                               // filling each element of the four-vector at the
                               // node with the same function
                               node = FourVector(inv_r, inv_r, inv_r, inv_r);
                             });

  // time step for the time derivative
  const double time_step = 0.1;
  // create the lattice storing the four-gradient
  RectangularLattice<std::array<FourVector, 4>> grad_lat(
      l, n, origin, periodicity, LatticeUpdate::EveryTimestep);
  // compute the four-gradient
  new_lat.compute_four_gradient_lattice(old_lat, time_step, grad_lat);
}

// If one of the dimensions is 1, then 3D gradient calculation is impossible
TEST_CATCH(fourgradient_impossible_lattice, std::runtime_error) {
  const std::array<double, 3> l = {9.0, 7.0, 13.0};
  const std::array<int, 3> n = {5, 42, 1};
  const std::array<double, 3> origin = {-5.2, -4.3, -6.7};
  bool periodicity = false;

  RectangularLattice<FourVector> old_lat(l, n, origin, periodicity,
                                         LatticeUpdate::EveryTimestep);
  RectangularLattice<FourVector> new_lat(l, n, origin, periodicity,
                                         LatticeUpdate::EveryTimestep);

  // time step for the time derivative
  const double time_step = 0.1;
  // create the lattice storing the four-gradient
  RectangularLattice<std::array<FourVector, 4>> grad_lat(
      l, n, origin, periodicity, LatticeUpdate::EveryTimestep);
  // compute the four-gradient
  new_lat.compute_four_gradient_lattice(old_lat, time_step, grad_lat);
}

TEST(iterate_in_cube) {
  // 1) Lattice is not periodic
  auto lattice = create_lattice(false);
  // Set all FourVectors to zeros
  lattice->reset();
  FourVector mark = FourVector(1.0, 2.0, 3.0, 4.0);
  ThreeVector r0 = ThreeVector(2.0, 2.0, 1.0), r;
  double r_cut = 1.0;
  lattice->iterate_in_cube(
      r0, r_cut, [&](FourVector &node, int, int, int) { node = mark; });
  /* Iterate all the lattice and check that marked nodes are within r_cut cube,
     while not marked nodes are out of r_cut cube */
  lattice->iterate_sublattice(
      {0, 0, 0}, lattice->n_cells(),
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
  lattice->iterate_in_cube(
      r0, r_cut, [&](FourVector &node, int, int, int) { node = mark; });
  double l;
  std::array<double, 3> d;
  // auxiliary structures to choose correct distance
  std::array<double, 3> d1;
  std::array<double, 3> d2;
  lattice->iterate_sublattice(
      {0, 0, 0}, lattice->n_cells(),
      [&](FourVector &node, int ix, int iy, int iz) {
        r = lattice->cell_center(ix, iy, iz);
        for (int i = 0; i < 3; i++) {
          d1[i] = r[i] - r0[i];
          l = lattice->lattice_sizes()[i];
          d1[i] -= static_cast<int>(d1[i] / l) * l;
          d1[i] = std::abs(d1[i]);

          d2[i] = r[i] - r0[i] - lattice->lattice_sizes()[i];
          l = lattice->lattice_sizes()[i];
          d2[i] -= static_cast<int>(d2[i] / l) * l;
          d2[i] = std::abs(d2[i]);

          d[i] = (d1[i] < d2[i]) ? d1[i] : d2[i];
        }
        if (d[0] <= r_cut && d[1] <= r_cut && d[2] <= r_cut) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else {
          COMPARE(node, FourVector()) << ix << " " << iy << " " << iz;
        }
      });
}

TEST(iterate_in_rectangle) {
  // 1) Lattice is not periodic
  auto lattice = create_lattice(false);
  // Set all FourVectors to zeros
  lattice->reset();
  // this fourvector will play a role of a mark
  FourVector mark = FourVector(1.0, 2.0, 3.0, 4.0);
  // position vector around which we iterate
  ThreeVector r0 = ThreeVector(2.0, 2.0, 1.0), r;
  // iteration range
  double range = 2.0;

  const std::array<double, 3> rectangle = {range * (lattice->cell_sizes())[0],
                                           range * (lattice->cell_sizes())[1],
                                           range * (lattice->cell_sizes())[2]};

  // "mark" the nodes within the iteration rectangle
  lattice->iterate_in_rectangle(
      r0, rectangle, [&](FourVector &node, int, int, int) { node = mark; });

  /* Iterate all the lattice and check that marked nodes are within the
     rectangle, while not marked nodes are out of the rectangle */
  lattice->iterate_sublattice({0, 0, 0}, lattice->n_cells(),
                              [&](FourVector &node, int ix, int iy, int iz) {
                                r = lattice->cell_center(ix, iy, iz);
                                if (std::abs(r[0] - r0[0]) <= rectangle[0] &&
                                    std::abs(r[1] - r0[1]) <= rectangle[1] &&
                                    std::abs(r[2] - r0[2]) <= rectangle[2]) {
                                  COMPARE(node, mark)
                                      << ix << " " << iy << " " << iz
                                      << "\t periodic = false";
                                } else {
                                  COMPARE(node, FourVector())
                                      << ix << " " << iy << " " << iz
                                      << "\t periodic = false";
                                }
                              });

  // 2) Lattice is periodic: here d(x1, x2) = shortest distance between points,
  // which may be "through the edges"
  lattice = create_lattice(true);
  lattice->reset();

  // "mark" nodes within the iteration rectangle
  lattice->iterate_in_rectangle(
      r0, rectangle, [&](FourVector &node, int, int, int) { node = mark; });
  double l;
  std::array<double, 3> d;
  // auxiliary structures to choose correct distance
  std::array<double, 3> d1;
  std::array<double, 3> d2;
  /* Iterate all the lattice and check that marked nodes are within the
     rectangle, while not marked nodes are out of the rectangle */
  lattice->iterate_sublattice(
      {0, 0, 0}, lattice->n_cells(),
      [&](FourVector &node, int ix, int iy, int iz) {
        r = lattice->cell_center(ix, iy, iz);
        for (int i = 0; i < 3; i++) {
          d1[i] = r[i] - r0[i];
          l = lattice->lattice_sizes()[i];
          d1[i] -= static_cast<int>(d1[i] / l) * l;
          d1[i] = std::abs(d1[i]);

          d2[i] = r[i] - r0[i] - lattice->lattice_sizes()[i];
          l = lattice->lattice_sizes()[i];
          d2[i] -= static_cast<int>(d2[i] / l) * l;
          d2[i] = std::abs(d2[i]);

          // d[i] = d1[i];
          d[i] = (d1[i] < d2[i]) ? d1[i] : d2[i];
        }
        if (d[0] <= rectangle[0] && d[1] <= rectangle[1] &&
            d[2] <= rectangle[2]) {
          COMPARE(node, mark)
              << ix << " " << iy << " " << iz << "\t periodic = true"
              << "\n d[0] = " << d[0] << "\td[1] = " << d[1]
              << "\td[2] = " << d[2];
        } else {
          COMPARE(node, FourVector())
              << ix << " " << iy << " " << iz << "\t periodic = true"
              << "\n   lattice n_cells = (" << lattice->n_cells()[0] << ", "
              << lattice->n_cells()[1] << ", " << lattice->n_cells()[2] << ")"
              << "\nlattice cell_sizes = (" << lattice->cell_sizes()[0] << ", "
              << lattice->cell_sizes()[1] << ", " << lattice->cell_sizes()[2]
              << ")"
              << "\n     lattice sizes = (" << lattice->lattice_sizes()[0]
              << ", " << lattice->lattice_sizes()[1] << ", "
              << lattice->lattice_sizes()[2] << ")"
              << "\n node = " << node << "\n   r0 = " << r0 << "\n    r = " << r
              << "\n   d1 = (" << d1[0] << ", " << d1[1] << ", " << d1[2] << ")"
              << "\n   d2 = (" << d2[0] << ", " << d2[1] << ", " << d2[2] << ")"
              << "\n    d = (" << d[0] << ", " << d[1] << ", " << d[2] << ")";
        }
      });
}

TEST(iterate_nearest_neighbors) {
  // 1) Lattice is not periodic
  auto lattice = create_lattice(false);
  // Set all FourVectors to zeros
  lattice->reset();
  // this fourvector will play a role of a mark
  FourVector mark = FourVector(1.0, 2.0, 3.0, 4.0);
  // position vector around which we iterate
  ThreeVector r0 = ThreeVector(2.0, 0.0, 1.0), r;

  // "mark" the nearest neighbors
  lattice->iterate_nearest_neighbors(
      r0, [&](FourVector &node, int, int) { node = mark; });

  std::array<double, 3> origin = lattice->origin();
  std::array<double, 3> cell_sizes = lattice->cell_sizes();

  // get the index of the cell to which r0 points
  int ix_0 = std::floor((r0.x1() - origin[0]) / cell_sizes[0]);
  int iy_0 = std::floor((r0.x2() - origin[1]) / cell_sizes[1]);
  int iz_0 = std::floor((r0.x3() - origin[2]) / cell_sizes[2]);
  int center_index = lattice->index1d(ix_0, iy_0, iz_0);
  int left_index = lattice->index_left(ix_0, iy_0, iz_0);
  int right_index = lattice->index_right(ix_0, iy_0, iz_0);
  int down_index = lattice->index_down(ix_0, iy_0, iz_0);
  int up_index = lattice->index_up(ix_0, iy_0, iz_0);
  int backward_index = lattice->index_backward(ix_0, iy_0, iz_0);
  int forward_index = lattice->index_forward(ix_0, iy_0, iz_0);

  int current_index = 0;
  // Iterate all the lattice and check that the center node and the neighboring
  // nodes are marked
  lattice->iterate_sublattice(
      {0, 0, 0}, lattice->n_cells(),
      [&](FourVector &node, int ix, int iy, int iz) {
        current_index = lattice->index1d(ix, iy, iz);
        if (current_index == center_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == left_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == right_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == down_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == up_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == backward_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == forward_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else {
          COMPARE(node, FourVector()) << ix << " " << iy << " " << iz;
        }
      });

  // 2) Lattice is periodic:
  lattice = create_lattice(true);
  lattice->reset();

  // "mark" the nearest neighbors
  lattice->iterate_nearest_neighbors(
      r0, [&](FourVector &node, int, int) { node = mark; });

  origin = lattice->origin();
  cell_sizes = lattice->cell_sizes();

  // get the index of the cell to which r0 points
  ix_0 = std::floor((r0.x1() - origin[0]) / cell_sizes[0]);
  iy_0 = std::floor((r0.x2() - origin[1]) / cell_sizes[1]);
  iz_0 = std::floor((r0.x3() - origin[2]) / cell_sizes[2]);
  center_index = lattice->index1d(ix_0, iy_0, iz_0);
  left_index = lattice->index_left(ix_0, iy_0, iz_0);
  right_index = lattice->index_right(ix_0, iy_0, iz_0);
  down_index = lattice->index_down(ix_0, iy_0, iz_0);
  up_index = lattice->index_up(ix_0, iy_0, iz_0);
  backward_index = lattice->index_backward(ix_0, iy_0, iz_0);
  forward_index = lattice->index_forward(ix_0, iy_0, iz_0);

  // Iterate all the lattice and check that the center node and the neighboring
  // nodes are marked
  lattice->iterate_sublattice(
      {0, 0, 0}, lattice->n_cells(),
      [&](FourVector &node, int ix, int iy, int iz) {
        current_index = lattice->index1d(ix, iy, iz);
        if (current_index == center_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == left_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == right_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == down_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == up_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == backward_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else if (current_index == forward_index) {
          COMPARE(node, mark) << ix << " " << iy << " " << iz;
        } else {
          COMPARE(node, FourVector()) << ix << " " << iy << " " << iz;
        }
      });
}

TEST(copy_constructor) {
  const std::array<double, 3> l = {10., 6., 2.};
  const std::array<int, 3> n = {4, 8, 3};
  const std::array<double, 3> origin = {0., 3., 7.};
  RectangularLattice<double> lattice(l, n, origin, false,
                                     LatticeUpdate::EveryTimestep);
  int i = 0;
  for (auto &node : lattice) {
    node = static_cast<double>(i++);
  }

  // Use copy constructor to create a new lattice
  RectangularLattice<double> lattice2 = lattice;

  // Make sure lattice2 is indeed identical to lattice
  VERIFY(lattice2.identical_to_lattice(&lattice));

  lattice.iterate_sublattice({0, 0, 0}, lattice.n_cells(),
                             [&](double &node, int ix, int iy, int iz) {
                               COMPARE(node, lattice2.node(ix, iy, iz));
                             });
}

static double integrand(ThreeVector pos, double &value, ThreeVector point) {
  return value / ((pos - point).abs());
}

TEST(integrate_volume) {
  int ncells = 20;
  const std::array<double, 3> l = {20., 20., 20.};
  const std::array<int, 3> n = {ncells, ncells, ncells};
  const std::array<double, 3> origin = {-10, -10, -10};
  RectangularLattice<double> lattice(l, n, origin, false,
                                     LatticeUpdate::EveryTimestep);
  const double radius = 8.;
  const double density = 2.0;
  const ThreeVector r0 = {1., 1., -1.};
  for (int i = 0; i < ncells * ncells * ncells; i++) {
    ThreeVector r = lattice.cell_center(i);
    if ((r - r0).abs() <= radius) {
      lattice[i] = density * ((r - r0).sqr());
    } else {
      lattice[i] = 0.;
    }
  }
  // The field is 2*|r-r0|^2. The integral of |r-r0|^2/|r-r0| over a
  // circle around r0 should be 2*pi*R^4
  double integral = 0;
  lattice.integrate_volume(integral, integrand, radius, r0);
  COMPARE_RELATIVE_ERROR(integral, 2 * M_PI * std::pow(radius, 4), 0.03);
}
