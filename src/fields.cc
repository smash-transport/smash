/*
 *
 *    Copyright (c) 2013-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/fields.h"

namespace smash {

void update_fields_lattice(
    RectangularLattice<FieldsOnLattice> *fields_lat,
    RectangularLattice<FourVector> *old_fields,
    RectangularLattice<FourVector> *new_fields,
    RectangularLattice<std::array<FourVector, 4>> *fields_four_grad_lattice,
    DensityLattice *jmuB_lat, const LatticeUpdate fields_lat_update,
    const Potentials &potentials, const double time_step) {
  // Do not proceed if lattice does not exists/update not required
  if (fields_lat == nullptr || fields_lat->when_update() != fields_lat_update) {
    return;
  }
  // get the number of nodes on the fields lattice
  const std::array<int, 3> lattice_n_cells = fields_lat->n_cells();
  const int number_of_nodes =
      lattice_n_cells[0] * lattice_n_cells[1] * lattice_n_cells[2];

  /*
   * Take the provided FieldsOnLattice lattice, fields_lat, and use the
   * information about the fields to populate the lattice of A^mu FourVectors at
   * t0, old_fields.
   */
  for (int i = 0; i < number_of_nodes; i++) {
    old_fields->assign_value(i, ((*fields_lat)[i]).A_mu());
  }

  /*
   * Update the fields lattice
   */
  fields_lat->reset();

  // Get the potential parameters
  const double rhoB_0 = potentials.saturation_density();

  // update the fields lattice
  for (int i = 0; i < number_of_nodes; i++) {
    // read values off the jmu_B lattice (which holds values at t0 + Delta t)
    double rhoB_at_i = ((*jmuB_lat)[i]).rho();
    FourVector jmuB_at_i = ((*jmuB_lat)[i]).jmu_net();

    double abs_rhoB_at_i = std::abs(rhoB_at_i);
    // this is to prevent nan expressions
    if (abs_rhoB_at_i < very_small_double) {
      abs_rhoB_at_i = very_small_double;
    }

    // this needs to be used in order to prevent trying to calculate something
    // an expression like (-rhoB)^{3.4}
    const int sgn = rhoB_at_i > 0 ? 1 : -1;

    // field contributions as obtained in the VDF model
    double field_contribution = 0.0;
    for (int j = 0; j < potentials.number_of_terms(); j++) {
      field_contribution +=
          sgn * potentials.coeffs()[j] *
          std::pow(abs_rhoB_at_i / rhoB_0, potentials.powers()[j] - 2.0) /
          rhoB_0;
    }
    FourVector field_at_i = field_contribution * jmuB_at_i;

    // fill the A_mu lattice
    ((*fields_lat)[i]).overwrite_A_mu(field_at_i);
  }

  /*
   * Use the updated fields lattice, fields_lat, to populate the lattice of A^mu
   * FourVectors at t0 + Delta t, new_fields.
   */
  for (int i = 0; i < number_of_nodes; i++) {
    new_fields->assign_value(i, ((*fields_lat)[i]).A_mu());
  }

  /*
   * Compute time derivatives and gradients of all components of A^mu
   */
  new_fields->compute_four_gradient_lattice(*old_fields, time_step,
                                            *fields_four_grad_lattice);

  // substitute new derivatives
  for (int i = 0; i < number_of_nodes; i++) {
    auto tmp = (*fields_four_grad_lattice)[i];
    ((*fields_lat)[i]).overwrite_dAmu_dxnu(tmp[0], tmp[1], tmp[2], tmp[3]);
  }
}  // void update_fields_lattice()

}  // namespace smash
