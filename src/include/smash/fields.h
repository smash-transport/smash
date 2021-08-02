/*
 *
 *    Copyright (c) 2013-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_SMASH_FIELDS_H_
#define SRC_INCLUDE_SMASH_FIELDS_H_

#include <iostream>
#include <typeinfo>
#include <utility>
#include <vector>

#include "density.h"
#include "experimentparameters.h"
//#include "forwarddeclarations.h"
#include "fourvector.h"
//#include "lattice.h"
#include "potentials.h"
#include "threevector.h"

namespace smash {

/**
 * A class for calculating the fields A^mu associated with the VDF potentials.
 * The structure of the class is heavily based on the DensityOnLattice class.
 * 
 * It holds the values of the A^mu FourVector as well as fourgradients of its
 * components.
 */
class FieldsOnLattice {
 public:
  /// Default constructor
  FieldsOnLattice()
      : A_mu_(FourVector()),
        dAmu_dxnu_({FourVector(), FourVector(), FourVector(), FourVector()}) {}


  /**
   * \return The field four-vector A^mu
   */
  FourVector A_mu() const { return A_mu_; }

  /**
   * \return The four-gradient of the field four-vector A^mu
   */
  std::array<FourVector, 4> dAmu_dxnu() const { return dAmu_dxnu_; }

  /**
   * \return The time derivative of \vec{A} on the local lattice
   */
  ThreeVector dvecA_dt() { return dAmu_dxnu_[0].threevec(); }

  /**
   * Compute the gradient of A^0 on the local lattice
   *
   * \return \f$\nabla A^0\f$
   */
  ThreeVector grad_A_0() {
    ThreeVector A_0_grad = ThreeVector();
    for (int i = 1; i < 4; i++) {
      A_0_grad[i - 1] = dAmu_dxnu_[i].x0();
    }
    return A_0_grad;
  }

  /**
   * Compute the curl of the field on the local lattice
   *
   * \return \f$\nabla\times\vec{A}\f$
   */
  ThreeVector curl_vec_A() {
    ThreeVector curl_vecA = ThreeVector();
    curl_vecA.set_x1(dAmu_dxnu_[2].x3() - dAmu_dxnu_[3].x2());
    curl_vecA.set_x2(dAmu_dxnu_[3].x1() - dAmu_dxnu_[1].x3());
    curl_vecA.set_x3(dAmu_dxnu_[1].x2() - dAmu_dxnu_[2].x1());
    return curl_vecA;
  }

  /**
   * Overwrite the value of the field on the local lattice
   *
   * \param[in] new_A_mu new value of the field
   */
  void overwrite_A_mu (FourVector new_A_mu)
  {
    A_mu_ = new_A_mu;
  }

  
  /**
   * Overwrite the time derivative of A^mu to zero
   */
  void overwrite_dAmu_dt_to_zero ()
  {
    FourVector tmp0 (0.0, 0.0, 0.0, 0.0);

    dAmu_dxnu_[0] = tmp0;
  }

  /**
   * Overwrite the four-gradient of A^mu on the local lattice, using the
   * provided values of its components.
   *
   * \param[in] dAmu_dt new value of the time derivative of A^mu
   * \param[in] dAmu_dx new value of the x-derivative of A^mu
   * \param[in] dAmu_dy new value of the y-derivative of A^mu
   * \param[in] dAmu_dz new value of the z-derivative of A^mu
   */
  void overwrite_dAmu_dxnu (FourVector dAmu_dt,FourVector dAmu_dx,
			    FourVector dAmu_dy, FourVector dAmu_dz) {
    dAmu_dxnu_[0] = dAmu_dt;
    dAmu_dxnu_[1] = dAmu_dx;
    dAmu_dxnu_[2] = dAmu_dy;
    dAmu_dxnu_[3] = dAmu_dz;
  }

  // TO DO: consider calculating the mean-field energy here
  //double calculate_mean_field_from_fields ();


 private:
  /// Four-vector density of the field
  FourVector A_mu_;
  /// Four-gradient of the four-vector density of the field
  std::array<FourVector, 4> dAmu_dxnu_;
};

/// Conveniency typedef for lattice of fields
typedef RectangularLattice<FieldsOnLattice> FieldsLattice;



/**
 * Updates the contents on the lattice of <FieldsOnLattice> type.
 *
 * \param[out] fields_lat The lattice of FieldsOnLattice type on which the
 *             content will be updated
 * \param[in] old_fields Auxiliary lattice, filled with field values at t0,
 *            needed for calculating time derivatives
 * \param[in] new_fields Auxiliary lattice, filled with field values at t0 + dt,
 *            needed for calculating time derivatives
 * \param[in] fields_four_grad_lattice Auxiliary lattice for calculating the
 *            fourgradient of the fields
 * \param[in] fields_lat_update Tells if called for update at printout or at
 *            timestep
 * \param[in] time_step Time step used in the simulation
 */
void update_lattice(
    RectangularLattice<FieldsOnLattice> *fields_lat,
    RectangularLattice<FourVector> *old_fields,
    RectangularLattice<FourVector> *new_fields,
    RectangularLattice<std::array<FourVector, 4>> *fields_four_grad_lattice,
    const LatticeUpdate fields_lat_update, const Potentials &potentials,
    const double time_step);
}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_DENSITY_H_
