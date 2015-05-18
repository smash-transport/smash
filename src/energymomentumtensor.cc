/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/energymomentumtensor.h"
#include "include/numerics.h"
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

namespace Smash {

FourVector EnergyMomentumTensor::landau_frame_4velocity() const {
  using namespace Eigen;
  /* We want to solve the generalized eigenvalue problem
     T^{\mu \nu} h_{nu} = \lambda g^{\mu \nu} h_{nu}, or in the other way
     T_{\mu}^{\nu} h_{nu} = \lambda h_{mu}. Eigenvector
     corresponding to the largest (and the only positive) eigenvalue is
     proportional to 4-velocity of the Landau frame.
     Denote T^{\mu}_{\nu} as A and g^{\mu \nu} as B.
     A is symmetric and positive semi-definite. B is symmetric.
     A x = \lambda B x. I have to solve generalized eigenvalue
     problem, because A can be not positively defined (e.g. if
     energy-momentum tensor is computed for particles with momenta lying
     in one plane). For positively defined A a more efficient solution
     is possible, but I (oliiny) do not consider it until it becomes
     important for SMASH performance.
     */
  Matrix4d A;
  // A = T_{\mu}^{\nu} = g_{\mu \mu'} T^{\mu' \nu}
  A <<  Tmn_[0],  Tmn_[1],  Tmn_[2],  Tmn_[3],
       -Tmn_[1], -Tmn_[4], -Tmn_[5], -Tmn_[6],
       -Tmn_[2], -Tmn_[5], -Tmn_[7], -Tmn_[8],
       -Tmn_[3], -Tmn_[6], -Tmn_[8], -Tmn_[9];

  // log.debug("Looking for Landau frame for T_{mu}^{nu} ", A);
  EigenSolver<Matrix4d> es(A);

  // Eigen values should be strictly real and non-negative.

  // Here and further I assume that eigenvalues are given in
  // descending order. TODO(oliiny): check Eigen documentation
  // to make sure this is always true.
  Vector4d eig_im = es.eigenvalues().imag();
  Vector4d eig_re = es.eigenvalues().real();
  for (size_t i = 0; i < 4; i++){
    assert(std::abs(eig_im(i)) < really_small);
    if (i == 0) {
      assert(eig_re(i) > -really_small);
    } else {
      assert(eig_re(i) < really_small);
    }
  }

  auto tmp = es.eigenvectors().col(0).real();
  if (tmp(0) < 0.0) {
    tmp = -tmp;
  }

  FourVector u(tmp(0), tmp(1), tmp(2), tmp(3));
  const double u_sqr = u.sqr();
  if (u_sqr > really_small) {
    u /= std::sqrt(u_sqr);
  } else {
/*    log.error("Landau frame is not defined.",
              " Eigen vector", u, " of ", A, " is not time-like and",
              " cannot be 4-velocity. This may happen if energy-momentum",
              " tensor was constructed for a massless particle.");*/
    u = FourVector(1., 0., 0., 0.);
  }
//  std::cout << "Eigen values: " << eig_re << std::endl;
//  std::cout << "Eigen vectors: " << es.eigenvectors().real() << std::endl;
  return u;
}

std::ostream &operator<<(std::ostream &out, const EnergyMomentumTensor &Tmn) {
  using namespace std;
  out.width(12);
  for (size_t mu = 0; mu < 4; mu++) {
    for (size_t nu = 0; nu < 4; nu++) {
      out << setprecision(3) << setw(12) << fixed <<
             Tmn[EnergyMomentumTensor::tmn_index(mu,nu)];
    }
    out << endl;
  }
  return out;
}

}  // namespace Smash
