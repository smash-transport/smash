/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/quantumnumbers.h"

#include <sstream>

#include "smash/numerics.h"

namespace smash {

std::string QuantumNumbers::report_deviations(const QuantumNumbers& rhs) const {
  if (rhs == *this) {
    return "";
  }
  std::stringstream error_msg;
  error_msg << "Conservation law violations detected (old vs. new)\n";
  if (momentum_ != rhs.momentum_) {
    error_msg << "Deviation in Four-Momentum:\n" << std::scientific;
  }
  /* programmer's note: here, I'd like to simultaneously loop over an
   * integer (for the output; so that we know which component is
   * faulty) and both the current and rhs's momentum four-vector. If
   * there is a better way to do this, feel free to implement.
   *
   * I chose mu < 4 as the breaking condition out of the vague feeling
   * that comparing integers may be faster than accessing the
   * iterators. */
  int mu = 0;
  for (auto here_iter = momentum_.cbegin(), rhs_iter = rhs.momentum_.cbegin();
       mu < 4; ++here_iter, ++rhs_iter, ++mu) {
    if (!almost_equal_physics(*here_iter, *rhs_iter)) {
      error_msg << " P_" << mu << ": " << *here_iter << " vs. " << *rhs_iter
                << "; Δ = " << (*here_iter - *rhs_iter) << "\n";
    }
  }
  if (charge_ != rhs.charge_) {
    error_msg << "Deviation in Charge:\n " << charge_ << " vs. " << rhs.charge_
              << "\n";
  }
  if (isospin3_ != rhs.isospin3_) {
    error_msg << "Deviation in Isospin 3:\n " << isospin3_ << " vs. "
              << rhs.isospin3_ << "\n";
  }
  if (strangeness_ != rhs.strangeness_) {
    error_msg << "Deviation in Strangeness:\n " << strangeness_ << " vs. "
              << rhs.strangeness_ << "\n";
  }
  if (charmness_ != rhs.charmness_) {
    error_msg << "Deviation in Charmness:\n " << charmness_ << " vs. "
              << rhs.charmness_ << "\n";
  }
  if (bottomness_ != rhs.bottomness_) {
    error_msg << "Deviation in Bottomness:\n " << bottomness_ << " vs. "
              << rhs.bottomness_ << "\n";
  }
  if (baryon_number_ != rhs.baryon_number_) {
    error_msg << "Deviation in Baryon Number:\n " << baryon_number_ << " vs. "
              << rhs.baryon_number_ << "\n";
  }
  return error_msg.str();
}

}  // namespace smash
