/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/parametrizations.h"

#include <cmath>

#include "include/kinematics.h"

namespace Smash {

/* pp elastic cross section parametrization.
 * Source: J. Weil, PhD thesis, eq. (44) */
float pp_elastic(double mandelstam_s) {
  double p_lab = plab_from_s_NN(mandelstam_s);
  if (p_lab < 0.435) {
    return 5.12 * mN / (mandelstam_s - 4 * mN * mN) + 1.67;
  } else if (p_lab < 0.8) {
    return 23.5 + 1000 * (p_lab - 0.7) * (p_lab - 0.7)
      * (p_lab - 0.7) *(p_lab - 0.7);
  } else if (p_lab < 2.0) {
    return 1250 / (p_lab + 50) - 4 * (p_lab - 1.3) * (p_lab - 1.3);
  } else if (p_lab < 2.776) {
    return 77 / (p_lab + 1.5);
  } else {
    const auto logp = std::log(p_lab);
    return 11.9 + 26.9 * std::pow(p_lab, -1.21) + 0.169 * logp * logp
           - 1.85 * logp;
  }
}

/* pp total cross section parametrization */
/* Sources:
 * low-p: J. Cugnon, D. L'Hote, J. Vandermeulen,
 * Nuclear Instruments and Methods at Physics Research B 111, 215 (1996)
 * highest-p:  O. Buss et al., Physics Reports 512, 1 (2012)
 */
float pp_total(double mandelstam_s) {
  double p_lab = plab_from_s_NN(mandelstam_s);
  if (p_lab < 0.4) {
    return 34 * std::pow(p_lab / 0.4, -2.104);
  } else if (p_lab < 0.8) {
    return 23.5 + 1000 * (p_lab - 0.7) * (p_lab - 0.7)
      * (p_lab - 0.7) *(p_lab - 0.7);
  } else if (p_lab < 1.5) {
    return 23.5 + 24.6 / (1 + std::exp(-(p_lab - 1.2) / 0.1));
  } else if (p_lab < 5.0) {
    return 41 + 60 * (p_lab - 0.9) * std::exp(-1.2 * p_lab);
  } else {
    const auto logp = std::log(p_lab);
    return 48.0 + 0.522 * logp * logp - 4.51 * logp;
  }
}

/* np elastic cross section parametrization.
 * Source: J. Weil, PhD thesis, eq. (45) */
float np_elastic(double mandelstam_s) {
  double p_lab = plab_from_s_NN(mandelstam_s);
  if (p_lab < 0.525) {
    return 17.05 * mN / (mandelstam_s - 4 * mN * mN) - 6.83;
  } else if (p_lab < 0.8) {
    return 33 + 196 * std::pow(fabs(p_lab - 0.95), 2.5);
  } else if (p_lab < 2.0) {
    return 31 / sqrt(p_lab);
  } else if (p_lab < 2.776) {
    return 77 / (p_lab + 1.5);
  } else {
    const auto logp = std::log(p_lab);
    return 11.9 + 26.9 * std::pow(p_lab, -1.21) + 0.169 * logp * logp
           - 1.85 * logp;
  }
}

/* np total cross section parametrization */
/* Sources:
 * low-p: J. Cugnon, D. L'Hote, J. Vandermeulen,
 * Nuclear Instruments and Methods at Physics Research B 111, 215 (1996)
 * highest-p:  O. Buss et al., Physics Reports 512, 1 (2012)
 */
float np_total(double mandelstam_s) {
  double p_lab = plab_from_s_NN(mandelstam_s);
  const auto logp = std::log(p_lab);
  if (p_lab < 0.4) {
    return 6.3555 * std::pow(p_lab, -3.2481) * std::exp(-0.377 * logp * logp);
  } else if (p_lab < 1.0) {
    return 33 + 196 * std::pow(fabs(p_lab - 0.95), 2.5);
  } else if (p_lab < 2.0) {
    return 24.2 + 8.9 * p_lab;
  } else if (p_lab < 5.0) {
    return 42;
  } else {
    return 48.0 + 0.522 * logp * logp - 4.51 * logp;
  }
}

/* ppbar elastic cross section parametrization */
/* Source: S. Bass et al., Prog.Part.Nucl.Phys. 41, 255 (1998) */
float ppbar_elastic(double mandelstam_s) {
  double p_lab = plab_from_s_NN(mandelstam_s);
  if (p_lab < 0.3) {
    return 78.6;
  } else if (p_lab < 5.0) {
    return 31.6 + 18.3 / p_lab - 1.1 / (p_lab * p_lab) - 3.8 * p_lab;
  } else {
    const auto logp = std::log(p_lab);
    return 10.2 + 52.7 * std::pow(p_lab, -1.16) + 0.125 * logp * logp
           - 1.28 * logp;
  }
}

/* ppbar total cross section parametrization */
/* Source: S. Bass et al., Prog.Part.Nucl.Phys. 41, 255 (1998) */
float ppbar_total(double mandelstam_s) {
  double p_lab = plab_from_s_NN(mandelstam_s);
  if (p_lab < 0.3) {
    return 271.6 * std::exp(-1.1 * p_lab * p_lab);
  } else if (p_lab < 5.0) {
    return 75.0 + 43.1 / p_lab + 2.6 / (p_lab * p_lab) - 3.9 * p_lab;
  } else {
    const auto logp = std::log(p_lab);
    return 38.4 + 77.6 * std::pow(p_lab, -0.64) + 0.26 * logp * logp
           - 1.2 * logp;
  }
}

}  // namespace Smash
