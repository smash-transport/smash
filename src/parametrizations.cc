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

/** pp elastic cross section parametrization.
 * Source: \iref{Weil:2013mya}, eq. (44) */
float pp_elastic(double mandelstam_s) {
  double p_lab = plab_from_s_NN(mandelstam_s);
  if (p_lab < 0.435) {
    return 5.12 * nucleon_mass
        / (mandelstam_s - 4 * nucleon_mass * nucleon_mass) + 1.67;
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

/** pp total cross section parametrization.
 * Sources:
 * low-p: \iref{Cugnon:1996kh}
 * highest-p: \iref{Buss:2011mx}
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

/** np elastic cross section parametrization.
 * Source: \iref{Weil:2013mya}, eq. (45) */
float np_elastic(double mandelstam_s) {
  double p_lab = plab_from_s_NN(mandelstam_s);
  if (p_lab < 0.525) {
    return 17.05 * nucleon_mass
        / (mandelstam_s - 4 * nucleon_mass * nucleon_mass) - 6.83;
  } else if (p_lab < 0.8) {
    return 33 + 196 * std::pow(std::abs(p_lab - 0.95), 2.5);
  } else if (p_lab < 2.0) {
    return 31 / std::sqrt(p_lab);
  } else if (p_lab < 2.776) {
    return 77 / (p_lab + 1.5);
  } else {
    const auto logp = std::log(p_lab);
    return 11.9 + 26.9 * std::pow(p_lab, -1.21) + 0.169 * logp * logp
           - 1.85 * logp;
  }
}

/** np total cross section parametrization.
 * Sources:
 * low-p: \iref{Cugnon:1996kh}
 * highest-p: \iref{Buss:2011mx}
 */
float np_total(double mandelstam_s) {
  double p_lab = plab_from_s_NN(mandelstam_s);
  const auto logp = std::log(p_lab);
  if (p_lab < 0.4) {
    return 6.3555 * std::pow(p_lab, -3.2481) * std::exp(-0.377 * logp * logp);
  } else if (p_lab < 1.0) {
    return 33 + 196 * std::pow(std::abs(p_lab - 0.95), 2.5);
  } else if (p_lab < 2.0) {
    return 24.2 + 8.9 * p_lab;
  } else if (p_lab < 5.0) {
    return 42;
  } else {
    return 48.0 + 0.522 * logp * logp - 4.51 * logp;
  }
}

/** ppbar elastic cross section parametrization.
 * Source: \iref{Bass:1998ca} */
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

/** ppbar total cross section parametrization.
 * Source: \iref{Bass:1998ca} */
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

/** K+ p elastic cross section parametrization.
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8 */
float kplusp_elastic(double mandelstam_s) {
  constexpr double a0 = 10.508;  // mb
  constexpr double a1 = -3.716;  // mb/GeV
  constexpr double a2 = 1.845;  // mb/GeV^2
  constexpr double a3 = -0.764;  // GeV^-1
  constexpr double a4 = 0.508;  // GeV^-2

  const double p_lab = plab_from_s_NN(mandelstam_s);
  const double p_lab2 = p_lab*p_lab;

  return (a0 + a1*p_lab + a2*p_lab2) / (1 + a3*p_lab + a4*p_lab2);
}

/** K+ n elastic cross section parametrization.
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8 */
float kplusn_elastic(double mandelstam_s) {
  return 0.5 * kplusp_elastic(mandelstam_s);
}


/** K- p elastic cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
float kminusp_elastic(double mandelstam_s) {
  if (mandelstam_s < 1.7) {
    // The parametrization here also works for anti-K0 n, Lambda pi0,
    // Sigma+ pi-, Sigma- pi+, Sigma0 pi0 with different parameters a0, a1, a2.
    constexpr double a0 = 150;  // mb GeV^2
    constexpr double a1 = 0.35;  // Gev
    constexpr int a2 = 2;

    const double p_lab = plab_from_s_NN(mandelstam_s);
    const double p_i = p_lab;
    const double p_f = p_lab;

    const double ratio = a1*a1 / (a1*a1 + p_f*p_f);
    return a0 * p_f / (p_i * mandelstam_s) * std::pow(ratio, a2);
  } else {
    // TODO: Use spline interpolation of experimental data.
    throw std::runtime_error("not implemented");
  }
}

/** K- n elastic cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
float kminusn_elastic(double) {
    return 4.0;
}

/** K0 p elastic cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
float k0p_elastic(double mandelstam_s) {
    // by isospin symmetry
    return kplusn_elastic(mandelstam_s);
}

/** K0 n elastic cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
float k0n_elastic(double mandelstam_s) {
    // by isospin symmetry
    return kplusp_elastic(mandelstam_s);
}

/** Kbar0 p elastic cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
float kbar0p_elastic(double mandelstam_s) {
    // by isospin symmetry
    return kminusn_elastic(mandelstam_s);
}

/** Kbar0 n elastic cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
float kbar0n_elastic(double mandelstam_s) {
    // by isospin symmetry
    return kminusp_elastic(mandelstam_s);
}

}  // namespace Smash
