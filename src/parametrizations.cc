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
#include <iostream>
#include <memory>
#include <vector>

#include "include/average.h"
#include "include/cxx14compat.h"
#include "include/kinematics.h"
#include "include/interpolation.h"
#include "include/lowess.h"

namespace Smash {

/** pp elastic cross section parametrization.
 * Source: \iref{Weil:2013mya}, eq. (44) */
float pp_elastic(double mandelstam_s) {
  double p_lab = plab_from_s(mandelstam_s);
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
  double p_lab = plab_from_s(mandelstam_s);
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
  double p_lab = plab_from_s(mandelstam_s);
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
  double p_lab = plab_from_s(mandelstam_s);
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
  double p_lab = plab_from_s(mandelstam_s);
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
  double p_lab = plab_from_s(mandelstam_s);
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

  const double p_lab = plab_from_s(mandelstam_s, kaon_mass, nucleon_mass);
  const double p_lab2 = p_lab*p_lab;

  return (a0 + a1*p_lab + a2*p_lab2) / (1 + a3*p_lab + a4*p_lab2);
}

/** K+ n elastic cross section parametrization.
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8 */
float kplusn_elastic(double mandelstam_s) {
  return 0.5 * kplusp_elastic(mandelstam_s);
}

/// PDG data on K- p elastic cross section: momentum in lab frame.
const std::initializer_list<double> KMINUSP_ELASTIC_P_LAB = {
  0.03000,   0.05000,   0.06287,   0.07000,  0.07044,  0.07300,  0.08730,
  0.09000,   0.11000,   0.11000,   0.11210,  0.11300,  0.12262,  0.13000,
  0.13000,   0.13722,   0.14552,   0.15000,  0.15000,  0.15300,  0.15909,
  0.16269,   0.17000,   0.17000,   0.17470,  0.18768,  0.18916,  0.19000,
  0.19000,   0.19300,   0.20534,   0.21000,  0.21000,  0.21554,  0.22049,
  0.22500,   0.23000,   0.23000,   0.23300,  0.23500,  0.23944,  0.24500,
  0.24623,   0.25000,   0.25500,   0.26154,  0.26500,  0.27000,  0.27500,
  0.27618,   0.28227,   0.28500,   0.29000,  0.29290,  0.29300,  0.29500,
  0.30000,   0.30191,   0.30500,   0.31500,  0.32500,  0.33500,  0.34500,
  0.35000,   0.35012,   0.35500,   0.36500,  0.37500,  0.38500,  0.38700,
  0.38986,   0.39200,   0.39500,   0.40000,  0.40500,  0.41500,  0.42000,
  0.42500,   0.43400,   0.43407,   0.43500,  0.43600,  0.44500,  0.45500,
  0.45500,   0.46500,   0.47500,   0.49500,  0.51300,  0.51300,  0.51400,
  0.52228,   0.53400,   0.55400,   0.57300,  0.59700,  0.61000,  0.61700,
  0.62004,   0.63700,   0.64200,   0.65800,  0.66000,  0.67200,  0.67700,
  0.69900,   0.70000,   0.70300,   0.70800,  0.71900,  0.72500,  0.73000,
  0.74000,   0.74000,   0.74100,   0.75800,  0.75996,  0.76100,  0.76800,
  0.77300,   0.77700,   0.77700,   0.78000,  0.78500,  0.79300,  0.80200,
  0.80600,   0.80600,   0.81000,   0.82000,  0.82000,  0.83300,  0.83800,
  0.83800,   0.84999,   0.85300,   0.85300,  0.85600,  0.86000,  0.87000,
  0.87400,   0.87400,   0.87600,   0.89400,  0.89400,  0.89900,  0.90000,
  0.90400,   0.90400,   0.90500,   0.91600,  0.91600,  0.92200,  0.92500,
  0.93500,   0.93500,   0.93500,   0.94000,  0.94300,  0.94500,  0.95400,
  0.95400,   0.95500,   0.96000,   0.96500,  0.97000,  0.97000,  0.98000,
  0.99000,   0.99100,   0.99100,   1.00500,  1.02000,  1.02200,  1.02200,
  1.04000,   1.04400,   1.04400,   1.04500,  1.06000,  1.06100,  1.06100,
  1.08000,   1.08000,   1.08500,   1.10000,  1.10200,  1.10200,  1.11700,
  1.11700,   1.12500,   1.13400,   1.13400,  1.13800,  1.14000,  1.15000,
  1.15000,   1.15300,   1.15300,   1.16100,  1.16500,  1.17400,  1.17400,
  1.17900,   1.18000,   1.18300,   1.18300,  1.20100,  1.20500,  1.22000,
  1.22600,   1.22600,   1.23300,   1.24500,  1.25300,  1.26000,  1.26000,
  1.26300,   1.26300,   1.27600,   1.28500,  1.29600,  1.30000,  1.31600,
  1.31600,   1.32000,   1.32800,   1.34000,  1.35500,  1.36800,  1.36800,
  1.38000,   1.38300,   1.41500,   1.41500,  1.42300,  1.43300,  1.46200,
  1.46500,   1.48300,   1.51300,   1.51400,  1.53000,  1.53400,  1.54500,
  1.54600,   1.58400,   1.60600,   1.60600,  1.63400,  1.65200,  1.65300,
  1.68000,   1.68400,   1.70500,   1.70500,  1.73400,  1.73900,  1.74100,
  1.78400,   1.80000,   1.80000,   1.81500,  1.84300,  1.84300,  1.88400,
  1.93400,   1.93400,   1.98400,   2.00000,  2.03100,  2.03400,  2.08400,
  2.13400,   2.13500,   2.17500,   2.23400,  2.24000,  2.28400,  2.32500,
  2.33100,   2.37400,   2.41200,   2.51600,  2.66000,  2.66000,  2.98500,
  3.00000,   3.00000,   3.46000,   3.59000,  3.65000,  4.10000,  4.20000,
  4.60000,   5.00000,   5.50000,   6.00000,  7.20020,  9.00010,  10.12000,
  14.30000,  14.30000,  25.20000,  32.10000, 40.10000, 50.00000, 70.00000,
  100.00000, 140.00000, 147.00000, 175.00000
};
/// PDG data on K- p elastic cross section: cross section.
const std::initializer_list<double> KMINUSP_ELASTIC_SIG = {
  313.50, 103.60, 113.00, 44.800, 58.500, 187.00, 92.000, 71.500, 92.800,
  87.290, 82.000, 105.00, 59.400, 40.400, 79.220, 82.000, 49.000, 41.400,
  69.610, 108.00, 53.900, 98.000, 75.760, 32.800, 45.000, 73.000, 66.000,
  59.090, 53.300, 68.000, 37.000, 53.300, 60.490, 48.000, 41.000, 65.410,
  62.900, 55.690, 50.000, 55.740, 41.200, 53.240, 37.000, 51.500, 49.220,
  43.600, 47.710, 58.060, 48.750, 30.000, 44.900, 39.420, 38.270, 47.800,
  48.200, 41.220, 44.500, 44.200, 40.360, 37.020, 40.280, 37.840, 37.260,
  34.000, 33.500, 34.770, 34.210, 36.670, 33.890, 31.900, 34.700, 34.000,
  38.470, 38.900, 32.060, 32.590, 48.400, 31.190, 30.600, 32.800, 26.830,
  25.800, 28.830, 23.800, 30.320, 31.990, 23.100, 21.500, 26.500, 27.600,
  21.700, 35.000, 19.300, 19.100, 17.500, 17.700, 17.660, 18.600, 16.000,
  16.000, 17.230, 16.400, 12.100, 16.220, 15.600, 15.200, 14.200, 15.220,
  13.500, 14.200, 11.500, 14.070, 15.900, 14.000, 11.600, 16.930, 16.700,
  17.300, 15.200, 18.600, 18.300, 18.700, 17.900, 19.370, 20.500, 19.300,
  20.000, 20.450, 20.670, 18.700, 19.300, 19.630, 19.800, 19.530, 22.400,
  19.100, 19.390, 19.780, 19.500, 21.360, 20.100, 20.310, 21.070, 21.600,
  21.660, 22.020, 21.500, 20.900, 20.710, 21.340, 20.850, 20.400, 22.340,
  20.940, 20.120, 20.100, 19.980, 21.800, 21.010, 19.330, 20.640, 20.700,
  20.610, 21.270, 20.340, 20.400, 20.720, 22.400, 21.220, 21.400, 22.150,
  21.750, 20.800, 22.100, 22.150, 23.300, 21.500, 21.460, 22.220, 21.200,
  20.600, 20.560, 18.700, 18.740, 19.830, 18.300, 18.600, 18.770, 17.600,
  17.820, 17.890, 17.000, 17.750, 15.700, 17.000, 18.300, 17.300, 17.200,
  17.230, 15.300, 15.390, 16.500, 16.460, 14.660, 16.800, 15.900, 15.890,
  12.680, 13.890, 15.700, 11.300, 11.810, 13.110, 12.320, 11.870, 15.200,
  14.000, 10.870, 10.870, 11.440, 10.960, 11.570, 12.000, 10.200, 11.200,
  10.260, 9.7400, 14.400, 9.5300, 10.300, 10.300, 16.600, 10.500, 8.8300,
  8.8300, 8.4200, 8.6000, 9.1100, 9.1100, 8.3000, 7.7000, 7.7000, 8.1500,
  8.8000, 8.0600, 8.0600, 8.7000, 8.8600, 8.8600, 9.2000, 8.4000, 8.4000,
  8.1900, 8.9000, 9.0800, 9.0800, 9.5000, 8.5100, 8.5100, 9.0000, 8.1300,
  8.1300, 6.9500, 7.8600, 7.8600, 9.0000, 7.9000, 7.4800, 7.8000, 7.4600,
  7.3100, 7.1000, 7.4000, 7.1000, 7.2600, 7.9000, 6.8700, 6.2000, 7.3000,
  6.5000, 6.7600, 6.6000, 6.3200, 6.2100, 5.7000, 6.0400, 4.9500, 5.0600,
  4.9500, 4.9400, 4.4000, 4.6000, 4.3000, 4.5000, 4.2000, 3.8400, 4.1000,
  3.6200, 4.2300, 3.9500, 3.2400, 2.9600, 3.0100, 2.4600, 2.5600, 2.3300,
  2.5400, 2.5300, 2.5100, 2.5200, 2.7400, 2.5900
};
static std::unique_ptr<InterpolateDataLinear<double>>
    kminusp_elastic_interpolation = nullptr;

/// Center-of-mass energy.
const std::initializer_list<double> KMINUSP_RES_SQRTS = {
  1.4325, 1.45, 1.47, 1.49, 1.51, 1.53, 1.55, 1.57, 1.59, 1.61, 1.63, 1.65,
  1.67,   1.69, 1.71, 1.73, 1.75, 1.77, 1.79, 1.81, 1.83, 1.85, 1.87, 1.89,
  1.91,   1.93, 1.95, 1.97, 1.99, 2.01, 2.03, 2.05, 2.07, 2.09, 2.11, 2.13,
  2.15,   2.17, 2.19, 2.21, 2.23, 2.25, 2.27, 2.29, 2.31, 2.33, 2.35, 2.37,
  2.39,   2.41, 2.43, 2.45, 2.47, 2.49, 2.51, 2.53,
};
/// Elastic K- p cross section contributions from decays.
///
/// These need to be subtracted from the interpolation of the PDG data on
/// elastic cross sections. This data was generated using the SMASH analysis
/// suite and should be updated when strange resonances are changed or added.
const std::initializer_list<double> KMINUSP_RES_SIG = {
    0.181508888856, 0.276824217055, 0.52709680347,   0.809170022877,
    5.4306649414,   2.73623096896,  1.00335277069,   0.831092759384,
    0.96487816871,  1.11429990074,  1.63535440366,   2.43199415481,
    3.35420549211,  3.5880964655,   3.55175458501,   4.17980540217,
    5.69214168171,  6.8367923342,   9.03411977658,   10.7000792054,
    10.4654818833,  9.42857422653,  10.2164594715,   10.7163541445,
    8.93164259792,  7.09269895697,  5.17386434884,   3.75480313343,
    2.62782096271,  2.07496543859,  1.64469989838,   1.17591535168,
    0.966825136065, 0.885320244975, 0.758643195654,  0.637513694965,
    0.467191033361, 0.462018802703, 0.371174185596,  0.366647793383,
    0.297093948546, 0.267488464682, 0.256312096717,  0.208579409826,
    0.190407843837, 0.175917906512, 0.178264019383,  0.149864315214,
    0.135615565435, 0.13676468557,  0.113707047901,  0.127858971204,
    0.111330153509, 0.105314540758, 0.0897379082466, 0.0810976351292,
};
static std::unique_ptr<InterpolateDataSpline> kminusp_elastic_res_interpolation
    = nullptr;

/** K- p elastic cross section parametrization, PDG data.
 *
 * The PDG data is being interpolated using cubic splines. If more than one
 * cross section was given for one p_lab value, the corresponding cross sections
 * are averaged.
 */
static float kminusp_elastic_pdg(double mandelstam_s) {
  if (kminusp_elastic_interpolation == nullptr) {
    std::vector<double> x = KMINUSP_ELASTIC_P_LAB;
    std::vector<double> y = KMINUSP_ELASTIC_SIG;
    std::vector<double> dedup_x;
    std::vector<double> dedup_y;
    std::tie(dedup_x, dedup_y) = dedup_avg(x, y);
    dedup_y = smooth(dedup_x, dedup_y, 0.1, 5);
    kminusp_elastic_interpolation =
        make_unique<InterpolateDataLinear<double>>(dedup_x, dedup_y);
    /*
    // Output interpolation for plotting.
    constexpr int N = 10000;
    const double x_min = 0.1;
    const double x_max = 100;
    std::cout << "\n-------------------\n";
    for (int i = 0; i < N; i++) {
        const double xi = x_min + (x_max - x_min) * (i /
    static_cast<double>(N));
        std::cout << xi << " " << (*kminusp_elastic_interpolation)(xi) << "\n";
    }
    std::cout << "-------------------" << std::endl;
    */
  }
  const double p_lab = plab_from_s(mandelstam_s, kaon_mass, nucleon_mass);
  return (*kminusp_elastic_interpolation)(p_lab);
}

/** K- p elastic cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
float kminusp_elastic(double mandelstam_s) {
  const double p_lab = plab_from_s(mandelstam_s, kaon_mass, nucleon_mass);
  double sigma;
  if (std::sqrt(mandelstam_s) < 1.68) {
    // The parametrization here also works for anti-K0 n, Lambda pi0,
    // Sigma+ pi-, Sigma- pi+, Sigma0 pi0 with different parameters a0, a1, a2.
    //
    // The values of the parameters are *not* taken from the source above,
    // they come from a fit to PDG data.
    constexpr double a0 = 186.03567644;  // mb GeV^2
    constexpr double a1 = 0.22002795;  // Gev
    constexpr double a2 = 0.64907116;

    const double p_i = p_lab;
    const double p_f = p_lab;

    const double ratio = a1*a1 / (a1*a1 + p_f*p_f);
    sigma = a0 * p_f / (p_i * mandelstam_s) * std::pow(ratio, a2);
  } else {
    sigma = kminusp_elastic_pdg(mandelstam_s);
  }
  // The elastic contributions from decays still need to be subtracted.
  if (kminusp_elastic_res_interpolation == nullptr) {
    std::vector<double> x = KMINUSP_RES_SQRTS;
    for (auto& i : x) {
      i = plab_from_s(i * i, kaon_mass, nucleon_mass);
    }
    std::vector<double> y = KMINUSP_RES_SIG;
    kminusp_elastic_res_interpolation =
        make_unique<InterpolateDataSpline>(x, y);
  }
  const auto old_sigma = sigma;
  sigma -= (*kminusp_elastic_res_interpolation)(p_lab);
  if (sigma < 0) {
    std::cout << "NEGATIVE SIGMA: sigma=" << sigma
              << ", sqrt(s)=" << std::sqrt(mandelstam_s)
              << ", sig_el_exp=" << old_sigma
              << ", sig_el_res=" << (*kminusp_elastic_res_interpolation)(p_lab)
              << std::endl;
  }
  assert(sigma >= 0);
  return sigma;
}

/** K- n elastic cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
float kminusn_elastic(double) { return 4.0; }

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
