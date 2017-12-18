/*
 *
 *    Copyright (c) 2013-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/parametrizations.h"
#include "include/parametrizations_data.h"

#include <cmath>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <vector>

#include "include/average.h"
#include "include/clebschgordan.h"
#include "include/cxx14compat.h"
#include "include/interpolation.h"
#include "include/kinematics.h"
#include "include/lowess.h"
#include "include/pow.h"

// All quantities in this file use they same units as the rest of SMASH.
// That is: GeV for energies and momenta, fm for distances and time, and mb for
// cross sections.

namespace smash {

/** total hadronic cross sections at high energies parametrized in the 2016 PDG
 *  book(http://pdg.lbl.gov/2016/reviews/rpp2016-rev-cross-section-plots.pdf) */
double xs_high_energy(double mandelstam_s, bool is_opposite_charge, double ma,
                      double mb, double P, double R1, double R2) {
  const double M = 2.1206;
  const double H = 0.272;
  const double eta1 = 0.4473;
  const double eta2 = 0.5486;
  const double s_sab = mandelstam_s / (ma + mb + M) / (ma + mb + M);
  double xs =
      H * std::log(s_sab) * std::log(s_sab) + P + R1 * std::pow(s_sab, -eta1);
  xs = is_opposite_charge ? xs + R2 * std::pow(s_sab, -eta2)
                          : xs - R2 * std::pow(s_sab, -eta2);
  return xs;
}

double pp_high_energy(double mandelstam_s) {
  return xs_high_energy(mandelstam_s, false, 0.939, 0.939, 34.41, 13.07, 7.394);
}

double ppbar_high_energy(double mandelstam_s) {
  return xs_high_energy(mandelstam_s, true, 0.939, 0.939, 34.41, 13.07, 7.394);
}

double np_high_energy(double mandelstam_s) {
  return xs_high_energy(mandelstam_s, false, 0.939, 0.939, 34.41, 12.52, 6.66);
}

double npbar_high_energy(double mandelstam_s) {
  return xs_high_energy(mandelstam_s, true, 0.939, 0.939, 34.41, 12.52, 6.66);
}

double piplusp_high_energy(double mandelstam_s) {
  return xs_high_energy(mandelstam_s, false, 0.939, 0.138, 18.75, 9.56, 1.767);
}

double piminusp_high_energy(double mandelstam_s) {
  return xs_high_energy(mandelstam_s, true, 0.939, 0.138, 18.75, 9.56, 1.767);
}

double xs_string_hard(double mandelstam_s,
                      double xs_0, double e_0, double lambda_pow) {
  const double sqrts = std::sqrt(mandelstam_s);
  if(sqrts < e_0) {
    return 0.;
  } else {
    const double xs = xs_0 * std::pow(std::log(sqrts / e_0), lambda_pow);
    return xs;
  }
}

double NN_string_hard(double mandelstam_s) {
  return xs_string_hard(mandelstam_s, 0.087, 4.1, 3.8);
}

double Npi_string_hard(double mandelstam_s) {
  return xs_string_hard(mandelstam_s, 0.042, 3.5, 4.2);
}

double pipi_string_hard(double mandelstam_s) {
  return xs_string_hard(mandelstam_s, 0.013, 2.3, 4.7);
}

/** pi+p elastic cross section parametrization.
 * Source: GiBUU:parametrizationBarMes_HighEnergy.f90
 * The parametrizations of the elastic pion+nucleon cross sections
 * are still under tuning. The parametrizaton is employed to give a
 * non-zero cross section at high energies. To make sure it
 * doesn't affect the cross section at the low energies, I truncate
 * the parametrization at p_lab = 8 GeV, which correspons to square
 * root of s equal to 4 GeV. */
double piplusp_elastic(double mandelstam_s) {
  double p_lab = plab_from_s(mandelstam_s, pion_mass, nucleon_mass);
  if (p_lab < 8.0 /*1.45*/) {
    return really_small;
    //} else if (p_lab < 2.0) {
    //    return 16.0 - (p_lab - 1.45) * 7.5 / 0.55;
  } else {
    const auto logp = std::log(p_lab);
    return 11.4 * std::pow(p_lab, -0.4) + 0.079 * logp * logp;
  }
}

/** pi-p elastic cross section parametrization.
 * Source: GiBUU:parametrizationBarMes_HighEnergy.f90 */
double piminusp_elastic(double mandelstam_s) {
  double p_lab = plab_from_s(mandelstam_s, pion_mass, nucleon_mass);
  if (p_lab < 8.0 /*2.0*/) {
    return really_small;
  } else {
    const auto logp = std::log(p_lab);
    return 1.76 + 11.2 * std::pow(p_lab, -0.64) + 0.043 * logp * logp;
  }
}

/** pp elastic cross section parametrization.
 * Source: \iref{Weil:2013mya}, eq. (44) */
double pp_elastic(double mandelstam_s) {
  double p_lab = plab_from_s(mandelstam_s);
  if (p_lab < 0.435) {
    return 5.12 * nucleon_mass /
               (mandelstam_s - 4 * nucleon_mass * nucleon_mass) +
           1.67;
  } else if (p_lab < 0.8) {
    return 23.5 +
           1000 * (p_lab - 0.7) * (p_lab - 0.7) * (p_lab - 0.7) * (p_lab - 0.7);
  } else if (p_lab < 2.0) {
    return 1250 / (p_lab + 50) - 4 * (p_lab - 1.3) * (p_lab - 1.3);
  } else if (p_lab < 2.776) {
    return 77 / (p_lab + 1.5);
  } else {
    const auto logp = std::log(p_lab);
    return 11.9 + 26.9 * std::pow(p_lab, -1.21) + 0.169 * logp * logp -
           1.85 * logp;
  }
}

/** pp total cross section parametrization.
 * Sources:
 * low-p: \iref{Cugnon:1996kh}
 * highest-p: \iref{Buss:2011mx}
 */
double pp_total(double mandelstam_s) {
  double p_lab = plab_from_s(mandelstam_s);
  if (p_lab < 0.4) {
    return 34 * std::pow(p_lab / 0.4, -2.104);
  } else if (p_lab < 0.8) {
    return 23.5 +
           1000 * (p_lab - 0.7) * (p_lab - 0.7) * (p_lab - 0.7) * (p_lab - 0.7);
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
double np_elastic(double mandelstam_s) {
  double p_lab = plab_from_s(mandelstam_s);
  if (p_lab < 0.525) {
    return 17.05 * nucleon_mass /
               (mandelstam_s - 4 * nucleon_mass * nucleon_mass) -
           6.83;
  } else if (p_lab < 0.8) {
    return 33 + 196 * std::pow(std::abs(p_lab - 0.95), 2.5);
  } else if (p_lab < 2.0) {
    return 31 / std::sqrt(p_lab);
  } else if (p_lab < 2.776) {
    return 77 / (p_lab + 1.5);
  } else {
    const auto logp = std::log(p_lab);
    return 11.9 + 26.9 * std::pow(p_lab, -1.21) + 0.169 * logp * logp -
           1.85 * logp;
  }
}

/** np total cross section parametrization.
 * Sources:
 * low-p: \iref{Cugnon:1996kh}
 * highest-p: \iref{Buss:2011mx}
 */
double np_total(double mandelstam_s) {
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
double ppbar_elastic(double mandelstam_s) {
  double p_lab = plab_from_s(mandelstam_s);
  if (p_lab < 0.3) {
    return 78.6;
  } else if (p_lab < 5.0) {
    return 31.6 + 18.3 / p_lab - 1.1 / (p_lab * p_lab) - 3.8 * p_lab;
  } else {
    const auto logp = std::log(p_lab);
    return 10.2 + 52.7 * std::pow(p_lab, -1.16) + 0.125 * logp * logp -
           1.28 * logp;
  }
}

/** ppbar total cross section parametrization.
 * Source: \iref{Bass:1998ca} */
double ppbar_total(double mandelstam_s) {
  double p_lab = plab_from_s(mandelstam_s);
  if (p_lab < 0.3) {
    return 271.6 * std::exp(-1.1 * p_lab * p_lab);
  } else if (p_lab < 5.0) {
    return 75.0 + 43.1 / p_lab + 2.6 / (p_lab * p_lab) - 3.9 * p_lab;
  } else {
    const auto logp = std::log(p_lab);
    return 38.4 + 77.6 * std::pow(p_lab, -0.64) + 0.26 * logp * logp -
           1.2 * logp;
  }
}

/** K+ p elastic background cross section parametrization.
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8 */
double kplusp_elastic_background(double mandelstam_s) {
  constexpr double a0 = 10.508;  // mb
  constexpr double a1 = -3.716;  // mb/GeV
  constexpr double a2 = 1.845;   // mb/GeV^2
  constexpr double a3 = -0.764;  // GeV^-1
  constexpr double a4 = 0.508;   // GeV^-2

  const double p_lab = plab_from_s(mandelstam_s, kaon_mass, nucleon_mass);
  const double p_lab2 = p_lab * p_lab;

  return (a0 + a1 * p_lab + a2 * p_lab2) / (1 + a3 * p_lab + a4 * p_lab2);
}

/** K+ n elastic background cross section parametrization.
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8 */
double kplusn_elastic_background(double mandelstam_s) {
  return 0.5 * kplusp_elastic_background(mandelstam_s);
}

/** K+ n charge exchange cross section parametrization.
 * sigma(K+n->K+n) = sigma(K+n->K0p) = 0.5 * sigma(K+p->K+p)
 * Source: \iref{Buss:2011mx}, B.3.8 */
double kplusn_k0p(double mandelstam_s) {
  return 0.5 * kplusp_elastic_background(mandelstam_s);
}

/** K- p elastic cross section parametrization, PDG data.
 *
 * The PDG data is smoothed using the LOWESS algorithm. If more than one
 * cross section was given for one p_lab value, the corresponding cross sections
 * are averaged.
 */
static double kminusp_elastic_pdg(double mandelstam_s) {
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

/** K- p elastic background cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
double kminusp_elastic_background(double mandelstam_s) {
  const double p_lab = plab_from_s(mandelstam_s, kaon_mass, nucleon_mass);
  double sigma;
  if (std::sqrt(mandelstam_s) < 1.68) {
    // The parametrization here also works for anti-K0 n, Lambda pi0,
    // Sigma+ pi-, Sigma- pi+, Sigma0 pi0 with different parameters a0, a1, a2.
    //
    // The values of the parameters are *not* taken from the source above,
    // they come from a fit to PDG data.
    constexpr double a0 = 186.03567644;  // mb GeV^2
    constexpr double a1 = 0.22002795;    // Gev
    constexpr double a2 = 0.64907116;

    const double p_i = p_lab;
    const double p_f = p_lab;

    const double ratio = a1 * a1 / (a1 * a1 + p_f * p_f);
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

/** K- n elastic background cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
double kminusn_elastic_background(double) { return 4.0; }

/** K0 p elastic background cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
double k0p_elastic_background(double mandelstam_s) {
  // by isospin symmetry
  return kplusn_elastic_background(mandelstam_s);
}

/** K0 n elastic background cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
double k0n_elastic_background(double mandelstam_s) {
  // by isospin symmetry
  return kplusp_elastic_background(mandelstam_s);
}

/** Kbar0 p elastic background cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
double kbar0p_elastic_background(double mandelstam_s) {
  // by isospin symmetry
  return kminusn_elastic_background(mandelstam_s);
}

/** Kbar0 n elastic background cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
double kbar0n_elastic_background(double mandelstam_s) {
  // by isospin symmetry
  return kminusp_elastic_background(mandelstam_s);
}

/** K+ p inelastic background cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.8 */
double kplusp_inelastic_background(double mandelstam_s) {
  if (kplusp_total_interpolation == nullptr) {
    std::vector<double> x = KPLUSP_TOT_PLAB;
    std::vector<double> y = KPLUSP_TOT_SIG;
    std::vector<double> dedup_x;
    std::vector<double> dedup_y;
    std::tie(dedup_x, dedup_y) = dedup_avg(x, y);
    dedup_y = smooth(dedup_x, dedup_y, 0.1, 5);
    kplusp_total_interpolation =
        make_unique<InterpolateDataLinear<double>>(dedup_x, dedup_y);
  }
  const double p_lab = plab_from_s(mandelstam_s, kaon_mass, nucleon_mass);
  return (*kplusp_total_interpolation)(p_lab)-kplusp_elastic_background(
      mandelstam_s);
}

/** K+ n inelastic background cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.8 */
double kplusn_inelastic_background(double mandelstam_s) {
  if (kplusn_total_interpolation == nullptr) {
    std::vector<double> x = KPLUSN_TOT_PLAB;
    std::vector<double> y = KPLUSN_TOT_SIG;
    std::vector<double> dedup_x;
    std::vector<double> dedup_y;
    std::tie(dedup_x, dedup_y) = dedup_avg(x, y);
    dedup_y = smooth(dedup_x, dedup_y, 0.05, 5);
    kplusn_total_interpolation =
        make_unique<InterpolateDataLinear<double>>(dedup_x, dedup_y);
  }
  const double p_lab = plab_from_s(mandelstam_s, kaon_mass, nucleon_mass);
  return (*kplusn_total_interpolation)(p_lab)-kplusn_elastic_background(
             mandelstam_s) -
         kplusn_k0p(mandelstam_s);
}

/// Calculate and store all isospin ratios for K+ N reactions.
static void initialize(std::unordered_map<std::pair<uint64_t, uint64_t>, double,
                                          pair_hash>& ratios) {
  const auto& type_p = ParticleType::find(pdg::p);
  const auto& type_n = ParticleType::find(pdg::n);
  const auto& type_K_p = ParticleType::find(pdg::K_p);
  const auto& type_K_z = ParticleType::find(pdg::K_z);
  const auto& type_Delta_pp = ParticleType::find(pdg::Delta_pp);
  const auto& type_Delta_p = ParticleType::find(pdg::Delta_p);
  const auto& type_Delta_z = ParticleType::find(pdg::Delta_z);
  const auto& type_Delta_m = ParticleType::find(pdg::Delta_m);

  auto add_to_ratios = [&](const ParticleType& a, const ParticleType& b,
                           const ParticleType& c, const ParticleType& d,
                           double weight_numerator, double weight_other) {
    assert(weight_numerator + weight_other != 0);
    const auto key =
        std::make_pair(pack(a.pdgcode().code(), b.pdgcode().code()),
                       pack(c.pdgcode().code(), d.pdgcode().code()));
    const double ratio = weight_numerator / (weight_numerator + weight_other);
    ratios[key] = ratio;
  };

  // All inelastic channels are K+ N -> K Delta -> K pi N or charge exchange,
  // with identical cross section, weighted by the isospin factor.
  {
    const auto weight1 = isospin_clebsch_gordan_sqr_2to2(
        type_p, type_K_p, type_K_z, type_Delta_pp);
    const auto weight2 = isospin_clebsch_gordan_sqr_2to2(
        type_p, type_K_p, type_K_p, type_Delta_p);

    add_to_ratios(type_p, type_K_p, type_K_z, type_Delta_pp, weight1, weight2);
    add_to_ratios(type_p, type_K_p, type_K_p, type_Delta_p, weight2, weight1);
  }
  {
    const auto weight1 = isospin_clebsch_gordan_sqr_2to2(
        type_n, type_K_p, type_K_z, type_Delta_p);
    const auto weight2 = isospin_clebsch_gordan_sqr_2to2(
        type_n, type_K_p, type_K_p, type_Delta_z);

    add_to_ratios(type_n, type_K_p, type_K_z, type_Delta_p, weight1, weight2);
    add_to_ratios(type_n, type_K_p, type_K_p, type_Delta_z, weight2, weight1);
  }
  {
    const auto weight1 =
        isospin_clebsch_gordan_sqr_2to2(type_n, type_K_p, type_K_z, type_p);
    const auto weight2 =
        isospin_clebsch_gordan_sqr_2to2(type_p, type_K_z, type_K_p, type_n);

    add_to_ratios(type_n, type_K_p, type_K_z, type_p, weight1, weight2);
    add_to_ratios(type_p, type_K_z, type_K_p, type_n, weight2, weight1);
  }

  // K+ and K0 have the same isospin projection, they are assumed to have
  // the same cross section here.
  {
    const auto weight1 = isospin_clebsch_gordan_sqr_2to2(
        type_p, type_K_z, type_K_z, type_Delta_p);
    const auto weight2 = isospin_clebsch_gordan_sqr_2to2(
        type_p, type_K_z, type_K_p, type_Delta_z);

    add_to_ratios(type_p, type_K_z, type_K_z, type_Delta_p, weight1, weight2);
    add_to_ratios(type_p, type_K_z, type_K_p, type_Delta_z, weight2, weight1);
  }
  {
    const auto weight1 = isospin_clebsch_gordan_sqr_2to2(
        type_n, type_K_z, type_K_z, type_Delta_z);
    const auto weight2 = isospin_clebsch_gordan_sqr_2to2(
        type_n, type_K_z, type_K_p, type_Delta_m);

    add_to_ratios(type_n, type_K_z, type_K_z, type_Delta_z, weight1, weight2);
    add_to_ratios(type_n, type_K_z, type_K_p, type_Delta_m, weight2, weight1);
  }
}

/// Return the isospin ratio of the given K+ N reaction's cross section.
double KplusNRatios::get_ratio(const ParticleType& a, const ParticleType& b,
                               const ParticleType& c,
                               const ParticleType& d) const {
  /* If this method is called with anti-nucleons, flip all particles to
   * anti-particles;
   * the ratio is equal */
  int flip = 0;
  for (const auto& p : {&a, &b, &c, &d}) {
    if (p->is_nucleon()) {
      if (flip == 0) {
        flip = p->antiparticle_sign();
      } else {
        assert(p->antiparticle_sign() == flip);
      }
    }
  }
  const auto key = std::make_pair(
      pack(a.pdgcode().code() * flip, b.pdgcode().code() * flip),
      pack(c.pdgcode().code() * flip, d.pdgcode().code() * flip));
  if (ratios_.empty()) {
    initialize(ratios_);
  }
  return ratios_.at(key);
}

/*thread_local (see #3075)*/ KplusNRatios kplusn_ratios;

/** K- p -> Kbar0 n cross section parametrization.
 * Source: \iref{Buss:2011mx}, B.3.9 */
double kminusp_kbar0n(double mandelstam_s) {
  constexpr double a0 = 100;   // mb GeV^2
  constexpr double a1 = 0.15;  // GeV
  constexpr unsigned a2 = 2;

  const double p_lab = plab_from_s(mandelstam_s, kaon_mass, nucleon_mass);
  const double p_i = p_lab;
  const double p_f = p_lab;

  return a0 * p_f / (p_i * mandelstam_s) *
         pow_int(a1 * a1 / (a1 * a1 + p_f * p_f), a2);
}

/* Parametrizations of strangeness exchange channels
 *
 * Taken from UrQMD (\iref{Graef:2014mra}).
 */

double kminusp_piminussigmaplus(double sqrts) {
  return 0.0788265 / smash::square(sqrts - 1.38841);
}

double kminusp_piplussigmaminus(double sqrts) {
  return 0.0196741 / smash::square(sqrts - 1.42318);
}

double kminusp_pi0sigma0(double sqrts) {
  // Fit to Landolt-Börnstein instead of UrQMD values
  return 0.0403364 / smash::square(sqrts - 1.39830305);
}

double kminusp_pi0lambda(double sqrts) {
  // Fit to Landolt-Börnstein instead of UrQMD values
  return 0.05932562 / smash::square(sqrts - 1.38786692);
}

// The other channels follow from the paramatriziation with the same strange
// product via isospin symmetry.

double kminusn_piminussigma0(double sqrts) {
  return 1. / 6 * 2 * kminusp_pi0sigma0(sqrts);
}

double kminusn_pi0sigmaminus(double sqrts) {
  return (0.25 + 1. / 6) * 2 * kminusp_piplussigmaminus(sqrts);
}

double kminusn_piminuslambda(double sqrts) {
  return 0.5 * kminusp_pi0lambda(sqrts);
}

// All K+ p and K+ n channels are forbidden by isospin.

// Two hyperon exchange, based on effective model by Feng Li,
// as in UrQMD (\iref{Graef:2014mra}).

double lambdalambda_ximinusp(double sqrts_sqrts0, double p_N, double p_lambda) {
  assert(p_lambda != 0);
  assert(sqrts_sqrts0 >= 0);
  return 37.15 / 2 * p_N / p_lambda * std::pow(sqrts_sqrts0, -0.16);
}

double lambdalambda_xi0n(double sqrts_sqrts0, double p_N, double p_lambda) {
  return lambdalambda_ximinusp(sqrts_sqrts0, p_N, p_lambda);
}

double lambdasigmaplus_xi0p(double sqrts_sqrts0) {
  assert(sqrts_sqrts0 >= 0);
  return 24.3781 * std::pow(sqrts_sqrts0, -0.479);
}

double lambdasigmaminus_ximinusn(double sqrts_sqrts0) {
  return lambdasigmaplus_xi0p(sqrts_sqrts0);
}

double lambdasigma0_ximinusp(double sqrts_sqrts0) {
  assert(sqrts_sqrts0 >= 0);
  if (sqrts_sqrts0 < 0.03336) {
    return 6.475 * std::pow(sqrts_sqrts0, -0.4167);
  } else {
    return 14.5054 * std::pow(sqrts_sqrts0, -0.1795);
  }
}

double lambdasigma0_xi0n(double sqrts_sqrts0) {
  return lambdasigma0_ximinusp(sqrts_sqrts0);
}

double sigma0sigma0_ximinusp(double sqrts_sqrts0) {
  assert(sqrts_sqrts0 >= 0);
  if (sqrts_sqrts0 < 0.09047) {
    return 5.625 * std::pow(sqrts_sqrts0, -0.318);
  } else {
    return 4.174 * std::pow(sqrts_sqrts0, -0.4421);
  }
}

// Note that there is a typo in the paper in equation (6):
// "Lambda Sigma0 -> Xi0 n" should be "Sigma0 Sigma0 -> Xi0 n".
double sigma0sigma0_xi0n(double sqrts_sqrts0) {
  return sigma0sigma0_ximinusp(sqrts_sqrts0);
}

double sigmaplussigmaminus_xi0p(double sqrts_sqrts0) {
  return 4 * sigma0sigma0_ximinusp(sqrts_sqrts0);
}

double sigma0sigmaminus_ximinusn(double sqrts_sqrts0) {
  return 4 * sigma0sigma0_ximinusp(sqrts_sqrts0);
}

double sigmaplussigmaminus_ximinusp(double sqrts_sqrts0) {
  return 14.194 * std::pow(sqrts_sqrts0, -0.442);
}

double sigmaplussigmaminus_xi0n(double sqrts_sqrts0) {
  return sigmaplussigmaminus_ximinusp(sqrts_sqrts0);
}

}  // namespace smash
