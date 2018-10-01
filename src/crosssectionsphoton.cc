/*
 *
 *    Copyright (c) 2016-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/crosssectionsphoton.h"
#include <memory>
#include "smash/constants.h"
#include "smash/logging.h"
#include "smash/particletype.h"
#include "smash/pdgcode.h"

namespace {

/* Necessary for the implementation of a non-stable rho meson in the
 * pi0 + rho0 -> omega -> pi0 + gamma and pi + pi0 -> rho + gamma channel.
 * For these specific scattering processes, there are s and t channels with
 * different thresholds. While the t-channel can always be performed, the
 * s-channel is kinematically only accessible if sqrt(s) >= mass of the exchange
 * particle. The corresponding s-channels need to be excluded from the
 * cross sections below their specific thresholds. */

/**
 * Heavyside step function
 *
 * \param[in] x value to be compared to zero
 * \return 0 if x is smaller than 0, else 1
 */
double HeavisideTheta(double x) {
  if (x >= 0.0) {
    return 1.0;
  } else {
    return 0.0;
  }
}

/**
 * Cross section after cut off
 *
 * Cross sections larger than a certain value are cut off in smash. Either the
 * cross section is returned or, if the cross section is larger than the cut
 * off, the cut off value is returned
 *
 * \param[in] sigma_mb cross section before cut off [mb]
 * \return Cross section after cut off [mb]
 */
double cut_off(const double sigma_mb) {
  return (sigma_mb > smash::maximum_cross_section)
             ? smash::maximum_cross_section
             : sigma_mb;
}

}  // anonymous namespace

/*
   The cross sections presented in this file are calculated applying an average
   over initial states and sum over final states. In transport simulations
   individual particles are propagated which have a specific degeneracy
   state such that an average over initial states is superflous. The cross
   sections need thus be divided by the spin degeneracy factor of the initial
   particles to account for these particle properties.
   The C-naming is in accordance to the matrix elements from
   (\iref{Turbide:2006})
*/

namespace smash {
// template class CrosssectionsPhoton<ComputationMethod::Analytic>; (really
// remove this?)

/*----------------------------------------------------------------------------*/
/*                               Pi + Rho -> Pi + Photon channels mediated by */
/*                                                      (Pi, Rho, a1)
 */
/*----------------------------------------------------------------------------*/

// C11
double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi_rho0_pi(
    const double s, const double m_rho) {
  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  auto t_mandelstam = get_t_range(sqrt(s), pion_mass, m_rho, pion_mass, 0.);
  const double &tmax = t_mandelstam[0];
  const double &tmin = t_mandelstam[1];

  const double spin_deg_factor = 3.0;

  // clang-format off
  const double xs =
      (pow(Const, 2) * pow(ghat, 4) *
       ((pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(a1_mass, 8) + pow(pion_mass, 8) -
               pow(pion_mass, 4) * pow(m_rho, 4) -
               2 * pow(a1_mass, 2) * pow(pion_mass, 2) *
                   (pow(pion_mass, 2) - pow(m_rho, 2)) * (pow(m_rho, 2) + s) +
               pow(a1_mass, 6) * (-4 * pow(pion_mass, 2) + 2 * s) +
               pow(a1_mass, 4) *
                   (4 * pow(pion_mass, 4) - pow(m_rho, 4) +
                    2 * pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) -
                    2 * pow(m_rho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(a1_mass, 8) +
               pow(pion_mass, 4) * pow(pow(pion_mass, 2) - pow(m_rho, 2), 2) +
               2 * pow(a1_mass, 6) *
                   (-2 * pow(pion_mass, 2) + pow(m_rho, 2) + s) +
               2 * pow(a1_mass, 2) * pow(pion_mass, 2) *
                   (-pow(m_rho, 4) +
                    pow(pion_mass, 2) * (2 * pow(m_rho, 2) - s) +
                    pow(m_rho, 2) * s) +
               pow(a1_mass, 4) *
                   (4 * pow(pion_mass, 4) + pow(m_rho, 4) -
                    2 * pow(m_rho, 2) * s + 2 * pow(s, 2) -
                    4 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
          pow(eta1, 2) *
              (pow(a1_mass, 8) + pow(pion_mass, 8) -
               2 * pow(pion_mass, 6) * pow(m_rho, 2) -
               2 * pow(a1_mass, 6) *
                   (2 * pow(pion_mass, 2) + pow(m_rho, 2) - s) -
               2 * pow(pion_mass, 2) * pow(m_rho, 4) * s +
               pow(pion_mass, 4) * (3 * pow(m_rho, 4) + 2 * pow(m_rho, 2) * s) +
               pow(a1_mass, 4) *
                   (4 * pow(pion_mass, 4) + pow(m_rho, 4) +
                    pow(pion_mass, 2) * (8 * pow(m_rho, 2) - 4 * s) -
                    4 * pow(m_rho, 2) * s + 2 * pow(s, 2)) -
               2 * pow(a1_mass, 2) *
                   (pow(m_rho, 2) * s * (-pow(m_rho, 2) + s) +
                    pow(pion_mass, 4) * (3 * pow(m_rho, 2) + s) +
                    pow(pion_mass, 2) *
                        (2 * pow(m_rho, 4) - 3 * pow(m_rho, 2) * s))))) /
            ((pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
             (pow(a1_mass, 2) - tmax)) +
        (8 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
         (4 * pow(pion_mass, 2) - pow(m_rho, 2))) /
            ((pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
             (pow(pion_mass, 2) - tmax)) -
        (8 * pow(-2 + delta, 2) * pow(pion_mass, 2) * tmax) /
            (pow(m_rho, 2) * pow(pow(pion_mass, 2) - s, 2)) -
        (8 * pow(-2 + delta, 2) * pow(pion_mass, 2) * tmax) /
            (pow(m_rho, 2) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (8 * (-2 + delta) *
         (-8 * C4 * pow(m_rho, 4) +
          pow(pion_mass, 2) * (2 + delta - 8 * C4 * pow(m_rho, 2)) -
          (2 + 3 * delta) * s +
          pow(m_rho, 2) * (-2 + 3 * delta + 16 * C4 * s)) *
         tmax) /
            (pow(m_rho, 2) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (eta2 * (pow(pion_mass, 2) + pow(m_rho, 2) - 2 * s) *
              (pow(pion_mass, 2) + s) +
          eta1 *
              (-2 * pow(pion_mass, 4) + pow(m_rho, 4) - 3 * pow(m_rho, 2) * s +
               2 * pow(s, 2) + pow(pion_mass, 2) * (pow(m_rho, 2) + s))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (3 * pow(a1_mass, 4) + 4 * pow(pion_mass, 4) + pow(m_rho, 4) +
               pow(pion_mass, 2) * (8 * pow(m_rho, 2) - 4 * s) -
               4 * pow(a1_mass, 2) *
                   (2 * pow(pion_mass, 2) + pow(m_rho, 2) - s) -
               4 * pow(m_rho, 2) * s + 2 * pow(s, 2)) +
          pow(eta2, 2) *
              (3 * pow(a1_mass, 4) + 4 * pow(pion_mass, 4) + pow(m_rho, 4) -
               2 * pow(m_rho, 2) * s + 2 * pow(s, 2) -
               4 * pow(pion_mass, 2) * (pow(m_rho, 2) + s) +
               4 * pow(a1_mass, 2) *
                   (-2 * pow(pion_mass, 2) + pow(m_rho, 2) + s)) -
          2 * eta1 * eta2 *
              (3 * pow(a1_mass, 4) + 4 * pow(pion_mass, 4) - pow(m_rho, 4) +
               2 * pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) -
               2 * pow(m_rho, 2) * s + 2 * pow(s, 2) +
               pow(a1_mass, 2) * (-8 * pow(pion_mass, 2) + 4 * s))) *
         tmax) /
            (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
             2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) +
        (8 *
         (pow(delta, 2) * (8 * pow(pion_mass, 4) + 3 * pow(m_rho, 4) +
                           4 * pow(pion_mass, 2) * (3 * pow(m_rho, 2) - 2 * s) -
                           6 * pow(m_rho, 2) * s + 2 * pow(s, 2)) +
          4 * pow(m_rho, 4) *
              (3 + 12 * C4 * (2 * pow(pion_mass, 2) - s) +
               8 * pow(C4, 2) * pow(-2 * pow(pion_mass, 2) + s, 2)) -
          4 * delta * pow(m_rho, 2) *
              (16 * C4 * pow(pion_mass, 4) +
               2 * pow(pion_mass, 2) *
                   (3 + 6 * C4 * pow(m_rho, 2) - 8 * C4 * s) +
               pow(m_rho, 2) * (3 - 6 * C4 * s) + s * (-3 + 4 * C4 * s))) *
         tmax) /
            (pow(m_rho, 4) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (8 * (-2 + delta) *
         (pow(pion_mass, 4) * (-2 + 3 * delta - 8 * C4 * pow(m_rho, 2)) +
          (pow(m_rho, 2) - s) * ((-2 + 3 * delta) * s +
                                 pow(m_rho, 2) * (-2 + delta - 8 * C4 * s)) +
          4 * pow(pion_mass, 2) *
              (2 * C4 * pow(m_rho, 4) + delta * s -
               pow(m_rho, 2) * (-1 + delta + 4 * C4 * s))) *
         tmax) /
            (pow(m_rho, 2) * (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (eta2 * (pow(pion_mass, 2) + s) *
              (pow(pion_mass, 4) - pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) +
               (pow(m_rho, 2) - s) * s) +
          eta1 * (-4 * pow(pion_mass, 6) + pow(pow(m_rho, 2) - s, 2) * s +
                  pow(pion_mass, 4) * (3 * pow(m_rho, 2) + s) -
                  pow(pion_mass, 2) *
                      (pow(m_rho, 4) - pow(m_rho, 2) * s + 2 * pow(s, 2)))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(pion_mass, 8) - pow(m_rho, 4) * pow(s, 2) + pow(s, 4) -
               pow(pion_mass, 4) *
                   (pow(m_rho, 4) + 2 * pow(m_rho, 2) * s - 4 * pow(s, 2)) +
               2 * pow(pion_mass, 2) * s *
                   (pow(m_rho, 4) + pow(m_rho, 2) * s - 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(s, 2) * pow(pow(m_rho, 2) + s, 2) +
               pow(pion_mass, 4) * pow(pow(m_rho, 2) + 2 * s, 2) -
               2 * pow(pion_mass, 2) * s *
                   (pow(m_rho, 4) + 2 * pow(m_rho, 2) * s + 2 * pow(s, 2))) +
          pow(eta1, 2) *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) -
               4 * pow(pion_mass, 2) * pow(pow(m_rho, 2) - s, 2) * s +
               pow(pow(m_rho, 2) - s, 2) * pow(s, 2) +
               pow(pion_mass, 4) * (3 * pow(m_rho, 4) - 6 * pow(m_rho, 2) * s +
                                    4 * pow(s, 2)))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) * (pow(a1_mass, 2) - s) *
         (pow(eta1, 2) *
              (pow(a1_mass, 4) * s +
               pow(pion_mass, 4) * (-3 * pow(m_rho, 2) + 2 * s) +
               s * (2 * pow(m_rho, 4) - 3 * pow(m_rho, 2) * s + pow(s, 2)) -
               2 * pow(pion_mass, 2) *
                   (pow(m_rho, 4) - 4 * pow(m_rho, 2) * s + 2 * pow(s, 2)) +
               pow(a1_mass, 2) *
                   (2 * pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) +
                    3 * s * (-pow(m_rho, 2) + s))) -
          2 * eta1 * eta2 *
              (pow(a1_mass, 4) * s +
               s * (2 * pow(pion_mass, 4) +
                    4 * pow(pion_mass, 2) * (pow(m_rho, 2) - s) +
                    s * (-2 * pow(m_rho, 2) + s)) +
               pow(a1_mass, 2) * (pow(pion_mass, 2) * (pow(m_rho, 2) - 4 * s) +
                                  s * (-2 * pow(m_rho, 2) + 3 * s))) +
          pow(eta2, 2) * (-4 * pow(pion_mass, 2) * s *
                              (pow(a1_mass, 2) + pow(m_rho, 2) + s) +
                          pow(pion_mass, 4) * (pow(m_rho, 2) + 2 * s) +
                          s * (pow(a1_mass, 4) + s * (pow(m_rho, 2) + s) +
                               pow(a1_mass, 2) * (pow(m_rho, 2) + 3 * s)))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (eta1 * (-4 * pow(pion_mass, 4) *
                      (6 * C4 * pow(m_rho, 4) + 2 * delta * s +
                       pow(m_rho, 2) * (1 - 2 * delta - 8 * C4 * s)) +
                  2 * pow(pion_mass, 2) *
                      (4 * delta * pow(s, 2) +
                       pow(m_rho, 2) * s * (6 - 7 * delta - 16 * C4 * s) +
                       2 * pow(m_rho, 4) * (-2 + delta + 8 * C4 * s)) -
                  (pow(m_rho, 2) - s) * s *
                      (-2 * delta * s +
                       pow(m_rho, 2) * (-6 + 3 * delta + 8 * C4 * s))) +
          eta2 * (delta *
                      (2 * pow(pion_mass, 4) * (pow(m_rho, 2) + 4 * s) +
                       pow(pion_mass, 2) * (2 * pow(m_rho, 4) +
                                            pow(m_rho, 2) * s - 8 * pow(s, 2)) +
                       s * (-2 * pow(m_rho, 4) - pow(m_rho, 2) * s +
                            2 * pow(s, 2))) -
                  2 * pow(m_rho, 2) *
                      (4 * C4 * pow(pion_mass, 4) * (pow(m_rho, 2) + 4 * s) +
                       pow(pion_mass, 2) * (s * (5 - 16 * C4 * s) +
                                            pow(m_rho, 2) * (2 - 8 * C4 * s)) +
                       s * (s * (-3 + 4 * C4 * s) +
                            pow(m_rho, 2) * (-2 + 4 * C4 * s))))) *
         tmax) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (8 * (eta1 - eta2) *
         (delta *
              (eta1 * (4 * pow(pion_mass, 6) +
                       pow(pion_mass, 4) * (7 * pow(m_rho, 2) - 8 * s) +
                       pow(a1_mass, 4) * (pow(pion_mass, 2) - s) -
                       pow(a1_mass, 2) *
                           (2 * pow(pion_mass, 2) + pow(m_rho, 2) - 2 * s) *
                           (2 * pow(pion_mass, 2) - s) +
                       pow(pion_mass, 2) * s * (-8 * pow(m_rho, 2) + 5 * s) +
                       s * (pow(m_rho, 4) + pow(m_rho, 2) * s - pow(s, 2))) +
               eta2 *
                   (-4 * pow(pion_mass, 6) -
                    pow(pion_mass, 4) * (pow(m_rho, 2) - 8 * s) +
                    pow(a1_mass, 4) * (-pow(pion_mass, 2) + s) +
                    pow(pion_mass, 2) * (2 * pow(m_rho, 4) - 5 * pow(s, 2)) +
                    s * (-2 * pow(m_rho, 4) + pow(m_rho, 2) * s + pow(s, 2)) +
                    pow(a1_mass, 2) *
                        (4 * pow(pion_mass, 4) - 6 * pow(pion_mass, 2) * s +
                         s * (pow(m_rho, 2) + 2 * s)))) -
          2 * pow(m_rho, 2) *
              (eta1 *
                   (8 * C4 * pow(pion_mass, 6) +
                    pow(pion_mass, 4) * (3 + 8 * C4 * (pow(m_rho, 2) - 2 * s)) +
                    2 * C4 * pow(a1_mass, 4) * (pow(pion_mass, 2) - s) +
                    2 * pow(pion_mass, 2) * s *
                        (-1 - 6 * C4 * pow(m_rho, 2) + 5 * C4 * s) -
                    pow(a1_mass, 2) *
                        (8 * C4 * pow(pion_mass, 4) +
                         pow(pion_mass, 2) *
                             (1 + 2 * C4 * (pow(m_rho, 2) - 6 * s)) +
                         2 * C4 * s * (-pow(m_rho, 2) + 2 * s)) +
                    s * (-(s * (1 + 2 * C4 * s)) +
                         pow(m_rho, 2) * (1 + 4 * C4 * s))) +
               eta2 * (2 * C4 * pow(a1_mass, 4) * (-pow(pion_mass, 2) + s) -
                       (pow(pion_mass, 2) - s) *
                           (8 * C4 * pow(pion_mass, 4) - 2 * pow(m_rho, 2) + s +
                            2 * C4 * pow(s, 2) +
                            pow(pion_mass, 2) *
                                (3 - 4 * C4 * (pow(m_rho, 2) + 2 * s))) +
                       pow(a1_mass, 2) *
                           (8 * C4 * pow(pion_mass, 4) +
                            2 * C4 * s * (pow(m_rho, 2) + 2 * s) +
                            pow(pion_mass, 2) *
                                (1 - 2 * C4 * (pow(m_rho, 2) + 6 * s)))))) *
         tmax) /
            (pow(m_rho, 2) * (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (8 * (-2 + delta) * (delta - 4 * C4 * pow(m_rho, 2)) * pow(tmax, 2)) /
            (pow(m_rho, 2) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (16 * (-2 + delta) * (delta - 4 * C4 * pow(m_rho, 2)) * s *
         pow(tmax, 2)) /
            (pow(m_rho, 2) * (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) * (pow(m_rho, 2) - s) + 2 * eta1 * eta2 * s -
          pow(eta2, 2) * s) *
         (pow(pion_mass, 4) - pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) +
          (pow(m_rho, 2) - s) * s) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (-(eta1 * (pow(pion_mass, 2) + 2 * pow(m_rho, 2) - 3 * s)) -
          eta2 * (pow(pion_mass, 2) + s)) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (pow(pion_mass, 2) + s) *
         (-2 * eta2 * s + eta1 * (pow(pion_mass, 2) - pow(m_rho, 2) + s)) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (pow(eta1 - eta2, 3) *
         (eta1 * (pow(a1_mass, 2) - 2 * pow(pion_mass, 2) - pow(m_rho, 2) + s) -
          eta2 *
              (pow(a1_mass, 2) - 2 * pow(pion_mass, 2) + pow(m_rho, 2) + s)) *
         pow(tmax, 2)) /
            (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
             2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) -
        (8 * (delta - 4 * C4 * pow(m_rho, 2)) *
         (delta * (4 * pow(pion_mass, 2) + 3 * pow(m_rho, 2) - 2 * s) -
          2 * pow(m_rho, 2) * (3 + 8 * C4 * pow(pion_mass, 2) - 4 * C4 * s)) *
         pow(tmax, 2)) /
            (pow(m_rho, 4) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (pow(eta1 - eta2, 2) * (pow(a1_mass, 2) - s) *
         (pow(eta2, 2) * s *
              (pow(a1_mass, 2) - 4 * pow(pion_mass, 2) + pow(m_rho, 2) +
               3 * s) +
          pow(eta1, 2) * (2 * pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) +
                          s * (pow(a1_mass, 2) - 3 * pow(m_rho, 2) + 3 * s)) -
          2 * eta1 * eta2 *
              (pow(pion_mass, 2) * (pow(m_rho, 2) - 4 * s) +
               s * (pow(a1_mass, 2) - 2 * pow(m_rho, 2) + 3 * s))) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (2 * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (eta1 * (4 * delta * pow(s, 2) -
                  2 * pow(m_rho, 2) * s * (-2 + 3 * delta + 8 * C4 * s) +
                  pow(m_rho, 4) * (-2 + delta + 16 * C4 * s) -
                  2 * pow(pion_mass, 2) *
                      (8 * C4 * pow(m_rho, 4) + 4 * delta * s +
                       pow(m_rho, 2) * (2 - 3 * delta - 16 * C4 * s))) +
          eta2 * (pow(pion_mass, 2) *
                      (8 * delta * s +
                       pow(m_rho, 2) * (-2 + delta - 32 * C4 * s)) +
                  s * (-4 * delta * s +
                       pow(m_rho, 2) * (-2 + delta + 16 * C4 * s)))) *
         pow(tmax, 2)) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (4 * (eta1 - eta2) *
         (delta * (eta1 * (pow(a1_mass, 2) * (pow(pion_mass, 2) - s) -
                           (2 * pow(pion_mass, 2) + pow(m_rho, 2) - 2 * s) *
                               (2 * pow(pion_mass, 2) - s)) +
                   eta2 * (4 * pow(pion_mass, 4) - 6 * pow(pion_mass, 2) * s +
                           pow(a1_mass, 2) * (-pow(pion_mass, 2) + s) +
                           s * (pow(m_rho, 2) + 2 * s))) +
          2 * pow(m_rho, 2) *
              (eta1 * (8 * C4 * pow(pion_mass, 4) +
                       2 * C4 * s * (pow(a1_mass, 2) - pow(m_rho, 2) + 2 * s) +
                       pow(pion_mass, 2) * (1 - 2 * C4 *
                                                    (pow(a1_mass, 2) -
                                                     pow(m_rho, 2) + 6 * s))) -
               eta2 * (8 * C4 * pow(pion_mass, 4) +
                       2 * C4 * s * (pow(a1_mass, 2) + pow(m_rho, 2) + 2 * s) -
                       pow(pion_mass, 2) *
                           (-1 +
                            2 * C4 *
                                (pow(a1_mass, 2) + pow(m_rho, 2) + 6 * s))))) *
         pow(tmax, 2)) /
            (pow(m_rho, 2) * (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (pow(eta1 - eta2, 4) * pow(tmax, 3)) /
            (3. * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                   2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (8 * pow(eta1 - eta2, 2) * (delta - 4 * C4 * pow(m_rho, 2)) *
         pow(tmax, 3)) /
            (3. * pow(m_rho, 2) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (16 * pow(delta - 4 * C4 * pow(m_rho, 2), 2) * pow(tmax, 3)) /
            (3. * pow(m_rho, 4) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (4 * (-2 + delta) * eta1 * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         pow(tmax, 3)) /
            (3. *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (2 * pow(eta1 - eta2, 4) * (pow(a1_mass, 2) - s) * s * pow(tmax, 3)) /
            (3. *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) * s *
         (-2 * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-pow(m_rho, 2) + s)) *
         pow(tmax, 3)) /
            (3. *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (2 * eta2 * (delta - 4 * C4 * pow(m_rho, 2)) * s +
          eta1 * (-2 * delta * s + pow(m_rho, 2) * (-2 + delta + 8 * C4 * s))) *
         pow(tmax, 3)) /
            (3. * pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(a1_mass, 8) + pow(pion_mass, 8) -
               pow(pion_mass, 4) * pow(m_rho, 4) -
               2 * pow(a1_mass, 2) * pow(pion_mass, 2) *
                   (pow(pion_mass, 2) - pow(m_rho, 2)) * (pow(m_rho, 2) + s) +
               pow(a1_mass, 6) * (-4 * pow(pion_mass, 2) + 2 * s) +
               pow(a1_mass, 4) *
                   (4 * pow(pion_mass, 4) - pow(m_rho, 4) +
                    2 * pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) -
                    2 * pow(m_rho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(a1_mass, 8) +
               pow(pion_mass, 4) * pow(pow(pion_mass, 2) - pow(m_rho, 2), 2) +
               2 * pow(a1_mass, 6) *
                   (-2 * pow(pion_mass, 2) + pow(m_rho, 2) + s) +
               2 * pow(a1_mass, 2) * pow(pion_mass, 2) *
                   (-pow(m_rho, 4) +
                    pow(pion_mass, 2) * (2 * pow(m_rho, 2) - s) +
                    pow(m_rho, 2) * s) +
               pow(a1_mass, 4) *
                   (4 * pow(pion_mass, 4) + pow(m_rho, 4) -
                    2 * pow(m_rho, 2) * s + 2 * pow(s, 2) -
                    4 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
          pow(eta1, 2) *
              (pow(a1_mass, 8) + pow(pion_mass, 8) -
               2 * pow(pion_mass, 6) * pow(m_rho, 2) -
               2 * pow(a1_mass, 6) *
                   (2 * pow(pion_mass, 2) + pow(m_rho, 2) - s) -
               2 * pow(pion_mass, 2) * pow(m_rho, 4) * s +
               pow(pion_mass, 4) * (3 * pow(m_rho, 4) + 2 * pow(m_rho, 2) * s) +
               pow(a1_mass, 4) *
                   (4 * pow(pion_mass, 4) + pow(m_rho, 4) +
                    pow(pion_mass, 2) * (8 * pow(m_rho, 2) - 4 * s) -
                    4 * pow(m_rho, 2) * s + 2 * pow(s, 2)) -
               2 * pow(a1_mass, 2) *
                   (pow(m_rho, 2) * s * (-pow(m_rho, 2) + s) +
                    pow(pion_mass, 4) * (3 * pow(m_rho, 2) + s) +
                    pow(pion_mass, 2) *
                        (2 * pow(m_rho, 4) - 3 * pow(m_rho, 2) * s))))) /
            ((pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
             (pow(a1_mass, 2) - tmin)) -
        (8 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
         (4 * pow(pion_mass, 2) - pow(m_rho, 2))) /
            ((pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
             (pow(pion_mass, 2) - tmin)) +
        (8 * pow(-2 + delta, 2) * pow(pion_mass, 2) * tmin) /
            (pow(m_rho, 2) * pow(pow(pion_mass, 2) - s, 2)) +
        (8 * pow(-2 + delta, 2) * pow(pion_mass, 2) * tmin) /
            (pow(m_rho, 2) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (8 * (-2 + delta) *
         (-8 * C4 * pow(m_rho, 4) +
          pow(pion_mass, 2) * (2 + delta - 8 * C4 * pow(m_rho, 2)) -
          (2 + 3 * delta) * s +
          pow(m_rho, 2) * (-2 + 3 * delta + 16 * C4 * s)) *
         tmin) /
            (pow(m_rho, 2) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (eta2 * (pow(pion_mass, 2) + pow(m_rho, 2) - 2 * s) *
              (pow(pion_mass, 2) + s) +
          eta1 *
              (-2 * pow(pion_mass, 4) + pow(m_rho, 4) - 3 * pow(m_rho, 2) * s +
               2 * pow(s, 2) + pow(pion_mass, 2) * (pow(m_rho, 2) + s))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (3 * pow(a1_mass, 4) + 4 * pow(pion_mass, 4) + pow(m_rho, 4) +
               pow(pion_mass, 2) * (8 * pow(m_rho, 2) - 4 * s) -
               4 * pow(a1_mass, 2) *
                   (2 * pow(pion_mass, 2) + pow(m_rho, 2) - s) -
               4 * pow(m_rho, 2) * s + 2 * pow(s, 2)) +
          pow(eta2, 2) *
              (3 * pow(a1_mass, 4) + 4 * pow(pion_mass, 4) + pow(m_rho, 4) -
               2 * pow(m_rho, 2) * s + 2 * pow(s, 2) -
               4 * pow(pion_mass, 2) * (pow(m_rho, 2) + s) +
               4 * pow(a1_mass, 2) *
                   (-2 * pow(pion_mass, 2) + pow(m_rho, 2) + s)) -
          2 * eta1 * eta2 *
              (3 * pow(a1_mass, 4) + 4 * pow(pion_mass, 4) - pow(m_rho, 4) +
               2 * pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) -
               2 * pow(m_rho, 2) * s + 2 * pow(s, 2) +
               pow(a1_mass, 2) * (-8 * pow(pion_mass, 2) + 4 * s))) *
         tmin) /
            (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
             2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) -
        (8 *
         (pow(delta, 2) * (8 * pow(pion_mass, 4) + 3 * pow(m_rho, 4) +
                           4 * pow(pion_mass, 2) * (3 * pow(m_rho, 2) - 2 * s) -
                           6 * pow(m_rho, 2) * s + 2 * pow(s, 2)) +
          4 * pow(m_rho, 4) *
              (3 + 12 * C4 * (2 * pow(pion_mass, 2) - s) +
               8 * pow(C4, 2) * pow(-2 * pow(pion_mass, 2) + s, 2)) -
          4 * delta * pow(m_rho, 2) *
              (16 * C4 * pow(pion_mass, 4) +
               2 * pow(pion_mass, 2) *
                   (3 + 6 * C4 * pow(m_rho, 2) - 8 * C4 * s) +
               pow(m_rho, 2) * (3 - 6 * C4 * s) + s * (-3 + 4 * C4 * s))) *
         tmin) /
            (pow(m_rho, 4) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (8 * (-2 + delta) *
         (pow(pion_mass, 4) * (-2 + 3 * delta - 8 * C4 * pow(m_rho, 2)) +
          (pow(m_rho, 2) - s) * ((-2 + 3 * delta) * s +
                                 pow(m_rho, 2) * (-2 + delta - 8 * C4 * s)) +
          4 * pow(pion_mass, 2) *
              (2 * C4 * pow(m_rho, 4) + delta * s -
               pow(m_rho, 2) * (-1 + delta + 4 * C4 * s))) *
         tmin) /
            (pow(m_rho, 2) * (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (eta2 * (pow(pion_mass, 2) + s) *
              (pow(pion_mass, 4) - pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) +
               (pow(m_rho, 2) - s) * s) +
          eta1 * (-4 * pow(pion_mass, 6) + pow(pow(m_rho, 2) - s, 2) * s +
                  pow(pion_mass, 4) * (3 * pow(m_rho, 2) + s) -
                  pow(pion_mass, 2) *
                      (pow(m_rho, 4) - pow(m_rho, 2) * s + 2 * pow(s, 2)))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(pion_mass, 8) - pow(m_rho, 4) * pow(s, 2) + pow(s, 4) -
               pow(pion_mass, 4) *
                   (pow(m_rho, 4) + 2 * pow(m_rho, 2) * s - 4 * pow(s, 2)) +
               2 * pow(pion_mass, 2) * s *
                   (pow(m_rho, 4) + pow(m_rho, 2) * s - 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(s, 2) * pow(pow(m_rho, 2) + s, 2) +
               pow(pion_mass, 4) * pow(pow(m_rho, 2) + 2 * s, 2) -
               2 * pow(pion_mass, 2) * s *
                   (pow(m_rho, 4) + 2 * pow(m_rho, 2) * s + 2 * pow(s, 2))) +
          pow(eta1, 2) *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) -
               4 * pow(pion_mass, 2) * pow(pow(m_rho, 2) - s, 2) * s +
               pow(pow(m_rho, 2) - s, 2) * pow(s, 2) +
               pow(pion_mass, 4) * (3 * pow(m_rho, 4) - 6 * pow(m_rho, 2) * s +
                                    4 * pow(s, 2)))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) * (pow(a1_mass, 2) - s) *
         (pow(eta1, 2) *
              (pow(a1_mass, 4) * s +
               pow(pion_mass, 4) * (-3 * pow(m_rho, 2) + 2 * s) +
               s * (2 * pow(m_rho, 4) - 3 * pow(m_rho, 2) * s + pow(s, 2)) -
               2 * pow(pion_mass, 2) *
                   (pow(m_rho, 4) - 4 * pow(m_rho, 2) * s + 2 * pow(s, 2)) +
               pow(a1_mass, 2) *
                   (2 * pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) +
                    3 * s * (-pow(m_rho, 2) + s))) -
          2 * eta1 * eta2 *
              (pow(a1_mass, 4) * s +
               s * (2 * pow(pion_mass, 4) +
                    4 * pow(pion_mass, 2) * (pow(m_rho, 2) - s) +
                    s * (-2 * pow(m_rho, 2) + s)) +
               pow(a1_mass, 2) * (pow(pion_mass, 2) * (pow(m_rho, 2) - 4 * s) +
                                  s * (-2 * pow(m_rho, 2) + 3 * s))) +
          pow(eta2, 2) * (-4 * pow(pion_mass, 2) * s *
                              (pow(a1_mass, 2) + pow(m_rho, 2) + s) +
                          pow(pion_mass, 4) * (pow(m_rho, 2) + 2 * s) +
                          s * (pow(a1_mass, 4) + s * (pow(m_rho, 2) + s) +
                               pow(a1_mass, 2) * (pow(m_rho, 2) + 3 * s)))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (eta1 * (4 * pow(pion_mass, 4) *
                      (6 * C4 * pow(m_rho, 4) + 2 * delta * s +
                       pow(m_rho, 2) * (1 - 2 * delta - 8 * C4 * s)) -
                  2 * pow(pion_mass, 2) *
                      (4 * delta * pow(s, 2) +
                       pow(m_rho, 2) * s * (6 - 7 * delta - 16 * C4 * s) +
                       2 * pow(m_rho, 4) * (-2 + delta + 8 * C4 * s)) +
                  (pow(m_rho, 2) - s) * s *
                      (-2 * delta * s +
                       pow(m_rho, 2) * (-6 + 3 * delta + 8 * C4 * s))) +
          eta2 * (-(delta *
                    (2 * pow(pion_mass, 4) * (pow(m_rho, 2) + 4 * s) +
                     pow(pion_mass, 2) * (2 * pow(m_rho, 4) +
                                          pow(m_rho, 2) * s - 8 * pow(s, 2)) +
                     s * (-2 * pow(m_rho, 4) - pow(m_rho, 2) * s +
                          2 * pow(s, 2)))) +
                  2 * pow(m_rho, 2) *
                      (4 * C4 * pow(pion_mass, 4) * (pow(m_rho, 2) + 4 * s) +
                       pow(pion_mass, 2) * (s * (5 - 16 * C4 * s) +
                                            pow(m_rho, 2) * (2 - 8 * C4 * s)) +
                       s * (s * (-3 + 4 * C4 * s) +
                            pow(m_rho, 2) * (-2 + 4 * C4 * s))))) *
         tmin) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (8 * (eta1 - eta2) *
         (delta *
              (eta1 * (4 * pow(pion_mass, 6) +
                       pow(pion_mass, 4) * (7 * pow(m_rho, 2) - 8 * s) +
                       pow(a1_mass, 4) * (pow(pion_mass, 2) - s) -
                       pow(a1_mass, 2) *
                           (2 * pow(pion_mass, 2) + pow(m_rho, 2) - 2 * s) *
                           (2 * pow(pion_mass, 2) - s) +
                       pow(pion_mass, 2) * s * (-8 * pow(m_rho, 2) + 5 * s) +
                       s * (pow(m_rho, 4) + pow(m_rho, 2) * s - pow(s, 2))) +
               eta2 *
                   (-4 * pow(pion_mass, 6) -
                    pow(pion_mass, 4) * (pow(m_rho, 2) - 8 * s) +
                    pow(a1_mass, 4) * (-pow(pion_mass, 2) + s) +
                    pow(pion_mass, 2) * (2 * pow(m_rho, 4) - 5 * pow(s, 2)) +
                    s * (-2 * pow(m_rho, 4) + pow(m_rho, 2) * s + pow(s, 2)) +
                    pow(a1_mass, 2) *
                        (4 * pow(pion_mass, 4) - 6 * pow(pion_mass, 2) * s +
                         s * (pow(m_rho, 2) + 2 * s)))) -
          2 * pow(m_rho, 2) *
              (eta1 *
                   (8 * C4 * pow(pion_mass, 6) +
                    pow(pion_mass, 4) * (3 + 8 * C4 * (pow(m_rho, 2) - 2 * s)) +
                    2 * C4 * pow(a1_mass, 4) * (pow(pion_mass, 2) - s) +
                    2 * pow(pion_mass, 2) * s *
                        (-1 - 6 * C4 * pow(m_rho, 2) + 5 * C4 * s) -
                    pow(a1_mass, 2) *
                        (8 * C4 * pow(pion_mass, 4) +
                         pow(pion_mass, 2) *
                             (1 + 2 * C4 * (pow(m_rho, 2) - 6 * s)) +
                         2 * C4 * s * (-pow(m_rho, 2) + 2 * s)) +
                    s * (-(s * (1 + 2 * C4 * s)) +
                         pow(m_rho, 2) * (1 + 4 * C4 * s))) +
               eta2 * (2 * C4 * pow(a1_mass, 4) * (-pow(pion_mass, 2) + s) -
                       (pow(pion_mass, 2) - s) *
                           (8 * C4 * pow(pion_mass, 4) - 2 * pow(m_rho, 2) + s +
                            2 * C4 * pow(s, 2) +
                            pow(pion_mass, 2) *
                                (3 - 4 * C4 * (pow(m_rho, 2) + 2 * s))) +
                       pow(a1_mass, 2) *
                           (8 * C4 * pow(pion_mass, 4) +
                            2 * C4 * s * (pow(m_rho, 2) + 2 * s) +
                            pow(pion_mass, 2) *
                                (1 - 2 * C4 * (pow(m_rho, 2) + 6 * s)))))) *
         tmin) /
            (pow(m_rho, 2) * (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (8 * (-2 + delta) * (delta - 4 * C4 * pow(m_rho, 2)) * pow(tmin, 2)) /
            (pow(m_rho, 2) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (16 * (-2 + delta) * (delta - 4 * C4 * pow(m_rho, 2)) * s *
         pow(tmin, 2)) /
            (pow(m_rho, 2) * (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) * (pow(m_rho, 2) - s) + 2 * eta1 * eta2 * s -
          pow(eta2, 2) * s) *
         (pow(pion_mass, 4) - pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) +
          (pow(m_rho, 2) - s) * s) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (-(eta1 * (pow(pion_mass, 2) + 2 * pow(m_rho, 2) - 3 * s)) -
          eta2 * (pow(pion_mass, 2) + s)) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (pow(pion_mass, 2) + s) *
         (-2 * eta2 * s + eta1 * (pow(pion_mass, 2) - pow(m_rho, 2) + s)) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (pow(eta1 - eta2, 3) *
         (-(eta1 *
            (pow(a1_mass, 2) - 2 * pow(pion_mass, 2) - pow(m_rho, 2) + s)) +
          eta2 *
              (pow(a1_mass, 2) - 2 * pow(pion_mass, 2) + pow(m_rho, 2) + s)) *
         pow(tmin, 2)) /
            (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
             2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) +
        (8 * (delta - 4 * C4 * pow(m_rho, 2)) *
         (delta * (4 * pow(pion_mass, 2) + 3 * pow(m_rho, 2) - 2 * s) -
          2 * pow(m_rho, 2) * (3 + 8 * C4 * pow(pion_mass, 2) - 4 * C4 * s)) *
         pow(tmin, 2)) /
            (pow(m_rho, 4) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (pow(eta1 - eta2, 2) * (pow(a1_mass, 2) - s) *
         (pow(eta2, 2) * s *
              (pow(a1_mass, 2) - 4 * pow(pion_mass, 2) + pow(m_rho, 2) +
               3 * s) +
          pow(eta1, 2) * (2 * pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) +
                          s * (pow(a1_mass, 2) - 3 * pow(m_rho, 2) + 3 * s)) -
          2 * eta1 * eta2 *
              (pow(pion_mass, 2) * (pow(m_rho, 2) - 4 * s) +
               s * (pow(a1_mass, 2) - 2 * pow(m_rho, 2) + 3 * s))) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (2 * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (eta1 * (4 * delta * pow(s, 2) -
                  2 * pow(m_rho, 2) * s * (-2 + 3 * delta + 8 * C4 * s) +
                  pow(m_rho, 4) * (-2 + delta + 16 * C4 * s) -
                  2 * pow(pion_mass, 2) *
                      (8 * C4 * pow(m_rho, 4) + 4 * delta * s +
                       pow(m_rho, 2) * (2 - 3 * delta - 16 * C4 * s))) +
          eta2 * (pow(pion_mass, 2) *
                      (8 * delta * s +
                       pow(m_rho, 2) * (-2 + delta - 32 * C4 * s)) +
                  s * (-4 * delta * s +
                       pow(m_rho, 2) * (-2 + delta + 16 * C4 * s)))) *
         pow(tmin, 2)) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (4 * (eta1 - eta2) *
         (delta * (eta1 * (pow(a1_mass, 2) * (pow(pion_mass, 2) - s) -
                           (2 * pow(pion_mass, 2) + pow(m_rho, 2) - 2 * s) *
                               (2 * pow(pion_mass, 2) - s)) +
                   eta2 * (4 * pow(pion_mass, 4) - 6 * pow(pion_mass, 2) * s +
                           pow(a1_mass, 2) * (-pow(pion_mass, 2) + s) +
                           s * (pow(m_rho, 2) + 2 * s))) +
          2 * pow(m_rho, 2) *
              (eta1 * (8 * C4 * pow(pion_mass, 4) +
                       2 * C4 * s * (pow(a1_mass, 2) - pow(m_rho, 2) + 2 * s) +
                       pow(pion_mass, 2) * (1 - 2 * C4 *
                                                    (pow(a1_mass, 2) -
                                                     pow(m_rho, 2) + 6 * s))) -
               eta2 * (8 * C4 * pow(pion_mass, 4) +
                       2 * C4 * s * (pow(a1_mass, 2) + pow(m_rho, 2) + 2 * s) -
                       pow(pion_mass, 2) *
                           (-1 +
                            2 * C4 *
                                (pow(a1_mass, 2) + pow(m_rho, 2) + 6 * s))))) *
         pow(tmin, 2)) /
            (pow(m_rho, 2) * (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (pow(eta1 - eta2, 4) * pow(tmin, 3)) /
            (3. * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                   2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (8 * pow(eta1 - eta2, 2) * (delta - 4 * C4 * pow(m_rho, 2)) *
         pow(tmin, 3)) /
            (3. * pow(m_rho, 2) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (16 * pow(delta - 4 * C4 * pow(m_rho, 2), 2) * pow(tmin, 3)) /
            (3. * pow(m_rho, 4) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (4 * (-2 + delta) * eta1 * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         pow(tmin, 3)) /
            (3. *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (2 * pow(eta1 - eta2, 4) * (pow(a1_mass, 2) - s) * s * pow(tmin, 3)) /
            (3. *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) * s *
         (-2 * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-pow(m_rho, 2) + s)) *
         pow(tmin, 3)) /
            (3. *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (4 * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (2 * eta2 * (delta - 4 * C4 * pow(m_rho, 2)) * s +
          eta1 * (-2 * delta * s + pow(m_rho, 2) * (-2 + delta + 8 * C4 * s))) *
         pow(tmin, 3)) /
            (3. * pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) *
         (pow(eta1, 2) * (2 * pow(a1_mass, 6) -
                          3 * pow(a1_mass, 4) *
                              (2 * pow(pion_mass, 2) + pow(m_rho, 2) - s) +
                          pow(m_rho, 2) * (pow(m_rho, 2) - s) * s -
                          pow(pion_mass, 4) * (3 * pow(m_rho, 2) + s) +
                          pow(pion_mass, 2) *
                              (-2 * pow(m_rho, 4) + 3 * pow(m_rho, 2) * s) +
                          pow(a1_mass, 2) *
                              (4 * pow(pion_mass, 4) + pow(m_rho, 4) +
                               pow(pion_mass, 2) * (8 * pow(m_rho, 2) - 4 * s) -
                               4 * pow(m_rho, 2) * s + 2 * pow(s, 2))) -
          2 * eta1 * eta2 *
              (2 * pow(a1_mass, 6) -
               pow(pion_mass, 2) * (pow(pion_mass, 2) - pow(m_rho, 2)) *
                   (pow(m_rho, 2) + s) +
               pow(a1_mass, 4) * (-6 * pow(pion_mass, 2) + 3 * s) +
               pow(a1_mass, 2) *
                   (4 * pow(pion_mass, 4) - pow(m_rho, 4) +
                    2 * pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) -
                    2 * pow(m_rho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) * (2 * pow(a1_mass, 6) +
                          3 * pow(a1_mass, 4) *
                              (-2 * pow(pion_mass, 2) + pow(m_rho, 2) + s) +
                          pow(pion_mass, 2) *
                              (-pow(m_rho, 4) +
                               pow(pion_mass, 2) * (2 * pow(m_rho, 2) - s) +
                               pow(m_rho, 2) * s) +
                          pow(a1_mass, 2) *
                              (4 * pow(pion_mass, 4) + pow(m_rho, 4) -
                               2 * pow(m_rho, 2) * s + 2 * pow(s, 2) -
                               4 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)))) *
         log(abs(-pow(a1_mass, 2) + tmax))) /
            (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
             2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) -
        (2 * pow(eta1 - eta2, 2) * (pow(a1_mass, 2) - s) *
         (-2 * eta1 * eta2 *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
               2 * pow(a1_mass, 2) * pow(pion_mass, 4) * s +
               pow(pion_mass, 2) *
                   (pow(a1_mass, 4) * (pow(m_rho, 2) - 4 * s) +
                    4 * pow(a1_mass, 2) * (pow(m_rho, 2) - s) * s +
                    pow(m_rho, 2) * pow(s, 2)) +
               pow(a1_mass, 2) * s *
                   (pow(a1_mass, 4) + s * (-2 * pow(m_rho, 2) + s) +
                    pow(a1_mass, 2) * (-2 * pow(m_rho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (pow(pion_mass, 8) -
               4 * pow(a1_mass, 2) * pow(pion_mass, 2) * s *
                   (pow(a1_mass, 2) + pow(m_rho, 2) + s) +
               pow(pion_mass, 4) * (pow(m_rho, 2) * s +
                                    pow(a1_mass, 2) * (pow(m_rho, 2) + 2 * s)) +
               pow(a1_mass, 2) * s *
                   (pow(a1_mass, 4) + s * (pow(m_rho, 2) + s) +
                    pow(a1_mass, 2) * (pow(m_rho, 2) + 3 * s))) +
          pow(eta1, 2) *
              (pow(pion_mass, 8) +
               pow(a1_mass, 2) * s *
                   (pow(a1_mass, 4) + 2 * pow(m_rho, 4) -
                    3 * pow(a1_mass, 2) * (pow(m_rho, 2) - s) -
                    3 * pow(m_rho, 2) * s + pow(s, 2)) +
               pow(pion_mass, 4) *
                   (2 * pow(m_rho, 4) - 3 * pow(m_rho, 2) * s +
                    pow(a1_mass, 2) * (-3 * pow(m_rho, 2) + 2 * s)) +
               2 * pow(pion_mass, 2) *
                   (pow(a1_mass, 4) * (pow(m_rho, 2) - 2 * s) +
                    pow(m_rho, 2) * s * (-pow(m_rho, 2) + s) -
                    pow(a1_mass, 2) * (pow(m_rho, 4) - 4 * pow(m_rho, 2) * s +
                                       2 * pow(s, 2))))) *
         log(abs(-pow(a1_mass, 2) + tmax))) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (8 * (eta1 - eta2) *
         (delta *
              (eta2 * (pow(pion_mass, 6) * pow(m_rho, 2) *
                           (2 * pow(pion_mass, 2) - s) +
                       pow(a1_mass, 8) * (-pow(pion_mass, 2) + s) +
                       pow(a1_mass, 6) *
                           (5 * pow(pion_mass, 4) - 7 * pow(pion_mass, 2) * s +
                            s * (pow(m_rho, 2) + 2 * s)) +
                       pow(a1_mass, 4) *
                           (-8 * pow(pion_mass, 6) -
                            pow(pion_mass, 4) * (pow(m_rho, 2) - 14 * s) +
                            pow(pion_mass, 2) *
                                (2 * pow(m_rho, 4) - pow(m_rho, 2) * s -
                                 7 * pow(s, 2)) +
                            s * (-2 * pow(m_rho, 4) + pow(m_rho, 2) * s +
                                 pow(s, 2))) +
                       pow(a1_mass, 2) * pow(pion_mass, 2) *
                           (4 * pow(pion_mass, 6) +
                            pow(pion_mass, 4) * (pow(m_rho, 2) - 8 * s) +
                            s * (2 * pow(m_rho, 4) + pow(m_rho, 2) * s -
                                 pow(s, 2)) +
                            pow(pion_mass, 2) *
                                (-2 * pow(m_rho, 4) - 3 * pow(m_rho, 2) * s +
                                 5 * pow(s, 2)))) +
               eta1 *
                   (pow(a1_mass, 8) * (pow(pion_mass, 2) - s) +
                    pow(a1_mass, 6) *
                        (-5 * pow(pion_mass, 4) + (pow(m_rho, 2) - 2 * s) * s +
                         pow(pion_mass, 2) * (-2 * pow(m_rho, 2) + 7 * s)) +
                    pow(pion_mass, 2) * pow(m_rho, 2) *
                        (2 * pow(pion_mass, 6) +
                         pow(pion_mass, 4) * (4 * pow(m_rho, 2) - 5 * s) +
                         pow(m_rho, 4) * s -
                         pow(pion_mass, 2) *
                             (pow(m_rho, 4) + 3 * pow(m_rho, 2) * s -
                              2 * pow(s, 2))) +
                    pow(a1_mass, 4) *
                        (8 * pow(pion_mass, 6) +
                         pow(pion_mass, 4) * (9 * pow(m_rho, 2) - 14 * s) +
                         pow(pion_mass, 2) * s * (-9 * pow(m_rho, 2) + 7 * s) +
                         s * (pow(m_rho, 4) + pow(m_rho, 2) * s - pow(s, 2))) +
                    pow(a1_mass, 2) *
                        (-4 * pow(pion_mass, 8) +
                         pow(m_rho, 4) * s * (-pow(m_rho, 2) + s) +
                         pow(pion_mass, 6) * (-11 * pow(m_rho, 2) + 8 * s) +
                         pow(pion_mass, 4) *
                             (-3 * pow(m_rho, 4) + 17 * pow(m_rho, 2) * s -
                              5 * pow(s, 2)) +
                         pow(pion_mass, 2) *
                             (pow(m_rho, 6) - 5 * pow(m_rho, 2) * pow(s, 2) +
                              pow(s, 3))))) -
          2 * pow(m_rho, 2) *
              (eta2 *
                   (pow(pion_mass, 8) * (1 + 2 * C4 * pow(m_rho, 2)) -
                    2 * C4 * pow(pion_mass, 6) * pow(m_rho, 2) * s +
                    2 * C4 * pow(a1_mass, 8) * (-pow(pion_mass, 2) + s) +
                    pow(a1_mass, 4) *
                        (-16 * C4 * pow(pion_mass, 6) +
                         pow(pion_mass, 4) *
                             (-4 + 6 * C4 * pow(m_rho, 2) + 28 * C4 * s) +
                         2 * pow(pion_mass, 2) *
                             (pow(m_rho, 2) + s - 3 * C4 * pow(m_rho, 2) * s -
                              7 * C4 * pow(s, 2)) +
                         s * (-2 * pow(m_rho, 2) + s + 2 * C4 * pow(s, 2))) +
                    pow(a1_mass, 6) *
                        (10 * C4 * pow(pion_mass, 4) +
                         2 * C4 * s * (pow(m_rho, 2) + 2 * s) +
                         pow(pion_mass, 2) *
                             (1 - 2 * C4 * (pow(m_rho, 2) + 7 * s))) +
                    pow(a1_mass, 2) * pow(pion_mass, 2) *
                        (8 * C4 * pow(pion_mass, 6) -
                         2 * pow(pion_mass, 4) *
                             (-2 + 3 * C4 * pow(m_rho, 2) + 8 * C4 * s) +
                         s * (2 * pow(m_rho, 2) + s - 2 * C4 * pow(s, 2)) +
                         2 * pow(pion_mass, 2) *
                             (pow(m_rho, 2) * (-1 + 3 * C4 * s) +
                              s * (-3 + 5 * C4 * s)))) +
               eta1 *
                   (pow(pion_mass, 8) * (-1 + 6 * C4 * pow(m_rho, 2)) +
                    2 * C4 * pow(a1_mass, 8) * (pow(pion_mass, 2) - s) +
                    pow(pion_mass, 2) * pow(m_rho, 4) * s +
                    2 * pow(pion_mass, 6) * pow(m_rho, 2) * (2 - 5 * C4 * s) -
                    pow(a1_mass, 6) *
                        (10 * C4 * pow(pion_mass, 4) +
                         pow(pion_mass, 2) *
                             (1 + 2 * C4 * (pow(m_rho, 2) - 7 * s)) +
                         2 * C4 * s * (-pow(m_rho, 2) + 2 * s)) -
                    pow(pion_mass, 4) * pow(m_rho, 2) *
                        (pow(m_rho, 2) + s * (3 - 4 * C4 * s)) +
                    pow(a1_mass, 4) *
                        (16 * C4 * pow(pion_mass, 6) +
                         2 * pow(pion_mass, 4) *
                             (2 + 5 * C4 * pow(m_rho, 2) - 14 * C4 * s) +
                         2 * pow(pion_mass, 2) * s *
                             (-1 - 7 * C4 * pow(m_rho, 2) + 7 * C4 * s) +
                         s * (-(s * (1 + 2 * C4 * s)) +
                              pow(m_rho, 2) * (1 + 4 * C4 * s))) -
                    pow(a1_mass, 2) *
                        (8 * C4 * pow(pion_mass, 8) +
                         pow(m_rho, 2) * (pow(m_rho, 2) - s) * s +
                         2 * pow(pion_mass, 6) *
                             (2 + 7 * C4 * pow(m_rho, 2) - 8 * C4 * s) +
                         pow(pion_mass, 2) *
                             (-pow(m_rho, 4) + pow(s, 2) +
                              8 * C4 * pow(m_rho, 2) * pow(s, 2) -
                              2 * C4 * pow(s, 3)) +
                         pow(pion_mass, 4) *
                             (pow(m_rho, 2) * (3 - 22 * C4 * s) +
                              2 * s * (-3 + 5 * C4 * s)))))) *
         log(abs(-pow(a1_mass, 2) + tmax))) /
            ((pow(a1_mass, 2) - pow(pion_mass, 2)) * pow(m_rho, 2) *
             (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (16 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
         log(abs(-pow(pion_mass, 2) + tmax))) /
            (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
             2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) -
        (8 * pow(-2 + delta, 2) *
         (3 * pow(pion_mass, 4) - 4 * pow(pion_mass, 2) * (pow(m_rho, 2) - s) +
          pow(pow(m_rho, 2) - s, 2)) *
         log(abs(-pow(pion_mass, 2) + tmax))) /
            ((pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (8 * (-2 + delta) * (eta1 - eta2) * pow(pion_mass, 2) *
         (2 * eta1 * pow(pion_mass, 2) - 2 * eta2 * pow(pion_mass, 2) -
          eta1 * pow(m_rho, 2)) *
         (pow(pion_mass, 2) - s) * log(abs(-pow(pion_mass, 2) + tmax))) /
            ((pow(a1_mass, 2) - pow(pion_mass, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (8 * (-2 + delta) * (eta1 - eta2) * pow(pion_mass, 2) *
         (pow(a1_mass, 2) - s) * (pow(pion_mass, 2) - s) *
         (-(eta2 * (pow(pion_mass, 2) + s)) +
          eta1 * (pow(pion_mass, 2) - pow(m_rho, 2) + s)) *
         log(abs(-pow(pion_mass, 2) + tmax))) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (8 * (-2 + delta) *
         (-(delta * (4 * pow(pion_mass, 2) - pow(m_rho, 2)) *
            (pow(pion_mass, 2) + pow(m_rho, 2) - s)) +
          2 * pow(m_rho, 2) *
              (8 * C4 * pow(pion_mass, 4) - pow(m_rho, 2) + s +
               pow(pion_mass, 2) * (3 - 8 * C4 * s))) *
         log(abs(-pow(pion_mass, 2) + tmax))) /
            (pow(m_rho, 2) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) *
         (pow(eta1, 2) * (2 * pow(a1_mass, 6) -
                          3 * pow(a1_mass, 4) *
                              (2 * pow(pion_mass, 2) + pow(m_rho, 2) - s) +
                          pow(m_rho, 2) * (pow(m_rho, 2) - s) * s -
                          pow(pion_mass, 4) * (3 * pow(m_rho, 2) + s) +
                          pow(pion_mass, 2) *
                              (-2 * pow(m_rho, 4) + 3 * pow(m_rho, 2) * s) +
                          pow(a1_mass, 2) *
                              (4 * pow(pion_mass, 4) + pow(m_rho, 4) +
                               pow(pion_mass, 2) * (8 * pow(m_rho, 2) - 4 * s) -
                               4 * pow(m_rho, 2) * s + 2 * pow(s, 2))) -
          2 * eta1 * eta2 *
              (2 * pow(a1_mass, 6) -
               pow(pion_mass, 2) * (pow(pion_mass, 2) - pow(m_rho, 2)) *
                   (pow(m_rho, 2) + s) +
               pow(a1_mass, 4) * (-6 * pow(pion_mass, 2) + 3 * s) +
               pow(a1_mass, 2) *
                   (4 * pow(pion_mass, 4) - pow(m_rho, 4) +
                    2 * pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) -
                    2 * pow(m_rho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) * (2 * pow(a1_mass, 6) +
                          3 * pow(a1_mass, 4) *
                              (-2 * pow(pion_mass, 2) + pow(m_rho, 2) + s) +
                          pow(pion_mass, 2) *
                              (-pow(m_rho, 4) +
                               pow(pion_mass, 2) * (2 * pow(m_rho, 2) - s) +
                               pow(m_rho, 2) * s) +
                          pow(a1_mass, 2) *
                              (4 * pow(pion_mass, 4) + pow(m_rho, 4) -
                               2 * pow(m_rho, 2) * s + 2 * pow(s, 2) -
                               4 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)))) *
         log(abs(-pow(a1_mass, 2) + tmin))) /
            (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
             2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) +
        (2 * pow(eta1 - eta2, 2) * (pow(a1_mass, 2) - s) *
         (-2 * eta1 * eta2 *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
               2 * pow(a1_mass, 2) * pow(pion_mass, 4) * s +
               pow(pion_mass, 2) *
                   (pow(a1_mass, 4) * (pow(m_rho, 2) - 4 * s) +
                    4 * pow(a1_mass, 2) * (pow(m_rho, 2) - s) * s +
                    pow(m_rho, 2) * pow(s, 2)) +
               pow(a1_mass, 2) * s *
                   (pow(a1_mass, 4) + s * (-2 * pow(m_rho, 2) + s) +
                    pow(a1_mass, 2) * (-2 * pow(m_rho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (pow(pion_mass, 8) -
               4 * pow(a1_mass, 2) * pow(pion_mass, 2) * s *
                   (pow(a1_mass, 2) + pow(m_rho, 2) + s) +
               pow(pion_mass, 4) * (pow(m_rho, 2) * s +
                                    pow(a1_mass, 2) * (pow(m_rho, 2) + 2 * s)) +
               pow(a1_mass, 2) * s *
                   (pow(a1_mass, 4) + s * (pow(m_rho, 2) + s) +
                    pow(a1_mass, 2) * (pow(m_rho, 2) + 3 * s))) +
          pow(eta1, 2) *
              (pow(pion_mass, 8) +
               pow(a1_mass, 2) * s *
                   (pow(a1_mass, 4) + 2 * pow(m_rho, 4) -
                    3 * pow(a1_mass, 2) * (pow(m_rho, 2) - s) -
                    3 * pow(m_rho, 2) * s + pow(s, 2)) +
               pow(pion_mass, 4) *
                   (2 * pow(m_rho, 4) - 3 * pow(m_rho, 2) * s +
                    pow(a1_mass, 2) * (-3 * pow(m_rho, 2) + 2 * s)) +
               2 * pow(pion_mass, 2) *
                   (pow(a1_mass, 4) * (pow(m_rho, 2) - 2 * s) +
                    pow(m_rho, 2) * s * (-pow(m_rho, 2) + s) -
                    pow(a1_mass, 2) * (pow(m_rho, 4) - 4 * pow(m_rho, 2) * s +
                                       2 * pow(s, 2))))) *
         log(abs(-pow(a1_mass, 2) + tmin))) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (8 * (eta1 - eta2) *
         (delta *
              (eta2 * (pow(pion_mass, 6) * pow(m_rho, 2) *
                           (2 * pow(pion_mass, 2) - s) +
                       pow(a1_mass, 8) * (-pow(pion_mass, 2) + s) +
                       pow(a1_mass, 6) *
                           (5 * pow(pion_mass, 4) - 7 * pow(pion_mass, 2) * s +
                            s * (pow(m_rho, 2) + 2 * s)) +
                       pow(a1_mass, 4) *
                           (-8 * pow(pion_mass, 6) -
                            pow(pion_mass, 4) * (pow(m_rho, 2) - 14 * s) +
                            pow(pion_mass, 2) *
                                (2 * pow(m_rho, 4) - pow(m_rho, 2) * s -
                                 7 * pow(s, 2)) +
                            s * (-2 * pow(m_rho, 4) + pow(m_rho, 2) * s +
                                 pow(s, 2))) +
                       pow(a1_mass, 2) * pow(pion_mass, 2) *
                           (4 * pow(pion_mass, 6) +
                            pow(pion_mass, 4) * (pow(m_rho, 2) - 8 * s) +
                            s * (2 * pow(m_rho, 4) + pow(m_rho, 2) * s -
                                 pow(s, 2)) +
                            pow(pion_mass, 2) *
                                (-2 * pow(m_rho, 4) - 3 * pow(m_rho, 2) * s +
                                 5 * pow(s, 2)))) +
               eta1 *
                   (pow(a1_mass, 8) * (pow(pion_mass, 2) - s) +
                    pow(a1_mass, 6) *
                        (-5 * pow(pion_mass, 4) + (pow(m_rho, 2) - 2 * s) * s +
                         pow(pion_mass, 2) * (-2 * pow(m_rho, 2) + 7 * s)) +
                    pow(pion_mass, 2) * pow(m_rho, 2) *
                        (2 * pow(pion_mass, 6) +
                         pow(pion_mass, 4) * (4 * pow(m_rho, 2) - 5 * s) +
                         pow(m_rho, 4) * s -
                         pow(pion_mass, 2) *
                             (pow(m_rho, 4) + 3 * pow(m_rho, 2) * s -
                              2 * pow(s, 2))) +
                    pow(a1_mass, 4) *
                        (8 * pow(pion_mass, 6) +
                         pow(pion_mass, 4) * (9 * pow(m_rho, 2) - 14 * s) +
                         pow(pion_mass, 2) * s * (-9 * pow(m_rho, 2) + 7 * s) +
                         s * (pow(m_rho, 4) + pow(m_rho, 2) * s - pow(s, 2))) +
                    pow(a1_mass, 2) *
                        (-4 * pow(pion_mass, 8) +
                         pow(m_rho, 4) * s * (-pow(m_rho, 2) + s) +
                         pow(pion_mass, 6) * (-11 * pow(m_rho, 2) + 8 * s) +
                         pow(pion_mass, 4) *
                             (-3 * pow(m_rho, 4) + 17 * pow(m_rho, 2) * s -
                              5 * pow(s, 2)) +
                         pow(pion_mass, 2) *
                             (pow(m_rho, 6) - 5 * pow(m_rho, 2) * pow(s, 2) +
                              pow(s, 3))))) -
          2 * pow(m_rho, 2) *
              (eta2 *
                   (pow(pion_mass, 8) * (1 + 2 * C4 * pow(m_rho, 2)) -
                    2 * C4 * pow(pion_mass, 6) * pow(m_rho, 2) * s +
                    2 * C4 * pow(a1_mass, 8) * (-pow(pion_mass, 2) + s) +
                    pow(a1_mass, 4) *
                        (-16 * C4 * pow(pion_mass, 6) +
                         pow(pion_mass, 4) *
                             (-4 + 6 * C4 * pow(m_rho, 2) + 28 * C4 * s) +
                         2 * pow(pion_mass, 2) *
                             (pow(m_rho, 2) + s - 3 * C4 * pow(m_rho, 2) * s -
                              7 * C4 * pow(s, 2)) +
                         s * (-2 * pow(m_rho, 2) + s + 2 * C4 * pow(s, 2))) +
                    pow(a1_mass, 6) *
                        (10 * C4 * pow(pion_mass, 4) +
                         2 * C4 * s * (pow(m_rho, 2) + 2 * s) +
                         pow(pion_mass, 2) *
                             (1 - 2 * C4 * (pow(m_rho, 2) + 7 * s))) +
                    pow(a1_mass, 2) * pow(pion_mass, 2) *
                        (8 * C4 * pow(pion_mass, 6) -
                         2 * pow(pion_mass, 4) *
                             (-2 + 3 * C4 * pow(m_rho, 2) + 8 * C4 * s) +
                         s * (2 * pow(m_rho, 2) + s - 2 * C4 * pow(s, 2)) +
                         2 * pow(pion_mass, 2) *
                             (pow(m_rho, 2) * (-1 + 3 * C4 * s) +
                              s * (-3 + 5 * C4 * s)))) +
               eta1 *
                   (pow(pion_mass, 8) * (-1 + 6 * C4 * pow(m_rho, 2)) +
                    2 * C4 * pow(a1_mass, 8) * (pow(pion_mass, 2) - s) +
                    pow(pion_mass, 2) * pow(m_rho, 4) * s +
                    2 * pow(pion_mass, 6) * pow(m_rho, 2) * (2 - 5 * C4 * s) -
                    pow(a1_mass, 6) *
                        (10 * C4 * pow(pion_mass, 4) +
                         pow(pion_mass, 2) *
                             (1 + 2 * C4 * (pow(m_rho, 2) - 7 * s)) +
                         2 * C4 * s * (-pow(m_rho, 2) + 2 * s)) -
                    pow(pion_mass, 4) * pow(m_rho, 2) *
                        (pow(m_rho, 2) + s * (3 - 4 * C4 * s)) +
                    pow(a1_mass, 4) *
                        (16 * C4 * pow(pion_mass, 6) +
                         2 * pow(pion_mass, 4) *
                             (2 + 5 * C4 * pow(m_rho, 2) - 14 * C4 * s) +
                         2 * pow(pion_mass, 2) * s *
                             (-1 - 7 * C4 * pow(m_rho, 2) + 7 * C4 * s) +
                         s * (-(s * (1 + 2 * C4 * s)) +
                              pow(m_rho, 2) * (1 + 4 * C4 * s))) -
                    pow(a1_mass, 2) *
                        (8 * C4 * pow(pion_mass, 8) +
                         pow(m_rho, 2) * (pow(m_rho, 2) - s) * s +
                         2 * pow(pion_mass, 6) *
                             (2 + 7 * C4 * pow(m_rho, 2) - 8 * C4 * s) +
                         pow(pion_mass, 2) *
                             (-pow(m_rho, 4) + pow(s, 2) +
                              8 * C4 * pow(m_rho, 2) * pow(s, 2) -
                              2 * C4 * pow(s, 3)) +
                         pow(pion_mass, 4) *
                             (pow(m_rho, 2) * (3 - 22 * C4 * s) +
                              2 * s * (-3 + 5 * C4 * s)))))) *
         log(abs(-pow(a1_mass, 2) + tmin))) /
            ((pow(a1_mass, 2) - pow(pion_mass, 2)) * pow(m_rho, 2) *
             (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (16 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
         log(abs(-pow(pion_mass, 2) + tmin))) /
            (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
             2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) +
        (8 * pow(-2 + delta, 2) *
         (3 * pow(pion_mass, 4) - 4 * pow(pion_mass, 2) * (pow(m_rho, 2) - s) +
          pow(pow(m_rho, 2) - s, 2)) *
         log(abs(-pow(pion_mass, 2) + tmin))) /
            ((pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (8 * (-2 + delta) * (eta1 - eta2) * pow(pion_mass, 2) *
         (2 * eta1 * pow(pion_mass, 2) - 2 * eta2 * pow(pion_mass, 2) -
          eta1 * pow(m_rho, 2)) *
         (pow(pion_mass, 2) - s) * log(abs(-pow(pion_mass, 2) + tmin))) /
            ((pow(a1_mass, 2) - pow(pion_mass, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
        (8 * (-2 + delta) * (eta1 - eta2) * pow(pion_mass, 2) *
         (pow(a1_mass, 2) - s) * (pow(pion_mass, 2) - s) *
         (-(eta2 * (pow(pion_mass, 2) + s)) +
          eta1 * (pow(pion_mass, 2) - pow(m_rho, 2) + s)) *
         log(abs(-pow(pion_mass, 2) + tmin))) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (8 * (-2 + delta) *
         (delta * (4 * pow(pion_mass, 2) - pow(m_rho, 2)) *
              (pow(pion_mass, 2) + pow(m_rho, 2) - s) -
          2 * pow(m_rho, 2) *
              (8 * C4 * pow(pion_mass, 4) - pow(m_rho, 2) + s +
               pow(pion_mass, 2) * (3 - 8 * C4 * s))) *
         log(abs(-pow(pion_mass, 2) + tmin))) /
            (pow(m_rho, 2) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))))) /
      (512. * Pi);
  // clang-format on

  return cut_off(xs * gev2_mb / spin_deg_factor);
}

double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi_rho0_pi(
    const double s, const double t, const double m_rho) {
  const double spin_deg_factor = 3.0;

  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  // clang-format off
  const double diff_xs =
      (pow(Const, 2) * pow(ghat, 4) *
       ((-8 * pow(-2 + delta, 2) * pow(pion_mass, 2)) /
            (pow(m_rho, 2) * pow(pow(pion_mass, 2) - s, 2)) -
        (8 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
         (pow(pion_mass, 4) + pow(pow(m_rho, 2) - t, 2) -
          2 * pow(pion_mass, 2) * (pow(m_rho, 2) + t))) /
            (pow(m_rho, 2) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
             pow(pow(pion_mass, 2) - t, 2)) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (-(eta2 * (pow(pion_mass, 2) + s)) + eta1 * (-pow(m_rho, 2) + s + t)) *
         (-pow(pion_mass, 4) + pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * t) +
          t * (-pow(m_rho, 2) + 2 * s + t))) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
             (pow(pion_mass, 2) - t)) -
        (8 * (-2 + delta) *
         (pow(pion_mass, 4) * (2 - 3 * delta + 8 * C4 * pow(m_rho, 2)) +
          pow(m_rho, 4) * (-2 + delta + 8 * C4 * t) +
          t * ((2 + 3 * delta) * s + 2 * delta * t) +
          pow(pion_mass, 2) *
              (-8 * C4 * pow(m_rho, 4) + (-2 + delta) * s -
               (2 + 3 * delta) * t + 4 * pow(m_rho, 2) * (1 + 4 * C4 * t)) -
          pow(m_rho, 2) * (t * (-2 + 3 * delta + 8 * C4 * t) +
                           s * (-2 + delta + 16 * C4 * t)))) /
            (pow(m_rho, 2) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
             (pow(pion_mass, 2) - t)) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(a1_mass, 2) - s) *
         (eta2 * (pow(pion_mass, 2) + s) *
              (pow(pion_mass, 4) - pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) +
               s * (pow(m_rho, 2) - s - 2 * t)) +
          eta1 * (-4 * pow(pion_mass, 6) +
                  s * (-pow(m_rho, 2) + s) * (-pow(m_rho, 2) + s + t) +
                  pow(pion_mass, 4) * (3 * pow(m_rho, 2) + s + t) -
                  pow(pion_mass, 2) * (pow(m_rho, 4) + 2 * s * (s - t) +
                                       pow(m_rho, 2) * (-s + t))))) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta2, 2) *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(pion_mass, 4) * (pow(pow(m_rho, 2) + 2 * s, 2) - 2 * s * t) +
               pow(s, 2) * (pow(pow(m_rho, 2) + s, 2) +
                            2 * (-pow(m_rho, 2) + s) * t + 2 * pow(t, 2)) -
               2 * pow(pion_mass, 2) * s *
                   (pow(m_rho, 4) + pow(m_rho, 2) * (2 * s - t) +
                    2 * s * (s + t))) +
          2 * eta1 * eta2 *
              (-pow(pion_mass, 8) +
               pow(pion_mass, 4) * (pow(m_rho, 4) + 2 * pow(m_rho, 2) * s +
                                    2 * s * (-2 * s + t)) -
               2 * pow(pion_mass, 2) * s *
                   (pow(m_rho, 4) + pow(m_rho, 2) * (s + t) - 2 * s * (s + t)) +
               pow(s, 2) * (pow(m_rho, 4) - pow(s, 2) + 2 * pow(m_rho, 2) * t -
                            2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(pion_mass, 4) * (3 * pow(m_rho, 4) + 2 * s * (2 * s - t) +
                                    2 * pow(m_rho, 2) * (-3 * s + t)) -
               2 * pow(pion_mass, 2) * (pow(m_rho, 2) - s) *
                   (-2 * s * (s + t) + pow(m_rho, 2) * (2 * s + t)) +
               s * (-pow(m_rho, 2) + s) *
                   (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                    pow(m_rho, 2) * (s + 2 * t))))) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(pion_mass, 8) -
               pow(pion_mass, 4) *
                   (pow(m_rho, 4) + 2 * (pow(m_rho, 2) + s) * t -
                    4 * pow(t, 2)) +
               pow(t, 2) * (-pow(m_rho, 4) - 2 * pow(m_rho, 2) * s +
                            2 * pow(s, 2) + 2 * s * t + pow(t, 2)) +
               2 * pow(pion_mass, 2) * t *
                   (pow(m_rho, 4) + pow(m_rho, 2) * (s + t) -
                    2 * t * (s + t))) +
          pow(eta2, 2) *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(pion_mass, 4) * (pow(m_rho, 4) + 4 * pow(m_rho, 2) * t -
                                    2 * (s - 2 * t) * t) +
               pow(t, 2) * (pow(m_rho, 4) + 2 * pow(s, 2) + 2 * s * t +
                            pow(t, 2) + 2 * pow(m_rho, 2) * (-s + t)) -
               2 * pow(pion_mass, 2) * t *
                   (pow(m_rho, 4) - pow(m_rho, 2) * (s - 2 * t) +
                    2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(pion_mass, 4) *
                   (3 * pow(m_rho, 4) + 2 * pow(m_rho, 2) * (s - 3 * t) -
                    2 * (s - 2 * t) * t) +
               t * (-pow(m_rho, 2) + t) *
                   (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                    pow(m_rho, 2) * (2 * s + t)) -
               2 * pow(pion_mass, 2) * (-pow(m_rho, 2) + t) *
                   (2 * t * (s + t) - pow(m_rho, 2) * (s + 2 * t))))) /
            ((pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
             pow(pow(a1_mass, 2) - t, 2)) +
        (8 * (-2 + delta) *
         ((-2 + delta) * pow(m_rho, 6) +
          pow(pion_mass, 6) * (-2 + 3 * delta - 8 * C4 * pow(m_rho, 2)) +
          s * t * ((-2 + 3 * delta) * s + 4 * delta * t) +
          pow(pion_mass, 4) *
              (8 * C4 * pow(m_rho, 4) + 4 * delta * s + 2 * t - 3 * delta * t -
               pow(m_rho, 2) * (2 + delta + 16 * C4 * s - 8 * C4 * t)) +
          pow(m_rho, 4) *
              (-((-2 + delta) * t) + s * (4 - 2 * delta + 8 * C4 * t)) +
          pow(m_rho, 2) * s *
              (s * (-2 + delta - 8 * C4 * t) - 2 * t * (delta + 8 * C4 * t)) +
          pow(pion_mass, 2) *
              (s * ((2 - 3 * delta) * s - 8 * delta * t) -
               pow(m_rho, 4) * (-6 + 3 * delta + 8 * C4 * (s + t)) +
               pow(m_rho, 2) * (8 * C4 * pow(s, 2) + 4 * (-1 + delta) * t +
                                s * (-8 + 6 * delta + 32 * C4 * t))))) /
            (pow(m_rho, 2) * (pow(pion_mass, 2) - s) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
             (pow(pion_mass, 2) - t)) +
        (2 * pow(eta1 - eta2, 2) * (pow(a1_mass, 2) - s) *
         (pow(eta1, 2) * (pow(pion_mass, 8) +
                          pow(pion_mass, 4) * (2 * pow(m_rho, 4) + 2 * s * t -
                                               3 * pow(m_rho, 2) * (s + t)) +
                          s * t *
                              (2 * pow(m_rho, 4) + pow(s, 2) + 3 * s * t +
                               pow(t, 2) - 3 * pow(m_rho, 2) * (s + t)) -
                          2 * pow(pion_mass, 2) * (pow(m_rho, 2) - s - t) *
                              (-2 * s * t + pow(m_rho, 2) * (s + t))) +
          pow(eta2, 2) *
              (pow(pion_mass, 8) -
               4 * pow(pion_mass, 2) * s * t * (pow(m_rho, 2) + s + t) +
               pow(pion_mass, 4) * (2 * s * t + pow(m_rho, 2) * (s + t)) +
               s * t *
                   (pow(s, 2) + 3 * s * t + pow(t, 2) +
                    pow(m_rho, 2) * (s + t))) +
          2 * eta1 * eta2 *
              (-pow(pion_mass, 8) + 2 * pow(pion_mass, 6) * pow(m_rho, 2) -
               2 * pow(pion_mass, 4) * s * t -
               s * t *
                   (pow(s, 2) + 3 * s * t + pow(t, 2) -
                    2 * pow(m_rho, 2) * (s + t)) -
               pow(pion_mass, 2) *
                   (-4 * s * t * (s + t) +
                    pow(m_rho, 2) * (pow(s, 2) + 4 * s * t + pow(t, 2)))))) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
             (pow(a1_mass, 2) - t)) +
        (8 * (pow(delta, 2) *
                  (8 * pow(pion_mass, 4) + 3 * pow(m_rho, 4) -
                   6 * pow(m_rho, 2) * (s + t) + 2 * pow(s + t, 2) +
                   4 * pow(pion_mass, 2) * (3 * pow(m_rho, 2) - 2 * (s + t))) -
              4 * delta * pow(m_rho, 2) *
                  (16 * C4 * pow(pion_mass, 4) +
                   pow(m_rho, 2) * (3 - 6 * C4 * (s + t)) +
                   (s + t) * (-3 + 4 * C4 * (s + t)) +
                   2 * pow(pion_mass, 2) *
                       (3 + C4 * (6 * pow(m_rho, 2) - 8 * (s + t)))) +
              4 * pow(m_rho, 4) *
                  (3 + 4 * C4 * (2 * pow(pion_mass, 2) - s - t) *
                           (3 + C4 * (4 * pow(pion_mass, 2) - 2 * (s + t)))))) /
            (pow(m_rho, 4) * (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
                              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (4 * (eta1 - eta2) * (-pow(a1_mass, 2) + s) *
         (eta2 * (-2 * pow(pion_mass, 4) * (delta - 4 * C4 * pow(m_rho, 2)) *
                      (pow(m_rho, 2) + 4 * s) +
                  pow(pion_mass, 2) *
                      (-2 * pow(m_rho, 4) * (-2 + delta + 8 * C4 * s) +
                       8 * delta * s * (s + t) -
                       pow(m_rho, 2) * ((-10 + delta) * s - (-2 + delta) * t +
                                        32 * C4 * s * (s + t))) +
                  s * (2 * pow(m_rho, 4) * (-2 + delta + 4 * C4 * s) -
                       2 * delta * pow(s + t, 2) +
                       pow(m_rho, 2) * ((-6 + delta) * s + (-2 + delta) * t +
                                        8 * C4 * pow(s + t, 2)))) +
          eta1 * (4 * pow(pion_mass, 4) *
                      (6 * C4 * pow(m_rho, 4) + 2 * delta * s +
                       pow(m_rho, 2) * (1 - 2 * delta - 8 * C4 * s)) +
                  2 * delta * s * pow(s + t, 2) -
                  pow(m_rho, 2) *
                      ((-6 + 5 * delta) * pow(s, 2) +
                       2 * (-2 + 3 * delta) * s * t + (-2 + delta) * pow(t, 2) +
                       8 * C4 * s * pow(s + t, 2)) +
                  pow(m_rho, 4) *
                      ((-2 + delta) * (3 * s + t) + 8 * C4 * s * (s + 2 * t)) -
                  2 * pow(pion_mass, 2) *
                      (4 * delta * s * (s + t) -
                       pow(m_rho, 2) * (-6 * s + 7 * delta * s - 2 * t +
                                        3 * delta * t + 16 * C4 * s * (s + t)) +
                       2 * pow(m_rho, 4) *
                           (-2 + delta + 4 * C4 * (2 * s + t)))))) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
              2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) +
        (4 * (eta1 - eta2) *
         (((-2 + delta) *
           (pow(pion_mass, 4) - pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) +
            s * (pow(m_rho, 2) - s - 2 * t)) *
           (eta1 * (pow(m_rho, 2) - s - t) + eta2 * (pow(pion_mass, 2) + t))) /
              ((pow(pion_mass, 2) - s) * (pow(a1_mass, 2) - t)) +
          ((-2 + delta) *
           (eta2 * (pow(pion_mass, 2) + t) *
                (pow(pion_mass, 4) -
                 pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * t) +
                 (pow(m_rho, 2) - 2 * s - t) * t) +
            eta1 * (-4 * pow(pion_mass, 6) +
                    (pow(m_rho, 2) - t) * (pow(m_rho, 2) - s - t) * t +
                    pow(pion_mass, 4) * (3 * pow(m_rho, 2) + s + t) -
                    pow(pion_mass, 2) *
                        (pow(m_rho, 4) + pow(m_rho, 2) * (s - t) +
                         2 * t * (-s + t))))) /
              ((-pow(a1_mass, 2) + t) * (-pow(pion_mass, 2) + t)) +
          (eta2 *
               (-2 * pow(pion_mass, 4) * (delta - 4 * C4 * pow(m_rho, 2)) *
                    (pow(m_rho, 2) + 4 * t) +
                pow(pion_mass, 2) *
                    (8 * delta * t * (s + t) -
                     2 * pow(m_rho, 4) * (-2 + delta + 8 * C4 * t) -
                     pow(m_rho, 2) * (-((-2 + delta) * s) + (-10 + delta) * t +
                                      32 * C4 * t * (s + t))) +
                t * (-2 * delta * pow(s + t, 2) +
                     2 * pow(m_rho, 4) * (-2 + delta + 4 * C4 * t) +
                     pow(m_rho, 2) * ((-2 + delta) * s + (-6 + delta) * t +
                                      8 * C4 * pow(s + t, 2)))) +
           eta1 *
               (2 * delta * t * pow(s + t, 2) -
                pow(m_rho, 2) *
                    ((-2 + delta) * pow(s, 2) + 2 * (-2 + 3 * delta) * s * t +
                     (-6 + 5 * delta) * pow(t, 2) +
                     8 * C4 * t * pow(s + t, 2)) +
                pow(m_rho, 4) *
                    (8 * C4 * t * (2 * s + t) + (-2 + delta) * (s + 3 * t)) +
                4 * pow(pion_mass, 4) *
                    (6 * C4 * pow(m_rho, 4) + 2 * delta * t +
                     pow(m_rho, 2) * (1 - 2 * delta - 8 * C4 * t)) -
                2 * pow(pion_mass, 2) *
                    (4 * delta * t * (s + t) -
                     pow(m_rho, 2) * (-2 * s + 3 * delta * s - 6 * t +
                                      7 * delta * t + 16 * C4 * t * (s + t)) +
                     2 * pow(m_rho, 4) *
                         (-2 + delta + 4 * C4 * (s + 2 * t))))) /
              (pow(m_rho, 2) * (-pow(a1_mass, 2) + t)))) /
            (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
             2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)))) /
      (512. * Pi);

  // clang-format on
  return cut_off(gev2_mb * diff_xs / spin_deg_factor);
}

// C12
double
CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi0_rho_pi_rho_mediated(
    const double s, const double m_rho) {
  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  auto t_mandelstam = get_t_range(sqrt(s), pion_mass, m_rho, pion_mass, 0.);
  const double &tmin = t_mandelstam[1];
  const double &tmax = t_mandelstam[0];
  const double spin_deg_factor = 3.0;

  // clang-format off
  const double xs =
      (pow(Const, 2) * pow(ghat, 4) *
       (0. -
        (0.25 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
         (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
          2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
         tmax) /
            (pow(m_rho, 2) * pow(pow(pion_mass, 2) - s, 2)) -
        (0.0625 *
         (0. + 8. * pow(pion_mass, 2) * pow(m_rho, 4) - 12. * pow(m_rho, 6) +
          4. * pow(m_rho, 4) * s +
          delta * pow(m_rho, 2) *
              (-16. * pow(pion_mass, 4) -
               16. * pow(pion_mass, 2) * pow(m_rho, 2) - 4. * pow(m_rho, 4) +
               16. * pow(m_rho, 2) * s + 4. * pow(s, 2)) +
          pow(delta, 2) *
              (8. * pow(pion_mass, 6) + 9. * pow(m_rho, 6) +
               pow(pion_mass, 4) * (4. * pow(m_rho, 2) - 4. * s) -
               13. * pow(m_rho, 4) * s - 5. * pow(m_rho, 2) * pow(s, 2) +
               1. * pow(s, 3) +
               pow(pion_mass, 2) * (-2. * pow(m_rho, 4) +
                                    4. * pow(m_rho, 2) * s - 2. * pow(s, 2)))) *
         tmax) /
            pow(m_rho, 6) -
        (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
         (pow(a1_mass, 2) - 1. * s) *
         (eta2 * (pow(pion_mass, 6) + pow(pion_mass, 2) * pow(s, 2) +
                  (pow(m_rho, 2) - 1. * s) * pow(s, 2) +
                  pow(pion_mass, 4) * (-1. * pow(m_rho, 2) + 3. * s)) +
          eta1 * (-4. * pow(pion_mass, 6) +
                  pow(pion_mass, 4) * (3. * pow(m_rho, 2) + s) +
                  pow(pion_mass, 2) * (-1. * pow(m_rho, 4) + pow(m_rho, 2) * s -
                                       2. * pow(s, 2)) +
                  s * (pow(m_rho, 4) - 2. * pow(m_rho, 2) * s + pow(s, 2)))) *
         tmax) /
            ((-1. * pow(pion_mass, 2) + s) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta1, 2) *
              (1. * pow(pion_mass, 8) - 2. * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(pion_mass, 2) * s *
                   (-4. * pow(m_rho, 4) + 8. * pow(m_rho, 2) * s -
                    4. * pow(s, 2)) +
               pow(s, 2) * (1. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s +
                            1. * pow(s, 2)) +
               pow(pion_mass, 4) * (3. * pow(m_rho, 4) -
                                    6. * pow(m_rho, 2) * s + 4. * pow(s, 2))) +
          pow(eta2, 2) *
              (1. * pow(pion_mass, 8) - 2. * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(pion_mass, 2) * s *
                   (-2. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s -
                    4. * pow(s, 2)) +
               pow(s, 2) * (1. * pow(m_rho, 4) + 2. * pow(m_rho, 2) * s +
                            1. * pow(s, 2)) +
               pow(pion_mass, 4) * (1. * pow(m_rho, 4) +
                                    4. * pow(m_rho, 2) * s + 4. * pow(s, 2))) +
          eta1 * eta2 *
              (-2. * pow(pion_mass, 8) + 2. * pow(m_rho, 4) * pow(s, 2) -
               2. * pow(s, 4) +
               pow(pion_mass, 4) * (2. * pow(m_rho, 4) +
                                    4. * pow(m_rho, 2) * s - 8. * pow(s, 2)) +
               pow(pion_mass, 2) * s *
                   (-4. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s +
                    8. * pow(s, 2)))) *
         tmax) /
            (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
             2. * pow(a1_mass, 2) * s + pow(s, 2)) +
        (0.5 *
         (0. - 4. * C4 * pow(m_rho, 8) - 0.5 * pow(m_rho, 4) * s +
          pow(m_rho, 6) * (2. + 2. * C4 * s) +
          pow(delta, 2) *
              (2. * pow(pion_mass, 6) + 2. * pow(m_rho, 6) -
               1.5 * pow(pion_mass, 4) * s - 2.375 * pow(m_rho, 4) * s -
               0.75 * pow(m_rho, 2) * pow(s, 2) + 0.125 * pow(s, 3) +
               pow(pion_mass, 2) *
                   (-1.5 * pow(m_rho, 4) + 0.5 * pow(m_rho, 2) * s)) +
          delta * pow(m_rho, 2) *
              (-2. * pow(pion_mass, 4) +
               pow(pion_mass, 2) *
                   (1. * s + pow(m_rho, 2) * (-1. - 2. * C4 * s)) +
               pow(m_rho, 2) * (2. * C4 * pow(m_rho, 4) +
                                pow(m_rho, 2) * (-3. + 1. * C4 * s) +
                                s * (2. + 1. * C4 * s)))) *
         tmax) /
            pow(m_rho, 6) -
        (0.5 *
         (pow(m_rho, 6) *
              (-1.5 + C4 * (-12. * pow(pion_mass, 2) + 6. * s) +
               pow(C4, 2) * (-16. * pow(pion_mass, 4) +
                             16. * pow(pion_mass, 2) * s - 4. * pow(s, 2))) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 6) + 0.125 * pow(m_rho, 6) +
               pow(pion_mass, 4) * (-2. * pow(m_rho, 2) - 1. * s) +
               0.25 * pow(m_rho, 4) * s - 0.625 * pow(m_rho, 2) * pow(s, 2) +
               pow(pion_mass, 2) *
                   (-2.5 * pow(m_rho, 4) + 1.75 * pow(m_rho, 2) * s +
                    0.25 * pow(s, 2))) +
          delta * pow(m_rho, 2) *
              (pow(pion_mass, 4) * (1. + 8. * C4 * pow(m_rho, 2)) +
               pow(pion_mass, 2) * (6. * C4 * pow(m_rho, 4) - 0.5 * s +
                                    pow(m_rho, 2) * (3. - 10. * C4 * s)) +
               pow(m_rho, 2) * (pow(m_rho, 2) * (1.5 - 1. * C4 * s) +
                                s * (-2.5 + 3. * C4 * s)))) *
         tmax) /
            pow(m_rho, 6) -
        (0.25 *
         (pow(delta, 2) *
              (1. * pow(pion_mass, 6) - 1. * pow(m_rho, 6) +
               pow(pion_mass, 4) *
                   (-2.499999999999999 * pow(m_rho, 2) - 2.5 * s) -
               1.5 * pow(m_rho, 4) * s + 2. * pow(m_rho, 2) * pow(s, 2) -
               0.5 * pow(s, 3) +
               pow(pion_mass, 2) *
                   (3.5 * pow(m_rho, 4) -
                    1.5000000000000004 * pow(m_rho, 2) * s + 2. * pow(s, 2))) +
          pow(m_rho, 2) *
              (pow(pion_mass, 4) * (-6. - 8. * C4 * pow(m_rho, 2)) +
               2. * pow(s, 2) + pow(m_rho, 4) * (-4. - 8. * C4 * s) +
               pow(m_rho, 2) * s * (-2. + 8. * C4 * s) +
               pow(pion_mass, 2) * (8. * C4 * pow(m_rho, 4) + 4. * s +
                                    pow(m_rho, 2) * (10. - 16. * C4 * s))) +
          delta *
              (-2. * pow(pion_mass, 6) - 5. * pow(m_rho, 2) * pow(s, 2) +
               1. * pow(s, 3) +
               pow(pion_mass, 4) *
                   (8. * pow(m_rho, 2) + 4. * C4 * pow(m_rho, 4) + 5. * s) +
               pow(m_rho, 4) * s * (4. - 4. * C4 * s) +
               pow(m_rho, 6) * (4. + 4. * C4 * s) +
               pow(pion_mass, 2) *
                   (-4. * C4 * pow(m_rho, 6) + 1. * pow(m_rho, 2) * s -
                    4. * pow(s, 2) + pow(m_rho, 4) * (-12. + 8. * C4 * s)))) *
         tmax) /
            (pow(m_rho, 4) * (pow(pion_mass, 2) - 1. * s)) +
        (0.0625 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (pow(m_rho, 2) *
              (eta2 * (4. * pow(pion_mass, 4) - 6. * pow(pion_mass, 2) * s +
                       s * (8. * pow(m_rho, 2) + 6. * s)) +
               eta1 * (-12. * pow(pion_mass, 4) + 4. * pow(m_rho, 4) +
                       2. * pow(m_rho, 2) * s - 6. * pow(s, 2) +
                       pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 8. * s))) +
          delta *
              (eta1 *
                   (8. * pow(pion_mass, 6) - 2. * pow(m_rho, 6) +
                    pow(pion_mass, 4) * (2. * pow(m_rho, 2) - 2. * s) -
                    3. * pow(m_rho, 4) * s + 4. * pow(m_rho, 2) * pow(s, 2) +
                    1. * pow(s, 3) +
                    pow(pion_mass, 2) * (2. * pow(m_rho, 4) - 2. * pow(s, 2))) +
               eta2 * (pow(pion_mass, 4) * (-2. * pow(m_rho, 2) - 4. * s) +
                       pow(pion_mass, 2) * s * (3. * pow(m_rho, 2) + 3. * s) +
                       s * (-4. * pow(m_rho, 4) - 7. * pow(m_rho, 2) * s -
                            1. * pow(s, 2))))) *
         tmax) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) -
        (0.1875 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (delta * (eta1 * (2.6666666666666665 * pow(pion_mass, 6) +
                           pow(pion_mass, 4) * (-4. * pow(m_rho, 2) + 2. * s) +
                           pow(pion_mass, 2) *
                               (-1.3333333333333333 * pow(m_rho, 4) +
                                6. * pow(m_rho, 2) * s -
                                3.3333333333333335 * pow(s, 2)) +
                           s * (0.3333333333333333 * pow(m_rho, 4) -
                                1.3333333333333333 * pow(m_rho, 2) * s +
                                1. * pow(s, 2))) +
                   eta2 * (pow(pion_mass, 4) *
                               (-0.6666666666666666 * pow(m_rho, 2) - 4. * s) +
                           s * (0.6666666666666666 * pow(m_rho, 4) -
                                1. * pow(m_rho, 2) * s - 1. * pow(s, 2)) +
                           pow(pion_mass, 2) *
                               (-0.6666666666666666 * pow(m_rho, 4) -
                                0.3333333333333333 * pow(m_rho, 2) * s +
                                3.6666666666666665 * pow(s, 2)))) +
          pow(m_rho, 2) *
              (eta2 * (C4 * pow(pion_mass, 4) *
                           (2.6666666666666665 * pow(m_rho, 2) +
                            10.666666666666666 * s) +
                       pow(pion_mass, 2) *
                           (s * (3.3333333333333335 -
                                 10.666666666666666 * C4 * s) +
                            pow(m_rho, 2) * (1.3333333333333333 -
                                             5.333333333333333 * C4 * s)) +
                       s * (s * (-2. + 2.6666666666666665 * C4 * s) +
                            pow(m_rho, 2) * (-1.3333333333333333 +
                                             2.6666666666666665 * C4 * s))) +
               eta1 *
                   (pow(pion_mass, 4) *
                        (1.3333333333333333 + 8. * C4 * pow(m_rho, 2) -
                         10.666666666666666 * C4 * s) +
                    s * (s * (2. - 2.6666666666666665 * C4 * s) +
                         pow(m_rho, 2) * (-2. + 2.6666666666666665 * C4 * s)) +
                    pow(pion_mass, 2) *
                        (pow(m_rho, 2) * (2.6666666666666665 -
                                          10.666666666666666 * C4 * s) +
                         s * (-4. + 10.666666666666666 * C4 * s))))) *
         tmax) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) -
        (0.0625 * (-2. + delta) * (eta1 - 1. * eta2) *
         (pow(a1_mass, 2) - 1. * s) * (pow(pion_mass, 2) + s) *
         (-2. * eta2 * s +
          eta1 * (pow(pion_mass, 2) - 1. * pow(m_rho, 2) + s)) *
         pow(tmax, 2)) /
            ((-1. * pow(pion_mass, 2) + s) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta1, 2) * (pow(m_rho, 2) - 1. * s) + 2. * eta1 * eta2 * s -
          1. * pow(eta2, 2) * s) *
         (pow(pion_mass, 4) + (pow(m_rho, 2) - 1. * s) * s +
          pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 2. * s)) *
         pow(tmax, 2)) /
            (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
             2. * pow(a1_mass, 2) * s + pow(s, 2)) -
        (0.125 *
         (-1. * pow(m_rho, 4) + 4. * C4 * pow(m_rho, 6) +
          delta * pow(m_rho, 2) *
              (2. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (2. - 4. * C4 * pow(m_rho, 2)) - 2. * s) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 4) + 0.25 * pow(m_rho, 4) -
               1.25 * pow(s, 2) +
               pow(pion_mass, 2) * (-3. * pow(m_rho, 2) + 2. * s))) *
         pow(tmax, 2)) /
            pow(m_rho, 6) +
        (0.03125 *
         (0. - 4. * pow(m_rho, 4) +
          delta * (16. * pow(m_rho, 4) - 8. * pow(m_rho, 2) * s) +
          pow(delta, 2) *
              (4. * pow(pion_mass, 4) - 3. * pow(m_rho, 4) +
               2. * pow(m_rho, 2) * s - 3. * pow(s, 2) +
               pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 4. * s))) *
         pow(tmax, 2)) /
            pow(m_rho, 6) +
        (0.0625 *
         (-32. * C4 * pow(m_rho, 4) * s +
          pow(delta, 2) * (1. * pow(pion_mass, 4) +
                           pow(pion_mass, 2) *
                               (-1.0000000000000009 * pow(m_rho, 2) - 2. * s) +
                           s * (-3. * pow(m_rho, 2) + 1. * s)) +
          delta *
              (-2. * pow(pion_mass, 4) +
               (6. * pow(m_rho, 2) + 16. * C4 * pow(m_rho, 4) - 2. * s) * s +
               pow(pion_mass, 2) * (2. * pow(m_rho, 2) + 4. * s))) *
         pow(tmax, 2)) /
            (pow(m_rho, 4) * (pow(pion_mass, 2) - 1. * s)) -
        (0.5625 *
         (C4 * pow(m_rho, 6) *
              (2.6666666666666665 + 7.111111111111112 * C4 * pow(pion_mass, 2) -
               3.555555555555556 * C4 * s) +
          pow(delta, 2) * (0.11111111111111112 * pow(m_rho, 4) +
                           pow(pion_mass, 2) *
                               (1. * pow(m_rho, 2) - 0.22222222222222224 * s) -
                           0.22222222222222224 * pow(m_rho, 2) * s +
                           0.11111111111111112 * pow(s, 2)) +
          delta * pow(m_rho, 2) *
              (-2.2222222222222223 * C4 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (-0.6666666666666666 -
                                    2.6666666666666665 * C4 * pow(m_rho, 2)) +
               0.22222222222222224 * s +
               pow(m_rho, 2) *
                   (-0.22222222222222224 + 1.777777777777778 * C4 * s))) *
         pow(tmax, 2)) /
            pow(m_rho, 6) +
        (0.03125 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (pow(m_rho, 2) * (-2. * eta2 * pow(pion_mass, 2) -
                           5.999999999999999 * eta1 * pow(m_rho, 2) +
                           8. * eta1 * s - 2. * eta2 * s) +
          delta *
              (eta1 * (-5.999999999999999 * pow(pion_mass, 4) +
                       5. * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) -
                       5.999999999999999 * pow(m_rho, 2) * s + 1. * pow(s, 2)) +
               eta2 * (4. * pow(pion_mass, 4) +
                       pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 2. * s) +
                       s * (5. * pow(m_rho, 2) + 2. * s)))) *
         pow(tmax, 2)) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) -
        (0.15625 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (delta *
              (eta1 * (-1.2 * pow(pion_mass, 4) + 0.6 * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 2.4 * s) -
                       1.6 * pow(m_rho, 2) * s + 1. * pow(s, 2)) +
               eta2 * (0.8 * pow(pion_mass, 4) +
                       (1. * pow(m_rho, 2) - 0.4 * s) * s +
                       pow(pion_mass, 2) * (0.2 * pow(m_rho, 2) + 1.2 * s))) +
          pow(m_rho, 2) *
              (eta2 * (pow(pion_mass, 2) * (-0.4 - 6.4 * C4 * s) +
                       s * (-0.4 + 3.2 * C4 * s)) +
               eta1 * (s * (0.8 - 3.2 * C4 * s) +
                       pow(m_rho, 2) * (-0.4 + 3.2 * C4 * s) +
                       pow(pion_mass, 2) *
                           (-0.8 - 3.2 * C4 * pow(m_rho, 2) + 6.4 * C4 * s)))) *
         pow(tmax, 2)) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) -
        (0.20833333333333331 * delta *
         (-0.8 * pow(m_rho, 2) + 0.8 * C4 * pow(m_rho, 4) +
          delta * (0.8 * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 0.7 * s)) *
         pow(tmax, 3)) /
            pow(m_rho, 6) +
        (0.125 *
         (5.333333333333333 * pow(C4, 2) * pow(m_rho, 6) +
          delta * (-0.6666666666666666 * pow(m_rho, 2) -
                   1.3333333333333333 * C4 * pow(m_rho, 4)) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 2) + 1.1666666666666667 * pow(m_rho, 2) -
               0.6666666666666666 * s)) *
         pow(tmax, 3)) /
            pow(m_rho, 6) +
        (0.10416666666666666 * delta *
         (-0.8 * pow(m_rho, 2) +
          delta * (0.4 * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 0.6 * s)) *
         pow(tmax, 3)) /
            pow(m_rho, 6) +
        (0.020833333333333332 * pow(eta1 - 1. * eta2, 2) * s *
         (-2. * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-1. * pow(m_rho, 2) + s)) *
         pow(tmax, 3)) /
            (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
             2. * pow(a1_mass, 2) * s + pow(s, 2)) +
        (0.10416666666666666 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (0.4 * eta1 * pow(m_rho, 2) +
          delta * (-0.2 * eta2 * pow(pion_mass, 2) - 0.2 * eta2 * s +
                   eta1 * (-0.4 * pow(pion_mass, 2) - 0.8 * pow(m_rho, 2) +
                           1. * s))) *
         pow(tmax, 3)) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) -
        (0.14583333333333331 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (delta * (-0.14285714285714285 * eta2 * pow(pion_mass, 2) -
                   0.42857142857142855 * eta2 * s +
                   eta1 * (-0.2857142857142857 * pow(pion_mass, 2) -
                           0.5714285714285714 * pow(m_rho, 2) + 1. * s)) +
          pow(m_rho, 2) *
              (1.1428571428571428 * C4 * eta2 * s +
               eta1 * (0.2857142857142857 - 1.1428571428571428 * C4 * s))) *
         pow(tmax, 3)) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) +
        (0. -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(pion_mass, 4) *
             pow(m_rho, 6) +
         0.25 * pow(2. - 1. * delta, 2) * pow(m_rho, 10) -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(m_rho, 6) *
             pow(s, 2) +
         pow(2. - 1. * delta, 2) * pow(pion_mass, 2) * pow(m_rho, 6) *
             (-1. * pow(m_rho, 2) + 1. * s)) /
            (pow(m_rho, 6) * (-2. * pow(pion_mass, 2) + 1. * s + 1. * tmax)) +
        (2. *
         (0. - 2. * pow(pion_mass, 4) * pow(m_rho, 4) - 0.5 * pow(m_rho, 8) +
          delta * pow(m_rho, 4) *
              (2. * pow(pion_mass, 4) + 0.5 * pow(m_rho, 4) +
               pow(pion_mass, 2) *
                   (-2. * pow(m_rho, 2) - 1.9999999999999998 * s)) +
          pow(pion_mass, 2) * (2. * pow(m_rho, 6) + 2. * pow(m_rho, 4) * s) +
          pow(delta, 2) * pow(m_rho, 2) *
              (-2.220446049250313e-16 * pow(pion_mass, 6) -
               0.125 * pow(m_rho, 6) +
               pow(pion_mass, 4) *
                   (-0.5 * pow(m_rho, 2) + 2.220446049250313e-16 * s) +
               pow(pion_mass, 2) *
                   (0.5 * pow(m_rho, 4) + 0.5 * pow(m_rho, 2) * s))) *
         log(abs(-1. * pow(pion_mass, 2) + 0.5 * s + 0.5 * tmax))) /
            (pow(m_rho, 4) * (pow(pion_mass, 2) - 1. * s)) -
        (0.25 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (eta2 * ((-2. + 1. * delta) * pow(pion_mass, 6) +
                  (6. - 3. * delta) * pow(pion_mass, 4) * s +
                  pow(s, 2) * ((4. - 2. * delta) * pow(m_rho, 2) +
                               (2. - 1. * delta) * s) +
                  pow(pion_mass, 2) * s *
                      ((-4. + 2. * delta) * pow(m_rho, 2) +
                       (-6. + 3. * delta) * s)) +
          eta1 * ((2. - 1. * delta) * pow(pion_mass, 6) +
                  (2. - 1. * delta) * pow(m_rho, 4) * s +
                  (-2. + 1. * delta) * pow(s, 3) +
                  pow(pion_mass, 4) * ((4. - 2. * delta) * pow(m_rho, 2) +
                                       (-6. + 3. * delta) * s) +
                  pow(pion_mass, 2) * ((-2. + 1. * delta) * pow(m_rho, 4) +
                                       (-4. + 2. * delta) * pow(m_rho, 2) * s +
                                       (6. - 3. * delta) * pow(s, 2)))) *
         log(abs(-2. * pow(pion_mass, 2) + s + tmax))) /
            (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
             2. * pow(a1_mass, 2) * s + pow(s, 2)) +
        (0.25 *
         (0. + 8. * pow(pion_mass, 2) * pow(m_rho, 6) + 4. * pow(m_rho, 8) -
          8. * pow(m_rho, 6) * s +
          delta * pow(m_rho, 4) *
              (8. * pow(pion_mass, 4) - 8. * pow(m_rho, 4) +
               pow(pion_mass, 2) * (8. * pow(m_rho, 2) - 16. * s) +
               8. * pow(m_rho, 2) * s + 8. * pow(s, 2)) +
          pow(delta, 2) * pow(m_rho, 4) *
              (-4. * pow(pion_mass, 4) + 3. * pow(m_rho, 4) -
               2. * pow(m_rho, 2) * s - 4. * pow(s, 2) +
               pow(pion_mass, 2) * (-6. * pow(m_rho, 2) + 8. * s))) *
         log(abs(-2. * pow(pion_mass, 2) + 1. * s + 1. * tmax))) /
            pow(m_rho, 6) +
        (0.5 *
         (0. +
          pow(pion_mass, 2) * (4. * pow(m_rho, 6) - 8. * C4 * pow(m_rho, 8)) -
          4. * pow(m_rho, 6) * s + 8. * C4 * pow(m_rho, 8) * s +
          pow(delta, 2) * pow(m_rho, 4) *
              (2. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
               pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) +
               2. * pow(s, 2)) +
          delta * pow(m_rho, 4) *
              (-4. * pow(pion_mass, 4) + 2. * pow(m_rho, 2) * s -
               4. * pow(s, 2) +
               pow(pion_mass, 2) *
                   (-10. * pow(m_rho, 2) + 4. * C4 * pow(m_rho, 4) + 8. * s) +
               pow(m_rho, 4) * (2. - 4. * C4 * s))) *
         log(abs(-2. * pow(pion_mass, 2) + 1. * s + 1. * tmax))) /
            pow(m_rho, 6))) /
          (16. * Pi *
           (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
            2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
      (pow(Const, 2) * pow(ghat, 4) *
       (0. -
        (0.25 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
         (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
          2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
         tmin) /
            (pow(m_rho, 2) * pow(pow(pion_mass, 2) - s, 2)) -
        (0.0625 *
         (0. + 8. * pow(pion_mass, 2) * pow(m_rho, 4) - 12. * pow(m_rho, 6) +
          4. * pow(m_rho, 4) * s +
          delta * pow(m_rho, 2) *
              (-16. * pow(pion_mass, 4) -
               16. * pow(pion_mass, 2) * pow(m_rho, 2) - 4. * pow(m_rho, 4) +
               16. * pow(m_rho, 2) * s + 4. * pow(s, 2)) +
          pow(delta, 2) *
              (8. * pow(pion_mass, 6) + 9. * pow(m_rho, 6) +
               pow(pion_mass, 4) * (4. * pow(m_rho, 2) - 4. * s) -
               13. * pow(m_rho, 4) * s - 5. * pow(m_rho, 2) * pow(s, 2) +
               1. * pow(s, 3) +
               pow(pion_mass, 2) * (-2. * pow(m_rho, 4) +
                                    4. * pow(m_rho, 2) * s - 2. * pow(s, 2)))) *
         tmin) /
            pow(m_rho, 6) -
        (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
         (pow(a1_mass, 2) - 1. * s) *
         (eta2 * (pow(pion_mass, 6) + pow(pion_mass, 2) * pow(s, 2) +
                  (pow(m_rho, 2) - 1. * s) * pow(s, 2) +
                  pow(pion_mass, 4) * (-1. * pow(m_rho, 2) + 3. * s)) +
          eta1 * (-4. * pow(pion_mass, 6) +
                  pow(pion_mass, 4) * (3. * pow(m_rho, 2) + s) +
                  pow(pion_mass, 2) * (-1. * pow(m_rho, 4) + pow(m_rho, 2) * s -
                                       2. * pow(s, 2)) +
                  s * (pow(m_rho, 4) - 2. * pow(m_rho, 2) * s + pow(s, 2)))) *
         tmin) /
            ((-1. * pow(pion_mass, 2) + s) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta1, 2) *
              (1. * pow(pion_mass, 8) - 2. * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(pion_mass, 2) * s *
                   (-4. * pow(m_rho, 4) + 8. * pow(m_rho, 2) * s -
                    4. * pow(s, 2)) +
               pow(s, 2) * (1. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s +
                            1. * pow(s, 2)) +
               pow(pion_mass, 4) * (3. * pow(m_rho, 4) -
                                    6. * pow(m_rho, 2) * s + 4. * pow(s, 2))) +
          pow(eta2, 2) *
              (1. * pow(pion_mass, 8) - 2. * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(pion_mass, 2) * s *
                   (-2. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s -
                    4. * pow(s, 2)) +
               pow(s, 2) * (1. * pow(m_rho, 4) + 2. * pow(m_rho, 2) * s +
                            1. * pow(s, 2)) +
               pow(pion_mass, 4) * (1. * pow(m_rho, 4) +
                                    4. * pow(m_rho, 2) * s + 4. * pow(s, 2))) +
          eta1 * eta2 *
              (-2. * pow(pion_mass, 8) + 2. * pow(m_rho, 4) * pow(s, 2) -
               2. * pow(s, 4) +
               pow(pion_mass, 4) * (2. * pow(m_rho, 4) +
                                    4. * pow(m_rho, 2) * s - 8. * pow(s, 2)) +
               pow(pion_mass, 2) * s *
                   (-4. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s +
                    8. * pow(s, 2)))) *
         tmin) /
            (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
             2. * pow(a1_mass, 2) * s + pow(s, 2)) +
        (0.5 *
         (0. - 4. * C4 * pow(m_rho, 8) - 0.5 * pow(m_rho, 4) * s +
          pow(m_rho, 6) * (2. + 2. * C4 * s) +
          pow(delta, 2) *
              (2. * pow(pion_mass, 6) + 2. * pow(m_rho, 6) -
               1.5 * pow(pion_mass, 4) * s - 2.375 * pow(m_rho, 4) * s -
               0.75 * pow(m_rho, 2) * pow(s, 2) + 0.125 * pow(s, 3) +
               pow(pion_mass, 2) *
                   (-1.5 * pow(m_rho, 4) + 0.5 * pow(m_rho, 2) * s)) +
          delta * pow(m_rho, 2) *
              (-2. * pow(pion_mass, 4) +
               pow(pion_mass, 2) *
                   (1. * s + pow(m_rho, 2) * (-1. - 2. * C4 * s)) +
               pow(m_rho, 2) * (2. * C4 * pow(m_rho, 4) +
                                pow(m_rho, 2) * (-3. + 1. * C4 * s) +
                                s * (2. + 1. * C4 * s)))) *
         tmin) /
            pow(m_rho, 6) -
        (0.5 *
         (pow(m_rho, 6) *
              (-1.5 + C4 * (-12. * pow(pion_mass, 2) + 6. * s) +
               pow(C4, 2) * (-16. * pow(pion_mass, 4) +
                             16. * pow(pion_mass, 2) * s - 4. * pow(s, 2))) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 6) + 0.125 * pow(m_rho, 6) +
               pow(pion_mass, 4) * (-2. * pow(m_rho, 2) - 1. * s) +
               0.25 * pow(m_rho, 4) * s - 0.625 * pow(m_rho, 2) * pow(s, 2) +
               pow(pion_mass, 2) *
                   (-2.5 * pow(m_rho, 4) + 1.75 * pow(m_rho, 2) * s +
                    0.25 * pow(s, 2))) +
          delta * pow(m_rho, 2) *
              (pow(pion_mass, 4) * (1. + 8. * C4 * pow(m_rho, 2)) +
               pow(pion_mass, 2) * (6. * C4 * pow(m_rho, 4) - 0.5 * s +
                                    pow(m_rho, 2) * (3. - 10. * C4 * s)) +
               pow(m_rho, 2) * (pow(m_rho, 2) * (1.5 - 1. * C4 * s) +
                                s * (-2.5 + 3. * C4 * s)))) *
         tmin) /
            pow(m_rho, 6) -
        (0.25 *
         (pow(delta, 2) *
              (1. * pow(pion_mass, 6) - 1. * pow(m_rho, 6) +
               pow(pion_mass, 4) *
                   (-2.499999999999999 * pow(m_rho, 2) - 2.5 * s) -
               1.5 * pow(m_rho, 4) * s + 2. * pow(m_rho, 2) * pow(s, 2) -
               0.5 * pow(s, 3) +
               pow(pion_mass, 2) *
                   (3.5 * pow(m_rho, 4) -
                    1.5000000000000004 * pow(m_rho, 2) * s + 2. * pow(s, 2))) +
          pow(m_rho, 2) *
              (pow(pion_mass, 4) * (-6. - 8. * C4 * pow(m_rho, 2)) +
               2. * pow(s, 2) + pow(m_rho, 4) * (-4. - 8. * C4 * s) +
               pow(m_rho, 2) * s * (-2. + 8. * C4 * s) +
               pow(pion_mass, 2) * (8. * C4 * pow(m_rho, 4) + 4. * s +
                                    pow(m_rho, 2) * (10. - 16. * C4 * s))) +
          delta *
              (-2. * pow(pion_mass, 6) - 5. * pow(m_rho, 2) * pow(s, 2) +
               1. * pow(s, 3) +
               pow(pion_mass, 4) *
                   (8. * pow(m_rho, 2) + 4. * C4 * pow(m_rho, 4) + 5. * s) +
               pow(m_rho, 4) * s * (4. - 4. * C4 * s) +
               pow(m_rho, 6) * (4. + 4. * C4 * s) +
               pow(pion_mass, 2) *
                   (-4. * C4 * pow(m_rho, 6) + 1. * pow(m_rho, 2) * s -
                    4. * pow(s, 2) + pow(m_rho, 4) * (-12. + 8. * C4 * s)))) *
         tmin) /
            (pow(m_rho, 4) * (pow(pion_mass, 2) - 1. * s)) +
        (0.0625 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (pow(m_rho, 2) *
              (eta2 * (4. * pow(pion_mass, 4) - 6. * pow(pion_mass, 2) * s +
                       s * (8. * pow(m_rho, 2) + 6. * s)) +
               eta1 * (-12. * pow(pion_mass, 4) + 4. * pow(m_rho, 4) +
                       2. * pow(m_rho, 2) * s - 6. * pow(s, 2) +
                       pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 8. * s))) +
          delta *
              (eta1 *
                   (8. * pow(pion_mass, 6) - 2. * pow(m_rho, 6) +
                    pow(pion_mass, 4) * (2. * pow(m_rho, 2) - 2. * s) -
                    3. * pow(m_rho, 4) * s + 4. * pow(m_rho, 2) * pow(s, 2) +
                    1. * pow(s, 3) +
                    pow(pion_mass, 2) * (2. * pow(m_rho, 4) - 2. * pow(s, 2))) +
               eta2 * (pow(pion_mass, 4) * (-2. * pow(m_rho, 2) - 4. * s) +
                       pow(pion_mass, 2) * s * (3. * pow(m_rho, 2) + 3. * s) +
                       s * (-4. * pow(m_rho, 4) - 7. * pow(m_rho, 2) * s -
                            1. * pow(s, 2))))) *
         tmin) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) -
        (0.1875 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (delta * (eta1 * (2.6666666666666665 * pow(pion_mass, 6) +
                           pow(pion_mass, 4) * (-4. * pow(m_rho, 2) + 2. * s) +
                           pow(pion_mass, 2) *
                               (-1.3333333333333333 * pow(m_rho, 4) +
                                6. * pow(m_rho, 2) * s -
                                3.3333333333333335 * pow(s, 2)) +
                           s * (0.3333333333333333 * pow(m_rho, 4) -
                                1.3333333333333333 * pow(m_rho, 2) * s +
                                1. * pow(s, 2))) +
                   eta2 * (pow(pion_mass, 4) *
                               (-0.6666666666666666 * pow(m_rho, 2) - 4. * s) +
                           s * (0.6666666666666666 * pow(m_rho, 4) -
                                1. * pow(m_rho, 2) * s - 1. * pow(s, 2)) +
                           pow(pion_mass, 2) *
                               (-0.6666666666666666 * pow(m_rho, 4) -
                                0.3333333333333333 * pow(m_rho, 2) * s +
                                3.6666666666666665 * pow(s, 2)))) +
          pow(m_rho, 2) *
              (eta2 * (C4 * pow(pion_mass, 4) *
                           (2.6666666666666665 * pow(m_rho, 2) +
                            10.666666666666666 * s) +
                       pow(pion_mass, 2) *
                           (s * (3.3333333333333335 -
                                 10.666666666666666 * C4 * s) +
                            pow(m_rho, 2) * (1.3333333333333333 -
                                             5.333333333333333 * C4 * s)) +
                       s * (s * (-2. + 2.6666666666666665 * C4 * s) +
                            pow(m_rho, 2) * (-1.3333333333333333 +
                                             2.6666666666666665 * C4 * s))) +
               eta1 *
                   (pow(pion_mass, 4) *
                        (1.3333333333333333 + 8. * C4 * pow(m_rho, 2) -
                         10.666666666666666 * C4 * s) +
                    s * (s * (2. - 2.6666666666666665 * C4 * s) +
                         pow(m_rho, 2) * (-2. + 2.6666666666666665 * C4 * s)) +
                    pow(pion_mass, 2) *
                        (pow(m_rho, 2) * (2.6666666666666665 -
                                          10.666666666666666 * C4 * s) +
                         s * (-4. + 10.666666666666666 * C4 * s))))) *
         tmin) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) -
        (0.0625 * (-2. + delta) * (eta1 - 1. * eta2) *
         (pow(a1_mass, 2) - 1. * s) * (pow(pion_mass, 2) + s) *
         (-2. * eta2 * s +
          eta1 * (pow(pion_mass, 2) - 1. * pow(m_rho, 2) + s)) *
         pow(tmin, 2)) /
            ((-1. * pow(pion_mass, 2) + s) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta1, 2) * (pow(m_rho, 2) - 1. * s) + 2. * eta1 * eta2 * s -
          1. * pow(eta2, 2) * s) *
         (pow(pion_mass, 4) + (pow(m_rho, 2) - 1. * s) * s +
          pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 2. * s)) *
         pow(tmin, 2)) /
            (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
             2. * pow(a1_mass, 2) * s + pow(s, 2)) -
        (0.125 *
         (-1. * pow(m_rho, 4) + 4. * C4 * pow(m_rho, 6) +
          delta * pow(m_rho, 2) *
              (2. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (2. - 4. * C4 * pow(m_rho, 2)) - 2. * s) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 4) + 0.25 * pow(m_rho, 4) -
               1.25 * pow(s, 2) +
               pow(pion_mass, 2) * (-3. * pow(m_rho, 2) + 2. * s))) *
         pow(tmin, 2)) /
            pow(m_rho, 6) +
        (0.03125 *
         (0. - 4. * pow(m_rho, 4) +
          delta * (16. * pow(m_rho, 4) - 8. * pow(m_rho, 2) * s) +
          pow(delta, 2) *
              (4. * pow(pion_mass, 4) - 3. * pow(m_rho, 4) +
               2. * pow(m_rho, 2) * s - 3. * pow(s, 2) +
               pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 4. * s))) *
         pow(tmin, 2)) /
            pow(m_rho, 6) +
        (0.0625 *
         (-32. * C4 * pow(m_rho, 4) * s +
          pow(delta, 2) * (1. * pow(pion_mass, 4) +
                           pow(pion_mass, 2) *
                               (-1.0000000000000009 * pow(m_rho, 2) - 2. * s) +
                           s * (-3. * pow(m_rho, 2) + 1. * s)) +
          delta *
              (-2. * pow(pion_mass, 4) +
               (6. * pow(m_rho, 2) + 16. * C4 * pow(m_rho, 4) - 2. * s) * s +
               pow(pion_mass, 2) * (2. * pow(m_rho, 2) + 4. * s))) *
         pow(tmin, 2)) /
            (pow(m_rho, 4) * (pow(pion_mass, 2) - 1. * s)) -
        (0.5625 *
         (C4 * pow(m_rho, 6) *
              (2.6666666666666665 + 7.111111111111112 * C4 * pow(pion_mass, 2) -
               3.555555555555556 * C4 * s) +
          pow(delta, 2) * (0.11111111111111112 * pow(m_rho, 4) +
                           pow(pion_mass, 2) *
                               (1. * pow(m_rho, 2) - 0.22222222222222224 * s) -
                           0.22222222222222224 * pow(m_rho, 2) * s +
                           0.11111111111111112 * pow(s, 2)) +
          delta * pow(m_rho, 2) *
              (-2.2222222222222223 * C4 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (-0.6666666666666666 -
                                    2.6666666666666665 * C4 * pow(m_rho, 2)) +
               0.22222222222222224 * s +
               pow(m_rho, 2) *
                   (-0.22222222222222224 + 1.777777777777778 * C4 * s))) *
         pow(tmin, 2)) /
            pow(m_rho, 6) +
        (0.03125 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (pow(m_rho, 2) * (-2. * eta2 * pow(pion_mass, 2) -
                           5.999999999999999 * eta1 * pow(m_rho, 2) +
                           8. * eta1 * s - 2. * eta2 * s) +
          delta *
              (eta1 * (-5.999999999999999 * pow(pion_mass, 4) +
                       5. * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) -
                       5.999999999999999 * pow(m_rho, 2) * s + 1. * pow(s, 2)) +
               eta2 * (4. * pow(pion_mass, 4) +
                       pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 2. * s) +
                       s * (5. * pow(m_rho, 2) + 2. * s)))) *
         pow(tmin, 2)) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) -
        (0.15625 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (delta *
              (eta1 * (-1.2 * pow(pion_mass, 4) + 0.6 * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 2.4 * s) -
                       1.6 * pow(m_rho, 2) * s + 1. * pow(s, 2)) +
               eta2 * (0.8 * pow(pion_mass, 4) +
                       (1. * pow(m_rho, 2) - 0.4 * s) * s +
                       pow(pion_mass, 2) * (0.2 * pow(m_rho, 2) + 1.2 * s))) +
          pow(m_rho, 2) *
              (eta2 * (pow(pion_mass, 2) * (-0.4 - 6.4 * C4 * s) +
                       s * (-0.4 + 3.2 * C4 * s)) +
               eta1 * (s * (0.8 - 3.2 * C4 * s) +
                       pow(m_rho, 2) * (-0.4 + 3.2 * C4 * s) +
                       pow(pion_mass, 2) *
                           (-0.8 - 3.2 * C4 * pow(m_rho, 2) + 6.4 * C4 * s)))) *
         pow(tmin, 2)) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) -
        (0.20833333333333331 * delta *
         (-0.8 * pow(m_rho, 2) + 0.8 * C4 * pow(m_rho, 4) +
          delta * (0.8 * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 0.7 * s)) *
         pow(tmin, 3)) /
            pow(m_rho, 6) +
        (0.125 *
         (5.333333333333333 * pow(C4, 2) * pow(m_rho, 6) +
          delta * (-0.6666666666666666 * pow(m_rho, 2) -
                   1.3333333333333333 * C4 * pow(m_rho, 4)) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 2) + 1.1666666666666667 * pow(m_rho, 2) -
               0.6666666666666666 * s)) *
         pow(tmin, 3)) /
            pow(m_rho, 6) +
        (0.10416666666666666 * delta *
         (-0.8 * pow(m_rho, 2) +
          delta * (0.4 * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 0.6 * s)) *
         pow(tmin, 3)) /
            pow(m_rho, 6) +
        (0.020833333333333332 * pow(eta1 - 1. * eta2, 2) * s *
         (-2. * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-1. * pow(m_rho, 2) + s)) *
         pow(tmin, 3)) /
            (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
             2. * pow(a1_mass, 2) * s + pow(s, 2)) +
        (0.10416666666666666 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (0.4 * eta1 * pow(m_rho, 2) +
          delta * (-0.2 * eta2 * pow(pion_mass, 2) - 0.2 * eta2 * s +
                   eta1 * (-0.4 * pow(pion_mass, 2) - 0.8 * pow(m_rho, 2) +
                           1. * s))) *
         pow(tmin, 3)) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) -
        (0.14583333333333331 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (delta * (-0.14285714285714285 * eta2 * pow(pion_mass, 2) -
                   0.42857142857142855 * eta2 * s +
                   eta1 * (-0.2857142857142857 * pow(pion_mass, 2) -
                           0.5714285714285714 * pow(m_rho, 2) + 1. * s)) +
          pow(m_rho, 2) *
              (1.1428571428571428 * C4 * eta2 * s +
               eta1 * (0.2857142857142857 - 1.1428571428571428 * C4 * s))) *
         pow(tmin, 3)) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * s + pow(s, 2))) +
        (0. -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(pion_mass, 4) *
             pow(m_rho, 6) +
         0.25 * pow(2. - 1. * delta, 2) * pow(m_rho, 10) -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(m_rho, 6) *
             pow(s, 2) +
         pow(2. - 1. * delta, 2) * pow(pion_mass, 2) * pow(m_rho, 6) *
             (-1. * pow(m_rho, 2) + 1. * s)) /
            (pow(m_rho, 6) * (-2. * pow(pion_mass, 2) + 1. * s + 1. * tmin)) +
        (2. *
         (0. - 2. * pow(pion_mass, 4) * pow(m_rho, 4) - 0.5 * pow(m_rho, 8) +
          delta * pow(m_rho, 4) *
              (2. * pow(pion_mass, 4) + 0.5 * pow(m_rho, 4) +
               pow(pion_mass, 2) *
                   (-2. * pow(m_rho, 2) - 1.9999999999999998 * s)) +
          pow(pion_mass, 2) * (2. * pow(m_rho, 6) + 2. * pow(m_rho, 4) * s) +
          pow(delta, 2) * pow(m_rho, 2) *
              (-2.220446049250313e-16 * pow(pion_mass, 6) -
               0.125 * pow(m_rho, 6) +
               pow(pion_mass, 4) *
                   (-0.5 * pow(m_rho, 2) + 2.220446049250313e-16 * s) +
               pow(pion_mass, 2) *
                   (0.5 * pow(m_rho, 4) + 0.5 * pow(m_rho, 2) * s))) *
         log(abs(-1. * pow(pion_mass, 2) + 0.5 * s + 0.5 * tmin))) /
            (pow(m_rho, 4) * (pow(pion_mass, 2) - 1. * s)) -
        (0.25 * (eta1 - 1. * eta2) * (pow(a1_mass, 2) - 1. * s) *
         (eta2 * ((-2. + 1. * delta) * pow(pion_mass, 6) +
                  (6. - 3. * delta) * pow(pion_mass, 4) * s +
                  pow(s, 2) * ((4. - 2. * delta) * pow(m_rho, 2) +
                               (2. - 1. * delta) * s) +
                  pow(pion_mass, 2) * s *
                      ((-4. + 2. * delta) * pow(m_rho, 2) +
                       (-6. + 3. * delta) * s)) +
          eta1 * ((2. - 1. * delta) * pow(pion_mass, 6) +
                  (2. - 1. * delta) * pow(m_rho, 4) * s +
                  (-2. + 1. * delta) * pow(s, 3) +
                  pow(pion_mass, 4) * ((4. - 2. * delta) * pow(m_rho, 2) +
                                       (-6. + 3. * delta) * s) +
                  pow(pion_mass, 2) * ((-2. + 1. * delta) * pow(m_rho, 4) +
                                       (-4. + 2. * delta) * pow(m_rho, 2) * s +
                                       (6. - 3. * delta) * pow(s, 2)))) *
         log(abs(-2. * pow(pion_mass, 2) + s + tmin))) /
            (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
             2. * pow(a1_mass, 2) * s + pow(s, 2)) +
        (0.25 *
         (0. + 8. * pow(pion_mass, 2) * pow(m_rho, 6) + 4. * pow(m_rho, 8) -
          8. * pow(m_rho, 6) * s +
          delta * pow(m_rho, 4) *
              (8. * pow(pion_mass, 4) - 8. * pow(m_rho, 4) +
               pow(pion_mass, 2) * (8. * pow(m_rho, 2) - 16. * s) +
               8. * pow(m_rho, 2) * s + 8. * pow(s, 2)) +
          pow(delta, 2) * pow(m_rho, 4) *
              (-4. * pow(pion_mass, 4) + 3. * pow(m_rho, 4) -
               2. * pow(m_rho, 2) * s - 4. * pow(s, 2) +
               pow(pion_mass, 2) * (-6. * pow(m_rho, 2) + 8. * s))) *
         log(abs(-2. * pow(pion_mass, 2) + 1. * s + 1. * tmin))) /
            pow(m_rho, 6) +
        (0.5 *
         (0. +
          pow(pion_mass, 2) * (4. * pow(m_rho, 6) - 8. * C4 * pow(m_rho, 8)) -
          4. * pow(m_rho, 6) * s + 8. * C4 * pow(m_rho, 8) * s +
          pow(delta, 2) * pow(m_rho, 4) *
              (2. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
               pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) +
               2. * pow(s, 2)) +
          delta * pow(m_rho, 4) *
              (-4. * pow(pion_mass, 4) + 2. * pow(m_rho, 2) * s -
               4. * pow(s, 2) +
               pow(pion_mass, 2) *
                   (-10. * pow(m_rho, 2) + 4. * C4 * pow(m_rho, 4) + 8. * s) +
               pow(m_rho, 4) * (2. - 4. * C4 * s))) *
         log(abs(-2. * pow(pion_mass, 2) + 1. * s + 1. * tmin))) /
            pow(m_rho, 6))) /
          (16. * Pi *
           (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
            2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)));

  // clang-format on
  return cut_off(xs * gev2_mb / spin_deg_factor);
}

double CrosssectionsPhoton<ComputationMethod::Analytic>::
    xs_diff_pi0_rho_pi_rho_mediated(const double s, const double t,
                                    const double m_rho) {
  const double spin_deg_factor = 3.0;

  // clang-format off
  const double diff_xs =
      (pow(Const, 2) * pow(ghat, 4) *
       ((-0.25 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
         (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
          2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) /
            (pow(m_rho, 2) * pow(pow(pion_mass, 2) - s, 2)) -
        (0.0625 * (eta1 - eta2) * (-pow(a1_mass, 2) + s) *
         (2 * pow(m_rho, 2) +
          delta * (-2 * pow(pion_mass, 2) - pow(m_rho, 2) + s + t)) *
         (-(eta2 * (s - t) *
            (4 * pow(pion_mass, 4) + s * (4 * pow(m_rho, 2) + s - t) -
             pow(pion_mass, 2) * (3 * s + t))) +
          eta1 *
              (8 * pow(pion_mass, 6) + pow(s, 3) + pow(s, 2) * t +
               5 * s * pow(t, 2) + pow(t, 3) + 2 * pow(m_rho, 4) * (-s + t) +
               pow(m_rho, 2) * (s - 3 * t) * (s + t) +
               2 * pow(pion_mass, 2) * (2 * pow(m_rho, 2) - s - t) * (s + t) -
               2 * pow(pion_mass, 4) * (2 * pow(m_rho, 2) + s + 3 * t)))) /
            (pow(m_rho, 2) *
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (-2 * pow(pion_mass, 2) + s + t)) -
        (0.0625 *
         pow(-2. * pow(m_rho, 2) + delta * (2. * pow(pion_mass, 2) +
                                            pow(m_rho, 2) - 1. * s - 1. * t),
             2) *
         (8. * pow(pion_mass, 6) + 4. * pow(m_rho, 6) + pow(s, 3) +
          pow(m_rho, 4) * (-4. * s - 4. * t) +
          pow(pion_mass, 4) * (-4. * pow(m_rho, 2) - 4. * s - 4. * t) +
          3. * pow(s, 2) * t + 3. * s * pow(t, 2) + pow(t, 3) +
          pow(m_rho, 2) * (-3. * pow(s, 2) + 2. * s * t - 3. * pow(t, 2)) +
          pow(pion_mass, 2) *
              (-8. * pow(m_rho, 4) - 2. * pow(s, 2) - 4. * s * t -
               2. * pow(t, 2) + pow(m_rho, 2) * (4. * s + 4. * t)))) /
            (pow(m_rho, 6) * pow(2. * pow(pion_mass, 2) - 1. * s - 1. * t, 2)) +
        (0.125 * (-2 + delta) * (eta1 - eta2) * (-pow(a1_mass, 2) + s) *
         (-(eta2 * (pow(pion_mass, 2) + s) *
            (-pow(pion_mass, 4) + pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s) +
             s * (-pow(m_rho, 2) + s + 2 * t))) +
          eta1 * (-4 * pow(pion_mass, 6) +
                  s * (-pow(m_rho, 2) + s) * (-pow(m_rho, 2) + s + t) +
                  pow(pion_mass, 4) * (3 * pow(m_rho, 2) + s + t) -
                  pow(pion_mass, 2) * (pow(m_rho, 4) + 2 * s * (s - t) +
                                       pow(m_rho, 2) * (-s + t))))) /
            ((pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) *
             (-pow(pion_mass, 2) + s)) +
        (0.03125 * pow(eta1 - eta2, 2) *
         (pow(eta2, 2) *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(pion_mass, 4) * (pow(pow(m_rho, 2) + 2 * s, 2) - 2 * s * t) +
               pow(s, 2) * (pow(pow(m_rho, 2) + s, 2) +
                            2 * (-pow(m_rho, 2) + s) * t + 2 * pow(t, 2)) -
               2 * pow(pion_mass, 2) * s *
                   (pow(m_rho, 4) + pow(m_rho, 2) * (2 * s - t) +
                    2 * s * (s + t))) -
          2 * eta1 * eta2 *
              (pow(pion_mass, 8) -
               pow(pion_mass, 4) * (pow(m_rho, 4) + 2 * pow(m_rho, 2) * s +
                                    2 * s * (-2 * s + t)) +
               2 * pow(pion_mass, 2) * s *
                   (pow(m_rho, 4) + pow(m_rho, 2) * (s + t) - 2 * s * (s + t)) +
               pow(s, 2) * (-pow(m_rho, 4) + pow(s, 2) - 2 * pow(m_rho, 2) * t +
                            2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
               pow(pion_mass, 4) * (3 * pow(m_rho, 4) + 2 * s * (2 * s - t) +
                                    2 * pow(m_rho, 2) * (-3 * s + t)) -
               2 * pow(pion_mass, 2) * (-pow(m_rho, 2) + s) *
                   (2 * s * (s + t) - pow(m_rho, 2) * (2 * s + t)) +
               s * (-pow(m_rho, 2) + s) *
                   (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                    pow(m_rho, 2) * (s + 2 * t))))) /
            (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(pow(a1_mass, 2) - s, 2)) +
        (0.5 *
         (-2. * pow(m_rho, 2) +
          delta * (2. * pow(pion_mass, 2) + pow(m_rho, 2) - 1. * s - 1. * t)) *
         (delta *
              (1. * pow(pion_mass, 6) + 0.5 * pow(m_rho, 6) +
               0.0625 * pow(s, 3) + pow(m_rho, 4) * (-0.5 * s - 0.5 * t) +
               pow(pion_mass, 4) *
                   (-0.5 * pow(m_rho, 2) - 0.75 * s - 0.25 * t) +
               0.3125 * pow(s, 2) * t + 0.4375 * s * pow(t, 2) +
               0.1875 * pow(t, 3) +
               pow(pion_mass, 2) * (-1. * pow(m_rho, 4) +
                                    pow(m_rho, 2) * (0.375 * s + 0.625 * t) +
                                    (-0.5 * s - 0.5 * t) * t) +
               pow(m_rho, 2) *
                   (-0.3125 * pow(s, 2) + 0.25 * s * t - 0.4375 * pow(t, 2))) +
          pow(m_rho, 2) *
              (-0.125 * pow(s, 2) + C4 * pow(m_rho, 4) * (1. * s - 1. * t) +
               0.125 * pow(t, 2) +
               pow(pion_mass, 2) * ((0.25 - 1. * C4 * pow(m_rho, 2)) * s +
                                    (-0.25 + 1. * C4 * pow(m_rho, 2)) * t) +
               pow(m_rho, 2) * (-0.5 * s + 0.5 * C4 * pow(s, 2) +
                                t * (0.5 - 0.5 * C4 * t))))) /
            (pow(m_rho, 6) * (1. * pow(pion_mass, 2) - 0.5 * s - 0.5 * t)) +
        (pow(delta, 2) *
             (-0.5 * pow(pion_mass, 6) - 0.0625 * pow(m_rho, 6) +
              pow(pion_mass, 4) * (1. * pow(m_rho, 2) + 0.5 * s) +
              pow(m_rho, 4) * (-0.125 * s - 0.125 * t) +
              t * (-0.125 * pow(s, 2) - 0.25 * s * t - 0.125 * pow(t, 2)) +
              pow(pion_mass, 2) * (1.25 * pow(m_rho, 4) - 0.125 * pow(s, 2) +
                                   pow(m_rho, 2) * (-0.875 * s - 1.125 * t) +
                                   0.25 * s * t + 0.375 * pow(t, 2)) +
              pow(m_rho, 2) *
                  (0.3125 * pow(s, 2) + 0.25 * s * t + 0.4375 * pow(t, 2))) +
         delta * pow(m_rho, 2) *
             (pow(pion_mass, 4) * (-0.5 - 4. * C4 * pow(m_rho, 2)) +
              (-0.25 * s - 0.25 * t) * t +
              pow(m_rho, 4) * (-0.75 + 0.5 * C4 * s + 2.5 * C4 * t) +
              pow(m_rho, 2) *
                  (-1.5 * C4 * pow(s, 2) + s * (1.25 - 2. * C4 * t) +
                   t * (0.25 - 0.5 * C4 * t)) +
              pow(pion_mass, 2) *
                  (-3. * C4 * pow(m_rho, 4) + 0.25 * s + 0.75 * t +
                   pow(m_rho, 2) * (-1.5 + 5. * C4 * s + 3. * C4 * t))) +
         pow(m_rho, 6) *
             (0.75 +
              C4 * (8. * C4 * pow(pion_mass, 4) + 2. * C4 * pow(s, 2) +
                    pow(pion_mass, 2) * (6. - 8. * C4 * s - 8. * C4 * t) +
                    t * (-3. + 2. * C4 * t) + s * (-3. + 4. * C4 * t)))) /
            pow(m_rho, 6) +
        (0.0625 * (eta1 - eta2) * (-pow(a1_mass, 2) + s) *
         (-(eta2 *
            (2 * pow(pion_mass, 4) *
                 (-4 * C4 * pow(m_rho, 2) * (pow(m_rho, 2) + 4 * s) +
                  delta * (pow(m_rho, 2) + 6 * s - 2 * t)) +
             pow(pion_mass, 2) *
                 (2 * pow(m_rho, 4) * (-2 + delta + 8 * C4 * s) +
                  delta * (-11 * pow(s, 2) - 6 * s * t + pow(t, 2)) +
                  pow(m_rho, 2) * ((-10 + delta) * s - (-2 + delta) * t +
                                   32 * C4 * s * (s + t))) +
             s * (-2 * pow(m_rho, 4) * (-2 + delta + 4 * C4 * s) +
                  delta * (3 * pow(s, 2) + 2 * s * t + 3 * pow(t, 2)) +
                  pow(m_rho, 2) * (3 * (2 + delta) * s + (2 - 5 * delta) * t -
                                   8 * C4 * pow(s + t, 2))))) +
          eta1 *
              (8 * delta * pow(pion_mass, 6) +
               2 * pow(pion_mass, 4) *
                   (12 * C4 * pow(m_rho, 4) -
                    2 * pow(m_rho, 2) * (-1 + 3 * delta + 8 * C4 * s) +
                    3 * delta * (s - t)) +
               delta * (3 * pow(s, 3) + 5 * pow(s, 2) * t + 7 * s * pow(t, 2) +
                        pow(t, 3)) -
               2 * pow(m_rho, 2) *
                   ((-3 + 2 * delta) * pow(s, 2) +
                    2 * (-1 + 2 * delta) * s * t +
                    (-1 + 2 * delta) * pow(t, 2) + 4 * C4 * s * pow(s + t, 2)) +
               pow(m_rho, 4) * ((-6 + delta) * s + (-2 + 3 * delta) * t +
                                8 * C4 * s * (s + 2 * t)) -
               2 * pow(pion_mass, 2) *
                   (delta * (s + t) * (5 * s + t) -
                    pow(m_rho, 2) * (-6 * s + 9 * delta * s - 2 * t +
                                     5 * delta * t + 16 * C4 * s * (s + t)) +
                    2 * pow(m_rho, 4) *
                        (-2 + delta + 4 * C4 * (2 * s + t)))))) /
            (pow(m_rho, 2) * (pow(Gammaa1, 2) * pow(a1_mass, 2) +
                              pow(pow(a1_mass, 2) - s, 2))) +
        (2 *
         ((0.0625 * (-2. + delta) *
           (-2. * pow(m_rho, 2) + delta * (2. * pow(pion_mass, 2) +
                                           pow(m_rho, 2) - 1. * s - 1. * t)) *
           (2. * pow(pion_mass, 6) + 1. * pow(m_rho, 6) +
            pow(pion_mass, 4) * (-3. * pow(m_rho, 2) - 2. * s) +
            pow(m_rho, 4) * (-1.5 * s - 1.5 * t) +
            pow(m_rho, 2) * (0.5 * s + 0.5 * t) * t +
            s * (0.5 * pow(s, 2) + 1. * s * t + 0.5 * pow(t, 2)) +
            pow(pion_mass, 2) *
                (-1. * pow(m_rho, 4) - 0.5 * pow(s, 2) - 1. * s * t -
                 0.5 * pow(t, 2) + pow(m_rho, 2) * (-0.5 * s + 2.5 * t)))) /
              ((pow(pion_mass, 2) - 1. * s) *
               (1. * pow(pion_mass, 2) - 0.5 * s - 0.5 * t)) +
          (0.0625 * (-2 + delta) *
           (delta * (6 * pow(pion_mass, 6) -
                     pow(pion_mass, 4) * (9 * (pow(m_rho, 2) + s) + t) -
                     pow(pion_mass, 2) *
                         (2 * pow(m_rho, 4) - 3 * pow(s, 2) +
                          pow(m_rho, 2) * (5 * s - 7 * t) + pow(t, 2)) +
                     (pow(m_rho, 2) - s - t) * (3 * pow(m_rho, 4) - s * t -
                                                pow(m_rho, 2) * (3 * s + t))) +
            2 * pow(m_rho, 2) *
                (pow(pion_mass, 4) * (1 + 4 * C4 * pow(m_rho, 2)) +
                 pow(m_rho, 4) * (-1 + 4 * C4 * s) + s * t -
                 pow(pion_mass, 2) *
                     (4 * C4 * pow(m_rho, 4) + s -
                      2 * pow(m_rho, 2) * (1 + 4 * C4 * s) + t) +
                 pow(m_rho, 2) * (t + s * (1 - 4 * C4 * (s + 2 * t)))))) /
              (-pow(pion_mass, 2) + s))) /
            pow(m_rho, 4))) /
      (16. * Pi *
       (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
        2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)));

  // clang-format on
  return cut_off(gev2_mb * diff_xs / spin_deg_factor);
}

// C13
double
CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi_rho_pi0_rho_mediated(
    const double s, const double m_rho) {
  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  auto t_mandelstam = get_t_range(sqrt(s), pion_mass, m_rho, pion_mass, 0.);
  const double &tmax = t_mandelstam[0];
  const double &tmin = t_mandelstam[1];
  const double spin_deg_factor = 3.0;

  // clang-format off
  const double xs =
      (pow(Const, 2) * pow(ghat, 4) *
       (0. +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta2, 2) *
              (1. * pow(a1_mass, 8) + 1. * pow(pion_mass, 8) -
               2. * pow(pion_mass, 6) * pow(m_rho, 2) +
               1. * pow(pion_mass, 4) * pow(m_rho, 4) +
               pow(a1_mass, 6) *
                   (-4. * pow(pion_mass, 2) + 2. * pow(m_rho, 2) + 2. * s) +
               pow(a1_mass, 2) * pow(pion_mass, 2) *
                   (-2. * pow(m_rho, 4) +
                    pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 2. * s) +
                    2. * pow(m_rho, 2) * s) +
               pow(a1_mass, 4) *
                   (4. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                    pow(pion_mass, 2) * (-4. * pow(m_rho, 2) - 4. * s) -
                    2. * pow(m_rho, 2) * s + 2. * pow(s, 2))) +
          eta1 * eta2 *
              (-2. * pow(a1_mass, 8) - 2. * pow(pion_mass, 8) +
               2. * pow(pion_mass, 4) * pow(m_rho, 4) +
               pow(a1_mass, 6) * (8. * pow(pion_mass, 2) - 4. * s) +
               pow(a1_mass, 2) * pow(pion_mass, 2) *
                   (-4. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s +
                    pow(pion_mass, 2) * (4. * pow(m_rho, 2) + 4. * s)) +
               pow(a1_mass, 4) *
                   (-8. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                    4. * pow(m_rho, 2) * s - 4. * pow(s, 2) +
                    pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 8. * s))) +
          pow(eta1, 2) *
              (1. * pow(a1_mass, 8) + 1. * pow(pion_mass, 8) -
               2. * pow(pion_mass, 6) * pow(m_rho, 2) -
               2. * pow(pion_mass, 2) * pow(m_rho, 4) * s +
               pow(a1_mass, 6) *
                   (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
               pow(pion_mass, 4) *
                   (3. * pow(m_rho, 4) + 2. * pow(m_rho, 2) * s) +
               pow(a1_mass, 4) *
                   (4. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                    pow(pion_mass, 2) * (8. * pow(m_rho, 2) - 4. * s) -
                    4. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
               pow(a1_mass, 2) *
                   (pow(pion_mass, 4) * (-6. * pow(m_rho, 2) - 2. * s) +
                    pow(m_rho, 2) * (2. * pow(m_rho, 2) - 2. * s) * s +
                    pow(pion_mass, 2) *
                        (-4. * pow(m_rho, 4) + 6. * pow(m_rho, 2) * s))))) /
            (1. * pow(a1_mass, 2) - 1. * tmax) +
        (1. * pow(-2. + delta, 2) * pow(pion_mass, 2) *
         (1. * pow(pion_mass, 2) - 0.25 * pow(m_rho, 2))) /
            (1. * pow(pion_mass, 2) - 1. * tmax) -
        (0.25 * pow(-2. + delta, 2) * pow(pion_mass, 2) * tmax) /
            pow(m_rho, 2) +
        0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
            (eta2 * (-1. * pow(a1_mass, 2) + pow(m_rho, 2) - 2. * s) +
             eta1 * (pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
                     2. * pow(m_rho, 2) + s)) *
            tmax +
        0.03125 * pow(eta1 - 1. * eta2, 2) *
            (pow(eta1, 2) * (3. * pow(a1_mass, 4) + 4. * pow(pion_mass, 4) +
                             pow(m_rho, 4) +
                             pow(pion_mass, 2) * (8. * pow(m_rho, 2) - 4. * s) -
                             4. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                             pow(a1_mass, 2) * (-8. * pow(pion_mass, 2) -
                                                4. * pow(m_rho, 2) + 4. * s)) +
             pow(eta2, 2) *
                 (3. * pow(a1_mass, 4) + 4. * pow(pion_mass, 4) +
                  pow(m_rho, 4) +
                  pow(pion_mass, 2) * (-4. * pow(m_rho, 2) - 4. * s) -
                  2. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                  pow(a1_mass, 2) *
                      (-8. * pow(pion_mass, 2) + 4. * pow(m_rho, 2) + 4. * s)) +
             eta1 * eta2 *
                 (-6. * pow(a1_mass, 4) - 8. * pow(pion_mass, 4) +
                  2. * pow(m_rho, 4) +
                  pow(a1_mass, 2) * (16. * pow(pion_mass, 2) - 8. * s) +
                  4. * pow(m_rho, 2) * s - 4. * pow(s, 2) +
                  pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 8. * s))) *
            tmax +
        (2. *
         (0. - 0.25 * pow(m_rho, 4) - 1. * C4 * pow(m_rho, 6) +
          pow(pion_mass, 2) * (0.75 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
          2. * C4 * pow(m_rho, 4) * s +
          pow(delta, 2) *
              (-0.125 * pow(pion_mass, 4) - 0.1875 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (0.0625 * pow(m_rho, 2) + 0.0625 * s) +
               0.1875 * pow(m_rho, 2) * s) +
          delta * (0.25 * pow(pion_mass, 4) + 0.5 * C4 * pow(m_rho, 6) +
                   pow(pion_mass, 2) * (-0.5 * pow(m_rho, 2) +
                                        0.5 * C4 * pow(m_rho, 4) - 0.125 * s) -
                   0.375 * pow(m_rho, 2) * s +
                   pow(m_rho, 4) * (0.5 - 1. * C4 * s))) *
         tmax) /
            pow(m_rho, 4) -
        (0.0625 *
         (0. + 8. * pow(pion_mass, 2) * pow(m_rho, 4) - 12. * pow(m_rho, 6) +
          4. * pow(m_rho, 4) * s +
          delta * pow(m_rho, 2) *
              (-16. * pow(pion_mass, 4) -
               16. * pow(pion_mass, 2) * pow(m_rho, 2) - 4. * pow(m_rho, 4) +
               16. * pow(m_rho, 2) * s + 4. * pow(s, 2)) +
          pow(delta, 2) *
              (8. * pow(pion_mass, 6) + 9. * pow(m_rho, 6) +
               pow(pion_mass, 4) * (4. * pow(m_rho, 2) - 4. * s) -
               13. * pow(m_rho, 4) * s - 5. * pow(m_rho, 2) * pow(s, 2) +
               1. * pow(s, 3) +
               pow(pion_mass, 2) * (-2. * pow(m_rho, 4) +
                                    4. * pow(m_rho, 2) * s - 2. * pow(s, 2)))) *
         tmax) /
            pow(m_rho, 6) -
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (pow(m_rho, 2) *
              (eta1 * (2. * pow(a1_mass, 2) + 2. * pow(m_rho, 2)) +
               eta2 * (-2. * pow(a1_mass, 2) + 2. * pow(pion_mass, 2) -
                       8. * pow(m_rho, 2) + 6. * s)) +
          delta *
              (eta1 * (1. * pow(a1_mass, 4) - 2. * pow(pion_mass, 4) -
                       3. * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) -
                       2. * pow(m_rho, 2) * s + 5.000000000000001 * pow(s, 2) +
                       pow(a1_mass, 2) * (-2. * pow(pion_mass, 2) + 1. * s)) +
               eta2 * (-1. * pow(a1_mass, 4) - 4. * pow(pion_mass, 4) +
                       4. * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (-1. * pow(m_rho, 2) - 2. * s) +
                       1. * pow(m_rho, 2) * s - 1. * pow(s, 2) +
                       pow(a1_mass, 2) * (3. * pow(pion_mass, 2) -
                                          3. * pow(m_rho, 2) + 2. * s)))) *
         tmax) /
            pow(m_rho, 2) -
        (0.5 *
         (pow(m_rho, 6) *
              (-1.5 + C4 * (-12. * pow(pion_mass, 2) + 6. * s) +
               pow(C4, 2) * (-16. * pow(pion_mass, 4) +
                             16. * pow(pion_mass, 2) * s - 4. * pow(s, 2))) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 6) - 2. * pow(pion_mass, 4) * pow(m_rho, 2) +
               0.125 * pow(m_rho, 6) + 0.25 * pow(m_rho, 4) * s -
               0.875 * pow(m_rho, 2) * pow(s, 2) + 0.25 * pow(s, 3) +
               pow(pion_mass, 2) *
                   (-2.5 * pow(m_rho, 4) + 2.25 * pow(m_rho, 2) * s -
                    0.75 * pow(s, 2))) +
          delta * pow(m_rho, 2) *
              (pow(pion_mass, 4) * (1. + 8. * C4 * pow(m_rho, 2)) +
               0.5 * pow(s, 2) + pow(m_rho, 4) * (1.5 - 5. * C4 * s) +
               pow(m_rho, 2) * s * (-0.5 + 1. * C4 * s) +
               pow(pion_mass, 2) * (6. * C4 * pow(m_rho, 4) - 1.5 * s +
                                    pow(m_rho, 2) * (3. - 6. * C4 * s)))) *
         tmax) /
            pow(m_rho, 6) -
        (0.5 *
         (0. - 4. * C4 * pow(m_rho, 8) - 0.5 * pow(m_rho, 4) * s +
          pow(m_rho, 6) * (2. + 2. * C4 * s) +
          pow(delta, 2) *
              (-2. * pow(pion_mass, 6) - 2. * pow(m_rho, 6) +
               0.5 * pow(pion_mass, 4) * s + 2.125 * pow(m_rho, 4) * s +
               1.25 * pow(m_rho, 2) * pow(s, 2) - 0.375 * pow(s, 3) +
               pow(pion_mass, 2) * (1.5 * pow(m_rho, 4) -
                                    1.5 * pow(m_rho, 2) * s + 1. * pow(s, 2))) +
          delta * pow(m_rho, 2) *
              (2. * pow(pion_mass, 4) + 2. * C4 * pow(m_rho, 6) -
               1. * pow(s, 2) + pow(m_rho, 2) * s * (-3. + 1. * C4 * s) +
               pow(m_rho, 4) * (1. + 1. * C4 * s) +
               pow(pion_mass, 2) *
                   (1. * s + pow(m_rho, 2) * (1. - 2. * C4 * s)))) *
         tmax) /
            pow(m_rho, 6) +
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (delta *
              (eta1 * (3. * pow(a1_mass, 4) + 6. * pow(pion_mass, 4) +
                       pow(m_rho, 4) +
                       pow(pion_mass, 2) * (18. * pow(m_rho, 2) - 12. * s) -
                       8. * pow(m_rho, 2) * s + 7. * pow(s, 2) +
                       pow(a1_mass, 2) * (-10. * pow(pion_mass, 2) -
                                          4. * pow(m_rho, 2) + 5. * s)) +
               eta2 * (-3. * pow(a1_mass, 4) - 12. * pow(pion_mass, 4) +
                       2. * pow(m_rho, 4) +
                       pow(a1_mass, 2) * (11. * pow(pion_mass, 2) -
                                          3. * pow(m_rho, 2) - 2. * s) +
                       5. * pow(m_rho, 2) * s - 3. * pow(s, 2) +
                       pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 6. * s))) +
          pow(m_rho, 2) *
              (eta1 * (-8. * C4 * pow(a1_mass, 4) -
                       32. * C4 * pow(pion_mass, 4) - 6. * pow(m_rho, 2) +
                       pow(a1_mass, 2) *
                           (6. + C4 * (32. * pow(pion_mass, 2) +
                                       8. * pow(m_rho, 2) - 16. * s)) +
                       4. * s + 16. * C4 * pow(m_rho, 2) * s -
                       8. * C4 * pow(s, 2) +
                       pow(pion_mass, 2) *
                           (-12. - 32. * C4 * pow(m_rho, 2) + 32. * C4 * s)) +
               eta2 *
                   (8. * C4 * pow(a1_mass, 4) + 32. * C4 * pow(pion_mass, 4) -
                    4. * pow(m_rho, 2) - 2. * s + 8. * C4 * pow(s, 2) +
                    pow(pion_mass, 2) *
                        (10. - 16. * C4 * pow(m_rho, 2) - 32. * C4 * s) +
                    pow(a1_mass, 2) *
                        (-6. + C4 * (-32. * pow(pion_mass, 2) +
                                     8. * pow(m_rho, 2) + 16. * s))))) *
         tmax) /
            pow(m_rho, 2) +
        0.0625 * (-2. + delta) * pow(eta1 - 1. * eta2, 2) * pow(tmax, 2) +
        (0.1875 *
         (1.3333333333333333 * pow(m_rho, 2) +
          5.333333333333333 * C4 * pow(m_rho, 4) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 2) + 1.3333333333333333 * pow(m_rho, 2) -
               0.3333333333333333 * s) +
          delta *
              (-2. * pow(pion_mass, 2) - 3.3333333333333335 * pow(m_rho, 2) -
               2.6666666666666665 * C4 * pow(m_rho, 4) +
               0.6666666666666666 * s)) *
         pow(tmax, 2)) /
            pow(m_rho, 4) +
        0.03125 * pow(eta1 - 1. * eta2, 3) *
            (eta2 * (-1. * pow(a1_mass, 2) + 2. * pow(pion_mass, 2) -
                     1. * pow(m_rho, 2) - 1. * s) +
             eta1 * (pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                     1. * pow(m_rho, 2) + s)) *
            pow(tmax, 2) -
        (0.375 *
         (0.3333333333333333 * pow(m_rho, 4) -
          1.3333333333333333 * C4 * pow(m_rho, 6) +
          delta * pow(m_rho, 2) *
              (1.3333333333333333 * pow(m_rho, 2) -
               0.6666666666666666 * C4 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (-0.6666666666666666 +
                                    1.3333333333333333 * C4 * pow(m_rho, 2)) -
               0.6666666666666666 * s) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 4) + 0.25 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (-0.3333333333333333 * pow(m_rho, 2) +
                                    0.6666666666666666 * s) -
               0.5833333333333334 * pow(s, 2))) *
         pow(tmax, 2)) /
            pow(m_rho, 6) -
        (0.03125 * (1. * eta1 - 1. * eta2) *
         ((2. * eta1 - 2. * eta2) * pow(m_rho, 2) +
          delta *
              (eta1 * (1. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) + 1. * s) +
               eta2 * (-1. * pow(a1_mass, 2) + 3. * pow(pion_mass, 2) -
                       3. * pow(m_rho, 2) + 2. * s))) *
         pow(tmax, 2)) /
            pow(m_rho, 2) +
        (0.03125 *
         (0. - 4. * pow(m_rho, 4) +
          delta * (16. * pow(m_rho, 4) - 8. * pow(m_rho, 2) * s) +
          pow(delta, 2) *
              (4. * pow(pion_mass, 4) - 3. * pow(m_rho, 4) +
               2. * pow(m_rho, 2) * s - 3. * pow(s, 2) +
               pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 4. * s))) *
         pow(tmax, 2)) /
            pow(m_rho, 6) +
        (0.25 *
         (C4 * pow(m_rho, 6) *
              (-6. - 16. * C4 * pow(pion_mass, 2) + 8. * C4 * s) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 4) - 0.25 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (-1.75 * pow(m_rho, 2) + 0.5 * s) +
               0.5 * pow(m_rho, 2) * s - 0.5 * pow(s, 2)) +
          delta * pow(m_rho, 2) *
              (1. * C4 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (0.5 + 10. * C4 * pow(m_rho, 2)) - 0.5 * s +
               pow(m_rho, 2) * (2.5 - 4. * C4 * s))) *
         pow(tmax, 2)) /
            pow(m_rho, 6) +
        (0.09375 * (1. * eta1 - 1. * eta2) *
         (delta * (eta2 * (-1. * pow(a1_mass, 2) +
                           3.6666666666666665 * pow(pion_mass, 2) -
                           1. * pow(m_rho, 2) - 0.6666666666666666 * s) +
                   eta1 * (1. * pow(a1_mass, 2) -
                           3.3333333333333335 * pow(pion_mass, 2) -
                           1.3333333333333333 * pow(m_rho, 2) +
                           1.6666666666666667 * s)) +
          pow(m_rho, 2) *
              (eta1 * (2. + C4 * (-2.6666666666666665 * pow(a1_mass, 2) +
                                  10.666666666666666 * pow(pion_mass, 2) +
                                  2.6666666666666665 * pow(m_rho, 2) -
                                  5.333333333333333 * s)) +
               eta2 * (-2. + C4 * (2.6666666666666665 * pow(a1_mass, 2) -
                                   10.666666666666666 * pow(pion_mass, 2) +
                                   2.6666666666666665 * pow(m_rho, 2) +
                                   5.333333333333333 * s)))) *
         pow(tmax, 2)) /
            pow(m_rho, 2) +
        0.010416666666666666 * pow(eta1 - 1. * eta2, 4) * pow(tmax, 3) -
        (0.041666666666666664 * delta * (-2. + 1. * delta) * pow(tmax, 3)) /
            pow(m_rho, 4) -
        (0.020833333333333332 * delta * pow(1. * eta1 - 1. * eta2, 2) *
         pow(tmax, 3)) /
            pow(m_rho, 2) -
        (0.16666666666666666 * pow(1. * eta1 - 1. * eta2, 2) *
         (-0.375 * delta + 1. * C4 * pow(m_rho, 2)) * pow(tmax, 3)) /
            pow(m_rho, 2) +
        (0.10416666666666666 * delta *
         (-0.8 * pow(m_rho, 2) +
          delta * (0.4 * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 0.6 * s)) *
         pow(tmax, 3)) /
            pow(m_rho, 6) +
        (0.16666666666666666 * delta *
         (0. - 0.75 * delta * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4) +
          0.625 * delta * s) *
         pow(tmax, 3)) /
            pow(m_rho, 6) -
        (0.041666666666666664 *
         (12. * C4 * delta * pow(m_rho, 4) - 16. * pow(C4, 2) * pow(m_rho, 6) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 2) - 2.5 * pow(m_rho, 2) + 1. * s)) *
         pow(tmax, 3)) /
            pow(m_rho, 6) +
        (0. -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(pion_mass, 4) *
             pow(m_rho, 6) +
         0.25 * pow(2. - 1. * delta, 2) * pow(m_rho, 10) -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(m_rho, 6) *
             pow(s, 2) +
         pow(2. - 1. * delta, 2) * pow(pion_mass, 2) * pow(m_rho, 6) *
             (-1. * pow(m_rho, 2) + 1. * s)) /
            (pow(m_rho, 6) * (-2. * pow(pion_mass, 2) + 1. * s + 1. * tmax)) -
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (2. * pow(m_rho, 2) +
          delta * (1. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                   1. * pow(m_rho, 2) + 1. * s)) *
         (eta2 * (-1. * pow(a1_mass, 6) +
                  pow(pion_mass, 2) * (4. * pow(pion_mass, 2) - 1. * s) * s +
                  pow(a1_mass, 4) *
                      (3. * pow(pion_mass, 2) - 4. * pow(m_rho, 2) + 2. * s) +
                  pow(a1_mass, 2) *
                      (-4. * pow(pion_mass, 4) - 2. * pow(pion_mass, 2) * s +
                       (4. * pow(m_rho, 2) - 1. * s) * s)) +
          eta1 * (1. * pow(a1_mass, 6) + 8. * pow(pion_mass, 6) +
                  pow(pion_mass, 4) * (-4. * pow(m_rho, 2) - 6. * s) +
                  pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 2. * s) * s +
                  pow(a1_mass, 4) *
                      (-2. * pow(pion_mass, 2) + 1. * pow(m_rho, 2) + 1. * s) +
                  s * (2. * pow(m_rho, 4) - 3. * pow(m_rho, 2) * s +
                       1. * pow(s, 2)) +
                  pow(a1_mass, 2) *
                      (-2. * pow(pion_mass, 4) - 2. * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) -
                       2. * pow(m_rho, 2) * s + 5. * pow(s, 2)))) *
         log(abs(-pow(a1_mass, 2) + tmax))) /
            (pow(m_rho, 2) * (pow(a1_mass, 2) - 2. * pow(pion_mass, 2) + s)) +
        (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
         (eta2 *
              (-1. * pow(a1_mass, 6) + pow(pion_mass, 6) -
               1. * pow(pion_mass, 4) * pow(m_rho, 2) +
               pow(a1_mass, 4) * (pow(pion_mass, 2) + pow(m_rho, 2) - 2. * s) +
               pow(a1_mass, 2) *
                   (3. * pow(pion_mass, 4) - 2. * pow(pion_mass, 2) * s)) +
          eta1 *
              (pow(a1_mass, 6) +
               pow(a1_mass, 4) *
                   (-2. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + s) +
               pow(pion_mass, 2) *
                   (-4. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) -
                    1. * pow(m_rho, 2) * s +
                    pow(pion_mass, 2) * (3. * pow(m_rho, 2) + s)) +
               pow(a1_mass, 2) *
                   (pow(pion_mass, 4) + pow(m_rho, 4) - 1. * pow(m_rho, 2) * s +
                    pow(pion_mass, 2) * (pow(m_rho, 2) + 2. * s)))) *
         log(abs(-pow(a1_mass, 2) + tmax))) /
            (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) +
        0.0625 * pow(eta1 - 1. * eta2, 2) *
            (pow(eta1, 2) *
                 (2. * pow(a1_mass, 6) +
                  pow(pion_mass, 4) * (-3. * pow(m_rho, 2) - 1. * s) +
                  pow(m_rho, 2) * (pow(m_rho, 2) - 1. * s) * s +
                  pow(a1_mass, 4) *
                      (-6. * pow(pion_mass, 2) - 3. * pow(m_rho, 2) + 3. * s) +
                  pow(pion_mass, 2) *
                      (-2. * pow(m_rho, 4) + 3. * pow(m_rho, 2) * s) +
                  pow(a1_mass, 2) *
                      (4. * pow(pion_mass, 4) + pow(m_rho, 4) +
                       pow(pion_mass, 2) * (8. * pow(m_rho, 2) - 4. * s) -
                       4. * pow(m_rho, 2) * s + 2. * pow(s, 2))) +
             pow(eta2, 2) *
                 (2. * pow(a1_mass, 6) +
                  pow(a1_mass, 4) *
                      (-6. * pow(pion_mass, 2) + 3. * pow(m_rho, 2) + 3. * s) +
                  pow(pion_mass, 2) *
                      (-1. * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 1. * s) +
                       pow(m_rho, 2) * s) +
                  pow(a1_mass, 2) *
                      (4. * pow(pion_mass, 4) + pow(m_rho, 4) +
                       pow(pion_mass, 2) * (-4. * pow(m_rho, 2) - 4. * s) -
                       2. * pow(m_rho, 2) * s + 2. * pow(s, 2))) +
             eta1 * eta2 *
                 (-4. * pow(a1_mass, 6) +
                  pow(a1_mass, 4) * (12. * pow(pion_mass, 2) - 6. * s) +
                  pow(pion_mass, 2) *
                      (-2. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s +
                       pow(pion_mass, 2) * (2. * pow(m_rho, 2) + 2. * s)) +
                  pow(a1_mass, 2) *
                      (-8. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                       4. * pow(m_rho, 2) * s - 4. * pow(s, 2) +
                       pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 8. * s)))) *
            log(abs(-pow(a1_mass, 2) + tmax)) +
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (delta *
              (eta1 *
                   (3. * pow(a1_mass, 6) + 8. * pow(pion_mass, 6) +
                    pow(pion_mass, 4) * (-12. * pow(m_rho, 2) - 6. * s) +
                    pow(a1_mass, 4) * (-10. * pow(pion_mass, 2) -
                                       4. * pow(m_rho, 2) + 5. * s) +
                    pow(pion_mass, 2) *
                        (-4. * pow(m_rho, 4) + 10. * pow(m_rho, 2) * s -
                         2. * pow(s, 2)) +
                    s * (3. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s +
                         pow(s, 2)) +
                    pow(a1_mass, 2) *
                        (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                         pow(pion_mass, 2) * (18. * pow(m_rho, 2) - 12. * s) -
                         8. * pow(m_rho, 2) * s + 7. * pow(s, 2))) +
               eta2 *
                   (-3. * pow(a1_mass, 6) +
                    pow(a1_mass, 4) * (11. * pow(pion_mass, 2) -
                                       3. * pow(m_rho, 2) - 2. * s) +
                    pow(pion_mass, 2) *
                        (-2. * pow(m_rho, 4) + pow(m_rho, 2) * s -
                         1. * pow(s, 2) +
                         pow(pion_mass, 2) * (-2. * pow(m_rho, 2) + 4. * s)) +
                    pow(a1_mass, 2) *
                        (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                         5. * pow(m_rho, 2) * s - 3. * pow(s, 2) +
                         pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 6. * s)))) +
          pow(m_rho, 2) *
              (eta2 *
                   (8. * C4 * pow(a1_mass, 6) +
                    pow(pion_mass, 2) *
                        ((4. + 8. * C4 * pow(pion_mass, 2)) * pow(m_rho, 2) -
                         2. * s) +
                    pow(a1_mass, 4) *
                        (-6. + C4 * (-32. * pow(pion_mass, 2) +
                                     8. * pow(m_rho, 2) + 16. * s)) +
                    pow(a1_mass, 2) *
                        (32. * C4 * pow(pion_mass, 4) - 4. * pow(m_rho, 2) +
                         pow(pion_mass, 2) *
                             (10. - 16. * C4 * pow(m_rho, 2) - 32. * C4 * s) +
                         s * (-2. + 8. * C4 * s))) +
               eta1 * (-8. * C4 * pow(a1_mass, 6) +
                       pow(pion_mass, 4) * (4. + 24. * C4 * pow(m_rho, 2)) +
                       pow(a1_mass, 4) *
                           (6. + C4 * (32. * pow(pion_mass, 2) +
                                       8. * pow(m_rho, 2) - 16. * s)) +
                       s * (-2. * pow(m_rho, 2) + 2. * s) +
                       pow(pion_mass, 2) *
                           (-4. * s + pow(m_rho, 2) * (8. - 16. * C4 * s)) +
                       pow(a1_mass, 2) * (-32. * C4 * pow(pion_mass, 4) +
                                          s * (4. - 8. * C4 * s) +
                                          pow(m_rho, 2) * (-6. + 16. * C4 * s) +
                                          pow(pion_mass, 2) *
                                              (-12. - 32. * C4 * pow(m_rho, 2) +
                                               32. * C4 * s))))) *
         log(abs(-pow(a1_mass, 2) + tmax))) /
            pow(m_rho, 2) +
        0.5 * pow(-2. + delta, 2) * pow(pion_mass, 2) *
            log(abs(-pow(pion_mass, 2) + tmax)) -
        (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
         (eta2 * (2. * pow(pion_mass, 6) - 2. * pow(pion_mass, 4) * s) +
          eta1 * (-2. * pow(pion_mass, 6) -
                  1. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                  pow(pion_mass, 4) * (pow(m_rho, 2) + 2. * s))) *
         log(abs(-pow(pion_mass, 2) + tmax))) /
            (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) -
        (0.125 *
         (0. - 32. * C4 * pow(pion_mass, 6) * pow(m_rho, 4) -
          8. * pow(m_rho, 8) + 8. * pow(m_rho, 6) * s +
          pow(pion_mass, 4) * pow(m_rho, 4) * (16. + 64. * C4 * s) +
          pow(pion_mass, 2) * pow(m_rho, 4) *
              (24. * pow(m_rho, 2) + s * (-16. - 32. * C4 * s)) +
          pow(delta, 2) * pow(m_rho, 2) *
              (-4. * pow(pion_mass, 6) - 2. * pow(m_rho, 6) +
               2. * pow(m_rho, 4) * s +
               pow(pion_mass, 4) * (4. * pow(m_rho, 2) + 8. * s) +
               pow(pion_mass, 2) *
                   (6. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s -
                    4.000000000000001 * pow(s, 2))) +
          delta * pow(m_rho, 2) *
              (8. * pow(m_rho, 6) +
               pow(pion_mass, 6) * (8. + 16. * C4 * pow(m_rho, 2)) -
               8. * pow(m_rho, 4) * s +
               pow(pion_mass, 4) *
                   (-16. * s +
                    pow(m_rho, 2) * (-15.999999999999996 - 32. * C4 * s)) +
               pow(pion_mass, 2) *
                   (-24. * pow(m_rho, 4) + 8. * pow(s, 2) +
                    pow(m_rho, 2) * s * (16. + 16. * C4 * s)))) *
         log(abs(-pow(pion_mass, 2) + tmax))) /
            (pow(m_rho, 4) * (-1. * pow(pion_mass, 2) + 1. * s)) -
        (0.25 * (1. * eta1 - 1. * eta2) *
         (eta2 * ((2. - 1. * delta) * pow(pion_mass, 6) +
                  pow(pion_mass, 2) * s *
                      ((-12. + 6. * delta) * pow(m_rho, 2) +
                       (6. - 3. * delta) * s) +
                  pow(s, 2) * ((4. - 2. * delta) * pow(m_rho, 2) +
                               (-2. + 1. * delta) * s) +
                  pow(pion_mass, 4) * ((8. - 4. * delta) * pow(m_rho, 2) +
                                       (-6. + 3. * delta) * s)) +
          eta1 * ((-2. + 1. * delta) * pow(pion_mass, 6) +
                  (-2. + 1. * delta) * pow(m_rho, 4) * s +
                  (2. - 1. * delta) * pow(s, 3) +
                  pow(pion_mass, 4) * ((-4. + 2. * delta) * pow(m_rho, 2) +
                                       (6. - 3. * delta) * s) +
                  pow(pion_mass, 2) * ((2. - 1. * delta) * pow(m_rho, 4) +
                                       (4. - 2. * delta) * pow(m_rho, 2) * s +
                                       (-6. + 3. * delta) * pow(s, 2)))) *
         log(abs(-2. * pow(pion_mass, 2) + s + tmax))) /
            (pow(a1_mass, 2) - 2. * pow(pion_mass, 2) + s) +
        (0.125 *
         (0. +
          (32. - 31.999999999999993 * delta + 8. * pow(delta, 2)) *
              pow(pion_mass, 4) * pow(m_rho, 4) -
          2.0000000000000004 * pow(2. - 1. * delta, 2) * pow(m_rho, 8) +
          pow(pion_mass, 2) * pow(m_rho, 4) *
              (8.000000000000002 * pow(2. - 1. * delta, 2) * pow(m_rho, 2) +
               (-32. + 31.999999999999996 * delta - 8. * pow(delta, 2)) * s)) *
         log(abs(-2. * pow(pion_mass, 2) + s + tmax))) /
            (pow(m_rho, 4) * (-1. * pow(pion_mass, 2) + 1. * s)) +
        (0.25 *
         (0. + 8. * pow(pion_mass, 2) * pow(m_rho, 6) + 4. * pow(m_rho, 8) -
          8. * pow(m_rho, 6) * s +
          delta * pow(m_rho, 4) *
              (8. * pow(pion_mass, 4) - 8. * pow(m_rho, 4) +
               pow(pion_mass, 2) * (8. * pow(m_rho, 2) - 16. * s) +
               8. * pow(m_rho, 2) * s + 8. * pow(s, 2)) +
          pow(delta, 2) * pow(m_rho, 4) *
              (-4. * pow(pion_mass, 4) + 3. * pow(m_rho, 4) -
               2. * pow(m_rho, 2) * s - 4. * pow(s, 2) +
               pow(pion_mass, 2) * (-6. * pow(m_rho, 2) + 8. * s))) *
         log(abs(-2. * pow(pion_mass, 2) + s + tmax))) /
            pow(m_rho, 6) -
        (0.5 *
         (0. +
          pow(pion_mass, 2) * (4. * pow(m_rho, 6) - 8. * C4 * pow(m_rho, 8)) -
          4. * pow(m_rho, 6) * s + 8. * C4 * pow(m_rho, 8) * s +
          pow(delta, 2) * pow(m_rho, 4) *
              (-2. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) - 2. * pow(s, 2) +
               pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 4. * s)) +
          delta * pow(m_rho, 4) *
              (4. * pow(pion_mass, 4) +
               pow(pion_mass, 2) *
                   (6. * pow(m_rho, 2) + 4. * C4 * pow(m_rho, 4) - 8. * s) +
               2. * pow(m_rho, 2) * s + 4. * pow(s, 2) +
               pow(m_rho, 4) * (-2. - 4. * C4 * s))) *
         log(abs(-2. * pow(pion_mass, 2) + s + tmax))) /
            pow(m_rho, 6))) /
          (16. * Pi *
           (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
            2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))) -
      (pow(Const, 2) * pow(ghat, 4) *
       (0. +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta2, 2) *
              (1. * pow(a1_mass, 8) + 1. * pow(pion_mass, 8) -
               2. * pow(pion_mass, 6) * pow(m_rho, 2) +
               1. * pow(pion_mass, 4) * pow(m_rho, 4) +
               pow(a1_mass, 6) *
                   (-4. * pow(pion_mass, 2) + 2. * pow(m_rho, 2) + 2. * s) +
               pow(a1_mass, 2) * pow(pion_mass, 2) *
                   (-2. * pow(m_rho, 4) +
                    pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 2. * s) +
                    2. * pow(m_rho, 2) * s) +
               pow(a1_mass, 4) *
                   (4. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                    pow(pion_mass, 2) * (-4. * pow(m_rho, 2) - 4. * s) -
                    2. * pow(m_rho, 2) * s + 2. * pow(s, 2))) +
          eta1 * eta2 *
              (-2. * pow(a1_mass, 8) - 2. * pow(pion_mass, 8) +
               2. * pow(pion_mass, 4) * pow(m_rho, 4) +
               pow(a1_mass, 6) * (8. * pow(pion_mass, 2) - 4. * s) +
               pow(a1_mass, 2) * pow(pion_mass, 2) *
                   (-4. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s +
                    pow(pion_mass, 2) * (4. * pow(m_rho, 2) + 4. * s)) +
               pow(a1_mass, 4) *
                   (-8. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                    4. * pow(m_rho, 2) * s - 4. * pow(s, 2) +
                    pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 8. * s))) +
          pow(eta1, 2) *
              (1. * pow(a1_mass, 8) + 1. * pow(pion_mass, 8) -
               2. * pow(pion_mass, 6) * pow(m_rho, 2) -
               2. * pow(pion_mass, 2) * pow(m_rho, 4) * s +
               pow(a1_mass, 6) *
                   (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
               pow(pion_mass, 4) *
                   (3. * pow(m_rho, 4) + 2. * pow(m_rho, 2) * s) +
               pow(a1_mass, 4) *
                   (4. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                    pow(pion_mass, 2) * (8. * pow(m_rho, 2) - 4. * s) -
                    4. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
               pow(a1_mass, 2) *
                   (pow(pion_mass, 4) * (-6. * pow(m_rho, 2) - 2. * s) +
                    pow(m_rho, 2) * (2. * pow(m_rho, 2) - 2. * s) * s +
                    pow(pion_mass, 2) *
                        (-4. * pow(m_rho, 4) + 6. * pow(m_rho, 2) * s))))) /
            (1. * pow(a1_mass, 2) - 1. * tmin) +
        (1. * pow(-2. + delta, 2) * pow(pion_mass, 2) *
         (1. * pow(pion_mass, 2) - 0.25 * pow(m_rho, 2))) /
            (1. * pow(pion_mass, 2) - 1. * tmin) -
        (0.25 * pow(-2. + delta, 2) * pow(pion_mass, 2) * tmin) /
            pow(m_rho, 2) +
        0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
            (eta2 * (-1. * pow(a1_mass, 2) + pow(m_rho, 2) - 2. * s) +
             eta1 * (pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
                     2. * pow(m_rho, 2) + s)) *
            tmin +
        0.03125 * pow(eta1 - 1. * eta2, 2) *
            (pow(eta1, 2) * (3. * pow(a1_mass, 4) + 4. * pow(pion_mass, 4) +
                             pow(m_rho, 4) +
                             pow(pion_mass, 2) * (8. * pow(m_rho, 2) - 4. * s) -
                             4. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                             pow(a1_mass, 2) * (-8. * pow(pion_mass, 2) -
                                                4. * pow(m_rho, 2) + 4. * s)) +
             pow(eta2, 2) *
                 (3. * pow(a1_mass, 4) + 4. * pow(pion_mass, 4) +
                  pow(m_rho, 4) +
                  pow(pion_mass, 2) * (-4. * pow(m_rho, 2) - 4. * s) -
                  2. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                  pow(a1_mass, 2) *
                      (-8. * pow(pion_mass, 2) + 4. * pow(m_rho, 2) + 4. * s)) +
             eta1 * eta2 *
                 (-6. * pow(a1_mass, 4) - 8. * pow(pion_mass, 4) +
                  2. * pow(m_rho, 4) +
                  pow(a1_mass, 2) * (16. * pow(pion_mass, 2) - 8. * s) +
                  4. * pow(m_rho, 2) * s - 4. * pow(s, 2) +
                  pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 8. * s))) *
            tmin +
        (2. *
         (0. - 0.25 * pow(m_rho, 4) - 1. * C4 * pow(m_rho, 6) +
          pow(pion_mass, 2) * (0.75 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
          2. * C4 * pow(m_rho, 4) * s +
          pow(delta, 2) *
              (-0.125 * pow(pion_mass, 4) - 0.1875 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (0.0625 * pow(m_rho, 2) + 0.0625 * s) +
               0.1875 * pow(m_rho, 2) * s) +
          delta * (0.25 * pow(pion_mass, 4) + 0.5 * C4 * pow(m_rho, 6) +
                   pow(pion_mass, 2) * (-0.5 * pow(m_rho, 2) +
                                        0.5 * C4 * pow(m_rho, 4) - 0.125 * s) -
                   0.375 * pow(m_rho, 2) * s +
                   pow(m_rho, 4) * (0.5 - 1. * C4 * s))) *
         tmin) /
            pow(m_rho, 4) -
        (0.0625 *
         (0. + 8. * pow(pion_mass, 2) * pow(m_rho, 4) - 12. * pow(m_rho, 6) +
          4. * pow(m_rho, 4) * s +
          delta * pow(m_rho, 2) *
              (-16. * pow(pion_mass, 4) -
               16. * pow(pion_mass, 2) * pow(m_rho, 2) - 4. * pow(m_rho, 4) +
               16. * pow(m_rho, 2) * s + 4. * pow(s, 2)) +
          pow(delta, 2) *
              (8. * pow(pion_mass, 6) + 9. * pow(m_rho, 6) +
               pow(pion_mass, 4) * (4. * pow(m_rho, 2) - 4. * s) -
               13. * pow(m_rho, 4) * s - 5. * pow(m_rho, 2) * pow(s, 2) +
               1. * pow(s, 3) +
               pow(pion_mass, 2) * (-2. * pow(m_rho, 4) +
                                    4. * pow(m_rho, 2) * s - 2. * pow(s, 2)))) *
         tmin) /
            pow(m_rho, 6) -
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (pow(m_rho, 2) *
              (eta1 * (2. * pow(a1_mass, 2) + 2. * pow(m_rho, 2)) +
               eta2 * (-2. * pow(a1_mass, 2) + 2. * pow(pion_mass, 2) -
                       8. * pow(m_rho, 2) + 6. * s)) +
          delta *
              (eta1 * (1. * pow(a1_mass, 4) - 2. * pow(pion_mass, 4) -
                       3. * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) -
                       2. * pow(m_rho, 2) * s + 5.000000000000001 * pow(s, 2) +
                       pow(a1_mass, 2) * (-2. * pow(pion_mass, 2) + 1. * s)) +
               eta2 * (-1. * pow(a1_mass, 4) - 4. * pow(pion_mass, 4) +
                       4. * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (-1. * pow(m_rho, 2) - 2. * s) +
                       1. * pow(m_rho, 2) * s - 1. * pow(s, 2) +
                       pow(a1_mass, 2) * (3. * pow(pion_mass, 2) -
                                          3. * pow(m_rho, 2) + 2. * s)))) *
         tmin) /
            pow(m_rho, 2) -
        (0.5 *
         (pow(m_rho, 6) *
              (-1.5 + C4 * (-12. * pow(pion_mass, 2) + 6. * s) +
               pow(C4, 2) * (-16. * pow(pion_mass, 4) +
                             16. * pow(pion_mass, 2) * s - 4. * pow(s, 2))) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 6) - 2. * pow(pion_mass, 4) * pow(m_rho, 2) +
               0.125 * pow(m_rho, 6) + 0.25 * pow(m_rho, 4) * s -
               0.875 * pow(m_rho, 2) * pow(s, 2) + 0.25 * pow(s, 3) +
               pow(pion_mass, 2) *
                   (-2.5 * pow(m_rho, 4) + 2.25 * pow(m_rho, 2) * s -
                    0.75 * pow(s, 2))) +
          delta * pow(m_rho, 2) *
              (pow(pion_mass, 4) * (1. + 8. * C4 * pow(m_rho, 2)) +
               0.5 * pow(s, 2) + pow(m_rho, 4) * (1.5 - 5. * C4 * s) +
               pow(m_rho, 2) * s * (-0.5 + 1. * C4 * s) +
               pow(pion_mass, 2) * (6. * C4 * pow(m_rho, 4) - 1.5 * s +
                                    pow(m_rho, 2) * (3. - 6. * C4 * s)))) *
         tmin) /
            pow(m_rho, 6) -
        (0.5 *
         (0. - 4. * C4 * pow(m_rho, 8) - 0.5 * pow(m_rho, 4) * s +
          pow(m_rho, 6) * (2. + 2. * C4 * s) +
          pow(delta, 2) *
              (-2. * pow(pion_mass, 6) - 2. * pow(m_rho, 6) +
               0.5 * pow(pion_mass, 4) * s + 2.125 * pow(m_rho, 4) * s +
               1.25 * pow(m_rho, 2) * pow(s, 2) - 0.375 * pow(s, 3) +
               pow(pion_mass, 2) * (1.5 * pow(m_rho, 4) -
                                    1.5 * pow(m_rho, 2) * s + 1. * pow(s, 2))) +
          delta * pow(m_rho, 2) *
              (2. * pow(pion_mass, 4) + 2. * C4 * pow(m_rho, 6) -
               1. * pow(s, 2) + pow(m_rho, 2) * s * (-3. + 1. * C4 * s) +
               pow(m_rho, 4) * (1. + 1. * C4 * s) +
               pow(pion_mass, 2) *
                   (1. * s + pow(m_rho, 2) * (1. - 2. * C4 * s)))) *
         tmin) /
            pow(m_rho, 6) +
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (delta *
              (eta1 * (3. * pow(a1_mass, 4) + 6. * pow(pion_mass, 4) +
                       pow(m_rho, 4) +
                       pow(pion_mass, 2) * (18. * pow(m_rho, 2) - 12. * s) -
                       8. * pow(m_rho, 2) * s + 7. * pow(s, 2) +
                       pow(a1_mass, 2) * (-10. * pow(pion_mass, 2) -
                                          4. * pow(m_rho, 2) + 5. * s)) +
               eta2 * (-3. * pow(a1_mass, 4) - 12. * pow(pion_mass, 4) +
                       2. * pow(m_rho, 4) +
                       pow(a1_mass, 2) * (11. * pow(pion_mass, 2) -
                                          3. * pow(m_rho, 2) - 2. * s) +
                       5. * pow(m_rho, 2) * s - 3. * pow(s, 2) +
                       pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 6. * s))) +
          pow(m_rho, 2) *
              (eta1 * (-8. * C4 * pow(a1_mass, 4) -
                       32. * C4 * pow(pion_mass, 4) - 6. * pow(m_rho, 2) +
                       pow(a1_mass, 2) *
                           (6. + C4 * (32. * pow(pion_mass, 2) +
                                       8. * pow(m_rho, 2) - 16. * s)) +
                       4. * s + 16. * C4 * pow(m_rho, 2) * s -
                       8. * C4 * pow(s, 2) +
                       pow(pion_mass, 2) *
                           (-12. - 32. * C4 * pow(m_rho, 2) + 32. * C4 * s)) +
               eta2 *
                   (8. * C4 * pow(a1_mass, 4) + 32. * C4 * pow(pion_mass, 4) -
                    4. * pow(m_rho, 2) - 2. * s + 8. * C4 * pow(s, 2) +
                    pow(pion_mass, 2) *
                        (10. - 16. * C4 * pow(m_rho, 2) - 32. * C4 * s) +
                    pow(a1_mass, 2) *
                        (-6. + C4 * (-32. * pow(pion_mass, 2) +
                                     8. * pow(m_rho, 2) + 16. * s))))) *
         tmin) /
            pow(m_rho, 2) +
        0.0625 * (-2. + delta) * pow(eta1 - 1. * eta2, 2) * pow(tmin, 2) +
        (0.1875 *
         (1.3333333333333333 * pow(m_rho, 2) +
          5.333333333333333 * C4 * pow(m_rho, 4) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 2) + 1.3333333333333333 * pow(m_rho, 2) -
               0.3333333333333333 * s) +
          delta *
              (-2. * pow(pion_mass, 2) - 3.3333333333333335 * pow(m_rho, 2) -
               2.6666666666666665 * C4 * pow(m_rho, 4) +
               0.6666666666666666 * s)) *
         pow(tmin, 2)) /
            pow(m_rho, 4) +
        0.03125 * pow(eta1 - 1. * eta2, 3) *
            (eta2 * (-1. * pow(a1_mass, 2) + 2. * pow(pion_mass, 2) -
                     1. * pow(m_rho, 2) - 1. * s) +
             eta1 * (pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                     1. * pow(m_rho, 2) + s)) *
            pow(tmin, 2) -
        (0.375 *
         (0.3333333333333333 * pow(m_rho, 4) -
          1.3333333333333333 * C4 * pow(m_rho, 6) +
          delta * pow(m_rho, 2) *
              (1.3333333333333333 * pow(m_rho, 2) -
               0.6666666666666666 * C4 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (-0.6666666666666666 +
                                    1.3333333333333333 * C4 * pow(m_rho, 2)) -
               0.6666666666666666 * s) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 4) + 0.25 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (-0.3333333333333333 * pow(m_rho, 2) +
                                    0.6666666666666666 * s) -
               0.5833333333333334 * pow(s, 2))) *
         pow(tmin, 2)) /
            pow(m_rho, 6) -
        (0.03125 * (1. * eta1 - 1. * eta2) *
         ((2. * eta1 - 2. * eta2) * pow(m_rho, 2) +
          delta *
              (eta1 * (1. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) + 1. * s) +
               eta2 * (-1. * pow(a1_mass, 2) + 3. * pow(pion_mass, 2) -
                       3. * pow(m_rho, 2) + 2. * s))) *
         pow(tmin, 2)) /
            pow(m_rho, 2) +
        (0.03125 *
         (0. - 4. * pow(m_rho, 4) +
          delta * (16. * pow(m_rho, 4) - 8. * pow(m_rho, 2) * s) +
          pow(delta, 2) *
              (4. * pow(pion_mass, 4) - 3. * pow(m_rho, 4) +
               2. * pow(m_rho, 2) * s - 3. * pow(s, 2) +
               pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 4. * s))) *
         pow(tmin, 2)) /
            pow(m_rho, 6) +
        (0.25 *
         (C4 * pow(m_rho, 6) *
              (-6. - 16. * C4 * pow(pion_mass, 2) + 8. * C4 * s) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 4) - 0.25 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (-1.75 * pow(m_rho, 2) + 0.5 * s) +
               0.5 * pow(m_rho, 2) * s - 0.5 * pow(s, 2)) +
          delta * pow(m_rho, 2) *
              (1. * C4 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (0.5 + 10. * C4 * pow(m_rho, 2)) - 0.5 * s +
               pow(m_rho, 2) * (2.5 - 4. * C4 * s))) *
         pow(tmin, 2)) /
            pow(m_rho, 6) +
        (0.09375 * (1. * eta1 - 1. * eta2) *
         (delta * (eta2 * (-1. * pow(a1_mass, 2) +
                           3.6666666666666665 * pow(pion_mass, 2) -
                           1. * pow(m_rho, 2) - 0.6666666666666666 * s) +
                   eta1 * (1. * pow(a1_mass, 2) -
                           3.3333333333333335 * pow(pion_mass, 2) -
                           1.3333333333333333 * pow(m_rho, 2) +
                           1.6666666666666667 * s)) +
          pow(m_rho, 2) *
              (eta1 * (2. + C4 * (-2.6666666666666665 * pow(a1_mass, 2) +
                                  10.666666666666666 * pow(pion_mass, 2) +
                                  2.6666666666666665 * pow(m_rho, 2) -
                                  5.333333333333333 * s)) +
               eta2 * (-2. + C4 * (2.6666666666666665 * pow(a1_mass, 2) -
                                   10.666666666666666 * pow(pion_mass, 2) +
                                   2.6666666666666665 * pow(m_rho, 2) +
                                   5.333333333333333 * s)))) *
         pow(tmin, 2)) /
            pow(m_rho, 2) +
        0.010416666666666666 * pow(eta1 - 1. * eta2, 4) * pow(tmin, 3) -
        (0.041666666666666664 * delta * (-2. + 1. * delta) * pow(tmin, 3)) /
            pow(m_rho, 4) -
        (0.020833333333333332 * delta * pow(1. * eta1 - 1. * eta2, 2) *
         pow(tmin, 3)) /
            pow(m_rho, 2) -
        (0.16666666666666666 * pow(1. * eta1 - 1. * eta2, 2) *
         (-0.375 * delta + 1. * C4 * pow(m_rho, 2)) * pow(tmin, 3)) /
            pow(m_rho, 2) +
        (0.10416666666666666 * delta *
         (-0.8 * pow(m_rho, 2) +
          delta * (0.4 * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 0.6 * s)) *
         pow(tmin, 3)) /
            pow(m_rho, 6) +
        (0.16666666666666666 * delta *
         (0. - 0.75 * delta * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4) +
          0.625 * delta * s) *
         pow(tmin, 3)) /
            pow(m_rho, 6) -
        (0.041666666666666664 *
         (12. * C4 * delta * pow(m_rho, 4) - 16. * pow(C4, 2) * pow(m_rho, 6) +
          pow(delta, 2) *
              (1. * pow(pion_mass, 2) - 2.5 * pow(m_rho, 2) + 1. * s)) *
         pow(tmin, 3)) /
            pow(m_rho, 6) +
        (0. -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(pion_mass, 4) *
             pow(m_rho, 6) +
         0.25 * pow(2. - 1. * delta, 2) * pow(m_rho, 10) -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(m_rho, 6) *
             pow(s, 2) +
         pow(2. - 1. * delta, 2) * pow(pion_mass, 2) * pow(m_rho, 6) *
             (-1. * pow(m_rho, 2) + 1. * s)) /
            (pow(m_rho, 6) * (-2. * pow(pion_mass, 2) + 1. * s + 1. * tmin)) -
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (2. * pow(m_rho, 2) +
          delta * (1. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                   1. * pow(m_rho, 2) + 1. * s)) *
         (eta2 * (-1. * pow(a1_mass, 6) +
                  pow(pion_mass, 2) * (4. * pow(pion_mass, 2) - 1. * s) * s +
                  pow(a1_mass, 4) *
                      (3. * pow(pion_mass, 2) - 4. * pow(m_rho, 2) + 2. * s) +
                  pow(a1_mass, 2) *
                      (-4. * pow(pion_mass, 4) - 2. * pow(pion_mass, 2) * s +
                       (4. * pow(m_rho, 2) - 1. * s) * s)) +
          eta1 * (1. * pow(a1_mass, 6) + 8. * pow(pion_mass, 6) +
                  pow(pion_mass, 4) * (-4. * pow(m_rho, 2) - 6. * s) +
                  pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 2. * s) * s +
                  pow(a1_mass, 4) *
                      (-2. * pow(pion_mass, 2) + 1. * pow(m_rho, 2) + 1. * s) +
                  s * (2. * pow(m_rho, 4) - 3. * pow(m_rho, 2) * s +
                       1. * pow(s, 2)) +
                  pow(a1_mass, 2) *
                      (-2. * pow(pion_mass, 4) - 2. * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) -
                       2. * pow(m_rho, 2) * s + 5. * pow(s, 2)))) *
         log(abs(-pow(a1_mass, 2) + tmin))) /
            (pow(m_rho, 2) * (pow(a1_mass, 2) - 2. * pow(pion_mass, 2) + s)) +
        (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
         (eta2 *
              (-1. * pow(a1_mass, 6) + pow(pion_mass, 6) -
               1. * pow(pion_mass, 4) * pow(m_rho, 2) +
               pow(a1_mass, 4) * (pow(pion_mass, 2) + pow(m_rho, 2) - 2. * s) +
               pow(a1_mass, 2) *
                   (3. * pow(pion_mass, 4) - 2. * pow(pion_mass, 2) * s)) +
          eta1 *
              (pow(a1_mass, 6) +
               pow(a1_mass, 4) *
                   (-2. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + s) +
               pow(pion_mass, 2) *
                   (-4. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) -
                    1. * pow(m_rho, 2) * s +
                    pow(pion_mass, 2) * (3. * pow(m_rho, 2) + s)) +
               pow(a1_mass, 2) *
                   (pow(pion_mass, 4) + pow(m_rho, 4) - 1. * pow(m_rho, 2) * s +
                    pow(pion_mass, 2) * (pow(m_rho, 2) + 2. * s)))) *
         log(abs(-pow(a1_mass, 2) + tmin))) /
            (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) +
        0.0625 * pow(eta1 - 1. * eta2, 2) *
            (pow(eta1, 2) *
                 (2. * pow(a1_mass, 6) +
                  pow(pion_mass, 4) * (-3. * pow(m_rho, 2) - 1. * s) +
                  pow(m_rho, 2) * (pow(m_rho, 2) - 1. * s) * s +
                  pow(a1_mass, 4) *
                      (-6. * pow(pion_mass, 2) - 3. * pow(m_rho, 2) + 3. * s) +
                  pow(pion_mass, 2) *
                      (-2. * pow(m_rho, 4) + 3. * pow(m_rho, 2) * s) +
                  pow(a1_mass, 2) *
                      (4. * pow(pion_mass, 4) + pow(m_rho, 4) +
                       pow(pion_mass, 2) * (8. * pow(m_rho, 2) - 4. * s) -
                       4. * pow(m_rho, 2) * s + 2. * pow(s, 2))) +
             pow(eta2, 2) *
                 (2. * pow(a1_mass, 6) +
                  pow(a1_mass, 4) *
                      (-6. * pow(pion_mass, 2) + 3. * pow(m_rho, 2) + 3. * s) +
                  pow(pion_mass, 2) *
                      (-1. * pow(m_rho, 4) +
                       pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 1. * s) +
                       pow(m_rho, 2) * s) +
                  pow(a1_mass, 2) *
                      (4. * pow(pion_mass, 4) + pow(m_rho, 4) +
                       pow(pion_mass, 2) * (-4. * pow(m_rho, 2) - 4. * s) -
                       2. * pow(m_rho, 2) * s + 2. * pow(s, 2))) +
             eta1 * eta2 *
                 (-4. * pow(a1_mass, 6) +
                  pow(a1_mass, 4) * (12. * pow(pion_mass, 2) - 6. * s) +
                  pow(pion_mass, 2) *
                      (-2. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s +
                       pow(pion_mass, 2) * (2. * pow(m_rho, 2) + 2. * s)) +
                  pow(a1_mass, 2) *
                      (-8. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                       4. * pow(m_rho, 2) * s - 4. * pow(s, 2) +
                       pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 8. * s)))) *
            log(abs(-pow(a1_mass, 2) + tmin)) +
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (delta *
              (eta1 *
                   (3. * pow(a1_mass, 6) + 8. * pow(pion_mass, 6) +
                    pow(pion_mass, 4) * (-12. * pow(m_rho, 2) - 6. * s) +
                    pow(a1_mass, 4) * (-10. * pow(pion_mass, 2) -
                                       4. * pow(m_rho, 2) + 5. * s) +
                    pow(pion_mass, 2) *
                        (-4. * pow(m_rho, 4) + 10. * pow(m_rho, 2) * s -
                         2. * pow(s, 2)) +
                    s * (3. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s +
                         pow(s, 2)) +
                    pow(a1_mass, 2) *
                        (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                         pow(pion_mass, 2) * (18. * pow(m_rho, 2) - 12. * s) -
                         8. * pow(m_rho, 2) * s + 7. * pow(s, 2))) +
               eta2 *
                   (-3. * pow(a1_mass, 6) +
                    pow(a1_mass, 4) * (11. * pow(pion_mass, 2) -
                                       3. * pow(m_rho, 2) - 2. * s) +
                    pow(pion_mass, 2) *
                        (-2. * pow(m_rho, 4) + pow(m_rho, 2) * s -
                         1. * pow(s, 2) +
                         pow(pion_mass, 2) * (-2. * pow(m_rho, 2) + 4. * s)) +
                    pow(a1_mass, 2) *
                        (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                         5. * pow(m_rho, 2) * s - 3. * pow(s, 2) +
                         pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 6. * s)))) +
          pow(m_rho, 2) *
              (eta2 *
                   (8. * C4 * pow(a1_mass, 6) +
                    pow(pion_mass, 2) *
                        ((4. + 8. * C4 * pow(pion_mass, 2)) * pow(m_rho, 2) -
                         2. * s) +
                    pow(a1_mass, 4) *
                        (-6. + C4 * (-32. * pow(pion_mass, 2) +
                                     8. * pow(m_rho, 2) + 16. * s)) +
                    pow(a1_mass, 2) *
                        (32. * C4 * pow(pion_mass, 4) - 4. * pow(m_rho, 2) +
                         pow(pion_mass, 2) *
                             (10. - 16. * C4 * pow(m_rho, 2) - 32. * C4 * s) +
                         s * (-2. + 8. * C4 * s))) +
               eta1 * (-8. * C4 * pow(a1_mass, 6) +
                       pow(pion_mass, 4) * (4. + 24. * C4 * pow(m_rho, 2)) +
                       pow(a1_mass, 4) *
                           (6. + C4 * (32. * pow(pion_mass, 2) +
                                       8. * pow(m_rho, 2) - 16. * s)) +
                       s * (-2. * pow(m_rho, 2) + 2. * s) +
                       pow(pion_mass, 2) *
                           (-4. * s + pow(m_rho, 2) * (8. - 16. * C4 * s)) +
                       pow(a1_mass, 2) * (-32. * C4 * pow(pion_mass, 4) +
                                          s * (4. - 8. * C4 * s) +
                                          pow(m_rho, 2) * (-6. + 16. * C4 * s) +
                                          pow(pion_mass, 2) *
                                              (-12. - 32. * C4 * pow(m_rho, 2) +
                                               32. * C4 * s))))) *
         log(abs(-pow(a1_mass, 2) + tmin))) /
            pow(m_rho, 2) +
        0.5 * pow(-2. + delta, 2) * pow(pion_mass, 2) *
            log(abs(-pow(pion_mass, 2) + tmin)) -
        (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
         (eta2 * (2. * pow(pion_mass, 6) - 2. * pow(pion_mass, 4) * s) +
          eta1 * (-2. * pow(pion_mass, 6) -
                  1. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                  pow(pion_mass, 4) * (pow(m_rho, 2) + 2. * s))) *
         log(abs(-pow(pion_mass, 2) + tmin))) /
            (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) -
        (0.125 *
         (0. - 32. * C4 * pow(pion_mass, 6) * pow(m_rho, 4) -
          8. * pow(m_rho, 8) + 8. * pow(m_rho, 6) * s +
          pow(pion_mass, 4) * pow(m_rho, 4) * (16. + 64. * C4 * s) +
          pow(pion_mass, 2) * pow(m_rho, 4) *
              (24. * pow(m_rho, 2) + s * (-16. - 32. * C4 * s)) +
          pow(delta, 2) * pow(m_rho, 2) *
              (-4. * pow(pion_mass, 6) - 2. * pow(m_rho, 6) +
               2. * pow(m_rho, 4) * s +
               pow(pion_mass, 4) * (4. * pow(m_rho, 2) + 8. * s) +
               pow(pion_mass, 2) *
                   (6. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s -
                    4.000000000000001 * pow(s, 2))) +
          delta * pow(m_rho, 2) *
              (8. * pow(m_rho, 6) +
               pow(pion_mass, 6) * (8. + 16. * C4 * pow(m_rho, 2)) -
               8. * pow(m_rho, 4) * s +
               pow(pion_mass, 4) *
                   (-16. * s +
                    pow(m_rho, 2) * (-15.999999999999996 - 32. * C4 * s)) +
               pow(pion_mass, 2) *
                   (-24. * pow(m_rho, 4) + 8. * pow(s, 2) +
                    pow(m_rho, 2) * s * (16. + 16. * C4 * s)))) *
         log(abs(-pow(pion_mass, 2) + tmin))) /
            (pow(m_rho, 4) * (-1. * pow(pion_mass, 2) + 1. * s)) -
        (0.25 * (1. * eta1 - 1. * eta2) *
         (eta2 * ((2. - 1. * delta) * pow(pion_mass, 6) +
                  pow(pion_mass, 2) * s *
                      ((-12. + 6. * delta) * pow(m_rho, 2) +
                       (6. - 3. * delta) * s) +
                  pow(s, 2) * ((4. - 2. * delta) * pow(m_rho, 2) +
                               (-2. + 1. * delta) * s) +
                  pow(pion_mass, 4) * ((8. - 4. * delta) * pow(m_rho, 2) +
                                       (-6. + 3. * delta) * s)) +
          eta1 * ((-2. + 1. * delta) * pow(pion_mass, 6) +
                  (-2. + 1. * delta) * pow(m_rho, 4) * s +
                  (2. - 1. * delta) * pow(s, 3) +
                  pow(pion_mass, 4) * ((-4. + 2. * delta) * pow(m_rho, 2) +
                                       (6. - 3. * delta) * s) +
                  pow(pion_mass, 2) * ((2. - 1. * delta) * pow(m_rho, 4) +
                                       (4. - 2. * delta) * pow(m_rho, 2) * s +
                                       (-6. + 3. * delta) * pow(s, 2)))) *
         log(abs(-2. * pow(pion_mass, 2) + s + tmin))) /
            (pow(a1_mass, 2) - 2. * pow(pion_mass, 2) + s) +
        (0.125 *
         (0. +
          (32. - 31.999999999999993 * delta + 8. * pow(delta, 2)) *
              pow(pion_mass, 4) * pow(m_rho, 4) -
          2.0000000000000004 * pow(2. - 1. * delta, 2) * pow(m_rho, 8) +
          pow(pion_mass, 2) * pow(m_rho, 4) *
              (8.000000000000002 * pow(2. - 1. * delta, 2) * pow(m_rho, 2) +
               (-32. + 31.999999999999996 * delta - 8. * pow(delta, 2)) * s)) *
         log(abs(-2. * pow(pion_mass, 2) + s + tmin))) /
            (pow(m_rho, 4) * (-1. * pow(pion_mass, 2) + 1. * s)) +
        (0.25 *
         (0. + 8. * pow(pion_mass, 2) * pow(m_rho, 6) + 4. * pow(m_rho, 8) -
          8. * pow(m_rho, 6) * s +
          delta * pow(m_rho, 4) *
              (8. * pow(pion_mass, 4) - 8. * pow(m_rho, 4) +
               pow(pion_mass, 2) * (8. * pow(m_rho, 2) - 16. * s) +
               8. * pow(m_rho, 2) * s + 8. * pow(s, 2)) +
          pow(delta, 2) * pow(m_rho, 4) *
              (-4. * pow(pion_mass, 4) + 3. * pow(m_rho, 4) -
               2. * pow(m_rho, 2) * s - 4. * pow(s, 2) +
               pow(pion_mass, 2) * (-6. * pow(m_rho, 2) + 8. * s))) *
         log(abs(-2. * pow(pion_mass, 2) + s + tmin))) /
            pow(m_rho, 6) -
        (0.5 *
         (0. +
          pow(pion_mass, 2) * (4. * pow(m_rho, 6) - 8. * C4 * pow(m_rho, 8)) -
          4. * pow(m_rho, 6) * s + 8. * C4 * pow(m_rho, 8) * s +
          pow(delta, 2) * pow(m_rho, 4) *
              (-2. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) - 2. * pow(s, 2) +
               pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 4. * s)) +
          delta * pow(m_rho, 4) *
              (4. * pow(pion_mass, 4) +
               pow(pion_mass, 2) *
                   (6. * pow(m_rho, 2) + 4. * C4 * pow(m_rho, 4) - 8. * s) +
               2. * pow(m_rho, 2) * s + 4. * pow(s, 2) +
               pow(m_rho, 4) * (-2. - 4. * C4 * s))) *
         log(abs(-2. * pow(pion_mass, 2) + s + tmin))) /
            pow(m_rho, 6))) /
          (16. * Pi *
           (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
            2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)));

  // clang-format on
  return cut_off(xs * gev2_mb / spin_deg_factor);
}

double CrosssectionsPhoton<ComputationMethod::Analytic>::
    xs_diff_pi_rho_pi0_rho_mediated(const double s, const double t,
                                    const double m_rho) {
  const double spin_deg_factor = 3.0;

  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  // clang-format off
  const double diff_xs =
      ((pow(Const, 2) * pow(ghat, 4) *
        ((-0.25 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
          (pow(pion_mass, 4) + pow(pow(m_rho, 2) - t, 2) -
           2 * pow(pion_mass, 2) * (pow(m_rho, 2) + t))) /
             (pow(m_rho, 2) * pow(pow(pion_mass, 2) - t, 2)) -
         (0.0625 * (eta1 - eta2) *
          (2 * pow(m_rho, 2) +
           delta * (-2 * pow(pion_mass, 2) - pow(m_rho, 2) + s + t)) *
          (eta1 *
               (8 * pow(pion_mass, 6) + pow(s, 3) +
                2 * pow(m_rho, 4) * (s - t) + 5 * pow(s, 2) * t +
                s * pow(t, 2) + pow(t, 3) +
                2 * pow(pion_mass, 2) * (2 * pow(m_rho, 2) - s - t) * (s + t) -
                pow(m_rho, 2) * (3 * s - t) * (s + t) -
                2 * pow(pion_mass, 4) * (2 * pow(m_rho, 2) + 3 * s + t)) +
           eta2 * (s - t) *
               (4 * pow(pion_mass, 4) + t * (4 * pow(m_rho, 2) - s + t) -
                pow(pion_mass, 2) * (s + 3 * t)))) /
             (pow(m_rho, 2) * (-pow(a1_mass, 2) + t) *
              (-2 * pow(pion_mass, 2) + s + t)) -
         (0.0625 *
          pow(-2. * pow(m_rho, 2) + delta * (2. * pow(pion_mass, 2) +
                                             pow(m_rho, 2) - 1. * s - 1. * t),
              2) *
          (8. * pow(pion_mass, 6) + 4. * pow(m_rho, 6) + pow(s, 3) +
           pow(m_rho, 4) * (-4. * s - 4. * t) +
           pow(pion_mass, 4) * (-4. * pow(m_rho, 2) - 4. * s - 4. * t) +
           3. * pow(s, 2) * t + 3. * s * pow(t, 2) + pow(t, 3) +
           pow(m_rho, 2) * (-3. * pow(s, 2) + 2. * s * t - 3. * pow(t, 2)) +
           pow(pion_mass, 2) *
               (-8. * pow(m_rho, 4) - 2. * pow(s, 2) - 4. * s * t -
                2. * pow(t, 2) + pow(m_rho, 2) * (4. * s + 4. * t)))) /
             (pow(m_rho, 6) *
              pow(2. * pow(pion_mass, 2) - 1. * s - 1. * t, 2)) +
         (0.125 * (-2 + delta) * (eta1 - eta2) *
          (eta2 * (pow(pion_mass, 2) + t) *
               (pow(pion_mass, 4) -
                pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * t) +
                (pow(m_rho, 2) - 2 * s - t) * t) +
           eta1 *
               (-4 * pow(pion_mass, 6) +
                (pow(m_rho, 2) - t) * (pow(m_rho, 2) - s - t) * t +
                pow(pion_mass, 4) * (3 * pow(m_rho, 2) + s + t) -
                pow(pion_mass, 2) * (pow(m_rho, 4) + pow(m_rho, 2) * (s - t) +
                                     2 * t * (-s + t))))) /
             ((-pow(a1_mass, 2) + t) * (-pow(pion_mass, 2) + t)) +
         (0.03125 * pow(eta1 - eta2, 2) *
          (-2 * eta1 * eta2 *
               (pow(pion_mass, 8) -
                pow(pion_mass, 4) *
                    (pow(m_rho, 4) + 2 * (pow(m_rho, 2) + s) * t -
                     4 * pow(t, 2)) +
                pow(t, 2) * (-pow(m_rho, 4) - 2 * pow(m_rho, 2) * s +
                             2 * pow(s, 2) + 2 * s * t + pow(t, 2)) +
                2 * pow(pion_mass, 2) * t *
                    (pow(m_rho, 4) + pow(m_rho, 2) * (s + t) -
                     2 * t * (s + t))) +
           pow(eta2, 2) *
               (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
                pow(pion_mass, 4) * (pow(m_rho, 4) + 4 * pow(m_rho, 2) * t -
                                     2 * (s - 2 * t) * t) +
                pow(t, 2) * (pow(m_rho, 4) + 2 * pow(s, 2) + 2 * s * t +
                             pow(t, 2) + 2 * pow(m_rho, 2) * (-s + t)) -
                2 * pow(pion_mass, 2) * t *
                    (pow(m_rho, 4) - pow(m_rho, 2) * (s - 2 * t) +
                     2 * t * (s + t))) +
           pow(eta1, 2) *
               (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
                pow(pion_mass, 4) *
                    (3 * pow(m_rho, 4) + 2 * pow(m_rho, 2) * (s - 3 * t) -
                     2 * (s - 2 * t) * t) +
                t * (-pow(m_rho, 2) + t) *
                    (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                     pow(m_rho, 2) * (2 * s + t)) -
                2 * pow(pion_mass, 2) * (-pow(m_rho, 2) + t) *
                    (2 * t * (s + t) - pow(m_rho, 2) * (s + 2 * t))))) /
             pow(pow(a1_mass, 2) - t, 2) -
         (0.5 *
          (-2. * pow(m_rho, 2) +
           delta * (2. * pow(pion_mass, 2) + pow(m_rho, 2) - 1. * s - 1. * t)) *
          (delta *
               (-1. * pow(pion_mass, 6) - 0.5 * pow(m_rho, 6) -
                0.1875 * pow(s, 3) +
                pow(pion_mass, 2) * (1. * pow(m_rho, 4) +
                                     pow(m_rho, 2) * (-0.625 * s - 0.375 * t) +
                                     s * (0.5 * s + 0.5 * t)) +
                pow(m_rho, 4) * (0.5 * s + 0.5 * t) +
                pow(pion_mass, 4) *
                    (0.5 * pow(m_rho, 2) + 0.25 * s + 0.75 * t) -
                0.4375 * pow(s, 2) * t - 0.3125 * s * pow(t, 2) -
                0.0625 * pow(t, 3) +
                pow(m_rho, 2) *
                    (0.4375 * pow(s, 2) - 0.25 * s * t + 0.3125 * pow(t, 2))) +
           pow(m_rho, 2) *
               (-0.125 * pow(s, 2) + C4 * pow(m_rho, 4) * (1. * s - 1. * t) +
                0.125 * pow(t, 2) +
                pow(pion_mass, 2) * ((0.25 - 1. * C4 * pow(m_rho, 2)) * s +
                                     (-0.25 + 1. * C4 * pow(m_rho, 2)) * t) +
                pow(m_rho, 2) * (-0.5 * s + 0.5 * C4 * pow(s, 2) +
                                 t * (0.5 - 0.5 * C4 * t))))) /
             (pow(m_rho, 6) * (1. * pow(pion_mass, 2) - 0.5 * s - 0.5 * t)) +
         (pow(delta, 2) *
              (-0.5 * pow(pion_mass, 6) - 0.0625 * pow(m_rho, 6) +
               pow(m_rho, 4) * (-0.125 * s - 0.125 * t) +
               pow(pion_mass, 4) * (1. * pow(m_rho, 2) + 0.5 * t) +
               s * (-0.125 * pow(s, 2) - 0.25 * s * t - 0.125 * pow(t, 2)) +
               pow(pion_mass, 2) * (1.25 * pow(m_rho, 4) + 0.375 * pow(s, 2) +
                                    pow(m_rho, 2) * (-1.125 * s - 0.875 * t) +
                                    0.25 * s * t - 0.125 * pow(t, 2)) +
               pow(m_rho, 2) *
                   (0.4375 * pow(s, 2) + 0.25 * s * t + 0.3125 * pow(t, 2))) +
          pow(m_rho, 6) *
              (0.75 +
               C4 * (8. * C4 * pow(pion_mass, 4) + 2. * C4 * pow(s, 2) +
                     pow(pion_mass, 2) * (6. - 8. * C4 * s - 8. * C4 * t) +
                     t * (-3. + 2. * C4 * t) + s * (-3. + 4. * C4 * t))) +
          delta * pow(m_rho, 2) *
              (pow(pion_mass, 4) * (-0.5 - 4. * C4 * pow(m_rho, 2)) +
               s * (-0.25 * s - 0.25 * t) +
               pow(m_rho, 4) * (-0.75 + 2.5 * C4 * s + 0.5 * C4 * t) +
               pow(m_rho, 2) *
                   (-0.5 * C4 * pow(s, 2) + s * (0.25 - 2. * C4 * t) +
                    t * (1.25 - 1.5 * C4 * t)) +
               pow(pion_mass, 2) *
                   (-3. * C4 * pow(m_rho, 4) + 0.75 * s + 0.25 * t +
                    pow(m_rho, 2) * (-1.5 + 3. * C4 * s + 5. * C4 * t)))) /
             pow(m_rho, 6) +
         (2 *
          ((0.0625 * (-2. + delta) *
            (-2. * pow(m_rho, 2) + delta * (2. * pow(pion_mass, 2) +
                                            pow(m_rho, 2) - 1. * s - 1. * t)) *
            (2. * pow(pion_mass, 6) + 1. * pow(m_rho, 6) +
             pow(pion_mass, 4) * (-3. * pow(m_rho, 2) - 2. * t) +
             pow(m_rho, 4) * (-1.5 * s - 1.5 * t) +
             pow(m_rho, 2) * s * (0.5 * s + 0.5 * t) +
             pow(pion_mass, 2) * (-1. * pow(m_rho, 4) - 0.5 * pow(s, 2) +
                                  pow(m_rho, 2) * (2.5 * s - 0.5 * t) -
                                  1. * s * t - 0.5 * pow(t, 2)) +
             t * (0.5 * pow(s, 2) + 1. * s * t + 0.5 * pow(t, 2)))) /
               ((pow(pion_mass, 2) - 1. * t) *
                (1. * pow(pion_mass, 2) - 0.5 * s - 0.5 * t)) +
           (0.0625 * (-2 + delta) *
            (6 * delta * pow(pion_mass, 6) + delta * s * t * (s + t) +
             pow(m_rho, 6) * (-2 + 3 * delta + 8 * C4 * t) -
             pow(pion_mass, 4) *
                 ((-2 + 9 * delta) * pow(m_rho, 2) - 8 * C4 * pow(m_rho, 4) +
                  delta * (s + 9 * t)) -
             2 * pow(m_rho, 4) *
                 (t * (-1 + 3 * delta + 4 * C4 * t) +
                  s * (-1 + 2 * delta + 8 * C4 * t)) -
             pow(pion_mass, 2) *
                 (8 * C4 * pow(m_rho, 6) +
                  2 * pow(m_rho, 4) * (-2 + delta - 8 * C4 * t) +
                  pow(m_rho, 2) * ((2 - 7 * delta) * s + (2 + 5 * delta) * t) +
                  delta * (pow(s, 2) - 3 * pow(t, 2))) +
             pow(m_rho, 2) * (2 * s * t + delta * (pow(s, 2) + 3 * s * t +
                                                   3 * pow(t, 2))))) /
               (-pow(pion_mass, 2) + t))) /
             pow(m_rho, 4) +
         (0.0625 * (eta1 - eta2) *
          (-(eta2 *
             (-2 * pow(pion_mass, 4) *
                  (4 * C4 * pow(m_rho, 2) * (pow(m_rho, 2) + 4 * t) -
                   delta * (pow(m_rho, 2) - 2 * s + 6 * t)) +
              pow(pion_mass, 2) *
                  (2 * pow(m_rho, 4) * (-2 + delta + 8 * C4 * t) +
                   delta * (pow(s, 2) - 6 * s * t - 11 * pow(t, 2)) +
                   pow(m_rho, 2) * (-((-2 + delta) * s) + (-10 + delta) * t +
                                    32 * C4 * t * (s + t))) +
              t * (-2 * pow(m_rho, 4) * (-2 + delta + 4 * C4 * t) +
                   delta * (3 * pow(s, 2) + 2 * s * t + 3 * pow(t, 2)) +
                   pow(m_rho, 2) * ((2 - 5 * delta) * s + 3 * (2 + delta) * t -
                                    8 * C4 * pow(s + t, 2))))) +
           eta1 *
               (8 * delta * pow(pion_mass, 6) +
                delta * (pow(s, 3) + 7 * pow(s, 2) * t + 5 * s * pow(t, 2) +
                         3 * pow(t, 3)) -
                2 * pow(m_rho, 2) *
                    ((-1 + 2 * delta) * pow(s, 2) +
                     2 * (-1 + 2 * delta) * s * t +
                     (-3 + 2 * delta) * pow(t, 2) +
                     4 * C4 * t * pow(s + t, 2)) +
                pow(pion_mass, 4) *
                    (24 * C4 * pow(m_rho, 4) + 6 * delta * (-s + t) -
                     4 * pow(m_rho, 2) * (-1 + 3 * delta + 8 * C4 * t)) +
                pow(m_rho, 4) * (t * (-6 + delta + 8 * C4 * t) +
                                 s * (-2 + 3 * delta + 16 * C4 * t)) -
                2 * pow(pion_mass, 2) *
                    (delta * (s + t) * (s + 5 * t) -
                     pow(m_rho, 2) * (-2 * s + 5 * delta * s - 6 * t +
                                      9 * delta * t + 16 * C4 * t * (s + t)) +
                     2 * pow(m_rho, 4) *
                         (-2 + delta + 4 * C4 * (s + 2 * t)))))) /
             (pow(m_rho, 2) * (-pow(a1_mass, 2) + t)))) /
       (16. * Pi *
        (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
         2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))));

  // clang-format on
  return cut_off(gev2_mb * diff_xs / spin_deg_factor);
}
/*----------------------------------------------------------------------------*/
/* 					Pi + Rho -> Pi + Photon channels
 * mediated by (omega) 						  */
/*----------------------------------------------------------------------------*/
// C14
double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi0_rho0_pi0(
    const double s, const double m_rho) {
  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  auto t_mandelstam = get_t_range(sqrt(s), pion_mass, m_rho, pion_mass, 0.0);
  const double tmin = t_mandelstam[1];
  const double tmax = t_mandelstam[0];
  const double spin_deg_factor = 3.0;

  // clang-format off
  const double xs =
      ((pow(Const, 2) * pow(g_POR, 4) *
        ((-0.125 * pow(omega_mass, 8) - 0.125 * pow(pion_mass, 8) +
          0.25 * pow(pion_mass, 6) * pow(m_rho, 2) -
          0.125 * pow(pion_mass, 4) * pow(m_rho, 4) +
          pow(omega_mass, 6) *
              (0.5 * pow(pion_mass, 2) + 0.25 * pow(m_rho, 2) - 0.25 * s) +
          pow(omega_mass, 2) * pow(pion_mass, 2) *
              (0.25 * pow(m_rho, 4) + 0.25 * pow(pion_mass, 2) * s -
               0.25 * pow(m_rho, 2) * s) +
          pow(omega_mass, 4) *
              (-0.5 * pow(pion_mass, 4) - 0.125 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (-0.5 * pow(m_rho, 2) + 0.5 * s) +
               0.25 * pow(m_rho, 2) * s - 0.25 * pow(s, 2))) /
             (1. * pow(omega_mass, 2) - 1. * tmin) -
         0.125 *
             (3. * pow(omega_mass, 4) + 4. * pow(pion_mass, 4) + pow(m_rho, 4) +
              pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) -
              2. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
              pow(omega_mass, 2) *
                  (-8. * pow(pion_mass, 2) - 4. * pow(m_rho, 2) + 4. * s)) *
             tmin -
         0.125 *
             (pow(omega_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
              s) *
             pow(tmin, 2) -
         0.041666666666666664 * pow(tmin, 3) +
         (0.125 * pow(omega_mass, 8) + 0.125 * pow(pion_mass, 8) -
          0.25 * pow(pion_mass, 6) * pow(m_rho, 2) +
          0.125 * pow(pion_mass, 4) * pow(m_rho, 4) +
          pow(omega_mass, 6) *
              (-0.5 * pow(pion_mass, 2) - 0.25 * pow(m_rho, 2) + 0.25 * s) +
          pow(omega_mass, 2) * pow(pion_mass, 2) *
              (-0.25 * pow(m_rho, 4) - 0.25 * pow(pion_mass, 2) * s +
               0.25 * pow(m_rho, 2) * s) +
          pow(omega_mass, 4) *
              (0.5 * pow(pion_mass, 4) + 0.125 * pow(m_rho, 4) +
               pow(pion_mass, 2) * (0.5 * pow(m_rho, 2) - 0.5 * s) -
               0.25 * pow(m_rho, 2) * s + 0.25 * pow(s, 2))) /
             (1. * pow(omega_mass, 2) - 1. * tmax) +
         0.125 *
             (3. * pow(omega_mass, 4) + 4. * pow(pion_mass, 4) + pow(m_rho, 4) +
              pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) -
              2. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
              pow(omega_mass, 2) *
                  (-8. * pow(pion_mass, 2) - 4. * pow(m_rho, 2) + 4. * s)) *
             tmax +
         0.125 *
             (pow(omega_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
              s) *
             pow(tmax, 2) +
         0.041666666666666664 * pow(tmax, 3) -
         0.25 *
             (2. * pow(omega_mass, 6) +
              pow(omega_mass, 4) *
                  (-6. * pow(pion_mass, 2) - 3. * pow(m_rho, 2) + 3. * s) +
              pow(pion_mass, 2) *
                  (-1. * pow(m_rho, 4) - 1. * pow(pion_mass, 2) * s +
                   pow(m_rho, 2) * s) +
              pow(omega_mass, 2) *
                  (4. * pow(pion_mass, 4) + pow(m_rho, 4) +
                   pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) -
                   2. * pow(m_rho, 2) * s + 2. * pow(s, 2))) *
             log(fabs(-1. * pow(omega_mass, 2) + tmin)) -
         (HeavisideTheta(-omega_mass + sqrt(s)) *
          (tmin *
               (0.125 * pow(pion_mass, 8) -
                0.25 * pow(pion_mass, 6) * pow(m_rho, 2) +
                pow(pion_mass, 4) *
                    (0.125 * pow(m_rho, 4) +
                     pow(omega_mass, 2) * (0.25 * pow(m_rho, 2) - 0.5 * s) -
                     0.25 * pow(m_rho, 2) * s + s * (1. * s - 0.125 * tmin)) +
                pow(pion_mass, 2) * s *
                    (1. * pow(omega_mass, 4) - 0.25 * pow(m_rho, 4) +
                     s * (-1.5 * s - 0.75 * tmin) +
                     pow(m_rho, 2) * (1.5 * s + 0.125 * tmin) +
                     pow(omega_mass, 2) *
                         (-1.0000000000000002 * pow(m_rho, 2) + 0.5 * tmin)) +
                s * (-0.25 * pow(omega_mass, 6) +
                     pow(omega_mass, 4) *
                         (0.25 * pow(m_rho, 2) - 0.5 * s - 0.125 * tmin) +
                     pow(omega_mass, 2) *
                         (0.5 * pow(s, 2) - 0.25 * s * tmin +
                          (0.125 * pow(m_rho, 2) - 0.08333333333333333 * tmin) *
                              tmin) +
                     s * (0.125 * pow(m_rho, 4) + 0.375 * pow(s, 2) +
                          pow(m_rho, 2) * (-0.5 * s - 0.25 * tmin) +
                          0.5 * s * tmin +
                          0.16666666666666666 * pow(tmin, 2)))) +
           (-0.25 * pow(omega_mass, 8) * s +
            pow(omega_mass, 6) *
                (1. * pow(pion_mass, 2) + 0.25 * pow(m_rho, 2) - 0.5 * s) * s +
            pow(pion_mass, 4) * s *
                (0.25 * pow(pion_mass, 4) - 0.25 * pow(m_rho, 2) * s) +
            pow(omega_mass, 4) *
                (pow(pion_mass, 4) * (0.25 * pow(m_rho, 2) - 0.5 * s) -
                 1.0000000000000002 * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                 0.5 * pow(s, 3)) +
            pow(omega_mass, 2) *
                (-0.25 * pow(pion_mass, 8) +
                 0.5 * pow(pion_mass, 4) * pow(s, 2) +
                 pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) * pow(s, 2) +
                 (-0.25 * pow(m_rho, 2) + 0.25 * s) * pow(s, 3))) *
               log(fabs(-1. * pow(omega_mass, 2) + tmin)))) /
             pow(pow(omega_mass, 2) - 1. * s, 2) +
         0.25 *
             (2. * pow(omega_mass, 6) +
              pow(omega_mass, 4) *
                  (-6. * pow(pion_mass, 2) - 3. * pow(m_rho, 2) + 3. * s) +
              pow(pion_mass, 2) *
                  (-1. * pow(m_rho, 4) - 1. * pow(pion_mass, 2) * s +
                   pow(m_rho, 2) * s) +
              pow(omega_mass, 2) *
                  (4. * pow(pion_mass, 4) + pow(m_rho, 4) +
                   pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) -
                   2. * pow(m_rho, 2) * s + 2. * pow(s, 2))) *
             log(fabs(-1. * pow(omega_mass, 2) + tmax)) +
         (HeavisideTheta(-omega_mass + sqrt(s)) *
          (tmax *
               (0.125 * pow(pion_mass, 8) -
                0.25 * pow(pion_mass, 6) * pow(m_rho, 2) +
                pow(pion_mass, 4) *
                    (0.125 * pow(m_rho, 4) +
                     pow(omega_mass, 2) * (0.25 * pow(m_rho, 2) - 0.5 * s) -
                     0.25 * pow(m_rho, 2) * s + s * (1. * s - 0.125 * tmax)) +
                pow(pion_mass, 2) * s *
                    (1. * pow(omega_mass, 4) - 0.25 * pow(m_rho, 4) +
                     s * (-1.5 * s - 0.75 * tmax) +
                     pow(m_rho, 2) * (1.5 * s + 0.125 * tmax) +
                     pow(omega_mass, 2) *
                         (-1.0000000000000002 * pow(m_rho, 2) + 0.5 * tmax)) +
                s * (-0.25 * pow(omega_mass, 6) +
                     pow(omega_mass, 4) *
                         (0.25 * pow(m_rho, 2) - 0.5 * s - 0.125 * tmax) +
                     pow(omega_mass, 2) *
                         (0.5 * pow(s, 2) - 0.25 * s * tmax +
                          (0.125 * pow(m_rho, 2) - 0.08333333333333333 * tmax) *
                              tmax) +
                     s * (0.125 * pow(m_rho, 4) + 0.375 * pow(s, 2) +
                          pow(m_rho, 2) * (-0.5 * s - 0.25 * tmax) +
                          0.5 * s * tmax +
                          0.16666666666666666 * pow(tmax, 2)))) +
           (-0.25 * pow(omega_mass, 8) * s +
            pow(omega_mass, 6) *
                (1. * pow(pion_mass, 2) + 0.25 * pow(m_rho, 2) - 0.5 * s) * s +
            pow(pion_mass, 4) * s *
                (0.25 * pow(pion_mass, 4) - 0.25 * pow(m_rho, 2) * s) +
            pow(omega_mass, 4) *
                (pow(pion_mass, 4) * (0.25 * pow(m_rho, 2) - 0.5 * s) -
                 1.0000000000000002 * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                 0.5 * pow(s, 3)) +
            pow(omega_mass, 2) *
                (-0.25 * pow(pion_mass, 8) +
                 0.5 * pow(pion_mass, 4) * pow(s, 2) +
                 pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) * pow(s, 2) +
                 (-0.25 * pow(m_rho, 2) + 0.25 * s) * pow(s, 3))) *
               log(fabs(-1. * pow(omega_mass, 2) + tmax)))) /
             pow(pow(omega_mass, 2) - 1. * s, 2))) /
       (16. * Pi *
        (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
         2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))));

  // clang-format on
  return cut_off(xs * gev2_mb / spin_deg_factor);
}

double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi0_rho0_pi0(
    const double s, const double t, const double m_rho) {
  const double spin_deg_factor = 3.0;

  using std::log;
  using std::pow;

  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  // clang-format off
  double diff_xs =
      ((pow(Const, 2) * pow(g_POR, 4) *
        ((0.125 * (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
                   pow(pion_mass, 4) * (pow(m_rho, 4) - 2 * (s - 2 * t) * t) +
                   pow(t, 2) * (pow(m_rho, 4) + 2 * pow(s, 2) + 2 * s * t +
                                pow(t, 2) - 2 * pow(m_rho, 2) * (s + t)) -
                   2 * pow(pion_mass, 2) * t *
                       (pow(m_rho, 4) + 2 * t * (s + t) -
                        pow(m_rho, 2) * (s + 2 * t)))) /
             pow(pow(omega_mass, 2) - t, 2) +
         (((0.25 * (-pow(omega_mass, 2) + s) *
            (pow(pion_mass, 8) -
             4 * pow(pion_mass, 2) * s * t * (-pow(m_rho, 2) + s + t) +
             pow(pion_mass, 4) * (2 * s * t - pow(m_rho, 2) * (s + t)) +
             s * t *
                 (pow(s, 2) + 3 * s * t + pow(t, 2) -
                  pow(m_rho, 2) * (s + t)))) /
               (-pow(omega_mass, 2) + t) +
           0.125 * (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
                    pow(pion_mass, 4) *
                        (pow(m_rho, 4) + 4 * pow(s, 2) - 2 * s * t) +
                    pow(s, 2) * (pow(m_rho, 4) + pow(s, 2) + 2 * s * t +
                                 2 * pow(t, 2) - 2 * pow(m_rho, 2) * (s + t)) -
                    2 * pow(pion_mass, 2) * s *
                        (pow(m_rho, 4) + 2 * s * (s + t) -
                         pow(m_rho, 2) * (2 * s + t)))) *
          HeavisideTheta(-omega_mass + sqrt(s))) /
             pow(pow(omega_mass, 2) - s, 2))) /
       (16. * Pi *
        (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
         2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))));

  // clang-format on
  return cut_off(gev2_mb * diff_xs / spin_deg_factor);
}

// C15
double
CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi_rho_pi0_omega_mediated(
    const double s, const double m_rho) {
  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  // we need this check here in case of using the summed up cross sections
  if (sqrt(s) < omega_mass)
    return 0.0;

  auto t_mandelstam = get_t_range(sqrt(s), pion_mass, m_rho, pion_mass, 0.);
  const double &tmax = t_mandelstam[0];
  const double &tmin = t_mandelstam[1];
  const double spin_deg_factor = 3.0;

  // clang-format off
  const double xs =
      (0.0024867959858108648 * pow(Const, 2) * pow(g_POR, 4) *
       (pow(pion_mass, 8) * (1. * tmax - 1. * tmin) +
        pow(pion_mass, 6) * pow(m_rho, 2) * (-2. * tmax + 2. * tmin) +
        pow(pion_mass, 4) * (pow(m_rho, 4) * (1. * tmax - 1. * tmin) +
                             s * (4. * s * tmax - 1. * pow(tmax, 2) -
                                  4. * s * tmin + 1. * pow(tmin, 2))) +
        pow(s, 2) *
            (1. * pow(s, 2) * tmax + 1. * s * pow(tmax, 2) +
             0.6666666666666666 * pow(tmax, 3) +
             pow(m_rho, 4) * (1. * tmax - 1. * tmin) - 1. * pow(s, 2) * tmin -
             1. * s * pow(tmin, 2) - 0.6666666666666666 * pow(tmin, 3) +
             pow(m_rho, 2) * (-2. * s * tmax - 1. * pow(tmax, 2) +
                              2. * s * tmin + 1. * pow(tmin, 2))) +
        pow(pion_mass, 2) * s *
            (pow(m_rho, 4) * (-2. * tmax + 2. * tmin) +
             pow(m_rho, 2) * (4. * s * tmax + 1. * pow(tmax, 2) -
                              4. * s * tmin - 1. * pow(tmin, 2)) +
             s * (-4. * s * tmax - 2. * pow(tmax, 2) + 4. * s * tmin +
                  2. * pow(tmin, 2))))) /
      ((pow(pow(omega_mass, 2) - 1. * s, 2) *
        (pow(pion_mass, 4) + pow(m_rho, 4) +
         pow(pion_mass, 2) * (-2. * pow(m_rho, 2) - 2. * s) -
         2. * pow(m_rho, 2) * s + pow(s, 2))));

  // clang-format on
  return cut_off(xs * gev2_mb / spin_deg_factor);
}

double CrosssectionsPhoton<ComputationMethod::Analytic>::
    xs_diff_pi_rho_pi0_omega_mediated(const double s, const double t,
                                      const double m_rho) {
  const double spin_deg_factor = 3.0;

  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  // we need this check here in case of using the summed up cross sections
  if (sqrt(s) < omega_mass)
    return 0.0;

  // clang-format off
  const double diff_xs =
      (0.0024867959858108648 * pow(Const, 2) * pow(g_POR, 4) *
       (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
        pow(pion_mass, 4) * (pow(m_rho, 4) + 4 * pow(s, 2) - 2 * s * t) +
        pow(s, 2) * (pow(m_rho, 4) + pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                     2 * pow(m_rho, 2) * (s + t)) -
        2 * pow(pion_mass, 2) * s *
            (pow(m_rho, 4) + 2 * s * (s + t) - pow(m_rho, 2) * (2 * s + t)))) /
      ((pow(pow(omega_mass, 2) - s, 2) *
        (pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
         2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s))));

  // clang-format on
  return cut_off(gev2_mb * diff_xs / spin_deg_factor);
}

// C16
double
CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi0_rho_pi_omega_mediated(
    const double s, const double m_rho) {
  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  auto t_mandelstam = get_t_range(sqrt(s), pion_mass, m_rho, pion_mass, 0.);
  const double &tmin = t_mandelstam[1];
  const double &tmax = t_mandelstam[0];
  const double spin_deg_factor = 3.0;

  // clang-format off
  const double xs =
      (0.0024868 * pow(Const, 2) * pow(g_POR, 4) *
       ((pow(omega_mass, 8) +
         pow(pion_mass, 4) * pow(pow(pion_mass, 2) - pow(m_rho, 2), 2) -
         2 * pow(omega_mass, 6) * (2 * pow(pion_mass, 2) + pow(m_rho, 2) - s) -
         2 * pow(omega_mass, 2) * pow(pion_mass, 2) *
             (pow(m_rho, 4) + pow(pion_mass, 2) * s - pow(m_rho, 2) * s) +
         pow(omega_mass, 4) * (4 * pow(pion_mass, 4) + pow(m_rho, 4) +
                               4 * pow(pion_mass, 2) * (pow(m_rho, 2) - s) -
                               2 * pow(m_rho, 2) * s + 2 * pow(s, 2))) /
            (pow(omega_mass, 2) - tmax) +
        3 * pow(omega_mass, 4) * tmax -
        8 * pow(omega_mass, 2) * pow(pion_mass, 2) * tmax +
        4 * pow(pion_mass, 4) * tmax -
        4 * pow(omega_mass, 2) * pow(m_rho, 2) * tmax +
        4 * pow(pion_mass, 2) * pow(m_rho, 2) * tmax + pow(m_rho, 4) * tmax +
        4 * pow(omega_mass, 2) * s * tmax - 4 * pow(pion_mass, 2) * s * tmax -
        2 * pow(m_rho, 2) * s * tmax + 2 * pow(s, 2) * tmax +
        pow(omega_mass, 2) * pow(tmax, 2) -
        2 * pow(pion_mass, 2) * pow(tmax, 2) - pow(m_rho, 2) * pow(tmax, 2) +
        s * pow(tmax, 2) + pow(tmax, 3) / 3. -
        (pow(omega_mass, 8) +
         pow(pion_mass, 4) * pow(pow(pion_mass, 2) - pow(m_rho, 2), 2) -
         2 * pow(omega_mass, 6) * (2 * pow(pion_mass, 2) + pow(m_rho, 2) - s) -
         2 * pow(omega_mass, 2) * pow(pion_mass, 2) *
             (pow(m_rho, 4) + pow(pion_mass, 2) * s - pow(m_rho, 2) * s) +
         pow(omega_mass, 4) * (4 * pow(pion_mass, 4) + pow(m_rho, 4) +
                               4 * pow(pion_mass, 2) * (pow(m_rho, 2) - s) -
                               2 * pow(m_rho, 2) * s + 2 * pow(s, 2))) /
            (pow(omega_mass, 2) - tmin) -
        3 * pow(omega_mass, 4) * tmin +
        8 * pow(omega_mass, 2) * pow(pion_mass, 2) * tmin -
        4 * pow(pion_mass, 4) * tmin +
        4 * pow(omega_mass, 2) * pow(m_rho, 2) * tmin -
        4 * pow(pion_mass, 2) * pow(m_rho, 2) * tmin - pow(m_rho, 4) * tmin -
        4 * pow(omega_mass, 2) * s * tmin + 4 * pow(pion_mass, 2) * s * tmin +
        2 * pow(m_rho, 2) * s * tmin - 2 * pow(s, 2) * tmin -
        pow(omega_mass, 2) * pow(tmin, 2) +
        2 * pow(pion_mass, 2) * pow(tmin, 2) + pow(m_rho, 2) * pow(tmin, 2) -
        s * pow(tmin, 2) - pow(tmin, 3) / 3. +
        2 *
            (2 * pow(omega_mass, 6) -
             3 * pow(omega_mass, 4) *
                 (2 * pow(pion_mass, 2) + pow(m_rho, 2) - s) -
             pow(pion_mass, 2) *
                 (pow(m_rho, 4) + pow(pion_mass, 2) * s - pow(m_rho, 2) * s) +
             pow(omega_mass, 2) * (4 * pow(pion_mass, 4) + pow(m_rho, 4) +
                                   4 * pow(pion_mass, 2) * (pow(m_rho, 2) - s) -
                                   2 * pow(m_rho, 2) * s + 2 * pow(s, 2))) *
            log(abs(-pow(omega_mass, 2) + tmax)) -
        2 *
            (2 * pow(omega_mass, 6) -
             3 * pow(omega_mass, 4) *
                 (2 * pow(pion_mass, 2) + pow(m_rho, 2) - s) -
             pow(pion_mass, 2) *
                 (pow(m_rho, 4) + pow(pion_mass, 2) * s - pow(m_rho, 2) * s) +
             pow(omega_mass, 2) * (4 * pow(pion_mass, 4) + pow(m_rho, 4) +
                                   4 * pow(pion_mass, 2) * (pow(m_rho, 2) - s) -
                                   2 * pow(m_rho, 2) * s + 2 * pow(s, 2))) *
            log(abs(-pow(omega_mass, 2) + tmin)))) /
      ((pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
        2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)));

  // clang-format on
  return cut_off(xs * gev2_mb / spin_deg_factor);
}

double CrosssectionsPhoton<ComputationMethod::Analytic>::
    xs_diff_pi0_rho_pi_omega_mediated(const double s, const double t,
                                      const double m_rho) {
  const double spin_deg_factor = 3.0;

  // clang-format off
  const double diff_xs =
      (0.0024867959858108648 * pow(Const, 2) * pow(g_POR, 4) *
       (pow(pion_mass, 8) - 2 * pow(pion_mass, 6) * pow(m_rho, 2) +
        pow(pion_mass, 4) * (pow(m_rho, 4) - 2 * (s - 2 * t) * t) +
        pow(t, 2) * (pow(m_rho, 4) + 2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                     2 * pow(m_rho, 2) * (s + t)) -
        2 * pow(pion_mass, 2) * t *
            (pow(m_rho, 4) + 2 * t * (s + t) - pow(m_rho, 2) * (s + 2 * t)))) /
      ((pow(pion_mass, 4) + pow(pow(m_rho, 2) - s, 2) -
        2 * pow(pion_mass, 2) * (pow(m_rho, 2) + s)) *
       pow(pow(omega_mass, 2) - t, 2));

  // clang-format on
  return cut_off(gev2_mb * diff_xs / spin_deg_factor);
}

/*----------------------------------------------------------------------------*/
/*				 Pi + Rho -> Pi + Photon channels mediated by
 * (Pi, Rho, a1)         */
/*  			 and Omega summed
 */
/*----------------------------------------------------------------------------*/

// C12 + C16
double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi0_rho_pi(
    const double s, const double m_rho) {
  return cut_off(xs_pi0_rho_pi_rho_mediated(s, m_rho) +
                 xs_pi0_rho_pi_omega_mediated(s, m_rho));
}

double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi0_rho_pi(
    const double s, const double t, const double m_rho) {
  return cut_off(xs_diff_pi0_rho_pi_rho_mediated(s, t, m_rho) +
                 xs_diff_pi0_rho_pi_omega_mediated(s, t, m_rho));
}

// C13 + C15
double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi_rho_pi0(
    const double s, const double m_rho) {
  return cut_off(xs_pi_rho_pi0_rho_mediated(s, m_rho) +
                 xs_pi_rho_pi0_omega_mediated(s, m_rho));
}

double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi_rho_pi0(
    const double s, const double t, const double m_rho) {
  return cut_off(xs_diff_pi_rho_pi0_rho_mediated(s, t, m_rho) +
                 xs_diff_pi_rho_pi0_omega_mediated(s, t, m_rho));
}

/*----------------------------------------------------------------------------*/
/* 					Pi + Pi -> Rho + Photon channels
 * mediated by (Pi, Rho, a1) 				*/
/*----------------------------------------------------------------------------*/
// C21
double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi_pi_rho0(
    const double s, const double m_rho) {
  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  const double s_sqrt = sqrt(s);
  const double spin_deg_factor = 1.0;

  auto mandelstam_t = get_t_range(s_sqrt, pion_mass, pion_mass, m_rho, 0.);
  double tmax = mandelstam_t[0];
  double tmin = mandelstam_t[1];

  // clang-format off
  const double xs =
      (-(pow(Const, 2) * pow(ghat, 4) *
         (0. +
          (0.03125 * pow(eta1 - 1. * eta2, 2) *
           (eta1 * eta2 *
                (-2. * pow(a1_mass, 8) - 2. * pow(pion_mass, 8) +
                 2. * pow(pion_mass, 4) * pow(m_rho, 4) +
                 pow(a1_mass, 6) * (8. * pow(pion_mass, 2) - 4. * s) +
                 pow(a1_mass, 2) * pow(pion_mass, 2) *
                     (8. * pow(pion_mass, 4) - 8. * pow(m_rho, 4) -
                      4. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s) +
                 pow(a1_mass, 4) *
                     (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                      8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                      4. * pow(s, 2))) +
            pow(eta2, 2) *
                (1. * pow(a1_mass, 8) + 1. * pow(pion_mass, 8) -
                 2. * pow(pion_mass, 6) * pow(m_rho, 2) +
                 1. * pow(pion_mass, 4) * pow(m_rho, 4) +
                 pow(a1_mass, 6) *
                     (-4. * pow(pion_mass, 2) + 2. * pow(m_rho, 2) + 2. * s) +
                 pow(a1_mass, 4) *
                     (6. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                      pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                      2. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                 pow(a1_mass, 2) *
                     (-4. * pow(pion_mass, 6) -
                      2. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                      pow(pion_mass, 4) * (6. * pow(m_rho, 2) + 2. * s))) +
            pow(eta1, 2) *
                (1. * pow(a1_mass, 8) +
                 pow(a1_mass, 6) *
                     (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
                 pow(a1_mass, 4) *
                     (6. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                      pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                      4. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                 pow(a1_mass, 2) *
                     (-4. * pow(pion_mass, 6) +
                      2. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                      pow(m_rho, 2) * (2. * pow(m_rho, 2) - 2. * s) * s +
                      pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 2. * s)) +
                 pow(pion_mass, 2) *
                     (1. * pow(pion_mass, 6) +
                      2. * pow(pion_mass, 4) * pow(m_rho, 2) -
                      2. * pow(m_rho, 6) + 2. * pow(m_rho, 4) * s +
                      pow(pion_mass, 2) *
                          (1. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s))))) /
              (1. * pow(a1_mass, 2) - 1. * tmin) +
          (1. * pow(-2. + delta, 2) * pow(pion_mass, 2) *
           (1. * pow(pion_mass, 2) - 0.25 * pow(m_rho, 2))) /
              (1. * pow(pion_mass, 2) - 1. * tmin) +
          (1. * pow(-2. + delta, 2) * pow(pion_mass, 2) *
           (1. * pow(pion_mass, 2) - 0.25 * pow(m_rho, 2))) /
              (1. * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 1. * s -
               1. * tmin) -
          (0.5 * pow(-2. + delta, 2) * pow(pion_mass, 2) * tmin) /
              pow(m_rho, 2) -
          0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
              (-0.5 * eta2 * pow(a1_mass, 2) + 1. * eta1 * pow(pion_mass, 2) +
               0.5 * eta2 * pow(m_rho, 2) + 0.5 * eta1 * s - 1. * eta2 * s) *
              tmin +
          (0.25 *
           (pow(pion_mass, 2) *
                (12. + 1. * pow(delta, 2) - 16. * C4 * pow(m_rho, 2) +
                 delta * (-8. + 8. * C4 * pow(m_rho, 2))) +
            (-4. - 3. * pow(delta, 2) - 16. * C4 * pow(m_rho, 2) +
             delta * (8. + 8. * C4 * pow(m_rho, 2))) *
                s) *
           tmin) /
              pow(m_rho, 2) -
          0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
              (eta2 * (pow(a1_mass, 2) - 1. * s) +
               eta1 * (-2. * pow(pion_mass, 2) + s)) *
              tmin -
          0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
              (eta2 *
                   (-1. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * s) +
               eta1 * (2. * pow(pion_mass, 2) + s)) *
              tmin +
          0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
              (eta1 * (1. * pow(pion_mass, 2) - 0.5 * s) +
               eta2 * (-0.5 * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
                       0.5 * pow(m_rho, 2) + 1. * s)) *
              tmin +
          (0.25 * (-2. + 1. * delta) *
           (-8. * C4 * pow(m_rho, 4) +
            pow(pion_mass, 2) * (2. + 1. * delta - 8. * C4 * pow(m_rho, 2)) +
            (-2. - 3. * delta) * s +
            pow(m_rho, 2) * (2. + 1. * delta + 16. * C4 * s)) *
           tmin) /
              pow(m_rho, 2) +
          (0.25 *
           (32 * pow(C4, 2) * pow(m_rho, 8) + 2 * pow(delta, 2) * pow(s, 2) +
            8 * C4 * pow(m_rho, 6) * (-6 + delta - 8 * C4 * s) +
            2 * delta * pow(m_rho, 2) * s * (-6 + delta - 8 * C4 * s) +
            pow(m_rho, 4) * (12 - pow(delta, 2) + 8 * C4 * (6 + delta) * s +
                             32 * pow(C4, 2) * pow(s, 2))) *
           tmin) /
              pow(m_rho, 4) -
          (1. * (1. * eta1 - 1. * eta2) *
           (eta2 *
                (0.75 * pow(m_rho, 4) - 0.125 * delta * pow(m_rho, 4) -
                 1. * C4 * pow(m_rho, 6) +
                 pow(a1_mass, 2) *
                     (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                 pow(pion_mass, 2) *
                     (-1. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4)) -
                 0.25 * pow(m_rho, 2) * s - 0.375 * delta * pow(m_rho, 2) * s +
                 2. * C4 * pow(m_rho, 4) * s + 0.25 * delta * pow(s, 2) -
                 1. * C4 * pow(m_rho, 2) * pow(s, 2)) +
            eta1 * (0.5 * pow(m_rho, 4) - 1. * C4 * pow(m_rho, 6) +
                    pow(pion_mass, 2) *
                        (1. * pow(m_rho, 2) - 2. * C4 * pow(m_rho, 4)) +
                    pow(a1_mass, 2) *
                        (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) -
                    0.25 * delta * pow(s, 2) +
                    1. * C4 * pow(m_rho, 2) * pow(s, 2))) *
           tmin) /
              pow(m_rho, 2) +
          0.0625 * pow(eta1 - 1. * eta2, 2) *
              (pow(eta2, 2) * (pow(Gammaa1, 2) * pow(a1_mass, 2) -
                               1. * pow(a1_mass, 4) - 2. * pow(pion_mass, 4) -
                               2. * pow(pion_mass, 2) * pow(m_rho, 2) +
                               2. * pow(m_rho, 4) +
                               pow(a1_mass, 2) * (2. * pow(pion_mass, 2) +
                                                  pow(m_rho, 2) - 1. * s) -
                               3. * pow(m_rho, 2) * s + pow(s, 2)) +
               pow(eta1, 2) * (pow(Gammaa1, 2) * pow(a1_mass, 2) -
                               1. * pow(a1_mass, 4) - 2. * pow(pion_mass, 4) -
                               2. * pow(pion_mass, 2) * pow(m_rho, 2) +
                               pow(a1_mass, 2) * (2. * pow(pion_mass, 2) +
                                                  pow(m_rho, 2) - 1. * s) +
                               pow(m_rho, 2) * s + pow(s, 2)) +
               eta1 * eta2 *
                   (-2. * pow(Gammaa1, 2) * pow(a1_mass, 2) +
                    2. * pow(a1_mass, 4) + 4. * pow(pion_mass, 4) +
                    4. * pow(pion_mass, 2) * pow(m_rho, 2) +
                    2. * pow(m_rho, 4) - 2. * pow(s, 2) +
                    pow(a1_mass, 2) * (-4. * pow(pion_mass, 2) -
                                       2. * pow(m_rho, 2) + 2. * s))) *
              tmin +
          0.03125 * pow(eta1 - 1. * eta2, 2) *
              (eta1 * eta2 *
                   (-6. * pow(a1_mass, 4) - 12. * pow(pion_mass, 4) +
                    2. * pow(m_rho, 4) +
                    pow(a1_mass, 2) * (16. * pow(pion_mass, 2) - 8. * s) +
                    8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                    4. * pow(s, 2)) +
               pow(eta1, 2) *
                   (3. * pow(a1_mass, 4) + 6. * pow(pion_mass, 4) +
                    pow(m_rho, 4) +
                    pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                    4. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                    pow(a1_mass, 2) * (-8. * pow(pion_mass, 2) -
                                       4. * pow(m_rho, 2) + 4. * s)) +
               pow(eta2, 2) *
                   (3. * pow(a1_mass, 4) + 6. * pow(pion_mass, 4) +
                    pow(m_rho, 4) +
                    pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                    2. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                    pow(a1_mass, 2) * (-8. * pow(pion_mass, 2) +
                                       4. * pow(m_rho, 2) + 4. * s))) *
              tmin -
          (0.125 * (-2. + 1. * delta) *
           (2. + 1. * delta - 8. * C4 * pow(m_rho, 2)) * pow(tmin, 2)) /
              pow(m_rho, 2) -
          0.5 * pow(1. * eta1 - 1. * eta2, 2) *
              (-0.5 + 1. * C4 * pow(m_rho, 2)) * pow(tmin, 2) -
          (1. *
           (0.5 - 0.125 * pow(delta, 2) - 2. * C4 * pow(m_rho, 2) +
            1. * C4 * delta * pow(m_rho, 2)) *
           pow(tmin, 2)) /
              pow(m_rho, 2) +
          0.0625 * pow(1. * eta1 - 1. * eta2, 4) *
              (1. * pow(pion_mass, 2) + 0.5 * pow(m_rho, 2) - 0.5 * s) *
              pow(tmin, 2) +
          0.03125 * pow(eta1 - 1. * eta2, 3) *
              (eta2 * (-1. * pow(a1_mass, 2) + 2. * pow(pion_mass, 2) -
                       1. * pow(m_rho, 2) - 1. * s) +
               eta1 * (pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                       1. * pow(m_rho, 2) + s)) *
              pow(tmin, 2) +
          0.010416666666666666 * pow(eta1 - 1. * eta2, 4) * pow(tmin, 3) -
          0.020833333333333332 * pow(1. * eta1 - 1. * eta2, 4) * pow(tmin, 3) +
          0.03125 * pow(eta1 - 1. * eta2, 2) *
              (eta1 * eta2 *
                   (2. * pow(Gammaa1, 2) * pow(a1_mass, 2) -
                    6. * pow(a1_mass, 4) - 4. * pow(pion_mass, 4) +
                    2. * pow(m_rho, 4) +
                    pow(a1_mass, 2) * (8. * pow(pion_mass, 2) - 8. * s) +
                    4. * pow(m_rho, 2) * s - 4. * pow(s, 2)) +
               pow(eta1, 2) *
                   (-1. * pow(Gammaa1, 2) * pow(a1_mass, 2) +
                    3. * pow(a1_mass, 4) + 2. * pow(pion_mass, 4) +
                    2. * pow(pion_mass, 2) * pow(m_rho, 2) + pow(m_rho, 4) -
                    4. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                    pow(a1_mass, 2) * (-4. * pow(pion_mass, 2) -
                                       4. * pow(m_rho, 2) + 4. * s)) +
               pow(eta2, 2) *
                   (-1. * pow(Gammaa1, 2) * pow(a1_mass, 2) +
                    3. * pow(a1_mass, 4) + 2. * pow(pion_mass, 4) -
                    2. * pow(pion_mass, 2) * pow(m_rho, 2) + pow(m_rho, 4) -
                    2. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                    pow(a1_mass, 2) * (-4. * pow(pion_mass, 2) +
                                       4. * pow(m_rho, 2) + 4. * s))) *
              (-1. * pow(m_rho, 2) + s + tmin) -
          0.03125 * pow(eta1 - 1. * eta2, 3) *
              (eta2 * (-1. * pow(a1_mass, 2) - 1. * pow(m_rho, 2) - 1. * s) +
               eta1 * (pow(a1_mass, 2) - 1. * pow(m_rho, 2) + s)) *
              pow(-1. * pow(m_rho, 2) + s + tmin, 2) +
          0.010416666666666666 * pow(eta1 - 1. * eta2, 4) *
              pow(-1. * pow(m_rho, 2) + s + tmin, 3) +
          0.25 * (eta1 - 1. * eta2) * (1. * eta1 - 1. * eta2) *
              (-1. + 2. * C4 * pow(m_rho, 2)) *
              pow(pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                      1. * pow(m_rho, 2) + s + tmin,
                  2) -
          (2. * (1. * eta1 - 1. * eta2) *
           (eta2 * (0.375 * pow(m_rho, 4) - 0.0625 * delta * pow(m_rho, 4) -
                    0.5 * C4 * pow(m_rho, 6) +
                    pow(a1_mass, 2) *
                        (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                    pow(pion_mass, 2) *
                        (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) -
                    0.125 * pow(m_rho, 2) * s -
                    0.1875 * delta * pow(m_rho, 2) * s +
                    1. * C4 * pow(m_rho, 4) * s + 0.125 * delta * pow(s, 2) -
                    0.5 * C4 * pow(m_rho, 2) * pow(s, 2)) +
            eta1 * (0.25 * pow(m_rho, 4) - 0.5 * C4 * pow(m_rho, 6) +
                    pow(pion_mass, 2) *
                        (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                    pow(a1_mass, 2) *
                        (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) -
                    0.125 * delta * pow(s, 2) +
                    0.5 * C4 * pow(m_rho, 2) * pow(s, 2))) *
           (1. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
            1. * s + 1. * tmin)) /
              pow(m_rho, 2) +
          (2. * (1. * eta1 - 1. * eta2) * Gammaa1 * a1_mass *
           (eta2 * (0.375 * pow(m_rho, 4) - 0.0625 * delta * pow(m_rho, 4) -
                    0.5 * C4 * pow(m_rho, 6) +
                    pow(a1_mass, 2) *
                        (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                    pow(pion_mass, 2) *
                        (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) -
                    0.125 * pow(m_rho, 2) * s -
                    0.1875 * delta * pow(m_rho, 2) * s +
                    1. * C4 * pow(m_rho, 4) * s + 0.125 * delta * pow(s, 2) -
                    0.5 * C4 * pow(m_rho, 2) * pow(s, 2)) +
            eta1 * (0.25 * pow(m_rho, 4) - 0.5 * C4 * pow(m_rho, 6) +
                    pow(pion_mass, 2) *
                        (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                    pow(a1_mass, 2) *
                        (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) -
                    0.125 * delta * pow(s, 2) +
                    0.5 * C4 * pow(m_rho, 2) * pow(s, 2))) *
           atan((pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
                 s + tmin) /
                (Gammaa1 * a1_mass))) /
              pow(m_rho, 2) +
          (0.25 * (-2. + delta) * (eta1 - 1. * eta2) * Gammaa1 * a1_mass *
           (eta2 *
                (-1. * pow(a1_mass, 6) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (-1. * pow(a1_mass, 2) + 0.5 * pow(m_rho, 2) - 1. * s) +
                 pow(a1_mass, 4) *
                     (2. * pow(pion_mass, 2) + 0.5 * pow(m_rho, 2) - 1. * s) +
                 pow(pion_mass, 4) * (-1.5 * pow(m_rho, 2) + 1. * s) +
                 pow(a1_mass, 2) * pow(pion_mass, 2) *
                     (-1. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + 2. * s)) +
            eta1 * (pow(Gammaa1, 2) * pow(a1_mass, 2) *
                        (1. * pow(pion_mass, 2) + 0.5 * s) +
                    pow(a1_mass, 4) * (1. * pow(pion_mass, 2) + 0.5 * s) +
                    pow(a1_mass, 2) *
                        (-2. * pow(pion_mass, 4) - 1. * pow(pion_mass, 2) * s) +
                    pow(pion_mass, 2) *
                        (1. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                         pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 1.5 * s) +
                         1. * pow(m_rho, 2) * s))) *
           atan((pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
                 s + tmin) /
                (Gammaa1 * a1_mass))) /
              (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
               2. * pow(a1_mass, 2) * pow(pion_mass, 2) + pow(pion_mass, 4)) -
          (0.25 * (-2. + delta) * (eta1 - 1. * eta2) * Gammaa1 * a1_mass *
           (eta2 *
                (-1. * pow(a1_mass, 6) -
                 2. * pow(pion_mass, 4) * pow(m_rho, 2) -
                 1. * pow(pion_mass, 2) * pow(m_rho, 4) +
                 pow(a1_mass, 4) *
                     (2. * pow(pion_mass, 2) + 2. * pow(m_rho, 2) - 1.5 * s) +
                 2.5 * pow(pion_mass, 4) * s +
                 3. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                 0.5 * pow(m_rho, 4) * s - 2. * pow(pion_mass, 2) * pow(s, 2) -
                 1. * pow(m_rho, 2) * pow(s, 2) + 0.5 * pow(s, 3) +
                 pow(Gammaa1, 2) *
                     (-1. * pow(a1_mass, 4) + 0.5 * pow(a1_mass, 2) * s) +
                 pow(a1_mass, 2) *
                     (-1. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                      1. * pow(m_rho, 2) * s +
                      pow(pion_mass, 2) * (-2. * pow(m_rho, 2) + 1. * s))) +
            eta1 *
                (1. * pow(pion_mass, 6) +
                 4. * pow(pion_mass, 4) * pow(m_rho, 2) +
                 1. * pow(pion_mass, 2) * pow(m_rho, 4) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (1. * pow(pion_mass, 2) - 0.5 * s) +
                 pow(a1_mass, 4) * (1. * pow(pion_mass, 2) - 0.5 * s) -
                 4.5 * pow(pion_mass, 4) * s -
                 4. * pow(pion_mass, 2) * pow(m_rho, 2) * s -
                 0.5 * pow(m_rho, 4) * s + 3. * pow(pion_mass, 2) * pow(s, 2) +
                 1. * pow(m_rho, 2) * pow(s, 2) - 0.5 * pow(s, 3) +
                 pow(a1_mass, 2) *
                     (-2. * pow(pion_mass, 4) +
                      (1. * pow(m_rho, 2) - 1. * s) * s +
                      pow(pion_mass, 2) * (-2. * pow(m_rho, 2) + 3. * s)))) *
           atan((pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
                 s + tmin) /
                (Gammaa1 * a1_mass))) /
              (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
               pow(pion_mass, 4) + 2. * pow(pion_mass, 2) * pow(m_rho, 2) +
               pow(m_rho, 4) - 2. * pow(pion_mass, 2) * s -
               2. * pow(m_rho, 2) * s + pow(s, 2) +
               pow(a1_mass, 2) *
                   (-2. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s)) +
          (0.03125 * pow(eta1 - 1. * eta2, 2) *
           (pow(eta2, 2) *
                (pow(Gammaa1, 4) * pow(a1_mass, 4) + pow(a1_mass, 8) +
                 pow(pion_mass, 8) - 2. * pow(pion_mass, 6) * pow(m_rho, 2) +
                 pow(pion_mass, 4) * pow(m_rho, 4) +
                 pow(a1_mass, 6) *
                     (-4. * pow(pion_mass, 2) + 2. * pow(m_rho, 2) + 2. * s) +
                 pow(a1_mass, 4) *
                     (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                      pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                      2. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                 pow(a1_mass, 2) *
                     (-4. * pow(pion_mass, 6) -
                      2. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                      pow(pion_mass, 4) * (6. * pow(m_rho, 2) + 2. * s)) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (-6. * pow(a1_mass, 4) - 6. * pow(pion_mass, 4) -
                      1. * pow(m_rho, 4) +
                      pow(a1_mass, 2) * (12. * pow(pion_mass, 2) -
                                         6. * pow(m_rho, 2) - 6. * s) +
                      2. * pow(m_rho, 2) * s - 2. * pow(s, 2) +
                      pow(pion_mass, 2) * (6. * pow(m_rho, 2) + 4. * s))) +
            eta1 * eta2 *
                (-2. * pow(Gammaa1, 4) * pow(a1_mass, 4) -
                 2. * pow(a1_mass, 8) - 2. * pow(pion_mass, 8) +
                 2. * pow(pion_mass, 4) * pow(m_rho, 4) +
                 pow(a1_mass, 6) * (8. * pow(pion_mass, 2) - 4. * s) +
                 pow(a1_mass, 2) * pow(pion_mass, 2) *
                     (8. * pow(pion_mass, 4) - 8. * pow(m_rho, 4) -
                      4. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s) +
                 pow(a1_mass, 4) *
                     (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                      8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                      4. * pow(s, 2)) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (12. * pow(a1_mass, 4) + 12. * pow(pion_mass, 4) -
                      2. * pow(m_rho, 4) - 8. * pow(pion_mass, 2) * s -
                      4. * pow(m_rho, 2) * s + 4. * pow(s, 2) +
                      pow(a1_mass, 2) * (-24. * pow(pion_mass, 2) + 12. * s))) +
            pow(eta1, 2) *
                (pow(Gammaa1, 4) * pow(a1_mass, 4) + pow(a1_mass, 8) +
                 pow(a1_mass, 6) *
                     (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
                 pow(a1_mass, 4) *
                     (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                      pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                      4. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                 pow(a1_mass, 2) *
                     (-4. * pow(pion_mass, 6) +
                      2. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                      pow(m_rho, 2) * (2. * pow(m_rho, 2) - 2. * s) * s +
                      pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 2. * s)) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (-6. * pow(a1_mass, 4) - 6. * pow(pion_mass, 4) -
                      1. * pow(m_rho, 4) +
                      pow(a1_mass, 2) * (12. * pow(pion_mass, 2) +
                                         6. * pow(m_rho, 2) - 6. * s) +
                      4. * pow(m_rho, 2) * s - 2. * pow(s, 2) +
                      pow(pion_mass, 2) * (-6. * pow(m_rho, 2) + 4. * s)) +
                 pow(pion_mass, 2) *
                     (pow(pion_mass, 6) +
                      2. * pow(pion_mass, 4) * pow(m_rho, 2) -
                      2. * pow(m_rho, 6) + 2. * pow(m_rho, 4) * s +
                      pow(pion_mass, 2) *
                          (pow(m_rho, 4) - 2. * pow(m_rho, 2) * s)))) *
           atan((pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
                 s + tmin) /
                (Gammaa1 * a1_mass))) /
              (Gammaa1 * a1_mass) -
          (0.0625 * pow(eta1 - 1. * eta2, 2) * Gammaa1 * a1_mass *
           (eta1 * eta2 *
                (-2. * pow(Gammaa1, 4) * pow(a1_mass, 4) +
                 14. * pow(a1_mass, 8) + 14. * pow(pion_mass, 8) +
                 28. * pow(pion_mass, 6) * pow(m_rho, 2) +
                 20. * pow(pion_mass, 4) * pow(m_rho, 4) +
                 10. * pow(pion_mass, 2) * pow(m_rho, 6) + 2. * pow(m_rho, 8) -
                 16. * pow(pion_mass, 6) * s -
                 16. * pow(pion_mass, 4) * pow(m_rho, 2) * s -
                 12. * pow(pion_mass, 2) * pow(m_rho, 4) * s -
                 4. * pow(m_rho, 6) * s - 4. * pow(pion_mass, 4) * pow(s, 2) -
                 6. * pow(pion_mass, 2) * pow(m_rho, 2) * pow(s, 2) +
                 8. * pow(pion_mass, 2) * pow(s, 3) +
                 4. * pow(m_rho, 2) * pow(s, 3) - 2. * pow(s, 4) +
                 pow(a1_mass, 6) * (-56. * pow(pion_mass, 2) -
                                    28. * pow(m_rho, 2) + 28. * s) +
                 pow(a1_mass, 4) *
                     (84. * pow(pion_mass, 4) + 24. * pow(m_rho, 4) +
                      pow(pion_mass, 2) * (84. * pow(m_rho, 2) - 72. * s) -
                      36. * pow(m_rho, 2) * s + 12. * pow(s, 2)) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (-4. * pow(a1_mass, 4) - 4. * pow(pion_mass, 4) +
                      pow(a1_mass, 2) * (8. * pow(pion_mass, 2) +
                                         4. * pow(m_rho, 2) - 4. * s) +
                      (4. * pow(m_rho, 2) - 4. * s) * s +
                      pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 8. * s)) +
                 pow(a1_mass, 2) *
                     (-56. * pow(pion_mass, 6) - 10. * pow(m_rho, 6) +
                      18. * pow(m_rho, 4) * s - 6. * pow(m_rho, 2) * pow(s, 2) -
                      2. * pow(s, 3) +
                      pow(pion_mass, 4) * (-84. * pow(m_rho, 2) + 60. * s) +
                      pow(pion_mass, 2) *
                          (-48. * pow(m_rho, 4) + 60. * pow(m_rho, 2) * s -
                           12. * pow(s, 2)))) +
            pow(eta1, 2) *
                (1. * pow(Gammaa1, 4) * pow(a1_mass, 4) - 7. * pow(a1_mass, 8) -
                 7. * pow(pion_mass, 8) -
                 14. * pow(pion_mass, 6) * pow(m_rho, 2) -
                 7. * pow(pion_mass, 4) * pow(m_rho, 4) -
                 2. * pow(pion_mass, 2) * pow(m_rho, 6) +
                 pow(a1_mass, 6) *
                     (28. * pow(pion_mass, 2) + 14. * pow(m_rho, 2) - 14. * s) +
                 8. * pow(pion_mass, 6) * s +
                 11. * pow(pion_mass, 4) * pow(m_rho, 2) * s +
                 6. * pow(pion_mass, 2) * pow(m_rho, 4) * s +
                 1. * pow(m_rho, 6) * s + 2. * pow(pion_mass, 4) * pow(s, 2) -
                 1. * pow(m_rho, 4) * pow(s, 2) -
                 4. * pow(pion_mass, 2) * pow(s, 3) -
                 1. * pow(m_rho, 2) * pow(s, 3) + 1. * pow(s, 4) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (2. * pow(a1_mass, 4) + 2. * pow(pion_mass, 4) +
                      1. * pow(m_rho, 4) +
                      pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 4. * s) -
                      1. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                      pow(a1_mass, 2) * (-4. * pow(pion_mass, 2) -
                                         2. * pow(m_rho, 2) + 2. * s)) +
                 pow(a1_mass, 4) *
                     (-42. * pow(pion_mass, 4) - 9. * pow(m_rho, 4) +
                      21. * pow(m_rho, 2) * s - 6. * pow(s, 2) +
                      pow(pion_mass, 2) * (-42. * pow(m_rho, 2) + 36. * s)) +
                 pow(a1_mass, 2) *
                     (28. * pow(pion_mass, 6) + 2. * pow(m_rho, 6) +
                      pow(pion_mass, 4) * (42. * pow(m_rho, 2) - 30. * s) -
                      9. * pow(m_rho, 4) * s + 6. * pow(m_rho, 2) * pow(s, 2) +
                      1. * pow(s, 3) +
                      pow(pion_mass, 2) *
                          (18. * pow(m_rho, 4) - 36. * pow(m_rho, 2) * s +
                           6. * pow(s, 2)))) +
            pow(eta2, 2) *
                (1. * pow(Gammaa1, 4) * pow(a1_mass, 4) - 7. * pow(a1_mass, 8) -
                 7. * pow(pion_mass, 8) -
                 14. * pow(pion_mass, 6) * pow(m_rho, 2) -
                 1. * pow(pion_mass, 4) * pow(m_rho, 4) +
                 6. * pow(pion_mass, 2) * pow(m_rho, 6) + 2. * pow(m_rho, 8) +
                 pow(a1_mass, 6) *
                     (28. * pow(pion_mass, 2) + 14. * pow(m_rho, 2) - 14. * s) +
                 8. * pow(pion_mass, 6) * s -
                 1. * pow(pion_mass, 4) * pow(m_rho, 2) * s -
                 16. * pow(pion_mass, 2) * pow(m_rho, 4) * s -
                 7. * pow(m_rho, 6) * s + 2. * pow(pion_mass, 4) * pow(s, 2) +
                 14. * pow(pion_mass, 2) * pow(m_rho, 2) * pow(s, 2) +
                 9. * pow(m_rho, 4) * pow(s, 2) -
                 4. * pow(pion_mass, 2) * pow(s, 3) -
                 5. * pow(m_rho, 2) * pow(s, 3) + 1. * pow(s, 4) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (2. * pow(a1_mass, 4) + 2. * pow(pion_mass, 4) +
                      3. * pow(m_rho, 4) +
                      pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 4. * s) -
                      5. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                      pow(a1_mass, 2) * (-4. * pow(pion_mass, 2) -
                                         2. * pow(m_rho, 2) + 2. * s)) +
                 pow(a1_mass, 4) *
                     (-42. * pow(pion_mass, 4) - 3. * pow(m_rho, 4) +
                      9. * pow(m_rho, 2) * s - 6. * pow(s, 2) +
                      pow(pion_mass, 2) * (-42. * pow(m_rho, 2) + 36. * s)) +
                 pow(a1_mass, 2) *
                     (28. * pow(pion_mass, 6) - 4. * pow(m_rho, 6) +
                      pow(pion_mass, 4) * (42. * pow(m_rho, 2) - 30. * s) +
                      9. * pow(m_rho, 4) * s - 6. * pow(m_rho, 2) * pow(s, 2) +
                      1. * pow(s, 3) +
                      pow(pion_mass, 2) *
                          (6. * pow(m_rho, 4) - 12. * pow(m_rho, 2) * s +
                           6. * pow(s, 2))))) *
           atan((pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
                 s + tmin) /
                (Gammaa1 * a1_mass))) /
              (pow(Gammaa1, 2) * pow(a1_mass, 2) + 4. * pow(a1_mass, 4) +
               4. * pow(pion_mass, 4) + 4. * pow(pion_mass, 2) * pow(m_rho, 2) +
               pow(m_rho, 4) - 4. * pow(pion_mass, 2) * s -
               2. * pow(m_rho, 2) * s + pow(s, 2) +
               pow(a1_mass, 2) *
                   (-8. * pow(pion_mass, 2) - 4. * pow(m_rho, 2) + 4. * s)) +
          0.0625 * pow(eta1 - 1. * eta2, 2) *
              (eta1 * eta2 *
                   (-4. * pow(a1_mass, 6) +
                    pow(a1_mass, 4) * (12. * pow(pion_mass, 2) - 6. * s) +
                    pow(pion_mass, 2) *
                        (4. * pow(pion_mass, 4) - 4. * pow(m_rho, 4) -
                         2. * pow(pion_mass, 2) * s + 2. * pow(m_rho, 2) * s) +
                    pow(a1_mass, 2) *
                        (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                         8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                         4. * pow(s, 2))) +
               pow(eta1, 2) *
                   (2. * pow(a1_mass, 6) - 2. * pow(pion_mass, 6) +
                    pow(pion_mass, 2) * pow(m_rho, 2) * s +
                    pow(m_rho, 2) * (pow(m_rho, 2) - 1. * s) * s +
                    pow(pion_mass, 4) * (-3. * pow(m_rho, 2) + s) +
                    pow(a1_mass, 4) * (-6. * pow(pion_mass, 2) -
                                       3. * pow(m_rho, 2) + 3. * s) +
                    pow(a1_mass, 2) *
                        (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                         pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                         4. * pow(m_rho, 2) * s + 2. * pow(s, 2))) +
               pow(eta2, 2) *
                   (2. * pow(a1_mass, 6) - 2. * pow(pion_mass, 6) -
                    1. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                    pow(pion_mass, 4) * (3. * pow(m_rho, 2) + s) +
                    pow(a1_mass, 4) * (-6. * pow(pion_mass, 2) +
                                       3. * pow(m_rho, 2) + 3. * s) +
                    pow(a1_mass, 2) *
                        (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                         pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                         2. * pow(m_rho, 2) * s + 2. * pow(s, 2)))) *
              log(abs(-1. * pow(a1_mass, 2) + tmin)) -
          (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
           (eta2 *
                (-0.5 * pow(a1_mass, 6) - 0.5 * pow(pion_mass, 6) +
                 0.5 * pow(pion_mass, 4) * pow(m_rho, 2) +
                 pow(a1_mass, 4) *
                     (0.5 * pow(pion_mass, 2) + 0.5 * pow(m_rho, 2) - 1. * s) +
                 pow(a1_mass, 2) * pow(pion_mass, 2) *
                     (0.5 * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 1. * s)) +
            eta1 * (pow(a1_mass, 4) * (1. * pow(pion_mass, 2) + 0.5 * s) +
                    pow(pion_mass, 2) *
                        (1. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                         pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 0.5 * s) -
                         0.5 * pow(m_rho, 2) * s) +
                    pow(a1_mass, 2) *
                        (-2. * pow(pion_mass, 4) - 0.5 * pow(m_rho, 2) * s +
                         pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s)))) *
           log(abs(-1. * pow(a1_mass, 2) + tmin))) /
              (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) +
          (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
           (eta2 * (-0.5 * pow(a1_mass, 6) - 0.5 * pow(pion_mass, 6) +
                    pow(a1_mass, 4) *
                        (0.5 * pow(pion_mass, 2) + 0.5 * pow(m_rho, 2)) +
                    pow(pion_mass, 4) * (0.5 * pow(m_rho, 2) - 1. * s) +
                    pow(pion_mass, 2) * (-0.5 * pow(m_rho, 2) + 0.5 * s) * s +
                    pow(a1_mass, 2) *
                        (0.5 * pow(pion_mass, 4) +
                         pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) +
                         (-0.5 * pow(m_rho, 2) + 0.5 * s) * s)) +
            eta1 * (1. * pow(pion_mass, 6) +
                    pow(a1_mass, 4) * (1. * pow(pion_mass, 2) - 0.5 * s) +
                    pow(pion_mass, 2) * (1.5 * pow(m_rho, 2) - 2. * s) * s +
                    (-0.5 * pow(m_rho, 2) + 0.5 * s) * pow(s, 2) +
                    pow(pion_mass, 4) * (-1. * pow(m_rho, 2) + 1.5 * s) +
                    pow(a1_mass, 2) *
                        (-2. * pow(pion_mass, 4) + 0.5 * pow(m_rho, 2) * s +
                         pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s)))) *
           log(abs(-1. * pow(a1_mass, 2) + tmin))) /
              (1. * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
               1. * pow(m_rho, 2) + 1. * s) -
          (0.03125 * pow(eta1 - 1. * eta2, 2) *
           (1. * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
            0.5 * pow(m_rho, 2) + 0.5 * s) *
           (eta1 * eta2 *
                (-2. * pow(a1_mass, 8) +
                 pow(a1_mass, 6) *
                     (8. * pow(pion_mass, 2) + 4. * pow(m_rho, 2) - 4. * s) +
                 pow(a1_mass, 4) *
                     (-12. * pow(pion_mass, 4) - 4. * pow(m_rho, 4) +
                      4. * pow(m_rho, 2) * s +
                      pow(pion_mass, 2) * (-12. * pow(m_rho, 2) + 8. * s)) +
                 pow(pion_mass, 2) *
                     (-2. * pow(pion_mass, 6) -
                      4. * pow(pion_mass, 4) * pow(m_rho, 2) -
                      2. * pow(m_rho, 6) + 4. * pow(m_rho, 4) * s -
                      2. * pow(m_rho, 2) * pow(s, 2) +
                      pow(pion_mass, 2) *
                          (-8. * pow(m_rho, 4) + 8. * pow(m_rho, 2) * s)) +
                 pow(a1_mass, 2) *
                     (8. * pow(pion_mass, 6) + 2. * pow(m_rho, 6) +
                      pow(pion_mass, 4) * (12. * pow(m_rho, 2) - 4. * s) -
                      2. * pow(m_rho, 4) * s - 2. * pow(m_rho, 2) * pow(s, 2) +
                      2. * pow(s, 3) +
                      pow(pion_mass, 2) *
                          (8. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s -
                           4. * pow(s, 2)))) +
            pow(eta2, 2) *
                (pow(a1_mass, 8) +
                 pow(a1_mass, 6) *
                     (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
                 pow(pion_mass, 4) * (pow(pion_mass, 4) +
                                      2. * pow(pion_mass, 2) * pow(m_rho, 2) +
                                      pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
                 pow(a1_mass, 4) *
                     (6. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                      pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) +
                      pow(m_rho, 2) * s) +
                 pow(a1_mass, 2) *
                     (-4. * pow(pion_mass, 6) + 2. * pow(m_rho, 6) -
                      5. * pow(m_rho, 4) * s + 4. * pow(m_rho, 2) * pow(s, 2) -
                      1. * pow(s, 3) +
                      pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 2. * s) +
                      pow(pion_mass, 2) *
                          (2. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s +
                           2. * pow(s, 2)))) +
            pow(eta1, 2) *
                (pow(a1_mass, 8) + pow(pion_mass, 8) +
                 2. * pow(pion_mass, 6) * pow(m_rho, 2) +
                 pow(pion_mass, 2) * pow(m_rho, 2) * s *
                     (-2. * pow(m_rho, 2) + 2. * s) +
                 pow(a1_mass, 6) *
                     (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
                 pow(pion_mass, 4) *
                     (3. * pow(m_rho, 4) - 5. * pow(m_rho, 2) * s) +
                 pow(a1_mass, 4) *
                     (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                      pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                      3. * pow(m_rho, 2) * s) +
                 pow(a1_mass, 2) *
                     (-4. * pow(pion_mass, 6) + pow(m_rho, 4) * s -
                      1. * pow(s, 3) +
                      pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 2. * s) +
                      pow(pion_mass, 2) *
                          (-2. * pow(m_rho, 4) + 4. * pow(m_rho, 2) * s +
                           2. * pow(s, 2))))) *
           log(abs(-1. * pow(a1_mass, 2) + tmin))) /
              (0.25 * pow(Gammaa1, 2) * pow(a1_mass, 2) + 1. * pow(a1_mass, 4) +
               1. * pow(pion_mass, 4) + 1. * pow(pion_mass, 2) * pow(m_rho, 2) +
               0.25 * pow(m_rho, 4) - 1. * pow(pion_mass, 2) * s -
               0.5 * pow(m_rho, 2) * s + 0.25 * pow(s, 2) +
               pow(a1_mass, 2) *
                   (-2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + 1. * s)) -
          (1. * (1. * eta1 - 1. * eta2) *
           (eta2 *
                (pow(a1_mass, 4) *
                     (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                 pow(pion_mass, 2) * pow(m_rho, 2) *
                     (pow(pion_mass, 2) * (0.5 - 1. * C4 * pow(m_rho, 2)) +
                      (-0.25 + 0.125 * delta) * (pow(m_rho, 2) + s)) +
                 pow(a1_mass, 2) *
                     (-1. * C4 * pow(m_rho, 6) +
                      pow(pion_mass, 2) *
                          (-1. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4)) +
                      0.25 * delta * pow(s, 2) +
                      pow(m_rho, 2) * s *
                          (-0.25 - 0.375 * delta - 1. * C4 * s) +
                      pow(m_rho, 4) * (0.75 - 0.125 * delta + 2. * C4 * s))) +
            eta1 *
                (pow(a1_mass, 4) *
                     (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) +
                 pow(a1_mass, 2) *
                     (0.5 * pow(m_rho, 4) - 1. * C4 * pow(m_rho, 6) +
                      pow(pion_mass, 2) *
                          (1. * pow(m_rho, 2) - 2. * C4 * pow(m_rho, 4)) -
                      0.25 * delta * pow(s, 2) +
                      1. * C4 * pow(m_rho, 2) * pow(s, 2)) +
                 pow(m_rho, 2) *
                     (pow(pion_mass, 4) * (-0.5 + 1. * C4 * pow(m_rho, 2)) +
                      s * ((0.25 - 0.125 * delta) * pow(m_rho, 2) +
                           (-0.25 + 0.125 * delta) * s) +
                      pow(pion_mass, 2) *
                          (2. * C4 * pow(m_rho, 4) + (0.5 + 0.25 * delta) * s +
                           pow(m_rho, 2) * (-1. - 2. * C4 * s))))) *
           log(abs(-1. * pow(a1_mass, 2) + tmin))) /
              pow(m_rho, 2) +
          0.5 * pow(-2. + delta, 2) * pow(pion_mass, 2) *
              log(abs(-1. * pow(pion_mass, 2) + tmin)) +
          (0.25 *
           (0. +
            8.000000000000002 * pow(2. - 1. * delta, 2) * pow(pion_mass, 4) *
                pow(m_rho, 2) -
            5.999999999999999 * pow(2. - 1. * delta, 2) * pow(pion_mass, 2) *
                pow(m_rho, 2) * s +
            1. * pow(2. - 1. * delta, 2) * pow(m_rho, 2) * pow(s, 2)) *
           log(abs(-1. * pow(pion_mass, 2) + tmin))) /
              (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
          (0.125 * (-2. + delta) * (eta1 - 1. * eta2) * pow(pion_mass, 2) *
           (0. + eta2 * pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) +
            eta1 * (2. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s +
                    pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 4. * s))) *
           log(abs(-1. * pow(pion_mass, 2) + tmin))) /
              (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) +
          (2. * (-2. + 1. * delta) *
           (0. + (-0.25 + 0.125 * delta) * pow(m_rho, 2) * s +
            pow(pion_mass, 2) * (-2. * C4 * pow(m_rho, 4) - 0.5 * delta * s +
                                 pow(m_rho, 2) * (1. + 2. * C4 * s))) *
           log(abs(-1. * pow(pion_mass, 2) + tmin))) /
              pow(m_rho, 2) -
          (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(pion_mass, 2) *
           (eta1 * (pow(pion_mass, 4) * (-1. * pow(m_rho, 2) + 1. * s) +
                    pow(a1_mass, 2) *
                        (pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) +
                         (-0.5 * pow(m_rho, 2) + 0.5 * s) * s) +
                    pow(pion_mass, 2) *
                        (-1. * pow(m_rho, 4) + 2.5 * pow(m_rho, 2) * s -
                         1.5 * pow(s, 2)) +
                    s * (0.5 * pow(m_rho, 4) - 1. * pow(m_rho, 2) * s +
                         0.5 * pow(s, 2))) +
            eta2 * (0.5 * pow(m_rho, 6) +
                    pow(pion_mass, 4) * (1. * pow(m_rho, 2) - 1. * s) -
                    1.5 * pow(m_rho, 4) * s + 1.5 * pow(m_rho, 2) * pow(s, 2) -
                    0.5 * pow(s, 3) +
                    pow(pion_mass, 2) *
                        (1.5 * pow(m_rho, 4) - 3. * pow(m_rho, 2) * s +
                         1.5 * pow(s, 2)) +
                    pow(a1_mass, 2) *
                        (-0.5 * pow(m_rho, 4) + 1. * pow(m_rho, 2) * s -
                         0.5 * pow(s, 2) +
                         pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s)))) *
           log(abs(-1. * pow(pion_mass, 2) + tmin))) /
              (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
               pow(pion_mass, 4) + 2. * pow(pion_mass, 2) * pow(m_rho, 2) +
               pow(m_rho, 4) - 2. * pow(pion_mass, 2) * s -
               2. * pow(m_rho, 2) * s + pow(s, 2) +
               pow(a1_mass, 2) *
                   (-2. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s)) -
          0.5 * pow(-2. + delta, 2) * pow(pion_mass, 2) *
              log(abs(-1. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + s +
                      tmin)) +
          (0.5 * (-2. + delta) * (eta1 - 1. * eta2) *
           (eta2 * pow(pion_mass, 6) * (1. * pow(m_rho, 2) - 1. * s) +
            eta2 * pow(a1_mass, 2) * pow(pion_mass, 4) *
                (-1. * pow(m_rho, 2) + 1. * s) +
            eta1 * pow(a1_mass, 2) * pow(pion_mass, 2) *
                (-0.5 * pow(m_rho, 4) +
                 pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) +
                 0.5 * pow(m_rho, 2) * s) +
            eta1 * pow(pion_mass, 4) *
                (0.5 * pow(m_rho, 4) - 0.5 * pow(m_rho, 2) * s +
                 pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s))) *
           log(abs(-1. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + s + tmin))) /
              (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
               2. * pow(a1_mass, 2) * pow(pion_mass, 2) + pow(pion_mass, 4)) +
          (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(pion_mass, 2) *
           (eta1 * (pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) +
                    (-0.5 * pow(m_rho, 2) + 0.5 * s) * s) +
            eta2 * (-0.5 * pow(m_rho, 4) + 1. * pow(m_rho, 2) * s -
                    0.5 * pow(s, 2) +
                    pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s))) *
           log(abs(-1. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + s + tmin))) /
              (1. * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
               1. * pow(m_rho, 2) + 1. * s) -
          (0.25 *
           (0. +
            8.000000000000002 * pow(2. - 1. * delta, 2) * pow(pion_mass, 4) *
                pow(m_rho, 2) +
            1. * pow(2. - 1. * delta, 2) * pow(m_rho, 4) * s +
            pow(pion_mass, 2) *
                (C4 * (32. - 16. * delta) * pow(m_rho, 6) +
                 delta * (-8. + 4. * delta) * pow(s, 2) +
                 pow(m_rho, 2) * s *
                     (-8. + 24. * delta - 10. * pow(delta, 2) + 32. * C4 * s -
                      16. * C4 * delta * s) +
                 pow(m_rho, 4) * (-16. + 8. * delta - 64. * C4 * s +
                                  32. * C4 * delta * s))) *
           log(abs(-1. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + s + tmin))) /
              (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
          0.03125 * pow(eta1 - 1. * eta2, 2) *
              (eta1 * eta2 *
                   (4. * pow(a1_mass, 6) +
                    pow(Gammaa1, 2) * pow(a1_mass, 2) *
                        (-4. * pow(a1_mass, 2) + 4. * pow(pion_mass, 2) -
                         2. * s) +
                    pow(a1_mass, 4) * (-12. * pow(pion_mass, 2) + 6. * s) +
                    pow(pion_mass, 2) *
                        (-4. * pow(pion_mass, 4) + 4. * pow(m_rho, 4) +
                         2. * pow(pion_mass, 2) * s - 2. * pow(m_rho, 2) * s) +
                    pow(a1_mass, 2) *
                        (12. * pow(pion_mass, 4) - 2. * pow(m_rho, 4) -
                         8. * pow(pion_mass, 2) * s - 4. * pow(m_rho, 2) * s +
                         4. * pow(s, 2))) +
               pow(eta1, 2) *
                   (-2. * pow(a1_mass, 6) + 2. * pow(pion_mass, 6) +
                    3. * pow(pion_mass, 4) * pow(m_rho, 2) +
                    pow(a1_mass, 4) *
                        (6. * pow(pion_mass, 2) + 3. * pow(m_rho, 2) - 3. * s) -
                    1. * pow(pion_mass, 4) * s -
                    1. * pow(pion_mass, 2) * pow(m_rho, 2) * s -
                    1. * pow(m_rho, 4) * s + pow(m_rho, 2) * pow(s, 2) +
                    pow(Gammaa1, 2) * pow(a1_mass, 2) *
                        (2. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                         1. * pow(m_rho, 2) + s) +
                    pow(a1_mass, 2) *
                        (-6. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                         4. * pow(m_rho, 2) * s - 2. * pow(s, 2) +
                         pow(pion_mass, 2) * (-6. * pow(m_rho, 2) + 4. * s))) +
               pow(eta2, 2) *
                   (-2. * pow(a1_mass, 6) + 2. * pow(pion_mass, 6) -
                    3. * pow(pion_mass, 4) * pow(m_rho, 2) +
                    pow(a1_mass, 4) *
                        (6. * pow(pion_mass, 2) - 3. * pow(m_rho, 2) - 3. * s) -
                    1. * pow(pion_mass, 4) * s +
                    pow(pion_mass, 2) * pow(m_rho, 2) * s +
                    pow(Gammaa1, 2) * pow(a1_mass, 2) *
                        (2. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) +
                         pow(m_rho, 2) + s) +
                    pow(a1_mass, 2) *
                        (-6. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                         2. * pow(m_rho, 2) * s - 2. * pow(s, 2) +
                         pow(pion_mass, 2) * (6. * pow(m_rho, 2) + 4. * s)))) *
              log(abs(pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
                      4. * pow(a1_mass, 2) * pow(pion_mass, 2) +
                      4. * pow(pion_mass, 4) +
                      2. * pow(a1_mass, 2) * (-1. * pow(m_rho, 2) + s + tmin) -
                      4. * pow(pion_mass, 2) *
                          (-1. * pow(m_rho, 2) + s + tmin) +
                      pow(-1. * pow(m_rho, 2) + s + tmin, 2))) -
          (0.5 * (1. * eta1 - 1. * eta2) *
           (eta2 * (pow(Gammaa1, 2) * pow(a1_mass, 2) * pow(m_rho, 2) *
                        (0.5 - 1. * C4 * pow(m_rho, 2)) +
                    pow(a1_mass, 4) *
                        (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) +
                    pow(pion_mass, 2) * pow(m_rho, 2) *
                        (pow(pion_mass, 2) * (-0.5 + 1. * C4 * pow(m_rho, 2)) +
                         (0.25 - 0.125 * delta) * (pow(m_rho, 2) + s)) +
                    pow(a1_mass, 2) *
                        (1. * C4 * pow(m_rho, 6) +
                         pow(pion_mass, 2) *
                             (1. * pow(m_rho, 2) - 2. * C4 * pow(m_rho, 4)) -
                         0.25 * delta * pow(s, 2) +
                         pow(m_rho, 4) * (-0.75 + 0.125 * delta - 2. * C4 * s) +
                         pow(m_rho, 2) * s *
                             (0.25 + 0.375 * delta + 1. * C4 * s))) +
            eta1 * (pow(Gammaa1, 2) * pow(a1_mass, 2) * pow(m_rho, 2) *
                        (-0.5 + 1. * C4 * pow(m_rho, 2)) +
                    pow(a1_mass, 4) *
                        (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                    pow(a1_mass, 2) *
                        (-0.5 * pow(m_rho, 4) + 1. * C4 * pow(m_rho, 6) +
                         pow(pion_mass, 2) *
                             (-1. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4)) +
                         0.25 * delta * pow(s, 2) -
                         1. * C4 * pow(m_rho, 2) * pow(s, 2)) +
                    pow(m_rho, 2) *
                        (pow(pion_mass, 4) * (0.5 - 1. * C4 * pow(m_rho, 2)) +
                         s * ((-0.25 + 0.125 * delta) * pow(m_rho, 2) +
                              (0.25 - 0.125 * delta) * s) +
                         pow(pion_mass, 2) *
                             (-2. * C4 * pow(m_rho, 4) +
                              (-0.5 - 0.25 * delta) * s +
                              pow(m_rho, 2) * (1. + 2. * C4 * s))))) *
           log(abs(pow(Gammaa1, 2) * pow(a1_mass, 2) +
                   pow(pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                           1. * pow(m_rho, 2) + s + tmin,
                       2)))) /
              pow(m_rho, 2) +
          (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
           (eta2 *
                (0.5 * pow(Gammaa1, 4) * pow(a1_mass, 4) -
                 0.5 * pow(a1_mass, 8) + 0.5 * pow(pion_mass, 8) +
                 0.5 * pow(a1_mass, 4) * pow(pion_mass, 2) * pow(m_rho, 2) -
                 0.5 * pow(pion_mass, 6) * pow(m_rho, 2) +
                 pow(Gammaa1, 2) *
                     (pow(a1_mass, 2) * pow(pion_mass, 2) *
                          (1. * pow(pion_mass, 2) + 1.5 * pow(m_rho, 2) -
                           2. * s) +
                      pow(a1_mass, 4) * (-1. * pow(pion_mass, 2) +
                                         0.5 * pow(m_rho, 2) - 1. * s)) +
                 pow(a1_mass, 6) *
                     (1. * pow(pion_mass, 2) + 0.5 * pow(m_rho, 2) - 1. * s) +
                 pow(a1_mass, 2) * pow(pion_mass, 4) *
                     (-1. * pow(pion_mass, 2) - 0.5 * pow(m_rho, 2) + 1. * s)) +
            eta1 *
                (pow(a1_mass, 6) * (1. * pow(pion_mass, 2) + 0.5 * s) +
                 pow(a1_mass, 2) * (3. * pow(pion_mass, 6) +
                                    1. * pow(pion_mass, 2) * pow(m_rho, 4) -
                                    0.5 * pow(pion_mass, 4) * s) +
                 pow(a1_mass, 4) *
                     (-3. * pow(pion_mass, 4) +
                      pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 0.5 * s) -
                      0.5 * pow(m_rho, 2) * s) +
                 pow(pion_mass, 4) *
                     (-1. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                      pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 0.5 * s) +
                      0.5 * pow(m_rho, 2) * s) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (-1. * pow(pion_mass, 4) +
                      pow(a1_mass, 2) * (1. * pow(pion_mass, 2) + 0.5 * s) -
                      0.5 * pow(m_rho, 2) * s +
                      pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1.5 * s)))) *
           log(abs(pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
                   4. * pow(pion_mass, 4) +
                   4. * pow(pion_mass, 2) * pow(m_rho, 2) + pow(m_rho, 4) -
                   4. * pow(pion_mass, 2) * s - 2. * pow(m_rho, 2) * s +
                   pow(s, 2) - 4. * pow(pion_mass, 2) * tmin -
                   2. * pow(m_rho, 2) * tmin + 2. * s * tmin + pow(tmin, 2) +
                   pow(a1_mass, 2) *
                       (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s +
                        2. * tmin)))) /
              (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
               2. * pow(a1_mass, 2) * pow(pion_mass, 2) + pow(pion_mass, 4)) -
          (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
           (eta1 * (-1. * pow(pion_mass, 8) +
                    1. * pow(pion_mass, 4) * pow(m_rho, 4) +
                    pow(a1_mass, 6) * (1. * pow(pion_mass, 2) - 0.5 * s) -
                    0.5 * pow(pion_mass, 6) * s -
                    4. * pow(pion_mass, 4) * pow(m_rho, 2) * s -
                    1.5 * pow(pion_mass, 2) * pow(m_rho, 4) * s +
                    3.5 * pow(pion_mass, 4) * pow(s, 2) +
                    4. * pow(pion_mass, 2) * pow(m_rho, 2) * pow(s, 2) +
                    0.5 * pow(m_rho, 4) * pow(s, 2) -
                    2.5 * pow(pion_mass, 2) * pow(s, 3) -
                    1. * pow(m_rho, 2) * pow(s, 3) + 0.5 * pow(s, 4) +
                    pow(Gammaa1, 2) * pow(a1_mass, 2) *
                        (-1. * pow(pion_mass, 4) +
                         pow(a1_mass, 2) * (1. * pow(pion_mass, 2) - 0.5 * s) -
                         0.5 * pow(pion_mass, 2) * s + 0.5 * pow(s, 2)) +
                    pow(a1_mass, 4) *
                        (-3. * pow(pion_mass, 4) +
                         (1. * pow(m_rho, 2) - 0.5 * s) * s +
                         pow(pion_mass, 2) * (-2. * pow(m_rho, 2) + 2.5 * s)) +
                    pow(a1_mass, 2) *
                        (3. * pow(pion_mass, 6) +
                         pow(pion_mass, 4) * (2. * pow(m_rho, 2) - 1.5 * s) -
                         0.5 * pow(m_rho, 4) * s + 0.5 * pow(s, 3) +
                         pow(pion_mass, 2) *
                             (1. * pow(m_rho, 4) - 1. * pow(m_rho, 2) * s -
                              1. * pow(s, 2)))) +
            eta2 *
                (0.5 * pow(Gammaa1, 4) * pow(a1_mass, 4) -
                 0.5 * pow(a1_mass, 8) + 0.5 * pow(pion_mass, 8) +
                 pow(a1_mass, 6) *
                     (1. * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 0.5 * s) +
                 0.5 * pow(pion_mass, 6) * s +
                 pow(a1_mass, 4) * (-0.5 * pow(m_rho, 4) +
                                    (-0.5 * pow(pion_mass, 2) + 0.5 * s) * s) +
                 pow(pion_mass, 4) *
                     (-0.5 * pow(m_rho, 4) + 2. * pow(m_rho, 2) * s -
                      1.5 * pow(s, 2)) +
                 pow(pion_mass, 2) * s *
                     (0.5 * pow(m_rho, 4) - 1. * pow(m_rho, 2) * s +
                      0.5 * pow(s, 2)) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (1. * pow(pion_mass, 4) + 0.5 * pow(m_rho, 4) +
                      pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 1.5 * s) -
                      1. * pow(m_rho, 2) * s + 0.5 * pow(s, 2) +
                      pow(a1_mass, 2) * (-1. * pow(pion_mass, 2) -
                                         1. * pow(m_rho, 2) + 1.5 * s)) +
                 pow(a1_mass, 2) *
                     (-1. * pow(pion_mass, 6) +
                      pow(pion_mass, 4) * (-1. * pow(m_rho, 2) + 0.5 * s) +
                      pow(pion_mass, 2) *
                          (-1. * pow(m_rho, 4) + 2. * pow(m_rho, 2) * s -
                           1. * pow(s, 2)) +
                      s * (0.5 * pow(m_rho, 4) - 1. * pow(m_rho, 2) * s +
                           0.5 * pow(s, 2))))) *
           log(abs(pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
                   4. * pow(pion_mass, 4) +
                   4. * pow(pion_mass, 2) * pow(m_rho, 2) + pow(m_rho, 4) -
                   4. * pow(pion_mass, 2) * s - 2. * pow(m_rho, 2) * s +
                   pow(s, 2) - 4. * pow(pion_mass, 2) * tmin -
                   2. * pow(m_rho, 2) * tmin + 2. * s * tmin + pow(tmin, 2) +
                   pow(a1_mass, 2) *
                       (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s +
                        2. * tmin)))) /
              (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
               pow(pion_mass, 4) + 2. * pow(pion_mass, 2) * pow(m_rho, 2) +
               pow(m_rho, 4) - 2. * pow(pion_mass, 2) * s -
               2. * pow(m_rho, 2) * s + pow(s, 2) +
               pow(a1_mass, 2) *
                   (-2. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s)) -
          (0.0625 * pow(eta1 - 1. * eta2, 2) *
           (pow(eta2, 2) *
                (-1. * pow(a1_mass, 10) +
                 pow(a1_mass, 8) *
                     (5. * pow(pion_mass, 2) + 2.5 * pow(m_rho, 2) - 2.5 * s) +
                 pow(Gammaa1, 4) * pow(a1_mass, 4) *
                     (1. * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
                      0.5 * pow(m_rho, 2) + 0.5 * s) +
                 pow(a1_mass, 4) *
                     (10. * pow(pion_mass, 6) - 2.5 * pow(m_rho, 6) +
                      pow(pion_mass, 4) * (15. * pow(m_rho, 2) - 9. * s) +
                      6. * pow(m_rho, 4) * s - 4.5 * pow(m_rho, 2) * pow(s, 2) +
                      1. * pow(s, 3)) +
                 pow(a1_mass, 6) *
                     (-10. * pow(pion_mass, 4) +
                      (1. * pow(m_rho, 2) - 1. * s) * s +
                      pow(pion_mass, 2) * (-10. * pow(m_rho, 2) + 8. * s)) +
                 pow(pion_mass, 4) *
                     (1. * pow(pion_mass, 6) + 0.5 * pow(m_rho, 6) +
                      pow(pion_mass, 4) * (2.5 * pow(m_rho, 2) - 0.5 * s) -
                      1. * pow(m_rho, 4) * s + 0.5 * pow(m_rho, 2) * pow(s, 2) +
                      pow(pion_mass, 2) *
                          (2. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s)) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (4. * pow(a1_mass, 6) - 4. * pow(pion_mass, 6) -
                      0.5 * pow(m_rho, 6) + 1.5 * pow(m_rho, 4) * s -
                      1.5 * pow(m_rho, 2) * pow(s, 2) + 0.5 * pow(s, 3) +
                      pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 6. * s) +
                      pow(a1_mass, 4) * (-12. * pow(pion_mass, 2) -
                                         6. * pow(m_rho, 2) + 6. * s) +
                      pow(pion_mass, 2) *
                          (-3. * pow(m_rho, 4) + 6. * pow(m_rho, 2) * s -
                           3. * pow(s, 2)) +
                      pow(a1_mass, 2) *
                          (12. * pow(pion_mass, 4) + 3. * pow(m_rho, 4) +
                           pow(pion_mass, 2) * (12. * pow(m_rho, 2) - 12. * s) -
                           6. * pow(m_rho, 2) * s + 3. * pow(s, 2))) +
                 pow(a1_mass, 2) *
                     (-5. * pow(pion_mass, 8) + 1. * pow(m_rho, 8) -
                      3.5 * pow(m_rho, 6) * s +
                      4.5 * pow(m_rho, 4) * pow(s, 2) -
                      2.5 * pow(m_rho, 2) * pow(s, 3) + 0.5 * pow(s, 4) +
                      pow(pion_mass, 6) * (-10. * pow(m_rho, 2) + 4. * s) +
                      pow(pion_mass, 4) *
                          (-2. * pow(m_rho, 4) + 1. * pow(m_rho, 2) * s +
                           1. * pow(s, 2)) +
                      pow(pion_mass, 2) *
                          (3. * pow(m_rho, 6) - 8. * pow(m_rho, 4) * s +
                           7. * pow(m_rho, 2) * pow(s, 2) - 2. * pow(s, 3)))) +
            pow(eta1, 2) *
                (-1. * pow(a1_mass, 10) +
                 pow(a1_mass, 8) *
                     (5. * pow(pion_mass, 2) + 2.5 * pow(m_rho, 2) - 2.5 * s) +
                 pow(Gammaa1, 4) * pow(a1_mass, 4) *
                     (1. * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
                      0.5 * pow(m_rho, 2) + 0.5 * s) +
                 pow(a1_mass, 6) *
                     (-10. * pow(pion_mass, 4) - 2. * pow(m_rho, 4) +
                      5. * pow(m_rho, 2) * s - 1. * pow(s, 2) +
                      pow(pion_mass, 2) * (-10. * pow(m_rho, 2) + 8. * s)) +
                 pow(a1_mass, 4) *
                     (10. * pow(pion_mass, 6) + 0.5 * pow(m_rho, 6) +
                      pow(pion_mass, 4) * (15. * pow(m_rho, 2) - 9. * s) -
                      3. * pow(m_rho, 4) * s + 1.5 * pow(m_rho, 2) * pow(s, 2) +
                      1. * pow(s, 3) +
                      pow(pion_mass, 2) *
                          (6. * pow(m_rho, 4) - 12. * pow(m_rho, 2) * s)) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (4. * pow(a1_mass, 6) - 4. * pow(pion_mass, 6) -
                      0.5 * pow(m_rho, 6) + 1.5 * pow(m_rho, 4) * s -
                      1.5 * pow(m_rho, 2) * pow(s, 2) + 0.5 * pow(s, 3) +
                      pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 6. * s) +
                      pow(a1_mass, 4) * (-12. * pow(pion_mass, 2) -
                                         6. * pow(m_rho, 2) + 6. * s) +
                      pow(pion_mass, 2) *
                          (-3. * pow(m_rho, 4) + 6. * pow(m_rho, 2) * s -
                           3. * pow(s, 2)) +
                      pow(a1_mass, 2) *
                          (12. * pow(pion_mass, 4) + 3. * pow(m_rho, 4) +
                           pow(pion_mass, 2) * (12. * pow(m_rho, 2) - 12. * s) -
                           6. * pow(m_rho, 2) * s + 3. * pow(s, 2))) +
                 pow(pion_mass, 2) *
                     (1. * pow(pion_mass, 8) +
                      pow(pion_mass, 6) * (2.5 * pow(m_rho, 2) - 0.5 * s) +
                      pow(pion_mass, 4) *
                          (4. * pow(m_rho, 4) - 6. * pow(m_rho, 2) * s) +
                      pow(m_rho, 2) * s *
                          (-1. * pow(m_rho, 4) + 2. * pow(m_rho, 2) * s -
                           1. * pow(s, 2)) +
                      pow(pion_mass, 2) *
                          (1.5 * pow(m_rho, 6) - 6. * pow(m_rho, 4) * s +
                           4.5 * pow(m_rho, 2) * pow(s, 2))) +
                 pow(a1_mass, 2) *
                     (-5. * pow(pion_mass, 8) +
                      pow(pion_mass, 6) * (-10. * pow(m_rho, 2) + 4. * s) +
                      pow(pion_mass, 4) *
                          (-8. * pow(m_rho, 4) + 13. * pow(m_rho, 2) * s +
                           1. * pow(s, 2)) +
                      pow(pion_mass, 2) *
                          (-1. * pow(m_rho, 6) + 6. * pow(m_rho, 4) * s -
                           3. * pow(m_rho, 2) * pow(s, 2) - 2. * pow(s, 3)) +
                      s * (0.5 * pow(m_rho, 6) - 0.5 * pow(m_rho, 4) * s -
                           0.5 * pow(m_rho, 2) * pow(s, 2) +
                           0.5 * pow(s, 3)))) +
            eta1 * eta2 *
                (2. * pow(a1_mass, 10) +
                 pow(Gammaa1, 4) * pow(a1_mass, 4) *
                     (-2. * pow(a1_mass, 2) + 2. * pow(pion_mass, 2) +
                      1. * pow(m_rho, 2) - 1. * s) +
                 pow(a1_mass, 8) *
                     (-10. * pow(pion_mass, 2) - 5. * pow(m_rho, 2) + 5. * s) +
                 pow(a1_mass, 6) *
                     (20. * pow(pion_mass, 4) + 6. * pow(m_rho, 4) +
                      pow(pion_mass, 2) * (20. * pow(m_rho, 2) - 16. * s) -
                      8. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                 pow(a1_mass, 4) *
                     (-20. * pow(pion_mass, 6) - 4. * pow(m_rho, 6) +
                      6. * pow(m_rho, 4) * s - 2. * pow(s, 3) +
                      pow(pion_mass, 4) * (-30. * pow(m_rho, 2) + 18. * s) +
                      pow(pion_mass, 2) *
                          (-18. * pow(m_rho, 4) + 18. * pow(m_rho, 2) * s)) +
                 pow(pion_mass, 2) *
                     (-2. * pow(pion_mass, 8) - 1. * pow(m_rho, 8) +
                      3. * pow(m_rho, 6) * s - 3. * pow(m_rho, 4) * pow(s, 2) +
                      1. * pow(m_rho, 2) * pow(s, 3) +
                      pow(pion_mass, 6) * (-5. * pow(m_rho, 2) + 1. * s) +
                      pow(pion_mass, 4) *
                          (-10. * pow(m_rho, 4) + 10. * pow(m_rho, 2) * s) +
                      pow(pion_mass, 2) *
                          (-6. * pow(m_rho, 6) + 12. * pow(m_rho, 4) * s -
                           6. * pow(m_rho, 2) * pow(s, 2))) +
                 pow(a1_mass, 2) *
                     (10. * pow(pion_mass, 8) + 1. * pow(m_rho, 8) +
                      pow(pion_mass, 6) * (20. * pow(m_rho, 2) - 8. * s) -
                      2. * pow(m_rho, 6) * s + 2. * pow(m_rho, 2) * pow(s, 3) -
                      1. * pow(s, 4) +
                      pow(pion_mass, 4) *
                          (22. * pow(m_rho, 4) - 20. * pow(m_rho, 2) * s -
                           2. * pow(s, 2)) +
                      pow(pion_mass, 2) *
                          (8. * pow(m_rho, 6) - 12. * pow(m_rho, 4) * s +
                           4. * pow(s, 3))) +
                 pow(Gammaa1, 2) * pow(a1_mass, 2) *
                     (-8. * pow(a1_mass, 6) + 8. * pow(pion_mass, 6) +
                      1. * pow(m_rho, 6) +
                      pow(pion_mass, 4) * (12. * pow(m_rho, 2) - 12. * s) +
                      pow(a1_mass, 4) * (24. * pow(pion_mass, 2) +
                                         12. * pow(m_rho, 2) - 12. * s) -
                      3. * pow(m_rho, 4) * s + 3. * pow(m_rho, 2) * pow(s, 2) -
                      1. * pow(s, 3) +
                      pow(pion_mass, 2) *
                          (6. * pow(m_rho, 4) - 12. * pow(m_rho, 2) * s +
                           6. * pow(s, 2)) +
                      pow(a1_mass, 2) *
                          (-24. * pow(pion_mass, 4) - 6. * pow(m_rho, 4) +
                           12. * pow(m_rho, 2) * s - 6. * pow(s, 2) +
                           pow(pion_mass, 2) *
                               (-24. * pow(m_rho, 2) + 24. * s))))) *
           log(abs(pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
                   4. * pow(pion_mass, 4) +
                   4. * pow(pion_mass, 2) * pow(m_rho, 2) + pow(m_rho, 4) -
                   4. * pow(pion_mass, 2) * s - 2. * pow(m_rho, 2) * s +
                   pow(s, 2) - 4. * pow(pion_mass, 2) * tmin -
                   2. * pow(m_rho, 2) * tmin + 2. * s * tmin + pow(tmin, 2) +
                   pow(a1_mass, 2) *
                       (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s +
                        2. * tmin)))) /
              (pow(Gammaa1, 2) * pow(a1_mass, 2) + 4. * pow(a1_mass, 4) +
               4. * pow(pion_mass, 4) + 4. * pow(pion_mass, 2) * pow(m_rho, 2) +
               pow(m_rho, 4) - 4. * pow(pion_mass, 2) * s -
               2. * pow(m_rho, 2) * s + pow(s, 2) +
               pow(a1_mass, 2) *
                   (-8. * pow(pion_mass, 2) - 4. * pow(m_rho, 2) + 4. * s)))) /
           (16. * Pi * s * (-4 * pow(pion_mass, 2) + s)) +
       (pow(Const, 2) * pow(ghat, 4) *
        (0. +
         (0.03125 * pow(eta1 - 1. * eta2, 2) *
          (eta1 * eta2 *
               (-2. * pow(a1_mass, 8) - 2. * pow(pion_mass, 8) +
                2. * pow(pion_mass, 4) * pow(m_rho, 4) +
                pow(a1_mass, 6) * (8. * pow(pion_mass, 2) - 4. * s) +
                pow(a1_mass, 2) * pow(pion_mass, 2) *
                    (8. * pow(pion_mass, 4) - 8. * pow(m_rho, 4) -
                     4. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s) +
                pow(a1_mass, 4) *
                    (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                     8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                     4. * pow(s, 2))) +
           pow(eta2, 2) *
               (1. * pow(a1_mass, 8) + 1. * pow(pion_mass, 8) -
                2. * pow(pion_mass, 6) * pow(m_rho, 2) +
                1. * pow(pion_mass, 4) * pow(m_rho, 4) +
                pow(a1_mass, 6) *
                    (-4. * pow(pion_mass, 2) + 2. * pow(m_rho, 2) + 2. * s) +
                pow(a1_mass, 4) *
                    (6. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                     pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                     2. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                pow(a1_mass, 2) *
                    (-4. * pow(pion_mass, 6) -
                     2. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                     pow(pion_mass, 4) * (6. * pow(m_rho, 2) + 2. * s))) +
           pow(eta1, 2) *
               (1. * pow(a1_mass, 8) +
                pow(a1_mass, 6) *
                    (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
                pow(a1_mass, 4) *
                    (6. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                     pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                     4. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                pow(a1_mass, 2) *
                    (-4. * pow(pion_mass, 6) +
                     2. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                     pow(m_rho, 2) * (2. * pow(m_rho, 2) - 2. * s) * s +
                     pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 2. * s)) +
                pow(pion_mass, 2) *
                    (1. * pow(pion_mass, 6) +
                     2. * pow(pion_mass, 4) * pow(m_rho, 2) -
                     2. * pow(m_rho, 6) + 2. * pow(m_rho, 4) * s +
                     pow(pion_mass, 2) *
                         (1. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s))))) /
             (1. * pow(a1_mass, 2) - 1. * tmax) +
         (1. * pow(-2. + delta, 2) * pow(pion_mass, 2) *
          (1. * pow(pion_mass, 2) - 0.25 * pow(m_rho, 2))) /
             (1. * pow(pion_mass, 2) - 1. * tmax) +
         (1. * pow(-2. + delta, 2) * pow(pion_mass, 2) *
          (1. * pow(pion_mass, 2) - 0.25 * pow(m_rho, 2))) /
             (1. * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 1. * s -
              1. * tmax) -
         (0.5 * pow(-2. + delta, 2) * pow(pion_mass, 2) * tmax) /
             pow(m_rho, 2) -
         0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
             (-0.5 * eta2 * pow(a1_mass, 2) + 1. * eta1 * pow(pion_mass, 2) +
              0.5 * eta2 * pow(m_rho, 2) + 0.5 * eta1 * s - 1. * eta2 * s) *
             tmax +
         (0.25 *
          (pow(pion_mass, 2) *
               (12. + 1. * pow(delta, 2) - 16. * C4 * pow(m_rho, 2) +
                delta * (-8. + 8. * C4 * pow(m_rho, 2))) +
           (-4. - 3. * pow(delta, 2) - 16. * C4 * pow(m_rho, 2) +
            delta * (8. + 8. * C4 * pow(m_rho, 2))) *
               s) *
          tmax) /
             pow(m_rho, 2) -
         0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta2 * (pow(a1_mass, 2) - 1. * s) +
              eta1 * (-2. * pow(pion_mass, 2) + s)) *
             tmax -
         0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta2 * (-1. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * s) +
              eta1 * (2. * pow(pion_mass, 2) + s)) *
             tmax +
         0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta1 * (1. * pow(pion_mass, 2) - 0.5 * s) +
              eta2 * (-0.5 * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
                      0.5 * pow(m_rho, 2) + 1. * s)) *
             tmax +
         (0.25 * (-2. + 1. * delta) *
          (-8. * C4 * pow(m_rho, 4) +
           pow(pion_mass, 2) * (2. + 1. * delta - 8. * C4 * pow(m_rho, 2)) +
           (-2. - 3. * delta) * s +
           pow(m_rho, 2) * (2. + 1. * delta + 16. * C4 * s)) *
          tmax) /
             pow(m_rho, 2) +
         (0.25 *
          (32 * pow(C4, 2) * pow(m_rho, 8) + 2 * pow(delta, 2) * pow(s, 2) +
           8 * C4 * pow(m_rho, 6) * (-6 + delta - 8 * C4 * s) +
           2 * delta * pow(m_rho, 2) * s * (-6 + delta - 8 * C4 * s) +
           pow(m_rho, 4) * (12 - pow(delta, 2) + 8 * C4 * (6 + delta) * s +
                            32 * pow(C4, 2) * pow(s, 2))) *
          tmax) /
             pow(m_rho, 4) -
         (1. * (1. * eta1 - 1. * eta2) *
          (eta2 *
               (0.75 * pow(m_rho, 4) - 0.125 * delta * pow(m_rho, 4) -
                1. * C4 * pow(m_rho, 6) +
                pow(a1_mass, 2) *
                    (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                pow(pion_mass, 2) *
                    (-1. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4)) -
                0.25 * pow(m_rho, 2) * s - 0.375 * delta * pow(m_rho, 2) * s +
                2. * C4 * pow(m_rho, 4) * s + 0.25 * delta * pow(s, 2) -
                1. * C4 * pow(m_rho, 2) * pow(s, 2)) +
           eta1 * (0.5 * pow(m_rho, 4) - 1. * C4 * pow(m_rho, 6) +
                   pow(pion_mass, 2) *
                       (1. * pow(m_rho, 2) - 2. * C4 * pow(m_rho, 4)) +
                   pow(a1_mass, 2) *
                       (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) -
                   0.25 * delta * pow(s, 2) +
                   1. * C4 * pow(m_rho, 2) * pow(s, 2))) *
          tmax) /
             pow(m_rho, 2) +
         0.0625 * pow(eta1 - 1. * eta2, 2) *
             (pow(eta2, 2) *
                  (pow(Gammaa1, 2) * pow(a1_mass, 2) - 1. * pow(a1_mass, 4) -
                   2. * pow(pion_mass, 4) -
                   2. * pow(pion_mass, 2) * pow(m_rho, 2) + 2. * pow(m_rho, 4) +
                   pow(a1_mass, 2) *
                       (2. * pow(pion_mass, 2) + pow(m_rho, 2) - 1. * s) -
                   3. * pow(m_rho, 2) * s + pow(s, 2)) +
              pow(eta1, 2) * (pow(Gammaa1, 2) * pow(a1_mass, 2) -
                              1. * pow(a1_mass, 4) - 2. * pow(pion_mass, 4) -
                              2. * pow(pion_mass, 2) * pow(m_rho, 2) +
                              pow(a1_mass, 2) * (2. * pow(pion_mass, 2) +
                                                 pow(m_rho, 2) - 1. * s) +
                              pow(m_rho, 2) * s + pow(s, 2)) +
              eta1 * eta2 *
                  (-2. * pow(Gammaa1, 2) * pow(a1_mass, 2) +
                   2. * pow(a1_mass, 4) + 4. * pow(pion_mass, 4) +
                   4. * pow(pion_mass, 2) * pow(m_rho, 2) + 2. * pow(m_rho, 4) -
                   2. * pow(s, 2) +
                   pow(a1_mass, 2) * (-4. * pow(pion_mass, 2) -
                                      2. * pow(m_rho, 2) + 2. * s))) *
             tmax +
         0.03125 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-6. * pow(a1_mass, 4) - 12. * pow(pion_mass, 4) +
                   2. * pow(m_rho, 4) +
                   pow(a1_mass, 2) * (16. * pow(pion_mass, 2) - 8. * s) +
                   8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                   4. * pow(s, 2)) +
              pow(eta1, 2) *
                  (3. * pow(a1_mass, 4) + 6. * pow(pion_mass, 4) +
                   pow(m_rho, 4) +
                   pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                   4. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                   pow(a1_mass, 2) * (-8. * pow(pion_mass, 2) -
                                      4. * pow(m_rho, 2) + 4. * s)) +
              pow(eta2, 2) *
                  (3. * pow(a1_mass, 4) + 6. * pow(pion_mass, 4) +
                   pow(m_rho, 4) +
                   pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                   2. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                   pow(a1_mass, 2) * (-8. * pow(pion_mass, 2) +
                                      4. * pow(m_rho, 2) + 4. * s))) *
             tmax -
         (0.125 * (-2. + 1. * delta) *
          (2. + 1. * delta - 8. * C4 * pow(m_rho, 2)) * pow(tmax, 2)) /
             pow(m_rho, 2) -
         0.5 * pow(1. * eta1 - 1. * eta2, 2) *
             (-0.5 + 1. * C4 * pow(m_rho, 2)) * pow(tmax, 2) -
         (1. *
          (0.5 - 0.125 * pow(delta, 2) - 2. * C4 * pow(m_rho, 2) +
           1. * C4 * delta * pow(m_rho, 2)) *
          pow(tmax, 2)) /
             pow(m_rho, 2) +
         0.0625 * pow(1. * eta1 - 1. * eta2, 4) *
             (1. * pow(pion_mass, 2) + 0.5 * pow(m_rho, 2) - 0.5 * s) *
             pow(tmax, 2) +
         0.03125 * pow(eta1 - 1. * eta2, 3) *
             (eta2 * (-1. * pow(a1_mass, 2) + 2. * pow(pion_mass, 2) -
                      1. * pow(m_rho, 2) - 1. * s) +
              eta1 * (pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                      1. * pow(m_rho, 2) + s)) *
             pow(tmax, 2) +
         0.010416666666666666 * pow(eta1 - 1. * eta2, 4) * pow(tmax, 3) -
         0.020833333333333332 * pow(1. * eta1 - 1. * eta2, 4) * pow(tmax, 3) +
         0.03125 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (2. * pow(Gammaa1, 2) * pow(a1_mass, 2) -
                   6. * pow(a1_mass, 4) - 4. * pow(pion_mass, 4) +
                   2. * pow(m_rho, 4) +
                   pow(a1_mass, 2) * (8. * pow(pion_mass, 2) - 8. * s) +
                   4. * pow(m_rho, 2) * s - 4. * pow(s, 2)) +
              pow(eta1, 2) *
                  (-1. * pow(Gammaa1, 2) * pow(a1_mass, 2) +
                   3. * pow(a1_mass, 4) + 2. * pow(pion_mass, 4) +
                   2. * pow(pion_mass, 2) * pow(m_rho, 2) + pow(m_rho, 4) -
                   4. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                   pow(a1_mass, 2) * (-4. * pow(pion_mass, 2) -
                                      4. * pow(m_rho, 2) + 4. * s)) +
              pow(eta2, 2) *
                  (-1. * pow(Gammaa1, 2) * pow(a1_mass, 2) +
                   3. * pow(a1_mass, 4) + 2. * pow(pion_mass, 4) -
                   2. * pow(pion_mass, 2) * pow(m_rho, 2) + pow(m_rho, 4) -
                   2. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                   pow(a1_mass, 2) * (-4. * pow(pion_mass, 2) +
                                      4. * pow(m_rho, 2) + 4. * s))) *
             (-1. * pow(m_rho, 2) + s + tmax) -
         0.03125 * pow(eta1 - 1. * eta2, 3) *
             (eta2 * (-1. * pow(a1_mass, 2) - 1. * pow(m_rho, 2) - 1. * s) +
              eta1 * (pow(a1_mass, 2) - 1. * pow(m_rho, 2) + s)) *
             pow(-1. * pow(m_rho, 2) + s + tmax, 2) +
         0.010416666666666666 * pow(eta1 - 1. * eta2, 4) *
             pow(-1. * pow(m_rho, 2) + s + tmax, 3) +
         0.25 * (eta1 - 1. * eta2) * (1. * eta1 - 1. * eta2) *
             (-1. + 2. * C4 * pow(m_rho, 2)) *
             pow(pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
                     s + tmax,
                 2) -
         (2. * (1. * eta1 - 1. * eta2) *
          (eta2 *
               (0.375 * pow(m_rho, 4) - 0.0625 * delta * pow(m_rho, 4) -
                0.5 * C4 * pow(m_rho, 6) +
                pow(a1_mass, 2) *
                    (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                pow(pion_mass, 2) *
                    (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) -
                0.125 * pow(m_rho, 2) * s - 0.1875 * delta * pow(m_rho, 2) * s +
                1. * C4 * pow(m_rho, 4) * s + 0.125 * delta * pow(s, 2) -
                0.5 * C4 * pow(m_rho, 2) * pow(s, 2)) +
           eta1 * (0.25 * pow(m_rho, 4) - 0.5 * C4 * pow(m_rho, 6) +
                   pow(pion_mass, 2) *
                       (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                   pow(a1_mass, 2) *
                       (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) -
                   0.125 * delta * pow(s, 2) +
                   0.5 * C4 * pow(m_rho, 2) * pow(s, 2))) *
          (1. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
           1. * s + 1. * tmax)) /
             pow(m_rho, 2) +
         (2. * (1. * eta1 - 1. * eta2) * Gammaa1 * a1_mass *
          (eta2 *
               (0.375 * pow(m_rho, 4) - 0.0625 * delta * pow(m_rho, 4) -
                0.5 * C4 * pow(m_rho, 6) +
                pow(a1_mass, 2) *
                    (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                pow(pion_mass, 2) *
                    (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) -
                0.125 * pow(m_rho, 2) * s - 0.1875 * delta * pow(m_rho, 2) * s +
                1. * C4 * pow(m_rho, 4) * s + 0.125 * delta * pow(s, 2) -
                0.5 * C4 * pow(m_rho, 2) * pow(s, 2)) +
           eta1 * (0.25 * pow(m_rho, 4) - 0.5 * C4 * pow(m_rho, 6) +
                   pow(pion_mass, 2) *
                       (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                   pow(a1_mass, 2) *
                       (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) -
                   0.125 * delta * pow(s, 2) +
                   0.5 * C4 * pow(m_rho, 2) * pow(s, 2))) *
          atan((pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
                s + tmax) /
               (Gammaa1 * a1_mass))) /
             pow(m_rho, 2) +
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) * Gammaa1 * a1_mass *
          (eta2 *
               (-1. * pow(a1_mass, 6) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (-1. * pow(a1_mass, 2) + 0.5 * pow(m_rho, 2) - 1. * s) +
                pow(a1_mass, 4) *
                    (2. * pow(pion_mass, 2) + 0.5 * pow(m_rho, 2) - 1. * s) +
                pow(pion_mass, 4) * (-1.5 * pow(m_rho, 2) + 1. * s) +
                pow(a1_mass, 2) * pow(pion_mass, 2) *
                    (-1. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + 2. * s)) +
           eta1 * (pow(Gammaa1, 2) * pow(a1_mass, 2) *
                       (1. * pow(pion_mass, 2) + 0.5 * s) +
                   pow(a1_mass, 4) * (1. * pow(pion_mass, 2) + 0.5 * s) +
                   pow(a1_mass, 2) *
                       (-2. * pow(pion_mass, 4) - 1. * pow(pion_mass, 2) * s) +
                   pow(pion_mass, 2) *
                       (1. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                        pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 1.5 * s) +
                        1. * pow(m_rho, 2) * s))) *
          atan((pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
                s + tmax) /
               (Gammaa1 * a1_mass))) /
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * pow(pion_mass, 2) + pow(pion_mass, 4)) -
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) * Gammaa1 * a1_mass *
          (eta2 *
               (-1. * pow(a1_mass, 6) - 2. * pow(pion_mass, 4) * pow(m_rho, 2) -
                1. * pow(pion_mass, 2) * pow(m_rho, 4) +
                pow(a1_mass, 4) *
                    (2. * pow(pion_mass, 2) + 2. * pow(m_rho, 2) - 1.5 * s) +
                2.5 * pow(pion_mass, 4) * s +
                3. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                0.5 * pow(m_rho, 4) * s - 2. * pow(pion_mass, 2) * pow(s, 2) -
                1. * pow(m_rho, 2) * pow(s, 2) + 0.5 * pow(s, 3) +
                pow(Gammaa1, 2) *
                    (-1. * pow(a1_mass, 4) + 0.5 * pow(a1_mass, 2) * s) +
                pow(a1_mass, 2) *
                    (-1. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                     1. * pow(m_rho, 2) * s +
                     pow(pion_mass, 2) * (-2. * pow(m_rho, 2) + 1. * s))) +
           eta1 *
               (1. * pow(pion_mass, 6) +
                4. * pow(pion_mass, 4) * pow(m_rho, 2) +
                1. * pow(pion_mass, 2) * pow(m_rho, 4) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (1. * pow(pion_mass, 2) - 0.5 * s) +
                pow(a1_mass, 4) * (1. * pow(pion_mass, 2) - 0.5 * s) -
                4.5 * pow(pion_mass, 4) * s -
                4. * pow(pion_mass, 2) * pow(m_rho, 2) * s -
                0.5 * pow(m_rho, 4) * s + 3. * pow(pion_mass, 2) * pow(s, 2) +
                1. * pow(m_rho, 2) * pow(s, 2) - 0.5 * pow(s, 3) +
                pow(a1_mass, 2) *
                    (-2. * pow(pion_mass, 4) +
                     (1. * pow(m_rho, 2) - 1. * s) * s +
                     pow(pion_mass, 2) * (-2. * pow(m_rho, 2) + 3. * s)))) *
          atan((pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
                s + tmax) /
               (Gammaa1 * a1_mass))) /
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
              pow(pion_mass, 4) + 2. * pow(pion_mass, 2) * pow(m_rho, 2) +
              pow(m_rho, 4) - 2. * pow(pion_mass, 2) * s -
              2. * pow(m_rho, 2) * s + pow(s, 2) +
              pow(a1_mass, 2) *
                  (-2. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s)) +
         (0.03125 * pow(eta1 - 1. * eta2, 2) *
          (pow(eta2, 2) *
               (pow(Gammaa1, 4) * pow(a1_mass, 4) + pow(a1_mass, 8) +
                pow(pion_mass, 8) - 2. * pow(pion_mass, 6) * pow(m_rho, 2) +
                pow(pion_mass, 4) * pow(m_rho, 4) +
                pow(a1_mass, 6) *
                    (-4. * pow(pion_mass, 2) + 2. * pow(m_rho, 2) + 2. * s) +
                pow(a1_mass, 4) *
                    (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                     pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                     2. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                pow(a1_mass, 2) *
                    (-4. * pow(pion_mass, 6) -
                     2. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                     pow(pion_mass, 4) * (6. * pow(m_rho, 2) + 2. * s)) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (-6. * pow(a1_mass, 4) - 6. * pow(pion_mass, 4) -
                     1. * pow(m_rho, 4) +
                     pow(a1_mass, 2) * (12. * pow(pion_mass, 2) -
                                        6. * pow(m_rho, 2) - 6. * s) +
                     2. * pow(m_rho, 2) * s - 2. * pow(s, 2) +
                     pow(pion_mass, 2) * (6. * pow(m_rho, 2) + 4. * s))) +
           eta1 * eta2 *
               (-2. * pow(Gammaa1, 4) * pow(a1_mass, 4) - 2. * pow(a1_mass, 8) -
                2. * pow(pion_mass, 8) +
                2. * pow(pion_mass, 4) * pow(m_rho, 4) +
                pow(a1_mass, 6) * (8. * pow(pion_mass, 2) - 4. * s) +
                pow(a1_mass, 2) * pow(pion_mass, 2) *
                    (8. * pow(pion_mass, 4) - 8. * pow(m_rho, 4) -
                     4. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s) +
                pow(a1_mass, 4) *
                    (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                     8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                     4. * pow(s, 2)) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (12. * pow(a1_mass, 4) + 12. * pow(pion_mass, 4) -
                     2. * pow(m_rho, 4) - 8. * pow(pion_mass, 2) * s -
                     4. * pow(m_rho, 2) * s + 4. * pow(s, 2) +
                     pow(a1_mass, 2) * (-24. * pow(pion_mass, 2) + 12. * s))) +
           pow(eta1, 2) *
               (pow(Gammaa1, 4) * pow(a1_mass, 4) + pow(a1_mass, 8) +
                pow(a1_mass, 6) *
                    (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
                pow(a1_mass, 4) *
                    (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                     pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                     4. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                pow(a1_mass, 2) *
                    (-4. * pow(pion_mass, 6) +
                     2. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                     pow(m_rho, 2) * (2. * pow(m_rho, 2) - 2. * s) * s +
                     pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 2. * s)) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (-6. * pow(a1_mass, 4) - 6. * pow(pion_mass, 4) -
                     1. * pow(m_rho, 4) +
                     pow(a1_mass, 2) * (12. * pow(pion_mass, 2) +
                                        6. * pow(m_rho, 2) - 6. * s) +
                     4. * pow(m_rho, 2) * s - 2. * pow(s, 2) +
                     pow(pion_mass, 2) * (-6. * pow(m_rho, 2) + 4. * s)) +
                pow(pion_mass, 2) *
                    (pow(pion_mass, 6) +
                     2. * pow(pion_mass, 4) * pow(m_rho, 2) -
                     2. * pow(m_rho, 6) + 2. * pow(m_rho, 4) * s +
                     pow(pion_mass, 2) *
                         (pow(m_rho, 4) - 2. * pow(m_rho, 2) * s)))) *
          atan((pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
                s + tmax) /
               (Gammaa1 * a1_mass))) /
             (Gammaa1 * a1_mass) -
         (0.0625 * pow(eta1 - 1. * eta2, 2) * Gammaa1 * a1_mass *
          (eta1 * eta2 *
               (-2. * pow(Gammaa1, 4) * pow(a1_mass, 4) +
                14. * pow(a1_mass, 8) + 14. * pow(pion_mass, 8) +
                28. * pow(pion_mass, 6) * pow(m_rho, 2) +
                20. * pow(pion_mass, 4) * pow(m_rho, 4) +
                10. * pow(pion_mass, 2) * pow(m_rho, 6) + 2. * pow(m_rho, 8) -
                16. * pow(pion_mass, 6) * s -
                16. * pow(pion_mass, 4) * pow(m_rho, 2) * s -
                12. * pow(pion_mass, 2) * pow(m_rho, 4) * s -
                4. * pow(m_rho, 6) * s - 4. * pow(pion_mass, 4) * pow(s, 2) -
                6. * pow(pion_mass, 2) * pow(m_rho, 2) * pow(s, 2) +
                8. * pow(pion_mass, 2) * pow(s, 3) +
                4. * pow(m_rho, 2) * pow(s, 3) - 2. * pow(s, 4) +
                pow(a1_mass, 6) *
                    (-56. * pow(pion_mass, 2) - 28. * pow(m_rho, 2) + 28. * s) +
                pow(a1_mass, 4) *
                    (84. * pow(pion_mass, 4) + 24. * pow(m_rho, 4) +
                     pow(pion_mass, 2) * (84. * pow(m_rho, 2) - 72. * s) -
                     36. * pow(m_rho, 2) * s + 12. * pow(s, 2)) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (-4. * pow(a1_mass, 4) - 4. * pow(pion_mass, 4) +
                     pow(a1_mass, 2) * (8. * pow(pion_mass, 2) +
                                        4. * pow(m_rho, 2) - 4. * s) +
                     (4. * pow(m_rho, 2) - 4. * s) * s +
                     pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 8. * s)) +
                pow(a1_mass, 2) *
                    (-56. * pow(pion_mass, 6) - 10. * pow(m_rho, 6) +
                     18. * pow(m_rho, 4) * s - 6. * pow(m_rho, 2) * pow(s, 2) -
                     2. * pow(s, 3) +
                     pow(pion_mass, 4) * (-84. * pow(m_rho, 2) + 60. * s) +
                     pow(pion_mass, 2) *
                         (-48. * pow(m_rho, 4) + 60. * pow(m_rho, 2) * s -
                          12. * pow(s, 2)))) +
           pow(eta1, 2) *
               (1. * pow(Gammaa1, 4) * pow(a1_mass, 4) - 7. * pow(a1_mass, 8) -
                7. * pow(pion_mass, 8) -
                14. * pow(pion_mass, 6) * pow(m_rho, 2) -
                7. * pow(pion_mass, 4) * pow(m_rho, 4) -
                2. * pow(pion_mass, 2) * pow(m_rho, 6) +
                pow(a1_mass, 6) *
                    (28. * pow(pion_mass, 2) + 14. * pow(m_rho, 2) - 14. * s) +
                8. * pow(pion_mass, 6) * s +
                11. * pow(pion_mass, 4) * pow(m_rho, 2) * s +
                6. * pow(pion_mass, 2) * pow(m_rho, 4) * s +
                1. * pow(m_rho, 6) * s + 2. * pow(pion_mass, 4) * pow(s, 2) -
                1. * pow(m_rho, 4) * pow(s, 2) -
                4. * pow(pion_mass, 2) * pow(s, 3) -
                1. * pow(m_rho, 2) * pow(s, 3) + 1. * pow(s, 4) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (2. * pow(a1_mass, 4) + 2. * pow(pion_mass, 4) +
                     1. * pow(m_rho, 4) +
                     pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 4. * s) -
                     1. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                     pow(a1_mass, 2) * (-4. * pow(pion_mass, 2) -
                                        2. * pow(m_rho, 2) + 2. * s)) +
                pow(a1_mass, 4) *
                    (-42. * pow(pion_mass, 4) - 9. * pow(m_rho, 4) +
                     21. * pow(m_rho, 2) * s - 6. * pow(s, 2) +
                     pow(pion_mass, 2) * (-42. * pow(m_rho, 2) + 36. * s)) +
                pow(a1_mass, 2) *
                    (28. * pow(pion_mass, 6) + 2. * pow(m_rho, 6) +
                     pow(pion_mass, 4) * (42. * pow(m_rho, 2) - 30. * s) -
                     9. * pow(m_rho, 4) * s + 6. * pow(m_rho, 2) * pow(s, 2) +
                     1. * pow(s, 3) +
                     pow(pion_mass, 2) *
                         (18. * pow(m_rho, 4) - 36. * pow(m_rho, 2) * s +
                          6. * pow(s, 2)))) +
           pow(eta2, 2) *
               (1. * pow(Gammaa1, 4) * pow(a1_mass, 4) - 7. * pow(a1_mass, 8) -
                7. * pow(pion_mass, 8) -
                14. * pow(pion_mass, 6) * pow(m_rho, 2) -
                1. * pow(pion_mass, 4) * pow(m_rho, 4) +
                6. * pow(pion_mass, 2) * pow(m_rho, 6) + 2. * pow(m_rho, 8) +
                pow(a1_mass, 6) *
                    (28. * pow(pion_mass, 2) + 14. * pow(m_rho, 2) - 14. * s) +
                8. * pow(pion_mass, 6) * s -
                1. * pow(pion_mass, 4) * pow(m_rho, 2) * s -
                16. * pow(pion_mass, 2) * pow(m_rho, 4) * s -
                7. * pow(m_rho, 6) * s + 2. * pow(pion_mass, 4) * pow(s, 2) +
                14. * pow(pion_mass, 2) * pow(m_rho, 2) * pow(s, 2) +
                9. * pow(m_rho, 4) * pow(s, 2) -
                4. * pow(pion_mass, 2) * pow(s, 3) -
                5. * pow(m_rho, 2) * pow(s, 3) + 1. * pow(s, 4) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (2. * pow(a1_mass, 4) + 2. * pow(pion_mass, 4) +
                     3. * pow(m_rho, 4) +
                     pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 4. * s) -
                     5. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                     pow(a1_mass, 2) * (-4. * pow(pion_mass, 2) -
                                        2. * pow(m_rho, 2) + 2. * s)) +
                pow(a1_mass, 4) *
                    (-42. * pow(pion_mass, 4) - 3. * pow(m_rho, 4) +
                     9. * pow(m_rho, 2) * s - 6. * pow(s, 2) +
                     pow(pion_mass, 2) * (-42. * pow(m_rho, 2) + 36. * s)) +
                pow(a1_mass, 2) *
                    (28. * pow(pion_mass, 6) - 4. * pow(m_rho, 6) +
                     pow(pion_mass, 4) * (42. * pow(m_rho, 2) - 30. * s) +
                     9. * pow(m_rho, 4) * s - 6. * pow(m_rho, 2) * pow(s, 2) +
                     1. * pow(s, 3) +
                     pow(pion_mass, 2) *
                         (6. * pow(m_rho, 4) - 12. * pow(m_rho, 2) * s +
                          6. * pow(s, 2))))) *
          atan((pow(a1_mass, 2) - 2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) +
                s + tmax) /
               (Gammaa1 * a1_mass))) /
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + 4. * pow(a1_mass, 4) +
              4. * pow(pion_mass, 4) + 4. * pow(pion_mass, 2) * pow(m_rho, 2) +
              pow(m_rho, 4) - 4. * pow(pion_mass, 2) * s -
              2. * pow(m_rho, 2) * s + pow(s, 2) +
              pow(a1_mass, 2) *
                  (-8. * pow(pion_mass, 2) - 4. * pow(m_rho, 2) + 4. * s)) +
         0.0625 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-4. * pow(a1_mass, 6) +
                   pow(a1_mass, 4) * (12. * pow(pion_mass, 2) - 6. * s) +
                   pow(pion_mass, 2) *
                       (4. * pow(pion_mass, 4) - 4. * pow(m_rho, 4) -
                        2. * pow(pion_mass, 2) * s + 2. * pow(m_rho, 2) * s) +
                   pow(a1_mass, 2) *
                       (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                        8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                        4. * pow(s, 2))) +
              pow(eta1, 2) *
                  (2. * pow(a1_mass, 6) - 2. * pow(pion_mass, 6) +
                   pow(pion_mass, 2) * pow(m_rho, 2) * s +
                   pow(m_rho, 2) * (pow(m_rho, 2) - 1. * s) * s +
                   pow(pion_mass, 4) * (-3. * pow(m_rho, 2) + s) +
                   pow(a1_mass, 4) *
                       (-6. * pow(pion_mass, 2) - 3. * pow(m_rho, 2) + 3. * s) +
                   pow(a1_mass, 2) *
                       (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                        pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                        4. * pow(m_rho, 2) * s + 2. * pow(s, 2))) +
              pow(eta2, 2) *
                  (2. * pow(a1_mass, 6) - 2. * pow(pion_mass, 6) -
                   1. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                   pow(pion_mass, 4) * (3. * pow(m_rho, 2) + s) +
                   pow(a1_mass, 4) *
                       (-6. * pow(pion_mass, 2) + 3. * pow(m_rho, 2) + 3. * s) +
                   pow(a1_mass, 2) *
                       (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                        pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                        2. * pow(m_rho, 2) * s + 2. * pow(s, 2)))) *
             log(abs(-1. * pow(a1_mass, 2) + tmax)) -
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 *
               (-0.5 * pow(a1_mass, 6) - 0.5 * pow(pion_mass, 6) +
                0.5 * pow(pion_mass, 4) * pow(m_rho, 2) +
                pow(a1_mass, 4) *
                    (0.5 * pow(pion_mass, 2) + 0.5 * pow(m_rho, 2) - 1. * s) +
                pow(a1_mass, 2) * pow(pion_mass, 2) *
                    (0.5 * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 1. * s)) +
           eta1 * (pow(a1_mass, 4) * (1. * pow(pion_mass, 2) + 0.5 * s) +
                   pow(pion_mass, 2) *
                       (1. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                        pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 0.5 * s) -
                        0.5 * pow(m_rho, 2) * s) +
                   pow(a1_mass, 2) *
                       (-2. * pow(pion_mass, 4) - 0.5 * pow(m_rho, 2) * s +
                        pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s)))) *
          log(abs(-1. * pow(a1_mass, 2) + tmax))) /
             (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) +
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * (-0.5 * pow(a1_mass, 6) - 0.5 * pow(pion_mass, 6) +
                   pow(a1_mass, 4) *
                       (0.5 * pow(pion_mass, 2) + 0.5 * pow(m_rho, 2)) +
                   pow(pion_mass, 4) * (0.5 * pow(m_rho, 2) - 1. * s) +
                   pow(pion_mass, 2) * (-0.5 * pow(m_rho, 2) + 0.5 * s) * s +
                   pow(a1_mass, 2) *
                       (0.5 * pow(pion_mass, 4) +
                        pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) +
                        (-0.5 * pow(m_rho, 2) + 0.5 * s) * s)) +
           eta1 * (1. * pow(pion_mass, 6) +
                   pow(a1_mass, 4) * (1. * pow(pion_mass, 2) - 0.5 * s) +
                   pow(pion_mass, 2) * (1.5 * pow(m_rho, 2) - 2. * s) * s +
                   (-0.5 * pow(m_rho, 2) + 0.5 * s) * pow(s, 2) +
                   pow(pion_mass, 4) * (-1. * pow(m_rho, 2) + 1.5 * s) +
                   pow(a1_mass, 2) *
                       (-2. * pow(pion_mass, 4) + 0.5 * pow(m_rho, 2) * s +
                        pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s)))) *
          log(abs(-1. * pow(a1_mass, 2) + tmax))) /
             (1. * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
              1. * pow(m_rho, 2) + 1. * s) -
         (0.03125 * pow(eta1 - 1. * eta2, 2) *
          (1. * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) - 0.5 * pow(m_rho, 2) +
           0.5 * s) *
          (eta1 * eta2 *
               (-2. * pow(a1_mass, 8) +
                pow(a1_mass, 6) *
                    (8. * pow(pion_mass, 2) + 4. * pow(m_rho, 2) - 4. * s) +
                pow(a1_mass, 4) *
                    (-12. * pow(pion_mass, 4) - 4. * pow(m_rho, 4) +
                     4. * pow(m_rho, 2) * s +
                     pow(pion_mass, 2) * (-12. * pow(m_rho, 2) + 8. * s)) +
                pow(pion_mass, 2) *
                    (-2. * pow(pion_mass, 6) -
                     4. * pow(pion_mass, 4) * pow(m_rho, 2) -
                     2. * pow(m_rho, 6) + 4. * pow(m_rho, 4) * s -
                     2. * pow(m_rho, 2) * pow(s, 2) +
                     pow(pion_mass, 2) *
                         (-8. * pow(m_rho, 4) + 8. * pow(m_rho, 2) * s)) +
                pow(a1_mass, 2) *
                    (8. * pow(pion_mass, 6) + 2. * pow(m_rho, 6) +
                     pow(pion_mass, 4) * (12. * pow(m_rho, 2) - 4. * s) -
                     2. * pow(m_rho, 4) * s - 2. * pow(m_rho, 2) * pow(s, 2) +
                     2. * pow(s, 3) +
                     pow(pion_mass, 2) *
                         (8. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s -
                          4. * pow(s, 2)))) +
           pow(eta2, 2) *
               (pow(a1_mass, 8) +
                pow(a1_mass, 6) *
                    (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
                pow(pion_mass, 4) * (pow(pion_mass, 4) +
                                     2. * pow(pion_mass, 2) * pow(m_rho, 2) +
                                     pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
                pow(a1_mass, 4) *
                    (6. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                     pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) +
                     pow(m_rho, 2) * s) +
                pow(a1_mass, 2) *
                    (-4. * pow(pion_mass, 6) + 2. * pow(m_rho, 6) -
                     5. * pow(m_rho, 4) * s + 4. * pow(m_rho, 2) * pow(s, 2) -
                     1. * pow(s, 3) +
                     pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 2. * s) +
                     pow(pion_mass, 2) *
                         (2. * pow(m_rho, 4) - 4. * pow(m_rho, 2) * s +
                          2. * pow(s, 2)))) +
           pow(eta1, 2) *
               (pow(a1_mass, 8) + pow(pion_mass, 8) +
                2. * pow(pion_mass, 6) * pow(m_rho, 2) +
                pow(pion_mass, 2) * pow(m_rho, 2) * s *
                    (-2. * pow(m_rho, 2) + 2. * s) +
                pow(a1_mass, 6) *
                    (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
                pow(pion_mass, 4) *
                    (3. * pow(m_rho, 4) - 5. * pow(m_rho, 2) * s) +
                pow(a1_mass, 4) *
                    (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                     pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                     3. * pow(m_rho, 2) * s) +
                pow(a1_mass, 2) *
                    (-4. * pow(pion_mass, 6) + pow(m_rho, 4) * s -
                     1. * pow(s, 3) +
                     pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 2. * s) +
                     pow(pion_mass, 2) *
                         (-2. * pow(m_rho, 4) + 4. * pow(m_rho, 2) * s +
                          2. * pow(s, 2))))) *
          log(abs(-1. * pow(a1_mass, 2) + tmax))) /
             (0.25 * pow(Gammaa1, 2) * pow(a1_mass, 2) + 1. * pow(a1_mass, 4) +
              1. * pow(pion_mass, 4) + 1. * pow(pion_mass, 2) * pow(m_rho, 2) +
              0.25 * pow(m_rho, 4) - 1. * pow(pion_mass, 2) * s -
              0.5 * pow(m_rho, 2) * s + 0.25 * pow(s, 2) +
              pow(a1_mass, 2) *
                  (-2. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + 1. * s)) -
         (1. * (1. * eta1 - 1. * eta2) *
          (eta2 *
               (pow(a1_mass, 4) *
                    (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                pow(pion_mass, 2) * pow(m_rho, 2) *
                    (pow(pion_mass, 2) * (0.5 - 1. * C4 * pow(m_rho, 2)) +
                     (-0.25 + 0.125 * delta) * (pow(m_rho, 2) + s)) +
                pow(a1_mass, 2) *
                    (-1. * C4 * pow(m_rho, 6) +
                     pow(pion_mass, 2) *
                         (-1. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4)) +
                     0.25 * delta * pow(s, 2) +
                     pow(m_rho, 2) * s * (-0.25 - 0.375 * delta - 1. * C4 * s) +
                     pow(m_rho, 4) * (0.75 - 0.125 * delta + 2. * C4 * s))) +
           eta1 *
               (pow(a1_mass, 4) *
                    (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) +
                pow(a1_mass, 2) *
                    (0.5 * pow(m_rho, 4) - 1. * C4 * pow(m_rho, 6) +
                     pow(pion_mass, 2) *
                         (1. * pow(m_rho, 2) - 2. * C4 * pow(m_rho, 4)) -
                     0.25 * delta * pow(s, 2) +
                     1. * C4 * pow(m_rho, 2) * pow(s, 2)) +
                pow(m_rho, 2) *
                    (pow(pion_mass, 4) * (-0.5 + 1. * C4 * pow(m_rho, 2)) +
                     s * ((0.25 - 0.125 * delta) * pow(m_rho, 2) +
                          (-0.25 + 0.125 * delta) * s) +
                     pow(pion_mass, 2) *
                         (2. * C4 * pow(m_rho, 4) + (0.5 + 0.25 * delta) * s +
                          pow(m_rho, 2) * (-1. - 2. * C4 * s))))) *
          log(abs(-1. * pow(a1_mass, 2) + tmax))) /
             pow(m_rho, 2) +
         0.5 * pow(-2. + delta, 2) * pow(pion_mass, 2) *
             log(abs(-1. * pow(pion_mass, 2) + tmax)) +
         (0.25 *
          (0. +
           8.000000000000002 * pow(2. - 1. * delta, 2) * pow(pion_mass, 4) *
               pow(m_rho, 2) -
           5.999999999999999 * pow(2. - 1. * delta, 2) * pow(pion_mass, 2) *
               pow(m_rho, 2) * s +
           1. * pow(2. - 1. * delta, 2) * pow(m_rho, 2) * pow(s, 2)) *
          log(abs(-1. * pow(pion_mass, 2) + tmax))) /
             (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
         (0.125 * (-2. + delta) * (eta1 - 1. * eta2) * pow(pion_mass, 2) *
          (0. + eta2 * pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 4. * s) +
           eta1 * (2. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s +
                   pow(pion_mass, 2) * (-4. * pow(m_rho, 2) + 4. * s))) *
          log(abs(-1. * pow(pion_mass, 2) + tmax))) /
             (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) +
         (2. * (-2. + 1. * delta) *
          (0. + (-0.25 + 0.125 * delta) * pow(m_rho, 2) * s +
           pow(pion_mass, 2) * (-2. * C4 * pow(m_rho, 4) - 0.5 * delta * s +
                                pow(m_rho, 2) * (1. + 2. * C4 * s))) *
          log(abs(-1. * pow(pion_mass, 2) + tmax))) /
             pow(m_rho, 2) -
         (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(pion_mass, 2) *
          (eta1 * (pow(pion_mass, 4) * (-1. * pow(m_rho, 2) + 1. * s) +
                   pow(a1_mass, 2) *
                       (pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) +
                        (-0.5 * pow(m_rho, 2) + 0.5 * s) * s) +
                   pow(pion_mass, 2) *
                       (-1. * pow(m_rho, 4) + 2.5 * pow(m_rho, 2) * s -
                        1.5 * pow(s, 2)) +
                   s * (0.5 * pow(m_rho, 4) - 1. * pow(m_rho, 2) * s +
                        0.5 * pow(s, 2))) +
           eta2 *
               (0.5 * pow(m_rho, 6) +
                pow(pion_mass, 4) * (1. * pow(m_rho, 2) - 1. * s) -
                1.5 * pow(m_rho, 4) * s + 1.5 * pow(m_rho, 2) * pow(s, 2) -
                0.5 * pow(s, 3) +
                pow(pion_mass, 2) * (1.5 * pow(m_rho, 4) -
                                     3. * pow(m_rho, 2) * s + 1.5 * pow(s, 2)) +
                pow(a1_mass, 2) *
                    (-0.5 * pow(m_rho, 4) + 1. * pow(m_rho, 2) * s -
                     0.5 * pow(s, 2) +
                     pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s)))) *
          log(abs(-1. * pow(pion_mass, 2) + tmax))) /
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
              pow(pion_mass, 4) + 2. * pow(pion_mass, 2) * pow(m_rho, 2) +
              pow(m_rho, 4) - 2. * pow(pion_mass, 2) * s -
              2. * pow(m_rho, 2) * s + pow(s, 2) +
              pow(a1_mass, 2) *
                  (-2. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s)) -
         0.5 * pow(-2. + delta, 2) * pow(pion_mass, 2) *
             log(abs(-1. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + s + tmax)) +
         (0.5 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * pow(pion_mass, 6) * (1. * pow(m_rho, 2) - 1. * s) +
           eta2 * pow(a1_mass, 2) * pow(pion_mass, 4) *
               (-1. * pow(m_rho, 2) + 1. * s) +
           eta1 * pow(a1_mass, 2) * pow(pion_mass, 2) *
               (-0.5 * pow(m_rho, 4) +
                pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) +
                0.5 * pow(m_rho, 2) * s) +
           eta1 * pow(pion_mass, 4) *
               (0.5 * pow(m_rho, 4) - 0.5 * pow(m_rho, 2) * s +
                pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s))) *
          log(abs(-1. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + s + tmax))) /
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * pow(pion_mass, 2) + pow(pion_mass, 4)) +
         (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(pion_mass, 2) *
          (eta1 * (pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) +
                   (-0.5 * pow(m_rho, 2) + 0.5 * s) * s) +
           eta2 * (-0.5 * pow(m_rho, 4) + 1. * pow(m_rho, 2) * s -
                   0.5 * pow(s, 2) +
                   pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s))) *
          log(abs(-1. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + s + tmax))) /
             (1. * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
              1. * pow(m_rho, 2) + 1. * s) -
         (0.25 *
          (0. +
           8.000000000000002 * pow(2. - 1. * delta, 2) * pow(pion_mass, 4) *
               pow(m_rho, 2) +
           1. * pow(2. - 1. * delta, 2) * pow(m_rho, 4) * s +
           pow(pion_mass, 2) *
               (C4 * (32. - 16. * delta) * pow(m_rho, 6) +
                delta * (-8. + 4. * delta) * pow(s, 2) +
                pow(m_rho, 2) * s *
                    (-8. + 24. * delta - 10. * pow(delta, 2) + 32. * C4 * s -
                     16. * C4 * delta * s) +
                pow(m_rho, 4) * (-16. + 8. * delta - 64. * C4 * s +
                                 32. * C4 * delta * s))) *
          log(abs(-1. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) + s + tmax))) /
             (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
         0.03125 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (4. * pow(a1_mass, 6) +
                   pow(Gammaa1, 2) * pow(a1_mass, 2) *
                       (-4. * pow(a1_mass, 2) + 4. * pow(pion_mass, 2) -
                        2. * s) +
                   pow(a1_mass, 4) * (-12. * pow(pion_mass, 2) + 6. * s) +
                   pow(pion_mass, 2) *
                       (-4. * pow(pion_mass, 4) + 4. * pow(m_rho, 4) +
                        2. * pow(pion_mass, 2) * s - 2. * pow(m_rho, 2) * s) +
                   pow(a1_mass, 2) *
                       (12. * pow(pion_mass, 4) - 2. * pow(m_rho, 4) -
                        8. * pow(pion_mass, 2) * s - 4. * pow(m_rho, 2) * s +
                        4. * pow(s, 2))) +
              pow(eta1, 2) *
                  (-2. * pow(a1_mass, 6) + 2. * pow(pion_mass, 6) +
                   3. * pow(pion_mass, 4) * pow(m_rho, 2) +
                   pow(a1_mass, 4) *
                       (6. * pow(pion_mass, 2) + 3. * pow(m_rho, 2) - 3. * s) -
                   1. * pow(pion_mass, 4) * s -
                   1. * pow(pion_mass, 2) * pow(m_rho, 2) * s -
                   1. * pow(m_rho, 4) * s + pow(m_rho, 2) * pow(s, 2) +
                   pow(Gammaa1, 2) * pow(a1_mass, 2) *
                       (2. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                        1. * pow(m_rho, 2) + s) +
                   pow(a1_mass, 2) *
                       (-6. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                        4. * pow(m_rho, 2) * s - 2. * pow(s, 2) +
                        pow(pion_mass, 2) * (-6. * pow(m_rho, 2) + 4. * s))) +
              pow(eta2, 2) *
                  (-2. * pow(a1_mass, 6) + 2. * pow(pion_mass, 6) -
                   3. * pow(pion_mass, 4) * pow(m_rho, 2) +
                   pow(a1_mass, 4) *
                       (6. * pow(pion_mass, 2) - 3. * pow(m_rho, 2) - 3. * s) -
                   1. * pow(pion_mass, 4) * s +
                   pow(pion_mass, 2) * pow(m_rho, 2) * s +
                   pow(Gammaa1, 2) * pow(a1_mass, 2) *
                       (2. * pow(a1_mass, 2) - 2. * pow(pion_mass, 2) +
                        pow(m_rho, 2) + s) +
                   pow(a1_mass, 2) *
                       (-6. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                        2. * pow(m_rho, 2) * s - 2. * pow(s, 2) +
                        pow(pion_mass, 2) * (6. * pow(m_rho, 2) + 4. * s)))) *
             log(abs(pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
                     4. * pow(a1_mass, 2) * pow(pion_mass, 2) +
                     4. * pow(pion_mass, 4) +
                     2. * pow(a1_mass, 2) * (-1. * pow(m_rho, 2) + s + tmax) -
                     4. * pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + s + tmax) +
                     pow(-1. * pow(m_rho, 2) + s + tmax, 2))) -
         (0.5 * (1. * eta1 - 1. * eta2) *
          (eta2 * (pow(Gammaa1, 2) * pow(a1_mass, 2) * pow(m_rho, 2) *
                       (0.5 - 1. * C4 * pow(m_rho, 2)) +
                   pow(a1_mass, 4) *
                       (-0.5 * pow(m_rho, 2) + 1. * C4 * pow(m_rho, 4)) +
                   pow(pion_mass, 2) * pow(m_rho, 2) *
                       (pow(pion_mass, 2) * (-0.5 + 1. * C4 * pow(m_rho, 2)) +
                        (0.25 - 0.125 * delta) * (pow(m_rho, 2) + s)) +
                   pow(a1_mass, 2) *
                       (1. * C4 * pow(m_rho, 6) +
                        pow(pion_mass, 2) *
                            (1. * pow(m_rho, 2) - 2. * C4 * pow(m_rho, 4)) -
                        0.25 * delta * pow(s, 2) +
                        pow(m_rho, 4) * (-0.75 + 0.125 * delta - 2. * C4 * s) +
                        pow(m_rho, 2) * s *
                            (0.25 + 0.375 * delta + 1. * C4 * s))) +
           eta1 *
               (pow(Gammaa1, 2) * pow(a1_mass, 2) * pow(m_rho, 2) *
                    (-0.5 + 1. * C4 * pow(m_rho, 2)) +
                pow(a1_mass, 4) *
                    (0.5 * pow(m_rho, 2) - 1. * C4 * pow(m_rho, 4)) +
                pow(a1_mass, 2) *
                    (-0.5 * pow(m_rho, 4) + 1. * C4 * pow(m_rho, 6) +
                     pow(pion_mass, 2) *
                         (-1. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4)) +
                     0.25 * delta * pow(s, 2) -
                     1. * C4 * pow(m_rho, 2) * pow(s, 2)) +
                pow(m_rho, 2) *
                    (pow(pion_mass, 4) * (0.5 - 1. * C4 * pow(m_rho, 2)) +
                     s * ((-0.25 + 0.125 * delta) * pow(m_rho, 2) +
                          (0.25 - 0.125 * delta) * s) +
                     pow(pion_mass, 2) *
                         (-2. * C4 * pow(m_rho, 4) + (-0.5 - 0.25 * delta) * s +
                          pow(m_rho, 2) * (1. + 2. * C4 * s))))) *
          log(abs(pow(Gammaa1, 2) * pow(a1_mass, 2) +
                  pow(pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                          1. * pow(m_rho, 2) + s + tmax,
                      2)))) /
             pow(m_rho, 2) +
         (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 *
               (0.5 * pow(Gammaa1, 4) * pow(a1_mass, 4) -
                0.5 * pow(a1_mass, 8) + 0.5 * pow(pion_mass, 8) +
                0.5 * pow(a1_mass, 4) * pow(pion_mass, 2) * pow(m_rho, 2) -
                0.5 * pow(pion_mass, 6) * pow(m_rho, 2) +
                pow(Gammaa1, 2) *
                    (pow(a1_mass, 2) * pow(pion_mass, 2) *
                         (1. * pow(pion_mass, 2) + 1.5 * pow(m_rho, 2) -
                          2. * s) +
                     pow(a1_mass, 4) * (-1. * pow(pion_mass, 2) +
                                        0.5 * pow(m_rho, 2) - 1. * s)) +
                pow(a1_mass, 6) *
                    (1. * pow(pion_mass, 2) + 0.5 * pow(m_rho, 2) - 1. * s) +
                pow(a1_mass, 2) * pow(pion_mass, 4) *
                    (-1. * pow(pion_mass, 2) - 0.5 * pow(m_rho, 2) + 1. * s)) +
           eta1 * (pow(a1_mass, 6) * (1. * pow(pion_mass, 2) + 0.5 * s) +
                   pow(a1_mass, 2) * (3. * pow(pion_mass, 6) +
                                      1. * pow(pion_mass, 2) * pow(m_rho, 4) -
                                      0.5 * pow(pion_mass, 4) * s) +
                   pow(a1_mass, 4) *
                       (-3. * pow(pion_mass, 4) +
                        pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 0.5 * s) -
                        0.5 * pow(m_rho, 2) * s) +
                   pow(pion_mass, 4) *
                       (-1. * pow(pion_mass, 4) - 1. * pow(m_rho, 4) +
                        pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 0.5 * s) +
                        0.5 * pow(m_rho, 2) * s) +
                   pow(Gammaa1, 2) * pow(a1_mass, 2) *
                       (-1. * pow(pion_mass, 4) +
                        pow(a1_mass, 2) * (1. * pow(pion_mass, 2) + 0.5 * s) -
                        0.5 * pow(m_rho, 2) * s +
                        pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1.5 * s)))) *
          log(abs(pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
                  4. * pow(pion_mass, 4) +
                  4. * pow(pion_mass, 2) * pow(m_rho, 2) + pow(m_rho, 4) -
                  4. * pow(pion_mass, 2) * s - 2. * pow(m_rho, 2) * s +
                  pow(s, 2) - 4. * pow(pion_mass, 2) * tmax -
                  2. * pow(m_rho, 2) * tmax + 2. * s * tmax + pow(tmax, 2) +
                  pow(a1_mass, 2) *
                      (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s +
                       2. * tmax)))) /
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) -
              2. * pow(a1_mass, 2) * pow(pion_mass, 2) + pow(pion_mass, 4)) -
         (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta1 * (-1. * pow(pion_mass, 8) +
                   1. * pow(pion_mass, 4) * pow(m_rho, 4) +
                   pow(a1_mass, 6) * (1. * pow(pion_mass, 2) - 0.5 * s) -
                   0.5 * pow(pion_mass, 6) * s -
                   4. * pow(pion_mass, 4) * pow(m_rho, 2) * s -
                   1.5 * pow(pion_mass, 2) * pow(m_rho, 4) * s +
                   3.5 * pow(pion_mass, 4) * pow(s, 2) +
                   4. * pow(pion_mass, 2) * pow(m_rho, 2) * pow(s, 2) +
                   0.5 * pow(m_rho, 4) * pow(s, 2) -
                   2.5 * pow(pion_mass, 2) * pow(s, 3) -
                   1. * pow(m_rho, 2) * pow(s, 3) + 0.5 * pow(s, 4) +
                   pow(Gammaa1, 2) * pow(a1_mass, 2) *
                       (-1. * pow(pion_mass, 4) +
                        pow(a1_mass, 2) * (1. * pow(pion_mass, 2) - 0.5 * s) -
                        0.5 * pow(pion_mass, 2) * s + 0.5 * pow(s, 2)) +
                   pow(a1_mass, 4) *
                       (-3. * pow(pion_mass, 4) +
                        (1. * pow(m_rho, 2) - 0.5 * s) * s +
                        pow(pion_mass, 2) * (-2. * pow(m_rho, 2) + 2.5 * s)) +
                   pow(a1_mass, 2) *
                       (3. * pow(pion_mass, 6) +
                        pow(pion_mass, 4) * (2. * pow(m_rho, 2) - 1.5 * s) -
                        0.5 * pow(m_rho, 4) * s + 0.5 * pow(s, 3) +
                        pow(pion_mass, 2) *
                            (1. * pow(m_rho, 4) - 1. * pow(m_rho, 2) * s -
                             1. * pow(s, 2)))) +
           eta2 *
               (0.5 * pow(Gammaa1, 4) * pow(a1_mass, 4) -
                0.5 * pow(a1_mass, 8) + 0.5 * pow(pion_mass, 8) +
                pow(a1_mass, 6) *
                    (1. * pow(pion_mass, 2) + 1. * pow(m_rho, 2) - 0.5 * s) +
                0.5 * pow(pion_mass, 6) * s +
                pow(a1_mass, 4) * (-0.5 * pow(m_rho, 4) +
                                   (-0.5 * pow(pion_mass, 2) + 0.5 * s) * s) +
                pow(pion_mass, 4) * (-0.5 * pow(m_rho, 4) +
                                     2. * pow(m_rho, 2) * s - 1.5 * pow(s, 2)) +
                pow(pion_mass, 2) * s *
                    (0.5 * pow(m_rho, 4) - 1. * pow(m_rho, 2) * s +
                     0.5 * pow(s, 2)) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (1. * pow(pion_mass, 4) + 0.5 * pow(m_rho, 4) +
                     pow(pion_mass, 2) * (2. * pow(m_rho, 2) - 1.5 * s) -
                     1. * pow(m_rho, 2) * s + 0.5 * pow(s, 2) +
                     pow(a1_mass, 2) * (-1. * pow(pion_mass, 2) -
                                        1. * pow(m_rho, 2) + 1.5 * s)) +
                pow(a1_mass, 2) *
                    (-1. * pow(pion_mass, 6) +
                     pow(pion_mass, 4) * (-1. * pow(m_rho, 2) + 0.5 * s) +
                     pow(pion_mass, 2) *
                         (-1. * pow(m_rho, 4) + 2. * pow(m_rho, 2) * s -
                          1. * pow(s, 2)) +
                     s * (0.5 * pow(m_rho, 4) - 1. * pow(m_rho, 2) * s +
                          0.5 * pow(s, 2))))) *
          log(abs(pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
                  4. * pow(pion_mass, 4) +
                  4. * pow(pion_mass, 2) * pow(m_rho, 2) + pow(m_rho, 4) -
                  4. * pow(pion_mass, 2) * s - 2. * pow(m_rho, 2) * s +
                  pow(s, 2) - 4. * pow(pion_mass, 2) * tmax -
                  2. * pow(m_rho, 2) * tmax + 2. * s * tmax + pow(tmax, 2) +
                  pow(a1_mass, 2) *
                      (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s +
                       2. * tmax)))) /
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
              pow(pion_mass, 4) + 2. * pow(pion_mass, 2) * pow(m_rho, 2) +
              pow(m_rho, 4) - 2. * pow(pion_mass, 2) * s -
              2. * pow(m_rho, 2) * s + pow(s, 2) +
              pow(a1_mass, 2) *
                  (-2. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s)) -
         (0.0625 * pow(eta1 - 1. * eta2, 2) *
          (pow(eta2, 2) *
               (-1. * pow(a1_mass, 10) +
                pow(a1_mass, 8) *
                    (5. * pow(pion_mass, 2) + 2.5 * pow(m_rho, 2) - 2.5 * s) +
                pow(Gammaa1, 4) * pow(a1_mass, 4) *
                    (1. * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
                     0.5 * pow(m_rho, 2) + 0.5 * s) +
                pow(a1_mass, 4) *
                    (10. * pow(pion_mass, 6) - 2.5 * pow(m_rho, 6) +
                     pow(pion_mass, 4) * (15. * pow(m_rho, 2) - 9. * s) +
                     6. * pow(m_rho, 4) * s - 4.5 * pow(m_rho, 2) * pow(s, 2) +
                     1. * pow(s, 3)) +
                pow(a1_mass, 6) *
                    (-10. * pow(pion_mass, 4) +
                     (1. * pow(m_rho, 2) - 1. * s) * s +
                     pow(pion_mass, 2) * (-10. * pow(m_rho, 2) + 8. * s)) +
                pow(pion_mass, 4) *
                    (1. * pow(pion_mass, 6) + 0.5 * pow(m_rho, 6) +
                     pow(pion_mass, 4) * (2.5 * pow(m_rho, 2) - 0.5 * s) -
                     1. * pow(m_rho, 4) * s + 0.5 * pow(m_rho, 2) * pow(s, 2) +
                     pow(pion_mass, 2) *
                         (2. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s)) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (4. * pow(a1_mass, 6) - 4. * pow(pion_mass, 6) -
                     0.5 * pow(m_rho, 6) + 1.5 * pow(m_rho, 4) * s -
                     1.5 * pow(m_rho, 2) * pow(s, 2) + 0.5 * pow(s, 3) +
                     pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 6. * s) +
                     pow(a1_mass, 4) * (-12. * pow(pion_mass, 2) -
                                        6. * pow(m_rho, 2) + 6. * s) +
                     pow(pion_mass, 2) *
                         (-3. * pow(m_rho, 4) + 6. * pow(m_rho, 2) * s -
                          3. * pow(s, 2)) +
                     pow(a1_mass, 2) *
                         (12. * pow(pion_mass, 4) + 3. * pow(m_rho, 4) +
                          pow(pion_mass, 2) * (12. * pow(m_rho, 2) - 12. * s) -
                          6. * pow(m_rho, 2) * s + 3. * pow(s, 2))) +
                pow(a1_mass, 2) *
                    (-5. * pow(pion_mass, 8) + 1. * pow(m_rho, 8) -
                     3.5 * pow(m_rho, 6) * s + 4.5 * pow(m_rho, 4) * pow(s, 2) -
                     2.5 * pow(m_rho, 2) * pow(s, 3) + 0.5 * pow(s, 4) +
                     pow(pion_mass, 6) * (-10. * pow(m_rho, 2) + 4. * s) +
                     pow(pion_mass, 4) *
                         (-2. * pow(m_rho, 4) + 1. * pow(m_rho, 2) * s +
                          1. * pow(s, 2)) +
                     pow(pion_mass, 2) *
                         (3. * pow(m_rho, 6) - 8. * pow(m_rho, 4) * s +
                          7. * pow(m_rho, 2) * pow(s, 2) - 2. * pow(s, 3)))) +
           pow(eta1, 2) *
               (-1. * pow(a1_mass, 10) +
                pow(a1_mass, 8) *
                    (5. * pow(pion_mass, 2) + 2.5 * pow(m_rho, 2) - 2.5 * s) +
                pow(Gammaa1, 4) * pow(a1_mass, 4) *
                    (1. * pow(a1_mass, 2) - 1. * pow(pion_mass, 2) -
                     0.5 * pow(m_rho, 2) + 0.5 * s) +
                pow(a1_mass, 6) *
                    (-10. * pow(pion_mass, 4) - 2. * pow(m_rho, 4) +
                     5. * pow(m_rho, 2) * s - 1. * pow(s, 2) +
                     pow(pion_mass, 2) * (-10. * pow(m_rho, 2) + 8. * s)) +
                pow(a1_mass, 4) *
                    (10. * pow(pion_mass, 6) + 0.5 * pow(m_rho, 6) +
                     pow(pion_mass, 4) * (15. * pow(m_rho, 2) - 9. * s) -
                     3. * pow(m_rho, 4) * s + 1.5 * pow(m_rho, 2) * pow(s, 2) +
                     1. * pow(s, 3) +
                     pow(pion_mass, 2) *
                         (6. * pow(m_rho, 4) - 12. * pow(m_rho, 2) * s)) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (4. * pow(a1_mass, 6) - 4. * pow(pion_mass, 6) -
                     0.5 * pow(m_rho, 6) + 1.5 * pow(m_rho, 4) * s -
                     1.5 * pow(m_rho, 2) * pow(s, 2) + 0.5 * pow(s, 3) +
                     pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 6. * s) +
                     pow(a1_mass, 4) * (-12. * pow(pion_mass, 2) -
                                        6. * pow(m_rho, 2) + 6. * s) +
                     pow(pion_mass, 2) *
                         (-3. * pow(m_rho, 4) + 6. * pow(m_rho, 2) * s -
                          3. * pow(s, 2)) +
                     pow(a1_mass, 2) *
                         (12. * pow(pion_mass, 4) + 3. * pow(m_rho, 4) +
                          pow(pion_mass, 2) * (12. * pow(m_rho, 2) - 12. * s) -
                          6. * pow(m_rho, 2) * s + 3. * pow(s, 2))) +
                pow(pion_mass, 2) *
                    (1. * pow(pion_mass, 8) +
                     pow(pion_mass, 6) * (2.5 * pow(m_rho, 2) - 0.5 * s) +
                     pow(pion_mass, 4) *
                         (4. * pow(m_rho, 4) - 6. * pow(m_rho, 2) * s) +
                     pow(m_rho, 2) * s *
                         (-1. * pow(m_rho, 4) + 2. * pow(m_rho, 2) * s -
                          1. * pow(s, 2)) +
                     pow(pion_mass, 2) *
                         (1.5 * pow(m_rho, 6) - 6. * pow(m_rho, 4) * s +
                          4.5 * pow(m_rho, 2) * pow(s, 2))) +
                pow(a1_mass, 2) *
                    (-5. * pow(pion_mass, 8) +
                     pow(pion_mass, 6) * (-10. * pow(m_rho, 2) + 4. * s) +
                     pow(pion_mass, 4) *
                         (-8. * pow(m_rho, 4) + 13. * pow(m_rho, 2) * s +
                          1. * pow(s, 2)) +
                     pow(pion_mass, 2) *
                         (-1. * pow(m_rho, 6) + 6. * pow(m_rho, 4) * s -
                          3. * pow(m_rho, 2) * pow(s, 2) - 2. * pow(s, 3)) +
                     s * (0.5 * pow(m_rho, 6) - 0.5 * pow(m_rho, 4) * s -
                          0.5 * pow(m_rho, 2) * pow(s, 2) + 0.5 * pow(s, 3)))) +
           eta1 * eta2 *
               (2. * pow(a1_mass, 10) +
                pow(Gammaa1, 4) * pow(a1_mass, 4) *
                    (-2. * pow(a1_mass, 2) + 2. * pow(pion_mass, 2) +
                     1. * pow(m_rho, 2) - 1. * s) +
                pow(a1_mass, 8) *
                    (-10. * pow(pion_mass, 2) - 5. * pow(m_rho, 2) + 5. * s) +
                pow(a1_mass, 6) *
                    (20. * pow(pion_mass, 4) + 6. * pow(m_rho, 4) +
                     pow(pion_mass, 2) * (20. * pow(m_rho, 2) - 16. * s) -
                     8. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                pow(a1_mass, 4) *
                    (-20. * pow(pion_mass, 6) - 4. * pow(m_rho, 6) +
                     6. * pow(m_rho, 4) * s - 2. * pow(s, 3) +
                     pow(pion_mass, 4) * (-30. * pow(m_rho, 2) + 18. * s) +
                     pow(pion_mass, 2) *
                         (-18. * pow(m_rho, 4) + 18. * pow(m_rho, 2) * s)) +
                pow(pion_mass, 2) *
                    (-2. * pow(pion_mass, 8) - 1. * pow(m_rho, 8) +
                     3. * pow(m_rho, 6) * s - 3. * pow(m_rho, 4) * pow(s, 2) +
                     1. * pow(m_rho, 2) * pow(s, 3) +
                     pow(pion_mass, 6) * (-5. * pow(m_rho, 2) + 1. * s) +
                     pow(pion_mass, 4) *
                         (-10. * pow(m_rho, 4) + 10. * pow(m_rho, 2) * s) +
                     pow(pion_mass, 2) *
                         (-6. * pow(m_rho, 6) + 12. * pow(m_rho, 4) * s -
                          6. * pow(m_rho, 2) * pow(s, 2))) +
                pow(a1_mass, 2) *
                    (10. * pow(pion_mass, 8) + 1. * pow(m_rho, 8) +
                     pow(pion_mass, 6) * (20. * pow(m_rho, 2) - 8. * s) -
                     2. * pow(m_rho, 6) * s + 2. * pow(m_rho, 2) * pow(s, 3) -
                     1. * pow(s, 4) +
                     pow(pion_mass, 4) *
                         (22. * pow(m_rho, 4) - 20. * pow(m_rho, 2) * s -
                          2. * pow(s, 2)) +
                     pow(pion_mass, 2) *
                         (8. * pow(m_rho, 6) - 12. * pow(m_rho, 4) * s +
                          4. * pow(s, 3))) +
                pow(Gammaa1, 2) * pow(a1_mass, 2) *
                    (-8. * pow(a1_mass, 6) + 8. * pow(pion_mass, 6) +
                     1. * pow(m_rho, 6) +
                     pow(pion_mass, 4) * (12. * pow(m_rho, 2) - 12. * s) +
                     pow(a1_mass, 4) * (24. * pow(pion_mass, 2) +
                                        12. * pow(m_rho, 2) - 12. * s) -
                     3. * pow(m_rho, 4) * s + 3. * pow(m_rho, 2) * pow(s, 2) -
                     1. * pow(s, 3) +
                     pow(pion_mass, 2) *
                         (6. * pow(m_rho, 4) - 12. * pow(m_rho, 2) * s +
                          6. * pow(s, 2)) +
                     pow(a1_mass, 2) *
                         (-24. * pow(pion_mass, 4) - 6. * pow(m_rho, 4) +
                          12. * pow(m_rho, 2) * s - 6. * pow(s, 2) +
                          pow(pion_mass, 2) *
                              (-24. * pow(m_rho, 2) + 24. * s))))) *
          log(abs(pow(Gammaa1, 2) * pow(a1_mass, 2) + pow(a1_mass, 4) +
                  4. * pow(pion_mass, 4) +
                  4. * pow(pion_mass, 2) * pow(m_rho, 2) + pow(m_rho, 4) -
                  4. * pow(pion_mass, 2) * s - 2. * pow(m_rho, 2) * s +
                  pow(s, 2) - 4. * pow(pion_mass, 2) * tmax -
                  2. * pow(m_rho, 2) * tmax + 2. * s * tmax + pow(tmax, 2) +
                  pow(a1_mass, 2) *
                      (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s +
                       2. * tmax)))) /
             (pow(Gammaa1, 2) * pow(a1_mass, 2) + 4. * pow(a1_mass, 4) +
              4. * pow(pion_mass, 4) + 4. * pow(pion_mass, 2) * pow(m_rho, 2) +
              pow(m_rho, 4) - 4. * pow(pion_mass, 2) * s -
              2. * pow(m_rho, 2) * s + pow(s, 2) +
              pow(a1_mass, 2) *
                  (-8. * pow(pion_mass, 2) - 4. * pow(m_rho, 2) + 4. * s)))) /
           (16. * Pi * s * (-4 * pow(pion_mass, 2) + s)));

  // clang-format on
  return cut_off(xs * gev2_mb / spin_deg_factor);
}

double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi_pi_rho0(
    const double s, const double t, const double m_rho) {
  const double spin_deg_factor = 1.0;

  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  // clang-format off
  const double diff_xs =
      ((pow(Const, 2) * pow(ghat, 4) *
        ((0.25 *
          (32 * pow(C4, 2) * pow(m_rho, 8) + 2 * pow(delta, 2) * pow(s, 2) +
           8 * C4 * pow(m_rho, 6) * (-6 + delta - 8 * C4 * s) +
           2 * delta * pow(m_rho, 2) * s * (-6 + delta - 8 * C4 * s) +
           pow(m_rho, 4) * (12 - pow(delta, 2) + 8 * C4 * (6 + delta) * s +
                            32 * pow(C4, 2) * pow(s, 2)))) /
             pow(m_rho, 4) -
         (0.25 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
          (pow(pion_mass, 4) + pow(pow(m_rho, 2) - t, 2) -
           2 * pow(pion_mass, 2) * (pow(m_rho, 2) + t))) /
             (pow(m_rho, 2) * pow(pow(pion_mass, 2) - t, 2)) -
         (0.25 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
          (pow(pion_mass, 4) + pow(s + t, 2) -
           2 * pow(pion_mass, 2) * (2 * pow(m_rho, 2) + s + t))) /
             (pow(m_rho, 2) *
              pow(pow(pion_mass, 2) + pow(m_rho, 2) - s - t, 2)) +
         (0.125 * (-2 + delta) * (eta1 - eta2) *
          (pow(a1_mass, 2) - 2 * pow(pion_mass, 2) - pow(m_rho, 2) + s + t) *
          (eta1 * (2 * pow(pion_mass, 2) - s) +
           eta2 * (-3 * pow(pion_mass, 2) - pow(m_rho, 2) + s + t)) *
          (pow(pion_mass, 4) + t * (-pow(m_rho, 2) + 2 * s + t) -
           pow(pion_mass, 2) * (pow(m_rho, 2) + 2 * t))) /
             ((-pow(pion_mass, 2) + t) *
              (pow(Gammaa1, 2) * pow(a1_mass, 2) +
               pow(pow(a1_mass, 2) - 2 * pow(pion_mass, 2) - pow(m_rho, 2) + s +
                       t,
                   2))) -
         (0.125 * (-2 + delta) * (eta1 - eta2) *
          (eta1 * (-2 * pow(pion_mass, 2) + s) +
           eta2 * (pow(pion_mass, 2) + t)) *
          (-pow(pion_mass, 4) + pow(s, 2) - pow(t, 2) +
           pow(m_rho, 2) * (-s + t) +
           pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s + 2 * t))) /
             ((pow(pion_mass, 2) + pow(m_rho, 2) - s - t) *
              (-pow(a1_mass, 2) + t)) +
         (0.25 * (-2. + delta) *
          (pow(pion_mass, 4) * (2. + delta - 8. * C4 * pow(m_rho, 2)) +
           8. * C4 * pow(m_rho, 4) * t +
           t * ((2. + 3. * delta) * s + (2. + delta) * t) +
           pow(m_rho, 2) * (s * (2. - 1. * delta - 16. * C4 * t) +
                            t * (-2. - 1. * delta - 8. * C4 * t)) +
           pow(pion_mass, 2) *
               (8. * C4 * pow(m_rho, 4) + (-2. + delta) * s +
                (-4. - 2. * delta) * t +
                pow(m_rho, 2) * (-6. + delta + 16. * C4 * t)))) /
             (pow(m_rho, 2) * (pow(pion_mass, 2) - 1. * t)) -
         (0.125 * (-2 + delta) * (eta1 - eta2) *
          (pow(a1_mass, 2) - 2 * pow(pion_mass, 2) - pow(m_rho, 2) + s + t) *
          (-(eta2 * (3 * pow(pion_mass, 2) + pow(m_rho, 2) - s - t) *
             (pow(pion_mass, 4) + (pow(m_rho, 2) - s - t) * (s - t) -
              pow(pion_mass, 2) * (pow(m_rho, 2) - 2 * s + 2 * t))) +
           eta1 * (2 * pow(pion_mass, 6) +
                   pow(pion_mass, 4) * (-2 * pow(m_rho, 2) + 5 * s - 4 * t) +
                   s * (s + t) * (-pow(m_rho, 2) + s + t) +
                   pow(pion_mass, 2) *
                       (2 * pow(m_rho, 4) + pow(m_rho, 2) * (s - 2 * t) -
                        2 * (2 * s - t) * (s + t))))) /
             ((-pow(pion_mass, 2) - pow(m_rho, 2) + s + t) *
              (pow(Gammaa1, 2) * pow(a1_mass, 2) +
               pow(pow(a1_mass, 2) - 2 * pow(pion_mass, 2) - pow(m_rho, 2) + s +
                       t,
                   2))) +
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * (-0.5 * pow(pion_mass, 6) +
                   pow(pion_mass, 4) * (0.5 * pow(m_rho, 2) + 0.5 * t) +
                   pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s + 0.5 * t) *
                       t +
                   (0.5 * pow(m_rho, 2) - 1. * s - 0.5 * t) * pow(t, 2)) +
           eta1 *
               (1. * pow(pion_mass, 6) +
                pow(pion_mass, 4) * (-1. * pow(m_rho, 2) + 0.5 * s - 2. * t) +
                s * (-0.5 * pow(m_rho, 2) + 0.5 * t) * t +
                pow(pion_mass, 2) *
                    (1. * pow(m_rho, 4) + pow(m_rho, 2) * (-0.5 * s - 1. * t) +
                     t * (1. * s + 1. * t))))) /
             ((pow(a1_mass, 2) - 1. * t) * (-1. * pow(pion_mass, 2) + t)) +
         (0.03125 * pow(eta1 - eta2, 2) *
          (-2 * eta1 * eta2 *
               (pow(pion_mass, 8) - 4 * pow(pion_mass, 6) * t +
                pow(t, 2) * (-pow(m_rho, 4) - 2 * pow(m_rho, 2) * s +
                             2 * pow(s, 2) + 2 * s * t + pow(t, 2)) -
                2 * pow(pion_mass, 2) * t *
                    (-2 * pow(m_rho, 4) + pow(m_rho, 2) * s + 2 * t * (s + t)) +
                pow(pion_mass, 4) * (-pow(m_rho, 4) + 2 * t * (s + 3 * t))) +
           pow(eta2, 2) *
               (pow(pion_mass, 8) -
                2 * pow(pion_mass, 6) * (pow(m_rho, 2) + 2 * t) +
                pow(t, 2) * (pow(m_rho, 4) + 2 * pow(s, 2) + 2 * s * t +
                             pow(t, 2) + 2 * pow(m_rho, 2) * (-s + t)) -
                2 * pow(pion_mass, 2) * t *
                    (2 * t * (s + t) + pow(m_rho, 2) * (s + 3 * t)) +
                pow(pion_mass, 4) * (pow(m_rho, 4) + 6 * pow(m_rho, 2) * t +
                                     2 * t * (s + 3 * t))) +
           pow(eta1, 2) *
               (pow(pion_mass, 8) +
                2 * pow(pion_mass, 6) * (pow(m_rho, 2) - 2 * t) -
                2 * pow(pion_mass, 2) * (pow(m_rho, 2) - s - t) *
                    (pow(m_rho, 4) + pow(m_rho, 2) * t - 2 * pow(t, 2)) +
                t * (-pow(m_rho, 2) + t) *
                    (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                     pow(m_rho, 2) * (2 * s + t)) +
                pow(pion_mass, 4) *
                    (pow(m_rho, 4) - 2 * pow(m_rho, 2) * (s + 3 * t) +
                     2 * t * (s + 3 * t))))) /
             pow(pow(a1_mass, 2) - t, 2) +
         ((0.25 * pow(-2 + delta, 2) * (2 * pow(pion_mass, 2) - s) *
           (pow(pion_mass, 4) + pow(m_rho, 2) * (s - t) + t * (s + t) -
            pow(pion_mass, 2) * (3 * pow(m_rho, 2) + s + 2 * t))) /
              ((pow(pion_mass, 2) - t) *
               (pow(pion_mass, 2) + pow(m_rho, 2) - s - t)) -
          (0.25 * (-2. + delta) *
           (pow(pion_mass, 4) * (2. + delta - 8. * C4 * pow(m_rho, 2)) -
            2. * delta * pow(s, 2) + 2. * s * t - 1. * delta * s * t +
            2. * pow(t, 2) + delta * pow(t, 2) +
            C4 * pow(m_rho, 4) * (-8. * s + 8. * t) +
            pow(m_rho, 2) * ((2. + delta) * s + 8. * C4 * pow(s, 2) +
                             t * (-2. - 1. * delta - 8. * C4 * t)) +
            pow(pion_mass, 2) * (8. * C4 * pow(m_rho, 4) - 2. * s +
                                 5. * delta * s - 4. * t - 2. * delta * t +
                                 pow(m_rho, 2) * (-6. + delta - 16. * C4 * s +
                                                  16. * C4 * t)))) /
              (pow(pion_mass, 2) + pow(m_rho, 2) - 1. * s - 1. * t)) /
             pow(m_rho, 2) +
         (0.03125 * pow(eta1 - eta2, 2) *
          (-2 * eta1 * eta2 *
               (pow(pion_mass, 8) +
                4 * pow(pion_mass, 6) * (pow(m_rho, 2) - t) +
                pow(-pow(m_rho, 2) + s + t, 2) *
                    (pow(s, 2) + pow(t, 2) - 2 * pow(m_rho, 2) * (s + t)) +
                pow(pion_mass, 4) *
                    (9 * pow(m_rho, 4) + 4 * pow(s, 2) + 2 * s * t +
                     6 * pow(t, 2) - 2 * pow(m_rho, 2) * (7 * s + 6 * t)) +
                2 * pow(pion_mass, 2) * (pow(m_rho, 2) - s - t) *
                    (2 * pow(m_rho, 4) - pow(m_rho, 2) * (5 * s + 4 * t) +
                     2 * (pow(s, 2) + pow(t, 2)))) +
           pow(eta2, 2) *
               (pow(pion_mass, 8) +
                pow(pion_mass, 6) * (6 * pow(m_rho, 2) - 4 * t) +
                pow(-pow(m_rho, 2) + s + t, 2) *
                    (4 * pow(m_rho, 4) + pow(s, 2) + pow(t, 2) -
                     4 * pow(m_rho, 2) * (s + t)) +
                pow(pion_mass, 4) *
                    (17 * pow(m_rho, 4) + 4 * pow(s, 2) + 2 * s * t +
                     6 * pow(t, 2) - 2 * pow(m_rho, 2) * (10 * s + 9 * t)) +
                2 * pow(pion_mass, 2) * (pow(m_rho, 2) - s - t) *
                    (7 * pow(m_rho, 4) - pow(m_rho, 2) * (8 * s + 7 * t) +
                     2 * (pow(s, 2) + pow(t, 2)))) +
           pow(eta1, 2) *
               (pow(pion_mass, 8) +
                2 * pow(pion_mass, 6) * (pow(m_rho, 2) - 2 * t) +
                (s + t) * (-pow(m_rho, 2) + s + t) *
                    (pow(s, 2) + pow(t, 2) - pow(m_rho, 2) * (s + t)) +
                pow(pion_mass, 4) *
                    (5 * pow(m_rho, 4) + 4 * pow(s, 2) + 2 * s * t +
                     6 * pow(t, 2) - 2 * pow(m_rho, 2) * (5 * s + 3 * t)) -
                2 * pow(pion_mass, 2) *
                    (2 * pow(m_rho, 4) * (s + t) +
                     2 * (s + t) * (pow(s, 2) + pow(t, 2)) -
                     pow(m_rho, 2) *
                         (4 * pow(s, 2) + 5 * s * t + 3 * pow(t, 2)))))) /
             (pow(Gammaa1, 2) * pow(a1_mass, 2) +
              pow(pow(a1_mass, 2) - 2 * pow(pion_mass, 2) - pow(m_rho, 2) + s +
                      t,
                  2)) +
         (0.0625 * pow(eta1 - eta2, 2) *
          (pow(a1_mass, 2) - 2 * pow(pion_mass, 2) - pow(m_rho, 2) + s + t) *
          (-(pow(eta2, 2) *
             (pow(pion_mass, 8) +
              2 * pow(pion_mass, 6) * (pow(m_rho, 2) - 2 * t) +
              2 * pow(pion_mass, 2) * t *
                  (pow(pow(m_rho, 2) - s, 2) + (3 * pow(m_rho, 2) - 2 * s) * t -
                   2 * pow(t, 2)) +
              (pow(m_rho, 2) - s - t) * t *
                  (2 * pow(m_rho, 4) + pow(s, 2) - s * t - pow(t, 2) +
                   pow(m_rho, 2) * (-3 * s + t)) +
              pow(pion_mass, 4) * (pow(m_rho, 4) + 2 * t * (s + 3 * t) -
                                   pow(m_rho, 2) * (s + 6 * t)))) -
           pow(eta1, 2) *
               (pow(pion_mass, 8) +
                2 * pow(pion_mass, 6) * (pow(m_rho, 2) - 2 * t) +
                (pow(m_rho, 2) - s - t) * t *
                    (pow(s, 2) - s * t - pow(t, 2) + pow(m_rho, 2) * (s + t)) +
                pow(pion_mass, 4) * (3 * pow(m_rho, 4) + 2 * t * (s + 3 * t) -
                                     pow(m_rho, 2) * (5 * s + 6 * t)) +
                2 * pow(pion_mass, 2) *
                    (-(pow(m_rho, 4) * (s + t)) +
                     t * (pow(s, 2) - 2 * s * t - 2 * pow(t, 2)) +
                     pow(m_rho, 2) * (pow(s, 2) + 2 * s * t + 3 * pow(t, 2)))) +
           2 * eta1 * eta2 *
               (pow(pion_mass, 8) +
                2 * pow(pion_mass, 6) * (pow(m_rho, 2) - 2 * t) -
                (pow(m_rho, 2) - s - t) * t *
                    (pow(m_rho, 4) - pow(s, 2) - pow(m_rho, 2) * t +
                     t * (s + t)) +
                2 * pow(pion_mass, 4) *
                    (2 * pow(m_rho, 4) + t * (s + 3 * t) -
                     pow(m_rho, 2) * (2 * s + 3 * t)) +
                pow(pion_mass, 2) *
                    (pow(m_rho, 6) - 2 * pow(m_rho, 4) * (s + 2 * t) +
                     2 * t * (pow(s, 2) - 2 * s * t - 2 * pow(t, 2)) +
                     pow(m_rho, 2) *
                         (pow(s, 2) + 2 * s * t + 6 * pow(t, 2)))))) /
             ((-pow(a1_mass, 2) + t) *
              (pow(Gammaa1, 2) * pow(a1_mass, 2) +
               pow(pow(a1_mass, 2) - 2 * pow(pion_mass, 2) - pow(m_rho, 2) + s +
                       t,
                   2))) +
         (0.125 * (eta1 - eta2) *
          (eta2 * (8 * C4 * pow(m_rho, 6) * t - 2 * delta * pow(s, 2) * t +
                   pow(m_rho, 2) *
                       (-4 * pow(pion_mass, 4) +
                        (s * (2 + 3 * delta + 8 * C4 * s) - 4 * t) * t +
                        pow(pion_mass, 2) * (-((-2 + delta) * s) + 8 * t)) +
                   pow(m_rho, 4) *
                       (8 * C4 * pow(pion_mass, 4) -
                        pow(pion_mass, 2) * (-2 + delta + 16 * C4 * t) +
                        t * (-6 + delta + 8 * C4 * (-2 * s + t)))) +
           eta1 * (2 * delta * pow(s, 2) * t +
                   8 * C4 * pow(m_rho, 6) * (-2 * pow(pion_mass, 2) + t) -
                   pow(m_rho, 2) *
                       (-4 * pow(pion_mass, 4) - 4 * pow(t, 2) +
                        2 * pow(pion_mass, 2) * ((2 + delta) * s + 4 * t) +
                        pow(s, 2) * (-2 + delta + 8 * C4 * t)) +
                   pow(m_rho, 4) *
                       (-8 * C4 * pow(pion_mass, 4) + (-2 + delta) * s -
                        4 * t * (1 + 2 * C4 * t) +
                        8 * pow(pion_mass, 2) * (1 + 2 * C4 * (s + t)))))) /
             (pow(m_rho, 2) * (-pow(a1_mass, 2) + t)) -
         (0.125 * (eta1 - eta2) *
          (pow(a1_mass, 2) - 2 * pow(pion_mass, 2) - pow(m_rho, 2) + s + t) *
          (eta1 *
               (pow(pion_mass, 4) *
                    (4 * pow(m_rho, 2) - 8 * C4 * pow(m_rho, 4)) +
                8 * C4 * pow(m_rho, 6) * (s + t) -
                2 * delta * pow(s, 2) * (s + t) +
                pow(m_rho, 2) * ((6 + delta) * pow(s, 2) + 8 * s * t +
                                 4 * pow(t, 2) + 8 * C4 * pow(s, 2) * (s + t)) -
                pow(m_rho, 4) *
                    (-((-6 + delta) * s) + 4 * t +
                     8 * C4 * (2 * pow(s, 2) + 2 * s * t + pow(t, 2))) +
                2 * pow(pion_mass, 2) *
                    (-8 * C4 * pow(m_rho, 6) + 2 * delta * pow(s, 2) -
                     pow(m_rho, 2) * (s * (6 + delta + 8 * C4 * s) + 4 * t) +
                     4 * pow(m_rho, 4) * (1 + 2 * C4 * (2 * s + t)))) +
           eta2 * (pow(pion_mass, 4) *
                       (-4 * pow(m_rho, 2) + 8 * C4 * pow(m_rho, 4)) -
                   (-pow(m_rho, 2) + s + t) *
                       (16 * C4 * pow(m_rho, 6) - 2 * delta * pow(s, 2) +
                        pow(m_rho, 2) *
                            (s * (6 + 3 * delta + 8 * C4 * s) + 4 * t) +
                        pow(m_rho, 4) * (-10 + delta - 8 * C4 * (3 * s + t))) +
                   pow(pion_mass, 2) *
                       (32 * C4 * pow(m_rho, 6) - 4 * delta * pow(s, 2) +
                        pow(m_rho, 2) *
                            (s * (14 + 5 * delta + 16 * C4 * s) + 8 * t) +
                        pow(m_rho, 4) *
                            (delta - 2 * (9 + 8 * C4 * (3 * s + t))))))) /
             (pow(m_rho, 2) * (pow(Gammaa1, 2) * pow(a1_mass, 2) +
                               pow(pow(a1_mass, 2) - 2 * pow(pion_mass, 2) -
                                       pow(m_rho, 2) + s + t,
                                   2))))) /
       (16. * Pi * s * (-4 * pow(pion_mass, 2) + s)));

  // clang-format on
  return cut_off(gev2_mb * diff_xs / spin_deg_factor);
}

// C22
double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_pi_pi0_rho(
    const double s, const double m_rho) {
  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  const double spin_deg_factor = 1.0;
  auto t_mandelstam = get_t_range(sqrt(s), pion_mass, pion_mass, m_rho, 0.0);

  const double tmin = t_mandelstam[1];
  const double tmax = t_mandelstam[0];

  // clang-format off
  const double
      xs =
          (pow(Const, 2) * pow(ghat, 4) *
           ((-0.03125 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-2. * pow(a1_mass, 8) - 2. * pow(pion_mass, 8) +
                   2. * pow(pion_mass, 4) * pow(m_rho, 4) +
                   pow(a1_mass, 6) * (8. * pow(pion_mass, 2) - 4. * s) +
                   pow(a1_mass, 2) * pow(pion_mass, 2) *
                       (8. * pow(pion_mass, 4) - 8. * pow(m_rho, 4) -
                        4. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s) +
                   pow(a1_mass, 4) *
                       (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                        8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                        4. * pow(s, 2))) +
              pow(eta2, 2) *
                  (1. * pow(a1_mass, 8) + 1. * pow(pion_mass, 8) -
                   2. * pow(pion_mass, 6) * pow(m_rho, 2) +
                   1. * pow(pion_mass, 4) * pow(m_rho, 4) +
                   pow(a1_mass, 6) *
                       (-4. * pow(pion_mass, 2) + 2. * pow(m_rho, 2) + 2. * s) +
                   pow(a1_mass, 4) *
                       (6. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                        pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                        2. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                   pow(a1_mass, 2) *
                       (-4. * pow(pion_mass, 6) -
                        2. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                        pow(pion_mass, 4) * (6. * pow(m_rho, 2) + 2. * s))) +
              pow(eta1, 2) *
                  (1. * pow(a1_mass, 8) +
                   pow(a1_mass, 6) *
                       (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
                   pow(a1_mass, 4) *
                       (6. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                        pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                        4. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                   pow(a1_mass, 2) *
                       (-4. * pow(pion_mass, 6) +
                        2. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                        pow(m_rho, 2) * (2. * pow(m_rho, 2) - 2. * s) * s +
                        pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 2. * s)) +
                   pow(pion_mass, 2) *
                       (1. * pow(pion_mass, 6) +
                        2. * pow(pion_mass, 4) * pow(m_rho, 2) -
                        2. * pow(m_rho, 6) + 2. * pow(m_rho, 4) * s +
                        pow(pion_mass, 2) *
                            (1. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s))))) /
                (1. * pow(a1_mass, 2) - 1. * tmin) -
            (1. * pow(-2. + delta, 2) * pow(pion_mass, 2) *
             (1. * pow(pion_mass, 2) - 0.25 * pow(m_rho, 2))) /
                (1. * pow(pion_mass, 2) - 1. * tmin) -
            0.75 * tmin +
            (0.25 * pow(-2. + delta, 2) * pow(pion_mass, 2) * tmin) /
                pow(m_rho, 2) +
            0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
                (eta2 * (-1. * pow(a1_mass, 2) + pow(m_rho, 2) - 2. * s) +
                 eta1 * (2. * pow(pion_mass, 2) + s)) *
                tmin -
            C4 *
                (2. * C4 * pow(m_rho, 4) + pow(m_rho, 2) * (-3. - 4. * C4 * s) +
                 s * (3. + 2. * C4 * s)) *
                tmin -
            (0.5 * pow(delta, 2) *
             (1. * pow(pion_mass, 4) * pow(m_rho, 2) + 0.25 * pow(m_rho, 6) -
              0.75 * pow(m_rho, 4) * s + 0.125 * pow(m_rho, 2) * pow(s, 2) +
              0.25 * pow(s, 3) +
              pow(pion_mass, 2) *
                  (2.5 * pow(m_rho, 4) + 0.25 * pow(m_rho, 2) * s -
                   0.75 * pow(s, 2))) *
             tmin) /
                pow(m_rho, 6) -
            0.03125 * pow(eta1 - 1. * eta2, 2) *
                (eta1 * eta2 *
                     (-6. * pow(a1_mass, 4) - 12. * pow(pion_mass, 4) +
                      2. * pow(m_rho, 4) +
                      pow(a1_mass, 2) * (16. * pow(pion_mass, 2) - 8. * s) +
                      8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                      4. * pow(s, 2)) +
                 pow(eta1, 2) *
                     (3. * pow(a1_mass, 4) + 6. * pow(pion_mass, 4) +
                      pow(m_rho, 4) +
                      pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                      4. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                      pow(a1_mass, 2) * (-8. * pow(pion_mass, 2) -
                                         4. * pow(m_rho, 2) + 4. * s)) +
                 pow(eta2, 2) *
                     (3. * pow(a1_mass, 4) + 6. * pow(pion_mass, 4) +
                      pow(m_rho, 4) +
                      pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                      2. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                      pow(a1_mass, 2) * (-8. * pow(pion_mass, 2) +
                                         4. * pow(m_rho, 2) + 4. * s))) *
                tmin -
            (3. * delta *
             (0.6666666666666666 * C4 * pow(m_rho, 6) -
              0.08333333333333333 * pow(s, 2) +
              pow(m_rho, 4) * (-0.25 - 0.5 * C4 * s) +
              pow(m_rho, 2) *
                  s * (0.08333333333333333 - 0.16666666666666666 * C4 * s) +
              pow(pion_mass, 2) *
                  (1. * C4 * pow(m_rho, 4) + 0.08333333333333333 * s +
                   pow(m_rho, 2) *
                       (-0.4166666666666667 - 0.3333333333333333 * C4 * s))) *
             tmin) /
                pow(m_rho, 4) -
            (0.25 * (1. * eta1 - 1. * eta2) *
             (pow(m_rho, 2) *
                  (eta1 * (-1. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4) +
                           pow(a1_mass, 2) * (1. - 2. * C4 * pow(m_rho, 2)) +
                           pow(pion_mass, 2) * (-2. + 4. * C4 * pow(m_rho, 2)) -
                           2. * C4 * pow(s, 2)) +
                   eta2 * (-1.5 * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4) +
                           pow(pion_mass, 2) * (2. - 4. * C4 * pow(m_rho, 2)) +
                           pow(a1_mass, 2) * (-1. + 2. * C4 * pow(m_rho, 2)) +
                           0.5 * s - 4. * C4 * pow(m_rho, 2) * s +
                           2. * C4 * pow(s, 2))) +
              delta * (eta2 *
                           (-1. * pow(a1_mass, 4) - 3. * pow(pion_mass, 4) +
                            1. * pow(m_rho, 4) +
                            pow(a1_mass, 2) * (3. * pow(pion_mass, 2) -
                                               1. * pow(m_rho, 2) - 1. * s) +
                            0.25 * pow(m_rho, 2) * s - 0.75 * pow(s, 2) +
                            pow(pion_mass, 2) * (1. * pow(m_rho, 2) + 1. * s)) +
                       eta1 *
                           (1. * pow(a1_mass, 4) + 3. * pow(pion_mass, 4) +
                            0.5 * pow(m_rho, 4) +
                            pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 2. * s) -
                            2. * pow(m_rho, 2) * s + 1. * pow(s, 2) +
                            pow(a1_mass, 2) *
                                (-3. * pow(pion_mass, 2) - 1.5 * pow(m_rho, 2) +
                                 1.5 * s)))) *
             tmin) /
                pow(m_rho, 2) -
            0.0625 * (-2. + delta) * (eta1 - 1. * eta2) * eta2 * pow(tmin, 2) +
            (0.25 * pow(delta, 2) *
             (2. * pow(pion_mass, 2) * pow(m_rho, 2) + 1. * pow(m_rho, 4) -
              0.75 * pow(m_rho, 2) * s - 0.25 * pow(s, 2)) *
             pow(tmin, 2)) /
                pow(m_rho, 6) -
            0.03125 * pow(eta1 - 1. * eta2, 3) *
                (eta2 * (-1. * pow(a1_mass, 2) + 2. * pow(pion_mass, 2) -
                         1. * pow(m_rho, 2) - 1. * s) +
                 eta1 * (pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                         1. * pow(m_rho, 2) + s)) *
                pow(tmin, 2) +
            (1.5 * delta *
             (1. * C4 * pow(m_rho, 4) + 0.08333333333333333 * s +
              pow(m_rho, 2) *
                  (-0.4166666666666667 - 0.3333333333333333 * C4 * s)) *
             pow(tmin, 2)) /
                pow(m_rho, 4) -
            (0.125 * (1. * eta1 - 1. * eta2) *
             (pow(m_rho, 2) * (eta1 * (1. - 2. * C4 * pow(m_rho, 2)) +
                               eta2 * (-1. + 2. * C4 * pow(m_rho, 2))) +
              delta * (eta2 * (-1. * pow(a1_mass, 2) + 3. * pow(pion_mass, 2) -
                               1. * pow(m_rho, 2) - 1. * s) +
                       eta1 * (1. * pow(a1_mass, 2) - 3. * pow(pion_mass, 2) -
                               1.5 * pow(m_rho, 2) + 1.5 * s))) *
             pow(tmin, 2)) /
                pow(m_rho, 2) -
            0.010416666666666666 * pow(eta1 - 1. * eta2, 4) * pow(tmin, 3) -
            (0.16666666666666666 * pow(delta, 2) * pow(tmin, 3)) /
                pow(m_rho, 4) -
            (0.08333333333333333 * delta * pow(1. * eta1 - 1. * eta2, 2) *
             pow(tmin, 3)) /
                pow(m_rho, 2) +
            (0.03125 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-2. * pow(a1_mass, 8) - 2. * pow(pion_mass, 8) +
                   2. * pow(pion_mass, 4) * pow(m_rho, 4) +
                   pow(a1_mass, 6) * (8. * pow(pion_mass, 2) - 4. * s) +
                   pow(a1_mass, 2) * pow(pion_mass, 2) *
                       (8. * pow(pion_mass, 4) - 8. * pow(m_rho, 4) -
                        4. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s) +
                   pow(a1_mass, 4) *
                       (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                        8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                        4. * pow(s, 2))) +
              pow(eta2, 2) *
                  (1. * pow(a1_mass, 8) + 1. * pow(pion_mass, 8) -
                   2. * pow(pion_mass, 6) * pow(m_rho, 2) +
                   1. * pow(pion_mass, 4) * pow(m_rho, 4) +
                   pow(a1_mass, 6) *
                       (-4. * pow(pion_mass, 2) + 2. * pow(m_rho, 2) + 2. * s) +
                   pow(a1_mass, 4) *
                       (6. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                        pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                        2. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                   pow(a1_mass, 2) *
                       (-4. * pow(pion_mass, 6) -
                        2. * pow(pion_mass, 2) *
                            pow(m_rho, 2) * s +
                        pow(pion_mass, 4) * (6. * pow(m_rho, 2) + 2. * s))) +
              pow(eta1, 2) *
                  (1. * pow(a1_mass, 8) +
                   pow(a1_mass, 6) *
                       (-4. * pow(pion_mass, 2) - 2. * pow(m_rho, 2) + 2. * s) +
                   pow(a1_mass, 4) *
                       (6. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                        pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                        4. * pow(m_rho, 2) * s + 2. * pow(s, 2)) +
                   pow(a1_mass, 2) *
                       (-4. * pow(pion_mass, 6) +
                        2. * pow(pion_mass, 2) *
                            pow(m_rho, 2) * s +
                        pow(m_rho, 2) * (2. * pow(m_rho, 2) - 2. * s) * s +
                        pow(pion_mass, 4) * (-6. * pow(m_rho, 2) + 2. * s)) +
                   pow(pion_mass, 2) *
                       (1. * pow(pion_mass, 6) +
                        2. * pow(pion_mass, 4) * pow(m_rho, 2) -
                        2. * pow(m_rho, 6) + 2. * pow(m_rho, 4) * s +
                        pow(pion_mass, 2) *
                            (1. * pow(m_rho, 4) - 2. * pow(m_rho, 2) * s))))) /
                (1. * pow(a1_mass, 2) - 1. * tmax) +
            (1. * pow(-2. + delta, 2) * pow(pion_mass, 2) *
             (1. * pow(pion_mass, 2) - 0.25 * pow(m_rho, 2))) /
                (1. * pow(pion_mass, 2) - 1. * tmax) +
            0.75 * tmax -
            (0.25 * pow(-2. + delta, 2) * pow(pion_mass, 2) * tmax) /
                pow(m_rho, 2) -
            0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
                (eta2 * (-1. * pow(a1_mass, 2) + pow(m_rho, 2) - 2. * s) +
                 eta1 * (2. * pow(pion_mass, 2) + s)) *
                tmax +
            C4 *
                (2. * C4 * pow(m_rho, 4) + pow(m_rho, 2) * (-3. - 4. * C4 * s) +
                 s * (3. + 2. * C4 * s)) *
                tmax +
            (0.5 * pow(delta, 2) *
             (1. * pow(pion_mass, 4) * pow(m_rho, 2) + 0.25 * pow(m_rho, 6) -
              0.75 * pow(m_rho, 4) * s + 0.125 * pow(m_rho, 2) * pow(s, 2) +
              0.25 * pow(s, 3) +
              pow(pion_mass, 2) *
                  (2.5 * pow(m_rho, 4) + 0.25 * pow(m_rho, 2) * s -
                   0.75 * pow(s, 2))) *
             tmax) /
                pow(m_rho, 6) +
            0.03125 * pow(eta1 - 1. * eta2, 2) *
                (eta1 * eta2 *
                     (-6. * pow(a1_mass, 4) - 12. * pow(pion_mass, 4) +
                      2. * pow(m_rho, 4) +
                      pow(a1_mass, 2) * (16. * pow(pion_mass, 2) - 8. * s) +
                      8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                      4. * pow(s, 2)) +
                 pow(eta1, 2) *
                     (3. * pow(a1_mass, 4) + 6. * pow(pion_mass, 4) +
                      pow(m_rho, 4) +
                      pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                      4. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                      pow(a1_mass, 2) * (-8. * pow(pion_mass, 2) -
                                         4. * pow(m_rho, 2) + 4. * s)) +
                 pow(eta2, 2) *
                     (3. * pow(a1_mass, 4) + 6. * pow(pion_mass, 4) +
                      pow(m_rho, 4) +
                      pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                      2. * pow(m_rho, 2) * s + 2. * pow(s, 2) +
                      pow(a1_mass, 2) * (-8. * pow(pion_mass, 2) +
                                         4. * pow(m_rho, 2) + 4. * s))) *
                tmax +
            (3. * delta *
             (0.6666666666666666 * C4 * pow(m_rho, 6) -
              0.08333333333333333 * pow(s, 2) +
              pow(m_rho, 4) * (-0.25 - 0.5 * C4 * s) +
              pow(m_rho, 2) *
                  s * (0.08333333333333333 - 0.16666666666666666 * C4 * s) +
              pow(pion_mass, 2) *
                  (1. * C4 * pow(m_rho, 4) + 0.08333333333333333 * s +
                   pow(m_rho, 2) *
                       (-0.4166666666666667 - 0.3333333333333333 * C4 * s))) *
             tmax) /
                pow(m_rho, 4) +
            (0.25 * (1. * eta1 - 1. * eta2) *
             (pow(m_rho, 2) *
                  (eta1 * (-1. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4) +
                           pow(a1_mass, 2) * (1. - 2. * C4 * pow(m_rho, 2)) +
                           pow(pion_mass, 2) * (-2. + 4. * C4 * pow(m_rho, 2)) -
                           2. * C4 * pow(s, 2)) +
                   eta2 * (-1.5 * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4) +
                           pow(pion_mass, 2) * (2. - 4. * C4 * pow(m_rho, 2)) +
                           pow(a1_mass, 2) * (-1. + 2. * C4 * pow(m_rho, 2)) +
                           0.5 * s - 4. * C4 * pow(m_rho, 2) * s +
                           2. * C4 * pow(s, 2))) +
              delta * (eta2 *
                           (-1. * pow(a1_mass, 4) - 3. * pow(pion_mass, 4) +
                            1. * pow(m_rho, 4) +
                            pow(a1_mass, 2) * (3. * pow(pion_mass, 2) -
                                               1. * pow(m_rho, 2) - 1. * s) +
                            0.25 * pow(m_rho, 2) * s - 0.75 * pow(s, 2) +
                            pow(pion_mass, 2) * (1. * pow(m_rho, 2) + 1. * s)) +
                       eta1 *
                           (1. * pow(a1_mass, 4) + 3. * pow(pion_mass, 4) +
                            0.5 * pow(m_rho, 4) +
                            pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 2. * s) -
                            2. * pow(m_rho, 2) * s + 1. * pow(s, 2) +
                            pow(a1_mass, 2) *
                                (-3. * pow(pion_mass, 2) - 1.5 * pow(m_rho, 2) +
                                 1.5 * s)))) *
             tmax) /
                pow(m_rho, 2) +
            0.0625 * (-2. + delta) * (eta1 - 1. * eta2) * eta2 * pow(tmax, 2) -
            (0.25 * pow(delta, 2) *
             (2. * pow(pion_mass, 2) * pow(m_rho, 2) + 1. * pow(m_rho, 4) -
              0.75 * pow(m_rho, 2) * s - 0.25 * pow(s, 2)) *
             pow(tmax, 2)) /
                pow(m_rho, 6) +
            0.03125 * pow(eta1 - 1. * eta2, 3) *
                (eta2 * (-1. * pow(a1_mass, 2) + 2. * pow(pion_mass, 2) -
                         1. * pow(m_rho, 2) - 1. * s) +
                 eta1 * (pow(a1_mass, 2) - 2. * pow(pion_mass, 2) -
                         1. * pow(m_rho, 2) + s)) *
                pow(tmax, 2) -
            (1.5 * delta *
             (1. * C4 * pow(m_rho, 4) + 0.08333333333333333 * s +
              pow(m_rho, 2) *
                  (-0.4166666666666667 - 0.3333333333333333 * C4 * s)) *
             pow(tmax, 2)) /
                pow(m_rho, 4) +
            (0.125 * (1. * eta1 - 1. * eta2) *
             (pow(m_rho, 2) * (eta1 * (1. - 2. * C4 * pow(m_rho, 2)) +
                               eta2 * (-1. + 2. * C4 * pow(m_rho, 2))) +
              delta * (eta2 * (-1. * pow(a1_mass, 2) + 3. * pow(pion_mass, 2) -
                               1. * pow(m_rho, 2) - 1. * s) +
                       eta1 * (1. * pow(a1_mass, 2) - 3. * pow(pion_mass, 2) -
                               1.5 * pow(m_rho, 2) + 1.5 * s))) *
             pow(tmax, 2)) /
                pow(m_rho, 2) +
            0.010416666666666666 * pow(eta1 - 1. * eta2, 4) * pow(tmax, 3) +
            (0.16666666666666666 * pow(delta, 2) * pow(tmax, 3)) /
                pow(m_rho, 4) +
            (0.08333333333333333 * delta * pow(1. * eta1 - 1. * eta2, 2) *
             pow(tmax, 3)) /
                pow(m_rho, 2) -
            (2. * pow(pion_mass, 4) * tmin * HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (5. * pow(pion_mass, 2) * pow(m_rho, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (0.5 * pow(m_rho, 4) * tmin * HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (2.499999 * pow(pion_mass, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (5.0000100000000005 * delta * pow(pion_mass, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (1.2500010000000001 * pow(m_rho, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (0.500001 * delta * pow(m_rho, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (6. * C4 * pow(pion_mass, 2) * pow(m_rho, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (3. * C4 * pow(m_rho, 4) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (5. * delta * pow(pion_mass, 2) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (1.75 * pow(m_rho, 2) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (0.5 * delta * pow(m_rho, 2) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (1.5 * s * tmin * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (1.0000005 * delta * s * tmin * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (0.2500005 * pow(delta, 2) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (2.000001 * C4 * pow(pion_mass, 2) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (3. * C4 * delta * pow(pion_mass, 2) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (3.9999899999999995 * C4 * pow(m_rho, 2) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (1.5 * C4 * delta * pow(m_rho, 2) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (1.75 * delta * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (0.125 * pow(delta, 2) * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (0.5 * pow(delta, 2) * pow(pion_mass, 4) * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) +
            (0.999999 * C4 * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (1.9999949999999997 * C4 * delta * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (1. * delta * pow(pion_mass, 2) * pow(s, 3) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) +
            (0.25 * pow(delta, 2) * pow(pion_mass, 2) * pow(s, 4) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) * pow(pow(m_rho, 2) - 1. * s, 2)) +
            (0.25 * delta * pow(s, 4) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) -
            (0.0625 * pow(delta, 2) * pow(s, 5) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) * pow(pow(m_rho, 2) - 1. * s, 2)) +
            (2. * delta * pow(pion_mass, 4) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) +
            (1. * pow(pion_mass, 2) * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) -
            (1.25 * pow(delta, 2) * pow(pion_mass, 2) * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) -
            (0.25 * pow(s, 3) * tmin * HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) +
            (0.4375 * pow(delta, 2) * pow(s, 3) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) +
            (2.000001 * delta * pow(pion_mass, 4) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (0.500001 * pow(pion_mass, 2) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (1.4999993999999999 * delta * pow(pion_mass, 2) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (2.5000050000000003 * pow(delta, 2) * pow(pion_mass, 2) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (0.2499999 * pow(s, 2) * tmin * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (0.9999998999999999 * delta * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (0.8125005000000001 * pow(delta, 2) * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (1.0000005 * C4 * delta * pow(pion_mass, 2) * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (0.4999995 * C4 * delta * pow(s, 3) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (1.0000005 * pow(delta, 2) * pow(pion_mass, 4) * s * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) -
            (1.5000015000000002 * delta * pow(pion_mass, 2) * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) -
            (0.12499995 * pow(delta, 2) * pow(pion_mass, 2) * pow(s, 2) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) +
            (0.49999994999999997 * delta * pow(s, 3) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) +
            (0.12499995 * pow(delta, 2) * pow(s, 3) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) +
            (0.6250005000000001 * pow(delta, 2) * pow(pion_mass, 2) *
             pow(s, 3) * tmin * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 8) - 1. * pow(m_rho, 6) * s) -
            (0.1875 * pow(delta, 2) * pow(s, 4) * tmin *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 8) - 1. * pow(m_rho, 6) * s) +
            (2. * pow(pion_mass, 2) * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (1. * pow(m_rho, 2) * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (1.2499995 * pow(tmin, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (1.0000005 * delta * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (3. * C4 * pow(m_rho, 2) * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (1. * s * pow(tmin, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (1. * delta * s * pow(tmin, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (1.0000005 * C4 * s * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (1.5 * C4 * delta * s * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (0.5 * pow(delta, 2) * pow(pion_mass, 2) * pow(s, 2) *
             pow(tmin, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) -
            (0.25 * pow(delta, 2) * pow(s, 3) * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) -
            (2. * delta * pow(pion_mass, 2) * s * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) +
            (1. * delta * pow(s, 2) * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) +
            (0.25 * pow(delta, 2) * pow(s, 2) * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) -
            (1.9999949999999997 * delta * pow(pion_mass, 2) * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (0.2500005 * s * pow(tmin, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (0.24999974999999997 * delta * s * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (0.50000025 * pow(delta, 2) * s * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (0.50000025 * C4 * delta * pow(s, 2) * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (0.9999974999999999 * pow(delta, 2) * pow(pion_mass, 2) *
             s * pow(tmin, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) +
            (0.2500002 * delta * pow(s, 2) * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) -
            (0.43749974999999997 * pow(delta, 2) * pow(s, 2) * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) -
            (0.062499975 * pow(delta, 2) * pow(s, 3) * pow(tmin, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 8) - 1. * pow(m_rho, 6) * s) -
            (0.6666666666666666 * pow(tmin, 3) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (0.16666666666666666 * pow(delta, 2) * pow(s, 2) * pow(tmin, 3) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) +
            (0.6666666666666666 * delta * s * pow(tmin, 3) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) +
            (0.666667 * delta * pow(tmin, 3) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (0.3333335 * pow(delta, 2) * s * pow(tmin, 3) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) +
            (2. * pow(pion_mass, 4) * tmax * HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (5. * pow(pion_mass, 2) * pow(m_rho, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (0.5 * pow(m_rho, 4) * tmax * HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (2.499999 * pow(pion_mass, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (5.0000100000000005 * delta * pow(pion_mass, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (1.2500010000000001 * pow(m_rho, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (0.500001 * delta * pow(m_rho, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (6. * C4 * pow(pion_mass, 2) * pow(m_rho, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (3. * C4 * pow(m_rho, 4) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (5. * delta * pow(pion_mass, 2) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (1.75 * pow(m_rho, 2) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (0.5 * delta * pow(m_rho, 2) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (1.5 * s * tmax * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (1.0000005 * delta * s * tmax * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (0.2500005 * pow(delta, 2) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (2.000001 * C4 * pow(pion_mass, 2) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (3. * C4 * delta * pow(pion_mass, 2) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (3.9999899999999995 * C4 * pow(m_rho, 2) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (1.5 * C4 * delta * pow(m_rho, 2) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (1.75 * delta * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (0.125 * pow(delta, 2) * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (0.5 * pow(delta, 2) * pow(pion_mass, 4) * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) -
            (0.999999 * C4 * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (1.9999949999999997 * C4 * delta * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (1. * delta * pow(pion_mass, 2) * pow(s, 3) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) -
            (0.25 * pow(delta, 2) * pow(pion_mass, 2) * pow(s, 4) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) * pow(pow(m_rho, 2) - 1. * s, 2)) -
            (0.25 * delta * pow(s, 4) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) +
            (0.0625 * pow(delta, 2) * pow(s, 5) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) * pow(pow(m_rho, 2) - 1. * s, 2)) -
            (2. * delta * pow(pion_mass, 4) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) -
            (1. * pow(pion_mass, 2) * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) +
            (1.25 * pow(delta, 2) * pow(pion_mass, 2) * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) +
            (0.25 * pow(s, 3) * tmax * HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) -
            (0.4375 * pow(delta, 2) * pow(s, 3) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) -
            (2.000001 * delta * pow(pion_mass, 4) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (0.500001 * pow(pion_mass, 2) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (1.4999993999999999 * delta * pow(pion_mass, 2) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (2.5000050000000003 * pow(delta, 2) * pow(pion_mass, 2) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (0.2499999 * pow(s, 2) * tmax * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (0.9999998999999999 * delta * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (0.8125005000000001 * pow(delta, 2) * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (1.0000005 * C4 * delta * pow(pion_mass, 2) * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (0.4999995 * C4 * delta * pow(s, 3) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (1.0000005 * pow(delta, 2) * pow(pion_mass, 4) * s * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) +
            (1.5000015000000002 * delta * pow(pion_mass, 2) * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) +
            (0.12499995 * pow(delta, 2) * pow(pion_mass, 2) * pow(s, 2) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) -
            (0.49999994999999997 * delta * pow(s, 3) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) -
            (0.12499995 * pow(delta, 2) * pow(s, 3) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) -
            (0.6250005000000001 * pow(delta, 2) * pow(pion_mass, 2) *
             pow(s, 3) * tmax * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 8) - 1. * pow(m_rho, 6) * s) +
            (0.1875 * pow(delta, 2) * pow(s, 4) * tmax *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 8) - 1. * pow(m_rho, 6) * s) -
            (2. * pow(pion_mass, 2) * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (1. * pow(m_rho, 2) * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (1.2499995 * pow(tmax, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (1.0000005 * delta * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (3. * C4 * pow(m_rho, 2) * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) +
            (1. * s * pow(tmax, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (1. * delta * s * pow(tmax, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) -
            (1.0000005 * C4 * s * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (1.5 * C4 * delta * s * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 2) - 1. * s) -
            (0.5 * pow(delta, 2) * pow(pion_mass, 2) * pow(s, 2) *
             pow(tmax, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) +
            (0.25 * pow(delta, 2) * pow(s, 3) * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) +
            (2. * delta * pow(pion_mass, 2) * s * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) -
            (1. * delta * pow(s, 2) * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) -
            (0.25 * pow(delta, 2) * pow(s, 2) * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) +
            (1.9999949999999997 * delta * pow(pion_mass, 2) * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (0.2500005 * s * pow(tmax, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (0.24999974999999997 * delta * s * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (0.50000025 * pow(delta, 2) * s * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (0.50000025 * C4 * delta * pow(s, 2) * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) -
            (0.9999974999999999 * pow(delta, 2) * pow(pion_mass, 2) *
             s * pow(tmax, 2) * HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) -
            (0.2500002 * delta * pow(s, 2) * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) +
            (0.43749974999999997 * pow(delta, 2) * pow(s, 2) * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) +
            (0.062499975 * pow(delta, 2) * pow(s, 3) * pow(tmax, 2) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 8) - 1. * pow(m_rho, 6) * s) +
            (0.6666666666666666 * pow(tmax, 3) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 2) - 1. * s, 2) +
            (0.16666666666666666 * pow(delta, 2) * pow(s, 2) * pow(tmax, 3) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) * pow(pow(m_rho, 2) - 1. * s, 2)) -
            (0.6666666666666666 * delta * s * pow(tmax, 3) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                pow(pow(m_rho, 3) - 1. * m_rho * s, 2) -
            (0.666667 * delta * pow(tmax, 3) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 4) - 1. * pow(m_rho, 2) * s) +
            (0.3333335 * pow(delta, 2) * s * pow(tmax, 3) *
             HeavisideTheta(-m_rho + sqrt(s))) /
                (pow(m_rho, 6) - 1. * pow(m_rho, 4) * s) -
            0.0625 * pow(eta1 - 1. * eta2, 2) *
                (eta1 * eta2 *
                     (-4. * pow(a1_mass, 6) +
                      pow(a1_mass, 4) * (12. * pow(pion_mass, 2) - 6. * s) +
                      pow(pion_mass, 2) *
                          (4. * pow(pion_mass, 4) - 4. * pow(m_rho, 4) -
                           2. * pow(pion_mass, 2) * s +
                           2. * pow(m_rho, 2) * s) +
                      pow(a1_mass, 2) *
                          (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                           8. * pow(pion_mass, 2) * s + 4. * pow(m_rho, 2) * s -
                           4. * pow(s, 2))) +
                 pow(eta1, 2) *
                     (2. * pow(a1_mass, 6) - 2. * pow(pion_mass, 6) +
                      pow(pion_mass, 2) * pow(m_rho, 2) * s +
                      pow(m_rho, 2) * (pow(m_rho, 2) - 1. * s) * s +
                      pow(pion_mass, 4) * (-3. * pow(m_rho, 2) + s) +
                      pow(a1_mass, 4) * (-6. * pow(pion_mass, 2) -
                                         3. * pow(m_rho, 2) + 3. * s) +
                      pow(a1_mass, 2) *
                          (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                           pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                           4. * pow(m_rho, 2) * s + 2. * pow(s, 2))) +
                 pow(eta2, 2) *
                     (2. * pow(a1_mass, 6) - 2. * pow(pion_mass, 6) -
                      1. * pow(pion_mass, 2) * pow(m_rho, 2) * s +
                      pow(pion_mass, 4) * (3. * pow(m_rho, 2) + s) +
                      pow(a1_mass, 4) * (-6. * pow(pion_mass, 2) +
                                         3. * pow(m_rho, 2) + 3. * s) +
                      pow(a1_mass, 2) *
                          (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                           pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                           2. * pow(m_rho, 2) * s + 2. * pow(s, 2)))) *
                log(fabs(-1. * pow(a1_mass, 2) + tmin)) +
            (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta2 * (-0.5 * pow(a1_mass, 6) - 0.5 * pow(pion_mass, 6) +
                      0.5 * pow(pion_mass, 4) * pow(m_rho, 2) +
                      pow(a1_mass, 4) * (0.5 * pow(pion_mass, 2) +
                                         0.5 * pow(m_rho, 2) - 1. * s) +
                      pow(a1_mass, 2) * pow(pion_mass, 2) *
                          (0.5 * pow(pion_mass, 2) + 1. * pow(m_rho, 2) -
                           1. * s)) +
              eta1 *
                  (pow(a1_mass, 4) * (1. * pow(pion_mass, 2) + 0.5 * s) +
                   pow(pion_mass, 2) *
                       (1. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                        pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 0.5 * s) -
                        0.5 * pow(m_rho, 2) * s) +
                   pow(a1_mass, 2) *
                       (-2. * pow(pion_mass, 4) - 0.5 * pow(m_rho, 2) * s +
                        pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s)))) *
             log(fabs(-1. * pow(a1_mass, 2) + tmin))) /
                (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) -
            (0.25 * (1. * eta1 - 1. * eta2) *
             (delta * (eta1 *
                           (1. * pow(a1_mass, 6) - 1. * pow(pion_mass, 6) +
                            pow(pion_mass, 4) *
                                (-2.5 * pow(m_rho, 2) + 0.5 * s) +
                            pow(pion_mass, 2) *
                                s * (-0.5 * pow(m_rho, 2) + 1. * s) +
                            pow(a1_mass, 4) * (-3. * pow(pion_mass, 2) -
                                               1.5 * pow(m_rho, 2) + 1.5 * s) +
                            s * (0.5 * pow(m_rho, 4) -
                                 0.25 * pow(m_rho, 2) * s - 0.25 * pow(s, 2)) +
                            pow(a1_mass, 2) *
                                (3. * pow(pion_mass, 4) + 0.5 * pow(m_rho, 4) +
                                 pow(pion_mass, 2) *
                                     (4. * pow(m_rho, 2) - 2. * s) -
                                 2. * pow(m_rho, 2) * s + 1. * pow(s, 2))) +
                       eta2 *
                           (-1. * pow(a1_mass, 6) +
                            pow(a1_mass, 4) * (3. * pow(pion_mass, 2) -
                                               1. * pow(m_rho, 2) - 1. * s) +
                            pow(pion_mass, 2) *
                                (1. * pow(pion_mass, 4) - 0.5 * pow(m_rho, 4) +
                                 0.25 * pow(m_rho, 2) * s - 0.25 * pow(s, 2)) +
                            pow(a1_mass, 2) *
                                (-3. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                                 0.25 * pow(m_rho, 2) * s - 0.75 * pow(s, 2) +
                                 pow(pion_mass, 2) *
                                     (1. * pow(m_rho, 2) + 1. * s)))) +
              pow(m_rho, 2) *
                  (eta2 * (pow(a1_mass, 4) * (-1. + 2. * C4 * pow(m_rho, 2)) +
                           pow(pion_mass, 2) *
                               (0.5 * pow(m_rho, 2) +
                                pow(pion_mass, 2) *
                                    (-1. + 2. * C4 * pow(m_rho, 2)) +
                                0.5 * s) +
                           pow(a1_mass, 2) *
                               (2. * C4 * pow(m_rho, 4) +
                                pow(pion_mass, 2) *
                                    (2. - 4. * C4 * pow(m_rho, 2)) +
                                pow(m_rho, 2) * (-1.5 - 4. * C4 * s) +
                                s * (0.5 + 2. * C4 * s))) +
                   eta1 *
                       (pow(a1_mass, 4) * (1. - 2. * C4 * pow(m_rho, 2)) +
                        pow(pion_mass, 4) * (1. - 2. * C4 * pow(m_rho, 2)) +
                        (-0.5 * pow(m_rho, 2) + 0.5 * s) * s +
                        pow(a1_mass, 2) *
                            (-1. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4) +
                             pow(pion_mass, 2) *
                                 (-2. + 4. * C4 * pow(m_rho, 2)) -
                             2. * C4 * pow(s, 2)) +
                        pow(pion_mass, 2) *
                            (-4. * C4 * pow(m_rho, 4) -
                             1. * s + pow(m_rho, 2) * (2. + 4. * C4 * s))))) *
             log(fabs(-1. * pow(a1_mass, 2) + tmin))) /
                pow(m_rho, 2) -
            (HeavisideTheta(-m_rho + sqrt(s)) *
             (-0.5 *
                  (eta1 * eta2 *
                       (0.25 * pow(m_rho, 6) + 1.5 * pow(m_rho, 4) * s -
                        0.125 * delta * pow(m_rho, 4) * s -
                        0.75 * pow(m_rho, 2) * pow(s, 2) -
                        0.75 * delta * pow(m_rho, 2) * pow(s, 2) +
                        0.375 * delta * pow(s, 3) +
                        pow(a1_mass, 4) *
                            (-2. * pow(m_rho, 2) + 1. * delta * s) +
                        pow(pion_mass, 4) *
                            (-6. * pow(m_rho, 2) + 3. * delta * s) +
                        pow(pion_mass, 2) *
                            (-3. * pow(m_rho, 4) +
                             (3. + 1.5 * delta) * pow(m_rho, 2) * s -
                             1.5 * delta * pow(s, 2)) +
                        pow(a1_mass, 2) *
                            (0.5 * pow(m_rho, 4) +
                             (-2.5 - 0.25 * delta) * pow(m_rho, 2) * s +
                             1.25 * delta * pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (6. * pow(m_rho, 2) - 3. * delta * s))) +
                   pow(eta1, 2) *
                       (0.5 * pow(m_rho, 6) - 2. * pow(m_rho, 4) * s -
                        0.25 * delta * pow(m_rho, 4) * s +
                        0.5 * pow(m_rho, 2) * pow(s, 2) +
                        1. * delta * pow(m_rho, 2) * pow(s, 2) -
                        0.25 * delta * pow(s, 3) +
                        pow(pion_mass, 4) *
                            (3. * pow(m_rho, 2) - 1.5 * delta * s) +
                        pow(a1_mass, 4) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 2) *
                            (4. * pow(m_rho, 4) +
                             (-2. - 2. * delta) * pow(m_rho, 2) * s +
                             1. * delta * pow(s, 2)) +
                        pow(a1_mass, 2) *
                            (-1.5 * pow(m_rho, 4) +
                             (1.5 + 0.75 * delta) * pow(m_rho, 2) * s -
                             0.75 * delta * pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (-3. * pow(m_rho, 2) + 1.5 * delta * s))) +
                   pow(eta2, 2) *
                       (-0.75 * pow(m_rho, 6) + 0.5 * pow(m_rho, 4) * s +
                        0.375 * delta * pow(m_rho, 4) * s +
                        0.25 * pow(m_rho, 2) * pow(s, 2) -
                        0.25 * delta * pow(m_rho, 2) * pow(s, 2) -
                        0.125 * delta * pow(s, 3) +
                        pow(pion_mass, 4) *
                            (3. * pow(m_rho, 2) - 1.5 * delta * s) +
                        pow(a1_mass, 4) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 2) *
                            (-1. * pow(m_rho, 4) +
                             (-1. + 0.5 * delta) * pow(m_rho, 2) * s +
                             0.5 * delta * pow(s, 2)) +
                        pow(a1_mass, 2) *
                            (1. * pow(m_rho, 4) +
                             (1. - 0.5 * delta) * pow(m_rho, 2) * s -
                             0.5 * delta *
                                 pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (-3. * pow(m_rho, 2) + 1.5 * delta * s)))) *
                  tmin -
              0.25 *
                  (eta1 * eta2 *
                       (0.5 * pow(m_rho, 4) - 2.5 * pow(m_rho, 2) * s -
                        0.25 * delta * pow(m_rho, 2) * s +
                        1.25 * delta * pow(s, 2) +
                        pow(pion_mass, 2) *
                            (6. * pow(m_rho, 2) - 3. * delta * s) +
                        pow(a1_mass, 2) *
                            (-2. * pow(m_rho, 2) + 1. * delta * s)) +
                   pow(eta1, 2) *
                       (-1.5 * pow(m_rho, 4) + 1.5 * pow(m_rho, 2) * s +
                        0.75 * delta * pow(m_rho, 2) * s -
                        0.75 * delta * pow(s, 2) +
                        pow(a1_mass, 2) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 2) *
                            (-3. * pow(m_rho, 2) + 1.5 * delta * s)) +
                   pow(eta2, 2) *
                       (1. * pow(m_rho, 4) + 1. * pow(m_rho, 2) * s -
                        0.5 * delta * pow(m_rho, 2) * s -
                        0.5 * delta * pow(s, 2) +
                        pow(a1_mass, 2) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 2) *
                            (-3. * pow(m_rho, 2) + 1.5 * delta * s))) *
                  pow(tmin, 2) -
              0.16666666666666666 *
                  (pow(eta1, 2) * (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                   pow(eta2, 2) * (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                   eta1 * eta2 * (-2. * pow(m_rho, 2) + 1. * delta * s)) *
                  pow(tmin, 3) -
              0.5 *
                  (eta1 * eta2 *
                       (pow(pion_mass, 6) *
                            (2. * pow(m_rho, 2) - 1. * delta * s) +
                        pow(a1_mass, 6) *
                            (-2. * pow(m_rho, 2) + 1. * delta * s) +
                        pow(pion_mass, 4) *
                            (2.5 * pow(m_rho, 4) +
                             (-0.5 - 1.25 * delta) * pow(m_rho, 2) * s +
                             0.25 * delta * pow(s, 2)) +
                        s * (-0.25 * pow(m_rho, 6) +
                             0.125 * delta * pow(m_rho, 4) * s +
                             0.25 * pow(m_rho, 2) * pow(s, 2) -
                             0.125 * delta * pow(s, 3)) +
                        pow(pion_mass, 2) *
                            (-0.25 * pow(m_rho, 6) +
                             (0.5 + 0.125 * delta) * pow(m_rho, 4) * s +
                             (-1.25 - 0.25 * delta) * pow(m_rho, 2) *
                                 pow(s, 2) +
                             0.625 * delta *
                                 pow(s, 3)) +
                        pow(a1_mass, 4) *
                            (0.5 * pow(m_rho, 4) +
                             (-2.5 - 0.25 * delta) * pow(m_rho, 2) * s +
                             1.25 * delta * pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (6. * pow(m_rho, 2) - 3. * delta * s)) +
                        pow(a1_mass, 2) *
                            (0.25 * pow(m_rho, 6) +
                             (1.5 - 0.125 * delta) * pow(m_rho, 4) * s +
                             (-0.75 - 0.75 * delta) * pow(m_rho, 2) *
                                 pow(s, 2) +
                             0.375 * delta * pow(s, 3) +
                             pow(pion_mass, 4) *
                                 (-6. * pow(m_rho, 2) + 3. * delta * s) +
                             pow(pion_mass, 2) * (-3. * pow(m_rho, 4) +
                                                  (3. + 1.5 * delta) *
                                                      pow(m_rho, 2) * s -
                                                  1.5 * delta * pow(s, 2)))) +
                   pow(eta2, 2) *
                       (pow(a1_mass, 6) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 2) *
                            (0.25 * pow(m_rho, 6) +
                             (-0.5 - 0.125 * delta) * pow(m_rho, 4) * s +
                             (0.25 + 0.25 * delta) * pow(m_rho, 2) * pow(s, 2) -
                             0.125 * delta *
                                 pow(s, 3) +
                             pow(pion_mass, 4) *
                                 (-1. * pow(m_rho, 2) + 0.5 * delta * s)) +
                        pow(a1_mass, 4) *
                            (1. * pow(m_rho, 4) +
                             (1. - 0.5 * delta) * pow(m_rho, 2) * s -
                             0.5 * delta *
                                 pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (-3. * pow(m_rho, 2) + 1.5 * delta * s)) +
                        pow(a1_mass, 2) *
                            (-0.75 * pow(m_rho, 6) +
                             (0.5 + 0.375 * delta) * pow(m_rho, 4) * s +
                             (0.25 - 0.25 * delta) * pow(m_rho, 2) * pow(s, 2) -
                             0.125 * delta * pow(s, 3) +
                             pow(pion_mass, 4) *
                                 (3. * pow(m_rho, 2) - 1.5 * delta * s) +
                             pow(pion_mass, 2) * (-1. * pow(m_rho, 4) +
                                                  (-1. + 0.5 * delta) *
                                                      pow(m_rho, 2) * s +
                                                  0.5 * delta * pow(s, 2)))) +
                   pow(eta1, 2) *
                       (pow(a1_mass, 6) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 2) * pow(s, 2) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 6) *
                            (-1. * pow(m_rho, 2) + 0.5 * delta * s) +
                        pow(pion_mass, 4) *
                            (-2.5 * pow(m_rho, 4) +
                             (0.5 + 1.25 * delta) * pow(m_rho, 2) * s -
                             0.25 * delta * pow(s, 2)) +
                        s * (0.25 * pow(m_rho, 6) -
                             0.125 * delta * pow(m_rho, 4) * s -
                             0.25 * pow(m_rho, 2) * pow(s, 2) +
                             0.125 * delta * pow(s, 3)) +
                        pow(a1_mass, 4) *
                            (-1.5 * pow(m_rho, 4) +
                             (1.5 + 0.75 * delta) * pow(m_rho, 2) * s -
                             0.75 * delta * pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (-3. * pow(m_rho, 2) + 1.5 * delta * s)) +
                        pow(a1_mass, 2) *
                            (0.5 * pow(m_rho, 6) +
                             (-2. - 0.25 * delta) * pow(m_rho, 4) * s +
                             (0.5 + 1. * delta) * pow(m_rho, 2) * pow(s, 2) -
                             0.25 * delta * pow(s, 3) +
                             pow(pion_mass, 4) *
                                 (3. * pow(m_rho, 2) - 1.5 * delta * s) +
                             pow(pion_mass, 2) *
                                 (4. * pow(m_rho, 4) +
                                  (-2. - 2. * delta) * pow(m_rho, 2) * s +
                                  1. * delta * pow(s, 2))))) *
                  log(fabs(-1. * pow(a1_mass, 2) + tmin)))) /
                (pow(m_rho, 2) * (pow(m_rho, 2) - 1. * s)) -
            0.5 * pow(-2. + delta, 2) * pow(pion_mass, 2) *
                log(fabs(-1. * pow(pion_mass, 2) + tmin)) +
            (0.5 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta2 * pow(pion_mass, 4) * (-1. * pow(m_rho, 2) + 1. * s) +
              eta1 * pow(pion_mass, 2) *
                  (-0.5 * pow(m_rho, 4) +
                   pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) +
                   0.5 * pow(m_rho, 2) * s)) *
             log(fabs(-1. * pow(pion_mass, 2) + tmin))) /
                (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) -
            (0.5 * (pow(m_rho, 2) - 1. * s) *
                 ((0.5 - 0.25 * delta) * pow(m_rho, 2) +
                  C4 * (-2. + 1. * delta) * pow(m_rho, 4) +
                  (-0.25 + 0.125 * delta) * delta * s) *
                 pow(tmin, 2) +
             tmin * (-0.5 * pow(m_rho, 6) + 0.25 * delta * pow(m_rho, 6) +
                     2. * C4 * pow(m_rho, 8) - 1. * C4 * delta * pow(m_rho, 8) +
                     1. * pow(m_rho, 4) * s + 0.25 * delta * pow(m_rho, 4) * s -
                     0.375 * pow(delta, 2) * pow(m_rho, 4) * s -
                     6. * C4 * pow(m_rho, 6) * s +
                     3. * C4 * delta * pow(m_rho, 6) * s -
                     0.5 * pow(m_rho, 2) * pow(s, 2) -
                     0.75 * delta * pow(m_rho, 2) * pow(s, 2) +
                     0.5 * pow(delta, 2) * pow(m_rho, 2) * pow(s, 2) +
                     4. * C4 * pow(m_rho, 4) * pow(s, 2) -
                     2. * C4 * delta * pow(m_rho, 4) * pow(s, 2) +
                     0.25 * delta * pow(s, 3) -
                     0.125 * pow(delta, 2) * pow(s, 3) +
                     pow(pion_mass, 2) *
                         (C4 * (2. - 1. * delta) * pow(m_rho, 6) +
                          (0.5 - 0.125 * pow(delta, 2)) * pow(m_rho, 2) * s +
                          (-1.25 + 0.625 * delta) * delta * pow(s, 2) +
                          pow(m_rho, 4) *
                              (-0.5 + 1.25 * delta - 0.5 * pow(delta, 2) -
                               2. * C4 * s + 1. * C4 * delta * s)) +
                     (pow(pion_mass, 2) *
                          ((-2. + 1. * delta) * pow(m_rho, 4) -
                           0.5000000000000001 * pow(2. - 1. * delta, 2) *
                               pow(m_rho, 2) * s +
                           (1. - 0.5 * delta) * delta * pow(s, 2)) +
                      s * ((-0.5 + 0.25 * delta) * pow(m_rho, 4) +
                           (0.5 - 0.125 * pow(delta, 2)) * pow(m_rho, 2) * s +
                           (-0.25 + 0.125 * delta) * delta * pow(s, 2))) *
                         HeavisideTheta(-m_rho + sqrt(s))) +
             pow(m_rho, 2) *
                 (s * ((0.5 - 0.75 * delta + 0.25 * pow(delta, 2)) *
                           pow(m_rho, 4) -
                       0.12500000000000003 * pow(2. - 1. * delta, 2) *
                           pow(m_rho, 2) * s +
                       (0.25 - 0.125 * delta) * delta * pow(s, 2)) +
                  pow(pion_mass, 2) *
                      (C4 * (4. - 2. * delta) * pow(m_rho, 6) +
                       delta * (-3. + 1.5 * delta) * pow(s, 2) +
                       pow(m_rho, 2) *
                           s * (2. - 1.5 * pow(delta, 2) + 4. * C4 * s + delta *
                           (2. - 2. * C4 * s)) +
                       pow(m_rho, 4) *
                           (-2. - 8. * C4 * s + delta * (1. + 4. * C4 * s))) +
                  s * ((0.5 - 0.25 * delta) * pow(m_rho, 4) +
                  0.12500000000000003 * pow(2. - 1. * delta, 2) * pow(m_rho, 2)
                  * s + (-0.25 + 0.125 * delta) * delta * pow(s, 2) +
                  pow(pion_mass, 2) * ((-4. + 2. * delta) * pow(m_rho, 2) +
                  (2. - 1. * delta) * delta * s)) *
                      HeavisideTheta(-m_rho + sqrt(s))) *
                 log(fabs(-1. * pow(pion_mass, 2) + tmin))) /
                (pow(m_rho, 4) * (pow(m_rho, 2) - 1. * s)) +
            0.0625 * pow(eta1 - 1. * eta2, 2) *
                (eta1 * eta2 *
                     (-4. * pow(a1_mass, 6) +
                      pow(a1_mass, 4) * (12. * pow(pion_mass, 2) - 6. * s) +
                      pow(pion_mass, 2) *
                          (4. * pow(pion_mass, 4) - 4. * pow(m_rho, 4) -
                           2. *
                               pow(pion_mass, 2) * s +
                           2. *
                               pow(m_rho, 2) * s) +
                      pow(a1_mass, 2) *
                          (-12. * pow(pion_mass, 4) + 2. * pow(m_rho, 4) +
                           8. *
                               pow(pion_mass, 2) * s +
                           4. *
                               pow(m_rho, 2) * s -
                           4. *
                               pow(s, 2))) +
                 pow(eta1, 2) *
                     (2. * pow(a1_mass, 6) - 2. * pow(pion_mass, 6) +
                      pow(pion_mass, 2) * pow(m_rho, 2) * s +
                      pow(m_rho, 2) *
                          (pow(m_rho, 2) - 1. * s) * s +
                      pow(pion_mass, 4) * (-3. * pow(m_rho, 2) + s) +
                      pow(a1_mass, 4) *
                          (-6. * pow(pion_mass, 2) - 3. * pow(m_rho, 2) +
                           3. *
                               s) +
                      pow(a1_mass, 2) *
                          (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                           pow(pion_mass, 2) * (6. * pow(m_rho, 2) - 4. * s) -
                           4. *
                               pow(m_rho, 2) * s +
                           2. *
                               pow(s, 2))) +
                 pow(eta2, 2) *
                     (2. * pow(a1_mass, 6) - 2. * pow(pion_mass, 6) -
                      1. *
                          pow(pion_mass, 2) * pow(m_rho, 2) * s +
                      pow(pion_mass, 4) * (3. * pow(m_rho, 2) + s) +
                      pow(a1_mass, 4) *
                          (-6. * pow(pion_mass, 2) + 3. * pow(m_rho, 2) +
                           3. *
                               s) +
                      pow(a1_mass, 2) *
                          (6. * pow(pion_mass, 4) + pow(m_rho, 4) +
                           pow(pion_mass, 2) * (-6. * pow(m_rho, 2) - 4. * s) -
                           2. *
                               pow(m_rho, 2) * s +
                           2. *
                               pow(s, 2)))) *
                log(fabs(-1. * pow(a1_mass, 2) + tmax)) -
            (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta2 * (-0.5 * pow(a1_mass, 6) - 0.5 * pow(pion_mass, 6) +
                      0.5 * pow(pion_mass, 4) * pow(m_rho, 2) +
                      pow(a1_mass, 4) * (0.5 * pow(pion_mass, 2) +
                                         0.5 * pow(m_rho, 2) - 1. * s) +
                      pow(a1_mass, 2) * pow(pion_mass, 2) *
                          (0.5 * pow(pion_mass, 2) + 1. * pow(m_rho, 2) -
                           1. *
                               s)) +
              eta1 *
                  (pow(a1_mass, 4) * (1. * pow(pion_mass, 2) + 0.5 * s) +
                   pow(pion_mass, 2) *
                       (1. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                        pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 0.5 * s) -
                        0.5 * pow(m_rho, 2) * s) +
                   pow(a1_mass, 2) *
                       (-2. * pow(pion_mass, 4) - 0.5 * pow(m_rho, 2) * s +
                        pow(pion_mass, 2) * (-1. * pow(m_rho, 2) + 1. * s)))) *
             log(fabs(-1. * pow(a1_mass, 2) + tmax))) /
                (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) +
            (0.25 * (1. * eta1 - 1. * eta2) *
             (delta *
                  (eta1 *
                       (1. * pow(a1_mass, 6) - 1. * pow(pion_mass, 6) +
                        pow(pion_mass, 4) * (-2.5 * pow(m_rho, 2) + 0.5 * s) +
                        pow(pion_mass, 2) * s *
                            (-0.5 * pow(m_rho, 2) + 1. * s) +
                        pow(a1_mass, 4) * (-3. * pow(pion_mass, 2) -
                                           1.5 * pow(m_rho, 2) + 1.5 * s) +
                        s * (0.5 * pow(m_rho, 4) - 0.25 * pow(m_rho, 2) * s -
                             0.25 * pow(s, 2)) +
                        pow(a1_mass, 2) *
                            (3. * pow(pion_mass, 4) + 0.5 * pow(m_rho, 4) +
                             pow(pion_mass, 2) * (4. * pow(m_rho, 2) - 2. * s) -
                             2. *
                                 pow(m_rho, 2) * s +
                             1. *
                                 pow(s, 2))) +
                   eta2 * (-1. * pow(a1_mass, 6) +
                           pow(a1_mass, 4) *
                               (3. * pow(pion_mass, 2) - 1. * pow(m_rho, 2) -
                                1. *
                                    s) +
                           pow(pion_mass, 2) *
                               (1. * pow(pion_mass, 4) - 0.5 * pow(m_rho, 4) +
                                0.25 * pow(m_rho, 2) * s - 0.25 * pow(s, 2)) +
                           pow(a1_mass, 2) *
                               (-3. * pow(pion_mass, 4) + 1. * pow(m_rho, 4) +
                                0.25 * pow(m_rho, 2) * s - 0.75 * pow(s, 2) +
                                pow(pion_mass, 2) *
                                    (1. * pow(m_rho, 2) + 1. * s)))) +
              pow(m_rho, 2) *
                  (eta2 *
                       (pow(a1_mass, 4) * (-1. + 2. * C4 * pow(m_rho, 2)) +
                        pow(pion_mass, 2) *
                            (0.5 * pow(m_rho, 2) +
                             pow(pion_mass, 2) *
                                 (-1. + 2. * C4 * pow(m_rho, 2)) +
                             0.5 * s) +
                        pow(a1_mass, 2) *
                            (2. * C4 * pow(m_rho, 4) +
                             pow(pion_mass, 2) *
                                 (2. - 4. * C4 * pow(m_rho, 2)) +
                             pow(m_rho, 2) * (-1.5 - 4. * C4 * s) +
                             s * (0.5 + 2. * C4 * s))) +
                   eta1 * (pow(a1_mass, 4) * (1. - 2. * C4 * pow(m_rho, 2)) +
                           pow(pion_mass, 4) * (1. - 2. * C4 * pow(m_rho, 2)) +
                           (-0.5 * pow(m_rho, 2) + 0.5 * s) * s +
                           pow(a1_mass, 2) *
                               (-1. * pow(m_rho, 2) + 2. * C4 * pow(m_rho, 4) +
                                pow(pion_mass, 2) *
                                    (-2. + 4. * C4 * pow(m_rho, 2)) -
                                2. *
                                    C4 * pow(s, 2)) +
                           pow(pion_mass, 2) *
                               (-4. * C4 * pow(m_rho, 4) - 1. * s +
                                pow(m_rho, 2) * (2. + 4. * C4 * s))))) *
             log(fabs(-1. * pow(a1_mass, 2) + tmax))) /
                pow(m_rho, 2) +
            (HeavisideTheta(-m_rho + sqrt(s)) *
             (-0.5 *
                  (eta1 * eta2 *
                       (0.25 * pow(m_rho, 6) + 1.5 * pow(m_rho, 4) * s -
                        0.125 * delta * pow(m_rho, 4) * s -
                        0.75 * pow(m_rho, 2) * pow(s, 2) -
                        0.75 * delta * pow(m_rho, 2) * pow(s, 2) +
                        0.375 * delta * pow(s, 3) +
                        pow(a1_mass, 4) *
                            (-2. * pow(m_rho, 2) + 1. * delta * s) +
                        pow(pion_mass, 4) *
                            (-6. * pow(m_rho, 2) + 3. * delta * s) +
                        pow(pion_mass, 2) * (-3. * pow(m_rho, 4) +
                                             (3. + 1.5 * delta) *
                                                 pow(m_rho, 2) * s -
                                             1.5 * delta * pow(s, 2)) +
                        pow(a1_mass, 2) *
                            (0.5 * pow(m_rho, 4) +
                             (-2.5 - 0.25 * delta) * pow(m_rho, 2) * s +
                             1.25 * delta * pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (6. * pow(m_rho, 2) - 3. * delta * s))) +
                   pow(eta1, 2) *
                       (0.5 * pow(m_rho, 6) - 2. * pow(m_rho, 4) * s -
                        0.25 * delta * pow(m_rho, 4) * s +
                        0.5 * pow(m_rho, 2) * pow(s, 2) +
                        1. *
                            delta * pow(m_rho, 2) * pow(s, 2) -
                        0.25 * delta * pow(s, 3) +
                        pow(pion_mass, 4) *
                            (3. * pow(m_rho, 2) - 1.5 * delta * s) +
                        pow(a1_mass, 4) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 2) * (4. * pow(m_rho, 4) +
                                             (-2. - 2. * delta) *
                                                 pow(m_rho, 2) * s +
                                             1. *
                                                 delta * pow(s, 2)) +
                        pow(a1_mass, 2) *
                            (-1.5 * pow(m_rho, 4) +
                             (1.5 + 0.75 * delta) * pow(m_rho, 2) * s -
                             0.75 * delta * pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (-3. * pow(m_rho, 2) + 1.5 * delta * s))) +
                   pow(eta2, 2) *
                       (-0.75 * pow(m_rho, 6) + 0.5 * pow(m_rho, 4) * s +
                        0.375 * delta *
                            pow(m_rho, 4) * s +
                        0.25 * pow(m_rho, 2) * pow(s, 2) -
                        0.25 * delta * pow(m_rho, 2) * pow(s, 2) -
                        0.125 * delta *
                            pow(s, 3) +
                        pow(pion_mass, 4) *
                            (3. * pow(m_rho, 2) - 1.5 * delta * s) +
                        pow(a1_mass, 4) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 2) *
                            (-1. * pow(m_rho, 4) +
                             (-1. + 0.5 * delta) * pow(m_rho, 2) * s +
                             0.5 * delta * pow(s, 2)) +
                        pow(a1_mass, 2) *
                            (1. * pow(m_rho, 4) +
                             (1. - 0.5 * delta) * pow(m_rho, 2) * s -
                             0.5 * delta * pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (-3. * pow(m_rho, 2) + 1.5 * delta * s)))) *
                  tmax -
              0.25 *
                  (eta1 * eta2 *
                       (0.5 * pow(m_rho, 4) - 2.5 * pow(m_rho, 2) * s -
                        0.25 * delta * pow(m_rho, 2) * s +
                        1.25 * delta * pow(s, 2) +
                        pow(pion_mass, 2) *
                            (6. * pow(m_rho, 2) - 3. * delta * s) +
                        pow(a1_mass, 2) *
                            (-2. * pow(m_rho, 2) + 1. * delta * s)) +
                   pow(eta1, 2) *
                       (-1.5 * pow(m_rho, 4) + 1.5 * pow(m_rho, 2) * s +
                        0.75 * delta * pow(m_rho, 2) * s -
                        0.75 * delta * pow(s, 2) +
                        pow(a1_mass, 2) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 2) *
                            (-3. * pow(m_rho, 2) + 1.5 * delta * s)) +
                   pow(eta2, 2) * (1. * pow(m_rho, 4) + 1. * pow(m_rho, 2) * s -
                                   0.5 * delta * pow(m_rho, 2) * s -
                                   0.5 * delta * pow(s, 2) +
                                   pow(a1_mass, 2) *
                                       (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                                   pow(pion_mass, 2) * (-3. * pow(m_rho, 2) +
                                                        1.5 * delta * s))) *
                  pow(tmax, 2) -
              0.16666666666666666 *
                  (pow(eta1, 2) * (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                   pow(eta2, 2) * (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                   eta1 * eta2 * (-2. * pow(m_rho, 2) + 1. * delta * s)) *
                  pow(tmax, 3) -
              0.5 *
                  (eta1 * eta2 *
                       (pow(pion_mass, 6) *
                            (2. * pow(m_rho, 2) - 1. * delta * s) +
                        pow(a1_mass, 6) *
                            (-2. * pow(m_rho, 2) + 1. * delta * s) +
                        pow(pion_mass, 4) *
                            (2.5 * pow(m_rho, 4) +
                             (-0.5 - 1.25 * delta) * pow(m_rho, 2) * s +
                             0.25 * delta * pow(s, 2)) +
                        s * (-0.25 * pow(m_rho, 6) +
                             0.125 * delta * pow(m_rho, 4) * s +
                             0.25 * pow(m_rho, 2) * pow(s, 2) -
                             0.125 * delta * pow(s, 3)) +
                        pow(pion_mass, 2) *
                            (-0.25 * pow(m_rho, 6) +
                             (0.5 + 0.125 * delta) * pow(m_rho, 4) * s +
                             (-1.25 - 0.25 * delta) * pow(m_rho, 2) *
                                 pow(s, 2) +
                             0.625 * delta * pow(s, 3)) +
                        pow(a1_mass, 4) *
                            (0.5 * pow(m_rho, 4) +
                             (-2.5 - 0.25 * delta) * pow(m_rho, 2) * s +
                             1.25 * delta * pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (6. * pow(m_rho, 2) - 3. * delta * s)) +
                        pow(a1_mass, 2) *
                            (0.25 * pow(m_rho, 6) +
                             (1.5 - 0.125 * delta) * pow(m_rho, 4) * s +
                             (-0.75 - 0.75 * delta) * pow(m_rho, 2) *
                                 pow(s, 2) +
                             0.375 * delta * pow(s, 3) +
                             pow(pion_mass, 4) *
                                 (-6. * pow(m_rho, 2) + 3. * delta * s) +
                             pow(pion_mass, 2) * (-3. * pow(m_rho, 4) +
                                                  (3. + 1.5 * delta) *
                                                      pow(m_rho, 2) * s -
                                                  1.5 * delta * pow(s, 2)))) +
                   pow(eta2, 2) *
                       (pow(a1_mass, 6) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 2) *
                            (0.25 * pow(m_rho, 6) +
                             (-0.5 - 0.125 * delta) * pow(m_rho, 4) * s +
                             (0.25 + 0.25 * delta) * pow(m_rho, 2) * pow(s, 2) -
                             0.125 * delta *
                                 pow(s, 3) +
                             pow(pion_mass, 4) *
                                 (-1. * pow(m_rho, 2) + 0.5 * delta * s)) +
                        pow(a1_mass, 4) *
                            (1. * pow(m_rho, 4) +
                             (1. - 0.5 * delta) * pow(m_rho, 2) * s -
                             0.5 * delta * pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (-3. * pow(m_rho, 2) + 1.5 * delta * s)) +
                        pow(a1_mass, 2) *
                            (-0.75 * pow(m_rho, 6) +
                             (0.5 + 0.375 * delta) * pow(m_rho, 4) * s +
                             (0.25 - 0.25 * delta) * pow(m_rho, 2) * pow(s, 2) -
                             0.125 * delta *
                                 pow(s, 3) +
                             pow(pion_mass, 4) *
                                 (3. * pow(m_rho, 2) - 1.5 * delta * s) +
                             pow(pion_mass, 2) *
                                 (-1. * pow(m_rho, 4) +
                                  (-1. + 0.5 * delta) * pow(m_rho, 2) * s +
                                  0.5 * delta * pow(s, 2)))) +
                   pow(eta1, 2) *
                       (pow(a1_mass, 6) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 2) * pow(s, 2) *
                            (1. * pow(m_rho, 2) - 0.5 * delta * s) +
                        pow(pion_mass, 6) *
                            (-1. * pow(m_rho, 2) + 0.5 * delta * s) +
                        pow(pion_mass, 4) *
                            (-2.5 * pow(m_rho, 4) +
                             (0.5 + 1.25 * delta) * pow(m_rho, 2) * s -
                             0.25 * delta * pow(s, 2)) +
                        s * (0.25 * pow(m_rho, 6) -
                             0.125 * delta *
                                 pow(m_rho, 4) * s -
                             0.25 * pow(m_rho, 2) * pow(s, 2) +
                             0.125 * delta *
                                 pow(s, 3)) +
                        pow(a1_mass, 4) *
                            (-1.5 * pow(m_rho, 4) +
                             (1.5 + 0.75 * delta) * pow(m_rho, 2) * s -
                             0.75 * delta * pow(s, 2) +
                             pow(pion_mass, 2) *
                                 (-3. * pow(m_rho, 2) + 1.5 * delta * s)) +
                        pow(a1_mass, 2) *
                            (0.5 * pow(m_rho, 6) +
                             (-2. - 0.25 * delta) * pow(m_rho, 4) * s +
                             (0.5 + 1. * delta) * pow(m_rho, 2) * pow(s, 2) -
                             0.25 * delta * pow(s, 3) +
                             pow(pion_mass, 4) *
                                 (3. * pow(m_rho, 2) - 1.5 * delta * s) +
                             pow(pion_mass, 2) *
                                 (4. * pow(m_rho, 4) +
                                  (-2. - 2. * delta) *
                                      pow(m_rho, 2) * s +
                                  1. *
                                      delta * pow(s, 2))))) *
                  log(fabs(-1. * pow(a1_mass, 2) + tmax)))) /
                (pow(m_rho, 2) * (pow(m_rho, 2) - 1. * s)) +
            0.5 * pow(-2. + delta, 2) * pow(pion_mass, 2) *
                log(fabs(-1. * pow(pion_mass, 2) + tmax)) -
            (0.5 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta2 * pow(pion_mass, 4) * (-1. * pow(m_rho, 2) + 1. * s) +
              eta1 * pow(pion_mass, 2) *
                  (-0.5 * pow(m_rho, 4) +
                   pow(pion_mass, 2) * (1. * pow(m_rho, 2) - 1. * s) +
                   0.5 * pow(m_rho, 2) * s)) *
             log(fabs(-1. * pow(pion_mass, 2) + tmax))) /
                (pow(a1_mass, 2) - 1. * pow(pion_mass, 2)) +
            (0.5 * (pow(m_rho, 2) - 1. * s) *
                 ((0.5 - 0.25 * delta) * pow(m_rho, 2) +
                  C4 *
                      (-2. + 1. * delta) * pow(m_rho, 4) +
                  (-0.25 + 0.125 * delta) * delta * s) *
                 pow(tmax, 2) +
             tmax * (-0.5 * pow(m_rho, 6) +
                     0.25 * delta * pow(m_rho, 6) + 2. * C4 * pow(m_rho, 8) -
                     1. *
                         C4 * delta * pow(m_rho, 8) +
                     1. *
                         pow(m_rho, 4) * s +
                     0.25 * delta * pow(m_rho, 4) * s -
                     0.375 * pow(delta, 2) * pow(m_rho, 4) * s -
                     6. *
                         C4 * pow(m_rho, 6) * s +
                     3. *
                         C4 * delta * pow(m_rho, 6) * s -
                     0.5 * pow(m_rho, 2) * pow(s, 2) -
                     0.75 * delta * pow(m_rho, 2) * pow(s, 2) +
                     0.5 * pow(delta, 2) * pow(m_rho, 2) * pow(s, 2) +
                     4. *
                         C4 * pow(m_rho, 4) * pow(s, 2) -
                     2. *
                         C4 * delta * pow(m_rho, 4) * pow(s, 2) +
                     0.25 * delta * pow(s, 3) -
                     0.125 * pow(delta, 2) * pow(s, 3) +
                     pow(pion_mass, 2) *
                         (C4 * (2. - 1. * delta) * pow(m_rho, 6) +
                          (0.5 - 0.125 * pow(delta, 2)) * pow(m_rho, 2) * s +
                          (-1.25 + 0.625 * delta) * delta * pow(s, 2) +
                          pow(m_rho, 4) *
                              (-0.5 + 1.25 * delta -
                               0.5 * pow(delta, 2) - 2. * C4 * s +
                               1. *
                                   C4 * delta *
                                   s)) +
                     (pow(pion_mass, 2) *
                          ((-2. + 1. * delta) * pow(m_rho, 4) -
                           0.5000000000000001 *
                               pow(2. - 1. * delta, 2) * pow(m_rho, 2) * s +
                           (1. - 0.5 * delta) *
                               delta * pow(s, 2)) +
                      s * ((-0.5 + 0.25 * delta) * pow(m_rho, 4) +
                           (0.5 - 0.125 * pow(delta, 2)) * pow(m_rho, 2) * s +
                           (-0.25 + 0.125 * delta) * delta * pow(s, 2))) *
                         HeavisideTheta(-m_rho + sqrt(s))) +
             pow(m_rho, 2) *
                 (s * ((0.5 - 0.75 * delta + 0.25 * pow(delta, 2)) *
                           pow(m_rho, 4) -
                       0.12500000000000003 * pow(2. - 1. * delta, 2) *
                           pow(m_rho, 2) * s +
                       (0.25 - 0.125 * delta) * delta * pow(s, 2)) +
                  pow(pion_mass, 2) *
                      (C4 * (4. - 2. * delta) * pow(m_rho, 6) +
                       delta * (-3. + 1.5 * delta) * pow(s, 2) +
                       pow(m_rho, 2) * s *
                           (2. - 1.5 * pow(delta, 2) + 4. * C4 * s +
                            delta * (2. - 2. * C4 * s)) +
                       pow(m_rho, 4) *
                           (-2. - 8. * C4 * s + delta * (1. + 4. * C4 * s))) +
                  s * ((0.5 - 0.25 * delta) * pow(m_rho, 4) +
                  0.12500000000000003 * pow(2. - 1. * delta, 2) * pow(m_rho, 2)
                  * s + (-0.25 + 0.125 * delta) * delta * pow(s, 2) +
                  pow(pion_mass, 2) * ((-4. + 2. * delta) * pow(m_rho, 2) +
                  (2. - 1. * delta) * delta * s)) *
                      HeavisideTheta(-m_rho + sqrt(s))) *
                 log(fabs(-1. * pow(pion_mass, 2) + tmax))) /
                (pow(m_rho, 4) * (pow(m_rho, 2) - 1. * s)))) /
          (16. * Pi * s * (-4 * pow(pion_mass, 2) + s));

  // clang-format on
  return cut_off(xs * gev2_mb / spin_deg_factor);
}

double CrosssectionsPhoton<ComputationMethod::Analytic>::xs_diff_pi_pi0_rho(
    const double s, const double t, const double m_rho) {
  const double spin_deg_factor = 1.0;

  using std::abs;
  using std::atan;
  using std::pow;
  using std::sqrt;

  // clang-format off
  const double diff_xs =
      (pow(Const, 2) * pow(ghat, 4) *
       (0.75 +
        C4 * (2. * C4 * pow(m_rho, 4) + pow(m_rho, 2) * (-3. - 4. * C4 * s) +
              s * (3. + 2. * C4 * s)) -
        (0.25 * pow(-2 + delta, 2) * pow(pion_mass, 2) *
         (pow(pion_mass, 4) + pow(pow(m_rho, 2) - t, 2) -
          2 * pow(pion_mass, 2) * (pow(m_rho, 2) + t))) /
            (pow(m_rho, 2) * pow(pow(pion_mass, 2) - t, 2)) +
        (pow(delta, 2) *
         (0.5 * pow(pion_mass, 4) * pow(m_rho, 2) + 0.125 * pow(m_rho, 6) +
          pow(pion_mass, 2) * (1.25 * pow(m_rho, 4) - 0.375 * pow(s, 2) +
                               pow(m_rho, 2) * (0.125 * s - 1. * t)) +
          pow(m_rho, 4) * (-0.375 * s - 0.5 * t) +
          pow(s, 2) * (0.125 * s + 0.125 * t) +
          pow(m_rho, 2) *
              (0.0625 * pow(s, 2) + 0.375 * s * t + 0.5 * pow(t, 2)))) /
            pow(m_rho, 6) +
        (delta * (2. * C4 * pow(m_rho, 6) +
                  pow(pion_mass, 2) * (3. * C4 * pow(m_rho, 4) + 0.25 * s +
                                       pow(m_rho, 2) * (-1.25 - 1. * C4 * s)) +
                  s * (-0.25 * s - 0.25 * t) +
                  pow(m_rho, 4) * (-0.75 - 1.5 * C4 * s - 3. * C4 * t) +
                  pow(m_rho, 2) *
                      (1.25 * t + s * (0.25 - 0.5 * C4 * s + 1. * C4 * t)))) /
            pow(m_rho, 4) -
        (0.125 * (-2 + delta) * (eta1 - eta2) *
         (-(eta2 * (pow(pion_mass, 2) + t) *
            (pow(pion_mass, 4) + t * (-pow(m_rho, 2) + 2 * s + t) -
             pow(pion_mass, 2) * (pow(m_rho, 2) + 2 * t))) +
          eta1 * (2 * pow(pion_mass, 6) +
                  pow(pion_mass, 4) * (-2 * pow(m_rho, 2) + s - 4 * t) +
                  s * t * (-pow(m_rho, 2) + t) +
                  pow(pion_mass, 2) * (2 * pow(m_rho, 4) + 2 * t * (s + t) -
                                       pow(m_rho, 2) * (s + 2 * t))))) /
            ((-pow(a1_mass, 2) + t) * (-pow(pion_mass, 2) + t)) +
        (0.03125 * pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(pion_mass, 8) - 4 * pow(pion_mass, 6) * t +
               pow(t, 2) * (-pow(m_rho, 4) - 2 * pow(m_rho, 2) * s +
                            2 * pow(s, 2) + 2 * s * t + pow(t, 2)) -
               2 * pow(pion_mass, 2) * t *
                   (-2 * pow(m_rho, 4) + pow(m_rho, 2) * s + 2 * t * (s + t)) +
               pow(pion_mass, 4) * (-pow(m_rho, 4) + 2 * t * (s + 3 * t))) +
          pow(eta2, 2) *
              (pow(pion_mass, 8) -
               2 * pow(pion_mass, 6) * (pow(m_rho, 2) + 2 * t) +
               pow(t, 2) * (pow(m_rho, 4) + 2 * pow(s, 2) + 2 * s * t +
                            pow(t, 2) + 2 * pow(m_rho, 2) * (-s + t)) -
               2 * pow(pion_mass, 2) * t *
                   (2 * t * (s + t) + pow(m_rho, 2) * (s + 3 * t)) +
               pow(pion_mass, 4) * (pow(m_rho, 4) + 6 * pow(m_rho, 2) * t +
                                    2 * t * (s + 3 * t))) +
          pow(eta1, 2) *
              (pow(pion_mass, 8) +
               2 * pow(pion_mass, 6) * (pow(m_rho, 2) - 2 * t) -
               2 * pow(pion_mass, 2) * (pow(m_rho, 2) - s - t) *
                   (pow(m_rho, 4) + pow(m_rho, 2) * t - 2 * pow(t, 2)) +
               t * (-pow(m_rho, 2) + t) *
                   (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                    pow(m_rho, 2) * (2 * s + t)) +
               pow(pion_mass, 4) *
                   (pow(m_rho, 4) - 2 * pow(m_rho, 2) * (s + 3 * t) +
                    2 * t * (s + 3 * t))))) /
            pow(pow(a1_mass, 2) - t, 2) -
        (0.0625 * (eta1 - eta2) *
         (eta2 *
              (-4 * delta * pow(pion_mass, 6) +
               4 * pow(pion_mass, 4) *
                   (pow(m_rho, 2) - 2 * C4 * pow(m_rho, 4) + 3 * delta * t) +
               pow(pion_mass, 2) *
                   (delta * (s - 6 * t) * (s + 2 * t) -
                    (2 + delta) * pow(m_rho, 2) * (s + 4 * t) +
                    2 * pow(m_rho, 4) * (-1 + delta + 8 * C4 * t)) +
               t * (-8 * C4 * pow(m_rho, 6) +
                    pow(m_rho, 4) * (6 - 4 * delta + 16 * C4 * s - 8 * C4 * t) +
                    pow(m_rho, 2) * (-(s * (2 + delta + 8 * C4 * s)) +
                                     4 * (1 + delta) * t) +
                    delta * (3 * pow(s, 2) + 4 * s * t + 4 * pow(t, 2)))) +
          eta1 * (4 * delta * pow(pion_mass, 6) - 8 * C4 * pow(m_rho, 6) * t +
                  delta * (pow(s, 3) - 4 * pow(s, 2) * t - 6 * s * pow(t, 2) -
                           4 * pow(t, 3)) -
                  2 * pow(pion_mass, 4) *
                      ((2 - 5 * delta) * pow(m_rho, 2) -
                       4 * C4 * pow(m_rho, 4) + delta * (s + 6 * t)) +
                  2 * pow(m_rho, 4) *
                      (s - delta * s + t * (2 - delta + 4 * C4 * t)) +
                  pow(m_rho, 2) *
                      (8 * delta * s * t + 2 * (-2 + 3 * delta) * pow(t, 2) +
                       pow(s, 2) * (-2 + delta + 8 * C4 * t)) -
                  2 * pow(pion_mass, 2) *
                      (-8 * C4 * pow(m_rho, 6) +
                       2 * delta * (s - 3 * t) * (s + t) -
                       pow(m_rho, 2) * ((2 + delta) * s + (4 - 8 * delta) * t) +
                       pow(m_rho, 4) * (4 + 8 * C4 * (s + t)))))) /
            (pow(m_rho, 2) * (-pow(a1_mass, 2) + t)) +
        (0.0625 * pow(-2 * pow(m_rho, 2) + delta * s, 2) *
         (8 * pow(pion_mass, 4) * pow(m_rho, 2) + 2 * pow(m_rho, 6) +
          pow(s, 3) + 8 * pow(m_rho, 2) * t * (s + t) -
          pow(m_rho, 4) * (7 * s + 8 * t) +
          4 * pow(pion_mass, 2) *
              (5 * pow(m_rho, 4) - pow(s, 2) - 4 * pow(m_rho, 2) * t)) *
         HeavisideTheta(-m_rho + sqrt(s))) /
            (pow(m_rho, 6) * pow(pow(m_rho, 2) - s, 2)) -
        (0.0625 * (eta1 - eta2) * (2 * pow(m_rho, 2) - delta * s) *
         (-(eta2 * (2 * pow(pion_mass, 2) + pow(m_rho, 2) - s - 2 * t) *
            (2 * pow(pion_mass, 4) +
             pow(pion_mass, 2) * (-pow(m_rho, 2) + s - 4 * t) +
             t * (3 * pow(m_rho, 2) + s + 2 * t))) +
          eta1 * (4 * pow(pion_mass, 6) - pow(m_rho, 4) * s + pow(s, 3) +
                  2 * pow(pion_mass, 4) * (5 * pow(m_rho, 2) - s - 6 * t) -
                  2 * (pow(m_rho, 4) - 4 * pow(m_rho, 2) * s + pow(s, 2)) * t +
                  6 * (pow(m_rho, 2) - s) * pow(t, 2) - 4 * pow(t, 3) -
                  4 * pow(pion_mass, 2) *
                      (4 * pow(m_rho, 2) * t + (s - 3 * t) * (s + t)))) *
         HeavisideTheta(-m_rho + sqrt(s))) /
            (pow(m_rho, 2) * (pow(m_rho, 2) - s) * (pow(a1_mass, 2) - t)) -
        (3. * (1. * pow(m_rho, 2) - 0.5 * delta * s) *
         (delta * (0.666667 * pow(pion_mass, 4) * pow(m_rho, 2) +
                   0.166667 * pow(m_rho, 6) +
                   pow(pion_mass, 2) *
                       (1.66667 * pow(m_rho, 4) - 0.416667 * pow(s, 2) +
                        pow(m_rho, 2) * (0.0833333 * s - 1.33333 * t)) +
                   pow(m_rho, 4) * (-0.541667 * s - 0.666667 * t) +
                   pow(s, 2) * (0.125 * s + 0.0833333 * t) +
                   pow(m_rho, 2) * (-0.0833333 * pow(s, 2) + 0.583333 * s * t +
                                    0.666667 * pow(t, 2))) +
          pow(m_rho, 2) *
              (1. * C4 * pow(m_rho, 6) +
               pow(pion_mass, 2) *
                   (2. * C4 * pow(m_rho, 4) + 0.166667 * s +
                    pow(m_rho, 2) * (-0.833333 - 0.666667 * C4 * s)) +
               s * (-0.0833333 * s - 0.166667 * t) +
               pow(m_rho, 4) * (-0.416667 - 1.33333 * C4 * s - 2. * C4 * t) +
               pow(m_rho, 2) * (0.833333 * t + s * (0.5 + 0.333333 * C4 * s +
                                                    0.666667 * C4 * t)))) *
         HeavisideTheta(-m_rho + sqrt(s))) /
            (pow(m_rho, 8) - 1. * pow(m_rho, 6) * s) +
        ((0.125 * (-2 + delta) *
          (pow(pion_mass, 4) * ((-2 + 4 * delta) * pow(m_rho, 2) +
                                8 * C4 * pow(m_rho, 4) + 5 * delta * s) -
           8 * C4 * pow(m_rho, 6) * t + delta * s * t * (s + t) +
           pow(m_rho, 2) * (delta * s * (s - 3 * t) - 2 * t * (s + t)) +
           2 * pow(m_rho, 4) *
               ((-1 + delta) * s + t + 4 * C4 * t * (2 * s + t)) -
           pow(pion_mass, 2) *
               (8 * C4 * pow(m_rho, 6) + delta * s * (s + 6 * t) +
                2 * pow(m_rho, 4) * (-3 + 8 * C4 * t) +
                pow(m_rho, 2) *
                    ((-2 + 9 * delta) * s + 4 * (-1 + delta) * t)))) /
             (-pow(pion_mass, 2) + t) -
         (0.125 * (-2. + delta) * (-2. * pow(m_rho, 2) + delta * s) *
          (pow(pion_mass, 4) * (4. * pow(m_rho, 2) + 4. * s) +
           pow(pion_mass, 2) *
               (pow(m_rho, 2) * (-7. * s - 4. * t) + s * (-1. * s - 4. * t)) +
           s * (pow(m_rho, 4) + pow(m_rho, 2) * (s - 1. * t) + s * t)) *
          HeavisideTheta(-m_rho + sqrt(s))) /
             ((pow(m_rho, 2) - 1. * s) * (pow(pion_mass, 2) - 1. * t))) /
            pow(m_rho, 4))) /
      (16. * Pi * s * (-4 * pow(pion_mass, 2) + s));

  // clang-format on
  return cut_off(gev2_mb * diff_xs / spin_deg_factor);
}

}  //  namespace smash
