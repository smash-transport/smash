/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/photoncrosssections.h"

typedef double (*Fun2D)(double, double);
typedef double (*Fun3D)(double, double, double);

template <Fun2D F>
double xs_wrapper(const double s, const double m) {
  if (std::sqrt(s) > m)
    return F(s, m);
  else
    return 0;
}

template <Fun3D F>
double xs_wrapper(const double s, const double t, const double m) {
  if (std::sqrt(s) > m)
    return F(s, t, m);
  else
    return 0;
}

/*
   The cross sections presented in this file are calculated applying an average
   over initial states and sum over final states. In transport simulations
         individual particles are propagated which have a specific degeneracy
   state such that an average over initial states is superflous. The cross
   sections need thus be divided by the spin degeneracy factor of the initial
   particles to account for these particle properties.
*/

using namespace Smash;

template class PhotonCrossSection<ComputationMethod::Analytic>;
template class PhotonCrossSection<ComputationMethod::Lookup>;

constexpr double PhotonCrossSection<ComputationMethod::Analytic>::m_pion_;

/*----------------------------------------------------------------------------*/
/*				 Pi + Rho -> Pi + Photon channels mediated by
 * (Pi, Rho, a1) 				*/
/*----------------------------------------------------------------------------*/

// C11
double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_rho0_pi(
    const double s, const double mrho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  // const double m_pi = m_pion_;
  const double &mpion = m_pion_;
  // const double &mrho = m_rho;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, mrho, m_pion_, 0.);
  const double &tmax = t_mandelstam[0];
  const double &tmin = t_mandelstam[1];
  const double &t2 = tmax;
  const double &t1 = tmin;
  const double spin_deg_factor = 3.0;

  const double xs =
      (pow(Const, 2) * pow(ghat, 4) *
       ((pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(ma1, 8) + pow(mpion, 8) - pow(mpion, 4) * pow(mrho, 4) -
               2 * pow(ma1, 2) * pow(mpion, 2) *
                   (pow(mpion, 2) - pow(mrho, 2)) * (pow(mrho, 2) + s) +
               pow(ma1, 6) * (-4 * pow(mpion, 2) + 2 * s) +
               pow(ma1, 4) * (4 * pow(mpion, 4) - pow(mrho, 4) +
                              2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(ma1, 8) +
               pow(mpion, 4) * pow(pow(mpion, 2) - pow(mrho, 2), 2) +
               2 * pow(ma1, 6) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s) +
               2 * pow(ma1, 2) * pow(mpion, 2) *
                   (-pow(mrho, 4) + pow(mpion, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 4) *
                   (4 * pow(mpion, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
          pow(eta1, 2) *
              (pow(ma1, 8) + pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) -
               2 * pow(ma1, 6) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
               2 * pow(mpion, 2) * pow(mrho, 4) * s +
               pow(mpion, 4) * (3 * pow(mrho, 4) + 2 * pow(mrho, 2) * s) +
               pow(ma1, 4) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                              pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2)) -
               2 * pow(ma1, 2) *
                   (pow(mrho, 2) * s * (-pow(mrho, 2) + s) +
                    pow(mpion, 4) * (3 * pow(mrho, 2) + s) +
                    pow(mpion, 2) *
                        (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s))))) /
            ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(ma1, 2) - tmax)) +
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (4 * pow(mpion, 2) - pow(mrho, 2))) /
            ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(mpion, 2) - tmax)) -
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) * tmax) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - s, 2)) -
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) * tmax) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) *
         (-8 * C4 * pow(mrho, 4) +
          pow(mpion, 2) * (2 + delta - 8 * C4 * pow(mrho, 2)) -
          (2 + 3 * delta) * s + pow(mrho, 2) * (-2 + 3 * delta + 16 * C4 * s)) *
         tmax) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(mpion, 2) + pow(mrho, 2) - 2 * s) * (pow(mpion, 2) + s) +
          eta1 * (-2 * pow(mpion, 4) + pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                  2 * pow(s, 2) + pow(mpion, 2) * (pow(mrho, 2) + s))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) + pow(mrho, 4) +
               pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
               4 * pow(ma1, 2) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
               4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          pow(eta2, 2) *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) + pow(mrho, 4) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) -
               4 * pow(mpion, 2) * (pow(mrho, 2) + s) +
               4 * pow(ma1, 2) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s)) -
          2 * eta1 * eta2 *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) - pow(mrho, 4) +
               2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) +
               pow(ma1, 2) * (-8 * pow(mpion, 2) + 4 * s))) *
         tmax) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) +
        (8 *
         (pow(delta, 2) * (8 * pow(mpion, 4) + 3 * pow(mrho, 4) +
                           4 * pow(mpion, 2) * (3 * pow(mrho, 2) - 2 * s) -
                           6 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          4 * pow(mrho, 4) *
              (3 + 12 * C4 * (2 * pow(mpion, 2) - s) +
               8 * pow(C4, 2) * pow(-2 * pow(mpion, 2) + s, 2)) -
          4 * delta * pow(mrho, 2) *
              (16 * C4 * pow(mpion, 4) +
               2 * pow(mpion, 2) * (3 + 6 * C4 * pow(mrho, 2) - 8 * C4 * s) +
               pow(mrho, 2) * (3 - 6 * C4 * s) + s * (-3 + 4 * C4 * s))) *
         tmax) /
            (pow(mrho, 4) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (pow(mpion, 4) * (-2 + 3 * delta - 8 * C4 * pow(mrho, 2)) +
          (pow(mrho, 2) - s) * ((-2 + 3 * delta) * s +
                                pow(mrho, 2) * (-2 + delta - 8 * C4 * s)) +
          4 * pow(mpion, 2) *
              (2 * C4 * pow(mrho, 4) + delta * s -
               pow(mrho, 2) * (-1 + delta + 4 * C4 * s))) *
         tmax) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(mpion, 2) + s) *
              (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
               (pow(mrho, 2) - s) * s) +
          eta1 * (-4 * pow(mpion, 6) + pow(pow(mrho, 2) - s, 2) * s +
                  pow(mpion, 4) * (3 * pow(mrho, 2) + s) -
                  pow(mpion, 2) *
                      (pow(mrho, 4) - pow(mrho, 2) * s + 2 * pow(s, 2)))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) - pow(mrho, 4) * pow(s, 2) + pow(s, 4) -
               pow(mpion, 4) *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s - 4 * pow(s, 2)) +
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * s - 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(s, 2) * pow(pow(mrho, 2) + s, 2) +
               pow(mpion, 4) * pow(pow(mrho, 2) + 2 * s, 2) -
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta1, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) -
               4 * pow(mpion, 2) * pow(pow(mrho, 2) - s, 2) * s +
               pow(pow(mrho, 2) - s, 2) * pow(s, 2) +
               pow(mpion, 4) *
                   (3 * pow(mrho, 4) - 6 * pow(mrho, 2) * s + 4 * pow(s, 2)))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta1, 2) *
              (pow(ma1, 4) * s + pow(mpion, 4) * (-3 * pow(mrho, 2) + 2 * s) +
               s * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s + pow(s, 2)) -
               2 * pow(mpion, 2) *
                   (pow(mrho, 4) - 4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
               pow(ma1, 2) * (2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
                              3 * s * (-pow(mrho, 2) + s))) -
          2 * eta1 * eta2 *
              (pow(ma1, 4) * s +
               s * (2 * pow(mpion, 4) + 4 * pow(mpion, 2) * (pow(mrho, 2) - s) +
                    s * (-2 * pow(mrho, 2) + s)) +
               pow(ma1, 2) * (pow(mpion, 2) * (pow(mrho, 2) - 4 * s) +
                              s * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (-4 * pow(mpion, 2) * s * (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(mpion, 4) * (pow(mrho, 2) + 2 * s) +
               s * (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s)))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (-4 * pow(mpion, 4) *
                      (6 * C4 * pow(mrho, 4) + 2 * delta * s +
                       pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * s)) +
                  2 * pow(mpion, 2) *
                      (4 * delta * pow(s, 2) +
                       pow(mrho, 2) * s * (6 - 7 * delta - 16 * C4 * s) +
                       2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * s)) -
                  (pow(mrho, 2) - s) * s *
                      (-2 * delta * s +
                       pow(mrho, 2) * (-6 + 3 * delta + 8 * C4 * s))) +
          eta2 * (delta * (2 * pow(mpion, 4) * (pow(mrho, 2) + 4 * s) +
                           pow(mpion, 2) * (2 * pow(mrho, 4) +
                                            pow(mrho, 2) * s - 8 * pow(s, 2)) +
                           s * (-2 * pow(mrho, 4) - pow(mrho, 2) * s +
                                2 * pow(s, 2))) -
                  2 * pow(mrho, 2) *
                      (4 * C4 * pow(mpion, 4) * (pow(mrho, 2) + 4 * s) +
                       pow(mpion, 2) * (s * (5 - 16 * C4 * s) +
                                        pow(mrho, 2) * (2 - 8 * C4 * s)) +
                       s * (s * (-3 + 4 * C4 * s) +
                            pow(mrho, 2) * (-2 + 4 * C4 * s))))) *
         tmax) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (eta1 - eta2) *
         (delta *
              (eta1 *
                   (4 * pow(mpion, 6) +
                    pow(mpion, 4) * (7 * pow(mrho, 2) - 8 * s) +
                    pow(ma1, 4) * (pow(mpion, 2) - s) -
                    pow(ma1, 2) * (2 * pow(mpion, 2) + pow(mrho, 2) - 2 * s) *
                        (2 * pow(mpion, 2) - s) +
                    pow(mpion, 2) * s * (-8 * pow(mrho, 2) + 5 * s) +
                    s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
               eta2 *
                   (-4 * pow(mpion, 6) -
                    pow(mpion, 4) * (pow(mrho, 2) - 8 * s) +
                    pow(ma1, 4) * (-pow(mpion, 2) + s) +
                    pow(mpion, 2) * (2 * pow(mrho, 4) - 5 * pow(s, 2)) +
                    s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s + pow(s, 2)) +
                    pow(ma1, 2) * (4 * pow(mpion, 4) - 6 * pow(mpion, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)))) -
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(mpion, 6) +
                    pow(mpion, 4) * (3 + 8 * C4 * (pow(mrho, 2) - 2 * s)) +
                    2 * C4 * pow(ma1, 4) * (pow(mpion, 2) - s) +
                    2 * pow(mpion, 2) * s *
                        (-1 - 6 * C4 * pow(mrho, 2) + 5 * C4 * s) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(mpion, 4) +
                         pow(mpion, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 6 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) +
                    s * (-(s * (1 + 2 * C4 * s)) +
                         pow(mrho, 2) * (1 + 4 * C4 * s))) +
               eta2 * (2 * C4 * pow(ma1, 4) * (-pow(mpion, 2) + s) -
                       (pow(mpion, 2) - s) *
                           (8 * C4 * pow(mpion, 4) - 2 * pow(mrho, 2) + s +
                            2 * C4 * pow(s, 2) +
                            pow(mpion, 2) *
                                (3 - 4 * C4 * (pow(mrho, 2) + 2 * s))) +
                       pow(ma1, 2) *
                           (8 * C4 * pow(mpion, 4) +
                            2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                            pow(mpion, 2) *
                                (1 - 2 * C4 * (pow(mrho, 2) + 6 * s)))))) *
         tmax) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * pow(tmax, 2)) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (16 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * s *
         pow(tmax, 2)) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) * (pow(mrho, 2) - s) + 2 * eta1 * eta2 * s -
          pow(eta2, 2) * s) *
         (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
          (pow(mrho, 2) - s) * s) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (-(eta1 * (pow(mpion, 2) + 2 * pow(mrho, 2) - 3 * s)) -
          eta2 * (pow(mpion, 2) + s)) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (pow(mpion, 2) + s) *
         (-2 * eta2 * s + eta1 * (pow(mpion, 2) - pow(mrho, 2) + s)) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 3) *
         (eta1 * (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s) -
          eta2 * (pow(ma1, 2) - 2 * pow(mpion, 2) + pow(mrho, 2) + s)) *
         pow(tmax, 2)) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) -
        (8 * (delta - 4 * C4 * pow(mrho, 2)) *
         (delta * (4 * pow(mpion, 2) + 3 * pow(mrho, 2) - 2 * s) -
          2 * pow(mrho, 2) * (3 + 8 * C4 * pow(mpion, 2) - 4 * C4 * s)) *
         pow(tmax, 2)) /
            (pow(mrho, 4) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta2, 2) * s *
              (pow(ma1, 2) - 4 * pow(mpion, 2) + pow(mrho, 2) + 3 * s) +
          pow(eta1, 2) * (2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
                          s * (pow(ma1, 2) - 3 * pow(mrho, 2) + 3 * s)) -
          2 * eta1 * eta2 *
              (pow(mpion, 2) * (pow(mrho, 2) - 4 * s) +
               s * (pow(ma1, 2) - 2 * pow(mrho, 2) + 3 * s))) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (4 * delta * pow(s, 2) -
                  2 * pow(mrho, 2) * s * (-2 + 3 * delta + 8 * C4 * s) +
                  pow(mrho, 4) * (-2 + delta + 16 * C4 * s) -
                  2 * pow(mpion, 2) *
                      (8 * C4 * pow(mrho, 4) + 4 * delta * s +
                       pow(mrho, 2) * (2 - 3 * delta - 16 * C4 * s))) +
          eta2 * (pow(mpion, 2) * (8 * delta * s +
                                   pow(mrho, 2) * (-2 + delta - 32 * C4 * s)) +
                  s * (-4 * delta * s +
                       pow(mrho, 2) * (-2 + delta + 16 * C4 * s)))) *
         pow(tmax, 2)) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) *
         (delta * (eta1 * (pow(ma1, 2) * (pow(mpion, 2) - s) -
                           (2 * pow(mpion, 2) + pow(mrho, 2) - 2 * s) *
                               (2 * pow(mpion, 2) - s)) +
                   eta2 * (4 * pow(mpion, 4) - 6 * pow(mpion, 2) * s +
                           pow(ma1, 2) * (-pow(mpion, 2) + s) +
                           s * (pow(mrho, 2) + 2 * s))) +
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(mpion, 4) +
                    2 * C4 * s * (pow(ma1, 2) - pow(mrho, 2) + 2 * s) +
                    pow(mpion, 2) *
                        (1 - 2 * C4 * (pow(ma1, 2) - pow(mrho, 2) + 6 * s))) -
               eta2 * (8 * C4 * pow(mpion, 4) +
                       2 * C4 * s * (pow(ma1, 2) + pow(mrho, 2) + 2 * s) -
                       pow(mpion, 2) *
                           (-1 +
                            2 * C4 * (pow(ma1, 2) + pow(mrho, 2) + 6 * s))))) *
         pow(tmax, 2)) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 4) * pow(tmax, 3)) /
            (3. * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                   2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * pow(eta1 - eta2, 2) * (delta - 4 * C4 * pow(mrho, 2)) *
         pow(tmax, 3)) /
            (3. * pow(mrho, 2) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (16 * pow(delta - 4 * C4 * pow(mrho, 2), 2) * pow(tmax, 3)) /
            (3. * pow(mrho, 4) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (4 * (-2 + delta) * eta1 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         pow(tmax, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 4) * (pow(ma1, 2) - s) * s * pow(tmax, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) * s *
         (-2 * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-pow(mrho, 2) + s)) *
         pow(tmax, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (2 * eta2 * (delta - 4 * C4 * pow(mrho, 2)) * s +
          eta1 * (-2 * delta * s + pow(mrho, 2) * (-2 + delta + 8 * C4 * s))) *
         pow(tmax, 3)) /
            (3. * pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(ma1, 8) + pow(mpion, 8) - pow(mpion, 4) * pow(mrho, 4) -
               2 * pow(ma1, 2) * pow(mpion, 2) *
                   (pow(mpion, 2) - pow(mrho, 2)) * (pow(mrho, 2) + s) +
               pow(ma1, 6) * (-4 * pow(mpion, 2) + 2 * s) +
               pow(ma1, 4) * (4 * pow(mpion, 4) - pow(mrho, 4) +
                              2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(ma1, 8) +
               pow(mpion, 4) * pow(pow(mpion, 2) - pow(mrho, 2), 2) +
               2 * pow(ma1, 6) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s) +
               2 * pow(ma1, 2) * pow(mpion, 2) *
                   (-pow(mrho, 4) + pow(mpion, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 4) *
                   (4 * pow(mpion, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
          pow(eta1, 2) *
              (pow(ma1, 8) + pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) -
               2 * pow(ma1, 6) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
               2 * pow(mpion, 2) * pow(mrho, 4) * s +
               pow(mpion, 4) * (3 * pow(mrho, 4) + 2 * pow(mrho, 2) * s) +
               pow(ma1, 4) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                              pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2)) -
               2 * pow(ma1, 2) *
                   (pow(mrho, 2) * s * (-pow(mrho, 2) + s) +
                    pow(mpion, 4) * (3 * pow(mrho, 2) + s) +
                    pow(mpion, 2) *
                        (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s))))) /
            ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(ma1, 2) - tmin)) -
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (4 * pow(mpion, 2) - pow(mrho, 2))) /
            ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(mpion, 2) - tmin)) +
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) * tmin) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - s, 2)) +
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) * tmin) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (-8 * C4 * pow(mrho, 4) +
          pow(mpion, 2) * (2 + delta - 8 * C4 * pow(mrho, 2)) -
          (2 + 3 * delta) * s + pow(mrho, 2) * (-2 + 3 * delta + 16 * C4 * s)) *
         tmin) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(mpion, 2) + pow(mrho, 2) - 2 * s) * (pow(mpion, 2) + s) +
          eta1 * (-2 * pow(mpion, 4) + pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                  2 * pow(s, 2) + pow(mpion, 2) * (pow(mrho, 2) + s))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) + pow(mrho, 4) +
               pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
               4 * pow(ma1, 2) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
               4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          pow(eta2, 2) *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) + pow(mrho, 4) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) -
               4 * pow(mpion, 2) * (pow(mrho, 2) + s) +
               4 * pow(ma1, 2) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s)) -
          2 * eta1 * eta2 *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) - pow(mrho, 4) +
               2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) +
               pow(ma1, 2) * (-8 * pow(mpion, 2) + 4 * s))) *
         tmin) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) -
        (8 *
         (pow(delta, 2) * (8 * pow(mpion, 4) + 3 * pow(mrho, 4) +
                           4 * pow(mpion, 2) * (3 * pow(mrho, 2) - 2 * s) -
                           6 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          4 * pow(mrho, 4) *
              (3 + 12 * C4 * (2 * pow(mpion, 2) - s) +
               8 * pow(C4, 2) * pow(-2 * pow(mpion, 2) + s, 2)) -
          4 * delta * pow(mrho, 2) *
              (16 * C4 * pow(mpion, 4) +
               2 * pow(mpion, 2) * (3 + 6 * C4 * pow(mrho, 2) - 8 * C4 * s) +
               pow(mrho, 2) * (3 - 6 * C4 * s) + s * (-3 + 4 * C4 * s))) *
         tmin) /
            (pow(mrho, 4) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) *
         (pow(mpion, 4) * (-2 + 3 * delta - 8 * C4 * pow(mrho, 2)) +
          (pow(mrho, 2) - s) * ((-2 + 3 * delta) * s +
                                pow(mrho, 2) * (-2 + delta - 8 * C4 * s)) +
          4 * pow(mpion, 2) *
              (2 * C4 * pow(mrho, 4) + delta * s -
               pow(mrho, 2) * (-1 + delta + 4 * C4 * s))) *
         tmin) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(mpion, 2) + s) *
              (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
               (pow(mrho, 2) - s) * s) +
          eta1 * (-4 * pow(mpion, 6) + pow(pow(mrho, 2) - s, 2) * s +
                  pow(mpion, 4) * (3 * pow(mrho, 2) + s) -
                  pow(mpion, 2) *
                      (pow(mrho, 4) - pow(mrho, 2) * s + 2 * pow(s, 2)))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) - pow(mrho, 4) * pow(s, 2) + pow(s, 4) -
               pow(mpion, 4) *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s - 4 * pow(s, 2)) +
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * s - 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(s, 2) * pow(pow(mrho, 2) + s, 2) +
               pow(mpion, 4) * pow(pow(mrho, 2) + 2 * s, 2) -
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta1, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) -
               4 * pow(mpion, 2) * pow(pow(mrho, 2) - s, 2) * s +
               pow(pow(mrho, 2) - s, 2) * pow(s, 2) +
               pow(mpion, 4) *
                   (3 * pow(mrho, 4) - 6 * pow(mrho, 2) * s + 4 * pow(s, 2)))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta1, 2) *
              (pow(ma1, 4) * s + pow(mpion, 4) * (-3 * pow(mrho, 2) + 2 * s) +
               s * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s + pow(s, 2)) -
               2 * pow(mpion, 2) *
                   (pow(mrho, 4) - 4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
               pow(ma1, 2) * (2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
                              3 * s * (-pow(mrho, 2) + s))) -
          2 * eta1 * eta2 *
              (pow(ma1, 4) * s +
               s * (2 * pow(mpion, 4) + 4 * pow(mpion, 2) * (pow(mrho, 2) - s) +
                    s * (-2 * pow(mrho, 2) + s)) +
               pow(ma1, 2) * (pow(mpion, 2) * (pow(mrho, 2) - 4 * s) +
                              s * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (-4 * pow(mpion, 2) * s * (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(mpion, 4) * (pow(mrho, 2) + 2 * s) +
               s * (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s)))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (4 * pow(mpion, 4) *
                      (6 * C4 * pow(mrho, 4) + 2 * delta * s +
                       pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * s)) -
                  2 * pow(mpion, 2) *
                      (4 * delta * pow(s, 2) +
                       pow(mrho, 2) * s * (6 - 7 * delta - 16 * C4 * s) +
                       2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * s)) +
                  (pow(mrho, 2) - s) * s *
                      (-2 * delta * s +
                       pow(mrho, 2) * (-6 + 3 * delta + 8 * C4 * s))) +
          eta2 *
              (-(delta *
                 (2 * pow(mpion, 4) * (pow(mrho, 2) + 4 * s) +
                  pow(mpion, 2) *
                      (2 * pow(mrho, 4) + pow(mrho, 2) * s - 8 * pow(s, 2)) +
                  s * (-2 * pow(mrho, 4) - pow(mrho, 2) * s + 2 * pow(s, 2)))) +
               2 * pow(mrho, 2) *
                   (4 * C4 * pow(mpion, 4) * (pow(mrho, 2) + 4 * s) +
                    pow(mpion, 2) * (s * (5 - 16 * C4 * s) +
                                     pow(mrho, 2) * (2 - 8 * C4 * s)) +
                    s * (s * (-3 + 4 * C4 * s) +
                         pow(mrho, 2) * (-2 + 4 * C4 * s))))) *
         tmin) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (eta1 - eta2) *
         (delta *
              (eta1 *
                   (4 * pow(mpion, 6) +
                    pow(mpion, 4) * (7 * pow(mrho, 2) - 8 * s) +
                    pow(ma1, 4) * (pow(mpion, 2) - s) -
                    pow(ma1, 2) * (2 * pow(mpion, 2) + pow(mrho, 2) - 2 * s) *
                        (2 * pow(mpion, 2) - s) +
                    pow(mpion, 2) * s * (-8 * pow(mrho, 2) + 5 * s) +
                    s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
               eta2 *
                   (-4 * pow(mpion, 6) -
                    pow(mpion, 4) * (pow(mrho, 2) - 8 * s) +
                    pow(ma1, 4) * (-pow(mpion, 2) + s) +
                    pow(mpion, 2) * (2 * pow(mrho, 4) - 5 * pow(s, 2)) +
                    s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s + pow(s, 2)) +
                    pow(ma1, 2) * (4 * pow(mpion, 4) - 6 * pow(mpion, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)))) -
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(mpion, 6) +
                    pow(mpion, 4) * (3 + 8 * C4 * (pow(mrho, 2) - 2 * s)) +
                    2 * C4 * pow(ma1, 4) * (pow(mpion, 2) - s) +
                    2 * pow(mpion, 2) * s *
                        (-1 - 6 * C4 * pow(mrho, 2) + 5 * C4 * s) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(mpion, 4) +
                         pow(mpion, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 6 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) +
                    s * (-(s * (1 + 2 * C4 * s)) +
                         pow(mrho, 2) * (1 + 4 * C4 * s))) +
               eta2 * (2 * C4 * pow(ma1, 4) * (-pow(mpion, 2) + s) -
                       (pow(mpion, 2) - s) *
                           (8 * C4 * pow(mpion, 4) - 2 * pow(mrho, 2) + s +
                            2 * C4 * pow(s, 2) +
                            pow(mpion, 2) *
                                (3 - 4 * C4 * (pow(mrho, 2) + 2 * s))) +
                       pow(ma1, 2) *
                           (8 * C4 * pow(mpion, 4) +
                            2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                            pow(mpion, 2) *
                                (1 - 2 * C4 * (pow(mrho, 2) + 6 * s)))))) *
         tmin) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * pow(tmin, 2)) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (16 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * s *
         pow(tmin, 2)) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) * (pow(mrho, 2) - s) + 2 * eta1 * eta2 * s -
          pow(eta2, 2) * s) *
         (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
          (pow(mrho, 2) - s) * s) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (-(eta1 * (pow(mpion, 2) + 2 * pow(mrho, 2) - 3 * s)) -
          eta2 * (pow(mpion, 2) + s)) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (pow(mpion, 2) + s) *
         (-2 * eta2 * s + eta1 * (pow(mpion, 2) - pow(mrho, 2) + s)) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 3) *
         (-(eta1 * (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s)) +
          eta2 * (pow(ma1, 2) - 2 * pow(mpion, 2) + pow(mrho, 2) + s)) *
         pow(tmin, 2)) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) +
        (8 * (delta - 4 * C4 * pow(mrho, 2)) *
         (delta * (4 * pow(mpion, 2) + 3 * pow(mrho, 2) - 2 * s) -
          2 * pow(mrho, 2) * (3 + 8 * C4 * pow(mpion, 2) - 4 * C4 * s)) *
         pow(tmin, 2)) /
            (pow(mrho, 4) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta2, 2) * s *
              (pow(ma1, 2) - 4 * pow(mpion, 2) + pow(mrho, 2) + 3 * s) +
          pow(eta1, 2) * (2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
                          s * (pow(ma1, 2) - 3 * pow(mrho, 2) + 3 * s)) -
          2 * eta1 * eta2 *
              (pow(mpion, 2) * (pow(mrho, 2) - 4 * s) +
               s * (pow(ma1, 2) - 2 * pow(mrho, 2) + 3 * s))) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (4 * delta * pow(s, 2) -
                  2 * pow(mrho, 2) * s * (-2 + 3 * delta + 8 * C4 * s) +
                  pow(mrho, 4) * (-2 + delta + 16 * C4 * s) -
                  2 * pow(mpion, 2) *
                      (8 * C4 * pow(mrho, 4) + 4 * delta * s +
                       pow(mrho, 2) * (2 - 3 * delta - 16 * C4 * s))) +
          eta2 * (pow(mpion, 2) * (8 * delta * s +
                                   pow(mrho, 2) * (-2 + delta - 32 * C4 * s)) +
                  s * (-4 * delta * s +
                       pow(mrho, 2) * (-2 + delta + 16 * C4 * s)))) *
         pow(tmin, 2)) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (4 * (eta1 - eta2) *
         (delta * (eta1 * (pow(ma1, 2) * (pow(mpion, 2) - s) -
                           (2 * pow(mpion, 2) + pow(mrho, 2) - 2 * s) *
                               (2 * pow(mpion, 2) - s)) +
                   eta2 * (4 * pow(mpion, 4) - 6 * pow(mpion, 2) * s +
                           pow(ma1, 2) * (-pow(mpion, 2) + s) +
                           s * (pow(mrho, 2) + 2 * s))) +
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(mpion, 4) +
                    2 * C4 * s * (pow(ma1, 2) - pow(mrho, 2) + 2 * s) +
                    pow(mpion, 2) *
                        (1 - 2 * C4 * (pow(ma1, 2) - pow(mrho, 2) + 6 * s))) -
               eta2 * (8 * C4 * pow(mpion, 4) +
                       2 * C4 * s * (pow(ma1, 2) + pow(mrho, 2) + 2 * s) -
                       pow(mpion, 2) *
                           (-1 +
                            2 * C4 * (pow(ma1, 2) + pow(mrho, 2) + 6 * s))))) *
         pow(tmin, 2)) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 4) * pow(tmin, 3)) /
            (3. * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                   2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * pow(eta1 - eta2, 2) * (delta - 4 * C4 * pow(mrho, 2)) *
         pow(tmin, 3)) /
            (3. * pow(mrho, 2) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (16 * pow(delta - 4 * C4 * pow(mrho, 2), 2) * pow(tmin, 3)) /
            (3. * pow(mrho, 4) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (-2 + delta) * eta1 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         pow(tmin, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 4) * (pow(ma1, 2) - s) * s * pow(tmin, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) * s *
         (-2 * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-pow(mrho, 2) + s)) *
         pow(tmin, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (2 * eta2 * (delta - 4 * C4 * pow(mrho, 2)) * s +
          eta1 * (-2 * delta * s + pow(mrho, 2) * (-2 + delta + 8 * C4 * s))) *
         pow(tmin, 3)) /
            (3. * pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (2 * pow(ma1, 6) -
               3 * pow(ma1, 4) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) +
               pow(mrho, 2) * (pow(mrho, 2) - s) * s -
               pow(mpion, 4) * (3 * pow(mrho, 2) + s) +
               pow(mpion, 2) * (-2 * pow(mrho, 4) + 3 * pow(mrho, 2) * s) +
               pow(ma1, 2) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                              pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2))) -
          2 * eta1 * eta2 *
              (2 * pow(ma1, 6) -
               pow(mpion, 2) * (pow(mpion, 2) - pow(mrho, 2)) *
                   (pow(mrho, 2) + s) +
               pow(ma1, 4) * (-6 * pow(mpion, 2) + 3 * s) +
               pow(ma1, 2) * (4 * pow(mpion, 4) - pow(mrho, 4) +
                              2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (2 * pow(ma1, 6) +
               3 * pow(ma1, 4) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s) +
               pow(mpion, 2) *
                   (-pow(mrho, 4) + pow(mpion, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 2) *
                   (4 * pow(mpion, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(mpion, 2) * (pow(mrho, 2) + s)))) *
         log(abs(-pow(ma1, 2) + tmax))) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) -
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               2 * pow(ma1, 2) * pow(mpion, 4) * s +
               pow(mpion, 2) * (pow(ma1, 4) * (pow(mrho, 2) - 4 * s) +
                                4 * pow(ma1, 2) * (pow(mrho, 2) - s) * s +
                                pow(mrho, 2) * pow(s, 2)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (-2 * pow(mrho, 2) + s) +
                    pow(ma1, 2) * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (pow(mpion, 8) -
               4 * pow(ma1, 2) * pow(mpion, 2) * s *
                   (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(mpion, 4) *
                   (pow(mrho, 2) * s + pow(ma1, 2) * (pow(mrho, 2) + 2 * s)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s))) +
          pow(eta1, 2) *
              (pow(mpion, 8) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + 2 * pow(mrho, 4) -
                    3 * pow(ma1, 2) * (pow(mrho, 2) - s) -
                    3 * pow(mrho, 2) * s + pow(s, 2)) +
               pow(mpion, 4) * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                                pow(ma1, 2) * (-3 * pow(mrho, 2) + 2 * s)) +
               2 * pow(mpion, 2) *
                   (pow(ma1, 4) * (pow(mrho, 2) - 2 * s) +
                    pow(mrho, 2) * s * (-pow(mrho, 2) + s) -
                    pow(ma1, 2) * (pow(mrho, 4) - 4 * pow(mrho, 2) * s +
                                   2 * pow(s, 2))))) *
         log(abs(-pow(ma1, 2) + tmax))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (eta1 - eta2) *
         (delta *
              (eta2 *
                   (pow(mpion, 6) * pow(mrho, 2) * (2 * pow(mpion, 2) - s) +
                    pow(ma1, 8) * (-pow(mpion, 2) + s) +
                    pow(ma1, 6) * (5 * pow(mpion, 4) - 7 * pow(mpion, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)) +
                    pow(ma1, 4) *
                        (-8 * pow(mpion, 6) -
                         pow(mpion, 4) * (pow(mrho, 2) - 14 * s) +
                         pow(mpion, 2) * (2 * pow(mrho, 4) - pow(mrho, 2) * s -
                                          7 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s +
                              pow(s, 2))) +
                    pow(ma1, 2) * pow(mpion, 2) *
                        (4 * pow(mpion, 6) +
                         pow(mpion, 4) * (pow(mrho, 2) - 8 * s) +
                         s * (2 * pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2)) +
                         pow(mpion, 2) *
                             (-2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                              5 * pow(s, 2)))) +
               eta1 *
                   (pow(ma1, 8) * (pow(mpion, 2) - s) +
                    pow(ma1, 6) *
                        (-5 * pow(mpion, 4) + (pow(mrho, 2) - 2 * s) * s +
                         pow(mpion, 2) * (-2 * pow(mrho, 2) + 7 * s)) +
                    pow(mpion, 2) * pow(mrho, 2) *
                        (2 * pow(mpion, 6) +
                         pow(mpion, 4) * (4 * pow(mrho, 2) - 5 * s) +
                         pow(mrho, 4) * s -
                         pow(mpion, 2) * (pow(mrho, 4) + 3 * pow(mrho, 2) * s -
                                          2 * pow(s, 2))) +
                    pow(ma1, 4) *
                        (8 * pow(mpion, 6) +
                         pow(mpion, 4) * (9 * pow(mrho, 2) - 14 * s) +
                         pow(mpion, 2) * s * (-9 * pow(mrho, 2) + 7 * s) +
                         s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
                    pow(ma1, 2) *
                        (-4 * pow(mpion, 8) +
                         pow(mrho, 4) * s * (-pow(mrho, 2) + s) +
                         pow(mpion, 6) * (-11 * pow(mrho, 2) + 8 * s) +
                         pow(mpion, 4) *
                             (-3 * pow(mrho, 4) + 17 * pow(mrho, 2) * s -
                              5 * pow(s, 2)) +
                         pow(mpion, 2) *
                             (pow(mrho, 6) - 5 * pow(mrho, 2) * pow(s, 2) +
                              pow(s, 3))))) -
          2 * pow(mrho, 2) *
              (eta2 *
                   (pow(mpion, 8) * (1 + 2 * C4 * pow(mrho, 2)) -
                    2 * C4 * pow(mpion, 6) * pow(mrho, 2) * s +
                    2 * C4 * pow(ma1, 8) * (-pow(mpion, 2) + s) +
                    pow(ma1, 4) *
                        (-16 * C4 * pow(mpion, 6) +
                         pow(mpion, 4) *
                             (-4 + 6 * C4 * pow(mrho, 2) + 28 * C4 * s) +
                         2 * pow(mpion, 2) *
                             (pow(mrho, 2) + s - 3 * C4 * pow(mrho, 2) * s -
                              7 * C4 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 2) + s + 2 * C4 * pow(s, 2))) +
                    pow(ma1, 6) * (10 * C4 * pow(mpion, 4) +
                                   2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                                   pow(mpion, 2) *
                                       (1 - 2 * C4 * (pow(mrho, 2) + 7 * s))) +
                    pow(ma1, 2) * pow(mpion, 2) *
                        (8 * C4 * pow(mpion, 6) -
                         2 * pow(mpion, 4) *
                             (-2 + 3 * C4 * pow(mrho, 2) + 8 * C4 * s) +
                         s * (2 * pow(mrho, 2) + s - 2 * C4 * pow(s, 2)) +
                         2 * pow(mpion, 2) *
                             (pow(mrho, 2) * (-1 + 3 * C4 * s) +
                              s * (-3 + 5 * C4 * s)))) +
               eta1 *
                   (pow(mpion, 8) * (-1 + 6 * C4 * pow(mrho, 2)) +
                    2 * C4 * pow(ma1, 8) * (pow(mpion, 2) - s) +
                    pow(mpion, 2) * pow(mrho, 4) * s +
                    2 * pow(mpion, 6) * pow(mrho, 2) * (2 - 5 * C4 * s) -
                    pow(ma1, 6) *
                        (10 * C4 * pow(mpion, 4) +
                         pow(mpion, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 7 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) -
                    pow(mpion, 4) * pow(mrho, 2) *
                        (pow(mrho, 2) + s * (3 - 4 * C4 * s)) +
                    pow(ma1, 4) *
                        (16 * C4 * pow(mpion, 6) +
                         2 * pow(mpion, 4) *
                             (2 + 5 * C4 * pow(mrho, 2) - 14 * C4 * s) +
                         2 * pow(mpion, 2) * s *
                             (-1 - 7 * C4 * pow(mrho, 2) + 7 * C4 * s) +
                         s * (-(s * (1 + 2 * C4 * s)) +
                              pow(mrho, 2) * (1 + 4 * C4 * s))) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(mpion, 8) +
                         pow(mrho, 2) * (pow(mrho, 2) - s) * s +
                         2 * pow(mpion, 6) *
                             (2 + 7 * C4 * pow(mrho, 2) - 8 * C4 * s) +
                         pow(mpion, 2) * (-pow(mrho, 4) + pow(s, 2) +
                                          8 * C4 * pow(mrho, 2) * pow(s, 2) -
                                          2 * C4 * pow(s, 3)) +
                         pow(mpion, 4) * (pow(mrho, 2) * (3 - 22 * C4 * s) +
                                          2 * s * (-3 + 5 * C4 * s)))))) *
         log(abs(-pow(ma1, 2) + tmax))) /
            ((pow(ma1, 2) - pow(mpion, 2)) * pow(mrho, 2) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (16 * pow(-2 + delta, 2) * pow(mpion, 2) *
         log(abs(-pow(mpion, 2) + tmax))) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) -
        (8 * pow(-2 + delta, 2) *
         (3 * pow(mpion, 4) - 4 * pow(mpion, 2) * (pow(mrho, 2) - s) +
          pow(pow(mrho, 2) - s, 2)) *
         log(abs(-pow(mpion, 2) + tmax))) /
            ((pow(mpion, 2) - s) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                                    2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) * (eta1 - eta2) * pow(mpion, 2) *
         (2 * eta1 * pow(mpion, 2) - 2 * eta2 * pow(mpion, 2) -
          eta1 * pow(mrho, 2)) *
         (pow(mpion, 2) - s) * log(abs(-pow(mpion, 2) + tmax))) /
            ((pow(ma1, 2) - pow(mpion, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) * (eta1 - eta2) * pow(mpion, 2) * (pow(ma1, 2) - s) *
         (pow(mpion, 2) - s) *
         (-(eta2 * (pow(mpion, 2) + s)) +
          eta1 * (pow(mpion, 2) - pow(mrho, 2) + s)) *
         log(abs(-pow(mpion, 2) + tmax))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (-(delta * (4 * pow(mpion, 2) - pow(mrho, 2)) *
            (pow(mpion, 2) + pow(mrho, 2) - s)) +
          2 * pow(mrho, 2) *
              (8 * C4 * pow(mpion, 4) - pow(mrho, 2) + s +
               pow(mpion, 2) * (3 - 8 * C4 * s))) *
         log(abs(-pow(mpion, 2) + tmax))) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (2 * pow(ma1, 6) -
               3 * pow(ma1, 4) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) +
               pow(mrho, 2) * (pow(mrho, 2) - s) * s -
               pow(mpion, 4) * (3 * pow(mrho, 2) + s) +
               pow(mpion, 2) * (-2 * pow(mrho, 4) + 3 * pow(mrho, 2) * s) +
               pow(ma1, 2) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                              pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2))) -
          2 * eta1 * eta2 *
              (2 * pow(ma1, 6) -
               pow(mpion, 2) * (pow(mpion, 2) - pow(mrho, 2)) *
                   (pow(mrho, 2) + s) +
               pow(ma1, 4) * (-6 * pow(mpion, 2) + 3 * s) +
               pow(ma1, 2) * (4 * pow(mpion, 4) - pow(mrho, 4) +
                              2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (2 * pow(ma1, 6) +
               3 * pow(ma1, 4) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s) +
               pow(mpion, 2) *
                   (-pow(mrho, 4) + pow(mpion, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 2) *
                   (4 * pow(mpion, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(mpion, 2) * (pow(mrho, 2) + s)))) *
         log(abs(-pow(ma1, 2) + tmin))) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) +
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               2 * pow(ma1, 2) * pow(mpion, 4) * s +
               pow(mpion, 2) * (pow(ma1, 4) * (pow(mrho, 2) - 4 * s) +
                                4 * pow(ma1, 2) * (pow(mrho, 2) - s) * s +
                                pow(mrho, 2) * pow(s, 2)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (-2 * pow(mrho, 2) + s) +
                    pow(ma1, 2) * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (pow(mpion, 8) -
               4 * pow(ma1, 2) * pow(mpion, 2) * s *
                   (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(mpion, 4) *
                   (pow(mrho, 2) * s + pow(ma1, 2) * (pow(mrho, 2) + 2 * s)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s))) +
          pow(eta1, 2) *
              (pow(mpion, 8) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + 2 * pow(mrho, 4) -
                    3 * pow(ma1, 2) * (pow(mrho, 2) - s) -
                    3 * pow(mrho, 2) * s + pow(s, 2)) +
               pow(mpion, 4) * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                                pow(ma1, 2) * (-3 * pow(mrho, 2) + 2 * s)) +
               2 * pow(mpion, 2) *
                   (pow(ma1, 4) * (pow(mrho, 2) - 2 * s) +
                    pow(mrho, 2) * s * (-pow(mrho, 2) + s) -
                    pow(ma1, 2) * (pow(mrho, 4) - 4 * pow(mrho, 2) * s +
                                   2 * pow(s, 2))))) *
         log(abs(-pow(ma1, 2) + tmin))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (eta1 - eta2) *
         (delta *
              (eta2 *
                   (pow(mpion, 6) * pow(mrho, 2) * (2 * pow(mpion, 2) - s) +
                    pow(ma1, 8) * (-pow(mpion, 2) + s) +
                    pow(ma1, 6) * (5 * pow(mpion, 4) - 7 * pow(mpion, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)) +
                    pow(ma1, 4) *
                        (-8 * pow(mpion, 6) -
                         pow(mpion, 4) * (pow(mrho, 2) - 14 * s) +
                         pow(mpion, 2) * (2 * pow(mrho, 4) - pow(mrho, 2) * s -
                                          7 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s +
                              pow(s, 2))) +
                    pow(ma1, 2) * pow(mpion, 2) *
                        (4 * pow(mpion, 6) +
                         pow(mpion, 4) * (pow(mrho, 2) - 8 * s) +
                         s * (2 * pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2)) +
                         pow(mpion, 2) *
                             (-2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                              5 * pow(s, 2)))) +
               eta1 *
                   (pow(ma1, 8) * (pow(mpion, 2) - s) +
                    pow(ma1, 6) *
                        (-5 * pow(mpion, 4) + (pow(mrho, 2) - 2 * s) * s +
                         pow(mpion, 2) * (-2 * pow(mrho, 2) + 7 * s)) +
                    pow(mpion, 2) * pow(mrho, 2) *
                        (2 * pow(mpion, 6) +
                         pow(mpion, 4) * (4 * pow(mrho, 2) - 5 * s) +
                         pow(mrho, 4) * s -
                         pow(mpion, 2) * (pow(mrho, 4) + 3 * pow(mrho, 2) * s -
                                          2 * pow(s, 2))) +
                    pow(ma1, 4) *
                        (8 * pow(mpion, 6) +
                         pow(mpion, 4) * (9 * pow(mrho, 2) - 14 * s) +
                         pow(mpion, 2) * s * (-9 * pow(mrho, 2) + 7 * s) +
                         s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
                    pow(ma1, 2) *
                        (-4 * pow(mpion, 8) +
                         pow(mrho, 4) * s * (-pow(mrho, 2) + s) +
                         pow(mpion, 6) * (-11 * pow(mrho, 2) + 8 * s) +
                         pow(mpion, 4) *
                             (-3 * pow(mrho, 4) + 17 * pow(mrho, 2) * s -
                              5 * pow(s, 2)) +
                         pow(mpion, 2) *
                             (pow(mrho, 6) - 5 * pow(mrho, 2) * pow(s, 2) +
                              pow(s, 3))))) -
          2 * pow(mrho, 2) *
              (eta2 *
                   (pow(mpion, 8) * (1 + 2 * C4 * pow(mrho, 2)) -
                    2 * C4 * pow(mpion, 6) * pow(mrho, 2) * s +
                    2 * C4 * pow(ma1, 8) * (-pow(mpion, 2) + s) +
                    pow(ma1, 4) *
                        (-16 * C4 * pow(mpion, 6) +
                         pow(mpion, 4) *
                             (-4 + 6 * C4 * pow(mrho, 2) + 28 * C4 * s) +
                         2 * pow(mpion, 2) *
                             (pow(mrho, 2) + s - 3 * C4 * pow(mrho, 2) * s -
                              7 * C4 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 2) + s + 2 * C4 * pow(s, 2))) +
                    pow(ma1, 6) * (10 * C4 * pow(mpion, 4) +
                                   2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                                   pow(mpion, 2) *
                                       (1 - 2 * C4 * (pow(mrho, 2) + 7 * s))) +
                    pow(ma1, 2) * pow(mpion, 2) *
                        (8 * C4 * pow(mpion, 6) -
                         2 * pow(mpion, 4) *
                             (-2 + 3 * C4 * pow(mrho, 2) + 8 * C4 * s) +
                         s * (2 * pow(mrho, 2) + s - 2 * C4 * pow(s, 2)) +
                         2 * pow(mpion, 2) *
                             (pow(mrho, 2) * (-1 + 3 * C4 * s) +
                              s * (-3 + 5 * C4 * s)))) +
               eta1 *
                   (pow(mpion, 8) * (-1 + 6 * C4 * pow(mrho, 2)) +
                    2 * C4 * pow(ma1, 8) * (pow(mpion, 2) - s) +
                    pow(mpion, 2) * pow(mrho, 4) * s +
                    2 * pow(mpion, 6) * pow(mrho, 2) * (2 - 5 * C4 * s) -
                    pow(ma1, 6) *
                        (10 * C4 * pow(mpion, 4) +
                         pow(mpion, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 7 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) -
                    pow(mpion, 4) * pow(mrho, 2) *
                        (pow(mrho, 2) + s * (3 - 4 * C4 * s)) +
                    pow(ma1, 4) *
                        (16 * C4 * pow(mpion, 6) +
                         2 * pow(mpion, 4) *
                             (2 + 5 * C4 * pow(mrho, 2) - 14 * C4 * s) +
                         2 * pow(mpion, 2) * s *
                             (-1 - 7 * C4 * pow(mrho, 2) + 7 * C4 * s) +
                         s * (-(s * (1 + 2 * C4 * s)) +
                              pow(mrho, 2) * (1 + 4 * C4 * s))) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(mpion, 8) +
                         pow(mrho, 2) * (pow(mrho, 2) - s) * s +
                         2 * pow(mpion, 6) *
                             (2 + 7 * C4 * pow(mrho, 2) - 8 * C4 * s) +
                         pow(mpion, 2) * (-pow(mrho, 4) + pow(s, 2) +
                                          8 * C4 * pow(mrho, 2) * pow(s, 2) -
                                          2 * C4 * pow(s, 3)) +
                         pow(mpion, 4) * (pow(mrho, 2) * (3 - 22 * C4 * s) +
                                          2 * s * (-3 + 5 * C4 * s)))))) *
         log(abs(-pow(ma1, 2) + tmin))) /
            ((pow(ma1, 2) - pow(mpion, 2)) * pow(mrho, 2) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (16 * pow(-2 + delta, 2) * pow(mpion, 2) *
         log(abs(-pow(mpion, 2) + tmin))) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) +
        (8 * pow(-2 + delta, 2) *
         (3 * pow(mpion, 4) - 4 * pow(mpion, 2) * (pow(mrho, 2) - s) +
          pow(pow(mrho, 2) - s, 2)) *
         log(abs(-pow(mpion, 2) + tmin))) /
            ((pow(mpion, 2) - s) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                                    2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) * (eta1 - eta2) * pow(mpion, 2) *
         (2 * eta1 * pow(mpion, 2) - 2 * eta2 * pow(mpion, 2) -
          eta1 * pow(mrho, 2)) *
         (pow(mpion, 2) - s) * log(abs(-pow(mpion, 2) + tmin))) /
            ((pow(ma1, 2) - pow(mpion, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) * (eta1 - eta2) * pow(mpion, 2) * (pow(ma1, 2) - s) *
         (pow(mpion, 2) - s) *
         (-(eta2 * (pow(mpion, 2) + s)) +
          eta1 * (pow(mpion, 2) - pow(mrho, 2) + s)) *
         log(abs(-pow(mpion, 2) + tmin))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (delta * (4 * pow(mpion, 2) - pow(mrho, 2)) *
              (pow(mpion, 2) + pow(mrho, 2) - s) -
          2 * pow(mrho, 2) *
              (8 * C4 * pow(mpion, 4) - pow(mrho, 2) + s +
               pow(mpion, 2) * (3 - 8 * C4 * s))) *
         log(abs(-pow(mpion, 2) + tmin))) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))))) /
      (512. * Pi);

  return xs * to_mb / spin_deg_factor;
}

double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_rho0_pi(
    const double s, const double t, const double m_rho) {
  const double &mrho = m_rho;
  const double &mpion = m_pion_;
  const double spin_deg_factor = 3.0;

  const double diff_xs =
      (pow(Const, 2) * pow(ghat, 4) *
       ((-8 * pow(-2 + delta, 2) * pow(mpion, 2)) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - s, 2)) -
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (pow(mpion, 4) + pow(pow(mrho, 2) - t, 2) -
          2 * pow(mpion, 2) * (pow(mrho, 2) + t))) /
            (pow(mrho, 2) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             pow(pow(mpion, 2) - t, 2)) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (-(eta2 * (pow(mpion, 2) + s)) + eta1 * (-pow(mrho, 2) + s + t)) *
         (-pow(mpion, 4) + pow(mpion, 2) * (pow(mrho, 2) - 2 * t) +
          t * (-pow(mrho, 2) + 2 * s + t))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(mpion, 2) - t)) -
        (8 * (-2 + delta) *
         (pow(mpion, 4) * (2 - 3 * delta + 8 * C4 * pow(mrho, 2)) +
          pow(mrho, 4) * (-2 + delta + 8 * C4 * t) +
          t * ((2 + 3 * delta) * s + 2 * delta * t) +
          pow(mpion, 2) *
              (-8 * C4 * pow(mrho, 4) + (-2 + delta) * s - (2 + 3 * delta) * t +
               4 * pow(mrho, 2) * (1 + 4 * C4 * t)) -
          pow(mrho, 2) * (t * (-2 + 3 * delta + 8 * C4 * t) +
                          s * (-2 + delta + 16 * C4 * t)))) /
            (pow(mrho, 2) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(mpion, 2) - t)) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(mpion, 2) + s) *
              (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
               s * (pow(mrho, 2) - s - 2 * t)) +
          eta1 * (-4 * pow(mpion, 6) +
                  s * (-pow(mrho, 2) + s) * (-pow(mrho, 2) + s + t) +
                  pow(mpion, 4) * (3 * pow(mrho, 2) + s + t) -
                  pow(mpion, 2) * (pow(mrho, 4) + 2 * s * (s - t) +
                                   pow(mrho, 2) * (-s + t))))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) * (pow(pow(mrho, 2) + 2 * s, 2) - 2 * s * t) +
               pow(s, 2) * (pow(pow(mrho, 2) + s, 2) +
                            2 * (-pow(mrho, 2) + s) * t + 2 * pow(t, 2)) -
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * (2 * s - t) +
                    2 * s * (s + t))) +
          2 * eta1 * eta2 *
              (-pow(mpion, 8) +
               pow(mpion, 4) * (pow(mrho, 4) + 2 * pow(mrho, 2) * s +
                                2 * s * (-2 * s + t)) -
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * (s + t) - 2 * s * (s + t)) +
               pow(s, 2) * (pow(mrho, 4) - pow(s, 2) + 2 * pow(mrho, 2) * t -
                            2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) * (3 * pow(mrho, 4) + 2 * s * (2 * s - t) +
                                2 * pow(mrho, 2) * (-3 * s + t)) -
               2 * pow(mpion, 2) * (pow(mrho, 2) - s) *
                   (-2 * s * (s + t) + pow(mrho, 2) * (2 * s + t)) +
               s * (-pow(mrho, 2) + s) *
                   (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                    pow(mrho, 2) * (s + 2 * t))))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) -
               pow(mpion, 4) *
                   (pow(mrho, 4) + 2 * (pow(mrho, 2) + s) * t - 4 * pow(t, 2)) +
               pow(t, 2) * (-pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                            2 * pow(s, 2) + 2 * s * t + pow(t, 2)) +
               2 * pow(mpion, 2) * t *
                   (pow(mrho, 4) + pow(mrho, 2) * (s + t) - 2 * t * (s + t))) +
          pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) *
                   (pow(mrho, 4) + 4 * pow(mrho, 2) * t - 2 * (s - 2 * t) * t) +
               pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t +
                            pow(t, 2) + 2 * pow(mrho, 2) * (-s + t)) -
               2 * pow(mpion, 2) * t *
                   (pow(mrho, 4) - pow(mrho, 2) * (s - 2 * t) +
                    2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) *
                   (3 * pow(mrho, 4) + 2 * pow(mrho, 2) * (s - 3 * t) -
                    2 * (s - 2 * t) * t) +
               t * (-pow(mrho, 2) + t) *
                   (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                    pow(mrho, 2) * (2 * s + t)) -
               2 * pow(mpion, 2) * (-pow(mrho, 2) + t) *
                   (2 * t * (s + t) - pow(mrho, 2) * (s + 2 * t))))) /
            ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             pow(pow(ma1, 2) - t, 2)) +
        (8 * (-2 + delta) *
         ((-2 + delta) * pow(mrho, 6) +
          pow(mpion, 6) * (-2 + 3 * delta - 8 * C4 * pow(mrho, 2)) +
          s * t * ((-2 + 3 * delta) * s + 4 * delta * t) +
          pow(mpion, 4) *
              (8 * C4 * pow(mrho, 4) + 4 * delta * s + 2 * t - 3 * delta * t -
               pow(mrho, 2) * (2 + delta + 16 * C4 * s - 8 * C4 * t)) +
          pow(mrho, 4) *
              (-((-2 + delta) * t) + s * (4 - 2 * delta + 8 * C4 * t)) +
          pow(mrho, 2) * s *
              (s * (-2 + delta - 8 * C4 * t) - 2 * t * (delta + 8 * C4 * t)) +
          pow(mpion, 2) *
              (s * ((2 - 3 * delta) * s - 8 * delta * t) -
               pow(mrho, 4) * (-6 + 3 * delta + 8 * C4 * (s + t)) +
               pow(mrho, 2) * (8 * C4 * pow(s, 2) + 4 * (-1 + delta) * t +
                               s * (-8 + 6 * delta + 32 * C4 * t))))) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(mpion, 2) - t)) +
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta1, 2) * (pow(mpion, 8) +
                          pow(mpion, 4) * (2 * pow(mrho, 4) + 2 * s * t -
                                           3 * pow(mrho, 2) * (s + t)) +
                          s * t *
                              (2 * pow(mrho, 4) + pow(s, 2) + 3 * s * t +
                               pow(t, 2) - 3 * pow(mrho, 2) * (s + t)) -
                          2 * pow(mpion, 2) * (pow(mrho, 2) - s - t) *
                              (-2 * s * t + pow(mrho, 2) * (s + t))) +
          pow(eta2, 2) * (pow(mpion, 8) -
                          4 * pow(mpion, 2) * s * t * (pow(mrho, 2) + s + t) +
                          pow(mpion, 4) * (2 * s * t + pow(mrho, 2) * (s + t)) +
                          s * t *
                              (pow(s, 2) + 3 * s * t + pow(t, 2) +
                               pow(mrho, 2) * (s + t))) +
          2 * eta1 * eta2 *
              (-pow(mpion, 8) + 2 * pow(mpion, 6) * pow(mrho, 2) -
               2 * pow(mpion, 4) * s * t -
               s * t *
                   (pow(s, 2) + 3 * s * t + pow(t, 2) -
                    2 * pow(mrho, 2) * (s + t)) -
               pow(mpion, 2) *
                   (-4 * s * t * (s + t) +
                    pow(mrho, 2) * (pow(s, 2) + 4 * s * t + pow(t, 2)))))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(ma1, 2) - t)) +
        (8 *
         (pow(delta, 2) * (8 * pow(mpion, 4) + 3 * pow(mrho, 4) -
                           6 * pow(mrho, 2) * (s + t) + 2 * pow(s + t, 2) +
                           4 * pow(mpion, 2) *
                               (3 * pow(mrho, 2) - 2 * (s + t))) -
          4 * delta * pow(mrho, 2) *
              (16 * C4 * pow(mpion, 4) + pow(mrho, 2) * (3 - 6 * C4 * (s + t)) +
               (s + t) * (-3 + 4 * C4 * (s + t)) +
               2 * pow(mpion, 2) *
                   (3 + C4 * (6 * pow(mrho, 2) - 8 * (s + t)))) +
          4 * pow(mrho, 4) *
              (3 + 4 * C4 * (2 * pow(mpion, 2) - s - t) *
                       (3 + C4 * (4 * pow(mpion, 2) - 2 * (s + t)))))) /
            (pow(mrho, 4) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (-pow(ma1, 2) + s) *
         (eta2 * (-2 * pow(mpion, 4) * (delta - 4 * C4 * pow(mrho, 2)) *
                      (pow(mrho, 2) + 4 * s) +
                  pow(mpion, 2) *
                      (-2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * s) +
                       8 * delta * s * (s + t) -
                       pow(mrho, 2) * ((-10 + delta) * s - (-2 + delta) * t +
                                       32 * C4 * s * (s + t))) +
                  s * (2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * s) -
                       2 * delta * pow(s + t, 2) +
                       pow(mrho, 2) * ((-6 + delta) * s + (-2 + delta) * t +
                                       8 * C4 * pow(s + t, 2)))) +
          eta1 *
              (4 * pow(mpion, 4) *
                   (6 * C4 * pow(mrho, 4) + 2 * delta * s +
                    pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * s)) +
               2 * delta * s * pow(s + t, 2) -
               pow(mrho, 2) *
                   ((-6 + 5 * delta) * pow(s, 2) +
                    2 * (-2 + 3 * delta) * s * t + (-2 + delta) * pow(t, 2) +
                    8 * C4 * s * pow(s + t, 2)) +
               pow(mrho, 4) *
                   ((-2 + delta) * (3 * s + t) + 8 * C4 * s * (s + 2 * t)) -
               2 * pow(mpion, 2) *
                   (4 * delta * s * (s + t) -
                    pow(mrho, 2) * (-6 * s + 7 * delta * s - 2 * t +
                                    3 * delta * t + 16 * C4 * s * (s + t)) +
                    2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * (2 * s + t)))))) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) *
         (((-2 + delta) *
           (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
            s * (pow(mrho, 2) - s - 2 * t)) *
           (eta1 * (pow(mrho, 2) - s - t) + eta2 * (pow(mpion, 2) + t))) /
              ((pow(mpion, 2) - s) * (pow(ma1, 2) - t)) +
          ((-2 + delta) *
           (eta2 * (pow(mpion, 2) + t) *
                (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * t) +
                 (pow(mrho, 2) - 2 * s - t) * t) +
            eta1 * (-4 * pow(mpion, 6) +
                    (pow(mrho, 2) - t) * (pow(mrho, 2) - s - t) * t +
                    pow(mpion, 4) * (3 * pow(mrho, 2) + s + t) -
                    pow(mpion, 2) * (pow(mrho, 4) + pow(mrho, 2) * (s - t) +
                                     2 * t * (-s + t))))) /
              ((-pow(ma1, 2) + t) * (-pow(mpion, 2) + t)) +
          (eta2 *
               (-2 * pow(mpion, 4) * (delta - 4 * C4 * pow(mrho, 2)) *
                    (pow(mrho, 2) + 4 * t) +
                pow(mpion, 2) * (8 * delta * t * (s + t) -
                                 2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * t) -
                                 pow(mrho, 2) *
                                     (-((-2 + delta) * s) + (-10 + delta) * t +
                                      32 * C4 * t * (s + t))) +
                t * (-2 * delta * pow(s + t, 2) +
                     2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * t) +
                     pow(mrho, 2) * ((-2 + delta) * s + (-6 + delta) * t +
                                     8 * C4 * pow(s + t, 2)))) +
           eta1 *
               (2 * delta * t * pow(s + t, 2) -
                pow(mrho, 2) *
                    ((-2 + delta) * pow(s, 2) + 2 * (-2 + 3 * delta) * s * t +
                     (-6 + 5 * delta) * pow(t, 2) +
                     8 * C4 * t * pow(s + t, 2)) +
                pow(mrho, 4) *
                    (8 * C4 * t * (2 * s + t) + (-2 + delta) * (s + 3 * t)) +
                4 * pow(mpion, 4) *
                    (6 * C4 * pow(mrho, 4) + 2 * delta * t +
                     pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * t)) -
                2 * pow(mpion, 2) *
                    (4 * delta * t * (s + t) -
                     pow(mrho, 2) * (-2 * s + 3 * delta * s - 6 * t +
                                     7 * delta * t + 16 * C4 * t * (s + t)) +
                     2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * (s + 2 * t))))) /
              (pow(mrho, 2) * (-pow(ma1, 2) + t)))) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)))) /
      (512. * Pi);

  return to_mb * diff_xs / spin_deg_factor;
}

// C12
double
PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho_pi_rho_mediated(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_rho, m_pion_, 0.);
  const double &tmin = t_mandelstam[1];
  const double &tmax = t_mandelstam[0];
  const double spin_deg_factor = 3.0;

  const double xs =
      (pow(Const, 2) * pow(ghat, 4) *
       (0. -
        (0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
          2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
         tmax) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - s, 2)) -
        (0.0625 *
         (0. + 8. * pow(mpion, 2) * pow(mrho, 4) - 12. * pow(mrho, 6) +
          4. * pow(mrho, 4) * s +
          delta * pow(mrho, 2) *
              (-16. * pow(mpion, 4) - 16. * pow(mpion, 2) * pow(mrho, 2) -
               4. * pow(mrho, 4) + 16. * pow(mrho, 2) * s + 4. * pow(s, 2)) +
          pow(delta, 2) *
              (8. * pow(mpion, 6) + 9. * pow(mrho, 6) +
               pow(mpion, 4) * (4. * pow(mrho, 2) - 4. * s) -
               13. * pow(mrho, 4) * s - 5. * pow(mrho, 2) * pow(s, 2) +
               1. * pow(s, 3) +
               pow(mpion, 2) * (-2. * pow(mrho, 4) + 4. * pow(mrho, 2) * s -
                                2. * pow(s, 2)))) *
         tmax) /
            pow(mrho, 6) -
        (0.125 * (-2. + delta) * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (eta2 * (pow(mpion, 6) + pow(mpion, 2) * pow(s, 2) +
                  (pow(mrho, 2) - 1. * s) * pow(s, 2) +
                  pow(mpion, 4) * (-1. * pow(mrho, 2) + 3. * s)) +
          eta1 *
              (-4. * pow(mpion, 6) + pow(mpion, 4) * (3. * pow(mrho, 2) + s) +
               pow(mpion, 2) *
                   (-1. * pow(mrho, 4) + pow(mrho, 2) * s - 2. * pow(s, 2)) +
               s * (pow(mrho, 4) - 2. * pow(mrho, 2) * s + pow(s, 2)))) *
         tmax) /
            ((-1. * pow(mpion, 2) + s) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
              2. * pow(ma1, 2) * s + pow(s, 2))) +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta1, 2) *
              (1. * pow(mpion, 8) - 2. * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 2) * s *
                   (-4. * pow(mrho, 4) + 8. * pow(mrho, 2) * s -
                    4. * pow(s, 2)) +
               pow(s, 2) * (1. * pow(mrho, 4) - 2. * pow(mrho, 2) * s +
                            1. * pow(s, 2)) +
               pow(mpion, 4) * (3. * pow(mrho, 4) - 6. * pow(mrho, 2) * s +
                                4. * pow(s, 2))) +
          pow(eta2, 2) *
              (1. * pow(mpion, 8) - 2. * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 2) * s *
                   (-2. * pow(mrho, 4) - 4. * pow(mrho, 2) * s -
                    4. * pow(s, 2)) +
               pow(s, 2) * (1. * pow(mrho, 4) + 2. * pow(mrho, 2) * s +
                            1. * pow(s, 2)) +
               pow(mpion, 4) * (1. * pow(mrho, 4) + 4. * pow(mrho, 2) * s +
                                4. * pow(s, 2))) +
          eta1 * eta2 *
              (-2. * pow(mpion, 8) + 2. * pow(mrho, 4) * pow(s, 2) -
               2. * pow(s, 4) +
               pow(mpion, 4) * (2. * pow(mrho, 4) + 4. * pow(mrho, 2) * s -
                                8. * pow(s, 2)) +
               pow(mpion, 2) * s *
                   (-4. * pow(mrho, 4) - 4. * pow(mrho, 2) * s +
                    8. * pow(s, 2)))) *
         tmax) /
            (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
             2. * pow(ma1, 2) * s + pow(s, 2)) +
        (0.5 *
         (0. - 4. * C4 * pow(mrho, 8) - 0.5 * pow(mrho, 4) * s +
          pow(mrho, 6) * (2. + 2. * C4 * s) +
          pow(delta, 2) *
              (2. * pow(mpion, 6) + 2. * pow(mrho, 6) -
               1.5 * pow(mpion, 4) * s - 2.375 * pow(mrho, 4) * s -
               0.75 * pow(mrho, 2) * pow(s, 2) + 0.125 * pow(s, 3) +
               pow(mpion, 2) * (-1.5 * pow(mrho, 4) + 0.5 * pow(mrho, 2) * s)) +
          delta * pow(mrho, 2) *
              (-2. * pow(mpion, 4) +
               pow(mpion, 2) * (1. * s + pow(mrho, 2) * (-1. - 2. * C4 * s)) +
               pow(mrho, 2) * (2. * C4 * pow(mrho, 4) +
                               pow(mrho, 2) * (-3. + 1. * C4 * s) +
                               s * (2. + 1. * C4 * s)))) *
         tmax) /
            pow(mrho, 6) -
        (0.5 *
         (pow(mrho, 6) *
              (-1.5 + C4 * (-12. * pow(mpion, 2) + 6. * s) +
               pow(C4, 2) * (-16. * pow(mpion, 4) + 16. * pow(mpion, 2) * s -
                             4. * pow(s, 2))) +
          pow(delta, 2) *
              (1. * pow(mpion, 6) + 0.125 * pow(mrho, 6) +
               pow(mpion, 4) * (-2. * pow(mrho, 2) - 1. * s) +
               0.25 * pow(mrho, 4) * s - 0.625 * pow(mrho, 2) * pow(s, 2) +
               pow(mpion, 2) * (-2.5 * pow(mrho, 4) + 1.75 * pow(mrho, 2) * s +
                                0.25 * pow(s, 2))) +
          delta * pow(mrho, 2) *
              (pow(mpion, 4) * (1. + 8. * C4 * pow(mrho, 2)) +
               pow(mpion, 2) * (6. * C4 * pow(mrho, 4) - 0.5 * s +
                                pow(mrho, 2) * (3. - 10. * C4 * s)) +
               pow(mrho, 2) * (pow(mrho, 2) * (1.5 - 1. * C4 * s) +
                               s * (-2.5 + 3. * C4 * s)))) *
         tmax) /
            pow(mrho, 6) -
        (0.25 *
         (pow(delta, 2) *
              (1. * pow(mpion, 6) - 1. * pow(mrho, 6) +
               pow(mpion, 4) * (-2.499999999999999 * pow(mrho, 2) - 2.5 * s) -
               1.5 * pow(mrho, 4) * s + 2. * pow(mrho, 2) * pow(s, 2) -
               0.5 * pow(s, 3) +
               pow(mpion, 2) *
                   (3.5 * pow(mrho, 4) - 1.5000000000000004 * pow(mrho, 2) * s +
                    2. * pow(s, 2))) +
          pow(mrho, 2) *
              (pow(mpion, 4) * (-6. - 8. * C4 * pow(mrho, 2)) + 2. * pow(s, 2) +
               pow(mrho, 4) * (-4. - 8. * C4 * s) +
               pow(mrho, 2) * s * (-2. + 8. * C4 * s) +
               pow(mpion, 2) * (8. * C4 * pow(mrho, 4) + 4. * s +
                                pow(mrho, 2) * (10. - 16. * C4 * s))) +
          delta * (-2. * pow(mpion, 6) - 5. * pow(mrho, 2) * pow(s, 2) +
                   1. * pow(s, 3) +
                   pow(mpion, 4) *
                       (8. * pow(mrho, 2) + 4. * C4 * pow(mrho, 4) + 5. * s) +
                   pow(mrho, 4) * s * (4. - 4. * C4 * s) +
                   pow(mrho, 6) * (4. + 4. * C4 * s) +
                   pow(mpion, 2) * (-4. * C4 * pow(mrho, 6) +
                                    1. * pow(mrho, 2) * s - 4. * pow(s, 2) +
                                    pow(mrho, 4) * (-12. + 8. * C4 * s)))) *
         tmax) /
            (pow(mrho, 4) * (pow(mpion, 2) - 1. * s)) +
        (0.0625 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (pow(mrho, 2) *
              (eta2 * (4. * pow(mpion, 4) - 6. * pow(mpion, 2) * s +
                       s * (8. * pow(mrho, 2) + 6. * s)) +
               eta1 * (-12. * pow(mpion, 4) + 4. * pow(mrho, 4) +
                       2. * pow(mrho, 2) * s - 6. * pow(s, 2) +
                       pow(mpion, 2) * (-4. * pow(mrho, 2) + 8. * s))) +
          delta *
              (eta1 * (8. * pow(mpion, 6) - 2. * pow(mrho, 6) +
                       pow(mpion, 4) * (2. * pow(mrho, 2) - 2. * s) -
                       3. * pow(mrho, 4) * s + 4. * pow(mrho, 2) * pow(s, 2) +
                       1. * pow(s, 3) +
                       pow(mpion, 2) * (2. * pow(mrho, 4) - 2. * pow(s, 2))) +
               eta2 * (pow(mpion, 4) * (-2. * pow(mrho, 2) - 4. * s) +
                       pow(mpion, 2) * s * (3. * pow(mrho, 2) + 3. * s) +
                       s * (-4. * pow(mrho, 4) - 7. * pow(mrho, 2) * s -
                            1. * pow(s, 2))))) *
         tmax) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) -
        (0.1875 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (delta *
              (eta1 * (2.6666666666666665 * pow(mpion, 6) +
                       pow(mpion, 4) * (-4. * pow(mrho, 2) + 2. * s) +
                       pow(mpion, 2) * (-1.3333333333333333 * pow(mrho, 4) +
                                        6. * pow(mrho, 2) * s -
                                        3.3333333333333335 * pow(s, 2)) +
                       s * (0.3333333333333333 * pow(mrho, 4) -
                            1.3333333333333333 * pow(mrho, 2) * s +
                            1. * pow(s, 2))) +
               eta2 * (pow(mpion, 4) *
                           (-0.6666666666666666 * pow(mrho, 2) - 4. * s) +
                       s * (0.6666666666666666 * pow(mrho, 4) -
                            1. * pow(mrho, 2) * s - 1. * pow(s, 2)) +
                       pow(mpion, 2) * (-0.6666666666666666 * pow(mrho, 4) -
                                        0.3333333333333333 * pow(mrho, 2) * s +
                                        3.6666666666666665 * pow(s, 2)))) +
          pow(mrho, 2) *
              (eta2 * (C4 * pow(mpion, 4) *
                           (2.6666666666666665 * pow(mrho, 2) +
                            10.666666666666666 * s) +
                       pow(mpion, 2) *
                           (s * (3.3333333333333335 -
                                 10.666666666666666 * C4 * s) +
                            pow(mrho, 2) * (1.3333333333333333 -
                                            5.333333333333333 * C4 * s)) +
                       s * (s * (-2. + 2.6666666666666665 * C4 * s) +
                            pow(mrho, 2) * (-1.3333333333333333 +
                                            2.6666666666666665 * C4 * s))) +
               eta1 *
                   (pow(mpion, 4) *
                        (1.3333333333333333 + 8. * C4 * pow(mrho, 2) -
                         10.666666666666666 * C4 * s) +
                    s * (s * (2. - 2.6666666666666665 * C4 * s) +
                         pow(mrho, 2) * (-2. + 2.6666666666666665 * C4 * s)) +
                    pow(mpion, 2) *
                        (pow(mrho, 2) * (2.6666666666666665 -
                                         10.666666666666666 * C4 * s) +
                         s * (-4. + 10.666666666666666 * C4 * s))))) *
         tmax) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) -
        (0.0625 * (-2. + delta) * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (pow(mpion, 2) + s) *
         (-2. * eta2 * s + eta1 * (pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
         pow(tmax, 2)) /
            ((-1. * pow(mpion, 2) + s) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
              2. * pow(ma1, 2) * s + pow(s, 2))) +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta1, 2) * (pow(mrho, 2) - 1. * s) + 2. * eta1 * eta2 * s -
          1. * pow(eta2, 2) * s) *
         (pow(mpion, 4) + (pow(mrho, 2) - 1. * s) * s +
          pow(mpion, 2) * (-1. * pow(mrho, 2) + 2. * s)) *
         pow(tmax, 2)) /
            (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
             2. * pow(ma1, 2) * s + pow(s, 2)) -
        (0.125 *
         (-1. * pow(mrho, 4) + 4. * C4 * pow(mrho, 6) +
          delta * pow(mrho, 2) *
              (2. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
               pow(mpion, 2) * (2. - 4. * C4 * pow(mrho, 2)) - 2. * s) +
          pow(delta, 2) *
              (1. * pow(mpion, 4) + 0.25 * pow(mrho, 4) - 1.25 * pow(s, 2) +
               pow(mpion, 2) * (-3. * pow(mrho, 2) + 2. * s))) *
         pow(tmax, 2)) /
            pow(mrho, 6) +
        (0.03125 *
         (0. - 4. * pow(mrho, 4) +
          delta * (16. * pow(mrho, 4) - 8. * pow(mrho, 2) * s) +
          pow(delta, 2) * (4. * pow(mpion, 4) - 3. * pow(mrho, 4) +
                           2. * pow(mrho, 2) * s - 3. * pow(s, 2) +
                           pow(mpion, 2) * (-4. * pow(mrho, 2) + 4. * s))) *
         pow(tmax, 2)) /
            pow(mrho, 6) +
        (0.0625 *
         (-32. * C4 * pow(mrho, 4) * s +
          pow(delta, 2) *
              (1. * pow(mpion, 4) +
               pow(mpion, 2) * (-1.0000000000000009 * pow(mrho, 2) - 2. * s) +
               s * (-3. * pow(mrho, 2) + 1. * s)) +
          delta * (-2. * pow(mpion, 4) +
                   (6. * pow(mrho, 2) + 16. * C4 * pow(mrho, 4) - 2. * s) * s +
                   pow(mpion, 2) * (2. * pow(mrho, 2) + 4. * s))) *
         pow(tmax, 2)) /
            (pow(mrho, 4) * (pow(mpion, 2) - 1. * s)) -
        (0.5625 *
         (C4 * pow(mrho, 6) *
              (2.6666666666666665 + 7.111111111111112 * C4 * pow(mpion, 2) -
               3.555555555555556 * C4 * s) +
          pow(delta, 2) *
              (0.11111111111111112 * pow(mrho, 4) +
               pow(mpion, 2) * (1. * pow(mrho, 2) - 0.22222222222222224 * s) -
               0.22222222222222224 * pow(mrho, 2) * s +
               0.11111111111111112 * pow(s, 2)) +
          delta * pow(mrho, 2) *
              (-2.2222222222222223 * C4 * pow(mrho, 4) +
               pow(mpion, 2) * (-0.6666666666666666 -
                                2.6666666666666665 * C4 * pow(mrho, 2)) +
               0.22222222222222224 * s +
               pow(mrho, 2) *
                   (-0.22222222222222224 + 1.777777777777778 * C4 * s))) *
         pow(tmax, 2)) /
            pow(mrho, 6) +
        (0.03125 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (pow(mrho, 2) * (-2. * eta2 * pow(mpion, 2) -
                          5.999999999999999 * eta1 * pow(mrho, 2) +
                          8. * eta1 * s - 2. * eta2 * s) +
          delta *
              (eta1 * (-5.999999999999999 * pow(mpion, 4) + 5. * pow(mrho, 4) +
                       pow(mpion, 2) * (4. * pow(mrho, 2) - 4. * s) -
                       5.999999999999999 * pow(mrho, 2) * s + 1. * pow(s, 2)) +
               eta2 * (4. * pow(mpion, 4) +
                       pow(mpion, 2) * (1. * pow(mrho, 2) - 2. * s) +
                       s * (5. * pow(mrho, 2) + 2. * s)))) *
         pow(tmax, 2)) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) -
        (0.15625 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (delta *
              (eta1 * (-1.2 * pow(mpion, 4) + 0.6 * pow(mrho, 4) +
                       pow(mpion, 2) * (2. * pow(mrho, 2) - 2.4 * s) -
                       1.6 * pow(mrho, 2) * s + 1. * pow(s, 2)) +
               eta2 * (0.8 * pow(mpion, 4) + (1. * pow(mrho, 2) - 0.4 * s) * s +
                       pow(mpion, 2) * (0.2 * pow(mrho, 2) + 1.2 * s))) +
          pow(mrho, 2) *
              (eta2 * (pow(mpion, 2) * (-0.4 - 6.4 * C4 * s) +
                       s * (-0.4 + 3.2 * C4 * s)) +
               eta1 * (s * (0.8 - 3.2 * C4 * s) +
                       pow(mrho, 2) * (-0.4 + 3.2 * C4 * s) +
                       pow(mpion, 2) *
                           (-0.8 - 3.2 * C4 * pow(mrho, 2) + 6.4 * C4 * s)))) *
         pow(tmax, 2)) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) -
        (0.20833333333333331 * delta *
         (-0.8 * pow(mrho, 2) + 0.8 * C4 * pow(mrho, 4) +
          delta * (0.8 * pow(mpion, 2) + 1. * pow(mrho, 2) - 0.7 * s)) *
         pow(tmax, 3)) /
            pow(mrho, 6) +
        (0.125 *
         (5.333333333333333 * pow(C4, 2) * pow(mrho, 6) +
          delta * (-0.6666666666666666 * pow(mrho, 2) -
                   1.3333333333333333 * C4 * pow(mrho, 4)) +
          pow(delta, 2) *
              (1. * pow(mpion, 2) + 1.1666666666666667 * pow(mrho, 2) -
               0.6666666666666666 * s)) *
         pow(tmax, 3)) /
            pow(mrho, 6) +
        (0.10416666666666666 * delta *
         (-0.8 * pow(mrho, 2) +
          delta * (0.4 * pow(mpion, 2) + 1. * pow(mrho, 2) - 0.6 * s)) *
         pow(tmax, 3)) /
            pow(mrho, 6) +
        (0.020833333333333332 * pow(eta1 - 1. * eta2, 2) * s *
         (-2. * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-1. * pow(mrho, 2) + s)) *
         pow(tmax, 3)) /
            (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
             2. * pow(ma1, 2) * s + pow(s, 2)) +
        (0.10416666666666666 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (0.4 * eta1 * pow(mrho, 2) +
          delta *
              (-0.2 * eta2 * pow(mpion, 2) - 0.2 * eta2 * s +
               eta1 * (-0.4 * pow(mpion, 2) - 0.8 * pow(mrho, 2) + 1. * s))) *
         pow(tmax, 3)) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) -
        (0.14583333333333331 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (delta * (-0.14285714285714285 * eta2 * pow(mpion, 2) -
                   0.42857142857142855 * eta2 * s +
                   eta1 * (-0.2857142857142857 * pow(mpion, 2) -
                           0.5714285714285714 * pow(mrho, 2) + 1. * s)) +
          pow(mrho, 2) *
              (1.1428571428571428 * C4 * eta2 * s +
               eta1 * (0.2857142857142857 - 1.1428571428571428 * C4 * s))) *
         pow(tmax, 3)) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) +
        (0. -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(mpion, 4) *
             pow(mrho, 6) +
         0.25 * pow(2. - 1. * delta, 2) * pow(mrho, 10) -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(mrho, 6) *
             pow(s, 2) +
         pow(2. - 1. * delta, 2) * pow(mpion, 2) * pow(mrho, 6) *
             (-1. * pow(mrho, 2) + 1. * s)) /
            (pow(mrho, 6) * (-2. * pow(mpion, 2) + 1. * s + 1. * tmax)) +
        (2. *
         (0. - 2. * pow(mpion, 4) * pow(mrho, 4) - 0.5 * pow(mrho, 8) +
          delta * pow(mrho, 4) *
              (2. * pow(mpion, 4) + 0.5 * pow(mrho, 4) +
               pow(mpion, 2) * (-2. * pow(mrho, 2) - 1.9999999999999998 * s)) +
          pow(mpion, 2) * (2. * pow(mrho, 6) + 2. * pow(mrho, 4) * s) +
          pow(delta, 2) * pow(mrho, 2) *
              (-2.220446049250313e-16 * pow(mpion, 6) - 0.125 * pow(mrho, 6) +
               pow(mpion, 4) *
                   (-0.5 * pow(mrho, 2) + 2.220446049250313e-16 * s) +
               pow(mpion, 2) * (0.5 * pow(mrho, 4) + 0.5 * pow(mrho, 2) * s))) *
         log(abs(-1. * pow(mpion, 2) + 0.5 * s + 0.5 * tmax))) /
            (pow(mrho, 4) * (pow(mpion, 2) - 1. * s)) -
        (0.25 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (eta2 * ((-2. + 1. * delta) * pow(mpion, 6) +
                  (6. - 3. * delta) * pow(mpion, 4) * s +
                  pow(s, 2) * ((4. - 2. * delta) * pow(mrho, 2) +
                               (2. - 1. * delta) * s) +
                  pow(mpion, 2) * s *
                      ((-4. + 2. * delta) * pow(mrho, 2) +
                       (-6. + 3. * delta) * s)) +
          eta1 * ((2. - 1. * delta) * pow(mpion, 6) +
                  (2. - 1. * delta) * pow(mrho, 4) * s +
                  (-2. + 1. * delta) * pow(s, 3) +
                  pow(mpion, 4) * ((4. - 2. * delta) * pow(mrho, 2) +
                                   (-6. + 3. * delta) * s) +
                  pow(mpion, 2) * ((-2. + 1. * delta) * pow(mrho, 4) +
                                   (-4. + 2. * delta) * pow(mrho, 2) * s +
                                   (6. - 3. * delta) * pow(s, 2)))) *
         log(abs(-2. * pow(mpion, 2) + s + tmax))) /
            (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
             2. * pow(ma1, 2) * s + pow(s, 2)) +
        (0.25 *
         (0. + 8. * pow(mpion, 2) * pow(mrho, 6) + 4. * pow(mrho, 8) -
          8. * pow(mrho, 6) * s +
          delta * pow(mrho, 4) *
              (8. * pow(mpion, 4) - 8. * pow(mrho, 4) +
               pow(mpion, 2) * (8. * pow(mrho, 2) - 16. * s) +
               8. * pow(mrho, 2) * s + 8. * pow(s, 2)) +
          pow(delta, 2) * pow(mrho, 4) *
              (-4. * pow(mpion, 4) + 3. * pow(mrho, 4) - 2. * pow(mrho, 2) * s -
               4. * pow(s, 2) +
               pow(mpion, 2) * (-6. * pow(mrho, 2) + 8. * s))) *
         log(abs(-2. * pow(mpion, 2) + 1. * s + 1. * tmax))) /
            pow(mrho, 6) +
        (0.5 *
         (0. + pow(mpion, 2) * (4. * pow(mrho, 6) - 8. * C4 * pow(mrho, 8)) -
          4. * pow(mrho, 6) * s + 8. * C4 * pow(mrho, 8) * s +
          pow(delta, 2) * pow(mrho, 4) *
              (2. * pow(mpion, 4) - 1. * pow(mrho, 4) +
               pow(mpion, 2) * (4. * pow(mrho, 2) - 4. * s) + 2. * pow(s, 2)) +
          delta * pow(mrho, 4) *
              (-4. * pow(mpion, 4) + 2. * pow(mrho, 2) * s - 4. * pow(s, 2) +
               pow(mpion, 2) *
                   (-10. * pow(mrho, 2) + 4. * C4 * pow(mrho, 4) + 8. * s) +
               pow(mrho, 4) * (2. - 4. * C4 * s))) *
         log(abs(-2. * pow(mpion, 2) + 1. * s + 1. * tmax))) /
            pow(mrho, 6))) /
          (16. * Pi *
           (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
            2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
      (pow(Const, 2) * pow(ghat, 4) *
       (0. -
        (0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
          2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
         tmin) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - s, 2)) -
        (0.0625 *
         (0. + 8. * pow(mpion, 2) * pow(mrho, 4) - 12. * pow(mrho, 6) +
          4. * pow(mrho, 4) * s +
          delta * pow(mrho, 2) *
              (-16. * pow(mpion, 4) - 16. * pow(mpion, 2) * pow(mrho, 2) -
               4. * pow(mrho, 4) + 16. * pow(mrho, 2) * s + 4. * pow(s, 2)) +
          pow(delta, 2) *
              (8. * pow(mpion, 6) + 9. * pow(mrho, 6) +
               pow(mpion, 4) * (4. * pow(mrho, 2) - 4. * s) -
               13. * pow(mrho, 4) * s - 5. * pow(mrho, 2) * pow(s, 2) +
               1. * pow(s, 3) +
               pow(mpion, 2) * (-2. * pow(mrho, 4) + 4. * pow(mrho, 2) * s -
                                2. * pow(s, 2)))) *
         tmin) /
            pow(mrho, 6) -
        (0.125 * (-2. + delta) * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (eta2 * (pow(mpion, 6) + pow(mpion, 2) * pow(s, 2) +
                  (pow(mrho, 2) - 1. * s) * pow(s, 2) +
                  pow(mpion, 4) * (-1. * pow(mrho, 2) + 3. * s)) +
          eta1 *
              (-4. * pow(mpion, 6) + pow(mpion, 4) * (3. * pow(mrho, 2) + s) +
               pow(mpion, 2) *
                   (-1. * pow(mrho, 4) + pow(mrho, 2) * s - 2. * pow(s, 2)) +
               s * (pow(mrho, 4) - 2. * pow(mrho, 2) * s + pow(s, 2)))) *
         tmin) /
            ((-1. * pow(mpion, 2) + s) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
              2. * pow(ma1, 2) * s + pow(s, 2))) +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta1, 2) *
              (1. * pow(mpion, 8) - 2. * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 2) * s *
                   (-4. * pow(mrho, 4) + 8. * pow(mrho, 2) * s -
                    4. * pow(s, 2)) +
               pow(s, 2) * (1. * pow(mrho, 4) - 2. * pow(mrho, 2) * s +
                            1. * pow(s, 2)) +
               pow(mpion, 4) * (3. * pow(mrho, 4) - 6. * pow(mrho, 2) * s +
                                4. * pow(s, 2))) +
          pow(eta2, 2) *
              (1. * pow(mpion, 8) - 2. * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 2) * s *
                   (-2. * pow(mrho, 4) - 4. * pow(mrho, 2) * s -
                    4. * pow(s, 2)) +
               pow(s, 2) * (1. * pow(mrho, 4) + 2. * pow(mrho, 2) * s +
                            1. * pow(s, 2)) +
               pow(mpion, 4) * (1. * pow(mrho, 4) + 4. * pow(mrho, 2) * s +
                                4. * pow(s, 2))) +
          eta1 * eta2 *
              (-2. * pow(mpion, 8) + 2. * pow(mrho, 4) * pow(s, 2) -
               2. * pow(s, 4) +
               pow(mpion, 4) * (2. * pow(mrho, 4) + 4. * pow(mrho, 2) * s -
                                8. * pow(s, 2)) +
               pow(mpion, 2) * s *
                   (-4. * pow(mrho, 4) - 4. * pow(mrho, 2) * s +
                    8. * pow(s, 2)))) *
         tmin) /
            (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
             2. * pow(ma1, 2) * s + pow(s, 2)) +
        (0.5 *
         (0. - 4. * C4 * pow(mrho, 8) - 0.5 * pow(mrho, 4) * s +
          pow(mrho, 6) * (2. + 2. * C4 * s) +
          pow(delta, 2) *
              (2. * pow(mpion, 6) + 2. * pow(mrho, 6) -
               1.5 * pow(mpion, 4) * s - 2.375 * pow(mrho, 4) * s -
               0.75 * pow(mrho, 2) * pow(s, 2) + 0.125 * pow(s, 3) +
               pow(mpion, 2) * (-1.5 * pow(mrho, 4) + 0.5 * pow(mrho, 2) * s)) +
          delta * pow(mrho, 2) *
              (-2. * pow(mpion, 4) +
               pow(mpion, 2) * (1. * s + pow(mrho, 2) * (-1. - 2. * C4 * s)) +
               pow(mrho, 2) * (2. * C4 * pow(mrho, 4) +
                               pow(mrho, 2) * (-3. + 1. * C4 * s) +
                               s * (2. + 1. * C4 * s)))) *
         tmin) /
            pow(mrho, 6) -
        (0.5 *
         (pow(mrho, 6) *
              (-1.5 + C4 * (-12. * pow(mpion, 2) + 6. * s) +
               pow(C4, 2) * (-16. * pow(mpion, 4) + 16. * pow(mpion, 2) * s -
                             4. * pow(s, 2))) +
          pow(delta, 2) *
              (1. * pow(mpion, 6) + 0.125 * pow(mrho, 6) +
               pow(mpion, 4) * (-2. * pow(mrho, 2) - 1. * s) +
               0.25 * pow(mrho, 4) * s - 0.625 * pow(mrho, 2) * pow(s, 2) +
               pow(mpion, 2) * (-2.5 * pow(mrho, 4) + 1.75 * pow(mrho, 2) * s +
                                0.25 * pow(s, 2))) +
          delta * pow(mrho, 2) *
              (pow(mpion, 4) * (1. + 8. * C4 * pow(mrho, 2)) +
               pow(mpion, 2) * (6. * C4 * pow(mrho, 4) - 0.5 * s +
                                pow(mrho, 2) * (3. - 10. * C4 * s)) +
               pow(mrho, 2) * (pow(mrho, 2) * (1.5 - 1. * C4 * s) +
                               s * (-2.5 + 3. * C4 * s)))) *
         tmin) /
            pow(mrho, 6) -
        (0.25 *
         (pow(delta, 2) *
              (1. * pow(mpion, 6) - 1. * pow(mrho, 6) +
               pow(mpion, 4) * (-2.499999999999999 * pow(mrho, 2) - 2.5 * s) -
               1.5 * pow(mrho, 4) * s + 2. * pow(mrho, 2) * pow(s, 2) -
               0.5 * pow(s, 3) +
               pow(mpion, 2) *
                   (3.5 * pow(mrho, 4) - 1.5000000000000004 * pow(mrho, 2) * s +
                    2. * pow(s, 2))) +
          pow(mrho, 2) *
              (pow(mpion, 4) * (-6. - 8. * C4 * pow(mrho, 2)) + 2. * pow(s, 2) +
               pow(mrho, 4) * (-4. - 8. * C4 * s) +
               pow(mrho, 2) * s * (-2. + 8. * C4 * s) +
               pow(mpion, 2) * (8. * C4 * pow(mrho, 4) + 4. * s +
                                pow(mrho, 2) * (10. - 16. * C4 * s))) +
          delta * (-2. * pow(mpion, 6) - 5. * pow(mrho, 2) * pow(s, 2) +
                   1. * pow(s, 3) +
                   pow(mpion, 4) *
                       (8. * pow(mrho, 2) + 4. * C4 * pow(mrho, 4) + 5. * s) +
                   pow(mrho, 4) * s * (4. - 4. * C4 * s) +
                   pow(mrho, 6) * (4. + 4. * C4 * s) +
                   pow(mpion, 2) * (-4. * C4 * pow(mrho, 6) +
                                    1. * pow(mrho, 2) * s - 4. * pow(s, 2) +
                                    pow(mrho, 4) * (-12. + 8. * C4 * s)))) *
         tmin) /
            (pow(mrho, 4) * (pow(mpion, 2) - 1. * s)) +
        (0.0625 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (pow(mrho, 2) *
              (eta2 * (4. * pow(mpion, 4) - 6. * pow(mpion, 2) * s +
                       s * (8. * pow(mrho, 2) + 6. * s)) +
               eta1 * (-12. * pow(mpion, 4) + 4. * pow(mrho, 4) +
                       2. * pow(mrho, 2) * s - 6. * pow(s, 2) +
                       pow(mpion, 2) * (-4. * pow(mrho, 2) + 8. * s))) +
          delta *
              (eta1 * (8. * pow(mpion, 6) - 2. * pow(mrho, 6) +
                       pow(mpion, 4) * (2. * pow(mrho, 2) - 2. * s) -
                       3. * pow(mrho, 4) * s + 4. * pow(mrho, 2) * pow(s, 2) +
                       1. * pow(s, 3) +
                       pow(mpion, 2) * (2. * pow(mrho, 4) - 2. * pow(s, 2))) +
               eta2 * (pow(mpion, 4) * (-2. * pow(mrho, 2) - 4. * s) +
                       pow(mpion, 2) * s * (3. * pow(mrho, 2) + 3. * s) +
                       s * (-4. * pow(mrho, 4) - 7. * pow(mrho, 2) * s -
                            1. * pow(s, 2))))) *
         tmin) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) -
        (0.1875 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (delta *
              (eta1 * (2.6666666666666665 * pow(mpion, 6) +
                       pow(mpion, 4) * (-4. * pow(mrho, 2) + 2. * s) +
                       pow(mpion, 2) * (-1.3333333333333333 * pow(mrho, 4) +
                                        6. * pow(mrho, 2) * s -
                                        3.3333333333333335 * pow(s, 2)) +
                       s * (0.3333333333333333 * pow(mrho, 4) -
                            1.3333333333333333 * pow(mrho, 2) * s +
                            1. * pow(s, 2))) +
               eta2 * (pow(mpion, 4) *
                           (-0.6666666666666666 * pow(mrho, 2) - 4. * s) +
                       s * (0.6666666666666666 * pow(mrho, 4) -
                            1. * pow(mrho, 2) * s - 1. * pow(s, 2)) +
                       pow(mpion, 2) * (-0.6666666666666666 * pow(mrho, 4) -
                                        0.3333333333333333 * pow(mrho, 2) * s +
                                        3.6666666666666665 * pow(s, 2)))) +
          pow(mrho, 2) *
              (eta2 * (C4 * pow(mpion, 4) *
                           (2.6666666666666665 * pow(mrho, 2) +
                            10.666666666666666 * s) +
                       pow(mpion, 2) *
                           (s * (3.3333333333333335 -
                                 10.666666666666666 * C4 * s) +
                            pow(mrho, 2) * (1.3333333333333333 -
                                            5.333333333333333 * C4 * s)) +
                       s * (s * (-2. + 2.6666666666666665 * C4 * s) +
                            pow(mrho, 2) * (-1.3333333333333333 +
                                            2.6666666666666665 * C4 * s))) +
               eta1 *
                   (pow(mpion, 4) *
                        (1.3333333333333333 + 8. * C4 * pow(mrho, 2) -
                         10.666666666666666 * C4 * s) +
                    s * (s * (2. - 2.6666666666666665 * C4 * s) +
                         pow(mrho, 2) * (-2. + 2.6666666666666665 * C4 * s)) +
                    pow(mpion, 2) *
                        (pow(mrho, 2) * (2.6666666666666665 -
                                         10.666666666666666 * C4 * s) +
                         s * (-4. + 10.666666666666666 * C4 * s))))) *
         tmin) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) -
        (0.0625 * (-2. + delta) * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (pow(mpion, 2) + s) *
         (-2. * eta2 * s + eta1 * (pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
         pow(tmin, 2)) /
            ((-1. * pow(mpion, 2) + s) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
              2. * pow(ma1, 2) * s + pow(s, 2))) +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta1, 2) * (pow(mrho, 2) - 1. * s) + 2. * eta1 * eta2 * s -
          1. * pow(eta2, 2) * s) *
         (pow(mpion, 4) + (pow(mrho, 2) - 1. * s) * s +
          pow(mpion, 2) * (-1. * pow(mrho, 2) + 2. * s)) *
         pow(tmin, 2)) /
            (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
             2. * pow(ma1, 2) * s + pow(s, 2)) -
        (0.125 *
         (-1. * pow(mrho, 4) + 4. * C4 * pow(mrho, 6) +
          delta * pow(mrho, 2) *
              (2. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
               pow(mpion, 2) * (2. - 4. * C4 * pow(mrho, 2)) - 2. * s) +
          pow(delta, 2) *
              (1. * pow(mpion, 4) + 0.25 * pow(mrho, 4) - 1.25 * pow(s, 2) +
               pow(mpion, 2) * (-3. * pow(mrho, 2) + 2. * s))) *
         pow(tmin, 2)) /
            pow(mrho, 6) +
        (0.03125 *
         (0. - 4. * pow(mrho, 4) +
          delta * (16. * pow(mrho, 4) - 8. * pow(mrho, 2) * s) +
          pow(delta, 2) * (4. * pow(mpion, 4) - 3. * pow(mrho, 4) +
                           2. * pow(mrho, 2) * s - 3. * pow(s, 2) +
                           pow(mpion, 2) * (-4. * pow(mrho, 2) + 4. * s))) *
         pow(tmin, 2)) /
            pow(mrho, 6) +
        (0.0625 *
         (-32. * C4 * pow(mrho, 4) * s +
          pow(delta, 2) *
              (1. * pow(mpion, 4) +
               pow(mpion, 2) * (-1.0000000000000009 * pow(mrho, 2) - 2. * s) +
               s * (-3. * pow(mrho, 2) + 1. * s)) +
          delta * (-2. * pow(mpion, 4) +
                   (6. * pow(mrho, 2) + 16. * C4 * pow(mrho, 4) - 2. * s) * s +
                   pow(mpion, 2) * (2. * pow(mrho, 2) + 4. * s))) *
         pow(tmin, 2)) /
            (pow(mrho, 4) * (pow(mpion, 2) - 1. * s)) -
        (0.5625 *
         (C4 * pow(mrho, 6) *
              (2.6666666666666665 + 7.111111111111112 * C4 * pow(mpion, 2) -
               3.555555555555556 * C4 * s) +
          pow(delta, 2) *
              (0.11111111111111112 * pow(mrho, 4) +
               pow(mpion, 2) * (1. * pow(mrho, 2) - 0.22222222222222224 * s) -
               0.22222222222222224 * pow(mrho, 2) * s +
               0.11111111111111112 * pow(s, 2)) +
          delta * pow(mrho, 2) *
              (-2.2222222222222223 * C4 * pow(mrho, 4) +
               pow(mpion, 2) * (-0.6666666666666666 -
                                2.6666666666666665 * C4 * pow(mrho, 2)) +
               0.22222222222222224 * s +
               pow(mrho, 2) *
                   (-0.22222222222222224 + 1.777777777777778 * C4 * s))) *
         pow(tmin, 2)) /
            pow(mrho, 6) +
        (0.03125 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (pow(mrho, 2) * (-2. * eta2 * pow(mpion, 2) -
                          5.999999999999999 * eta1 * pow(mrho, 2) +
                          8. * eta1 * s - 2. * eta2 * s) +
          delta *
              (eta1 * (-5.999999999999999 * pow(mpion, 4) + 5. * pow(mrho, 4) +
                       pow(mpion, 2) * (4. * pow(mrho, 2) - 4. * s) -
                       5.999999999999999 * pow(mrho, 2) * s + 1. * pow(s, 2)) +
               eta2 * (4. * pow(mpion, 4) +
                       pow(mpion, 2) * (1. * pow(mrho, 2) - 2. * s) +
                       s * (5. * pow(mrho, 2) + 2. * s)))) *
         pow(tmin, 2)) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) -
        (0.15625 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (delta *
              (eta1 * (-1.2 * pow(mpion, 4) + 0.6 * pow(mrho, 4) +
                       pow(mpion, 2) * (2. * pow(mrho, 2) - 2.4 * s) -
                       1.6 * pow(mrho, 2) * s + 1. * pow(s, 2)) +
               eta2 * (0.8 * pow(mpion, 4) + (1. * pow(mrho, 2) - 0.4 * s) * s +
                       pow(mpion, 2) * (0.2 * pow(mrho, 2) + 1.2 * s))) +
          pow(mrho, 2) *
              (eta2 * (pow(mpion, 2) * (-0.4 - 6.4 * C4 * s) +
                       s * (-0.4 + 3.2 * C4 * s)) +
               eta1 * (s * (0.8 - 3.2 * C4 * s) +
                       pow(mrho, 2) * (-0.4 + 3.2 * C4 * s) +
                       pow(mpion, 2) *
                           (-0.8 - 3.2 * C4 * pow(mrho, 2) + 6.4 * C4 * s)))) *
         pow(tmin, 2)) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) -
        (0.20833333333333331 * delta *
         (-0.8 * pow(mrho, 2) + 0.8 * C4 * pow(mrho, 4) +
          delta * (0.8 * pow(mpion, 2) + 1. * pow(mrho, 2) - 0.7 * s)) *
         pow(tmin, 3)) /
            pow(mrho, 6) +
        (0.125 *
         (5.333333333333333 * pow(C4, 2) * pow(mrho, 6) +
          delta * (-0.6666666666666666 * pow(mrho, 2) -
                   1.3333333333333333 * C4 * pow(mrho, 4)) +
          pow(delta, 2) *
              (1. * pow(mpion, 2) + 1.1666666666666667 * pow(mrho, 2) -
               0.6666666666666666 * s)) *
         pow(tmin, 3)) /
            pow(mrho, 6) +
        (0.10416666666666666 * delta *
         (-0.8 * pow(mrho, 2) +
          delta * (0.4 * pow(mpion, 2) + 1. * pow(mrho, 2) - 0.6 * s)) *
         pow(tmin, 3)) /
            pow(mrho, 6) +
        (0.020833333333333332 * pow(eta1 - 1. * eta2, 2) * s *
         (-2. * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-1. * pow(mrho, 2) + s)) *
         pow(tmin, 3)) /
            (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
             2. * pow(ma1, 2) * s + pow(s, 2)) +
        (0.10416666666666666 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (0.4 * eta1 * pow(mrho, 2) +
          delta *
              (-0.2 * eta2 * pow(mpion, 2) - 0.2 * eta2 * s +
               eta1 * (-0.4 * pow(mpion, 2) - 0.8 * pow(mrho, 2) + 1. * s))) *
         pow(tmin, 3)) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) -
        (0.14583333333333331 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (delta * (-0.14285714285714285 * eta2 * pow(mpion, 2) -
                   0.42857142857142855 * eta2 * s +
                   eta1 * (-0.2857142857142857 * pow(mpion, 2) -
                           0.5714285714285714 * pow(mrho, 2) + 1. * s)) +
          pow(mrho, 2) *
              (1.1428571428571428 * C4 * eta2 * s +
               eta1 * (0.2857142857142857 - 1.1428571428571428 * C4 * s))) *
         pow(tmin, 3)) /
            (pow(mrho, 2) * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                             2. * pow(ma1, 2) * s + pow(s, 2))) +
        (0. -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(mpion, 4) *
             pow(mrho, 6) +
         0.25 * pow(2. - 1. * delta, 2) * pow(mrho, 10) -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(mrho, 6) *
             pow(s, 2) +
         pow(2. - 1. * delta, 2) * pow(mpion, 2) * pow(mrho, 6) *
             (-1. * pow(mrho, 2) + 1. * s)) /
            (pow(mrho, 6) * (-2. * pow(mpion, 2) + 1. * s + 1. * tmin)) +
        (2. *
         (0. - 2. * pow(mpion, 4) * pow(mrho, 4) - 0.5 * pow(mrho, 8) +
          delta * pow(mrho, 4) *
              (2. * pow(mpion, 4) + 0.5 * pow(mrho, 4) +
               pow(mpion, 2) * (-2. * pow(mrho, 2) - 1.9999999999999998 * s)) +
          pow(mpion, 2) * (2. * pow(mrho, 6) + 2. * pow(mrho, 4) * s) +
          pow(delta, 2) * pow(mrho, 2) *
              (-2.220446049250313e-16 * pow(mpion, 6) - 0.125 * pow(mrho, 6) +
               pow(mpion, 4) *
                   (-0.5 * pow(mrho, 2) + 2.220446049250313e-16 * s) +
               pow(mpion, 2) * (0.5 * pow(mrho, 4) + 0.5 * pow(mrho, 2) * s))) *
         log(abs(-1. * pow(mpion, 2) + 0.5 * s + 0.5 * tmin))) /
            (pow(mrho, 4) * (pow(mpion, 2) - 1. * s)) -
        (0.25 * (eta1 - 1. * eta2) * (pow(ma1, 2) - 1. * s) *
         (eta2 * ((-2. + 1. * delta) * pow(mpion, 6) +
                  (6. - 3. * delta) * pow(mpion, 4) * s +
                  pow(s, 2) * ((4. - 2. * delta) * pow(mrho, 2) +
                               (2. - 1. * delta) * s) +
                  pow(mpion, 2) * s *
                      ((-4. + 2. * delta) * pow(mrho, 2) +
                       (-6. + 3. * delta) * s)) +
          eta1 * ((2. - 1. * delta) * pow(mpion, 6) +
                  (2. - 1. * delta) * pow(mrho, 4) * s +
                  (-2. + 1. * delta) * pow(s, 3) +
                  pow(mpion, 4) * ((4. - 2. * delta) * pow(mrho, 2) +
                                   (-6. + 3. * delta) * s) +
                  pow(mpion, 2) * ((-2. + 1. * delta) * pow(mrho, 4) +
                                   (-4. + 2. * delta) * pow(mrho, 2) * s +
                                   (6. - 3. * delta) * pow(s, 2)))) *
         log(abs(-2. * pow(mpion, 2) + s + tmin))) /
            (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
             2. * pow(ma1, 2) * s + pow(s, 2)) +
        (0.25 *
         (0. + 8. * pow(mpion, 2) * pow(mrho, 6) + 4. * pow(mrho, 8) -
          8. * pow(mrho, 6) * s +
          delta * pow(mrho, 4) *
              (8. * pow(mpion, 4) - 8. * pow(mrho, 4) +
               pow(mpion, 2) * (8. * pow(mrho, 2) - 16. * s) +
               8. * pow(mrho, 2) * s + 8. * pow(s, 2)) +
          pow(delta, 2) * pow(mrho, 4) *
              (-4. * pow(mpion, 4) + 3. * pow(mrho, 4) - 2. * pow(mrho, 2) * s -
               4. * pow(s, 2) +
               pow(mpion, 2) * (-6. * pow(mrho, 2) + 8. * s))) *
         log(abs(-2. * pow(mpion, 2) + 1. * s + 1. * tmin))) /
            pow(mrho, 6) +
        (0.5 *
         (0. + pow(mpion, 2) * (4. * pow(mrho, 6) - 8. * C4 * pow(mrho, 8)) -
          4. * pow(mrho, 6) * s + 8. * C4 * pow(mrho, 8) * s +
          pow(delta, 2) * pow(mrho, 4) *
              (2. * pow(mpion, 4) - 1. * pow(mrho, 4) +
               pow(mpion, 2) * (4. * pow(mrho, 2) - 4. * s) + 2. * pow(s, 2)) +
          delta * pow(mrho, 4) *
              (-4. * pow(mpion, 4) + 2. * pow(mrho, 2) * s - 4. * pow(s, 2) +
               pow(mpion, 2) *
                   (-10. * pow(mrho, 2) + 4. * C4 * pow(mrho, 4) + 8. * s) +
               pow(mrho, 4) * (2. - 4. * C4 * s))) *
         log(abs(-2. * pow(mpion, 2) + 1. * s + 1. * tmin))) /
            pow(mrho, 6))) /
          (16. * Pi *
           (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
            2 * pow(mpion, 2) * (pow(mrho, 2) + s)));

  return xs * to_mb / spin_deg_factor;
}

double PhotonCrossSection<ComputationMethod::Analytic>::
    xs_diff_pi0_rho_pi_rho_mediated(const double s, const double t,
                                    const double m_rho) {
  const double mpion = m_pion_;
  const double mrho = m_rho;
  const double spin_deg_factor = 3.0;

  const double diff_xs =
      (pow(Const, 2) * pow(ghat, 4) *
       ((-0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
          2 * pow(mpion, 2) * (pow(mrho, 2) + s))) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - s, 2)) -
        (0.0625 * (eta1 - eta2) * (-pow(ma1, 2) + s) *
         (2 * pow(mrho, 2) +
          delta * (-2 * pow(mpion, 2) - pow(mrho, 2) + s + t)) *
         (-(eta2 * (s - t) *
            (4 * pow(mpion, 4) + s * (4 * pow(mrho, 2) + s - t) -
             pow(mpion, 2) * (3 * s + t))) +
          eta1 * (8 * pow(mpion, 6) + pow(s, 3) + pow(s, 2) * t +
                  5 * s * pow(t, 2) + pow(t, 3) + 2 * pow(mrho, 4) * (-s + t) +
                  pow(mrho, 2) * (s - 3 * t) * (s + t) +
                  2 * pow(mpion, 2) * (2 * pow(mrho, 2) - s - t) * (s + t) -
                  2 * pow(mpion, 4) * (2 * pow(mrho, 2) + s + 3 * t)))) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (-2 * pow(mpion, 2) + s + t)) -
        (0.0625 *
         pow(-2. * pow(mrho, 2) +
                 delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t),
             2) *
         (8. * pow(mpion, 6) + 4. * pow(mrho, 6) + pow(s, 3) +
          pow(mrho, 4) * (-4. * s - 4. * t) +
          pow(mpion, 4) * (-4. * pow(mrho, 2) - 4. * s - 4. * t) +
          3. * pow(s, 2) * t + 3. * s * pow(t, 2) + pow(t, 3) +
          pow(mrho, 2) * (-3. * pow(s, 2) + 2. * s * t - 3. * pow(t, 2)) +
          pow(mpion, 2) *
              (-8. * pow(mrho, 4) - 2. * pow(s, 2) - 4. * s * t -
               2. * pow(t, 2) + pow(mrho, 2) * (4. * s + 4. * t)))) /
            (pow(mrho, 6) * pow(2. * pow(mpion, 2) - 1. * s - 1. * t, 2)) +
        (0.125 * (-2 + delta) * (eta1 - eta2) * (-pow(ma1, 2) + s) *
         (-(eta2 * (pow(mpion, 2) + s) *
            (-pow(mpion, 4) + pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
             s * (-pow(mrho, 2) + s + 2 * t))) +
          eta1 * (-4 * pow(mpion, 6) +
                  s * (-pow(mrho, 2) + s) * (-pow(mrho, 2) + s + t) +
                  pow(mpion, 4) * (3 * pow(mrho, 2) + s + t) -
                  pow(mpion, 2) * (pow(mrho, 4) + 2 * s * (s - t) +
                                   pow(mrho, 2) * (-s + t))))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (-pow(mpion, 2) + s)) +
        (0.03125 * pow(eta1 - eta2, 2) *
         (pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) * (pow(pow(mrho, 2) + 2 * s, 2) - 2 * s * t) +
               pow(s, 2) * (pow(pow(mrho, 2) + s, 2) +
                            2 * (-pow(mrho, 2) + s) * t + 2 * pow(t, 2)) -
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * (2 * s - t) +
                    2 * s * (s + t))) -
          2 * eta1 * eta2 *
              (pow(mpion, 8) -
               pow(mpion, 4) * (pow(mrho, 4) + 2 * pow(mrho, 2) * s +
                                2 * s * (-2 * s + t)) +
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * (s + t) - 2 * s * (s + t)) +
               pow(s, 2) * (-pow(mrho, 4) + pow(s, 2) - 2 * pow(mrho, 2) * t +
                            2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) * (3 * pow(mrho, 4) + 2 * s * (2 * s - t) +
                                2 * pow(mrho, 2) * (-3 * s + t)) -
               2 * pow(mpion, 2) * (-pow(mrho, 2) + s) *
                   (2 * s * (s + t) - pow(mrho, 2) * (2 * s + t)) +
               s * (-pow(mrho, 2) + s) *
                   (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                    pow(mrho, 2) * (s + 2 * t))))) /
            (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) +
        (0.5 *
         (-2. * pow(mrho, 2) +
          delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t)) *
         (delta * (1. * pow(mpion, 6) + 0.5 * pow(mrho, 6) +
                   0.0625 * pow(s, 3) + pow(mrho, 4) * (-0.5 * s - 0.5 * t) +
                   pow(mpion, 4) * (-0.5 * pow(mrho, 2) - 0.75 * s - 0.25 * t) +
                   0.3125 * pow(s, 2) * t + 0.4375 * s * pow(t, 2) +
                   0.1875 * pow(t, 3) +
                   pow(mpion, 2) * (-1. * pow(mrho, 4) +
                                    pow(mrho, 2) * (0.375 * s + 0.625 * t) +
                                    (-0.5 * s - 0.5 * t) * t) +
                   pow(mrho, 2) * (-0.3125 * pow(s, 2) + 0.25 * s * t -
                                   0.4375 * pow(t, 2))) +
          pow(mrho, 2) *
              (-0.125 * pow(s, 2) + C4 * pow(mrho, 4) * (1. * s - 1. * t) +
               0.125 * pow(t, 2) +
               pow(mpion, 2) * ((0.25 - 1. * C4 * pow(mrho, 2)) * s +
                                (-0.25 + 1. * C4 * pow(mrho, 2)) * t) +
               pow(mrho, 2) * (-0.5 * s + 0.5 * C4 * pow(s, 2) +
                               t * (0.5 - 0.5 * C4 * t))))) /
            (pow(mrho, 6) * (1. * pow(mpion, 2) - 0.5 * s - 0.5 * t)) +
        (pow(delta, 2) *
             (-0.5 * pow(mpion, 6) - 0.0625 * pow(mrho, 6) +
              pow(mpion, 4) * (1. * pow(mrho, 2) + 0.5 * s) +
              pow(mrho, 4) * (-0.125 * s - 0.125 * t) +
              t * (-0.125 * pow(s, 2) - 0.25 * s * t - 0.125 * pow(t, 2)) +
              pow(mpion, 2) * (1.25 * pow(mrho, 4) - 0.125 * pow(s, 2) +
                               pow(mrho, 2) * (-0.875 * s - 1.125 * t) +
                               0.25 * s * t + 0.375 * pow(t, 2)) +
              pow(mrho, 2) *
                  (0.3125 * pow(s, 2) + 0.25 * s * t + 0.4375 * pow(t, 2))) +
         delta * pow(mrho, 2) *
             (pow(mpion, 4) * (-0.5 - 4. * C4 * pow(mrho, 2)) +
              (-0.25 * s - 0.25 * t) * t +
              pow(mrho, 4) * (-0.75 + 0.5 * C4 * s + 2.5 * C4 * t) +
              pow(mrho, 2) * (-1.5 * C4 * pow(s, 2) + s * (1.25 - 2. * C4 * t) +
                              t * (0.25 - 0.5 * C4 * t)) +
              pow(mpion, 2) *
                  (-3. * C4 * pow(mrho, 4) + 0.25 * s + 0.75 * t +
                   pow(mrho, 2) * (-1.5 + 5. * C4 * s + 3. * C4 * t))) +
         pow(mrho, 6) *
             (0.75 +
              C4 * (8. * C4 * pow(mpion, 4) + 2. * C4 * pow(s, 2) +
                    pow(mpion, 2) * (6. - 8. * C4 * s - 8. * C4 * t) +
                    t * (-3. + 2. * C4 * t) + s * (-3. + 4. * C4 * t)))) /
            pow(mrho, 6) +
        (0.0625 * (eta1 - eta2) * (-pow(ma1, 2) + s) *
         (-(eta2 *
            (2 * pow(mpion, 4) *
                 (-4 * C4 * pow(mrho, 2) * (pow(mrho, 2) + 4 * s) +
                  delta * (pow(mrho, 2) + 6 * s - 2 * t)) +
             pow(mpion, 2) *
                 (2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * s) +
                  delta * (-11 * pow(s, 2) - 6 * s * t + pow(t, 2)) +
                  pow(mrho, 2) * ((-10 + delta) * s - (-2 + delta) * t +
                                  32 * C4 * s * (s + t))) +
             s * (-2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * s) +
                  delta * (3 * pow(s, 2) + 2 * s * t + 3 * pow(t, 2)) +
                  pow(mrho, 2) * (3 * (2 + delta) * s + (2 - 5 * delta) * t -
                                  8 * C4 * pow(s + t, 2))))) +
          eta1 *
              (8 * delta * pow(mpion, 6) +
               2 * pow(mpion, 4) *
                   (12 * C4 * pow(mrho, 4) -
                    2 * pow(mrho, 2) * (-1 + 3 * delta + 8 * C4 * s) +
                    3 * delta * (s - t)) +
               delta * (3 * pow(s, 3) + 5 * pow(s, 2) * t + 7 * s * pow(t, 2) +
                        pow(t, 3)) -
               2 * pow(mrho, 2) *
                   ((-3 + 2 * delta) * pow(s, 2) +
                    2 * (-1 + 2 * delta) * s * t +
                    (-1 + 2 * delta) * pow(t, 2) + 4 * C4 * s * pow(s + t, 2)) +
               pow(mrho, 4) * ((-6 + delta) * s + (-2 + 3 * delta) * t +
                               8 * C4 * s * (s + 2 * t)) -
               2 * pow(mpion, 2) *
                   (delta * (s + t) * (5 * s + t) -
                    pow(mrho, 2) * (-6 * s + 9 * delta * s - 2 * t +
                                    5 * delta * t + 16 * C4 * s * (s + t)) +
                    2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * (2 * s + t)))))) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2))) +
        (2 *
         ((0.0625 * (-2. + delta) *
           (-2. * pow(mrho, 2) +
            delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t)) *
           (2. * pow(mpion, 6) + 1. * pow(mrho, 6) +
            pow(mpion, 4) * (-3. * pow(mrho, 2) - 2. * s) +
            pow(mrho, 4) * (-1.5 * s - 1.5 * t) +
            pow(mrho, 2) * (0.5 * s + 0.5 * t) * t +
            s * (0.5 * pow(s, 2) + 1. * s * t + 0.5 * pow(t, 2)) +
            pow(mpion, 2) *
                (-1. * pow(mrho, 4) - 0.5 * pow(s, 2) - 1. * s * t -
                 0.5 * pow(t, 2) + pow(mrho, 2) * (-0.5 * s + 2.5 * t)))) /
              ((pow(mpion, 2) - 1. * s) *
               (1. * pow(mpion, 2) - 0.5 * s - 0.5 * t)) +
          (0.0625 * (-2 + delta) *
           (delta *
                (6 * pow(mpion, 6) -
                 pow(mpion, 4) * (9 * (pow(mrho, 2) + s) + t) -
                 pow(mpion, 2) * (2 * pow(mrho, 4) - 3 * pow(s, 2) +
                                  pow(mrho, 2) * (5 * s - 7 * t) + pow(t, 2)) +
                 (pow(mrho, 2) - s - t) *
                     (3 * pow(mrho, 4) - s * t - pow(mrho, 2) * (3 * s + t))) +
            2 * pow(mrho, 2) *
                (pow(mpion, 4) * (1 + 4 * C4 * pow(mrho, 2)) +
                 pow(mrho, 4) * (-1 + 4 * C4 * s) + s * t -
                 pow(mpion, 2) * (4 * C4 * pow(mrho, 4) + s -
                                  2 * pow(mrho, 2) * (1 + 4 * C4 * s) + t) +
                 pow(mrho, 2) * (t + s * (1 - 4 * C4 * (s + 2 * t)))))) /
              (-pow(mpion, 2) + s))) /
            pow(mrho, 4))) /
      (16. * Pi *
       (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
        2 * pow(mpion, 2) * (pow(mrho, 2) + s)));

  return to_mb * diff_xs / spin_deg_factor;
}

// C13
double
PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_rho_pi0_rho_mediated(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_rho, m_pion_, 0.);
  const double &tmax = t_mandelstam[0];
  const double &tmin = t_mandelstam[1];
  const double spin_deg_factor = 3.0;

  const double xs =
      (pow(Const, 2) * pow(ghat, 4) *
       (0. +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta2, 2) *
              (1. * pow(ma1, 8) + 1. * pow(mpion, 8) -
               2. * pow(mpion, 6) * pow(mrho, 2) +
               1. * pow(mpion, 4) * pow(mrho, 4) +
               pow(ma1, 6) *
                   (-4. * pow(mpion, 2) + 2. * pow(mrho, 2) + 2. * s) +
               pow(ma1, 2) * pow(mpion, 2) *
                   (-2. * pow(mrho, 4) +
                    pow(mpion, 2) * (4. * pow(mrho, 2) - 2. * s) +
                    2. * pow(mrho, 2) * s) +
               pow(ma1, 4) * (4. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                              pow(mpion, 2) * (-4. * pow(mrho, 2) - 4. * s) -
                              2. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
          eta1 * eta2 *
              (-2. * pow(ma1, 8) - 2. * pow(mpion, 8) +
               2. * pow(mpion, 4) * pow(mrho, 4) +
               pow(ma1, 6) * (8. * pow(mpion, 2) - 4. * s) +
               pow(ma1, 2) * pow(mpion, 2) *
                   (-4. * pow(mrho, 4) - 4. * pow(mrho, 2) * s +
                    pow(mpion, 2) * (4. * pow(mrho, 2) + 4. * s)) +
               pow(ma1, 4) * (-8. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                              4. * pow(mrho, 2) * s - 4. * pow(s, 2) +
                              pow(mpion, 2) * (-4. * pow(mrho, 2) + 8. * s))) +
          pow(eta1, 2) *
              (1. * pow(ma1, 8) + 1. * pow(mpion, 8) -
               2. * pow(mpion, 6) * pow(mrho, 2) -
               2. * pow(mpion, 2) * pow(mrho, 4) * s +
               pow(ma1, 6) *
                   (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
               pow(mpion, 4) * (3. * pow(mrho, 4) + 2. * pow(mrho, 2) * s) +
               pow(ma1, 4) * (4. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                              pow(mpion, 2) * (8. * pow(mrho, 2) - 4. * s) -
                              4. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
               pow(ma1, 2) * (pow(mpion, 4) * (-6. * pow(mrho, 2) - 2. * s) +
                              pow(mrho, 2) * (2. * pow(mrho, 2) - 2. * s) * s +
                              pow(mpion, 2) * (-4. * pow(mrho, 4) +
                                               6. * pow(mrho, 2) * s))))) /
            (1. * pow(ma1, 2) - 1. * tmax) +
        (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
         (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
            (1. * pow(mpion, 2) - 1. * tmax) -
        (0.25 * pow(-2. + delta, 2) * pow(mpion, 2) * tmax) / pow(mrho, 2) +
        0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
            (eta2 * (-1. * pow(ma1, 2) + pow(mrho, 2) - 2. * s) +
             eta1 *
                 (pow(ma1, 2) - 1. * pow(mpion, 2) - 2. * pow(mrho, 2) + s)) *
            tmax +
        0.03125 * pow(eta1 - 1. * eta2, 2) *
            (pow(eta1, 2) *
                 (3. * pow(ma1, 4) + 4. * pow(mpion, 4) + pow(mrho, 4) +
                  pow(mpion, 2) * (8. * pow(mrho, 2) - 4. * s) -
                  4. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                  pow(ma1, 2) *
                      (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) +
             pow(eta2, 2) *
                 (3. * pow(ma1, 4) + 4. * pow(mpion, 4) + pow(mrho, 4) +
                  pow(mpion, 2) * (-4. * pow(mrho, 2) - 4. * s) -
                  2. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                  pow(ma1, 2) *
                      (-8. * pow(mpion, 2) + 4. * pow(mrho, 2) + 4. * s)) +
             eta1 * eta2 *
                 (-6. * pow(ma1, 4) - 8. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                  pow(ma1, 2) * (16. * pow(mpion, 2) - 8. * s) +
                  4. * pow(mrho, 2) * s - 4. * pow(s, 2) +
                  pow(mpion, 2) * (-4. * pow(mrho, 2) + 8. * s))) *
            tmax +
        (2. *
         (0. - 0.25 * pow(mrho, 4) - 1. * C4 * pow(mrho, 6) +
          pow(mpion, 2) * (0.75 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
          2. * C4 * pow(mrho, 4) * s +
          pow(delta, 2) *
              (-0.125 * pow(mpion, 4) - 0.1875 * pow(mrho, 4) +
               pow(mpion, 2) * (0.0625 * pow(mrho, 2) + 0.0625 * s) +
               0.1875 * pow(mrho, 2) * s) +
          delta *
              (0.25 * pow(mpion, 4) + 0.5 * C4 * pow(mrho, 6) +
               pow(mpion, 2) *
                   (-0.5 * pow(mrho, 2) + 0.5 * C4 * pow(mrho, 4) - 0.125 * s) -
               0.375 * pow(mrho, 2) * s + pow(mrho, 4) * (0.5 - 1. * C4 * s))) *
         tmax) /
            pow(mrho, 4) -
        (0.0625 *
         (0. + 8. * pow(mpion, 2) * pow(mrho, 4) - 12. * pow(mrho, 6) +
          4. * pow(mrho, 4) * s +
          delta * pow(mrho, 2) *
              (-16. * pow(mpion, 4) - 16. * pow(mpion, 2) * pow(mrho, 2) -
               4. * pow(mrho, 4) + 16. * pow(mrho, 2) * s + 4. * pow(s, 2)) +
          pow(delta, 2) *
              (8. * pow(mpion, 6) + 9. * pow(mrho, 6) +
               pow(mpion, 4) * (4. * pow(mrho, 2) - 4. * s) -
               13. * pow(mrho, 4) * s - 5. * pow(mrho, 2) * pow(s, 2) +
               1. * pow(s, 3) +
               pow(mpion, 2) * (-2. * pow(mrho, 4) + 4. * pow(mrho, 2) * s -
                                2. * pow(s, 2)))) *
         tmax) /
            pow(mrho, 6) -
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (pow(mrho, 2) * (eta1 * (2. * pow(ma1, 2) + 2. * pow(mrho, 2)) +
                          eta2 * (-2. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                                  8. * pow(mrho, 2) + 6. * s)) +
          delta *
              (eta1 *
                   (1. * pow(ma1, 4) - 2. * pow(mpion, 4) - 3. * pow(mrho, 4) +
                    pow(mpion, 2) * (4. * pow(mrho, 2) - 4. * s) -
                    2. * pow(mrho, 2) * s + 5.000000000000001 * pow(s, 2) +
                    pow(ma1, 2) * (-2. * pow(mpion, 2) + 1. * s)) +
               eta2 *
                   (-1. * pow(ma1, 4) - 4. * pow(mpion, 4) + 4. * pow(mrho, 4) +
                    pow(mpion, 2) * (-1. * pow(mrho, 2) - 2. * s) +
                    1. * pow(mrho, 2) * s - 1. * pow(s, 2) +
                    pow(ma1, 2) *
                        (3. * pow(mpion, 2) - 3. * pow(mrho, 2) + 2. * s)))) *
         tmax) /
            pow(mrho, 2) -
        (0.5 *
         (pow(mrho, 6) *
              (-1.5 + C4 * (-12. * pow(mpion, 2) + 6. * s) +
               pow(C4, 2) * (-16. * pow(mpion, 4) + 16. * pow(mpion, 2) * s -
                             4. * pow(s, 2))) +
          pow(delta, 2) *
              (1. * pow(mpion, 6) - 2. * pow(mpion, 4) * pow(mrho, 2) +
               0.125 * pow(mrho, 6) + 0.25 * pow(mrho, 4) * s -
               0.875 * pow(mrho, 2) * pow(s, 2) + 0.25 * pow(s, 3) +
               pow(mpion, 2) * (-2.5 * pow(mrho, 4) + 2.25 * pow(mrho, 2) * s -
                                0.75 * pow(s, 2))) +
          delta * pow(mrho, 2) *
              (pow(mpion, 4) * (1. + 8. * C4 * pow(mrho, 2)) + 0.5 * pow(s, 2) +
               pow(mrho, 4) * (1.5 - 5. * C4 * s) +
               pow(mrho, 2) * s * (-0.5 + 1. * C4 * s) +
               pow(mpion, 2) * (6. * C4 * pow(mrho, 4) - 1.5 * s +
                                pow(mrho, 2) * (3. - 6. * C4 * s)))) *
         tmax) /
            pow(mrho, 6) -
        (0.5 *
         (0. - 4. * C4 * pow(mrho, 8) - 0.5 * pow(mrho, 4) * s +
          pow(mrho, 6) * (2. + 2. * C4 * s) +
          pow(delta, 2) *
              (-2. * pow(mpion, 6) - 2. * pow(mrho, 6) +
               0.5 * pow(mpion, 4) * s + 2.125 * pow(mrho, 4) * s +
               1.25 * pow(mrho, 2) * pow(s, 2) - 0.375 * pow(s, 3) +
               pow(mpion, 2) * (1.5 * pow(mrho, 4) - 1.5 * pow(mrho, 2) * s +
                                1. * pow(s, 2))) +
          delta * pow(mrho, 2) *
              (2. * pow(mpion, 4) + 2. * C4 * pow(mrho, 6) - 1. * pow(s, 2) +
               pow(mrho, 2) * s * (-3. + 1. * C4 * s) +
               pow(mrho, 4) * (1. + 1. * C4 * s) +
               pow(mpion, 2) * (1. * s + pow(mrho, 2) * (1. - 2. * C4 * s)))) *
         tmax) /
            pow(mrho, 6) +
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (delta *
              (eta1 * (3. * pow(ma1, 4) + 6. * pow(mpion, 4) + pow(mrho, 4) +
                       pow(mpion, 2) * (18. * pow(mrho, 2) - 12. * s) -
                       8. * pow(mrho, 2) * s + 7. * pow(s, 2) +
                       pow(ma1, 2) * (-10. * pow(mpion, 2) - 4. * pow(mrho, 2) +
                                      5. * s)) +
               eta2 * (-3. * pow(ma1, 4) - 12. * pow(mpion, 4) +
                       2. * pow(mrho, 4) +
                       pow(ma1, 2) *
                           (11. * pow(mpion, 2) - 3. * pow(mrho, 2) - 2. * s) +
                       5. * pow(mrho, 2) * s - 3. * pow(s, 2) +
                       pow(mpion, 2) * (-1. * pow(mrho, 2) + 6. * s))) +
          pow(mrho, 2) *
              (eta1 *
                   (-8. * C4 * pow(ma1, 4) - 32. * C4 * pow(mpion, 4) -
                    6. * pow(mrho, 2) +
                    pow(ma1, 2) * (6. + C4 * (32. * pow(mpion, 2) +
                                              8. * pow(mrho, 2) - 16. * s)) +
                    4. * s + 16. * C4 * pow(mrho, 2) * s - 8. * C4 * pow(s, 2) +
                    pow(mpion, 2) *
                        (-12. - 32. * C4 * pow(mrho, 2) + 32. * C4 * s)) +
               eta2 * (8. * C4 * pow(ma1, 4) + 32. * C4 * pow(mpion, 4) -
                       4. * pow(mrho, 2) - 2. * s + 8. * C4 * pow(s, 2) +
                       pow(mpion, 2) *
                           (10. - 16. * C4 * pow(mrho, 2) - 32. * C4 * s) +
                       pow(ma1, 2) *
                           (-6. + C4 * (-32. * pow(mpion, 2) +
                                        8. * pow(mrho, 2) + 16. * s))))) *
         tmax) /
            pow(mrho, 2) +
        0.0625 * (-2. + delta) * pow(eta1 - 1. * eta2, 2) * pow(tmax, 2) +
        (0.1875 *
         (1.3333333333333333 * pow(mrho, 2) +
          5.333333333333333 * C4 * pow(mrho, 4) +
          pow(delta, 2) *
              (1. * pow(mpion, 2) + 1.3333333333333333 * pow(mrho, 2) -
               0.3333333333333333 * s) +
          delta * (-2. * pow(mpion, 2) - 3.3333333333333335 * pow(mrho, 2) -
                   2.6666666666666665 * C4 * pow(mrho, 4) +
                   0.6666666666666666 * s)) *
         pow(tmax, 2)) /
            pow(mrho, 4) +
        0.03125 * pow(eta1 - 1. * eta2, 3) *
            (eta2 * (-1. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                     1. * pow(mrho, 2) - 1. * s) +
             eta1 *
                 (pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
            pow(tmax, 2) -
        (0.375 *
         (0.3333333333333333 * pow(mrho, 4) -
          1.3333333333333333 * C4 * pow(mrho, 6) +
          delta * pow(mrho, 2) *
              (1.3333333333333333 * pow(mrho, 2) -
               0.6666666666666666 * C4 * pow(mrho, 4) +
               pow(mpion, 2) * (-0.6666666666666666 +
                                1.3333333333333333 * C4 * pow(mrho, 2)) -
               0.6666666666666666 * s) +
          pow(delta, 2) * (1. * pow(mpion, 4) + 0.25 * pow(mrho, 4) +
                           pow(mpion, 2) * (-0.3333333333333333 * pow(mrho, 2) +
                                            0.6666666666666666 * s) -
                           0.5833333333333334 * pow(s, 2))) *
         pow(tmax, 2)) /
            pow(mrho, 6) -
        (0.03125 * (1. * eta1 - 1. * eta2) *
         ((2. * eta1 - 2. * eta2) * pow(mrho, 2) +
          delta * (eta1 * (1. * pow(ma1, 2) - 2. * pow(mpion, 2) + 1. * s) +
                   eta2 * (-1. * pow(ma1, 2) + 3. * pow(mpion, 2) -
                           3. * pow(mrho, 2) + 2. * s))) *
         pow(tmax, 2)) /
            pow(mrho, 2) +
        (0.03125 *
         (0. - 4. * pow(mrho, 4) +
          delta * (16. * pow(mrho, 4) - 8. * pow(mrho, 2) * s) +
          pow(delta, 2) * (4. * pow(mpion, 4) - 3. * pow(mrho, 4) +
                           2. * pow(mrho, 2) * s - 3. * pow(s, 2) +
                           pow(mpion, 2) * (-4. * pow(mrho, 2) + 4. * s))) *
         pow(tmax, 2)) /
            pow(mrho, 6) +
        (0.25 *
         (C4 * pow(mrho, 6) * (-6. - 16. * C4 * pow(mpion, 2) + 8. * C4 * s) +
          pow(delta, 2) * (1. * pow(mpion, 4) - 0.25 * pow(mrho, 4) +
                           pow(mpion, 2) * (-1.75 * pow(mrho, 2) + 0.5 * s) +
                           0.5 * pow(mrho, 2) * s - 0.5 * pow(s, 2)) +
          delta * pow(mrho, 2) *
              (1. * C4 * pow(mrho, 4) +
               pow(mpion, 2) * (0.5 + 10. * C4 * pow(mrho, 2)) - 0.5 * s +
               pow(mrho, 2) * (2.5 - 4. * C4 * s))) *
         pow(tmax, 2)) /
            pow(mrho, 6) +
        (0.09375 * (1. * eta1 - 1. * eta2) *
         (delta *
              (eta2 * (-1. * pow(ma1, 2) + 3.6666666666666665 * pow(mpion, 2) -
                       1. * pow(mrho, 2) - 0.6666666666666666 * s) +
               eta1 * (1. * pow(ma1, 2) - 3.3333333333333335 * pow(mpion, 2) -
                       1.3333333333333333 * pow(mrho, 2) +
                       1.6666666666666667 * s)) +
          pow(mrho, 2) *
              (eta1 * (2. + C4 * (-2.6666666666666665 * pow(ma1, 2) +
                                  10.666666666666666 * pow(mpion, 2) +
                                  2.6666666666666665 * pow(mrho, 2) -
                                  5.333333333333333 * s)) +
               eta2 * (-2. + C4 * (2.6666666666666665 * pow(ma1, 2) -
                                   10.666666666666666 * pow(mpion, 2) +
                                   2.6666666666666665 * pow(mrho, 2) +
                                   5.333333333333333 * s)))) *
         pow(tmax, 2)) /
            pow(mrho, 2) +
        0.010416666666666666 * pow(eta1 - 1. * eta2, 4) * pow(tmax, 3) -
        (0.041666666666666664 * delta * (-2. + 1. * delta) * pow(tmax, 3)) /
            pow(mrho, 4) -
        (0.020833333333333332 * delta * pow(1. * eta1 - 1. * eta2, 2) *
         pow(tmax, 3)) /
            pow(mrho, 2) -
        (0.16666666666666666 * pow(1. * eta1 - 1. * eta2, 2) *
         (-0.375 * delta + 1. * C4 * pow(mrho, 2)) * pow(tmax, 3)) /
            pow(mrho, 2) +
        (0.10416666666666666 * delta *
         (-0.8 * pow(mrho, 2) +
          delta * (0.4 * pow(mpion, 2) + 1. * pow(mrho, 2) - 0.6 * s)) *
         pow(tmax, 3)) /
            pow(mrho, 6) +
        (0.16666666666666666 * delta *
         (0. - 0.75 * delta * pow(mrho, 2) + 1. * C4 * pow(mrho, 4) +
          0.625 * delta * s) *
         pow(tmax, 3)) /
            pow(mrho, 6) -
        (0.041666666666666664 *
         (12. * C4 * delta * pow(mrho, 4) - 16. * pow(C4, 2) * pow(mrho, 6) +
          pow(delta, 2) * (1. * pow(mpion, 2) - 2.5 * pow(mrho, 2) + 1. * s)) *
         pow(tmax, 3)) /
            pow(mrho, 6) +
        (0. -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(mpion, 4) *
             pow(mrho, 6) +
         0.25 * pow(2. - 1. * delta, 2) * pow(mrho, 10) -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(mrho, 6) *
             pow(s, 2) +
         pow(2. - 1. * delta, 2) * pow(mpion, 2) * pow(mrho, 6) *
             (-1. * pow(mrho, 2) + 1. * s)) /
            (pow(mrho, 6) * (-2. * pow(mpion, 2) + 1. * s + 1. * tmax)) -
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (2. * pow(mrho, 2) + delta * (1. * pow(ma1, 2) - 2. * pow(mpion, 2) -
                                       1. * pow(mrho, 2) + 1. * s)) *
         (eta2 *
              (-1. * pow(ma1, 6) +
               pow(mpion, 2) * (4. * pow(mpion, 2) - 1. * s) * s +
               pow(ma1, 4) * (3. * pow(mpion, 2) - 4. * pow(mrho, 2) + 2. * s) +
               pow(ma1, 2) * (-4. * pow(mpion, 4) - 2. * pow(mpion, 2) * s +
                              (4. * pow(mrho, 2) - 1. * s) * s)) +
          eta1 * (1. * pow(ma1, 6) + 8. * pow(mpion, 6) +
                  pow(mpion, 4) * (-4. * pow(mrho, 2) - 6. * s) +
                  pow(mpion, 2) * (4. * pow(mrho, 2) - 2. * s) * s +
                  pow(ma1, 4) *
                      (-2. * pow(mpion, 2) + 1. * pow(mrho, 2) + 1. * s) +
                  s * (2. * pow(mrho, 4) - 3. * pow(mrho, 2) * s +
                       1. * pow(s, 2)) +
                  pow(ma1, 2) * (-2. * pow(mpion, 4) - 2. * pow(mrho, 4) +
                                 pow(mpion, 2) * (4. * pow(mrho, 2) - 4. * s) -
                                 2. * pow(mrho, 2) * s + 5. * pow(s, 2)))) *
         log(abs(-pow(ma1, 2) + tmax))) /
            (pow(mrho, 2) * (pow(ma1, 2) - 2. * pow(mpion, 2) + s)) +
        (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
         (eta2 * (-1. * pow(ma1, 6) + pow(mpion, 6) -
                  1. * pow(mpion, 4) * pow(mrho, 2) +
                  pow(ma1, 4) * (pow(mpion, 2) + pow(mrho, 2) - 2. * s) +
                  pow(ma1, 2) * (3. * pow(mpion, 4) - 2. * pow(mpion, 2) * s)) +
          eta1 * (pow(ma1, 6) +
                  pow(ma1, 4) * (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + s) +
                  pow(mpion, 2) * (-4. * pow(mpion, 4) - 1. * pow(mrho, 4) -
                                   1. * pow(mrho, 2) * s +
                                   pow(mpion, 2) * (3. * pow(mrho, 2) + s)) +
                  pow(ma1, 2) *
                      (pow(mpion, 4) + pow(mrho, 4) - 1. * pow(mrho, 2) * s +
                       pow(mpion, 2) * (pow(mrho, 2) + 2. * s)))) *
         log(abs(-pow(ma1, 2) + tmax))) /
            (pow(ma1, 2) - 1. * pow(mpion, 2)) +
        0.0625 * pow(eta1 - 1. * eta2, 2) *
            (pow(eta1, 2) *
                 (2. * pow(ma1, 6) +
                  pow(mpion, 4) * (-3. * pow(mrho, 2) - 1. * s) +
                  pow(mrho, 2) * (pow(mrho, 2) - 1. * s) * s +
                  pow(ma1, 4) *
                      (-6. * pow(mpion, 2) - 3. * pow(mrho, 2) + 3. * s) +
                  pow(mpion, 2) * (-2. * pow(mrho, 4) + 3. * pow(mrho, 2) * s) +
                  pow(ma1, 2) * (4. * pow(mpion, 4) + pow(mrho, 4) +
                                 pow(mpion, 2) * (8. * pow(mrho, 2) - 4. * s) -
                                 4. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
             pow(eta2, 2) *
                 (2. * pow(ma1, 6) +
                  pow(ma1, 4) *
                      (-6. * pow(mpion, 2) + 3. * pow(mrho, 2) + 3. * s) +
                  pow(mpion, 2) *
                      (-1. * pow(mrho, 4) +
                       pow(mpion, 2) * (2. * pow(mrho, 2) - 1. * s) +
                       pow(mrho, 2) * s) +
                  pow(ma1, 2) * (4. * pow(mpion, 4) + pow(mrho, 4) +
                                 pow(mpion, 2) * (-4. * pow(mrho, 2) - 4. * s) -
                                 2. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
             eta1 * eta2 *
                 (-4. * pow(ma1, 6) +
                  pow(ma1, 4) * (12. * pow(mpion, 2) - 6. * s) +
                  pow(mpion, 2) *
                      (-2. * pow(mrho, 4) - 2. * pow(mrho, 2) * s +
                       pow(mpion, 2) * (2. * pow(mrho, 2) + 2. * s)) +
                  pow(ma1, 2) *
                      (-8. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                       4. * pow(mrho, 2) * s - 4. * pow(s, 2) +
                       pow(mpion, 2) * (-4. * pow(mrho, 2) + 8. * s)))) *
            log(abs(-pow(ma1, 2) + tmax)) +
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (delta *
              (eta1 *
                   (3. * pow(ma1, 6) + 8. * pow(mpion, 6) +
                    pow(mpion, 4) * (-12. * pow(mrho, 2) - 6. * s) +
                    pow(ma1, 4) *
                        (-10. * pow(mpion, 2) - 4. * pow(mrho, 2) + 5. * s) +
                    pow(mpion, 2) * (-4. * pow(mrho, 4) +
                                     10. * pow(mrho, 2) * s - 2. * pow(s, 2)) +
                    s * (3. * pow(mrho, 4) - 4. * pow(mrho, 2) * s +
                         pow(s, 2)) +
                    pow(ma1, 2) *
                        (6. * pow(mpion, 4) + pow(mrho, 4) +
                         pow(mpion, 2) * (18. * pow(mrho, 2) - 12. * s) -
                         8. * pow(mrho, 2) * s + 7. * pow(s, 2))) +
               eta2 * (-3. * pow(ma1, 6) +
                       pow(ma1, 4) *
                           (11. * pow(mpion, 2) - 3. * pow(mrho, 2) - 2. * s) +
                       pow(mpion, 2) *
                           (-2. * pow(mrho, 4) + pow(mrho, 2) * s -
                            1. * pow(s, 2) +
                            pow(mpion, 2) * (-2. * pow(mrho, 2) + 4. * s)) +
                       pow(ma1, 2) *
                           (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                            5. * pow(mrho, 2) * s - 3. * pow(s, 2) +
                            pow(mpion, 2) * (-1. * pow(mrho, 2) + 6. * s)))) +
          pow(mrho, 2) *
              (eta2 *
                   (8. * C4 * pow(ma1, 6) +
                    pow(mpion, 2) *
                        ((4. + 8. * C4 * pow(mpion, 2)) * pow(mrho, 2) -
                         2. * s) +
                    pow(ma1, 4) * (-6. + C4 * (-32. * pow(mpion, 2) +
                                               8. * pow(mrho, 2) + 16. * s)) +
                    pow(ma1, 2) *
                        (32. * C4 * pow(mpion, 4) - 4. * pow(mrho, 2) +
                         pow(mpion, 2) *
                             (10. - 16. * C4 * pow(mrho, 2) - 32. * C4 * s) +
                         s * (-2. + 8. * C4 * s))) +
               eta1 * (-8. * C4 * pow(ma1, 6) +
                       pow(mpion, 4) * (4. + 24. * C4 * pow(mrho, 2)) +
                       pow(ma1, 4) * (6. + C4 * (32. * pow(mpion, 2) +
                                                 8. * pow(mrho, 2) - 16. * s)) +
                       s * (-2. * pow(mrho, 2) + 2. * s) +
                       pow(mpion, 2) *
                           (-4. * s + pow(mrho, 2) * (8. - 16. * C4 * s)) +
                       pow(ma1, 2) *
                           (-32. * C4 * pow(mpion, 4) + s * (4. - 8. * C4 * s) +
                            pow(mrho, 2) * (-6. + 16. * C4 * s) +
                            pow(mpion, 2) * (-12. - 32. * C4 * pow(mrho, 2) +
                                             32. * C4 * s))))) *
         log(abs(-pow(ma1, 2) + tmax))) /
            pow(mrho, 2) +
        0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
            log(abs(-pow(mpion, 2) + tmax)) -
        (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
         (eta2 * (2. * pow(mpion, 6) - 2. * pow(mpion, 4) * s) +
          eta1 * (-2. * pow(mpion, 6) - 1. * pow(mpion, 2) * pow(mrho, 2) * s +
                  pow(mpion, 4) * (pow(mrho, 2) + 2. * s))) *
         log(abs(-pow(mpion, 2) + tmax))) /
            (pow(ma1, 2) - 1. * pow(mpion, 2)) -
        (0.125 *
         (0. - 32. * C4 * pow(mpion, 6) * pow(mrho, 4) - 8. * pow(mrho, 8) +
          8. * pow(mrho, 6) * s +
          pow(mpion, 4) * pow(mrho, 4) * (16. + 64. * C4 * s) +
          pow(mpion, 2) * pow(mrho, 4) *
              (24. * pow(mrho, 2) + s * (-16. - 32. * C4 * s)) +
          pow(delta, 2) * pow(mrho, 2) *
              (-4. * pow(mpion, 6) - 2. * pow(mrho, 6) + 2. * pow(mrho, 4) * s +
               pow(mpion, 4) * (4. * pow(mrho, 2) + 8. * s) +
               pow(mpion, 2) * (6. * pow(mrho, 4) - 4. * pow(mrho, 2) * s -
                                4.000000000000001 * pow(s, 2))) +
          delta * pow(mrho, 2) *
              (8. * pow(mrho, 6) +
               pow(mpion, 6) * (8. + 16. * C4 * pow(mrho, 2)) -
               8. * pow(mrho, 4) * s +
               pow(mpion, 4) * (-16. * s + pow(mrho, 2) * (-15.999999999999996 -
                                                           32. * C4 * s)) +
               pow(mpion, 2) * (-24. * pow(mrho, 4) + 8. * pow(s, 2) +
                                pow(mrho, 2) * s * (16. + 16. * C4 * s)))) *
         log(abs(-pow(mpion, 2) + tmax))) /
            (pow(mrho, 4) * (-1. * pow(mpion, 2) + 1. * s)) -
        (0.25 * (1. * eta1 - 1. * eta2) *
         (eta2 * ((2. - 1. * delta) * pow(mpion, 6) +
                  pow(mpion, 2) * s *
                      ((-12. + 6. * delta) * pow(mrho, 2) +
                       (6. - 3. * delta) * s) +
                  pow(s, 2) * ((4. - 2. * delta) * pow(mrho, 2) +
                               (-2. + 1. * delta) * s) +
                  pow(mpion, 4) * ((8. - 4. * delta) * pow(mrho, 2) +
                                   (-6. + 3. * delta) * s)) +
          eta1 * ((-2. + 1. * delta) * pow(mpion, 6) +
                  (-2. + 1. * delta) * pow(mrho, 4) * s +
                  (2. - 1. * delta) * pow(s, 3) +
                  pow(mpion, 4) * ((-4. + 2. * delta) * pow(mrho, 2) +
                                   (6. - 3. * delta) * s) +
                  pow(mpion, 2) * ((2. - 1. * delta) * pow(mrho, 4) +
                                   (4. - 2. * delta) * pow(mrho, 2) * s +
                                   (-6. + 3. * delta) * pow(s, 2)))) *
         log(abs(-2. * pow(mpion, 2) + s + tmax))) /
            (pow(ma1, 2) - 2. * pow(mpion, 2) + s) +
        (0.125 *
         (0. +
          (32. - 31.999999999999993 * delta + 8. * pow(delta, 2)) *
              pow(mpion, 4) * pow(mrho, 4) -
          2.0000000000000004 * pow(2. - 1. * delta, 2) * pow(mrho, 8) +
          pow(mpion, 2) * pow(mrho, 4) *
              (8.000000000000002 * pow(2. - 1. * delta, 2) * pow(mrho, 2) +
               (-32. + 31.999999999999996 * delta - 8. * pow(delta, 2)) * s)) *
         log(abs(-2. * pow(mpion, 2) + s + tmax))) /
            (pow(mrho, 4) * (-1. * pow(mpion, 2) + 1. * s)) +
        (0.25 *
         (0. + 8. * pow(mpion, 2) * pow(mrho, 6) + 4. * pow(mrho, 8) -
          8. * pow(mrho, 6) * s +
          delta * pow(mrho, 4) *
              (8. * pow(mpion, 4) - 8. * pow(mrho, 4) +
               pow(mpion, 2) * (8. * pow(mrho, 2) - 16. * s) +
               8. * pow(mrho, 2) * s + 8. * pow(s, 2)) +
          pow(delta, 2) * pow(mrho, 4) *
              (-4. * pow(mpion, 4) + 3. * pow(mrho, 4) - 2. * pow(mrho, 2) * s -
               4. * pow(s, 2) +
               pow(mpion, 2) * (-6. * pow(mrho, 2) + 8. * s))) *
         log(abs(-2. * pow(mpion, 2) + s + tmax))) /
            pow(mrho, 6) -
        (0.5 *
         (0. + pow(mpion, 2) * (4. * pow(mrho, 6) - 8. * C4 * pow(mrho, 8)) -
          4. * pow(mrho, 6) * s + 8. * C4 * pow(mrho, 8) * s +
          pow(delta, 2) * pow(mrho, 4) *
              (-2. * pow(mpion, 4) + 1. * pow(mrho, 4) - 2. * pow(s, 2) +
               pow(mpion, 2) * (-4. * pow(mrho, 2) + 4. * s)) +
          delta * pow(mrho, 4) *
              (4. * pow(mpion, 4) +
               pow(mpion, 2) *
                   (6. * pow(mrho, 2) + 4. * C4 * pow(mrho, 4) - 8. * s) +
               2. * pow(mrho, 2) * s + 4. * pow(s, 2) +
               pow(mrho, 4) * (-2. - 4. * C4 * s))) *
         log(abs(-2. * pow(mpion, 2) + s + tmax))) /
            pow(mrho, 6))) /
          (16. * Pi *
           (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
            2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
      (pow(Const, 2) * pow(ghat, 4) *
       (0. +
        (0.03125 * pow(eta1 - 1. * eta2, 2) *
         (pow(eta2, 2) *
              (1. * pow(ma1, 8) + 1. * pow(mpion, 8) -
               2. * pow(mpion, 6) * pow(mrho, 2) +
               1. * pow(mpion, 4) * pow(mrho, 4) +
               pow(ma1, 6) *
                   (-4. * pow(mpion, 2) + 2. * pow(mrho, 2) + 2. * s) +
               pow(ma1, 2) * pow(mpion, 2) *
                   (-2. * pow(mrho, 4) +
                    pow(mpion, 2) * (4. * pow(mrho, 2) - 2. * s) +
                    2. * pow(mrho, 2) * s) +
               pow(ma1, 4) * (4. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                              pow(mpion, 2) * (-4. * pow(mrho, 2) - 4. * s) -
                              2. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
          eta1 * eta2 *
              (-2. * pow(ma1, 8) - 2. * pow(mpion, 8) +
               2. * pow(mpion, 4) * pow(mrho, 4) +
               pow(ma1, 6) * (8. * pow(mpion, 2) - 4. * s) +
               pow(ma1, 2) * pow(mpion, 2) *
                   (-4. * pow(mrho, 4) - 4. * pow(mrho, 2) * s +
                    pow(mpion, 2) * (4. * pow(mrho, 2) + 4. * s)) +
               pow(ma1, 4) * (-8. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                              4. * pow(mrho, 2) * s - 4. * pow(s, 2) +
                              pow(mpion, 2) * (-4. * pow(mrho, 2) + 8. * s))) +
          pow(eta1, 2) *
              (1. * pow(ma1, 8) + 1. * pow(mpion, 8) -
               2. * pow(mpion, 6) * pow(mrho, 2) -
               2. * pow(mpion, 2) * pow(mrho, 4) * s +
               pow(ma1, 6) *
                   (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
               pow(mpion, 4) * (3. * pow(mrho, 4) + 2. * pow(mrho, 2) * s) +
               pow(ma1, 4) * (4. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                              pow(mpion, 2) * (8. * pow(mrho, 2) - 4. * s) -
                              4. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
               pow(ma1, 2) * (pow(mpion, 4) * (-6. * pow(mrho, 2) - 2. * s) +
                              pow(mrho, 2) * (2. * pow(mrho, 2) - 2. * s) * s +
                              pow(mpion, 2) * (-4. * pow(mrho, 4) +
                                               6. * pow(mrho, 2) * s))))) /
            (1. * pow(ma1, 2) - 1. * tmin) +
        (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
         (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
            (1. * pow(mpion, 2) - 1. * tmin) -
        (0.25 * pow(-2. + delta, 2) * pow(mpion, 2) * tmin) / pow(mrho, 2) +
        0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
            (eta2 * (-1. * pow(ma1, 2) + pow(mrho, 2) - 2. * s) +
             eta1 *
                 (pow(ma1, 2) - 1. * pow(mpion, 2) - 2. * pow(mrho, 2) + s)) *
            tmin +
        0.03125 * pow(eta1 - 1. * eta2, 2) *
            (pow(eta1, 2) *
                 (3. * pow(ma1, 4) + 4. * pow(mpion, 4) + pow(mrho, 4) +
                  pow(mpion, 2) * (8. * pow(mrho, 2) - 4. * s) -
                  4. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                  pow(ma1, 2) *
                      (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) +
             pow(eta2, 2) *
                 (3. * pow(ma1, 4) + 4. * pow(mpion, 4) + pow(mrho, 4) +
                  pow(mpion, 2) * (-4. * pow(mrho, 2) - 4. * s) -
                  2. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                  pow(ma1, 2) *
                      (-8. * pow(mpion, 2) + 4. * pow(mrho, 2) + 4. * s)) +
             eta1 * eta2 *
                 (-6. * pow(ma1, 4) - 8. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                  pow(ma1, 2) * (16. * pow(mpion, 2) - 8. * s) +
                  4. * pow(mrho, 2) * s - 4. * pow(s, 2) +
                  pow(mpion, 2) * (-4. * pow(mrho, 2) + 8. * s))) *
            tmin +
        (2. *
         (0. - 0.25 * pow(mrho, 4) - 1. * C4 * pow(mrho, 6) +
          pow(mpion, 2) * (0.75 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
          2. * C4 * pow(mrho, 4) * s +
          pow(delta, 2) *
              (-0.125 * pow(mpion, 4) - 0.1875 * pow(mrho, 4) +
               pow(mpion, 2) * (0.0625 * pow(mrho, 2) + 0.0625 * s) +
               0.1875 * pow(mrho, 2) * s) +
          delta *
              (0.25 * pow(mpion, 4) + 0.5 * C4 * pow(mrho, 6) +
               pow(mpion, 2) *
                   (-0.5 * pow(mrho, 2) + 0.5 * C4 * pow(mrho, 4) - 0.125 * s) -
               0.375 * pow(mrho, 2) * s + pow(mrho, 4) * (0.5 - 1. * C4 * s))) *
         tmin) /
            pow(mrho, 4) -
        (0.0625 *
         (0. + 8. * pow(mpion, 2) * pow(mrho, 4) - 12. * pow(mrho, 6) +
          4. * pow(mrho, 4) * s +
          delta * pow(mrho, 2) *
              (-16. * pow(mpion, 4) - 16. * pow(mpion, 2) * pow(mrho, 2) -
               4. * pow(mrho, 4) + 16. * pow(mrho, 2) * s + 4. * pow(s, 2)) +
          pow(delta, 2) *
              (8. * pow(mpion, 6) + 9. * pow(mrho, 6) +
               pow(mpion, 4) * (4. * pow(mrho, 2) - 4. * s) -
               13. * pow(mrho, 4) * s - 5. * pow(mrho, 2) * pow(s, 2) +
               1. * pow(s, 3) +
               pow(mpion, 2) * (-2. * pow(mrho, 4) + 4. * pow(mrho, 2) * s -
                                2. * pow(s, 2)))) *
         tmin) /
            pow(mrho, 6) -
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (pow(mrho, 2) * (eta1 * (2. * pow(ma1, 2) + 2. * pow(mrho, 2)) +
                          eta2 * (-2. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                                  8. * pow(mrho, 2) + 6. * s)) +
          delta *
              (eta1 *
                   (1. * pow(ma1, 4) - 2. * pow(mpion, 4) - 3. * pow(mrho, 4) +
                    pow(mpion, 2) * (4. * pow(mrho, 2) - 4. * s) -
                    2. * pow(mrho, 2) * s + 5.000000000000001 * pow(s, 2) +
                    pow(ma1, 2) * (-2. * pow(mpion, 2) + 1. * s)) +
               eta2 *
                   (-1. * pow(ma1, 4) - 4. * pow(mpion, 4) + 4. * pow(mrho, 4) +
                    pow(mpion, 2) * (-1. * pow(mrho, 2) - 2. * s) +
                    1. * pow(mrho, 2) * s - 1. * pow(s, 2) +
                    pow(ma1, 2) *
                        (3. * pow(mpion, 2) - 3. * pow(mrho, 2) + 2. * s)))) *
         tmin) /
            pow(mrho, 2) -
        (0.5 *
         (pow(mrho, 6) *
              (-1.5 + C4 * (-12. * pow(mpion, 2) + 6. * s) +
               pow(C4, 2) * (-16. * pow(mpion, 4) + 16. * pow(mpion, 2) * s -
                             4. * pow(s, 2))) +
          pow(delta, 2) *
              (1. * pow(mpion, 6) - 2. * pow(mpion, 4) * pow(mrho, 2) +
               0.125 * pow(mrho, 6) + 0.25 * pow(mrho, 4) * s -
               0.875 * pow(mrho, 2) * pow(s, 2) + 0.25 * pow(s, 3) +
               pow(mpion, 2) * (-2.5 * pow(mrho, 4) + 2.25 * pow(mrho, 2) * s -
                                0.75 * pow(s, 2))) +
          delta * pow(mrho, 2) *
              (pow(mpion, 4) * (1. + 8. * C4 * pow(mrho, 2)) + 0.5 * pow(s, 2) +
               pow(mrho, 4) * (1.5 - 5. * C4 * s) +
               pow(mrho, 2) * s * (-0.5 + 1. * C4 * s) +
               pow(mpion, 2) * (6. * C4 * pow(mrho, 4) - 1.5 * s +
                                pow(mrho, 2) * (3. - 6. * C4 * s)))) *
         tmin) /
            pow(mrho, 6) -
        (0.5 *
         (0. - 4. * C4 * pow(mrho, 8) - 0.5 * pow(mrho, 4) * s +
          pow(mrho, 6) * (2. + 2. * C4 * s) +
          pow(delta, 2) *
              (-2. * pow(mpion, 6) - 2. * pow(mrho, 6) +
               0.5 * pow(mpion, 4) * s + 2.125 * pow(mrho, 4) * s +
               1.25 * pow(mrho, 2) * pow(s, 2) - 0.375 * pow(s, 3) +
               pow(mpion, 2) * (1.5 * pow(mrho, 4) - 1.5 * pow(mrho, 2) * s +
                                1. * pow(s, 2))) +
          delta * pow(mrho, 2) *
              (2. * pow(mpion, 4) + 2. * C4 * pow(mrho, 6) - 1. * pow(s, 2) +
               pow(mrho, 2) * s * (-3. + 1. * C4 * s) +
               pow(mrho, 4) * (1. + 1. * C4 * s) +
               pow(mpion, 2) * (1. * s + pow(mrho, 2) * (1. - 2. * C4 * s)))) *
         tmin) /
            pow(mrho, 6) +
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (delta *
              (eta1 * (3. * pow(ma1, 4) + 6. * pow(mpion, 4) + pow(mrho, 4) +
                       pow(mpion, 2) * (18. * pow(mrho, 2) - 12. * s) -
                       8. * pow(mrho, 2) * s + 7. * pow(s, 2) +
                       pow(ma1, 2) * (-10. * pow(mpion, 2) - 4. * pow(mrho, 2) +
                                      5. * s)) +
               eta2 * (-3. * pow(ma1, 4) - 12. * pow(mpion, 4) +
                       2. * pow(mrho, 4) +
                       pow(ma1, 2) *
                           (11. * pow(mpion, 2) - 3. * pow(mrho, 2) - 2. * s) +
                       5. * pow(mrho, 2) * s - 3. * pow(s, 2) +
                       pow(mpion, 2) * (-1. * pow(mrho, 2) + 6. * s))) +
          pow(mrho, 2) *
              (eta1 *
                   (-8. * C4 * pow(ma1, 4) - 32. * C4 * pow(mpion, 4) -
                    6. * pow(mrho, 2) +
                    pow(ma1, 2) * (6. + C4 * (32. * pow(mpion, 2) +
                                              8. * pow(mrho, 2) - 16. * s)) +
                    4. * s + 16. * C4 * pow(mrho, 2) * s - 8. * C4 * pow(s, 2) +
                    pow(mpion, 2) *
                        (-12. - 32. * C4 * pow(mrho, 2) + 32. * C4 * s)) +
               eta2 * (8. * C4 * pow(ma1, 4) + 32. * C4 * pow(mpion, 4) -
                       4. * pow(mrho, 2) - 2. * s + 8. * C4 * pow(s, 2) +
                       pow(mpion, 2) *
                           (10. - 16. * C4 * pow(mrho, 2) - 32. * C4 * s) +
                       pow(ma1, 2) *
                           (-6. + C4 * (-32. * pow(mpion, 2) +
                                        8. * pow(mrho, 2) + 16. * s))))) *
         tmin) /
            pow(mrho, 2) +
        0.0625 * (-2. + delta) * pow(eta1 - 1. * eta2, 2) * pow(tmin, 2) +
        (0.1875 *
         (1.3333333333333333 * pow(mrho, 2) +
          5.333333333333333 * C4 * pow(mrho, 4) +
          pow(delta, 2) *
              (1. * pow(mpion, 2) + 1.3333333333333333 * pow(mrho, 2) -
               0.3333333333333333 * s) +
          delta * (-2. * pow(mpion, 2) - 3.3333333333333335 * pow(mrho, 2) -
                   2.6666666666666665 * C4 * pow(mrho, 4) +
                   0.6666666666666666 * s)) *
         pow(tmin, 2)) /
            pow(mrho, 4) +
        0.03125 * pow(eta1 - 1. * eta2, 3) *
            (eta2 * (-1. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                     1. * pow(mrho, 2) - 1. * s) +
             eta1 *
                 (pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
            pow(tmin, 2) -
        (0.375 *
         (0.3333333333333333 * pow(mrho, 4) -
          1.3333333333333333 * C4 * pow(mrho, 6) +
          delta * pow(mrho, 2) *
              (1.3333333333333333 * pow(mrho, 2) -
               0.6666666666666666 * C4 * pow(mrho, 4) +
               pow(mpion, 2) * (-0.6666666666666666 +
                                1.3333333333333333 * C4 * pow(mrho, 2)) -
               0.6666666666666666 * s) +
          pow(delta, 2) * (1. * pow(mpion, 4) + 0.25 * pow(mrho, 4) +
                           pow(mpion, 2) * (-0.3333333333333333 * pow(mrho, 2) +
                                            0.6666666666666666 * s) -
                           0.5833333333333334 * pow(s, 2))) *
         pow(tmin, 2)) /
            pow(mrho, 6) -
        (0.03125 * (1. * eta1 - 1. * eta2) *
         ((2. * eta1 - 2. * eta2) * pow(mrho, 2) +
          delta * (eta1 * (1. * pow(ma1, 2) - 2. * pow(mpion, 2) + 1. * s) +
                   eta2 * (-1. * pow(ma1, 2) + 3. * pow(mpion, 2) -
                           3. * pow(mrho, 2) + 2. * s))) *
         pow(tmin, 2)) /
            pow(mrho, 2) +
        (0.03125 *
         (0. - 4. * pow(mrho, 4) +
          delta * (16. * pow(mrho, 4) - 8. * pow(mrho, 2) * s) +
          pow(delta, 2) * (4. * pow(mpion, 4) - 3. * pow(mrho, 4) +
                           2. * pow(mrho, 2) * s - 3. * pow(s, 2) +
                           pow(mpion, 2) * (-4. * pow(mrho, 2) + 4. * s))) *
         pow(tmin, 2)) /
            pow(mrho, 6) +
        (0.25 *
         (C4 * pow(mrho, 6) * (-6. - 16. * C4 * pow(mpion, 2) + 8. * C4 * s) +
          pow(delta, 2) * (1. * pow(mpion, 4) - 0.25 * pow(mrho, 4) +
                           pow(mpion, 2) * (-1.75 * pow(mrho, 2) + 0.5 * s) +
                           0.5 * pow(mrho, 2) * s - 0.5 * pow(s, 2)) +
          delta * pow(mrho, 2) *
              (1. * C4 * pow(mrho, 4) +
               pow(mpion, 2) * (0.5 + 10. * C4 * pow(mrho, 2)) - 0.5 * s +
               pow(mrho, 2) * (2.5 - 4. * C4 * s))) *
         pow(tmin, 2)) /
            pow(mrho, 6) +
        (0.09375 * (1. * eta1 - 1. * eta2) *
         (delta *
              (eta2 * (-1. * pow(ma1, 2) + 3.6666666666666665 * pow(mpion, 2) -
                       1. * pow(mrho, 2) - 0.6666666666666666 * s) +
               eta1 * (1. * pow(ma1, 2) - 3.3333333333333335 * pow(mpion, 2) -
                       1.3333333333333333 * pow(mrho, 2) +
                       1.6666666666666667 * s)) +
          pow(mrho, 2) *
              (eta1 * (2. + C4 * (-2.6666666666666665 * pow(ma1, 2) +
                                  10.666666666666666 * pow(mpion, 2) +
                                  2.6666666666666665 * pow(mrho, 2) -
                                  5.333333333333333 * s)) +
               eta2 * (-2. + C4 * (2.6666666666666665 * pow(ma1, 2) -
                                   10.666666666666666 * pow(mpion, 2) +
                                   2.6666666666666665 * pow(mrho, 2) +
                                   5.333333333333333 * s)))) *
         pow(tmin, 2)) /
            pow(mrho, 2) +
        0.010416666666666666 * pow(eta1 - 1. * eta2, 4) * pow(tmin, 3) -
        (0.041666666666666664 * delta * (-2. + 1. * delta) * pow(tmin, 3)) /
            pow(mrho, 4) -
        (0.020833333333333332 * delta * pow(1. * eta1 - 1. * eta2, 2) *
         pow(tmin, 3)) /
            pow(mrho, 2) -
        (0.16666666666666666 * pow(1. * eta1 - 1. * eta2, 2) *
         (-0.375 * delta + 1. * C4 * pow(mrho, 2)) * pow(tmin, 3)) /
            pow(mrho, 2) +
        (0.10416666666666666 * delta *
         (-0.8 * pow(mrho, 2) +
          delta * (0.4 * pow(mpion, 2) + 1. * pow(mrho, 2) - 0.6 * s)) *
         pow(tmin, 3)) /
            pow(mrho, 6) +
        (0.16666666666666666 * delta *
         (0. - 0.75 * delta * pow(mrho, 2) + 1. * C4 * pow(mrho, 4) +
          0.625 * delta * s) *
         pow(tmin, 3)) /
            pow(mrho, 6) -
        (0.041666666666666664 *
         (12. * C4 * delta * pow(mrho, 4) - 16. * pow(C4, 2) * pow(mrho, 6) +
          pow(delta, 2) * (1. * pow(mpion, 2) - 2.5 * pow(mrho, 2) + 1. * s)) *
         pow(tmin, 3)) /
            pow(mrho, 6) +
        (0. -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(mpion, 4) *
             pow(mrho, 6) +
         0.25 * pow(2. - 1. * delta, 2) * pow(mrho, 10) -
         0.5000000000000001 * pow(2. - 1. * delta, 2) * pow(mrho, 6) *
             pow(s, 2) +
         pow(2. - 1. * delta, 2) * pow(mpion, 2) * pow(mrho, 6) *
             (-1. * pow(mrho, 2) + 1. * s)) /
            (pow(mrho, 6) * (-2. * pow(mpion, 2) + 1. * s + 1. * tmin)) -
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (2. * pow(mrho, 2) + delta * (1. * pow(ma1, 2) - 2. * pow(mpion, 2) -
                                       1. * pow(mrho, 2) + 1. * s)) *
         (eta2 *
              (-1. * pow(ma1, 6) +
               pow(mpion, 2) * (4. * pow(mpion, 2) - 1. * s) * s +
               pow(ma1, 4) * (3. * pow(mpion, 2) - 4. * pow(mrho, 2) + 2. * s) +
               pow(ma1, 2) * (-4. * pow(mpion, 4) - 2. * pow(mpion, 2) * s +
                              (4. * pow(mrho, 2) - 1. * s) * s)) +
          eta1 * (1. * pow(ma1, 6) + 8. * pow(mpion, 6) +
                  pow(mpion, 4) * (-4. * pow(mrho, 2) - 6. * s) +
                  pow(mpion, 2) * (4. * pow(mrho, 2) - 2. * s) * s +
                  pow(ma1, 4) *
                      (-2. * pow(mpion, 2) + 1. * pow(mrho, 2) + 1. * s) +
                  s * (2. * pow(mrho, 4) - 3. * pow(mrho, 2) * s +
                       1. * pow(s, 2)) +
                  pow(ma1, 2) * (-2. * pow(mpion, 4) - 2. * pow(mrho, 4) +
                                 pow(mpion, 2) * (4. * pow(mrho, 2) - 4. * s) -
                                 2. * pow(mrho, 2) * s + 5. * pow(s, 2)))) *
         log(abs(-pow(ma1, 2) + tmin))) /
            (pow(mrho, 2) * (pow(ma1, 2) - 2. * pow(mpion, 2) + s)) +
        (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
         (eta2 * (-1. * pow(ma1, 6) + pow(mpion, 6) -
                  1. * pow(mpion, 4) * pow(mrho, 2) +
                  pow(ma1, 4) * (pow(mpion, 2) + pow(mrho, 2) - 2. * s) +
                  pow(ma1, 2) * (3. * pow(mpion, 4) - 2. * pow(mpion, 2) * s)) +
          eta1 * (pow(ma1, 6) +
                  pow(ma1, 4) * (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + s) +
                  pow(mpion, 2) * (-4. * pow(mpion, 4) - 1. * pow(mrho, 4) -
                                   1. * pow(mrho, 2) * s +
                                   pow(mpion, 2) * (3. * pow(mrho, 2) + s)) +
                  pow(ma1, 2) *
                      (pow(mpion, 4) + pow(mrho, 4) - 1. * pow(mrho, 2) * s +
                       pow(mpion, 2) * (pow(mrho, 2) + 2. * s)))) *
         log(abs(-pow(ma1, 2) + tmin))) /
            (pow(ma1, 2) - 1. * pow(mpion, 2)) +
        0.0625 * pow(eta1 - 1. * eta2, 2) *
            (pow(eta1, 2) *
                 (2. * pow(ma1, 6) +
                  pow(mpion, 4) * (-3. * pow(mrho, 2) - 1. * s) +
                  pow(mrho, 2) * (pow(mrho, 2) - 1. * s) * s +
                  pow(ma1, 4) *
                      (-6. * pow(mpion, 2) - 3. * pow(mrho, 2) + 3. * s) +
                  pow(mpion, 2) * (-2. * pow(mrho, 4) + 3. * pow(mrho, 2) * s) +
                  pow(ma1, 2) * (4. * pow(mpion, 4) + pow(mrho, 4) +
                                 pow(mpion, 2) * (8. * pow(mrho, 2) - 4. * s) -
                                 4. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
             pow(eta2, 2) *
                 (2. * pow(ma1, 6) +
                  pow(ma1, 4) *
                      (-6. * pow(mpion, 2) + 3. * pow(mrho, 2) + 3. * s) +
                  pow(mpion, 2) *
                      (-1. * pow(mrho, 4) +
                       pow(mpion, 2) * (2. * pow(mrho, 2) - 1. * s) +
                       pow(mrho, 2) * s) +
                  pow(ma1, 2) * (4. * pow(mpion, 4) + pow(mrho, 4) +
                                 pow(mpion, 2) * (-4. * pow(mrho, 2) - 4. * s) -
                                 2. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
             eta1 * eta2 *
                 (-4. * pow(ma1, 6) +
                  pow(ma1, 4) * (12. * pow(mpion, 2) - 6. * s) +
                  pow(mpion, 2) *
                      (-2. * pow(mrho, 4) - 2. * pow(mrho, 2) * s +
                       pow(mpion, 2) * (2. * pow(mrho, 2) + 2. * s)) +
                  pow(ma1, 2) *
                      (-8. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                       4. * pow(mrho, 2) * s - 4. * pow(s, 2) +
                       pow(mpion, 2) * (-4. * pow(mrho, 2) + 8. * s)))) *
            log(abs(-pow(ma1, 2) + tmin)) +
        (0.0625 * (1. * eta1 - 1. * eta2) *
         (delta *
              (eta1 *
                   (3. * pow(ma1, 6) + 8. * pow(mpion, 6) +
                    pow(mpion, 4) * (-12. * pow(mrho, 2) - 6. * s) +
                    pow(ma1, 4) *
                        (-10. * pow(mpion, 2) - 4. * pow(mrho, 2) + 5. * s) +
                    pow(mpion, 2) * (-4. * pow(mrho, 4) +
                                     10. * pow(mrho, 2) * s - 2. * pow(s, 2)) +
                    s * (3. * pow(mrho, 4) - 4. * pow(mrho, 2) * s +
                         pow(s, 2)) +
                    pow(ma1, 2) *
                        (6. * pow(mpion, 4) + pow(mrho, 4) +
                         pow(mpion, 2) * (18. * pow(mrho, 2) - 12. * s) -
                         8. * pow(mrho, 2) * s + 7. * pow(s, 2))) +
               eta2 * (-3. * pow(ma1, 6) +
                       pow(ma1, 4) *
                           (11. * pow(mpion, 2) - 3. * pow(mrho, 2) - 2. * s) +
                       pow(mpion, 2) *
                           (-2. * pow(mrho, 4) + pow(mrho, 2) * s -
                            1. * pow(s, 2) +
                            pow(mpion, 2) * (-2. * pow(mrho, 2) + 4. * s)) +
                       pow(ma1, 2) *
                           (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                            5. * pow(mrho, 2) * s - 3. * pow(s, 2) +
                            pow(mpion, 2) * (-1. * pow(mrho, 2) + 6. * s)))) +
          pow(mrho, 2) *
              (eta2 *
                   (8. * C4 * pow(ma1, 6) +
                    pow(mpion, 2) *
                        ((4. + 8. * C4 * pow(mpion, 2)) * pow(mrho, 2) -
                         2. * s) +
                    pow(ma1, 4) * (-6. + C4 * (-32. * pow(mpion, 2) +
                                               8. * pow(mrho, 2) + 16. * s)) +
                    pow(ma1, 2) *
                        (32. * C4 * pow(mpion, 4) - 4. * pow(mrho, 2) +
                         pow(mpion, 2) *
                             (10. - 16. * C4 * pow(mrho, 2) - 32. * C4 * s) +
                         s * (-2. + 8. * C4 * s))) +
               eta1 * (-8. * C4 * pow(ma1, 6) +
                       pow(mpion, 4) * (4. + 24. * C4 * pow(mrho, 2)) +
                       pow(ma1, 4) * (6. + C4 * (32. * pow(mpion, 2) +
                                                 8. * pow(mrho, 2) - 16. * s)) +
                       s * (-2. * pow(mrho, 2) + 2. * s) +
                       pow(mpion, 2) *
                           (-4. * s + pow(mrho, 2) * (8. - 16. * C4 * s)) +
                       pow(ma1, 2) *
                           (-32. * C4 * pow(mpion, 4) + s * (4. - 8. * C4 * s) +
                            pow(mrho, 2) * (-6. + 16. * C4 * s) +
                            pow(mpion, 2) * (-12. - 32. * C4 * pow(mrho, 2) +
                                             32. * C4 * s))))) *
         log(abs(-pow(ma1, 2) + tmin))) /
            pow(mrho, 2) +
        0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
            log(abs(-pow(mpion, 2) + tmin)) -
        (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
         (eta2 * (2. * pow(mpion, 6) - 2. * pow(mpion, 4) * s) +
          eta1 * (-2. * pow(mpion, 6) - 1. * pow(mpion, 2) * pow(mrho, 2) * s +
                  pow(mpion, 4) * (pow(mrho, 2) + 2. * s))) *
         log(abs(-pow(mpion, 2) + tmin))) /
            (pow(ma1, 2) - 1. * pow(mpion, 2)) -
        (0.125 *
         (0. - 32. * C4 * pow(mpion, 6) * pow(mrho, 4) - 8. * pow(mrho, 8) +
          8. * pow(mrho, 6) * s +
          pow(mpion, 4) * pow(mrho, 4) * (16. + 64. * C4 * s) +
          pow(mpion, 2) * pow(mrho, 4) *
              (24. * pow(mrho, 2) + s * (-16. - 32. * C4 * s)) +
          pow(delta, 2) * pow(mrho, 2) *
              (-4. * pow(mpion, 6) - 2. * pow(mrho, 6) + 2. * pow(mrho, 4) * s +
               pow(mpion, 4) * (4. * pow(mrho, 2) + 8. * s) +
               pow(mpion, 2) * (6. * pow(mrho, 4) - 4. * pow(mrho, 2) * s -
                                4.000000000000001 * pow(s, 2))) +
          delta * pow(mrho, 2) *
              (8. * pow(mrho, 6) +
               pow(mpion, 6) * (8. + 16. * C4 * pow(mrho, 2)) -
               8. * pow(mrho, 4) * s +
               pow(mpion, 4) * (-16. * s + pow(mrho, 2) * (-15.999999999999996 -
                                                           32. * C4 * s)) +
               pow(mpion, 2) * (-24. * pow(mrho, 4) + 8. * pow(s, 2) +
                                pow(mrho, 2) * s * (16. + 16. * C4 * s)))) *
         log(abs(-pow(mpion, 2) + tmin))) /
            (pow(mrho, 4) * (-1. * pow(mpion, 2) + 1. * s)) -
        (0.25 * (1. * eta1 - 1. * eta2) *
         (eta2 * ((2. - 1. * delta) * pow(mpion, 6) +
                  pow(mpion, 2) * s *
                      ((-12. + 6. * delta) * pow(mrho, 2) +
                       (6. - 3. * delta) * s) +
                  pow(s, 2) * ((4. - 2. * delta) * pow(mrho, 2) +
                               (-2. + 1. * delta) * s) +
                  pow(mpion, 4) * ((8. - 4. * delta) * pow(mrho, 2) +
                                   (-6. + 3. * delta) * s)) +
          eta1 * ((-2. + 1. * delta) * pow(mpion, 6) +
                  (-2. + 1. * delta) * pow(mrho, 4) * s +
                  (2. - 1. * delta) * pow(s, 3) +
                  pow(mpion, 4) * ((-4. + 2. * delta) * pow(mrho, 2) +
                                   (6. - 3. * delta) * s) +
                  pow(mpion, 2) * ((2. - 1. * delta) * pow(mrho, 4) +
                                   (4. - 2. * delta) * pow(mrho, 2) * s +
                                   (-6. + 3. * delta) * pow(s, 2)))) *
         log(abs(-2. * pow(mpion, 2) + s + tmin))) /
            (pow(ma1, 2) - 2. * pow(mpion, 2) + s) +
        (0.125 *
         (0. +
          (32. - 31.999999999999993 * delta + 8. * pow(delta, 2)) *
              pow(mpion, 4) * pow(mrho, 4) -
          2.0000000000000004 * pow(2. - 1. * delta, 2) * pow(mrho, 8) +
          pow(mpion, 2) * pow(mrho, 4) *
              (8.000000000000002 * pow(2. - 1. * delta, 2) * pow(mrho, 2) +
               (-32. + 31.999999999999996 * delta - 8. * pow(delta, 2)) * s)) *
         log(abs(-2. * pow(mpion, 2) + s + tmin))) /
            (pow(mrho, 4) * (-1. * pow(mpion, 2) + 1. * s)) +
        (0.25 *
         (0. + 8. * pow(mpion, 2) * pow(mrho, 6) + 4. * pow(mrho, 8) -
          8. * pow(mrho, 6) * s +
          delta * pow(mrho, 4) *
              (8. * pow(mpion, 4) - 8. * pow(mrho, 4) +
               pow(mpion, 2) * (8. * pow(mrho, 2) - 16. * s) +
               8. * pow(mrho, 2) * s + 8. * pow(s, 2)) +
          pow(delta, 2) * pow(mrho, 4) *
              (-4. * pow(mpion, 4) + 3. * pow(mrho, 4) - 2. * pow(mrho, 2) * s -
               4. * pow(s, 2) +
               pow(mpion, 2) * (-6. * pow(mrho, 2) + 8. * s))) *
         log(abs(-2. * pow(mpion, 2) + s + tmin))) /
            pow(mrho, 6) -
        (0.5 *
         (0. + pow(mpion, 2) * (4. * pow(mrho, 6) - 8. * C4 * pow(mrho, 8)) -
          4. * pow(mrho, 6) * s + 8. * C4 * pow(mrho, 8) * s +
          pow(delta, 2) * pow(mrho, 4) *
              (-2. * pow(mpion, 4) + 1. * pow(mrho, 4) - 2. * pow(s, 2) +
               pow(mpion, 2) * (-4. * pow(mrho, 2) + 4. * s)) +
          delta * pow(mrho, 4) *
              (4. * pow(mpion, 4) +
               pow(mpion, 2) *
                   (6. * pow(mrho, 2) + 4. * C4 * pow(mrho, 4) - 8. * s) +
               2. * pow(mrho, 2) * s + 4. * pow(s, 2) +
               pow(mrho, 4) * (-2. - 4. * C4 * s))) *
         log(abs(-2. * pow(mpion, 2) + s + tmin))) /
            pow(mrho, 6))) /
          (16. * Pi *
           (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
            2 * pow(mpion, 2) * (pow(mrho, 2) + s)));

  return xs * to_mb / spin_deg_factor;
}

double PhotonCrossSection<ComputationMethod::Analytic>::
    xs_diff_pi_rho_pi0_rho_mediated(const double s, const double t,
                                    const double m_rho) {
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  const double spin_deg_factor = 3.0;

  const double diff_xs =
      ((pow(Const, 2) * pow(ghat, 4) *
        ((-0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
          (pow(mpion, 4) + pow(pow(mrho, 2) - t, 2) -
           2 * pow(mpion, 2) * (pow(mrho, 2) + t))) /
             (pow(mrho, 2) * pow(pow(mpion, 2) - t, 2)) -
         (0.0625 * (eta1 - eta2) *
          (2 * pow(mrho, 2) +
           delta * (-2 * pow(mpion, 2) - pow(mrho, 2) + s + t)) *
          (eta1 * (8 * pow(mpion, 6) + pow(s, 3) + 2 * pow(mrho, 4) * (s - t) +
                   5 * pow(s, 2) * t + s * pow(t, 2) + pow(t, 3) +
                   2 * pow(mpion, 2) * (2 * pow(mrho, 2) - s - t) * (s + t) -
                   pow(mrho, 2) * (3 * s - t) * (s + t) -
                   2 * pow(mpion, 4) * (2 * pow(mrho, 2) + 3 * s + t)) +
           eta2 * (s - t) *
               (4 * pow(mpion, 4) + t * (4 * pow(mrho, 2) - s + t) -
                pow(mpion, 2) * (s + 3 * t)))) /
             (pow(mrho, 2) * (-pow(ma1, 2) + t) *
              (-2 * pow(mpion, 2) + s + t)) -
         (0.0625 *
          pow(-2. * pow(mrho, 2) +
                  delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t),
              2) *
          (8. * pow(mpion, 6) + 4. * pow(mrho, 6) + pow(s, 3) +
           pow(mrho, 4) * (-4. * s - 4. * t) +
           pow(mpion, 4) * (-4. * pow(mrho, 2) - 4. * s - 4. * t) +
           3. * pow(s, 2) * t + 3. * s * pow(t, 2) + pow(t, 3) +
           pow(mrho, 2) * (-3. * pow(s, 2) + 2. * s * t - 3. * pow(t, 2)) +
           pow(mpion, 2) *
               (-8. * pow(mrho, 4) - 2. * pow(s, 2) - 4. * s * t -
                2. * pow(t, 2) + pow(mrho, 2) * (4. * s + 4. * t)))) /
             (pow(mrho, 6) * pow(2. * pow(mpion, 2) - 1. * s - 1. * t, 2)) +
         (0.125 * (-2 + delta) * (eta1 - eta2) *
          (eta2 * (pow(mpion, 2) + t) *
               (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * t) +
                (pow(mrho, 2) - 2 * s - t) * t) +
           eta1 * (-4 * pow(mpion, 6) +
                   (pow(mrho, 2) - t) * (pow(mrho, 2) - s - t) * t +
                   pow(mpion, 4) * (3 * pow(mrho, 2) + s + t) -
                   pow(mpion, 2) * (pow(mrho, 4) + pow(mrho, 2) * (s - t) +
                                    2 * t * (-s + t))))) /
             ((-pow(ma1, 2) + t) * (-pow(mpion, 2) + t)) +
         (0.03125 * pow(eta1 - eta2, 2) *
          (-2 * eta1 * eta2 *
               (pow(mpion, 8) -
                pow(mpion, 4) * (pow(mrho, 4) + 2 * (pow(mrho, 2) + s) * t -
                                 4 * pow(t, 2)) +
                pow(t, 2) * (-pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                             2 * pow(s, 2) + 2 * s * t + pow(t, 2)) +
                2 * pow(mpion, 2) * t *
                    (pow(mrho, 4) + pow(mrho, 2) * (s + t) - 2 * t * (s + t))) +
           pow(eta2, 2) *
               (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
                pow(mpion, 4) * (pow(mrho, 4) + 4 * pow(mrho, 2) * t -
                                 2 * (s - 2 * t) * t) +
                pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t +
                             pow(t, 2) + 2 * pow(mrho, 2) * (-s + t)) -
                2 * pow(mpion, 2) * t *
                    (pow(mrho, 4) - pow(mrho, 2) * (s - 2 * t) +
                     2 * t * (s + t))) +
           pow(eta1, 2) *
               (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
                pow(mpion, 4) *
                    (3 * pow(mrho, 4) + 2 * pow(mrho, 2) * (s - 3 * t) -
                     2 * (s - 2 * t) * t) +
                t * (-pow(mrho, 2) + t) *
                    (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                     pow(mrho, 2) * (2 * s + t)) -
                2 * pow(mpion, 2) * (-pow(mrho, 2) + t) *
                    (2 * t * (s + t) - pow(mrho, 2) * (s + 2 * t))))) /
             pow(pow(ma1, 2) - t, 2) -
         (0.5 *
          (-2. * pow(mrho, 2) +
           delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t)) *
          (delta *
               (-1. * pow(mpion, 6) - 0.5 * pow(mrho, 6) - 0.1875 * pow(s, 3) +
                pow(mpion, 2) * (1. * pow(mrho, 4) +
                                 pow(mrho, 2) * (-0.625 * s - 0.375 * t) +
                                 s * (0.5 * s + 0.5 * t)) +
                pow(mrho, 4) * (0.5 * s + 0.5 * t) +
                pow(mpion, 4) * (0.5 * pow(mrho, 2) + 0.25 * s + 0.75 * t) -
                0.4375 * pow(s, 2) * t - 0.3125 * s * pow(t, 2) -
                0.0625 * pow(t, 3) +
                pow(mrho, 2) *
                    (0.4375 * pow(s, 2) - 0.25 * s * t + 0.3125 * pow(t, 2))) +
           pow(mrho, 2) *
               (-0.125 * pow(s, 2) + C4 * pow(mrho, 4) * (1. * s - 1. * t) +
                0.125 * pow(t, 2) +
                pow(mpion, 2) * ((0.25 - 1. * C4 * pow(mrho, 2)) * s +
                                 (-0.25 + 1. * C4 * pow(mrho, 2)) * t) +
                pow(mrho, 2) * (-0.5 * s + 0.5 * C4 * pow(s, 2) +
                                t * (0.5 - 0.5 * C4 * t))))) /
             (pow(mrho, 6) * (1. * pow(mpion, 2) - 0.5 * s - 0.5 * t)) +
         (pow(delta, 2) *
              (-0.5 * pow(mpion, 6) - 0.0625 * pow(mrho, 6) +
               pow(mrho, 4) * (-0.125 * s - 0.125 * t) +
               pow(mpion, 4) * (1. * pow(mrho, 2) + 0.5 * t) +
               s * (-0.125 * pow(s, 2) - 0.25 * s * t - 0.125 * pow(t, 2)) +
               pow(mpion, 2) * (1.25 * pow(mrho, 4) + 0.375 * pow(s, 2) +
                                pow(mrho, 2) * (-1.125 * s - 0.875 * t) +
                                0.25 * s * t - 0.125 * pow(t, 2)) +
               pow(mrho, 2) *
                   (0.4375 * pow(s, 2) + 0.25 * s * t + 0.3125 * pow(t, 2))) +
          pow(mrho, 6) *
              (0.75 +
               C4 * (8. * C4 * pow(mpion, 4) + 2. * C4 * pow(s, 2) +
                     pow(mpion, 2) * (6. - 8. * C4 * s - 8. * C4 * t) +
                     t * (-3. + 2. * C4 * t) + s * (-3. + 4. * C4 * t))) +
          delta * pow(mrho, 2) *
              (pow(mpion, 4) * (-0.5 - 4. * C4 * pow(mrho, 2)) +
               s * (-0.25 * s - 0.25 * t) +
               pow(mrho, 4) * (-0.75 + 2.5 * C4 * s + 0.5 * C4 * t) +
               pow(mrho, 2) *
                   (-0.5 * C4 * pow(s, 2) + s * (0.25 - 2. * C4 * t) +
                    t * (1.25 - 1.5 * C4 * t)) +
               pow(mpion, 2) *
                   (-3. * C4 * pow(mrho, 4) + 0.75 * s + 0.25 * t +
                    pow(mrho, 2) * (-1.5 + 3. * C4 * s + 5. * C4 * t)))) /
             pow(mrho, 6) +
         (2 *
          ((0.0625 * (-2. + delta) *
            (-2. * pow(mrho, 2) +
             delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t)) *
            (2. * pow(mpion, 6) + 1. * pow(mrho, 6) +
             pow(mpion, 4) * (-3. * pow(mrho, 2) - 2. * t) +
             pow(mrho, 4) * (-1.5 * s - 1.5 * t) +
             pow(mrho, 2) * s * (0.5 * s + 0.5 * t) +
             pow(mpion, 2) * (-1. * pow(mrho, 4) - 0.5 * pow(s, 2) +
                              pow(mrho, 2) * (2.5 * s - 0.5 * t) - 1. * s * t -
                              0.5 * pow(t, 2)) +
             t * (0.5 * pow(s, 2) + 1. * s * t + 0.5 * pow(t, 2)))) /
               ((pow(mpion, 2) - 1. * t) *
                (1. * pow(mpion, 2) - 0.5 * s - 0.5 * t)) +
           (0.0625 * (-2 + delta) *
            (6 * delta * pow(mpion, 6) + delta * s * t * (s + t) +
             pow(mrho, 6) * (-2 + 3 * delta + 8 * C4 * t) -
             pow(mpion, 4) * ((-2 + 9 * delta) * pow(mrho, 2) -
                              8 * C4 * pow(mrho, 4) + delta * (s + 9 * t)) -
             2 * pow(mrho, 4) *
                 (t * (-1 + 3 * delta + 4 * C4 * t) +
                  s * (-1 + 2 * delta + 8 * C4 * t)) -
             pow(mpion, 2) *
                 (8 * C4 * pow(mrho, 6) +
                  2 * pow(mrho, 4) * (-2 + delta - 8 * C4 * t) +
                  pow(mrho, 2) * ((2 - 7 * delta) * s + (2 + 5 * delta) * t) +
                  delta * (pow(s, 2) - 3 * pow(t, 2))) +
             pow(mrho, 2) * (2 * s * t + delta * (pow(s, 2) + 3 * s * t +
                                                  3 * pow(t, 2))))) /
               (-pow(mpion, 2) + t))) /
             pow(mrho, 4) +
         (0.0625 * (eta1 - eta2) *
          (-(eta2 *
             (-2 * pow(mpion, 4) *
                  (4 * C4 * pow(mrho, 2) * (pow(mrho, 2) + 4 * t) -
                   delta * (pow(mrho, 2) - 2 * s + 6 * t)) +
              pow(mpion, 2) *
                  (2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * t) +
                   delta * (pow(s, 2) - 6 * s * t - 11 * pow(t, 2)) +
                   pow(mrho, 2) * (-((-2 + delta) * s) + (-10 + delta) * t +
                                   32 * C4 * t * (s + t))) +
              t * (-2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * t) +
                   delta * (3 * pow(s, 2) + 2 * s * t + 3 * pow(t, 2)) +
                   pow(mrho, 2) * ((2 - 5 * delta) * s + 3 * (2 + delta) * t -
                                   8 * C4 * pow(s + t, 2))))) +
           eta1 * (8 * delta * pow(mpion, 6) +
                   delta * (pow(s, 3) + 7 * pow(s, 2) * t + 5 * s * pow(t, 2) +
                            3 * pow(t, 3)) -
                   2 * pow(mrho, 2) *
                       ((-1 + 2 * delta) * pow(s, 2) +
                        2 * (-1 + 2 * delta) * s * t +
                        (-3 + 2 * delta) * pow(t, 2) +
                        4 * C4 * t * pow(s + t, 2)) +
                   pow(mpion, 4) *
                       (24 * C4 * pow(mrho, 4) + 6 * delta * (-s + t) -
                        4 * pow(mrho, 2) * (-1 + 3 * delta + 8 * C4 * t)) +
                   pow(mrho, 4) * (t * (-6 + delta + 8 * C4 * t) +
                                   s * (-2 + 3 * delta + 16 * C4 * t)) -
                   2 * pow(mpion, 2) *
                       (delta * (s + t) * (s + 5 * t) -
                        pow(mrho, 2) * (-2 * s + 5 * delta * s - 6 * t +
                                        9 * delta * t + 16 * C4 * t * (s + t)) +
                        2 * pow(mrho, 4) *
                            (-2 + delta + 4 * C4 * (s + 2 * t)))))) /
             (pow(mrho, 2) * (-pow(ma1, 2) + t)))) /
       (16. * Pi *
        (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
         2 * pow(mpion, 2) * (pow(mrho, 2) + s))));

  return to_mb * diff_xs / spin_deg_factor;
}

/*----------------------------------------------------------------------------*/
/* 					Pi + Rho -> Pi + Photon channels mediated
 * by (omega) 						  */
/*----------------------------------------------------------------------------*/
// C14
double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho0_pi0(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_rho, m_pion_, 0.0);
  const double tmin = t_mandelstam[1];
  const double tmax = t_mandelstam[0];
  const double spin_deg_factor = 3.0;

  const double xs =
      (pow(Const, 2) * pow(g_POR, 4) *
       ((pow(pow(momega, 2) - s, 2) *
         (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
          pow(mpion, 4) *
              (pow(mrho, 4) + 4 * pow(momega, 4) - 2 * pow(momega, 2) * s) +
          pow(momega, 4) *
              (pow(mrho, 4) + pow(momega, 4) + 2 * pow(momega, 2) * s +
               2 * pow(s, 2) - 2 * pow(mrho, 2) * (pow(momega, 2) + s)) -
          2 * pow(mpion, 2) * pow(momega, 2) *
              (pow(mrho, 4) + 2 * pow(momega, 2) * (pow(momega, 2) + s) -
               pow(mrho, 2) * (2 * pow(momega, 2) + s)))) /
            (pow(momega, 2) - tmax) +
        (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) + 3 * pow(momega, 8) -
         4 * pow(momega, 6) * s - 7 * pow(momega, 4) * pow(s, 2) +
         4 * pow(momega, 2) * pow(s, 3) + 5 * pow(s, 4) +
         pow(mrho, 4) *
             (pow(momega, 4) - 2 * pow(momega, 2) * s + 2 * pow(s, 2)) +
         pow(mrho, 2) *
             (-4 * pow(momega, 6) + 8 * pow(momega, 4) * s - 6 * pow(s, 3)) -
         2 * pow(mpion, 2) *
             (4 * pow(momega, 6) -
              2 * pow(mrho, 2) * pow(pow(momega, 2) - 2 * s, 2) +
              pow(mrho, 4) * s - 10 * pow(momega, 4) * s + 8 * pow(s, 3)) +
         pow(mpion, 4) *
             (pow(mrho, 4) + 2 * pow(mrho, 2) * (pow(momega, 2) - s) +
              4 * (pow(momega, 4) - 3 * pow(momega, 2) * s + 3 * pow(s, 2)))) *
            tmax -
        2 * pow(mpion, 2) * pow(momega, 4) * pow(tmax, 2) -
        pow(mrho, 2) * pow(momega, 4) * pow(tmax, 2) +
        pow(momega, 6) * pow(tmax, 2) - pow(mpion, 4) * s * pow(tmax, 2) +
        pow(mpion, 2) * pow(mrho, 2) * s * pow(tmax, 2) +
        8 * pow(mpion, 2) * pow(momega, 2) * s * pow(tmax, 2) +
        3 * pow(mrho, 2) * pow(momega, 2) * s * pow(tmax, 2) -
        2 * pow(momega, 4) * s * pow(tmax, 2) -
        8 * pow(mpion, 2) * pow(s, 2) * pow(tmax, 2) -
        3 * pow(mrho, 2) * pow(s, 2) * pow(tmax, 2) -
        3 * pow(momega, 2) * pow(s, 2) * pow(tmax, 2) +
        5 * pow(s, 3) * pow(tmax, 2) +
        ((pow(momega, 4) - 4 * pow(momega, 2) * s + 5 * pow(s, 2)) *
         pow(tmax, 3)) /
            3. -
        (pow(pow(momega, 2) - s, 2) *
         (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
          pow(mpion, 4) *
              (pow(mrho, 4) + 4 * pow(momega, 4) - 2 * pow(momega, 2) * s) +
          pow(momega, 4) *
              (pow(mrho, 4) + pow(momega, 4) + 2 * pow(momega, 2) * s +
               2 * pow(s, 2) - 2 * pow(mrho, 2) * (pow(momega, 2) + s)) -
          2 * pow(mpion, 2) * pow(momega, 2) *
              (pow(mrho, 4) + 2 * pow(momega, 2) * (pow(momega, 2) + s) -
               pow(mrho, 2) * (2 * pow(momega, 2) + s)))) /
            (pow(momega, 2) - tmin) -
        (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) + 3 * pow(momega, 8) -
         4 * pow(momega, 6) * s - 7 * pow(momega, 4) * pow(s, 2) +
         4 * pow(momega, 2) * pow(s, 3) + 5 * pow(s, 4) +
         pow(mrho, 4) *
             (pow(momega, 4) - 2 * pow(momega, 2) * s + 2 * pow(s, 2)) +
         pow(mrho, 2) *
             (-4 * pow(momega, 6) + 8 * pow(momega, 4) * s - 6 * pow(s, 3)) -
         2 * pow(mpion, 2) *
             (4 * pow(momega, 6) -
              2 * pow(mrho, 2) * pow(pow(momega, 2) - 2 * s, 2) +
              pow(mrho, 4) * s - 10 * pow(momega, 4) * s + 8 * pow(s, 3)) +
         pow(mpion, 4) *
             (pow(mrho, 4) + 2 * pow(mrho, 2) * (pow(momega, 2) - s) +
              4 * (pow(momega, 4) - 3 * pow(momega, 2) * s + 3 * pow(s, 2)))) *
            tmin +
        2 * pow(mpion, 2) * pow(momega, 4) * pow(tmin, 2) +
        pow(mrho, 2) * pow(momega, 4) * pow(tmin, 2) -
        pow(momega, 6) * pow(tmin, 2) + pow(mpion, 4) * s * pow(tmin, 2) -
        pow(mpion, 2) * pow(mrho, 2) * s * pow(tmin, 2) -
        8 * pow(mpion, 2) * pow(momega, 2) * s * pow(tmin, 2) -
        3 * pow(mrho, 2) * pow(momega, 2) * s * pow(tmin, 2) +
        2 * pow(momega, 4) * s * pow(tmin, 2) +
        8 * pow(mpion, 2) * pow(s, 2) * pow(tmin, 2) +
        3 * pow(mrho, 2) * pow(s, 2) * pow(tmin, 2) +
        3 * pow(momega, 2) * pow(s, 2) * pow(tmin, 2) -
        5 * pow(s, 3) * pow(tmin, 2) -
        ((pow(momega, 4) - 4 * pow(momega, 2) * s + 5 * pow(s, 2)) *
         pow(tmin, 3)) /
            3. +
        2 * (pow(momega, 2) - s) *
            (-pow(mpion, 8) +
             pow(mpion, 4) * (4 * pow(momega, 4) - 7 * pow(momega, 2) * s +
                              pow(s, 2) + pow(mrho, 2) * (pow(momega, 2) + s)) +
             pow(mpion, 2) *
                 (-6 * pow(momega, 6) + 6 * pow(momega, 4) * s +
                  8 * pow(momega, 2) * pow(s, 2) +
                  pow(mrho, 4) * (-pow(momega, 2) + s) +
                  pow(mrho, 2) * (4 * pow(momega, 4) - 7 * pow(momega, 2) * s -
                                  pow(s, 2))) +
             pow(momega, 2) *
                 (2 * pow(momega, 6) + pow(mrho, 4) * (pow(momega, 2) - s) -
                  4 * pow(momega, 2) * pow(s, 2) - 3 * pow(s, 3) +
                  pow(mrho, 2) * (-3 * pow(momega, 4) + 2 * pow(momega, 2) * s +
                                  3 * pow(s, 2)))) *
            log((-pow(momega, 2) + tmax) / (-pow(momega, 2) + tmin)))) /
      (128. * Pi * pow(pow(momega, 2) - s, 2) *
       (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
        2 * pow(mpion, 2) * (pow(mrho, 2) + s)));

  return xs * to_mb / spin_deg_factor;
}

double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi0_rho0_pi0(
    const double s, const double t, const double m_rho) {
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  const double spin_deg_factor = 3.0;

  double diff_xs =
      (pow(Const, 2) * pow(g_POR, 4) *
       (pow(momega, 4) * pow(s, 4) + 4 * pow(momega, 4) * pow(s, 3) * t -
        4 * pow(momega, 2) * pow(s, 4) * t +
        10 * pow(momega, 4) * pow(s, 2) * pow(t, 2) -
        16 * pow(momega, 2) * pow(s, 3) * pow(t, 2) +
        5 * pow(s, 4) * pow(t, 2) + 4 * pow(momega, 4) * s * pow(t, 3) -
        16 * pow(momega, 2) * pow(s, 2) * pow(t, 3) +
        10 * pow(s, 3) * pow(t, 3) + pow(momega, 4) * pow(t, 4) -
        4 * pow(momega, 2) * s * pow(t, 4) + 5 * pow(s, 2) * pow(t, 4) +
        pow(mpion, 8) * pow(-2 * pow(momega, 2) + s + t, 2) -
        2 * pow(mpion, 6) * pow(mrho, 2) *
            (2 * pow(momega, 4) + pow(s, 2) + pow(t, 2) -
             2 * pow(momega, 2) * (s + t)) +
        pow(mrho, 4) *
            (2 * pow(s, 2) * pow(t, 2) - 2 * pow(momega, 2) * s * t * (s + t) +
             pow(momega, 4) * (pow(s, 2) + pow(t, 2))) -
        2 * pow(mrho, 2) *
            (3 * pow(s, 2) * pow(t, 2) * (s + t) -
             3 * pow(momega, 2) * s * t * pow(s + t, 2) +
             pow(momega, 4) * (pow(s, 3) + 2 * pow(s, 2) * t +
                               2 * s * pow(t, 2) + pow(t, 3))) +
        pow(mpion, 4) *
            (-2 * pow(mrho, 2) * (pow(momega, 2) - s) * (pow(momega, 2) - t) *
                 (s + t) -
             8 * pow(momega, 2) * s * t * (s + t) +
             4 * pow(momega, 4) * (pow(s, 2) + pow(t, 2)) -
             2 * s * t * (pow(s, 2) - 6 * s * t + pow(t, 2)) +
             pow(mrho, 4) * (2 * pow(momega, 4) + pow(s, 2) + pow(t, 2) -
                             2 * pow(momega, 2) * (s + t))) -
        2 * pow(mpion, 2) *
            (2 * (s + t) * pow(-2 * s * t + pow(momega, 2) * (s + t), 2) +
             pow(mrho, 4) * (-4 * pow(momega, 2) * s * t +
                             pow(momega, 4) * (s + t) + s * t * (s + t)) -
             pow(mrho, 2) *
                 (-10 * pow(momega, 2) * s * t * (s + t) +
                  2 * pow(momega, 4) * (pow(s, 2) + 3 * s * t + pow(t, 2)) +
                  s * t * (pow(s, 2) + 8 * s * t + pow(t, 2)))))) /
      (128. * Pi * pow(pow(momega, 2) - s, 2) *
       (pow(pow(mpion, 2) - pow(mrho, 2), 2) -
        2 * (pow(mpion, 2) + pow(mrho, 2)) * s + pow(s, 2)) *
       pow(pow(momega, 2) - t, 2));

  return to_mb * diff_xs / spin_deg_factor;
}

// C15
double
PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_rho_pi0_omega_mediated(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_rho, m_pion_, 0.);
  const double &tmax = t_mandelstam[0];
  const double &tmin = t_mandelstam[1];
  const double spin_deg_factor = 3.0;

  const double xs =
      (0.0024867959858108648 * pow(Const, 2) * pow(g_POR, 4) *
       (pow(mpion, 8) * (1. * tmax - 1. * tmin) +
        pow(mpion, 6) * pow(mrho, 2) * (-2. * tmax + 2. * tmin) +
        pow(mpion, 4) * (pow(mrho, 4) * (1. * tmax - 1. * tmin) +
                         s * (4. * s * tmax - 1. * pow(tmax, 2) -
                              4. * s * tmin + 1. * pow(tmin, 2))) +
        pow(s, 2) *
            (1. * pow(s, 2) * tmax + 1. * s * pow(tmax, 2) +
             0.6666666666666666 * pow(tmax, 3) +
             pow(mrho, 4) * (1. * tmax - 1. * tmin) - 1. * pow(s, 2) * tmin -
             1. * s * pow(tmin, 2) - 0.6666666666666666 * pow(tmin, 3) +
             pow(mrho, 2) * (-2. * s * tmax - 1. * pow(tmax, 2) +
                             2. * s * tmin + 1. * pow(tmin, 2))) +
        pow(mpion, 2) * s *
            (pow(mrho, 4) * (-2. * tmax + 2. * tmin) +
             pow(mrho, 2) * (4. * s * tmax + 1. * pow(tmax, 2) - 4. * s * tmin -
                             1. * pow(tmin, 2)) +
             s * (-4. * s * tmax - 2. * pow(tmax, 2) + 4. * s * tmin +
                  2. * pow(tmin, 2))))) /
      ((pow(pow(momega, 2) - 1. * s, 2) *
        (pow(mpion, 4) + pow(mrho, 4) +
         pow(mpion, 2) * (-2. * pow(mrho, 2) - 2. * s) - 2. * pow(mrho, 2) * s +
         pow(s, 2))));

  return xs * to_mb / spin_deg_factor;
}

double PhotonCrossSection<ComputationMethod::Analytic>::
    xs_diff_pi_rho_pi0_omega_mediated(const double s, const double t,
                                      const double m_rho) {
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  const double spin_deg_factor = 3.0;
  const double diff_xs =
      (0.0024867959858108648 * pow(Const, 2) * pow(g_POR, 4) *
       (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
        pow(mpion, 4) * (pow(mrho, 4) + 4 * pow(s, 2) - 2 * s * t) +
        pow(s, 2) * (pow(mrho, 4) + pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                     2 * pow(mrho, 2) * (s + t)) -
        2 * pow(mpion, 2) * s *
            (pow(mrho, 4) + 2 * s * (s + t) - pow(mrho, 2) * (2 * s + t)))) /
      ((pow(pow(momega, 2) - s, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                                      2 * pow(mpion, 2) * (pow(mrho, 2) + s))));

  return to_mb * diff_xs / spin_deg_factor;
}

// C16
double
PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho_pi_omega_mediated(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_rho, m_pion_, 0.);
  const double &tmin = t_mandelstam[1];
  const double &tmax = t_mandelstam[0];
  const double spin_deg_factor = 3.0;

  const double xs =
      (0.0024868 * pow(Const, 2) * pow(g_POR, 4) *
       ((pow(momega, 8) + pow(mpion, 4) * pow(pow(mpion, 2) - pow(mrho, 2), 2) -
         2 * pow(momega, 6) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
         2 * pow(momega, 2) * pow(mpion, 2) *
             (pow(mrho, 4) + pow(mpion, 2) * s - pow(mrho, 2) * s) +
         pow(momega, 4) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                           4 * pow(mpion, 2) * (pow(mrho, 2) - s) -
                           2 * pow(mrho, 2) * s + 2 * pow(s, 2))) /
            (pow(momega, 2) - tmax) +
        3 * pow(momega, 4) * tmax - 8 * pow(momega, 2) * pow(mpion, 2) * tmax +
        4 * pow(mpion, 4) * tmax - 4 * pow(momega, 2) * pow(mrho, 2) * tmax +
        4 * pow(mpion, 2) * pow(mrho, 2) * tmax + pow(mrho, 4) * tmax +
        4 * pow(momega, 2) * s * tmax - 4 * pow(mpion, 2) * s * tmax -
        2 * pow(mrho, 2) * s * tmax + 2 * pow(s, 2) * tmax +
        pow(momega, 2) * pow(tmax, 2) - 2 * pow(mpion, 2) * pow(tmax, 2) -
        pow(mrho, 2) * pow(tmax, 2) + s * pow(tmax, 2) + pow(tmax, 3) / 3. -
        (pow(momega, 8) + pow(mpion, 4) * pow(pow(mpion, 2) - pow(mrho, 2), 2) -
         2 * pow(momega, 6) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
         2 * pow(momega, 2) * pow(mpion, 2) *
             (pow(mrho, 4) + pow(mpion, 2) * s - pow(mrho, 2) * s) +
         pow(momega, 4) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                           4 * pow(mpion, 2) * (pow(mrho, 2) - s) -
                           2 * pow(mrho, 2) * s + 2 * pow(s, 2))) /
            (pow(momega, 2) - tmin) -
        3 * pow(momega, 4) * tmin + 8 * pow(momega, 2) * pow(mpion, 2) * tmin -
        4 * pow(mpion, 4) * tmin + 4 * pow(momega, 2) * pow(mrho, 2) * tmin -
        4 * pow(mpion, 2) * pow(mrho, 2) * tmin - pow(mrho, 4) * tmin -
        4 * pow(momega, 2) * s * tmin + 4 * pow(mpion, 2) * s * tmin +
        2 * pow(mrho, 2) * s * tmin - 2 * pow(s, 2) * tmin -
        pow(momega, 2) * pow(tmin, 2) + 2 * pow(mpion, 2) * pow(tmin, 2) +
        pow(mrho, 2) * pow(tmin, 2) - s * pow(tmin, 2) - pow(tmin, 3) / 3. +
        2 *
            (2 * pow(momega, 6) -
             3 * pow(momega, 4) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
             pow(mpion, 2) *
                 (pow(mrho, 4) + pow(mpion, 2) * s - pow(mrho, 2) * s) +
             pow(momega, 2) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                               4 * pow(mpion, 2) * (pow(mrho, 2) - s) -
                               2 * pow(mrho, 2) * s + 2 * pow(s, 2))) *
            log(fabs(-pow(momega, 2) + tmax)) -
        2 *
            (2 * pow(momega, 6) -
             3 * pow(momega, 4) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
             pow(mpion, 2) *
                 (pow(mrho, 4) + pow(mpion, 2) * s - pow(mrho, 2) * s) +
             pow(momega, 2) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                               4 * pow(mpion, 2) * (pow(mrho, 2) - s) -
                               2 * pow(mrho, 2) * s + 2 * pow(s, 2))) *
            log(fabs(-pow(momega, 2) + tmin)))) /
      ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
        2 * pow(mpion, 2) * (pow(mrho, 2) + s)));

  return xs * to_mb / spin_deg_factor;
}

double PhotonCrossSection<ComputationMethod::Analytic>::
    xs_diff_pi0_rho_pi_omega_mediated(const double s, const double t,
                                      const double m_rho) {
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  const double spin_deg_factor = 3.0;

  const double diff_xs =
      (0.0024867959858108648 * pow(Const, 2) * pow(g_POR, 4) *
       (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
        pow(mpion, 4) * (pow(mrho, 4) - 2 * (s - 2 * t) * t) +
        pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                     2 * pow(mrho, 2) * (s + t)) -
        2 * pow(mpion, 2) * t *
            (pow(mrho, 4) + 2 * t * (s + t) - pow(mrho, 2) * (s + 2 * t)))) /
      ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
        2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
       pow(pow(momega, 2) - t, 2));

  return to_mb * diff_xs / spin_deg_factor;
}

/*----------------------------------------------------------------------------*/
/*				 Pi + Rho -> Pi + Photon channels mediated by
 * (Pi, Rho, a1)         */
/*  			 and Omega summed
 */
/*----------------------------------------------------------------------------*/

// C12 + C16
double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho_pi(
    const double s, const double m_rho) {
  return xs_pi0_rho_pi_rho_mediated(s, m_rho) +
         xs_pi0_rho_pi_omega_mediated(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi0_rho_pi(
    const double s, const double t, const double m_rho) {
  return xs_diff_pi0_rho_pi_rho_mediated(s, t, m_rho) +
         xs_diff_pi0_rho_pi_omega_mediated(s, t, m_rho);
}

// C13 + C15
double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_rho_pi0(
    const double s, const double m_rho) {
  return xs_pi_rho_pi0_rho_mediated(s, m_rho) +
         xs_pi_rho_pi0_omega_mediated(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_rho_pi0(
    const double s, const double t, const double m_rho) {
  return xs_diff_pi_rho_pi0_rho_mediated(s, t, m_rho) +
         xs_diff_pi_rho_pi0_omega_mediated(s, t, m_rho);
}

/*----------------------------------------------------------------------------*/
/* 					Pi + Pi -> Rho + Photon channels mediated
 * by (Pi, Rho, a1) 				*/
/*----------------------------------------------------------------------------*/
// C21
double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_pi_rho0(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const double s_sqrt = sqrt(s);
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  const double spin_deg_factor = 1.0;

  auto mandelstam_t = get_t_range(s_sqrt, m_pion_, m_pion_, m_rho, 0.);
  double tmax = mandelstam_t[0];
  double tmin = mandelstam_t[1];

  const double xs =
      (-(pow(Const, 2) * pow(ghat, 4) *
         (0. +
          (0.03125 * pow(eta1 - 1. * eta2, 2) *
           (eta1 * eta2 *
                (-2. * pow(ma1, 8) - 2. * pow(mpion, 8) +
                 2. * pow(mpion, 4) * pow(mrho, 4) +
                 pow(ma1, 6) * (8. * pow(mpion, 2) - 4. * s) +
                 pow(ma1, 2) * pow(mpion, 2) *
                     (8. * pow(mpion, 4) - 8. * pow(mrho, 4) -
                      4. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s) +
                 pow(ma1, 4) * (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                                8. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s -
                                4. * pow(s, 2))) +
            pow(eta2, 2) *
                (1. * pow(ma1, 8) + 1. * pow(mpion, 8) -
                 2. * pow(mpion, 6) * pow(mrho, 2) +
                 1. * pow(mpion, 4) * pow(mrho, 4) +
                 pow(ma1, 6) *
                     (-4. * pow(mpion, 2) + 2. * pow(mrho, 2) + 2. * s) +
                 pow(ma1, 4) * (6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                                pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                                2. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                 pow(ma1, 2) * (-4. * pow(mpion, 6) -
                                2. * pow(mpion, 2) * pow(mrho, 2) * s +
                                pow(mpion, 4) * (6. * pow(mrho, 2) + 2. * s))) +
            pow(eta1, 2) *
                (1. * pow(ma1, 8) +
                 pow(ma1, 6) *
                     (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                 pow(ma1, 4) * (6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                                pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                                4. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                 pow(ma1, 2) *
                     (-4. * pow(mpion, 6) +
                      2. * pow(mpion, 2) * pow(mrho, 2) * s +
                      pow(mrho, 2) * (2. * pow(mrho, 2) - 2. * s) * s +
                      pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s)) +
                 pow(mpion, 2) *
                     (1. * pow(mpion, 6) + 2. * pow(mpion, 4) * pow(mrho, 2) -
                      2. * pow(mrho, 6) + 2. * pow(mrho, 4) * s +
                      pow(mpion, 2) *
                          (1. * pow(mrho, 4) - 2. * pow(mrho, 2) * s))))) /
              (1. * pow(ma1, 2) - 1. * tmin) +
          (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
           (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
              (1. * pow(mpion, 2) - 1. * tmin) +
          (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
           (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
              (1. * pow(mpion, 2) + 1. * pow(mrho, 2) - 1. * s - 1. * tmin) -
          (0.5 * pow(-2. + delta, 2) * pow(mpion, 2) * tmin) / pow(mrho, 2) -
          0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
              (-0.5 * eta2 * pow(ma1, 2) + 1. * eta1 * pow(mpion, 2) +
               0.5 * eta2 * pow(mrho, 2) + 0.5 * eta1 * s - 1. * eta2 * s) *
              tmin +
          (0.25 *
           (pow(mpion, 2) *
                (12. + 1. * pow(delta, 2) - 16. * C4 * pow(mrho, 2) +
                 delta * (-8. + 8. * C4 * pow(mrho, 2))) +
            (-4. - 3. * pow(delta, 2) - 16. * C4 * pow(mrho, 2) +
             delta * (8. + 8. * C4 * pow(mrho, 2))) *
                s) *
           tmin) /
              pow(mrho, 2) -
          0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
              (eta2 * (pow(ma1, 2) - 1. * s) +
               eta1 * (-2. * pow(mpion, 2) + s)) *
              tmin -
          0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
              (eta2 * (-1. * pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * s) +
               eta1 * (2. * pow(mpion, 2) + s)) *
              tmin +
          0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
              (eta1 * (1. * pow(mpion, 2) - 0.5 * s) +
               eta2 * (-0.5 * pow(ma1, 2) - 1. * pow(mpion, 2) -
                       0.5 * pow(mrho, 2) + 1. * s)) *
              tmin +
          (0.25 * (-2. + 1. * delta) *
           (-8. * C4 * pow(mrho, 4) +
            pow(mpion, 2) * (2. + 1. * delta - 8. * C4 * pow(mrho, 2)) +
            (-2. - 3. * delta) * s +
            pow(mrho, 2) * (2. + 1. * delta + 16. * C4 * s)) *
           tmin) /
              pow(mrho, 2) +
          (0.25 *
           (32 * pow(C4, 2) * pow(mrho, 8) + 2 * pow(delta, 2) * pow(s, 2) +
            8 * C4 * pow(mrho, 6) * (-6 + delta - 8 * C4 * s) +
            2 * delta * pow(mrho, 2) * s * (-6 + delta - 8 * C4 * s) +
            pow(mrho, 4) * (12 - pow(delta, 2) + 8 * C4 * (6 + delta) * s +
                            32 * pow(C4, 2) * pow(s, 2))) *
           tmin) /
              pow(mrho, 4) -
          (1. * (1. * eta1 - 1. * eta2) *
           (eta2 *
                (0.75 * pow(mrho, 4) - 0.125 * delta * pow(mrho, 4) -
                 1. * C4 * pow(mrho, 6) +
                 pow(ma1, 2) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                 pow(mpion, 2) * (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4)) -
                 0.25 * pow(mrho, 2) * s - 0.375 * delta * pow(mrho, 2) * s +
                 2. * C4 * pow(mrho, 4) * s + 0.25 * delta * pow(s, 2) -
                 1. * C4 * pow(mrho, 2) * pow(s, 2)) +
            eta1 *
                (0.5 * pow(mrho, 4) - 1. * C4 * pow(mrho, 6) +
                 pow(mpion, 2) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                 pow(ma1, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                 0.25 * delta * pow(s, 2) +
                 1. * C4 * pow(mrho, 2) * pow(s, 2))) *
           tmin) /
              pow(mrho, 2) +
          0.0625 * pow(eta1 - 1. * eta2, 2) *
              (pow(eta2, 2) *
                   (pow(Gammaa1, 2) * pow(ma1, 2) - 1. * pow(ma1, 4) -
                    2. * pow(mpion, 4) - 2. * pow(mpion, 2) * pow(mrho, 2) +
                    2. * pow(mrho, 4) +
                    pow(ma1, 2) * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s) -
                    3. * pow(mrho, 2) * s + pow(s, 2)) +
               pow(eta1, 2) *
                   (pow(Gammaa1, 2) * pow(ma1, 2) - 1. * pow(ma1, 4) -
                    2. * pow(mpion, 4) - 2. * pow(mpion, 2) * pow(mrho, 2) +
                    pow(ma1, 2) * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s) +
                    pow(mrho, 2) * s + pow(s, 2)) +
               eta1 * eta2 *
                   (-2. * pow(Gammaa1, 2) * pow(ma1, 2) + 2. * pow(ma1, 4) +
                    4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                    2. * pow(mrho, 4) - 2. * pow(s, 2) +
                    pow(ma1, 2) *
                        (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s))) *
              tmin +
          0.03125 * pow(eta1 - 1. * eta2, 2) *
              (eta1 * eta2 *
                   (-6. * pow(ma1, 4) - 12. * pow(mpion, 4) +
                    2. * pow(mrho, 4) +
                    pow(ma1, 2) * (16. * pow(mpion, 2) - 8. * s) +
                    8. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s -
                    4. * pow(s, 2)) +
               pow(eta1, 2) *
                   (3. * pow(ma1, 4) + 6. * pow(mpion, 4) + pow(mrho, 4) +
                    pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                    4. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                    pow(ma1, 2) *
                        (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) +
               pow(eta2, 2) *
                   (3. * pow(ma1, 4) + 6. * pow(mpion, 4) + pow(mrho, 4) +
                    pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                    2. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                    pow(ma1, 2) *
                        (-8. * pow(mpion, 2) + 4. * pow(mrho, 2) + 4. * s))) *
              tmin -
          (0.125 * (-2. + 1. * delta) *
           (2. + 1. * delta - 8. * C4 * pow(mrho, 2)) * pow(tmin, 2)) /
              pow(mrho, 2) -
          0.5 * pow(1. * eta1 - 1. * eta2, 2) *
              (-0.5 + 1. * C4 * pow(mrho, 2)) * pow(tmin, 2) -
          (1. *
           (0.5 - 0.125 * pow(delta, 2) - 2. * C4 * pow(mrho, 2) +
            1. * C4 * delta * pow(mrho, 2)) *
           pow(tmin, 2)) /
              pow(mrho, 2) +
          0.0625 * pow(1. * eta1 - 1. * eta2, 4) *
              (1. * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 0.5 * s) *
              pow(tmin, 2) +
          0.03125 * pow(eta1 - 1. * eta2, 3) *
              (eta2 * (-1. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                       1. * pow(mrho, 2) - 1. * s) +
               eta1 *
                   (pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
              pow(tmin, 2) +
          0.010416666666666666 * pow(eta1 - 1. * eta2, 4) * pow(tmin, 3) -
          0.020833333333333332 * pow(1. * eta1 - 1. * eta2, 4) * pow(tmin, 3) +
          0.03125 * pow(eta1 - 1. * eta2, 2) *
              (eta1 * eta2 *
                   (2. * pow(Gammaa1, 2) * pow(ma1, 2) - 6. * pow(ma1, 4) -
                    4. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                    pow(ma1, 2) * (8. * pow(mpion, 2) - 8. * s) +
                    4. * pow(mrho, 2) * s - 4. * pow(s, 2)) +
               pow(eta1, 2) *
                   (-1. * pow(Gammaa1, 2) * pow(ma1, 2) + 3. * pow(ma1, 4) +
                    2. * pow(mpion, 4) + 2. * pow(mpion, 2) * pow(mrho, 2) +
                    pow(mrho, 4) - 4. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                    pow(ma1, 2) *
                        (-4. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) +
               pow(eta2, 2) *
                   (-1. * pow(Gammaa1, 2) * pow(ma1, 2) + 3. * pow(ma1, 4) +
                    2. * pow(mpion, 4) - 2. * pow(mpion, 2) * pow(mrho, 2) +
                    pow(mrho, 4) - 2. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                    pow(ma1, 2) *
                        (-4. * pow(mpion, 2) + 4. * pow(mrho, 2) + 4. * s))) *
              (-1. * pow(mrho, 2) + s + tmin) -
          0.03125 * pow(eta1 - 1. * eta2, 3) *
              (eta2 * (-1. * pow(ma1, 2) - 1. * pow(mrho, 2) - 1. * s) +
               eta1 * (pow(ma1, 2) - 1. * pow(mrho, 2) + s)) *
              pow(-1. * pow(mrho, 2) + s + tmin, 2) +
          0.010416666666666666 * pow(eta1 - 1. * eta2, 4) *
              pow(-1. * pow(mrho, 2) + s + tmin, 3) +
          0.25 * (eta1 - 1. * eta2) * (1. * eta1 - 1. * eta2) *
              (-1. + 2. * C4 * pow(mrho, 2)) *
              pow(pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                      tmin,
                  2) -
          (2. * (1. * eta1 - 1. * eta2) *
           (eta2 *
                (0.375 * pow(mrho, 4) - 0.0625 * delta * pow(mrho, 4) -
                 0.5 * C4 * pow(mrho, 6) +
                 pow(ma1, 2) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                 pow(mpion, 2) *
                     (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                 0.125 * pow(mrho, 2) * s - 0.1875 * delta * pow(mrho, 2) * s +
                 1. * C4 * pow(mrho, 4) * s + 0.125 * delta * pow(s, 2) -
                 0.5 * C4 * pow(mrho, 2) * pow(s, 2)) +
            eta1 *
                (0.25 * pow(mrho, 4) - 0.5 * C4 * pow(mrho, 6) +
                 pow(mpion, 2) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                 pow(ma1, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                 0.125 * delta * pow(s, 2) +
                 0.5 * C4 * pow(mrho, 2) * pow(s, 2))) *
           (1. * pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + 1. * s +
            1. * tmin)) /
              pow(mrho, 2) +
          (2. * (1. * eta1 - 1. * eta2) * Gammaa1 * ma1 *
           (eta2 *
                (0.375 * pow(mrho, 4) - 0.0625 * delta * pow(mrho, 4) -
                 0.5 * C4 * pow(mrho, 6) +
                 pow(ma1, 2) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                 pow(mpion, 2) *
                     (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                 0.125 * pow(mrho, 2) * s - 0.1875 * delta * pow(mrho, 2) * s +
                 1. * C4 * pow(mrho, 4) * s + 0.125 * delta * pow(s, 2) -
                 0.5 * C4 * pow(mrho, 2) * pow(s, 2)) +
            eta1 *
                (0.25 * pow(mrho, 4) - 0.5 * C4 * pow(mrho, 6) +
                 pow(mpion, 2) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                 pow(ma1, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                 0.125 * delta * pow(s, 2) +
                 0.5 * C4 * pow(mrho, 2) * pow(s, 2))) *
           atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                 tmin) /
                (Gammaa1 * ma1))) /
              pow(mrho, 2) +
          (0.25 * (-2. + delta) * (eta1 - 1. * eta2) * Gammaa1 * ma1 *
           (eta2 * (-1. * pow(ma1, 6) +
                    pow(Gammaa1, 2) * pow(ma1, 2) *
                        (-1. * pow(ma1, 2) + 0.5 * pow(mrho, 2) - 1. * s) +
                    pow(ma1, 4) *
                        (2. * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 1. * s) +
                    pow(mpion, 4) * (-1.5 * pow(mrho, 2) + 1. * s) +
                    pow(ma1, 2) * pow(mpion, 2) *
                        (-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + 2. * s)) +
            eta1 *
                (pow(Gammaa1, 2) * pow(ma1, 2) *
                     (1. * pow(mpion, 2) + 0.5 * s) +
                 pow(ma1, 4) * (1. * pow(mpion, 2) + 0.5 * s) +
                 pow(ma1, 2) * (-2. * pow(mpion, 4) - 1. * pow(mpion, 2) * s) +
                 pow(mpion, 2) *
                     (1. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                      pow(mpion, 2) * (2. * pow(mrho, 2) - 1.5 * s) +
                      1. * pow(mrho, 2) * s))) *
           atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                 tmin) /
                (Gammaa1 * ma1))) /
              (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
               2. * pow(ma1, 2) * pow(mpion, 2) + pow(mpion, 4)) -
          (0.25 * (-2. + delta) * (eta1 - 1. * eta2) * Gammaa1 * ma1 *
           (eta2 *
                (-1. * pow(ma1, 6) - 2. * pow(mpion, 4) * pow(mrho, 2) -
                 1. * pow(mpion, 2) * pow(mrho, 4) +
                 pow(ma1, 4) *
                     (2. * pow(mpion, 2) + 2. * pow(mrho, 2) - 1.5 * s) +
                 2.5 * pow(mpion, 4) * s +
                 3. * pow(mpion, 2) * pow(mrho, 2) * s +
                 0.5 * pow(mrho, 4) * s - 2. * pow(mpion, 2) * pow(s, 2) -
                 1. * pow(mrho, 2) * pow(s, 2) + 0.5 * pow(s, 3) +
                 pow(Gammaa1, 2) * (-1. * pow(ma1, 4) + 0.5 * pow(ma1, 2) * s) +
                 pow(ma1, 2) *
                     (-1. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                      1. * pow(mrho, 2) * s +
                      pow(mpion, 2) * (-2. * pow(mrho, 2) + 1. * s))) +
            eta1 *
                (1. * pow(mpion, 6) + 4. * pow(mpion, 4) * pow(mrho, 2) +
                 1. * pow(mpion, 2) * pow(mrho, 4) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (1. * pow(mpion, 2) - 0.5 * s) +
                 pow(ma1, 4) * (1. * pow(mpion, 2) - 0.5 * s) -
                 4.5 * pow(mpion, 4) * s -
                 4. * pow(mpion, 2) * pow(mrho, 2) * s -
                 0.5 * pow(mrho, 4) * s + 3. * pow(mpion, 2) * pow(s, 2) +
                 1. * pow(mrho, 2) * pow(s, 2) - 0.5 * pow(s, 3) +
                 pow(ma1, 2) *
                     (-2. * pow(mpion, 4) + (1. * pow(mrho, 2) - 1. * s) * s +
                      pow(mpion, 2) * (-2. * pow(mrho, 2) + 3. * s)))) *
           atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                 tmin) /
                (Gammaa1 * ma1))) /
              (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) + pow(mpion, 4) +
               2. * pow(mpion, 2) * pow(mrho, 2) + pow(mrho, 4) -
               2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s + pow(s, 2) +
               pow(ma1, 2) *
                   (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) +
          (0.03125 * pow(eta1 - 1. * eta2, 2) *
           (pow(eta2, 2) *
                (pow(Gammaa1, 4) * pow(ma1, 4) + pow(ma1, 8) + pow(mpion, 8) -
                 2. * pow(mpion, 6) * pow(mrho, 2) +
                 pow(mpion, 4) * pow(mrho, 4) +
                 pow(ma1, 6) *
                     (-4. * pow(mpion, 2) + 2. * pow(mrho, 2) + 2. * s) +
                 pow(ma1, 4) * (6. * pow(mpion, 4) + pow(mrho, 4) +
                                pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                                2. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                 pow(ma1, 2) * (-4. * pow(mpion, 6) -
                                2. * pow(mpion, 2) * pow(mrho, 2) * s +
                                pow(mpion, 4) * (6. * pow(mrho, 2) + 2. * s)) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (-6. * pow(ma1, 4) - 6. * pow(mpion, 4) -
                      1. * pow(mrho, 4) +
                      pow(ma1, 2) *
                          (12. * pow(mpion, 2) - 6. * pow(mrho, 2) - 6. * s) +
                      2. * pow(mrho, 2) * s - 2. * pow(s, 2) +
                      pow(mpion, 2) * (6. * pow(mrho, 2) + 4. * s))) +
            eta1 * eta2 *
                (-2. * pow(Gammaa1, 4) * pow(ma1, 4) - 2. * pow(ma1, 8) -
                 2. * pow(mpion, 8) + 2. * pow(mpion, 4) * pow(mrho, 4) +
                 pow(ma1, 6) * (8. * pow(mpion, 2) - 4. * s) +
                 pow(ma1, 2) * pow(mpion, 2) *
                     (8. * pow(mpion, 4) - 8. * pow(mrho, 4) -
                      4. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s) +
                 pow(ma1, 4) * (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                                8. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s -
                                4. * pow(s, 2)) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (12. * pow(ma1, 4) + 12. * pow(mpion, 4) -
                      2. * pow(mrho, 4) - 8. * pow(mpion, 2) * s -
                      4. * pow(mrho, 2) * s + 4. * pow(s, 2) +
                      pow(ma1, 2) * (-24. * pow(mpion, 2) + 12. * s))) +
            pow(eta1, 2) *
                (pow(Gammaa1, 4) * pow(ma1, 4) + pow(ma1, 8) +
                 pow(ma1, 6) *
                     (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                 pow(ma1, 4) * (6. * pow(mpion, 4) + pow(mrho, 4) +
                                pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                                4. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                 pow(ma1, 2) *
                     (-4. * pow(mpion, 6) +
                      2. * pow(mpion, 2) * pow(mrho, 2) * s +
                      pow(mrho, 2) * (2. * pow(mrho, 2) - 2. * s) * s +
                      pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s)) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (-6. * pow(ma1, 4) - 6. * pow(mpion, 4) -
                      1. * pow(mrho, 4) +
                      pow(ma1, 2) *
                          (12. * pow(mpion, 2) + 6. * pow(mrho, 2) - 6. * s) +
                      4. * pow(mrho, 2) * s - 2. * pow(s, 2) +
                      pow(mpion, 2) * (-6. * pow(mrho, 2) + 4. * s)) +
                 pow(mpion, 2) *
                     (pow(mpion, 6) + 2. * pow(mpion, 4) * pow(mrho, 2) -
                      2. * pow(mrho, 6) + 2. * pow(mrho, 4) * s +
                      pow(mpion, 2) *
                          (pow(mrho, 4) - 2. * pow(mrho, 2) * s)))) *
           atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                 tmin) /
                (Gammaa1 * ma1))) /
              (Gammaa1 * ma1) -
          (0.0625 * pow(eta1 - 1. * eta2, 2) * Gammaa1 * ma1 *
           (eta1 * eta2 *
                (-2. * pow(Gammaa1, 4) * pow(ma1, 4) + 14. * pow(ma1, 8) +
                 14. * pow(mpion, 8) + 28. * pow(mpion, 6) * pow(mrho, 2) +
                 20. * pow(mpion, 4) * pow(mrho, 4) +
                 10. * pow(mpion, 2) * pow(mrho, 6) + 2. * pow(mrho, 8) -
                 16. * pow(mpion, 6) * s -
                 16. * pow(mpion, 4) * pow(mrho, 2) * s -
                 12. * pow(mpion, 2) * pow(mrho, 4) * s -
                 4. * pow(mrho, 6) * s - 4. * pow(mpion, 4) * pow(s, 2) -
                 6. * pow(mpion, 2) * pow(mrho, 2) * pow(s, 2) +
                 8. * pow(mpion, 2) * pow(s, 3) +
                 4. * pow(mrho, 2) * pow(s, 3) - 2. * pow(s, 4) +
                 pow(ma1, 6) *
                     (-56. * pow(mpion, 2) - 28. * pow(mrho, 2) + 28. * s) +
                 pow(ma1, 4) * (84. * pow(mpion, 4) + 24. * pow(mrho, 4) +
                                pow(mpion, 2) * (84. * pow(mrho, 2) - 72. * s) -
                                36. * pow(mrho, 2) * s + 12. * pow(s, 2)) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (-4. * pow(ma1, 4) - 4. * pow(mpion, 4) +
                      pow(ma1, 2) *
                          (8. * pow(mpion, 2) + 4. * pow(mrho, 2) - 4. * s) +
                      (4. * pow(mrho, 2) - 4. * s) * s +
                      pow(mpion, 2) * (-4. * pow(mrho, 2) + 8. * s)) +
                 pow(ma1, 2) *
                     (-56. * pow(mpion, 6) - 10. * pow(mrho, 6) +
                      18. * pow(mrho, 4) * s - 6. * pow(mrho, 2) * pow(s, 2) -
                      2. * pow(s, 3) +
                      pow(mpion, 4) * (-84. * pow(mrho, 2) + 60. * s) +
                      pow(mpion, 2) *
                          (-48. * pow(mrho, 4) + 60. * pow(mrho, 2) * s -
                           12. * pow(s, 2)))) +
            pow(eta1, 2) *
                (1. * pow(Gammaa1, 4) * pow(ma1, 4) - 7. * pow(ma1, 8) -
                 7. * pow(mpion, 8) - 14. * pow(mpion, 6) * pow(mrho, 2) -
                 7. * pow(mpion, 4) * pow(mrho, 4) -
                 2. * pow(mpion, 2) * pow(mrho, 6) +
                 pow(ma1, 6) *
                     (28. * pow(mpion, 2) + 14. * pow(mrho, 2) - 14. * s) +
                 8. * pow(mpion, 6) * s +
                 11. * pow(mpion, 4) * pow(mrho, 2) * s +
                 6. * pow(mpion, 2) * pow(mrho, 4) * s + 1. * pow(mrho, 6) * s +
                 2. * pow(mpion, 4) * pow(s, 2) -
                 1. * pow(mrho, 4) * pow(s, 2) -
                 4. * pow(mpion, 2) * pow(s, 3) -
                 1. * pow(mrho, 2) * pow(s, 3) + 1. * pow(s, 4) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (2. * pow(ma1, 4) + 2. * pow(mpion, 4) +
                      1. * pow(mrho, 4) +
                      pow(mpion, 2) * (2. * pow(mrho, 2) - 4. * s) -
                      1. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                      pow(ma1, 2) *
                          (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) +
                 pow(ma1, 4) *
                     (-42. * pow(mpion, 4) - 9. * pow(mrho, 4) +
                      21. * pow(mrho, 2) * s - 6. * pow(s, 2) +
                      pow(mpion, 2) * (-42. * pow(mrho, 2) + 36. * s)) +
                 pow(ma1, 2) * (28. * pow(mpion, 6) + 2. * pow(mrho, 6) +
                                pow(mpion, 4) * (42. * pow(mrho, 2) - 30. * s) -
                                9. * pow(mrho, 4) * s +
                                6. * pow(mrho, 2) * pow(s, 2) + 1. * pow(s, 3) +
                                pow(mpion, 2) * (18. * pow(mrho, 4) -
                                                 36. * pow(mrho, 2) * s +
                                                 6. * pow(s, 2)))) +
            pow(eta2, 2) *
                (1. * pow(Gammaa1, 4) * pow(ma1, 4) - 7. * pow(ma1, 8) -
                 7. * pow(mpion, 8) - 14. * pow(mpion, 6) * pow(mrho, 2) -
                 1. * pow(mpion, 4) * pow(mrho, 4) +
                 6. * pow(mpion, 2) * pow(mrho, 6) + 2. * pow(mrho, 8) +
                 pow(ma1, 6) *
                     (28. * pow(mpion, 2) + 14. * pow(mrho, 2) - 14. * s) +
                 8. * pow(mpion, 6) * s -
                 1. * pow(mpion, 4) * pow(mrho, 2) * s -
                 16. * pow(mpion, 2) * pow(mrho, 4) * s -
                 7. * pow(mrho, 6) * s + 2. * pow(mpion, 4) * pow(s, 2) +
                 14. * pow(mpion, 2) * pow(mrho, 2) * pow(s, 2) +
                 9. * pow(mrho, 4) * pow(s, 2) -
                 4. * pow(mpion, 2) * pow(s, 3) -
                 5. * pow(mrho, 2) * pow(s, 3) + 1. * pow(s, 4) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (2. * pow(ma1, 4) + 2. * pow(mpion, 4) +
                      3. * pow(mrho, 4) +
                      pow(mpion, 2) * (2. * pow(mrho, 2) - 4. * s) -
                      5. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                      pow(ma1, 2) *
                          (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) +
                 pow(ma1, 4) *
                     (-42. * pow(mpion, 4) - 3. * pow(mrho, 4) +
                      9. * pow(mrho, 2) * s - 6. * pow(s, 2) +
                      pow(mpion, 2) * (-42. * pow(mrho, 2) + 36. * s)) +
                 pow(ma1, 2) * (28. * pow(mpion, 6) - 4. * pow(mrho, 6) +
                                pow(mpion, 4) * (42. * pow(mrho, 2) - 30. * s) +
                                9. * pow(mrho, 4) * s -
                                6. * pow(mrho, 2) * pow(s, 2) + 1. * pow(s, 3) +
                                pow(mpion, 2) * (6. * pow(mrho, 4) -
                                                 12. * pow(mrho, 2) * s +
                                                 6. * pow(s, 2))))) *
           atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                 tmin) /
                (Gammaa1 * ma1))) /
              (pow(Gammaa1, 2) * pow(ma1, 2) + 4. * pow(ma1, 4) +
               4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
               pow(mrho, 4) - 4. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s +
               pow(s, 2) +
               pow(ma1, 2) *
                   (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) +
          0.0625 * pow(eta1 - 1. * eta2, 2) *
              (eta1 * eta2 *
                   (-4. * pow(ma1, 6) +
                    pow(ma1, 4) * (12. * pow(mpion, 2) - 6. * s) +
                    pow(mpion, 2) *
                        (4. * pow(mpion, 4) - 4. * pow(mrho, 4) -
                         2. * pow(mpion, 2) * s + 2. * pow(mrho, 2) * s) +
                    pow(ma1, 2) * (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                                   8. * pow(mpion, 2) * s +
                                   4. * pow(mrho, 2) * s - 4. * pow(s, 2))) +
               pow(eta1, 2) *
                   (2. * pow(ma1, 6) - 2. * pow(mpion, 6) +
                    pow(mpion, 2) * pow(mrho, 2) * s +
                    pow(mrho, 2) * (pow(mrho, 2) - 1. * s) * s +
                    pow(mpion, 4) * (-3. * pow(mrho, 2) + s) +
                    pow(ma1, 4) *
                        (-6. * pow(mpion, 2) - 3. * pow(mrho, 2) + 3. * s) +
                    pow(ma1, 2) *
                        (6. * pow(mpion, 4) + pow(mrho, 4) +
                         pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                         4. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
               pow(eta2, 2) *
                   (2. * pow(ma1, 6) - 2. * pow(mpion, 6) -
                    1. * pow(mpion, 2) * pow(mrho, 2) * s +
                    pow(mpion, 4) * (3. * pow(mrho, 2) + s) +
                    pow(ma1, 4) *
                        (-6. * pow(mpion, 2) + 3. * pow(mrho, 2) + 3. * s) +
                    pow(ma1, 2) *
                        (6. * pow(mpion, 4) + pow(mrho, 4) +
                         pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                         2. * pow(mrho, 2) * s + 2. * pow(s, 2)))) *
              log(fabs(-1. * pow(ma1, 2) + tmin)) -
          (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
           (eta2 * (-0.5 * pow(ma1, 6) - 0.5 * pow(mpion, 6) +
                    0.5 * pow(mpion, 4) * pow(mrho, 2) +
                    pow(ma1, 4) *
                        (0.5 * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 1. * s) +
                    pow(ma1, 2) * pow(mpion, 2) *
                        (0.5 * pow(mpion, 2) + 1. * pow(mrho, 2) - 1. * s)) +
            eta1 * (pow(ma1, 4) * (1. * pow(mpion, 2) + 0.5 * s) +
                    pow(mpion, 2) *
                        (1. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                         pow(mpion, 2) * (-1. * pow(mrho, 2) + 0.5 * s) -
                         0.5 * pow(mrho, 2) * s) +
                    pow(ma1, 2) *
                        (-2. * pow(mpion, 4) - 0.5 * pow(mrho, 2) * s +
                         pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s)))) *
           log(fabs(-1. * pow(ma1, 2) + tmin))) /
              (pow(ma1, 2) - 1. * pow(mpion, 2)) +
          (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
           (eta2 *
                (-0.5 * pow(ma1, 6) - 0.5 * pow(mpion, 6) +
                 pow(ma1, 4) * (0.5 * pow(mpion, 2) + 0.5 * pow(mrho, 2)) +
                 pow(mpion, 4) * (0.5 * pow(mrho, 2) - 1. * s) +
                 pow(mpion, 2) * (-0.5 * pow(mrho, 2) + 0.5 * s) * s +
                 pow(ma1, 2) * (0.5 * pow(mpion, 4) +
                                pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                                (-0.5 * pow(mrho, 2) + 0.5 * s) * s)) +
            eta1 * (1. * pow(mpion, 6) +
                    pow(ma1, 4) * (1. * pow(mpion, 2) - 0.5 * s) +
                    pow(mpion, 2) * (1.5 * pow(mrho, 2) - 2. * s) * s +
                    (-0.5 * pow(mrho, 2) + 0.5 * s) * pow(s, 2) +
                    pow(mpion, 4) * (-1. * pow(mrho, 2) + 1.5 * s) +
                    pow(ma1, 2) *
                        (-2. * pow(mpion, 4) + 0.5 * pow(mrho, 2) * s +
                         pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s)))) *
           log(fabs(-1. * pow(ma1, 2) + tmin))) /
              (1. * pow(ma1, 2) - 1. * pow(mpion, 2) - 1. * pow(mrho, 2) +
               1. * s) -
          (0.03125 * pow(eta1 - 1. * eta2, 2) *
           (1. * pow(ma1, 2) - 1. * pow(mpion, 2) - 0.5 * pow(mrho, 2) +
            0.5 * s) *
           (eta1 * eta2 *
                (-2. * pow(ma1, 8) +
                 pow(ma1, 6) *
                     (8. * pow(mpion, 2) + 4. * pow(mrho, 2) - 4. * s) +
                 pow(ma1, 4) *
                     (-12. * pow(mpion, 4) - 4. * pow(mrho, 4) +
                      4. * pow(mrho, 2) * s +
                      pow(mpion, 2) * (-12. * pow(mrho, 2) + 8. * s)) +
                 pow(mpion, 2) *
                     (-2. * pow(mpion, 6) - 4. * pow(mpion, 4) * pow(mrho, 2) -
                      2. * pow(mrho, 6) + 4. * pow(mrho, 4) * s -
                      2. * pow(mrho, 2) * pow(s, 2) +
                      pow(mpion, 2) *
                          (-8. * pow(mrho, 4) + 8. * pow(mrho, 2) * s)) +
                 pow(ma1, 2) * (8. * pow(mpion, 6) + 2. * pow(mrho, 6) +
                                pow(mpion, 4) * (12. * pow(mrho, 2) - 4. * s) -
                                2. * pow(mrho, 4) * s -
                                2. * pow(mrho, 2) * pow(s, 2) + 2. * pow(s, 3) +
                                pow(mpion, 2) *
                                    (8. * pow(mrho, 4) - 4. * pow(mrho, 2) * s -
                                     4. * pow(s, 2)))) +
            pow(eta2, 2) *
                (pow(ma1, 8) +
                 pow(ma1, 6) *
                     (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                 pow(mpion, 4) *
                     (pow(mpion, 4) + 2. * pow(mpion, 2) * pow(mrho, 2) +
                      pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
                 pow(ma1, 4) * (6. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                                pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) +
                                pow(mrho, 2) * s) +
                 pow(ma1, 2) * (-4. * pow(mpion, 6) + 2. * pow(mrho, 6) -
                                5. * pow(mrho, 4) * s +
                                4. * pow(mrho, 2) * pow(s, 2) - 1. * pow(s, 3) +
                                pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s) +
                                pow(mpion, 2) *
                                    (2. * pow(mrho, 4) - 4. * pow(mrho, 2) * s +
                                     2. * pow(s, 2)))) +
            pow(eta1, 2) *
                (pow(ma1, 8) + pow(mpion, 8) +
                 2. * pow(mpion, 6) * pow(mrho, 2) +
                 pow(mpion, 2) * pow(mrho, 2) * s *
                     (-2. * pow(mrho, 2) + 2. * s) +
                 pow(ma1, 6) *
                     (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                 pow(mpion, 4) * (3. * pow(mrho, 4) - 5. * pow(mrho, 2) * s) +
                 pow(ma1, 4) * (6. * pow(mpion, 4) + pow(mrho, 4) +
                                pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                                3. * pow(mrho, 2) * s) +
                 pow(ma1, 2) *
                     (-4. * pow(mpion, 6) + pow(mrho, 4) * s - 1. * pow(s, 3) +
                      pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s) +
                      pow(mpion, 2) *
                          (-2. * pow(mrho, 4) + 4. * pow(mrho, 2) * s +
                           2. * pow(s, 2))))) *
           log(fabs(-1. * pow(ma1, 2) + tmin))) /
              (0.25 * pow(Gammaa1, 2) * pow(ma1, 2) + 1. * pow(ma1, 4) +
               1. * pow(mpion, 4) + 1. * pow(mpion, 2) * pow(mrho, 2) +
               0.25 * pow(mrho, 4) - 1. * pow(mpion, 2) * s -
               0.5 * pow(mrho, 2) * s + 0.25 * pow(s, 2) +
               pow(ma1, 2) *
                   (-2. * pow(mpion, 2) - 1. * pow(mrho, 2) + 1. * s)) -
          (1. * (1. * eta1 - 1. * eta2) *
           (eta2 *
                (pow(ma1, 4) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                 pow(mpion, 2) * pow(mrho, 2) *
                     (pow(mpion, 2) * (0.5 - 1. * C4 * pow(mrho, 2)) +
                      (-0.25 + 0.125 * delta) * (pow(mrho, 2) + s)) +
                 pow(ma1, 2) *
                     (-1. * C4 * pow(mrho, 6) +
                      pow(mpion, 2) *
                          (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4)) +
                      0.25 * delta * pow(s, 2) +
                      pow(mrho, 2) * s * (-0.25 - 0.375 * delta - 1. * C4 * s) +
                      pow(mrho, 4) * (0.75 - 0.125 * delta + 2. * C4 * s))) +
            eta1 *
                (pow(ma1, 4) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) +
                 pow(ma1, 2) * (0.5 * pow(mrho, 4) - 1. * C4 * pow(mrho, 6) +
                                pow(mpion, 2) * (1. * pow(mrho, 2) -
                                                 2. * C4 * pow(mrho, 4)) -
                                0.25 * delta * pow(s, 2) +
                                1. * C4 * pow(mrho, 2) * pow(s, 2)) +
                 pow(mrho, 2) *
                     (pow(mpion, 4) * (-0.5 + 1. * C4 * pow(mrho, 2)) +
                      s * ((0.25 - 0.125 * delta) * pow(mrho, 2) +
                           (-0.25 + 0.125 * delta) * s) +
                      pow(mpion, 2) *
                          (2. * C4 * pow(mrho, 4) + (0.5 + 0.25 * delta) * s +
                           pow(mrho, 2) * (-1. - 2. * C4 * s))))) *
           log(fabs(-1. * pow(ma1, 2) + tmin))) /
              pow(mrho, 2) +
          0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
              log(fabs(-1. * pow(mpion, 2) + tmin)) +
          (0.25 *
           (0. +
            8.000000000000002 * pow(2. - 1. * delta, 2) * pow(mpion, 4) *
                pow(mrho, 2) -
            5.999999999999999 * pow(2. - 1. * delta, 2) * pow(mpion, 2) *
                pow(mrho, 2) * s +
            1. * pow(2. - 1. * delta, 2) * pow(mrho, 2) * pow(s, 2)) *
           log(fabs(-1. * pow(mpion, 2) + tmin))) /
              (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
          (0.125 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
           (0. + eta2 * pow(mpion, 2) * (4. * pow(mrho, 2) - 4. * s) +
            eta1 * (2. * pow(mrho, 4) - 2. * pow(mrho, 2) * s +
                    pow(mpion, 2) * (-4. * pow(mrho, 2) + 4. * s))) *
           log(fabs(-1. * pow(mpion, 2) + tmin))) /
              (pow(ma1, 2) - 1. * pow(mpion, 2)) +
          (2. * (-2. + 1. * delta) *
           (0. + (-0.25 + 0.125 * delta) * pow(mrho, 2) * s +
            pow(mpion, 2) * (-2. * C4 * pow(mrho, 4) - 0.5 * delta * s +
                             pow(mrho, 2) * (1. + 2. * C4 * s))) *
           log(fabs(-1. * pow(mpion, 2) + tmin))) /
              pow(mrho, 2) -
          (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
           (eta1 *
                (pow(mpion, 4) * (-1. * pow(mrho, 2) + 1. * s) +
                 pow(ma1, 2) * (pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                                (-0.5 * pow(mrho, 2) + 0.5 * s) * s) +
                 pow(mpion, 2) * (-1. * pow(mrho, 4) + 2.5 * pow(mrho, 2) * s -
                                  1.5 * pow(s, 2)) +
                 s * (0.5 * pow(mrho, 4) - 1. * pow(mrho, 2) * s +
                      0.5 * pow(s, 2))) +
            eta2 * (0.5 * pow(mrho, 6) +
                    pow(mpion, 4) * (1. * pow(mrho, 2) - 1. * s) -
                    1.5 * pow(mrho, 4) * s + 1.5 * pow(mrho, 2) * pow(s, 2) -
                    0.5 * pow(s, 3) +
                    pow(mpion, 2) * (1.5 * pow(mrho, 4) -
                                     3. * pow(mrho, 2) * s + 1.5 * pow(s, 2)) +
                    pow(ma1, 2) *
                        (-0.5 * pow(mrho, 4) + 1. * pow(mrho, 2) * s -
                         0.5 * pow(s, 2) +
                         pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s)))) *
           log(fabs(-1. * pow(mpion, 2) + tmin))) /
              (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) + pow(mpion, 4) +
               2. * pow(mpion, 2) * pow(mrho, 2) + pow(mrho, 4) -
               2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s + pow(s, 2) +
               pow(ma1, 2) *
                   (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) -
          0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
              log(fabs(-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + tmin)) +
          (0.5 * (-2. + delta) * (eta1 - 1. * eta2) *
           (eta2 * pow(mpion, 6) * (1. * pow(mrho, 2) - 1. * s) +
            eta2 * pow(ma1, 2) * pow(mpion, 4) * (-1. * pow(mrho, 2) + 1. * s) +
            eta1 * pow(ma1, 2) * pow(mpion, 2) *
                (-0.5 * pow(mrho, 4) +
                 pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                 0.5 * pow(mrho, 2) * s) +
            eta1 * pow(mpion, 4) *
                (0.5 * pow(mrho, 4) - 0.5 * pow(mrho, 2) * s +
                 pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s))) *
           log(fabs(-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + tmin))) /
              (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
               2. * pow(ma1, 2) * pow(mpion, 2) + pow(mpion, 4)) +
          (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
           (eta1 * (pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                    (-0.5 * pow(mrho, 2) + 0.5 * s) * s) +
            eta2 *
                (-0.5 * pow(mrho, 4) + 1. * pow(mrho, 2) * s - 0.5 * pow(s, 2) +
                 pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s))) *
           log(fabs(-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + tmin))) /
              (1. * pow(ma1, 2) - 1. * pow(mpion, 2) - 1. * pow(mrho, 2) +
               1. * s) -
          (0.25 *
           (0. +
            8.000000000000002 * pow(2. - 1. * delta, 2) * pow(mpion, 4) *
                pow(mrho, 2) +
            1. * pow(2. - 1. * delta, 2) * pow(mrho, 4) * s +
            pow(mpion, 2) * (C4 * (32. - 16. * delta) * pow(mrho, 6) +
                             delta * (-8. + 4. * delta) * pow(s, 2) +
                             pow(mrho, 2) * s *
                                 (-8. + 24. * delta - 10. * pow(delta, 2) +
                                  32. * C4 * s - 16. * C4 * delta * s) +
                             pow(mrho, 4) * (-16. + 8. * delta - 64. * C4 * s +
                                             32. * C4 * delta * s))) *
           log(fabs(-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + tmin))) /
              (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
          0.03125 * pow(eta1 - 1. * eta2, 2) *
              (eta1 * eta2 *
                   (4. * pow(ma1, 6) +
                    pow(Gammaa1, 2) * pow(ma1, 2) *
                        (-4. * pow(ma1, 2) + 4. * pow(mpion, 2) - 2. * s) +
                    pow(ma1, 4) * (-12. * pow(mpion, 2) + 6. * s) +
                    pow(mpion, 2) *
                        (-4. * pow(mpion, 4) + 4. * pow(mrho, 4) +
                         2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s) +
                    pow(ma1, 2) * (12. * pow(mpion, 4) - 2. * pow(mrho, 4) -
                                   8. * pow(mpion, 2) * s -
                                   4. * pow(mrho, 2) * s + 4. * pow(s, 2))) +
               pow(eta1, 2) *
                   (-2. * pow(ma1, 6) + 2. * pow(mpion, 6) +
                    3. * pow(mpion, 4) * pow(mrho, 2) +
                    pow(ma1, 4) *
                        (6. * pow(mpion, 2) + 3. * pow(mrho, 2) - 3. * s) -
                    1. * pow(mpion, 4) * s -
                    1. * pow(mpion, 2) * pow(mrho, 2) * s -
                    1. * pow(mrho, 4) * s + pow(mrho, 2) * pow(s, 2) +
                    pow(Gammaa1, 2) * pow(ma1, 2) *
                        (2. * pow(ma1, 2) - 2. * pow(mpion, 2) -
                         1. * pow(mrho, 2) + s) +
                    pow(ma1, 2) *
                        (-6. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                         4. * pow(mrho, 2) * s - 2. * pow(s, 2) +
                         pow(mpion, 2) * (-6. * pow(mrho, 2) + 4. * s))) +
               pow(eta2, 2) *
                   (-2. * pow(ma1, 6) + 2. * pow(mpion, 6) -
                    3. * pow(mpion, 4) * pow(mrho, 2) +
                    pow(ma1, 4) *
                        (6. * pow(mpion, 2) - 3. * pow(mrho, 2) - 3. * s) -
                    1. * pow(mpion, 4) * s + pow(mpion, 2) * pow(mrho, 2) * s +
                    pow(Gammaa1, 2) * pow(ma1, 2) *
                        (2. * pow(ma1, 2) - 2. * pow(mpion, 2) + pow(mrho, 2) +
                         s) +
                    pow(ma1, 2) *
                        (-6. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                         2. * pow(mrho, 2) * s - 2. * pow(s, 2) +
                         pow(mpion, 2) * (6. * pow(mrho, 2) + 4. * s)))) *
              log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                       4. * pow(ma1, 2) * pow(mpion, 2) + 4. * pow(mpion, 4) +
                       2. * pow(ma1, 2) * (-1. * pow(mrho, 2) + s + tmin) -
                       4. * pow(mpion, 2) * (-1. * pow(mrho, 2) + s + tmin) +
                       pow(-1. * pow(mrho, 2) + s + tmin, 2))) -
          (0.5 * (1. * eta1 - 1. * eta2) *
           (eta2 *
                (pow(Gammaa1, 2) * pow(ma1, 2) * pow(mrho, 2) *
                     (0.5 - 1. * C4 * pow(mrho, 2)) +
                 pow(ma1, 4) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) +
                 pow(mpion, 2) * pow(mrho, 2) *
                     (pow(mpion, 2) * (-0.5 + 1. * C4 * pow(mrho, 2)) +
                      (0.25 - 0.125 * delta) * (pow(mrho, 2) + s)) +
                 pow(ma1, 2) *
                     (1. * C4 * pow(mrho, 6) +
                      pow(mpion, 2) *
                          (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) -
                      0.25 * delta * pow(s, 2) +
                      pow(mrho, 4) * (-0.75 + 0.125 * delta - 2. * C4 * s) +
                      pow(mrho, 2) * s *
                          (0.25 + 0.375 * delta + 1. * C4 * s))) +
            eta1 *
                (pow(Gammaa1, 2) * pow(ma1, 2) * pow(mrho, 2) *
                     (-0.5 + 1. * C4 * pow(mrho, 2)) +
                 pow(ma1, 4) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                 pow(ma1, 2) * (-0.5 * pow(mrho, 4) + 1. * C4 * pow(mrho, 6) +
                                pow(mpion, 2) * (-1. * pow(mrho, 2) +
                                                 2. * C4 * pow(mrho, 4)) +
                                0.25 * delta * pow(s, 2) -
                                1. * C4 * pow(mrho, 2) * pow(s, 2)) +
                 pow(mrho, 2) *
                     (pow(mpion, 4) * (0.5 - 1. * C4 * pow(mrho, 2)) +
                      s * ((-0.25 + 0.125 * delta) * pow(mrho, 2) +
                           (0.25 - 0.125 * delta) * s) +
                      pow(mpion, 2) *
                          (-2. * C4 * pow(mrho, 4) + (-0.5 - 0.25 * delta) * s +
                           pow(mrho, 2) * (1. + 2. * C4 * s))))) *
           log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) +
                    pow(pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) +
                            s + tmin,
                        2)))) /
              pow(mrho, 2) +
          (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
           (eta2 *
                (0.5 * pow(Gammaa1, 4) * pow(ma1, 4) - 0.5 * pow(ma1, 8) +
                 0.5 * pow(mpion, 8) +
                 0.5 * pow(ma1, 4) * pow(mpion, 2) * pow(mrho, 2) -
                 0.5 * pow(mpion, 6) * pow(mrho, 2) +
                 pow(Gammaa1, 2) *
                     (pow(ma1, 2) * pow(mpion, 2) *
                          (1. * pow(mpion, 2) + 1.5 * pow(mrho, 2) - 2. * s) +
                      pow(ma1, 4) *
                          (-1. * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 1. * s)) +
                 pow(ma1, 6) *
                     (1. * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 1. * s) +
                 pow(ma1, 2) * pow(mpion, 4) *
                     (-1. * pow(mpion, 2) - 0.5 * pow(mrho, 2) + 1. * s)) +
            eta1 *
                (pow(ma1, 6) * (1. * pow(mpion, 2) + 0.5 * s) +
                 pow(ma1, 2) *
                     (3. * pow(mpion, 6) + 1. * pow(mpion, 2) * pow(mrho, 4) -
                      0.5 * pow(mpion, 4) * s) +
                 pow(ma1, 4) * (-3. * pow(mpion, 4) +
                                pow(mpion, 2) * (-1. * pow(mrho, 2) + 0.5 * s) -
                                0.5 * pow(mrho, 2) * s) +
                 pow(mpion, 4) *
                     (-1. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                      pow(mpion, 2) * (1. * pow(mrho, 2) - 0.5 * s) +
                      0.5 * pow(mrho, 2) * s) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (-1. * pow(mpion, 4) +
                      pow(ma1, 2) * (1. * pow(mpion, 2) + 0.5 * s) -
                      0.5 * pow(mrho, 2) * s +
                      pow(mpion, 2) * (-1. * pow(mrho, 2) + 1.5 * s)))) *
           log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) +
                    4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                    pow(mrho, 4) - 4. * pow(mpion, 2) * s -
                    2. * pow(mrho, 2) * s + pow(s, 2) -
                    4. * pow(mpion, 2) * tmin - 2. * pow(mrho, 2) * tmin +
                    2. * s * tmin + pow(tmin, 2) +
                    pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                   2. * s + 2. * tmin)))) /
              (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
               2. * pow(ma1, 2) * pow(mpion, 2) + pow(mpion, 4)) -
          (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
           (eta1 *
                (-1. * pow(mpion, 8) + 1. * pow(mpion, 4) * pow(mrho, 4) +
                 pow(ma1, 6) * (1. * pow(mpion, 2) - 0.5 * s) -
                 0.5 * pow(mpion, 6) * s -
                 4. * pow(mpion, 4) * pow(mrho, 2) * s -
                 1.5 * pow(mpion, 2) * pow(mrho, 4) * s +
                 3.5 * pow(mpion, 4) * pow(s, 2) +
                 4. * pow(mpion, 2) * pow(mrho, 2) * pow(s, 2) +
                 0.5 * pow(mrho, 4) * pow(s, 2) -
                 2.5 * pow(mpion, 2) * pow(s, 3) -
                 1. * pow(mrho, 2) * pow(s, 3) + 0.5 * pow(s, 4) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (-1. * pow(mpion, 4) +
                      pow(ma1, 2) * (1. * pow(mpion, 2) - 0.5 * s) -
                      0.5 * pow(mpion, 2) * s + 0.5 * pow(s, 2)) +
                 pow(ma1, 4) *
                     (-3. * pow(mpion, 4) + (1. * pow(mrho, 2) - 0.5 * s) * s +
                      pow(mpion, 2) * (-2. * pow(mrho, 2) + 2.5 * s)) +
                 pow(ma1, 2) * (3. * pow(mpion, 6) +
                                pow(mpion, 4) * (2. * pow(mrho, 2) - 1.5 * s) -
                                0.5 * pow(mrho, 4) * s + 0.5 * pow(s, 3) +
                                pow(mpion, 2) *
                                    (1. * pow(mrho, 4) - 1. * pow(mrho, 2) * s -
                                     1. * pow(s, 2)))) +
            eta2 *
                (0.5 * pow(Gammaa1, 4) * pow(ma1, 4) - 0.5 * pow(ma1, 8) +
                 0.5 * pow(mpion, 8) +
                 pow(ma1, 6) *
                     (1. * pow(mpion, 2) + 1. * pow(mrho, 2) - 0.5 * s) +
                 0.5 * pow(mpion, 6) * s +
                 pow(ma1, 4) * (-0.5 * pow(mrho, 4) +
                                (-0.5 * pow(mpion, 2) + 0.5 * s) * s) +
                 pow(mpion, 4) * (-0.5 * pow(mrho, 4) + 2. * pow(mrho, 2) * s -
                                  1.5 * pow(s, 2)) +
                 pow(mpion, 2) * s *
                     (0.5 * pow(mrho, 4) - 1. * pow(mrho, 2) * s +
                      0.5 * pow(s, 2)) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (1. * pow(mpion, 4) + 0.5 * pow(mrho, 4) +
                      pow(mpion, 2) * (2. * pow(mrho, 2) - 1.5 * s) -
                      1. * pow(mrho, 2) * s + 0.5 * pow(s, 2) +
                      pow(ma1, 2) *
                          (-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + 1.5 * s)) +
                 pow(ma1, 2) *
                     (-1. * pow(mpion, 6) +
                      pow(mpion, 4) * (-1. * pow(mrho, 2) + 0.5 * s) +
                      pow(mpion, 2) * (-1. * pow(mrho, 4) +
                                       2. * pow(mrho, 2) * s - 1. * pow(s, 2)) +
                      s * (0.5 * pow(mrho, 4) - 1. * pow(mrho, 2) * s +
                           0.5 * pow(s, 2))))) *
           log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) +
                    4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                    pow(mrho, 4) - 4. * pow(mpion, 2) * s -
                    2. * pow(mrho, 2) * s + pow(s, 2) -
                    4. * pow(mpion, 2) * tmin - 2. * pow(mrho, 2) * tmin +
                    2. * s * tmin + pow(tmin, 2) +
                    pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                   2. * s + 2. * tmin)))) /
              (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) + pow(mpion, 4) +
               2. * pow(mpion, 2) * pow(mrho, 2) + pow(mrho, 4) -
               2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s + pow(s, 2) +
               pow(ma1, 2) *
                   (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) -
          (0.0625 * pow(eta1 - 1. * eta2, 2) *
           (pow(eta2, 2) *
                (-1. * pow(ma1, 10) +
                 pow(ma1, 8) *
                     (5. * pow(mpion, 2) + 2.5 * pow(mrho, 2) - 2.5 * s) +
                 pow(Gammaa1, 4) * pow(ma1, 4) *
                     (1. * pow(ma1, 2) - 1. * pow(mpion, 2) -
                      0.5 * pow(mrho, 2) + 0.5 * s) +
                 pow(ma1, 4) *
                     (10. * pow(mpion, 6) - 2.5 * pow(mrho, 6) +
                      pow(mpion, 4) * (15. * pow(mrho, 2) - 9. * s) +
                      6. * pow(mrho, 4) * s - 4.5 * pow(mrho, 2) * pow(s, 2) +
                      1. * pow(s, 3)) +
                 pow(ma1, 6) *
                     (-10. * pow(mpion, 4) + (1. * pow(mrho, 2) - 1. * s) * s +
                      pow(mpion, 2) * (-10. * pow(mrho, 2) + 8. * s)) +
                 pow(mpion, 4) *
                     (1. * pow(mpion, 6) + 0.5 * pow(mrho, 6) +
                      pow(mpion, 4) * (2.5 * pow(mrho, 2) - 0.5 * s) -
                      1. * pow(mrho, 4) * s + 0.5 * pow(mrho, 2) * pow(s, 2) +
                      pow(mpion, 2) *
                          (2. * pow(mrho, 4) - 2. * pow(mrho, 2) * s)) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (4. * pow(ma1, 6) - 4. * pow(mpion, 6) -
                      0.5 * pow(mrho, 6) + 1.5 * pow(mrho, 4) * s -
                      1.5 * pow(mrho, 2) * pow(s, 2) + 0.5 * pow(s, 3) +
                      pow(mpion, 4) * (-6. * pow(mrho, 2) + 6. * s) +
                      pow(ma1, 4) *
                          (-12. * pow(mpion, 2) - 6. * pow(mrho, 2) + 6. * s) +
                      pow(mpion, 2) * (-3. * pow(mrho, 4) +
                                       6. * pow(mrho, 2) * s - 3. * pow(s, 2)) +
                      pow(ma1, 2) *
                          (12. * pow(mpion, 4) + 3. * pow(mrho, 4) +
                           pow(mpion, 2) * (12. * pow(mrho, 2) - 12. * s) -
                           6. * pow(mrho, 2) * s + 3. * pow(s, 2))) +
                 pow(ma1, 2) *
                     (-5. * pow(mpion, 8) + 1. * pow(mrho, 8) -
                      3.5 * pow(mrho, 6) * s + 4.5 * pow(mrho, 4) * pow(s, 2) -
                      2.5 * pow(mrho, 2) * pow(s, 3) + 0.5 * pow(s, 4) +
                      pow(mpion, 6) * (-10. * pow(mrho, 2) + 4. * s) +
                      pow(mpion, 4) * (-2. * pow(mrho, 4) +
                                       1. * pow(mrho, 2) * s + 1. * pow(s, 2)) +
                      pow(mpion, 2) *
                          (3. * pow(mrho, 6) - 8. * pow(mrho, 4) * s +
                           7. * pow(mrho, 2) * pow(s, 2) - 2. * pow(s, 3)))) +
            pow(eta1, 2) *
                (-1. * pow(ma1, 10) +
                 pow(ma1, 8) *
                     (5. * pow(mpion, 2) + 2.5 * pow(mrho, 2) - 2.5 * s) +
                 pow(Gammaa1, 4) * pow(ma1, 4) *
                     (1. * pow(ma1, 2) - 1. * pow(mpion, 2) -
                      0.5 * pow(mrho, 2) + 0.5 * s) +
                 pow(ma1, 6) *
                     (-10. * pow(mpion, 4) - 2. * pow(mrho, 4) +
                      5. * pow(mrho, 2) * s - 1. * pow(s, 2) +
                      pow(mpion, 2) * (-10. * pow(mrho, 2) + 8. * s)) +
                 pow(ma1, 4) *
                     (10. * pow(mpion, 6) + 0.5 * pow(mrho, 6) +
                      pow(mpion, 4) * (15. * pow(mrho, 2) - 9. * s) -
                      3. * pow(mrho, 4) * s + 1.5 * pow(mrho, 2) * pow(s, 2) +
                      1. * pow(s, 3) +
                      pow(mpion, 2) *
                          (6. * pow(mrho, 4) - 12. * pow(mrho, 2) * s)) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (4. * pow(ma1, 6) - 4. * pow(mpion, 6) -
                      0.5 * pow(mrho, 6) + 1.5 * pow(mrho, 4) * s -
                      1.5 * pow(mrho, 2) * pow(s, 2) + 0.5 * pow(s, 3) +
                      pow(mpion, 4) * (-6. * pow(mrho, 2) + 6. * s) +
                      pow(ma1, 4) *
                          (-12. * pow(mpion, 2) - 6. * pow(mrho, 2) + 6. * s) +
                      pow(mpion, 2) * (-3. * pow(mrho, 4) +
                                       6. * pow(mrho, 2) * s - 3. * pow(s, 2)) +
                      pow(ma1, 2) *
                          (12. * pow(mpion, 4) + 3. * pow(mrho, 4) +
                           pow(mpion, 2) * (12. * pow(mrho, 2) - 12. * s) -
                           6. * pow(mrho, 2) * s + 3. * pow(s, 2))) +
                 pow(mpion, 2) *
                     (1. * pow(mpion, 8) +
                      pow(mpion, 6) * (2.5 * pow(mrho, 2) - 0.5 * s) +
                      pow(mpion, 4) *
                          (4. * pow(mrho, 4) - 6. * pow(mrho, 2) * s) +
                      pow(mrho, 2) * s *
                          (-1. * pow(mrho, 4) + 2. * pow(mrho, 2) * s -
                           1. * pow(s, 2)) +
                      pow(mpion, 2) *
                          (1.5 * pow(mrho, 6) - 6. * pow(mrho, 4) * s +
                           4.5 * pow(mrho, 2) * pow(s, 2))) +
                 pow(ma1, 2) *
                     (-5. * pow(mpion, 8) +
                      pow(mpion, 6) * (-10. * pow(mrho, 2) + 4. * s) +
                      pow(mpion, 4) *
                          (-8. * pow(mrho, 4) + 13. * pow(mrho, 2) * s +
                           1. * pow(s, 2)) +
                      pow(mpion, 2) *
                          (-1. * pow(mrho, 6) + 6. * pow(mrho, 4) * s -
                           3. * pow(mrho, 2) * pow(s, 2) - 2. * pow(s, 3)) +
                      s * (0.5 * pow(mrho, 6) - 0.5 * pow(mrho, 4) * s -
                           0.5 * pow(mrho, 2) * pow(s, 2) + 0.5 * pow(s, 3)))) +
            eta1 * eta2 *
                (2. * pow(ma1, 10) +
                 pow(Gammaa1, 4) * pow(ma1, 4) *
                     (-2. * pow(ma1, 2) + 2. * pow(mpion, 2) +
                      1. * pow(mrho, 2) - 1. * s) +
                 pow(ma1, 8) *
                     (-10. * pow(mpion, 2) - 5. * pow(mrho, 2) + 5. * s) +
                 pow(ma1, 6) * (20. * pow(mpion, 4) + 6. * pow(mrho, 4) +
                                pow(mpion, 2) * (20. * pow(mrho, 2) - 16. * s) -
                                8. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                 pow(ma1, 4) *
                     (-20. * pow(mpion, 6) - 4. * pow(mrho, 6) +
                      6. * pow(mrho, 4) * s - 2. * pow(s, 3) +
                      pow(mpion, 4) * (-30. * pow(mrho, 2) + 18. * s) +
                      pow(mpion, 2) *
                          (-18. * pow(mrho, 4) + 18. * pow(mrho, 2) * s)) +
                 pow(mpion, 2) *
                     (-2. * pow(mpion, 8) - 1. * pow(mrho, 8) +
                      3. * pow(mrho, 6) * s - 3. * pow(mrho, 4) * pow(s, 2) +
                      1. * pow(mrho, 2) * pow(s, 3) +
                      pow(mpion, 6) * (-5. * pow(mrho, 2) + 1. * s) +
                      pow(mpion, 4) *
                          (-10. * pow(mrho, 4) + 10. * pow(mrho, 2) * s) +
                      pow(mpion, 2) *
                          (-6. * pow(mrho, 6) + 12. * pow(mrho, 4) * s -
                           6. * pow(mrho, 2) * pow(s, 2))) +
                 pow(ma1, 2) * (10. * pow(mpion, 8) + 1. * pow(mrho, 8) +
                                pow(mpion, 6) * (20. * pow(mrho, 2) - 8. * s) -
                                2. * pow(mrho, 6) * s +
                                2. * pow(mrho, 2) * pow(s, 3) - 1. * pow(s, 4) +
                                pow(mpion, 4) *
                                    (22. * pow(mrho, 4) -
                                     20. * pow(mrho, 2) * s - 2. * pow(s, 2)) +
                                pow(mpion, 2) *
                                    (8. * pow(mrho, 6) -
                                     12. * pow(mrho, 4) * s + 4. * pow(s, 3))) +
                 pow(Gammaa1, 2) * pow(ma1, 2) *
                     (-8. * pow(ma1, 6) + 8. * pow(mpion, 6) +
                      1. * pow(mrho, 6) +
                      pow(mpion, 4) * (12. * pow(mrho, 2) - 12. * s) +
                      pow(ma1, 4) *
                          (24. * pow(mpion, 2) + 12. * pow(mrho, 2) - 12. * s) -
                      3. * pow(mrho, 4) * s + 3. * pow(mrho, 2) * pow(s, 2) -
                      1. * pow(s, 3) +
                      pow(mpion, 2) *
                          (6. * pow(mrho, 4) - 12. * pow(mrho, 2) * s +
                           6. * pow(s, 2)) +
                      pow(ma1, 2) *
                          (-24. * pow(mpion, 4) - 6. * pow(mrho, 4) +
                           12. * pow(mrho, 2) * s - 6. * pow(s, 2) +
                           pow(mpion, 2) * (-24. * pow(mrho, 2) + 24. * s))))) *
           log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) +
                    4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                    pow(mrho, 4) - 4. * pow(mpion, 2) * s -
                    2. * pow(mrho, 2) * s + pow(s, 2) -
                    4. * pow(mpion, 2) * tmin - 2. * pow(mrho, 2) * tmin +
                    2. * s * tmin + pow(tmin, 2) +
                    pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                   2. * s + 2. * tmin)))) /
              (pow(Gammaa1, 2) * pow(ma1, 2) + 4. * pow(ma1, 4) +
               4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
               pow(mrho, 4) - 4. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s +
               pow(s, 2) +
               pow(ma1, 2) *
                   (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)))) /
           (16. * Pi * s * (-4 * pow(mpion, 2) + s)) +
       (pow(Const, 2) * pow(ghat, 4) *
        (0. +
         (0.03125 * pow(eta1 - 1. * eta2, 2) *
          (eta1 * eta2 *
               (-2. * pow(ma1, 8) - 2. * pow(mpion, 8) +
                2. * pow(mpion, 4) * pow(mrho, 4) +
                pow(ma1, 6) * (8. * pow(mpion, 2) - 4. * s) +
                pow(ma1, 2) * pow(mpion, 2) *
                    (8. * pow(mpion, 4) - 8. * pow(mrho, 4) -
                     4. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s) +
                pow(ma1, 4) * (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                               8. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s -
                               4. * pow(s, 2))) +
           pow(eta2, 2) *
               (1. * pow(ma1, 8) + 1. * pow(mpion, 8) -
                2. * pow(mpion, 6) * pow(mrho, 2) +
                1. * pow(mpion, 4) * pow(mrho, 4) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) + 2. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                               pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                               2. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                pow(ma1, 2) * (-4. * pow(mpion, 6) -
                               2. * pow(mpion, 2) * pow(mrho, 2) * s +
                               pow(mpion, 4) * (6. * pow(mrho, 2) + 2. * s))) +
           pow(eta1, 2) *
               (1. * pow(ma1, 8) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                               pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                               4. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                pow(ma1, 2) * (-4. * pow(mpion, 6) +
                               2. * pow(mpion, 2) * pow(mrho, 2) * s +
                               pow(mrho, 2) * (2. * pow(mrho, 2) - 2. * s) * s +
                               pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s)) +
                pow(mpion, 2) *
                    (1. * pow(mpion, 6) + 2. * pow(mpion, 4) * pow(mrho, 2) -
                     2. * pow(mrho, 6) + 2. * pow(mrho, 4) * s +
                     pow(mpion, 2) *
                         (1. * pow(mrho, 4) - 2. * pow(mrho, 2) * s))))) /
             (1. * pow(ma1, 2) - 1. * tmax) +
         (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
          (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
             (1. * pow(mpion, 2) - 1. * tmax) +
         (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
          (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
             (1. * pow(mpion, 2) + 1. * pow(mrho, 2) - 1. * s - 1. * tmax) -
         (0.5 * pow(-2. + delta, 2) * pow(mpion, 2) * tmax) / pow(mrho, 2) -
         0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
             (-0.5 * eta2 * pow(ma1, 2) + 1. * eta1 * pow(mpion, 2) +
              0.5 * eta2 * pow(mrho, 2) + 0.5 * eta1 * s - 1. * eta2 * s) *
             tmax +
         (0.25 *
          (pow(mpion, 2) * (12. + 1. * pow(delta, 2) - 16. * C4 * pow(mrho, 2) +
                            delta * (-8. + 8. * C4 * pow(mrho, 2))) +
           (-4. - 3. * pow(delta, 2) - 16. * C4 * pow(mrho, 2) +
            delta * (8. + 8. * C4 * pow(mrho, 2))) *
               s) *
          tmax) /
             pow(mrho, 2) -
         0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta2 * (pow(ma1, 2) - 1. * s) +
              eta1 * (-2. * pow(mpion, 2) + s)) *
             tmax -
         0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta2 * (-1. * pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * s) +
              eta1 * (2. * pow(mpion, 2) + s)) *
             tmax +
         0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta1 * (1. * pow(mpion, 2) - 0.5 * s) +
              eta2 * (-0.5 * pow(ma1, 2) - 1. * pow(mpion, 2) -
                      0.5 * pow(mrho, 2) + 1. * s)) *
             tmax +
         (0.25 * (-2. + 1. * delta) *
          (-8. * C4 * pow(mrho, 4) +
           pow(mpion, 2) * (2. + 1. * delta - 8. * C4 * pow(mrho, 2)) +
           (-2. - 3. * delta) * s +
           pow(mrho, 2) * (2. + 1. * delta + 16. * C4 * s)) *
          tmax) /
             pow(mrho, 2) +
         (0.25 *
          (32 * pow(C4, 2) * pow(mrho, 8) + 2 * pow(delta, 2) * pow(s, 2) +
           8 * C4 * pow(mrho, 6) * (-6 + delta - 8 * C4 * s) +
           2 * delta * pow(mrho, 2) * s * (-6 + delta - 8 * C4 * s) +
           pow(mrho, 4) * (12 - pow(delta, 2) + 8 * C4 * (6 + delta) * s +
                           32 * pow(C4, 2) * pow(s, 2))) *
          tmax) /
             pow(mrho, 4) -
         (1. * (1. * eta1 - 1. * eta2) *
          (eta2 *
               (0.75 * pow(mrho, 4) - 0.125 * delta * pow(mrho, 4) -
                1. * C4 * pow(mrho, 6) +
                pow(ma1, 2) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                pow(mpion, 2) * (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4)) -
                0.25 * pow(mrho, 2) * s - 0.375 * delta * pow(mrho, 2) * s +
                2. * C4 * pow(mrho, 4) * s + 0.25 * delta * pow(s, 2) -
                1. * C4 * pow(mrho, 2) * pow(s, 2)) +
           eta1 *
               (0.5 * pow(mrho, 4) - 1. * C4 * pow(mrho, 6) +
                pow(mpion, 2) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                pow(ma1, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                0.25 * delta * pow(s, 2) +
                1. * C4 * pow(mrho, 2) * pow(s, 2))) *
          tmax) /
             pow(mrho, 2) +
         0.0625 * pow(eta1 - 1. * eta2, 2) *
             (pow(eta2, 2) *
                  (pow(Gammaa1, 2) * pow(ma1, 2) - 1. * pow(ma1, 4) -
                   2. * pow(mpion, 4) - 2. * pow(mpion, 2) * pow(mrho, 2) +
                   2. * pow(mrho, 4) +
                   pow(ma1, 2) * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s) -
                   3. * pow(mrho, 2) * s + pow(s, 2)) +
              pow(eta1, 2) *
                  (pow(Gammaa1, 2) * pow(ma1, 2) - 1. * pow(ma1, 4) -
                   2. * pow(mpion, 4) - 2. * pow(mpion, 2) * pow(mrho, 2) +
                   pow(ma1, 2) * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s) +
                   pow(mrho, 2) * s + pow(s, 2)) +
              eta1 * eta2 *
                  (-2. * pow(Gammaa1, 2) * pow(ma1, 2) + 2. * pow(ma1, 4) +
                   4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                   2. * pow(mrho, 4) - 2. * pow(s, 2) +
                   pow(ma1, 2) *
                       (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s))) *
             tmax +
         0.03125 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-6. * pow(ma1, 4) - 12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                   pow(ma1, 2) * (16. * pow(mpion, 2) - 8. * s) +
                   8. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s -
                   4. * pow(s, 2)) +
              pow(eta1, 2) *
                  (3. * pow(ma1, 4) + 6. * pow(mpion, 4) + pow(mrho, 4) +
                   pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                   4. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                   pow(ma1, 2) *
                       (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) +
              pow(eta2, 2) *
                  (3. * pow(ma1, 4) + 6. * pow(mpion, 4) + pow(mrho, 4) +
                   pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                   2. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                   pow(ma1, 2) *
                       (-8. * pow(mpion, 2) + 4. * pow(mrho, 2) + 4. * s))) *
             tmax -
         (0.125 * (-2. + 1. * delta) *
          (2. + 1. * delta - 8. * C4 * pow(mrho, 2)) * pow(tmax, 2)) /
             pow(mrho, 2) -
         0.5 * pow(1. * eta1 - 1. * eta2, 2) * (-0.5 + 1. * C4 * pow(mrho, 2)) *
             pow(tmax, 2) -
         (1. *
          (0.5 - 0.125 * pow(delta, 2) - 2. * C4 * pow(mrho, 2) +
           1. * C4 * delta * pow(mrho, 2)) *
          pow(tmax, 2)) /
             pow(mrho, 2) +
         0.0625 * pow(1. * eta1 - 1. * eta2, 4) *
             (1. * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 0.5 * s) *
             pow(tmax, 2) +
         0.03125 * pow(eta1 - 1. * eta2, 3) *
             (eta2 * (-1. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                      1. * pow(mrho, 2) - 1. * s) +
              eta1 *
                  (pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
             pow(tmax, 2) +
         0.010416666666666666 * pow(eta1 - 1. * eta2, 4) * pow(tmax, 3) -
         0.020833333333333332 * pow(1. * eta1 - 1. * eta2, 4) * pow(tmax, 3) +
         0.03125 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (2. * pow(Gammaa1, 2) * pow(ma1, 2) - 6. * pow(ma1, 4) -
                   4. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                   pow(ma1, 2) * (8. * pow(mpion, 2) - 8. * s) +
                   4. * pow(mrho, 2) * s - 4. * pow(s, 2)) +
              pow(eta1, 2) *
                  (-1. * pow(Gammaa1, 2) * pow(ma1, 2) + 3. * pow(ma1, 4) +
                   2. * pow(mpion, 4) + 2. * pow(mpion, 2) * pow(mrho, 2) +
                   pow(mrho, 4) - 4. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                   pow(ma1, 2) *
                       (-4. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) +
              pow(eta2, 2) *
                  (-1. * pow(Gammaa1, 2) * pow(ma1, 2) + 3. * pow(ma1, 4) +
                   2. * pow(mpion, 4) - 2. * pow(mpion, 2) * pow(mrho, 2) +
                   pow(mrho, 4) - 2. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                   pow(ma1, 2) *
                       (-4. * pow(mpion, 2) + 4. * pow(mrho, 2) + 4. * s))) *
             (-1. * pow(mrho, 2) + s + tmax) -
         0.03125 * pow(eta1 - 1. * eta2, 3) *
             (eta2 * (-1. * pow(ma1, 2) - 1. * pow(mrho, 2) - 1. * s) +
              eta1 * (pow(ma1, 2) - 1. * pow(mrho, 2) + s)) *
             pow(-1. * pow(mrho, 2) + s + tmax, 2) +
         0.010416666666666666 * pow(eta1 - 1. * eta2, 4) *
             pow(-1. * pow(mrho, 2) + s + tmax, 3) +
         0.25 * (eta1 - 1. * eta2) * (1. * eta1 - 1. * eta2) *
             (-1. + 2. * C4 * pow(mrho, 2)) *
             pow(pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                     tmax,
                 2) -
         (2. * (1. * eta1 - 1. * eta2) *
          (eta2 *
               (0.375 * pow(mrho, 4) - 0.0625 * delta * pow(mrho, 4) -
                0.5 * C4 * pow(mrho, 6) +
                pow(ma1, 2) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                pow(mpion, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                0.125 * pow(mrho, 2) * s - 0.1875 * delta * pow(mrho, 2) * s +
                1. * C4 * pow(mrho, 4) * s + 0.125 * delta * pow(s, 2) -
                0.5 * C4 * pow(mrho, 2) * pow(s, 2)) +
           eta1 *
               (0.25 * pow(mrho, 4) - 0.5 * C4 * pow(mrho, 6) +
                pow(mpion, 2) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                pow(ma1, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                0.125 * delta * pow(s, 2) +
                0.5 * C4 * pow(mrho, 2) * pow(s, 2))) *
          (1. * pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + 1. * s +
           1. * tmax)) /
             pow(mrho, 2) +
         (2. * (1. * eta1 - 1. * eta2) * Gammaa1 * ma1 *
          (eta2 *
               (0.375 * pow(mrho, 4) - 0.0625 * delta * pow(mrho, 4) -
                0.5 * C4 * pow(mrho, 6) +
                pow(ma1, 2) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                pow(mpion, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                0.125 * pow(mrho, 2) * s - 0.1875 * delta * pow(mrho, 2) * s +
                1. * C4 * pow(mrho, 4) * s + 0.125 * delta * pow(s, 2) -
                0.5 * C4 * pow(mrho, 2) * pow(s, 2)) +
           eta1 *
               (0.25 * pow(mrho, 4) - 0.5 * C4 * pow(mrho, 6) +
                pow(mpion, 2) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                pow(ma1, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                0.125 * delta * pow(s, 2) +
                0.5 * C4 * pow(mrho, 2) * pow(s, 2))) *
          atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                tmax) /
               (Gammaa1 * ma1))) /
             pow(mrho, 2) +
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) * Gammaa1 * ma1 *
          (eta2 * (-1. * pow(ma1, 6) +
                   pow(Gammaa1, 2) * pow(ma1, 2) *
                       (-1. * pow(ma1, 2) + 0.5 * pow(mrho, 2) - 1. * s) +
                   pow(ma1, 4) *
                       (2. * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 1. * s) +
                   pow(mpion, 4) * (-1.5 * pow(mrho, 2) + 1. * s) +
                   pow(ma1, 2) * pow(mpion, 2) *
                       (-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + 2. * s)) +
           eta1 *
               (pow(Gammaa1, 2) * pow(ma1, 2) * (1. * pow(mpion, 2) + 0.5 * s) +
                pow(ma1, 4) * (1. * pow(mpion, 2) + 0.5 * s) +
                pow(ma1, 2) * (-2. * pow(mpion, 4) - 1. * pow(mpion, 2) * s) +
                pow(mpion, 2) * (1. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                                 pow(mpion, 2) * (2. * pow(mrho, 2) - 1.5 * s) +
                                 1. * pow(mrho, 2) * s))) *
          atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                tmax) /
               (Gammaa1 * ma1))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
              2. * pow(ma1, 2) * pow(mpion, 2) + pow(mpion, 4)) -
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) * Gammaa1 * ma1 *
          (eta2 *
               (-1. * pow(ma1, 6) - 2. * pow(mpion, 4) * pow(mrho, 2) -
                1. * pow(mpion, 2) * pow(mrho, 4) +
                pow(ma1, 4) *
                    (2. * pow(mpion, 2) + 2. * pow(mrho, 2) - 1.5 * s) +
                2.5 * pow(mpion, 4) * s +
                3. * pow(mpion, 2) * pow(mrho, 2) * s + 0.5 * pow(mrho, 4) * s -
                2. * pow(mpion, 2) * pow(s, 2) - 1. * pow(mrho, 2) * pow(s, 2) +
                0.5 * pow(s, 3) +
                pow(Gammaa1, 2) * (-1. * pow(ma1, 4) + 0.5 * pow(ma1, 2) * s) +
                pow(ma1, 2) * (-1. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                               1. * pow(mrho, 2) * s +
                               pow(mpion, 2) * (-2. * pow(mrho, 2) + 1. * s))) +
           eta1 *
               (1. * pow(mpion, 6) + 4. * pow(mpion, 4) * pow(mrho, 2) +
                1. * pow(mpion, 2) * pow(mrho, 4) +
                pow(Gammaa1, 2) * pow(ma1, 2) * (1. * pow(mpion, 2) - 0.5 * s) +
                pow(ma1, 4) * (1. * pow(mpion, 2) - 0.5 * s) -
                4.5 * pow(mpion, 4) * s -
                4. * pow(mpion, 2) * pow(mrho, 2) * s - 0.5 * pow(mrho, 4) * s +
                3. * pow(mpion, 2) * pow(s, 2) + 1. * pow(mrho, 2) * pow(s, 2) -
                0.5 * pow(s, 3) +
                pow(ma1, 2) *
                    (-2. * pow(mpion, 4) + (1. * pow(mrho, 2) - 1. * s) * s +
                     pow(mpion, 2) * (-2. * pow(mrho, 2) + 3. * s)))) *
          atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                tmax) /
               (Gammaa1 * ma1))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) + pow(mpion, 4) +
              2. * pow(mpion, 2) * pow(mrho, 2) + pow(mrho, 4) -
              2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s + pow(s, 2) +
              pow(ma1, 2) *
                  (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) +
         (0.03125 * pow(eta1 - 1. * eta2, 2) *
          (pow(eta2, 2) *
               (pow(Gammaa1, 4) * pow(ma1, 4) + pow(ma1, 8) + pow(mpion, 8) -
                2. * pow(mpion, 6) * pow(mrho, 2) +
                pow(mpion, 4) * pow(mrho, 4) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) + 2. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) + pow(mrho, 4) +
                               pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                               2. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                pow(ma1, 2) * (-4. * pow(mpion, 6) -
                               2. * pow(mpion, 2) * pow(mrho, 2) * s +
                               pow(mpion, 4) * (6. * pow(mrho, 2) + 2. * s)) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (-6. * pow(ma1, 4) - 6. * pow(mpion, 4) -
                     1. * pow(mrho, 4) +
                     pow(ma1, 2) *
                         (12. * pow(mpion, 2) - 6. * pow(mrho, 2) - 6. * s) +
                     2. * pow(mrho, 2) * s - 2. * pow(s, 2) +
                     pow(mpion, 2) * (6. * pow(mrho, 2) + 4. * s))) +
           eta1 * eta2 *
               (-2. * pow(Gammaa1, 4) * pow(ma1, 4) - 2. * pow(ma1, 8) -
                2. * pow(mpion, 8) + 2. * pow(mpion, 4) * pow(mrho, 4) +
                pow(ma1, 6) * (8. * pow(mpion, 2) - 4. * s) +
                pow(ma1, 2) * pow(mpion, 2) *
                    (8. * pow(mpion, 4) - 8. * pow(mrho, 4) -
                     4. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s) +
                pow(ma1, 4) * (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                               8. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s -
                               4. * pow(s, 2)) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (12. * pow(ma1, 4) + 12. * pow(mpion, 4) -
                     2. * pow(mrho, 4) - 8. * pow(mpion, 2) * s -
                     4. * pow(mrho, 2) * s + 4. * pow(s, 2) +
                     pow(ma1, 2) * (-24. * pow(mpion, 2) + 12. * s))) +
           pow(eta1, 2) *
               (pow(Gammaa1, 4) * pow(ma1, 4) + pow(ma1, 8) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) + pow(mrho, 4) +
                               pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                               4. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                pow(ma1, 2) * (-4. * pow(mpion, 6) +
                               2. * pow(mpion, 2) * pow(mrho, 2) * s +
                               pow(mrho, 2) * (2. * pow(mrho, 2) - 2. * s) * s +
                               pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s)) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (-6. * pow(ma1, 4) - 6. * pow(mpion, 4) -
                     1. * pow(mrho, 4) +
                     pow(ma1, 2) *
                         (12. * pow(mpion, 2) + 6. * pow(mrho, 2) - 6. * s) +
                     4. * pow(mrho, 2) * s - 2. * pow(s, 2) +
                     pow(mpion, 2) * (-6. * pow(mrho, 2) + 4. * s)) +
                pow(mpion, 2) *
                    (pow(mpion, 6) + 2. * pow(mpion, 4) * pow(mrho, 2) -
                     2. * pow(mrho, 6) + 2. * pow(mrho, 4) * s +
                     pow(mpion, 2) * (pow(mrho, 4) - 2. * pow(mrho, 2) * s)))) *
          atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                tmax) /
               (Gammaa1 * ma1))) /
             (Gammaa1 * ma1) -
         (0.0625 * pow(eta1 - 1. * eta2, 2) * Gammaa1 * ma1 *
          (eta1 * eta2 *
               (-2. * pow(Gammaa1, 4) * pow(ma1, 4) + 14. * pow(ma1, 8) +
                14. * pow(mpion, 8) + 28. * pow(mpion, 6) * pow(mrho, 2) +
                20. * pow(mpion, 4) * pow(mrho, 4) +
                10. * pow(mpion, 2) * pow(mrho, 6) + 2. * pow(mrho, 8) -
                16. * pow(mpion, 6) * s -
                16. * pow(mpion, 4) * pow(mrho, 2) * s -
                12. * pow(mpion, 2) * pow(mrho, 4) * s - 4. * pow(mrho, 6) * s -
                4. * pow(mpion, 4) * pow(s, 2) -
                6. * pow(mpion, 2) * pow(mrho, 2) * pow(s, 2) +
                8. * pow(mpion, 2) * pow(s, 3) + 4. * pow(mrho, 2) * pow(s, 3) -
                2. * pow(s, 4) +
                pow(ma1, 6) *
                    (-56. * pow(mpion, 2) - 28. * pow(mrho, 2) + 28. * s) +
                pow(ma1, 4) * (84. * pow(mpion, 4) + 24. * pow(mrho, 4) +
                               pow(mpion, 2) * (84. * pow(mrho, 2) - 72. * s) -
                               36. * pow(mrho, 2) * s + 12. * pow(s, 2)) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (-4. * pow(ma1, 4) - 4. * pow(mpion, 4) +
                     pow(ma1, 2) *
                         (8. * pow(mpion, 2) + 4. * pow(mrho, 2) - 4. * s) +
                     (4. * pow(mrho, 2) - 4. * s) * s +
                     pow(mpion, 2) * (-4. * pow(mrho, 2) + 8. * s)) +
                pow(ma1, 2) * (-56. * pow(mpion, 6) - 10. * pow(mrho, 6) +
                               18. * pow(mrho, 4) * s -
                               6. * pow(mrho, 2) * pow(s, 2) - 2. * pow(s, 3) +
                               pow(mpion, 4) * (-84. * pow(mrho, 2) + 60. * s) +
                               pow(mpion, 2) * (-48. * pow(mrho, 4) +
                                                60. * pow(mrho, 2) * s -
                                                12. * pow(s, 2)))) +
           pow(eta1, 2) *
               (1. * pow(Gammaa1, 4) * pow(ma1, 4) - 7. * pow(ma1, 8) -
                7. * pow(mpion, 8) - 14. * pow(mpion, 6) * pow(mrho, 2) -
                7. * pow(mpion, 4) * pow(mrho, 4) -
                2. * pow(mpion, 2) * pow(mrho, 6) +
                pow(ma1, 6) *
                    (28. * pow(mpion, 2) + 14. * pow(mrho, 2) - 14. * s) +
                8. * pow(mpion, 6) * s +
                11. * pow(mpion, 4) * pow(mrho, 2) * s +
                6. * pow(mpion, 2) * pow(mrho, 4) * s + 1. * pow(mrho, 6) * s +
                2. * pow(mpion, 4) * pow(s, 2) - 1. * pow(mrho, 4) * pow(s, 2) -
                4. * pow(mpion, 2) * pow(s, 3) - 1. * pow(mrho, 2) * pow(s, 3) +
                1. * pow(s, 4) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (2. * pow(ma1, 4) + 2. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                     pow(mpion, 2) * (2. * pow(mrho, 2) - 4. * s) -
                     1. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                     pow(ma1, 2) *
                         (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) +
                pow(ma1, 4) *
                    (-42. * pow(mpion, 4) - 9. * pow(mrho, 4) +
                     21. * pow(mrho, 2) * s - 6. * pow(s, 2) +
                     pow(mpion, 2) * (-42. * pow(mrho, 2) + 36. * s)) +
                pow(ma1, 2) * (28. * pow(mpion, 6) + 2. * pow(mrho, 6) +
                               pow(mpion, 4) * (42. * pow(mrho, 2) - 30. * s) -
                               9. * pow(mrho, 4) * s +
                               6. * pow(mrho, 2) * pow(s, 2) + 1. * pow(s, 3) +
                               pow(mpion, 2) *
                                   (18. * pow(mrho, 4) -
                                    36. * pow(mrho, 2) * s + 6. * pow(s, 2)))) +
           pow(eta2, 2) *
               (1. * pow(Gammaa1, 4) * pow(ma1, 4) - 7. * pow(ma1, 8) -
                7. * pow(mpion, 8) - 14. * pow(mpion, 6) * pow(mrho, 2) -
                1. * pow(mpion, 4) * pow(mrho, 4) +
                6. * pow(mpion, 2) * pow(mrho, 6) + 2. * pow(mrho, 8) +
                pow(ma1, 6) *
                    (28. * pow(mpion, 2) + 14. * pow(mrho, 2) - 14. * s) +
                8. * pow(mpion, 6) * s - 1. * pow(mpion, 4) * pow(mrho, 2) * s -
                16. * pow(mpion, 2) * pow(mrho, 4) * s - 7. * pow(mrho, 6) * s +
                2. * pow(mpion, 4) * pow(s, 2) +
                14. * pow(mpion, 2) * pow(mrho, 2) * pow(s, 2) +
                9. * pow(mrho, 4) * pow(s, 2) - 4. * pow(mpion, 2) * pow(s, 3) -
                5. * pow(mrho, 2) * pow(s, 3) + 1. * pow(s, 4) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (2. * pow(ma1, 4) + 2. * pow(mpion, 4) + 3. * pow(mrho, 4) +
                     pow(mpion, 2) * (2. * pow(mrho, 2) - 4. * s) -
                     5. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                     pow(ma1, 2) *
                         (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) +
                pow(ma1, 4) *
                    (-42. * pow(mpion, 4) - 3. * pow(mrho, 4) +
                     9. * pow(mrho, 2) * s - 6. * pow(s, 2) +
                     pow(mpion, 2) * (-42. * pow(mrho, 2) + 36. * s)) +
                pow(ma1, 2) * (28. * pow(mpion, 6) - 4. * pow(mrho, 6) +
                               pow(mpion, 4) * (42. * pow(mrho, 2) - 30. * s) +
                               9. * pow(mrho, 4) * s -
                               6. * pow(mrho, 2) * pow(s, 2) + 1. * pow(s, 3) +
                               pow(mpion, 2) *
                                   (6. * pow(mrho, 4) - 12. * pow(mrho, 2) * s +
                                    6. * pow(s, 2))))) *
          atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                tmax) /
               (Gammaa1 * ma1))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + 4. * pow(ma1, 4) +
              4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
              pow(mrho, 4) - 4. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s +
              pow(s, 2) +
              pow(ma1, 2) *
                  (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) +
         0.0625 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-4. * pow(ma1, 6) +
                   pow(ma1, 4) * (12. * pow(mpion, 2) - 6. * s) +
                   pow(mpion, 2) *
                       (4. * pow(mpion, 4) - 4. * pow(mrho, 4) -
                        2. * pow(mpion, 2) * s + 2. * pow(mrho, 2) * s) +
                   pow(ma1, 2) * (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                                  8. * pow(mpion, 2) * s +
                                  4. * pow(mrho, 2) * s - 4. * pow(s, 2))) +
              pow(eta1, 2) *
                  (2. * pow(ma1, 6) - 2. * pow(mpion, 6) +
                   pow(mpion, 2) * pow(mrho, 2) * s +
                   pow(mrho, 2) * (pow(mrho, 2) - 1. * s) * s +
                   pow(mpion, 4) * (-3. * pow(mrho, 2) + s) +
                   pow(ma1, 4) *
                       (-6. * pow(mpion, 2) - 3. * pow(mrho, 2) + 3. * s) +
                   pow(ma1, 2) * (6. * pow(mpion, 4) + pow(mrho, 4) +
                                  pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                                  4. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
              pow(eta2, 2) *
                  (2. * pow(ma1, 6) - 2. * pow(mpion, 6) -
                   1. * pow(mpion, 2) * pow(mrho, 2) * s +
                   pow(mpion, 4) * (3. * pow(mrho, 2) + s) +
                   pow(ma1, 4) *
                       (-6. * pow(mpion, 2) + 3. * pow(mrho, 2) + 3. * s) +
                   pow(ma1, 2) *
                       (6. * pow(mpion, 4) + pow(mrho, 4) +
                        pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                        2. * pow(mrho, 2) * s + 2. * pow(s, 2)))) *
             log(fabs(-1. * pow(ma1, 2) + tmax)) -
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * (-0.5 * pow(ma1, 6) - 0.5 * pow(mpion, 6) +
                   0.5 * pow(mpion, 4) * pow(mrho, 2) +
                   pow(ma1, 4) *
                       (0.5 * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 1. * s) +
                   pow(ma1, 2) * pow(mpion, 2) *
                       (0.5 * pow(mpion, 2) + 1. * pow(mrho, 2) - 1. * s)) +
           eta1 * (pow(ma1, 4) * (1. * pow(mpion, 2) + 0.5 * s) +
                   pow(mpion, 2) *
                       (1. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                        pow(mpion, 2) * (-1. * pow(mrho, 2) + 0.5 * s) -
                        0.5 * pow(mrho, 2) * s) +
                   pow(ma1, 2) *
                       (-2. * pow(mpion, 4) - 0.5 * pow(mrho, 2) * s +
                        pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s)))) *
          log(fabs(-1. * pow(ma1, 2) + tmax))) /
             (pow(ma1, 2) - 1. * pow(mpion, 2)) +
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * (-0.5 * pow(ma1, 6) - 0.5 * pow(mpion, 6) +
                   pow(ma1, 4) * (0.5 * pow(mpion, 2) + 0.5 * pow(mrho, 2)) +
                   pow(mpion, 4) * (0.5 * pow(mrho, 2) - 1. * s) +
                   pow(mpion, 2) * (-0.5 * pow(mrho, 2) + 0.5 * s) * s +
                   pow(ma1, 2) * (0.5 * pow(mpion, 4) +
                                  pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                                  (-0.5 * pow(mrho, 2) + 0.5 * s) * s)) +
           eta1 * (1. * pow(mpion, 6) +
                   pow(ma1, 4) * (1. * pow(mpion, 2) - 0.5 * s) +
                   pow(mpion, 2) * (1.5 * pow(mrho, 2) - 2. * s) * s +
                   (-0.5 * pow(mrho, 2) + 0.5 * s) * pow(s, 2) +
                   pow(mpion, 4) * (-1. * pow(mrho, 2) + 1.5 * s) +
                   pow(ma1, 2) *
                       (-2. * pow(mpion, 4) + 0.5 * pow(mrho, 2) * s +
                        pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s)))) *
          log(fabs(-1. * pow(ma1, 2) + tmax))) /
             (1. * pow(ma1, 2) - 1. * pow(mpion, 2) - 1. * pow(mrho, 2) +
              1. * s) -
         (0.03125 * pow(eta1 - 1. * eta2, 2) *
          (1. * pow(ma1, 2) - 1. * pow(mpion, 2) - 0.5 * pow(mrho, 2) +
           0.5 * s) *
          (eta1 * eta2 *
               (-2. * pow(ma1, 8) +
                pow(ma1, 6) *
                    (8. * pow(mpion, 2) + 4. * pow(mrho, 2) - 4. * s) +
                pow(ma1, 4) * (-12. * pow(mpion, 4) - 4. * pow(mrho, 4) +
                               4. * pow(mrho, 2) * s +
                               pow(mpion, 2) * (-12. * pow(mrho, 2) + 8. * s)) +
                pow(mpion, 2) *
                    (-2. * pow(mpion, 6) - 4. * pow(mpion, 4) * pow(mrho, 2) -
                     2. * pow(mrho, 6) + 4. * pow(mrho, 4) * s -
                     2. * pow(mrho, 2) * pow(s, 2) +
                     pow(mpion, 2) *
                         (-8. * pow(mrho, 4) + 8. * pow(mrho, 2) * s)) +
                pow(ma1, 2) * (8. * pow(mpion, 6) + 2. * pow(mrho, 6) +
                               pow(mpion, 4) * (12. * pow(mrho, 2) - 4. * s) -
                               2. * pow(mrho, 4) * s -
                               2. * pow(mrho, 2) * pow(s, 2) + 2. * pow(s, 3) +
                               pow(mpion, 2) *
                                   (8. * pow(mrho, 4) - 4. * pow(mrho, 2) * s -
                                    4. * pow(s, 2)))) +
           pow(eta2, 2) *
               (pow(ma1, 8) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                pow(mpion, 4) *
                    (pow(mpion, 4) + 2. * pow(mpion, 2) * pow(mrho, 2) +
                     pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                               pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) +
                               pow(mrho, 2) * s) +
                pow(ma1, 2) * (-4. * pow(mpion, 6) + 2. * pow(mrho, 6) -
                               5. * pow(mrho, 4) * s +
                               4. * pow(mrho, 2) * pow(s, 2) - 1. * pow(s, 3) +
                               pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s) +
                               pow(mpion, 2) *
                                   (2. * pow(mrho, 4) - 4. * pow(mrho, 2) * s +
                                    2. * pow(s, 2)))) +
           pow(eta1, 2) *
               (pow(ma1, 8) + pow(mpion, 8) +
                2. * pow(mpion, 6) * pow(mrho, 2) +
                pow(mpion, 2) * pow(mrho, 2) * s *
                    (-2. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                pow(mpion, 4) * (3. * pow(mrho, 4) - 5. * pow(mrho, 2) * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) + pow(mrho, 4) +
                               pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                               3. * pow(mrho, 2) * s) +
                pow(ma1, 2) *
                    (-4. * pow(mpion, 6) + pow(mrho, 4) * s - 1. * pow(s, 3) +
                     pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s) +
                     pow(mpion, 2) *
                         (-2. * pow(mrho, 4) + 4. * pow(mrho, 2) * s +
                          2. * pow(s, 2))))) *
          log(fabs(-1. * pow(ma1, 2) + tmax))) /
             (0.25 * pow(Gammaa1, 2) * pow(ma1, 2) + 1. * pow(ma1, 4) +
              1. * pow(mpion, 4) + 1. * pow(mpion, 2) * pow(mrho, 2) +
              0.25 * pow(mrho, 4) - 1. * pow(mpion, 2) * s -
              0.5 * pow(mrho, 2) * s + 0.25 * pow(s, 2) +
              pow(ma1, 2) *
                  (-2. * pow(mpion, 2) - 1. * pow(mrho, 2) + 1. * s)) -
         (1. * (1. * eta1 - 1. * eta2) *
          (eta2 *
               (pow(ma1, 4) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                pow(mpion, 2) * pow(mrho, 2) *
                    (pow(mpion, 2) * (0.5 - 1. * C4 * pow(mrho, 2)) +
                     (-0.25 + 0.125 * delta) * (pow(mrho, 2) + s)) +
                pow(ma1, 2) *
                    (-1. * C4 * pow(mrho, 6) +
                     pow(mpion, 2) *
                         (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4)) +
                     0.25 * delta * pow(s, 2) +
                     pow(mrho, 2) * s * (-0.25 - 0.375 * delta - 1. * C4 * s) +
                     pow(mrho, 4) * (0.75 - 0.125 * delta + 2. * C4 * s))) +
           eta1 *
               (pow(ma1, 4) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) +
                pow(ma1, 2) * (0.5 * pow(mrho, 4) - 1. * C4 * pow(mrho, 6) +
                               pow(mpion, 2) * (1. * pow(mrho, 2) -
                                                2. * C4 * pow(mrho, 4)) -
                               0.25 * delta * pow(s, 2) +
                               1. * C4 * pow(mrho, 2) * pow(s, 2)) +
                pow(mrho, 2) *
                    (pow(mpion, 4) * (-0.5 + 1. * C4 * pow(mrho, 2)) +
                     s * ((0.25 - 0.125 * delta) * pow(mrho, 2) +
                          (-0.25 + 0.125 * delta) * s) +
                     pow(mpion, 2) *
                         (2. * C4 * pow(mrho, 4) + (0.5 + 0.25 * delta) * s +
                          pow(mrho, 2) * (-1. - 2. * C4 * s))))) *
          log(fabs(-1. * pow(ma1, 2) + tmax))) /
             pow(mrho, 2) +
         0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
             log(fabs(-1. * pow(mpion, 2) + tmax)) +
         (0.25 *
          (0. +
           8.000000000000002 * pow(2. - 1. * delta, 2) * pow(mpion, 4) *
               pow(mrho, 2) -
           5.999999999999999 * pow(2. - 1. * delta, 2) * pow(mpion, 2) *
               pow(mrho, 2) * s +
           1. * pow(2. - 1. * delta, 2) * pow(mrho, 2) * pow(s, 2)) *
          log(fabs(-1. * pow(mpion, 2) + tmax))) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
         (0.125 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
          (0. + eta2 * pow(mpion, 2) * (4. * pow(mrho, 2) - 4. * s) +
           eta1 * (2. * pow(mrho, 4) - 2. * pow(mrho, 2) * s +
                   pow(mpion, 2) * (-4. * pow(mrho, 2) + 4. * s))) *
          log(fabs(-1. * pow(mpion, 2) + tmax))) /
             (pow(ma1, 2) - 1. * pow(mpion, 2)) +
         (2. * (-2. + 1. * delta) *
          (0. + (-0.25 + 0.125 * delta) * pow(mrho, 2) * s +
           pow(mpion, 2) * (-2. * C4 * pow(mrho, 4) - 0.5 * delta * s +
                            pow(mrho, 2) * (1. + 2. * C4 * s))) *
          log(fabs(-1. * pow(mpion, 2) + tmax))) /
             pow(mrho, 2) -
         (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
          (eta1 * (pow(mpion, 4) * (-1. * pow(mrho, 2) + 1. * s) +
                   pow(ma1, 2) * (pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                                  (-0.5 * pow(mrho, 2) + 0.5 * s) * s) +
                   pow(mpion, 2) * (-1. * pow(mrho, 4) +
                                    2.5 * pow(mrho, 2) * s - 1.5 * pow(s, 2)) +
                   s * (0.5 * pow(mrho, 4) - 1. * pow(mrho, 2) * s +
                        0.5 * pow(s, 2))) +
           eta2 * (0.5 * pow(mrho, 6) +
                   pow(mpion, 4) * (1. * pow(mrho, 2) - 1. * s) -
                   1.5 * pow(mrho, 4) * s + 1.5 * pow(mrho, 2) * pow(s, 2) -
                   0.5 * pow(s, 3) +
                   pow(mpion, 2) * (1.5 * pow(mrho, 4) - 3. * pow(mrho, 2) * s +
                                    1.5 * pow(s, 2)) +
                   pow(ma1, 2) *
                       (-0.5 * pow(mrho, 4) + 1. * pow(mrho, 2) * s -
                        0.5 * pow(s, 2) +
                        pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s)))) *
          log(fabs(-1. * pow(mpion, 2) + tmax))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) + pow(mpion, 4) +
              2. * pow(mpion, 2) * pow(mrho, 2) + pow(mrho, 4) -
              2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s + pow(s, 2) +
              pow(ma1, 2) *
                  (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) -
         0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
             log(fabs(-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + tmax)) +
         (0.5 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * pow(mpion, 6) * (1. * pow(mrho, 2) - 1. * s) +
           eta2 * pow(ma1, 2) * pow(mpion, 4) * (-1. * pow(mrho, 2) + 1. * s) +
           eta1 * pow(ma1, 2) * pow(mpion, 2) *
               (-0.5 * pow(mrho, 4) +
                pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                0.5 * pow(mrho, 2) * s) +
           eta1 * pow(mpion, 4) *
               (0.5 * pow(mrho, 4) - 0.5 * pow(mrho, 2) * s +
                pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s))) *
          log(fabs(-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + tmax))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
              2. * pow(ma1, 2) * pow(mpion, 2) + pow(mpion, 4)) +
         (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
          (eta1 * (pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                   (-0.5 * pow(mrho, 2) + 0.5 * s) * s) +
           eta2 *
               (-0.5 * pow(mrho, 4) + 1. * pow(mrho, 2) * s - 0.5 * pow(s, 2) +
                pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s))) *
          log(fabs(-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + tmax))) /
             (1. * pow(ma1, 2) - 1. * pow(mpion, 2) - 1. * pow(mrho, 2) +
              1. * s) -
         (0.25 *
          (0. +
           8.000000000000002 * pow(2. - 1. * delta, 2) * pow(mpion, 4) *
               pow(mrho, 2) +
           1. * pow(2. - 1. * delta, 2) * pow(mrho, 4) * s +
           pow(mpion, 2) * (C4 * (32. - 16. * delta) * pow(mrho, 6) +
                            delta * (-8. + 4. * delta) * pow(s, 2) +
                            pow(mrho, 2) * s *
                                (-8. + 24. * delta - 10. * pow(delta, 2) +
                                 32. * C4 * s - 16. * C4 * delta * s) +
                            pow(mrho, 4) * (-16. + 8. * delta - 64. * C4 * s +
                                            32. * C4 * delta * s))) *
          log(fabs(-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + tmax))) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
         0.03125 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (4. * pow(ma1, 6) +
                   pow(Gammaa1, 2) * pow(ma1, 2) *
                       (-4. * pow(ma1, 2) + 4. * pow(mpion, 2) - 2. * s) +
                   pow(ma1, 4) * (-12. * pow(mpion, 2) + 6. * s) +
                   pow(mpion, 2) *
                       (-4. * pow(mpion, 4) + 4. * pow(mrho, 4) +
                        2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s) +
                   pow(ma1, 2) * (12. * pow(mpion, 4) - 2. * pow(mrho, 4) -
                                  8. * pow(mpion, 2) * s -
                                  4. * pow(mrho, 2) * s + 4. * pow(s, 2))) +
              pow(eta1, 2) *
                  (-2. * pow(ma1, 6) + 2. * pow(mpion, 6) +
                   3. * pow(mpion, 4) * pow(mrho, 2) +
                   pow(ma1, 4) *
                       (6. * pow(mpion, 2) + 3. * pow(mrho, 2) - 3. * s) -
                   1. * pow(mpion, 4) * s -
                   1. * pow(mpion, 2) * pow(mrho, 2) * s -
                   1. * pow(mrho, 4) * s + pow(mrho, 2) * pow(s, 2) +
                   pow(Gammaa1, 2) * pow(ma1, 2) *
                       (2. * pow(ma1, 2) - 2. * pow(mpion, 2) -
                        1. * pow(mrho, 2) + s) +
                   pow(ma1, 2) *
                       (-6. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                        4. * pow(mrho, 2) * s - 2. * pow(s, 2) +
                        pow(mpion, 2) * (-6. * pow(mrho, 2) + 4. * s))) +
              pow(eta2, 2) *
                  (-2. * pow(ma1, 6) + 2. * pow(mpion, 6) -
                   3. * pow(mpion, 4) * pow(mrho, 2) +
                   pow(ma1, 4) *
                       (6. * pow(mpion, 2) - 3. * pow(mrho, 2) - 3. * s) -
                   1. * pow(mpion, 4) * s + pow(mpion, 2) * pow(mrho, 2) * s +
                   pow(Gammaa1, 2) * pow(ma1, 2) *
                       (2. * pow(ma1, 2) - 2. * pow(mpion, 2) + pow(mrho, 2) +
                        s) +
                   pow(ma1, 2) *
                       (-6. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                        2. * pow(mrho, 2) * s - 2. * pow(s, 2) +
                        pow(mpion, 2) * (6. * pow(mrho, 2) + 4. * s)))) *
             log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
                      4. * pow(ma1, 2) * pow(mpion, 2) + 4. * pow(mpion, 4) +
                      2. * pow(ma1, 2) * (-1. * pow(mrho, 2) + s + tmax) -
                      4. * pow(mpion, 2) * (-1. * pow(mrho, 2) + s + tmax) +
                      pow(-1. * pow(mrho, 2) + s + tmax, 2))) -
         (0.5 * (1. * eta1 - 1. * eta2) *
          (eta2 *
               (pow(Gammaa1, 2) * pow(ma1, 2) * pow(mrho, 2) *
                    (0.5 - 1. * C4 * pow(mrho, 2)) +
                pow(ma1, 4) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) +
                pow(mpion, 2) * pow(mrho, 2) *
                    (pow(mpion, 2) * (-0.5 + 1. * C4 * pow(mrho, 2)) +
                     (0.25 - 0.125 * delta) * (pow(mrho, 2) + s)) +
                pow(ma1, 2) *
                    (1. * C4 * pow(mrho, 6) +
                     pow(mpion, 2) *
                         (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) -
                     0.25 * delta * pow(s, 2) +
                     pow(mrho, 4) * (-0.75 + 0.125 * delta - 2. * C4 * s) +
                     pow(mrho, 2) * s * (0.25 + 0.375 * delta + 1. * C4 * s))) +
           eta1 * (pow(Gammaa1, 2) * pow(ma1, 2) * pow(mrho, 2) *
                       (-0.5 + 1. * C4 * pow(mrho, 2)) +
                   pow(ma1, 4) * (0.5 * pow(mrho, 2) - 1. * C4 * pow(mrho, 4)) +
                   pow(ma1, 2) * (-0.5 * pow(mrho, 4) + 1. * C4 * pow(mrho, 6) +
                                  pow(mpion, 2) * (-1. * pow(mrho, 2) +
                                                   2. * C4 * pow(mrho, 4)) +
                                  0.25 * delta * pow(s, 2) -
                                  1. * C4 * pow(mrho, 2) * pow(s, 2)) +
                   pow(mrho, 2) *
                       (pow(mpion, 4) * (0.5 - 1. * C4 * pow(mrho, 2)) +
                        s * ((-0.25 + 0.125 * delta) * pow(mrho, 2) +
                             (0.25 - 0.125 * delta) * s) +
                        pow(mpion, 2) * (-2. * C4 * pow(mrho, 4) +
                                         (-0.5 - 0.25 * delta) * s +
                                         pow(mrho, 2) * (1. + 2. * C4 * s))))) *
          log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) +
                   pow(pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) +
                           s + tmax,
                       2)))) /
             pow(mrho, 2) +
         (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * (0.5 * pow(Gammaa1, 4) * pow(ma1, 4) - 0.5 * pow(ma1, 8) +
                   0.5 * pow(mpion, 8) +
                   0.5 * pow(ma1, 4) * pow(mpion, 2) * pow(mrho, 2) -
                   0.5 * pow(mpion, 6) * pow(mrho, 2) +
                   pow(Gammaa1, 2) *
                       (pow(ma1, 2) * pow(mpion, 2) *
                            (1. * pow(mpion, 2) + 1.5 * pow(mrho, 2) - 2. * s) +
                        pow(ma1, 4) * (-1. * pow(mpion, 2) +
                                       0.5 * pow(mrho, 2) - 1. * s)) +
                   pow(ma1, 6) *
                       (1. * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 1. * s) +
                   pow(ma1, 2) * pow(mpion, 4) *
                       (-1. * pow(mpion, 2) - 0.5 * pow(mrho, 2) + 1. * s)) +
           eta1 *
               (pow(ma1, 6) * (1. * pow(mpion, 2) + 0.5 * s) +
                pow(ma1, 2) *
                    (3. * pow(mpion, 6) + 1. * pow(mpion, 2) * pow(mrho, 4) -
                     0.5 * pow(mpion, 4) * s) +
                pow(ma1, 4) * (-3. * pow(mpion, 4) +
                               pow(mpion, 2) * (-1. * pow(mrho, 2) + 0.5 * s) -
                               0.5 * pow(mrho, 2) * s) +
                pow(mpion, 4) * (-1. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                                 pow(mpion, 2) * (1. * pow(mrho, 2) - 0.5 * s) +
                                 0.5 * pow(mrho, 2) * s) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (-1. * pow(mpion, 4) +
                     pow(ma1, 2) * (1. * pow(mpion, 2) + 0.5 * s) -
                     0.5 * pow(mrho, 2) * s +
                     pow(mpion, 2) * (-1. * pow(mrho, 2) + 1.5 * s)))) *
          log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) +
                   4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                   pow(mrho, 4) - 4. * pow(mpion, 2) * s -
                   2. * pow(mrho, 2) * s + pow(s, 2) -
                   4. * pow(mpion, 2) * tmax - 2. * pow(mrho, 2) * tmax +
                   2. * s * tmax + pow(tmax, 2) +
                   pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                  2. * s + 2. * tmax)))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
              2. * pow(ma1, 2) * pow(mpion, 2) + pow(mpion, 4)) -
         (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta1 *
               (-1. * pow(mpion, 8) + 1. * pow(mpion, 4) * pow(mrho, 4) +
                pow(ma1, 6) * (1. * pow(mpion, 2) - 0.5 * s) -
                0.5 * pow(mpion, 6) * s -
                4. * pow(mpion, 4) * pow(mrho, 2) * s -
                1.5 * pow(mpion, 2) * pow(mrho, 4) * s +
                3.5 * pow(mpion, 4) * pow(s, 2) +
                4. * pow(mpion, 2) * pow(mrho, 2) * pow(s, 2) +
                0.5 * pow(mrho, 4) * pow(s, 2) -
                2.5 * pow(mpion, 2) * pow(s, 3) -
                1. * pow(mrho, 2) * pow(s, 3) + 0.5 * pow(s, 4) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (-1. * pow(mpion, 4) +
                     pow(ma1, 2) * (1. * pow(mpion, 2) - 0.5 * s) -
                     0.5 * pow(mpion, 2) * s + 0.5 * pow(s, 2)) +
                pow(ma1, 4) *
                    (-3. * pow(mpion, 4) + (1. * pow(mrho, 2) - 0.5 * s) * s +
                     pow(mpion, 2) * (-2. * pow(mrho, 2) + 2.5 * s)) +
                pow(ma1, 2) * (3. * pow(mpion, 6) +
                               pow(mpion, 4) * (2. * pow(mrho, 2) - 1.5 * s) -
                               0.5 * pow(mrho, 4) * s + 0.5 * pow(s, 3) +
                               pow(mpion, 2) *
                                   (1. * pow(mrho, 4) - 1. * pow(mrho, 2) * s -
                                    1. * pow(s, 2)))) +
           eta2 *
               (0.5 * pow(Gammaa1, 4) * pow(ma1, 4) - 0.5 * pow(ma1, 8) +
                0.5 * pow(mpion, 8) +
                pow(ma1, 6) *
                    (1. * pow(mpion, 2) + 1. * pow(mrho, 2) - 0.5 * s) +
                0.5 * pow(mpion, 6) * s +
                pow(ma1, 4) * (-0.5 * pow(mrho, 4) +
                               (-0.5 * pow(mpion, 2) + 0.5 * s) * s) +
                pow(mpion, 4) * (-0.5 * pow(mrho, 4) + 2. * pow(mrho, 2) * s -
                                 1.5 * pow(s, 2)) +
                pow(mpion, 2) * s *
                    (0.5 * pow(mrho, 4) - 1. * pow(mrho, 2) * s +
                     0.5 * pow(s, 2)) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (1. * pow(mpion, 4) + 0.5 * pow(mrho, 4) +
                     pow(mpion, 2) * (2. * pow(mrho, 2) - 1.5 * s) -
                     1. * pow(mrho, 2) * s + 0.5 * pow(s, 2) +
                     pow(ma1, 2) *
                         (-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + 1.5 * s)) +
                pow(ma1, 2) *
                    (-1. * pow(mpion, 6) +
                     pow(mpion, 4) * (-1. * pow(mrho, 2) + 0.5 * s) +
                     pow(mpion, 2) * (-1. * pow(mrho, 4) +
                                      2. * pow(mrho, 2) * s - 1. * pow(s, 2)) +
                     s * (0.5 * pow(mrho, 4) - 1. * pow(mrho, 2) * s +
                          0.5 * pow(s, 2))))) *
          log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) +
                   4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                   pow(mrho, 4) - 4. * pow(mpion, 2) * s -
                   2. * pow(mrho, 2) * s + pow(s, 2) -
                   4. * pow(mpion, 2) * tmax - 2. * pow(mrho, 2) * tmax +
                   2. * s * tmax + pow(tmax, 2) +
                   pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                  2. * s + 2. * tmax)))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) + pow(mpion, 4) +
              2. * pow(mpion, 2) * pow(mrho, 2) + pow(mrho, 4) -
              2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s + pow(s, 2) +
              pow(ma1, 2) *
                  (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) -
         (0.0625 * pow(eta1 - 1. * eta2, 2) *
          (pow(eta2, 2) *
               (-1. * pow(ma1, 10) +
                pow(ma1, 8) *
                    (5. * pow(mpion, 2) + 2.5 * pow(mrho, 2) - 2.5 * s) +
                pow(Gammaa1, 4) * pow(ma1, 4) *
                    (1. * pow(ma1, 2) - 1. * pow(mpion, 2) -
                     0.5 * pow(mrho, 2) + 0.5 * s) +
                pow(ma1, 4) *
                    (10. * pow(mpion, 6) - 2.5 * pow(mrho, 6) +
                     pow(mpion, 4) * (15. * pow(mrho, 2) - 9. * s) +
                     6. * pow(mrho, 4) * s - 4.5 * pow(mrho, 2) * pow(s, 2) +
                     1. * pow(s, 3)) +
                pow(ma1, 6) *
                    (-10. * pow(mpion, 4) + (1. * pow(mrho, 2) - 1. * s) * s +
                     pow(mpion, 2) * (-10. * pow(mrho, 2) + 8. * s)) +
                pow(mpion, 4) *
                    (1. * pow(mpion, 6) + 0.5 * pow(mrho, 6) +
                     pow(mpion, 4) * (2.5 * pow(mrho, 2) - 0.5 * s) -
                     1. * pow(mrho, 4) * s + 0.5 * pow(mrho, 2) * pow(s, 2) +
                     pow(mpion, 2) *
                         (2. * pow(mrho, 4) - 2. * pow(mrho, 2) * s)) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (4. * pow(ma1, 6) - 4. * pow(mpion, 6) -
                     0.5 * pow(mrho, 6) + 1.5 * pow(mrho, 4) * s -
                     1.5 * pow(mrho, 2) * pow(s, 2) + 0.5 * pow(s, 3) +
                     pow(mpion, 4) * (-6. * pow(mrho, 2) + 6. * s) +
                     pow(ma1, 4) *
                         (-12. * pow(mpion, 2) - 6. * pow(mrho, 2) + 6. * s) +
                     pow(mpion, 2) * (-3. * pow(mrho, 4) +
                                      6. * pow(mrho, 2) * s - 3. * pow(s, 2)) +
                     pow(ma1, 2) *
                         (12. * pow(mpion, 4) + 3. * pow(mrho, 4) +
                          pow(mpion, 2) * (12. * pow(mrho, 2) - 12. * s) -
                          6. * pow(mrho, 2) * s + 3. * pow(s, 2))) +
                pow(ma1, 2) *
                    (-5. * pow(mpion, 8) + 1. * pow(mrho, 8) -
                     3.5 * pow(mrho, 6) * s + 4.5 * pow(mrho, 4) * pow(s, 2) -
                     2.5 * pow(mrho, 2) * pow(s, 3) + 0.5 * pow(s, 4) +
                     pow(mpion, 6) * (-10. * pow(mrho, 2) + 4. * s) +
                     pow(mpion, 4) * (-2. * pow(mrho, 4) +
                                      1. * pow(mrho, 2) * s + 1. * pow(s, 2)) +
                     pow(mpion, 2) *
                         (3. * pow(mrho, 6) - 8. * pow(mrho, 4) * s +
                          7. * pow(mrho, 2) * pow(s, 2) - 2. * pow(s, 3)))) +
           pow(eta1, 2) *
               (-1. * pow(ma1, 10) +
                pow(ma1, 8) *
                    (5. * pow(mpion, 2) + 2.5 * pow(mrho, 2) - 2.5 * s) +
                pow(Gammaa1, 4) * pow(ma1, 4) *
                    (1. * pow(ma1, 2) - 1. * pow(mpion, 2) -
                     0.5 * pow(mrho, 2) + 0.5 * s) +
                pow(ma1, 6) * (-10. * pow(mpion, 4) - 2. * pow(mrho, 4) +
                               5. * pow(mrho, 2) * s - 1. * pow(s, 2) +
                               pow(mpion, 2) * (-10. * pow(mrho, 2) + 8. * s)) +
                pow(ma1, 4) * (10. * pow(mpion, 6) + 0.5 * pow(mrho, 6) +
                               pow(mpion, 4) * (15. * pow(mrho, 2) - 9. * s) -
                               3. * pow(mrho, 4) * s +
                               1.5 * pow(mrho, 2) * pow(s, 2) + 1. * pow(s, 3) +
                               pow(mpion, 2) * (6. * pow(mrho, 4) -
                                                12. * pow(mrho, 2) * s)) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (4. * pow(ma1, 6) - 4. * pow(mpion, 6) -
                     0.5 * pow(mrho, 6) + 1.5 * pow(mrho, 4) * s -
                     1.5 * pow(mrho, 2) * pow(s, 2) + 0.5 * pow(s, 3) +
                     pow(mpion, 4) * (-6. * pow(mrho, 2) + 6. * s) +
                     pow(ma1, 4) *
                         (-12. * pow(mpion, 2) - 6. * pow(mrho, 2) + 6. * s) +
                     pow(mpion, 2) * (-3. * pow(mrho, 4) +
                                      6. * pow(mrho, 2) * s - 3. * pow(s, 2)) +
                     pow(ma1, 2) *
                         (12. * pow(mpion, 4) + 3. * pow(mrho, 4) +
                          pow(mpion, 2) * (12. * pow(mrho, 2) - 12. * s) -
                          6. * pow(mrho, 2) * s + 3. * pow(s, 2))) +
                pow(mpion, 2) *
                    (1. * pow(mpion, 8) +
                     pow(mpion, 6) * (2.5 * pow(mrho, 2) - 0.5 * s) +
                     pow(mpion, 4) *
                         (4. * pow(mrho, 4) - 6. * pow(mrho, 2) * s) +
                     pow(mrho, 2) * s *
                         (-1. * pow(mrho, 4) + 2. * pow(mrho, 2) * s -
                          1. * pow(s, 2)) +
                     pow(mpion, 2) *
                         (1.5 * pow(mrho, 6) - 6. * pow(mrho, 4) * s +
                          4.5 * pow(mrho, 2) * pow(s, 2))) +
                pow(ma1, 2) *
                    (-5. * pow(mpion, 8) +
                     pow(mpion, 6) * (-10. * pow(mrho, 2) + 4. * s) +
                     pow(mpion, 4) * (-8. * pow(mrho, 4) +
                                      13. * pow(mrho, 2) * s + 1. * pow(s, 2)) +
                     pow(mpion, 2) *
                         (-1. * pow(mrho, 6) + 6. * pow(mrho, 4) * s -
                          3. * pow(mrho, 2) * pow(s, 2) - 2. * pow(s, 3)) +
                     s * (0.5 * pow(mrho, 6) - 0.5 * pow(mrho, 4) * s -
                          0.5 * pow(mrho, 2) * pow(s, 2) + 0.5 * pow(s, 3)))) +
           eta1 * eta2 *
               (2. * pow(ma1, 10) +
                pow(Gammaa1, 4) * pow(ma1, 4) *
                    (-2. * pow(ma1, 2) + 2. * pow(mpion, 2) +
                     1. * pow(mrho, 2) - 1. * s) +
                pow(ma1, 8) *
                    (-10. * pow(mpion, 2) - 5. * pow(mrho, 2) + 5. * s) +
                pow(ma1, 6) * (20. * pow(mpion, 4) + 6. * pow(mrho, 4) +
                               pow(mpion, 2) * (20. * pow(mrho, 2) - 16. * s) -
                               8. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                pow(ma1, 4) * (-20. * pow(mpion, 6) - 4. * pow(mrho, 6) +
                               6. * pow(mrho, 4) * s - 2. * pow(s, 3) +
                               pow(mpion, 4) * (-30. * pow(mrho, 2) + 18. * s) +
                               pow(mpion, 2) * (-18. * pow(mrho, 4) +
                                                18. * pow(mrho, 2) * s)) +
                pow(mpion, 2) *
                    (-2. * pow(mpion, 8) - 1. * pow(mrho, 8) +
                     3. * pow(mrho, 6) * s - 3. * pow(mrho, 4) * pow(s, 2) +
                     1. * pow(mrho, 2) * pow(s, 3) +
                     pow(mpion, 6) * (-5. * pow(mrho, 2) + 1. * s) +
                     pow(mpion, 4) *
                         (-10. * pow(mrho, 4) + 10. * pow(mrho, 2) * s) +
                     pow(mpion, 2) *
                         (-6. * pow(mrho, 6) + 12. * pow(mrho, 4) * s -
                          6. * pow(mrho, 2) * pow(s, 2))) +
                pow(ma1, 2) *
                    (10. * pow(mpion, 8) + 1. * pow(mrho, 8) +
                     pow(mpion, 6) * (20. * pow(mrho, 2) - 8. * s) -
                     2. * pow(mrho, 6) * s + 2. * pow(mrho, 2) * pow(s, 3) -
                     1. * pow(s, 4) +
                     pow(mpion, 4) * (22. * pow(mrho, 4) -
                                      20. * pow(mrho, 2) * s - 2. * pow(s, 2)) +
                     pow(mpion, 2) *
                         (8. * pow(mrho, 6) - 12. * pow(mrho, 4) * s +
                          4. * pow(s, 3))) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (-8. * pow(ma1, 6) + 8. * pow(mpion, 6) +
                     1. * pow(mrho, 6) +
                     pow(mpion, 4) * (12. * pow(mrho, 2) - 12. * s) +
                     pow(ma1, 4) *
                         (24. * pow(mpion, 2) + 12. * pow(mrho, 2) - 12. * s) -
                     3. * pow(mrho, 4) * s + 3. * pow(mrho, 2) * pow(s, 2) -
                     1. * pow(s, 3) +
                     pow(mpion, 2) * (6. * pow(mrho, 4) -
                                      12. * pow(mrho, 2) * s + 6. * pow(s, 2)) +
                     pow(ma1, 2) *
                         (-24. * pow(mpion, 4) - 6. * pow(mrho, 4) +
                          12. * pow(mrho, 2) * s - 6. * pow(s, 2) +
                          pow(mpion, 2) * (-24. * pow(mrho, 2) + 24. * s))))) *
          log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) +
                   4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                   pow(mrho, 4) - 4. * pow(mpion, 2) * s -
                   2. * pow(mrho, 2) * s + pow(s, 2) -
                   4. * pow(mpion, 2) * tmax - 2. * pow(mrho, 2) * tmax +
                   2. * s * tmax + pow(tmax, 2) +
                   pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                  2. * s + 2. * tmax)))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + 4. * pow(ma1, 4) +
              4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
              pow(mrho, 4) - 4. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s +
              pow(s, 2) +
              pow(ma1, 2) *
                  (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)))) /
           (16. * Pi * s * (-4 * pow(mpion, 2) + s)));

  return xs * to_mb / spin_deg_factor;
}

double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_pi_rho0(
    const double s, const double t, const double m_rho) {
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  const double spin_deg_factor = 1.0;

  const double diff_xs =
      ((pow(Const, 2) * pow(ghat, 4) *
        ((0.25 *
          (32 * pow(C4, 2) * pow(mrho, 8) + 2 * pow(delta, 2) * pow(s, 2) +
           8 * C4 * pow(mrho, 6) * (-6 + delta - 8 * C4 * s) +
           2 * delta * pow(mrho, 2) * s * (-6 + delta - 8 * C4 * s) +
           pow(mrho, 4) * (12 - pow(delta, 2) + 8 * C4 * (6 + delta) * s +
                           32 * pow(C4, 2) * pow(s, 2)))) /
             pow(mrho, 4) -
         (0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
          (pow(mpion, 4) + pow(pow(mrho, 2) - t, 2) -
           2 * pow(mpion, 2) * (pow(mrho, 2) + t))) /
             (pow(mrho, 2) * pow(pow(mpion, 2) - t, 2)) -
         (0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
          (pow(mpion, 4) + pow(s + t, 2) -
           2 * pow(mpion, 2) * (2 * pow(mrho, 2) + s + t))) /
             (pow(mrho, 2) * pow(pow(mpion, 2) + pow(mrho, 2) - s - t, 2)) +
         (0.125 * (-2 + delta) * (eta1 - eta2) *
          (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t) *
          (eta1 * (2 * pow(mpion, 2) - s) +
           eta2 * (-3 * pow(mpion, 2) - pow(mrho, 2) + s + t)) *
          (pow(mpion, 4) + t * (-pow(mrho, 2) + 2 * s + t) -
           pow(mpion, 2) * (pow(mrho, 2) + 2 * t))) /
             ((-pow(mpion, 2) + t) *
              (pow(Gammaa1, 2) * pow(ma1, 2) +
               pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t,
                   2))) -
         (0.125 * (-2 + delta) * (eta1 - eta2) *
          (eta1 * (-2 * pow(mpion, 2) + s) + eta2 * (pow(mpion, 2) + t)) *
          (-pow(mpion, 4) + pow(s, 2) - pow(t, 2) + pow(mrho, 2) * (-s + t) +
           pow(mpion, 2) * (pow(mrho, 2) - 2 * s + 2 * t))) /
             ((pow(mpion, 2) + pow(mrho, 2) - s - t) * (-pow(ma1, 2) + t)) +
         (0.25 * (-2. + delta) *
          (pow(mpion, 4) * (2. + delta - 8. * C4 * pow(mrho, 2)) +
           8. * C4 * pow(mrho, 4) * t +
           t * ((2. + 3. * delta) * s + (2. + delta) * t) +
           pow(mrho, 2) * (s * (2. - 1. * delta - 16. * C4 * t) +
                           t * (-2. - 1. * delta - 8. * C4 * t)) +
           pow(mpion, 2) * (8. * C4 * pow(mrho, 4) + (-2. + delta) * s +
                            (-4. - 2. * delta) * t +
                            pow(mrho, 2) * (-6. + delta + 16. * C4 * t)))) /
             (pow(mrho, 2) * (pow(mpion, 2) - 1. * t)) -
         (0.125 * (-2 + delta) * (eta1 - eta2) *
          (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t) *
          (-(eta2 * (3 * pow(mpion, 2) + pow(mrho, 2) - s - t) *
             (pow(mpion, 4) + (pow(mrho, 2) - s - t) * (s - t) -
              pow(mpion, 2) * (pow(mrho, 2) - 2 * s + 2 * t))) +
           eta1 *
               (2 * pow(mpion, 6) +
                pow(mpion, 4) * (-2 * pow(mrho, 2) + 5 * s - 4 * t) +
                s * (s + t) * (-pow(mrho, 2) + s + t) +
                pow(mpion, 2) * (2 * pow(mrho, 4) + pow(mrho, 2) * (s - 2 * t) -
                                 2 * (2 * s - t) * (s + t))))) /
             ((-pow(mpion, 2) - pow(mrho, 2) + s + t) *
              (pow(Gammaa1, 2) * pow(ma1, 2) +
               pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t,
                   2))) +
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * (-0.5 * pow(mpion, 6) +
                   pow(mpion, 4) * (0.5 * pow(mrho, 2) + 0.5 * t) +
                   pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s + 0.5 * t) * t +
                   (0.5 * pow(mrho, 2) - 1. * s - 0.5 * t) * pow(t, 2)) +
           eta1 * (1. * pow(mpion, 6) +
                   pow(mpion, 4) * (-1. * pow(mrho, 2) + 0.5 * s - 2. * t) +
                   s * (-0.5 * pow(mrho, 2) + 0.5 * t) * t +
                   pow(mpion, 2) *
                       (1. * pow(mrho, 4) + pow(mrho, 2) * (-0.5 * s - 1. * t) +
                        t * (1. * s + 1. * t))))) /
             ((pow(ma1, 2) - 1. * t) * (-1. * pow(mpion, 2) + t)) +
         (0.03125 * pow(eta1 - eta2, 2) *
          (-2 * eta1 * eta2 *
               (pow(mpion, 8) - 4 * pow(mpion, 6) * t +
                pow(t, 2) * (-pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                             2 * pow(s, 2) + 2 * s * t + pow(t, 2)) -
                2 * pow(mpion, 2) * t *
                    (-2 * pow(mrho, 4) + pow(mrho, 2) * s + 2 * t * (s + t)) +
                pow(mpion, 4) * (-pow(mrho, 4) + 2 * t * (s + 3 * t))) +
           pow(eta2, 2) *
               (pow(mpion, 8) - 2 * pow(mpion, 6) * (pow(mrho, 2) + 2 * t) +
                pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t +
                             pow(t, 2) + 2 * pow(mrho, 2) * (-s + t)) -
                2 * pow(mpion, 2) * t *
                    (2 * t * (s + t) + pow(mrho, 2) * (s + 3 * t)) +
                pow(mpion, 4) * (pow(mrho, 4) + 6 * pow(mrho, 2) * t +
                                 2 * t * (s + 3 * t))) +
           pow(eta1, 2) *
               (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) -
                2 * pow(mpion, 2) * (pow(mrho, 2) - s - t) *
                    (pow(mrho, 4) + pow(mrho, 2) * t - 2 * pow(t, 2)) +
                t * (-pow(mrho, 2) + t) *
                    (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                     pow(mrho, 2) * (2 * s + t)) +
                pow(mpion, 4) * (pow(mrho, 4) - 2 * pow(mrho, 2) * (s + 3 * t) +
                                 2 * t * (s + 3 * t))))) /
             pow(pow(ma1, 2) - t, 2) +
         ((0.25 * pow(-2 + delta, 2) * (2 * pow(mpion, 2) - s) *
           (pow(mpion, 4) + pow(mrho, 2) * (s - t) + t * (s + t) -
            pow(mpion, 2) * (3 * pow(mrho, 2) + s + 2 * t))) /
              ((pow(mpion, 2) - t) * (pow(mpion, 2) + pow(mrho, 2) - s - t)) -
          (0.25 * (-2. + delta) *
           (pow(mpion, 4) * (2. + delta - 8. * C4 * pow(mrho, 2)) -
            2. * delta * pow(s, 2) + 2. * s * t - 1. * delta * s * t +
            2. * pow(t, 2) + delta * pow(t, 2) +
            C4 * pow(mrho, 4) * (-8. * s + 8. * t) +
            pow(mrho, 2) * ((2. + delta) * s + 8. * C4 * pow(s, 2) +
                            t * (-2. - 1. * delta - 8. * C4 * t)) +
            pow(mpion, 2) *
                (8. * C4 * pow(mrho, 4) - 2. * s + 5. * delta * s - 4. * t -
                 2. * delta * t +
                 pow(mrho, 2) * (-6. + delta - 16. * C4 * s + 16. * C4 * t)))) /
              (pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t)) /
             pow(mrho, 2) +
         (0.03125 * pow(eta1 - eta2, 2) *
          (-2 * eta1 * eta2 *
               (pow(mpion, 8) + 4 * pow(mpion, 6) * (pow(mrho, 2) - t) +
                pow(-pow(mrho, 2) + s + t, 2) *
                    (pow(s, 2) + pow(t, 2) - 2 * pow(mrho, 2) * (s + t)) +
                pow(mpion, 4) *
                    (9 * pow(mrho, 4) + 4 * pow(s, 2) + 2 * s * t +
                     6 * pow(t, 2) - 2 * pow(mrho, 2) * (7 * s + 6 * t)) +
                2 * pow(mpion, 2) * (pow(mrho, 2) - s - t) *
                    (2 * pow(mrho, 4) - pow(mrho, 2) * (5 * s + 4 * t) +
                     2 * (pow(s, 2) + pow(t, 2)))) +
           pow(eta2, 2) *
               (pow(mpion, 8) + pow(mpion, 6) * (6 * pow(mrho, 2) - 4 * t) +
                pow(-pow(mrho, 2) + s + t, 2) *
                    (4 * pow(mrho, 4) + pow(s, 2) + pow(t, 2) -
                     4 * pow(mrho, 2) * (s + t)) +
                pow(mpion, 4) *
                    (17 * pow(mrho, 4) + 4 * pow(s, 2) + 2 * s * t +
                     6 * pow(t, 2) - 2 * pow(mrho, 2) * (10 * s + 9 * t)) +
                2 * pow(mpion, 2) * (pow(mrho, 2) - s - t) *
                    (7 * pow(mrho, 4) - pow(mrho, 2) * (8 * s + 7 * t) +
                     2 * (pow(s, 2) + pow(t, 2)))) +
           pow(eta1, 2) *
               (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) +
                (s + t) * (-pow(mrho, 2) + s + t) *
                    (pow(s, 2) + pow(t, 2) - pow(mrho, 2) * (s + t)) +
                pow(mpion, 4) *
                    (5 * pow(mrho, 4) + 4 * pow(s, 2) + 2 * s * t +
                     6 * pow(t, 2) - 2 * pow(mrho, 2) * (5 * s + 3 * t)) -
                2 * pow(mpion, 2) *
                    (2 * pow(mrho, 4) * (s + t) +
                     2 * (s + t) * (pow(s, 2) + pow(t, 2)) -
                     pow(mrho, 2) *
                         (4 * pow(s, 2) + 5 * s * t + 3 * pow(t, 2)))))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) +
              pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t, 2)) +
         (0.0625 * pow(eta1 - eta2, 2) *
          (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t) *
          (-(pow(eta2, 2) *
             (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) +
              2 * pow(mpion, 2) * t *
                  (pow(pow(mrho, 2) - s, 2) + (3 * pow(mrho, 2) - 2 * s) * t -
                   2 * pow(t, 2)) +
              (pow(mrho, 2) - s - t) * t *
                  (2 * pow(mrho, 4) + pow(s, 2) - s * t - pow(t, 2) +
                   pow(mrho, 2) * (-3 * s + t)) +
              pow(mpion, 4) * (pow(mrho, 4) + 2 * t * (s + 3 * t) -
                               pow(mrho, 2) * (s + 6 * t)))) -
           pow(eta1, 2) *
               (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) +
                (pow(mrho, 2) - s - t) * t *
                    (pow(s, 2) - s * t - pow(t, 2) + pow(mrho, 2) * (s + t)) +
                pow(mpion, 4) * (3 * pow(mrho, 4) + 2 * t * (s + 3 * t) -
                                 pow(mrho, 2) * (5 * s + 6 * t)) +
                2 * pow(mpion, 2) *
                    (-(pow(mrho, 4) * (s + t)) +
                     t * (pow(s, 2) - 2 * s * t - 2 * pow(t, 2)) +
                     pow(mrho, 2) * (pow(s, 2) + 2 * s * t + 3 * pow(t, 2)))) +
           2 * eta1 * eta2 *
               (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) -
                (pow(mrho, 2) - s - t) * t *
                    (pow(mrho, 4) - pow(s, 2) - pow(mrho, 2) * t +
                     t * (s + t)) +
                2 * pow(mpion, 4) *
                    (2 * pow(mrho, 4) + t * (s + 3 * t) -
                     pow(mrho, 2) * (2 * s + 3 * t)) +
                pow(mpion, 2) *
                    (pow(mrho, 6) - 2 * pow(mrho, 4) * (s + 2 * t) +
                     2 * t * (pow(s, 2) - 2 * s * t - 2 * pow(t, 2)) +
                     pow(mrho, 2) *
                         (pow(s, 2) + 2 * s * t + 6 * pow(t, 2)))))) /
             ((-pow(ma1, 2) + t) *
              (pow(Gammaa1, 2) * pow(ma1, 2) +
               pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t,
                   2))) +
         (0.125 * (eta1 - eta2) *
          (eta2 *
               (8 * C4 * pow(mrho, 6) * t - 2 * delta * pow(s, 2) * t +
                pow(mrho, 2) * (-4 * pow(mpion, 4) +
                                (s * (2 + 3 * delta + 8 * C4 * s) - 4 * t) * t +
                                pow(mpion, 2) * (-((-2 + delta) * s) + 8 * t)) +
                pow(mrho, 4) * (8 * C4 * pow(mpion, 4) -
                                pow(mpion, 2) * (-2 + delta + 16 * C4 * t) +
                                t * (-6 + delta + 8 * C4 * (-2 * s + t)))) +
           eta1 *
               (2 * delta * pow(s, 2) * t +
                8 * C4 * pow(mrho, 6) * (-2 * pow(mpion, 2) + t) -
                pow(mrho, 2) * (-4 * pow(mpion, 4) - 4 * pow(t, 2) +
                                2 * pow(mpion, 2) * ((2 + delta) * s + 4 * t) +
                                pow(s, 2) * (-2 + delta + 8 * C4 * t)) +
                pow(mrho, 4) * (-8 * C4 * pow(mpion, 4) + (-2 + delta) * s -
                                4 * t * (1 + 2 * C4 * t) +
                                8 * pow(mpion, 2) * (1 + 2 * C4 * (s + t)))))) /
             (pow(mrho, 2) * (-pow(ma1, 2) + t)) -
         (0.125 * (eta1 - eta2) *
          (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t) *
          (eta1 *
               (pow(mpion, 4) * (4 * pow(mrho, 2) - 8 * C4 * pow(mrho, 4)) +
                8 * C4 * pow(mrho, 6) * (s + t) -
                2 * delta * pow(s, 2) * (s + t) +
                pow(mrho, 2) * ((6 + delta) * pow(s, 2) + 8 * s * t +
                                4 * pow(t, 2) + 8 * C4 * pow(s, 2) * (s + t)) -
                pow(mrho, 4) *
                    (-((-6 + delta) * s) + 4 * t +
                     8 * C4 * (2 * pow(s, 2) + 2 * s * t + pow(t, 2))) +
                2 * pow(mpion, 2) *
                    (-8 * C4 * pow(mrho, 6) + 2 * delta * pow(s, 2) -
                     pow(mrho, 2) * (s * (6 + delta + 8 * C4 * s) + 4 * t) +
                     4 * pow(mrho, 4) * (1 + 2 * C4 * (2 * s + t)))) +
           eta2 *
               (pow(mpion, 4) * (-4 * pow(mrho, 2) + 8 * C4 * pow(mrho, 4)) -
                (-pow(mrho, 2) + s + t) *
                    (16 * C4 * pow(mrho, 6) - 2 * delta * pow(s, 2) +
                     pow(mrho, 2) * (s * (6 + 3 * delta + 8 * C4 * s) + 4 * t) +
                     pow(mrho, 4) * (-10 + delta - 8 * C4 * (3 * s + t))) +
                pow(mpion, 2) *
                    (32 * C4 * pow(mrho, 6) - 4 * delta * pow(s, 2) +
                     pow(mrho, 2) *
                         (s * (14 + 5 * delta + 16 * C4 * s) + 8 * t) +
                     pow(mrho, 4) *
                         (delta - 2 * (9 + 8 * C4 * (3 * s + t))))))) /
             (pow(mrho, 2) *
              (pow(Gammaa1, 2) * pow(ma1, 2) +
               pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t,
                   2))))) /
       (16. * Pi * s * (-4 * pow(mpion, 2) + s)));

  return to_mb * diff_xs / spin_deg_factor;
}

// C22
double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_pi0_rho(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  const double spin_deg_factor = 1.0;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_pion_, m_rho, 0.0);

  const double tmin = t_mandelstam[1];
  const double tmax = t_mandelstam[0];

  const double xs =
      -(pow(Const, 2) * pow(ghat, 4) *
        ((0.03125 * pow(eta1 - 1. * eta2, 2) *
          (eta1 * eta2 *
               (-2. * pow(ma1, 8) - 2. * pow(mpion, 8) +
                2. * pow(mpion, 4) * pow(mrho, 4) +
                pow(ma1, 6) * (8. * pow(mpion, 2) - 4. * s) +
                pow(ma1, 2) * pow(mpion, 2) *
                    (8. * pow(mpion, 4) - 8. * pow(mrho, 4) -
                     4. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s) +
                pow(ma1, 4) * (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                               8. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s -
                               4. * pow(s, 2))) +
           pow(eta2, 2) *
               (1. * pow(ma1, 8) + 1. * pow(mpion, 8) -
                2. * pow(mpion, 6) * pow(mrho, 2) +
                1. * pow(mpion, 4) * pow(mrho, 4) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) + 2. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                               pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                               2. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                pow(ma1, 2) * (-4. * pow(mpion, 6) -
                               2. * pow(mpion, 2) * pow(mrho, 2) * s +
                               pow(mpion, 4) * (6. * pow(mrho, 2) + 2. * s))) +
           pow(eta1, 2) *
               (1. * pow(ma1, 8) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                               pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                               4. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                pow(ma1, 2) * (-4. * pow(mpion, 6) +
                               2. * pow(mpion, 2) * pow(mrho, 2) * s +
                               pow(mrho, 2) * (2. * pow(mrho, 2) - 2. * s) * s +
                               pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s)) +
                pow(mpion, 2) *
                    (1. * pow(mpion, 6) + 2. * pow(mpion, 4) * pow(mrho, 2) -
                     2. * pow(mrho, 6) + 2. * pow(mrho, 4) * s +
                     pow(mpion, 2) *
                         (1. * pow(mrho, 4) - 2. * pow(mrho, 2) * s))))) /
             (1. * pow(ma1, 2) - 1. * tmax) +
         (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
          (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
             (1. * pow(mpion, 2) - 1. * tmax) -
         (0.25 * pow(-2. + delta, 2) * pow(mpion, 2) * tmax) / pow(mrho, 2) -
         0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta2 * (-1. * pow(ma1, 2) + pow(mrho, 2) - 2. * s) +
              eta1 * (2. * pow(mpion, 2) + s)) *
             tmax +
         (0.5 * pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) *
          (4. * pow(mpion, 4) * pow(mrho, 2) + 1. * pow(mrho, 6) -
           3.5 * pow(mrho, 4) * s + 0.5 * pow(s, 3) +
           pow(mpion, 2) * (10. * pow(mrho, 4) - 2. * pow(s, 2))) *
          tmax) /
             (pow(mrho, 6) * pow(pow(mrho, 2) - 1. * s, 2)) -
         (0.25 * (eta1 - 1. * eta2) * (1. * pow(mrho, 2) - 0.5 * delta * s) *
          (eta2 * (-2. * pow(ma1, 4) - 6. * pow(mpion, 4) + 1.5 * pow(mrho, 4) +
                   pow(ma1, 2) *
                       (6. * pow(mpion, 2) - 2. * pow(mrho, 2) - 2. * s) -
                   1. * pow(mrho, 2) * s - 0.5 * pow(s, 2) +
                   pow(mpion, 2) * (2. * pow(mrho, 2) + 2. * s)) +
           eta1 * (2. * pow(ma1, 4) + 6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                   pow(mpion, 2) * (8. * pow(mrho, 2) - 4. * s) -
                   4. * pow(mrho, 2) * s + 1. * pow(s, 2) +
                   pow(ma1, 2) *
                       (-6. * pow(mpion, 2) - 3. * pow(mrho, 2) + 3. * s))) *
          tmax) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
         0.03125 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-6. * pow(ma1, 4) - 12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                   pow(ma1, 2) * (16. * pow(mpion, 2) - 8. * s) +
                   8. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s -
                   4. * pow(s, 2)) +
              pow(eta1, 2) *
                  (3. * pow(ma1, 4) + 6. * pow(mpion, 4) + pow(mrho, 4) +
                   pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                   4. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                   pow(ma1, 2) *
                       (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) +
              pow(eta2, 2) *
                  (3. * pow(ma1, 4) + 6. * pow(mpion, 4) + pow(mrho, 4) +
                   pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                   2. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                   pow(ma1, 2) *
                       (-8. * pow(mpion, 2) + 4. * pow(mrho, 2) + 4. * s))) *
             tmax -
         (1. *
          (pow(mpion, 2) *
               (C4 * (-2. + 1. * delta) * pow(mrho, 6) +
                (1.5 - 2. * delta + 0.625 * pow(delta, 2)) * pow(mrho, 2) * s +
                (0.25 - 0.125 * delta) * delta * pow(s, 2) +
                pow(mrho, 4) * (2.5 - 2.25 * delta + 0.5 * pow(delta, 2) +
                                2. * C4 * s - 1. * C4 * delta * s)) +
           pow(mrho, 2) * (C4 * (-2. + 1. * delta) * pow(mrho, 6) +
                           (0.75 - 0.375 * delta) * delta * pow(s, 2) +
                           pow(mrho, 4) * (0.5 - 0.25 * delta + 6. * C4 * s -
                                           3. * C4 * delta * s) +
                           pow(mrho, 2) * s *
                               (-0.5 - 0.5 * delta + 0.375 * pow(delta, 2) -
                                4. * C4 * s + 2. * C4 * delta * s))) *
          tmax) /
             (pow(mrho, 6) - 1. * pow(mrho, 4) * s) +
         (0.25 * (1. * eta1 - 1. * eta2) *
          (pow(mrho, 2) *
               (eta1 * (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
                        pow(ma1, 2) * (1. - 2. * C4 * pow(mrho, 2)) +
                        pow(mpion, 2) * (-2. + 4. * C4 * pow(mrho, 2)) -
                        2. * C4 * pow(s, 2)) +
                eta2 * (-1.5 * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
                        pow(mpion, 2) * (2. - 4. * C4 * pow(mrho, 2)) +
                        pow(ma1, 2) * (-1. + 2. * C4 * pow(mrho, 2)) + 0.5 * s -
                        4. * C4 * pow(mrho, 2) * s + 2. * C4 * pow(s, 2))) +
           delta * (eta2 * (-1. * pow(ma1, 4) - 3. * pow(mpion, 4) +
                            1. * pow(mrho, 4) +
                            pow(ma1, 2) * (3. * pow(mpion, 2) -
                                           1. * pow(mrho, 2) - 1. * s) +
                            0.25 * pow(mrho, 2) * s - 0.75 * pow(s, 2) +
                            pow(mpion, 2) * (1. * pow(mrho, 2) + 1. * s)) +
                    eta1 * (1. * pow(ma1, 4) + 3. * pow(mpion, 4) +
                            0.5 * pow(mrho, 4) +
                            pow(mpion, 2) * (4. * pow(mrho, 2) - 2. * s) -
                            2. * pow(mrho, 2) * s + 1. * pow(s, 2) +
                            pow(ma1, 2) * (-3. * pow(mpion, 2) -
                                           1.5 * pow(mrho, 2) + 1.5 * s)))) *
          tmax) /
             pow(mrho, 2) +
         (0.5 *
          (pow(delta, 2) *
               (1. * pow(mpion, 4) * pow(mrho, 2) + 0.25 * pow(mrho, 6) -
                0.75 * pow(mrho, 4) * s + 0.125 * pow(mrho, 2) * pow(s, 2) +
                0.25 * pow(s, 3) +
                pow(mpion, 2) * (2.5 * pow(mrho, 4) + 0.25 * pow(mrho, 2) * s -
                                 0.75 * pow(s, 2))) +
           pow(mrho, 6) *
               (1.5 + C4 * (-6. * pow(mrho, 2) + 6. * s) +
                pow(C4, 2) * (4. * pow(mrho, 4) - 8. * pow(mrho, 2) * s +
                              4. * pow(s, 2))) +
           delta * pow(mrho, 2) *
               (4. * C4 * pow(mrho, 6) - 0.5 * pow(s, 2) +
                pow(mrho, 4) * (-1.5 - 3. * C4 * s) +
                pow(mrho, 2) * s * (0.5 - 1. * C4 * s) +
                pow(mpion, 2) * (6. * C4 * pow(mrho, 4) + 0.5 * s +
                                 pow(mrho, 2) * (-2.5 - 2. * C4 * s)))) *
          tmax) /
             pow(mrho, 6) -
         (3. * (1. * pow(mrho, 2) - 0.5 * delta * s) *
          (delta * (0.666667 * pow(mpion, 4) * pow(mrho, 2) +
                    0.166667 * pow(mrho, 6) - 0.541667 * pow(mrho, 4) * s -
                    0.0833333 * pow(mrho, 2) * pow(s, 2) + 0.125 * pow(s, 3) +
                    pow(mpion, 2) *
                        (1.66667 * pow(mrho, 4) + 0.0833333 * pow(mrho, 2) * s -
                         0.416667 * pow(s, 2))) +
           pow(mrho, 2) *
               (1. * C4 * pow(mrho, 6) - 0.0833333 * pow(s, 2) +
                pow(mrho, 4) * (-0.416667 - 1.33333 * C4 * s) +
                pow(mrho, 2) * s * (0.5 + 0.333333 * C4 * s) +
                pow(mpion, 2) *
                    (2. * C4 * pow(mrho, 4) + 0.166667 * s +
                     pow(mrho, 2) * (-0.833333 - 0.666667 * C4 * s)))) *
          tmax) /
             (pow(mrho, 8) - 1. * pow(mrho, 6) * s) -
         1. * C4 * pow(tmax, 2) - 1. * C4 * delta * pow(tmax, 2) +
         0.0625 * (-2. + delta) * (eta1 - 1. * eta2) * eta2 * pow(tmax, 2) -
         (0.5 * pow(delta, 2) * pow(mpion, 2) * pow(tmax, 2)) / pow(mrho, 4) +
         (0.25 * pow(tmax, 2)) / pow(mrho, 2) +
         (0.5 * delta * pow(tmax, 2)) / pow(mrho, 2) -
         (0.25 * pow(delta, 2) * pow(tmax, 2)) / pow(mrho, 2) -
         (0.25 * delta * s * pow(tmax, 2)) / pow(mrho, 4) +
         (0.25 * pow(delta, 2) * s * pow(tmax, 2)) / pow(mrho, 4) +
         (0.5 * C4 * delta * s * pow(tmax, 2)) / pow(mrho, 2) +
         (0.0625 * pow(delta, 2) * pow(s, 2) * pow(tmax, 2)) / pow(mrho, 6) -
         (1. * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s) *
          pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) * pow(tmax, 2)) /
             (pow(mrho, 4) * pow(pow(mrho, 2) - 1. * s, 2)) +
         (0.375 * (eta1 - 1. * eta2) *
          (eta1 * (-0.6666666666666666 * pow(ma1, 2) + 2. * pow(mpion, 2) +
                   1. * pow(mrho, 2) - 1. * s) +
           eta2 *
               (0.6666666666666666 * pow(ma1, 2) - 2. * pow(mpion, 2) +
                0.6666666666666666 * pow(mrho, 2) + 0.6666666666666666 * s)) *
          (1. * pow(mrho, 2) - 0.5 * delta * s) * pow(tmax, 2)) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
         0.03125 * pow(eta1 - 1. * eta2, 3) *
             (eta2 * (-1. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                      1. * pow(mrho, 2) - 1. * s) +
              eta1 *
                  (pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
             pow(tmax, 2) +
         (3. * (1. * pow(mrho, 2) - 0.5 * delta * s) *
          (1. * C4 * pow(mrho, 6) + 0.0833335 * pow(mrho, 2) * s +
           pow(mrho, 4) * (-0.416667 - 0.333334 * C4 * s) +
           delta * (0.666665 * pow(mpion, 2) * pow(mrho, 2) +
                    0.333334 * pow(mrho, 4) - 0.291667 * pow(mrho, 2) * s -
                    0.0416667 * pow(s, 2))) *
          pow(tmax, 2)) /
             (pow(mrho, 8) - 1. * pow(mrho, 6) * s) +
         (0.125 * (1. * eta1 - 1. * eta2) *
          (pow(mrho, 2) * (eta1 * (1. - 2. * C4 * pow(mrho, 2)) +
                           eta2 * (-1. + 2. * C4 * pow(mrho, 2))) +
           delta * (eta2 * (-1. * pow(ma1, 2) + 3. * pow(mpion, 2) -
                            1. * pow(mrho, 2) - 1. * s) +
                    eta1 * (1. * pow(ma1, 2) - 3. * pow(mpion, 2) -
                            1.5 * pow(mrho, 2) + 1.5 * s))) *
          pow(tmax, 2)) /
             pow(mrho, 2) +
         0.0104167 * pow(eta1 - 1. * eta2, 4) * pow(tmax, 3) +
         (0.166667 * pow(delta, 2) * pow(tmax, 3)) / pow(mrho, 4) +
         (0.0833333 * delta * pow(1. * eta1 - 1. * eta2, 2) * pow(tmax, 3)) /
             pow(mrho, 2) +
         (0.666667 * pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) *
          pow(tmax, 3)) /
             (pow(mrho, 4) * pow(pow(mrho, 2) - 1. * s, 2)) -
         (0.166667 * pow(1. * eta1 - 1. * eta2, 2) *
          (1. * pow(mrho, 2) - 0.5 * delta * s) * pow(tmax, 3)) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
         (0.333334 * delta * (-2. * pow(mrho, 2) + 1. * delta * s) *
          pow(tmax, 3)) /
             (pow(mrho, 6) - 1. * pow(mrho, 4) * s) -
         (0.03125 * pow(eta1 - 1. * eta2, 2) *
          (eta1 * eta2 *
               (-2. * pow(ma1, 8) - 2. * pow(mpion, 8) +
                2. * pow(mpion, 4) * pow(mrho, 4) +
                pow(ma1, 6) * (8. * pow(mpion, 2) - 4. * s) +
                pow(ma1, 2) * pow(mpion, 2) *
                    (8. * pow(mpion, 4) - 8. * pow(mrho, 4) -
                     4. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s) +
                pow(ma1, 4) * (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                               8. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s -
                               4. * pow(s, 2))) +
           pow(eta2, 2) *
               (1. * pow(ma1, 8) + 1. * pow(mpion, 8) -
                2. * pow(mpion, 6) * pow(mrho, 2) +
                1. * pow(mpion, 4) * pow(mrho, 4) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) + 2. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                               pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                               2. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                pow(ma1, 2) * (-4. * pow(mpion, 6) -
                               2. * pow(mpion, 2) * pow(mrho, 2) * s +
                               pow(mpion, 4) * (6. * pow(mrho, 2) + 2. * s))) +
           pow(eta1, 2) *
               (1. * pow(ma1, 8) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                               pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                               4. * pow(mrho, 2) * s + 2. * pow(s, 2)) +
                pow(ma1, 2) * (-4. * pow(mpion, 6) +
                               2. * pow(mpion, 2) * pow(mrho, 2) * s +
                               pow(mrho, 2) * (2. * pow(mrho, 2) - 2. * s) * s +
                               pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s)) +
                pow(mpion, 2) *
                    (1. * pow(mpion, 6) + 2. * pow(mpion, 4) * pow(mrho, 2) -
                     2. * pow(mrho, 6) + 2. * pow(mrho, 4) * s +
                     pow(mpion, 2) *
                         (1. * pow(mrho, 4) - 2. * pow(mrho, 2) * s))))) /
             (1. * pow(ma1, 2) - 1. * tmin) -
         (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
          (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
             (1. * pow(mpion, 2) - 1. * tmin) +
         (0.25 * pow(-2. + delta, 2) * pow(mpion, 2) * tmin) / pow(mrho, 2) +
         0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta2 * (-1. * pow(ma1, 2) + pow(mrho, 2) - 2. * s) +
              eta1 * (2. * pow(mpion, 2) + s)) *
             tmin -
         (0.5 * pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) *
          (4. * pow(mpion, 4) * pow(mrho, 2) + 1. * pow(mrho, 6) -
           3.5 * pow(mrho, 4) * s + 0.5 * pow(s, 3) +
           pow(mpion, 2) * (10. * pow(mrho, 4) - 2. * pow(s, 2))) *
          tmin) /
             (pow(mrho, 6) * pow(pow(mrho, 2) - 1. * s, 2)) +
         (0.25 * (eta1 - 1. * eta2) * (1. * pow(mrho, 2) - 0.5 * delta * s) *
          (eta2 * (-2. * pow(ma1, 4) - 6. * pow(mpion, 4) + 1.5 * pow(mrho, 4) +
                   pow(ma1, 2) *
                       (6. * pow(mpion, 2) - 2. * pow(mrho, 2) - 2. * s) -
                   1. * pow(mrho, 2) * s - 0.5 * pow(s, 2) +
                   pow(mpion, 2) * (2. * pow(mrho, 2) + 2. * s)) +
           eta1 * (2. * pow(ma1, 4) + 6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                   pow(mpion, 2) * (8. * pow(mrho, 2) - 4. * s) -
                   4. * pow(mrho, 2) * s + 1. * pow(s, 2) +
                   pow(ma1, 2) *
                       (-6. * pow(mpion, 2) - 3. * pow(mrho, 2) + 3. * s))) *
          tmin) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) -
         0.03125 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-6. * pow(ma1, 4) - 12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                   pow(ma1, 2) * (16. * pow(mpion, 2) - 8. * s) +
                   8. * pow(mpion, 2) * s + 4. * pow(mrho, 2) * s -
                   4. * pow(s, 2)) +
              pow(eta1, 2) *
                  (3. * pow(ma1, 4) + 6. * pow(mpion, 4) + pow(mrho, 4) +
                   pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                   4. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                   pow(ma1, 2) *
                       (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) +
              pow(eta2, 2) *
                  (3. * pow(ma1, 4) + 6. * pow(mpion, 4) + pow(mrho, 4) +
                   pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                   2. * pow(mrho, 2) * s + 2. * pow(s, 2) +
                   pow(ma1, 2) *
                       (-8. * pow(mpion, 2) + 4. * pow(mrho, 2) + 4. * s))) *
             tmin +
         (1. *
          (pow(mpion, 2) *
               (C4 * (-2. + 1. * delta) * pow(mrho, 6) +
                (1.5 - 2. * delta + 0.625 * pow(delta, 2)) * pow(mrho, 2) * s +
                (0.25 - 0.125 * delta) * delta * pow(s, 2) +
                pow(mrho, 4) * (2.5 - 2.25 * delta + 0.5 * pow(delta, 2) +
                                2. * C4 * s - 1. * C4 * delta * s)) +
           pow(mrho, 2) * (C4 * (-2. + 1. * delta) * pow(mrho, 6) +
                           (0.75 - 0.375 * delta) * delta * pow(s, 2) +
                           pow(mrho, 4) * (0.5 - 0.25 * delta + 6. * C4 * s -
                                           3. * C4 * delta * s) +
                           pow(mrho, 2) * s *
                               (-0.5 - 0.5 * delta + 0.375 * pow(delta, 2) -
                                4. * C4 * s + 2. * C4 * delta * s))) *
          tmin) /
             (pow(mrho, 6) - 1. * pow(mrho, 4) * s) -
         (0.25 * (1. * eta1 - 1. * eta2) *
          (pow(mrho, 2) *
               (eta1 * (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
                        pow(ma1, 2) * (1. - 2. * C4 * pow(mrho, 2)) +
                        pow(mpion, 2) * (-2. + 4. * C4 * pow(mrho, 2)) -
                        2. * C4 * pow(s, 2)) +
                eta2 * (-1.5 * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
                        pow(mpion, 2) * (2. - 4. * C4 * pow(mrho, 2)) +
                        pow(ma1, 2) * (-1. + 2. * C4 * pow(mrho, 2)) + 0.5 * s -
                        4. * C4 * pow(mrho, 2) * s + 2. * C4 * pow(s, 2))) +
           delta * (eta2 * (-1. * pow(ma1, 4) - 3. * pow(mpion, 4) +
                            1. * pow(mrho, 4) +
                            pow(ma1, 2) * (3. * pow(mpion, 2) -
                                           1. * pow(mrho, 2) - 1. * s) +
                            0.25 * pow(mrho, 2) * s - 0.75 * pow(s, 2) +
                            pow(mpion, 2) * (1. * pow(mrho, 2) + 1. * s)) +
                    eta1 * (1. * pow(ma1, 4) + 3. * pow(mpion, 4) +
                            0.5 * pow(mrho, 4) +
                            pow(mpion, 2) * (4. * pow(mrho, 2) - 2. * s) -
                            2. * pow(mrho, 2) * s + 1. * pow(s, 2) +
                            pow(ma1, 2) * (-3. * pow(mpion, 2) -
                                           1.5 * pow(mrho, 2) + 1.5 * s)))) *
          tmin) /
             pow(mrho, 2) -
         (0.5 *
          (pow(delta, 2) *
               (1. * pow(mpion, 4) * pow(mrho, 2) + 0.25 * pow(mrho, 6) -
                0.75 * pow(mrho, 4) * s + 0.125 * pow(mrho, 2) * pow(s, 2) +
                0.25 * pow(s, 3) +
                pow(mpion, 2) * (2.5 * pow(mrho, 4) + 0.25 * pow(mrho, 2) * s -
                                 0.75 * pow(s, 2))) +
           pow(mrho, 6) *
               (1.5 + C4 * (-6. * pow(mrho, 2) + 6. * s) +
                pow(C4, 2) * (4. * pow(mrho, 4) - 8. * pow(mrho, 2) * s +
                              4. * pow(s, 2))) +
           delta * pow(mrho, 2) *
               (4. * C4 * pow(mrho, 6) - 0.5 * pow(s, 2) +
                pow(mrho, 4) * (-1.5 - 3. * C4 * s) +
                pow(mrho, 2) * s * (0.5 - 1. * C4 * s) +
                pow(mpion, 2) * (6. * C4 * pow(mrho, 4) + 0.5 * s +
                                 pow(mrho, 2) * (-2.5 - 2. * C4 * s)))) *
          tmin) /
             pow(mrho, 6) +
         (3. * (1. * pow(mrho, 2) - 0.5 * delta * s) *
          (delta * (0.666667 * pow(mpion, 4) * pow(mrho, 2) +
                    0.166667 * pow(mrho, 6) - 0.541667 * pow(mrho, 4) * s -
                    0.0833333 * pow(mrho, 2) * pow(s, 2) + 0.125 * pow(s, 3) +
                    pow(mpion, 2) *
                        (1.66667 * pow(mrho, 4) + 0.0833333 * pow(mrho, 2) * s -
                         0.416667 * pow(s, 2))) +
           pow(mrho, 2) *
               (1. * C4 * pow(mrho, 6) - 0.0833333 * pow(s, 2) +
                pow(mrho, 4) * (-0.416667 - 1.33333 * C4 * s) +
                pow(mrho, 2) * s * (0.5 + 0.333333 * C4 * s) +
                pow(mpion, 2) *
                    (2. * C4 * pow(mrho, 4) + 0.166667 * s +
                     pow(mrho, 2) * (-0.833333 - 0.666667 * C4 * s)))) *
          tmin) /
             (pow(mrho, 8) - 1. * pow(mrho, 6) * s) +
         1. * C4 * pow(tmin, 2) + 1. * C4 * delta * pow(tmin, 2) -
         0.0625 * (-2. + delta) * (eta1 - 1. * eta2) * eta2 * pow(tmin, 2) +
         (0.5 * pow(delta, 2) * pow(mpion, 2) * pow(tmin, 2)) / pow(mrho, 4) -
         (0.25 * pow(tmin, 2)) / pow(mrho, 2) -
         (0.5 * delta * pow(tmin, 2)) / pow(mrho, 2) +
         (0.25 * pow(delta, 2) * pow(tmin, 2)) / pow(mrho, 2) +
         (0.25 * delta * s * pow(tmin, 2)) / pow(mrho, 4) -
         (0.25 * pow(delta, 2) * s * pow(tmin, 2)) / pow(mrho, 4) -
         (0.5 * C4 * delta * s * pow(tmin, 2)) / pow(mrho, 2) -
         (0.0625 * pow(delta, 2) * pow(s, 2) * pow(tmin, 2)) / pow(mrho, 6) +
         (1. * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s) *
          pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) * pow(tmin, 2)) /
             (pow(mrho, 4) * pow(pow(mrho, 2) - 1. * s, 2)) -
         (0.375 * (eta1 - 1. * eta2) *
          (eta1 * (-0.6666666666666666 * pow(ma1, 2) + 2. * pow(mpion, 2) +
                   1. * pow(mrho, 2) - 1. * s) +
           eta2 *
               (0.6666666666666666 * pow(ma1, 2) - 2. * pow(mpion, 2) +
                0.6666666666666666 * pow(mrho, 2) + 0.6666666666666666 * s)) *
          (1. * pow(mrho, 2) - 0.5 * delta * s) * pow(tmin, 2)) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) -
         0.03125 * pow(eta1 - 1. * eta2, 3) *
             (eta2 * (-1. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                      1. * pow(mrho, 2) - 1. * s) +
              eta1 *
                  (pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
             pow(tmin, 2) -
         (3. * (1. * pow(mrho, 2) - 0.5 * delta * s) *
          (1. * C4 * pow(mrho, 6) + 0.0833335 * pow(mrho, 2) * s +
           pow(mrho, 4) * (-0.416667 - 0.333334 * C4 * s) +
           delta * (0.666665 * pow(mpion, 2) * pow(mrho, 2) +
                    0.333334 * pow(mrho, 4) - 0.291667 * pow(mrho, 2) * s -
                    0.0416667 * pow(s, 2))) *
          pow(tmin, 2)) /
             (pow(mrho, 8) - 1. * pow(mrho, 6) * s) -
         (0.125 * (1. * eta1 - 1. * eta2) *
          (pow(mrho, 2) * (eta1 * (1. - 2. * C4 * pow(mrho, 2)) +
                           eta2 * (-1. + 2. * C4 * pow(mrho, 2))) +
           delta * (eta2 * (-1. * pow(ma1, 2) + 3. * pow(mpion, 2) -
                            1. * pow(mrho, 2) - 1. * s) +
                    eta1 * (1. * pow(ma1, 2) - 3. * pow(mpion, 2) -
                            1.5 * pow(mrho, 2) + 1.5 * s))) *
          pow(tmin, 2)) /
             pow(mrho, 2) -
         0.0104167 * pow(eta1 - 1. * eta2, 4) * pow(tmin, 3) -
         (0.166667 * pow(delta, 2) * pow(tmin, 3)) / pow(mrho, 4) -
         (0.0833333 * delta * pow(1. * eta1 - 1. * eta2, 2) * pow(tmin, 3)) /
             pow(mrho, 2) -
         (0.666667 * pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) *
          pow(tmin, 3)) /
             (pow(mrho, 4) * pow(pow(mrho, 2) - 1. * s, 2)) +
         (0.166667 * pow(1. * eta1 - 1. * eta2, 2) *
          (1. * pow(mrho, 2) - 0.5 * delta * s) * pow(tmin, 3)) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) -
         (0.333334 * delta * (-2. * pow(mrho, 2) + 1. * delta * s) *
          pow(tmin, 3)) /
             (pow(mrho, 6) - 1. * pow(mrho, 4) * s) +
         0.0625 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-4. * pow(ma1, 6) +
                   pow(ma1, 4) * (12. * pow(mpion, 2) - 6. * s) +
                   pow(mpion, 2) *
                       (4. * pow(mpion, 4) - 4. * pow(mrho, 4) -
                        2. * pow(mpion, 2) * s + 2. * pow(mrho, 2) * s) +
                   pow(ma1, 2) * (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                                  8. * pow(mpion, 2) * s +
                                  4. * pow(mrho, 2) * s - 4. * pow(s, 2))) +
              pow(eta1, 2) *
                  (2. * pow(ma1, 6) - 2. * pow(mpion, 6) +
                   pow(mpion, 2) * pow(mrho, 2) * s +
                   pow(mrho, 2) * (pow(mrho, 2) - 1. * s) * s +
                   pow(mpion, 4) * (-3. * pow(mrho, 2) + s) +
                   pow(ma1, 4) *
                       (-6. * pow(mpion, 2) - 3. * pow(mrho, 2) + 3. * s) +
                   pow(ma1, 2) * (6. * pow(mpion, 4) + pow(mrho, 4) +
                                  pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                                  4. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
              pow(eta2, 2) *
                  (2. * pow(ma1, 6) - 2. * pow(mpion, 6) -
                   1. * pow(mpion, 2) * pow(mrho, 2) * s +
                   pow(mpion, 4) * (3. * pow(mrho, 2) + s) +
                   pow(ma1, 4) *
                       (-6. * pow(mpion, 2) + 3. * pow(mrho, 2) + 3. * s) +
                   pow(ma1, 2) *
                       (6. * pow(mpion, 4) + pow(mrho, 4) +
                        pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                        2. * pow(mrho, 2) * s + 2. * pow(s, 2)))) *
             log(fabs(-pow(ma1, 2) + tmax)) -
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * (-0.5 * pow(ma1, 6) - 0.5 * pow(mpion, 6) +
                   0.5 * pow(mpion, 4) * pow(mrho, 2) +
                   pow(ma1, 4) *
                       (0.5 * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 1. * s) +
                   pow(ma1, 2) * pow(mpion, 2) *
                       (0.5 * pow(mpion, 2) + 1. * pow(mrho, 2) - 1. * s)) +
           eta1 * (pow(ma1, 4) * (1. * pow(mpion, 2) + 0.5 * s) +
                   pow(mpion, 2) *
                       (1. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                        pow(mpion, 2) * (-1. * pow(mrho, 2) + 0.5 * s) -
                        0.5 * pow(mrho, 2) * s) +
                   pow(ma1, 2) *
                       (-2. * pow(mpion, 4) - 0.5 * pow(mrho, 2) * s +
                        pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s)))) *
          log(fabs(-pow(ma1, 2) + tmax))) /
             (pow(ma1, 2) - 1. * pow(mpion, 2)) -
         (0.125 * (eta1 - 1. * eta2) * (1. * pow(mrho, 2) - 0.5 * delta * s) *
          (eta1 *
               (4. * pow(ma1, 6) - 4. * pow(mpion, 6) + pow(mrho, 4) * s +
                4. * pow(mpion, 2) * pow(s, 2) - 1. * pow(s, 3) +
                pow(mpion, 4) * (-10. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) *
                    (-12. * pow(mpion, 2) - 6. * pow(mrho, 2) + 6. * s) +
                pow(ma1, 2) * (12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                               pow(mpion, 2) * (16. * pow(mrho, 2) - 8. * s) -
                               8. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
           eta2 *
               (-4. * pow(ma1, 6) +
                pow(ma1, 4) *
                    (12. * pow(mpion, 2) - 4. * pow(mrho, 2) - 4. * s) +
                pow(mpion, 2) * (4. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                                 2. * pow(mrho, 2) * s - 1. * pow(s, 2)) +
                pow(ma1, 2) * (-12. * pow(mpion, 4) + 3. * pow(mrho, 4) -
                               2. * pow(mrho, 2) * s - 1. * pow(s, 2) +
                               pow(mpion, 2) * (4. * pow(mrho, 2) + 4. * s)))) *
          log(fabs(-pow(ma1, 2) + tmax))) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
         (0.25 * (1. * eta1 - 1. * eta2) *
          (delta *
               (eta1 * (1. * pow(ma1, 6) - 1. * pow(mpion, 6) +
                        pow(mpion, 4) * (-2.5 * pow(mrho, 2) + 0.5 * s) +
                        pow(mpion, 2) * s * (-0.5 * pow(mrho, 2) + 1. * s) +
                        pow(ma1, 4) * (-3. * pow(mpion, 2) -
                                       1.5 * pow(mrho, 2) + 1.5 * s) +
                        s * (0.5 * pow(mrho, 4) - 0.25 * pow(mrho, 2) * s -
                             0.25 * pow(s, 2)) +
                        pow(ma1, 2) *
                            (3. * pow(mpion, 4) + 0.5 * pow(mrho, 4) +
                             pow(mpion, 2) * (4. * pow(mrho, 2) - 2. * s) -
                             2. * pow(mrho, 2) * s + 1. * pow(s, 2))) +
                eta2 * (-1. * pow(ma1, 6) +
                        pow(ma1, 4) *
                            (3. * pow(mpion, 2) - 1. * pow(mrho, 2) - 1. * s) +
                        pow(mpion, 2) *
                            (1. * pow(mpion, 4) - 0.5 * pow(mrho, 4) +
                             0.25 * pow(mrho, 2) * s - 0.25 * pow(s, 2)) +
                        pow(ma1, 2) *
                            (-3. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                             0.25 * pow(mrho, 2) * s - 0.75 * pow(s, 2) +
                             pow(mpion, 2) * (1. * pow(mrho, 2) + 1. * s)))) +
           pow(mrho, 2) *
               (eta2 * (pow(ma1, 4) * (-1. + 2. * C4 * pow(mrho, 2)) +
                        pow(mpion, 2) *
                            (0.5 * pow(mrho, 2) +
                             pow(mpion, 2) * (-1. + 2. * C4 * pow(mrho, 2)) +
                             0.5 * s) +
                        pow(ma1, 2) *
                            (2. * C4 * pow(mrho, 4) +
                             pow(mpion, 2) * (2. - 4. * C4 * pow(mrho, 2)) +
                             pow(mrho, 2) * (-1.5 - 4. * C4 * s) +
                             s * (0.5 + 2. * C4 * s))) +
                eta1 * (pow(ma1, 4) * (1. - 2. * C4 * pow(mrho, 2)) +
                        pow(mpion, 4) * (1. - 2. * C4 * pow(mrho, 2)) +
                        (-0.5 * pow(mrho, 2) + 0.5 * s) * s +
                        pow(ma1, 2) *
                            (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
                             pow(mpion, 2) * (-2. + 4. * C4 * pow(mrho, 2)) -
                             2. * C4 * pow(s, 2)) +
                        pow(mpion, 2) * (-4. * C4 * pow(mrho, 4) - 1. * s +
                                         pow(mrho, 2) * (2. + 4. * C4 * s))))) *
          log(fabs(-pow(ma1, 2) + tmax))) /
             pow(mrho, 2) +
         0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
             log(fabs(-pow(mpion, 2) + tmax)) +
         (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
          (eta2 * pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s) +
           eta1 * (-0.5 * pow(mrho, 4) +
                   pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                   0.5 * pow(mrho, 2) * s)) *
          log(fabs(-pow(mpion, 2) + tmax))) /
             (-1. * pow(ma1, 2) + 1. * pow(mpion, 2)) -
         (2. *
          (-0.12500000000000003 * pow(2. - 1. * delta, 2) * pow(mrho, 4) * s +
           pow(mpion, 2) * (C4 * (-2. + 1. * delta) * pow(mrho, 6) +
                            (0.5 - 0.25 * delta) * delta * pow(s, 2) +
                            pow(mrho, 4) * (1. - 0.5 * delta + 4. * C4 * s -
                                            2. * C4 * delta * s) +
                            pow(mrho, 2) * s *
                                (1. - 2. * delta + 0.75 * pow(delta, 2) -
                                 2. * C4 * s + 1. * C4 * delta * s))) *
          log(fabs(-pow(mpion, 2) + tmax))) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) -
         0.0625 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-4. * pow(ma1, 6) +
                   pow(ma1, 4) * (12. * pow(mpion, 2) - 6. * s) +
                   pow(mpion, 2) *
                       (4. * pow(mpion, 4) - 4. * pow(mrho, 4) -
                        2. * pow(mpion, 2) * s + 2. * pow(mrho, 2) * s) +
                   pow(ma1, 2) * (-12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                                  8. * pow(mpion, 2) * s +
                                  4. * pow(mrho, 2) * s - 4. * pow(s, 2))) +
              pow(eta1, 2) *
                  (2. * pow(ma1, 6) - 2. * pow(mpion, 6) +
                   pow(mpion, 2) * pow(mrho, 2) * s +
                   pow(mrho, 2) * (pow(mrho, 2) - 1. * s) * s +
                   pow(mpion, 4) * (-3. * pow(mrho, 2) + s) +
                   pow(ma1, 4) *
                       (-6. * pow(mpion, 2) - 3. * pow(mrho, 2) + 3. * s) +
                   pow(ma1, 2) * (6. * pow(mpion, 4) + pow(mrho, 4) +
                                  pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                                  4. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
              pow(eta2, 2) *
                  (2. * pow(ma1, 6) - 2. * pow(mpion, 6) -
                   1. * pow(mpion, 2) * pow(mrho, 2) * s +
                   pow(mpion, 4) * (3. * pow(mrho, 2) + s) +
                   pow(ma1, 4) *
                       (-6. * pow(mpion, 2) + 3. * pow(mrho, 2) + 3. * s) +
                   pow(ma1, 2) *
                       (6. * pow(mpion, 4) + pow(mrho, 4) +
                        pow(mpion, 2) * (-6. * pow(mrho, 2) - 4. * s) -
                        2. * pow(mrho, 2) * s + 2. * pow(s, 2)))) *
             log(fabs(-pow(ma1, 2) + tmin)) +
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * (-0.5 * pow(ma1, 6) - 0.5 * pow(mpion, 6) +
                   0.5 * pow(mpion, 4) * pow(mrho, 2) +
                   pow(ma1, 4) *
                       (0.5 * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 1. * s) +
                   pow(ma1, 2) * pow(mpion, 2) *
                       (0.5 * pow(mpion, 2) + 1. * pow(mrho, 2) - 1. * s)) +
           eta1 * (pow(ma1, 4) * (1. * pow(mpion, 2) + 0.5 * s) +
                   pow(mpion, 2) *
                       (1. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                        pow(mpion, 2) * (-1. * pow(mrho, 2) + 0.5 * s) -
                        0.5 * pow(mrho, 2) * s) +
                   pow(ma1, 2) *
                       (-2. * pow(mpion, 4) - 0.5 * pow(mrho, 2) * s +
                        pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s)))) *
          log(fabs(-pow(ma1, 2) + tmin))) /
             (pow(ma1, 2) - 1. * pow(mpion, 2)) +
         (0.125 * (eta1 - 1. * eta2) * (1. * pow(mrho, 2) - 0.5 * delta * s) *
          (eta1 *
               (4. * pow(ma1, 6) - 4. * pow(mpion, 6) + pow(mrho, 4) * s +
                4. * pow(mpion, 2) * pow(s, 2) - 1. * pow(s, 3) +
                pow(mpion, 4) * (-10. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) *
                    (-12. * pow(mpion, 2) - 6. * pow(mrho, 2) + 6. * s) +
                pow(ma1, 2) * (12. * pow(mpion, 4) + 2. * pow(mrho, 4) +
                               pow(mpion, 2) * (16. * pow(mrho, 2) - 8. * s) -
                               8. * pow(mrho, 2) * s + 2. * pow(s, 2))) +
           eta2 *
               (-4. * pow(ma1, 6) +
                pow(ma1, 4) *
                    (12. * pow(mpion, 2) - 4. * pow(mrho, 2) - 4. * s) +
                pow(mpion, 2) * (4. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                                 2. * pow(mrho, 2) * s - 1. * pow(s, 2)) +
                pow(ma1, 2) * (-12. * pow(mpion, 4) + 3. * pow(mrho, 4) -
                               2. * pow(mrho, 2) * s - 1. * pow(s, 2) +
                               pow(mpion, 2) * (4. * pow(mrho, 2) + 4. * s)))) *
          log(fabs(-pow(ma1, 2) + tmin))) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) -
         (0.25 * (1. * eta1 - 1. * eta2) *
          (delta *
               (eta1 * (1. * pow(ma1, 6) - 1. * pow(mpion, 6) +
                        pow(mpion, 4) * (-2.5 * pow(mrho, 2) + 0.5 * s) +
                        pow(mpion, 2) * s * (-0.5 * pow(mrho, 2) + 1. * s) +
                        pow(ma1, 4) * (-3. * pow(mpion, 2) -
                                       1.5 * pow(mrho, 2) + 1.5 * s) +
                        s * (0.5 * pow(mrho, 4) - 0.25 * pow(mrho, 2) * s -
                             0.25 * pow(s, 2)) +
                        pow(ma1, 2) *
                            (3. * pow(mpion, 4) + 0.5 * pow(mrho, 4) +
                             pow(mpion, 2) * (4. * pow(mrho, 2) - 2. * s) -
                             2. * pow(mrho, 2) * s + 1. * pow(s, 2))) +
                eta2 * (-1. * pow(ma1, 6) +
                        pow(ma1, 4) *
                            (3. * pow(mpion, 2) - 1. * pow(mrho, 2) - 1. * s) +
                        pow(mpion, 2) *
                            (1. * pow(mpion, 4) - 0.5 * pow(mrho, 4) +
                             0.25 * pow(mrho, 2) * s - 0.25 * pow(s, 2)) +
                        pow(ma1, 2) *
                            (-3. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                             0.25 * pow(mrho, 2) * s - 0.75 * pow(s, 2) +
                             pow(mpion, 2) * (1. * pow(mrho, 2) + 1. * s)))) +
           pow(mrho, 2) *
               (eta2 * (pow(ma1, 4) * (-1. + 2. * C4 * pow(mrho, 2)) +
                        pow(mpion, 2) *
                            (0.5 * pow(mrho, 2) +
                             pow(mpion, 2) * (-1. + 2. * C4 * pow(mrho, 2)) +
                             0.5 * s) +
                        pow(ma1, 2) *
                            (2. * C4 * pow(mrho, 4) +
                             pow(mpion, 2) * (2. - 4. * C4 * pow(mrho, 2)) +
                             pow(mrho, 2) * (-1.5 - 4. * C4 * s) +
                             s * (0.5 + 2. * C4 * s))) +
                eta1 * (pow(ma1, 4) * (1. - 2. * C4 * pow(mrho, 2)) +
                        pow(mpion, 4) * (1. - 2. * C4 * pow(mrho, 2)) +
                        (-0.5 * pow(mrho, 2) + 0.5 * s) * s +
                        pow(ma1, 2) *
                            (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
                             pow(mpion, 2) * (-2. + 4. * C4 * pow(mrho, 2)) -
                             2. * C4 * pow(s, 2)) +
                        pow(mpion, 2) * (-4. * C4 * pow(mrho, 4) - 1. * s +
                                         pow(mrho, 2) * (2. + 4. * C4 * s))))) *
          log(fabs(-pow(ma1, 2) + tmin))) /
             pow(mrho, 2) -
         0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
             log(fabs(-pow(mpion, 2) + tmin)) -
         (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
          (eta2 * pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s) +
           eta1 * (-0.5 * pow(mrho, 4) +
                   pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                   0.5 * pow(mrho, 2) * s)) *
          log(fabs(-pow(mpion, 2) + tmin))) /
             (-1. * pow(ma1, 2) + 1. * pow(mpion, 2)) +
         (2. *
          (-0.12500000000000003 * pow(2. - 1. * delta, 2) * pow(mrho, 4) * s +
           pow(mpion, 2) * (C4 * (-2. + 1. * delta) * pow(mrho, 6) +
                            (0.5 - 0.25 * delta) * delta * pow(s, 2) +
                            pow(mrho, 4) * (1. - 0.5 * delta + 4. * C4 * s -
                                            2. * C4 * delta * s) +
                            pow(mrho, 2) * s *
                                (1. - 2. * delta + 0.75 * pow(delta, 2) -
                                 2. * C4 * s + 1. * C4 * delta * s))) *
          log(fabs(-pow(mpion, 2) + tmin))) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s))) /
      (16. * Pi * (4 * pow(mpion, 2) - s) * s);

  return xs * to_mb / spin_deg_factor;
}

double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_pi0_rho(
    const double s, const double t, const double m_rho) {
  const double &mpion = m_pion_;
  const double &mrho = m_rho;
  const double spin_deg_factor = 1.0;

  const double diff_xs =
      (pow(Const, 2) * pow(ghat, 4) *
       ((-0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (pow(mpion, 4) + pow(pow(mrho, 2) - t, 2) -
          2 * pow(mpion, 2) * (pow(mrho, 2) + t))) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - t, 2)) +
        (0.0625 * pow(-2 * pow(mrho, 2) + delta * s, 2) *
         (8 * pow(mpion, 4) * pow(mrho, 2) + 2 * pow(mrho, 6) + pow(s, 3) +
          8 * pow(mrho, 2) * t * (s + t) - pow(mrho, 4) * (7 * s + 8 * t) +
          4 * pow(mpion, 2) *
              (5 * pow(mrho, 4) - pow(s, 2) - 4 * pow(mrho, 2) * t))) /
            (pow(mrho, 6) * pow(pow(mrho, 2) - s, 2)) -
        (0.0625 * (eta1 - eta2) * (2 * pow(mrho, 2) - delta * s) *
         (-(eta2 * (2 * pow(mpion, 2) + pow(mrho, 2) - s - 2 * t) *
            (2 * pow(mpion, 4) + pow(mpion, 2) * (-pow(mrho, 2) + s - 4 * t) +
             t * (3 * pow(mrho, 2) + s + 2 * t))) +
          eta1 * (4 * pow(mpion, 6) - pow(mrho, 4) * s + pow(s, 3) +
                  2 * pow(mpion, 4) * (5 * pow(mrho, 2) - s - 6 * t) -
                  2 * (pow(mrho, 4) - 4 * pow(mrho, 2) * s + pow(s, 2)) * t +
                  6 * (pow(mrho, 2) - s) * pow(t, 2) - 4 * pow(t, 3) -
                  4 * pow(mpion, 2) *
                      (4 * pow(mrho, 2) * t + (s - 3 * t) * (s + t))))) /
            (pow(mrho, 2) * (pow(mrho, 2) - s) * (pow(ma1, 2) - t)) -
        (0.125 * (-2 + delta) * (eta1 - eta2) *
         (-(eta2 * (pow(mpion, 2) + t) *
            (pow(mpion, 4) + t * (-pow(mrho, 2) + 2 * s + t) -
             pow(mpion, 2) * (pow(mrho, 2) + 2 * t))) +
          eta1 * (2 * pow(mpion, 6) +
                  pow(mpion, 4) * (-2 * pow(mrho, 2) + s - 4 * t) +
                  s * t * (-pow(mrho, 2) + t) +
                  pow(mpion, 2) * (2 * pow(mrho, 4) + 2 * t * (s + t) -
                                   pow(mrho, 2) * (s + 2 * t))))) /
            ((-pow(ma1, 2) + t) * (-pow(mpion, 2) + t)) +
        (0.03125 * pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) - 4 * pow(mpion, 6) * t +
               pow(t, 2) * (-pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                            2 * pow(s, 2) + 2 * s * t + pow(t, 2)) -
               2 * pow(mpion, 2) * t *
                   (-2 * pow(mrho, 4) + pow(mrho, 2) * s + 2 * t * (s + t)) +
               pow(mpion, 4) * (-pow(mrho, 4) + 2 * t * (s + 3 * t))) +
          pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * (pow(mrho, 2) + 2 * t) +
               pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t +
                            pow(t, 2) + 2 * pow(mrho, 2) * (-s + t)) -
               2 * pow(mpion, 2) * t *
                   (2 * t * (s + t) + pow(mrho, 2) * (s + 3 * t)) +
               pow(mpion, 4) * (pow(mrho, 4) + 6 * pow(mrho, 2) * t +
                                2 * t * (s + 3 * t))) +
          pow(eta1, 2) *
              (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) -
               2 * pow(mpion, 2) * (pow(mrho, 2) - s - t) *
                   (pow(mrho, 4) + pow(mrho, 2) * t - 2 * pow(t, 2)) +
               t * (-pow(mrho, 2) + t) *
                   (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                    pow(mrho, 2) * (2 * s + t)) +
               pow(mpion, 4) * (pow(mrho, 4) - 2 * pow(mrho, 2) * (s + 3 * t) +
                                2 * t * (s + 3 * t))))) /
            pow(pow(ma1, 2) - t, 2) -
        (3. * (1. * pow(mrho, 2) - 0.5 * delta * s) *
         (delta *
              (0.666667 * pow(mpion, 4) * pow(mrho, 2) +
               0.166667 * pow(mrho, 6) +
               pow(mpion, 2) * (1.66667 * pow(mrho, 4) - 0.416667 * pow(s, 2) +
                                pow(mrho, 2) * (0.0833333 * s - 1.33333 * t)) +
               pow(mrho, 4) * (-0.541667 * s - 0.666667 * t) +
               pow(s, 2) * (0.125 * s + 0.0833333 * t) +
               pow(mrho, 2) * (-0.0833333 * pow(s, 2) + 0.583333 * s * t +
                               0.666667 * pow(t, 2))) +
          pow(mrho, 2) *
              (1. * C4 * pow(mrho, 6) +
               pow(mpion, 2) *
                   (2. * C4 * pow(mrho, 4) + 0.166667 * s +
                    pow(mrho, 2) * (-0.833333 - 0.666667 * C4 * s)) +
               s * (-0.0833333 * s - 0.166667 * t) +
               pow(mrho, 4) * (-0.416667 - 1.33333 * C4 * s - 2. * C4 * t) +
               pow(mrho, 2) * (0.833333 * t + s * (0.5 + 0.333333 * C4 * s +
                                                   0.666667 * C4 * t))))) /
            (pow(mrho, 8) - 1. * pow(mrho, 6) * s) +
        (pow(mrho, 6) * (0.75 + C4 * (2. * C4 * pow(mrho, 4) +
                                      pow(mrho, 2) * (-3. - 4. * C4 * s) +
                                      s * (3. + 2. * C4 * s))) +
         pow(delta, 2) *
             (0.5 * pow(mpion, 4) * pow(mrho, 2) + 0.125 * pow(mrho, 6) +
              pow(mpion, 2) * (1.25 * pow(mrho, 4) - 0.375 * pow(s, 2) +
                               pow(mrho, 2) * (0.125 * s - 1. * t)) +
              pow(mrho, 4) * (-0.375 * s - 0.5 * t) +
              pow(s, 2) * (0.125 * s + 0.125 * t) +
              pow(mrho, 2) *
                  (0.0625 * pow(s, 2) + 0.375 * s * t + 0.5 * pow(t, 2))) +
         delta * pow(mrho, 2) *
             (2. * C4 * pow(mrho, 6) +
              pow(mpion, 2) * (3. * C4 * pow(mrho, 4) + 0.25 * s +
                               pow(mrho, 2) * (-1.25 - 1. * C4 * s)) +
              s * (-0.25 * s - 0.25 * t) +
              pow(mrho, 4) * (-0.75 - 1.5 * C4 * s - 3. * C4 * t) +
              pow(mrho, 2) *
                  (1.25 * t + s * (0.25 - 0.5 * C4 * s + 1. * C4 * t)))) /
            pow(mrho, 6) +
        (2 *
         ((-0.0625 * (-2. + delta) * (-2. * pow(mrho, 2) + delta * s) *
           (pow(mpion, 4) * (4. * pow(mrho, 2) + 4. * s) +
            pow(mpion, 2) *
                (pow(mrho, 2) * (-7. * s - 4. * t) + s * (-1. * s - 4. * t)) +
            s * (pow(mrho, 4) + pow(mrho, 2) * (s - 1. * t) + s * t))) /
              ((pow(mrho, 2) - 1. * s) * (pow(mpion, 2) - 1. * t)) +
          (0.0625 * (-2 + delta) *
           (pow(mpion, 4) * ((-2 + 4 * delta) * pow(mrho, 2) +
                             8 * C4 * pow(mrho, 4) + 5 * delta * s) -
            8 * C4 * pow(mrho, 6) * t + delta * s * t * (s + t) +
            pow(mrho, 2) * (delta * s * (s - 3 * t) - 2 * t * (s + t)) +
            2 * pow(mrho, 4) *
                ((-1 + delta) * s + t + 4 * C4 * t * (2 * s + t)) -
            pow(mpion, 2) * (8 * C4 * pow(mrho, 6) + delta * s * (s + 6 * t) +
                             2 * pow(mrho, 4) * (-3 + 8 * C4 * t) +
                             pow(mrho, 2) * ((-2 + 9 * delta) * s +
                                             4 * (-1 + delta) * t)))) /
              (-pow(mpion, 2) + t))) /
            pow(mrho, 4) -
        (0.0625 * (eta1 - eta2) *
         (eta2 *
              (-4 * delta * pow(mpion, 6) +
               4 * pow(mpion, 4) *
                   (pow(mrho, 2) - 2 * C4 * pow(mrho, 4) + 3 * delta * t) +
               pow(mpion, 2) * (delta * (s - 6 * t) * (s + 2 * t) -
                                (2 + delta) * pow(mrho, 2) * (s + 4 * t) +
                                2 * pow(mrho, 4) * (-1 + delta + 8 * C4 * t)) +
               t * (-8 * C4 * pow(mrho, 6) +
                    pow(mrho, 4) * (6 - 4 * delta + 16 * C4 * s - 8 * C4 * t) +
                    pow(mrho, 2) * (-(s * (2 + delta + 8 * C4 * s)) +
                                    4 * (1 + delta) * t) +
                    delta * (3 * pow(s, 2) + 4 * s * t + 4 * pow(t, 2)))) +
          eta1 *
              (4 * delta * pow(mpion, 6) - 8 * C4 * pow(mrho, 6) * t +
               delta * (pow(s, 3) - 4 * pow(s, 2) * t - 6 * s * pow(t, 2) -
                        4 * pow(t, 3)) -
               2 * pow(mpion, 4) *
                   ((2 - 5 * delta) * pow(mrho, 2) - 4 * C4 * pow(mrho, 4) +
                    delta * (s + 6 * t)) +
               2 * pow(mrho, 4) *
                   (s - delta * s + t * (2 - delta + 4 * C4 * t)) +
               pow(mrho, 2) *
                   (8 * delta * s * t + 2 * (-2 + 3 * delta) * pow(t, 2) +
                    pow(s, 2) * (-2 + delta + 8 * C4 * t)) -
               2 * pow(mpion, 2) *
                   (-8 * C4 * pow(mrho, 6) + 2 * delta * (s - 3 * t) * (s + t) -
                    pow(mrho, 2) * ((2 + delta) * s + (4 - 8 * delta) * t) +
                    pow(mrho, 4) * (4 + 8 * C4 * (s + t)))))) /
            (pow(mrho, 2) * (-pow(ma1, 2) + t)))) /
      (16. * Pi * s * (-4 * pow(mpion, 2) + s));

  return to_mb * diff_xs / spin_deg_factor;
}

/*----------------------------------------------------------------------------*/
/*																																						*/
/* 																Tabulation
 */
/*																																						*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*				 Pi + Rho -> Pi + Photon channels mediated by
 * (Pi, Rho, a1) 				*/
/*----------------------------------------------------------------------------*/

// C11

double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi_rho0_pi(
    const double s, const double m_rho) {
  if (tab_pi_rho0_pi_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi_rho0_pi\n");
    tab_pi_rho0_pi_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        xs_wrapper<
            PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_rho0_pi>);
  }

  return tab_pi_rho0_pi_->get_linear(s, m_rho);
}
double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi_rho0_pi(
    const double s, const double t, const double m_rho) {
  if (tab_pi_rho0_pi_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi_rho0_pi_diff\n";

    tab_pi_rho0_pi_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff,
        dm,
        xs_wrapper<PhotonCrossSection<
            ComputationMethod::Analytic>::xs_diff_pi_rho0_pi>);
  }

  return tab_pi_rho0_pi_diff_->get_linear(s, t, m_rho);
}

// C12
double
PhotonCrossSection<ComputationMethod::Lookup>::xs_pi0_rho_pi_rho_mediated(
    const double s, const double m_rho) {
  if (tab_pi0_rho_pi_rho_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi0_rho_pi\n");
    tab_pi0_rho_pi_rho_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        xs_wrapper<PhotonCrossSection<
            ComputationMethod::Analytic>::xs_pi0_rho_pi_rho_mediated>);
  }

  return tab_pi0_rho_pi_rho_->get_linear(s, m_rho);
}

double
PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi0_rho_pi_rho_mediated(
    const double s, const double t, const double m_rho) {
  if (tab_pi0_rho_pi_rho_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi0_rho_pi\n";
    tab_pi0_rho_pi_rho_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff,
        dm,
        xs_wrapper<PhotonCrossSection<
            ComputationMethod::Analytic>::xs_diff_pi0_rho_pi_rho_mediated>);
  }

  return tab_pi0_rho_pi_rho_diff_->get_linear(s, t, m_rho);
}

// C13
double
PhotonCrossSection<ComputationMethod::Lookup>::xs_pi_rho_pi0_rho_mediated(
    const double s, const double m_rho) {
  if (tab_pi_rho_pi0_rho_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi_rho_pi0\n");
    tab_pi_rho_pi0_rho_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        xs_wrapper<PhotonCrossSection<
            ComputationMethod::Analytic>::xs_pi_rho_pi0_rho_mediated>);
  }

  return tab_pi_rho_pi0_rho_->get_linear(s, m_rho);
}

double
PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi_rho_pi0_rho_mediated(
    const double s, const double t, const double m_rho) {
  if (tab_pi_rho_pi0_rho_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi_rho_pi0\n";
    tab_pi_rho_pi0_rho_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff,
        dm,
        xs_wrapper<PhotonCrossSection<
            ComputationMethod::Analytic>::xs_diff_pi_rho_pi0>);
  }

  return tab_pi_rho_pi0_rho_diff_->get_linear(s, t, m_rho);
}

/*----------------------------------------------------------------------------*/
/* 					Pi + Rho -> Pi + Photon channels mediated
 * by (omega) 						  */
/*----------------------------------------------------------------------------*/

// C14
double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi0_rho0_pi0(
    const double s, const double m_rho) {
  if (tab_pi0_rho0_pi0_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi0_rho0_pi0\n");
    tab_pi0_rho0_pi0_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        xs_wrapper<
            PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho0_pi0>);
  }

  return tab_pi0_rho0_pi0_->get_linear(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi0_rho0_pi0(
    const double s, const double t, const double m_rho) {
  if (tab_pi0_rho0_pi0_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi0_rho0_pi0_diff\n";
    tab_pi0_rho0_pi0_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff,
        dm,
        xs_wrapper<PhotonCrossSection<
            ComputationMethod::Analytic>::xs_diff_pi0_rho0_pi0>);
  }

  return tab_pi0_rho0_pi0_diff_->get_linear(s, t, m_rho);
}

// C15
double
PhotonCrossSection<ComputationMethod::Lookup>::xs_pi_rho_pi0_omega_mediated(
    const double s, const double m_rho) {
  if (tab_pi_rho_pi0_omega_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi_rho_pi0\n");
    tab_pi_rho_pi0_omega_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        xs_wrapper<PhotonCrossSection<
            ComputationMethod::Analytic>::xs_pi_rho_pi0_omega_mediated>);
  }

  return tab_pi_rho_pi0_omega_->get_linear(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::
    xs_diff_pi_rho_pi0_omega_mediated(const double s, const double t,
                                      const double m_rho) {
  if (tab_pi_rho_pi0_omega_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi_rho_pi0\n";
    tab_pi_rho_pi0_omega_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff,
        dm,
        xs_wrapper<PhotonCrossSection<
            ComputationMethod::Analytic>::xs_diff_pi_rho_pi0_omega_mediated>);
  }

  return tab_pi_rho_pi0_omega_diff_->get_linear(s, t, m_rho);
}

// C16
double
PhotonCrossSection<ComputationMethod::Lookup>::xs_pi0_rho_pi_omega_mediated(
    const double s, const double m_rho) {
  if (tab_pi0_rho_pi_omega_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi0_rho_pi\n");
    tab_pi0_rho_pi_omega_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        xs_wrapper<
            PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho_pi>);
  }

  return tab_pi0_rho_pi_omega_->get_linear(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::
    xs_diff_pi0_rho_pi_omega_mediated(const double s, const double t,
                                      const double m_rho) {
  if (tab_pi0_rho_pi_omega_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi0_rho_pi\n";
    tab_pi0_rho_pi_omega_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff,
        dm,
        xs_wrapper<PhotonCrossSection<
            ComputationMethod::Analytic>::xs_diff_pi0_rho_pi_omega_mediated>);
  }

  return tab_pi0_rho_pi_omega_diff_->get_linear(s, t, m_rho);
}

/*----------------------------------------------------------------------------*/
/* 					Pi + Pi -> Rho + Photon channels mediated
 * by (Pi, Rho, a1) 				*/
/*----------------------------------------------------------------------------*/
// C21
double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi_pi_rho0(
    const double s, const double m_rho) {
  if (tab_pi_pi_rho0_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi_pi_rho0\n");
    tab_pi_pi_rho0_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        xs_wrapper<
            PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_pi_rho0>);
  }

  return tab_pi_pi_rho0_->get_linear(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi_pi_rho0(
    const double s, const double t, const double m_rho) {
  if (tab_pi_pi_rho0_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi_pi_rho0_diff\n";
    tab_pi_pi_rho0_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff,
        dm,
        xs_wrapper<PhotonCrossSection<
            ComputationMethod::Analytic>::xs_diff_pi_pi_rho0>);
  }

  return tab_pi_pi_rho0_diff_->get_linear(s, t, m_rho);
}

// C22
double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi_pi0_rho(
    const double s, const double m_rho) {
  if (tab_pi_pi0_rho_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi_pi0_rho\n");
    tab_pi_pi0_rho_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        xs_wrapper<
            PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_pi0_rho>);
  }

  return tab_pi_pi0_rho_->get_linear(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi_pi0_rho(
    const double s, const double t, const double m_rho) {
  if (tab_pi_pi0_rho_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi_pi0_rho_diff\n";
    tab_pi_pi0_rho_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff,
        dm,
        xs_wrapper<PhotonCrossSection<
            ComputationMethod::Analytic>::xs_diff_pi_pi0_rho>);
  }

  return tab_pi_pi0_rho_diff_->get_linear(s, t, m_rho);
}

/*----------------------------------------------------------------------------*/
/*				 Pi + Rho -> Pi + Photon channels mediated by
 * (Pi, Rho, a1)         */
/*  			 and (Omega) summed
 */
/*----------------------------------------------------------------------------*/

// C12 + C16
double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi0_rho_pi(
    const double s, const double m_rho) {
  return xs_pi0_rho_pi_rho_mediated(s, m_rho) +
         xs_pi0_rho_pi_omega_mediated(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi0_rho_pi(
    const double s, const double t, const double m_rho) {
  return xs_diff_pi0_rho_pi_rho_mediated(s, t, m_rho) +
         xs_diff_pi0_rho_pi_omega_mediated(s, t, m_rho);
}

// C13 + C15
double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi_rho_pi0(
    const double s, const double m_rho) {
  return xs_pi_rho_pi0_rho_mediated(s, m_rho) +
         xs_pi_rho_pi0_omega_mediated(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi_rho_pi0(
    const double s, const double t, const double m_rho) {
  return xs_diff_pi_rho_pi0_rho_mediated(s, t, m_rho) +
         xs_diff_pi_rho_pi0_omega_mediated(s, t, m_rho);
}

// definition of differential xs-getters

// definition of static lookup tables. fine-tuning for parameters needed.
//
constexpr double s0 = 0.1, s1 = 18.4, t0 = -10., t1 = -0.001, ds = 0.01;
constexpr double dt = 0.01;

std::unique_ptr<TabulationND<2>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_pi_rho0_ = nullptr;

std::unique_ptr<TabulationND<2>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_pi0_rho_ = nullptr;

std::unique_ptr<TabulationND<2>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi0_rho0_pi0_ = nullptr;

std::unique_ptr<TabulationND<2>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_rho0_pi_ = nullptr;

std::unique_ptr<TabulationND<2>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi0_rho_pi_omega_ =
        nullptr;

std::unique_ptr<TabulationND<2>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi0_rho_pi_rho_ =
        nullptr;

std::unique_ptr<TabulationND<2>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_rho_pi0_omega_ =
        nullptr;

std::unique_ptr<TabulationND<2>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_rho_pi0_rho_ =
        nullptr;

std::unique_ptr<TabulationND<3>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_pi_rho0_diff_ =
        nullptr;
std::unique_ptr<TabulationND<3>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_pi0_rho_diff_ =
        nullptr;
std::unique_ptr<TabulationND<3>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi0_rho0_pi0_diff_ =
        nullptr;
std::unique_ptr<TabulationND<3>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_rho_pi0_omega_diff_ =
        nullptr;

std::unique_ptr<TabulationND<3>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_rho_pi0_rho_diff_ =
        nullptr;
std::unique_ptr<TabulationND<3>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi0_rho_pi_omega_diff_ =
        nullptr;

std::unique_ptr<TabulationND<3>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi0_rho_pi_rho_diff_ =
        nullptr;
std::unique_ptr<TabulationND<3>>

    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_rho0_pi_diff_ =
        nullptr;
