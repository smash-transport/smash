/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/photoncrosssections.h"

// using float_t = float;

float_t PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho0_pi0(
    const float_t m1, const float_t m2, const float_t m3, const float_t t1,
    const float_t t2, const float_t s, const float_t mpion,
    const float_t mrho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const float_t xs =
      to_mb * 1 / 3.0 *
      (pow(Const, 2) * pow(g_POR, 4) *
       ((pow(pow(m_omega, 2) - s, 2) *
         (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
          pow(mpion, 4) *
              (pow(mrho, 4) + 4 * pow(m_omega, 4) - 2 * pow(m_omega, 2) * s) +
          pow(m_omega, 4) *
              (pow(mrho, 4) + pow(m_omega, 4) + 2 * pow(m_omega, 2) * s +
               2 * pow(s, 2) - 2 * pow(mrho, 2) * (pow(m_omega, 2) + s)) -
          2 * pow(mpion, 2) * pow(m_omega, 2) *
              (pow(mrho, 4) + 2 * pow(m_omega, 2) * (pow(m_omega, 2) + s) -
               pow(mrho, 2) * (2 * pow(m_omega, 2) + s)))) /
            (pow(m_omega, 2) - t2) +
        (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
         3 * pow(m_omega, 8) - 4 * pow(m_omega, 6) * s -
         7 * pow(m_omega, 4) * pow(s, 2) + 4 * pow(m_omega, 2) * pow(s, 3) +
         5 * pow(s, 4) +
         pow(mrho, 4) *
             (pow(m_omega, 4) - 2 * pow(m_omega, 2) * s + 2 * pow(s, 2)) +
         pow(mrho, 2) *
             (-4 * pow(m_omega, 6) + 8 * pow(m_omega, 4) * s - 6 * pow(s, 3)) -
         2 * pow(mpion, 2) *
             (4 * pow(m_omega, 6) -
              2 * pow(mrho, 2) * pow(pow(m_omega, 2) - 2 * s, 2) +
              pow(mrho, 4) * s - 10 * pow(m_omega, 4) * s + 8 * pow(s, 3)) +
         pow(mpion, 4) *
             (pow(mrho, 4) + 2 * pow(mrho, 2) * (pow(m_omega, 2) - s) +
              4 * (pow(m_omega, 4) - 3 * pow(m_omega, 2) * s +
                   3 * pow(s, 2)))) *
            t2 -
        2 * pow(mpion, 2) * pow(m_omega, 4) * pow(t2, 2) -
        pow(mrho, 2) * pow(m_omega, 4) * pow(t2, 2) +
        pow(m_omega, 6) * pow(t2, 2) - pow(mpion, 4) * s * pow(t2, 2) +
        pow(mpion, 2) * pow(mrho, 2) * s * pow(t2, 2) +
        8 * pow(mpion, 2) * pow(m_omega, 2) * s * pow(t2, 2) +
        3 * pow(mrho, 2) * pow(m_omega, 2) * s * pow(t2, 2) -
        2 * pow(m_omega, 4) * s * pow(t2, 2) -
        8 * pow(mpion, 2) * pow(s, 2) * pow(t2, 2) -
        3 * pow(mrho, 2) * pow(s, 2) * pow(t2, 2) -
        3 * pow(m_omega, 2) * pow(s, 2) * pow(t2, 2) +
        5 * pow(s, 3) * pow(t2, 2) +
        ((pow(m_omega, 4) - 4 * pow(m_omega, 2) * s + 5 * pow(s, 2)) *
         pow(t2, 3)) /
            3. -
        (pow(pow(m_omega, 2) - s, 2) *
         (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
          pow(mpion, 4) *
              (pow(mrho, 4) + 4 * pow(m_omega, 4) - 2 * pow(m_omega, 2) * s) +
          pow(m_omega, 4) *
              (pow(mrho, 4) + pow(m_omega, 4) + 2 * pow(m_omega, 2) * s +
               2 * pow(s, 2) - 2 * pow(mrho, 2) * (pow(m_omega, 2) + s)) -
          2 * pow(mpion, 2) * pow(m_omega, 2) *
              (pow(mrho, 4) + 2 * pow(m_omega, 2) * (pow(m_omega, 2) + s) -
               pow(mrho, 2) * (2 * pow(m_omega, 2) + s)))) /
            (pow(m_omega, 2) - t1) -
        (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
         3 * pow(m_omega, 8) - 4 * pow(m_omega, 6) * s -
         7 * pow(m_omega, 4) * pow(s, 2) + 4 * pow(m_omega, 2) * pow(s, 3) +
         5 * pow(s, 4) +
         pow(mrho, 4) *
             (pow(m_omega, 4) - 2 * pow(m_omega, 2) * s + 2 * pow(s, 2)) +
         pow(mrho, 2) *
             (-4 * pow(m_omega, 6) + 8 * pow(m_omega, 4) * s - 6 * pow(s, 3)) -
         2 * pow(mpion, 2) *
             (4 * pow(m_omega, 6) -
              2 * pow(mrho, 2) * pow(pow(m_omega, 2) - 2 * s, 2) +
              pow(mrho, 4) * s - 10 * pow(m_omega, 4) * s + 8 * pow(s, 3)) +
         pow(mpion, 4) *
             (pow(mrho, 4) + 2 * pow(mrho, 2) * (pow(m_omega, 2) - s) +
              4 * (pow(m_omega, 4) - 3 * pow(m_omega, 2) * s +
                   3 * pow(s, 2)))) *
            t1 +
        2 * pow(mpion, 2) * pow(m_omega, 4) * pow(t1, 2) +
        pow(mrho, 2) * pow(m_omega, 4) * pow(t1, 2) -
        pow(m_omega, 6) * pow(t1, 2) + pow(mpion, 4) * s * pow(t1, 2) -
        pow(mpion, 2) * pow(mrho, 2) * s * pow(t1, 2) -
        8 * pow(mpion, 2) * pow(m_omega, 2) * s * pow(t1, 2) -
        3 * pow(mrho, 2) * pow(m_omega, 2) * s * pow(t1, 2) +
        2 * pow(m_omega, 4) * s * pow(t1, 2) +
        8 * pow(mpion, 2) * pow(s, 2) * pow(t1, 2) +
        3 * pow(mrho, 2) * pow(s, 2) * pow(t1, 2) +
        3 * pow(m_omega, 2) * pow(s, 2) * pow(t1, 2) -
        5 * pow(s, 3) * pow(t1, 2) -
        ((pow(m_omega, 4) - 4 * pow(m_omega, 2) * s + 5 * pow(s, 2)) *
         pow(t1, 3)) /
            3. +
        2 * (pow(m_omega, 2) - s) *
            (-pow(mpion, 8) +
             pow(mpion, 4) *
                 (4 * pow(m_omega, 4) - 7 * pow(m_omega, 2) * s + pow(s, 2) +
                  pow(mrho, 2) * (pow(m_omega, 2) + s)) +
             pow(mpion, 2) *
                 (-6 * pow(m_omega, 6) + 6 * pow(m_omega, 4) * s +
                  8 * pow(m_omega, 2) * pow(s, 2) +
                  pow(mrho, 4) * (-pow(m_omega, 2) + s) +
                  pow(mrho, 2) * (4 * pow(m_omega, 4) -
                                  7 * pow(m_omega, 2) * s - pow(s, 2))) +
             pow(m_omega, 2) *
                 (2 * pow(m_omega, 6) + pow(mrho, 4) * (pow(m_omega, 2) - s) -
                  4 * pow(m_omega, 2) * pow(s, 2) - 3 * pow(s, 3) +
                  pow(mrho, 2) * (-3 * pow(m_omega, 4) +
                                  2 * pow(m_omega, 2) * s + 3 * pow(s, 2)))) *
            log((-pow(m_omega, 2) + t2) / (-pow(m_omega, 2) + t1)))) /
      (128.0 * M_PI * pow(pow(m_omega, 2) - s, 2) *
       (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
        2 * pow(mpion, 2) * (pow(mrho, 2) + s)));

  return xs;
}

float_t PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho_pi(
    const float_t m1, const float_t m2, const float_t m3, const float_t t1,
    const float_t t2, const float_t s, const float_t mpion,
    const float_t mrho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const float_t xs =
      to_mb * 1 / 3.0 *
      (0.0024868 * pow(Const, 2) * pow(g_POR, 4) *
       ((pow(m_omega, 8) +
         pow(mpion, 4) * pow(pow(mpion, 2) - pow(mrho, 2), 2) -
         2 * pow(m_omega, 6) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
         2 * pow(m_omega, 2) * pow(mpion, 2) *
             (pow(mrho, 4) + pow(mpion, 2) * s - pow(mrho, 2) * s) +
         pow(m_omega, 4) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                            4 * pow(mpion, 2) * (pow(mrho, 2) - s) -
                            2 * pow(mrho, 2) * s + 2 * pow(s, 2))) /
            (pow(m_omega, 2) - t2) +
        3 * pow(m_omega, 4) * t2 - 8 * pow(m_omega, 2) * pow(mpion, 2) * t2 +
        4 * pow(mpion, 4) * t2 - 4 * pow(m_omega, 2) * pow(mrho, 2) * t2 +
        4 * pow(mpion, 2) * pow(mrho, 2) * t2 + pow(mrho, 4) * t2 +
        4 * pow(m_omega, 2) * s * t2 - 4 * pow(mpion, 2) * s * t2 -
        2 * pow(mrho, 2) * s * t2 + 2 * pow(s, 2) * t2 +
        pow(m_omega, 2) * pow(t2, 2) - 2 * pow(mpion, 2) * pow(t2, 2) -
        pow(mrho, 2) * pow(t2, 2) + s * pow(t2, 2) + pow(t2, 3) / 3. -
        (pow(m_omega, 8) +
         pow(mpion, 4) * pow(pow(mpion, 2) - pow(mrho, 2), 2) -
         2 * pow(m_omega, 6) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
         2 * pow(m_omega, 2) * pow(mpion, 2) *
             (pow(mrho, 4) + pow(mpion, 2) * s - pow(mrho, 2) * s) +
         pow(m_omega, 4) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                            4 * pow(mpion, 2) * (pow(mrho, 2) - s) -
                            2 * pow(mrho, 2) * s + 2 * pow(s, 2))) /
            (pow(m_omega, 2) - t1) -
        3 * pow(m_omega, 4) * t1 + 8 * pow(m_omega, 2) * pow(mpion, 2) * t1 -
        4 * pow(mpion, 4) * t1 + 4 * pow(m_omega, 2) * pow(mrho, 2) * t1 -
        4 * pow(mpion, 2) * pow(mrho, 2) * t1 - pow(mrho, 4) * t1 -
        4 * pow(m_omega, 2) * s * t1 + 4 * pow(mpion, 2) * s * t1 +
        2 * pow(mrho, 2) * s * t1 - 2 * pow(s, 2) * t1 -
        pow(m_omega, 2) * pow(t1, 2) + 2 * pow(mpion, 2) * pow(t1, 2) +
        pow(mrho, 2) * pow(t1, 2) - s * pow(t1, 2) - pow(t1, 3) / 3. +
        2 *
            (2 * pow(m_omega, 6) -
             3 * pow(m_omega, 4) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
             pow(mpion, 2) *
                 (pow(mrho, 4) + pow(mpion, 2) * s - pow(mrho, 2) * s) +
             pow(m_omega, 2) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                                4 * pow(mpion, 2) * (pow(mrho, 2) - s) -
                                2 * pow(mrho, 2) * s + 2 * pow(s, 2))) *
            log(fabs(-pow(m_omega, 2) + t2)) -
        2 *
            (2 * pow(m_omega, 6) -
             3 * pow(m_omega, 4) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
             pow(mpion, 2) *
                 (pow(mrho, 4) + pow(mpion, 2) * s - pow(mrho, 2) * s) +
             pow(m_omega, 2) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                                4 * pow(mpion, 2) * (pow(mrho, 2) - s) -
                                2 * pow(mrho, 2) * s + 2 * pow(s, 2))) *
            log(fabs(-pow(m_omega, 2) + t1)))) /
      (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
       2 * pow(mpion, 2) * (pow(mrho, 2) + s));

  return xs;
}

float_t PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_rho_pi0(
    const float_t m1, const float_t m2, const float_t m3, const float_t t1,
    const float_t t2, const float_t s, const float_t mpion,
    const float_t mrho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const float_t xs =
      to_mb * 1 / 3.0 *
      (0.0024867959858108648 * pow(Const, 2) * pow(g_POR, 4) *
       (pow(mpion, 8) * (t2 - t1) +
        pow(mpion, 6) * pow(mrho, 2) * (-2. * t2 + 2. * t1) +
        pow(mpion, 4) *
            (pow(mrho, 4) * (t2 - t1) +
             s * (4. * s * t2 - pow(t2, 2) - 4. * s * t1 + pow(t1, 2))) +
        pow(s, 2) *
            (pow(s, 2) * t2 + s * pow(t2, 2) + 0.6666666666666666 * pow(t2, 3) +
             pow(mrho, 4) * (t2 - t1) - pow(s, 2) * t1 - s * pow(t1, 2) -
             0.6666666666666666 * pow(t1, 3) +
             pow(mrho, 2) *
                 (-2. * s * t2 - pow(t2, 2) + 2. * s * t1 + pow(t1, 2))) +
        pow(mpion, 2) * s *
            (pow(mrho, 4) * (-2. * t2 + 2. * t1) +
             pow(mrho, 2) *
                 (4. * s * t2 + pow(t2, 2) - 4. * s * t1 - pow(t1, 2)) +
             s * (-4. * s * t2 - 2. * pow(t2, 2) + 4. * s * t1 +
                  2. * pow(t1, 2))))) /
      (pow(pow(m_omega, 2) - s, 2) *
       (pow(mpion, 4) + pow(mrho, 4) +
        pow(mpion, 2) * (-2. * pow(mrho, 2) - 2. * s) - 2. * pow(mrho, 2) * s +
        pow(s, 2)));

  return xs;
}

float_t PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_rho0_pi(
    const float_t m1, const float_t m2, const float_t m3, const float_t t1,
    const float_t t2, const float_t s, const float_t mpion,
    const float_t mrho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const float_t m_pi = mpion;

  const float_t xs =
      to_mb * 1 / 3.0 *
      (pow(Const, 2) * pow(ghat, 4) *
       ((pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(ma1, 8) + pow(m_pi, 8) - pow(m_pi, 4) * pow(mrho, 4) -
               2 * pow(ma1, 2) * pow(m_pi, 2) * (pow(m_pi, 2) - pow(mrho, 2)) *
                   (pow(mrho, 2) + s) +
               pow(ma1, 6) * (-4 * pow(m_pi, 2) + 2 * s) +
               pow(ma1, 4) * (4 * pow(m_pi, 4) - pow(mrho, 4) +
                              2 * pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(ma1, 8) +
               pow(m_pi, 4) * pow(pow(m_pi, 2) - pow(mrho, 2), 2) +
               2 * pow(ma1, 6) * (-2 * pow(m_pi, 2) + pow(mrho, 2) + s) +
               2 * pow(ma1, 2) * pow(m_pi, 2) *
                   (-pow(mrho, 4) + pow(m_pi, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 4) *
                   (4 * pow(m_pi, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
          pow(eta1, 2) *
              (pow(ma1, 8) + pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) -
               2 * pow(ma1, 6) * (2 * pow(m_pi, 2) + pow(mrho, 2) - s) -
               2 * pow(m_pi, 2) * pow(mrho, 4) * s +
               pow(m_pi, 4) * (3 * pow(mrho, 4) + 2 * pow(mrho, 2) * s) +
               pow(ma1, 4) * (4 * pow(m_pi, 4) + pow(mrho, 4) +
                              pow(m_pi, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2)) -
               2 * pow(ma1, 2) *
                   (pow(mrho, 2) * s * (-pow(mrho, 2) + s) +
                    pow(m_pi, 4) * (3 * pow(mrho, 2) + s) +
                    pow(m_pi, 2) *
                        (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s))))) /
            ((pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) *
             (pow(ma1, 2) - t2)) +
        (8 * pow(-2 + delta, 2) * pow(m_pi, 2) *
         (4 * pow(m_pi, 2) - pow(mrho, 2))) /
            ((pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) *
             (pow(m_pi, 2) - t2)) -
        (8 * pow(-2 + delta, 2) * pow(m_pi, 2) * t2) /
            (pow(mrho, 2) * pow(pow(m_pi, 2) - s, 2)) -
        (8 * pow(-2 + delta, 2) * pow(m_pi, 2) * t2) /
            (pow(mrho, 2) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) *
         (-8 * C4 * pow(mrho, 4) +
          pow(m_pi, 2) * (2 + delta - 8 * C4 * pow(mrho, 2)) -
          (2 + 3 * delta) * s + pow(mrho, 2) * (-2 + 3 * delta + 16 * C4 * s)) *
         t2) /
            (pow(mrho, 2) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(m_pi, 2) + pow(mrho, 2) - 2 * s) * (pow(m_pi, 2) + s) +
          eta1 * (-2 * pow(m_pi, 4) + pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                  2 * pow(s, 2) + pow(m_pi, 2) * (pow(mrho, 2) + s))) *
         t2) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (3 * pow(ma1, 4) + 4 * pow(m_pi, 4) + pow(mrho, 4) +
               pow(m_pi, 2) * (8 * pow(mrho, 2) - 4 * s) -
               4 * pow(ma1, 2) * (2 * pow(m_pi, 2) + pow(mrho, 2) - s) -
               4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          pow(eta2, 2) *
              (3 * pow(ma1, 4) + 4 * pow(m_pi, 4) + pow(mrho, 4) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) -
               4 * pow(m_pi, 2) * (pow(mrho, 2) + s) +
               4 * pow(ma1, 2) * (-2 * pow(m_pi, 2) + pow(mrho, 2) + s)) -
          2 * eta1 * eta2 *
              (3 * pow(ma1, 4) + 4 * pow(m_pi, 4) - pow(mrho, 4) +
               2 * pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) +
               pow(ma1, 2) * (-8 * pow(m_pi, 2) + 4 * s))) *
         t2) /
            (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) +
        (8 *
         (pow(delta, 2) * (8 * pow(m_pi, 4) + 3 * pow(mrho, 4) +
                           4 * pow(m_pi, 2) * (3 * pow(mrho, 2) - 2 * s) -
                           6 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          4 * pow(mrho, 4) *
              (3 + 12 * C4 * (2 * pow(m_pi, 2) - s) +
               8 * pow(C4, 2) * pow(-2 * pow(m_pi, 2) + s, 2)) -
          4 * delta * pow(mrho, 2) *
              (16 * C4 * pow(m_pi, 4) +
               2 * pow(m_pi, 2) * (3 + 6 * C4 * pow(mrho, 2) - 8 * C4 * s) +
               pow(mrho, 2) * (3 - 6 * C4 * s) + s * (-3 + 4 * C4 * s))) *
         t2) /
            (pow(mrho, 4) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (pow(m_pi, 4) * (-2 + 3 * delta - 8 * C4 * pow(mrho, 2)) +
          (pow(mrho, 2) - s) * ((-2 + 3 * delta) * s +
                                pow(mrho, 2) * (-2 + delta - 8 * C4 * s)) +
          4 * pow(m_pi, 2) *
              (2 * C4 * pow(mrho, 4) + delta * s -
               pow(mrho, 2) * (-1 + delta + 4 * C4 * s))) *
         t2) /
            (pow(mrho, 2) * (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(m_pi, 2) + s) *
              (pow(m_pi, 4) - pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) +
               (pow(mrho, 2) - s) * s) +
          eta1 * (-4 * pow(m_pi, 6) + pow(pow(mrho, 2) - s, 2) * s +
                  pow(m_pi, 4) * (3 * pow(mrho, 2) + s) -
                  pow(m_pi, 2) *
                      (pow(mrho, 4) - pow(mrho, 2) * s + 2 * pow(s, 2)))) *
         t2) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(m_pi, 8) - pow(mrho, 4) * pow(s, 2) + pow(s, 4) -
               pow(m_pi, 4) *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s - 4 * pow(s, 2)) +
               2 * pow(m_pi, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * s - 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) +
               pow(s, 2) * pow(pow(mrho, 2) + s, 2) +
               pow(m_pi, 4) * pow(pow(mrho, 2) + 2 * s, 2) -
               2 * pow(m_pi, 2) * s *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta1, 2) *
              (pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) -
               4 * pow(m_pi, 2) * pow(pow(mrho, 2) - s, 2) * s +
               pow(pow(mrho, 2) - s, 2) * pow(s, 2) +
               pow(m_pi, 4) *
                   (3 * pow(mrho, 4) - 6 * pow(mrho, 2) * s + 4 * pow(s, 2)))) *
         t2) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta1, 2) *
              (pow(ma1, 4) * s + pow(m_pi, 4) * (-3 * pow(mrho, 2) + 2 * s) +
               s * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s + pow(s, 2)) -
               2 * pow(m_pi, 2) *
                   (pow(mrho, 4) - 4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
               pow(ma1, 2) * (2 * pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) +
                              3 * s * (-pow(mrho, 2) + s))) -
          2 * eta1 * eta2 *
              (pow(ma1, 4) * s +
               s * (2 * pow(m_pi, 4) + 4 * pow(m_pi, 2) * (pow(mrho, 2) - s) +
                    s * (-2 * pow(mrho, 2) + s)) +
               pow(ma1, 2) * (pow(m_pi, 2) * (pow(mrho, 2) - 4 * s) +
                              s * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (-4 * pow(m_pi, 2) * s * (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(m_pi, 4) * (pow(mrho, 2) + 2 * s) +
               s * (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s)))) *
         t2) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (-4 * pow(m_pi, 4) *
                      (6 * C4 * pow(mrho, 4) + 2 * delta * s +
                       pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * s)) +
                  2 * pow(m_pi, 2) *
                      (4 * delta * pow(s, 2) +
                       pow(mrho, 2) * s * (6 - 7 * delta - 16 * C4 * s) +
                       2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * s)) -
                  (pow(mrho, 2) - s) * s *
                      (-2 * delta * s +
                       pow(mrho, 2) * (-6 + 3 * delta + 8 * C4 * s))) +
          eta2 * (delta * (2 * pow(m_pi, 4) * (pow(mrho, 2) + 4 * s) +
                           pow(m_pi, 2) * (2 * pow(mrho, 4) + pow(mrho, 2) * s -
                                           8 * pow(s, 2)) +
                           s * (-2 * pow(mrho, 4) - pow(mrho, 2) * s +
                                2 * pow(s, 2))) -
                  2 * pow(mrho, 2) *
                      (4 * C4 * pow(m_pi, 4) * (pow(mrho, 2) + 4 * s) +
                       pow(m_pi, 2) * (s * (5 - 16 * C4 * s) +
                                       pow(mrho, 2) * (2 - 8 * C4 * s)) +
                       s * (s * (-3 + 4 * C4 * s) +
                            pow(mrho, 2) * (-2 + 4 * C4 * s))))) *
         t2) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (8 * (eta1 - eta2) *
         (delta *
              (eta1 * (4 * pow(m_pi, 6) +
                       pow(m_pi, 4) * (7 * pow(mrho, 2) - 8 * s) +
                       pow(ma1, 4) * (pow(m_pi, 2) - s) -
                       pow(ma1, 2) * (2 * pow(m_pi, 2) + pow(mrho, 2) - 2 * s) *
                           (2 * pow(m_pi, 2) - s) +
                       pow(m_pi, 2) * s * (-8 * pow(mrho, 2) + 5 * s) +
                       s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
               eta2 *
                   (-4 * pow(m_pi, 6) - pow(m_pi, 4) * (pow(mrho, 2) - 8 * s) +
                    pow(ma1, 4) * (-pow(m_pi, 2) + s) +
                    pow(m_pi, 2) * (2 * pow(mrho, 4) - 5 * pow(s, 2)) +
                    s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s + pow(s, 2)) +
                    pow(ma1, 2) * (4 * pow(m_pi, 4) - 6 * pow(m_pi, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)))) -
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(m_pi, 6) +
                    pow(m_pi, 4) * (3 + 8 * C4 * (pow(mrho, 2) - 2 * s)) +
                    2 * C4 * pow(ma1, 4) * (pow(m_pi, 2) - s) +
                    2 * pow(m_pi, 2) * s *
                        (-1 - 6 * C4 * pow(mrho, 2) + 5 * C4 * s) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(m_pi, 4) +
                         pow(m_pi, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 6 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) +
                    s * (-(s * (1 + 2 * C4 * s)) +
                         pow(mrho, 2) * (1 + 4 * C4 * s))) +
               eta2 *
                   (2 * C4 * pow(ma1, 4) * (-pow(m_pi, 2) + s) -
                    (pow(m_pi, 2) - s) *
                        (8 * C4 * pow(m_pi, 4) - 2 * pow(mrho, 2) + s +
                         2 * C4 * pow(s, 2) +
                         pow(m_pi, 2) * (3 - 4 * C4 * (pow(mrho, 2) + 2 * s))) +
                    pow(ma1, 2) *
                        (8 * C4 * pow(m_pi, 4) +
                         2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                         pow(m_pi, 2) *
                             (1 - 2 * C4 * (pow(mrho, 2) + 6 * s)))))) *
         t2) /
            (pow(mrho, 2) * (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * pow(t2, 2)) /
            (pow(mrho, 2) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (16 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * s * pow(t2, 2)) /
            (pow(mrho, 2) * (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) * (pow(mrho, 2) - s) + 2 * eta1 * eta2 * s -
          pow(eta2, 2) * s) *
         (pow(m_pi, 4) - pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) +
          (pow(mrho, 2) - s) * s) *
         pow(t2, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (-(eta1 * (pow(m_pi, 2) + 2 * pow(mrho, 2) - 3 * s)) -
          eta2 * (pow(m_pi, 2) + s)) *
         pow(t2, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (pow(m_pi, 2) + s) *
         (-2 * eta2 * s + eta1 * (pow(m_pi, 2) - pow(mrho, 2) + s)) *
         pow(t2, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 3) *
         (eta1 * (pow(ma1, 2) - 2 * pow(m_pi, 2) - pow(mrho, 2) + s) -
          eta2 * (pow(ma1, 2) - 2 * pow(m_pi, 2) + pow(mrho, 2) + s)) *
         pow(t2, 2)) /
            (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) -
        (8 * (delta - 4 * C4 * pow(mrho, 2)) *
         (delta * (4 * pow(m_pi, 2) + 3 * pow(mrho, 2) - 2 * s) -
          2 * pow(mrho, 2) * (3 + 8 * C4 * pow(m_pi, 2) - 4 * C4 * s)) *
         pow(t2, 2)) /
            (pow(mrho, 4) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta2, 2) * s *
              (pow(ma1, 2) - 4 * pow(m_pi, 2) + pow(mrho, 2) + 3 * s) +
          pow(eta1, 2) * (2 * pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) +
                          s * (pow(ma1, 2) - 3 * pow(mrho, 2) + 3 * s)) -
          2 * eta1 * eta2 *
              (pow(m_pi, 2) * (pow(mrho, 2) - 4 * s) +
               s * (pow(ma1, 2) - 2 * pow(mrho, 2) + 3 * s))) *
         pow(t2, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (2 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (4 * delta * pow(s, 2) -
                  2 * pow(mrho, 2) * s * (-2 + 3 * delta + 8 * C4 * s) +
                  pow(mrho, 4) * (-2 + delta + 16 * C4 * s) -
                  2 * pow(m_pi, 2) *
                      (8 * C4 * pow(mrho, 4) + 4 * delta * s +
                       pow(mrho, 2) * (2 - 3 * delta - 16 * C4 * s))) +
          eta2 * (pow(m_pi, 2) * (8 * delta * s +
                                  pow(mrho, 2) * (-2 + delta - 32 * C4 * s)) +
                  s * (-4 * delta * s +
                       pow(mrho, 2) * (-2 + delta + 16 * C4 * s)))) *
         pow(t2, 2)) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) *
         (delta * (eta1 * (pow(ma1, 2) * (pow(m_pi, 2) - s) -
                           (2 * pow(m_pi, 2) + pow(mrho, 2) - 2 * s) *
                               (2 * pow(m_pi, 2) - s)) +
                   eta2 * (4 * pow(m_pi, 4) - 6 * pow(m_pi, 2) * s +
                           pow(ma1, 2) * (-pow(m_pi, 2) + s) +
                           s * (pow(mrho, 2) + 2 * s))) +
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(m_pi, 4) +
                    2 * C4 * s * (pow(ma1, 2) - pow(mrho, 2) + 2 * s) +
                    pow(m_pi, 2) *
                        (1 - 2 * C4 * (pow(ma1, 2) - pow(mrho, 2) + 6 * s))) -
               eta2 * (8 * C4 * pow(m_pi, 4) +
                       2 * C4 * s * (pow(ma1, 2) + pow(mrho, 2) + 2 * s) -
                       pow(m_pi, 2) *
                           (-1 +
                            2 * C4 * (pow(ma1, 2) + pow(mrho, 2) + 6 * s))))) *
         pow(t2, 2)) /
            (pow(mrho, 2) * (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 4) * pow(t2, 3)) /
            (3. * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                   2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (8 * pow(eta1 - eta2, 2) * (delta - 4 * C4 * pow(mrho, 2)) *
         pow(t2, 3)) /
            (3. * pow(mrho, 2) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (16 * pow(delta - 4 * C4 * pow(mrho, 2), 2) * pow(t2, 3)) /
            (3. * pow(mrho, 4) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (4 * (-2 + delta) * eta1 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         pow(t2, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 4) * (pow(ma1, 2) - s) * s * pow(t2, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) * s *
         (-2 * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-pow(mrho, 2) + s)) *
         pow(t2, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (2 * eta2 * (delta - 4 * C4 * pow(mrho, 2)) * s +
          eta1 * (-2 * delta * s + pow(mrho, 2) * (-2 + delta + 8 * C4 * s))) *
         pow(t2, 3)) /
            (3. * pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(ma1, 8) + pow(m_pi, 8) - pow(m_pi, 4) * pow(mrho, 4) -
               2 * pow(ma1, 2) * pow(m_pi, 2) * (pow(m_pi, 2) - pow(mrho, 2)) *
                   (pow(mrho, 2) + s) +
               pow(ma1, 6) * (-4 * pow(m_pi, 2) + 2 * s) +
               pow(ma1, 4) * (4 * pow(m_pi, 4) - pow(mrho, 4) +
                              2 * pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(ma1, 8) +
               pow(m_pi, 4) * pow(pow(m_pi, 2) - pow(mrho, 2), 2) +
               2 * pow(ma1, 6) * (-2 * pow(m_pi, 2) + pow(mrho, 2) + s) +
               2 * pow(ma1, 2) * pow(m_pi, 2) *
                   (-pow(mrho, 4) + pow(m_pi, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 4) *
                   (4 * pow(m_pi, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
          pow(eta1, 2) *
              (pow(ma1, 8) + pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) -
               2 * pow(ma1, 6) * (2 * pow(m_pi, 2) + pow(mrho, 2) - s) -
               2 * pow(m_pi, 2) * pow(mrho, 4) * s +
               pow(m_pi, 4) * (3 * pow(mrho, 4) + 2 * pow(mrho, 2) * s) +
               pow(ma1, 4) * (4 * pow(m_pi, 4) + pow(mrho, 4) +
                              pow(m_pi, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2)) -
               2 * pow(ma1, 2) *
                   (pow(mrho, 2) * s * (-pow(mrho, 2) + s) +
                    pow(m_pi, 4) * (3 * pow(mrho, 2) + s) +
                    pow(m_pi, 2) *
                        (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s))))) /
            ((pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) *
             (pow(ma1, 2) - t1)) -
        (8 * pow(-2 + delta, 2) * pow(m_pi, 2) *
         (4 * pow(m_pi, 2) - pow(mrho, 2))) /
            ((pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) *
             (pow(m_pi, 2) - t1)) +
        (8 * pow(-2 + delta, 2) * pow(m_pi, 2) * t1) /
            (pow(mrho, 2) * pow(pow(m_pi, 2) - s, 2)) +
        (8 * pow(-2 + delta, 2) * pow(m_pi, 2) * t1) /
            (pow(mrho, 2) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (-8 * C4 * pow(mrho, 4) +
          pow(m_pi, 2) * (2 + delta - 8 * C4 * pow(mrho, 2)) -
          (2 + 3 * delta) * s + pow(mrho, 2) * (-2 + 3 * delta + 16 * C4 * s)) *
         t1) /
            (pow(mrho, 2) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(m_pi, 2) + pow(mrho, 2) - 2 * s) * (pow(m_pi, 2) + s) +
          eta1 * (-2 * pow(m_pi, 4) + pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                  2 * pow(s, 2) + pow(m_pi, 2) * (pow(mrho, 2) + s))) *
         t1) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (3 * pow(ma1, 4) + 4 * pow(m_pi, 4) + pow(mrho, 4) +
               pow(m_pi, 2) * (8 * pow(mrho, 2) - 4 * s) -
               4 * pow(ma1, 2) * (2 * pow(m_pi, 2) + pow(mrho, 2) - s) -
               4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          pow(eta2, 2) *
              (3 * pow(ma1, 4) + 4 * pow(m_pi, 4) + pow(mrho, 4) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) -
               4 * pow(m_pi, 2) * (pow(mrho, 2) + s) +
               4 * pow(ma1, 2) * (-2 * pow(m_pi, 2) + pow(mrho, 2) + s)) -
          2 * eta1 * eta2 *
              (3 * pow(ma1, 4) + 4 * pow(m_pi, 4) - pow(mrho, 4) +
               2 * pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) +
               pow(ma1, 2) * (-8 * pow(m_pi, 2) + 4 * s))) *
         t1) /
            (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) -
        (8 *
         (pow(delta, 2) * (8 * pow(m_pi, 4) + 3 * pow(mrho, 4) +
                           4 * pow(m_pi, 2) * (3 * pow(mrho, 2) - 2 * s) -
                           6 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          4 * pow(mrho, 4) *
              (3 + 12 * C4 * (2 * pow(m_pi, 2) - s) +
               8 * pow(C4, 2) * pow(-2 * pow(m_pi, 2) + s, 2)) -
          4 * delta * pow(mrho, 2) *
              (16 * C4 * pow(m_pi, 4) +
               2 * pow(m_pi, 2) * (3 + 6 * C4 * pow(mrho, 2) - 8 * C4 * s) +
               pow(mrho, 2) * (3 - 6 * C4 * s) + s * (-3 + 4 * C4 * s))) *
         t1) /
            (pow(mrho, 4) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) *
         (pow(m_pi, 4) * (-2 + 3 * delta - 8 * C4 * pow(mrho, 2)) +
          (pow(mrho, 2) - s) * ((-2 + 3 * delta) * s +
                                pow(mrho, 2) * (-2 + delta - 8 * C4 * s)) +
          4 * pow(m_pi, 2) *
              (2 * C4 * pow(mrho, 4) + delta * s -
               pow(mrho, 2) * (-1 + delta + 4 * C4 * s))) *
         t1) /
            (pow(mrho, 2) * (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(m_pi, 2) + s) *
              (pow(m_pi, 4) - pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) +
               (pow(mrho, 2) - s) * s) +
          eta1 * (-4 * pow(m_pi, 6) + pow(pow(mrho, 2) - s, 2) * s +
                  pow(m_pi, 4) * (3 * pow(mrho, 2) + s) -
                  pow(m_pi, 2) *
                      (pow(mrho, 4) - pow(mrho, 2) * s + 2 * pow(s, 2)))) *
         t1) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(m_pi, 8) - pow(mrho, 4) * pow(s, 2) + pow(s, 4) -
               pow(m_pi, 4) *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s - 4 * pow(s, 2)) +
               2 * pow(m_pi, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * s - 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) +
               pow(s, 2) * pow(pow(mrho, 2) + s, 2) +
               pow(m_pi, 4) * pow(pow(mrho, 2) + 2 * s, 2) -
               2 * pow(m_pi, 2) * s *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta1, 2) *
              (pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) -
               4 * pow(m_pi, 2) * pow(pow(mrho, 2) - s, 2) * s +
               pow(pow(mrho, 2) - s, 2) * pow(s, 2) +
               pow(m_pi, 4) *
                   (3 * pow(mrho, 4) - 6 * pow(mrho, 2) * s + 4 * pow(s, 2)))) *
         t1) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta1, 2) *
              (pow(ma1, 4) * s + pow(m_pi, 4) * (-3 * pow(mrho, 2) + 2 * s) +
               s * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s + pow(s, 2)) -
               2 * pow(m_pi, 2) *
                   (pow(mrho, 4) - 4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
               pow(ma1, 2) * (2 * pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) +
                              3 * s * (-pow(mrho, 2) + s))) -
          2 * eta1 * eta2 *
              (pow(ma1, 4) * s +
               s * (2 * pow(m_pi, 4) + 4 * pow(m_pi, 2) * (pow(mrho, 2) - s) +
                    s * (-2 * pow(mrho, 2) + s)) +
               pow(ma1, 2) * (pow(m_pi, 2) * (pow(mrho, 2) - 4 * s) +
                              s * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (-4 * pow(m_pi, 2) * s * (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(m_pi, 4) * (pow(mrho, 2) + 2 * s) +
               s * (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s)))) *
         t1) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (4 * pow(m_pi, 4) *
                      (6 * C4 * pow(mrho, 4) + 2 * delta * s +
                       pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * s)) -
                  2 * pow(m_pi, 2) *
                      (4 * delta * pow(s, 2) +
                       pow(mrho, 2) * s * (6 - 7 * delta - 16 * C4 * s) +
                       2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * s)) +
                  (pow(mrho, 2) - s) * s *
                      (-2 * delta * s +
                       pow(mrho, 2) * (-6 + 3 * delta + 8 * C4 * s))) +
          eta2 * (-(delta * (2 * pow(m_pi, 4) * (pow(mrho, 2) + 4 * s) +
                             pow(m_pi, 2) * (2 * pow(mrho, 4) +
                                             pow(mrho, 2) * s - 8 * pow(s, 2)) +
                             s * (-2 * pow(mrho, 4) - pow(mrho, 2) * s +
                                  2 * pow(s, 2)))) +
                  2 * pow(mrho, 2) *
                      (4 * C4 * pow(m_pi, 4) * (pow(mrho, 2) + 4 * s) +
                       pow(m_pi, 2) * (s * (5 - 16 * C4 * s) +
                                       pow(mrho, 2) * (2 - 8 * C4 * s)) +
                       s * (s * (-3 + 4 * C4 * s) +
                            pow(mrho, 2) * (-2 + 4 * C4 * s))))) *
         t1) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (8 * (eta1 - eta2) *
         (delta *
              (eta1 * (4 * pow(m_pi, 6) +
                       pow(m_pi, 4) * (7 * pow(mrho, 2) - 8 * s) +
                       pow(ma1, 4) * (pow(m_pi, 2) - s) -
                       pow(ma1, 2) * (2 * pow(m_pi, 2) + pow(mrho, 2) - 2 * s) *
                           (2 * pow(m_pi, 2) - s) +
                       pow(m_pi, 2) * s * (-8 * pow(mrho, 2) + 5 * s) +
                       s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
               eta2 *
                   (-4 * pow(m_pi, 6) - pow(m_pi, 4) * (pow(mrho, 2) - 8 * s) +
                    pow(ma1, 4) * (-pow(m_pi, 2) + s) +
                    pow(m_pi, 2) * (2 * pow(mrho, 4) - 5 * pow(s, 2)) +
                    s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s + pow(s, 2)) +
                    pow(ma1, 2) * (4 * pow(m_pi, 4) - 6 * pow(m_pi, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)))) -
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(m_pi, 6) +
                    pow(m_pi, 4) * (3 + 8 * C4 * (pow(mrho, 2) - 2 * s)) +
                    2 * C4 * pow(ma1, 4) * (pow(m_pi, 2) - s) +
                    2 * pow(m_pi, 2) * s *
                        (-1 - 6 * C4 * pow(mrho, 2) + 5 * C4 * s) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(m_pi, 4) +
                         pow(m_pi, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 6 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) +
                    s * (-(s * (1 + 2 * C4 * s)) +
                         pow(mrho, 2) * (1 + 4 * C4 * s))) +
               eta2 *
                   (2 * C4 * pow(ma1, 4) * (-pow(m_pi, 2) + s) -
                    (pow(m_pi, 2) - s) *
                        (8 * C4 * pow(m_pi, 4) - 2 * pow(mrho, 2) + s +
                         2 * C4 * pow(s, 2) +
                         pow(m_pi, 2) * (3 - 4 * C4 * (pow(mrho, 2) + 2 * s))) +
                    pow(ma1, 2) *
                        (8 * C4 * pow(m_pi, 4) +
                         2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                         pow(m_pi, 2) *
                             (1 - 2 * C4 * (pow(mrho, 2) + 6 * s)))))) *
         t1) /
            (pow(mrho, 2) * (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * pow(t1, 2)) /
            (pow(mrho, 2) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (16 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * s * pow(t1, 2)) /
            (pow(mrho, 2) * (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) * (pow(mrho, 2) - s) + 2 * eta1 * eta2 * s -
          pow(eta2, 2) * s) *
         (pow(m_pi, 4) - pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) +
          (pow(mrho, 2) - s) * s) *
         pow(t1, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (-(eta1 * (pow(m_pi, 2) + 2 * pow(mrho, 2) - 3 * s)) -
          eta2 * (pow(m_pi, 2) + s)) *
         pow(t1, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (pow(m_pi, 2) + s) *
         (-2 * eta2 * s + eta1 * (pow(m_pi, 2) - pow(mrho, 2) + s)) *
         pow(t1, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 3) *
         (-(eta1 * (pow(ma1, 2) - 2 * pow(m_pi, 2) - pow(mrho, 2) + s)) +
          eta2 * (pow(ma1, 2) - 2 * pow(m_pi, 2) + pow(mrho, 2) + s)) *
         pow(t1, 2)) /
            (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) +
        (8 * (delta - 4 * C4 * pow(mrho, 2)) *
         (delta * (4 * pow(m_pi, 2) + 3 * pow(mrho, 2) - 2 * s) -
          2 * pow(mrho, 2) * (3 + 8 * C4 * pow(m_pi, 2) - 4 * C4 * s)) *
         pow(t1, 2)) /
            (pow(mrho, 4) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta2, 2) * s *
              (pow(ma1, 2) - 4 * pow(m_pi, 2) + pow(mrho, 2) + 3 * s) +
          pow(eta1, 2) * (2 * pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) +
                          s * (pow(ma1, 2) - 3 * pow(mrho, 2) + 3 * s)) -
          2 * eta1 * eta2 *
              (pow(m_pi, 2) * (pow(mrho, 2) - 4 * s) +
               s * (pow(ma1, 2) - 2 * pow(mrho, 2) + 3 * s))) *
         pow(t1, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (2 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (4 * delta * pow(s, 2) -
                  2 * pow(mrho, 2) * s * (-2 + 3 * delta + 8 * C4 * s) +
                  pow(mrho, 4) * (-2 + delta + 16 * C4 * s) -
                  2 * pow(m_pi, 2) *
                      (8 * C4 * pow(mrho, 4) + 4 * delta * s +
                       pow(mrho, 2) * (2 - 3 * delta - 16 * C4 * s))) +
          eta2 * (pow(m_pi, 2) * (8 * delta * s +
                                  pow(mrho, 2) * (-2 + delta - 32 * C4 * s)) +
                  s * (-4 * delta * s +
                       pow(mrho, 2) * (-2 + delta + 16 * C4 * s)))) *
         pow(t1, 2)) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (4 * (eta1 - eta2) *
         (delta * (eta1 * (pow(ma1, 2) * (pow(m_pi, 2) - s) -
                           (2 * pow(m_pi, 2) + pow(mrho, 2) - 2 * s) *
                               (2 * pow(m_pi, 2) - s)) +
                   eta2 * (4 * pow(m_pi, 4) - 6 * pow(m_pi, 2) * s +
                           pow(ma1, 2) * (-pow(m_pi, 2) + s) +
                           s * (pow(mrho, 2) + 2 * s))) +
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(m_pi, 4) +
                    2 * C4 * s * (pow(ma1, 2) - pow(mrho, 2) + 2 * s) +
                    pow(m_pi, 2) *
                        (1 - 2 * C4 * (pow(ma1, 2) - pow(mrho, 2) + 6 * s))) -
               eta2 * (8 * C4 * pow(m_pi, 4) +
                       2 * C4 * s * (pow(ma1, 2) + pow(mrho, 2) + 2 * s) -
                       pow(m_pi, 2) *
                           (-1 +
                            2 * C4 * (pow(ma1, 2) + pow(mrho, 2) + 6 * s))))) *
         pow(t1, 2)) /
            (pow(mrho, 2) * (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 4) * pow(t1, 3)) /
            (3. * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                   2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (8 * pow(eta1 - eta2, 2) * (delta - 4 * C4 * pow(mrho, 2)) *
         pow(t1, 3)) /
            (3. * pow(mrho, 2) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (16 * pow(delta - 4 * C4 * pow(mrho, 2), 2) * pow(t1, 3)) /
            (3. * pow(mrho, 4) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (4 * (-2 + delta) * eta1 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         pow(t1, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 4) * (pow(ma1, 2) - s) * s * pow(t1, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) * s *
         (-2 * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-pow(mrho, 2) + s)) *
         pow(t1, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (2 * eta2 * (delta - 4 * C4 * pow(mrho, 2)) * s +
          eta1 * (-2 * delta * s + pow(mrho, 2) * (-2 + delta + 8 * C4 * s))) *
         pow(t1, 3)) /
            (3. * pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (2 * pow(ma1, 6) -
               3 * pow(ma1, 4) * (2 * pow(m_pi, 2) + pow(mrho, 2) - s) +
               pow(mrho, 2) * (pow(mrho, 2) - s) * s -
               pow(m_pi, 4) * (3 * pow(mrho, 2) + s) +
               pow(m_pi, 2) * (-2 * pow(mrho, 4) + 3 * pow(mrho, 2) * s) +
               pow(ma1, 2) * (4 * pow(m_pi, 4) + pow(mrho, 4) +
                              pow(m_pi, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2))) -
          2 * eta1 * eta2 *
              (2 * pow(ma1, 6) -
               pow(m_pi, 2) * (pow(m_pi, 2) - pow(mrho, 2)) *
                   (pow(mrho, 2) + s) +
               pow(ma1, 4) * (-6 * pow(m_pi, 2) + 3 * s) +
               pow(ma1, 2) * (4 * pow(m_pi, 4) - pow(mrho, 4) +
                              2 * pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (2 * pow(ma1, 6) +
               3 * pow(ma1, 4) * (-2 * pow(m_pi, 2) + pow(mrho, 2) + s) +
               pow(m_pi, 2) *
                   (-pow(mrho, 4) + pow(m_pi, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 2) *
                   (4 * pow(m_pi, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(m_pi, 2) * (pow(mrho, 2) + s)))) *
         log(fabs(-pow(ma1, 2) + t2))) /
            (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) -
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (-2 * eta1 * eta2 *
              (pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) +
               2 * pow(ma1, 2) * pow(m_pi, 4) * s +
               pow(m_pi, 2) * (pow(ma1, 4) * (pow(mrho, 2) - 4 * s) +
                               4 * pow(ma1, 2) * (pow(mrho, 2) - s) * s +
                               pow(mrho, 2) * pow(s, 2)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (-2 * pow(mrho, 2) + s) +
                    pow(ma1, 2) * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (pow(m_pi, 8) -
               4 * pow(ma1, 2) * pow(m_pi, 2) * s *
                   (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(m_pi, 4) *
                   (pow(mrho, 2) * s + pow(ma1, 2) * (pow(mrho, 2) + 2 * s)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s))) +
          pow(eta1, 2) *
              (pow(m_pi, 8) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + 2 * pow(mrho, 4) -
                    3 * pow(ma1, 2) * (pow(mrho, 2) - s) -
                    3 * pow(mrho, 2) * s + pow(s, 2)) +
               pow(m_pi, 4) * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                               pow(ma1, 2) * (-3 * pow(mrho, 2) + 2 * s)) +
               2 * pow(m_pi, 2) *
                   (pow(ma1, 4) * (pow(mrho, 2) - 2 * s) +
                    pow(mrho, 2) * s * (-pow(mrho, 2) + s) -
                    pow(ma1, 2) * (pow(mrho, 4) - 4 * pow(mrho, 2) * s +
                                   2 * pow(s, 2))))) *
         log(fabs(-pow(ma1, 2) + t2))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (8 * (eta1 - eta2) *
         (delta *
              (eta2 *
                   (pow(m_pi, 6) * pow(mrho, 2) * (2 * pow(m_pi, 2) - s) +
                    pow(ma1, 8) * (-pow(m_pi, 2) + s) +
                    pow(ma1, 6) * (5 * pow(m_pi, 4) - 7 * pow(m_pi, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)) +
                    pow(ma1, 4) *
                        (-8 * pow(m_pi, 6) -
                         pow(m_pi, 4) * (pow(mrho, 2) - 14 * s) +
                         pow(m_pi, 2) * (2 * pow(mrho, 4) - pow(mrho, 2) * s -
                                         7 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s +
                              pow(s, 2))) +
                    pow(ma1, 2) * pow(m_pi, 2) *
                        (4 * pow(m_pi, 6) +
                         pow(m_pi, 4) * (pow(mrho, 2) - 8 * s) +
                         s * (2 * pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2)) +
                         pow(m_pi, 2) *
                             (-2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                              5 * pow(s, 2)))) +
               eta1 *
                   (pow(ma1, 8) * (pow(m_pi, 2) - s) +
                    pow(ma1, 6) *
                        (-5 * pow(m_pi, 4) + (pow(mrho, 2) - 2 * s) * s +
                         pow(m_pi, 2) * (-2 * pow(mrho, 2) + 7 * s)) +
                    pow(m_pi, 2) * pow(mrho, 2) *
                        (2 * pow(m_pi, 6) +
                         pow(m_pi, 4) * (4 * pow(mrho, 2) - 5 * s) +
                         pow(mrho, 4) * s -
                         pow(m_pi, 2) * (pow(mrho, 4) + 3 * pow(mrho, 2) * s -
                                         2 * pow(s, 2))) +
                    pow(ma1, 4) *
                        (8 * pow(m_pi, 6) +
                         pow(m_pi, 4) * (9 * pow(mrho, 2) - 14 * s) +
                         pow(m_pi, 2) * s * (-9 * pow(mrho, 2) + 7 * s) +
                         s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
                    pow(ma1, 2) *
                        (-4 * pow(m_pi, 8) +
                         pow(mrho, 4) * s * (-pow(mrho, 2) + s) +
                         pow(m_pi, 6) * (-11 * pow(mrho, 2) + 8 * s) +
                         pow(m_pi, 4) *
                             (-3 * pow(mrho, 4) + 17 * pow(mrho, 2) * s -
                              5 * pow(s, 2)) +
                         pow(m_pi, 2) *
                             (pow(mrho, 6) - 5 * pow(mrho, 2) * pow(s, 2) +
                              pow(s, 3))))) -
          2 * pow(mrho, 2) *
              (eta2 *
                   (pow(m_pi, 8) * (1 + 2 * C4 * pow(mrho, 2)) -
                    2 * C4 * pow(m_pi, 6) * pow(mrho, 2) * s +
                    2 * C4 * pow(ma1, 8) * (-pow(m_pi, 2) + s) +
                    pow(ma1, 4) *
                        (-16 * C4 * pow(m_pi, 6) +
                         pow(m_pi, 4) *
                             (-4 + 6 * C4 * pow(mrho, 2) + 28 * C4 * s) +
                         2 * pow(m_pi, 2) *
                             (pow(mrho, 2) + s - 3 * C4 * pow(mrho, 2) * s -
                              7 * C4 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 2) + s + 2 * C4 * pow(s, 2))) +
                    pow(ma1, 6) *
                        (10 * C4 * pow(m_pi, 4) +
                         2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                         pow(m_pi, 2) * (1 - 2 * C4 * (pow(mrho, 2) + 7 * s))) +
                    pow(ma1, 2) * pow(m_pi, 2) *
                        (8 * C4 * pow(m_pi, 6) -
                         2 * pow(m_pi, 4) *
                             (-2 + 3 * C4 * pow(mrho, 2) + 8 * C4 * s) +
                         s * (2 * pow(mrho, 2) + s - 2 * C4 * pow(s, 2)) +
                         2 * pow(m_pi, 2) *
                             (pow(mrho, 2) * (-1 + 3 * C4 * s) +
                              s * (-3 + 5 * C4 * s)))) +
               eta1 *
                   (pow(m_pi, 8) * (-1 + 6 * C4 * pow(mrho, 2)) +
                    2 * C4 * pow(ma1, 8) * (pow(m_pi, 2) - s) +
                    pow(m_pi, 2) * pow(mrho, 4) * s +
                    2 * pow(m_pi, 6) * pow(mrho, 2) * (2 - 5 * C4 * s) -
                    pow(ma1, 6) *
                        (10 * C4 * pow(m_pi, 4) +
                         pow(m_pi, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 7 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) -
                    pow(m_pi, 4) * pow(mrho, 2) *
                        (pow(mrho, 2) + s * (3 - 4 * C4 * s)) +
                    pow(ma1, 4) *
                        (16 * C4 * pow(m_pi, 6) +
                         2 * pow(m_pi, 4) *
                             (2 + 5 * C4 * pow(mrho, 2) - 14 * C4 * s) +
                         2 * pow(m_pi, 2) * s *
                             (-1 - 7 * C4 * pow(mrho, 2) + 7 * C4 * s) +
                         s * (-(s * (1 + 2 * C4 * s)) +
                              pow(mrho, 2) * (1 + 4 * C4 * s))) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(m_pi, 8) +
                         pow(mrho, 2) * (pow(mrho, 2) - s) * s +
                         2 * pow(m_pi, 6) *
                             (2 + 7 * C4 * pow(mrho, 2) - 8 * C4 * s) +
                         pow(m_pi, 2) * (-pow(mrho, 4) + pow(s, 2) +
                                         8 * C4 * pow(mrho, 2) * pow(s, 2) -
                                         2 * C4 * pow(s, 3)) +
                         pow(m_pi, 4) * (pow(mrho, 2) * (3 - 22 * C4 * s) +
                                         2 * s * (-3 + 5 * C4 * s)))))) *
         log(fabs(-pow(ma1, 2) + t2))) /
            ((pow(ma1, 2) - pow(m_pi, 2)) * pow(mrho, 2) * (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (16 * pow(-2 + delta, 2) * pow(m_pi, 2) *
         log(fabs(-pow(m_pi, 2) + t2))) /
            (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) -
        (8 * pow(-2 + delta, 2) *
         (3 * pow(m_pi, 4) - 4 * pow(m_pi, 2) * (pow(mrho, 2) - s) +
          pow(pow(mrho, 2) - s, 2)) *
         log(fabs(-pow(m_pi, 2) + t2))) /
            ((pow(m_pi, 2) - s) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                                   2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) * (eta1 - eta2) * pow(m_pi, 2) *
         (2 * eta1 * pow(m_pi, 2) - 2 * eta2 * pow(m_pi, 2) -
          eta1 * pow(mrho, 2)) *
         (pow(m_pi, 2) - s) * log(fabs(-pow(m_pi, 2) + t2))) /
            ((pow(ma1, 2) - pow(m_pi, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) * (eta1 - eta2) * pow(m_pi, 2) * (pow(ma1, 2) - s) *
         (pow(m_pi, 2) - s) *
         (-(eta2 * (pow(m_pi, 2) + s)) +
          eta1 * (pow(m_pi, 2) - pow(mrho, 2) + s)) *
         log(fabs(-pow(m_pi, 2) + t2))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (-(delta * (4 * pow(m_pi, 2) - pow(mrho, 2)) *
            (pow(m_pi, 2) + pow(mrho, 2) - s)) +
          2 * pow(mrho, 2) *
              (8 * C4 * pow(m_pi, 4) - pow(mrho, 2) + s +
               pow(m_pi, 2) * (3 - 8 * C4 * s))) *
         log(fabs(-pow(m_pi, 2) + t2))) /
            (pow(mrho, 2) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (2 * pow(ma1, 6) -
               3 * pow(ma1, 4) * (2 * pow(m_pi, 2) + pow(mrho, 2) - s) +
               pow(mrho, 2) * (pow(mrho, 2) - s) * s -
               pow(m_pi, 4) * (3 * pow(mrho, 2) + s) +
               pow(m_pi, 2) * (-2 * pow(mrho, 4) + 3 * pow(mrho, 2) * s) +
               pow(ma1, 2) * (4 * pow(m_pi, 4) + pow(mrho, 4) +
                              pow(m_pi, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2))) -
          2 * eta1 * eta2 *
              (2 * pow(ma1, 6) -
               pow(m_pi, 2) * (pow(m_pi, 2) - pow(mrho, 2)) *
                   (pow(mrho, 2) + s) +
               pow(ma1, 4) * (-6 * pow(m_pi, 2) + 3 * s) +
               pow(ma1, 2) * (4 * pow(m_pi, 4) - pow(mrho, 4) +
                              2 * pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (2 * pow(ma1, 6) +
               3 * pow(ma1, 4) * (-2 * pow(m_pi, 2) + pow(mrho, 2) + s) +
               pow(m_pi, 2) *
                   (-pow(mrho, 4) + pow(m_pi, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 2) *
                   (4 * pow(m_pi, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(m_pi, 2) * (pow(mrho, 2) + s)))) *
         log(fabs(-pow(ma1, 2) + t1))) /
            (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) +
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (-2 * eta1 * eta2 *
              (pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) +
               2 * pow(ma1, 2) * pow(m_pi, 4) * s +
               pow(m_pi, 2) * (pow(ma1, 4) * (pow(mrho, 2) - 4 * s) +
                               4 * pow(ma1, 2) * (pow(mrho, 2) - s) * s +
                               pow(mrho, 2) * pow(s, 2)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (-2 * pow(mrho, 2) + s) +
                    pow(ma1, 2) * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (pow(m_pi, 8) -
               4 * pow(ma1, 2) * pow(m_pi, 2) * s *
                   (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(m_pi, 4) *
                   (pow(mrho, 2) * s + pow(ma1, 2) * (pow(mrho, 2) + 2 * s)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s))) +
          pow(eta1, 2) *
              (pow(m_pi, 8) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + 2 * pow(mrho, 4) -
                    3 * pow(ma1, 2) * (pow(mrho, 2) - s) -
                    3 * pow(mrho, 2) * s + pow(s, 2)) +
               pow(m_pi, 4) * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                               pow(ma1, 2) * (-3 * pow(mrho, 2) + 2 * s)) +
               2 * pow(m_pi, 2) *
                   (pow(ma1, 4) * (pow(mrho, 2) - 2 * s) +
                    pow(mrho, 2) * s * (-pow(mrho, 2) + s) -
                    pow(ma1, 2) * (pow(mrho, 4) - 4 * pow(mrho, 2) * s +
                                   2 * pow(s, 2))))) *
         log(fabs(-pow(ma1, 2) + t1))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (8 * (eta1 - eta2) *
         (delta *
              (eta2 *
                   (pow(m_pi, 6) * pow(mrho, 2) * (2 * pow(m_pi, 2) - s) +
                    pow(ma1, 8) * (-pow(m_pi, 2) + s) +
                    pow(ma1, 6) * (5 * pow(m_pi, 4) - 7 * pow(m_pi, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)) +
                    pow(ma1, 4) *
                        (-8 * pow(m_pi, 6) -
                         pow(m_pi, 4) * (pow(mrho, 2) - 14 * s) +
                         pow(m_pi, 2) * (2 * pow(mrho, 4) - pow(mrho, 2) * s -
                                         7 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s +
                              pow(s, 2))) +
                    pow(ma1, 2) * pow(m_pi, 2) *
                        (4 * pow(m_pi, 6) +
                         pow(m_pi, 4) * (pow(mrho, 2) - 8 * s) +
                         s * (2 * pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2)) +
                         pow(m_pi, 2) *
                             (-2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                              5 * pow(s, 2)))) +
               eta1 *
                   (pow(ma1, 8) * (pow(m_pi, 2) - s) +
                    pow(ma1, 6) *
                        (-5 * pow(m_pi, 4) + (pow(mrho, 2) - 2 * s) * s +
                         pow(m_pi, 2) * (-2 * pow(mrho, 2) + 7 * s)) +
                    pow(m_pi, 2) * pow(mrho, 2) *
                        (2 * pow(m_pi, 6) +
                         pow(m_pi, 4) * (4 * pow(mrho, 2) - 5 * s) +
                         pow(mrho, 4) * s -
                         pow(m_pi, 2) * (pow(mrho, 4) + 3 * pow(mrho, 2) * s -
                                         2 * pow(s, 2))) +
                    pow(ma1, 4) *
                        (8 * pow(m_pi, 6) +
                         pow(m_pi, 4) * (9 * pow(mrho, 2) - 14 * s) +
                         pow(m_pi, 2) * s * (-9 * pow(mrho, 2) + 7 * s) +
                         s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
                    pow(ma1, 2) *
                        (-4 * pow(m_pi, 8) +
                         pow(mrho, 4) * s * (-pow(mrho, 2) + s) +
                         pow(m_pi, 6) * (-11 * pow(mrho, 2) + 8 * s) +
                         pow(m_pi, 4) *
                             (-3 * pow(mrho, 4) + 17 * pow(mrho, 2) * s -
                              5 * pow(s, 2)) +
                         pow(m_pi, 2) *
                             (pow(mrho, 6) - 5 * pow(mrho, 2) * pow(s, 2) +
                              pow(s, 3))))) -
          2 * pow(mrho, 2) *
              (eta2 *
                   (pow(m_pi, 8) * (1 + 2 * C4 * pow(mrho, 2)) -
                    2 * C4 * pow(m_pi, 6) * pow(mrho, 2) * s +
                    2 * C4 * pow(ma1, 8) * (-pow(m_pi, 2) + s) +
                    pow(ma1, 4) *
                        (-16 * C4 * pow(m_pi, 6) +
                         pow(m_pi, 4) *
                             (-4 + 6 * C4 * pow(mrho, 2) + 28 * C4 * s) +
                         2 * pow(m_pi, 2) *
                             (pow(mrho, 2) + s - 3 * C4 * pow(mrho, 2) * s -
                              7 * C4 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 2) + s + 2 * C4 * pow(s, 2))) +
                    pow(ma1, 6) *
                        (10 * C4 * pow(m_pi, 4) +
                         2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                         pow(m_pi, 2) * (1 - 2 * C4 * (pow(mrho, 2) + 7 * s))) +
                    pow(ma1, 2) * pow(m_pi, 2) *
                        (8 * C4 * pow(m_pi, 6) -
                         2 * pow(m_pi, 4) *
                             (-2 + 3 * C4 * pow(mrho, 2) + 8 * C4 * s) +
                         s * (2 * pow(mrho, 2) + s - 2 * C4 * pow(s, 2)) +
                         2 * pow(m_pi, 2) *
                             (pow(mrho, 2) * (-1 + 3 * C4 * s) +
                              s * (-3 + 5 * C4 * s)))) +
               eta1 *
                   (pow(m_pi, 8) * (-1 + 6 * C4 * pow(mrho, 2)) +
                    2 * C4 * pow(ma1, 8) * (pow(m_pi, 2) - s) +
                    pow(m_pi, 2) * pow(mrho, 4) * s +
                    2 * pow(m_pi, 6) * pow(mrho, 2) * (2 - 5 * C4 * s) -
                    pow(ma1, 6) *
                        (10 * C4 * pow(m_pi, 4) +
                         pow(m_pi, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 7 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) -
                    pow(m_pi, 4) * pow(mrho, 2) *
                        (pow(mrho, 2) + s * (3 - 4 * C4 * s)) +
                    pow(ma1, 4) *
                        (16 * C4 * pow(m_pi, 6) +
                         2 * pow(m_pi, 4) *
                             (2 + 5 * C4 * pow(mrho, 2) - 14 * C4 * s) +
                         2 * pow(m_pi, 2) * s *
                             (-1 - 7 * C4 * pow(mrho, 2) + 7 * C4 * s) +
                         s * (-(s * (1 + 2 * C4 * s)) +
                              pow(mrho, 2) * (1 + 4 * C4 * s))) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(m_pi, 8) +
                         pow(mrho, 2) * (pow(mrho, 2) - s) * s +
                         2 * pow(m_pi, 6) *
                             (2 + 7 * C4 * pow(mrho, 2) - 8 * C4 * s) +
                         pow(m_pi, 2) * (-pow(mrho, 4) + pow(s, 2) +
                                         8 * C4 * pow(mrho, 2) * pow(s, 2) -
                                         2 * C4 * pow(s, 3)) +
                         pow(m_pi, 4) * (pow(mrho, 2) * (3 - 22 * C4 * s) +
                                         2 * s * (-3 + 5 * C4 * s)))))) *
         log(fabs(-pow(ma1, 2) + t1))) /
            ((pow(ma1, 2) - pow(m_pi, 2)) * pow(mrho, 2) * (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (16 * pow(-2 + delta, 2) * pow(m_pi, 2) *
         log(fabs(-pow(m_pi, 2) + t1))) /
            (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) +
        (8 * pow(-2 + delta, 2) *
         (3 * pow(m_pi, 4) - 4 * pow(m_pi, 2) * (pow(mrho, 2) - s) +
          pow(pow(mrho, 2) - s, 2)) *
         log(fabs(-pow(m_pi, 2) + t1))) /
            ((pow(m_pi, 2) - s) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                                   2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) * (eta1 - eta2) * pow(m_pi, 2) *
         (2 * eta1 * pow(m_pi, 2) - 2 * eta2 * pow(m_pi, 2) -
          eta1 * pow(mrho, 2)) *
         (pow(m_pi, 2) - s) * log(fabs(-pow(m_pi, 2) + t1))) /
            ((pow(ma1, 2) - pow(m_pi, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) * (eta1 - eta2) * pow(m_pi, 2) * (pow(ma1, 2) - s) *
         (pow(m_pi, 2) - s) *
         (-(eta2 * (pow(m_pi, 2) + s)) +
          eta1 * (pow(m_pi, 2) - pow(mrho, 2) + s)) *
         log(fabs(-pow(m_pi, 2) + t1))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (delta * (4 * pow(m_pi, 2) - pow(mrho, 2)) *
              (pow(m_pi, 2) + pow(mrho, 2) - s) -
          2 * pow(mrho, 2) *
              (8 * C4 * pow(m_pi, 4) - pow(mrho, 2) + s +
               pow(m_pi, 2) * (3 - 8 * C4 * s))) *
         log(fabs(-pow(m_pi, 2) + t1))) /
            (pow(mrho, 2) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))))) /
      (512. * Pi);

  return xs;
}

float_t PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_pi_rho0(
    const float_t m1, const float_t m2, const float_t m3, const float_t t1,
    const float_t t2, const float_t s, const float_t mpion,
    const float_t mrho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const float_t xs =
      to_mb *
      ((pow(Const, 2) * pow(ghat, 4) *
        ((0.5 * pow(-2. + delta, 2) * pow(mpion, 2) * (0. - 1. * t2)) /
             pow(mrho, 2) +
         0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
             (-0.5 * eta2 * pow(ma1, 2) + 1. * eta1 * pow(mpion, 2) -
              1. * eta2 * pow(mpion, 2) + 0.5 * eta1 * s - 0.5 * eta2 * s) *
             (0. - 1. * t2) -
         (2. *
          (pow(mpion, 2) *
               (1.5 + 0.125 * pow(delta, 2) - 2. * C4 * pow(mrho, 2) +
                delta * (-1. + 1. * C4 * pow(mrho, 2))) +
           (-0.5 - 0.375 * pow(delta, 2) - 2. * C4 * pow(mrho, 2) +
            delta * (1. + 1. * C4 * pow(mrho, 2))) *
               s) *
          (0. - 1. * t2)) /
             pow(mrho, 2) -
         0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta1 * (1. * pow(mpion, 2) - 0.5 * s) +
              eta2 * (-0.5 * pow(ma1, 2) - 1. * pow(mpion, 2) -
                      0.5 * pow(mrho, 2) + 1. * s)) *
             (0. - 1. * t2) -
         (0.25 * (-2. + 1. * delta) *
          (-8. * C4 * pow(mrho, 4) +
           pow(mpion, 2) * (2. + 1. * delta - 8. * C4 * pow(mrho, 2)) +
           (-2. - 3. * delta) * s +
           pow(mrho, 2) * (2. + 1. * delta + 16. * C4 * s)) *
          (0. - 1. * t2)) /
             pow(mrho, 2) -
         0.09375 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-2. * pow(ma1, 4) - 4. * pow(mpion, 4) +
                   0.6666666666666666 * pow(mrho, 4) +
                   pow(ma1, 2) * (5.333333333333333 * pow(mpion, 2) -
                                  2.6666666666666665 * s) +
                   2.6666666666666665 * pow(mpion, 2) * s +
                   1.3333333333333333 * pow(mrho, 2) * s -
                   1.3333333333333333 * pow(s, 2)) +
              pow(eta1, 2) *
                  (1. * pow(ma1, 4) + 2. * pow(mpion, 4) +
                   0.3333333333333333 * pow(mrho, 4) +
                   pow(mpion, 2) *
                       (2. * pow(mrho, 2) - 1.3333333333333333 * s) -
                   1.3333333333333333 * pow(mrho, 2) * s +
                   0.6666666666666666 * pow(s, 2) +
                   pow(ma1, 2) * (-2.6666666666666665 * pow(mpion, 2) -
                                  1.3333333333333333 * pow(mrho, 2) +
                                  1.3333333333333333 * s)) +
              pow(eta2, 2) *
                  (1. * pow(ma1, 4) + 2. * pow(mpion, 4) +
                   0.3333333333333333 * pow(mrho, 4) +
                   pow(mpion, 2) *
                       (-2. * pow(mrho, 2) - 1.3333333333333333 * s) -
                   0.6666666666666666 * pow(mrho, 2) * s +
                   0.6666666666666666 * pow(s, 2) +
                   pow(ma1, 2) * (-2.6666666666666665 * pow(mpion, 2) +
                                  1.3333333333333333 * pow(mrho, 2) +
                                  1.3333333333333333 * s))) *
             (0. - 1. * t2) -
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
             (0. - 1. * t2) +
         (1. *
          (pow(eta1, 2) *
               (0.5 * pow(mrho, 4) - 1. * C4 * pow(mrho, 6) +
                pow(mpion, 2) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                pow(ma1, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                0.5 * pow(mrho, 2) * s + 0.25 * delta * pow(mrho, 2) * s -
                0.25 * delta * pow(s, 2) + 1. * C4 * pow(mrho, 2) * pow(s, 2)) +
           pow(eta2, 2) *
               (-0.5 * pow(mrho, 4) + 1. * C4 * pow(mrho, 6) +
                pow(mpion, 2) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                pow(ma1, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                0.5 * pow(mrho, 2) * s + 0.75 * delta * pow(mrho, 2) * s -
                2. * C4 * pow(mrho, 4) * s - 0.25 * delta * pow(s, 2) +
                1. * C4 * pow(mrho, 2) * pow(s, 2)) +
           eta1 * eta2 *
               (pow(ma1, 2) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                pow(mpion, 2) * (-2. * pow(mrho, 2) + 4. * C4 * pow(mrho, 4)) +
                s * (2. * C4 * pow(mrho, 4) + 0.5 * delta * s +
                     pow(mrho, 2) * (1. - 1. * delta - 2. * C4 * s)))) *
          (0. - 1. * t2)) /
             pow(mrho, 2) +
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
             (0. + 1. * pow(ma1, 2) - 1. * t2) +
         (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
          (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
             (0. + 1. * pow(mpion, 2) - 1. * t2) +
         0.03125 * pow(eta1 - 1. * eta2, 2) *
             (pow(eta2, 2) *
                  (1. * pow(Gammaa1, 2) * pow(ma1, 2) - 3. * pow(ma1, 4) -
                   2. * pow(mpion, 4) + 2. * pow(mpion, 2) * pow(mrho, 2) -
                   1. * pow(mrho, 4) +
                   pow(ma1, 2) *
                       (4. * pow(mpion, 2) - 4. * pow(mrho, 2) - 4. * s) +
                   2. * pow(mrho, 2) * s - 2. * pow(s, 2)) +
              pow(eta1, 2) *
                  (1. * pow(Gammaa1, 2) * pow(ma1, 2) - 3. * pow(ma1, 4) -
                   2. * pow(mpion, 4) - 2. * pow(mpion, 2) * pow(mrho, 2) -
                   1. * pow(mrho, 4) +
                   pow(ma1, 2) *
                       (4. * pow(mpion, 2) + 4. * pow(mrho, 2) - 4. * s) +
                   4. * pow(mrho, 2) * s - 2. * pow(s, 2)) +
              eta1 * eta2 *
                  (-2. * pow(Gammaa1, 2) * pow(ma1, 2) + 6. * pow(ma1, 4) +
                   4. * pow(mpion, 4) - 2. * pow(mrho, 4) -
                   4. * pow(mrho, 2) * s + 4. * pow(s, 2) +
                   pow(ma1, 2) * (-8. * pow(mpion, 2) + 8. * s))) *
             (0. + 1. * pow(mrho, 2) - 1. * s - 1. * t2) -
         0.03125 * pow(eta1 - 1. * eta2, 3) *
             (eta2 * (-1. * pow(ma1, 2) - 1. * pow(mrho, 2) - 1. * s) +
              eta1 * (pow(ma1, 2) - 1. * pow(mrho, 2) + s)) *
             pow(0. + 1. * pow(mrho, 2) - 1. * s - 1. * t2, 2) -
         0.010416666666666666 * pow(eta1 - 1. * eta2, 4) *
             pow(0. + 1. * pow(mrho, 2) - 1. * s - 1. * t2, 3) +
         (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
          (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
             (0. + 1. * pow(mpion, 2) + 1. * pow(mrho, 2) - 1. * s - 1. * t2) -
         (1. * (1. * eta1 - 1. * eta2) *
          (eta1 *
               (-0.5 * pow(mrho, 4) + 1. * C4 * pow(mrho, 6) +
                pow(ma1, 2) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                pow(mpion, 2) * (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4)) +
                0.25 * delta * pow(s, 2) - 1. * C4 * pow(mrho, 2) * pow(s, 2)) +
           eta2 *
               (-0.75 * pow(mrho, 4) + 0.125 * delta * pow(mrho, 4) +
                1. * C4 * pow(mrho, 6) +
                pow(mpion, 2) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                pow(ma1, 2) * (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4)) +
                0.25 * pow(mrho, 2) * s + 0.375 * delta * pow(mrho, 2) * s -
                2. * C4 * pow(mrho, 4) * s - 0.25 * delta * pow(s, 2) +
                1. * C4 * pow(mrho, 2) * pow(s, 2))) *
          (0. - 1. * pow(ma1, 2) + 2. * pow(mpion, 2) + 1. * pow(mrho, 2) -
           1. * s - 1. * t2)) /
             pow(mrho, 2) +
         0.5 * pow(1. * eta1 - 1. * eta2, 2) * (-0.5 + 1. * C4 * pow(mrho, 2)) *
             pow(0. - 1. * pow(ma1, 2) + 2. * pow(mpion, 2) +
                     1. * pow(mrho, 2) - 1. * s - 1. * t2,
                 2) +
         (0.25 *
          (32. * pow(C4, 2) * pow(mrho, 8) + 2. * pow(delta, 2) * pow(s, 2) +
           8. * C4 * pow(mrho, 6) * (-6. + delta - 8. * C4 * s) +
           2. * delta * pow(mrho, 2) * s * (-6. + delta - 8. * C4 * s) +
           pow(mrho, 4) *
               (12. - 1. * pow(delta, 2) + 8. * C4 * (6. + delta) * s +
                32. * pow(C4, 2) * pow(s, 2))) *
          (0. + t2)) /
             pow(mrho, 4) -
         (0.125 * (-2. + 1. * delta) *
          (2. + 1. * delta - 8. * C4 * pow(mrho, 2)) * (0. + 1. * pow(t2, 2))) /
             pow(mrho, 2) -
         (1. *
          (0.5 - 0.125 * pow(delta, 2) - 2. * C4 * pow(mrho, 2) +
           1. * C4 * delta * pow(mrho, 2)) *
          (0. + 1. * pow(t2, 2))) /
             pow(mrho, 2) -
         0.5 *
             (eta1 * eta2 * (1. - 2. * C4 * pow(mrho, 2)) +
              pow(eta1, 2) * (-0.5 + 1. * C4 * pow(mrho, 2)) +
              pow(eta2, 2) * (-0.5 + 1. * C4 * pow(mrho, 2))) *
             (0. + 1. * pow(t2, 2)) +
         0.0625 * pow(1. * eta1 - 1. * eta2, 4) *
             (1. * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 0.5 * s) *
             (0. + 1. * pow(t2, 2)) +
         0.03125 * pow(eta1 - 1. * eta2, 3) *
             (eta2 * (-1. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                      1. * pow(mrho, 2) - 1. * s) +
              eta1 *
                  (pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
             (0. + 1. * pow(t2, 2)) +
         0.010416666666666666 * pow(eta1 - 1. * eta2, 4) *
             (0. + 1. * pow(t2, 3)) -
         0.020833333333333332 * pow(1. * eta1 - 1. * eta2, 4) *
             (0. + 1. * pow(t2, 3)) -
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
             (1. * pow(ma1, 2) - 1. * t1) -
         (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
          (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
             (1. * pow(mpion, 2) - 1. * t1) -
         0.03125 * pow(eta1 - 1. * eta2, 2) *
             (pow(eta2, 2) *
                  (1. * pow(Gammaa1, 2) * pow(ma1, 2) - 3. * pow(ma1, 4) -
                   2. * pow(mpion, 4) + 2. * pow(mpion, 2) * pow(mrho, 2) -
                   1. * pow(mrho, 4) +
                   pow(ma1, 2) *
                       (4. * pow(mpion, 2) - 4. * pow(mrho, 2) - 4. * s) +
                   2. * pow(mrho, 2) * s - 2. * pow(s, 2)) +
              pow(eta1, 2) *
                  (1. * pow(Gammaa1, 2) * pow(ma1, 2) - 3. * pow(ma1, 4) -
                   2. * pow(mpion, 4) - 2. * pow(mpion, 2) * pow(mrho, 2) -
                   1. * pow(mrho, 4) +
                   pow(ma1, 2) *
                       (4. * pow(mpion, 2) + 4. * pow(mrho, 2) - 4. * s) +
                   4. * pow(mrho, 2) * s - 2. * pow(s, 2)) +
              eta1 * eta2 *
                  (-2. * pow(Gammaa1, 2) * pow(ma1, 2) + 6. * pow(ma1, 4) +
                   4. * pow(mpion, 4) - 2. * pow(mrho, 4) -
                   4. * pow(mrho, 2) * s + 4. * pow(s, 2) +
                   pow(ma1, 2) * (-8. * pow(mpion, 2) + 8. * s))) *
             (1. * pow(mrho, 2) - 1. * s - 1. * t1) +
         0.03125 * pow(eta1 - 1. * eta2, 3) *
             (eta2 * (-1. * pow(ma1, 2) - 1. * pow(mrho, 2) - 1. * s) +
              eta1 * (pow(ma1, 2) - 1. * pow(mrho, 2) + s)) *
             pow(1. * pow(mrho, 2) - 1. * s - 1. * t1, 2) +
         0.010416666666666666 * pow(eta1 - 1. * eta2, 4) *
             pow(1. * pow(mrho, 2) - 1. * s - 1. * t1, 3) -
         (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
          (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
             (1. * pow(mpion, 2) + 1. * pow(mrho, 2) - 1. * s - 1. * t1) +
         (0.5 * pow(-2. + delta, 2) * pow(mpion, 2) * t1) / pow(mrho, 2) +
         0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
             (-0.5 * eta2 * pow(ma1, 2) + 1. * eta1 * pow(mpion, 2) -
              1. * eta2 * pow(mpion, 2) + 0.5 * eta1 * s - 0.5 * eta2 * s) *
             t1 -
         (0.25 *
          (pow(mpion, 2) * (12. + 1. * pow(delta, 2) - 16. * C4 * pow(mrho, 2) +
                            delta * (-8. + 8. * C4 * pow(mrho, 2))) +
           (-4. - 3. * pow(delta, 2) - 16. * C4 * pow(mrho, 2) +
            delta * (8. + 8. * C4 * pow(mrho, 2))) *
               s) *
          t1) /
             pow(mrho, 2) -
         0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
             (eta1 * (1. * pow(mpion, 2) - 0.5 * s) +
              eta2 * (-0.5 * pow(ma1, 2) - 1. * pow(mpion, 2) -
                      0.5 * pow(mrho, 2) + 1. * s)) *
             t1 -
         (0.25 * (-2. + 1. * delta) *
          (-8. * C4 * pow(mrho, 4) +
           pow(mpion, 2) * (2. + 1. * delta - 8. * C4 * pow(mrho, 2)) +
           (-2. - 3. * delta) * s +
           pow(mrho, 2) * (2. + 1. * delta + 16. * C4 * s)) *
          t1) /
             pow(mrho, 2) -
         (0.25 *
          (32. * pow(C4, 2) * pow(mrho, 8) + 2. * pow(delta, 2) * pow(s, 2) +
           8. * C4 * pow(mrho, 6) * (-6. + delta - 8. * C4 * s) +
           2. * delta * pow(mrho, 2) * s * (-6. + delta - 8. * C4 * s) +
           pow(mrho, 4) *
               (12. - 1. * pow(delta, 2) + 8. * C4 * (6. + delta) * s +
                32. * pow(C4, 2) * pow(s, 2))) *
          t1) /
             pow(mrho, 4) -
         0.09375 * pow(eta1 - 1. * eta2, 2) *
             (eta1 * eta2 *
                  (-2. * pow(ma1, 4) - 4. * pow(mpion, 4) +
                   0.6666666666666666 * pow(mrho, 4) +
                   pow(ma1, 2) * (5.333333333333333 * pow(mpion, 2) -
                                  2.6666666666666665 * s) +
                   2.6666666666666665 * pow(mpion, 2) * s +
                   1.3333333333333333 * pow(mrho, 2) * s -
                   1.3333333333333333 * pow(s, 2)) +
              pow(eta1, 2) *
                  (1. * pow(ma1, 4) + 2. * pow(mpion, 4) +
                   0.3333333333333333 * pow(mrho, 4) +
                   pow(mpion, 2) *
                       (2. * pow(mrho, 2) - 1.3333333333333333 * s) -
                   1.3333333333333333 * pow(mrho, 2) * s +
                   0.6666666666666666 * pow(s, 2) +
                   pow(ma1, 2) * (-2.6666666666666665 * pow(mpion, 2) -
                                  1.3333333333333333 * pow(mrho, 2) +
                                  1.3333333333333333 * s)) +
              pow(eta2, 2) *
                  (1. * pow(ma1, 4) + 2. * pow(mpion, 4) +
                   0.3333333333333333 * pow(mrho, 4) +
                   pow(mpion, 2) *
                       (-2. * pow(mrho, 2) - 1.3333333333333333 * s) -
                   0.6666666666666666 * pow(mrho, 2) * s +
                   0.6666666666666666 * pow(s, 2) +
                   pow(ma1, 2) * (-2.6666666666666665 * pow(mpion, 2) +
                                  1.3333333333333333 * pow(mrho, 2) +
                                  1.3333333333333333 * s))) *
             t1 -
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
             t1 +
         (1. *
          (pow(eta1, 2) *
               (0.5 * pow(mrho, 4) - 1. * C4 * pow(mrho, 6) +
                pow(mpion, 2) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                pow(ma1, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                0.5 * pow(mrho, 2) * s + 0.25 * delta * pow(mrho, 2) * s -
                0.25 * delta * pow(s, 2) + 1. * C4 * pow(mrho, 2) * pow(s, 2)) +
           pow(eta2, 2) *
               (-0.5 * pow(mrho, 4) + 1. * C4 * pow(mrho, 6) +
                pow(mpion, 2) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                pow(ma1, 2) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) -
                0.5 * pow(mrho, 2) * s + 0.75 * delta * pow(mrho, 2) * s -
                2. * C4 * pow(mrho, 4) * s - 0.25 * delta * pow(s, 2) +
                1. * C4 * pow(mrho, 2) * pow(s, 2)) +
           eta1 * eta2 *
               (pow(ma1, 2) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                pow(mpion, 2) * (-2. * pow(mrho, 2) + 4. * C4 * pow(mrho, 4)) +
                s * (2. * C4 * pow(mrho, 4) + 0.5 * delta * s +
                     pow(mrho, 2) * (1. - 1. * delta - 2. * C4 * s)))) *
          t1) /
             pow(mrho, 2) +
         (0.125 * (-2. + 1. * delta) *
          (2. + 1. * delta - 8. * C4 * pow(mrho, 2)) * pow(t1, 2)) /
             pow(mrho, 2) +
         (1. *
          (0.5 - 0.125 * pow(delta, 2) - 2. * C4 * pow(mrho, 2) +
           1. * C4 * delta * pow(mrho, 2)) *
          pow(t1, 2)) /
             pow(mrho, 2) +
         0.5 *
             (eta1 * eta2 * (1. - 2. * C4 * pow(mrho, 2)) +
              pow(eta1, 2) * (-0.5 + 1. * C4 * pow(mrho, 2)) +
              pow(eta2, 2) * (-0.5 + 1. * C4 * pow(mrho, 2))) *
             pow(t1, 2) -
         0.0625 * pow(1. * eta1 - 1. * eta2, 4) *
             (1. * pow(mpion, 2) + 0.5 * pow(mrho, 2) - 0.5 * s) * pow(t1, 2) -
         0.03125 * pow(eta1 - 1. * eta2, 3) *
             (eta2 * (-1. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                      1. * pow(mrho, 2) - 1. * s) +
              eta1 *
                  (pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
             pow(t1, 2) -
         0.010416666666666666 * pow(eta1 - 1. * eta2, 4) * pow(t1, 3) +
         0.020833333333333332 * pow(1. * eta1 - 1. * eta2, 4) * pow(t1, 3) +
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
           1. * t1)) /
             pow(mrho, 2) -
         0.5 * pow(1. * eta1 - 1. * eta2, 2) * (-0.5 + 1. * C4 * pow(mrho, 2)) *
             pow(1. * pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) +
                     1. * s + 1. * t1,
                 2) +
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
          atan((0. + pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                t2) /
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
          atan((0. + pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                t2) /
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
          atan((0. + pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                t2) /
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
          atan((0. + pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                t2) /
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
          atan((0. + pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s +
                t2) /
               (Gammaa1 * ma1))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + 4. * pow(ma1, 4) +
              4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
              pow(mrho, 4) - 4. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s +
              pow(s, 2) +
              pow(ma1, 2) *
                  (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) -
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
          atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + t1) /
               (Gammaa1 * ma1))) /
             pow(mrho, 2) -
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
          atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + t1) /
               (Gammaa1 * ma1))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
              2. * pow(ma1, 2) * pow(mpion, 2) + pow(mpion, 4)) +
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
          atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + t1) /
               (Gammaa1 * ma1))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) + pow(mpion, 4) +
              2. * pow(mpion, 2) * pow(mrho, 2) + pow(mrho, 4) -
              2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s + pow(s, 2) +
              pow(ma1, 2) *
                  (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) -
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
          atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + t1) /
               (Gammaa1 * ma1))) /
             (Gammaa1 * ma1) +
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
          atan((pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s + t1) /
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
             log(fabs(-pow(ma1, 2) + t2)) -
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
               (pow(ma1, 8) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) + pow(mrho, 4) +
                               pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                               3. * pow(mrho, 2) * s) +
                pow(mpion, 2) *
                    (pow(mpion, 6) + 2. * pow(mpion, 4) * pow(mrho, 2) +
                     pow(mrho, 2) * s * (-2. * pow(mrho, 2) + 2. * s) +
                     pow(mpion, 2) *
                         (3. * pow(mrho, 4) - 5. * pow(mrho, 2) * s)) +
                pow(ma1, 2) *
                    (-4. * pow(mpion, 6) + pow(mrho, 4) * s - 1. * pow(s, 3) +
                     pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s) +
                     pow(mpion, 2) *
                         (-2. * pow(mrho, 4) + 4. * pow(mrho, 2) * s +
                          2. * pow(s, 2))))) *
          log(fabs(-pow(ma1, 2) + t2))) /
             (0.25 * pow(Gammaa1, 2) * pow(ma1, 2) + 1. * pow(ma1, 4) +
              1. * pow(mpion, 4) + 1. * pow(mpion, 2) * pow(mrho, 2) +
              0.25 * pow(mrho, 4) - 1. * pow(mpion, 2) * s -
              0.5 * pow(mrho, 2) * s + 0.25 * pow(s, 2) +
              pow(ma1, 2) *
                  (-2. * pow(mpion, 2) - 1. * pow(mrho, 2) + 1. * s)) -
         (1. *
          (eta1 * eta2 *
               (pow(ma1, 8) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                pow(ma1, 6) *
                    (-1. * pow(mrho, 4) + 2. * C4 * pow(mrho, 6) +
                     pow(mpion, 2) *
                         (-4. * pow(mrho, 2) + 8. * C4 * pow(mrho, 4)) +
                     0.5 * delta * pow(s, 2) +
                     pow(mrho, 2) * s * (2. - 1. * delta - 2. * C4 * s)) +
                pow(ma1, 2) *
                    (pow(mpion, 6) *
                         (-4. * pow(mrho, 2) + 8. * C4 * pow(mrho, 4)) +
                     pow(mrho, 4) * s *
                         ((0.5 - 0.25 * delta) * pow(mrho, 2) +
                          (-0.5 + 0.25 * delta) * s) +
                     pow(mpion, 4) *
                         (10. * C4 * pow(mrho, 6) + 0.5 * delta * pow(s, 2) +
                          pow(mrho, 4) * (-3. - 1. * delta - 8. * C4 * s) +
                          pow(mrho, 2) * s * (2. + 1. * delta - 2. * C4 * s)) +
                     pow(mpion, 2) *
                         (2. * C4 * pow(mrho, 8) - 0.5 * delta * pow(s, 3) +
                          pow(mrho, 6) * (1. - 1. * delta - 2. * C4 * s) +
                          pow(mrho, 4) * s * (-1. + 1. * delta - 2. * C4 * s) +
                          pow(mrho, 2) * pow(s, 2) * (1. + 2. * C4 * s))) +
                pow(mpion, 2) * pow(mrho, 2) *
                    (pow(mpion, 6) * (1. - 2. * C4 * pow(mrho, 2)) +
                     pow(mrho, 4) * ((-0.5 + 0.25 * delta) * pow(mrho, 2) +
                                     (0.5 - 0.25 * delta) * s) +
                     pow(mpion, 2) *
                         (-2. * C4 * pow(mrho, 6) + 0.5 * delta * pow(s, 2) +
                          pow(mrho, 2) * s *
                              (-1.5 - 0.25 * delta - 2. * C4 * s) +
                          pow(mrho, 4) * (1. + 4. * C4 * s)) +
                     pow(mpion, 4) *
                         (-4. * C4 * pow(mrho, 4) - 1. * delta * s +
                          pow(mrho, 2) * (1. + 0.5 * delta + 4. * C4 * s))) +
                pow(ma1, 4) *
                    (pow(mpion, 4) *
                         (6. * pow(mrho, 2) - 12. * C4 * pow(mrho, 4)) +
                     pow(mpion, 2) *
                         (-8. * C4 * pow(mrho, 6) - 1. * delta * pow(s, 2) +
                          pow(mrho, 4) * (3. + 0.5 * delta + 4. * C4 * s) +
                          pow(mrho, 2) * s * (-4. + 1. * delta + 4. * C4 * s)) +
                     s * (-2. * C4 * pow(mrho, 6) + 0.5 * delta * pow(s, 2) +
                          pow(mrho, 2) * s * (1. - 1.5 * delta - 2. * C4 * s) +
                          pow(mrho, 4) *
                              (-1.5 + 1.25 * delta + 4. * C4 * s)))) +
           pow(eta1, 2) *
               (pow(ma1, 8) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) +
                pow(ma1, 6) *
                    (-2. * C4 * pow(mrho, 6) +
                     pow(mpion, 2) *
                         (2. * pow(mrho, 2) - 4. * C4 * pow(mrho, 4)) -
                     0.25 * delta * pow(s, 2) +
                     pow(mrho, 4) * (1. + 1. * C4 * s) +
                     pow(mrho, 2) * s * (-1. + 0.25 * delta + 1. * C4 * s)) +
                pow(ma1, 4) *
                    (1. * C4 * pow(mrho, 8) +
                     pow(mpion, 4) *
                         (-3. * pow(mrho, 2) + 6. * C4 * pow(mrho, 4)) -
                     0.25 * delta * pow(s, 3) +
                     pow(mrho, 6) * (-0.5 - 1. * C4 * s) +
                     pow(mrho, 4) * s * (1.5 - 0.5 * delta - 1. * C4 * s) +
                     pow(mrho, 2) * pow(s, 2) *
                         (-0.5 + 0.5 * delta + 1. * C4 * s) +
                     pow(mpion, 2) *
                         (7. * C4 * pow(mrho, 6) + 0.5 * delta * pow(s, 2) +
                          pow(mrho, 4) * (-3. - 0.25 * delta - 5. * C4 * s) +
                          pow(mrho, 2) * s *
                              (2. + 0.25 * delta - 2. * C4 * s))) +
                pow(mpion, 2) * pow(mrho, 2) *
                    (pow(mpion, 6) * (-0.5 + 1. * C4 * pow(mrho, 2)) +
                     pow(mrho, 4) * ((0.5 - 0.25 * delta) * pow(mrho, 2) +
                                     (-0.5 + 0.25 * delta) * s) +
                     pow(mpion, 4) *
                         (3. * C4 * pow(mrho, 4) + 0.75 * delta * s +
                          pow(mrho, 2) * (-1. - 0.25 * delta - 3. * C4 * s)) +
                     pow(mpion, 2) *
                         (2. * C4 * pow(mrho, 6) - 0.5 * delta * pow(s, 2) +
                          pow(mrho, 4) * (-1. - 4. * C4 * s) +
                          pow(mrho, 2) * s *
                              (1.5 + 0.25 * delta + 2. * C4 * s))) +
                pow(ma1, 2) *
                    (pow(mpion, 6) *
                         (2. * pow(mrho, 2) - 4. * C4 * pow(mrho, 4)) +
                     pow(mrho, 4) * s *
                         ((-0.5 + 0.25 * delta) * pow(mrho, 2) +
                          (0.5 - 0.25 * delta) * s) +
                     pow(mpion, 2) *
                         (-3. * C4 * pow(mrho, 8) + 0.25 * delta * pow(s, 3) +
                          pow(mrho, 4) * s *
                              (-1. - 0.75 * delta - 1. * C4 * s) +
                          pow(mrho, 2) * pow(s, 2) *
                              (-0.5 + 0.5 * delta - 1. * C4 * s) +
                          pow(mrho, 6) * (0.5 + 0.5 * delta + 5. * C4 * s)) +
                     pow(mpion, 4) *
                         (-8. * C4 * pow(mrho, 6) - 0.25 * delta * pow(s, 2) +
                          pow(mrho, 2) * s *
                              (-1. - 1.25 * delta + 1. * C4 * s) +
                          pow(mrho, 4) * (3. + 0.5 * delta + 7. * C4 * s)))) +
           pow(eta2, 2) *
               (pow(ma1, 8) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) +
                pow(mpion, 6) * pow(mrho, 2) *
                    (-0.25 * delta * pow(mrho, 2) + 1. * C4 * pow(mrho, 4) +
                     pow(mpion, 2) * (-0.5 + 1. * C4 * pow(mrho, 2)) +
                     0.25 * delta * s - 1. * C4 * pow(mrho, 2) * s) +
                pow(ma1, 6) *
                    (pow(mpion, 2) *
                         (2. * pow(mrho, 2) - 4. * C4 * pow(mrho, 4)) +
                     s * (-1. * C4 * pow(mrho, 4) - 0.25 * delta * s +
                          pow(mrho, 2) * (-1. + 0.75 * delta + 1. * C4 * s))) +
                pow(ma1, 4) *
                    (-1. * C4 * pow(mrho, 8) +
                     pow(mpion, 4) *
                         (-3. * pow(mrho, 2) + 6. * C4 * pow(mrho, 4)) -
                     0.25 * delta * pow(s, 3) +
                     pow(mrho, 4) * s * (-0.75 * delta - 3. * C4 * s) +
                     pow(mrho, 2) * pow(s, 2) *
                         (-0.5 + 1. * delta + 1. * C4 * s) +
                     pow(mrho, 6) * (0.5 + 3. * C4 * s) +
                     pow(mpion, 2) *
                         (delta * (-0.25 * pow(mrho, 4) -
                                   1.25 * pow(mrho, 2) * s + 0.5 * pow(s, 2)) +
                          pow(mrho, 2) * (1. * C4 * pow(mrho, 4) + 2. * s +
                                          1. * C4 * pow(mrho, 2) * s -
                                          2. * C4 * pow(s, 2)))) +
                pow(ma1, 2) * pow(mpion, 2) *
                    (1. * C4 * pow(mrho, 8) +
                     pow(mpion, 4) *
                         (2. * pow(mrho, 2) - 4. * C4 * pow(mrho, 4)) +
                     0.25 * delta * pow(s, 3) +
                     pow(mrho, 6) * (-1.5 + 0.5 * delta - 3. * C4 * s) +
                     pow(mrho, 2) * pow(s, 2) *
                         (-0.5 - 0.5 * delta - 1. * C4 * s) +
                     pow(mrho, 4) * s * (2. - 0.25 * delta + 3. * C4 * s) +
                     pow(mpion, 2) *
                         (delta * (0.5 * pow(mrho, 4) +
                                   0.25 * pow(mrho, 2) * s - 0.25 * pow(s, 2)) +
                          pow(mrho, 2) * (-2. * C4 * pow(mrho, 4) - 1. * s +
                                          1. * C4 * pow(mrho, 2) * s +
                                          1. * C4 * pow(s, 2)))))) *
          log(fabs(-pow(ma1, 2) + t2))) /
             ((pow(ma1, 2) - 1. * pow(mpion, 2)) * pow(mrho, 2) *
              (1. * pow(ma1, 2) - 1. * pow(mpion, 2) - 1. * pow(mrho, 2) +
               1. * s)) +
         (0.5 * pow(mpion, 2) *
          (pow(eta2, 2) * pow(mpion, 2) *
               ((-2. + 1. * delta) * pow(mrho, 4) +
                (4. - 2. * delta) * pow(mrho, 2) * s +
                (-2. + 1. * delta) * pow(s, 2)) +
           eta1 * eta2 *
               (pow(mpion, 2) * ((4. - 2. * delta) * pow(mrho, 4) +
                                 (-8. + 4. * delta) * pow(mrho, 2) * s +
                                 (4. - 2. * delta) * pow(s, 2)) +
                pow(mrho, 2) * ((-1. + 0.5 * delta) * pow(mrho, 4) +
                                (2. - 1. * delta) * pow(mrho, 2) * s +
                                (-1. + 0.5 * delta) * pow(s, 2))) +
           pow(eta1, 2) *
               (pow(mrho, 2) * ((1. - 0.5 * delta) * pow(mrho, 4) +
                                (-2. + 1. * delta) * pow(mrho, 2) * s +
                                (1. - 0.5 * delta) * pow(s, 2)) +
                pow(mpion, 2) * ((-2. + 1. * delta) * pow(mrho, 4) +
                                 (4. - 2. * delta) * pow(mrho, 2) * s +
                                 (-2. + 1. * delta) * pow(s, 2)))) *
          log(fabs(-pow(mpion, 2) + t2))) /
             ((-1. * pow(ma1, 2) + 1. * pow(mpion, 2)) *
              (pow(mrho, 2) - 1. * s)) -
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
          log(fabs(-pow(mpion, 2) + t2))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) + pow(mpion, 4) +
              2. * pow(mpion, 2) * pow(mrho, 2) + pow(mrho, 4) -
              2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s + pow(s, 2) +
              pow(ma1, 2) *
                  (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) -
         (0.5 * pow(mpion, 2) *
          (eta1 * eta2 *
               ((1. - 0.5 * delta) * pow(mrho, 6) +
                (-4. + 2. * delta) * pow(mrho, 4) * s +
                (5. - 2.5 * delta) * pow(mrho, 2) * pow(s, 2) +
                (-2. + 1. * delta) * pow(s, 3) +
                pow(mpion, 2) * ((4. - 2. * delta) * pow(mrho, 4) +
                                 (-8. + 4. * delta) * pow(mrho, 2) * s +
                                 (4. - 2. * delta) * pow(s, 2))) +
           pow(eta2, 2) *
               ((-1. + 0.5 * delta) * pow(mrho, 6) +
                (3. - 1.5 * delta) * pow(mrho, 4) * s +
                (-3. + 1.5 * delta) * pow(mrho, 2) * pow(s, 2) +
                (1. - 0.5 * delta) * pow(s, 3) +
                pow(mpion, 2) * ((-2. + 1. * delta) * pow(mrho, 4) +
                                 (4. - 2. * delta) * pow(mrho, 2) * s +
                                 (-2. + 1. * delta) * pow(s, 2))) +
           pow(eta1, 2) *
               (s * ((1. - 0.5 * delta) * pow(mrho, 4) +
                     (-2. + 1. * delta) * pow(mrho, 2) * s +
                     (1. - 0.5 * delta) * pow(s, 2)) +
                pow(mpion, 2) * ((-2. + 1. * delta) * pow(mrho, 4) +
                                 (4. - 2. * delta) * pow(mrho, 2) * s +
                                 (-2. + 1. * delta) * pow(s, 2)))) *
          log(fabs(-pow(mpion, 2) - pow(mrho, 2) + s + t2))) /
             ((pow(mrho, 2) - 1. * s) *
              (-1. * pow(ma1, 2) + pow(mpion, 2) + pow(mrho, 2) - 1. * s)) +
         0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
             log(fabs(-pow(mpion, 2) - pow(mrho, 2) + s + t2)) +
         (0.25000000000000006 * pow(2. - 1. * delta, 2) *
          (7.999999999999998 * pow(mpion, 4) -
           5.999999999999998 * pow(mpion, 2) * s + 1. * pow(s, 2)) *
          log(fabs(-pow(mpion, 2) - pow(mrho, 2) + s + t2))) /
             (pow(mrho, 2) - 1. * s) -
         (1. * (-2. + 1. * delta) *
          ((0.5 - 0.25 * delta) * pow(mrho, 2) * s +
           pow(mpion, 2) * (4. * C4 * pow(mrho, 4) + 1. * delta * s +
                            pow(mrho, 2) * (-2. - 4. * C4 * s))) *
          log(fabs(-pow(mpion, 2) - pow(mrho, 2) + s + t2))) /
             pow(mrho, 2) -
         0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
             log(fabs(-pow(mpion, 2) - pow(mrho, 2) + s + t2)) -
         (2. *
          (1. * pow(2. - 1. * delta, 2) * pow(mpion, 4) * pow(mrho, 2) +
           0.12500000000000003 * pow(2. - 1. * delta, 2) * pow(mrho, 4) * s +
           pow(mpion, 2) * (C4 * (4. - 2. * delta) * pow(mrho, 6) +
                            (-1. + 0.5 * delta) * delta * pow(s, 2) +
                            pow(mrho, 2) * s *
                                (-1. + 3. * delta - 1.25 * pow(delta, 2) +
                                 4. * C4 * s - 2. * C4 * delta * s) +
                            pow(mrho, 4) * (-2. + 1. * delta - 8. * C4 * s +
                                            4. * C4 * delta * s))) *
          log(fabs(-pow(mpion, 2) - pow(mrho, 2) + s + t2))) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
         (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
          (eta2 * pow(mpion, 2) *
               (pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                pow(ma1, 2) * (-1. * pow(mrho, 2) + 1. * s)) +
           eta1 * (pow(ma1, 2) * (-0.5 * pow(mrho, 4) +
                                  pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                                  0.5 * pow(mrho, 2) * s) +
                   pow(mpion, 2) *
                       (0.5 * pow(mrho, 4) - 0.5 * pow(mrho, 2) * s +
                        pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s)))) *
          log(fabs(-pow(mpion, 2) - pow(mrho, 2) + s + t2))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
              2. * pow(ma1, 2) * pow(mpion, 2) + pow(mpion, 4)) +
         (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * (0.5 * pow(Gammaa1, 4) * pow(ma1, 4) - 0.5 * pow(ma1, 8) +
                   0.5 * pow(mpion, 8) +
                   0.5 * pow(ma1, 4) * pow(mpion, 2) * pow(mrho, 2) -
                   0.5 * pow(mpion, 6) * pow(mrho, 2) +
                   pow(Gammaa1, 2) * pow(ma1, 2) *
                       (pow(mpion, 2) *
                            (1. * pow(mpion, 2) + 1.5 * pow(mrho, 2) - 2. * s) +
                        pow(ma1, 2) * (-1. * pow(mpion, 2) +
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
                   2. * pow(mrho, 2) * s + pow(s, 2) - 4. * pow(mpion, 2) * t2 -
                   2. * pow(mrho, 2) * t2 + 2. * s * t2 + pow(t2, 2) +
                   pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                  2. * s + 2. * t2)))) /
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
                pow(ma1, 6) *
                    (1. * pow(mpion, 2) + 1. * pow(mrho, 2) - 0.5 * s) +
                pow(ma1, 4) * (-0.5 * pow(mrho, 4) +
                               (-0.5 * pow(mpion, 2) + 0.5 * s) * s) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (1. * pow(mpion, 4) + 0.5 * pow(mrho, 4) +
                     pow(mpion, 2) * (2. * pow(mrho, 2) - 1.5 * s) -
                     1. * pow(mrho, 2) * s + 0.5 * pow(s, 2) +
                     pow(ma1, 2) *
                         (-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + 1.5 * s)) +
                pow(mpion, 2) *
                    (0.5 * pow(mpion, 6) + 0.5 * pow(mpion, 4) * s +
                     pow(mpion, 2) * (-0.5 * pow(mrho, 4) +
                                      2. * pow(mrho, 2) * s - 1.5 * pow(s, 2)) +
                     s * (0.5 * pow(mrho, 4) - 1. * pow(mrho, 2) * s +
                          0.5 * pow(s, 2))) +
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
                   2. * pow(mrho, 2) * s + pow(s, 2) - 4. * pow(mpion, 2) * t2 -
                   2. * pow(mrho, 2) * t2 + 2. * s * t2 + pow(t2, 2) +
                   pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                  2. * s + 2. * t2)))) /
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
                   2. * pow(mrho, 2) * s + pow(s, 2) - 4. * pow(mpion, 2) * t2 -
                   2. * pow(mrho, 2) * t2 + 2. * s * t2 + pow(t2, 2) +
                   pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                  2. * s + 2. * t2)))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + 4. * pow(ma1, 4) +
              4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
              pow(mrho, 4) - 4. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s +
              pow(s, 2) +
              pow(ma1, 2) *
                  (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) +
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
             log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) +
                      4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                      1. * pow(mrho, 4) - 4. * pow(mpion, 2) * s -
                      2. * pow(mrho, 2) * s + pow(s, 2) -
                      4. * pow(mpion, 2) * t2 - 2. * pow(mrho, 2) * t2 +
                      2. * s * t2 + pow(t2, 2) +
                      pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                     2. * s + 2. * t2))) -
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
          log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) +
                   4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                   1. * pow(mrho, 4) - 4. * pow(mpion, 2) * s -
                   2. * pow(mrho, 2) * s + pow(s, 2) - 4. * pow(mpion, 2) * t2 -
                   2. * pow(mrho, 2) * t2 + 2. * s * t2 + pow(t2, 2) +
                   pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                  2. * s + 2. * t2)))) /
             pow(mrho, 2) -
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
             log(fabs(-pow(ma1, 2) + t1)) +
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
               (pow(ma1, 8) +
                pow(ma1, 6) *
                    (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s) +
                pow(ma1, 4) * (6. * pow(mpion, 4) + pow(mrho, 4) +
                               pow(mpion, 2) * (6. * pow(mrho, 2) - 4. * s) -
                               3. * pow(mrho, 2) * s) +
                pow(mpion, 2) *
                    (pow(mpion, 6) + 2. * pow(mpion, 4) * pow(mrho, 2) +
                     pow(mrho, 2) * s * (-2. * pow(mrho, 2) + 2. * s) +
                     pow(mpion, 2) *
                         (3. * pow(mrho, 4) - 5. * pow(mrho, 2) * s)) +
                pow(ma1, 2) *
                    (-4. * pow(mpion, 6) + pow(mrho, 4) * s - 1. * pow(s, 3) +
                     pow(mpion, 4) * (-6. * pow(mrho, 2) + 2. * s) +
                     pow(mpion, 2) *
                         (-2. * pow(mrho, 4) + 4. * pow(mrho, 2) * s +
                          2. * pow(s, 2))))) *
          log(fabs(-pow(ma1, 2) + t1))) /
             (0.25 * pow(Gammaa1, 2) * pow(ma1, 2) + 1. * pow(ma1, 4) +
              1. * pow(mpion, 4) + 1. * pow(mpion, 2) * pow(mrho, 2) +
              0.25 * pow(mrho, 4) - 1. * pow(mpion, 2) * s -
              0.5 * pow(mrho, 2) * s + 0.25 * pow(s, 2) +
              pow(ma1, 2) *
                  (-2. * pow(mpion, 2) - 1. * pow(mrho, 2) + 1. * s)) +
         (1. *
          (eta1 * eta2 *
               (pow(ma1, 8) * (1. * pow(mrho, 2) - 2. * C4 * pow(mrho, 4)) +
                pow(ma1, 6) *
                    (-1. * pow(mrho, 4) + 2. * C4 * pow(mrho, 6) +
                     pow(mpion, 2) *
                         (-4. * pow(mrho, 2) + 8. * C4 * pow(mrho, 4)) +
                     0.5 * delta * pow(s, 2) +
                     pow(mrho, 2) * s * (2. - 1. * delta - 2. * C4 * s)) +
                pow(ma1, 2) *
                    (pow(mpion, 6) *
                         (-4. * pow(mrho, 2) + 8. * C4 * pow(mrho, 4)) +
                     pow(mrho, 4) * s *
                         ((0.5 - 0.25 * delta) * pow(mrho, 2) +
                          (-0.5 + 0.25 * delta) * s) +
                     pow(mpion, 4) *
                         (10. * C4 * pow(mrho, 6) + 0.5 * delta * pow(s, 2) +
                          pow(mrho, 4) * (-3. - 1. * delta - 8. * C4 * s) +
                          pow(mrho, 2) * s * (2. + 1. * delta - 2. * C4 * s)) +
                     pow(mpion, 2) *
                         (2. * C4 * pow(mrho, 8) - 0.5 * delta * pow(s, 3) +
                          pow(mrho, 6) * (1. - 1. * delta - 2. * C4 * s) +
                          pow(mrho, 4) * s * (-1. + 1. * delta - 2. * C4 * s) +
                          pow(mrho, 2) * pow(s, 2) * (1. + 2. * C4 * s))) +
                pow(mpion, 2) * pow(mrho, 2) *
                    (pow(mpion, 6) * (1. - 2. * C4 * pow(mrho, 2)) +
                     pow(mrho, 4) * ((-0.5 + 0.25 * delta) * pow(mrho, 2) +
                                     (0.5 - 0.25 * delta) * s) +
                     pow(mpion, 2) *
                         (-2. * C4 * pow(mrho, 6) + 0.5 * delta * pow(s, 2) +
                          pow(mrho, 2) * s *
                              (-1.5 - 0.25 * delta - 2. * C4 * s) +
                          pow(mrho, 4) * (1. + 4. * C4 * s)) +
                     pow(mpion, 4) *
                         (-4. * C4 * pow(mrho, 4) - 1. * delta * s +
                          pow(mrho, 2) * (1. + 0.5 * delta + 4. * C4 * s))) +
                pow(ma1, 4) *
                    (pow(mpion, 4) *
                         (6. * pow(mrho, 2) - 12. * C4 * pow(mrho, 4)) +
                     pow(mpion, 2) *
                         (-8. * C4 * pow(mrho, 6) - 1. * delta * pow(s, 2) +
                          pow(mrho, 4) * (3. + 0.5 * delta + 4. * C4 * s) +
                          pow(mrho, 2) * s * (-4. + 1. * delta + 4. * C4 * s)) +
                     s * (-2. * C4 * pow(mrho, 6) + 0.5 * delta * pow(s, 2) +
                          pow(mrho, 2) * s * (1. - 1.5 * delta - 2. * C4 * s) +
                          pow(mrho, 4) *
                              (-1.5 + 1.25 * delta + 4. * C4 * s)))) +
           pow(eta1, 2) *
               (pow(ma1, 8) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) +
                pow(ma1, 6) *
                    (-2. * C4 * pow(mrho, 6) +
                     pow(mpion, 2) *
                         (2. * pow(mrho, 2) - 4. * C4 * pow(mrho, 4)) -
                     0.25 * delta * pow(s, 2) +
                     pow(mrho, 4) * (1. + 1. * C4 * s) +
                     pow(mrho, 2) * s * (-1. + 0.25 * delta + 1. * C4 * s)) +
                pow(ma1, 4) *
                    (1. * C4 * pow(mrho, 8) +
                     pow(mpion, 4) *
                         (-3. * pow(mrho, 2) + 6. * C4 * pow(mrho, 4)) -
                     0.25 * delta * pow(s, 3) +
                     pow(mrho, 6) * (-0.5 - 1. * C4 * s) +
                     pow(mrho, 4) * s * (1.5 - 0.5 * delta - 1. * C4 * s) +
                     pow(mrho, 2) * pow(s, 2) *
                         (-0.5 + 0.5 * delta + 1. * C4 * s) +
                     pow(mpion, 2) *
                         (7. * C4 * pow(mrho, 6) + 0.5 * delta * pow(s, 2) +
                          pow(mrho, 4) * (-3. - 0.25 * delta - 5. * C4 * s) +
                          pow(mrho, 2) * s *
                              (2. + 0.25 * delta - 2. * C4 * s))) +
                pow(mpion, 2) * pow(mrho, 2) *
                    (pow(mpion, 6) * (-0.5 + 1. * C4 * pow(mrho, 2)) +
                     pow(mrho, 4) * ((0.5 - 0.25 * delta) * pow(mrho, 2) +
                                     (-0.5 + 0.25 * delta) * s) +
                     pow(mpion, 4) *
                         (3. * C4 * pow(mrho, 4) + 0.75 * delta * s +
                          pow(mrho, 2) * (-1. - 0.25 * delta - 3. * C4 * s)) +
                     pow(mpion, 2) *
                         (2. * C4 * pow(mrho, 6) - 0.5 * delta * pow(s, 2) +
                          pow(mrho, 4) * (-1. - 4. * C4 * s) +
                          pow(mrho, 2) * s *
                              (1.5 + 0.25 * delta + 2. * C4 * s))) +
                pow(ma1, 2) *
                    (pow(mpion, 6) *
                         (2. * pow(mrho, 2) - 4. * C4 * pow(mrho, 4)) +
                     pow(mrho, 4) * s *
                         ((-0.5 + 0.25 * delta) * pow(mrho, 2) +
                          (0.5 - 0.25 * delta) * s) +
                     pow(mpion, 2) *
                         (-3. * C4 * pow(mrho, 8) + 0.25 * delta * pow(s, 3) +
                          pow(mrho, 4) * s *
                              (-1. - 0.75 * delta - 1. * C4 * s) +
                          pow(mrho, 2) * pow(s, 2) *
                              (-0.5 + 0.5 * delta - 1. * C4 * s) +
                          pow(mrho, 6) * (0.5 + 0.5 * delta + 5. * C4 * s)) +
                     pow(mpion, 4) *
                         (-8. * C4 * pow(mrho, 6) - 0.25 * delta * pow(s, 2) +
                          pow(mrho, 2) * s *
                              (-1. - 1.25 * delta + 1. * C4 * s) +
                          pow(mrho, 4) * (3. + 0.5 * delta + 7. * C4 * s)))) +
           pow(eta2, 2) *
               (pow(ma1, 8) * (-0.5 * pow(mrho, 2) + 1. * C4 * pow(mrho, 4)) +
                pow(mpion, 6) * pow(mrho, 2) *
                    (-0.25 * delta * pow(mrho, 2) + 1. * C4 * pow(mrho, 4) +
                     pow(mpion, 2) * (-0.5 + 1. * C4 * pow(mrho, 2)) +
                     0.25 * delta * s - 1. * C4 * pow(mrho, 2) * s) +
                pow(ma1, 6) *
                    (pow(mpion, 2) *
                         (2. * pow(mrho, 2) - 4. * C4 * pow(mrho, 4)) +
                     s * (-1. * C4 * pow(mrho, 4) - 0.25 * delta * s +
                          pow(mrho, 2) * (-1. + 0.75 * delta + 1. * C4 * s))) +
                pow(ma1, 4) *
                    (-1. * C4 * pow(mrho, 8) +
                     pow(mpion, 4) *
                         (-3. * pow(mrho, 2) + 6. * C4 * pow(mrho, 4)) -
                     0.25 * delta * pow(s, 3) +
                     pow(mrho, 4) * s * (-0.75 * delta - 3. * C4 * s) +
                     pow(mrho, 2) * pow(s, 2) *
                         (-0.5 + 1. * delta + 1. * C4 * s) +
                     pow(mrho, 6) * (0.5 + 3. * C4 * s) +
                     pow(mpion, 2) *
                         (delta * (-0.25 * pow(mrho, 4) -
                                   1.25 * pow(mrho, 2) * s + 0.5 * pow(s, 2)) +
                          pow(mrho, 2) * (1. * C4 * pow(mrho, 4) + 2. * s +
                                          1. * C4 * pow(mrho, 2) * s -
                                          2. * C4 * pow(s, 2)))) +
                pow(ma1, 2) * pow(mpion, 2) *
                    (1. * C4 * pow(mrho, 8) +
                     pow(mpion, 4) *
                         (2. * pow(mrho, 2) - 4. * C4 * pow(mrho, 4)) +
                     0.25 * delta * pow(s, 3) +
                     pow(mrho, 6) * (-1.5 + 0.5 * delta - 3. * C4 * s) +
                     pow(mrho, 2) * pow(s, 2) *
                         (-0.5 - 0.5 * delta - 1. * C4 * s) +
                     pow(mrho, 4) * s * (2. - 0.25 * delta + 3. * C4 * s) +
                     pow(mpion, 2) *
                         (delta * (0.5 * pow(mrho, 4) +
                                   0.25 * pow(mrho, 2) * s - 0.25 * pow(s, 2)) +
                          pow(mrho, 2) * (-2. * C4 * pow(mrho, 4) - 1. * s +
                                          1. * C4 * pow(mrho, 2) * s +
                                          1. * C4 * pow(s, 2)))))) *
          log(fabs(-pow(ma1, 2) + t1))) /
             ((pow(ma1, 2) - 1. * pow(mpion, 2)) * pow(mrho, 2) *
              (1. * pow(ma1, 2) - 1. * pow(mpion, 2) - 1. * pow(mrho, 2) +
               1. * s)) -
         (0.5 * pow(mpion, 2) *
          (pow(eta2, 2) * pow(mpion, 2) *
               ((-2. + 1. * delta) * pow(mrho, 4) +
                (4. - 2. * delta) * pow(mrho, 2) * s +
                (-2. + 1. * delta) * pow(s, 2)) +
           eta1 * eta2 *
               (pow(mpion, 2) * ((4. - 2. * delta) * pow(mrho, 4) +
                                 (-8. + 4. * delta) * pow(mrho, 2) * s +
                                 (4. - 2. * delta) * pow(s, 2)) +
                pow(mrho, 2) * ((-1. + 0.5 * delta) * pow(mrho, 4) +
                                (2. - 1. * delta) * pow(mrho, 2) * s +
                                (-1. + 0.5 * delta) * pow(s, 2))) +
           pow(eta1, 2) *
               (pow(mrho, 2) * ((1. - 0.5 * delta) * pow(mrho, 4) +
                                (-2. + 1. * delta) * pow(mrho, 2) * s +
                                (1. - 0.5 * delta) * pow(s, 2)) +
                pow(mpion, 2) * ((-2. + 1. * delta) * pow(mrho, 4) +
                                 (4. - 2. * delta) * pow(mrho, 2) * s +
                                 (-2. + 1. * delta) * pow(s, 2)))) *
          log(fabs(-pow(mpion, 2) + t1))) /
             ((-1. * pow(ma1, 2) + 1. * pow(mpion, 2)) *
              (pow(mrho, 2) - 1. * s)) +
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
          log(fabs(-pow(mpion, 2) + t1))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) + pow(mpion, 4) +
              2. * pow(mpion, 2) * pow(mrho, 2) + pow(mrho, 4) -
              2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s + pow(s, 2) +
              pow(ma1, 2) *
                  (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) +
         (0.5 * pow(mpion, 2) *
          (eta1 * eta2 *
               ((1. - 0.5 * delta) * pow(mrho, 6) +
                (-4. + 2. * delta) * pow(mrho, 4) * s +
                (5. - 2.5 * delta) * pow(mrho, 2) * pow(s, 2) +
                (-2. + 1. * delta) * pow(s, 3) +
                pow(mpion, 2) * ((4. - 2. * delta) * pow(mrho, 4) +
                                 (-8. + 4. * delta) * pow(mrho, 2) * s +
                                 (4. - 2. * delta) * pow(s, 2))) +
           pow(eta2, 2) *
               ((-1. + 0.5 * delta) * pow(mrho, 6) +
                (3. - 1.5 * delta) * pow(mrho, 4) * s +
                (-3. + 1.5 * delta) * pow(mrho, 2) * pow(s, 2) +
                (1. - 0.5 * delta) * pow(s, 3) +
                pow(mpion, 2) * ((-2. + 1. * delta) * pow(mrho, 4) +
                                 (4. - 2. * delta) * pow(mrho, 2) * s +
                                 (-2. + 1. * delta) * pow(s, 2))) +
           pow(eta1, 2) *
               (s * ((1. - 0.5 * delta) * pow(mrho, 4) +
                     (-2. + 1. * delta) * pow(mrho, 2) * s +
                     (1. - 0.5 * delta) * pow(s, 2)) +
                pow(mpion, 2) * ((-2. + 1. * delta) * pow(mrho, 4) +
                                 (4. - 2. * delta) * pow(mrho, 2) * s +
                                 (-2. + 1. * delta) * pow(s, 2)))) *
          log(fabs(-pow(mpion, 2) - pow(mrho, 2) + s + t1))) /
             ((pow(mrho, 2) - 1. * s) *
              (-1. * pow(ma1, 2) + pow(mpion, 2) + pow(mrho, 2) - 1. * s)) -
         0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
             log(fabs(-pow(mpion, 2) + t1)) -
         (0.25000000000000006 * pow(2. - 1. * delta, 2) *
          (7.999999999999998 * pow(mpion, 4) -
           5.999999999999998 * pow(mpion, 2) * s + 1. * pow(s, 2)) *
          log(fabs(-pow(mpion, 2) + t1))) /
             (pow(mrho, 2) - 1. * s) +
         (1. * (-2. + 1. * delta) *
          ((0.5 - 0.25 * delta) * pow(mrho, 2) * s +
           pow(mpion, 2) * (4. * C4 * pow(mrho, 4) + 1. * delta * s +
                            pow(mrho, 2) * (-2. - 4. * C4 * s))) *
          log(fabs(-pow(mpion, 2) + t1))) /
             pow(mrho, 2) +
         0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
             log(fabs(-pow(mpion, 2) - pow(mrho, 2) + s + t1)) +
         (2. *
          (1. * pow(2. - 1. * delta, 2) * pow(mpion, 4) * pow(mrho, 2) +
           0.12500000000000003 * pow(2. - 1. * delta, 2) * pow(mrho, 4) * s +
           pow(mpion, 2) * (C4 * (4. - 2. * delta) * pow(mrho, 6) +
                            (-1. + 0.5 * delta) * delta * pow(s, 2) +
                            pow(mrho, 2) * s *
                                (-1. + 3. * delta - 1.25 * pow(delta, 2) +
                                 4. * C4 * s - 2. * C4 * delta * s) +
                            pow(mrho, 4) * (-2. + 1. * delta - 8. * C4 * s +
                                            4. * C4 * delta * s))) *
          log(fabs(-pow(mpion, 2) - pow(mrho, 2) + s + t1))) /
             (pow(mrho, 4) - 1. * pow(mrho, 2) * s) -
         (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
          (eta2 * pow(mpion, 2) *
               (pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                pow(ma1, 2) * (-1. * pow(mrho, 2) + 1. * s)) +
           eta1 * (pow(ma1, 2) * (-0.5 * pow(mrho, 4) +
                                  pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                                  0.5 * pow(mrho, 2) * s) +
                   pow(mpion, 2) *
                       (0.5 * pow(mrho, 4) - 0.5 * pow(mrho, 2) * s +
                        pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s)))) *
          log(fabs(-pow(mpion, 2) - pow(mrho, 2) + s + t1))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
              2. * pow(ma1, 2) * pow(mpion, 2) + pow(mpion, 4)) -
         (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * (0.5 * pow(Gammaa1, 4) * pow(ma1, 4) - 0.5 * pow(ma1, 8) +
                   0.5 * pow(mpion, 8) +
                   0.5 * pow(ma1, 4) * pow(mpion, 2) * pow(mrho, 2) -
                   0.5 * pow(mpion, 6) * pow(mrho, 2) +
                   pow(Gammaa1, 2) * pow(ma1, 2) *
                       (pow(mpion, 2) *
                            (1. * pow(mpion, 2) + 1.5 * pow(mrho, 2) - 2. * s) +
                        pow(ma1, 2) * (-1. * pow(mpion, 2) +
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
                   2. * pow(mrho, 2) * s + pow(s, 2) - 4. * pow(mpion, 2) * t1 -
                   2. * pow(mrho, 2) * t1 + 2. * s * t1 + pow(t1, 2) +
                   pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                  2. * s + 2. * t1)))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) -
              2. * pow(ma1, 2) * pow(mpion, 2) + pow(mpion, 4)) +
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
                pow(ma1, 6) *
                    (1. * pow(mpion, 2) + 1. * pow(mrho, 2) - 0.5 * s) +
                pow(ma1, 4) * (-0.5 * pow(mrho, 4) +
                               (-0.5 * pow(mpion, 2) + 0.5 * s) * s) +
                pow(Gammaa1, 2) * pow(ma1, 2) *
                    (1. * pow(mpion, 4) + 0.5 * pow(mrho, 4) +
                     pow(mpion, 2) * (2. * pow(mrho, 2) - 1.5 * s) -
                     1. * pow(mrho, 2) * s + 0.5 * pow(s, 2) +
                     pow(ma1, 2) *
                         (-1. * pow(mpion, 2) - 1. * pow(mrho, 2) + 1.5 * s)) +
                pow(mpion, 2) *
                    (0.5 * pow(mpion, 6) + 0.5 * pow(mpion, 4) * s +
                     pow(mpion, 2) * (-0.5 * pow(mrho, 4) +
                                      2. * pow(mrho, 2) * s - 1.5 * pow(s, 2)) +
                     s * (0.5 * pow(mrho, 4) - 1. * pow(mrho, 2) * s +
                          0.5 * pow(s, 2))) +
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
                   2. * pow(mrho, 2) * s + pow(s, 2) - 4. * pow(mpion, 2) * t1 -
                   2. * pow(mrho, 2) * t1 + 2. * s * t1 + pow(t1, 2) +
                   pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                  2. * s + 2. * t1)))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) + pow(mpion, 4) +
              2. * pow(mpion, 2) * pow(mrho, 2) + pow(mrho, 4) -
              2. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s + pow(s, 2) +
              pow(ma1, 2) *
                  (-2. * pow(mpion, 2) - 2. * pow(mrho, 2) + 2. * s)) +
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
                   2. * pow(mrho, 2) * s + pow(s, 2) - 4. * pow(mpion, 2) * t1 -
                   2. * pow(mrho, 2) * t1 + 2. * s * t1 + pow(t1, 2) +
                   pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                  2. * s + 2. * t1)))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) + 4. * pow(ma1, 4) +
              4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
              pow(mrho, 4) - 4. * pow(mpion, 2) * s - 2. * pow(mrho, 2) * s +
              pow(s, 2) +
              pow(ma1, 2) *
                  (-8. * pow(mpion, 2) - 4. * pow(mrho, 2) + 4. * s)) -
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
             log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) +
                      4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                      1. * pow(mrho, 4) - 4. * pow(mpion, 2) * s -
                      2. * pow(mrho, 2) * s + pow(s, 2) -
                      4. * pow(mpion, 2) * t1 - 2. * pow(mrho, 2) * t1 +
                      2. * s * t1 + pow(t1, 2) +
                      pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                     2. * s + 2. * t1))) +
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
          log(fabs(pow(Gammaa1, 2) * pow(ma1, 2) + pow(ma1, 4) +
                   4. * pow(mpion, 4) + 4. * pow(mpion, 2) * pow(mrho, 2) +
                   1. * pow(mrho, 4) - 4. * pow(mpion, 2) * s -
                   2. * pow(mrho, 2) * s + pow(s, 2) - 4. * pow(mpion, 2) * t1 -
                   2. * pow(mrho, 2) * t1 + 2. * s * t1 + pow(t1, 2) +
                   pow(ma1, 2) * (-4. * pow(mpion, 2) - 2. * pow(mrho, 2) +
                                  2. * s + 2. * t1)))) /
             pow(mrho, 2)))) /
      (16. * M_PI * s * (-4 * pow(mpion, 2) + s));

  return xs;
}

float_t PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_pi_rho(
    const float_t m1, const float_t m2, const float_t m3, const float_t t1,
    const float_t t2, const float_t s, const float_t mpion,
    const float_t mrho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const float_t xs =
      to_mb *
      (-(pow(Const, 2) * pow(ghat, 4) *
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
              (1. * pow(ma1, 2) - 1. * t2) +
          (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
           (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
              (1. * pow(mpion, 2) - 1. * t2) -
          (0.25 * pow(-2. + delta, 2) * pow(mpion, 2) * t2) / pow(mrho, 2) -
          0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
              (eta2 * (-1. * pow(ma1, 2) + pow(mrho, 2) - 2. * s) +
               eta1 * (2. * pow(mpion, 2) + s)) *
              t2 +
          (0.5 * pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) *
           (4. * pow(mpion, 4) * pow(mrho, 2) + 1. * pow(mrho, 6) -
            3.5 * pow(mrho, 4) * s + 0.5 * pow(s, 3) +
            pow(mpion, 2) * (10. * pow(mrho, 4) - 2. * pow(s, 2))) *
           t2) /
              (pow(mrho, 6) * pow(pow(mrho, 2) - 1. * s, 2)) -
          (0.25 * (eta1 - 1. * eta2) * (1. * pow(mrho, 2) - 0.5 * delta * s) *
           (eta2 *
                (-2. * pow(ma1, 4) - 6. * pow(mpion, 4) + 1.5 * pow(mrho, 4) +
                 pow(ma1, 2) *
                     (6. * pow(mpion, 2) - 2. * pow(mrho, 2) - 2. * s) -
                 1. * pow(mrho, 2) * s - 0.5 * pow(s, 2) +
                 pow(mpion, 2) * (2. * pow(mrho, 2) + 2. * s)) +
            eta1 * (2. * pow(ma1, 4) + 6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                    pow(mpion, 2) * (8. * pow(mrho, 2) - 4. * s) -
                    4. * pow(mrho, 2) * s + 1. * pow(s, 2) +
                    pow(ma1, 2) *
                        (-6. * pow(mpion, 2) - 3. * pow(mrho, 2) + 3. * s))) *
           t2) /
              (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
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
              t2 -
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
           t2) /
              (pow(mrho, 6) - 1. * pow(mrho, 4) * s) +
          (0.25 * (1. * eta1 - 1. * eta2) *
           (pow(mrho, 2) *
                (eta1 * (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
                         pow(ma1, 2) * (1. - 2. * C4 * pow(mrho, 2)) +
                         pow(mpion, 2) * (-2. + 4. * C4 * pow(mrho, 2)) -
                         2. * C4 * pow(s, 2)) +
                 eta2 *
                     (-1.5 * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
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
           t2) /
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
           t2) /
              pow(mrho, 6) -
          (3. * (1. * pow(mrho, 2) - 0.5 * delta * s) *
           (delta * (0.666667 * pow(mpion, 4) * pow(mrho, 2) +
                     0.166667 * pow(mrho, 6) - 0.541667 * pow(mrho, 4) * s -
                     0.0833333 * pow(mrho, 2) * pow(s, 2) + 0.125 * pow(s, 3) +
                     pow(mpion, 2) * (1.66667 * pow(mrho, 4) +
                                      0.0833333 * pow(mrho, 2) * s -
                                      0.416667 * pow(s, 2))) +
            pow(mrho, 2) *
                (1. * C4 * pow(mrho, 6) - 0.0833333 * pow(s, 2) +
                 pow(mrho, 4) * (-0.416667 - 1.33333 * C4 * s) +
                 pow(mrho, 2) * s * (0.5 + 0.333333 * C4 * s) +
                 pow(mpion, 2) *
                     (2. * C4 * pow(mrho, 4) + 0.166667 * s +
                      pow(mrho, 2) * (-0.833333 - 0.666667 * C4 * s)))) *
           t2) /
              (pow(mrho, 8) - 1. * pow(mrho, 6) * s) -
          1. * C4 * pow(t2, 2) - 1. * C4 * delta * pow(t2, 2) +
          0.0625 * (-2. + delta) * (eta1 - 1. * eta2) * eta2 * pow(t2, 2) -
          (0.5 * pow(delta, 2) * pow(mpion, 2) * pow(t2, 2)) / pow(mrho, 4) +
          (0.25 * pow(t2, 2)) / pow(mrho, 2) +
          (0.5 * delta * pow(t2, 2)) / pow(mrho, 2) -
          (0.25 * pow(delta, 2) * pow(t2, 2)) / pow(mrho, 2) -
          (0.25 * delta * s * pow(t2, 2)) / pow(mrho, 4) +
          (0.25 * pow(delta, 2) * s * pow(t2, 2)) / pow(mrho, 4) +
          (0.5 * C4 * delta * s * pow(t2, 2)) / pow(mrho, 2) +
          (0.0625 * pow(delta, 2) * pow(s, 2) * pow(t2, 2)) / pow(mrho, 6) -
          (1. * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s) *
           pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) * pow(t2, 2)) /
              (pow(mrho, 4) * pow(pow(mrho, 2) - 1. * s, 2)) +
          (0.375 * (eta1 - 1. * eta2) *
           (eta1 * (-0.6666666666666666 * pow(ma1, 2) + 2. * pow(mpion, 2) +
                    1. * pow(mrho, 2) - 1. * s) +
            eta2 *
                (0.6666666666666666 * pow(ma1, 2) - 2. * pow(mpion, 2) +
                 0.6666666666666666 * pow(mrho, 2) + 0.6666666666666666 * s)) *
           (1. * pow(mrho, 2) - 0.5 * delta * s) * pow(t2, 2)) /
              (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
          0.03125 * pow(eta1 - 1. * eta2, 3) *
              (eta2 * (-1. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                       1. * pow(mrho, 2) - 1. * s) +
               eta1 *
                   (pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
              pow(t2, 2) +
          (3. * (1. * pow(mrho, 2) - 0.5 * delta * s) *
           (1. * C4 * pow(mrho, 6) + 0.0833335 * pow(mrho, 2) * s +
            pow(mrho, 4) * (-0.416667 - 0.333334 * C4 * s) +
            delta * (0.666665 * pow(mpion, 2) * pow(mrho, 2) +
                     0.333334 * pow(mrho, 4) - 0.291667 * pow(mrho, 2) * s -
                     0.0416667 * pow(s, 2))) *
           pow(t2, 2)) /
              (pow(mrho, 8) - 1. * pow(mrho, 6) * s) +
          (0.125 * (1. * eta1 - 1. * eta2) *
           (pow(mrho, 2) * (eta1 * (1. - 2. * C4 * pow(mrho, 2)) +
                            eta2 * (-1. + 2. * C4 * pow(mrho, 2))) +
            delta * (eta2 * (-1. * pow(ma1, 2) + 3. * pow(mpion, 2) -
                             1. * pow(mrho, 2) - 1. * s) +
                     eta1 * (1. * pow(ma1, 2) - 3. * pow(mpion, 2) -
                             1.5 * pow(mrho, 2) + 1.5 * s))) *
           pow(t2, 2)) /
              pow(mrho, 2) +
          0.0104167 * pow(eta1 - 1. * eta2, 4) * pow(t2, 3) +
          (0.166667 * pow(delta, 2) * pow(t2, 3)) / pow(mrho, 4) +
          (0.0833333 * delta * pow(1. * eta1 - 1. * eta2, 2) * pow(t2, 3)) /
              pow(mrho, 2) +
          (0.666667 * pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) *
           pow(t2, 3)) /
              (pow(mrho, 4) * pow(pow(mrho, 2) - 1. * s, 2)) -
          (0.166667 * pow(1. * eta1 - 1. * eta2, 2) *
           (1. * pow(mrho, 2) - 0.5 * delta * s) * pow(t2, 3)) /
              (pow(mrho, 4) - 1. * pow(mrho, 2) * s) +
          (0.333334 * delta * (-2. * pow(mrho, 2) + 1. * delta * s) *
           pow(t2, 3)) /
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
              (1. * pow(ma1, 2) - 1. * t1) -
          (1. * pow(-2. + delta, 2) * pow(mpion, 2) *
           (1. * pow(mpion, 2) - 0.25 * pow(mrho, 2))) /
              (1. * pow(mpion, 2) - 1. * t1) +
          (0.25 * pow(-2. + delta, 2) * pow(mpion, 2) * t1) / pow(mrho, 2) +
          0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
              (eta2 * (-1. * pow(ma1, 2) + pow(mrho, 2) - 2. * s) +
               eta1 * (2. * pow(mpion, 2) + s)) *
              t1 -
          (0.5 * pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) *
           (4. * pow(mpion, 4) * pow(mrho, 2) + 1. * pow(mrho, 6) -
            3.5 * pow(mrho, 4) * s + 0.5 * pow(s, 3) +
            pow(mpion, 2) * (10. * pow(mrho, 4) - 2. * pow(s, 2))) *
           t1) /
              (pow(mrho, 6) * pow(pow(mrho, 2) - 1. * s, 2)) +
          (0.25 * (eta1 - 1. * eta2) * (1. * pow(mrho, 2) - 0.5 * delta * s) *
           (eta2 *
                (-2. * pow(ma1, 4) - 6. * pow(mpion, 4) + 1.5 * pow(mrho, 4) +
                 pow(ma1, 2) *
                     (6. * pow(mpion, 2) - 2. * pow(mrho, 2) - 2. * s) -
                 1. * pow(mrho, 2) * s - 0.5 * pow(s, 2) +
                 pow(mpion, 2) * (2. * pow(mrho, 2) + 2. * s)) +
            eta1 * (2. * pow(ma1, 4) + 6. * pow(mpion, 4) + 1. * pow(mrho, 4) +
                    pow(mpion, 2) * (8. * pow(mrho, 2) - 4. * s) -
                    4. * pow(mrho, 2) * s + 1. * pow(s, 2) +
                    pow(ma1, 2) *
                        (-6. * pow(mpion, 2) - 3. * pow(mrho, 2) + 3. * s))) *
           t1) /
              (pow(mrho, 4) - 1. * pow(mrho, 2) * s) -
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
              t1 +
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
           t1) /
              (pow(mrho, 6) - 1. * pow(mrho, 4) * s) -
          (0.25 * (1. * eta1 - 1. * eta2) *
           (pow(mrho, 2) *
                (eta1 * (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
                         pow(ma1, 2) * (1. - 2. * C4 * pow(mrho, 2)) +
                         pow(mpion, 2) * (-2. + 4. * C4 * pow(mrho, 2)) -
                         2. * C4 * pow(s, 2)) +
                 eta2 *
                     (-1.5 * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
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
           t1) /
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
           t1) /
              pow(mrho, 6) +
          (3. * (1. * pow(mrho, 2) - 0.5 * delta * s) *
           (delta * (0.666667 * pow(mpion, 4) * pow(mrho, 2) +
                     0.166667 * pow(mrho, 6) - 0.541667 * pow(mrho, 4) * s -
                     0.0833333 * pow(mrho, 2) * pow(s, 2) + 0.125 * pow(s, 3) +
                     pow(mpion, 2) * (1.66667 * pow(mrho, 4) +
                                      0.0833333 * pow(mrho, 2) * s -
                                      0.416667 * pow(s, 2))) +
            pow(mrho, 2) *
                (1. * C4 * pow(mrho, 6) - 0.0833333 * pow(s, 2) +
                 pow(mrho, 4) * (-0.416667 - 1.33333 * C4 * s) +
                 pow(mrho, 2) * s * (0.5 + 0.333333 * C4 * s) +
                 pow(mpion, 2) *
                     (2. * C4 * pow(mrho, 4) + 0.166667 * s +
                      pow(mrho, 2) * (-0.833333 - 0.666667 * C4 * s)))) *
           t1) /
              (pow(mrho, 8) - 1. * pow(mrho, 6) * s) +
          1. * C4 * pow(t1, 2) + 1. * C4 * delta * pow(t1, 2) -
          0.0625 * (-2. + delta) * (eta1 - 1. * eta2) * eta2 * pow(t1, 2) +
          (0.5 * pow(delta, 2) * pow(mpion, 2) * pow(t1, 2)) / pow(mrho, 4) -
          (0.25 * pow(t1, 2)) / pow(mrho, 2) -
          (0.5 * delta * pow(t1, 2)) / pow(mrho, 2) +
          (0.25 * pow(delta, 2) * pow(t1, 2)) / pow(mrho, 2) +
          (0.25 * delta * s * pow(t1, 2)) / pow(mrho, 4) -
          (0.25 * pow(delta, 2) * s * pow(t1, 2)) / pow(mrho, 4) -
          (0.5 * C4 * delta * s * pow(t1, 2)) / pow(mrho, 2) -
          (0.0625 * pow(delta, 2) * pow(s, 2) * pow(t1, 2)) / pow(mrho, 6) +
          (1. * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s) *
           pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) * pow(t1, 2)) /
              (pow(mrho, 4) * pow(pow(mrho, 2) - 1. * s, 2)) -
          (0.375 * (eta1 - 1. * eta2) *
           (eta1 * (-0.6666666666666666 * pow(ma1, 2) + 2. * pow(mpion, 2) +
                    1. * pow(mrho, 2) - 1. * s) +
            eta2 *
                (0.6666666666666666 * pow(ma1, 2) - 2. * pow(mpion, 2) +
                 0.6666666666666666 * pow(mrho, 2) + 0.6666666666666666 * s)) *
           (1. * pow(mrho, 2) - 0.5 * delta * s) * pow(t1, 2)) /
              (pow(mrho, 4) - 1. * pow(mrho, 2) * s) -
          0.03125 * pow(eta1 - 1. * eta2, 3) *
              (eta2 * (-1. * pow(ma1, 2) + 2. * pow(mpion, 2) -
                       1. * pow(mrho, 2) - 1. * s) +
               eta1 *
                   (pow(ma1, 2) - 2. * pow(mpion, 2) - 1. * pow(mrho, 2) + s)) *
              pow(t1, 2) -
          (3. * (1. * pow(mrho, 2) - 0.5 * delta * s) *
           (1. * C4 * pow(mrho, 6) + 0.0833335 * pow(mrho, 2) * s +
            pow(mrho, 4) * (-0.416667 - 0.333334 * C4 * s) +
            delta * (0.666665 * pow(mpion, 2) * pow(mrho, 2) +
                     0.333334 * pow(mrho, 4) - 0.291667 * pow(mrho, 2) * s -
                     0.0416667 * pow(s, 2))) *
           pow(t1, 2)) /
              (pow(mrho, 8) - 1. * pow(mrho, 6) * s) -
          (0.125 * (1. * eta1 - 1. * eta2) *
           (pow(mrho, 2) * (eta1 * (1. - 2. * C4 * pow(mrho, 2)) +
                            eta2 * (-1. + 2. * C4 * pow(mrho, 2))) +
            delta * (eta2 * (-1. * pow(ma1, 2) + 3. * pow(mpion, 2) -
                             1. * pow(mrho, 2) - 1. * s) +
                     eta1 * (1. * pow(ma1, 2) - 3. * pow(mpion, 2) -
                             1.5 * pow(mrho, 2) + 1.5 * s))) *
           pow(t1, 2)) /
              pow(mrho, 2) -
          0.0104167 * pow(eta1 - 1. * eta2, 4) * pow(t1, 3) -
          (0.166667 * pow(delta, 2) * pow(t1, 3)) / pow(mrho, 4) -
          (0.0833333 * delta * pow(1. * eta1 - 1. * eta2, 2) * pow(t1, 3)) /
              pow(mrho, 2) -
          (0.666667 * pow(1. * pow(mrho, 2) - 0.5 * delta * s, 2) *
           pow(t1, 3)) /
              (pow(mrho, 4) * pow(pow(mrho, 2) - 1. * s, 2)) +
          (0.166667 * pow(1. * eta1 - 1. * eta2, 2) *
           (1. * pow(mrho, 2) - 0.5 * delta * s) * pow(t1, 3)) /
              (pow(mrho, 4) - 1. * pow(mrho, 2) * s) -
          (0.333334 * delta * (-2. * pow(mrho, 2) + 1. * delta * s) *
           pow(t1, 3)) /
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
              log(fabs(-pow(ma1, 2) + t2)) -
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
           log(fabs(-pow(ma1, 2) + t2))) /
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
            eta2 * (-4. * pow(ma1, 6) +
                    pow(ma1, 4) *
                        (12. * pow(mpion, 2) - 4. * pow(mrho, 2) - 4. * s) +
                    pow(mpion, 2) * (4. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                                     2. * pow(mrho, 2) * s - 1. * pow(s, 2)) +
                    pow(ma1, 2) *
                        (-12. * pow(mpion, 4) + 3. * pow(mrho, 4) -
                         2. * pow(mrho, 2) * s - 1. * pow(s, 2) +
                         pow(mpion, 2) * (4. * pow(mrho, 2) + 4. * s)))) *
           log(fabs(-pow(ma1, 2) + t2))) /
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
                 eta1 *
                     (pow(ma1, 4) * (1. - 2. * C4 * pow(mrho, 2)) +
                      pow(mpion, 4) * (1. - 2. * C4 * pow(mrho, 2)) +
                      (-0.5 * pow(mrho, 2) + 0.5 * s) * s +
                      pow(ma1, 2) *
                          (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
                           pow(mpion, 2) * (-2. + 4. * C4 * pow(mrho, 2)) -
                           2. * C4 * pow(s, 2)) +
                      pow(mpion, 2) * (-4. * C4 * pow(mrho, 4) - 1. * s +
                                       pow(mrho, 2) * (2. + 4. * C4 * s))))) *
           log(fabs(-pow(ma1, 2) + t2))) /
              pow(mrho, 2) +
          0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
              log(fabs(-pow(mpion, 2) + t2)) +
          (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
           (eta2 * pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s) +
            eta1 * (-0.5 * pow(mrho, 4) +
                    pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                    0.5 * pow(mrho, 2) * s)) *
           log(fabs(-pow(mpion, 2) + t2))) /
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
           log(fabs(-pow(mpion, 2) + t2))) /
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
              log(fabs(-pow(ma1, 2) + t1)) +
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
           log(fabs(-pow(ma1, 2) + t1))) /
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
            eta2 * (-4. * pow(ma1, 6) +
                    pow(ma1, 4) *
                        (12. * pow(mpion, 2) - 4. * pow(mrho, 2) - 4. * s) +
                    pow(mpion, 2) * (4. * pow(mpion, 4) - 1. * pow(mrho, 4) +
                                     2. * pow(mrho, 2) * s - 1. * pow(s, 2)) +
                    pow(ma1, 2) *
                        (-12. * pow(mpion, 4) + 3. * pow(mrho, 4) -
                         2. * pow(mrho, 2) * s - 1. * pow(s, 2) +
                         pow(mpion, 2) * (4. * pow(mrho, 2) + 4. * s)))) *
           log(fabs(-pow(ma1, 2) + t1))) /
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
                 eta1 *
                     (pow(ma1, 4) * (1. - 2. * C4 * pow(mrho, 2)) +
                      pow(mpion, 4) * (1. - 2. * C4 * pow(mrho, 2)) +
                      (-0.5 * pow(mrho, 2) + 0.5 * s) * s +
                      pow(ma1, 2) *
                          (-1. * pow(mrho, 2) + 2. * C4 * pow(mrho, 4) +
                           pow(mpion, 2) * (-2. + 4. * C4 * pow(mrho, 2)) -
                           2. * C4 * pow(s, 2)) +
                      pow(mpion, 2) * (-4. * C4 * pow(mrho, 4) - 1. * s +
                                       pow(mrho, 2) * (2. + 4. * C4 * s))))) *
           log(fabs(-pow(ma1, 2) + t1))) /
              pow(mrho, 2) -
          0.5 * pow(-2. + delta, 2) * pow(mpion, 2) *
              log(fabs(-pow(mpion, 2) + t1)) -
          (0.5 * (-2. + delta) * (eta1 - 1. * eta2) * pow(mpion, 2) *
           (eta2 * pow(mpion, 2) * (-1. * pow(mrho, 2) + 1. * s) +
            eta1 * (-0.5 * pow(mrho, 4) +
                    pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s) +
                    0.5 * pow(mrho, 2) * s)) *
           log(fabs(-pow(mpion, 2) + t1))) /
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
           log(fabs(-pow(mpion, 2) + t1))) /
              (pow(mrho, 4) - 1. * pow(mrho, 2) * s))) /
       (16. * Pi * (4 * pow(mpion, 2) - s) * s));

  return xs;
}
