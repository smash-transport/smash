/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/photoncrosssections.h"

using namespace Smash;

// template class PhotonCrossSection<ComputationMethod::Analytic>;
// template class PhotonCrossSection<ComputationMethod::Lookup>;

constexpr double PhotonCrossSection<ComputationMethod::Analytic>::m_pion_;


double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho0_pi0(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;

  const double mpion = m_pion_, mrho = m_rho;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_rho, m_pion_, 0.0);

  const double t1 = t_mandelstam[1];
  const double t2 = t_mandelstam[0];

  const double xs =
          to_mb*(pow(Const,2)*pow(g_POR,4)*((pow(pow(momega,2) - s,2)*(pow(mpion,8) - 2*pow(mpion,6)*pow(mrho,2) + pow(mpion,4)*(pow(mrho,4) + 4*pow(momega,4) - 2*pow(momega,2)*s) +
                               pow(momega,4)*(pow(mrho,4) + pow(momega,4) + 2*pow(momega,2)*s + 2*pow(s,2) - 2*pow(mrho,2)*(pow(momega,2) + s)) -
                               2*pow(mpion,2)*pow(momega,2)*(pow(mrho,4) + 2*pow(momega,2)*(pow(momega,2) + s) - pow(mrho,2)*(2*pow(momega,2) + s))))/(pow(momega,2) - t2) +
                          (pow(mpion,8) - 2*pow(mpion,6)*pow(mrho,2) + 3*pow(momega,8) - 4*pow(momega,6)*s - 7*pow(momega,4)*pow(s,2) + 4*pow(momega,2)*pow(s,3) + 5*pow(s,4) +
                             pow(mrho,4)*(pow(momega,4) - 2*pow(momega,2)*s + 2*pow(s,2)) + pow(mrho,2)*(-4*pow(momega,6) + 8*pow(momega,4)*s - 6*pow(s,3)) -
                             2*pow(mpion,2)*(4*pow(momega,6) - 2*pow(mrho,2)*pow(pow(momega,2) - 2*s,2) + pow(mrho,4)*s - 10*pow(momega,4)*s + 8*pow(s,3)) +
                             pow(mpion,4)*(pow(mrho,4) + 2*pow(mrho,2)*(pow(momega,2) - s) + 4*(pow(momega,4) - 3*pow(momega,2)*s + 3*pow(s,2))))*t2 -
                          2*pow(mpion,2)*pow(momega,4)*pow(t2,2) - pow(mrho,2)*pow(momega,4)*pow(t2,2) + pow(momega,6)*pow(t2,2) - pow(mpion,4)*s*pow(t2,2) +
                          pow(mpion,2)*pow(mrho,2)*s*pow(t2,2) + 8*pow(mpion,2)*pow(momega,2)*s*pow(t2,2) + 3*pow(mrho,2)*pow(momega,2)*s*pow(t2,2) -
                          2*pow(momega,4)*s*pow(t2,2) - 8*pow(mpion,2)*pow(s,2)*pow(t2,2) - 3*pow(mrho,2)*pow(s,2)*pow(t2,2) - 3*pow(momega,2)*pow(s,2)*pow(t2,2) +
                          5*pow(s,3)*pow(t2,2) + ((pow(momega,4) - 4*pow(momega,2)*s + 5*pow(s,2))*pow(t2,3))/3. -
                          (pow(pow(momega,2) - s,2)*(pow(mpion,8) - 2*pow(mpion,6)*pow(mrho,2) + pow(mpion,4)*(pow(mrho,4) + 4*pow(momega,4) - 2*pow(momega,2)*s) +
                               pow(momega,4)*(pow(mrho,4) + pow(momega,4) + 2*pow(momega,2)*s + 2*pow(s,2) - 2*pow(mrho,2)*(pow(momega,2) + s)) -
                               2*pow(mpion,2)*pow(momega,2)*(pow(mrho,4) + 2*pow(momega,2)*(pow(momega,2) + s) - pow(mrho,2)*(2*pow(momega,2) + s))))/(pow(momega,2) - t1) -
                          (pow(mpion,8) - 2*pow(mpion,6)*pow(mrho,2) + 3*pow(momega,8) - 4*pow(momega,6)*s - 7*pow(momega,4)*pow(s,2) + 4*pow(momega,2)*pow(s,3) + 5*pow(s,4) +
                             pow(mrho,4)*(pow(momega,4) - 2*pow(momega,2)*s + 2*pow(s,2)) + pow(mrho,2)*(-4*pow(momega,6) + 8*pow(momega,4)*s - 6*pow(s,3)) -
                             2*pow(mpion,2)*(4*pow(momega,6) - 2*pow(mrho,2)*pow(pow(momega,2) - 2*s,2) + pow(mrho,4)*s - 10*pow(momega,4)*s + 8*pow(s,3)) +
                             pow(mpion,4)*(pow(mrho,4) + 2*pow(mrho,2)*(pow(momega,2) - s) + 4*(pow(momega,4) - 3*pow(momega,2)*s + 3*pow(s,2))))*t1 +
                          2*pow(mpion,2)*pow(momega,4)*pow(t1,2) + pow(mrho,2)*pow(momega,4)*pow(t1,2) - pow(momega,6)*pow(t1,2) + pow(mpion,4)*s*pow(t1,2) -
                          pow(mpion,2)*pow(mrho,2)*s*pow(t1,2) - 8*pow(mpion,2)*pow(momega,2)*s*pow(t1,2) - 3*pow(mrho,2)*pow(momega,2)*s*pow(t1,2) +
                          2*pow(momega,4)*s*pow(t1,2) + 8*pow(mpion,2)*pow(s,2)*pow(t1,2) + 3*pow(mrho,2)*pow(s,2)*pow(t1,2) + 3*pow(momega,2)*pow(s,2)*pow(t1,2) -
                          5*pow(s,3)*pow(t1,2) - ((pow(momega,4) - 4*pow(momega,2)*s + 5*pow(s,2))*pow(t1,3))/3. +
                          2*(pow(momega,2) - s)*(-pow(mpion,8) + pow(mpion,4)*(4*pow(momega,4) - 7*pow(momega,2)*s + pow(s,2) + pow(mrho,2)*(pow(momega,2) + s)) +
                             pow(mpion,2)*(-6*pow(momega,6) + 6*pow(momega,4)*s + 8*pow(momega,2)*pow(s,2) + pow(mrho,4)*(-pow(momega,2) + s) +
                                pow(mrho,2)*(4*pow(momega,4) - 7*pow(momega,2)*s - pow(s,2))) +
                             pow(momega,2)*(2*pow(momega,6) + pow(mrho,4)*(pow(momega,2) - s) - 4*pow(momega,2)*pow(s,2) - 3*pow(s,3) +
                                pow(mrho,2)*(-3*pow(momega,4) + 2*pow(momega,2)*s + 3*pow(s,2))))*log((-pow(momega,2) + t2)/(-pow(momega,2) + t1))))/
                      (3.0*128.*Pi*pow(pow(momega,2) - s,2)*(pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s)));

  return xs;
}



double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho_pi(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const double &mpion = m_pion_, &mrho = m_rho;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_rho, m_pion_, 0.);
  const double &t1 = t_mandelstam[1];
  const double &t2 = t_mandelstam[0];

  const double xs =
          to_mb*1/3.0*(0.0024868*pow(Const,2)*pow(g_POR,4)*((pow(momega,8) + pow(mpion,4)*pow(pow(mpion,2) - pow(mrho,2),2) -
                           2*pow(momega,6)*(2*pow(mpion,2) + pow(mrho,2) - s) -
                           2*pow(momega,2)*pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
                           pow(momega,4)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))/
                         (pow(momega,2) - t2) + 3*pow(momega,4)*t2 - 8*pow(momega,2)*pow(mpion,2)*t2 + 4*pow(mpion,4)*t2 -
                        4*pow(momega,2)*pow(mrho,2)*t2 + 4*pow(mpion,2)*pow(mrho,2)*t2 + pow(mrho,4)*t2 + 4*pow(momega,2)*s*t2 -
                        4*pow(mpion,2)*s*t2 - 2*pow(mrho,2)*s*t2 + 2*pow(s,2)*t2 + pow(momega,2)*pow(t2,2) - 2*pow(mpion,2)*pow(t2,2) -
                        pow(mrho,2)*pow(t2,2) + s*pow(t2,2) + pow(t2,3)/3. -
                        (pow(momega,8) + pow(mpion,4)*pow(pow(mpion,2) - pow(mrho,2),2) - 2*pow(momega,6)*(2*pow(mpion,2) + pow(mrho,2) - s) -
                           2*pow(momega,2)*pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
                           pow(momega,4)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))/
                         (pow(momega,2) - t1) - 3*pow(momega,4)*t1 + 8*pow(momega,2)*pow(mpion,2)*t1 - 4*pow(mpion,4)*t1 +
                        4*pow(momega,2)*pow(mrho,2)*t1 - 4*pow(mpion,2)*pow(mrho,2)*t1 - pow(mrho,4)*t1 - 4*pow(momega,2)*s*t1 +
                        4*pow(mpion,2)*s*t1 + 2*pow(mrho,2)*s*t1 - 2*pow(s,2)*t1 - pow(momega,2)*pow(t1,2) + 2*pow(mpion,2)*pow(t1,2) +
                        pow(mrho,2)*pow(t1,2) - s*pow(t1,2) - pow(t1,3)/3. +
                        2*(2*pow(momega,6) - 3*pow(momega,4)*(2*pow(mpion,2) + pow(mrho,2) - s) -
                           pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
                           pow(momega,2)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))*
                         log(fabs(-pow(momega,2) + t2)) - 2*(2*pow(momega,6) - 3*pow(momega,4)*(2*pow(mpion,2) + pow(mrho,2) - s) -
                           pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
                           pow(momega,2)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))*
                         log(fabs(-pow(momega,2) + t1))))/(pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s));

  return xs;
}



double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_rho_pi0(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;

  const double &mpion = m_pion_, &mrho = m_rho;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_rho, m_pion_, 0.);
  const double &t1 = t_mandelstam[1];
  const double &t2 = t_mandelstam[0];

  const double xs =
          to_mb*1/3.0*(0.0024867959858108648*pow(Const,2)*pow(g_POR,4)*(pow(mpion,8)*(t2 - t1) + pow(mpion,6)*pow(mrho,2)*(-2.*t2 + 2.*t1) +
                          pow(mpion,4)*(pow(mrho,4)*(t2 - t1) + s*(4.*s*t2 - pow(t2,2) - 4.*s*t1 + pow(t1,2))) +
                          pow(s,2)*(pow(s,2)*t2 + s*pow(t2,2) + 0.6666666666666666*pow(t2,3) + pow(mrho,4)*(t2 - t1) -
                             pow(s,2)*t1 - s*pow(t1,2) - 0.6666666666666666*pow(t1,3) +
                             pow(mrho,2)*(-2.*s*t2 - pow(t2,2) + 2.*s*t1 + pow(t1,2))) +
                          pow(mpion,2)*s*(pow(mrho,4)*(-2.*t2 + 2.*t1) + pow(mrho,2)*(4.*s*t2 + pow(t2,2) - 4.*s*t1 - pow(t1,2)) +
                             s*(-4.*s*t2 - 2.*pow(t2,2) + 4.*s*t1 + 2.*pow(t1,2)))))/
                      (pow(pow(momega,2) - s,2)*(pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-2.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s + pow(s,2)));
  return xs;
}





double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_rho0_pi(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const double m_pi = m_pion_;
  const double &mpion = m_pion_, &mrho = m_rho;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_rho, m_pion_, 0.);
  const double &t1 = t_mandelstam[1];
  const double &t2 = t_mandelstam[0];

  const double xs =
          to_mb*1/3.0*(pow(Const,2)*pow(ghat,4)*((pow(eta1 - eta2,2)*(-2*eta1*eta2*
                                (pow(ma1,8) + pow(m_pi,8) - pow(m_pi,4)*pow(mrho,4) - 2*pow(ma1,2)*pow(m_pi,2)*(pow(m_pi,2) - pow(mrho,2))*(pow(mrho,2) + s) +
                                  pow(ma1,6)*(-4*pow(m_pi,2) + 2*s) + pow(ma1,4)*
                                   (4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2))) +
                               pow(eta2,2)*(pow(ma1,8) + pow(m_pi,4)*pow(pow(m_pi,2) - pow(mrho,2),2) + 2*pow(ma1,6)*(-2*pow(m_pi,2) + pow(mrho,2) + s) +
                                  2*pow(ma1,2)*pow(m_pi,2)*(-pow(mrho,4) + pow(m_pi,2)*(2*pow(mrho,2) - s) + pow(mrho,2)*s) +
                                  pow(ma1,4)*(4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s))) +
                               pow(eta1,2)*(pow(ma1,8) + pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) - 2*pow(ma1,6)*(2*pow(m_pi,2) + pow(mrho,2) - s) -
                                  2*pow(m_pi,2)*pow(mrho,4)*s + pow(m_pi,4)*(3*pow(mrho,4) + 2*pow(mrho,2)*s) +
                                  pow(ma1,4)*(4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) - 4*pow(mrho,2)*s + 2*pow(s,2)) -
                                  2*pow(ma1,2)*(pow(mrho,2)*s*(-pow(mrho,2) + s) + pow(m_pi,4)*(3*pow(mrho,2) + s) + pow(m_pi,2)*(2*pow(mrho,4) - 3*pow(mrho,2)*s)))))
                            /((pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*(pow(ma1,2) - t2)) +
                          (8*pow(-2 + delta,2)*pow(m_pi,2)*(4*pow(m_pi,2) - pow(mrho,2)))/
                           ((pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*(pow(m_pi,2) - t2)) -
                          (8*pow(-2 + delta,2)*pow(m_pi,2)*t2)/(pow(mrho,2)*pow(pow(m_pi,2) - s,2)) -
                          (8*pow(-2 + delta,2)*pow(m_pi,2)*t2)/(pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (8*(-2 + delta)*(-8*C4*pow(mrho,4) + pow(m_pi,2)*(2 + delta - 8*C4*pow(mrho,2)) - (2 + 3*delta)*s + pow(mrho,2)*(-2 + 3*delta + 16*C4*s))*t2)/
                           (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (4*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(eta2*(pow(m_pi,2) + pow(mrho,2) - 2*s)*(pow(m_pi,2) + s) +
                               eta1*(-2*pow(m_pi,4) + pow(mrho,4) - 3*pow(mrho,2)*s + 2*pow(s,2) + pow(m_pi,2)*(pow(mrho,2) + s)))*t2)/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (pow(eta1 - eta2,2)*(pow(eta1,2)*(3*pow(ma1,4) + 4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) -
                                  4*pow(ma1,2)*(2*pow(m_pi,2) + pow(mrho,2) - s) - 4*pow(mrho,2)*s + 2*pow(s,2)) +
                               pow(eta2,2)*(3*pow(ma1,4) + 4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s) +
                                  4*pow(ma1,2)*(-2*pow(m_pi,2) + pow(mrho,2) + s)) -
                               2*eta1*eta2*(3*pow(ma1,4) + 4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2) +
                                  pow(ma1,2)*(-8*pow(m_pi,2) + 4*s)))*t2)/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) +
                          (8*(pow(delta,2)*(8*pow(m_pi,4) + 3*pow(mrho,4) + 4*pow(m_pi,2)*(3*pow(mrho,2) - 2*s) - 6*pow(mrho,2)*s + 2*pow(s,2)) +
                               4*pow(mrho,4)*(3 + 12*C4*(2*pow(m_pi,2) - s) + 8*pow(C4,2)*pow(-2*pow(m_pi,2) + s,2)) -
                               4*delta*pow(mrho,2)*(16*C4*pow(m_pi,4) + 2*pow(m_pi,2)*(3 + 6*C4*pow(mrho,2) - 8*C4*s) + pow(mrho,2)*(3 - 6*C4*s) + s*(-3 + 4*C4*s)))*t2)/
                           (pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (8*(-2 + delta)*(pow(m_pi,4)*(-2 + 3*delta - 8*C4*pow(mrho,2)) + (pow(mrho,2) - s)*((-2 + 3*delta)*s + pow(mrho,2)*(-2 + delta - 8*C4*s)) +
                               4*pow(m_pi,2)*(2*C4*pow(mrho,4) + delta*s - pow(mrho,2)*(-1 + delta + 4*C4*s)))*t2)/
                           (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (4*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(eta2*(pow(m_pi,2) + s)*(pow(m_pi,4) - pow(m_pi,2)*(pow(mrho,2) - 2*s) + (pow(mrho,2) - s)*s) +
                               eta1*(-4*pow(m_pi,6) + pow(pow(mrho,2) - s,2)*s + pow(m_pi,4)*(3*pow(mrho,2) + s) -
                                  pow(m_pi,2)*(pow(mrho,4) - pow(mrho,2)*s + 2*pow(s,2))))*t2)/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,2) - s)*
                             (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (pow(eta1 - eta2,2)*(-2*eta1*eta2*(pow(m_pi,8) - pow(mrho,4)*pow(s,2) + pow(s,4) -
                                  pow(m_pi,4)*(pow(mrho,4) + 2*pow(mrho,2)*s - 4*pow(s,2)) + 2*pow(m_pi,2)*s*(pow(mrho,4) + pow(mrho,2)*s - 2*pow(s,2))) +
                               pow(eta2,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + pow(s,2)*pow(pow(mrho,2) + s,2) + pow(m_pi,4)*pow(pow(mrho,2) + 2*s,2) -
                                  2*pow(m_pi,2)*s*(pow(mrho,4) + 2*pow(mrho,2)*s + 2*pow(s,2))) +
                               pow(eta1,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) - 4*pow(m_pi,2)*pow(pow(mrho,2) - s,2)*s + pow(pow(mrho,2) - s,2)*pow(s,2) +
                                  pow(m_pi,4)*(3*pow(mrho,4) - 6*pow(mrho,2)*s + 4*pow(s,2))))*t2)/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (2*pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(pow(eta1,2)*(pow(ma1,4)*s + pow(m_pi,4)*(-3*pow(mrho,2) + 2*s) +
                                  s*(2*pow(mrho,4) - 3*pow(mrho,2)*s + pow(s,2)) - 2*pow(m_pi,2)*(pow(mrho,4) - 4*pow(mrho,2)*s + 2*pow(s,2)) +
                                  pow(ma1,2)*(2*pow(m_pi,2)*(pow(mrho,2) - 2*s) + 3*s*(-pow(mrho,2) + s))) -
                               2*eta1*eta2*(pow(ma1,4)*s + s*(2*pow(m_pi,4) + 4*pow(m_pi,2)*(pow(mrho,2) - s) + s*(-2*pow(mrho,2) + s)) +
                                  pow(ma1,2)*(pow(m_pi,2)*(pow(mrho,2) - 4*s) + s*(-2*pow(mrho,2) + 3*s))) +
                               pow(eta2,2)*(-4*pow(m_pi,2)*s*(pow(ma1,2) + pow(mrho,2) + s) + pow(m_pi,4)*(pow(mrho,2) + 2*s) +
                                  s*(pow(ma1,4) + s*(pow(mrho,2) + s) + pow(ma1,2)*(pow(mrho,2) + 3*s))))*t2)/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (4*(eta1 - eta2)*(pow(ma1,2) - s)*(eta1*(-4*pow(m_pi,4)*(6*C4*pow(mrho,4) + 2*delta*s + pow(mrho,2)*(1 - 2*delta - 8*C4*s)) +
                                  2*pow(m_pi,2)*(4*delta*pow(s,2) + pow(mrho,2)*s*(6 - 7*delta - 16*C4*s) + 2*pow(mrho,4)*(-2 + delta + 8*C4*s)) -
                                  (pow(mrho,2) - s)*s*(-2*delta*s + pow(mrho,2)*(-6 + 3*delta + 8*C4*s))) +
                               eta2*(delta*(2*pow(m_pi,4)*(pow(mrho,2) + 4*s) + pow(m_pi,2)*(2*pow(mrho,4) + pow(mrho,2)*s - 8*pow(s,2)) +
                                     s*(-2*pow(mrho,4) - pow(mrho,2)*s + 2*pow(s,2))) -
                                  2*pow(mrho,2)*(4*C4*pow(m_pi,4)*(pow(mrho,2) + 4*s) + pow(m_pi,2)*(s*(5 - 16*C4*s) + pow(mrho,2)*(2 - 8*C4*s)) +
                                     s*(s*(-3 + 4*C4*s) + pow(mrho,2)*(-2 + 4*C4*s)))))*t2)/
                           (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                             (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (8*(eta1 - eta2)*(delta*(eta1*(4*pow(m_pi,6) + pow(m_pi,4)*(7*pow(mrho,2) - 8*s) + pow(ma1,4)*(pow(m_pi,2) - s) -
                                     pow(ma1,2)*(2*pow(m_pi,2) + pow(mrho,2) - 2*s)*(2*pow(m_pi,2) - s) + pow(m_pi,2)*s*(-8*pow(mrho,2) + 5*s) +
                                     s*(pow(mrho,4) + pow(mrho,2)*s - pow(s,2))) +
                                  eta2*(-4*pow(m_pi,6) - pow(m_pi,4)*(pow(mrho,2) - 8*s) + pow(ma1,4)*(-pow(m_pi,2) + s) +
                                     pow(m_pi,2)*(2*pow(mrho,4) - 5*pow(s,2)) + s*(-2*pow(mrho,4) + pow(mrho,2)*s + pow(s,2)) +
                                     pow(ma1,2)*(4*pow(m_pi,4) - 6*pow(m_pi,2)*s + s*(pow(mrho,2) + 2*s)))) -
                               2*pow(mrho,2)*(eta1*(8*C4*pow(m_pi,6) + pow(m_pi,4)*(3 + 8*C4*(pow(mrho,2) - 2*s)) + 2*C4*pow(ma1,4)*(pow(m_pi,2) - s) +
                                     2*pow(m_pi,2)*s*(-1 - 6*C4*pow(mrho,2) + 5*C4*s) -
                                     pow(ma1,2)*(8*C4*pow(m_pi,4) + pow(m_pi,2)*(1 + 2*C4*(pow(mrho,2) - 6*s)) + 2*C4*s*(-pow(mrho,2) + 2*s)) +
                                     s*(-(s*(1 + 2*C4*s)) + pow(mrho,2)*(1 + 4*C4*s))) +
                                  eta2*(2*C4*pow(ma1,4)*(-pow(m_pi,2) + s) - (pow(m_pi,2) - s)*
                                      (8*C4*pow(m_pi,4) - 2*pow(mrho,2) + s + 2*C4*pow(s,2) + pow(m_pi,2)*(3 - 4*C4*(pow(mrho,2) + 2*s))) +
                                     pow(ma1,2)*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(mrho,2) + 2*s) + pow(m_pi,2)*(1 - 2*C4*(pow(mrho,2) + 6*s))))))*t2)/
                           (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (8*(-2 + delta)*(delta - 4*C4*pow(mrho,2))*pow(t2,2))/
                           (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (16*(-2 + delta)*(delta - 4*C4*pow(mrho,2))*s*pow(t2,2))/
                           (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (pow(eta1 - eta2,2)*(pow(eta1,2)*(pow(mrho,2) - s) + 2*eta1*eta2*s - pow(eta2,2)*s)*
                             (pow(m_pi,4) - pow(m_pi,2)*(pow(mrho,2) - 2*s) + (pow(mrho,2) - s)*s)*pow(t2,2))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (2*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(-(eta1*(pow(m_pi,2) + 2*pow(mrho,2) - 3*s)) - eta2*(pow(m_pi,2) + s))*pow(t2,2))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (2*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(pow(m_pi,2) + s)*(-2*eta2*s + eta1*(pow(m_pi,2) - pow(mrho,2) + s))*pow(t2,2))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,2) - s)*
                             (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (pow(eta1 - eta2,3)*(eta1*(pow(ma1,2) - 2*pow(m_pi,2) - pow(mrho,2) + s) - eta2*(pow(ma1,2) - 2*pow(m_pi,2) + pow(mrho,2) + s))*
                             pow(t2,2))/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) -
                          (8*(delta - 4*C4*pow(mrho,2))*(delta*(4*pow(m_pi,2) + 3*pow(mrho,2) - 2*s) - 2*pow(mrho,2)*(3 + 8*C4*pow(m_pi,2) - 4*C4*s))*pow(t2,2))/
                           (pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(pow(eta2,2)*s*(pow(ma1,2) - 4*pow(m_pi,2) + pow(mrho,2) + 3*s) +
                               pow(eta1,2)*(2*pow(m_pi,2)*(pow(mrho,2) - 2*s) + s*(pow(ma1,2) - 3*pow(mrho,2) + 3*s)) -
                               2*eta1*eta2*(pow(m_pi,2)*(pow(mrho,2) - 4*s) + s*(pow(ma1,2) - 2*pow(mrho,2) + 3*s)))*pow(t2,2))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (2*(eta1 - eta2)*(pow(ma1,2) - s)*(eta1*(4*delta*pow(s,2) - 2*pow(mrho,2)*s*(-2 + 3*delta + 8*C4*s) + pow(mrho,4)*(-2 + delta + 16*C4*s) -
                                  2*pow(m_pi,2)*(8*C4*pow(mrho,4) + 4*delta*s + pow(mrho,2)*(2 - 3*delta - 16*C4*s))) +
                               eta2*(pow(m_pi,2)*(8*delta*s + pow(mrho,2)*(-2 + delta - 32*C4*s)) + s*(-4*delta*s + pow(mrho,2)*(-2 + delta + 16*C4*s))))*pow(t2,2))/
                           (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                             (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (4*(eta1 - eta2)*(delta*(eta1*(pow(ma1,2)*(pow(m_pi,2) - s) - (2*pow(m_pi,2) + pow(mrho,2) - 2*s)*(2*pow(m_pi,2) - s)) +
                                  eta2*(4*pow(m_pi,4) - 6*pow(m_pi,2)*s + pow(ma1,2)*(-pow(m_pi,2) + s) + s*(pow(mrho,2) + 2*s))) +
                               2*pow(mrho,2)*(eta1*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(ma1,2) - pow(mrho,2) + 2*s) +
                                     pow(m_pi,2)*(1 - 2*C4*(pow(ma1,2) - pow(mrho,2) + 6*s))) -
                                  eta2*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(ma1,2) + pow(mrho,2) + 2*s) - pow(m_pi,2)*(-1 + 2*C4*(pow(ma1,2) + pow(mrho,2) + 6*s)))))*
                             pow(t2,2))/(pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (pow(eta1 - eta2,4)*pow(t2,3))/(3.*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (8*pow(eta1 - eta2,2)*(delta - 4*C4*pow(mrho,2))*pow(t2,3))/
                           (3.*pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (16*pow(delta - 4*C4*pow(mrho,2),2)*pow(t2,3))/
                           (3.*pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (4*(-2 + delta)*eta1*(eta1 - eta2)*(pow(ma1,2) - s)*pow(t2,3))/
                           (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (2*pow(eta1 - eta2,4)*(pow(ma1,2) - s)*s*pow(t2,3))/
                           (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (2*pow(eta1 - eta2,2)*s*(-2*eta1*eta2*s + pow(eta2,2)*s + pow(eta1,2)*(-pow(mrho,2) + s))*pow(t2,3))/
                           (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (4*(eta1 - eta2)*(pow(ma1,2) - s)*(2*eta2*(delta - 4*C4*pow(mrho,2))*s + eta1*(-2*delta*s + pow(mrho,2)*(-2 + delta + 8*C4*s)))*pow(t2,3))/
                           (3.*pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                             (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (pow(eta1 - eta2,2)*(-2*eta1*eta2*(pow(ma1,8) + pow(m_pi,8) - pow(m_pi,4)*pow(mrho,4) -
                                  2*pow(ma1,2)*pow(m_pi,2)*(pow(m_pi,2) - pow(mrho,2))*(pow(mrho,2) + s) + pow(ma1,6)*(-4*pow(m_pi,2) + 2*s) +
                                  pow(ma1,4)*(4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2))) +
                               pow(eta2,2)*(pow(ma1,8) + pow(m_pi,4)*pow(pow(m_pi,2) - pow(mrho,2),2) + 2*pow(ma1,6)*(-2*pow(m_pi,2) + pow(mrho,2) + s) +
                                  2*pow(ma1,2)*pow(m_pi,2)*(-pow(mrho,4) + pow(m_pi,2)*(2*pow(mrho,2) - s) + pow(mrho,2)*s) +
                                  pow(ma1,4)*(4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s))) +
                               pow(eta1,2)*(pow(ma1,8) + pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) - 2*pow(ma1,6)*(2*pow(m_pi,2) + pow(mrho,2) - s) -
                                  2*pow(m_pi,2)*pow(mrho,4)*s + pow(m_pi,4)*(3*pow(mrho,4) + 2*pow(mrho,2)*s) +
                                  pow(ma1,4)*(4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) - 4*pow(mrho,2)*s + 2*pow(s,2)) -
                                  2*pow(ma1,2)*(pow(mrho,2)*s*(-pow(mrho,2) + s) + pow(m_pi,4)*(3*pow(mrho,2) + s) + pow(m_pi,2)*(2*pow(mrho,4) - 3*pow(mrho,2)*s)))))
                            /((pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*(pow(ma1,2) - t1)) -
                          (8*pow(-2 + delta,2)*pow(m_pi,2)*(4*pow(m_pi,2) - pow(mrho,2)))/
                           ((pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))*(pow(m_pi,2) - t1)) +
                          (8*pow(-2 + delta,2)*pow(m_pi,2)*t1)/(pow(mrho,2)*pow(pow(m_pi,2) - s,2)) +
                          (8*pow(-2 + delta,2)*pow(m_pi,2)*t1)/(pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (8*(-2 + delta)*(-8*C4*pow(mrho,4) + pow(m_pi,2)*(2 + delta - 8*C4*pow(mrho,2)) - (2 + 3*delta)*s + pow(mrho,2)*(-2 + 3*delta + 16*C4*s))*t1)/
                           (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (4*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(eta2*(pow(m_pi,2) + pow(mrho,2) - 2*s)*(pow(m_pi,2) + s) +
                               eta1*(-2*pow(m_pi,4) + pow(mrho,4) - 3*pow(mrho,2)*s + 2*pow(s,2) + pow(m_pi,2)*(pow(mrho,2) + s)))*t1)/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (pow(eta1 - eta2,2)*(pow(eta1,2)*(3*pow(ma1,4) + 4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) -
                                  4*pow(ma1,2)*(2*pow(m_pi,2) + pow(mrho,2) - s) - 4*pow(mrho,2)*s + 2*pow(s,2)) +
                               pow(eta2,2)*(3*pow(ma1,4) + 4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s) +
                                  4*pow(ma1,2)*(-2*pow(m_pi,2) + pow(mrho,2) + s)) -
                               2*eta1*eta2*(3*pow(ma1,4) + 4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2) +
                                  pow(ma1,2)*(-8*pow(m_pi,2) + 4*s)))*t1)/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) -
                          (8*(pow(delta,2)*(8*pow(m_pi,4) + 3*pow(mrho,4) + 4*pow(m_pi,2)*(3*pow(mrho,2) - 2*s) - 6*pow(mrho,2)*s + 2*pow(s,2)) +
                               4*pow(mrho,4)*(3 + 12*C4*(2*pow(m_pi,2) - s) + 8*pow(C4,2)*pow(-2*pow(m_pi,2) + s,2)) -
                               4*delta*pow(mrho,2)*(16*C4*pow(m_pi,4) + 2*pow(m_pi,2)*(3 + 6*C4*pow(mrho,2) - 8*C4*s) + pow(mrho,2)*(3 - 6*C4*s) + s*(-3 + 4*C4*s)))*t1)/
                           (pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (8*(-2 + delta)*(pow(m_pi,4)*(-2 + 3*delta - 8*C4*pow(mrho,2)) + (pow(mrho,2) - s)*((-2 + 3*delta)*s + pow(mrho,2)*(-2 + delta - 8*C4*s)) +
                               4*pow(m_pi,2)*(2*C4*pow(mrho,4) + delta*s - pow(mrho,2)*(-1 + delta + 4*C4*s)))*t1)/
                           (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (4*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(eta2*(pow(m_pi,2) + s)*(pow(m_pi,4) - pow(m_pi,2)*(pow(mrho,2) - 2*s) + (pow(mrho,2) - s)*s) +
                               eta1*(-4*pow(m_pi,6) + pow(pow(mrho,2) - s,2)*s + pow(m_pi,4)*(3*pow(mrho,2) + s) -
                                  pow(m_pi,2)*(pow(mrho,4) - pow(mrho,2)*s + 2*pow(s,2))))*t1)/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,2) - s)*
                             (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (pow(eta1 - eta2,2)*(-2*eta1*eta2*(pow(m_pi,8) - pow(mrho,4)*pow(s,2) + pow(s,4) -
                                  pow(m_pi,4)*(pow(mrho,4) + 2*pow(mrho,2)*s - 4*pow(s,2)) + 2*pow(m_pi,2)*s*(pow(mrho,4) + pow(mrho,2)*s - 2*pow(s,2))) +
                               pow(eta2,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + pow(s,2)*pow(pow(mrho,2) + s,2) + pow(m_pi,4)*pow(pow(mrho,2) + 2*s,2) -
                                  2*pow(m_pi,2)*s*(pow(mrho,4) + 2*pow(mrho,2)*s + 2*pow(s,2))) +
                               pow(eta1,2)*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) - 4*pow(m_pi,2)*pow(pow(mrho,2) - s,2)*s + pow(pow(mrho,2) - s,2)*pow(s,2) +
                                  pow(m_pi,4)*(3*pow(mrho,4) - 6*pow(mrho,2)*s + 4*pow(s,2))))*t1)/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (2*pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(pow(eta1,2)*(pow(ma1,4)*s + pow(m_pi,4)*(-3*pow(mrho,2) + 2*s) +
                                  s*(2*pow(mrho,4) - 3*pow(mrho,2)*s + pow(s,2)) - 2*pow(m_pi,2)*(pow(mrho,4) - 4*pow(mrho,2)*s + 2*pow(s,2)) +
                                  pow(ma1,2)*(2*pow(m_pi,2)*(pow(mrho,2) - 2*s) + 3*s*(-pow(mrho,2) + s))) -
                               2*eta1*eta2*(pow(ma1,4)*s + s*(2*pow(m_pi,4) + 4*pow(m_pi,2)*(pow(mrho,2) - s) + s*(-2*pow(mrho,2) + s)) +
                                  pow(ma1,2)*(pow(m_pi,2)*(pow(mrho,2) - 4*s) + s*(-2*pow(mrho,2) + 3*s))) +
                               pow(eta2,2)*(-4*pow(m_pi,2)*s*(pow(ma1,2) + pow(mrho,2) + s) + pow(m_pi,4)*(pow(mrho,2) + 2*s) +
                                  s*(pow(ma1,4) + s*(pow(mrho,2) + s) + pow(ma1,2)*(pow(mrho,2) + 3*s))))*t1)/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (4*(eta1 - eta2)*(pow(ma1,2) - s)*(eta1*(4*pow(m_pi,4)*(6*C4*pow(mrho,4) + 2*delta*s + pow(mrho,2)*(1 - 2*delta - 8*C4*s)) -
                                  2*pow(m_pi,2)*(4*delta*pow(s,2) + pow(mrho,2)*s*(6 - 7*delta - 16*C4*s) + 2*pow(mrho,4)*(-2 + delta + 8*C4*s)) +
                                  (pow(mrho,2) - s)*s*(-2*delta*s + pow(mrho,2)*(-6 + 3*delta + 8*C4*s))) +
                               eta2*(-(delta*(2*pow(m_pi,4)*(pow(mrho,2) + 4*s) + pow(m_pi,2)*(2*pow(mrho,4) + pow(mrho,2)*s - 8*pow(s,2)) +
                                       s*(-2*pow(mrho,4) - pow(mrho,2)*s + 2*pow(s,2)))) +
                                  2*pow(mrho,2)*(4*C4*pow(m_pi,4)*(pow(mrho,2) + 4*s) + pow(m_pi,2)*(s*(5 - 16*C4*s) + pow(mrho,2)*(2 - 8*C4*s)) +
                                     s*(s*(-3 + 4*C4*s) + pow(mrho,2)*(-2 + 4*C4*s)))))*t1)/
                           (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                             (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (8*(eta1 - eta2)*(delta*(eta1*(4*pow(m_pi,6) + pow(m_pi,4)*(7*pow(mrho,2) - 8*s) + pow(ma1,4)*(pow(m_pi,2) - s) -
                                     pow(ma1,2)*(2*pow(m_pi,2) + pow(mrho,2) - 2*s)*(2*pow(m_pi,2) - s) + pow(m_pi,2)*s*(-8*pow(mrho,2) + 5*s) +
                                     s*(pow(mrho,4) + pow(mrho,2)*s - pow(s,2))) +
                                  eta2*(-4*pow(m_pi,6) - pow(m_pi,4)*(pow(mrho,2) - 8*s) + pow(ma1,4)*(-pow(m_pi,2) + s) +
                                     pow(m_pi,2)*(2*pow(mrho,4) - 5*pow(s,2)) + s*(-2*pow(mrho,4) + pow(mrho,2)*s + pow(s,2)) +
                                     pow(ma1,2)*(4*pow(m_pi,4) - 6*pow(m_pi,2)*s + s*(pow(mrho,2) + 2*s)))) -
                               2*pow(mrho,2)*(eta1*(8*C4*pow(m_pi,6) + pow(m_pi,4)*(3 + 8*C4*(pow(mrho,2) - 2*s)) + 2*C4*pow(ma1,4)*(pow(m_pi,2) - s) +
                                     2*pow(m_pi,2)*s*(-1 - 6*C4*pow(mrho,2) + 5*C4*s) -
                                     pow(ma1,2)*(8*C4*pow(m_pi,4) + pow(m_pi,2)*(1 + 2*C4*(pow(mrho,2) - 6*s)) + 2*C4*s*(-pow(mrho,2) + 2*s)) +
                                     s*(-(s*(1 + 2*C4*s)) + pow(mrho,2)*(1 + 4*C4*s))) +
                                  eta2*(2*C4*pow(ma1,4)*(-pow(m_pi,2) + s) - (pow(m_pi,2) - s)*
                                      (8*C4*pow(m_pi,4) - 2*pow(mrho,2) + s + 2*C4*pow(s,2) + pow(m_pi,2)*(3 - 4*C4*(pow(mrho,2) + 2*s))) +
                                     pow(ma1,2)*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(mrho,2) + 2*s) + pow(m_pi,2)*(1 - 2*C4*(pow(mrho,2) + 6*s))))))*t1)/
                           (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (8*(-2 + delta)*(delta - 4*C4*pow(mrho,2))*pow(t1,2))/
                           (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (16*(-2 + delta)*(delta - 4*C4*pow(mrho,2))*s*pow(t1,2))/
                           (pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (pow(eta1 - eta2,2)*(pow(eta1,2)*(pow(mrho,2) - s) + 2*eta1*eta2*s - pow(eta2,2)*s)*
                             (pow(m_pi,4) - pow(m_pi,2)*(pow(mrho,2) - 2*s) + (pow(mrho,2) - s)*s)*pow(t1,2))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (2*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(-(eta1*(pow(m_pi,2) + 2*pow(mrho,2) - 3*s)) - eta2*(pow(m_pi,2) + s))*pow(t1,2))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (2*(-2 + delta)*(eta1 - eta2)*(pow(ma1,2) - s)*(pow(m_pi,2) + s)*(-2*eta2*s + eta1*(pow(m_pi,2) - pow(mrho,2) + s))*pow(t1,2))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,2) - s)*
                             (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (pow(eta1 - eta2,3)*(-(eta1*(pow(ma1,2) - 2*pow(m_pi,2) - pow(mrho,2) + s)) + eta2*(pow(ma1,2) - 2*pow(m_pi,2) + pow(mrho,2) + s))*
                             pow(t1,2))/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) +
                          (8*(delta - 4*C4*pow(mrho,2))*(delta*(4*pow(m_pi,2) + 3*pow(mrho,2) - 2*s) - 2*pow(mrho,2)*(3 + 8*C4*pow(m_pi,2) - 4*C4*s))*pow(t1,2))/
                           (pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(pow(eta2,2)*s*(pow(ma1,2) - 4*pow(m_pi,2) + pow(mrho,2) + 3*s) +
                               pow(eta1,2)*(2*pow(m_pi,2)*(pow(mrho,2) - 2*s) + s*(pow(ma1,2) - 3*pow(mrho,2) + 3*s)) -
                               2*eta1*eta2*(pow(m_pi,2)*(pow(mrho,2) - 4*s) + s*(pow(ma1,2) - 2*pow(mrho,2) + 3*s)))*pow(t1,2))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (2*(eta1 - eta2)*(pow(ma1,2) - s)*(eta1*(4*delta*pow(s,2) - 2*pow(mrho,2)*s*(-2 + 3*delta + 8*C4*s) + pow(mrho,4)*(-2 + delta + 16*C4*s) -
                                  2*pow(m_pi,2)*(8*C4*pow(mrho,4) + 4*delta*s + pow(mrho,2)*(2 - 3*delta - 16*C4*s))) +
                               eta2*(pow(m_pi,2)*(8*delta*s + pow(mrho,2)*(-2 + delta - 32*C4*s)) + s*(-4*delta*s + pow(mrho,2)*(-2 + delta + 16*C4*s))))*pow(t1,2))/
                           (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                             (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (4*(eta1 - eta2)*(delta*(eta1*(pow(ma1,2)*(pow(m_pi,2) - s) - (2*pow(m_pi,2) + pow(mrho,2) - 2*s)*(2*pow(m_pi,2) - s)) +
                                  eta2*(4*pow(m_pi,4) - 6*pow(m_pi,2)*s + pow(ma1,2)*(-pow(m_pi,2) + s) + s*(pow(mrho,2) + 2*s))) +
                               2*pow(mrho,2)*(eta1*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(ma1,2) - pow(mrho,2) + 2*s) +
                                     pow(m_pi,2)*(1 - 2*C4*(pow(ma1,2) - pow(mrho,2) + 6*s))) -
                                  eta2*(8*C4*pow(m_pi,4) + 2*C4*s*(pow(ma1,2) + pow(mrho,2) + 2*s) - pow(m_pi,2)*(-1 + 2*C4*(pow(ma1,2) + pow(mrho,2) + 6*s)))))*
                             pow(t1,2))/(pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (pow(eta1 - eta2,4)*pow(t1,3))/(3.*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (8*pow(eta1 - eta2,2)*(delta - 4*C4*pow(mrho,2))*pow(t1,3))/
                           (3.*pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (16*pow(delta - 4*C4*pow(mrho,2),2)*pow(t1,3))/
                           (3.*pow(mrho,4)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (4*(-2 + delta)*eta1*(eta1 - eta2)*(pow(ma1,2) - s)*pow(t1,3))/
                           (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (2*pow(eta1 - eta2,4)*(pow(ma1,2) - s)*s*pow(t1,3))/
                           (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (2*pow(eta1 - eta2,2)*s*(-2*eta1*eta2*s + pow(eta2,2)*s + pow(eta1,2)*(-pow(mrho,2) + s))*pow(t1,3))/
                           (3.*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (4*(eta1 - eta2)*(pow(ma1,2) - s)*(2*eta2*(delta - 4*C4*pow(mrho,2))*s + eta1*(-2*delta*s + pow(mrho,2)*(-2 + delta + 8*C4*s)))*pow(t1,3))/
                           (3.*pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*
                             (pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (2*pow(eta1 - eta2,2)*(pow(eta1,2)*(2*pow(ma1,6) - 3*pow(ma1,4)*(2*pow(m_pi,2) + pow(mrho,2) - s) + pow(mrho,2)*(pow(mrho,2) - s)*s -
                                  pow(m_pi,4)*(3*pow(mrho,2) + s) + pow(m_pi,2)*(-2*pow(mrho,4) + 3*pow(mrho,2)*s) +
                                  pow(ma1,2)*(4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) - 4*pow(mrho,2)*s + 2*pow(s,2))) -
                               2*eta1*eta2*(2*pow(ma1,6) - pow(m_pi,2)*(pow(m_pi,2) - pow(mrho,2))*(pow(mrho,2) + s) + pow(ma1,4)*(-6*pow(m_pi,2) + 3*s) +
                                  pow(ma1,2)*(4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2))) +
                               pow(eta2,2)*(2*pow(ma1,6) + 3*pow(ma1,4)*(-2*pow(m_pi,2) + pow(mrho,2) + s) +
                                  pow(m_pi,2)*(-pow(mrho,4) + pow(m_pi,2)*(2*pow(mrho,2) - s) + pow(mrho,2)*s) +
                                  pow(ma1,2)*(4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s))))*log(fabs(-pow(ma1,2) + t2)))
                            /(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) -
                          (2*pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(-2*eta1*eta2*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + 2*pow(ma1,2)*pow(m_pi,4)*s +
                                  pow(m_pi,2)*(pow(ma1,4)*(pow(mrho,2) - 4*s) + 4*pow(ma1,2)*(pow(mrho,2) - s)*s + pow(mrho,2)*pow(s,2)) +
                                  pow(ma1,2)*s*(pow(ma1,4) + s*(-2*pow(mrho,2) + s) + pow(ma1,2)*(-2*pow(mrho,2) + 3*s))) +
                               pow(eta2,2)*(pow(m_pi,8) - 4*pow(ma1,2)*pow(m_pi,2)*s*(pow(ma1,2) + pow(mrho,2) + s) +
                                  pow(m_pi,4)*(pow(mrho,2)*s + pow(ma1,2)*(pow(mrho,2) + 2*s)) +
                                  pow(ma1,2)*s*(pow(ma1,4) + s*(pow(mrho,2) + s) + pow(ma1,2)*(pow(mrho,2) + 3*s))) +
                               pow(eta1,2)*(pow(m_pi,8) + pow(ma1,2)*s*(pow(ma1,4) + 2*pow(mrho,4) - 3*pow(ma1,2)*(pow(mrho,2) - s) - 3*pow(mrho,2)*s + pow(s,2)) +
                                  pow(m_pi,4)*(2*pow(mrho,4) - 3*pow(mrho,2)*s + pow(ma1,2)*(-3*pow(mrho,2) + 2*s)) +
                                  2*pow(m_pi,2)*(pow(ma1,4)*(pow(mrho,2) - 2*s) + pow(mrho,2)*s*(-pow(mrho,2) + s) -
                                     pow(ma1,2)*(pow(mrho,4) - 4*pow(mrho,2)*s + 2*pow(s,2)))))*log(fabs(-pow(ma1,2) + t2)))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (8*(eta1 - eta2)*(delta*(eta2*(pow(m_pi,6)*pow(mrho,2)*(2*pow(m_pi,2) - s) + pow(ma1,8)*(-pow(m_pi,2) + s) +
                                     pow(ma1,6)*(5*pow(m_pi,4) - 7*pow(m_pi,2)*s + s*(pow(mrho,2) + 2*s)) +
                                     pow(ma1,4)*(-8*pow(m_pi,6) - pow(m_pi,4)*(pow(mrho,2) - 14*s) + pow(m_pi,2)*(2*pow(mrho,4) - pow(mrho,2)*s - 7*pow(s,2)) +
                                        s*(-2*pow(mrho,4) + pow(mrho,2)*s + pow(s,2))) +
                                     pow(ma1,2)*pow(m_pi,2)*(4*pow(m_pi,6) + pow(m_pi,4)*(pow(mrho,2) - 8*s) + s*(2*pow(mrho,4) + pow(mrho,2)*s - pow(s,2)) +
                                        pow(m_pi,2)*(-2*pow(mrho,4) - 3*pow(mrho,2)*s + 5*pow(s,2)))) +
                                  eta1*(pow(ma1,8)*(pow(m_pi,2) - s) + pow(ma1,6)*(-5*pow(m_pi,4) + (pow(mrho,2) - 2*s)*s + pow(m_pi,2)*(-2*pow(mrho,2) + 7*s)) +
                                     pow(m_pi,2)*pow(mrho,2)*(2*pow(m_pi,6) + pow(m_pi,4)*(4*pow(mrho,2) - 5*s) + pow(mrho,4)*s -
                                        pow(m_pi,2)*(pow(mrho,4) + 3*pow(mrho,2)*s - 2*pow(s,2))) +
                                     pow(ma1,4)*(8*pow(m_pi,6) + pow(m_pi,4)*(9*pow(mrho,2) - 14*s) + pow(m_pi,2)*s*(-9*pow(mrho,2) + 7*s) +
                                        s*(pow(mrho,4) + pow(mrho,2)*s - pow(s,2))) +
                                     pow(ma1,2)*(-4*pow(m_pi,8) + pow(mrho,4)*s*(-pow(mrho,2) + s) + pow(m_pi,6)*(-11*pow(mrho,2) + 8*s) +
                                        pow(m_pi,4)*(-3*pow(mrho,4) + 17*pow(mrho,2)*s - 5*pow(s,2)) + pow(m_pi,2)*(pow(mrho,6) - 5*pow(mrho,2)*pow(s,2) + pow(s,3))
                                        ))) - 2*pow(mrho,2)*(eta2*(pow(m_pi,8)*(1 + 2*C4*pow(mrho,2)) - 2*C4*pow(m_pi,6)*pow(mrho,2)*s +
                                     2*C4*pow(ma1,8)*(-pow(m_pi,2) + s) + pow(ma1,4)*
                                      (-16*C4*pow(m_pi,6) + pow(m_pi,4)*(-4 + 6*C4*pow(mrho,2) + 28*C4*s) +
                                        2*pow(m_pi,2)*(pow(mrho,2) + s - 3*C4*pow(mrho,2)*s - 7*C4*pow(s,2)) + s*(-2*pow(mrho,2) + s + 2*C4*pow(s,2))) +
                                     pow(ma1,6)*(10*C4*pow(m_pi,4) + 2*C4*s*(pow(mrho,2) + 2*s) + pow(m_pi,2)*(1 - 2*C4*(pow(mrho,2) + 7*s))) +
                                     pow(ma1,2)*pow(m_pi,2)*(8*C4*pow(m_pi,6) - 2*pow(m_pi,4)*(-2 + 3*C4*pow(mrho,2) + 8*C4*s) +
                                        s*(2*pow(mrho,2) + s - 2*C4*pow(s,2)) + 2*pow(m_pi,2)*(pow(mrho,2)*(-1 + 3*C4*s) + s*(-3 + 5*C4*s)))) +
                                  eta1*(pow(m_pi,8)*(-1 + 6*C4*pow(mrho,2)) + 2*C4*pow(ma1,8)*(pow(m_pi,2) - s) + pow(m_pi,2)*pow(mrho,4)*s +
                                     2*pow(m_pi,6)*pow(mrho,2)*(2 - 5*C4*s) - pow(ma1,6)*
                                      (10*C4*pow(m_pi,4) + pow(m_pi,2)*(1 + 2*C4*(pow(mrho,2) - 7*s)) + 2*C4*s*(-pow(mrho,2) + 2*s)) -
                                     pow(m_pi,4)*pow(mrho,2)*(pow(mrho,2) + s*(3 - 4*C4*s)) +
                                     pow(ma1,4)*(16*C4*pow(m_pi,6) + 2*pow(m_pi,4)*(2 + 5*C4*pow(mrho,2) - 14*C4*s) + 2*pow(m_pi,2)*s*(-1 - 7*C4*pow(mrho,2) + 7*C4*s) +
                                        s*(-(s*(1 + 2*C4*s)) + pow(mrho,2)*(1 + 4*C4*s))) -
                                     pow(ma1,2)*(8*C4*pow(m_pi,8) + pow(mrho,2)*(pow(mrho,2) - s)*s + 2*pow(m_pi,6)*(2 + 7*C4*pow(mrho,2) - 8*C4*s) +
                                        pow(m_pi,2)*(-pow(mrho,4) + pow(s,2) + 8*C4*pow(mrho,2)*pow(s,2) - 2*C4*pow(s,3)) +
                                        pow(m_pi,4)*(pow(mrho,2)*(3 - 22*C4*s) + 2*s*(-3 + 5*C4*s))))))*log(fabs(-pow(ma1,2) + t2)))/
                           ((pow(ma1,2) - pow(m_pi,2))*pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (16*pow(-2 + delta,2)*pow(m_pi,2)*log(fabs(-pow(m_pi,2) + t2)))/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) -
                          (8*pow(-2 + delta,2)*(3*pow(m_pi,4) - 4*pow(m_pi,2)*(pow(mrho,2) - s) + pow(pow(mrho,2) - s,2))*log(fabs(-pow(m_pi,2) + t2)))/
                           ((pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (8*(-2 + delta)*(eta1 - eta2)*pow(m_pi,2)*(2*eta1*pow(m_pi,2) - 2*eta2*pow(m_pi,2) - eta1*pow(mrho,2))*(pow(m_pi,2) - s)*
                             log(fabs(-pow(m_pi,2) + t2)))/((pow(ma1,2) - pow(m_pi,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (8*(-2 + delta)*(eta1 - eta2)*pow(m_pi,2)*(pow(ma1,2) - s)*(pow(m_pi,2) - s)*
                             (-(eta2*(pow(m_pi,2) + s)) + eta1*(pow(m_pi,2) - pow(mrho,2) + s))*log(fabs(-pow(m_pi,2) + t2)))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (8*(-2 + delta)*(-(delta*(4*pow(m_pi,2) - pow(mrho,2))*(pow(m_pi,2) + pow(mrho,2) - s)) +
                               2*pow(mrho,2)*(8*C4*pow(m_pi,4) - pow(mrho,2) + s + pow(m_pi,2)*(3 - 8*C4*s)))*log(fabs(-pow(m_pi,2) + t2)))/
                           (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (2*pow(eta1 - eta2,2)*(pow(eta1,2)*(2*pow(ma1,6) - 3*pow(ma1,4)*(2*pow(m_pi,2) + pow(mrho,2) - s) + pow(mrho,2)*(pow(mrho,2) - s)*s -
                                  pow(m_pi,4)*(3*pow(mrho,2) + s) + pow(m_pi,2)*(-2*pow(mrho,4) + 3*pow(mrho,2)*s) +
                                  pow(ma1,2)*(4*pow(m_pi,4) + pow(mrho,4) + pow(m_pi,2)*(8*pow(mrho,2) - 4*s) - 4*pow(mrho,2)*s + 2*pow(s,2))) -
                               2*eta1*eta2*(2*pow(ma1,6) - pow(m_pi,2)*(pow(m_pi,2) - pow(mrho,2))*(pow(mrho,2) + s) + pow(ma1,4)*(-6*pow(m_pi,2) + 3*s) +
                                  pow(ma1,2)*(4*pow(m_pi,4) - pow(mrho,4) + 2*pow(m_pi,2)*(pow(mrho,2) - 2*s) - 2*pow(mrho,2)*s + 2*pow(s,2))) +
                               pow(eta2,2)*(2*pow(ma1,6) + 3*pow(ma1,4)*(-2*pow(m_pi,2) + pow(mrho,2) + s) +
                                  pow(m_pi,2)*(-pow(mrho,4) + pow(m_pi,2)*(2*pow(mrho,2) - s) + pow(mrho,2)*s) +
                                  pow(ma1,2)*(4*pow(m_pi,4) + pow(mrho,4) - 2*pow(mrho,2)*s + 2*pow(s,2) - 4*pow(m_pi,2)*(pow(mrho,2) + s))))*log(fabs(-pow(ma1,2) + t1)))
                            /(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) +
                          (2*pow(eta1 - eta2,2)*(pow(ma1,2) - s)*(-2*eta1*eta2*(pow(m_pi,8) - 2*pow(m_pi,6)*pow(mrho,2) + 2*pow(ma1,2)*pow(m_pi,4)*s +
                                  pow(m_pi,2)*(pow(ma1,4)*(pow(mrho,2) - 4*s) + 4*pow(ma1,2)*(pow(mrho,2) - s)*s + pow(mrho,2)*pow(s,2)) +
                                  pow(ma1,2)*s*(pow(ma1,4) + s*(-2*pow(mrho,2) + s) + pow(ma1,2)*(-2*pow(mrho,2) + 3*s))) +
                               pow(eta2,2)*(pow(m_pi,8) - 4*pow(ma1,2)*pow(m_pi,2)*s*(pow(ma1,2) + pow(mrho,2) + s) +
                                  pow(m_pi,4)*(pow(mrho,2)*s + pow(ma1,2)*(pow(mrho,2) + 2*s)) +
                                  pow(ma1,2)*s*(pow(ma1,4) + s*(pow(mrho,2) + s) + pow(ma1,2)*(pow(mrho,2) + 3*s))) +
                               pow(eta1,2)*(pow(m_pi,8) + pow(ma1,2)*s*(pow(ma1,4) + 2*pow(mrho,4) - 3*pow(ma1,2)*(pow(mrho,2) - s) - 3*pow(mrho,2)*s + pow(s,2)) +
                                  pow(m_pi,4)*(2*pow(mrho,4) - 3*pow(mrho,2)*s + pow(ma1,2)*(-3*pow(mrho,2) + 2*s)) +
                                  2*pow(m_pi,2)*(pow(ma1,4)*(pow(mrho,2) - 2*s) + pow(mrho,2)*s*(-pow(mrho,2) + s) -
                                     pow(ma1,2)*(pow(mrho,4) - 4*pow(mrho,2)*s + 2*pow(s,2)))))*log(fabs(-pow(ma1,2) + t1)))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (8*(eta1 - eta2)*(delta*(eta2*(pow(m_pi,6)*pow(mrho,2)*(2*pow(m_pi,2) - s) + pow(ma1,8)*(-pow(m_pi,2) + s) +
                                     pow(ma1,6)*(5*pow(m_pi,4) - 7*pow(m_pi,2)*s + s*(pow(mrho,2) + 2*s)) +
                                     pow(ma1,4)*(-8*pow(m_pi,6) - pow(m_pi,4)*(pow(mrho,2) - 14*s) + pow(m_pi,2)*(2*pow(mrho,4) - pow(mrho,2)*s - 7*pow(s,2)) +
                                        s*(-2*pow(mrho,4) + pow(mrho,2)*s + pow(s,2))) +
                                     pow(ma1,2)*pow(m_pi,2)*(4*pow(m_pi,6) + pow(m_pi,4)*(pow(mrho,2) - 8*s) + s*(2*pow(mrho,4) + pow(mrho,2)*s - pow(s,2)) +
                                        pow(m_pi,2)*(-2*pow(mrho,4) - 3*pow(mrho,2)*s + 5*pow(s,2)))) +
                                  eta1*(pow(ma1,8)*(pow(m_pi,2) - s) + pow(ma1,6)*(-5*pow(m_pi,4) + (pow(mrho,2) - 2*s)*s + pow(m_pi,2)*(-2*pow(mrho,2) + 7*s)) +
                                     pow(m_pi,2)*pow(mrho,2)*(2*pow(m_pi,6) + pow(m_pi,4)*(4*pow(mrho,2) - 5*s) + pow(mrho,4)*s -
                                        pow(m_pi,2)*(pow(mrho,4) + 3*pow(mrho,2)*s - 2*pow(s,2))) +
                                     pow(ma1,4)*(8*pow(m_pi,6) + pow(m_pi,4)*(9*pow(mrho,2) - 14*s) + pow(m_pi,2)*s*(-9*pow(mrho,2) + 7*s) +
                                        s*(pow(mrho,4) + pow(mrho,2)*s - pow(s,2))) +
                                     pow(ma1,2)*(-4*pow(m_pi,8) + pow(mrho,4)*s*(-pow(mrho,2) + s) + pow(m_pi,6)*(-11*pow(mrho,2) + 8*s) +
                                        pow(m_pi,4)*(-3*pow(mrho,4) + 17*pow(mrho,2)*s - 5*pow(s,2)) + pow(m_pi,2)*(pow(mrho,6) - 5*pow(mrho,2)*pow(s,2) + pow(s,3))
                                        ))) - 2*pow(mrho,2)*(eta2*(pow(m_pi,8)*(1 + 2*C4*pow(mrho,2)) - 2*C4*pow(m_pi,6)*pow(mrho,2)*s +
                                     2*C4*pow(ma1,8)*(-pow(m_pi,2) + s) + pow(ma1,4)*
                                      (-16*C4*pow(m_pi,6) + pow(m_pi,4)*(-4 + 6*C4*pow(mrho,2) + 28*C4*s) +
                                        2*pow(m_pi,2)*(pow(mrho,2) + s - 3*C4*pow(mrho,2)*s - 7*C4*pow(s,2)) + s*(-2*pow(mrho,2) + s + 2*C4*pow(s,2))) +
                                     pow(ma1,6)*(10*C4*pow(m_pi,4) + 2*C4*s*(pow(mrho,2) + 2*s) + pow(m_pi,2)*(1 - 2*C4*(pow(mrho,2) + 7*s))) +
                                     pow(ma1,2)*pow(m_pi,2)*(8*C4*pow(m_pi,6) - 2*pow(m_pi,4)*(-2 + 3*C4*pow(mrho,2) + 8*C4*s) +
                                        s*(2*pow(mrho,2) + s - 2*C4*pow(s,2)) + 2*pow(m_pi,2)*(pow(mrho,2)*(-1 + 3*C4*s) + s*(-3 + 5*C4*s)))) +
                                  eta1*(pow(m_pi,8)*(-1 + 6*C4*pow(mrho,2)) + 2*C4*pow(ma1,8)*(pow(m_pi,2) - s) + pow(m_pi,2)*pow(mrho,4)*s +
                                     2*pow(m_pi,6)*pow(mrho,2)*(2 - 5*C4*s) - pow(ma1,6)*
                                      (10*C4*pow(m_pi,4) + pow(m_pi,2)*(1 + 2*C4*(pow(mrho,2) - 7*s)) + 2*C4*s*(-pow(mrho,2) + 2*s)) -
                                     pow(m_pi,4)*pow(mrho,2)*(pow(mrho,2) + s*(3 - 4*C4*s)) +
                                     pow(ma1,4)*(16*C4*pow(m_pi,6) + 2*pow(m_pi,4)*(2 + 5*C4*pow(mrho,2) - 14*C4*s) + 2*pow(m_pi,2)*s*(-1 - 7*C4*pow(mrho,2) + 7*C4*s) +
                                        s*(-(s*(1 + 2*C4*s)) + pow(mrho,2)*(1 + 4*C4*s))) -
                                     pow(ma1,2)*(8*C4*pow(m_pi,8) + pow(mrho,2)*(pow(mrho,2) - s)*s + 2*pow(m_pi,6)*(2 + 7*C4*pow(mrho,2) - 8*C4*s) +
                                        pow(m_pi,2)*(-pow(mrho,4) + pow(s,2) + 8*C4*pow(mrho,2)*pow(s,2) - 2*C4*pow(s,3)) +
                                        pow(m_pi,4)*(pow(mrho,2)*(3 - 22*C4*s) + 2*s*(-3 + 5*C4*s))))))*log(fabs(-pow(ma1,2) + t1)))/
                           ((pow(ma1,2) - pow(m_pi,2))*pow(mrho,2)*(pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (16*pow(-2 + delta,2)*pow(m_pi,2)*log(fabs(-pow(m_pi,2) + t1)))/(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)) +
                          (8*pow(-2 + delta,2)*(3*pow(m_pi,4) - 4*pow(m_pi,2)*(pow(mrho,2) - s) + pow(pow(mrho,2) - s,2))*log(fabs(-pow(m_pi,2) + t1)))/
                           ((pow(m_pi,2) - s)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (8*(-2 + delta)*(eta1 - eta2)*pow(m_pi,2)*(2*eta1*pow(m_pi,2) - 2*eta2*pow(m_pi,2) - eta1*pow(mrho,2))*(pow(m_pi,2) - s)*
                             log(fabs(-pow(m_pi,2) + t1)))/((pow(ma1,2) - pow(m_pi,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) -
                          (8*(-2 + delta)*(eta1 - eta2)*pow(m_pi,2)*(pow(ma1,2) - s)*(pow(m_pi,2) - s)*
                             (-(eta2*(pow(m_pi,2) + s)) + eta1*(pow(m_pi,2) - pow(mrho,2) + s))*log(fabs(-pow(m_pi,2) + t1)))/
                           ((pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) - s,2))*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s))) +
                          (8*(-2 + delta)*(delta*(4*pow(m_pi,2) - pow(mrho,2))*(pow(m_pi,2) + pow(mrho,2) - s) -
                               2*pow(mrho,2)*(8*C4*pow(m_pi,4) - pow(mrho,2) + s + pow(m_pi,2)*(3 - 8*C4*s)))*log(fabs(-pow(m_pi,2) + t1)))/
                           (pow(mrho,2)*(pow(m_pi,4) + pow(pow(mrho,2) - s,2) - 2*pow(m_pi,2)*(pow(mrho,2) + s)))))/(512.*Pi);
  return xs;
}


double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_pi_rho0(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const double s_sqrt = sqrt(s);
  const double mpion = m_pion_, mrho = m_rho;
  
  auto mandelstam_t = get_t_range(s_sqrt, m_pion_, m_pion_, mrho, 0.);

  double t1 = mandelstam_t[1];
  double t2 = mandelstam_t[0];

  double &tmin = t1, &tmax = t2;


  const double xs =
          to_mb*(pow(Const,2)*pow(ghat,4)*((0.5*pow(-2. + delta,2)*pow(mpion,2)*(0. - 1.*tmax))/pow(mrho,2) +
        0.25*(-2. + delta)*(eta1 - 1.*eta2)*(-0.5*eta2*pow(ma1,2) + 1.*eta1*pow(mpion,2) - 1.*eta2*pow(mpion,2) + 0.5*eta1*s -
            0.5*eta2*s)*(0. - 1.*tmax) - (2.*(pow(mpion,2)*
              (1.5 + 0.125*pow(delta,2) - 2.*C4*pow(mrho,2) + delta*(-1. + 1.*C4*pow(mrho,2))) +
              (-0.5 - 0.375*pow(delta,2) - 2.*C4*pow(mrho,2) + delta*(1. + 1.*C4*pow(mrho,2)))*s)*(0. - 1.*tmax))/pow(mrho,2) -
        0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(1.*pow(mpion,2) - 0.5*s) +
            eta2*(-0.5*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 1.*s))*(0. - 1.*tmax) -
        (0.25*(-2. + 1.*delta)*(-8.*C4*pow(mrho,4) + pow(mpion,2)*(2. + 1.*delta - 8.*C4*pow(mrho,2)) + (-2. - 3.*delta)*s +
              pow(mrho,2)*(2. + 1.*delta + 16.*C4*s))*(0. - 1.*tmax))/pow(mrho,2) -
        0.09375*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-2.*pow(ma1,4) - 4.*pow(mpion,4) + 0.6666666666666666*pow(mrho,4) +
              pow(ma1,2)*(5.333333333333333*pow(mpion,2) - 2.6666666666666665*s) + 2.6666666666666665*pow(mpion,2)*s +
              1.3333333333333333*pow(mrho,2)*s - 1.3333333333333333*pow(s,2)) +
            pow(eta1,2)*(1.*pow(ma1,4) + 2.*pow(mpion,4) + 0.3333333333333333*pow(mrho,4) +
              pow(mpion,2)*(2.*pow(mrho,2) - 1.3333333333333333*s) - 1.3333333333333333*pow(mrho,2)*s + 0.6666666666666666*pow(s,2) +
              pow(ma1,2)*(-2.6666666666666665*pow(mpion,2) - 1.3333333333333333*pow(mrho,2) + 1.3333333333333333*s)) +
            pow(eta2,2)*(1.*pow(ma1,4) + 2.*pow(mpion,4) + 0.3333333333333333*pow(mrho,4) +
              pow(mpion,2)*(-2.*pow(mrho,2) - 1.3333333333333333*s) - 0.6666666666666666*pow(mrho,2)*s + 0.6666666666666666*pow(s,2) +
              pow(ma1,2)*(-2.6666666666666665*pow(mpion,2) + 1.3333333333333333*pow(mrho,2) + 1.3333333333333333*s)))*(0. - 1.*tmax) -
        0.0625*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(pow(Gammaa1,2)*pow(ma1,2) - 1.*pow(ma1,4) - 2.*pow(mpion,4) -
              2.*pow(mpion,2)*pow(mrho,2) + 2.*pow(mrho,4) + pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) -
              3.*pow(mrho,2)*s + pow(s,2)) + pow(eta1,2)*
            (pow(Gammaa1,2)*pow(ma1,2) - 1.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) +
              pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) + pow(mrho,2)*s + pow(s,2)) +
            eta1*eta2*(-2.*pow(Gammaa1,2)*pow(ma1,2) + 2.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) +
              2.*pow(mrho,4) - 2.*pow(s,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)))*(0. - 1.*tmax) +
        (1.*(pow(eta1,2)*(0.5*pow(mrho,4) - 1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) - 0.5*pow(mrho,2)*s + 0.25*delta*pow(mrho,2)*s -
                0.25*delta*pow(s,2) + 1.*C4*pow(mrho,2)*pow(s,2)) +
              pow(eta2,2)*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) - 0.5*pow(mrho,2)*s + 0.75*delta*pow(mrho,2)*s -
                2.*C4*pow(mrho,4)*s - 0.25*delta*pow(s,2) + 1.*C4*pow(mrho,2)*pow(s,2)) +
              eta1*eta2*(pow(ma1,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) + pow(mpion,2)*(-2.*pow(mrho,2) + 4.*C4*pow(mrho,4)) +
                s*(2.*C4*pow(mrho,4) + 0.5*delta*s + pow(mrho,2)*(1. - 1.*delta - 2.*C4*s))))*(0. - 1.*tmax))/pow(mrho,2) +
        (0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
                pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) + pow(ma1,2)*pow(mpion,2)*
                  (8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
                pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
              pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
                pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                    2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s +
                    pow(mpion,4)*(6.*pow(mrho,2) + 2.*s))) +
              pow(eta1,2)*(1.*pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                    2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
                    pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
                pow(mpion,2)*(1.*pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                    pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/(0. + 1.*pow(ma1,2) - 1.*tmax) +
        (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(0. + 1.*pow(mpion,2) - 1.*tmax) +
        0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(1.*pow(Gammaa1,2)*pow(ma1,2) - 3.*pow(ma1,4) - 2.*pow(mpion,4) +
              2.*pow(mpion,2)*pow(mrho,2) - 1.*pow(mrho,4) + pow(ma1,2)*(4.*pow(mpion,2) - 4.*pow(mrho,2) - 4.*s) +
              2.*pow(mrho,2)*s - 2.*pow(s,2)) + pow(eta1,2)*
            (1.*pow(Gammaa1,2)*pow(ma1,2) - 3.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) - 1.*pow(mrho,4) +
              pow(ma1,2)*(4.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) + 4.*pow(mrho,2)*s - 2.*pow(s,2)) +
            eta1*eta2*(-2.*pow(Gammaa1,2)*pow(ma1,2) + 6.*pow(ma1,4) + 4.*pow(mpion,4) - 2.*pow(mrho,4) - 4.*pow(mrho,2)*s +
              4.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) + 8.*s)))*(0. + 1.*pow(mrho,2) - 1.*s - 1.*tmax) -
        0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) - 1.*pow(mrho,2) - 1.*s) + eta1*(pow(ma1,2) - 1.*pow(mrho,2) + s))*
          pow(0. + 1.*pow(mrho,2) - 1.*s - 1.*tmax,2) -
        0.010416666666666666*pow(eta1 - 1.*eta2,4)*pow(0. + 1.*pow(mrho,2) - 1.*s - 1.*tmax,3) +
        (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/
          (0. + 1.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s - 1.*tmax) -
        (1.*(1.*eta1 - 1.*eta2)*(eta1*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6) + pow(ma1,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                pow(mpion,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) + 0.25*delta*pow(s,2) - 1.*C4*pow(mrho,2)*pow(s,2)) +
              eta2*(-0.75*pow(mrho,4) + 0.125*delta*pow(mrho,4) + 1.*C4*pow(mrho,6) +
                pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) + pow(ma1,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) +
                0.25*pow(mrho,2)*s + 0.375*delta*pow(mrho,2)*s - 2.*C4*pow(mrho,4)*s - 0.25*delta*pow(s,2) +
                1.*C4*pow(mrho,2)*pow(s,2)))*(0. - 1.*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s - 1.*tmax))/pow(mrho,2)\
          + 0.5*pow(1.*eta1 - 1.*eta2,2)*(-0.5 + 1.*C4*pow(mrho,2))*
          pow(0. - 1.*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s - 1.*tmax,2) +
        (0.25*(32.*pow(C4,2)*pow(mrho,8) + 2.*pow(delta,2)*pow(s,2) + 8.*C4*pow(mrho,6)*(-6. + delta - 8.*C4*s) +
              2.*delta*pow(mrho,2)*s*(-6. + delta - 8.*C4*s) +
              pow(mrho,4)*(12. - 1.*pow(delta,2) + 8.*C4*(6. + delta)*s + 32.*pow(C4,2)*pow(s,2)))*(0. + tmax))/pow(mrho,4) -
        (0.125*(-2. + 1.*delta)*(2. + 1.*delta - 8.*C4*pow(mrho,2))*(0. + 1.*pow(tmax,2)))/pow(mrho,2) -
        (1.*(0.5 - 0.125*pow(delta,2) - 2.*C4*pow(mrho,2) + 1.*C4*delta*pow(mrho,2))*(0. + 1.*pow(tmax,2)))/pow(mrho,2) -
        0.5*(eta1*eta2*(1. - 2.*C4*pow(mrho,2)) + pow(eta1,2)*(-0.5 + 1.*C4*pow(mrho,2)) + pow(eta2,2)*(-0.5 + 1.*C4*pow(mrho,2)))*
          (0. + 1.*pow(tmax,2)) + 0.0625*pow(1.*eta1 - 1.*eta2,4)*(1.*pow(mpion,2) + 0.5*pow(mrho,2) - 0.5*s)*(0. + 1.*pow(tmax,2)) +
        0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
            eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*(0. + 1.*pow(tmax,2)) +
        0.010416666666666666*pow(eta1 - 1.*eta2,4)*(0. + 1.*pow(tmax,3)) -
        0.020833333333333332*pow(1.*eta1 - 1.*eta2,4)*(0. + 1.*pow(tmax,3)) -
        (0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
                pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) + pow(ma1,2)*pow(mpion,2)*
                  (8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
                pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
              pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
                pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                    2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s +
                    pow(mpion,4)*(6.*pow(mrho,2) + 2.*s))) +
              pow(eta1,2)*(1.*pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                    2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
                    pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
                pow(mpion,2)*(1.*pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                    pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/(1.*pow(ma1,2) - 1.*tmin) -
        (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*tmin) -
        0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(1.*pow(Gammaa1,2)*pow(ma1,2) - 3.*pow(ma1,4) - 2.*pow(mpion,4) +
              2.*pow(mpion,2)*pow(mrho,2) - 1.*pow(mrho,4) + pow(ma1,2)*(4.*pow(mpion,2) - 4.*pow(mrho,2) - 4.*s) +
              2.*pow(mrho,2)*s - 2.*pow(s,2)) + pow(eta1,2)*
            (1.*pow(Gammaa1,2)*pow(ma1,2) - 3.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) - 1.*pow(mrho,4) +
              pow(ma1,2)*(4.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) + 4.*pow(mrho,2)*s - 2.*pow(s,2)) +
            eta1*eta2*(-2.*pow(Gammaa1,2)*pow(ma1,2) + 6.*pow(ma1,4) + 4.*pow(mpion,4) - 2.*pow(mrho,4) - 4.*pow(mrho,2)*s +
              4.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) + 8.*s)))*(1.*pow(mrho,2) - 1.*s - 1.*tmin) +
        0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) - 1.*pow(mrho,2) - 1.*s) + eta1*(pow(ma1,2) - 1.*pow(mrho,2) + s))*
          pow(1.*pow(mrho,2) - 1.*s - 1.*tmin,2) + 0.010416666666666666*pow(eta1 - 1.*eta2,4)*pow(1.*pow(mrho,2) - 1.*s - 1.*tmin,3) -
        (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/
          (1.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s - 1.*tmin) + (0.5*pow(-2. + delta,2)*pow(mpion,2)*tmin)/pow(mrho,2) +
        0.25*(-2. + delta)*(eta1 - 1.*eta2)*(-0.5*eta2*pow(ma1,2) + 1.*eta1*pow(mpion,2) - 1.*eta2*pow(mpion,2) + 0.5*eta1*s -
            0.5*eta2*s)*tmin - (0.25*(pow(mpion,2)*(12. + 1.*pow(delta,2) - 16.*C4*pow(mrho,2) + delta*(-8. + 8.*C4*pow(mrho,2))) +
              (-4. - 3.*pow(delta,2) - 16.*C4*pow(mrho,2) + delta*(8. + 8.*C4*pow(mrho,2)))*s)*tmin)/pow(mrho,2) -
        0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(1.*pow(mpion,2) - 0.5*s) +
            eta2*(-0.5*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 1.*s))*tmin -
        (0.25*(-2. + 1.*delta)*(-8.*C4*pow(mrho,4) + pow(mpion,2)*(2. + 1.*delta - 8.*C4*pow(mrho,2)) + (-2. - 3.*delta)*s +
              pow(mrho,2)*(2. + 1.*delta + 16.*C4*s))*tmin)/pow(mrho,2) -
        (0.25*(32.*pow(C4,2)*pow(mrho,8) + 2.*pow(delta,2)*pow(s,2) + 8.*C4*pow(mrho,6)*(-6. + delta - 8.*C4*s) +
              2.*delta*pow(mrho,2)*s*(-6. + delta - 8.*C4*s) +
              pow(mrho,4)*(12. - 1.*pow(delta,2) + 8.*C4*(6. + delta)*s + 32.*pow(C4,2)*pow(s,2)))*tmin)/pow(mrho,4) -
        0.09375*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-2.*pow(ma1,4) - 4.*pow(mpion,4) + 0.6666666666666666*pow(mrho,4) +
              pow(ma1,2)*(5.333333333333333*pow(mpion,2) - 2.6666666666666665*s) + 2.6666666666666665*pow(mpion,2)*s +
              1.3333333333333333*pow(mrho,2)*s - 1.3333333333333333*pow(s,2)) +
            pow(eta1,2)*(1.*pow(ma1,4) + 2.*pow(mpion,4) + 0.3333333333333333*pow(mrho,4) +
              pow(mpion,2)*(2.*pow(mrho,2) - 1.3333333333333333*s) - 1.3333333333333333*pow(mrho,2)*s + 0.6666666666666666*pow(s,2) +
              pow(ma1,2)*(-2.6666666666666665*pow(mpion,2) - 1.3333333333333333*pow(mrho,2) + 1.3333333333333333*s)) +
            pow(eta2,2)*(1.*pow(ma1,4) + 2.*pow(mpion,4) + 0.3333333333333333*pow(mrho,4) +
              pow(mpion,2)*(-2.*pow(mrho,2) - 1.3333333333333333*s) - 0.6666666666666666*pow(mrho,2)*s + 0.6666666666666666*pow(s,2) +
              pow(ma1,2)*(-2.6666666666666665*pow(mpion,2) + 1.3333333333333333*pow(mrho,2) + 1.3333333333333333*s)))*tmin -
        0.0625*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(pow(Gammaa1,2)*pow(ma1,2) - 1.*pow(ma1,4) - 2.*pow(mpion,4) -
              2.*pow(mpion,2)*pow(mrho,2) + 2.*pow(mrho,4) + pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) -
              3.*pow(mrho,2)*s + pow(s,2)) + pow(eta1,2)*
            (pow(Gammaa1,2)*pow(ma1,2) - 1.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) +
              pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) + pow(mrho,2)*s + pow(s,2)) +
            eta1*eta2*(-2.*pow(Gammaa1,2)*pow(ma1,2) + 2.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) +
              2.*pow(mrho,4) - 2.*pow(s,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)))*tmin +
        (1.*(pow(eta1,2)*(0.5*pow(mrho,4) - 1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) - 0.5*pow(mrho,2)*s + 0.25*delta*pow(mrho,2)*s -
                0.25*delta*pow(s,2) + 1.*C4*pow(mrho,2)*pow(s,2)) +
              pow(eta2,2)*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) - 0.5*pow(mrho,2)*s + 0.75*delta*pow(mrho,2)*s -
                2.*C4*pow(mrho,4)*s - 0.25*delta*pow(s,2) + 1.*C4*pow(mrho,2)*pow(s,2)) +
              eta1*eta2*(pow(ma1,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) + pow(mpion,2)*(-2.*pow(mrho,2) + 4.*C4*pow(mrho,4)) +
                s*(2.*C4*pow(mrho,4) + 0.5*delta*s + pow(mrho,2)*(1. - 1.*delta - 2.*C4*s))))*tmin)/pow(mrho,2) +
        (0.125*(-2. + 1.*delta)*(2. + 1.*delta - 8.*C4*pow(mrho,2))*pow(tmin,2))/pow(mrho,2) +
        (1.*(0.5 - 0.125*pow(delta,2) - 2.*C4*pow(mrho,2) + 1.*C4*delta*pow(mrho,2))*pow(tmin,2))/pow(mrho,2) +
        0.5*(eta1*eta2*(1. - 2.*C4*pow(mrho,2)) + pow(eta1,2)*(-0.5 + 1.*C4*pow(mrho,2)) + pow(eta2,2)*(-0.5 + 1.*C4*pow(mrho,2)))*
          pow(tmin,2) - 0.0625*pow(1.*eta1 - 1.*eta2,4)*(1.*pow(mpion,2) + 0.5*pow(mrho,2) - 0.5*s)*pow(tmin,2) -
        0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
            eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(tmin,2) -
        0.010416666666666666*pow(eta1 - 1.*eta2,4)*pow(tmin,3) + 0.020833333333333332*pow(1.*eta1 - 1.*eta2,4)*pow(tmin,3) +
        (2.*(1.*eta1 - 1.*eta2)*(eta2*(0.375*pow(mrho,4) - 0.0625*delta*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(mpion,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                0.125*pow(mrho,2)*s - 0.1875*delta*pow(mrho,2)*s + 1.*C4*pow(mrho,4)*s + 0.125*delta*pow(s,2) -
                0.5*C4*pow(mrho,2)*pow(s,2)) + eta1*(0.25*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                pow(mpion,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                0.125*delta*pow(s,2) + 0.5*C4*pow(mrho,2)*pow(s,2)))*
            (1.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s + 1.*tmin))/pow(mrho,2) -
        0.5*pow(1.*eta1 - 1.*eta2,2)*(-0.5 + 1.*C4*pow(mrho,2))*
          pow(1.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s + 1.*tmin,2) +
        (2.*(1.*eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*(0.375*pow(mrho,4) - 0.0625*delta*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(mpion,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                0.125*pow(mrho,2)*s - 0.1875*delta*pow(mrho,2)*s + 1.*C4*pow(mrho,4)*s + 0.125*delta*pow(s,2) -
                0.5*C4*pow(mrho,2)*pow(s,2)) + eta1*(0.25*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                pow(mpion,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                0.125*delta*pow(s,2) + 0.5*C4*pow(mrho,2)*pow(s,2)))*
            atan((0. + pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmax)/(Gammaa1*ma1)))/pow(mrho,2) +
        (0.25*(-2. + delta)*(eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*
              (-1.*pow(ma1,6) + pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(ma1,2) + 0.5*pow(mrho,2) - 1.*s) +
                pow(ma1,4)*(2.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) + pow(mpion,4)*(-1.5*pow(mrho,2) + 1.*s) +
                pow(ma1,2)*pow(mpion,2)*(-1.*pow(mpion,2) - 1.*pow(mrho,2) + 2.*s)) +
              eta1*(pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) + 0.5*s) + pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
                pow(ma1,2)*(-2.*pow(mpion,4) - 1.*pow(mpion,2)*s) +
                pow(mpion,2)*(1.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) + 1.*pow(mrho,2)*s)))*
            atan((0. + pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmax)/(Gammaa1*ma1)))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) -
        (0.25*(-2. + delta)*(eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*
              (-1.*pow(ma1,6) - 2.*pow(mpion,4)*pow(mrho,2) - 1.*pow(mpion,2)*pow(mrho,4) +
                pow(ma1,4)*(2.*pow(mpion,2) + 2.*pow(mrho,2) - 1.5*s) + 2.5*pow(mpion,4)*s + 3.*pow(mpion,2)*pow(mrho,2)*s +
                0.5*pow(mrho,4)*s - 2.*pow(mpion,2)*pow(s,2) - 1.*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) +
                pow(Gammaa1,2)*(-1.*pow(ma1,4) + 0.5*pow(ma1,2)*s) +
                pow(ma1,2)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) + 1.*pow(mrho,2)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 1.*s))) +
              eta1*(1.*pow(mpion,6) + 4.*pow(mpion,4)*pow(mrho,2) + 1.*pow(mpion,2)*pow(mrho,4) +
                pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) + pow(ma1,4)*(1.*pow(mpion,2) - 0.5*s) - 4.5*pow(mpion,4)*s -
                4.*pow(mpion,2)*pow(mrho,2)*s - 0.5*pow(mrho,4)*s + 3.*pow(mpion,2)*pow(s,2) + 1.*pow(mrho,2)*pow(s,2) -
                0.5*pow(s,3) + pow(ma1,2)*(-2.*pow(mpion,4) + (1.*pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 3.*s))))*
            atan((0. + pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmax)/(Gammaa1*ma1)))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
            2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
        (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8) + pow(mpion,8) -
                2.*pow(mpion,6)*pow(mrho,2) + pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                    2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s +
                    pow(mpion,4)*(6.*pow(mrho,2) + 2.*s)) +
                pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) - 6.*pow(mpion,4) - 1.*pow(mrho,4) +
                    pow(ma1,2)*(12.*pow(mpion,2) - 6.*pow(mrho,2) - 6.*s) + 2.*pow(mrho,2)*s - 2.*pow(s,2) +
                    pow(mpion,2)*(6.*pow(mrho,2) + 4.*s))) +
              eta1*eta2*(-2.*pow(Gammaa1,4)*pow(ma1,4) - 2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
                pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) + pow(ma1,2)*pow(mpion,2)*
                  (8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
                pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
                pow(Gammaa1,2)*pow(ma1,2)*(12.*pow(ma1,4) + 12.*pow(mpion,4) - 2.*pow(mrho,4) - 8.*pow(mpion,2)*s -
                    4.*pow(mrho,2)*s + 4.*pow(s,2) + pow(ma1,2)*(-24.*pow(mpion,2) + 12.*s))) +
              pow(eta1,2)*(pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                    2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
                    pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
                pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) - 6.*pow(mpion,4) - 1.*pow(mrho,4) +
                    pow(ma1,2)*(12.*pow(mpion,2) + 6.*pow(mrho,2) - 6.*s) + 4.*pow(mrho,2)*s - 2.*pow(s,2) +
                    pow(mpion,2)*(-6.*pow(mrho,2) + 4.*s)) +
                pow(mpion,2)*(pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                    pow(mpion,2)*(pow(mrho,4) - 2.*pow(mrho,2)*s))))*
            atan((0. + pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmax)/(Gammaa1*ma1)))/(Gammaa1*ma1) -
        (0.0625*pow(eta1 - 1.*eta2,2)*Gammaa1*ma1*(eta1*eta2*
              (-2.*pow(Gammaa1,4)*pow(ma1,4) + 14.*pow(ma1,8) + 14.*pow(mpion,8) + 28.*pow(mpion,6)*pow(mrho,2) +
                20.*pow(mpion,4)*pow(mrho,4) + 10.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) - 16.*pow(mpion,6)*s -
                16.*pow(mpion,4)*pow(mrho,2)*s - 12.*pow(mpion,2)*pow(mrho,4)*s - 4.*pow(mrho,6)*s - 4.*pow(mpion,4)*pow(s,2) -
                6.*pow(mpion,2)*pow(mrho,2)*pow(s,2) + 8.*pow(mpion,2)*pow(s,3) + 4.*pow(mrho,2)*pow(s,3) - 2.*pow(s,4) +
                pow(ma1,6)*(-56.*pow(mpion,2) - 28.*pow(mrho,2) + 28.*s) +
                pow(ma1,4)*(84.*pow(mpion,4) + 24.*pow(mrho,4) + pow(mpion,2)*(84.*pow(mrho,2) - 72.*s) - 36.*pow(mrho,2)*s +
                    12.*pow(s,2)) + pow(Gammaa1,2)*pow(ma1,2)*
                  (-4.*pow(ma1,4) - 4.*pow(mpion,4) + pow(ma1,2)*(8.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) +
                    (4.*pow(mrho,2) - 4.*s)*s + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)) +
                pow(ma1,2)*(-56.*pow(mpion,6) - 10.*pow(mrho,6) + 18.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3) +
                    pow(mpion,4)*(-84.*pow(mrho,2) + 60.*s) + pow(mpion,2)*(-48.*pow(mrho,4) + 60.*pow(mrho,2)*s - 12.*pow(s,2)))) +
              pow(eta1,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) - 7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
                7.*pow(mpion,4)*pow(mrho,4) - 2.*pow(mpion,2)*pow(mrho,6) +
                pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) + 8.*pow(mpion,6)*s + 11.*pow(mpion,4)*pow(mrho,2)*s +
                6.*pow(mpion,2)*pow(mrho,4)*s + 1.*pow(mrho,6)*s + 2.*pow(mpion,4)*pow(s,2) - 1.*pow(mrho,4)*pow(s,2) -
                4.*pow(mpion,2)*pow(s,3) - 1.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
                pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,4) + 2.*pow(mpion,4) + 1.*pow(mrho,4) +
                    pow(mpion,2)*(2.*pow(mrho,2) - 4.*s) - 1.*pow(mrho,2)*s + 2.*pow(s,2) +
                    pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
                pow(ma1,4)*(-42.*pow(mpion,4) - 9.*pow(mrho,4) + 21.*pow(mrho,2)*s - 6.*pow(s,2) +
                    pow(mpion,2)*(-42.*pow(mrho,2) + 36.*s)) +
                pow(ma1,2)*(28.*pow(mpion,6) + 2.*pow(mrho,6) + pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) - 9.*pow(mrho,4)*s +
                    6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(18.*pow(mrho,4) - 36.*pow(mrho,2)*s + 6.*pow(s,2)))) +
              pow(eta2,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) - 7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
                1.*pow(mpion,4)*pow(mrho,4) + 6.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) +
                pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) + 8.*pow(mpion,6)*s - 1.*pow(mpion,4)*pow(mrho,2)*s -
                16.*pow(mpion,2)*pow(mrho,4)*s - 7.*pow(mrho,6)*s + 2.*pow(mpion,4)*pow(s,2) +
                14.*pow(mpion,2)*pow(mrho,2)*pow(s,2) + 9.*pow(mrho,4)*pow(s,2) - 4.*pow(mpion,2)*pow(s,3) -
                5.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
                pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,4) + 2.*pow(mpion,4) + 3.*pow(mrho,4) +
                    pow(mpion,2)*(2.*pow(mrho,2) - 4.*s) - 5.*pow(mrho,2)*s + 2.*pow(s,2) +
                    pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
                pow(ma1,4)*(-42.*pow(mpion,4) - 3.*pow(mrho,4) + 9.*pow(mrho,2)*s - 6.*pow(s,2) +
                    pow(mpion,2)*(-42.*pow(mrho,2) + 36.*s)) +
                pow(ma1,2)*(28.*pow(mpion,6) - 4.*pow(mrho,6) + pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) + 9.*pow(mrho,4)*s -
                    6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s + 6.*pow(s,2)))))*
            atan((0. + pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmax)/(Gammaa1*ma1)))/
          (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
            4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) -
        (2.*(1.*eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*(0.375*pow(mrho,4) - 0.0625*delta*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(mpion,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                0.125*pow(mrho,2)*s - 0.1875*delta*pow(mrho,2)*s + 1.*C4*pow(mrho,4)*s + 0.125*delta*pow(s,2) -
                0.5*C4*pow(mrho,2)*pow(s,2)) + eta1*(0.25*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
                pow(mpion,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) + pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) -
                0.125*delta*pow(s,2) + 0.5*C4*pow(mrho,2)*pow(s,2)))*
            atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmin)/(Gammaa1*ma1)))/pow(mrho,2) -
        (0.25*(-2. + delta)*(eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*
              (-1.*pow(ma1,6) + pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(ma1,2) + 0.5*pow(mrho,2) - 1.*s) +
                pow(ma1,4)*(2.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) + pow(mpion,4)*(-1.5*pow(mrho,2) + 1.*s) +
                pow(ma1,2)*pow(mpion,2)*(-1.*pow(mpion,2) - 1.*pow(mrho,2) + 2.*s)) +
              eta1*(pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) + 0.5*s) + pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
                pow(ma1,2)*(-2.*pow(mpion,4) - 1.*pow(mpion,2)*s) +
                pow(mpion,2)*(1.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) + 1.*pow(mrho,2)*s)))*
            atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmin)/(Gammaa1*ma1)))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) +
        (0.25*(-2. + delta)*(eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*
              (-1.*pow(ma1,6) - 2.*pow(mpion,4)*pow(mrho,2) - 1.*pow(mpion,2)*pow(mrho,4) +
                pow(ma1,4)*(2.*pow(mpion,2) + 2.*pow(mrho,2) - 1.5*s) + 2.5*pow(mpion,4)*s + 3.*pow(mpion,2)*pow(mrho,2)*s +
                0.5*pow(mrho,4)*s - 2.*pow(mpion,2)*pow(s,2) - 1.*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) +
                pow(Gammaa1,2)*(-1.*pow(ma1,4) + 0.5*pow(ma1,2)*s) +
                pow(ma1,2)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) + 1.*pow(mrho,2)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 1.*s))) +
              eta1*(1.*pow(mpion,6) + 4.*pow(mpion,4)*pow(mrho,2) + 1.*pow(mpion,2)*pow(mrho,4) +
                pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) + pow(ma1,4)*(1.*pow(mpion,2) - 0.5*s) - 4.5*pow(mpion,4)*s -
                4.*pow(mpion,2)*pow(mrho,2)*s - 0.5*pow(mrho,4)*s + 3.*pow(mpion,2)*pow(s,2) + 1.*pow(mrho,2)*pow(s,2) -
                0.5*pow(s,3) + pow(ma1,2)*(-2.*pow(mpion,4) + (1.*pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 3.*s))))*
            atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmin)/(Gammaa1*ma1)))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
            2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) -
        (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8) + pow(mpion,8) -
                2.*pow(mpion,6)*pow(mrho,2) + pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                    2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s +
                    pow(mpion,4)*(6.*pow(mrho,2) + 2.*s)) +
                pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) - 6.*pow(mpion,4) - 1.*pow(mrho,4) +
                    pow(ma1,2)*(12.*pow(mpion,2) - 6.*pow(mrho,2) - 6.*s) + 2.*pow(mrho,2)*s - 2.*pow(s,2) +
                    pow(mpion,2)*(6.*pow(mrho,2) + 4.*s))) +
              eta1*eta2*(-2.*pow(Gammaa1,4)*pow(ma1,4) - 2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
                pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) + pow(ma1,2)*pow(mpion,2)*
                  (8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
                pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
                pow(Gammaa1,2)*pow(ma1,2)*(12.*pow(ma1,4) + 12.*pow(mpion,4) - 2.*pow(mrho,4) - 8.*pow(mpion,2)*s -
                    4.*pow(mrho,2)*s + 4.*pow(s,2) + pow(ma1,2)*(-24.*pow(mpion,2) + 12.*s))) +
              pow(eta1,2)*(pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                    2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
                    pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
                pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) - 6.*pow(mpion,4) - 1.*pow(mrho,4) +
                    pow(ma1,2)*(12.*pow(mpion,2) + 6.*pow(mrho,2) - 6.*s) + 4.*pow(mrho,2)*s - 2.*pow(s,2) +
                    pow(mpion,2)*(-6.*pow(mrho,2) + 4.*s)) +
                pow(mpion,2)*(pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                    pow(mpion,2)*(pow(mrho,4) - 2.*pow(mrho,2)*s))))*
            atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmin)/(Gammaa1*ma1)))/(Gammaa1*ma1) +
        (0.0625*pow(eta1 - 1.*eta2,2)*Gammaa1*ma1*(eta1*eta2*
              (-2.*pow(Gammaa1,4)*pow(ma1,4) + 14.*pow(ma1,8) + 14.*pow(mpion,8) + 28.*pow(mpion,6)*pow(mrho,2) +
                20.*pow(mpion,4)*pow(mrho,4) + 10.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) - 16.*pow(mpion,6)*s -
                16.*pow(mpion,4)*pow(mrho,2)*s - 12.*pow(mpion,2)*pow(mrho,4)*s - 4.*pow(mrho,6)*s - 4.*pow(mpion,4)*pow(s,2) -
                6.*pow(mpion,2)*pow(mrho,2)*pow(s,2) + 8.*pow(mpion,2)*pow(s,3) + 4.*pow(mrho,2)*pow(s,3) - 2.*pow(s,4) +
                pow(ma1,6)*(-56.*pow(mpion,2) - 28.*pow(mrho,2) + 28.*s) +
                pow(ma1,4)*(84.*pow(mpion,4) + 24.*pow(mrho,4) + pow(mpion,2)*(84.*pow(mrho,2) - 72.*s) - 36.*pow(mrho,2)*s +
                    12.*pow(s,2)) + pow(Gammaa1,2)*pow(ma1,2)*
                  (-4.*pow(ma1,4) - 4.*pow(mpion,4) + pow(ma1,2)*(8.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) +
                    (4.*pow(mrho,2) - 4.*s)*s + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)) +
                pow(ma1,2)*(-56.*pow(mpion,6) - 10.*pow(mrho,6) + 18.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3) +
                    pow(mpion,4)*(-84.*pow(mrho,2) + 60.*s) + pow(mpion,2)*(-48.*pow(mrho,4) + 60.*pow(mrho,2)*s - 12.*pow(s,2)))) +
              pow(eta1,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) - 7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
                7.*pow(mpion,4)*pow(mrho,4) - 2.*pow(mpion,2)*pow(mrho,6) +
                pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) + 8.*pow(mpion,6)*s + 11.*pow(mpion,4)*pow(mrho,2)*s +
                6.*pow(mpion,2)*pow(mrho,4)*s + 1.*pow(mrho,6)*s + 2.*pow(mpion,4)*pow(s,2) - 1.*pow(mrho,4)*pow(s,2) -
                4.*pow(mpion,2)*pow(s,3) - 1.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
                pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,4) + 2.*pow(mpion,4) + 1.*pow(mrho,4) +
                    pow(mpion,2)*(2.*pow(mrho,2) - 4.*s) - 1.*pow(mrho,2)*s + 2.*pow(s,2) +
                    pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
                pow(ma1,4)*(-42.*pow(mpion,4) - 9.*pow(mrho,4) + 21.*pow(mrho,2)*s - 6.*pow(s,2) +
                    pow(mpion,2)*(-42.*pow(mrho,2) + 36.*s)) +
                pow(ma1,2)*(28.*pow(mpion,6) + 2.*pow(mrho,6) + pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) - 9.*pow(mrho,4)*s +
                    6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(18.*pow(mrho,4) - 36.*pow(mrho,2)*s + 6.*pow(s,2)))) +
              pow(eta2,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) - 7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
                1.*pow(mpion,4)*pow(mrho,4) + 6.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) +
                pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) + 8.*pow(mpion,6)*s - 1.*pow(mpion,4)*pow(mrho,2)*s -
                16.*pow(mpion,2)*pow(mrho,4)*s - 7.*pow(mrho,6)*s + 2.*pow(mpion,4)*pow(s,2) +
                14.*pow(mpion,2)*pow(mrho,2)*pow(s,2) + 9.*pow(mrho,4)*pow(s,2) - 4.*pow(mpion,2)*pow(s,3) -
                5.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
                pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,4) + 2.*pow(mpion,4) + 3.*pow(mrho,4) +
                    pow(mpion,2)*(2.*pow(mrho,2) - 4.*s) - 5.*pow(mrho,2)*s + 2.*pow(s,2) +
                    pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
                pow(ma1,4)*(-42.*pow(mpion,4) - 3.*pow(mrho,4) + 9.*pow(mrho,2)*s - 6.*pow(s,2) +
                    pow(mpion,2)*(-42.*pow(mrho,2) + 36.*s)) +
                pow(ma1,2)*(28.*pow(mpion,6) - 4.*pow(mrho,6) + pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) + 9.*pow(mrho,4)*s -
                    6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s + 6.*pow(s,2)))))*
            atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmin)/(Gammaa1*ma1)))/
          (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
            4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
        0.0625*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
              pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) - 2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
              pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
            pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) + pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s +
              pow(mpion,4)*(-3.*pow(mrho,2) + s) + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
              pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2)))\
            + pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
              pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
              pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2)))
            )*log(fabs(-pow(ma1,2) + tmax)) - (0.03125*pow(eta1 - 1.*eta2,2)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s)*
            (eta1*eta2*(-2.*pow(ma1,8) + pow(ma1,6)*(8.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) +
                pow(ma1,4)*(-12.*pow(mpion,4) - 4.*pow(mrho,4) + 4.*pow(mrho,2)*s + pow(mpion,2)*(-12.*pow(mrho,2) + 8.*s)) +
                pow(mpion,2)*(-2.*pow(mpion,6) - 4.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 4.*pow(mrho,4)*s -
                    2.*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(-8.*pow(mrho,4) + 8.*pow(mrho,2)*s)) +
                pow(ma1,2)*(8.*pow(mpion,6) + 2.*pow(mrho,6) + pow(mpion,4)*(12.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,4)*s -
                    2.*pow(mrho,2)*pow(s,2) + 2.*pow(s,3) + pow(mpion,2)*(8.*pow(mrho,4) - 4.*pow(mrho,2)*s - 4.*pow(s,2)))) +
              pow(eta2,2)*(pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                pow(mpion,4)*(pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 1.*pow(mrho,2)*s) +
                pow(ma1,4)*(6.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) + pow(mrho,2)*s) +
                pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mrho,6) - 5.*pow(mrho,4)*s + 4.*pow(mrho,2)*pow(s,2) - 1.*pow(s,3) +
                    pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) + pow(mpion,2)*(2.*pow(mrho,4) - 4.*pow(mrho,2)*s + 2.*pow(s,2)))) +
              pow(eta1,2)*(pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 3.*pow(mrho,2)*s) +
                pow(mpion,2)*(pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) + pow(mrho,2)*s*(-2.*pow(mrho,2) + 2.*s) +
                    pow(mpion,2)*(3.*pow(mrho,4) - 5.*pow(mrho,2)*s)) +
                pow(ma1,2)*(-4.*pow(mpion,6) + pow(mrho,4)*s - 1.*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) +
                    pow(mpion,2)*(-2.*pow(mrho,4) + 4.*pow(mrho,2)*s + 2.*pow(s,2)))))*log(fabs(-pow(ma1,2) + tmax)))/
          (0.25*pow(Gammaa1,2)*pow(ma1,2) + 1.*pow(ma1,4) + 1.*pow(mpion,4) + 1.*pow(mpion,2)*pow(mrho,2) + 0.25*pow(mrho,4) -
            1.*pow(mpion,2)*s - 0.5*pow(mrho,2)*s + 0.25*pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s)) -
        (1.*(eta1*eta2*(pow(ma1,8)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                pow(ma1,6)*(-1.*pow(mrho,4) + 2.*C4*pow(mrho,6) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*C4*pow(mrho,4)) +
                    0.5*delta*pow(s,2) + pow(mrho,2)*s*(2. - 1.*delta - 2.*C4*s)) +
                pow(ma1,2)*(pow(mpion,6)*(-4.*pow(mrho,2) + 8.*C4*pow(mrho,4)) +
                    pow(mrho,4)*s*((0.5 - 0.25*delta)*pow(mrho,2) + (-0.5 + 0.25*delta)*s) +
                    pow(mpion,4)*(10.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,4)*(-3. - 1.*delta - 8.*C4*s) +
                      pow(mrho,2)*s*(2. + 1.*delta - 2.*C4*s)) +
                    pow(mpion,2)*(2.*C4*pow(mrho,8) - 0.5*delta*pow(s,3) + pow(mrho,6)*(1. - 1.*delta - 2.*C4*s) +
                      pow(mrho,4)*s*(-1. + 1.*delta - 2.*C4*s) + pow(mrho,2)*pow(s,2)*(1. + 2.*C4*s))) +
                pow(mpion,2)*pow(mrho,2)*(pow(mpion,6)*(1. - 2.*C4*pow(mrho,2)) +
                    pow(mrho,4)*((-0.5 + 0.25*delta)*pow(mrho,2) + (0.5 - 0.25*delta)*s) +
                    pow(mpion,2)*(-2.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,2)*s*(-1.5 - 0.25*delta - 2.*C4*s) +
                      pow(mrho,4)*(1. + 4.*C4*s)) + pow(mpion,4)*
                    (-4.*C4*pow(mrho,4) - 1.*delta*s + pow(mrho,2)*(1. + 0.5*delta + 4.*C4*s))) +
                pow(ma1,4)*(pow(mpion,4)*(6.*pow(mrho,2) - 12.*C4*pow(mrho,4)) +
                    pow(mpion,2)*(-8.*C4*pow(mrho,6) - 1.*delta*pow(s,2) + pow(mrho,4)*(3. + 0.5*delta + 4.*C4*s) +
                      pow(mrho,2)*s*(-4. + 1.*delta + 4.*C4*s)) +
                    s*(-2.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,2)*s*(1. - 1.5*delta - 2.*C4*s) +
                      pow(mrho,4)*(-1.5 + 1.25*delta + 4.*C4*s)))) +
              pow(eta1,2)*(pow(ma1,8)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                pow(ma1,6)*(-2.*C4*pow(mrho,6) + pow(mpion,2)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) - 0.25*delta*pow(s,2) +
                    pow(mrho,4)*(1. + 1.*C4*s) + pow(mrho,2)*s*(-1. + 0.25*delta + 1.*C4*s)) +
                pow(ma1,4)*(1.*C4*pow(mrho,8) + pow(mpion,4)*(-3.*pow(mrho,2) + 6.*C4*pow(mrho,4)) - 0.25*delta*pow(s,3) +
                    pow(mrho,6)*(-0.5 - 1.*C4*s) + pow(mrho,4)*s*(1.5 - 0.5*delta - 1.*C4*s) +
                    pow(mrho,2)*pow(s,2)*(-0.5 + 0.5*delta + 1.*C4*s) +
                    pow(mpion,2)*(7.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,4)*(-3. - 0.25*delta - 5.*C4*s) +
                      pow(mrho,2)*s*(2. + 0.25*delta - 2.*C4*s))) +
                pow(mpion,2)*pow(mrho,2)*(pow(mpion,6)*(-0.5 + 1.*C4*pow(mrho,2)) +
                    pow(mrho,4)*((0.5 - 0.25*delta)*pow(mrho,2) + (-0.5 + 0.25*delta)*s) +
                    pow(mpion,4)*(3.*C4*pow(mrho,4) + 0.75*delta*s + pow(mrho,2)*(-1. - 0.25*delta - 3.*C4*s)) +
                    pow(mpion,2)*(2.*C4*pow(mrho,6) - 0.5*delta*pow(s,2) + pow(mrho,4)*(-1. - 4.*C4*s) +
                      pow(mrho,2)*s*(1.5 + 0.25*delta + 2.*C4*s))) +
                pow(ma1,2)*(pow(mpion,6)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                    pow(mrho,4)*s*((-0.5 + 0.25*delta)*pow(mrho,2) + (0.5 - 0.25*delta)*s) +
                    pow(mpion,2)*(-3.*C4*pow(mrho,8) + 0.25*delta*pow(s,3) + pow(mrho,4)*s*(-1. - 0.75*delta - 1.*C4*s) +
                      pow(mrho,2)*pow(s,2)*(-0.5 + 0.5*delta - 1.*C4*s) + pow(mrho,6)*(0.5 + 0.5*delta + 5.*C4*s)) +
                    pow(mpion,4)*(-8.*C4*pow(mrho,6) - 0.25*delta*pow(s,2) + pow(mrho,2)*s*(-1. - 1.25*delta + 1.*C4*s) +
                      pow(mrho,4)*(3. + 0.5*delta + 7.*C4*s)))) +
              pow(eta2,2)*(pow(ma1,8)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                pow(mpion,6)*pow(mrho,2)*(-0.25*delta*pow(mrho,2) + 1.*C4*pow(mrho,4) + pow(mpion,2)*(-0.5 + 1.*C4*pow(mrho,2)) +
                    0.25*delta*s - 1.*C4*pow(mrho,2)*s) +
                pow(ma1,6)*(pow(mpion,2)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                    s*(-1.*C4*pow(mrho,4) - 0.25*delta*s + pow(mrho,2)*(-1. + 0.75*delta + 1.*C4*s))) +
                pow(ma1,4)*(-1.*C4*pow(mrho,8) + pow(mpion,4)*(-3.*pow(mrho,2) + 6.*C4*pow(mrho,4)) - 0.25*delta*pow(s,3) +
                    pow(mrho,4)*s*(-0.75*delta - 3.*C4*s) + pow(mrho,2)*pow(s,2)*(-0.5 + 1.*delta + 1.*C4*s) +
                    pow(mrho,6)*(0.5 + 3.*C4*s) + pow(mpion,2)*
                    (delta*(-0.25*pow(mrho,4) - 1.25*pow(mrho,2)*s + 0.5*pow(s,2)) +
                      pow(mrho,2)*(1.*C4*pow(mrho,4) + 2.*s + 1.*C4*pow(mrho,2)*s - 2.*C4*pow(s,2)))) +
                pow(ma1,2)*pow(mpion,2)*(1.*C4*pow(mrho,8) + pow(mpion,4)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                    0.25*delta*pow(s,3) + pow(mrho,6)*(-1.5 + 0.5*delta - 3.*C4*s) + pow(mrho,2)*pow(s,2)*(-0.5 - 0.5*delta - 1.*C4*s) +
                    pow(mrho,4)*s*(2. - 0.25*delta + 3.*C4*s) +
                    pow(mpion,2)*(delta*(0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
                      pow(mrho,2)*(-2.*C4*pow(mrho,4) - 1.*s + 1.*C4*pow(mrho,2)*s + 1.*C4*pow(s,2))))))*log(fabs(-pow(ma1,2) + tmax)))/
          ((pow(ma1,2) - 1.*pow(mpion,2))*pow(mrho,2)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s)) +
        (0.5*pow(mpion,2)*(pow(eta2,2)*pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s +
                (-2. + 1.*delta)*pow(s,2)) + eta1*eta2*(pow(mpion,2)*
                  ((4. - 2.*delta)*pow(mrho,4) + (-8. + 4.*delta)*pow(mrho,2)*s + (4. - 2.*delta)*pow(s,2)) +
                pow(mrho,2)*((-1. + 0.5*delta)*pow(mrho,4) + (2. - 1.*delta)*pow(mrho,2)*s + (-1. + 0.5*delta)*pow(s,2))) +
              pow(eta1,2)*(pow(mrho,2)*((1. - 0.5*delta)*pow(mrho,4) + (-2. + 1.*delta)*pow(mrho,2)*s + (1. - 0.5*delta)*pow(s,2)) +
                pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))))*
            log(fabs(-pow(mpion,2) + tmax)))/((-1.*pow(ma1,2) + 1.*pow(mpion,2))*(pow(mrho,2) - 1.*s)) -
        (0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*
            (eta1*(pow(mpion,4)*(-1.*pow(mrho,2) + 1.*s) +
                pow(ma1,2)*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + (-0.5*pow(mrho,2) + 0.5*s)*s) +
                pow(mpion,2)*(-1.*pow(mrho,4) + 2.5*pow(mrho,2)*s - 1.5*pow(s,2)) +
                s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2))) +
              eta2*(0.5*pow(mrho,6) + pow(mpion,4)*(1.*pow(mrho,2) - 1.*s) - 1.5*pow(mrho,4)*s + 1.5*pow(mrho,2)*pow(s,2) -
                0.5*pow(s,3) + pow(mpion,2)*(1.5*pow(mrho,4) - 3.*pow(mrho,2)*s + 1.5*pow(s,2)) +
                pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*pow(mrho,2)*s - 0.5*pow(s,2) + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*
            log(fabs(-pow(mpion,2) + tmax)))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
            2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) -
        (0.5*pow(mpion,2)*(eta1*eta2*((1. - 0.5*delta)*pow(mrho,6) + (-4. + 2.*delta)*pow(mrho,4)*s +
                (5. - 2.5*delta)*pow(mrho,2)*pow(s,2) + (-2. + 1.*delta)*pow(s,3) +
                pow(mpion,2)*((4. - 2.*delta)*pow(mrho,4) + (-8. + 4.*delta)*pow(mrho,2)*s + (4. - 2.*delta)*pow(s,2))) +
              pow(eta2,2)*((-1. + 0.5*delta)*pow(mrho,6) + (3. - 1.5*delta)*pow(mrho,4)*s + (-3. + 1.5*delta)*pow(mrho,2)*pow(s,2) +
                (1. - 0.5*delta)*pow(s,3) + pow(mpion,2)*
                  ((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))) +
              pow(eta1,2)*(s*((1. - 0.5*delta)*pow(mrho,4) + (-2. + 1.*delta)*pow(mrho,2)*s + (1. - 0.5*delta)*pow(s,2)) +
                pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))))*
            log(fabs(-pow(mpion,2) - pow(mrho,2) + s + tmax)))/
          ((pow(mrho,2) - 1.*s)*(-1.*pow(ma1,2) + pow(mpion,2) + pow(mrho,2) - 1.*s)) +
        0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + tmax)) +
        (0.25000000000000006*pow(2. - 1.*delta,2)*(7.999999999999998*pow(mpion,4) - 5.999999999999998*pow(mpion,2)*s + 1.*pow(s,2))*
            log(fabs(-pow(mpion,2) - pow(mrho,2) + s + tmax)))/(pow(mrho,2) - 1.*s) -
        (1.*(-2. + 1.*delta)*((0.5 - 0.25*delta)*pow(mrho,2)*s +
              pow(mpion,2)*(4.*C4*pow(mrho,4) + 1.*delta*s + pow(mrho,2)*(-2. - 4.*C4*s)))*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + tmax)))/pow(mrho,2)
          - 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + tmax)) -
        (2.*(1.*pow(2. - 1.*delta,2)*pow(mpion,4)*pow(mrho,2) + 0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,4)*s +
              pow(mpion,2)*(C4*(4. - 2.*delta)*pow(mrho,6) + (-1. + 0.5*delta)*delta*pow(s,2) +
                pow(mrho,2)*s*(-1. + 3.*delta - 1.25*pow(delta,2) + 4.*C4*s - 2.*C4*delta*s) +
                pow(mrho,4)*(-2. + 1.*delta - 8.*C4*s + 4.*C4*delta*s)))*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + tmax)))/
          (pow(mrho,4) - 1.*pow(mrho,2)*s) + (0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*
            (eta2*pow(mpion,2)*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + pow(ma1,2)*(-1.*pow(mrho,2) + 1.*s)) +
              eta1*(pow(ma1,2)*(-0.5*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s) +
                pow(mpion,2)*(0.5*pow(mrho,4) - 0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*
            log(fabs(-pow(mpion,2) - pow(mrho,2) + s + tmax)))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) +
        (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8) + 0.5*pow(mpion,8) +
                0.5*pow(ma1,4)*pow(mpion,2)*pow(mrho,2) - 0.5*pow(mpion,6)*pow(mrho,2) +
                pow(Gammaa1,2)*pow(ma1,2)*(pow(mpion,2)*(1.*pow(mpion,2) + 1.5*pow(mrho,2) - 2.*s) +
                    pow(ma1,2)*(-1.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s)) +
                pow(ma1,6)*(1.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) +
                pow(ma1,2)*pow(mpion,4)*(-1.*pow(mpion,2) - 0.5*pow(mrho,2) + 1.*s)) +
              eta1*(pow(ma1,6)*(1.*pow(mpion,2) + 0.5*s) +
                pow(ma1,2)*(3.*pow(mpion,6) + 1.*pow(mpion,2)*pow(mrho,4) - 0.5*pow(mpion,4)*s) +
                pow(ma1,4)*(-3.*pow(mpion,4) + pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
                pow(mpion,4)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 0.5*s) + 0.5*pow(mrho,2)*s) +
                pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) + pow(ma1,2)*(1.*pow(mpion,2) + 0.5*s) - 0.5*pow(mrho,2)*s +
                    pow(mpion,2)*(-1.*pow(mrho,2) + 1.5*s))))*
            log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*tmax - 2.*pow(mrho,2)*tmax + 2.*s*tmax + pow(tmax,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*tmax))))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) -
        (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(-1.*pow(mpion,8) + 1.*pow(mpion,4)*pow(mrho,4) +
                pow(ma1,6)*(1.*pow(mpion,2) - 0.5*s) - 0.5*pow(mpion,6)*s - 4.*pow(mpion,4)*pow(mrho,2)*s -
                1.5*pow(mpion,2)*pow(mrho,4)*s + 3.5*pow(mpion,4)*pow(s,2) + 4.*pow(mpion,2)*pow(mrho,2)*pow(s,2) +
                0.5*pow(mrho,4)*pow(s,2) - 2.5*pow(mpion,2)*pow(s,3) - 1.*pow(mrho,2)*pow(s,3) + 0.5*pow(s,4) +
                pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) + pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) - 0.5*pow(mpion,2)*s +
                    0.5*pow(s,2)) + pow(ma1,4)*(-3.*pow(mpion,4) + (1.*pow(mrho,2) - 0.5*s)*s +
                    pow(mpion,2)*(-2.*pow(mrho,2) + 2.5*s)) +
                pow(ma1,2)*(3.*pow(mpion,6) + pow(mpion,4)*(2.*pow(mrho,2) - 1.5*s) - 0.5*pow(mrho,4)*s + 0.5*pow(s,3) +
                    pow(mpion,2)*(1.*pow(mrho,4) - 1.*pow(mrho,2)*s - 1.*pow(s,2)))) +
              eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8) + pow(ma1,6)*(1.*pow(mpion,2) + 1.*pow(mrho,2) - 0.5*s) +
                pow(ma1,4)*(-0.5*pow(mrho,4) + (-0.5*pow(mpion,2) + 0.5*s)*s) +
                pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) -
                    1.*pow(mrho,2)*s + 0.5*pow(s,2) + pow(ma1,2)*(-1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.5*s)) +
                pow(mpion,2)*(0.5*pow(mpion,6) + 0.5*pow(mpion,4)*s +
                    pow(mpion,2)*(-0.5*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.5*pow(s,2)) +
                    s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2))) +
                pow(ma1,2)*(-1.*pow(mpion,6) + pow(mpion,4)*(-1.*pow(mrho,2) + 0.5*s) +
                    pow(mpion,2)*(-1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
                    s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2)))))*
            log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*tmax - 2.*pow(mrho,2)*tmax + 2.*s*tmax + pow(tmax,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*tmax))))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
            2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) -
        (0.0625*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(-1.*pow(ma1,10) + pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
                pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
                pow(ma1,4)*(10.*pow(mpion,6) - 2.5*pow(mrho,6) + pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) + 6.*pow(mrho,4)*s -
                    4.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3)) +
                pow(ma1,6)*(-10.*pow(mpion,4) + (1.*pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-10.*pow(mrho,2) + 8.*s)) +
                pow(mpion,4)*(1.*pow(mpion,6) + 0.5*pow(mrho,6) + pow(mpion,4)*(2.5*pow(mrho,2) - 0.5*s) - 1.*pow(mrho,4)*s +
                    0.5*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(2.*pow(mrho,4) - 2.*pow(mrho,2)*s)) +
                pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) - 4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
                    1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
                    pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
                    pow(mpion,2)*(-3.*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.*pow(s,2)) +
                    pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4) + pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
                      3.*pow(s,2))) + pow(ma1,2)*(-5.*pow(mpion,8) + 1.*pow(mrho,8) - 3.5*pow(mrho,6)*s +
                    4.5*pow(mrho,4)*pow(s,2) - 2.5*pow(mrho,2)*pow(s,3) + 0.5*pow(s,4) + pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
                    pow(mpion,4)*(-2.*pow(mrho,4) + 1.*pow(mrho,2)*s + 1.*pow(s,2)) +
                    pow(mpion,2)*(3.*pow(mrho,6) - 8.*pow(mrho,4)*s + 7.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)))) +
              pow(eta1,2)*(-1.*pow(ma1,10) + pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
                pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
                pow(ma1,6)*(-10.*pow(mpion,4) - 2.*pow(mrho,4) + 5.*pow(mrho,2)*s - 1.*pow(s,2) +
                    pow(mpion,2)*(-10.*pow(mrho,2) + 8.*s)) +
                pow(ma1,4)*(10.*pow(mpion,6) + 0.5*pow(mrho,6) + pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) - 3.*pow(mrho,4)*s +
                    1.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s)) +
                pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) - 4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
                    1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
                    pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
                    pow(mpion,2)*(-3.*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.*pow(s,2)) +
                    pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4) + pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
                      3.*pow(s,2))) + pow(mpion,2)*(1.*pow(mpion,8) + pow(mpion,6)*(2.5*pow(mrho,2) - 0.5*s) +
                    pow(mpion,4)*(4.*pow(mrho,4) - 6.*pow(mrho,2)*s) +
                    pow(mrho,2)*s*(-1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
                    pow(mpion,2)*(1.5*pow(mrho,6) - 6.*pow(mrho,4)*s + 4.5*pow(mrho,2)*pow(s,2))) +
                pow(ma1,2)*(-5.*pow(mpion,8) + pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
                    pow(mpion,4)*(-8.*pow(mrho,4) + 13.*pow(mrho,2)*s + 1.*pow(s,2)) +
                    pow(mpion,2)*(-1.*pow(mrho,6) + 6.*pow(mrho,4)*s - 3.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)) +
                    s*(0.5*pow(mrho,6) - 0.5*pow(mrho,4)*s - 0.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3)))) +
              eta1*eta2*(2.*pow(ma1,10) + pow(Gammaa1,4)*pow(ma1,4)*(-2.*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s) +
                pow(ma1,8)*(-10.*pow(mpion,2) - 5.*pow(mrho,2) + 5.*s) +
                pow(ma1,6)*(20.*pow(mpion,4) + 6.*pow(mrho,4) + pow(mpion,2)*(20.*pow(mrho,2) - 16.*s) - 8.*pow(mrho,2)*s +
                    2.*pow(s,2)) + pow(ma1,4)*(-20.*pow(mpion,6) - 4.*pow(mrho,6) + 6.*pow(mrho,4)*s - 2.*pow(s,3) +
                    pow(mpion,4)*(-30.*pow(mrho,2) + 18.*s) + pow(mpion,2)*(-18.*pow(mrho,4) + 18.*pow(mrho,2)*s)) +
                pow(mpion,2)*(-2.*pow(mpion,8) - 1.*pow(mrho,8) + 3.*pow(mrho,6)*s - 3.*pow(mrho,4)*pow(s,2) +
                    1.*pow(mrho,2)*pow(s,3) + pow(mpion,6)*(-5.*pow(mrho,2) + 1.*s) +
                    pow(mpion,4)*(-10.*pow(mrho,4) + 10.*pow(mrho,2)*s) +
                    pow(mpion,2)*(-6.*pow(mrho,6) + 12.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2))) +
                pow(ma1,2)*(10.*pow(mpion,8) + 1.*pow(mrho,8) + pow(mpion,6)*(20.*pow(mrho,2) - 8.*s) - 2.*pow(mrho,6)*s +
                    2.*pow(mrho,2)*pow(s,3) - 1.*pow(s,4) + pow(mpion,4)*(22.*pow(mrho,4) - 20.*pow(mrho,2)*s - 2.*pow(s,2)) +
                    pow(mpion,2)*(8.*pow(mrho,6) - 12.*pow(mrho,4)*s + 4.*pow(s,3))) +
                pow(Gammaa1,2)*pow(ma1,2)*(-8.*pow(ma1,6) + 8.*pow(mpion,6) + 1.*pow(mrho,6) +
                    pow(mpion,4)*(12.*pow(mrho,2) - 12.*s) + pow(ma1,4)*(24.*pow(mpion,2) + 12.*pow(mrho,2) - 12.*s) -
                    3.*pow(mrho,4)*s + 3.*pow(mrho,2)*pow(s,2) - 1.*pow(s,3) +
                    pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s + 6.*pow(s,2)) +
                    pow(ma1,2)*(-24.*pow(mpion,4) - 6.*pow(mrho,4) + 12.*pow(mrho,2)*s - 6.*pow(s,2) +
                      pow(mpion,2)*(-24.*pow(mrho,2) + 24.*s)))))*
            log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*tmax - 2.*pow(mrho,2)*tmax + 2.*s*tmax + pow(tmax,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*tmax))))/
          (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
            4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
        0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(4.*pow(ma1,6) +
              pow(Gammaa1,2)*pow(ma1,2)*(-4.*pow(ma1,2) + 4.*pow(mpion,2) - 2.*s) + pow(ma1,4)*(-12.*pow(mpion,2) + 6.*s) +
              pow(mpion,2)*(-4.*pow(mpion,4) + 4.*pow(mrho,4) + 2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s) +
              pow(ma1,2)*(12.*pow(mpion,4) - 2.*pow(mrho,4) - 8.*pow(mpion,2)*s - 4.*pow(mrho,2)*s + 4.*pow(s,2))) +
            pow(eta1,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) + 3.*pow(mpion,4)*pow(mrho,2) +
              pow(ma1,4)*(6.*pow(mpion,2) + 3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s - 1.*pow(mpion,2)*pow(mrho,2)*s -
              1.*pow(mrho,4)*s + pow(mrho,2)*pow(s,2) +
              pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s) +
              pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) + 4.*pow(mrho,2)*s - 2.*pow(s,2) +
                  pow(mpion,2)*(-6.*pow(mrho,2) + 4.*s))) +
            pow(eta2,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) - 3.*pow(mpion,4)*pow(mrho,2) +
              pow(ma1,4)*(6.*pow(mpion,2) - 3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s + pow(mpion,2)*pow(mrho,2)*s +
              pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) - 2.*pow(mpion,2) + pow(mrho,2) + s) +
              pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 2.*pow(s,2) +
                  pow(mpion,2)*(6.*pow(mrho,2) + 4.*s))))*
          log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + 1.*pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*tmax - 2.*pow(mrho,2)*tmax + 2.*s*tmax + pow(tmax,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*tmax))) -
        (0.5*(1.*eta1 - 1.*eta2)*(eta2*(pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(0.5 - 1.*C4*pow(mrho,2)) +
                pow(ma1,4)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                pow(mpion,2)*pow(mrho,2)*(pow(mpion,2)*(-0.5 + 1.*C4*pow(mrho,2)) + (0.25 - 0.125*delta)*(pow(mrho,2) + s)) +
                pow(ma1,2)*(1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) - 0.25*delta*pow(s,2) +
                    pow(mrho,4)*(-0.75 + 0.125*delta - 2.*C4*s) + pow(mrho,2)*s*(0.25 + 0.375*delta + 1.*C4*s))) +
              eta1*(pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(-0.5 + 1.*C4*pow(mrho,2)) +
                pow(ma1,4)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
                pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6) + pow(mpion,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) +
                    0.25*delta*pow(s,2) - 1.*C4*pow(mrho,2)*pow(s,2)) +
                pow(mrho,2)*(pow(mpion,4)*(0.5 - 1.*C4*pow(mrho,2)) + s*((-0.25 + 0.125*delta)*pow(mrho,2) + (0.25 - 0.125*delta)*s) +
                    pow(mpion,2)*(-2.*C4*pow(mrho,4) + (-0.5 - 0.25*delta)*s + pow(mrho,2)*(1. + 2.*C4*s)))))*
            log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + 1.*pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*tmax - 2.*pow(mrho,2)*tmax + 2.*s*tmax + pow(tmax,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*tmax))))/pow(mrho,2) -
        0.0625*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
              pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) - 2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
              pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
            pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) + pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s +
              pow(mpion,4)*(-3.*pow(mrho,2) + s) + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
              pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2)))\
            + pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
              pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
              pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2)))
            )*log(fabs(-pow(ma1,2) + tmin)) + (0.03125*pow(eta1 - 1.*eta2,2)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s)*
            (eta1*eta2*(-2.*pow(ma1,8) + pow(ma1,6)*(8.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) +
                pow(ma1,4)*(-12.*pow(mpion,4) - 4.*pow(mrho,4) + 4.*pow(mrho,2)*s + pow(mpion,2)*(-12.*pow(mrho,2) + 8.*s)) +
                pow(mpion,2)*(-2.*pow(mpion,6) - 4.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 4.*pow(mrho,4)*s -
                    2.*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(-8.*pow(mrho,4) + 8.*pow(mrho,2)*s)) +
                pow(ma1,2)*(8.*pow(mpion,6) + 2.*pow(mrho,6) + pow(mpion,4)*(12.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,4)*s -
                    2.*pow(mrho,2)*pow(s,2) + 2.*pow(s,3) + pow(mpion,2)*(8.*pow(mrho,4) - 4.*pow(mrho,2)*s - 4.*pow(s,2)))) +
              pow(eta2,2)*(pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                pow(mpion,4)*(pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 1.*pow(mrho,2)*s) +
                pow(ma1,4)*(6.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) + pow(mrho,2)*s) +
                pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mrho,6) - 5.*pow(mrho,4)*s + 4.*pow(mrho,2)*pow(s,2) - 1.*pow(s,3) +
                    pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) + pow(mpion,2)*(2.*pow(mrho,4) - 4.*pow(mrho,2)*s + 2.*pow(s,2)))) +
              pow(eta1,2)*(pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 3.*pow(mrho,2)*s) +
                pow(mpion,2)*(pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) + pow(mrho,2)*s*(-2.*pow(mrho,2) + 2.*s) +
                    pow(mpion,2)*(3.*pow(mrho,4) - 5.*pow(mrho,2)*s)) +
                pow(ma1,2)*(-4.*pow(mpion,6) + pow(mrho,4)*s - 1.*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) +
                    pow(mpion,2)*(-2.*pow(mrho,4) + 4.*pow(mrho,2)*s + 2.*pow(s,2)))))*log(fabs(-pow(ma1,2) + tmin)))/
          (0.25*pow(Gammaa1,2)*pow(ma1,2) + 1.*pow(ma1,4) + 1.*pow(mpion,4) + 1.*pow(mpion,2)*pow(mrho,2) + 0.25*pow(mrho,4) -
            1.*pow(mpion,2)*s - 0.5*pow(mrho,2)*s + 0.25*pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s)) +
        (1.*(eta1*eta2*(pow(ma1,8)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
                pow(ma1,6)*(-1.*pow(mrho,4) + 2.*C4*pow(mrho,6) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*C4*pow(mrho,4)) +
                    0.5*delta*pow(s,2) + pow(mrho,2)*s*(2. - 1.*delta - 2.*C4*s)) +
                pow(ma1,2)*(pow(mpion,6)*(-4.*pow(mrho,2) + 8.*C4*pow(mrho,4)) +
                    pow(mrho,4)*s*((0.5 - 0.25*delta)*pow(mrho,2) + (-0.5 + 0.25*delta)*s) +
                    pow(mpion,4)*(10.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,4)*(-3. - 1.*delta - 8.*C4*s) +
                      pow(mrho,2)*s*(2. + 1.*delta - 2.*C4*s)) +
                    pow(mpion,2)*(2.*C4*pow(mrho,8) - 0.5*delta*pow(s,3) + pow(mrho,6)*(1. - 1.*delta - 2.*C4*s) +
                      pow(mrho,4)*s*(-1. + 1.*delta - 2.*C4*s) + pow(mrho,2)*pow(s,2)*(1. + 2.*C4*s))) +
                pow(mpion,2)*pow(mrho,2)*(pow(mpion,6)*(1. - 2.*C4*pow(mrho,2)) +
                    pow(mrho,4)*((-0.5 + 0.25*delta)*pow(mrho,2) + (0.5 - 0.25*delta)*s) +
                    pow(mpion,2)*(-2.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,2)*s*(-1.5 - 0.25*delta - 2.*C4*s) +
                      pow(mrho,4)*(1. + 4.*C4*s)) + pow(mpion,4)*
                    (-4.*C4*pow(mrho,4) - 1.*delta*s + pow(mrho,2)*(1. + 0.5*delta + 4.*C4*s))) +
                pow(ma1,4)*(pow(mpion,4)*(6.*pow(mrho,2) - 12.*C4*pow(mrho,4)) +
                    pow(mpion,2)*(-8.*C4*pow(mrho,6) - 1.*delta*pow(s,2) + pow(mrho,4)*(3. + 0.5*delta + 4.*C4*s) +
                      pow(mrho,2)*s*(-4. + 1.*delta + 4.*C4*s)) +
                    s*(-2.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,2)*s*(1. - 1.5*delta - 2.*C4*s) +
                      pow(mrho,4)*(-1.5 + 1.25*delta + 4.*C4*s)))) +
              pow(eta1,2)*(pow(ma1,8)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                pow(ma1,6)*(-2.*C4*pow(mrho,6) + pow(mpion,2)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) - 0.25*delta*pow(s,2) +
                    pow(mrho,4)*(1. + 1.*C4*s) + pow(mrho,2)*s*(-1. + 0.25*delta + 1.*C4*s)) +
                pow(ma1,4)*(1.*C4*pow(mrho,8) + pow(mpion,4)*(-3.*pow(mrho,2) + 6.*C4*pow(mrho,4)) - 0.25*delta*pow(s,3) +
                    pow(mrho,6)*(-0.5 - 1.*C4*s) + pow(mrho,4)*s*(1.5 - 0.5*delta - 1.*C4*s) +
                    pow(mrho,2)*pow(s,2)*(-0.5 + 0.5*delta + 1.*C4*s) +
                    pow(mpion,2)*(7.*C4*pow(mrho,6) + 0.5*delta*pow(s,2) + pow(mrho,4)*(-3. - 0.25*delta - 5.*C4*s) +
                      pow(mrho,2)*s*(2. + 0.25*delta - 2.*C4*s))) +
                pow(mpion,2)*pow(mrho,2)*(pow(mpion,6)*(-0.5 + 1.*C4*pow(mrho,2)) +
                    pow(mrho,4)*((0.5 - 0.25*delta)*pow(mrho,2) + (-0.5 + 0.25*delta)*s) +
                    pow(mpion,4)*(3.*C4*pow(mrho,4) + 0.75*delta*s + pow(mrho,2)*(-1. - 0.25*delta - 3.*C4*s)) +
                    pow(mpion,2)*(2.*C4*pow(mrho,6) - 0.5*delta*pow(s,2) + pow(mrho,4)*(-1. - 4.*C4*s) +
                      pow(mrho,2)*s*(1.5 + 0.25*delta + 2.*C4*s))) +
                pow(ma1,2)*(pow(mpion,6)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                    pow(mrho,4)*s*((-0.5 + 0.25*delta)*pow(mrho,2) + (0.5 - 0.25*delta)*s) +
                    pow(mpion,2)*(-3.*C4*pow(mrho,8) + 0.25*delta*pow(s,3) + pow(mrho,4)*s*(-1. - 0.75*delta - 1.*C4*s) +
                      pow(mrho,2)*pow(s,2)*(-0.5 + 0.5*delta - 1.*C4*s) + pow(mrho,6)*(0.5 + 0.5*delta + 5.*C4*s)) +
                    pow(mpion,4)*(-8.*C4*pow(mrho,6) - 0.25*delta*pow(s,2) + pow(mrho,2)*s*(-1. - 1.25*delta + 1.*C4*s) +
                      pow(mrho,4)*(3. + 0.5*delta + 7.*C4*s)))) +
              pow(eta2,2)*(pow(ma1,8)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                pow(mpion,6)*pow(mrho,2)*(-0.25*delta*pow(mrho,2) + 1.*C4*pow(mrho,4) + pow(mpion,2)*(-0.5 + 1.*C4*pow(mrho,2)) +
                    0.25*delta*s - 1.*C4*pow(mrho,2)*s) +
                pow(ma1,6)*(pow(mpion,2)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                    s*(-1.*C4*pow(mrho,4) - 0.25*delta*s + pow(mrho,2)*(-1. + 0.75*delta + 1.*C4*s))) +
                pow(ma1,4)*(-1.*C4*pow(mrho,8) + pow(mpion,4)*(-3.*pow(mrho,2) + 6.*C4*pow(mrho,4)) - 0.25*delta*pow(s,3) +
                    pow(mrho,4)*s*(-0.75*delta - 3.*C4*s) + pow(mrho,2)*pow(s,2)*(-0.5 + 1.*delta + 1.*C4*s) +
                    pow(mrho,6)*(0.5 + 3.*C4*s) + pow(mpion,2)*
                    (delta*(-0.25*pow(mrho,4) - 1.25*pow(mrho,2)*s + 0.5*pow(s,2)) +
                      pow(mrho,2)*(1.*C4*pow(mrho,4) + 2.*s + 1.*C4*pow(mrho,2)*s - 2.*C4*pow(s,2)))) +
                pow(ma1,2)*pow(mpion,2)*(1.*C4*pow(mrho,8) + pow(mpion,4)*(2.*pow(mrho,2) - 4.*C4*pow(mrho,4)) +
                    0.25*delta*pow(s,3) + pow(mrho,6)*(-1.5 + 0.5*delta - 3.*C4*s) + pow(mrho,2)*pow(s,2)*(-0.5 - 0.5*delta - 1.*C4*s) +
                    pow(mrho,4)*s*(2. - 0.25*delta + 3.*C4*s) +
                    pow(mpion,2)*(delta*(0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
                      pow(mrho,2)*(-2.*C4*pow(mrho,4) - 1.*s + 1.*C4*pow(mrho,2)*s + 1.*C4*pow(s,2))))))*log(fabs(-pow(ma1,2) + tmin)))/
          ((pow(ma1,2) - 1.*pow(mpion,2))*pow(mrho,2)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s)) -
        (0.5*pow(mpion,2)*(pow(eta2,2)*pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s +
                (-2. + 1.*delta)*pow(s,2)) + eta1*eta2*(pow(mpion,2)*
                  ((4. - 2.*delta)*pow(mrho,4) + (-8. + 4.*delta)*pow(mrho,2)*s + (4. - 2.*delta)*pow(s,2)) +
                pow(mrho,2)*((-1. + 0.5*delta)*pow(mrho,4) + (2. - 1.*delta)*pow(mrho,2)*s + (-1. + 0.5*delta)*pow(s,2))) +
              pow(eta1,2)*(pow(mrho,2)*((1. - 0.5*delta)*pow(mrho,4) + (-2. + 1.*delta)*pow(mrho,2)*s + (1. - 0.5*delta)*pow(s,2)) +
                pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))))*
            log(fabs(-pow(mpion,2) + tmin)))/((-1.*pow(ma1,2) + 1.*pow(mpion,2))*(pow(mrho,2) - 1.*s)) +
        (0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*
            (eta1*(pow(mpion,4)*(-1.*pow(mrho,2) + 1.*s) +
                pow(ma1,2)*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + (-0.5*pow(mrho,2) + 0.5*s)*s) +
                pow(mpion,2)*(-1.*pow(mrho,4) + 2.5*pow(mrho,2)*s - 1.5*pow(s,2)) +
                s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2))) +
              eta2*(0.5*pow(mrho,6) + pow(mpion,4)*(1.*pow(mrho,2) - 1.*s) - 1.5*pow(mrho,4)*s + 1.5*pow(mrho,2)*pow(s,2) -
                0.5*pow(s,3) + pow(mpion,2)*(1.5*pow(mrho,4) - 3.*pow(mrho,2)*s + 1.5*pow(s,2)) +
                pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*pow(mrho,2)*s - 0.5*pow(s,2) + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*
            log(fabs(-pow(mpion,2) + tmin)))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
            2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
        (0.5*pow(mpion,2)*(eta1*eta2*((1. - 0.5*delta)*pow(mrho,6) + (-4. + 2.*delta)*pow(mrho,4)*s +
                (5. - 2.5*delta)*pow(mrho,2)*pow(s,2) + (-2. + 1.*delta)*pow(s,3) +
                pow(mpion,2)*((4. - 2.*delta)*pow(mrho,4) + (-8. + 4.*delta)*pow(mrho,2)*s + (4. - 2.*delta)*pow(s,2))) +
              pow(eta2,2)*((-1. + 0.5*delta)*pow(mrho,6) + (3. - 1.5*delta)*pow(mrho,4)*s + (-3. + 1.5*delta)*pow(mrho,2)*pow(s,2) +
                (1. - 0.5*delta)*pow(s,3) + pow(mpion,2)*
                  ((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))) +
              pow(eta1,2)*(s*((1. - 0.5*delta)*pow(mrho,4) + (-2. + 1.*delta)*pow(mrho,2)*s + (1. - 0.5*delta)*pow(s,2)) +
                pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s + (-2. + 1.*delta)*pow(s,2))))*
            log(fabs(-pow(mpion,2) - pow(mrho,2) + s + tmin)))/
          ((pow(mrho,2) - 1.*s)*(-1.*pow(ma1,2) + pow(mpion,2) + pow(mrho,2) - 1.*s)) -
        0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) + tmin)) -
        (0.25000000000000006*pow(2. - 1.*delta,2)*(7.999999999999998*pow(mpion,4) - 5.999999999999998*pow(mpion,2)*s + 1.*pow(s,2))*
            log(fabs(-pow(mpion,2) + tmin)))/(pow(mrho,2) - 1.*s) +
        (1.*(-2. + 1.*delta)*((0.5 - 0.25*delta)*pow(mrho,2)*s +
              pow(mpion,2)*(4.*C4*pow(mrho,4) + 1.*delta*s + pow(mrho,2)*(-2. - 4.*C4*s)))*log(fabs(-pow(mpion,2) + tmin)))/pow(mrho,2)
          + 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + tmin)) +
        (2.*(1.*pow(2. - 1.*delta,2)*pow(mpion,4)*pow(mrho,2) + 0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,4)*s +
              pow(mpion,2)*(C4*(4. - 2.*delta)*pow(mrho,6) + (-1. + 0.5*delta)*delta*pow(s,2) +
                pow(mrho,2)*s*(-1. + 3.*delta - 1.25*pow(delta,2) + 4.*C4*s - 2.*C4*delta*s) +
                pow(mrho,4)*(-2. + 1.*delta - 8.*C4*s + 4.*C4*delta*s)))*log(fabs(-pow(mpion,2) - pow(mrho,2) + s + tmin)))/
          (pow(mrho,4) - 1.*pow(mrho,2)*s) - (0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*
            (eta2*pow(mpion,2)*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + pow(ma1,2)*(-1.*pow(mrho,2) + 1.*s)) +
              eta1*(pow(ma1,2)*(-0.5*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s) +
                pow(mpion,2)*(0.5*pow(mrho,4) - 0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*
            log(fabs(-pow(mpion,2) - pow(mrho,2) + s + tmin)))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) -
        (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8) + 0.5*pow(mpion,8) +
                0.5*pow(ma1,4)*pow(mpion,2)*pow(mrho,2) - 0.5*pow(mpion,6)*pow(mrho,2) +
                pow(Gammaa1,2)*pow(ma1,2)*(pow(mpion,2)*(1.*pow(mpion,2) + 1.5*pow(mrho,2) - 2.*s) +
                    pow(ma1,2)*(-1.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s)) +
                pow(ma1,6)*(1.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) +
                pow(ma1,2)*pow(mpion,4)*(-1.*pow(mpion,2) - 0.5*pow(mrho,2) + 1.*s)) +
              eta1*(pow(ma1,6)*(1.*pow(mpion,2) + 0.5*s) +
                pow(ma1,2)*(3.*pow(mpion,6) + 1.*pow(mpion,2)*pow(mrho,4) - 0.5*pow(mpion,4)*s) +
                pow(ma1,4)*(-3.*pow(mpion,4) + pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
                pow(mpion,4)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 0.5*s) + 0.5*pow(mrho,2)*s) +
                pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) + pow(ma1,2)*(1.*pow(mpion,2) + 0.5*s) - 0.5*pow(mrho,2)*s +
                    pow(mpion,2)*(-1.*pow(mrho,2) + 1.5*s))))*
            log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*tmin - 2.*pow(mrho,2)*tmin + 2.*s*tmin + pow(tmin,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*tmin))))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) +
        (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(-1.*pow(mpion,8) + 1.*pow(mpion,4)*pow(mrho,4) +
                pow(ma1,6)*(1.*pow(mpion,2) - 0.5*s) - 0.5*pow(mpion,6)*s - 4.*pow(mpion,4)*pow(mrho,2)*s -
                1.5*pow(mpion,2)*pow(mrho,4)*s + 3.5*pow(mpion,4)*pow(s,2) + 4.*pow(mpion,2)*pow(mrho,2)*pow(s,2) +
                0.5*pow(mrho,4)*pow(s,2) - 2.5*pow(mpion,2)*pow(s,3) - 1.*pow(mrho,2)*pow(s,3) + 0.5*pow(s,4) +
                pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) + pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) - 0.5*pow(mpion,2)*s +
                    0.5*pow(s,2)) + pow(ma1,4)*(-3.*pow(mpion,4) + (1.*pow(mrho,2) - 0.5*s)*s +
                    pow(mpion,2)*(-2.*pow(mrho,2) + 2.5*s)) +
                pow(ma1,2)*(3.*pow(mpion,6) + pow(mpion,4)*(2.*pow(mrho,2) - 1.5*s) - 0.5*pow(mrho,4)*s + 0.5*pow(s,3) +
                    pow(mpion,2)*(1.*pow(mrho,4) - 1.*pow(mrho,2)*s - 1.*pow(s,2)))) +
              eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8) + pow(ma1,6)*(1.*pow(mpion,2) + 1.*pow(mrho,2) - 0.5*s) +
                pow(ma1,4)*(-0.5*pow(mrho,4) + (-0.5*pow(mpion,2) + 0.5*s)*s) +
                pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) -
                    1.*pow(mrho,2)*s + 0.5*pow(s,2) + pow(ma1,2)*(-1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.5*s)) +
                pow(mpion,2)*(0.5*pow(mpion,6) + 0.5*pow(mpion,4)*s +
                    pow(mpion,2)*(-0.5*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.5*pow(s,2)) +
                    s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2))) +
                pow(ma1,2)*(-1.*pow(mpion,6) + pow(mpion,4)*(-1.*pow(mrho,2) + 0.5*s) +
                    pow(mpion,2)*(-1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
                    s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s + 0.5*pow(s,2)))))*
            log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*tmin - 2.*pow(mrho,2)*tmin + 2.*s*tmin + pow(tmin,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*tmin))))/
          (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
            2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
        (0.0625*pow(eta1 - 1.*eta2,2)*(pow(eta2,2)*(-1.*pow(ma1,10) + pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
                pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
                pow(ma1,4)*(10.*pow(mpion,6) - 2.5*pow(mrho,6) + pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) + 6.*pow(mrho,4)*s -
                    4.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3)) +
                pow(ma1,6)*(-10.*pow(mpion,4) + (1.*pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-10.*pow(mrho,2) + 8.*s)) +
                pow(mpion,4)*(1.*pow(mpion,6) + 0.5*pow(mrho,6) + pow(mpion,4)*(2.5*pow(mrho,2) - 0.5*s) - 1.*pow(mrho,4)*s +
                    0.5*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(2.*pow(mrho,4) - 2.*pow(mrho,2)*s)) +
                pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) - 4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
                    1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
                    pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
                    pow(mpion,2)*(-3.*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.*pow(s,2)) +
                    pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4) + pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
                      3.*pow(s,2))) + pow(ma1,2)*(-5.*pow(mpion,8) + 1.*pow(mrho,8) - 3.5*pow(mrho,6)*s +
                    4.5*pow(mrho,4)*pow(s,2) - 2.5*pow(mrho,2)*pow(s,3) + 0.5*pow(s,4) + pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
                    pow(mpion,4)*(-2.*pow(mrho,4) + 1.*pow(mrho,2)*s + 1.*pow(s,2)) +
                    pow(mpion,2)*(3.*pow(mrho,6) - 8.*pow(mrho,4)*s + 7.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)))) +
              pow(eta1,2)*(-1.*pow(ma1,10) + pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
                pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
                pow(ma1,6)*(-10.*pow(mpion,4) - 2.*pow(mrho,4) + 5.*pow(mrho,2)*s - 1.*pow(s,2) +
                    pow(mpion,2)*(-10.*pow(mrho,2) + 8.*s)) +
                pow(ma1,4)*(10.*pow(mpion,6) + 0.5*pow(mrho,6) + pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) - 3.*pow(mrho,4)*s +
                    1.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s)) +
                pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) - 4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
                    1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
                    pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
                    pow(mpion,2)*(-3.*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.*pow(s,2)) +
                    pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4) + pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
                      3.*pow(s,2))) + pow(mpion,2)*(1.*pow(mpion,8) + pow(mpion,6)*(2.5*pow(mrho,2) - 0.5*s) +
                    pow(mpion,4)*(4.*pow(mrho,4) - 6.*pow(mrho,2)*s) +
                    pow(mrho,2)*s*(-1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
                    pow(mpion,2)*(1.5*pow(mrho,6) - 6.*pow(mrho,4)*s + 4.5*pow(mrho,2)*pow(s,2))) +
                pow(ma1,2)*(-5.*pow(mpion,8) + pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
                    pow(mpion,4)*(-8.*pow(mrho,4) + 13.*pow(mrho,2)*s + 1.*pow(s,2)) +
                    pow(mpion,2)*(-1.*pow(mrho,6) + 6.*pow(mrho,4)*s - 3.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)) +
                    s*(0.5*pow(mrho,6) - 0.5*pow(mrho,4)*s - 0.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3)))) +
              eta1*eta2*(2.*pow(ma1,10) + pow(Gammaa1,4)*pow(ma1,4)*(-2.*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s) +
                pow(ma1,8)*(-10.*pow(mpion,2) - 5.*pow(mrho,2) + 5.*s) +
                pow(ma1,6)*(20.*pow(mpion,4) + 6.*pow(mrho,4) + pow(mpion,2)*(20.*pow(mrho,2) - 16.*s) - 8.*pow(mrho,2)*s +
                    2.*pow(s,2)) + pow(ma1,4)*(-20.*pow(mpion,6) - 4.*pow(mrho,6) + 6.*pow(mrho,4)*s - 2.*pow(s,3) +
                    pow(mpion,4)*(-30.*pow(mrho,2) + 18.*s) + pow(mpion,2)*(-18.*pow(mrho,4) + 18.*pow(mrho,2)*s)) +
                pow(mpion,2)*(-2.*pow(mpion,8) - 1.*pow(mrho,8) + 3.*pow(mrho,6)*s - 3.*pow(mrho,4)*pow(s,2) +
                    1.*pow(mrho,2)*pow(s,3) + pow(mpion,6)*(-5.*pow(mrho,2) + 1.*s) +
                    pow(mpion,4)*(-10.*pow(mrho,4) + 10.*pow(mrho,2)*s) +
                    pow(mpion,2)*(-6.*pow(mrho,6) + 12.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2))) +
                pow(ma1,2)*(10.*pow(mpion,8) + 1.*pow(mrho,8) + pow(mpion,6)*(20.*pow(mrho,2) - 8.*s) - 2.*pow(mrho,6)*s +
                    2.*pow(mrho,2)*pow(s,3) - 1.*pow(s,4) + pow(mpion,4)*(22.*pow(mrho,4) - 20.*pow(mrho,2)*s - 2.*pow(s,2)) +
                    pow(mpion,2)*(8.*pow(mrho,6) - 12.*pow(mrho,4)*s + 4.*pow(s,3))) +
                pow(Gammaa1,2)*pow(ma1,2)*(-8.*pow(ma1,6) + 8.*pow(mpion,6) + 1.*pow(mrho,6) +
                    pow(mpion,4)*(12.*pow(mrho,2) - 12.*s) + pow(ma1,4)*(24.*pow(mpion,2) + 12.*pow(mrho,2) - 12.*s) -
                    3.*pow(mrho,4)*s + 3.*pow(mrho,2)*pow(s,2) - 1.*pow(s,3) +
                    pow(mpion,2)*(6.*pow(mrho,4) - 12.*pow(mrho,2)*s + 6.*pow(s,2)) +
                    pow(ma1,2)*(-24.*pow(mpion,4) - 6.*pow(mrho,4) + 12.*pow(mrho,2)*s - 6.*pow(s,2) +
                      pow(mpion,2)*(-24.*pow(mrho,2) + 24.*s)))))*
            log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*tmin - 2.*pow(mrho,2)*tmin + 2.*s*tmin + pow(tmin,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*tmin))))/
          (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
            4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) -
        0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(4.*pow(ma1,6) +
              pow(Gammaa1,2)*pow(ma1,2)*(-4.*pow(ma1,2) + 4.*pow(mpion,2) - 2.*s) + pow(ma1,4)*(-12.*pow(mpion,2) + 6.*s) +
              pow(mpion,2)*(-4.*pow(mpion,4) + 4.*pow(mrho,4) + 2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s) +
              pow(ma1,2)*(12.*pow(mpion,4) - 2.*pow(mrho,4) - 8.*pow(mpion,2)*s - 4.*pow(mrho,2)*s + 4.*pow(s,2))) +
            pow(eta1,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) + 3.*pow(mpion,4)*pow(mrho,2) +
              pow(ma1,4)*(6.*pow(mpion,2) + 3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s - 1.*pow(mpion,2)*pow(mrho,2)*s -
              1.*pow(mrho,4)*s + pow(mrho,2)*pow(s,2) +
              pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s) +
              pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) + 4.*pow(mrho,2)*s - 2.*pow(s,2) +
                  pow(mpion,2)*(-6.*pow(mrho,2) + 4.*s))) +
            pow(eta2,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) - 3.*pow(mpion,4)*pow(mrho,2) +
              pow(ma1,4)*(6.*pow(mpion,2) - 3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s + pow(mpion,2)*pow(mrho,2)*s +
              pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) - 2.*pow(mpion,2) + pow(mrho,2) + s) +
              pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 2.*pow(s,2) +
                  pow(mpion,2)*(6.*pow(mrho,2) + 4.*s))))*
          log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + 1.*pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*tmin - 2.*pow(mrho,2)*tmin + 2.*s*tmin + pow(tmin,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*tmin))) +
        (0.5*(1.*eta1 - 1.*eta2)*(eta2*(pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(0.5 - 1.*C4*pow(mrho,2)) +
                pow(ma1,4)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
                pow(mpion,2)*pow(mrho,2)*(pow(mpion,2)*(-0.5 + 1.*C4*pow(mrho,2)) + (0.25 - 0.125*delta)*(pow(mrho,2) + s)) +
                pow(ma1,2)*(1.*C4*pow(mrho,6) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) - 0.25*delta*pow(s,2) +
                    pow(mrho,4)*(-0.75 + 0.125*delta - 2.*C4*s) + pow(mrho,2)*s*(0.25 + 0.375*delta + 1.*C4*s))) +
              eta1*(pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(-0.5 + 1.*C4*pow(mrho,2)) +
                pow(ma1,4)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
                pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6) + pow(mpion,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) +
                    0.25*delta*pow(s,2) - 1.*C4*pow(mrho,2)*pow(s,2)) +
                pow(mrho,2)*(pow(mpion,4)*(0.5 - 1.*C4*pow(mrho,2)) + s*((-0.25 + 0.125*delta)*pow(mrho,2) + (0.25 - 0.125*delta)*s) +
                    pow(mpion,2)*(-2.*C4*pow(mrho,4) + (-0.5 - 0.25*delta)*s + pow(mrho,2)*(1. + 2.*C4*s)))))*
            log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + 1.*pow(mrho,4) - 4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) - 4.*pow(mpion,2)*tmin - 2.*pow(mrho,2)*tmin + 2.*s*tmin + pow(tmin,2) + pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s + 2.*tmin))))/pow(mrho,2)))/
    (16.*M_PI*s*(-4*pow(mpion,2) + s));
    return xs;
}


double PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_pi0_rho(
    const double s, const double m_rho) {
  using std::atan;
  using std::pow;
  using std::sqrt;
  const double mpion = m_pion_, mrho = m_rho;
  auto t_mandelstam = get_t_range(sqrt(s), m_pion_, m_pion_, m_rho, 0.0);

  const double t1 = t_mandelstam[1];
  const double t2 = t_mandelstam[0];
  const double xs = to_mb*(-(pow(Const,2)*pow(ghat,4)*((0.03125*pow(eta1 - 1.*eta2,2)*
                                     (eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
                                          pow(ma1,2)*pow(mpion,2)*(8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
                                          pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
                                       pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
                                          pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                                          pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                                             2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(6.*pow(mrho,2) + 2.*s)))
                                         + pow(eta1,2)*(1.*pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                                          pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                                             2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
                                             pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
                                          pow(mpion,2)*(1.*pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                                             pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/(1.*pow(ma1,2) - 1.*t2) +
                                  (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*t2) -
                                  (0.25*pow(-2. + delta,2)*pow(mpion,2)*t2)/pow(mrho,2) -
                                  0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,2) + pow(mrho,2) - 2.*s) + eta1*(2.*pow(mpion,2) + s))*t2 +
                                  (0.5*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*(4.*pow(mpion,4)*pow(mrho,2) + 1.*pow(mrho,6) - 3.5*pow(mrho,4)*s + 0.5*pow(s,3) +
                                       pow(mpion,2)*(10.*pow(mrho,4) - 2.*pow(s,2)))*t2)/(pow(mrho,6)*pow(pow(mrho,2) - 1.*s,2)) -
                                  (0.25*(eta1 - 1.*eta2)*(1.*pow(mrho,2) - 0.5*delta*s)*
                                     (eta2*(-2.*pow(ma1,4) - 6.*pow(mpion,4) + 1.5*pow(mrho,4) + pow(ma1,2)*(6.*pow(mpion,2) - 2.*pow(mrho,2) - 2.*s) -
                                          1.*pow(mrho,2)*s - 0.5*pow(s,2) + pow(mpion,2)*(2.*pow(mrho,2) + 2.*s)) +
                                       eta1*(2.*pow(ma1,4) + 6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                                          1.*pow(s,2) + pow(ma1,2)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s)))*t2)/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
                                  0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-6.*pow(ma1,4) - 12.*pow(mpion,4) + 2.*pow(mrho,4) +
                                        pow(ma1,2)*(16.*pow(mpion,2) - 8.*s) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
                                     pow(eta1,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                                        2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
                                     pow(eta2,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                                        2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) + 4.*s)))*t2 -
                                  (1.*(pow(mpion,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (1.5 - 2.*delta + 0.625*pow(delta,2))*pow(mrho,2)*s +
                                          (0.25 - 0.125*delta)*delta*pow(s,2) + pow(mrho,4)*(2.5 - 2.25*delta + 0.5*pow(delta,2) + 2.*C4*s - 1.*C4*delta*s)) +
                                       pow(mrho,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (0.75 - 0.375*delta)*delta*pow(s,2) +
                                          pow(mrho,4)*(0.5 - 0.25*delta + 6.*C4*s - 3.*C4*delta*s) +
                                          pow(mrho,2)*s*(-0.5 - 0.5*delta + 0.375*pow(delta,2) - 4.*C4*s + 2.*C4*delta*s)))*t2)/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
                                  (0.25*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*(eta1*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(ma1,2)*(1. - 2.*C4*pow(mrho,2)) +
                                             pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) - 2.*C4*pow(s,2)) +
                                          eta2*(-1.5*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) +
                                             pow(ma1,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s - 4.*C4*pow(mrho,2)*s + 2.*C4*pow(s,2))) +
                                       delta*(eta2*(-1.*pow(ma1,4) - 3.*pow(mpion,4) + 1.*pow(mrho,4) + pow(ma1,2)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                                             0.25*pow(mrho,2)*s - 0.75*pow(s,2) + pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)) +
                                          eta1*(1.*pow(ma1,4) + 3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s +
                                             1.*pow(s,2) + pow(ma1,2)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s))))*t2)/pow(mrho,2) +
                                  (0.5*(pow(delta,2)*(1.*pow(mpion,4)*pow(mrho,2) + 0.25*pow(mrho,6) - 0.75*pow(mrho,4)*s + 0.125*pow(mrho,2)*pow(s,2) +
                                          0.25*pow(s,3) + pow(mpion,2)*(2.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2))) +
                                       pow(mrho,6)*(1.5 + C4*(-6.*pow(mrho,2) + 6.*s) + pow(C4,2)*(4.*pow(mrho,4) - 8.*pow(mrho,2)*s + 4.*pow(s,2))) +
                                       delta*pow(mrho,2)*(4.*C4*pow(mrho,6) - 0.5*pow(s,2) + pow(mrho,4)*(-1.5 - 3.*C4*s) + pow(mrho,2)*s*(0.5 - 1.*C4*s) +
                                          pow(mpion,2)*(6.*C4*pow(mrho,4) + 0.5*s + pow(mrho,2)*(-2.5 - 2.*C4*s))))*t2)/pow(mrho,6) -
                                  (3.*(1.*pow(mrho,2) - 0.5*delta*s)*(delta*(0.666667*pow(mpion,4)*pow(mrho,2) + 0.166667*pow(mrho,6) - 0.541667*pow(mrho,4)*s -
                                          0.0833333*pow(mrho,2)*pow(s,2) + 0.125*pow(s,3) +
                                          pow(mpion,2)*(1.66667*pow(mrho,4) + 0.0833333*pow(mrho,2)*s - 0.416667*pow(s,2))) +
                                       pow(mrho,2)*(1.*C4*pow(mrho,6) - 0.0833333*pow(s,2) + pow(mrho,4)*(-0.416667 - 1.33333*C4*s) +
                                          pow(mrho,2)*s*(0.5 + 0.333333*C4*s) + pow(mpion,2)*(2.*C4*pow(mrho,4) + 0.166667*s + pow(mrho,2)*(-0.833333 - 0.666667*C4*s))
                                          ))*t2)/(pow(mrho,8) - 1.*pow(mrho,6)*s) - 1.*C4*pow(t2,2) - 1.*C4*delta*pow(t2,2) +
                                  0.0625*(-2. + delta)*(eta1 - 1.*eta2)*eta2*pow(t2,2) - (0.5*pow(delta,2)*pow(mpion,2)*pow(t2,2))/pow(mrho,4) +
                                  (0.25*pow(t2,2))/pow(mrho,2) + (0.5*delta*pow(t2,2))/pow(mrho,2) - (0.25*pow(delta,2)*pow(t2,2))/pow(mrho,2) -
                                  (0.25*delta*s*pow(t2,2))/pow(mrho,4) + (0.25*pow(delta,2)*s*pow(t2,2))/pow(mrho,4) +
                                  (0.5*C4*delta*s*pow(t2,2))/pow(mrho,2) + (0.0625*pow(delta,2)*pow(s,2)*pow(t2,2))/pow(mrho,6) -
                                  (1.*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s)*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*pow(t2,2))/
                                   (pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) + (0.375*(eta1 - 1.*eta2)*
                                     (eta1*(-0.6666666666666666*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s) +
                                       eta2*(0.6666666666666666*pow(ma1,2) - 2.*pow(mpion,2) + 0.6666666666666666*pow(mrho,2) + 0.6666666666666666*s))*
                                     (1.*pow(mrho,2) - 0.5*delta*s)*pow(t2,2))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
                                  0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                                     eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(t2,2) +
                                  (3.*(1.*pow(mrho,2) - 0.5*delta*s)*(1.*C4*pow(mrho,6) + 0.0833335*pow(mrho,2)*s + pow(mrho,4)*(-0.416667 - 0.333334*C4*s) +
                                       delta*(0.666665*pow(mpion,2)*pow(mrho,2) + 0.333334*pow(mrho,4) - 0.291667*pow(mrho,2)*s - 0.0416667*pow(s,2)))*pow(t2,2))
                                    /(pow(mrho,8) - 1.*pow(mrho,6)*s) + (0.125*(1.*eta1 - 1.*eta2)*
                                     (pow(mrho,2)*(eta1*(1. - 2.*C4*pow(mrho,2)) + eta2*(-1. + 2.*C4*pow(mrho,2))) +
                                       delta*(eta2*(-1.*pow(ma1,2) + 3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                                          eta1*(1.*pow(ma1,2) - 3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s)))*pow(t2,2))/pow(mrho,2) +
                                  0.0104167*pow(eta1 - 1.*eta2,4)*pow(t2,3) + (0.166667*pow(delta,2)*pow(t2,3))/pow(mrho,4) +
                                  (0.0833333*delta*pow(1.*eta1 - 1.*eta2,2)*pow(t2,3))/pow(mrho,2) +
                                  (0.666667*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*pow(t2,3))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) -
                                  (0.166667*pow(1.*eta1 - 1.*eta2,2)*(1.*pow(mrho,2) - 0.5*delta*s)*pow(t2,3))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
                                  (0.333334*delta*(-2.*pow(mrho,2) + 1.*delta*s)*pow(t2,3))/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
                                  (0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
                                          pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) + pow(ma1,2)*pow(mpion,2)*
                                           (8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
                                          pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
                                       pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
                                          pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                                          pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                                             2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(6.*pow(mrho,2) + 2.*s)))
                                         + pow(eta1,2)*(1.*pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                                          pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                                             2.*pow(s,2)) + pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s +
                                             pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
                                          pow(mpion,2)*(1.*pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                                             pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/(1.*pow(ma1,2) - 1.*t1) -
                                  (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*t1) +
                                  (0.25*pow(-2. + delta,2)*pow(mpion,2)*t1)/pow(mrho,2) +
                                  0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,2) + pow(mrho,2) - 2.*s) + eta1*(2.*pow(mpion,2) + s))*t1 -
                                  (0.5*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*(4.*pow(mpion,4)*pow(mrho,2) + 1.*pow(mrho,6) - 3.5*pow(mrho,4)*s + 0.5*pow(s,3) +
                                       pow(mpion,2)*(10.*pow(mrho,4) - 2.*pow(s,2)))*t1)/(pow(mrho,6)*pow(pow(mrho,2) - 1.*s,2)) +
                                  (0.25*(eta1 - 1.*eta2)*(1.*pow(mrho,2) - 0.5*delta*s)*
                                     (eta2*(-2.*pow(ma1,4) - 6.*pow(mpion,4) + 1.5*pow(mrho,4) + pow(ma1,2)*(6.*pow(mpion,2) - 2.*pow(mrho,2) - 2.*s) -
                                          1.*pow(mrho,2)*s - 0.5*pow(s,2) + pow(mpion,2)*(2.*pow(mrho,2) + 2.*s)) +
                                       eta1*(2.*pow(ma1,4) + 6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                                          1.*pow(s,2) + pow(ma1,2)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s)))*t1)/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
                                  0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-6.*pow(ma1,4) - 12.*pow(mpion,4) + 2.*pow(mrho,4) +
                                        pow(ma1,2)*(16.*pow(mpion,2) - 8.*s) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
                                     pow(eta1,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
                                        2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
                                     pow(eta2,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
                                        2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) + 4.*s)))*t1 +
                                  (1.*(pow(mpion,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (1.5 - 2.*delta + 0.625*pow(delta,2))*pow(mrho,2)*s +
                                          (0.25 - 0.125*delta)*delta*pow(s,2) + pow(mrho,4)*(2.5 - 2.25*delta + 0.5*pow(delta,2) + 2.*C4*s - 1.*C4*delta*s)) +
                                       pow(mrho,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (0.75 - 0.375*delta)*delta*pow(s,2) +
                                          pow(mrho,4)*(0.5 - 0.25*delta + 6.*C4*s - 3.*C4*delta*s) +
                                          pow(mrho,2)*s*(-0.5 - 0.5*delta + 0.375*pow(delta,2) - 4.*C4*s + 2.*C4*delta*s)))*t1)/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
                                  (0.25*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*(eta1*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(ma1,2)*(1. - 2.*C4*pow(mrho,2)) +
                                             pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) - 2.*C4*pow(s,2)) +
                                          eta2*(-1.5*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) +
                                             pow(ma1,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s - 4.*C4*pow(mrho,2)*s + 2.*C4*pow(s,2))) +
                                       delta*(eta2*(-1.*pow(ma1,4) - 3.*pow(mpion,4) + 1.*pow(mrho,4) + pow(ma1,2)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                                             0.25*pow(mrho,2)*s - 0.75*pow(s,2) + pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)) +
                                          eta1*(1.*pow(ma1,4) + 3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s +
                                             1.*pow(s,2) + pow(ma1,2)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s))))*t1)/pow(mrho,2) -
                                  (0.5*(pow(delta,2)*(1.*pow(mpion,4)*pow(mrho,2) + 0.25*pow(mrho,6) - 0.75*pow(mrho,4)*s + 0.125*pow(mrho,2)*pow(s,2) +
                                          0.25*pow(s,3) + pow(mpion,2)*(2.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2))) +
                                       pow(mrho,6)*(1.5 + C4*(-6.*pow(mrho,2) + 6.*s) + pow(C4,2)*(4.*pow(mrho,4) - 8.*pow(mrho,2)*s + 4.*pow(s,2))) +
                                       delta*pow(mrho,2)*(4.*C4*pow(mrho,6) - 0.5*pow(s,2) + pow(mrho,4)*(-1.5 - 3.*C4*s) + pow(mrho,2)*s*(0.5 - 1.*C4*s) +
                                          pow(mpion,2)*(6.*C4*pow(mrho,4) + 0.5*s + pow(mrho,2)*(-2.5 - 2.*C4*s))))*t1)/pow(mrho,6) +
                                  (3.*(1.*pow(mrho,2) - 0.5*delta*s)*(delta*(0.666667*pow(mpion,4)*pow(mrho,2) + 0.166667*pow(mrho,6) - 0.541667*pow(mrho,4)*s -
                                          0.0833333*pow(mrho,2)*pow(s,2) + 0.125*pow(s,3) +
                                          pow(mpion,2)*(1.66667*pow(mrho,4) + 0.0833333*pow(mrho,2)*s - 0.416667*pow(s,2))) +
                                       pow(mrho,2)*(1.*C4*pow(mrho,6) - 0.0833333*pow(s,2) + pow(mrho,4)*(-0.416667 - 1.33333*C4*s) +
                                          pow(mrho,2)*s*(0.5 + 0.333333*C4*s) + pow(mpion,2)*(2.*C4*pow(mrho,4) + 0.166667*s + pow(mrho,2)*(-0.833333 - 0.666667*C4*s))
                                          ))*t1)/(pow(mrho,8) - 1.*pow(mrho,6)*s) + 1.*C4*pow(t1,2) + 1.*C4*delta*pow(t1,2) -
                                  0.0625*(-2. + delta)*(eta1 - 1.*eta2)*eta2*pow(t1,2) + (0.5*pow(delta,2)*pow(mpion,2)*pow(t1,2))/pow(mrho,4) -
                                  (0.25*pow(t1,2))/pow(mrho,2) - (0.5*delta*pow(t1,2))/pow(mrho,2) + (0.25*pow(delta,2)*pow(t1,2))/pow(mrho,2) +
                                  (0.25*delta*s*pow(t1,2))/pow(mrho,4) - (0.25*pow(delta,2)*s*pow(t1,2))/pow(mrho,4) -
                                  (0.5*C4*delta*s*pow(t1,2))/pow(mrho,2) - (0.0625*pow(delta,2)*pow(s,2)*pow(t1,2))/pow(mrho,6) +
                                  (1.*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s)*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*pow(t1,2))/
                                   (pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) - (0.375*(eta1 - 1.*eta2)*
                                     (eta1*(-0.6666666666666666*pow(ma1,2) + 2.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s) +
                                       eta2*(0.6666666666666666*pow(ma1,2) - 2.*pow(mpion,2) + 0.6666666666666666*pow(mrho,2) + 0.6666666666666666*s))*
                                     (1.*pow(mrho,2) - 0.5*delta*s)*pow(t1,2))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
                                  0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                                     eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(t1,2) -
                                  (3.*(1.*pow(mrho,2) - 0.5*delta*s)*(1.*C4*pow(mrho,6) + 0.0833335*pow(mrho,2)*s + pow(mrho,4)*(-0.416667 - 0.333334*C4*s) +
                                       delta*(0.666665*pow(mpion,2)*pow(mrho,2) + 0.333334*pow(mrho,4) - 0.291667*pow(mrho,2)*s - 0.0416667*pow(s,2)))*pow(t1,2))
                                    /(pow(mrho,8) - 1.*pow(mrho,6)*s) - (0.125*(1.*eta1 - 1.*eta2)*
                                     (pow(mrho,2)*(eta1*(1. - 2.*C4*pow(mrho,2)) + eta2*(-1. + 2.*C4*pow(mrho,2))) +
                                       delta*(eta2*(-1.*pow(ma1,2) + 3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                                          eta1*(1.*pow(ma1,2) - 3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s)))*pow(t1,2))/pow(mrho,2) -
                                  0.0104167*pow(eta1 - 1.*eta2,4)*pow(t1,3) - (0.166667*pow(delta,2)*pow(t1,3))/pow(mrho,4) -
                                  (0.0833333*delta*pow(1.*eta1 - 1.*eta2,2)*pow(t1,3))/pow(mrho,2) -
                                  (0.666667*pow(1.*pow(mrho,2) - 0.5*delta*s,2)*pow(t1,3))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) +
                                  (0.166667*pow(1.*eta1 - 1.*eta2,2)*(1.*pow(mrho,2) - 0.5*delta*s)*pow(t1,3))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
                                  (0.333334*delta*(-2.*pow(mrho,2) + 1.*delta*s)*pow(t1,3))/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
                                  0.0625*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
                                        pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) - 2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
                                        pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
                                     pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) + pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s +
                                        pow(mpion,4)*(-3.*pow(mrho,2) + s) + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
                                        pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2))) +
                                     pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
                                        pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
                                        pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2))))*
                                   log(fabs(-pow(ma1,2) + t2)) - (0.25*(-2. + delta)*(eta1 - 1.*eta2)*
                                     (eta2*(-0.5*pow(ma1,6) - 0.5*pow(mpion,6) + 0.5*pow(mpion,4)*pow(mrho,2) +
                                          pow(ma1,4)*(0.5*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) +
                                          pow(ma1,2)*pow(mpion,2)*(0.5*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s)) +
                                       eta1*(pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
                                          pow(mpion,2)*(1.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
                                          pow(ma1,2)*(-2.*pow(mpion,4) - 0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*log(fabs(-pow(ma1,2) + t2)))/
                                   (pow(ma1,2) - 1.*pow(mpion,2)) - (0.125*(eta1 - 1.*eta2)*(1.*pow(mrho,2) - 0.5*delta*s)*
                                     (eta1*(4.*pow(ma1,6) - 4.*pow(mpion,6) + pow(mrho,4)*s + 4.*pow(mpion,2)*pow(s,2) - 1.*pow(s,3) +
                                          pow(mpion,4)*(-10.*pow(mrho,2) + 2.*s) + pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
                                          pow(ma1,2)*(12.*pow(mpion,4) + 2.*pow(mrho,4) + pow(mpion,2)*(16.*pow(mrho,2) - 8.*s) - 8.*pow(mrho,2)*s + 2.*pow(s,2))
                                          ) + eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 4.*pow(mrho,2) - 4.*s) +
                                          pow(mpion,2)*(4.*pow(mpion,4) - 1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
                                          pow(ma1,2)*(-12.*pow(mpion,4) + 3.*pow(mrho,4) - 2.*pow(mrho,2)*s - 1.*pow(s,2) + pow(mpion,2)*(4.*pow(mrho,2) + 4.*s))
                                          ))*log(fabs(-pow(ma1,2) + t2)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
                                  (0.25*(1.*eta1 - 1.*eta2)*(delta*(eta1*(1.*pow(ma1,6) - 1.*pow(mpion,6) + pow(mpion,4)*(-2.5*pow(mrho,2) + 0.5*s) +
                                             pow(mpion,2)*s*(-0.5*pow(mrho,2) + 1.*s) + pow(ma1,4)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s) +
                                             s*(0.5*pow(mrho,4) - 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
                                             pow(ma1,2)*(3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s +
                                                1.*pow(s,2))) + eta2*(-1.*pow(ma1,6) + pow(ma1,4)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                                             pow(mpion,2)*(1.*pow(mpion,4) - 0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
                                             pow(ma1,2)*(-3.*pow(mpion,4) + 1.*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2) +
                                                pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)))) +
                                       pow(mrho,2)*(eta2*(pow(ma1,4)*(-1. + 2.*C4*pow(mrho,2)) +
                                             pow(mpion,2)*(0.5*pow(mrho,2) + pow(mpion,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s) +
                                             pow(ma1,2)*(2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) + pow(mrho,2)*(-1.5 - 4.*C4*s) +
                                                s*(0.5 + 2.*C4*s))) + eta1*(pow(ma1,4)*(1. - 2.*C4*pow(mrho,2)) + pow(mpion,4)*(1. - 2.*C4*pow(mrho,2)) +
                                             (-0.5*pow(mrho,2) + 0.5*s)*s + pow(ma1,2)*
                                              (-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) - 2.*C4*pow(s,2)) +
                                             pow(mpion,2)*(-4.*C4*pow(mrho,4) - 1.*s + pow(mrho,2)*(2. + 4.*C4*s)))))*log(fabs(-pow(ma1,2) + t2)))/pow(mrho,2) +
                                  0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) + t2)) +
                                  (0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*(eta2*pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s) +
                                       eta1*(-0.5*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s))*log(fabs(-pow(mpion,2) + t2)))/
                                   (-1.*pow(ma1,2) + 1.*pow(mpion,2)) - (2.*(-0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,4)*s +
                                       pow(mpion,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (0.5 - 0.25*delta)*delta*pow(s,2) +
                                          pow(mrho,4)*(1. - 0.5*delta + 4.*C4*s - 2.*C4*delta*s) +
                                          pow(mrho,2)*s*(1. - 2.*delta + 0.75*pow(delta,2) - 2.*C4*s + 1.*C4*delta*s)))*log(fabs(-pow(mpion,2) + t2)))/
                                   (pow(mrho,4) - 1.*pow(mrho,2)*s) - 0.0625*pow(eta1 - 1.*eta2,2)*
                                   (eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
                                        pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) - 2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
                                        pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
                                     pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) + pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s +
                                        pow(mpion,4)*(-3.*pow(mrho,2) + s) + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
                                        pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2))) +
                                     pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
                                        pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
                                        pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2))))*
                                   log(fabs(-pow(ma1,2) + t1)) + (0.25*(-2. + delta)*(eta1 - 1.*eta2)*
                                     (eta2*(-0.5*pow(ma1,6) - 0.5*pow(mpion,6) + 0.5*pow(mpion,4)*pow(mrho,2) +
                                          pow(ma1,4)*(0.5*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) +
                                          pow(ma1,2)*pow(mpion,2)*(0.5*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s)) +
                                       eta1*(pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
                                          pow(mpion,2)*(1.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
                                          pow(ma1,2)*(-2.*pow(mpion,4) - 0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*log(fabs(-pow(ma1,2) + t1)))/
                                   (pow(ma1,2) - 1.*pow(mpion,2)) + (0.125*(eta1 - 1.*eta2)*(1.*pow(mrho,2) - 0.5*delta*s)*
                                     (eta1*(4.*pow(ma1,6) - 4.*pow(mpion,6) + pow(mrho,4)*s + 4.*pow(mpion,2)*pow(s,2) - 1.*pow(s,3) +
                                          pow(mpion,4)*(-10.*pow(mrho,2) + 2.*s) + pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
                                          pow(ma1,2)*(12.*pow(mpion,4) + 2.*pow(mrho,4) + pow(mpion,2)*(16.*pow(mrho,2) - 8.*s) - 8.*pow(mrho,2)*s + 2.*pow(s,2))
                                          ) + eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 4.*pow(mrho,2) - 4.*s) +
                                          pow(mpion,2)*(4.*pow(mpion,4) - 1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
                                          pow(ma1,2)*(-12.*pow(mpion,4) + 3.*pow(mrho,4) - 2.*pow(mrho,2)*s - 1.*pow(s,2) + pow(mpion,2)*(4.*pow(mrho,2) + 4.*s))
                                          ))*log(fabs(-pow(ma1,2) + t1)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
                                  (0.25*(1.*eta1 - 1.*eta2)*(delta*(eta1*(1.*pow(ma1,6) - 1.*pow(mpion,6) + pow(mpion,4)*(-2.5*pow(mrho,2) + 0.5*s) +
                                             pow(mpion,2)*s*(-0.5*pow(mrho,2) + 1.*s) + pow(ma1,4)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s) +
                                             s*(0.5*pow(mrho,4) - 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
                                             pow(ma1,2)*(3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s +
                                                1.*pow(s,2))) + eta2*(-1.*pow(ma1,6) + pow(ma1,4)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
                                             pow(mpion,2)*(1.*pow(mpion,4) - 0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
                                             pow(ma1,2)*(-3.*pow(mpion,4) + 1.*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2) +
                                                pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)))) +
                                       pow(mrho,2)*(eta2*(pow(ma1,4)*(-1. + 2.*C4*pow(mrho,2)) +
                                             pow(mpion,2)*(0.5*pow(mrho,2) + pow(mpion,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s) +
                                             pow(ma1,2)*(2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) + pow(mrho,2)*(-1.5 - 4.*C4*s) +
                                                s*(0.5 + 2.*C4*s))) + eta1*(pow(ma1,4)*(1. - 2.*C4*pow(mrho,2)) + pow(mpion,4)*(1. - 2.*C4*pow(mrho,2)) +
                                             (-0.5*pow(mrho,2) + 0.5*s)*s + pow(ma1,2)*
                                              (-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) - 2.*C4*pow(s,2)) +
                                             pow(mpion,2)*(-4.*C4*pow(mrho,4) - 1.*s + pow(mrho,2)*(2. + 4.*C4*s)))))*log(fabs(-pow(ma1,2) + t1)))/pow(mrho,2) -
                                  0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-pow(mpion,2) + t1)) -
                                  (0.5*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*(eta2*pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s) +
                                       eta1*(-0.5*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s))*log(fabs(-pow(mpion,2) + t1)))/
                                   (-1.*pow(ma1,2) + 1.*pow(mpion,2)) + (2.*(-0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,4)*s +
                                       pow(mpion,2)*(C4*(-2. + 1.*delta)*pow(mrho,6) + (0.5 - 0.25*delta)*delta*pow(s,2) +
                                          pow(mrho,4)*(1. - 0.5*delta + 4.*C4*s - 2.*C4*delta*s) +
                                          pow(mrho,2)*s*(1. - 2.*delta + 0.75*pow(delta,2) - 2.*C4*s + 1.*C4*delta*s)))*log(fabs(-pow(mpion,2) + t1)))/
                                   (pow(mrho,4) - 1.*pow(mrho,2)*s)))/(16.*Pi*(4*pow(mpion,2) - s)*s));

  return xs;
}

double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi0_rho0_pi0(
    const double s, const double t, const double m_rho) {
  const double m_pi = m_pion_;
  double diff_xsection =
      1 / 3.0 *
      (pow(Const, 2) * pow(g_POR, 4) *
       (pow(m_omega, 4) * pow(s, 4) + 4 * pow(m_omega, 4) * pow(s, 3) * t -
        4 * pow(m_omega, 2) * pow(s, 4) * t +
        10 * pow(m_omega, 4) * pow(s, 2) * pow(t, 2) -
        16 * pow(m_omega, 2) * pow(s, 3) * pow(t, 2) +
        5 * pow(s, 4) * pow(t, 2) + 4 * pow(m_omega, 4) * s * pow(t, 3) -
        16 * pow(m_omega, 2) * pow(s, 2) * pow(t, 3) +
        10 * pow(s, 3) * pow(t, 3) + pow(m_omega, 4) * pow(t, 4) -
        4 * pow(m_omega, 2) * s * pow(t, 4) + 5 * pow(s, 2) * pow(t, 4) +
        pow(m_pi, 8) * pow(-2 * pow(m_omega, 2) + s + t, 2) -
        2 * pow(m_pi, 6) * pow(m_rho, 2) *
            (2 * pow(m_omega, 4) + pow(s, 2) + pow(t, 2) -
             2 * pow(m_omega, 2) * (s + t)) +
        pow(m_rho, 4) *
            (2 * pow(s, 2) * pow(t, 2) - 2 * pow(m_omega, 2) * s * t * (s + t) +
             pow(m_omega, 4) * (pow(s, 2) + pow(t, 2))) -
        2 * pow(m_rho, 2) *
            (3 * pow(s, 2) * pow(t, 2) * (s + t) -
             3 * pow(m_omega, 2) * s * t * pow(s + t, 2) +
             pow(m_omega, 4) * (pow(s, 3) + 2 * pow(s, 2) * t +
                                2 * s * pow(t, 2) + pow(t, 3))) +
        pow(m_pi, 4) *
            (-2 * pow(m_rho, 2) * (pow(m_omega, 2) - s) *
                 (pow(m_omega, 2) - t) * (s + t) -
             8 * pow(m_omega, 2) * s * t * (s + t) +
             4 * pow(m_omega, 4) * (pow(s, 2) + pow(t, 2)) -
             2 * s * t * (pow(s, 2) - 6 * s * t + pow(t, 2)) +
             pow(m_rho, 4) * (2 * pow(m_omega, 4) + pow(s, 2) + pow(t, 2) -
                              2 * pow(m_omega, 2) * (s + t))) -
        2 * pow(m_pi, 2) *
            (2 * (s + t) * pow(-2 * s * t + pow(m_omega, 2) * (s + t), 2) +
             pow(m_rho, 4) * (-4 * pow(m_omega, 2) * s * t +
                              pow(m_omega, 4) * (s + t) + s * t * (s + t)) -
             pow(m_rho, 2) *
                 (-10 * pow(m_omega, 2) * s * t * (s + t) +
                  2 * pow(m_omega, 4) * (pow(s, 2) + 3 * s * t + pow(t, 2)) +
                  s * t * (pow(s, 2) + 8 * s * t + pow(t, 2)))))) /
      (128. * Pi * pow(pow(m_omega, 2) - s, 2) *
       (pow(pow(m_pi, 2) - pow(m_rho, 2), 2) -
        2 * (pow(m_pi, 2) + pow(m_rho, 2)) * s + pow(s, 2)) *
       pow(pow(m_omega, 2) - t, 2));

  return to_mb * diff_xsection;
}

double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_pi_rho0(
    const double s, const double t, const double m_rho) {
  const double mpion = m_pion_, mrho = m_rho;
  const double diff_xsection =
      (pow(Const, 2) * pow(ghat, 4) *
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
              pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t, 2))) +
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
              pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t, 2))) +
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
        (2 *
         ((0.125 * pow(-2 + delta, 2) * (2 * pow(mpion, 2) - s) *
           (pow(mpion, 4) + pow(mrho, 2) * (s - t) + t * (s + t) -
            pow(mpion, 2) * (3 * pow(mrho, 2) + s + 2 * t))) /
              ((pow(mpion, 2) - t) * (pow(mpion, 2) + pow(mrho, 2) - s - t)) -
          (0.125 * (-2. + delta) *
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
              (pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t))) /
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
                   (pow(mrho, 4) - pow(s, 2) - pow(mrho, 2) * t + t * (s + t)) +
               2 * pow(mpion, 4) *
                   (2 * pow(mrho, 4) + t * (s + 3 * t) -
                    pow(mrho, 2) * (2 * s + 3 * t)) +
               pow(mpion, 2) *
                   (pow(mrho, 6) - 2 * pow(mrho, 4) * (s + 2 * t) +
                    2 * t * (pow(s, 2) - 2 * s * t - 2 * pow(t, 2)) +
                    pow(mrho, 2) * (pow(s, 2) + 2 * s * t + 6 * pow(t, 2)))))) /
            ((-pow(ma1, 2) + t) *
             (pow(Gammaa1, 2) * pow(ma1, 2) +
              pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t, 2))) +
        2 * ((-0.0625 * (-2 + delta) * (eta1 - eta2) *
              (eta1 * (-2 * pow(mpion, 2) + s) + eta2 * (pow(mpion, 2) + t)) *
              (-pow(mpion, 4) + pow(s, 2) - pow(t, 2) +
               pow(mrho, 2) * (-s + t) +
               pow(mpion, 2) * (pow(mrho, 2) - 2 * s + 2 * t))) /
                 ((pow(mpion, 2) + pow(mrho, 2) - s - t) * (-pow(ma1, 2) + t)) +
             (0.125 * (-2. + delta) * (eta1 - 1. * eta2) *
              (eta2 *
                   (-0.5 * pow(mpion, 6) +
                    pow(mpion, 4) * (0.5 * pow(mrho, 2) + 0.5 * t) +
                    pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s + 0.5 * t) * t +
                    (0.5 * pow(mrho, 2) - 1. * s - 0.5 * t) * pow(t, 2)) +
               eta1 * (1. * pow(mpion, 6) +
                       pow(mpion, 4) * (-1. * pow(mrho, 2) + 0.5 * s - 2. * t) +
                       s * (-0.5 * pow(mrho, 2) + 0.5 * t) * t +
                       pow(mpion, 2) * (1. * pow(mrho, 4) +
                                        pow(mrho, 2) * (-0.5 * s - 1. * t) +
                                        t * (1. * s + 1. * t))))) /
                 ((pow(ma1, 2) - 1. * t) * (-1. * pow(mpion, 2) + t)) +
             (0.0625 * (eta1 - eta2) *
              (eta2 *
                   (8 * C4 * pow(mrho, 6) * t - 2 * delta * pow(s, 2) * t +
                    pow(mrho, 2) *
                        (-4 * pow(mpion, 4) +
                         (s * (2 + 3 * delta + 8 * C4 * s) - 4 * t) * t +
                         pow(mpion, 2) * (-((-2 + delta) * s) + 8 * t)) +
                    pow(mrho, 4) * (8 * C4 * pow(mpion, 4) -
                                    pow(mpion, 2) * (-2 + delta + 16 * C4 * t) +
                                    t * (-6 + delta + 8 * C4 * (-2 * s + t)))) +
               eta1 * (2 * delta * pow(s, 2) * t +
                       8 * C4 * pow(mrho, 6) * (-2 * pow(mpion, 2) + t) -
                       pow(mrho, 2) *
                           (-4 * pow(mpion, 4) - 4 * pow(t, 2) +
                            2 * pow(mpion, 2) * ((2 + delta) * s + 4 * t) +
                            pow(s, 2) * (-2 + delta + 8 * C4 * t)) +
                       pow(mrho, 4) *
                           (-8 * C4 * pow(mpion, 4) + (-2 + delta) * s -
                            4 * t * (1 + 2 * C4 * t) +
                            8 * pow(mpion, 2) * (1 + 2 * C4 * (s + t)))))) /
                 (pow(mrho, 2) * (-pow(ma1, 2) + t))) -
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
      (16. * M_PI * s * (-4 * pow(mpion, 2) + s));

  return to_mb * diff_xsection;
}

double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_pi0_rho(
    const double s, const double t, const double m_rho) {
  const double mpion = m_pion_, mrho = m_rho;
  const double diff_xsection =
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

  return to_mb * diff_xsection;
}

double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_rho_pi0(
    const double s, const double t, const double m_rho) {
  // omega:
  const double mpion = m_pion_, mrho = m_rho;
  const double diff_xsection =
      1 / 3.0 *
      (0.0024867959858108648 * pow(Const, 2) * pow(g_POR, 4) *
       (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
        pow(mpion, 4) * (pow(mrho, 4) + 4 * pow(s, 2) - 2 * s * t) +
        pow(s, 2) * (pow(mrho, 4) + pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                     2 * pow(mrho, 2) * (s + t)) -
        2 * pow(mpion, 2) * s *
            (pow(mrho, 4) + 2 * s * (s + t) - pow(mrho, 2) * (2 * s + t)))) /
      (pow(pow(m_omega, 2) - s, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                                      2 * pow(mpion, 2) * (pow(mrho, 2) + s)));

  return diff_xsection * to_mb;
}

double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_rho0_pi(
    const double s, const double t, const double m_rho) {
  const double mrho = m_rho;
  const double m_pi = m_pion_;

  const double diff_xsection =
      1 / 3.0 *
      (pow(Const, 2) * pow(ghat, 4) *
       ((-8 * pow(-2 + delta, 2) * pow(m_pi, 2)) /
            (pow(mrho, 2) * pow(pow(m_pi, 2) - s, 2)) -
        (8 * pow(-2 + delta, 2) * pow(m_pi, 2) *
         (pow(m_pi, 4) + pow(pow(mrho, 2) - t, 2) -
          2 * pow(m_pi, 2) * (pow(mrho, 2) + t))) /
            (pow(mrho, 2) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) *
             pow(pow(m_pi, 2) - t, 2)) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (-(eta2 * (pow(m_pi, 2) + s)) + eta1 * (-pow(mrho, 2) + s + t)) *
         (-pow(m_pi, 4) + pow(m_pi, 2) * (pow(mrho, 2) - 2 * t) +
          t * (-pow(mrho, 2) + 2 * s + t))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) *
             (pow(m_pi, 2) - t)) -
        (8 * (-2 + delta) *
         (pow(m_pi, 4) * (2 - 3 * delta + 8 * C4 * pow(mrho, 2)) +
          pow(mrho, 4) * (-2 + delta + 8 * C4 * t) +
          t * ((2 + 3 * delta) * s + 2 * delta * t) +
          pow(m_pi, 2) *
              (-8 * C4 * pow(mrho, 4) + (-2 + delta) * s - (2 + 3 * delta) * t +
               4 * pow(mrho, 2) * (1 + 4 * C4 * t)) -
          pow(mrho, 2) * (t * (-2 + 3 * delta + 8 * C4 * t) +
                          s * (-2 + delta + 16 * C4 * t)))) /
            (pow(mrho, 2) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) *
             (pow(m_pi, 2) - t)) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(m_pi, 2) + s) *
              (pow(m_pi, 4) - pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) +
               s * (pow(mrho, 2) - s - 2 * t)) +
          eta1 * (-4 * pow(m_pi, 6) +
                  s * (-pow(mrho, 2) + s) * (-pow(mrho, 2) + s + t) +
                  pow(m_pi, 4) * (3 * pow(mrho, 2) + s + t) -
                  pow(m_pi, 2) * (pow(mrho, 4) + 2 * s * (s - t) +
                                  pow(mrho, 2) * (-s + t))))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta2, 2) *
              (pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) +
               pow(m_pi, 4) * (pow(pow(mrho, 2) + 2 * s, 2) - 2 * s * t) +
               pow(s, 2) * (pow(pow(mrho, 2) + s, 2) +
                            2 * (-pow(mrho, 2) + s) * t + 2 * pow(t, 2)) -
               2 * pow(m_pi, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * (2 * s - t) +
                    2 * s * (s + t))) +
          2 * eta1 * eta2 *
              (-pow(m_pi, 8) +
               pow(m_pi, 4) * (pow(mrho, 4) + 2 * pow(mrho, 2) * s +
                               2 * s * (-2 * s + t)) -
               2 * pow(m_pi, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * (s + t) - 2 * s * (s + t)) +
               pow(s, 2) * (pow(mrho, 4) - pow(s, 2) + 2 * pow(mrho, 2) * t -
                            2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) +
               pow(m_pi, 4) * (3 * pow(mrho, 4) + 2 * s * (2 * s - t) +
                               2 * pow(mrho, 2) * (-3 * s + t)) -
               2 * pow(m_pi, 2) * (pow(mrho, 2) - s) *
                   (-2 * s * (s + t) + pow(mrho, 2) * (2 * s + t)) +
               s * (-pow(mrho, 2) + s) *
                   (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                    pow(mrho, 2) * (s + 2 * t))))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(m_pi, 8) -
               pow(m_pi, 4) *
                   (pow(mrho, 4) + 2 * (pow(mrho, 2) + s) * t - 4 * pow(t, 2)) +
               pow(t, 2) * (-pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                            2 * pow(s, 2) + 2 * s * t + pow(t, 2)) +
               2 * pow(m_pi, 2) * t *
                   (pow(mrho, 4) + pow(mrho, 2) * (s + t) - 2 * t * (s + t))) +
          pow(eta2, 2) *
              (pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) +
               pow(m_pi, 4) *
                   (pow(mrho, 4) + 4 * pow(mrho, 2) * t - 2 * (s - 2 * t) * t) +
               pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t +
                            pow(t, 2) + 2 * pow(mrho, 2) * (-s + t)) -
               2 * pow(m_pi, 2) * t *
                   (pow(mrho, 4) - pow(mrho, 2) * (s - 2 * t) +
                    2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(m_pi, 8) - 2 * pow(m_pi, 6) * pow(mrho, 2) +
               pow(m_pi, 4) *
                   (3 * pow(mrho, 4) + 2 * pow(mrho, 2) * (s - 3 * t) -
                    2 * (s - 2 * t) * t) +
               t * (-pow(mrho, 2) + t) *
                   (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                    pow(mrho, 2) * (2 * s + t)) -
               2 * pow(m_pi, 2) * (-pow(mrho, 2) + t) *
                   (2 * t * (s + t) - pow(mrho, 2) * (s + 2 * t))))) /
            ((pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) *
             pow(pow(ma1, 2) - t, 2)) +
        (8 * (-2 + delta) *
         ((-2 + delta) * pow(mrho, 6) +
          pow(m_pi, 6) * (-2 + 3 * delta - 8 * C4 * pow(mrho, 2)) +
          s * t * ((-2 + 3 * delta) * s + 4 * delta * t) +
          pow(m_pi, 4) *
              (8 * C4 * pow(mrho, 4) + 4 * delta * s + 2 * t - 3 * delta * t -
               pow(mrho, 2) * (2 + delta + 16 * C4 * s - 8 * C4 * t)) +
          pow(mrho, 4) *
              (-((-2 + delta) * t) + s * (4 - 2 * delta + 8 * C4 * t)) +
          pow(mrho, 2) * s *
              (s * (-2 + delta - 8 * C4 * t) - 2 * t * (delta + 8 * C4 * t)) +
          pow(m_pi, 2) *
              (s * ((2 - 3 * delta) * s - 8 * delta * t) -
               pow(mrho, 4) * (-6 + 3 * delta + 8 * C4 * (s + t)) +
               pow(mrho, 2) * (8 * C4 * pow(s, 2) + 4 * (-1 + delta) * t +
                               s * (-8 + 6 * delta + 32 * C4 * t))))) /
            (pow(mrho, 2) * (pow(m_pi, 2) - s) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) *
             (pow(m_pi, 2) - t)) +
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta1, 2) * (pow(m_pi, 8) +
                          pow(m_pi, 4) * (2 * pow(mrho, 4) + 2 * s * t -
                                          3 * pow(mrho, 2) * (s + t)) +
                          s * t *
                              (2 * pow(mrho, 4) + pow(s, 2) + 3 * s * t +
                               pow(t, 2) - 3 * pow(mrho, 2) * (s + t)) -
                          2 * pow(m_pi, 2) * (pow(mrho, 2) - s - t) *
                              (-2 * s * t + pow(mrho, 2) * (s + t))) +
          pow(eta2, 2) * (pow(m_pi, 8) -
                          4 * pow(m_pi, 2) * s * t * (pow(mrho, 2) + s + t) +
                          pow(m_pi, 4) * (2 * s * t + pow(mrho, 2) * (s + t)) +
                          s * t *
                              (pow(s, 2) + 3 * s * t + pow(t, 2) +
                               pow(mrho, 2) * (s + t))) +
          2 * eta1 * eta2 *
              (-pow(m_pi, 8) + 2 * pow(m_pi, 6) * pow(mrho, 2) -
               2 * pow(m_pi, 4) * s * t -
               s * t *
                   (pow(s, 2) + 3 * s * t + pow(t, 2) -
                    2 * pow(mrho, 2) * (s + t)) -
               pow(m_pi, 2) *
                   (-4 * s * t * (s + t) +
                    pow(mrho, 2) * (pow(s, 2) + 4 * s * t + pow(t, 2)))))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s)) *
             (pow(ma1, 2) - t)) +
        (8 *
         (pow(delta, 2) *
              (8 * pow(m_pi, 4) + 3 * pow(mrho, 4) -
               6 * pow(mrho, 2) * (s + t) + 2 * pow(s + t, 2) +
               4 * pow(m_pi, 2) * (3 * pow(mrho, 2) - 2 * (s + t))) -
          4 * delta * pow(mrho, 2) *
              (16 * C4 * pow(m_pi, 4) + pow(mrho, 2) * (3 - 6 * C4 * (s + t)) +
               (s + t) * (-3 + 4 * C4 * (s + t)) +
               2 * pow(m_pi, 2) * (3 + C4 * (6 * pow(mrho, 2) - 8 * (s + t)))) +
          4 * pow(mrho, 4) *
              (3 + 4 * C4 * (2 * pow(m_pi, 2) - s - t) *
                       (3 + C4 * (4 * pow(m_pi, 2) - 2 * (s + t)))))) /
            (pow(mrho, 4) * (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (-pow(ma1, 2) + s) *
         (eta2 * (-2 * pow(m_pi, 4) * (delta - 4 * C4 * pow(mrho, 2)) *
                      (pow(mrho, 2) + 4 * s) +
                  pow(m_pi, 2) *
                      (-2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * s) +
                       8 * delta * s * (s + t) -
                       pow(mrho, 2) * ((-10 + delta) * s - (-2 + delta) * t +
                                       32 * C4 * s * (s + t))) +
                  s * (2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * s) -
                       2 * delta * pow(s + t, 2) +
                       pow(mrho, 2) * ((-6 + delta) * s + (-2 + delta) * t +
                                       8 * C4 * pow(s + t, 2)))) +
          eta1 *
              (4 * pow(m_pi, 4) *
                   (6 * C4 * pow(mrho, 4) + 2 * delta * s +
                    pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * s)) +
               2 * delta * s * pow(s + t, 2) -
               pow(mrho, 2) *
                   ((-6 + 5 * delta) * pow(s, 2) +
                    2 * (-2 + 3 * delta) * s * t + (-2 + delta) * pow(t, 2) +
                    8 * C4 * s * pow(s + t, 2)) +
               pow(mrho, 4) *
                   ((-2 + delta) * (3 * s + t) + 8 * C4 * s * (s + 2 * t)) -
               2 * pow(m_pi, 2) *
                   (4 * delta * s * (s + t) -
                    pow(mrho, 2) * (-6 * s + 7 * delta * s - 2 * t +
                                    3 * delta * t + 16 * C4 * s * (s + t)) +
                    2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * (2 * s + t)))))) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(m_pi, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) *
         (((-2 + delta) *
           (pow(m_pi, 4) - pow(m_pi, 2) * (pow(mrho, 2) - 2 * s) +
            s * (pow(mrho, 2) - s - 2 * t)) *
           (eta1 * (pow(mrho, 2) - s - t) + eta2 * (pow(m_pi, 2) + t))) /
              ((pow(m_pi, 2) - s) * (pow(ma1, 2) - t)) +
          ((-2 + delta) *
           (eta2 * (pow(m_pi, 2) + t) *
                (pow(m_pi, 4) - pow(m_pi, 2) * (pow(mrho, 2) - 2 * t) +
                 (pow(mrho, 2) - 2 * s - t) * t) +
            eta1 * (-4 * pow(m_pi, 6) +
                    (pow(mrho, 2) - t) * (pow(mrho, 2) - s - t) * t +
                    pow(m_pi, 4) * (3 * pow(mrho, 2) + s + t) -
                    pow(m_pi, 2) * (pow(mrho, 4) + pow(mrho, 2) * (s - t) +
                                    2 * t * (-s + t))))) /
              ((-pow(ma1, 2) + t) * (-pow(m_pi, 2) + t)) +
          (eta2 *
               (-2 * pow(m_pi, 4) * (delta - 4 * C4 * pow(mrho, 2)) *
                    (pow(mrho, 2) + 4 * t) +
                pow(m_pi, 2) *
                    (8 * delta * t * (s + t) -
                     2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * t) -
                     pow(mrho, 2) * (-((-2 + delta) * s) + (-10 + delta) * t +
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
                4 * pow(m_pi, 4) *
                    (6 * C4 * pow(mrho, 4) + 2 * delta * t +
                     pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * t)) -
                2 * pow(m_pi, 2) *
                    (4 * delta * t * (s + t) -
                     pow(mrho, 2) * (-2 * s + 3 * delta * s - 6 * t +
                                     7 * delta * t + 16 * C4 * t * (s + t)) +
                     2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * (s + 2 * t))))) /
              (pow(mrho, 2) * (-pow(ma1, 2) + t)))) /
            (pow(m_pi, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(m_pi, 2) * (pow(mrho, 2) + s)))) /
      (512. * Pi);

  return to_mb * diff_xsection;
}
double PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi0_rho_pi(
    const double s, const double t, const double m_rho) {
  const double mpion = m_pion_, mrho = m_rho;

  const double diff_xsection =
      1 / 3.0 *
      (0.0024867959858108648 * pow(Const, 2) * pow(g_POR, 4) *
       (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
        pow(mpion, 4) * (pow(mrho, 4) - 2 * (s - 2 * t) * t) +
        pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                     2 * pow(mrho, 2) * (s + t)) -
        2 * pow(mpion, 2) * t *
            (pow(mrho, 4) + 2 * t * (s + t) - pow(mrho, 2) * (s + 2 * t)))) /
      ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
        2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
       pow(pow(m_omega, 2) - t, 2));

  return to_mb * diff_xsection;
}

// definition of total-xs getters
double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi_pi_rho0(
    const double s, const double m_rho) {
  if (tab_pi_pi_rho0_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi_pi_rho0\n");
    tab_pi_pi_rho0_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_pi_rho0);
  }

  return tab_pi_pi_rho0_->get_linear(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi_pi0_rho(
    const double s, const double m_rho) {
  if (tab_pi_pi0_rho_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi_pi0_rho\n");
    tab_pi_pi0_rho_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_pi0_rho);
  }

  return tab_pi_pi0_rho_->get_linear(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi0_rho0_pi0(
    const double s, const double m_rho) {
  if (tab_pi0_rho0_pi0_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi0_rho0_pi0\n");
    tab_pi0_rho0_pi0_ =make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho0_pi0);
  }

  return tab_pi0_rho0_pi0_->get_linear(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi_rho0_pi(
    const double s, const double m_rho) {
  if (tab_pi_rho0_pi_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi_rho0_pi\n");
    tab_pi_rho0_pi_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_rho0_pi);
  }

  return tab_pi_rho0_pi_->get_linear(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi0_rho_pi(
    const double s, const double m_rho) {
  if (tab_pi0_rho_pi_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi0_rho_pi\n");
    tab_pi0_rho_pi_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_pi0_rho_pi);
  }

  return tab_pi0_rho_pi_->get_linear(s, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_pi_rho_pi0(
    const double s, const double m_rho) {
  if (tab_pi_rho_pi0_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info("init table pi_rho_pi0\n");
    tab_pi_rho_pi0_ = make_unique<TabulationND<2>>(
        s0_tot, s1_tot, m_rho_0, m_rho_1, ds_tot, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_pi_rho_pi0);
  }

  return tab_pi_rho_pi0_->get_linear(s, m_rho);
}

// definition of differential xs-getters
double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi_pi_rho0(
    const double s, const double t, const double m_rho) {
  if (tab_pi_pi_rho0_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi_pi_rho0_diff\n";
    tab_pi_pi_rho0_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_pi_rho0);
  }

  return tab_pi_pi_rho0_diff_->get_linear(s, t, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi_pi0_rho(
    const double s, const double t, const double m_rho) {
  if (tab_pi_pi0_rho_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi_pi0_rho_diff\n";
    tab_pi_pi0_rho_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_pi0_rho);
  }

  return tab_pi_pi0_rho_diff_->get_linear(s, t, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi0_rho0_pi0(
    const double s, const double t, const double m_rho) {
  if (tab_pi0_rho0_pi0_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi0_rho0_pi0_diff\n";
    tab_pi0_rho0_pi0_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi0_rho0_pi0);
  }

  return tab_pi0_rho0_pi0_diff_->get_linear(s, t, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi_rho_pi0(
    const double s, const double t, const double m_rho) {
  if (tab_pi_rho_pi0_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi_rho_pi0\n";
    tab_pi_rho_pi0_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_rho_pi0);
  }

  return tab_pi_rho_pi0_diff_->get_linear(s, t, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi0_rho_pi(
    const double s, const double t, const double m_rho) {
  if (tab_pi0_rho_pi_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi0_rho_pi\n";
    tab_pi0_rho_pi_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_rho0_pi);
  }

  return tab_pi0_rho_pi_diff_->get_linear(s, t, m_rho);
}

double PhotonCrossSection<ComputationMethod::Lookup>::xs_diff_pi_rho0_pi(
    const double s, const double t, const double m_rho) {
  if (tab_pi_rho0_pi_diff_ == nullptr) {
    const auto &log = logger<LogArea::Experiment>();
    log.info() << "init table pi_rho0_pi_diff\n";
    tab_pi_rho0_pi_diff_ = make_unique<TabulationND<3>>(
        s0_diff, s1_diff, t0_diff, t1_diff, m_rho_0, m_rho_1, ds_diff, dt_diff, dm,
        PhotonCrossSection<ComputationMethod::Analytic>::xs_diff_pi_rho0_pi);
  }

  return tab_pi_rho0_pi_diff_->get_linear(s, t, m_rho);
}

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
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi0_rho_pi_ = nullptr;

std::unique_ptr<TabulationND<2>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_rho_pi0_ = nullptr;

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
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_rho_pi0_diff_ =
        nullptr;
std::unique_ptr<TabulationND<3>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi0_rho_pi_diff_ =
        nullptr;
std::unique_ptr<TabulationND<3>>
    PhotonCrossSection<ComputationMethod::Lookup>::tab_pi_rho0_pi_diff_ =
        nullptr;
