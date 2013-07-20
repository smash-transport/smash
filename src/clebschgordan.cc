/*
 * Class for calculating Clebsch-Gordan coefficients
 * Based on Caswell, R. S. and  Maximon, L. C.,
 * "Fortran programs for the calculation of Wigner 3j, 6j, and 9j
 * coefficients for angular momenta greater than or equal to 80",
 * NBS Technical Note 409 (1966).
 * http://archive.org/details/fortranprogramsf409casw
 *
 * Copyright (c) 2013
 *    Kai Gallmeister <gallmei@th.physik.uni-frankfurt.de>
 *    Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *  GNU General Public License (GPLv3)
 */

#include <math.h>
#include <stdlib.h>
#include <algorithm>

#include "include/ClebschGordan.h"

ClebschGordan::ClebschGordan() : factorial_log_max(4 * 50 + 2) {
  /* factorial log = log[(N-1)!] */
  factorial_log = new double[factorial_log_max + 1];
  factorial_log[0] = factorial_log[1] = factorial_log[2] = 0.0;
  double factorial_n = 1.0;
  for (unsigned int i = 3; i <= factorial_log_max; i++) {
    factorial_n += 1.0;
    /* log[n!] = log[(n-1)! * n] = log[(n-1)!] + log[n] */
    factorial_log[i] = factorial_log[i-1] + log(factorial_n);
  }
}

ClebschGordan::~ClebschGordan() {
  delete[] factorial_log;
}

/* ClebschGordan - calculate the (iso)spin combination coefficients
 * Note that both j and m -inputs are double the actual spin values!
 */
double ClebschGordan::operator()(int j1, int j2, int j3, int m1, int m2,
                                 int m3) {
  /* No negative spin numbers */
  if (j1 < 0 || j2 < 0 || j3 < 0)
    return 0.0;
  /* z-component quantum number m must be between -j and j */
  if (abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3)
    return 0.0;
  /* For total spin, |j1 + j2| < j3 < j1 + j2 */
  if (j3 > j1 + j2 || j3 < abs(j1 - j2))
    return 0.0;
  /* z-components must add up */
  if (m1 + m2 != m3)
    return 0.0;
  /* actual spin values must add to an integer */
  if (((j1 + j2 + j3) % 2) != 0)
    return 0.0;
  /* the actual total spin + z-component values
   * must also add to an integer for each particle
   * individually
   */
  if (((j1 + m1) % 2) != 0 || ((j2 + m2) % 2) != 0 || ((j3 + m3) % 2) != 0)
    return 0.0;
  /* If either of initial spins is zero, combination is trivial */
  if (j1 == 0 || j2 == 0)
    return 1.0;

  /* Return the Clebsch-Gordan coefficient */
  return pow(-1, (j1 - j2 + m3) / 2) * sqrt(static_cast<double>(j3 + 1))
    * f3j(j1, j2, j3, m1, m2, -m3);
}

/* f3j - Calculation of Wigner 3j symbol <j1 m1 j2 m2 | j3 m3>
 * This is basically equation (6) in Caswell & Maximon (see also Eq. (1'))
 * Input is twice the actual value
 */
double ClebschGordan::f3j(int j1, int j2, int j3, int m1, int m2, int m3) {
  int spinfactor[9];

  /* We'll need the factorials of these spin terms */
  spinfactor[0] = (j1 + j2 - j3) / 2;
  spinfactor[1] = (j1 - j2 + j3) / 2;
  spinfactor[2] = (-j1 + j2 + j3) / 2;
  spinfactor[3] = (j1 + m1) / 2;
  spinfactor[4] = (j1 - m1) / 2;
  spinfactor[5] = (j2 + m2) / 2;
  spinfactor[6] = (j2 - m2) / 2;
  spinfactor[7] = (j3 + m3) / 2;
  spinfactor[8] = (j3 - m3) / 2;

  /* Define sum limits */
  /* lower limit is the absolute value of the most negative
   * of 0, j3 - j2 + m1 and j3 - j1 - m2
   */
  int kmin = std::max(-j3 + j2 - m1, std::max(-j3 + j1 + m2, 0)) / 2;

  /* Upper limit is the smallest of j1 + j2 - j3, j1 - m1 and j2 + m2 */
  int kmax = (j2 - j3 + m1 < 0) ? j1 + j2 - j3 : j1 - m1;
  kmax = std::min(j2 + m2, kmax) / 2;

  /* k-dependent spin factors we'll be taking factorials of */
  int k_spinfactor1 = spinfactor[0] - kmin + 1;
  int k_spinfactor2  = spinfactor[4] - kmin + 1;
  int k_spinfactor3 = spinfactor[5] - kmin + 1;
  int k_spinfactor4 = (j3 - j2 + m1) / 2 + kmin;
  int k_spinfactor5 = (j3 - j1 - m2) / 2 + kmin;

  /* sum series in double precision */
  double sumterm_k = 1.0e-10;
  double sum = 1.0e-10;
  int ncut = 0;
  kmax -= kmin;
  if (kmax > 0) {
    for (int k = 1; k <= kmax; k++) {
      sumterm_k = -sumterm_k * (static_cast<double>((k_spinfactor1 - k)
                  * (k_spinfactor2 - k) * (k_spinfactor3 - k)))
        / (static_cast<double>((kmin + k) * (k_spinfactor4 + k)
                               * (k_spinfactor5 + k)));
      /* If terms are getting very large, renormalize */
      if (fabs(sumterm_k) >= 1.0e30) {
        sumterm_k *= 1.0e-10;
        sum *= 1.0e-10;
        ncut++;
      }
      /* If terms are getting very small, cut the sum */
      if (fabs(sumterm_k) < 1.0e-20)
        break;
      sum += sumterm_k;
    }
  }

  /* Calculate the sum prefactor */
  double spinfactor_log = 0.0;
  for (unsigned int i = 0; i < 9; i++)
    spinfactor_log += factorial_log[spinfactor[i] + 1];

  int numerator = (j1 + j2 + j3) / 2 + 2;
  spinfactor_log = (spinfactor_log - factorial_log[numerator]) / 2;
  double sum_min_log = -factorial_log[kmin + 1] - factorial_log[k_spinfactor1]
    - factorial_log[k_spinfactor2] - factorial_log[k_spinfactor3]
    - factorial_log[k_spinfactor4 + 1] - factorial_log[k_spinfactor5 + 1];
  double prefactor_log = spinfactor_log + sum_min_log;

  /* Compute the final result */
  double result;
  if ((prefactor_log < -80.0) || (ncut > 0)) {
    double sign = (sum > 0) ? 1 : ((sum < 0) ? -1 : 0);
    sum = fabs(sum);
    double sum_log = log(sum) + (static_cast<double>(ncut + 1)) * log(1.0e10);
    result = sign * exp(sum_log + prefactor_log);
  } else {
    sum *= 1.0e10;
    result = exp(prefactor_log) * sum;
  }
  numerator = kmin + (j1 - j2 - m3) / 2;
  if ((numerator % 2) != 0)
    result = -result;

  return result;
}

/* Check if it's possible to combine j1 and j2 to j3 */
bool ClebschGordan::MayBranch(int j1, int j2, int j3) {
  for (int m1 = -j1; m1 <= j1; m1 += 2) {
    for (int m2 = -j2; m2 <= j2; m2 += 2) {
      if (this->operator()(j1, j2, j3, m1, m2, m1+m2) > 0)
        return true;
    }
  }
  return false;
}
