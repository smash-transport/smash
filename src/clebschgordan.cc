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

ClebschGordan::ClebschGordan() : flmax(4 * 50 + 2) {
  fl = new double[flmax + 1];
  fl[0] = fl[1] = fl[2] = 0.0;
  double fn = 1.0;
  for (unsigned int i = 3; i <= flmax; i++) {
    fn += 1.0;
    fl[i] = fl[i-1] + log(fn);
  }
}

ClebschGordan::~ClebschGordan() {
  delete[] fl;
}

/* ClebschGordan - we calculate funky coefficients */
double ClebschGordan::operator()(int j1, int j2, int j3, int m1, int m2,
                                 int m3) {
  /* XXX: document what is happening */
  if (j1 < 0 || j2 < 0 || j3 < 0)
    return 0.0;
  if (abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3)
    return 0.0;
  if (j3 > j1 + j2 || j3 < abs(j1 - j2))
    return 0.0;
  if (m1 + m2 != m3)
    return 0.0;
  if (((j1 + j2 + j3) % 2) != 0)
    return 0.0;
  if (((j1 + m1) % 2) != 0 || ((j2 + m2) % 2) != 0 || ((j3 + m3) % 2) != 0)
    return 0.0;
  if (j1 == 0 || j2 == 0)
    return 1.0;

  return pow(-1, (j1 - j2 + m3) / 2) * sqrt(static_cast<double>(j3 + 1))
    * f3j(j1, j2, j3, m1, m2, -m3);
}

double ClebschGordan::f3j(int j1, int j2, int j3, int m1, int m2, int m3) {
  int mtri[10];

  /* XXX: c and c++ start to count at zero */
  mtri[1] = (j1 + j2 - j3) / 2;
  mtri[2] = (j1 - j2 + j3) / 2;
  mtri[3] = (-j1 + j2 + j3) / 2;
  mtri[4] = (j1 + m1) / 2;
  mtri[5] = (j1 - m1) / 2;
  mtri[6] = (j2 + m2) / 2;
  mtri[7] = (j2 - m2) / 2;
  mtri[8] = (j3 + m3) / 2;
  mtri[9] = (j3 - m3) / 2;

  int kmin, kmax;

  kmin = std::max(-j3 + j2 - m1, std::max(-j3 + j1 + m2, 0)) / 2;
  kmax = (j2 - j3 + m1 < 0) ? j1 + j2 - j3 : j1 - m1;
  kmax = std::min(j2 + m2, kmax) / 2;

  int mini, min2, min3, min4, min5;

  mini = mtri[1] - kmin + 1;
  min2 = mtri[5] - kmin + 1;
  min3 = mtri[6] - kmin + 1;
  min4 = (j3 - j2 + m1) / 2 + kmin;
  min5 = (j3 - j1 - m2) / 2 + kmin;


  /* sum series in double precision */
  double uk = 1.0e-10;
  double s = 1.0e-10;
  int ncut = 0;
  kmax -= kmin;
  if (kmax > 0) {
    for (int k = 1; k <= kmax; k++) {
      uk = -uk * (static_cast<double>((mini - k) * (min2 - k) * (min3 - k)))
        / (static_cast<double>((kmin + k) * (min4 + k) * (min5 + k)));
      if (fabs(uk) >= 1.0e30) {
        uk *= 1.0e-10;
        s *= 1.0e-10;
        ncut++;
      }
      if (fabs(uk) < 1.0e-20) break;
      s += uk;
    }
  }

  /* calculate delta functions */
  double delog = 0.0;
  for (unsigned int i = 1; i < 10; i++) delog += fl[mtri[i] + 1];

  int num = (j1 + j2 + j3) / 2 + 2;
  delog = (delog - fl[num]) / 2;
  double ulog = -fl[kmin + 1] - fl[mini] - fl[min2] - fl[min3]
    - fl[min4 + 1] - fl[min5 + 1];
  double plog = delog + ulog;
  double r;

  if ((plog < -80.0) || (ncut > 0)) {
    double sig = (s > 0) ? 1 : ((s < 0) ? -1 : 0);
    s = fabs(s);
    double slog = log(s) + (static_cast<double>(ncut + 1)) * log(1.0e10);
    r = sig * exp(slog + plog);
  } else {
    s *= 1.0e10;
    r = exp(plog) * s;
  }
  num = kmin + (j1 - j2 - m3) / 2;
  if ((num % 2) != 0)
    r = -r;

  return r;
}

bool ClebschGordan::MayBranch(int j1, int j2, int j3) {
  for (int m1 = -j1; m1 <= j1; m1 += 2) {
    for (int m2 = -j2; m2 <= j2; m2 += 2) {
      if (this->operator()(j1, j2, j3, m1, m2, m1+m2) > 0)
        return true;
    }
  }
  return false;
}
