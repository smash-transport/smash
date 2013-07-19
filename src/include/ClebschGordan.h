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

#ifndef SRC_INCLUDE_CLEBSCHGORDAN_H_
#define SRC_INCLUDE_CLEBSCHGORDAN_H_

class ClebschGordan {
 protected:
  unsigned int flmax;
  double * fl;
  double f3j(int j1, int j2, int j3, int m1, int m2, int m3);

 public:
  ClebschGordan();
  ~ClebschGordan();

  double operator()(int j1, int j2, int j3, int m1, int m2, int m3);
  bool MayBranch(int j1, int j2, int j3);
};

#endif  // SRC_INCLUDE_CLEBSCHGORDAN_H_
