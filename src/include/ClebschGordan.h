/*
 * Class for calculating Clebsch-Gordan coefficients
 * Based on Caswell, R. S. and  Maximon, L. C.,
 * "Fortran programs for the calculation of Wigner 3j, 6j, and 9j
 * coefficients for angular momenta greater than or equal to 80",
 * NBS Technical Note 409 (1966).
 * http://archive.org/details/fortranprogramsf409casw
 *
 * Copyright (c) 2013
 *      SMASH Team
 * 
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_CLEBSCHGORDAN_H_
#define SRC_INCLUDE_CLEBSCHGORDAN_H_

class ClebschGordan {
 protected:
  unsigned int factorial_log_max;
  double * factorial_log;
  /* Wigner 3j symbol <j1 m1 j2 m2 | j3 m3> */
  double f3j(int j1, int j2, int j3, int m1, int m2, int m3);

 public:
  ClebschGordan();
  ClebschGordan(const ClebschGordan& original);
  ~ClebschGordan();

  /* Calculate the Clebsch-Gordan coefficient */
  double operator()(int j1, int j2, int j3, int m1, int m2, int m3);
  /* Check if it's possible to combine j1 and j2 to j3 */
  bool MayBranch(int j1, int j2, int j3);
};

#endif  // SRC_INCLUDE_CLEBSCHGORDAN_H_
