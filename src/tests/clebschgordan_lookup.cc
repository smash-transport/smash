/*
 *
 *    Copyright (c) 2013-2018,2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/clebschgordan_lookup.h"

#include <vector>

#include "setup.h"
#include "smash/clebschgordan.h"

using namespace smash;

TEST(coefficient) {
  /* spins are two times the actual values,
   * so j = 1 for spin-1/2 particle, 2 for spin-1 particle, etc.
   * Ordering of spins in array:
   *  0: j1, 1: j2, 2: j3, 3: m1, 4: m2, 5: m3
   */
  int spin[7][3];
  int spinz[7][3];
  double correct_coefficient[7];

  spin[0][0] = 1;
  spin[0][1] = 1;
  spin[0][2] = 2;
  spinz[0][0] = 1;
  spinz[0][1] = 1;
  spinz[0][2] = 2;
  correct_coefficient[0] = 1.;

  spin[1][0] = 1;
  spin[1][1] = 1;
  spin[1][2] = 2;
  spinz[1][0] = 1;
  spinz[1][1] = -1;
  spinz[1][2] = 0;
  correct_coefficient[1] = 1 / std::sqrt(2.);

  spin[2][0] = 2;
  spin[2][1] = 1;
  spin[2][2] = 1;
  spinz[2][0] = 2;
  spinz[2][1] = -1;
  spinz[2][2] = 1;
  correct_coefficient[2] = std::sqrt(2. / 3.);

  spin[3][0] = 2;
  spin[3][1] = 1;
  spin[3][2] = 3;
  spinz[3][0] = -2;
  spinz[3][1] = 1;
  spinz[3][2] = -1;
  correct_coefficient[3] = std::sqrt(1. / 3.);

  spin[4][0] = 2;
  spin[4][1] = 2;
  spin[4][2] = 2;
  spinz[4][0] = 0;
  spinz[4][1] = 2;
  spinz[4][2] = 2;
  correct_coefficient[4] = -1 / std::sqrt(2.);

  spin[5][0] = 2;
  spin[5][1] = 2;
  spin[5][2] = 2;
  spinz[5][0] = 0;
  spinz[5][1] = 0;
  spinz[5][2] = 0;
  correct_coefficient[5] = 0.;

  spin[6][0] = 2;
  spin[6][1] = 2;
  spin[6][2] = 4;
  spinz[6][0] = 2;
  spinz[6][1] = -2;
  spinz[6][2] = 0;
  correct_coefficient[6] = 1 / std::sqrt(6.);
  for (int i = 0; i < 7; i++) {
    double cg =
        ClebschGordan::coefficient(spin[i][0], spin[i][1], spin[i][2],
                                   spinz[i][0], spinz[i][1], spinz[i][2]);
    FUZZY_COMPARE(cg, correct_coefficient[i])
        << '\n'  // Using double quotes here produces an error(?!)
        << "J1: " << spin[i][0] << " Jz1: " << spinz[i][0] << "\n"
        << "J2: " << spin[i][1] << " Jz2: " << spinz[i][1] << "\n"
        << "J3: " << spin[i][2] << " Jz3: " << spinz[i][2] << "\n"
        << "CG: " << cg << " Correct: " << correct_coefficient[i];
  }
}

[[maybe_unused]] static std::ostream &operator<<(
    std::ostream &out, const ClebschGordan::ThreeSpins &v) {
  out.put('{');
  out << v.j1 << ", " << v.j2 << ", " << v.j3 << ',' << std::showpos
      << field<2> << v.m1 << ',' << field<2> << v.m2 << ',' << field<2> << v.m3
      << std::noshowpos;
  return out << '}';
}

/*
 * This test is not really supposed to be a unit test, but rather a way to
 * include somewhere the calculation of the tabulated Clebsch-Gordan
 * coefficients. Here we use a couple of vectors in order to then print
 * the coefficients ordered w.r.t. to j1, j2, j3, m1, m2, m3, respectively.
 * Note that physics symmetries and properties could be used to reduce the
 * number of iterations, since j1+j2+j3 should be even (here all j and m
 * are twice the physical quantity, so even means that its half is integer),
 * m1+m2=m3 (since m3 here enters the CG formula as -m3) and |m3|≤j3. Instead,
 * these properties are used to test the obtained non-zero results.
 */
TEST(tabulate) {
  std::vector<ClebschGordan::ThreeSpins> keys{};
  std::vector<double> coefficients{};
  const int L = 3;
  for (auto l1 = 0; l1 <= L; l1++) {
    for (auto l2 = 0; l2 <= L; l2++) {
      for (auto l3 = std::abs(l1 - l2); l3 <= l1 + l2; l3++) {
        for (auto m1 = -l1; m1 <= l1; m1 += 2) {
          for (auto m2 = -l2; m2 <= l2; m2 += 2) {
            for (auto m3 = -l3; m3 <= l3; m3 += 2) {
              coefficients.push_back(
                  ClebschGordan::coefficient(l1, l2, l3, m1, m2, m3));
              if (coefficients.back() != 0.0)
                keys.push_back({l1, l2, l3, m1, m2, m3});
              else
                coefficients.pop_back();
            }
          }
        }
      }
    }
  }
  for (auto i = 0u; i < keys.size(); i++) {
    VERIFY((keys[i].j1 + keys[i].j2 + keys[i].j3) % 2 == 0);
    VERIFY(keys[i].m1 + keys[i].m2 == keys[i].m3);
  }
  /*
   * The printing is left out, in order not to pollute normal runs of tests.
   * Change the #if 0 to e.g. #if 1 to include this snippet in compilation in
   * case of need. Printing to standard error facilitates redirection to file,
   * especially in case L should be increased in the future (it is then easier
   * to copy from a file than from the terminal).
   */
#if 0
  std::cerr << "\nstd::unordered_map<ThreeSpins, double, "
               "smash::ThreeSpinHash> CG_coefficients =\n{\n";
  const auto double_digits = std::numeric_limits<double>::max_digits10;
  for (auto i = 0u; i < keys.size(); i++) {
    std::cerr << "  {" << keys[i] << ", " << std::setprecision(double_digits)
              << std::setw(double_digits + 3) << coefficients[i] << "},\n";
  }
  std::cerr << "};\n\n";
#endif
}
