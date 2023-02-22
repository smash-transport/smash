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
using namespace smash::clebsch_gordan;

[[maybe_unused]] static std::ostream &operator<<(std::ostream &out, const ThreeSpins &v) {
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
 */
TEST(tabulate) {
  std::vector<ThreeSpins> keys{};
  std::vector<double> coefficients{};
  const int L = 3;
  for (auto l1 = 0; l1 <= L; l1++) {
    for (auto l2 = 0; l2 <= L; l2++) {
      for (auto l3 = std::abs(l1 - l2); l3 <= l1 + l2; l3++) {
        for (auto m1 = -l1; m1 <= l1; m1++) {
          for (auto m2 = -l2; m2 <= l2; m2++) {
            for (auto m3 = -l3; m3 <= l3; m3++) {
              coefficients.push_back(
                  clebsch_gordan_calculation(l1, l2, l3, m1, m2, m3));
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
  VERIFY(keys.size() == 189);
  VERIFY(coefficients.size() == 189);
  /*
   * The printing is commented out, in order not to pollute normal runs of tests.
   * Uncomment and run in case of need. Printing to standard error facilitates
   * redirection to file, especially in case L should be increased in the future.
   */
  /*
  std::cerr << "\nstd::unordered_map<ThreeSpins, double, "
               "smash::ThreeSpinHash> CG_coefficients =\n{\n";
  const auto double_digits = std::numeric_limits<double>::max_digits10;
  for (auto i = 0u; i < keys.size(); i++) {
    std::cerr << "  {" << keys[i] << ", " << std::setprecision(double_digits)
              << std::setw(double_digits + 3) << coefficients[i] << "},\n";
  }
  std::cerr << "};\n\n";
  */
}
