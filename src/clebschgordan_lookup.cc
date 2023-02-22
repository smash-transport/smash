/*
 *    Copyright (c) 2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/clebschgordan_lookup.h"

#include <numeric>
#include <unordered_map>

#include "gsl/gsl_sf_coupling.h"

#include "smash/constants.h"
#include "smash/iomanipulators.h"
#include "smash/logging.h"

namespace smash {
static constexpr int LResonances = LogArea::Resonances::id;

double ClebschGordan::calculate_coefficient(const int j_a, const int j_b,
                                            const int j_c, const int m_a,
                                            const int m_b, const int m_c) {
  const double wigner_3j = gsl_sf_coupling_3j(j_a, j_b, j_c, m_a, m_b, -m_c);
  if (std::abs(wigner_3j) < really_small) {
    return 0.;
  }
  assert((j_a - j_b + m_c) % 2 == 0);
  const int j = (j_a - j_b + m_c) / 2;
  double result = std::sqrt(j_c + 1) * wigner_3j;
  result *= (j % 2 == 0) * 2 - 1;  // == (-1)**j

  logg[LResonances].debug("CG: ", result, " I1: ", j_a, " I2: ", j_b,
                          " IR: ", j_c, " iz1: ", m_a, " iz2: ", m_b,
                          " izR: ", m_c);

  return result;
}

double ClebschGordan::coefficient(const int j_a, const int j_b, const int j_c,
                                  const int m_a, const int m_b, const int m_c) {
  const ThreeSpins spin_information = {j_a, j_b, j_c, m_a, m_b, m_c};
  if (auto search = lookup_table.find(spin_information);
      search != lookup_table.end()) {
    return search->second;
  } else {
    double result = calculate_coefficient(j_a, j_b, j_c, m_a, m_b, m_c);
    lookup_table.insert({spin_information, result});
    return result;
  }
}

}  // namespace smash
