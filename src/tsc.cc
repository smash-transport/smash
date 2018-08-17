/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/tsc.h"

#include <cmath>
#include <iomanip>
#include <iostream>

namespace smash {

std::ostream &operator<<(std::ostream &out, const TimeStampCounter &tsc) {
  auto c = tsc.cycles();
  int blocks[10];
  int n = 0;
  for (int digits = std::log10(c); digits > 0; digits -= 3) {
    blocks[n++] = c % 1000;
    c /= 1000;
  }
  if (n == 0) {
    return out;
  }
  const auto lastFill = out.fill('0');
  out << blocks[--n];
  while (n > 0) {
    out << '\'' << std::setw(3) << blocks[--n];
  }
  out.fill(lastFill);
  return out << " Cycles";
}

}  // namespace smash
