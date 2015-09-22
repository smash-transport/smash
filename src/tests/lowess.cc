/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/lowess.h"

#include <numeric>

TEST(smooth_trivial) {
    std::vector<double> x(20);
    std::iota(x.begin(), x.end(), 0);
    std::vector<double> y(20);
    for (size_t i = 0; i < 10; i++) {
        y[i] = 0;
    }
    for (size_t i = 10; i < 20; i++) {
        y[i] = 1;
    }
    const auto result = Smash::smooth(x, y, 0.4, 3, 0);
    for (size_t i = 0; i < 20; i++) {
        COMPARE_ABSOLUTE_ERROR(result[i], y[i], 1e-7);
    }
}
