/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"

#include "../include/average.h"

TEST(average) {
    Average<double> avg;
    for (int i = 1; i <= 5; i++)
        avg.add(i);
    COMPARE(avg.average(), 3.0);
}

void test_dedup_avg(
        std::initializer_list<double>&& a,
        std::initializer_list<double>&& b,
        std::initializer_list<double>&& expected_a,
        std::initializer_list<double>&& expected_b) {
    const std::vector<double> x = a;
    const std::vector<double> y = b;
    std::vector<double> dedup_x;
    std::vector<double> dedup_y;
    std::tie(dedup_x, dedup_y) = dedup_avg(x, y);
    const std::vector<double> expected_x = expected_a;
    const std::vector<double> expected_y = expected_b;
    COMPARE(dedup_x.size(), dedup_y.size());
    COMPARE(dedup_x.size(), expected_x.size());
    COMPARE(dedup_y.size(), expected_y.size());
    for (size_t i = 0; i < dedup_x.size(); i++) {
        COMPARE(dedup_x[i], expected_x[i]);
        COMPARE(dedup_y[i], expected_y[i]);
    }
}

TEST(dedup_avg) {
    test_dedup_avg({1, 2, 2, 4, 5, 6, 7, 7, 9}, {1, 2, 3, 4, 5, 6, 7, 7, 9},
                   {1, 2, 4, 5, 6, 7, 9}, {1, 2.5, 4, 5, 6, 7, 9});
    test_dedup_avg({1, 2, 2, 4, 5, 6, 7, 7}, {1, 2, 3, 4, 5, 6, 7, 7},
                   {1, 2, 4, 5, 6, 7}, {1, 2.5, 4, 5, 6, 7});
    test_dedup_avg({1}, {2}, {1}, {2});
}
