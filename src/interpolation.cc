/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/interpolation.h"

#include <cassert>
#include <cstddef>

InterpolateLinear::InterpolateLinear(double x0, double y0, double x1, double y1) {
    assert(x0 != x1);
    m = (y1 - y0) / (x1 - x0);
    n = y0 - m*x0;
}

double InterpolateLinear::operator()(double x) const {
    return m*x + n;
}

InterpolateData::InterpolateData(const std::vector<double>& x, const std::vector<double>& y) {
    assert(x.size() == y.size());
    const size_t n = x.size();
    this->x = x;
    f.reserve(n - 1);
    for (size_t i = 0; i < n - 1; i++) {
        f.emplace_back( InterpolateLinear(x[i], y[i], x[i + 1], y[i + 1]) );
    }
}

double InterpolateData::operator()(double x0) const {
    size_t i = 0;
    while (x[i] > x0 && i < x.size() - 2) {
        i++;
    }
    return f[i](x0);
}
