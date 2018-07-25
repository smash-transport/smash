/*
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 *    This code has been adapted from ROOT's TGraphSmooth.cxx.
 *
 *    Copyright (c) 2006 Rene Brun and Fons Rademakers.
 *    Copyright (c) 2001 Christian Stratowa
 *    Copyright (c) 1999-2001 Robert Gentleman, Ross Ihaka and the R
 *                            Development Core Team
 */
#ifndef SRC_INCLUDE_LOWESS_H_
#define SRC_INCLUDE_LOWESS_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace smash {

namespace lowess {
/**
 * Fit value at x[i]
 * Based on R function lowest: Translated to C++ by C. Stratowa
 * (R source file: lowess.c by R Development Core Team (C) 1999-2001)
 */
template <typename T>
void lowest(const T *x, const T *y, size_t n, T xs, T &ys, size_t nleft,
            size_t nright, T *w, bool userw, T *rw, bool &ok) {
  // indices start at 1
  x--;
  y--;
  w--;
  rw--;

  const auto range = x[n] - x[1];
  const auto h = std::max(xs - x[nleft], x[nright] - xs);
  const auto h9 = 0.999 * h;
  const auto h1 = 0.001 * h;

  // sum of weights
  T a = 0.;
  auto j = nleft;
  while (j <= n) {
    // compute weights (pick up all ties on right)
    w[j] = 0.;
    const auto r = std::abs(x[j] - xs);
    if (r <= h9) {
      if (r <= h1) {
        w[j] = 1.;
      } else {
        const auto d = (r / h) * (r / h) * (r / h);
        w[j] = (1. - d) * (1. - d) * (1. - d);
      }
      if (userw)
        w[j] *= rw[j];
      a += w[j];
    } else if (x[j] > xs) {
      break;
    }
    j += 1;
  }

  // rightmost pt (may be greater than nright because of ties)
  const auto nrt = j - 1;
  if (a <= 0.) {
    ok = false;
  } else {
    ok = true;
    // weighted least squares: make sum of w[j] == 1
    for (j = nleft; j <= nrt; j++)
      w[j] /= a;
    if (h > 0.) {
      a = 0.;
      // use linear fit weighted center of x values
      for (j = nleft; j <= nrt; j++)
        a += w[j] * x[j];
      auto b = xs - a;
      T c = 0.;
      for (j = nleft; j <= nrt; j++)
        c += w[j] * (x[j] - a) * (x[j] - a);
      if (std::sqrt(c) > 0.001 * range) {
        b /= c;
        // points are spread out enough to compute slope
        for (j = nleft; j <= nrt; j++)
          w[j] *= (b * (x[j] - a) + 1.);
      }
    }
    ys = 0.;
    for (j = nleft; j <= nrt; j++)
      ys += w[j] * y[j];
  }
}

/**
 * Partial sort.
 * based on R function rPsort: adapted to C++ by Christian Stratowa
 * (R source file: R_sort.c by R Development Core Team (C) 1999-2001)
 */
template <typename T>
void psort(T *x, size_t n, size_t k) {
  for (size_t pL = 0, pR = n - 1; pL < pR;) {
    const auto v = x[k];
    size_t i;
    size_t j;
    for (i = pL, j = pR; i <= j;) {
      while (x[i] < v) {
        i++;
      }
      while (v < x[j]) {
        j--;
      }
      if (i <= j) {
        const auto w = x[i];
        x[i++] = x[j];
        x[j--] = w;
      }
    }
    if (j < k) {
      pL = i;
    }
    if (k < i) {
      pR = j;
    }
  }
}

/**
 * Lowess regression smoother.
 * Based on R function clowess: Translated to C++ by C. Stratowa
 * (R source file: lowess.c by R Development Core Team (C) 1999-2001)
 */
template <typename T>
void lowess(const T *x, const T *y, size_t n, T *ys, T span, size_t iter,
            T delta, T *rw, T *res) {
  if (n < 2) {
    ys[0] = y[0];
    return;
  }

  // nleft, nright, last, etc. must all be shifted to get rid of these:
  x--;
  y--;
  ys--;

  // at least two, at most n points
  constexpr size_t two = 2;
  const auto ns =
      std::max(two, std::min(n, static_cast<size_t>(span * n + 1e-7)));

  // robustness iterations
  size_t iiter = 1;
  while (iiter <= iter + 1) {
    size_t nleft = 1;
    size_t nright = ns;
    size_t last = 0;  // index of prev estimated point
    size_t i = 1;     // index of current point

    for (;;) {
      if (nright < n) {
        // move nleft,  nright to right if radius decreases
        const auto d1 = x[i] - x[nleft];
        const auto d2 = x[nright + 1] - x[i];

        // if d1 <= d2 with x[nright+1] == x[nright], lowest fixes
        if (d1 > d2) {
          // radius will not decrease by move right
          nleft++;
          nright++;
          continue;
        }
      }

      // fitted value at x[i]
      const bool iterg1 = iiter > 1;
      bool ok;
      lowest(&x[1], &y[1], n, x[i], ys[i], nleft, nright, res, iterg1, rw, ok);
      if (!ok)
        ys[i] = y[i];

      // all weights zero copy over value (all rw==0)
      if (last < i - 1) {
        const auto denom = x[i] - x[last];

        // skipped points -- interpolate non-zero - proof?
        for (auto j = last + 1; j < i; j++) {
          const auto alpha = (x[j] - x[last]) / denom;
          ys[j] = alpha * ys[i] + (1. - alpha) * ys[last];
        }
      }

      // last point actually estimated
      last = i;

      // x coord of close points
      const auto cut = x[last] + delta;
      for (i = last + 1; i <= n; i++) {
        if (x[i] > cut)
          break;
        if (x[i] == x[last]) {
          ys[i] = ys[last];
          last = i;
        }
      }
      i = std::max(last + 1, i - 1);
      if (last >= n)
        break;
    }

    // residuals
    for (i = 0; i < n; i++)
      res[i] = y[i + 1] - ys[i + 1];

    // compute robustness weights except last time
    if (iiter > iter)
      break;
    for (i = 0; i < n; i++)
      rw[i] = std::abs(res[i]);

    // compute cmad := 6 * median(rw[], n)
    const auto m1 = n / 2;
    // partial sort, for m1 & m2
    // TODO(steinberg): consider replacing psort with std::partial_sort
    psort(rw, n, m1);
    T cmad;
    if (n % 2 == 0) {
      const auto m2 = n - m1 - 1;
      psort(rw, n, m2);
      cmad = 3. * (rw[m1] + rw[m2]);
    } else { /* n odd */
      cmad = 6. * rw[m1];
    }

    const auto c9 = 0.999 * cmad;
    const auto c1 = 0.001 * cmad;
    for (i = 0; i < n; i++) {
      if (cmad == 0.) {
        // In this case, `r` cannot really be smaller than `c1` or `c2`, so we
        // would set `rw[i] = 0` anyway. To avoid divisions by zero, we do this
        // already here.
        rw[i] = 0;
        continue;
      }
      const auto r = std::abs(res[i]);

      if (r <= c1)
        rw[i] = 1.;
      else if (r <= c9)
        rw[i] = (1. - (r / cmad) * (r / cmad)) * (1. - (r / cmad) * (r / cmad));
      else
        rw[i] = 0.;
    }
    iiter++;
  }
}

}  // namespace lowess

/**
 * Apply the LOWESS smoother (see the reference below) to the given data
 * (x, y).
 *
 * \param x x-values.
 * \param y y-values.
 * \param span The smoother span. This gives the proportion of points in
 *     the plot which influence the smoothness at each value. Larger values
 *     give more smoothness.
 * \param iter The number of robustifying iterations which should be
 *     performed. Using smaller values of iter will make lowess run faster.
 * \param delta Values of x which lie within delta of each other replaced
 *     by a single value in the output from lowess.
 *     For delta = 0, delta will be calculated.
 * \return Smoothed y-values.
 *
 * References:
 *
 * - Cleveland, W. S. (1979) Robust locally weighted regression and smoothing
 *        scatterplots. J. Amer. Statist. Assoc. 74, 829-836.
 * - Cleveland, W. S. (1981) LOWESS: A program for smoothing scatterplots
 *        by robust locally weighted regression.
 *        The American Statistician, 35, 54.
 */
template <typename T>
std::vector<T> smooth(const std::vector<T> &x, const std::vector<T> &y,
                      T span = 2. / 3, size_t iter = 3, T delta = 0) {
  assert(x.size() == y.size());
  std::vector<T> result;
  result.resize(x.size());
  std::vector<T> rw;
  rw.resize(x.size());
  std::vector<T> res;
  res.resize(x.size());
  lowess::lowess(&x.front(), &y.front(), x.size(), &result.front(), span, iter,
                 delta, &rw.front(), &res.front());
  return std::move(result);
}

}  // namespace smash

#endif  // SRC_INCLUDE_LOWESS_H_
