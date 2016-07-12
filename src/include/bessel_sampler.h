/*
 *
 *    Copyright (c) 2016-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_BESSEL_SAMPLER_H_
#define SRC_INCLUDE_BESSEL_SAMPLER_H_

#include "random.h"

namespace Smash {

class BesselSampler {

 public:
  BesselSampler(const double poisson_mean1,
                const double poisson_mean2,
                const int fixed_difference) :
    a_(2.0*std::sqrt(poisson_mean1*poisson_mean2)),
    N_(fixed_difference) {
    m_ = 0.5 * (std::sqrt(a_*a_ + N_*N_) - N_);
    if (m_ >= m_switch_method_) {
      std::cout << "Using gaussian method" << std::endl;
      mu_ = 0.5 * a_ * r(N_, a_);
      std::cout << "mean = " << mu_ << std::endl;
      const double mean_sqr = mu_* (1.0 + 0.5 * a_ * r(N_+1, a_));
      sigma_ = std::sqrt(mean_sqr - mu_*mu_);
      std::cout << "sigma = " << sigma_ << std::endl;
    } else {
      std::cout << "Using direct sampling method" << std::endl;
      std::vector<double> probabilities;
      double wi = 1.0, sum = 0.0;
      int i = 0;
      do {
        sum += wi;
        probabilities.push_back(wi);
        wi *= 0.25*a_*a_ / (i + 1) / (N_ + i + 1);
        i++;
      } while (wi > negligible_probability_);
      i = 0;
      for (double p : probabilities) {
        p /= sum;
        std::cout << "probability (" << i << ") = " << p << std::endl;
        i++;
      }
      dist_.reset_weights(probabilities);
    }
  }

  int sample() {
    return (m_ >= m_switch_method_) ?
             std::round(Random::normal(mu_, sigma_)) :
             dist_();
  };

 private:
   double r(int n, double a) {
     // gsl_sf_bessel_In_scaled(n+1,a)/gsl_sf_bessel_In_scaled(n,a);
     const double a_inv = 1.0/a;
     double res = 0.0;
     int i = 200;
     for (; i > 0; i--) {
       res = 1.0 / (a_inv * 2 * (n + i) + res);
     }
     std::cout << a/(std::sqrt(a*a + (n+1)*(n+1)) + n + 1) << " <= " << res << " <= " <<
                  a/(std::sqrt(a*a + n*n) + n) << std::endl;
     return res;
   }
   Random::discrete_dist<double> dist_;
   const double a_;
   const int N_;
   static constexpr double m_switch_method_ = 6.0;
   static constexpr double negligible_probability_ = 1.e-12;
   double m_;
   double mu_;
   double sigma_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_BESSEL_SAMPLER_H_
