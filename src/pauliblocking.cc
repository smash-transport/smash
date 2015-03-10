/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
*/

#include "include/constants.h"
#include "include/logging.h"
#include "include/pauliblocking.h"

namespace Smash {

/**
 * PauliBlocker constructor. Gets parameters from configuration.
 * Tabulates necessary integrals.
 */
PauliBlocker::PauliBlocker(Configuration conf, const ExperimentParameters &param)
    : sig_(param.gaussian_sigma),
      rc_(conf.take({"Gaussian_Cutoff"}, 2.2)),
      rr_(conf.take({"Spatial_Averaging_Radius"}, 1.86)),
      rp_(conf.take({"Momentum_Averaging_Radius"}, 0.08)),
      ntest_(param.testparticles) {

  /*!\Userguide
   * \page potentials PauliBlocker
   *
   * \key Spatial_Averaging_Radius (float, optional, default = 1.86): \n
   * Radius [fm] of sphere for averaging in the coordinate space
   *
   * \key Gaussian_Cutoff (float, optional, default = 2.2): \n
   * Radius [fm] at which gaussians used for smoothing are cut
   *
   * \key Momentum_Averaging_Radius (float, optional, default = 0.08): \n
   * Radius [GeV/c] of sphere for averaging in the momentum space
   */

  const auto &log = logger<LogArea::PauliBlocking>();

  if (ntest_ < 20) {
    log.error("Phase-space density calculation in Pauli blocking"
              " will not work reasonably for small number of testparticles."
              " Recommended number of testparticles is 200.");
  }

  if (rc_ < rr_ || rr_ < 0.0 || rp_ < 0) {
    log.error("Please choose reasonable parameters for Pauli blocking:"
              "All radii have to be positive and Gaussian_Cutoff should"
              "be larger than Spatial_Averaging_Radius");
  }

  init_weights_analytical();
}

PauliBlocker::~PauliBlocker() {
}

float PauliBlocker::phasespace_dens(const ThreeVector r, const ThreeVector p,
                        const Particles &particles, const PdgCode pdg) const {

  float f = 0.0;
  float rdist_sqr, pdist_sqr;
  int index;

  for (const auto &part : particles.data()) {
    // Only consider identical particles
    if (part.pdgcode() != pdg) {
      continue;
    }
    // Only consider momenta in sphere of radius rp_ with center at p
    pdist_sqr = (part.momentum().threevec() - p).sqr();
    if (pdist_sqr > rp_*rp_) {
      continue;
    }
    rdist_sqr = (part.position().threevec() - r).sqr();
    // Only consider coordinates in sphere of radius rr_+rc_ with center at r
    if (rdist_sqr > (rr_+rc_)*(rr_+rc_)) {
      continue;
    }
    // 0th order interpolation using tabulated values
    index = std::round(std::sqrt(rdist_sqr) / (rr_ + rc_) * weights_.size());
    f += weights_[index];
  }
  return f / ntest_;
};

void PauliBlocker::init_weights() {
  const auto &log = logger<LogArea::PauliBlocking>();

  // Number of integration steps for radial integration
  const int steps_r = 2000;
  // Number of integration steps for theta integration
  const int steps_theta = 2000;
  // Volume of the phase-space area; Factor 2 stands for spin.
  const float phase_volume = 2. * (4./3.*M_PI*rr_*rr_*rr_) *
                                  (4./3.*M_PI*rp_*rp_*rp_) /
                                  ((2*M_PI*hbarc)*(2*M_PI*hbarc)*(2*M_PI*hbarc));
  // Analytical expression for integral in denominator
  const float norm = std::erf(rc_/std::sqrt(2)/sig_) -
                     rc_* std::sqrt(2./M_PI) / sig_ * std::exp(-rc_*rc_/2./sig_/sig_);

  // Define integration steps
  const float dr = rr_ / static_cast<float>(steps_r);
  const float d_cos_theta = 1./ static_cast<float>(steps_theta);
  // Step of the table for tabulated integral
  const float d_pos = (rr_ + rc_) / static_cast<float>(weights_.size());

  // Running variables of integration
  float rdist, r, cos_theta;
  // |vec(r) - vec(rdist)|^2
  float distance_sqr;
  // integration result
  float integral;

  for (size_t k = 0; k < weights_.size(); k++) {
    // rdist = 0 ... rc_ (gauss cut) + rr_ (position cut)
    rdist = d_pos * k;
    integral = 0.;
    for (auto n = -steps_theta; n < steps_theta; n++) {
      // integral over cos(theta) = -1 ... 1
      cos_theta = (n + 0.5) * d_cos_theta;
      for (auto m = 0; m < steps_r; m++) {
        // integral over r = 0 ... rr_ (position cut)
        r = (m + 0.5) * dr;
        distance_sqr = r*r + rdist*rdist - 2.*r*rdist*cos_theta;
        if (distance_sqr < rc_*rc_) {
          integral += r*r * std::exp(- distance_sqr/2./sig_/sig_);
        }
      }
    }
    integral *= 2.*M_PI*d_cos_theta*dr / std::pow(2.*M_PI*sig_*sig_, 3./2.);
    weights_[k] = integral / norm / phase_volume;
    log.debug("Numerical weights[",k,"] = ", weights_[k]);
  }
}

void PauliBlocker::init_weights_analytical() {
  const auto &log = logger<LogArea::PauliBlocking>();

  // Volume of the phase-space area; Factor 2 stands for spin.
  const float phase_volume = 2. * (4./3.*M_PI*rr_*rr_*rr_) *
                                  (4./3.*M_PI*rp_*rp_*rp_) /
                                  ((2*M_PI*hbarc)*(2*M_PI*hbarc)*(2*M_PI*hbarc));
  // Analytical expression for integral in denominator
  const float norm = std::erf(rc_/std::sqrt(2)/sig_) -
                     rc_* std::sqrt(2./M_PI) / sig_ * std::exp(-rc_*rc_/2./sig_/sig_);

  float rj, integral;
  // Step of the table for tabulated integral
  const float d_pos = (rr_ + rc_) / static_cast<float>(weights_.size());

  for (size_t k = 0; k < weights_.size(); k++) {
    // rdist = 0 ... rc_ (gauss cut) + rr_ (position cut)
    rj = d_pos * k;
    if (rj < really_small) {
      // Assuming rc_ > rr_
      const float A = rr_/std::sqrt(2.0)/sig_;
      integral = std::sqrt(2.0*M_PI)*sig_*std::erf(A) - 2.0*rr_*std::exp(-A*A);
      integral *= sig_*sig_;
    } else if (rc_ > rj + rr_) {
      const float A = (rj+rr_)/std::sqrt(2.0)/sig_;
      const float B = (rj-rr_)/std::sqrt(2.0)/sig_;
      integral = sig_ / rj * (std::exp(-A*A) - std::exp(-B*B)) +
                 std::sqrt(M_PI/2.0)*(std::erf(A) - std::erf(B));
      integral *= sig_*sig_*sig_;
    } else {
      const float A = rc_/std::sqrt(2.0)/sig_;
      const float B = (rj-rr_)/std::sqrt(2.0)/sig_;
      const float C = (rc_-rj)*(rc_-rj) - rr_*rr_ + 2.0*sig_*sig_;
      integral = (0.5*std::exp(-A*A)*C - sig_*sig_*std::exp(-B*B)) / rj +
                 std::sqrt(M_PI/2.0)*sig_*(std::erf(A) - std::erf(B));
      integral *= sig_*sig_;
    }
    integral *= 2.*M_PI / std::pow(2.*M_PI*sig_*sig_, 3./2.);
    weights_[k] = integral / norm / phase_volume;
    log.debug("Analytical weights[",k,"] = ", weights_[k]);
  }
}

}  // namespace Smash
