/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/pauliblocking.h"
#include "smash/constants.h"
#include "smash/logging.h"

namespace smash {

PauliBlocker::PauliBlocker(Configuration conf,
                           const ExperimentParameters &param)
    : sig_(param.gaussian_sigma),
      rc_(conf.take({"Gaussian_Cutoff"}, 2.2)),
      rr_(conf.take({"Spatial_Averaging_Radius"}, 1.86)),
      rp_(conf.take({"Momentum_Averaging_Radius"}, 0.08)),
      ntest_(param.testparticles) {
  /*!\Userguide
   * \page pauliblocker Pauli_Blocking
   *
   * \key Spatial_Averaging_Radius (double, optional, default = 1.86): \n
   * Radius [fm] of sphere for averaging in the coordinate space
   *
   * \key Gaussian_Cutoff (double, optional, default = 2.2): \n
   * Radius [fm] at which gaussians used for smoothing are cut
   *
   * \key Momentum_Averaging_Radius (double, optional, default = 0.08): \n
   * Radius [GeV/c] of sphere for averaging in the momentum space
   */

  const auto &log = logger<LogArea::PauliBlocking>();

  if (ntest_ < 20) {
    log.warn(
        "Phase-space density calculation in Pauli blocking"
        " will not work reasonably for a small number of testparticles."
        " The recommended number of testparticles is 20.");
  }

  if (rc_ < rr_ || rr_ < 0.0 || rp_ < 0) {
    log.error(
        "Please choose reasonable parameters for Pauli blocking:"
        "All radii have to be positive and Gaussian_Cutoff should"
        "be larger than Spatial_Averaging_Radius");
  }

  init_weights_analytical();
}

PauliBlocker::~PauliBlocker() {}

double PauliBlocker::phasespace_dens(const ThreeVector &r, const ThreeVector &p,
                                     const Particles &particles,
                                     const PdgCode pdg,
                                     const ParticleList &disregard) const {
  double f = 0.0;

  /* TODO(oliiny): looping over all particles is inefficient,
   * I need only particles within rp_ radius in momentum and
   * within rr_+rc_ in coordinate space. Some search algorithm might help. */
  for (const auto &part : particles) {
    // Only consider identical particles
    if (part.pdgcode() != pdg) {
      continue;
    }
    // Only consider momenta in sphere of radius rp_ with center at p
    const double pdist_sqr = (part.momentum().threevec() - p).sqr();
    if (pdist_sqr > rp_ * rp_) {
      continue;
    }
    const double rdist_sqr = (part.position().threevec() - r).sqr();
    // Only consider coordinates in sphere of radius rr_+rc_ with center at r
    if (rdist_sqr >= (rr_ + rc_) * (rr_ + rc_)) {
      continue;
    }
    // Do not count particles that should be disregarded.
    bool to_disregard = false;
    for (const auto &disregard_part : disregard) {
      if (part.id() == disregard_part.id()) {
        to_disregard = true;
      }
    }
    if (to_disregard) {
      continue;
    }
    // 1st order interpolation using tabulated values
    const double i_real = std::sqrt(rdist_sqr) / (rr_ + rc_) * weights_.size();
    const size_t i = std::floor(i_real);
    const double rest = i_real - i;
    if (likely(i + 1 < weights_.size())) {
      f += weights_[i] * rest + weights_[i + 1] * (1. - rest);
    }
  }
  return f / ntest_;
}

void PauliBlocker::init_weights_analytical() {
  const auto &log = logger<LogArea::PauliBlocking>();

  const double pi = M_PI;
  const double sqrt2 = std::sqrt(2.);
  const double sqrt_2pi = std::sqrt(2. * pi);
  // Volume of the phase-space area; Factor 2 stands for spin.
  const double phase_volume =
      2 * (4. / 3. * pi * rr_ * rr_ * rr_) * (4. / 3. * pi * rp_ * rp_ * rp_) /
      ((2 * pi * hbarc) * (2 * pi * hbarc) * (2 * pi * hbarc));
  // Analytical expression for integral in denominator
  const double norm =
      std::erf(rc_ / sqrt2 / sig_) -
      rc_ * 2 / sqrt_2pi / sig_ * std::exp(-0.5 * rc_ * rc_ / sig_ / sig_);

  double integral;
  // Step of the table for tabulated integral
  const double d_pos = (rr_ + rc_) / static_cast<double>(weights_.size());

  for (size_t k = 0; k < weights_.size(); k++) {
    // rdist = 0 ... rc_ (gauss cut) + rr_ (position cut)
    const double rj = d_pos * k;
    if (rj < really_small) {
      // Assuming rc_ > rr_
      const double A = rr_ / sqrt2 / sig_;
      integral = sqrt_2pi * sig_ * std::erf(A) - 2 * rr_ * std::exp(-A * A);
      integral *= sig_ * sig_;
    } else if (rc_ > rj + rr_) {
      const double A = (rj + rr_) / sqrt2 / sig_;
      const double B = (rj - rr_) / sqrt2 / sig_;
      integral = sig_ / rj * (std::exp(-A * A) - std::exp(-B * B)) +
                 0.5 * sqrt_2pi * (std::erf(A) - std::erf(B));
      integral *= sig_ * sig_ * sig_;
    } else {
      const double A = rc_ / sqrt2 / sig_;
      const double B = (rj - rr_) / sqrt2 / sig_;
      const double C = (rc_ - rj) * (rc_ - rj) - rr_ * rr_ + 2 * sig_ * sig_;
      integral =
          (0.5 * std::exp(-A * A) * C - sig_ * sig_ * std::exp(-B * B)) / rj +
          0.5 * sqrt_2pi * sig_ * (std::erf(A) - std::erf(B));
      integral *= sig_ * sig_;
    }
    integral *= 2 * pi / std::pow(2 * pi * sig_ * sig_, 1.5);
    weights_[k] = integral / norm / phase_volume;
    log.debug("Analytical weights[", k, "] = ", weights_[k]);
  }
}

}  // namespace smash
