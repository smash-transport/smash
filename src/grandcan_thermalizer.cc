/*
 *    Copyright (c) 2016-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/grandcan_thermalizer.h"

namespace Smash {

ThermLatticeNode::ThermLatticeNode() :
  Tmu0_(FourVector()),
  nb_(0.0),
  ns_(0.0) {}

void ThermLatticeNode::add_particle(const ParticleData& part, double factor) {
  Tmu0_ += part.momentum() * factor;
  nb_ += static_cast<double>(part.type().baryon_number()) * factor;
  ns_ += static_cast<double>(part.type().strangeness()) * factor;
}

void ThermLatticeNode::compute_rest_frame_quantities(HadgasEos& eos) {
  const int max_iter = 1000;
  v_ = ThreeVector(0.0, 0.0, 0.0);
  double e_previous_step = 0.0;
  const double tolerance = 1.e-8;
  for(int iter = 0; iter < max_iter; iter++) {
    e_previous_step = e_;
    e_ = Tmu0_.x0() - Tmu0_.threevec()*v_;
    if (std::abs(e_ - e_previous_step) < tolerance) {
      break;
    }
    const double gamma_inv = std::sqrt(1.0 - v_.sqr());
    const std::array<double, 3> T_mub_mus = eos.solve_hadgas_eos(e_,
                                                 gamma_inv*nb_, gamma_inv*ns_);
    T_ = T_mub_mus[0];
    mub_ = T_mub_mus[1];
    mus_ = T_mub_mus[2];
    p_ = HadgasEos::hadgas_pressure(T_, mub_, mus_);
    v_ = Tmu0_.threevec()/(Tmu0_.x0() + p_);
  }
}

std::ostream &operator<<(std::ostream &out, const ThermLatticeNode &node) {
  return out
         << "T[mu,0]: " << node.Tmu0()
         << ", nb: " << node.nb()
         << ", ns: " << node.ns()
         << ", v: " << node.v()
         << ", e: " << node.e()
         << ", p: " << node.p()
         << ", T: " << node.T()
         << ", mub: " << node.mub()
         << ", mus: " << node.mus();
}

}  // namespace Smash
