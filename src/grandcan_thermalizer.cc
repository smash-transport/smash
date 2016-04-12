/*
 *    Copyright (c) 2016-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/grandcan_thermalizer.h"

#include "include/cxx14compat.h"
#include "include/particles.h"

namespace Smash {

GrandCanThermalizer::GrandCanThermalizer(const std::array<float, 3> cell_sizes,
                                         const std::array<int, 3> n_cells,
                                         float e_critical) :
  e_crit_(e_critical) {
  const std::array<float, 3> l = {cell_sizes[0]*n_cells[0],
                                  cell_sizes[1]*n_cells[1],
                                  cell_sizes[2]*n_cells[2]};
  /* Lattice is placed such that the center is 0,0,0.
     If one wants to have a central cell with center at 0,0,0 then
     number of cells should be odd (2k+1) in each direction.
   */
  const std::array<float, 3> origin = {0.5f*l[0], 0.5f*l[1], 0.5f*l[2]};
  const bool periodicity = false;
  const LatticeUpdate upd = LatticeUpdate::EveryFixedInterval;
  lat_ = make_unique<RectangularLattice<ThermLatticeNode>>(l,
                                                           n_cells,
                                                           origin,
                                                           periodicity,
                                                           upd);
}

void GrandCanThermalizer::update_lattice(const Particles& particles,
                                         const DensityParameters& dens_par) {
  const DensityType dens_type = DensityType::Hadron;
  // ToDo(oliiny): fix this stuff with LatticeUpdate
  const LatticeUpdate update = LatticeUpdate::EveryFixedInterval;
  update_general_lattice(lat_.get(), update, dens_type, dens_par, particles);
}

void GrandCanThermalizer::thermalize(Particles& particles) {
  // 1. Loop over particles, remove those which lie in the cells with e > e_crit_
  // 2. Loop over cells, sample particles according isochronous Cooper-Frye

  ThermLatticeNode node;
  for (auto &particle : particles) {
    const bool is_on_lattice = lat_->value_at(particle.position().threevec(),
                                              node);
    if (is_on_lattice && node.e() > e_crit_) {
      particles.remove(particle);
    }
  }

}

ThermLatticeNode::ThermLatticeNode() :
  Tmu0_(FourVector()),
  nb_(0.0),
  ns_(0.0) {}

void ThermLatticeNode::add_particle(const ParticleData& part, double factor) {
  Tmu0_ += part.momentum() * factor;
  nb_ += static_cast<double>(part.type().baryon_number()) * factor;
  ns_ += static_cast<double>(part.type().strangeness()) * factor;
}

void ThermLatticeNode::compute_rest_frame_quantities(HadronGasEos& eos) {
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
    const std::array<double, 3> T_mub_mus = eos.solve_eos(e_,
                                                 gamma_inv*nb_, gamma_inv*ns_);
    T_ = T_mub_mus[0];
    mub_ = T_mub_mus[1];
    mus_ = T_mub_mus[2];
    p_ = HadronGasEos::pressure(T_, mub_, mus_);
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
