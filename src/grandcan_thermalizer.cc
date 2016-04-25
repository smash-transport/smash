/*
 *    Copyright (c) 2016-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/grandcan_thermalizer.h"

#include "include/angles.h"
#include "include/cxx14compat.h"
#include "include/distributions.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"
#include "include/quantumnumbers.h"
#include "include/random.h"

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
  const std::array<float, 3> origin = {-0.5f*l[0], -0.5f*l[1], -0.5f*l[2]};
  const bool periodicity = false;
  const LatticeUpdate upd = LatticeUpdate::EveryFixedInterval;
  lat_ = make_unique<RectangularLattice<ThermLatticeNode>>(l,
                                                           n_cells,
                                                           origin,
                                                           periodicity,
                                                           upd);
  cell_volume_ = cell_sizes[0] * cell_sizes[1] * cell_sizes[2];
  cells_to_sample_.resize(50000);
}

ThreeVector GrandCanThermalizer::uniform_in_cell() const {
  return ThreeVector(Random::uniform(-0.5 * static_cast<double>(lat_->cell_sizes()[0]),
                       +0.5 * static_cast<double>(lat_->cell_sizes()[0])),
       Random::uniform(-0.5 * static_cast<double>(lat_->cell_sizes()[1]),
                       +0.5 * static_cast<double>(lat_->cell_sizes()[1])),
       Random::uniform(-0.5 * static_cast<double>(lat_->cell_sizes()[2]),
                       +0.5 * static_cast<double>(lat_->cell_sizes()[2])));
}

void GrandCanThermalizer::sample_in_random_cell(ParticleList& plist,
                                                const double time) {
  plist.clear();
  // Choose random cell
  int cells_to_sample_size = cells_to_sample_.size();
  const int cell_index =
              cells_to_sample_[Random::uniform_int(0, cells_to_sample_size-1)];
  const ThermLatticeNode cell = (*lat_)[cell_index];
  const ThreeVector cell_center = lat_->cell_center(cell_index);
  const double gamma = 1.0 / std::sqrt(1.0 - cell.v().sqr());
  // Loop over all existing hadrons (no leptons/quarks/etc)
  for (const ParticleType &ptype : ParticleType::list_all()) {
    if (!ptype.is_hadron()) {
      continue;
    }
    // First find out how many particles of this kind to sample
    // N = n u^mu dsigma_mu = (isochronous hypersurface) n * V * gamma
    const double N_average = cell_volume_ * gamma *
      HadronGasEos::partial_density(ptype, cell.T(), cell.mub(), cell.mus());
    const int N_to_sample = Random::poisson(N_average);
    // Sample particles: uniform in the cell, momenta from Cooper-Frye
    for (int i = 0; i < N_to_sample; i++) {
      ParticleData particle(ptype);
      // Note: it's pole mass for resonances!
      const double m = static_cast<double>(ptype.mass());
      // Position
      particle.set_4position(FourVector(time, cell_center + uniform_in_cell()));
      // Momentum
      double momentum_radial = sample_momenta_from_thermal(cell.T(), m);
      Angles phitheta;
      phitheta.distribute_isotropically();
      particle.set_4momentum(m, phitheta.threevec() * momentum_radial);
      particle.boost_momentum(cell.v());
      plist.push_back(particle);
    }
  }
}

void GrandCanThermalizer::update_lattice(const Particles& particles,
                                         const DensityParameters& dens_par) {
  const DensityType dens_type = DensityType::Hadron;
  // ToDo(oliiny): fix this stuff with LatticeUpdate
  const LatticeUpdate update = LatticeUpdate::EveryFixedInterval;
  update_general_lattice(lat_.get(), update, dens_type, dens_par, particles);
  for (auto &node : *lat_) {
    /* If energy density is definitely below e_crit -
       no need to find T, mu, etc. So if e = T00 - T0i*vi <=
       T00 + sum abs(T0i) < e_crit, no efforts are necessary. */
    if (node.Tmu0().x0() +
        std::abs(node.Tmu0().x1()) +
        std::abs(node.Tmu0().x2()) +
        std::abs(node.Tmu0().x3()) >= e_crit_) {
      node.compute_rest_frame_quantities(eos_);
    }
  }
}

void GrandCanThermalizer::thermalize(Particles& particles, double time) {
  // Remove particles from the cells with e > e_crit_,
  // sum up their conserved quantities
  QuantumNumbers conserved_initial   = QuantumNumbers(),
                 conserved_remaining = QuantumNumbers();
  ThermLatticeNode node;
  int particles_removed = 0;
  for (auto &particle : particles) {
    const bool is_on_lattice = lat_->value_at(particle.position().threevec(),
                                              node);
    if (is_on_lattice && node.e() > e_crit_) {
      conserved_initial.add_values(particle);
      particles.remove(particle);
      particles_removed ++;
    }
  }
  std::cout << "Removed " << particles_removed << " particles." << std::endl;

  // Exit if there is nothing to thermalize
  if (conserved_initial == QuantumNumbers()) {
    return;
  }

  // Save the indices of cells inside the volume with e > e_crit_
  const size_t lattice_total_cells = lat_->size();
  for (size_t i = 0; i < lattice_total_cells; i++) {
    if ((*lat_)[i].e() > e_crit_) {
      cells_to_sample_.push_back(i);
    }
  }

  ParticleList mode_list, sampled_list;
  double energy = 0.0;
  int S_plus = 0, S_minus = 0,
      B_plus = 0, B_minus = 0,
      E_plus = 0, E_minus = 0;
  // Mode 1: sample until energy is conserved, take only strangeness < 0
  while (conserved_initial.momentum().x0() > energy ||
         S_plus < conserved_initial.strangeness()) {
    sample_in_random_cell(mode_list, time);
    for (auto &particle : mode_list) {
      energy += particle.momentum().x0();
      if (particle.pdgcode().strangeness() > 0) {
        sampled_list.push_back(particle);
        S_plus += particle.pdgcode().strangeness();
      }
    }
  }

  // Mode 2: sample until strangeness is conserved
  while (S_plus + S_minus > conserved_initial.strangeness()) {
    sample_in_random_cell(mode_list, time);
    for (auto &particle : mode_list) {
      if (particle.pdgcode().strangeness() < 0) {
        sampled_list.push_back(particle);
        S_minus += particle.pdgcode().strangeness();
      }
    }
  }

  // Mode 3: sample non-strange baryons
  conserved_remaining = conserved_initial - QuantumNumbers(sampled_list);
  energy = 0.0;
  while (conserved_remaining.momentum().x0() > energy ||
         B_plus < conserved_remaining.baryon_number()) {
    sample_in_random_cell(mode_list, time);
    for (auto &particle : mode_list) {
      energy += particle.momentum().x0();
      if (particle.pdgcode().strangeness() == 0 &&
          particle.pdgcode().baryon_number() > 0) {
        sampled_list.push_back(particle);
        B_plus += particle.pdgcode().baryon_number();
      }
    }
  }

  // Mode 4: sample non-strange anti-baryons
  while (B_plus + B_minus > conserved_remaining.baryon_number()) {
    sample_in_random_cell(mode_list, time);
    for (auto &particle : mode_list) {
      if (particle.pdgcode().strangeness() == 0 &&
          particle.pdgcode().baryon_number() < 0) {
        sampled_list.push_back(particle);
        B_minus += particle.pdgcode().baryon_number();
      }
    }
  }

  // Mode 5: sample non_strange mesons, but take only with charge > 0
  conserved_remaining = conserved_initial - QuantumNumbers(sampled_list);
  energy = 0.0;
  while (conserved_remaining.momentum().x0() > energy ||
         E_plus < conserved_remaining.charge()) {
    sample_in_random_cell(mode_list, time);
    for (auto &particle : mode_list) {
      energy += particle.momentum().x0();
      if (particle.pdgcode().strangeness() == 0 &&
          particle.pdgcode().baryon_number() == 0 &&
          particle.pdgcode().charge() > 0) {
        sampled_list.push_back(particle);
        E_plus += particle.pdgcode().charge();
      }
    }
  }

  // Mode 6: sample non_strange mesons to conserve charge
  while (E_plus + E_minus > conserved_remaining.charge()) {
    sample_in_random_cell(mode_list, time);
    for (auto &particle : mode_list) {
      if (particle.pdgcode().strangeness() == 0 &&
          particle.pdgcode().baryon_number() == 0 &&
          particle.pdgcode().charge() < 0) {
        sampled_list.push_back(particle);
        E_minus += particle.pdgcode().charge();
      }
    }
  }

  // Mode 7: sample neutral non-strange mesons to conserve energy
  conserved_remaining = conserved_initial - QuantumNumbers(sampled_list);
  energy = 0.0;
  while (conserved_remaining.momentum().x0() > energy) {
    sample_in_random_cell(mode_list, time);
    for (auto &particle : mode_list) {
      if (particle.pdgcode().strangeness() == 0 &&
          particle.pdgcode().baryon_number() == 0 &&
          particle.pdgcode().charge() == 0) {
        sampled_list.push_back(particle);
        energy += particle.momentum().x0();
      }
    }
  }
  std::cout << "Sampled " << sampled_list.size() << " particles." << std::endl;

  // Report conservation laws:
  QuantumNumbers conserved_final = QuantumNumbers(sampled_list);
  std::cout << conserved_initial.report_deviations(conserved_final) << std::endl;

  // Add sampled particles to particles
  for (auto &particle : sampled_list) {
    particles.insert(particle);
  }
}

ThermLatticeNode::ThermLatticeNode() :
  Tmu0_(FourVector()),
  nb_(0.0),
  ns_(0.0),
  e_(0.0),
  p_(0.0),
  v_(ThreeVector()),
  T_(0.0),
  mub_(0.0),
  mus_(0.0) {}

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
    auto tabulated = eos.from_table(e_, gamma_inv*nb_);
    if (tabulated == nullptr) {
      auto T_mub_mus = eos.solve_eos(e_, gamma_inv*nb_, gamma_inv*ns_);
      T_   = T_mub_mus[0];
      mub_ = T_mub_mus[1];
      mus_ = T_mub_mus[2];
      p_ = HadronGasEos::pressure(T_, mub_, mus_);
    } else {
      p_ = tabulated->p;
      T_ = tabulated->T;
      mub_ = tabulated->mub;
      mus_ = tabulated->mus;
    }
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
