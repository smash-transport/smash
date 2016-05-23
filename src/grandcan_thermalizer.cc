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
#include "include/logging.h"
#include "include/particles.h"
#include "include/quantumnumbers.h"
#include "include/random.h"

namespace Smash {

GrandCanThermalizer::GrandCanThermalizer(const std::array<float, 3> cell_sizes,
                                         const std::array<int, 3> n_cells,
                                         float e_critical,
                                         float t_start,
                                         float delta_t) :
  e_crit_(e_critical),
  t_start_(t_start),
  period_(delta_t) {
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
                                         const DensityParameters& dens_par,
                                         bool ignore_cells_under_treshold) {
  const DensityType dens_type = DensityType::Hadron;
  const LatticeUpdate update = LatticeUpdate::EveryFixedInterval;
  update_general_lattice(lat_.get(), update, dens_type, dens_par, particles);
  for (auto &node : *lat_) {
    /* If energy density is definitely below e_crit -
       no need to find T, mu, etc. So if e = T00 - T0i*vi <=
       T00 + sum abs(T0i) < e_crit, no efforts are necessary. */
    if (!ignore_cells_under_treshold ||
        node.Tmu0().x0() +
        std::abs(node.Tmu0().x1()) +
        std::abs(node.Tmu0().x2()) +
        std::abs(node.Tmu0().x3()) >= e_crit_) {
      node.compute_rest_frame_quantities(eos_);
    } else {
      node = ThermLatticeNode();
    }
  }
}

void GrandCanThermalizer::thermalize(Particles& particles, double time) {
  const auto &log = logger<LogArea::GrandcanThermalizer>();
  log.info("Starting forced thermalization, time ", time, " fm/c");
  // Remove particles from the cells with e > e_crit_,
  // sum up their conserved quantities
  QuantumNumbers conserved_initial   = QuantumNumbers(),
                 conserved_remaining = QuantumNumbers();
  ThermLatticeNode node;
  ParticleList to_remove;
  for (auto &particle : particles) {
    const bool is_on_lattice = lat_->value_at(particle.position().threevec(),
                                              node);
    if (is_on_lattice && node.e() > e_crit_) {
      to_remove.push_back(particle);
    }
  }
  // Do not thermalize too small number of particles
  if (to_remove.size() > 30) {
    for (auto &particle : to_remove) {
      conserved_initial.add_values(particle);
      particles.remove(particle);
    }
  } else {
    to_remove.clear();
    conserved_initial = QuantumNumbers();
  }
  log.info("Removed ", to_remove.size(), " particles.");

  // Exit if there is nothing to thermalize
  if (conserved_initial == QuantumNumbers()) {
    return;
  }
  // Save the indices of cells inside the volume with e > e_crit_
  cells_to_sample_.clear();
  const size_t lattice_total_cells = lat_->size();
  for (size_t i = 0; i < lattice_total_cells; i++) {
    if ((*lat_)[i].e() > e_crit_) {
      cells_to_sample_.push_back(i);
    }
  }
  log.info("Number of cells in the thermalization region = ",
           cells_to_sample_.size(), ", its total volume [fm^3]: ",
           cells_to_sample_.size()*cell_volume_, ", in \% of lattice: ",
           100.0*cells_to_sample_.size()*cell_volume_/lattice_total_cells);

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
      const int s_part = particle.pdgcode().strangeness();
      if (s_part < 0 &&
          // Do not allow particles with S = -2 or -3 spoil the total sum
          S_plus + S_minus + s_part >= conserved_initial.strangeness()) {
        sampled_list.push_back(particle);
        S_minus += s_part;
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
      const int bar = particle.pdgcode().baryon_number();
      if (particle.pdgcode().strangeness() == 0 && bar < 0 &&
          B_plus + B_minus + bar >= conserved_remaining.baryon_number()) {
        sampled_list.push_back(particle);
        B_minus += bar;
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
      const int charge = particle.pdgcode().charge();
      if (particle.pdgcode().strangeness() == 0 &&
          particle.pdgcode().baryon_number() == 0 &&
          charge < 0 &&
          E_plus + E_minus + charge >= conserved_remaining.charge()) {
        sampled_list.push_back(particle);
        E_minus += charge;
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
  log.info("Sampled ", sampled_list.size(), " particles.");

  // Centralize momenta
  QuantumNumbers conserved_final = QuantumNumbers(sampled_list);
  log.info("Initial particles' 4-momentum: ", conserved_initial.momentum());
  log.info("Samples particles' 4-momentum: ", conserved_final.momentum());
  const QuantumNumbers deviation = conserved_initial - conserved_final;
  const ThreeVector mom_to_add = deviation.momentum().threevec() /
                                 sampled_list.size();
  log.info("Adjusting momenta by ", mom_to_add);
  for (auto &particle : sampled_list) {
    particle.set_4momentum(particle.type().mass(),
                           particle.momentum().threevec() + mom_to_add);
  }

  // Boost every particle to the common center of mass frame
  conserved_final = QuantumNumbers(sampled_list);
  const ThreeVector beta_CM_generated = conserved_final.momentum().velocity();
  const ThreeVector beta_CM_initial = conserved_initial.momentum().velocity();

  double E = 0.0;
  double E_expected = conserved_initial.momentum().abs();
  for (auto &particle : sampled_list) {
    particle.boost_momentum(beta_CM_generated);
    E += particle.momentum().x0();
  }
  // Renorm. momenta by factor (1+a) to get the right energy, binary search
  const double tolerance = really_small;
  double a, a_min, a_max, er;
  const int max_iter = 50;
  int iter = 0;
  if (E_expected >= E) {
    a_min = 0.0;
    a_max = 0.5;
  } else {
    a_min = -0.5;
    a_max = 0.0;
  }
  do {
    a = 0.5 * (a_min + a_max);
    E = 0.0;
    for (const auto &particle : sampled_list) {
      const double p2 = particle.momentum().threevec().sqr();
      const double E2 = particle.momentum().x0() * particle.momentum().x0();
      E += std::sqrt(E2 + a*(a + 2.0) * p2);
    }
    er = E - E_expected;
    if (er >= 0.0) {
      a_max = a;
    } else {
      a_min = a;
    }
    log.debug("Iteration ", iter, ": a = ", a, ", Δ = ", er);
    iter++;
  } while (std::abs(er) > tolerance && iter < max_iter);

  log.info("Renormalizing momenta by factor 1+a, a = ", a);
  for (auto &particle : sampled_list) {
    particle.set_4momentum(particle.type().mass(),
                           (1+a)*particle.momentum().threevec());
    particle.boost_momentum(-beta_CM_initial);
  }

  // Add sampled particles to particles
  for (auto &particle : sampled_list) {
    particles.insert(particle);
  }
}

void GrandCanThermalizer::print_statistics(const Clock& clock) const {
  struct to_average {
    double T;
    double mub;
    double mus;
    double nb;
    double ns;
  };
  struct to_average on_lattice = {0.0, 0.0, 0.0, 0.0, 0.0};
  struct to_average in_therm_reg = {0.0, 0.0, 0.0, 0.0, 0.0};
  double e_sum_on_lattice = 0.0, e_sum_in_therm_reg = 0.0;
  int node_counter = 0;
  for (const auto &node : *lat_) {
    const double e = node.e();
    on_lattice.T   += node.T()   * e;
    on_lattice.mub += node.mub() * e;
    on_lattice.mus += node.mus() * e;
    on_lattice.nb  += node.nb()  * e;
    on_lattice.ns  += node.ns()  * e;
    e_sum_on_lattice += e;
    if (e >= e_crit_) {
      in_therm_reg.T   += node.T()   * e;
      in_therm_reg.mub += node.mub() * e;
      in_therm_reg.mus += node.mus() * e;
      in_therm_reg.nb  += node.nb()  * e;
      in_therm_reg.ns  += node.ns()  * e;
      e_sum_in_therm_reg += e;
      node_counter++;
    }
  }
  if (e_sum_on_lattice > really_small) {
    on_lattice.T   /= e_sum_on_lattice;
    on_lattice.mub /= e_sum_on_lattice;
    on_lattice.mus /= e_sum_on_lattice;
    on_lattice.nb  /= e_sum_on_lattice;
    on_lattice.ns  /= e_sum_on_lattice;
  }
  if (e_sum_in_therm_reg > really_small) {
    in_therm_reg.T   /= e_sum_in_therm_reg;
    in_therm_reg.mub /= e_sum_in_therm_reg;
    in_therm_reg.mus /= e_sum_in_therm_reg;
    in_therm_reg.nb  /= e_sum_in_therm_reg;
    in_therm_reg.ns  /= e_sum_in_therm_reg;
  }

  std::cout << "Current time [fm/c]: " << clock.current_time() << std::endl;
  std::cout << "Averages on the lattice - T[GeV], mub[GeV], mus[GeV], "
            << "nb[fm^-3], ns[fm^-3]: "
            << on_lattice.T << " " << on_lattice.mub << " " << on_lattice.mus
            << " " << on_lattice.nb << " " << on_lattice.ns << std::endl;
  std::cout << "Averages in therm. region - T[GeV], mub[GeV], mus[GeV], "
            << "nb[fm^-3], ns[fm^-3]: "
            << in_therm_reg.T << " " << in_therm_reg.mub << " "
            << in_therm_reg.mus << " " << in_therm_reg.nb << " "
            << in_therm_reg.ns << std::endl;
  std::cout << "Volume with e > e_crit [fm^3]: "
            << cell_volume_ * node_counter << std::endl;
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
  // ToDo(oliiny): use Newton's method instead of these iterations
  const int max_iter = 50;
  v_ = ThreeVector(0.0, 0.0, 0.0);
  double e_previous_step = 0.0;
  const double tolerance = 5.e-4;
  int iter;
  for(iter = 0; iter < max_iter; iter++) {
    e_previous_step = e_;
    e_ = Tmu0_.x0() - Tmu0_.threevec()*v_;
    if (std::abs(e_ - e_previous_step) < tolerance) {
      break;
    }
    const double gamma_inv = std::sqrt(1.0 - v_.sqr());
    auto tabulated = eos.from_table(e_, gamma_inv*nb_);
    if (!eos.is_tabulated() || tabulated.p < 0.0) {
      auto T_mub_mus = eos.solve_eos(e_, gamma_inv*nb_, gamma_inv*ns_);
      T_   = T_mub_mus[0];
      mub_ = T_mub_mus[1];
      mus_ = T_mub_mus[2];
      p_ = HadronGasEos::pressure(T_, mub_, mus_);
    } else {
      p_ = tabulated.p;
      T_ = tabulated.T;
      mub_ = tabulated.mub;
      mus_ = tabulated.mus;
    }
    v_ = Tmu0_.threevec()/(Tmu0_.x0() + p_);
  }
  if (iter == max_iter) {
    std::cout << "Warning from solver: max iterations exceeded." <<
                 " Accuracy: " << std::abs(e_ - e_previous_step) <<
                 " is less than tolerance " << tolerance << std::endl;
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
