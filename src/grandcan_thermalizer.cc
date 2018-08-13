/*
 *    Copyright (c) 2016-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/grandcan_thermalizer.h"

#include <time.h>

#include "smash/angles.h"
#include "smash/bessel_sampler.h"
#include "smash/cxx14compat.h"
#include "smash/distributions.h"
#include "smash/forwarddeclarations.h"
#include "smash/logging.h"
#include "smash/particles.h"
#include "smash/quantumnumbers.h"
#include "smash/random.h"

namespace smash {

ThermLatticeNode::ThermLatticeNode()
    : Tmu0_(FourVector()),
      nb_(0.0),
      ns_(0.0),
      e_(0.0),
      p_(0.0),
      v_(ThreeVector()),
      T_(0.0),
      mub_(0.0),
      mus_(0.0) {}

void ThermLatticeNode::add_particle(const ParticleData &part, double factor) {
  Tmu0_ += part.momentum() * factor;
  nb_ += static_cast<double>(part.type().baryon_number()) * factor;
  ns_ += static_cast<double>(part.type().strangeness()) * factor;
}

void ThermLatticeNode::compute_rest_frame_quantities(HadronGasEos &eos) {
  /// \todo(oliiny): use Newton's method instead of these iterations
  const int max_iter = 50;
  v_ = ThreeVector(0.0, 0.0, 0.0);
  double e_previous_step = 0.0;
  const double tolerance = 5.e-4;
  int iter;
  for (iter = 0; iter < max_iter; iter++) {
    e_previous_step = e_;
    e_ = Tmu0_.x0() - Tmu0_.threevec() * v_;
    if (std::abs(e_ - e_previous_step) < tolerance) {
      break;
    }
    const double gamma_inv = std::sqrt(1.0 - v_.sqr());
    EosTable::table_element tabulated;
    eos.from_table(tabulated, e_, gamma_inv * nb_);
    if (!eos.is_tabulated() || tabulated.p < 0.0) {
      auto T_mub_mus = eos.solve_eos(e_, gamma_inv * nb_, gamma_inv * ns_);
      T_ = T_mub_mus[0];
      mub_ = T_mub_mus[1];
      mus_ = T_mub_mus[2];
      p_ = HadronGasEos::pressure(T_, mub_, mus_);
    } else {
      p_ = tabulated.p;
      T_ = tabulated.T;
      mub_ = tabulated.mub;
      mus_ = tabulated.mus;
    }
    v_ = Tmu0_.threevec() / (Tmu0_.x0() + p_);
  }
  if (iter == max_iter) {
    std::cout << "Warning from solver: max iterations exceeded."
              << " Accuracy: " << std::abs(e_ - e_previous_step)
              << " is less than tolerance " << tolerance << std::endl;
  }
}

void ThermLatticeNode::set_rest_frame_quantities(double T0, double mub0,
                                                 double mus0,
                                                 const ThreeVector v0) {
  T_ = T0;
  mub_ = mub0;
  mus_ = mus0;
  v_ = v0;
  e_ = HadronGasEos::energy_density(T_, mub_, mus_);
  p_ = HadronGasEos::pressure(T_, mub_, mus_);
  nb_ = HadronGasEos::net_baryon_density(T_, mub_, mus_);
  ns_ = HadronGasEos::net_strange_density(T_, mub_, mus_);
}

std::ostream &operator<<(std::ostream &out, const ThermLatticeNode &node) {
  return out << "T[mu,0]: " << node.Tmu0() << ", nb: " << node.nb()
             << ", ns: " << node.ns() << ", v: " << node.v()
             << ", e: " << node.e() << ", p: " << node.p()
             << ", T: " << node.T() << ", mub: " << node.mub()
             << ", mus: " << node.mus();
}

GrandCanThermalizer::GrandCanThermalizer(const std::array<double, 3> lat_sizes,
                                         const std::array<int, 3> n_cells,
                                         const std::array<double, 3> origin,
                                         bool periodicity, double e_critical,
                                         double t_start, double delta_t,
                                         ThermalizationAlgorithm algo)
    : eos_typelist_(list_eos_particles()),
      N_sorts_(eos_typelist_.size()),
      e_crit_(e_critical),
      t_start_(t_start),
      period_(delta_t),
      algorithm_(algo) {
  const LatticeUpdate upd = LatticeUpdate::EveryFixedInterval;
  lat_ = make_unique<RectangularLattice<ThermLatticeNode>>(
      lat_sizes, n_cells, origin, periodicity, upd);
  const std::array<double, 3> abc = lat_->cell_sizes();
  cell_volume_ = abc[0] * abc[1] * abc[2];
  cells_to_sample_.resize(50000);
  mult_sort_.resize(N_sorts_);
  mult_int_.resize(N_sorts_);
}

void GrandCanThermalizer::update_lattice(const Particles &particles,
                                         const DensityParameters &dens_par,
                                         bool ignore_cells_under_treshold) {
  const DensityType dens_type = DensityType::Hadron;
  const LatticeUpdate update = LatticeUpdate::EveryFixedInterval;
  update_general_lattice(lat_.get(), update, dens_type, dens_par, particles);
  for (auto &node : *lat_) {
    /* If energy density is definitely below e_crit -
       no need to find T, mu, etc. So if e = T00 - T0i*vi <=
       T00 + sum abs(T0i) < e_crit, no efforts are necessary. */
    if (!ignore_cells_under_treshold ||
        node.Tmu0().x0() + std::abs(node.Tmu0().x1()) +
                std::abs(node.Tmu0().x2()) + std::abs(node.Tmu0().x3()) >=
            e_crit_) {
      node.compute_rest_frame_quantities(eos_);
    } else {
      node = ThermLatticeNode();
    }
  }
}

ThreeVector GrandCanThermalizer::uniform_in_cell() const {
  return ThreeVector(random::uniform(-0.5 * lat_->cell_sizes()[0],
                                     +0.5 * lat_->cell_sizes()[0]),
                     random::uniform(-0.5 * lat_->cell_sizes()[1],
                                     +0.5 * lat_->cell_sizes()[1]),
                     random::uniform(-0.5 * lat_->cell_sizes()[2],
                                     +0.5 * lat_->cell_sizes()[2]));
}

void GrandCanThermalizer::renormalize_momenta(
    ParticleList &plist, const FourVector required_total_momentum) {
  const auto &log = logger<LogArea::GrandcanThermalizer>();

  // Centralize momenta
  QuantumNumbers conserved = QuantumNumbers(plist);
  log.info("Required 4-momentum: ", required_total_momentum);
  log.info("Sampled 4-momentum: ", conserved.momentum());
  const ThreeVector mom_to_add =
      (required_total_momentum.threevec() - conserved.momentum().threevec()) /
      plist.size();
  log.info("Adjusting momenta by ", mom_to_add);
  for (auto &particle : plist) {
    particle.set_4momentum(particle.type().mass(),
                           particle.momentum().threevec() + mom_to_add);
  }

  // Boost every particle to the common center of mass frame
  conserved = QuantumNumbers(plist);
  const ThreeVector beta_CM_generated = conserved.momentum().velocity();
  const ThreeVector beta_CM_required = required_total_momentum.velocity();

  double E = 0.0;
  double E_expected = required_total_momentum.abs();
  for (auto &particle : plist) {
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
    a_max = 1.0;
  } else {
    a_min = -1.0;
    a_max = 0.0;
  }
  do {
    a = 0.5 * (a_min + a_max);
    E = 0.0;
    for (const auto &particle : plist) {
      const double p2 = particle.momentum().threevec().sqr();
      const double E2 = particle.momentum().x0() * particle.momentum().x0();
      E += std::sqrt(E2 + a * (a + 2.0) * p2);
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
  for (auto &particle : plist) {
    particle.set_4momentum(particle.type().mass(),
                           (1 + a) * particle.momentum().threevec());
    particle.boost_momentum(-beta_CM_required);
  }
}

void GrandCanThermalizer::sample_multinomial(HadronClass particle_class,
                                             int N_to_sample) {
  /* The array mult_sort_ contains real numbers \f$ a_i \f$. The numbers \f$
   * n_i \f$ are saved in the mult_int_ array. Only particles of class
   * particle_class are sampled, where particle_class is defined by the
   * get_class function. */
  double sum = mult_class(particle_class);
  for (size_t i_type = 0; (i_type < N_sorts_) && (N_to_sample > 0); i_type++) {
    if (get_class(i_type) != particle_class) {
      continue;
    }
    const double p = mult_sort_[i_type] / sum;
    mult_int_[i_type] = random::binomial(N_to_sample, p);
    /// \todo (oliiny) what to do with this output?
    /*std::cout << eos_typelist_[i_type]->name() <<
             ": mult_sort = " << mult_sort_[i_type] <<
             ", sum = " << sum <<
             ", p = " << p <<
             ", N to sample = " << N_to_sample <<
             ", mult_int_ = " << mult_int_[i_type] << std::endl;*/
    sum -= mult_sort_[i_type];
    N_to_sample -= mult_int_[i_type];
  }
}

void GrandCanThermalizer::sample_in_random_cell_BF_algo(ParticleList &plist,
                                                        const double time,
                                                        size_t type_index) {
  N_in_cells_.clear();
  N_total_in_cells_ = 0.0;
  for (auto cell_index : cells_to_sample_) {
    const ThermLatticeNode cell = (*lat_)[cell_index];
    const double gamma = 1.0 / std::sqrt(1.0 - cell.v().sqr());
    const double N_this_cell =
        cell_volume_ * gamma *
        HadronGasEos::partial_density(*eos_typelist_[type_index], cell.T(),
                                      cell.mub(), cell.mus());
    N_in_cells_.push_back(N_this_cell);
    N_total_in_cells_ += N_this_cell;
  }

  for (int i = 0; i < mult_int_[type_index]; i++) {
    // Choose random cell, probability = N_in_cell/N_total
    double r = random::uniform(0.0, N_total_in_cells_);
    double partial_sum = 0.0;
    int index_only_thermalized = -1;
    while (partial_sum < r) {
      index_only_thermalized++;
      partial_sum += N_in_cells_[index_only_thermalized];
    }
    const int cell_index = cells_to_sample_[index_only_thermalized];
    const ThermLatticeNode cell = (*lat_)[cell_index];
    const ThreeVector cell_center = lat_->cell_center(cell_index);

    ParticleData particle(*eos_typelist_[type_index]);
    // Note: it's pole mass for resonances!
    const double m = eos_typelist_[type_index]->mass();
    // Position
    particle.set_4position(FourVector(time, cell_center + uniform_in_cell()));
    // Momentum
    double momentum_radial = sample_momenta_from_thermal(cell.T(), m);
    Angles phitheta;
    phitheta.distribute_isotropically();
    particle.set_4momentum(m, phitheta.threevec() * momentum_radial);
    particle.boost_momentum(-cell.v());
    particle.set_formation_time(time);

    plist.push_back(particle);
  }
}

void GrandCanThermalizer::thermalize_BF_algo(QuantumNumbers &conserved_initial,
                                             double time, int ntest) {
  const auto &log = logger<LogArea::GrandcanThermalizer>();

  std::fill(mult_sort_.begin(), mult_sort_.end(), 0.0);
  for (auto cell_index : cells_to_sample_) {
    const ThermLatticeNode cell = (*lat_)[cell_index];
    const double gamma = 1.0 / std::sqrt(1.0 - cell.v().sqr());
    for (size_t i = 0; i < N_sorts_; i++) {
      // N_i = n u^mu dsigma_mu = (isochronous hypersurface) n * V * gamma
      mult_sort_[i] += cell_volume_ * gamma * ntest *
                       HadronGasEos::partial_density(
                           *eos_typelist_[i], cell.T(), cell.mub(), cell.mus());
    }
  }

  std::fill(mult_classes_.begin(), mult_classes_.end(), 0.0);
  for (size_t i = 0; i < N_sorts_; i++) {
    mult_classes_[static_cast<size_t>(get_class(i))] += mult_sort_[i];
  }

  BesselSampler bessel_sampler_B(mult_class(HadronClass::Baryon),
                                 mult_class(HadronClass::Antibaryon),
                                 conserved_initial.baryon_number());

  while (true) {
    sampled_list_.clear();
    std::fill(mult_int_.begin(), mult_int_.end(), 0);
    const auto Nbar_antibar = bessel_sampler_B.sample();

    sample_multinomial(HadronClass::Baryon, Nbar_antibar.first);
    sample_multinomial(HadronClass::Antibaryon, Nbar_antibar.second);

    // Count strangeness of the sampled particles
    int S_sampled = 0;
    for (size_t i = 0; i < N_sorts_; i++) {
      S_sampled += eos_typelist_[i]->strangeness() * mult_int_[i];
    }

    std::pair<int, int> NS_antiS;
    if (algorithm_ == ThermalizationAlgorithm::BiasedBF) {
      BesselSampler bessel_sampler_S(
          mult_class(HadronClass::PositiveSMeson),
          mult_class(HadronClass::NegativeSMeson),
          conserved_initial.strangeness() - S_sampled);
      NS_antiS = bessel_sampler_S.sample();
    } else if (algorithm_ == ThermalizationAlgorithm::UnbiasedBF) {
      NS_antiS = std::make_pair(
          random::poisson(mult_class(HadronClass::PositiveSMeson)),
          random::poisson(mult_class(HadronClass::NegativeSMeson)));
      if (NS_antiS.first - NS_antiS.second !=
          conserved_initial.strangeness() - S_sampled) {
        continue;
      }
    }

    sample_multinomial(HadronClass::PositiveSMeson, NS_antiS.first);
    sample_multinomial(HadronClass::NegativeSMeson, NS_antiS.second);
    // Count charge of the sampled particles
    int ch_sampled = 0;
    for (size_t i = 0; i < N_sorts_; i++) {
      ch_sampled += eos_typelist_[i]->charge() * mult_int_[i];
    }

    std::pair<int, int> NC_antiC;
    if (algorithm_ == ThermalizationAlgorithm::BiasedBF) {
      BesselSampler bessel_sampler_C(
          mult_class(HadronClass::PositiveQZeroSMeson),
          mult_class(HadronClass::NegativeQZeroSMeson),
          conserved_initial.charge() - ch_sampled);
      NC_antiC = bessel_sampler_C.sample();
    } else if (algorithm_ == ThermalizationAlgorithm::UnbiasedBF) {
      NC_antiC = std::make_pair(
          random::poisson(mult_class(HadronClass::PositiveQZeroSMeson)),
          random::poisson(mult_class(HadronClass::NegativeQZeroSMeson)));
      if (NC_antiC.first - NC_antiC.second !=
          conserved_initial.charge() - ch_sampled) {
        continue;
      }
    }

    sample_multinomial(HadronClass::PositiveQZeroSMeson, NC_antiC.first);
    sample_multinomial(HadronClass::NegativeQZeroSMeson, NC_antiC.second);
    sample_multinomial(
        HadronClass::ZeroQZeroSMeson,
        random::poisson(mult_class(HadronClass::ZeroQZeroSMeson)));

    for (size_t itype = 0; itype < N_sorts_; itype++) {
      sample_in_random_cell_BF_algo(sampled_list_, time, itype);
    }
    double e_tot;
    const double e_init = conserved_initial.momentum().x0();
    e_tot = 0.0;
    for (auto &particle : sampled_list_) {
      e_tot += particle.momentum().x0();
    }
    if (std::abs(e_tot - e_init) > 0.01 * e_init) {
      log.debug("Rejecting: energy ", e_tot, " too far from e_init = ", e_init);
      continue;
    }
    break;
  }
}

void GrandCanThermalizer::thermalize_mode_algo(
    QuantumNumbers &conserved_initial, double time) {
  double energy = 0.0;
  int S_plus = 0, S_minus = 0, B_plus = 0, B_minus = 0, E_plus = 0, E_minus = 0;
  // Mode 1: sample until energy is conserved, take only strangeness < 0
  auto condition1 = [](int, int, int) { return true; };
  compute_N_in_cells_mode_algo(condition1);
  while (conserved_initial.momentum().x0() > energy ||
         S_plus < conserved_initial.strangeness()) {
    ParticleData p = sample_in_random_cell_mode_algo(time, condition1);
    energy += p.momentum().x0();
    if (p.pdgcode().strangeness() > 0) {
      sampled_list_.push_back(p);
      S_plus += p.pdgcode().strangeness();
    }
  }

  // Mode 2: sample until strangeness is conserved
  auto condition2 = [](int S, int, int) { return (S < 0); };
  compute_N_in_cells_mode_algo(condition2);
  while (S_plus + S_minus > conserved_initial.strangeness()) {
    ParticleData p = sample_in_random_cell_mode_algo(time, condition2);
    const int s_part = p.pdgcode().strangeness();
    // Do not allow particles with S = -2 or -3 spoil the total sum
    if (S_plus + S_minus + s_part >= conserved_initial.strangeness()) {
      sampled_list_.push_back(p);
      S_minus += s_part;
    }
  }

  // Mode 3: sample non-strange baryons
  auto condition3 = [](int S, int, int) { return (S == 0); };
  QuantumNumbers conserved_remaining =
      conserved_initial - QuantumNumbers(sampled_list_);
  energy = 0.0;
  compute_N_in_cells_mode_algo(condition3);
  while (conserved_remaining.momentum().x0() > energy ||
         B_plus < conserved_remaining.baryon_number()) {
    ParticleData p = sample_in_random_cell_mode_algo(time, condition3);
    energy += p.momentum().x0();
    if (p.pdgcode().baryon_number() > 0) {
      sampled_list_.push_back(p);
      B_plus += p.pdgcode().baryon_number();
    }
  }

  // Mode 4: sample non-strange anti-baryons
  auto condition4 = [](int S, int B, int) { return (S == 0) && (B < 0); };
  compute_N_in_cells_mode_algo(condition4);
  while (B_plus + B_minus > conserved_remaining.baryon_number()) {
    ParticleData p = sample_in_random_cell_mode_algo(time, condition4);
    const int bar = p.pdgcode().baryon_number();
    if (B_plus + B_minus + bar >= conserved_remaining.baryon_number()) {
      sampled_list_.push_back(p);
      B_minus += bar;
    }
  }

  // Mode 5: sample non_strange mesons, but take only with charge > 0
  auto condition5 = [](int S, int B, int) { return (S == 0) && (B == 0); };
  conserved_remaining = conserved_initial - QuantumNumbers(sampled_list_);
  energy = 0.0;
  compute_N_in_cells_mode_algo(condition5);
  while (conserved_remaining.momentum().x0() > energy ||
         E_plus < conserved_remaining.charge()) {
    ParticleData p = sample_in_random_cell_mode_algo(time, condition5);
    energy += p.momentum().x0();
    if (p.pdgcode().charge() > 0) {
      sampled_list_.push_back(p);
      E_plus += p.pdgcode().charge();
    }
  }

  // Mode 6: sample non_strange mesons to conserve charge
  auto condition6 = [](int S, int B, int C) {
    return (S == 0) && (B == 0) && (C < 0);
  };
  compute_N_in_cells_mode_algo(condition6);
  while (E_plus + E_minus > conserved_remaining.charge()) {
    ParticleData p = sample_in_random_cell_mode_algo(time, condition6);
    const int charge = p.pdgcode().charge();
    if (E_plus + E_minus + charge >= conserved_remaining.charge()) {
      sampled_list_.push_back(p);
      E_minus += charge;
    }
  }

  // Mode 7: sample neutral non-strange mesons to conserve energy
  auto condition7 = [](int S, int B, int C) {
    return (S == 0) && (B == 0) && (C == 0);
  };
  conserved_remaining = conserved_initial - QuantumNumbers(sampled_list_);
  energy = 0.0;
  compute_N_in_cells_mode_algo(condition7);
  while (conserved_remaining.momentum().x0() > energy) {
    ParticleData p = sample_in_random_cell_mode_algo(time, condition7);
    sampled_list_.push_back(p);
    energy += p.momentum().x0();
  }
}

void GrandCanThermalizer::thermalize(const Particles &particles, double time,
                                     int ntest) {
  const auto &log = logger<LogArea::GrandcanThermalizer>();
  log.info("Starting forced thermalization, time ", time, " fm/c");
  to_remove_.clear();
  sampled_list_.clear();
  /* Remove particles from the cells with e > e_crit_,
   * sum up their conserved quantities */
  QuantumNumbers conserved_initial = QuantumNumbers();
  ThermLatticeNode node;
  for (auto &particle : particles) {
    const bool is_on_lattice =
        lat_->value_at(particle.position().threevec(), node);
    if (is_on_lattice && node.e() > e_crit_) {
      to_remove_.push_back(particle);
    }
  }
  /* Do not thermalize too small number of particles: for the number
   * of particles < 30 the algorithm tends to hang or crash too often. */
  if (to_remove_.size() > 30) {
    for (auto &particle : to_remove_) {
      conserved_initial.add_values(particle);
    }
  } else {
    to_remove_.clear();
    conserved_initial = QuantumNumbers();
  }
  log.info("Removed ", to_remove_.size(), " particles.");

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
           cells_to_sample_.size() * cell_volume_, ", in % of lattice: ",
           100.0 * cells_to_sample_.size() / lattice_total_cells);

  switch (algorithm_) {
    case ThermalizationAlgorithm::BiasedBF:
    case ThermalizationAlgorithm::UnbiasedBF:
      thermalize_BF_algo(conserved_initial, time, ntest);
      break;
    case ThermalizationAlgorithm::ModeSampling:
      thermalize_mode_algo(conserved_initial, time);
      break;
    default:
      throw std::invalid_argument(
          "This thermalization algorithm is"
          " not yet implemented");
  }
  log.info("Sampled ", sampled_list_.size(), " particles.");

  // Adjust momenta
  renormalize_momenta(sampled_list_, conserved_initial.momentum());
}

void GrandCanThermalizer::print_statistics(const Clock &clock) const {
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
    on_lattice.T += node.T() * e;
    on_lattice.mub += node.mub() * e;
    on_lattice.mus += node.mus() * e;
    on_lattice.nb += node.nb() * e;
    on_lattice.ns += node.ns() * e;
    e_sum_on_lattice += e;
    if (e >= e_crit_) {
      in_therm_reg.T += node.T() * e;
      in_therm_reg.mub += node.mub() * e;
      in_therm_reg.mus += node.mus() * e;
      in_therm_reg.nb += node.nb() * e;
      in_therm_reg.ns += node.ns() * e;
      e_sum_in_therm_reg += e;
      node_counter++;
    }
  }
  if (e_sum_on_lattice > really_small) {
    on_lattice.T /= e_sum_on_lattice;
    on_lattice.mub /= e_sum_on_lattice;
    on_lattice.mus /= e_sum_on_lattice;
    on_lattice.nb /= e_sum_on_lattice;
    on_lattice.ns /= e_sum_on_lattice;
  }
  if (e_sum_in_therm_reg > really_small) {
    in_therm_reg.T /= e_sum_in_therm_reg;
    in_therm_reg.mub /= e_sum_in_therm_reg;
    in_therm_reg.mus /= e_sum_in_therm_reg;
    in_therm_reg.nb /= e_sum_in_therm_reg;
    in_therm_reg.ns /= e_sum_in_therm_reg;
  }

  std::cout << "Current time [fm/c]: " << clock.current_time() << std::endl;
  std::cout << "Averages on the lattice - T[GeV], mub[GeV], mus[GeV], "
            << "nb[fm^-3], ns[fm^-3]: " << on_lattice.T << " " << on_lattice.mub
            << " " << on_lattice.mus << " " << on_lattice.nb << " "
            << on_lattice.ns << std::endl;
  std::cout << "Averages in therm. region - T[GeV], mub[GeV], mus[GeV], "
            << "nb[fm^-3], ns[fm^-3]: " << in_therm_reg.T << " "
            << in_therm_reg.mub << " " << in_therm_reg.mus << " "
            << in_therm_reg.nb << " " << in_therm_reg.ns << std::endl;
  std::cout << "Volume with e > e_crit [fm^3]: " << cell_volume_ * node_counter
            << std::endl;
}

}  // namespace smash
