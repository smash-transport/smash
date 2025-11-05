/*
 *
 *    Copyright (c) 2014-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/action.h"

#include <assert.h>

#include <algorithm>
#include <sstream>

#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/logging.h"
#include "smash/pauliblocking.h"
#include "smash/potential_globals.h"
#include "smash/quantumnumbers.h"

namespace smash {
/// Destructor
Action::~Action() = default;
static constexpr int LPauliBlocking = LogArea::PauliBlocking::id;

bool Action::is_valid(const Particles &particles) const {
  return std::all_of(
      incoming_particles_.begin(), incoming_particles_.end(),
      [&particles](const ParticleData &p) { return particles.is_valid(p); });
}

bool Action::is_pauli_blocked(const std::vector<Particles> &ensembles,
                              const PauliBlocker &p_bl) const {
  // Wall-crossing actions should never be blocked: currently
  // if the action is blocked, a particle continues to propagate in a straight
  // line. This would simply bring it out of the box.
  if (process_type_ == ProcessType::Wall) {
    return false;
  }
  for (const auto &p : outgoing_particles_) {
    if (p.is_baryon()) {
      const auto f =
          p_bl.phasespace_dens(p.position().threevec(), p.momentum().threevec(),
                               ensembles, p.pdgcode(), incoming_particles_);
      if (f > random::uniform(0., 1.)) {
        logg[LPauliBlocking].debug("Action ", *this,
                                   " is pauli-blocked with f = ", f);
        return true;
      }
    }
  }
  return false;
}

const ParticleList &Action::incoming_particles() const {
  return incoming_particles_;
}

void Action::update_incoming(const Particles &particles) {
  for (auto &p : incoming_particles_) {
    p = particles.lookup(p);
  }
}

FourVector Action::get_interaction_point() const {
  // Estimate for the interaction point in the calculational frame
  ThreeVector interaction_point = ThreeVector(0., 0., 0.);
  std::vector<ThreeVector> propagated_positions;
  for (const auto &part : incoming_particles_) {
    ThreeVector propagated_position =
        part.position().threevec() +
        part.velocity() * (time_of_execution_ - part.position().x0());
    propagated_positions.push_back(propagated_position);
    interaction_point += propagated_position;
  }
  interaction_point /= incoming_particles_.size();
  /*
   * In case of periodic boundaries interaction point is not necessarily
   * (x1 + x2)/2. Consider only one dimension, e.g. x, the rest are analogous.
   * Instead of x, there can be x + k * L, where k is any integer and L
   * is period.Interaction point is either. Therefore, interaction point is
   * (x1 + k * L + x2 + m * L) / 2  = (x1 + x2) / 2 + n * L / 2. We need
   * this interaction point to be with [0, L], so n can be {-1, 0, 1}.
   * Which n to choose? Our guiding principle is that n should be such that
   * interaction point is closest to interacting particles.
   */
  if (box_length_ > 0 && stochastic_position_idx_ < 0) {
    assert(incoming_particles_.size() == 2);
    const ThreeVector r = propagated_positions[0] - propagated_positions[1];
    for (int i = 0; i < 3; i++) {
      const double d = std::abs(r[i]);
      if (d > 0.5 * box_length_) {
        if (interaction_point[i] >= 0.5 * box_length_) {
          interaction_point[i] -= 0.5 * box_length_;
        } else {
          interaction_point[i] += 0.5 * box_length_;
        }
      }
    }
  }
  /* In case of scatterings via the stochastic criterion, use postion of random
   * incoming particle to prevent density hotspots in grid cell centers. */
  if (stochastic_position_idx_ >= 0) {
    return incoming_particles_[stochastic_position_idx_].position();
  }
  return FourVector(time_of_execution_, interaction_point);
}

std::pair<FourVector, FourVector> Action::get_potential_at_interaction_point()
    const {
  const ThreeVector r = get_interaction_point().threevec();
  FourVector UB = FourVector();
  FourVector UI3 = FourVector();
  /* Check:
   * Lattice is turned on. */
  if (UB_lat_pointer != nullptr) {
    UB_lat_pointer->value_at(r, UB);
  }
  if (UI3_lat_pointer != nullptr) {
    UI3_lat_pointer->value_at(r, UI3);
  }
  return std::make_pair(UB, UI3);
}

double Action::perform(Particles *particles, uint32_t id_process) {
  assert(id_process != 0);
  double energy_violation = 0.;
  for (ParticleData &p : outgoing_particles_) {
    /* Store the history info. Wall crossing and fluidization don't change the
     * last collision a particle went through. */
    if ((process_type_ != ProcessType::Wall) &&
        (process_type_ != ProcessType::FluidizationNoRemoval)) {
      p.set_history(p.get_history().collisions_per_particle + 1, id_process,
                    process_type_, time_of_execution_, incoming_particles_);
    }
  }

  /* For elastic collisions and box wall crossings it is not necessary to remove
   * particles from the list and insert new ones, it is enough to update their
   * properties. */
  const bool replace = (process_type_ != ProcessType::Elastic) &&
                       (process_type_ != ProcessType::Wall) &&
                       (process_type_ != ProcessType::FluidizationNoRemoval);
  particles->update(incoming_particles_, outgoing_particles_, replace);

  logg[LAction].debug("Particle map now has ", particles->size(), " elements.");

  /* Check the conservation laws if the modifications of the total kinetic
   * energy of the outgoing particles by the mean field potentials are not
   * taken into account. */
  if (UB_lat_pointer == nullptr && UI3_lat_pointer == nullptr) {
    energy_violation = check_conservation(id_process);
  }
  return energy_violation;
}

FourVector Action::total_momentum_of_outgoing_particles() const {
  const auto potentials = get_potential_at_interaction_point();
  /* scale_B returns the difference of the total force scales of the skyrme
   * potential between the initial and final states. */
  double scale_B = 0.0;
  /* scale_I3 returns the difference of the total force scales of the symmetry
   * potential between the initial and final states. */
  double scale_I3 = 0.0;
  for (const auto &p_in : incoming_particles_) {
    // Get the force scale of the incoming particle.
    const auto scale =
        ((pot_pointer != nullptr) ? pot_pointer->force_scale(p_in.type())
                                  : std::make_pair(0.0, 0));
    scale_B += scale.first;
    scale_I3 += scale.second * p_in.type().isospin3_rel();
  }
  for (const auto &p_out : outgoing_particles_) {
    // Get the force scale of the outgoing particle.
    const auto scale = ((pot_pointer != nullptr)
                            ? pot_pointer->force_scale(type_of_pout(p_out))
                            : std::make_pair(0.0, 0));
    scale_B -= scale.first;
    scale_I3 -= scale.second * type_of_pout(p_out).isospin3_rel();
  }
  /* Rescale to get the potential difference between the
   * initial and final state, and thus get the total momentum
   * of the outgoing particles*/
  return total_momentum() + potentials.first * scale_B +
         potentials.second * scale_I3;
}

void Action::assign_formation_time_to_outgoing_particles() {
  /* Find incoming particle with largest formation time i.e. the last formed
   * incoming particle. If all particles form at the same time, take the one
   * with the lowest cross section scaling factor */
  ParticleList::iterator last_formed_in_part;
  bool all_incoming_same_formation_time =
      std::all_of(incoming_particles_.begin() + 1, incoming_particles_.end(),
                  [&](const ParticleData &data_comp) {
                    return std::abs(incoming_particles_[0].formation_time() -
                                    data_comp.formation_time()) < really_small;
                  });
  if (all_incoming_same_formation_time) {
    last_formed_in_part =
        std::min_element(incoming_particles_.begin(), incoming_particles_.end(),
                         [](const ParticleData &a, const ParticleData &b) {
                           return a.initial_xsec_scaling_factor() <
                                  b.initial_xsec_scaling_factor();
                         });
  } else {
    last_formed_in_part =
        std::max_element(incoming_particles_.begin(), incoming_particles_.end(),
                         [](const ParticleData &a, const ParticleData &b) {
                           return a.formation_time() < b.formation_time();
                         });
  }

  const double form_time_begin = last_formed_in_part->begin_formation_time();
  const double sc = last_formed_in_part->initial_xsec_scaling_factor();

  if (last_formed_in_part->formation_time() > time_of_execution_) {
    for (ParticleData &new_particle : outgoing_particles_) {
      if (new_particle.initial_xsec_scaling_factor() < 1.0) {
        /* The new cross section scaling factor will be the product of the
         * cross section scaling factor of the ingoing particles and of the
         * outgoing ones (since the outgoing ones are also string fragments
         * and thus take time to form). */
        double sc_out = new_particle.initial_xsec_scaling_factor();
        new_particle.set_cross_section_scaling_factor(sc * sc_out);
        if (last_formed_in_part->formation_time() >
            new_particle.formation_time()) {
          /* If the unformed incoming particles' formation time is larger than
           * the current outgoing particle's formation time, then the latter
           * is overwritten by the former*/
          new_particle.set_slow_formation_times(
              time_of_execution_, last_formed_in_part->formation_time());
        }
      } else {
        // not a string product
        new_particle.set_slow_formation_times(
            form_time_begin, last_formed_in_part->formation_time());
        new_particle.set_cross_section_scaling_factor(sc);
      }
    }
  } else {
    for (ParticleData &new_particle : outgoing_particles_) {
      if (new_particle.initial_xsec_scaling_factor() == 1.0) {
        new_particle.set_formation_time(time_of_execution_);
      }
    }
  }
}

std::pair<double, double> Action::sample_masses(
    const double kinetic_energy_cm) const {
  const ParticleType &t_a = outgoing_particles_[0].type();
  const ParticleType &t_b = outgoing_particles_[1].type();
  // start with pole masses
  std::pair<double, double> masses = {t_a.mass(), t_b.mass()};

  if (kinetic_energy_cm < t_a.min_mass_kinematic() + t_b.min_mass_kinematic()) {
    const std::string reaction = incoming_particles_[0].type().name() +
                                 incoming_particles_[1].type().name() + "→" +
                                 t_a.name() + t_b.name();
    throw InvalidResonanceFormation(
        reaction + ": not enough energy, " + std::to_string(kinetic_energy_cm) +
        " < " + std::to_string(t_a.min_mass_kinematic()) + " + " +
        std::to_string(t_b.min_mass_kinematic()));
  }

  /* If one of the particles is a resonance, sample its mass. */
  if (!t_a.is_stable() && t_b.is_stable()) {
    masses.first = t_a.sample_resonance_mass(t_b.mass(), kinetic_energy_cm);
  } else if (!t_b.is_stable() && t_a.is_stable()) {
    masses.second = t_b.sample_resonance_mass(t_a.mass(), kinetic_energy_cm);
  } else if (!t_a.is_stable() && !t_b.is_stable()) {
    // two resonances in final state
    masses = t_a.sample_resonance_masses(t_b, kinetic_energy_cm);
  }
  return masses;
}

void Action::sample_angles(std::pair<double, double> masses,
                           const double kinetic_energy_cm) {
  ParticleData *p_a = &outgoing_particles_[0];
  ParticleData *p_b = &outgoing_particles_[1];

  const double pcm = pCM(kinetic_energy_cm, masses.first, masses.second);
  if (!(pcm > 0.0)) {
    logg[LAction].warn("Particle: ", p_a->pdgcode(), " radial momentum: ", pcm);
    logg[LAction].warn("Ektot: ", kinetic_energy_cm, " m_a: ", masses.first,
                       " m_b: ", masses.second);
  }
  /* Here we assume an isotropic angular distribution. */
  Angles phitheta;
  phitheta.distribute_isotropically();

  p_a->set_4momentum(masses.first, phitheta.threevec() * pcm);
  p_b->set_4momentum(masses.second, -phitheta.threevec() * pcm);
  /* Debug message is printed before boost, so that p_a and p_b are
   * the momenta in the center of mass frame and thus opposite to
   * each other.*/
  logg[LAction].debug("p_a: ", *p_a, "\np_b: ", *p_b);
}

void Action::sample_2body_phasespace() {
  /* This function only operates on 2-particle final states. */
  assert(outgoing_particles_.size() == 2);
  const FourVector p_tot = total_momentum_of_outgoing_particles();
  const double cm_kin_energy = p_tot.abs();
  // first sample the masses
  const std::pair<double, double> masses = sample_masses(cm_kin_energy);
  // after the masses are fixed (and thus also pcm), sample the angles
  sample_angles(masses, cm_kin_energy);
}

void Action::sample_manybody_phasespace_impl(
    double sqrts, const std::vector<double> &m,
    std::vector<FourVector> &sampled_momenta) {
  /**
   *  Using the M-method from CERN-68-15 report, paragraph 9.6
   *  1) Generate invariant masses M12, M123, M1234, etc from
   *      distribution dM12 x dM123 x dM1234 x ...
   *      This is not trivial because of the integration limits.
   *      Here the idea is to change variables to T12 = M12 - (m1 + m2),
   *      T123 = M123 - (m1 + m2 + m3), etc. Then we need to generate
   *      uniform T such that 0 <= T12 <= T123 <= T1234 <= ... <= sqrts - sum
   * (m_i). For the latter there is a trick: generate values uniformly in [0,
   * sqrts - sum (m_i)] and then sort the values. 2) accept or reject this
   * combination of invariant masses with weight proportional to R2(sqrt,
   * M_{n-1}, m_n) x R2(M_{n-1}, M_{n-2}, m_{n-1}) x
   *     ... x R2(M2, m1, m2) x (prod M_i). Maximum weight is estmated
   * heuristically, here I'm using an idea by Scott Pratt that maximum is close
   * to T12 = T123 = T1234 = ... = (sqrts - sum (m_i)) / (n - 1)
   */
  const size_t n = m.size();
  assert(n > 1);
  sampled_momenta.resize(n);

  // Arrange a convenient vector of m1, m1 + m2, m1 + m2 + m3, ...
  std::vector<double> msum(n);
  std::partial_sum(m.begin(), m.end(), msum.begin());
  const double msum_all = msum[n - 1];
  int rejection_counter = -1;
  if (sqrts <= msum_all) {
    logg[LAction].error()
        << "An interaction requiring " << sqrts
        << "GeV was attempted below the minimum energy threshold" << msum_all
        << " GeV, but was ignored.\nThis is a known internal error which does "
           "not significantly affect physical results, and will be fixed in a "
           "near-future release.";
    throw StochasticBelowEnergyThreshold("Ignoring this action.");
  }

  double w, r01;
  std::vector<double> Minv(n);

  double weight_sqr_max = 1;
  const double Ekin_share = (sqrts - msum_all) / (n - 1);
  for (size_t i = 1; i < n; i++) {
    // This maximum estimate idea is due Scott Pratt: maximum should be
    // roughly at equal kinetic energies
    weight_sqr_max *= pCM_sqr(i * Ekin_share + msum[i],
                              (i - 1) * Ekin_share + msum[i - 1], m[i]);
  }
  // Maximum estimate is rough and can be wrong. We multiply it by additional
  // factor to be on the safer side.
  const double safety_factor = 1.1 + (n - 2) * 0.2;
  weight_sqr_max *= (safety_factor * safety_factor);
  bool first_warning = true;

  do {
    // Generate invariant masses of 1, 12, 123, 1243, etc.
    // Minv = {m1, M12, M123, ..., M123n-1, sqrts}
    Minv[0] = 0.0;
    Minv[n - 1] = sqrts - msum_all;
    for (size_t i = 1; i < n - 1; i++) {
      Minv[i] = random::uniform(0.0, sqrts - msum_all);
    }
    std::sort(Minv.begin(), Minv.end());
    for (size_t i = 0; i < n; i++) {
      Minv[i] += msum[i];
    }

    double weight_sqr = 1;
    for (size_t i = 1; i < n; i++) {
      weight_sqr *= pCM_sqr(Minv[i], Minv[i - 1], m[i]);
    }

    rejection_counter++;
    r01 = random::canonical();
    w = weight_sqr / weight_sqr_max;
    if (w > 1.0) {
      logg[LAction].warn()
          << "sample_manybody_phasespace_impl: alarm, weight > 1, w^2 = " << w
          << ". Increase safety factor." << std::endl;
    }
    if (rejection_counter > 20 && first_warning) {
      logg[LAction].warn() << "sample_manybody_phasespace_impl: "
                           << "likely hanging, way too many rejections,"
                           << " n = " << n << ", sqrts = " << sqrts
                           << ", msum = " << msum_all;
      first_warning = false;
    }
  } while (w < r01 * r01);

  // Boost particles to the right frame
  std::vector<ThreeVector> beta(n);
  for (size_t i = n - 1; i > 0; i--) {
    const double pcm = pCM(Minv[i], Minv[i - 1], m[i]);
    Angles phitheta;
    phitheta.distribute_isotropically();
    const ThreeVector isotropic_unitvector = phitheta.threevec();
    sampled_momenta[i] = FourVector(std::sqrt(m[i] * m[i] + pcm * pcm),
                                    pcm * isotropic_unitvector);
    if (i >= 2) {
      beta[i - 2] = pcm * isotropic_unitvector /
                    std::sqrt(pcm * pcm + Minv[i - 1] * Minv[i - 1]);
    }
    if (i == 1) {
      sampled_momenta[0] = FourVector(std::sqrt(m[0] * m[0] + pcm * pcm),
                                      -pcm * isotropic_unitvector);
    }
  }

  for (size_t i = 0; i < n - 2; i++) {
    // After each boost except the last one the sum of 3-momenta should be 0
    FourVector ptot = FourVector(0.0, 0.0, 0.0, 0.0);
    for (size_t j = 0; j <= i + 1; j++) {
      ptot += sampled_momenta[j];
    }
    logg[LAction].debug() << "Total momentum of 0.." << i + 1 << " = "
                          << ptot.threevec() << " and should be (0, 0, 0). "
                          << std::endl;

    // Boost the first i+1 particles to the next CM frame
    for (size_t j = 0; j <= i + 1; j++) {
      sampled_momenta[j] = sampled_momenta[j].lorentz_boost(beta[i]);
    }
  }

  FourVector ptot_all = FourVector(0.0, 0.0, 0.0, 0.0);
  for (size_t j = 0; j < n; j++) {
    ptot_all += sampled_momenta[j];
  }
  logg[LAction].debug() << "Total 4-momentum = " << ptot_all << ", should be ("
                        << sqrts << ", 0, 0, 0)" << std::endl;
}

void Action::sample_manybody_phasespace() {
  const size_t n = outgoing_particles_.size();
  if (n < 3) {
    throw std::invalid_argument(
        "sample_manybody_phasespace: number of outgoing particles should be 3 "
        "or more");
  }
  bool all_stable = true;
  for (size_t i = 0; i < n; i++) {
    all_stable = all_stable && outgoing_particles_[i].type().is_stable();
  }
  if (!all_stable) {
    throw std::invalid_argument(
        "sample_manybody_phasespace: Found resonance in to be sampled outgoing "
        "particles, but assumes stable particles.");
  }

  std::vector<double> m(n);
  for (size_t i = 0; i < n; i++) {
    m[i] = outgoing_particles_[i].type().mass();
  }
  std::vector<FourVector> p(n);

  sample_manybody_phasespace_impl(sqrt_s(), m, p);
  for (size_t i = 0; i < n; i++) {
    outgoing_particles_[i].set_4momentum(p[i]);
  }
}

void Action::assign_unpolarized_spin_vector_to_outgoing_particles() {
  for (ParticleData &p : outgoing_particles_) {
    p.set_unpolarized_spin_vector();
  }
}

double Action::check_conservation(const uint32_t id_process) const {
  QuantumNumbers before(incoming_particles_);
  QuantumNumbers after(outgoing_particles_);
  double energy_violation = 0.;
  if (before != after) {
    std::stringstream particle_names;
    for (const auto &p : incoming_particles_) {
      particle_names << p.type().name();
    }
    particle_names << " vs. ";
    for (const auto &p : outgoing_particles_) {
      particle_names << p.type().name();
    }
    particle_names << "\n";
    std::string err_msg = before.report_deviations(after);
    /* Pythia does not conserve energy and momentum at high energy, so we just
     * print the warning and continue. */
    if ((is_string_soft_process(process_type_)) ||
        (process_type_ == ProcessType::StringHard)) {
      logg[LAction].warn() << "Conservation law violations due to Pythia\n"
                           << particle_names.str() << err_msg;
      energy_violation = after.momentum()[0] - before.momentum()[0];
      return energy_violation;
    }
    /* We allow decay of particles stable under the strong interaction to decay
     * at the end, so just warn about such a "weak" process violating
     * conservation laws */
    if (process_type_ == ProcessType::Decay &&
        incoming_particles_[0].type().is_stable()) {
      logg[LAction].warn()
          << "Conservation law violations of strong interaction in weak or "
             "e.m. decay\n"
          << particle_names.str() << err_msg;
      return energy_violation;
    }
    /* If particles are added or removed, it is not surprising that conservation
     * laws are potentially violated. Do not warn the user but print some
     * information for debug */
    if (process_type_ == ProcessType::Freeforall) {
      logg[LAction].debug()
          << "Conservation law violation, but we want it (Freeforall Action).\n"
          << particle_names.str() << err_msg;
      return energy_violation;
    }
    logg[LAction].error() << "Conservation law violations detected\n"
                          << particle_names.str() << err_msg;
    if (id_process == ID_PROCESS_PHOTON) {
      throw std::runtime_error("Conservation laws violated in photon process");
    } else {
      throw std::runtime_error("Conservation laws violated in process " +
                               std::to_string(id_process));
    }
  }
  return energy_violation;
}

std::ostream &operator<<(std::ostream &out, const ActionList &actions) {
  out << "ActionList {\n";
  for (const auto &a : actions) {
    out << "- " << a << '\n';
  }
  return out << '}';
}

}  // namespace smash
