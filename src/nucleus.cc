/*
 *    Copyright (c) 2014-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include "smash/nucleus.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>

#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/fourvector.h"
#include "smash/input_keys.h"
#include "smash/logging.h"
#include "smash/numerics.h"
#include "smash/particles.h"
#include "smash/random.h"
#include "smash/threevector.h"

namespace smash {
static constexpr int LNucleus = LogArea::Nucleus::id;

Nucleus::Nucleus(const std::map<PdgCode, int> &particle_list, int nTest) {
  fill_from_list(particle_list, nTest);
  make_nucleus_unpolarized();
  set_parameters_automatic();
  set_saturation_density(calculate_saturation_density());
}

Nucleus::Nucleus(Configuration &config, int nTest) {
  assert(has_projectile_or_target(config));
  const bool is_projectile = is_about_projectile(config);
  const auto &[particles_key, diffusiveness_key, radius_key,
               saturation_key] = [&is_projectile]() {
    return is_projectile
               ? std::make_tuple(
                     InputKeys::modi_collider_projectile_particles,
                     InputKeys::modi_collider_projectile_diffusiveness,
                     InputKeys::modi_collider_projectile_radius,
                     InputKeys::modi_collider_projectile_saturationDensity)
               : std::make_tuple(
                     InputKeys::modi_collider_target_particles,
                     InputKeys::modi_collider_target_diffusiveness,
                     InputKeys::modi_collider_target_radius,
                     InputKeys::modi_collider_target_saturationDensity);
  }();
  // Fill nuclei with particles.
  std::map<PdgCode, int> part = config.take(particles_key);
  fill_from_list(part, nTest);
  make_nucleus_unpolarized();
  // Look for user-defined values or take the default parameters.
  const bool is_diffusiveness_given = config.has_value(diffusiveness_key),
             is_radius_given = config.has_value(radius_key),
             is_saturation_given = config.has_value(saturation_key);
  if (is_diffusiveness_given && is_radius_given && is_saturation_given) {
    diffusiveness_ = config.take(diffusiveness_key);
    nuclear_radius_ = config.take(radius_key);
    saturation_density_ = config.take(saturation_key);
  } else if (!is_diffusiveness_given && !is_radius_given &&
             !is_saturation_given) {
    set_parameters_automatic();
    set_saturation_density(calculate_saturation_density());
  } else {
    throw std::invalid_argument(
        "Diffusiveness, Radius and Saturation_Density required to manually "
        "configure the Woods-Saxon distribution. Only one or two were provided."
        "\nProviding none of the above mentioned parameters automatically "
        "configures the distribution based on the atomic number.");
  }
}

double Nucleus::mass() const {
  double total_mass = 0.;
  for (auto i = cbegin(); i != cend(); i++) {
    total_mass += i->momentum().abs();
  }
  return total_mass / (testparticles_ + 0.0);
}

/**
 * <h1> Woods-Saxon-distribution </h1>
 *
 * <h2> The distribution </h2>
 *
 * Nucleons in nuclei are distributed according to a
 * Woods-Saxon-distribution (see \iref{Woods:1954zz})
 *
 * \f[\frac{dN}{d^3r} = \frac{\rho_0}{\exp\left(\frac{r-r_0}{d}\right)
 * +1},\f]
 *
 * where \f$d\f$ is the \em diffusiveness of the nucleus. For \f$d=0\f$,
 * the nucleus is a hard sphere.  \f$\rho_0\f$ and \f$r_0\f$ are, in
 * this limit, the nuclear ground state density and
 * nuclear radius, respectively. For small \f$d\f$, this is still
 * approximately true.
 *
 * This distribution is obviously spherically symmetric, hence we can
 * rewrite \f$d^3r = 4\pi r^2 dr\f$ and obtain
 *
 * \f[\frac{dN}{4\pi\rho_0dr} =
 * \frac{r^2}{\exp\left(\frac{r-r_0}{d}\right) + 1}.\f]
 *
 * Let us rewrite that in units of \f$d\f$ (that's the diffusiveness)
 * and drop any constraints on normalization (since in the end we only
 * care about relative probabilities: we create as many nucleons as we
 * need). Now, \f$p(B)\f$ is the un-normalized probability to obtain a
 * point at \f$r = Bd\f$ (with \f$R = r_0/d\f$):
 *
 * \f[p(B) = \frac{B^2}{\exp(B-R) + 1}.\f]
 *
 * <h2> Splitting it up in two regimes </h2>
 *
 * We shift the distribution so that \f$B-R\f$ is 0 at \f$t = 0\f$: \f$t
 * = B-R\f$:
 *
 * \f[p^{(1)}(t)= \frac{(t+R)^2}{\exp(t)+1}\f]
 *
 * and observe
 *
 * \f[\frac{1}{\exp(x)+1} = \frac{e^{-x}}{e^{-x}e^{x}+e^{-x}} =
 * \frac{e^{-x}}{e^{-x}+1}.\f]
 *
 * The distribution function can now be split into two cases. For
 * negative t (first case), \f$-|t| = t\f$, and for positive t (second
 * case), \f$-|t| = -t\f$:
 *
 * \f[p^{(1)}(t) = \frac{1}{e^{-|t|}+1} \cdot (t+R)^2 \cdot
 * \begin{cases}
 * 1 & -R \le t < 0 \\
 * e^{-t} & t \ge 0
 * \end{cases}.\f]
 *
 * Apart from the first term, all that remains here can easily and
 * exactly be generated from unrejected uniform random numbers (see
 * below). The first term itself - \f$(1+e^{-|t|})^{-1}\f$ - is a number
 * between 1/2 and 1.
 *
 * If we now have a variable \f$t\f$ distributed according to the
 * remainder, \f$p^{(2)}(t)\f$, and reject \f$t\f$ with a probability
 * \f$p^{(rej)}(t) = 1 - p^{(survive)}(t) = 1 - (1+e^{-|t|})^{-1}\f$,
 * the resulting distribution is \f$p^{(combined)}(t) = p^{(2)}(t) \cdot
 * p^{(survive)}(t)\f$. Hence, we need to generate \f$p^{(2)}(t)\f$,
 * which we can normalize to
 *
 * \f[\tilde{p}^{(2)}(t) = \frac{1}{1+3/R+6/R^2+6/R^3} \cdot \begin{cases}
 * \frac{3}{R^3} (t+R)^2 & -R \le t < 0 \\
 * e^{-t} \left( \frac{3}{R}+\frac{6}{R^2}t+\frac{6}{R^3}\frac{1}{2}t^2 \right)
 * & t \ge 0 \end{cases}.\f]
 *
 * (the tilde \f$\tilde{p}\f$ means that this is normalized).
 *
 * <h2> Four parts inside the rejection </h2>
 *
 * Let \f$c_1 = 1+3/R+6/R^2+6/R^3\f$. The above means:
 *
 * \f[\mbox{Choose: } \begin{cases}
 * \tilde p^{({\rm I})} = \frac{3}{R^3}(t+R)^2 \Theta(-t) \Theta(t+R) \\
 * \tilde p^{({\rm II})}= e^{-t}\Theta(t) \\
 * \tilde p^{({\rm III})}=e^{-t}\Theta(t) t \\
 * \tilde p^{({\rm IV})} =e^{-t}\Theta(t) \frac{1}{2} t^2
 * \end{cases} \mbox{ with a probability of }\begin{cases}
 * \frac{1}{c_1} \cdot 1 \\
 * \frac{1}{c_1} \cdot \frac{3}{R} \\
 * \frac{1}{c_1} \cdot \frac{6}{R^2} \\
 * \frac{1}{c_1} \cdot \frac{6}{R^3}
 * \end{cases}.\f]
 *
 * Let us see how those are generated. \f$\chi_i\f$ are uniformly
 * distributed numbers between 0 and 1.
 *
 * \f[p(\chi_i) = \Theta(\chi_i)\Theta(1-\chi_i)\f]
 *
 * For simple distributions (only one \f$\chi\f$ involved), we invert
 * \f$t(\chi)\f$, derive it w.r.t. \f$t\f$ and normalize.
 *
 * <h3> Case I: \f$p^{({\rm I})}\f$ </h3>
 *
 * Simply from one random number:
 *
 * \f[t = R\left( \sqrt[ 3 ]{\chi} - 1 \right)\f]
 * \f[\tilde p^{({\rm I})} = \frac{3}{R^3}(t+R)^2 \mbox{ for } -R \le t
 * \le 0\f]
 *
 * <h3> Case II: \f$p^{({\rm II})}\f$ </h3>
 *
 * Again, from one only:
 *
 * \f[t = -\log(\chi)\f]
 * \f[p(t) = \frac{d\chi}{dt}\f]
 * \f[p^{({\rm II})} = e^{-t} \mbox{ for } t > 0\f]
 *
 * <h3> Case III: \f$p^{({\rm III})}\f$ </h3>
 *
 * Here, we need two variables:
 *
 * \f[t = -\log{\chi_1} -\log{\chi_2}\f]
 *
 * \f$p^{({\rm III})}\f$ is now the folding of \f$p^{({\rm II})}\f$ with
 *itself[1]:
 *
 * \f[p^{({\rm III})} = \int_{-\infty}^{\infty} d\tau e^{-\tau} e^{-(t-\tau)}
 * \Theta(\tau) \Theta(t-\tau) = t e^{-t} \mbox{ for } t > 0\f]
 *
 * <h3> Case IV: \f$p^{({\rm IV})}\f$ </h3>
 *
 * Three variables needed:
 *
 * \f[t = -\log{\chi_1} -\log{\chi_2} -\log{\chi_3}\f]
 *
 * \f$p^{({\rm IV})}\f$ is now the folding of \f$p^{({\rm II})}\f$ with
 * \f$p^{({\rm III})}\f$:
 *
 * \f[p^{({\rm IV})} = \int_{ - \infty}^{\infty} d\tau e^{- \tau} \left(
 * t - \tau \right) e^{ - (t - \tau)} \Theta(\tau) \Theta(t - \tau) =
 * \frac{1}{2} t^2 e^{ -t} \mbox{ for } t > 0\f]
 *
 * [1]: This is [the probability to find a \f$\tau\f$] times [the
 * probability to find the value \f$\tau_2 = t-\tau\f$ that added to
 * \f$\tau\f$ yields \f$t\f$], integrated over all possible combinations
 * that have that property.
 *
 * <h2> From the beginning </h2>
 *
 *  So, the algorithm needs to do all this from the end:
 *
 */
ThreeVector Nucleus::distribute_nucleon() {
  // Get the solid angle of the nucleon.
  Angles dir;
  dir.distribute_isotropically();
  // diffusiveness_ zero or negative? Use hard sphere.
  if (almost_equal(diffusiveness_, 0.)) {
    return dir.threevec() * nuclear_radius_ * std::cbrt(random::canonical());
  }
  if (almost_equal(nuclear_radius_, 0.)) {
    return smash::ThreeVector();
  }
  double radius_scaled = nuclear_radius_ / diffusiveness_;
  double prob_range1 = 1.0;
  double prob_range2 = 3. / radius_scaled;
  double prob_range3 = 2. * prob_range2 / radius_scaled;
  double prob_range4 = 1. * prob_range3 / radius_scaled;
  double ranges234 = prob_range2 + prob_range3 + prob_range4;
  double t;
  /// \li Decide which branch \f$\tilde p^{({\rm I - IV})}\f$ to go into
  do {
    double which_range = random::uniform(-prob_range1, ranges234);
    if (which_range < 0.0) {
      t = radius_scaled * (std::cbrt(random::canonical()) - 1.);
    } else {
      t = -std::log(random::canonical());
      if (which_range >= prob_range2) {
        t -= std::log(random::canonical());
        if (which_range >= prob_range2 + prob_range3) {
          t -= std::log(random::canonical());
        }
      }
    }
    /**
     * \li Generate \f$t\f$ from the distribution in the respective
     * branches
     * \li \a Reject that number with a probability
     * \f$1-(1+\exp(-|t|))^{-1}\f$ (the efficiency of this should be
     * \f$\gg \frac{1}{2}\f$)
     */
  } while (random::canonical() > 1. / (1. + std::exp(-std::abs(t))));
  /// \li Shift and rescale \f$t\f$ to \f$r = d\cdot t + r_0\f$
  double position_scaled = t + radius_scaled;
  double position = position_scaled * diffusiveness_;
  return dir.threevec() * position;
}

double Nucleus::woods_saxon(double r) {
  return r * r / (std::exp((r - nuclear_radius_) / diffusiveness_) + 1);
}

void Nucleus::arrange_nucleons() {
  for (auto i = begin(); i != end(); i++) {
    // Initialize momentum
    i->set_4momentum(i->pole_mass(), 0.0, 0.0, 0.0);
    /* Sampling the Woods-Saxon, get the radial
     * position and solid angle for the nucleon. */
    ThreeVector pos = distribute_nucleon();

    // Set the position of the nucleon.
    i->set_4position(FourVector(0.0, pos));
  }

  // Recenter and rotate
  align_center();
  rotate();
}

void Nucleus::set_parameters_automatic() {
  int A = Nucleus::number_of_particles();
  int Z = Nucleus::number_of_protons();
  if (A == 1) {  // single particle
    /* In case of testparticles, an infinite reaction loop will be
     * avoided by a small finite spread according to a single particles
     * 'nucleus'. The proper solution will be to introduce parallel
     * ensembles. */
    set_nuclear_radius(
        testparticles_ == 1 ? 0. : 1. - std::exp(-(testparticles_ - 1.) * 0.1));
    set_diffusiveness(testparticles_ == 1 ? -1. : 0.02);
  } else if ((A == 238) && (Z == 92)) {  // Uranium
    // Default values.
    set_diffusiveness(0.556);
    set_nuclear_radius(6.86);
  } else if ((A == 208) && (Z == 82)) {  // Lead
    // Default values.
    set_diffusiveness(0.54);
    set_nuclear_radius(6.67);
  } else if ((A == 197) && (Z == 79)) {  // Gold
    // Default values from \iref{Schopper:2004qco}
    set_diffusiveness(0.523);
    set_nuclear_radius(6.55);
  } else if ((A == 129) && (Z == 54)) {  // Xenon
    // Default values.
    set_diffusiveness(0.59);
    set_nuclear_radius(5.36);
  } else if ((A == 63) && (Z == 29)) {  // Copper
    // Default values.
    set_diffusiveness(0.5977);
    set_nuclear_radius(4.20641);
  } else if (A == 96) {
    if (Z == 40) {  // Zirconium
      // Default values.
      set_diffusiveness(0.46);
      set_nuclear_radius(5.02);
    } else if (Z == 44) {  // Ruthenium
      // Default values.
      set_diffusiveness(0.46);
      set_nuclear_radius(5.085);
    } else {
      // radius and diffusiveness taken from \iref{Rybczynski:2013yba}
      set_diffusiveness(0.54);
      set_nuclear_radius(1.12 * std::pow(A, 1.0 / 3.0) -
                         0.86 * std::pow(A, -1.0 / 3.0));
    }
  } else {
    // saturation density already has reasonable default
    set_nuclear_radius(default_nuclear_radius());
    if (A <= 16) {
      set_diffusiveness(0.545);
    } else {
      // diffusiveness taken from \iref{Rybczynski:2013yba}
      set_diffusiveness(0.54);
    }
  }
}

void Nucleus::set_orientation_from_config(Configuration &config) {
  const bool is_projectile =
      has_projectile_or_target(config) ? is_about_projectile(config) : true;
  const auto &[rotation_key, theta_key, phi_key, psi_key] = [&is_projectile]() {
    return is_projectile
               ? std::make_tuple(
                     InputKeys::modi_collider_projectile_orientation_randRot,
                     InputKeys::modi_collider_projectile_orientation_theta,
                     InputKeys::modi_collider_projectile_orientation_phi,
                     InputKeys::modi_collider_projectile_orientation_psi)
               : std::make_tuple(
                     InputKeys::modi_collider_target_orientation_randRot,
                     InputKeys::modi_collider_target_orientation_theta,
                     InputKeys::modi_collider_target_orientation_phi,
                     InputKeys::modi_collider_target_orientation_psi);
  }();
  const bool was_any_angle_provided = config.has_value(theta_key) ||
                                      config.has_value(phi_key) ||
                                      config.has_value(psi_key);
  random_rotation_ = config.take(rotation_key);
  if (random_rotation_ && was_any_angle_provided) {
    throw std::domain_error(
        "The random rotation of nuclei has been requested, but some specific "
        "rotation angle is provided, too. Please specify only either of them.");
  } else {
    euler_theta_ = config.take(theta_key);
    euler_phi_ = config.take(phi_key);
    euler_psi_ = config.take(psi_key);
  }
}

void Nucleus::rotate() {
  if (random_rotation_) {
    // Randomly generate euler angles for theta and phi. Psi needs not be
    // assigned, as the nucleus objects are symmetric with respect to psi.
    random_euler_angles();
  }
  if (euler_phi_ != 0.0 || euler_theta_ != 0.0 || euler_psi_ != 0.0) {
    for (auto &particle : *this) {
      /* Rotate every vector by the euler angles phi, theta and psi.
       * This means applying the matrix for a rotation of phi around the z-axis,
       * followed by the matrix for a rotation of theta around the rotated
       * x-axis and the matrix for a rotation of psi around the rotated z-axis.
       */
      ThreeVector three_pos = particle.position().threevec();
      three_pos.rotate(euler_phi_, euler_theta_, euler_psi_);
      particle.set_3position(three_pos);
    }
  }
}

void Nucleus::generate_fermi_momenta() {
  const int N_n = std::count_if(begin(), end(), [](const ParticleData i) {
    return i.pdgcode() == pdg::n;
  });
  const int N_p = std::count_if(begin(), end(), [](const ParticleData i) {
    return i.pdgcode() == pdg::p;
  });
  const FourVector nucleus_center = center();
  const int A = N_n + N_p;
  constexpr double pi2_3 = 3.0 * M_PI * M_PI;
  logg[LNucleus].debug() << N_n << " neutrons, " << N_p << " protons.";

  ThreeVector ptot = ThreeVector(0.0, 0.0, 0.0);
  for (auto i = begin(); i != end(); i++) {
    // Only protons and neutrons get Fermi momenta
    if (i->pdgcode() != pdg::p && i->pdgcode() != pdg::n) {
      if (i->is_baryon()) {
        logg[LNucleus].warn() << "No rule to calculate Fermi momentum "
                              << "for particle " << i->pdgcode();
      }
      continue;
    }
    const double r = (i->position() - nucleus_center).abs3();
    const double theta = (i->position().threevec().get_theta());
    const double phi = (i->position().threevec().get_phi());
    double rho = nucleon_density(r, std::cos(theta), phi);

    if (i->pdgcode() == pdg::p) {
      rho = rho * N_p / A;
    }
    if (i->pdgcode() == pdg::n) {
      rho = rho * N_n / A;
    }
    const double p =
        hbarc * std::pow(pi2_3 * rho * random::uniform(0.0, 1.0), 1.0 / 3.0);
    Angles phitheta;
    phitheta.distribute_isotropically();
    const ThreeVector ith_3momentum = phitheta.threevec() * p;
    ptot += ith_3momentum;
    i->set_3momentum(ith_3momentum);
    logg[LNucleus].debug() << "Particle: " << *i << ", pF[GeV]: "
                           << hbarc * std::pow(pi2_3 * rho, 1.0 / 3.0)
                           << " r[fm]: " << r
                           << " Nuclear radius[fm]: " << nuclear_radius_;
  }
  if (A == 0) {
    // No Fermi momenta should be assigned
    assert(ptot.x1() == 0.0 && ptot.x2() == 0.0 && ptot.x3() == 0.0);
  } else {
    /* Ensure zero total momentum of nucleus - redistribute ptot equally
     * among protons and neutrons */
    const ThreeVector centralizer = ptot / A;
    for (auto i = begin(); i != end(); i++) {
      if (i->pdgcode() == pdg::p || i->pdgcode() == pdg::n) {
        i->set_4momentum(i->pole_mass(),
                         i->momentum().threevec() - centralizer);
      }
    }
  }
}

void Nucleus::boost(double beta_scalar) {
  double beta_squared = beta_scalar * beta_scalar;
  double one_over_gamma = std::sqrt(1.0 - beta_squared);
  double gamma = 1.0 / one_over_gamma;
  /* We are talking about a /passive/ lorentz transformation here, as
   * far as I can see, so we need to boost in the direction opposite to
   * where we want to go
   *     ( The vector we transform - p - stays unchanged, but we go into
   *       a system that moves with -beta. Now in this frame, it seems
   *       like p has been accelerated with +beta.
   *     ) */
  for (auto i = begin(); i != end(); i++) {
    /* a real Lorentz Transformation would leave the particles at
     * different times here, which we would then have to propagate back
     * to equal times. Since we know the result, we can simply multiply
     * the z-value with 1/gamma. */
    FourVector this_position = i->position();
    this_position.set_x3(this_position.x3() * one_over_gamma);
    i->set_4position(this_position);
    /* The simple Lorentz transformation of momenta does not take into account
     * that nucleus has binding energy. Here we apply the method used
     * in the JAM code \iref{Nara:1999dz}: p' = p_beam + gamma*p_F.
     * This formula is derived under assumption that all nucleons have
     * the same binding energy. */
    FourVector mom_i = i->momentum();
    i->set_4momentum(i->pole_mass(), mom_i.x1(), mom_i.x2(),
                     gamma * (beta_scalar * mom_i.x0() + mom_i.x3()));
  }
}

void Nucleus::fill_from_list(const std::map<PdgCode, int> &particle_list,
                             int testparticles) {
  testparticles_ = testparticles;
  for (auto n = particle_list.cbegin(); n != particle_list.cend(); ++n) {
    const ParticleType &current_type = ParticleType::find(n->first);
    double current_mass = current_type.mass();
    for (unsigned int i = 0; i < n->second * testparticles_; i++) {
      // append particle to list and set its PDG code.
      particles_.emplace_back(current_type);
      particles_.back().set_4momentum(current_mass, 0.0, 0.0, 0.0);
    }
  }
}

void Nucleus::shift(double z_offset, double x_offset, double simulation_time) {
  // Move the nucleus in z and x directions, and set the time.
  for (auto i = begin(); i != end(); i++) {
    FourVector this_position = i->position();
    this_position.set_x3(this_position.x3() + z_offset);
    this_position.set_x1(this_position.x1() + x_offset);
    this_position.set_x0(simulation_time);
    i->set_4position(this_position);
    i->set_formation_time(simulation_time);
  }
}

void Nucleus::copy_particles(Particles *external_particles) {
  for (auto p = begin(); p != end(); p++) {
    external_particles->insert(*p);
  }
}

FourVector Nucleus::center() const {
  FourVector centerpoint(0.0, 0.0, 0.0, 0.0);
  for (auto p = cbegin(); p != cend(); p++) {
    centerpoint += p->position();
  }
  centerpoint /= size();
  return centerpoint;
}

void Nucleus::random_euler_angles() {
  // Sample euler_theta_ such that cos(theta) is uniform
  euler_phi_ = twopi * random::uniform(0., 1.);
  euler_theta_ = std::acos(2 * random::uniform(0., 1.) - 1);
  euler_psi_ = twopi * random::uniform(0., 1.);
}

double Nucleus::nucleon_density(double r, double, double) const {
  return get_saturation_density() /
         (std::exp((r - nuclear_radius_) / diffusiveness_) + 1.);
}

double Nucleus::nucleon_density_unnormalized(double r, double, double) const {
  return 1.0 / (std::exp((r - nuclear_radius_) / diffusiveness_) + 1.);
}

double Nucleus::calculate_saturation_density() const {
  Integrator2d integrate;
  // Transform integral from (0, oo) to (0, 1) via r = (1 - t) / t.
  // To prevent overflow, the integration is only performed to t = 0.01 which
  // corresponds to r = 99fm. Additionally the precision settings in the
  // Integrator2d scheme are equally important. However both these point affect
  // the result only after the seventh digit which should not be relevant here.
  const auto result = integrate(0.01, 1, -1, 1, [&](double t, double cosx) {
    const double r = (1 - t) / t;
    return twopi * std::pow(r, 2.0) *
           nucleon_density_unnormalized(r, cosx, 0.0) / std::pow(t, 2.0);
  });
  const auto rho0 = number_of_particles() / result.value();
  return rho0;
}

std::ostream &operator<<(std::ostream &out, const Nucleus &n) {
  return out << "  #particles   #testparticles   mass [GeV]   "
                "radius [fm]  diffusiveness [fm]\n"
             << format(n.number_of_particles(), nullptr, 12)
             << format(n.size(), nullptr, 17) << format(n.mass(), nullptr, 13)
             << format(n.get_nuclear_radius(), nullptr, 14)
             << format(n.get_diffusiveness(), nullptr, 20);
}

bool has_projectile_or_target(const Configuration &config) {
  const bool is_projectile = config.has_section(InputSections::m_c_projectile);
  const bool is_target = config.has_section(InputSections::m_c_target);
  return is_projectile || is_target;
}

bool is_about_projectile(const Configuration &config) {
  const bool is_projectile = config.has_section(InputSections::m_c_projectile);
  const bool is_target = config.has_section(InputSections::m_c_target);
  if (is_projectile == is_target) {
    throw std::logic_error(
        "Error parsing configuration of EITHER projectile OR target.\n"
        "Configuration tested for it contains the following:\n------------\n" +
        config.to_string() + "\n------------\n");
  }
  return is_projectile;
}

}  // namespace smash
