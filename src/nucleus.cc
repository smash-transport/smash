/*
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include "include/nucleus.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>

#include "include/angles.h"
#include "include/constants.h"
#include "include/logging.h"
#include "include/numerics.h"
#include "include/particles.h"
#include "include/pdgcode.h"
#include "include/threevector.h"

namespace Smash {

Nucleus::Nucleus(const std::map<PdgCode, int>& particle_list, int nTest) {
  fill_from_list(particle_list, nTest);
}

Nucleus::Nucleus(Configuration &config, int nTest) {
  // Fill nuclei with particles.
  std::map<PdgCode, int> part = config.take({"Particles"});
  fill_from_list(part, nTest);
  // Ask to construct nuclei based on atomic number;
  // otherwise, look for the user-defined values or take the default parameters.
  if (config.has_value({"Automatic"}) && config.take({"Automatic"})) {
    set_parameters_automatic();
  } else {
    set_parameters_from_config(config);
  }
}

float Nucleus::mass() const {
  float total_mass = 0.f;
  for (auto i = cbegin(); i != cend(); i++) {
    total_mass += i->momentum().abs();
  }
  return total_mass/(testparticles_+0.0);
}

/**
 * Woods-Saxon-distribution
 * ========================
 *
 * The distribution
 * ----------------
 *
 *
 * Nucleons in nuclei are distributed according to a
 * Woods-Saxon-distribution (see \iref{Woods:1954zz})
 *
 * \f[\frac{dN}{d^3r} = \frac{\rho_0}{\exp\left(\frac{r-r_0}{d}\right)
 * +1},\f]
 *
 * where \f$d\f$ is the \em diffuseness of the nucleus. For \f$d=0\f$,
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
 * care about relative probabilities: we create as many nuclei as we
 * need). Now, \f$p(B)\f$ is the un-normalized probability to obtain a
 * point at \f$r = Bd\f$ (with \f$R = r_0/d\f$):
 *
 * \f[p(B) = \frac{B^2}{\exp(B-R) + 1}.\f]
 *
 * Splitting it up in two regimes
 * ------------------------------
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
 * Four parts inside the rejection
 * -------------------------------
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
 * \end{cases}.\]
 *
 * Let us see how those are generated. \f$\chi_i\f$ are uniformly
 * distributed numbers between 0 and 1.
 *
 * \f[p(\chi_i) = \Theta(\chi_i)\Theta(1-\chi_i)\f]
 *
 * For simple distributions (only one \f$\chi\f$ involved), we invert
 * \f$t(\chi)\f$, derive it w.r.t. \f$t\f$ and normalize.
 *
 * ### Case I: \f$p^{({\rm I})}\f$
 *
 * Simply from one random number:
 *
 * \f[t = R\left( \sqrt[ 3 ]{\chi} - 1 \right)\f]
 * \f[\tilde p^{({\rm I})} = \frac{3}{R^3}(t+R)^2 \mbox{ for } -R \le t
 * \le 0\f]
 *
 * ### Case II: \f$p^{({\rm II})}\f$
 *
 * Again, from one only:
 *
 * \f[t = -\log(\chi)\f]
 * \f[p(t) = \frac{d\chi}{dt}\f]
 * \f[p^{({\rm II})} = e^{-t} \mbox{ for } t > 0\f]
 *
 * ### Case III: \f$p^{({\rm III})}\f$
 *
 * Here, we need two variables:
 *
 * \f[t = -\log{\chi_1} -\log{\chi_2}\f]
 *
 * \f$p^{({\rm III})}\f$ is now the folding of \f$p^{({\rm II})}\f$ with itself[1]:
 *
 * \f[p^{({\rm III})} = \int_{-\infty}^{\infty} d\tau e^{-\tau} e^{-(t-\tau)}
 * \Theta(\tau) \Theta(t-\tau) = t e^{-t} \mbox{ for } t > 0\f]
 *
 * ### Case IV: \f$p^{({\rm IV})}\f$
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
 * \fpPrecision
 * Why does the `do-while` loop require double-precision?
 *
 * From the beginning
 * ------------------
 *
 *  So, the algorithm needs to do all this from the end:
 *
 **/
ThreeVector Nucleus::distribute_nucleon() const {
  // Get the solid angle of the nucleon.
  Angles dir;
  dir.distribute_isotropically();
  // diffusiveness_ zero or negative? Use hard sphere.
  if (almost_equal(diffusiveness_, 0.f)) {
    return dir.threevec() * nuclear_radius_ * std::cbrt(Random::canonical());
  }
  if (almost_equal(nuclear_radius_, 0.f)) {
    return Smash::ThreeVector();
  }
  float radius_scaled = nuclear_radius_/diffusiveness_;
  float prob_range1 = 1.0;
  float prob_range2 = 3. / radius_scaled;
  float prob_range3 = 2. * prob_range2 / radius_scaled;
  float prob_range4 = 1. * prob_range3 / radius_scaled;
  float ranges234 = prob_range2 + prob_range3 + prob_range4;
  float t;
  /// \li Decide which branch \f$\tilde p^{({\rm I - IV})}\f$ to go into
  do {
    float which_range = Random::uniform(-prob_range1, ranges234);
    if (which_range < 0.0) {
      t = radius_scaled * (std::cbrt(Random::canonical()) - 1.);
    } else {
      t = -std::log(Random::canonical());
      if (which_range >= prob_range2) {
        t -= std::log(Random::canonical());
        if (which_range >= prob_range2 + prob_range3) {
          t -= std::log(Random::canonical());
        }
      }
    }
    /** \li Generate \f$t\f$ from the distribution in the respective
     * branches
     * \li \a reject that number with a probability
     * \f$1-(1+\exp(-|t|))^{-1}\f$ (the efficiency of this should be
     * \f$\gg \frac{1}{2}\f$)
     **/
  } while (Random::canonical() > 1. / (1. + std::exp(-std::abs(t))));
  /// \li shift and rescale \f$t\f$ to \f$r = d\cdot t + r_0\f$
  float position_scaled = t + radius_scaled;
  float position = position_scaled * diffusiveness_;
  return dir.threevec() * position;
}

float Nucleus::woods_saxon(float r) {
  return r * r / (std::exp((r - nuclear_radius_) / diffusiveness_) + 1);
}

void Nucleus::arrange_nucleons() {
  for (auto i = begin(); i != end(); i++) {
    // Initialize momentum
    i->set_4momentum(i->pole_mass(), 0.0, 0.0, 0.0);
    // Sampling the W.S., get the radial
    // position and solid angle for the nucleon.
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
  switch (A) {
    case 1:  // single particle
      set_nuclear_radius(0.);
      set_diffusiveness(-1.);
      break;
    case 238:  // Uranium
      // Default values.
      set_diffusiveness(0.556);
      set_nuclear_radius(6.86);
      set_saturation_density(0.166);
      // Hirano, Huovinen, Nara - Corrections.
      // set_diffusiveness(0.44);
      // set_nuclear_radius(6.86);
      break;
    case 208:  // Lead
      // Default values.
      set_diffusiveness(0.54);
      set_nuclear_radius(6.67);
      set_saturation_density(0.161);
      break;
    case 197:  // Gold
      // Default values.
      set_diffusiveness(0.535);
      set_nuclear_radius(6.38);
      set_saturation_density(0.1695);
      // Hirano, Nara - Corrections.
      // set_diffusiveness(0.44);
      // set_nuclear_radius(6.42);
      break;
    case 63:  // Copper
      // Default values.
      set_diffusiveness(0.597);
      set_nuclear_radius(4.20641);
      set_saturation_density(0.1686);
      // Hirano, Nara - Corrections.
      // set_diffusiveness(0.50);
      // set_nuclear_radius(4.28);
      break;
    default:
      // radius: rough guess for all nuclei not listed explicitly
      set_nuclear_radius(1.2*std::cbrt(A));
      // diffusiveness and saturation density already have reasonable defaults
  }
}

void Nucleus::set_parameters_from_config(Configuration &config) {
  // Diffusiveness
  if (config.has_value({"Diffusiveness"})) {
    set_diffusiveness(static_cast<float>(config.take({"Diffusiveness"})));
  }
  // Radius
  if (config.has_value({"Radius"})) {
    set_nuclear_radius(static_cast<float>(config.take({"Radius"})));
  } else {
    set_nuclear_radius(default_nuclear_radius());
  }
}

void Nucleus::generate_fermi_momenta() {
  double r, rho, p;
  const int N_n = std::count_if(begin(), end(),
                  [](const ParticleData i) {return i.pdgcode() == 0x2112;});
  const int N_p = std::count_if(begin(), end(),
                  [](const ParticleData i) {return i.pdgcode() == 0x2212;});
  const FourVector nucleus_center = center();
  const int A = N_n + N_p;
  const double pi2_3 = 3.0 * M_PI * M_PI;
  Angles phitheta;
  const auto &log = logger<LogArea::Nucleus>();

  log.debug() << N_n << " neutrons, " << N_p << " protons.";

  ThreeVector ptot = ThreeVector(0.0, 0.0, 0.0);
  for (auto i = begin(); i != end(); i++) {
    // Only protons and neutrons get Fermi momenta
    if (i->pdgcode() != 0x2212 && i->pdgcode() != 0x2112) {
      if (i->is_baryon()) {
        log.error() << "No rule to calculate Fermi momentum " <<
                       "for particle " << i->pdgcode();
      }
      continue;
    }
    r = (i->position() - nucleus_center).abs3();
    rho = nuclear_density
          / (std::exp((r - nuclear_radius_)/diffusiveness_) + 1.);
    if (i->pdgcode() == 0x2212) {  // proton
      rho = rho * N_p / A;
    }
    if (i->pdgcode() == 0x2112) {  // neutron
      rho = rho * N_n / A;
    }
    p = hbarc * std::pow(pi2_3 * rho * Random::uniform(0.0, 1.0), 1.0/3.0);
    phitheta.distribute_isotropically();
    const ThreeVector ith_3momentum = phitheta.threevec() * p;
    ptot += ith_3momentum;
    i->set_3momentum(ith_3momentum);
    log.debug() << "Particle: " << *i <<
               ", pF[GeV]: " << hbarc * std::pow(pi2_3 * rho, 1.0/3.0) <<
               " r[fm]: " << r <<
               " Nuclear radius[fm]: " << nuclear_radius_;
  }
  if (A == 0) {
    // No Fermi momenta should be assigned
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wfloat-equal"
    assert(ptot.x1() == 0.0 && ptot.x2() == 0.0 && ptot.x3() == 0.0);
    #pragma GCC diagnostic pop
  } else {
    // Make sure that total momentum is zero - redistribute ptot equally among
    // protons and neutrons
    const ThreeVector centralizer = ptot/A;
    for (auto i = begin(); i != end(); i++) {
      if (i->pdgcode() == 0x2212 || i->pdgcode() == 0x2112) {
        i->set_4momentum(i->pole_mass(),
                         i->momentum().threevec() - centralizer);
      }
    }
  }
}

void Nucleus::boost(double beta_scalar) {
  double beta_squared = beta_scalar * beta_scalar;
  double one_over_gamma = std::sqrt(1.0 - beta_squared);
  /*double gamma = 1.0/one_over_gamma;
    double gammabeta = sign*sqrt(beta_squared)*gamma;
   */
  // We are talking about a /passive/ lorentz transformation here, as
  // far as I can see, so we need to boost in the direction opposite to
  // where we want to go
  //     ( The vector we transform - p - stays unchanged, but we go into
  //       a system that moves with -beta. Now in this frame, it seems
  //       like p has been accelerated with +beta.
  //     )
  ThreeVector beta(0., 0., -beta_scalar);
  for (auto i = begin(); i != end(); i++) {
    // a real Lorentz Transformation would leave the particles at
    // different times here, which we would then have to propagate back
    // to equal times. Since we know the result, we can simply multiply
    // the z-value with 1/gamma.
    FourVector this_position = i->position();
    this_position.set_x3(this_position.x3() * one_over_gamma);
    i->set_4position(this_position);
    // for momenta, though, we CAN do normal Lorentz Boosts, since we
    // *do* want to transform the zero-component (i.e., the energy).
    i->boost_momentum(beta);
  }
}

void Nucleus::fill_from_list(const std::map<PdgCode, int>& particle_list,
                             int testparticles) {
  testparticles_ = testparticles;
  for (auto n = particle_list.cbegin(); n != particle_list.cend(); ++n) {
    const ParticleType &current_type = ParticleType::find(n->first);
    float current_mass = current_type.mass();
    for (unsigned int i = 0; i < n->second*testparticles_; i++) {
      // append particle to list and set its PDG code.
      particles_.emplace_back(current_type);
      particles_.back().set_4momentum(current_mass, 0.0, 0.0, 0.0);
    }
  }
}

void Nucleus::shift(double z_offset,
                    double x_offset, float simulation_time) {
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

void Nucleus::copy_particles(Particles* external_particles) {
  for (auto p = begin(); p != end(); p++) {
    external_particles->insert(*p);
  }
}

/*void Nucleus::print_nucleus(const char * file_name) const {
  for (auto i = cbegin(); i != cend(); i++) {
    FourVector this_position = i->position();
    std::ofstream a_file;
    a_file.open(file_name, std::ios::app);
    a_file << std::to_string(this_position.x1()) + " " +
              std::to_string(this_position.x2()) + " " +
              std::to_string(this_position.x3()) << std::endl;
    a_file.close();
  }
}*/

FourVector Nucleus::center() const {
  FourVector centerpoint(0.0, 0.0, 0.0, 0.0);
  for (auto p = cbegin(); p != cend(); p++) {
    centerpoint += p->position();
  }
  centerpoint /= size();
  return centerpoint;
}

std::ostream &operator<<(std::ostream &out, const Nucleus &n) {
  return out << "  #particles   #testparticles   mass [GeV]   "
                "radius [fm]  diffusiveness [fm]\n"
             << format(n.number_of_particles(), nullptr, 12)
             << format(n.size(), nullptr, 17)
             << format(n.mass(), nullptr, 13)
             << format(n.get_nuclear_radius(), nullptr, 14)
             << format(n.get_diffusiveness(), nullptr, 20);
}

}  // namespace Smash
