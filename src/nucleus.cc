/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include <limits>

#include "include/nucleus.h"
#include "include/angles.h"

Nucleus::Nucleus() {}

float Nucleus::mass() const {
  float total_mass = 0.f;
  for (auto i = cbegin(); i != cend(); i++) {
    total_mass += sqrt(i->momentum().Dot());
  }
  return total_mass;
}

/**************************************************************************//**
 *
 * Woods-Saxon-distribution
 * ========================
 *
 * The distribution
 * ----------------
 * 
 *
 * Nucleons in nuclei are distributed according to a
 * Woods-Saxon-distribution[citation needed]
 *
 * \f[\frac{dN}{d^3r} = \frac{\rho_0}{\exp\left(\frac{r-r_0}{d}\right)
 * +1}.\f]
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
 * \f[\frac{1}{\exp(x)+1} = \frac{e^{-x}e^x}{e^{-x}e^{x}+e^{-x}} =
 * \frac{e^{-x}}{e^{-x}+1}.\f]
 *
 * The distribution function can now be split into two cases. For
 * negative t (first case), \f$-|t| = t\f$, and for positive t (second
 * case), \f$-|t| = -t\f$:
 *
 * \f[p^{(1)}(t) = \frac{1}{e^{-|t|}+1} \cdot (t+R)^2\begin{cases}
 * 1 & R \le t < 0 \\
 * e^{-t} & t \ge 0
 * \end{cases}.\f]
 *
 * Apart from the first term, all that remains here can easily and
 * exactly be generated from unrejected uniform random numbers (see
 * below). The first term itself - \f$(1+e^{-|t|})^{-1}\f$ - is a number
 * between 1/2 and 1.?
 *
 * If we now have a variable $t$ distributed according to
 * \f$p^{(2)}(t)\f$ and reject \f$t\f$ with a probability
 * \f$p^{(rej)}(t) = (1+e^{-|t|})^{-1}\f$, the resulting distribution is
 * \f$p^{(combined)}(t) = p^{(2)}(t) \cdot p^{(rej)}(t)\f$. Hence, what
 * we need to generate is (the tilde \f$\tilde p\f$ means that this is
 * normalized):
 *
 * \f[\tilde{p}^{(2)}(t) = \frac{1}{1+3/R+6/R^2+6/R^3} \cdot \begin{cases}
 * \frac{3}{R^3} (t+R)^2 & R \le t < 0 \\
 * e^{-x} \left( \frac{3}{R}+\frac{6}{R^2}t+\frac{6}{R^3}t^2 \right) & t
 * \ge 0 \end{cases}.\f]
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
 * From the beginning
 * ------------------
 *
 *  So, the algorithm needs to do all this from the end:
 *
 **/
float Nucleus::distribution_nucleons() const {
  // diffusiveness_ zero or negative? Use hard sphere.
  if (diffusiveness_ < std::numeric_limits<float>::min()) {
    return nuclear_radius()*(pow(drand48(), 1./3.));
  }
  float radius_scaled = nuclear_radius()/diffusiveness_;
  float prob_range1 = 1.0;
  float prob_range2 = 3. / radius_scaled;
  float prob_range3 = 2. * prob_range2 / radius_scaled;
  float prob_range4 = 1. * prob_range3 / radius_scaled;
  float all_ranges = prob_range1 + prob_range2 + prob_range3 + prob_range4;
  float t;
  /// \li Decide which branch \f$\tilde p^{({\rm I - IV})}\f$ to go into
  do {
    float which_range = drand48() * all_ranges - prob_range1;
    if (which_range < 0.0) {
      t = radius_scaled * (pow(drand48(), 1./3.) - 1.);
    } else {
      t = -log(drand48());
      if (which_range >= prob_range2) {
        t += -log(drand48());
        if (which_range >= prob_range2 + prob_range3) {
          t += -log(drand48());
        }
      }
    }
    /** \li Generate \f$t\f$ from the distribution in the respective
     * branches
     * \li \a reject that number with a probability
     * \f$1-(1+\exp(-|t|))^{-1}\f$ (the efficiency of this should be
     * \f$\gg \frac{1}{2}\f$)
     **/
  } while (drand48() > 1./(1. + exp(-fabs(t)) ) );
  /// \li shift and rescale \f$t\f$ to \f$r = d\cdot t + r_0\f$
  float position_scaled = t + radius_scaled;
  float position = position_scaled * diffusiveness_;
  return position;
  /// \li (choose direction; this is done outside of this routine,
  /// though).
}

void Nucleus::arrange_nucleons() {
  for (auto i = begin(); i != end(); i++) {
    // get radial position of current nucleon:
    float r = distribution_nucleons();
    // get solid angle for current nucleon:
    Angles dir;
    dir.distribute_isotropically();
    double z = r*dir.z();
    double x = r*dir.x();
    // set position of current nucleon:
    i->set_position(0.0, x, r*dir.y(), z);
    // update maximal and minimal z values
    z_max_ = (z > z_max_) ? z : z_max_;
    z_min_ = (z < z_min_) ? z : z_min_;
    // update maximal and minimal x values
    x_max_ = (x > x_max_) ? x : x_max_;
    x_min_ = (x < x_min_) ? x : x_min_;
  }
}

void Nucleus::boost(const double& beta_squared_with_sign) {
  // we use the sign of the squared velocity to get information about
  // the sign of the velocity itself.
  double sign = beta_squared_with_sign >= 0 ? 1 : -1;
  double beta_squared = std::abs(beta_squared_with_sign);
  double one_over_gamma = sqrt(1.0 - beta_squared);
  /*double gamma = 1.0/one_over_gamma;
    double gammabeta = sign*sqrt(beta_squared)*gamma;
   */
  //despite the type name and the interface, the lorentz-boost does NOT
  // take a FourVector in the physical sense, i.e., one that has
  // u_mu*u^mu=1.
  // We are talking about a /passive/ lorentz transformation here, as
  // far as I can see, so we need to boost in the direction opposite to
  // where we want to go
  //     ( The vector we transform - p - stays unchanged, but we go into
  //       a system that moves with -beta. Now in this frame, it seems
  //       like p has been accelerated with +beta.
  //     )
  FourVector u_mu(1, 0.0, 0.0, -sign*sqrt(beta_squared));
  for (auto i = begin(); i != end(); i++) {
    // a real Lorentz Transformation would leave the particles at
    // different times here, which we would then have to propagate back
    // to equal times. Since we know the result, we can simply multiply
    // the z-value with 1/gamma.
    FourVector this_position = i->position();
    this_position.set_x3(this_position.x3() * one_over_gamma);
    i->set_position(this_position);
    // for momenta, though, we CAN do normal Lorentz Boosts, since we
    // *do* want to transform the zero-component (i.e., the energy).
    FourVector this_momentum = i->momentum();
    this_momentum = this_momentum.LorentzBoost(u_mu);
    i->set_momentum(this_momentum);
  }
  // we also need to update z_max_ and z_min:
  z_max_ *= one_over_gamma;
  z_min_ *= one_over_gamma;
  return;
}

void Nucleus::fill_from_list(const std::map<int, int>& particle_list) {
  for (auto n = particle_list.cbegin(); n != particle_list.cend(); n++) {
    for (int i = 0; i < n->second; i++) {
      // append particle to list
      particles.push_back({});
      // set this particle's PDG code.
      particles.back().set_pdgcode(n->first);
    }
  }
}

void Nucleus::set_diffusiveness(const float& diffuse) {
  diffusiveness_ = diffuse;
}

void Nucleus::auto_set_masses(const Particles *external_particles) {
  for (auto p = begin(); p != end(); p++) {
    p->set_momentum(
      external_particles->particle_type(p->pdgcode()).mass(), 0.0, 0.0, 0.0);
  }
}

void Nucleus::shift(const bool is_projectile,
                    const double& initial_z_displacement,
                    const double& x_offset,
                    const double& simulation_time) {
  // amount to shift z value. If is_projectile, we shift to -z_max_,
  // else we shift to -z_min_ (z_min_ itself should be negative).
  double z_offset = is_projectile ? -z_max_ : -z_min_;
  // now, the nuclei would touch. We want them to be a little apart, so
  // we need a bigger offset.
  z_offset += initial_z_displacement;
  for (auto i = begin(); i != end(); i++) {
    FourVector this_position = i->position();
    this_position.set_x3(this_position.x3() + z_offset);
    this_position.set_x1(this_position.x1() + x_offset);
    this_position.set_x0(simulation_time);
    i->set_position(this_position);
  }
}

void Nucleus::copy_particles(Particles* external_particles) {
  for (auto p = begin(); p != end(); p++) {
    external_particles->add_data(*p);
  }
}
