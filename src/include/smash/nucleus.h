/*
 *    Copyright (c) 2014-2022,2024-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_NUCLEUS_H_
#define SRC_INCLUDE_SMASH_NUCLEUS_H_

#include <map>
#include <stdexcept>
#include <vector>

#include "configuration.h"
#include "constants.h"
#include "forwarddeclarations.h"
#include "fourvector.h"
#include "particledata.h"
#include "threevector.h"

namespace smash {

/**
 * A nucleus is a collection of particles that are initialized,
 * before the beginning of the simulation and all have the same velocity.
 */
class Nucleus {
 public:
  /// default constructor
  Nucleus() = default;

  /**
   * Constructor for Nucleus, that needs the configuration parameters from
   * the inputfile and the number of testparticles
   *
   * \param[in] config contains the parameters from the inputfile on the
   * numbers of particles with a certain PDG code
   * \param[in] nTest number of testparticles
   */
  Nucleus(Configuration &config, int nTest);

  /**
   * Constructor which directly initializes the Nucleus with particles
   * and respective counts.
   * Only used for testing.
   *
   * \param[in] particle_list std::map, which maps PdgCode and count
   * of this particle.
   * \param[in] nTest Number of test particles.
   * \param[in] spin_interaction_type whether to use spin interactions.
   */
  Nucleus(const std::map<PdgCode, int> &particle_list, int nTest,
          SpinInteractionType spin_interaction_type = SpinInteractionType::Off);

  virtual ~Nucleus() = default;

  /**
   * \return Mass of the nucleus [GeV].
   * It needs to be double to allow for calculations at LHC energies.
   */
  double mass() const;

  /**
   * The distribution of return values from this function is according to a
   * spherically symmetric Woods-Saxon distribution suitable for this nucleus.
   * \f$\frac{dN}{dr} = \frac{r^2}{\exp\left(\frac{r-r_0}{d}\right) +
   * 1}\f$ where \f$d\f$ is the diffusiveness_ parameter and \f$r_0\f$ is
   * nuclear_radius_.
   *
   * \return  Woods-Saxon distributed position.
   */
  virtual ThreeVector distribute_nucleon();

  /**
   * Woods-Saxon distribution
   * \param[in] x the position at which to evaluate the function
   * \return un-normalized Woods-saxon probability
   */
  double woods_saxon(double x);

  /// Sets the positions of the nucleons inside a nucleus.
  virtual void arrange_nucleons();

  /**
   * Sets the deformation parameters of the Woods-Saxon distribution
   * according to the current mass number.
   * The values are taken from \iref{DeVries:1987atn} and
   * \iref{Loizides:2014vua}. They are in agreement with MC-Glauber models such
   * as GLISSANDO (see \iref{Rybczynski:2013yba}) and TGlauber MC (see
   * \iref{Loizides:2017ack}).
   */
  virtual void set_parameters_automatic();

  /**
   * Generates momenta according to Fermi motion for the nucleons.
   * For neutrons and protons Fermi momenta are calculated as
   * \f$ p_{F} = (3 \pi^2 \rho)^{1/3}\f$, where \f$ rho \f$ is
   * neutron density for neutrons and proton density for protons.
   * The actual momenta \f$p_x\f$, \f$p_y\f$, \f$p_z\f$ are
   * uniformly distributed in the sphere with radius \f$p_F\f$.
   */
  virtual void generate_fermi_momenta();

  /**
   * Boosts the nuclei into the computational frame, such that
   * the nucleons have the appropriate momentum and the
   * nuclei are lorentz-contracted. Note that the usual boost cannot be
   * applied for nuclei, since the particles would end up with different
   * times and the binding energy needs to be taken into account.
   *
   * \param[in] beta_scalar velocity in z-direction used for boost.
   */
  void boost(double beta_scalar);

  /**
   * Adds particles from a map PDG code =>
   * Number_of_particles_with_that_PDG_code to the nucleus. E.g., the map [2212:
   * 6, 2112: 7] initializes C-13 (6 protons and 7 neutrons). The particles are
   * only created, no position or momenta are yet assigned. It is also possible
   * to use any other PDG code, in addition to nucleons.
   *
   * \param[out] particle_list The particle slots that are created.
   * \param[in] testparticles Number of test particles to use.
   */
  void fill_from_list(const std::map<PdgCode, int> &particle_list,
                      int testparticles);

  /**
   * Shifts the nucleus to correct impact parameter and z displacement.
   *
   * \param[in] z_offset is the shift in z-direction
   * \param[in] x_offset is the shift in x-direction
   * \param[in] simulation_time set the time and formation_time of each
   * particle to this value.
   */
  void shift(double z_offset, double x_offset, double simulation_time);

  /**
   * Rotates the nucleus using the three euler angles phi, theta and psi.
   */
  virtual void rotate();

  /**
   * Copies the particles from this nucleus into the particle list.
   *
   * \param[out] particles Particle list with all constituents of a nucleus
   */
  void copy_particles(Particles *particles);

  /// Number of numerical (=test-)particles in the nucleus:
  inline size_t size() const { return particles_.size(); }

  /**
   * Number of physical particles in the nucleus:
   *
   * \throw TestparticleConfusion if the number of the nucleons is not a
   *        multiple of testparticles_.
   */
  inline size_t number_of_particles() const {
    size_t nop = particles_.size() / testparticles_;
    /* If size() is not a multiple of testparticles_, this will throw an
     * error. */
    if (nop * testparticles_ != particles_.size()) {
      throw TestparticleConfusion(
          "Number of test particles and test particles"
          "per particle are incompatible.");
    }
    return nop;
  }

  /**
   * Number of physical protons in the nucleus:
   *
   * \return number of protons
   * \throw Testparticleconfusion if the number of the protons is not a
   *        multiple of testparticles_.
   */
  inline size_t number_of_protons() const {
    size_t proton_counter = 0;
    /* If n_protons is not a multiple of testparticles_, this will throw an
     * error. */
    for (auto &particle : particles_) {
      if (particle.type().pdgcode() == pdg::p) {
        proton_counter++;
      }
    }

    size_t n_protons = proton_counter / testparticles_;

    if (n_protons * testparticles_ != proton_counter) {
      throw TestparticleConfusion(
          "Number of test protons and test particles"
          "per proton are incompatible.");
    }

    return n_protons;
  }

  /**
   * Calculate geometrical center of the nucleus
   * \return \f$\mathbf{r}_s = \frac{1}{N} \sum_{i=1}^N \mathbf{r}_i\f$ (for a
   * nucleus with N particles that are at the positions \f$\mathbf{r}_i\f$).
   */
  FourVector center() const;

  /// Sets target / projectile labels on nucleons
  void set_label(BelongsTo label) {
    for (ParticleData &data : particles_) {
      data.set_belongs_to(label);
    }
  }

  /**
   * Shifts the nucleus so that its center is at (0,0,0)
   * \see center()
   */
  void align_center() {
    FourVector centerpoint = center();
    for (auto p = particles_.begin(); p != particles_.end(); ++p) {
      p->set_4position(p->position() - centerpoint);
    }
  }

  /**
   * Return the Woods-Saxon probability density for the given position. This
   * corresponds to the nuclear density at the very same position.
   *
   * \param[in] r The radius at which to sample
   * \return The Woods-Saxon density
   */
  // This function as well as nucleon_density_unnormalized could in principle
  // be defined without the second argument
  virtual double nucleon_density(double r, double, double) const;
  /**
   * Return the unnormalized Woods-Saxon distribution for the given position
   * without deformation.
   *
   * \param[in] r The radius
   * \return The unnormalized Woods-Saxon distribution
   */
  virtual double nucleon_density_unnormalized(double r, double, double) const;
  /**
   * \return the normalized ground state density for the corresponding
   * Woods-Saxon parameter. This is done by integrating the Woods-Saxon
   * distribution and setting the normalization such that the integral of the
   * Woods-Saxon distribution yields the number of particles in the nucleus
   * \f$\int\rho(r)d^3r = N_{particles}\f$.
   *
   */
  virtual double calculate_saturation_density() const;
  /**
   * Sets the saturation density of the nucleus
   * \see saturation_density_
   */
  virtual void set_saturation_density(double density) {
    saturation_density_ = density;
  }

  /// \ingroup exception
  struct TestparticleConfusion : public std::length_error {
    using std::length_error::length_error;
  };

 private:
  /**
   * Diffusiveness of Woods-Saxon distribution of this nucleus in fm
   * (for diffusiveness_ == 0, we obtain a hard sphere).
   */
  double diffusiveness_;
  /// Nuclear radius of this nucleus
  double nuclear_radius_;
  /**
   * Single proton radius in fm
   * \see default_nuclear_radius
   */
  double proton_radius_ = 1.2;
  /// Number of testparticles per physical particle
  size_t testparticles_ = 1;

  /// Set unpolarized spin vectors for all particles in the nucleus
  void make_nucleus_unpolarized() {
    for (auto &particle : particles_) {
      particle.set_unpolarized_spin_vector();
    }
  }

 protected:
  /// Particles associated with this nucleus.
  std::vector<ParticleData> particles_;

  /// Saturation density of this nucleus.
  // Needed as public member for inheritance to deformed nuclei
  double saturation_density_ = nuclear_density;

  /**
   * Randomly generate Euler angles. Necessary for rotation of deformed and
   * custom nuclei, whenever a new nucleus of this kind is initialized.
   */
  void random_euler_angles();

  /**
   * The Euler angle phi of the three Euler angles used to apply rotations to
   * the nucleus. We do not use the \c Angles class here to keep a clear
   * distinction between spherical coordinates and angles for rotations.
   */
  double euler_phi_ = 0.0;
  /// Euler angle theta
  double euler_theta_ = 0.0;
  /// Euler angle psi
  double euler_psi_ = 0.0;
  /// Whether the nucleus should be rotated randomly.
  bool random_rotation_ = false;

 public:
  /// For iterators over the particle list:
  inline std::vector<ParticleData>::iterator begin() {
    return particles_.begin();
  }
  /// For iterators over the particle list:
  inline std::vector<ParticleData>::iterator end() { return particles_.end(); }
  /// For const iterators over the particle list:
  inline std::vector<ParticleData>::const_iterator cbegin() const {
    return particles_.cbegin();
  }
  /// For const iterators over the particle list:
  inline std::vector<ParticleData>::const_iterator cend() const {
    return particles_.cend();
  }
  /**
   * Sets the diffusiveness of the nucleus
   * \see diffusiveness_
   */
  inline void set_diffusiveness(double diffuse) { diffusiveness_ = diffuse; }
  /**
   * \return the diffusiveness of the nucleus
   * \see diffusiveness_
   */
  inline double get_diffusiveness() const { return diffusiveness_; }
  /**
   * \return the saturation density of the nucleus
   * \see saturation_density_
   */
  inline double get_saturation_density() const { return saturation_density_; }
  /**
   * Default nuclear radius calculated as:
   * \li \f$ r = r_\mathrm{proton} \ A^{1/3} \qquad \qquad \qquad \ \f$ for A <=
   * 16 \li \f$ r = 1.12 \ A^{1/3} - 0.86 \ A^{-1/3} \qquad \f$ for A > 16
   *
   * \return default radius for the nucleus in fm\n
   */
  inline double default_nuclear_radius() {
    int A = number_of_particles();

    if (A <= 16) {
      // radius: rough guess for all nuclei not listed explicitly with A <= 16
      return (proton_radius_ * std::cbrt(A));
    } else {
      // radius taken from \iref{Rybczynski:2013yba}
      return (1.12 * std::pow(A, 1.0 / 3.0) - 0.86 * std::pow(A, -1.0 / 3.0));
    }
  }
  /**
   * Sets the nuclear radius
   * \see nuclear_radius
   */
  inline void set_nuclear_radius(double rad) { nuclear_radius_ = rad; }
  /**
   * \return the nuclear radius
   * \see nuclear_radius
   */
  inline double get_nuclear_radius() const { return nuclear_radius_; }
  /**
   * Set angles for rotation of the nucleus from config file.
   * \param[in] orientation_config The configuration for the rotation of this
   * nucleus (projectile or target).
   */
  void set_orientation_from_config(Configuration &orientation_config);
  /**
   * \ingroup logging
   * Writes the state of the Nucleus object to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &, const Nucleus &);
};

/**
 * Find out whether a configuration has a projectile or a target sub-section.
 *
 * \param config The configuration to be checked.
 */
bool has_projectile_or_target(const Configuration &config);

/**
 * Find out whether a configuration is about projectile or target.
 *
 * \param config The configuration to be checked.
 *
 * \throw An \c std::logic_error if there is neither a projectile nor a target
 * subsection or if both are present.
 */
bool is_about_projectile(const Configuration &config);

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_NUCLEUS_H_
