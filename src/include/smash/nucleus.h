/*
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_NUCLEUS_H_
#define SRC_INCLUDE_NUCLEUS_H_

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
  /// \todo unused or default constructor?
  Nucleus(const std::map<PdgCode, int> &particle_list, int nTest);
  /**
   * Constructor for Nucleus, that needs the configuration parameters from
   * the inputfile and the number of testparticles
   *
   * \param[in] config contains the parameters from the inputfile on the
   * numbers of particles with a certain PDG code
   * \param[in] nTest number of testparticles
   */
  Nucleus(Configuration &config, int nTest);

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
  virtual ThreeVector distribute_nucleon() const;

  /**
   * Woods-Saxon distribution
   * \param[in] x the position at which to evaluate the function
   * \return un-normalized Woods-saxon probability
   */
  double woods_saxon(double x);

  /// Sets the positions of the nucleons inside a nucleus.
  void arrange_nucleons();

  /**
   * Sets the deformation parameters of the Woods-Saxon distribution
   * according to the current mass number.
   * Ref. for nuclear radii is \iref{DeJager:1987qc}.
   * \todo For diffusiveness and saturation density, see [insert reference].
   * \todo Issue #4743 covers the update of this part with references; Also
   * the Hirano-Nara correction should be an option
   */
  virtual void set_parameters_automatic();

  /**
   * Sets the parameters of the Woods-Saxon according to
   * manually added values in the configuration file.
   *
   * \param config The configuration for this nucleus (projectile or target).
   */
  virtual void set_parameters_from_config(Configuration &config);

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
   * Adds a particle to the nucleus
   *
   * \param pdgcode PDG code of the particle.
   * \todo unused?
   */
  void add_particle(int pdgcode);

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
   * Rotates the nucleus. (Due to spherical symmetry of nondeformed nuclei,
   * there is nothing to do.)
   */
  virtual void rotate() {}

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
    /* if size() is not a multiple of testparticles_, this will throw an
     * error. */
    if (nop * testparticles_ != particles_.size()) {
      throw TestparticleConfusion(
          "Number of test particles and test particles"
          "per particle are incompatible");
    }
    return nop;
  }

  /**
   * Calculate geometrical center of the nucleus
   * \return \f$\vec r_s = \frac{1}{N} \sum_{i=1}^N \vec r_i\f$ (for a
   * nucleus with N particles that are at the positions \f$\vec r_i\f$).
   */
  FourVector center() const;

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

  /// \ingroup exception
  struct TestparticleConfusion : public std::length_error {
    using std::length_error::length_error;
  };

 private:
  /**
   * Diffusiveness of Woods-Saxon distribution of this nucleus in fm
   * (for diffusiveness_ == 0, we obtain a hard sphere).
   */
  double diffusiveness_ = .545;
  /// Saturation density of this nucleus.
  double saturation_density_ = nuclear_density;
  /// Nuclear radius of this nucleus
  double nuclear_radius_;
  /**
   * Single proton radius in fm
   * \see default_nuclear_radius
   */
  double proton_radius_ = 1.2;
  /// Number of testparticles per physical particle
  size_t testparticles_ = 1;
  /// Particles associated with this nucleus.
  std::vector<ParticleData> particles_;

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
   * Sets the saturation density of the nucleus
   * \see saturation_density_
   */
  inline void set_saturation_density(double density) {
    saturation_density_ = density;
  }
  /**
   * \return the saturation density of the nucleus
   * \see saturation_density_
   */
  inline double get_saturation_density() const { return saturation_density_; }
  /**
   * \return a default radius for the nucleus
   * Nuclear radius is calculated with the proton radius times the third
   * root of the number of nucleons.
   */
  inline double default_nuclear_radius() {
    return proton_radius_ * std::pow(number_of_particles(), 1. / 3.);
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
   * \ingroup logging
   * Writes the state of the Nucleus object to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &, const Nucleus &);
};

}  // namespace smash

#endif  // SRC_INCLUDE_NUCLEUS_H_
