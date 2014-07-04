/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_NUCLEUS_H_
#define SRC_INCLUDE_NUCLEUS_H_

#include "configuration.h"
#include "forwarddeclarations.h"
#include "fourvector.h"
#include "particledata.h"

#include <map>
#include <stdexcept>
#include <vector>

namespace Smash {

/// A Nucleus is a collection of Particles (ParticleData thingys) that
/// are initialized before the beginning of the simulation and all have
/// the same velocity (and spatial proximity).
class Nucleus {
 public:
  Nucleus();

  /// returns the mass of the nucleus
  float mass() const;

  /** Returns a Woods-Saxon distributed length.
   *
   * the distribution of return values from this function is according to a 
   * spherically symmetric Woods-Saxon distribution suitable for this nucleus.
   * \f$\frac{dN}{dr} = \frac{r^2}{\exp\left(\frac{r-R}{d}\right) +
   * 1}\f$ where \f$d\f$ is the diffusiveness_ parameter and \f$R\f$ is
   * nuclear_radius_. */
  float distribute_nucleon() const;

  /** returns the Woods-Saxon distribution directly
   *
   * @param x the position at which to evaluate the function
   * @return the
   **/
  float woods_saxon(float x);

  /// sets the positions of the nuclei inside nucleus A.
  virtual void arrange_nucleons();

  // Sets the parameters of the Woods-Saxon distribution
  // according to the current mass number.
  virtual size_t determine_nucleus();

  // Sets the parameters of the nucleus according to
  // manually added values in the configuration file.
  // Returns the type of nucleus (Projectile or Target).
  virtual void manual_nucleus(bool is_projectile, Configuration &config);

  /**
   * Boosts the nuclei so that the nucleons have the appropriate
   * momentum and the nuclei are lorentz-contracted.
   *
   * @param beta_squared_with_sign velocity used for boosting,
   * interpreted as z-value. Note that the sign of this variable is used
   * to determine the sign of the velocity, i.e., \f$\beta_z =
   * \mathop{sign}(\beta^2)\cdot\sqrt{|\beta^2|}\f$.
   **/
  void boost(double beta_squared_with_sign);

  /** Adds a particle to the nucleus
   *
   * @param pdgcode PDG code of the particle. */
  void add_particle(int pdgcode);

  /**
   * Adds particles from a map PDG_ID => Number_of_particles_with_that_PDG_ID
   * to the nucleus.
   *
   * If the map is, e.g., [2212: 6, 2112: 7] initializes C-13 (6 protons
   * and 7 neutrons). The particles are only created, no position or
   * momenta are yet assigned.
   *
   * \param particle_list The particles that are added.
   * \param testparticles Number of test particles to use.
   *
   **/
  void fill_from_list(const std::map<PdgCode, int>& particle_list,
                      int testparticles);

  /**
   * Shifts the nucleus to correct impact parameter and z displacement.
   *
   * @param is_projectile switches if the projectile is shifted to
   * -z_max_ or -z_min_ (the projetcile is shifted to -z_max_, so that
   *  the particle at highest z is at z = 0, and the target is shifted
   *  to -z_min_, so that the leftmost particle is at z = 0.
   *
   * @param initial_z_displacement is the additional shift in z
   * direction, so that two nuclei do not touch each other at the
   * beginning.
   *
   * @param x_offset is the shift in x-direction (for impact parameter
   * setting).
   *
   * @param simulation_time set the time of each particle to this value.
   **/
  virtual void shift(bool is_projectile, double initial_z_displacement,
                     double x_offset, float simulation_time);

  /// copies the particles from this nucleus into the particle list.
  void copy_particles(Particles* particles);

  /// Number of numerical (=test-)particles in the nucleus:
  inline size_t size() const {
    return particles_.size();
  }

  /// Number of physical particles in the nucleus:
  inline size_t number_of_particles() const {
    int nop = particles_.size()/testparticles_;
    // if size() is not a multiple of testparticles_, this will throw an
    // error.
    if (nop * testparticles_ != particles_.size()) {
      throw TestparticleConfusion("Number of test particles and test particles"
            "per particle are incompatible");
    }
    return nop;
  }

  /** returns the geometrical center of the nucleus.
   *
   * \return \f$\vec r_s = \frac{1}{N} \sum_{i=1}^N \vec r_i\f$ (for a
   * nucleus with N particles that are at the positions \f$\vec r_i\f$).
   */
  FourVector center() const;
  
  /** shifts the nucleus so that its center is at (0,0,0)
   *
   * \see center()
   */
  void align_center() {
    FourVector centerpoint = center();
    for (auto p = particles_.begin(); p != particles_.end(); ++p) {
      p->set_position(p->position()-centerpoint);
    }
  }

  struct TestparticleConfusion : public std::length_error {
    using std::length_error::length_error;
  };

 private:
  /** diffusiveness of Woods-Saxon-distribution in this nucleus in fm
   * (for diffusiveness_ == 0, we obtain a hard sphere. */
  float diffusiveness_ = .545f;
  // Saturation density of nuclear matter
  float saturation_density_ = 1.f;
  // Nuclear radius of this nucleus
  float nuclear_radius_;
  /** single-proton-radius
   *
   * \see default_nuclear_radius
   * */
  float proton_radius_ = 1.2f;
  /// z (beam direction-) coordinate of the outermost particle (highest
  /// z)
  float z_max_ = 0.f;
  /// z (beam direction-) coordinate of the outermost particle (lowest
  /// z)
  float z_min_ = 0.f;
  /// x (impact parameter direction-) coordinate of the outermost
  /// particle (highest x)
  float x_max_ = 0.f;
  /// x (impact parameter direction-) coordinate of the outermost
  /// particle (lowest x)
  float x_min_ = 0.f;
  /// Number of testparticles per physical particle
  size_t testparticles_ = 1;
  /// particles associated with this nucleus.
  std::vector<ParticleData> particles_;
 public:
  /// for iterators over the particle list:
  inline std::vector<ParticleData>::iterator begin() {
    return particles_.begin();
  }
  /// for iterators over the particle list:
  inline std::vector<ParticleData>::iterator end() {
    return particles_.end();
  }
  /// for iterators over the particle list:
  inline std::vector<ParticleData>::const_iterator cbegin() const {
    return particles_.cbegin();
  }
  /// for iterators over the particle list:
  inline std::vector<ParticleData>::const_iterator cend() const {
    return particles_.cend();
  }
  /// sets the diffusiveness of the nucleus
  ///
  /// \see diffusiveness_.
  inline void set_diffusiveness(float diffuse) {
    diffusiveness_ = diffuse;
  }
  /// gets the diffusiveness of the nucleus
  ///
  /// \see diffusiveness
  inline float get_diffusiveness() const {
    return diffusiveness_;
  }
  /// sets the saturation density of the nucleus
  ///
  /// \see saturation_density
  inline void set_saturation_density(float density) {
    saturation_density_ = density;
  }
  inline float get_saturation_density() const {
    return saturation_density_;
  }
  /** sets and returns a default radius for the nucleus
   *
   * Nuclear radius is calculated with the proton radius times the third
   * root of the number of nucleons.
   **/
  inline float default_nuclear_radius() {
    nuclear_radius_ = proton_radius_*pow(number_of_particles(), 1./3.);
    return nuclear_radius_;
  }
  /// sets the nuclear radius
  ///
  /// \see nuclear_radius 
  inline void set_nuclear_radius(float rad) {
    nuclear_radius_ = rad;
  }
  /// gets the nuclear radius
  ///
  /// \see nuclear_radius 
  inline float get_nuclear_radius() const {
    return nuclear_radius_;
  }
};

}  // namespace Smash

#endif  // SRC_INCLUDE_NUCLEUS_H_
