/*
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_NUCLEUS_H_
#define SRC_INCLUDE_NUCLEUS_H_

#include<map>
#include<vector>
#include "include/particledata.h"
#include "include/particles.h"

/// A Nucleus is a collection of Particles (ParticleData thingys) that
/// are initialized before the beginning of the simulation and all have
/// the same velocity (and spatial proximity).
/// This class inherits from Particles, which is the collection of all
/// particles in the simulation and contains special functions for the
/// initialization of nuclei.
class Nucleus {
 public:
  Nucleus();

  /// returns the mass of the nucleus
  float mass() const;
  /** returns the radius of the nucleus
   *
   * Nuclear radius is calculated with the proton radius times the third
   * root of the number of nucleons.
   **/
  inline float nuclear_radius() const {
    return proton_radius_*pow(size(), 1./3.);
  }

  /** returns a Woods-Saxon distributed length 
   *
   * the distribution of return values from this function is according
   * to a Woods-Saxon distribution suitable for this nucleus.
   * \f$\frac{dN}{dr} = \frac{r^2}{\exp\left(\frac{r-R}{d}\right) +
   * 1}\f$ where \f$d\f$ is the softness_ parameter and \f$R\f$ is
   * nuclear_radius(). */
  float distribution_nucleons() const;
  /// sets the positions of the nuclei inside nucleus A.
  void arrange_nucleons();
  /**
   * Boosts the nuclei so that the nucleons have the appropriate
   * momentum and the nuclei are lorentz-contracted.
   *
   * @param beta_squared_with_sign velocity used for boosting,
   * interpreted as z-value. Note that the sign of this variable is used
   * to determine the sign of the velocity, i.e., \f$\beta_z =
   * \mathop{sign}(\beta^2)\cdot\sqrt{|\beta^2|}\f$. 
   **/
  void boost(const double& beta_squared_with_sign);
  /** Adds a particle to the nucleus
   *
   * @param pdgcode PDG code of the particle. */
  void add_particle(const int pdgcode);
  /**
   * Adds particles from a map PDG_ID => Number_of_particles_with_that_PDG_ID
   * to the nucleus.
   *
   * If the map is, e.g., [2212: 6, 2112: 7] initializes C-13 (6 protons
   * and 7 neutrons). The particles are only created, no position or
   * momenta are yet assigned. */
  void fill_from_list(const std::map<int, int>& particle_list);
  /// sets the softness of the nucleus
  ///
  /// \see softness_.
  void set_softness(const float& soft);
  /**
   * sets the masses of all nucleons automatically from the PDG info in
   * particles.
   *
   * @param particles is an object of Particles which has all the PDG
   * types read in.
   **/
  void auto_set_masses(const Particles *particles);
  /**
   * shifts the nucleus to correct impact parameter and z displacement.
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
  void shift(const bool is_projectile,
             const double& initial_z_displacement,
             const double& x_offset,
             const double& simulation_time);
  /// copies the particles from this nucleus into the particle list.
  void copy_particles(Particles* particles);

 private:
  /** softness of Woods-Saxon-distribution in this nucleus im fm
   * (for softness_ == 0, we obtain a hard sphere. */
  float softness_ = .545f;
  /** single-proton-radius 
   *
   * \see nuclear_radius
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
  /// particles associated with this nucleus.
  std::vector<ParticleData> particles;
  /// for iterators over the particle list:
  inline std::vector<ParticleData>::iterator begin() {
    return particles.begin();
  }
  /// for iterators over the particle list:
  inline std::vector<ParticleData>::iterator end() {
    return particles.end();
  }
  /// for iterators over the particle list:
  inline std::vector<ParticleData>::const_iterator cbegin() const {
    return particles.cbegin();
  }
  /// for iterators over the particle list:
  inline std::vector<ParticleData>::const_iterator cend() const {
    return particles.cend();
  }
  /// Number of particles in the list:
  inline size_t size() const {
    return particles.size();
  }
};

#endif  // SRC_INCLUDE_NUCLEUS_H_
