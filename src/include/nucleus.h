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
#include "threevector.h"

#include <map>
#include <stdexcept>
#include <vector>

namespace Smash {

/** Nucleus: a nucleus is a collection of Particles (ParticleData thingys) that
 * are initialized before the beginning of the simulation and all have
 * the same velocity (and spatial proximity).
 *
 * Options added by Nucleus go in the "Modi"→"Nucleus"->"a nucleus" section of the
 * configuration, where "a nucleus" is either projectile or target.
 *
 * The following nucleus directives are understood:
 * -------------
 */
// !!USER:Input
/**
 * \if user
 * \page input_modi_nucleus_ Input Section Modi:Nucleus
 * \endif
 *
 * \li `AUTOMATIC:` Sets all necessary parameters based on the atomic number
 * of the input nucleus (true=automatic, false=manual, see additional directives).
 * \li `DIFFUSIVENESS` The woods-saxon parameter controlling the width of the
 * nucleus. If unspecified, the default value is 0.545.
 * \li `RADIUS` The radius of the nucleus. If not specified, the default is
 * proton_radius * A ^ (1/3).
 **/
 // !!/USER:Input
class Nucleus {
 public:
  Nucleus();

  /// returns the mass of the nucleus
  float mass() const;

  /** Returns a Woods-Saxon distributed position.
   * The distribution of return values from this function is according to a
   * spherically symmetric Woods-Saxon distribution suitable for this nucleus.
   * \f$\frac{dN}{dr} = \frac{r^2}{\exp\left(\frac{r-R}{d}\right) +
   * 1}\f$ where \f$d\f$ is the diffusiveness_ parameter and \f$R\f$ is
   * nuclear_radius_.
   **/
  virtual ThreeVector distribute_nucleon() const;

  /** returns the Woods-Saxon distribution directly
   *
   * @param x the position at which to evaluate the function
   * @return un-normalized woods-saxon probability for @param x
   **/
  float woods_saxon(float x);

  /// sets the positions of the nuclei inside nucleus A.
  void arrange_nucleons();

 /** Sets the deformation parameters of the Woods-Saxon distribution
  * according to the current mass number.
  *
  * Ref. for nuclear radii is Atom. Data Nucl. Data Tabl. by H. De Vries et. al.
  * For diffusiveness and saturation density, see [insert reference].
  */
   virtual void set_parameters_automatic();

  /** Sets the parameters of the Woods-Saxon according to
   * manually added values in the configuration file.
   *
   * @param nucleus_type A string determining the nucleus type
   *                     (either "Projectile" or "Target")
   * @param config The configuration file located at node Nucleus
   **/
  virtual void set_parameters_from_config(const char *nucleus_type,
                                          Configuration &config);

  /**
   * Boosts the nuclei so that the nucleons have the appropriate
   * momentum and the nuclei are lorentz-contracted.
   *
   * @param beta_scalar magnitude (with sign) of the z
   * component velocity used for boosting.
   **/
  void boost(double beta_scalar);

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
   * @param particle_list The particles that are added.
   * @param testparticles Number of test particles to use.
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
  void shift(bool is_projectile, double initial_z_displacement,
                     double x_offset, float simulation_time);

  /** Rotates the nucleus. (Spherical symmetry of nondeformed nuclei
   * means there is nothing to do.)
   **/
  virtual void rotate() {};

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
   * @return \f$\vec r_s = \frac{1}{N} \sum_{i=1}^N \vec r_i\f$ (for a
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
      p->set_4position(p->position()-centerpoint);
    }
  }

  /// Write the nucleon positions to a text file.
  void print_nucleus(const char * file_name) const;

  /// \ingroup exception
  struct TestparticleConfusion : public std::length_error {
    using std::length_error::length_error;
  };

 private:
  /** diffusiveness of Woods-Saxon-distribution in this nucleus in fm
   * (for diffusiveness_ == 0, we obtain a hard sphere.
   **/
  float diffusiveness_ = .545f;
  /// Saturation density of this nucleus.
  float saturation_density_ = .168f;
  /// Nuclear radius of this nucleus
  float nuclear_radius_;
  /** single-proton-radius
   *
   * \see default_nuclear_radius
   * */
  float proton_radius_ = 1.2f;
  /// Coordinate of the outermost particle.
  float r_max_ = 0.f;
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
  /// \see diffusiveness_
  inline void set_diffusiveness(float diffuse) {
    diffusiveness_ = diffuse;
  }
  /// gets the diffusiveness of the nucleus
  ///
  /// \see diffusiveness_
  inline float get_diffusiveness() const {
    return diffusiveness_;
  }
  /// sets the saturation density of the nucleus
  ///
  /// \see saturation_density_
  inline void set_saturation_density(float density) {
    saturation_density_ = density;
  }
  /// gets the saturation density of the nucleus
  ///
  /// \see saturation_density_
  inline float get_saturation_density() const {
    return saturation_density_;
  }
  /** returns a default radius for the nucleus
   *
   * Nuclear radius is calculated with the proton radius times the third
   * root of the number of nucleons.
   **/
  inline float default_nuclear_radius() {
    return proton_radius_*pow(number_of_particles(), 1./3.);
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

  /**\ingroup logging
   * Writes the state of the Nucleus object to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &, const Nucleus &);
};

}  // namespace Smash

#endif  // SRC_INCLUDE_NUCLEUS_H_
