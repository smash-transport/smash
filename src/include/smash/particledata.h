/*
 *    Copyright (c) 2012-2020,2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_PARTICLEDATA_H_
#define SRC_INCLUDE_SMASH_PARTICLEDATA_H_

#include <limits>
#include <utility>

#include "forwarddeclarations.h"
#include "fourvector.h"
#include "particletype.h"
#include "pdgcode.h"
#include "processbranch.h"

namespace smash {

enum class BelongsTo : uint8_t {
  Nothing = 0,
  Projectile = 1,
  Target = 2,
};

/**
 * A structure to hold information about the history of the particle,
 * e.g. the last interaction etc.
 */
struct HistoryData {
  /// Collision counter per particle, zero only for initially present particles
  int32_t collisions_per_particle = 0;
  /// id of the last action
  int32_t id_process = 0;
  /// type of the last action
  ProcessType process_type = ProcessType::None;
  /**
   * Time of the last action (excluding walls), time of kinetic freeze_out
   * for HBT analysis this time should be larger or equal to the formation
   * time of the particle, since only formed particles can freeze out
   * The full coordinate space 4-vector can be obtained by back-propagation
   */
  double time_last_collision = smash_NaN<double>;
  /// PdgCode of the first parent particles
  PdgCode p1 = 0x0;
  /// PdgCode of the second parent particles
  PdgCode p2 = 0x0;
};

/**
 * \ingroup data
 *
 * ParticleData contains the dynamic information of a certain particle.
 *
 * Each particle has its momentum, position and other relevant physical
 * data entry.
 */
class ParticleData {
 public:
  /**
   * Create a new particle with the given \p particle_type and optionally a
   * specific \p unique_id.
   *
   * All other values are initialized to unphysical values.
   *
   * \param[in] particle_type Type of particle to be created
   * \param[in] unique_id id of particle to be created
   */
  explicit ParticleData(const ParticleType &particle_type, int unique_id = -1)
      : id_(unique_id), type_(&particle_type) {}

  /**
   * Get the id of the particle
   * \return particle id
   */
  int32_t id() const { return id_; }
  /**
   * Set id of the particle
   * \param[in] i id to be assigned to the particle
   */
  void set_id(int i) { id_ = i; }

  /**
   * Get the pdgcode of the particle
   * \return pdgcode of the particle
   */
  PdgCode pdgcode() const { return type_->pdgcode(); }

  // Convenience accessors to PdgCode:
  /// \copydoc PdgCode::is_hadron
  bool is_hadron() const { return type_->is_hadron(); }

  /// \copydoc PdgCode::is_baryon
  bool is_baryon() const { return pdgcode().is_baryon(); }

  /// \copydoc PdgCode::is_nucleus
  bool is_nucleus() const { return pdgcode().is_nucleus(); }

  /// \copydoc PdgCode::is_rho
  bool is_rho() const { return type_->is_rho(); }

  /// \copydoc PdgCode::is_proton
  bool is_proton() const { return pdgcode().is_proton(); }

  /// \copydoc PdgCode::is_neutron
  bool is_neutron() const { return pdgcode().is_neutron(); }

  /// \copydoc PdgCode::is_pion
  bool is_pion() const { return pdgcode().is_pion(); }

  /// \copydoc PdgCode::is_sigmastar
  bool is_sigmastar() const { return pdgcode().is_Sigmastar(); }

  /**
   * Get the particle's pole mass ("on-shell").
   * \return pole mass of the particle [GeV]
   */
  double pole_mass() const { return type_->mass(); }
  /**
   * Get the particle's effective mass
   *
   * Determined from the 4-momentum \f$m=\sqrt{p_\mu p^\mu}\f$.
   * Possibly "off-shell".
   * \return Effective mass [GeV]
   */
  double effective_mass() const;
  /**
   * Get the type of the particle
   * \return ParticleType object associated to this particle.
   */
  const ParticleType &type() const { return *type_; }

  /**
   * Get the id of the last action
   * \return id of particle's latest collision
   */
  uint32_t id_process() const { return history_.id_process; }
  /**
   * Get history information
   * \return particle history struct
   */
  HistoryData get_history() const { return history_; }
  /**
   * Set history_ from rvalue reference. Meant to be used only in
   * special situations e.g. in the ListModus, where a temporary HistoryData
   * is constructed from the user input.
   *
   * \param[in] history object to be moved from.
   */
  void set_history(HistoryData &&history) { history_ = std::move(history); }

  /**
   * Store history information
   *
   * The history contains the type of process and possibly the
   * PdgCodes of the parent particles (\p plist). Note that  history is not set
   * for dileptons and photons.
   * \param[in] ncoll particle's number of collisions
   * \param[in] pid id of the particle's latest process
   * \param[in] pt process type of the particle's latest process
   * \param[in] time_last_coll time of latest collision [fm]
   * \param[in] plist list of parent particles */
  void set_history(int ncoll, uint32_t pid, ProcessType pt,
                   double time_last_coll, const ParticleList &plist);

  /**
   * Get the particle's 4-momentum
   * \return particle's 4-momentum [GeV]
   */
  const FourVector &momentum() const { return momentum_; }

  /**
   * Set the particle's 4-momentum directly
   * \param[in] momentum_vector 4-vector \f$p^\mu = (E,\mathbf{p})^T\f$
   */
  void set_4momentum(const FourVector &momentum_vector) {
    momentum_ = momentum_vector;
  }

  /**
   * Set the momentum of the particle given its mass and momentum three-vector.
   *
   * \param[in] mass the mass of the particle (without E_kin contribution) [GeV]
   * \param[in] mom the three-momentum of the particle [GeV]
   */
  void set_4momentum(double mass, const ThreeVector &mom) {
    momentum_ = FourVector(std::sqrt(mass * mass + mom * mom), mom);
  }

  /**
   * Set the momentum of the particle.
   *
   * \param[in] mass the mass of the particle (without E_kin contribution) [GeV]
   * \param[in] px x-component of the momentum [GeV]
   * \param[in] py y-component of the momentum [GeV]
   * \param[in] pz z-component of the momentum [GeV]
   */
  void set_4momentum(double mass, double px, double py, double pz) {
    momentum_ = FourVector(std::sqrt(mass * mass + px * px + py * py + pz * pz),
                           px, py, pz);
  }
  /**
   * Set the momentum of the particle without modifying the energy.
   *
   * WARNING: Mass gets modified.
   * \param[in] mom momentum 3-vector [GeV]
   */
  void set_3momentum(const ThreeVector &mom) {
    momentum_ = FourVector(momentum_.x0(), mom);
  }

  /**
   * Get the particle's position in Minkowski space
   * \return particle's position 4-vector
   */
  const FourVector &position() const { return position_; }
  /**
   * Set the particle's 4-position directly
   * \param[in] pos position 4-vector
   */
  void set_4position(const FourVector &pos) { position_ = pos; }
  /**
   * Set particle's 3-position
   *
   * The time component is not changed
   * \param[in] pos position 3-vector
   */
  void set_3position(const ThreeVector &pos) {
    position_ = FourVector(position_.x0(), pos);
  }

  /**
   * Translate the particle position
   * \param[in] delta 3-vector by which the particle is translated [fm]
   */
  ParticleData translated(const ThreeVector &delta) const {
    ParticleData p = *this;
    p.position_[1] += delta[0];
    p.position_[2] += delta[1];
    p.position_[3] += delta[2];
    return p;
  }

  /**
   * Get the absolute formation time of the particle
   * \return particle's formation time
   */
  double formation_time() const { return formation_time_; }
  /**
   * Get the absolute time, where the cross section scaling factor slowly
   * starts increasing from the given scaling factor to 1
   * \return time, when scaling factor starts increasing
   */
  double begin_formation_time() const { return begin_formation_time_; }

  /**
   * Set the absolute formation time
   *
   * The particle's cross section scaling factor will be a Heavyside fuction
   * of time.
   * \param[in] form_time absolute formation time
   */
  void set_formation_time(double form_time) {
    formation_time_ = form_time;
    // cross section scaling factor will be a step function in time
    begin_formation_time_ = form_time;
    // if time of the last collision is NaN set it to the formation time
    if (std::isnan(history_.time_last_collision)) {
      history_.time_last_collision = form_time;
    }
  }
  /**
   * Set the time, when the cross section scaling factor begins, and finishes
   * to increase from the given cross section scaling factor to 1.
   *
   * The cross section will only grow slowly, if the option is used.
   *
   * \param[in] begin_form_time time when the cross section starts to increase
   * \param[in] form_time time when the cross section reaches 1
   */
  void set_slow_formation_times(double begin_form_time, double form_time) {
    begin_formation_time_ = begin_form_time;
    formation_time_ = form_time;
    if (std::isnan(history_.time_last_collision)) {
      history_.time_last_collision = form_time;
    }
  }

  /**
   * Get the initially assigned cross section scaling factor.
   *
   * Depending on the config, the cross section scaling factor might change
   * with time, while this value will not be updated.
   *
   * \return particle's initially assigned cross section scaling factor
   */
  const double &initial_xsec_scaling_factor() const {
    return initial_xsec_scaling_factor_;
  }
  /**
   * Set the particle's initial cross_section_scaling_factor
   *
   * All cross sections of this particle are scaled down by this factor until
   * the formation time is over.
   *
   * If the particle formation power is set to be positive, this will only be
   * the initial scaling factor, while the actual scaling factor grows with
   * time.
   *
   * \param[in] xsec_scal cross section scaling factor
   */
  void set_cross_section_scaling_factor(const double &xsec_scal) {
    initial_xsec_scaling_factor_ = xsec_scal;
  }

  /**
   * Get the velocity 3-vector
   * \return 3-velocity of the particle
   */
  ThreeVector velocity() const { return momentum_.velocity(); }

  /**
   * Get the inverse of the gamma factor from the current velocity of the
   * particle.
   *
   * \f[\frac{1}{\gamma}=\sqrt{1-v^2}\f]
   *
   * This functions is more efficient than calculating the gamma factor from
   * \ref velocity, since the \ref velocity function must execute three
   * divisions (for every space component of the momentum vector).
   *
   * \returns inverse gamma factor
   */
  double inverse_gamma() const {
    return std::sqrt(1. - momentum_.sqr3() / (momentum_.x0() * momentum_.x0()));
  }

  /**
   * Apply a full Lorentz boost of momentum and position
   * \param[in] v boost 3-velocity
   */
  void boost(const ThreeVector &v) {
    set_4momentum(momentum_.lorentz_boost(v));
    set_4position(position_.lorentz_boost(v));
  }

  /**
   * Apply a Lorentz-boost to only the momentum
   * \param[in] v boost 3-veloctity
   */
  void boost_momentum(const ThreeVector &v) {
    set_4momentum(momentum_.lorentz_boost(v));
  }
  /**
   * Get the (maximum positive) spin s of a particle in multiples of 1/2.
   * E.g. for a spin-1 particle s=2.
   * \return particle's spin in multiples of 1/2
   */
  int spin() const { return pdgcode().spin(); }
  /**
   * Get the mean spin 4-vector (Pauli–Lubanski vector) of the particle (const
   * reference, no copy). \return particle's mean spin 4-vector
   */
  const FourVector &spin_vector() const { return spin_vector_; }
  /**
   * Get the mean spin 4-vector (Pauli–Lubanski vector) of the particle (non
   * const reference). \return particle's mean spin 4-vector
   */
  FourVector &spin_vector() { return spin_vector_; }
  /**
   * Set the mean spin 4-vector (Pauli–Lubanski vector) of the particle.
   * \param[in] s particle's mean spin 4-vector
   */
  void set_spin_vector(const FourVector &s) { spin_vector_ = s; }
  /**
   * Set a single component of the mean spin 4-vector (Pauli-Lubanski vector).
   * \param[in] index component index (0-3)
   * \param[in] value component value
   */
  void set_spin_vector_component(int index, double value) {
    if (index < 0 || index > 3) {
      throw std::out_of_range("Invalid spin vector component index");
    }
    spin_vector_[index] = value;
  }

  /**
   * Set the 4 components of the spin vector such that the particle is
   * unpolarized. This function is used only to initialize the spin vector of
   * particles at creation.
   */
  void set_unpolarized_spin_vector();
  /// Setter for belongs_to label
  void set_belongs_to(BelongsTo label) { belongs_to_ = label; }
  /// Getter for belongs_to label
  BelongsTo belongs_to() const { return belongs_to_; }

  /// Fluidize the particle
  void fluidize() { core_ = true; }
  /// Check whether the particle is core
  bool is_core() const { return core_; }
  /// Particle \f$tau\f$ (hyperbolic time)
  double hyperbolic_time() const { return position_.tau(); }
  /// Particle spacetime rapidity \f$\eta_s\f$
  double spatial_rapidity() const { return position_.eta(); }
  /// Particle \f$m_T\f$
  double transverse_mass() const { return momentum_.tau(); }
  /// Particle momentum rapidity \f$y_\mathrm{rap}\f$
  double rapidity() const { return momentum_.eta(); }

  /**
   * Check whether two particles have the same id
   * \param[in] a particle to compare to
   * \return whether the particles have the same id
   */
  bool operator==(const ParticleData &a) const { return this->id_ == a.id_; }
  /**
   * Check if this particle has a smaller id than another particle
   * \param[in] a particle to compare to
   * \return whether this particle has a smaller id than other particle
   */
  bool operator<(const ParticleData &a) const { return this->id_ < a.id_; }

  /**
   * Check if the particle has a given id
   * \param[in] id_a id to compare to
   * \return whether the particle has the given id
   */
  bool operator==(int id_a) const { return this->id_ == id_a; }
  /**
   * Check whether the particle's id is smaller than the given id
   * \param[in] id_a number to compare particle's id to
   */
  bool operator<(int id_a) const { return this->id_ < id_a; }

  /**
   * Construct a particle with the given type, id and index in Particles.
   *
   * This constructor may only be called (directly or indirectly) from
   * Particles. This constructor should be private, but can't be in order to
   * support vector::emplace_back.
   * \param[in] ptype Type of the particle to be constructed
   * \param[in] uid id of the particle to be constructed
   * \param[in] index index of the particle to be constructed
   */
  ParticleData(const ParticleType &ptype, int uid, int index)
      : id_(uid), index_(index), type_(&ptype) {}

  /**
   * Return the cross section scaling factor at a given time.
   *
   * \param[in] delta_time scaling factor at current time plus
   *            this time will be returned.
   * \return the cross section scaling factor at a specified time.
   */
  double xsec_scaling_factor(double delta_time = 0.) const;

  /// Power with which the cross section scaling factor grows in time.
  static double formation_power_;

 private:
  friend class Particles;
  /// Default constructor.
  ParticleData() = default;

  /**
   * Copies some information of the particle to the given particle \p dst.
   *
   * Specifically it avoids to copy id_, index_, and type_.
   * \param[in] dst particle values are copied to
   */
  void copy_to(ParticleData &dst) const {
    dst.history_ = history_;
    dst.momentum_ = momentum_;
    dst.position_ = position_;
    dst.spin_vector_ = spin_vector_;
    dst.formation_time_ = formation_time_;
    dst.initial_xsec_scaling_factor_ = initial_xsec_scaling_factor_;
    dst.begin_formation_time_ = begin_formation_time_;
    dst.belongs_to_ = belongs_to_;
    dst.core_ = core_;
  }

  /**
   * Each particle has a unique identifier. This identifier is used for
   * identifying the particle in the output files. It is specifically not used
   * for searching for ParticleData objects in lists of particles, though it may
   * be used to identify two ParticleData objects as referencing the same
   * particle. This is why the comparison operators depend only on the id_
   * member.
   */
  int32_t id_ = -1;

  /**
   * Internal index in the \ref Particles list. This number is used to find the
   * Experiment-wide original of this copy.
   *
   * The value is read and written from the Particles class.
   *
   * \see Particles::data_
   */
  unsigned index_ = std::numeric_limits<unsigned>::max();

  /**
   * A reference to the ParticleType object for this particle (this contains
   * all the static information). Default-initialized with an invalid index.
   */
  ParticleTypePtr type_;

  // this leaves us two Bytes padding to use for "free"
  static_assert(sizeof(ParticleTypePtr) == 2, "");
  // make sure we don't exceed that space
  static_assert(sizeof(bool) <= 2, "");
  /**
   * If \c true, the object is an entry in Particles::data_ and does not hold
   * valid particle data. Specifically iterations over Particles must skip
   * objects with `hole_ == true`. All other ParticleData instances should set
   * this member to \c false.
   *
   * \see Particles::data_
   */
  bool hole_ = false;

  /// If the particle is part of a pseudofluid.
  // A particle cannot be un-core, and any children it produces inherits
  // this trait.
  bool core_ = false;

  /// momenta of the particle: x0, x1, x2, x3 as E, px, py, pz
  FourVector momentum_;
  /// position in space: x0, x1, x2, x3 as t, x, y, z
  FourVector position_;
  /**
   * Pauli-Lubanski vector (mean spin 4-vector) of the particle. Each
   * component is initialized with NaN (double) to indicate that the spin vector
   * has not been set.
   */
  FourVector spin_vector_ = FourVector(smash_NaN<double>, smash_NaN<double>,
                                       smash_NaN<double>, smash_NaN<double>);
  /** Formation time at which the particle is fully formed
   *  given as an absolute value in the computational frame
   */
  double formation_time_ = smash_NaN<double>;
  /// time when the cross section scaling factor starts to increase to 1
  double begin_formation_time_ = smash_NaN<double>;
  /**
   * Initial cross section scaling factor.
   * 1 by default, since a particle is fully formed in this case.
   */
  double initial_xsec_scaling_factor_ = 1.0;
  /// history information
  HistoryData history_;
  /// is it part of projectile or target nuclei?
  BelongsTo belongs_to_ = BelongsTo::Nothing;
};

/**
 * \ingroup logging
 * Writes the state of the particle to the output stream.
 */
std::ostream &operator<<(std::ostream &s, const ParticleData &p);

/**
 * \ingroup logging
 * Writes a compact overview over the particles in the \p particle_list argument
 * to the stream.
 */
std::ostream &operator<<(std::ostream &out, const ParticleList &particle_list);

/**
 * \ingroup logging
 * \internal
 * Helper type to attach the request for detailed printing to the type.
 */
struct PrintParticleListDetailed {
  /// Particle list
  const ParticleList &list;
};
/**
 * \ingroup logging
 * Request the ParticleList to be printed in full detail (i.e. one full
 * ParticleData printout per line).
 */
inline PrintParticleListDetailed detailed(const ParticleList &list) {
  return {list};
}

/**
 * \ingroup logging
 * Writes a detailed overview over the particles in the \p particle_list
 * argument
 * to the stream. This overload is selected via the function detailed.
 */
std::ostream &operator<<(std::ostream &out,
                         const PrintParticleListDetailed &particle_list);

/**
 * This function creates a SMASH particle validating the provided information.
 *
 *  - The input position and momentum is checked for nan values.
 *  - A particle is first created using the given PDG code, setting its
 *    4-momentum to the given one.
 *  - Afterwards, if stable, its mass is compared to the given one and, if these
 *    do not match (up to numeric rounding), the internal SMASH value (i.e. that
 *    from the particles file) is used to put the particle on the SMASH mass
 *    shell.
 *  - Finally, for unstable particles or if the previous mass check passed, the
 *    particle is checked to be on its mass shell and, if not, its energy is
 *    adjusted to put the particle on its mass shell.
 *
 * This function possibly warns the user, if requested.
 *
 * \param[in] pdgcode PdgCode  of the particle which is supposed to be checked
 * \param[in] mass Mass of the new particle
 * \param[in] four_position Position four vector of the new particle
 * \param[in] four_momentum Momentum four vector of the new particle
 * \param[in] log_area Logging area for the warning
 * \param[inout] mass_warning Whether to warn about mass discrepancies
 * \param[inout] on_shell_warning Whether to warn about off-shell particles
 * \return Valid SMASH particle matching all given input quantities
 *
 * \note The boolean flags are passed by reference, since we want to allow
 * client code to warn the user only once per flag. Hence, this function is
 * turning the flags to \c false after having warned the user.
 */
ParticleData create_valid_smash_particle_matching_provided_quantities(
    PdgCode pdgcode, double mass, const FourVector &four_position,
    const FourVector &four_momentum, int log_area, bool &mass_warning,
    bool &on_shell_warning);

/**
 * Utility function to compare two \c ParticleData instances with respect to
 * their PDG code, 4-position and 4-momenta.
 *
 * The particles are propagated to the given time before being compared.
 * 4-vectors are compared using the \c FourVector::operator== overload.
 *
 * @param p1 The first particle
 * @param p2 The second particle
 * @param time The time at which the comparison should take place
 * @return \c true if the two particles have the same PDG codes, 4-position and
 * 4-momentum;
 * @return \c false otherwise.
 */
bool are_particles_identical_at_given_time(const ParticleData &p1,
                                           const ParticleData &p2, double time);

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_PARTICLEDATA_H_
