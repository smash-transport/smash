/*
 *    Copyright (c) 2012-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARTICLEDATA_H_
#define SRC_INCLUDE_PARTICLEDATA_H_

#include "forwarddeclarations.h"
#include "fourvector.h"
#include "particletype.h"
#include "pdgcode.h"

namespace Smash {

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
   * All other values are initialized to improbable values.
   */
  explicit ParticleData(const ParticleType &particle_type, int unique_id = -1)
      : id_(unique_id), type_(&particle_type) {}

  /// look up the id of the particle
  int id() const { return id_; }
  /// set id of the particle
  void set_id(int i) { id_ = i; }

  /// look up the pdgcode of the particle
  PdgCode pdgcode() const { return type_->pdgcode(); }

  // convenience accessors to PdgCode:
  /// \copydoc PdgCode::is_hadron
  bool is_hadron() const { return type_->is_hadron(); }

  /// \copydoc PdgCode::is_baryon
  bool is_baryon() const { return pdgcode().is_baryon(); }

  /** Returns the particle's pole mass ("on-shell"). */
  float pole_mass() const { return type_->mass(); }
  /** Returns the particle's effective mass
   * (as determined from the 4-momentum, possibly "off-shell"). */
  float effective_mass() const { return momentum().abs(); }

  /**
   * Return the ParticleType object associated to this particle.
   */
  const ParticleType &type() const { return *type_; }

  /// look up the id of the collision process
  int id_process() const { return id_process_; }
  /// set the id of the collision process
  void set_id_process(int i) { id_process_ = i; }

  /// return the particle's 4-momentum
  const FourVector &momentum() const { return momentum_; }
  /// set the particle's 4-momentum directly
  void set_4momentum(const FourVector &momentum_vector) {
    momentum_ = momentum_vector;
  }
  /**
   * Set the momentum of the particle.
   *
   * \param[in] mass the mass of the particle (without E_kin contribution)
   * \param[in] mom the three-momentum of the particle
   *
   * \fpPrecision The momentum FourVector requires double-precision.
   */
  void set_4momentum(double mass, const ThreeVector &mom) {
    momentum_ = FourVector(std::sqrt(mass * mass + mom * mom), mom);
  }
  /**
   * Set the momentum of the particle.
   *
   * \param[in] mass the mass of the particle (without E_kin contribution)
   * \param[in] px x-component of the momentum
   * \param[in] py y-component of the momentum
   * \param[in] pz z-component of the momentum
   *
   * \fpPrecision The momentum FourVector requires double-precision.
   */
  void set_4momentum(double mass, double px, double py, double pz) {
    momentum_ = FourVector(std::sqrt(mass * mass + px * px + py * py + pz * pz),
                           px, py, pz);
  }
  /**
   * Set the momentum of the particle without modifying the currently set mass.
   */
  void set_3momentum(const ThreeVector &mom) {
    momentum_ = FourVector(momentum_.x0(), mom);
  }

  /// The particle's position in Minkowski space
  const FourVector &position() const { return position_; }
  /// Set the particle's position directly
  void set_4position(const FourVector &pos) { position_ = pos; }
  /// Set the particle 3-position only (the time component is not changed)
  void set_3position(const ThreeVector &pos) {
    position_ = FourVector(position_.x0(), pos);
  }

  /// get the velocity 3-vector
  ThreeVector velocity() const { return momentum_.threevec() / momentum_.x0(); }

  /**
   * Returns the inverse of the gamma factor from the current velocity of the
   * particle.
   *
   * \f[\frac{1}{\gamma}=\sqrt{1-v^2}\f]
   *
   * This functions is more efficient than calculating the gamma factor from
   * \ref velocity, since the \ref velocity function must execute three
   * divisions (for every space component of the momentum vector).
   *
   * \fpPrecision This function must use double-precision for the calculation of
   * \f$ \beta \f$ and \f$ 1-\beta \f$ as the latter results in a value close to
   * zero and thus exhibits catastrophic cancelation.
   */
  double inverse_gamma() const {
    return std::sqrt(1. -
                     (momentum_.x1() * momentum_.x1() +
                      momentum_.x2() * momentum_.x2() +
                      momentum_.x3() * momentum_.x3()) /
                         (momentum_.x0() * momentum_.x0()));
  }

  /// Apply a full Lorentz boost of momentum and position
  void boost(const ThreeVector &v) {
    set_4momentum(momentum_.LorentzBoost(v));
    set_4position(position_.LorentzBoost(v));
  }

  /// Apply a Lorentz-boost of only the momentum
  void boost_momentum(const ThreeVector &v) {
    set_4momentum(momentum_.LorentzBoost(v));
  }

  /// Returns whether the particles are identical
  bool operator==(const ParticleData &a) const { return this->id_ == a.id_; }
  /// Defines a total order of particles according to their id.
  bool operator<(const ParticleData &a) const { return this->id_ < a.id_; }

  /// check if the particles are identical to a given id
  bool operator==(int id_a) const { return this->id_ == id_a; }
  /// sort particles along to given id
  bool operator<(int id_a) const { return this->id_ < id_a; }

 private:
  /// Each particle has a unique identifier
  int id_ = -1;

  /**
   * A reference to the ParticleType object for this particle (this contains
   * all the static information).
   */
  ParticleTypePtr type_ = nullptr;

  /// counter of the last collision/decay
  int id_process_ = -1;
  /// momenta of the particle: x0, x1, x2, x3 as E, px, py, pz
  FourVector momentum_;
  /// position in space: x0, x1, x2, x3 as t, x, y, z
  FourVector position_;
};

/**\ingroup logging
 * Writes the state of the particle to the output stream.
 */
std::ostream &operator<<(std::ostream &s, const ParticleData &p);

/** \ingroup logging
 * Writes a compact overview over the particles in the \p particle_list argument
 * to the stream.
 */
std::ostream &operator<<(std::ostream &out, const ParticleList &particle_list);

}  // namespace Smash

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
