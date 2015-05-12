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
#include <limits>

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
  ThreeVector velocity() const { return momentum_.velocity(); }

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

  /**
   * Construct a particle with the given type, id, and index in Particles.
   * This constructor may only be called (directly or indirectly) from
   * Particles. This constructor should be private, but can't be in order to
   * support vector::emplace_back.
   */
  ParticleData(const ParticleType &ptype, int uid, int index)
      : id_(uid), index_(index), type_(&ptype) {}

 private:
  friend class Particles;
  ParticleData() = default;

  /**
   * Each particle has a unique identifier. This identifier is used for
   * identifying the particle in the output files. It is specifically not used
   * for searching for ParticleData objects in lists of particles, though it may
   * be used to identify two ParticleData objects as referencing the same
   * particle. This is why the comparison operators depend only on the id_
   * member.
   */
  int id_ = -1;

  /**
   * Internal index in the \ref Particles list. This number is used to find the
   * Experiment-wide original of this copy.
   *
   * The value is read and written from the Particles class.
   *
   * \see Particles::data_
   */
  unsigned index_ = std::numeric_limits<unsigned>::max();

  /// counter of the last collision/decay
  int id_process_ = -1;

  /**
   * A reference to the ParticleType object for this particle (this contains
   * all the static information).
   */
  ParticleTypePtr type_ = nullptr;

  static_assert(sizeof(ParticleTypePtr) == 2,
                "");  // this leaves us two Bytes padding to use for "free"
  static_assert(sizeof(bool) <= 2, "");  // make sure we don't exceed that space
  /**
   * If \c true, the object is an entry in Particles::data_ and does not hold
   * valid particle data. Specifically iterations over Particles must skip
   * objects with `hole_ == true`. All other ParticleData instances should set
   * this member to \c false.
   *
   * \see Particles::data_
   */
  bool hole_ = false;

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

/**\ingroup logging
 * \internal
 * Helper type to attach the request for detailed printing to the type.
 */
struct PrintParticleListDetailed {
  const ParticleList &list;
};
/**\ingroup loggingk
 * Request the ParticleList to be printed in full detail (i.e. one full
 * ParticleData printout per line).
 */
inline PrintParticleListDetailed detailed(const ParticleList &list) {
  return {list};
}

/** \ingroup logging
 * Writes a detailed overview over the particles in the \p particle_list argument
 * to the stream. This overload is selected via the function \ref detailed.
 */
std::ostream &operator<<(std::ostream &out,
                         const PrintParticleListDetailed &particle_list);

}  // namespace Smash

#endif  // SRC_INCLUDE_PARTICLEDATA_H_
