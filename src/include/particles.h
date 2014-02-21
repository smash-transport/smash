/*
 *    Copyright (c) 2012-2013
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARTICLES_H_
#define SRC_INCLUDE_PARTICLES_H_

#include <map>
#include <utility>

#include "include/decaymodes.h"
#include "include/particledata.h"
#include "include/particletype.h"

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

/**
 * Particles contains both the available particle list and the types
 *
 * The particle data contains the current particle list, which
 * has the changing attributes of a particle.
 * The particle type lists the static property of a certain particle.
 */
class Particles {
 public:
  /// Use improbable values for default constructor
  Particles() :id_max_(-1) {}
  /// Return the specific data of a particle according to its id
  inline const ParticleData &data(int id);
  /// Return the specific datapointer of a particle according to its id
  inline ParticleData * data_pointer(int id);
  /// Return the type of a specific particle given its id
  inline ParticleType type(const int id) const;
  /// Return the type for a specific pdgcode
  inline ParticleType particle_type(const int pdgcode) const;
  /// Return decay modes of this particle type
  inline DecayModes decay_modes(int pdg);
  /// return the highest used id
  inline int id_max(void);
  /// inserts a new particle and returns its id
  inline int add_data(const ParticleData &particle_data);
  /// inserts a new particle type
  inline void add_type(const ParticleType &particle_type, int pdg);
  /// adds decay modes for a particle type
  inline void add_decaymodes(const DecayModes &new_decay_modes, int pdg);
  /// add a range of particles
  inline void create(size_t number, int pdg);
  /* add one particle and return pointer to it */
  inline ParticleData& create(const int pdg);
  /// remove a specific particle
  inline void remove(int id);
  /// size() of the ParticleData map
  inline size_t size(void) const;
  /// empty() check of the ParticleData map
  inline bool empty(void) const;
  /// check the existence of an element in the ParticleData map
  inline size_t count(int id) const;
  /// size() check of the ParticleType map
  inline size_t types_size(void) const;
  /// empty() check of the ParticleType map
  inline bool types_empty(void) const;
  /// return time of the computational frame
  inline double time(void) const;
  /* iterators */
  inline std::map<int, ParticleData>::iterator begin(void);
  inline std::map<int, ParticleData>::iterator end(void);
  inline std::map<int, ParticleData>::const_iterator cbegin(void) const;
  inline std::map<int, ParticleData>::const_iterator cend(void) const;
  inline std::map<int, ParticleType>::const_iterator types_cbegin(void) const;
  inline std::map<int, ParticleType>::const_iterator types_cend(void) const;

 private:
  /// Highest id of a given particle
  int id_max_;
  /**
   * dynamic data of the particles a map between its id and data
   *
   * A map structure is used as particles decay and hence this
   * a swiss cheese over the runtime of SMASH. Also we want direct
   * lookup of the corresponding particle id with its data.
   */
  std::map<int, ParticleData> data_;
  /**
   * a map between pdg and correspoding static data of the particles
   *
   * PDG ids are scattered in a large range of values, hence it is a map.
   */
  std::map<int, ParticleType> types_;
  /// a map between pdg and corresponding decay modes
  std::map<int, DecayModes> all_decay_modes_;
  /// google style recommendation
  DISALLOW_COPY_AND_ASSIGN(Particles);
};

/* return the data of a specific particle */
inline const ParticleData &Particles::data(int particle_id) {
  return data_.at(particle_id);
}

/* return the pointer to the data of a specific particle */
inline ParticleData* Particles::data_pointer(int particle_id) {
  return &data_.at(particle_id);
}

/* return the type of a specific particle */
inline ParticleType Particles::type(const int particle_id) const {
  return types_.at(data_.at(particle_id).pdgcode());
}

/* return a specific type */
inline ParticleType Particles::particle_type(const int pdgcode) const {
  return types_.at(pdgcode);
}

/* return the decay modes of specific type */
inline DecayModes Particles::decay_modes(int pdg) {
  return all_decay_modes_.at(pdg);
}

/* add a new particle data and return the id of the new particle */
inline int Particles::add_data(ParticleData const &particle_data) {
  id_max_++;
  data_.insert(std::pair<int, ParticleData>(id_max_, particle_data));
  data_.at(id_max_).set_id(id_max_);
  return id_max_;
}

/* create a bunch of particles */
inline void Particles::create(size_t number, int pdgcode) {
  ParticleData particle;
  /* fixed pdgcode and no collision yet */
  particle.set_pdgcode(pdgcode);
  particle.set_collision(-1, 0, -1);
  for (size_t i = 0; i < number; i++) {
    id_max_++;
    particle.set_id(id_max_);
    data_.insert(std::pair<int, ParticleData>(id_max_, particle));
  }
}

/* create a bunch of particles */
inline ParticleData& Particles::create(int pdgcode) {
  ParticleData particle;
  /* fixed pdgcode and no collision yet */
  particle.set_pdgcode(pdgcode);
  particle.set_collision(-1, 0, -1);
  id_max_++;
  particle.set_id(id_max_);
  data_.insert(std::pair<int, ParticleData>(id_max_, particle));
  return data_[id_max_];
}

/* return the highest used id */
inline int Particles::id_max() {
  return id_max_;
}

/* add a new particle type */
inline void Particles::add_type(ParticleType const &TYPE, int pdg) {
  types_.insert(std::pair<int, ParticleType>(pdg, TYPE));
}

/* add decay modes for a particle type */
inline void Particles::add_decaymodes(const DecayModes &new_decay_modes,
                                      int pdg) {
  all_decay_modes_.insert(std::pair<int, DecayModes>(pdg, new_decay_modes));
}

/* remove a particle */
inline void Particles::remove(int id) {
  data_.erase(id);
}

/* total number of particles */
inline size_t Particles::size() const {
  return data_.size();
}

/* check if we have particles */
inline bool Particles::empty() const {
  return data_.empty();
}

/* total number particle types */
inline size_t Particles::types_size() const {
  return types_.size();
}

/* check if we have particle types */
inline bool Particles::types_empty() const {
  return types_.empty();
}

inline std::map<int, ParticleData>::iterator Particles::begin() {
  return data_.begin();
}

inline std::map<int, ParticleData>::iterator Particles::end() {
  return data_.end();
}

inline std::map<int, ParticleData>::const_iterator Particles::cbegin() const {
  return data_.begin();
}

inline std::map<int, ParticleData>::const_iterator Particles::cend() const {
  return data_.end();
}

/* we only provide const ParticleType iterators as this shouldn't change */
inline std::map<int, ParticleType>::const_iterator Particles::types_cbegin()
  const {
  return types_.begin();
}

inline std::map<int, ParticleType>::const_iterator Particles::types_cend()
  const {
  return types_.end();
}

/* check the existence of an element in the ParticleData map */
inline size_t Particles::count(int id) const {
  return data_.count(id);
}

/* return computation time which is reduced by the start up time */
inline double Particles::time() const {
  return data_.begin()->second.position().x0() - 1.0;
}

/* boost_CM - boost to center of momentum */
void boost_CM(ParticleData *particle1, ParticleData *particle2,
  FourVector *velocity);

/* boost_from_CM - boost back from center of momentum */
void boost_back_CM(ParticleData *particle1, ParticleData *particle2,
  FourVector *velocity_orig);

/* particle_distance - measure distance between two particles */
double particle_distance(ParticleData *particle_orig1,
  ParticleData *particle_orig2);

/* time_collision - measure collision time of two particles */
double collision_time(const ParticleData &particle1,
  const ParticleData &particle2);

/* momenta_exchange - soft scattering */
void momenta_exchange(ParticleData *particle1, ParticleData *particle2);

/* Sample final state momenta in general 2->2 process */
void sample_cms_momenta(ParticleData *particle1, ParticleData *particle2,
  const double cms_energy, const double mass1, const double mass2);

#endif  // SRC_INCLUDE_PARTICLES_H_
