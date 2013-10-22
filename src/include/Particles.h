/*
 *    Copyright (c) 2012-2013
 *      maximilian attems <attems@fias.uni-frankfurt.de>
 *      Jussi Auvinen <auvinen@fias.uni-frankfurt.de>
 *
 *    GNU General Public License (GPLv3)
 */
#ifndef SRC_INCLUDE_PARTICLES_H_
#define SRC_INCLUDE_PARTICLES_H_

#include <map>
#include <utility>

#include "../include/DecayModes.h"
#include "../include/ParticleData.h"
#include "../include/ParticleType.h"

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

class Particles {
  public:
  /* Use improbable values for default constructor */
  Particles() :id_max_(-1) {}
  /* pass out the specific data of a particle as needed all across board */
  inline const ParticleData &data(int id);
  inline ParticleData * data_pointer(int id);
  /* pass out the type of a specific particle */
  ParticleType type(int id);
  /* pass out the specific type */
  inline ParticleType particle_type(int id);
  /* pass out decay modes of this particle type */
  inline DecayModes decay_modes(int pdg);
  /* inserts new data or type or decay modes */
  inline int id_max(void);
  inline int add_data(const ParticleData &particle_data);
  inline void add_type(const ParticleType &particle_type, int pdg);
  inline void add_decaymodes(const DecayModes &new_decay_modes, int pdg);
  /* add a range of particles */
  inline void create(size_t number, int pdg);
  /* remove the particle */
  inline void remove(int id);
  /* map methods that directly apply on the ParticleData */
  inline size_t size(void) const;
  inline bool empty(void) const;
  inline size_t count(int i) const;
  /* map methods that directly apply on the ParticleType */
  inline size_t types_size(void) const;
  inline bool types_empty(void) const;
  /* return time of the computational frame */
  inline double time(void) const;
  /* iterators */
  inline std::map<int, ParticleData>::iterator begin(void);
  inline std::map<int, ParticleData>::iterator end(void);
  inline std::map<int, ParticleData>::const_iterator cbegin(void) const;
  inline std::map<int, ParticleData>::const_iterator cend(void) const;
  inline std::map<int, ParticleType>::const_iterator types_cbegin(void) const;
  inline std::map<int, ParticleType>::const_iterator types_cend(void) const;

  private:
    /* Highest id of a given particle */
    int id_max_;
    /* dynamic data of the particles a map between it's id and data */
    std::map<int, ParticleData> data_;
    /* a map between pdg and correspoding static data of the particles */
    std::map<int, ParticleType> types_;
    /* a map between pdg and corresponding decay modes */
    std::map<int, DecayModes> all_decay_modes_;
    /* google style recommendation */
    DISALLOW_COPY_AND_ASSIGN(Particles);
};

/* return the data of a specific particle */
inline const ParticleData &Particles::data(int particle_id) {
  return data_[particle_id];
}

/* return the pointer to the data of a specific particle */
inline ParticleData* Particles::data_pointer(int particle_id) {
  return &data_[particle_id];
}

/* return the type of a specific particle */
inline ParticleType Particles::type(int particle_id) {
  return types_[data_[particle_id].pdgcode()];
}

/* return a specific type */
inline ParticleType Particles::particle_type(int pdgcode) {
  return types_[pdgcode];
}

/* return the decay modes of specific type */
inline DecayModes Particles::decay_modes(int pdg) {
  return all_decay_modes_[pdg];
}

/* add a new particle data */
inline int Particles::add_data(ParticleData const &particle_data) {
  id_max_++;
  data_.insert(std::pair<int, ParticleData>(id_max_, particle_data));
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

inline size_t Particles::count(int i) const {
  return data_.count(i);
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

#endif  // SRC_INCLUDE_PARTICLES_H_
