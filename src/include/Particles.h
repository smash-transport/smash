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

#include "../include/ParticleData.h"
#include "../include/ParticleType.h"

class Particles {
  public:
  /* Use improbable values for default constructor */
  Particles() :id_max_(-1) {}
  /* direct access to both data and types */
  inline std::map<int, ParticleData> data(void) const;
  inline std::map<int, ParticleType> types(void) const;
  /* pass out the specific data of a particle as needed all across board */
  inline ParticleData data(int id);
  inline ParticleData * data_pointer(int id);
  /* pass out the type of a specific particle */
  inline ParticleType type(int id);
  /* pass out the specific type */
  inline ParticleType particle_type(int id);
  /* inserts new data or type */
  inline void add_data(void);
  inline int new_data(void);
  inline void add_data(const ParticleData &particle_data);
  inline void add_type(const ParticleType &particle_type, int pdg_code);
  /* remove the particle */
  inline void remove(int id);
  /* return number of particles */
  inline size_t size(void) const;
  /* return time of the computanional frame */
  inline double time(void) const;

  private:
    /* Highest id of a given particle */
    int id_max_;
    /* dynamic data of the particles a map between it's id and data */
    std::map<int, ParticleData> data_;
    /* a map between pdg and correspoding static data of the particles */
    std::map<int, ParticleType> types_;
};

/* returns the particle data */
inline std::map<int, ParticleData> Particles::data(void) const {
  return data_;
}

/* returns the particle types */
inline std::map<int, ParticleType> Particles::types(void) const {
  return types_;
}

/* return the data of a specific particle */
inline ParticleData Particles::data(int particle_id) {
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
inline ParticleType Particles::particle_type(int type_id) {
  return types_[type_id];
}

/* add a new particle data */
inline void Particles::add_data(ParticleData const &particle_data) {
  id_max_++;
  data_.insert(std::pair<int, ParticleData>(id_max_, particle_data));
}
inline void Particles::add_data(void) {
  id_max_++;
  ParticleData new_particle(id_max_);
  data_.insert(std::pair<int, ParticleData>(id_max_, new_particle));
}
inline int Particles::new_data(void) {
  id_max_++;
  ParticleData new_particle(id_max_);
  data_.insert(std::pair<int, ParticleData>(id_max_, new_particle));
  return id_max_;
}


/* add a new particle type */
inline void Particles::add_type(ParticleType const &TYPE, int pdg) {
  types_.insert(std::pair<int, ParticleType>(pdg, TYPE));
}

/* remove a particle */
inline void Particles::remove(int id) {
  data_.erase(id);
}

/* total number of particles */
inline size_t Particles::size() const {
  return data_.size();
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
