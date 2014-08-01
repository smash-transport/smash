/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARTICLES_H_
#define SRC_INCLUDE_PARTICLES_H_

#include "macros.h"
#include "particledata.h"
#include "particletype.h"
#include "pdgcode.h"

#include <map>

namespace Smash {

/**
 * \ingroup data
 *
 * The Particles class abstracts the storage and manipulation of particles.
 *
 * There is one Particles object per Experiment. It stores
 * the data about all existing particles in the experiment (ParticleData).
 *
 * \note
 * The Particles object cannot be copied, because it does not make sense
 * semantically. Move semantics make sense and can be implemented when needed.
 */
class Particles {
  /**
   * Adapter class for making iterations over particle data and types easier to
   * use.
   *
   * This is an implementation detail and need not be understood for use of the
   * Particles class.
   */
  template <typename T>
  class MapIterationAdapter {
    /// the type T::begin() returns
    using map_iterator = decltype(std::declval<T &>().begin());
    /// the type T::cbegin() returns
    using const_map_iterator = decltype(std::declval<T &>().cbegin());
    /// the value type of T (with const added if T is const)
    using mapped_type = typename std::conditional<
        std::is_const<T>::value, const typename T::mapped_type,
        typename T::mapped_type>::type;

   public:
    MapIterationAdapter(T *map) : map_(map) {}  // NOLINT(runtime/explicit)

    /**
     * Adapter class that modifies the normal map iterator to return only the
     * value instead of the key/value pair.
     */
    class iterator : public map_iterator {
     public:
      iterator(map_iterator it)  // NOLINT(runtime/explicit)
          : map_iterator(it) {}

      /**
       * overwritten dereference operator to return only the value instead of
       * the
       * key/value pair.
       */
      mapped_type &operator*() { return map_iterator::operator*().second; }
      /// const overload of the above
      const mapped_type &operator*() const {
        return map_iterator::operator*().second;
      }

      mapped_type *operator->() { return &map_iterator::operator*().second; }
      const mapped_type *operator->() const {
        return &map_iterator::operator*().second;
      }
    };

    /**
     * Adapter class that modifies the normal map const_iterator to return only
     * the value instead of the key/value pair.
     */
    class const_iterator : public const_map_iterator {
     public:
      const_iterator(const_map_iterator it)  // NOLINT(runtime/explicit)
          : const_map_iterator(it) {}

      /**
       * overwritten dereference operator to return only the value instead of
       * the key/value pair.
       */
      const mapped_type &operator*() const {
        return const_map_iterator::operator*().second;
      }

      const mapped_type *operator->() const {
        return &const_map_iterator::operator*().second;
      }
    };


    /// returns an adapted iterator to the begin iterator of the map
    iterator begin() const { return map_->begin(); }
    /// returns an adapted iterator to the end iterator of the map
    iterator end() const { return map_->end(); }

    /// returns an adapted iterator to the begin iterator of the map
    const_iterator cbegin() const { return map_->cbegin(); }
    /// returns an adapted iterator to the end iterator of the map
    const_iterator cend() const { return map_->cend(); }

   private:
    /// points to the map the adapter class wraps; needed for calls to begin/end
    T *map_;
  };

  /** ParticleDataMap is the prime accessor for ParticleData
   *
   * It maps the unique Particle ID to its volatile data.
   */
  using ParticleDataMap = std::map<int, ParticleData>;

 public:
  /**
   * Set up the Particles object.
   *
   * This initializes all the members. The object is ready for usage right after
   * construction.
   */
  Particles();

  /// Cannot be copied
  Particles(const Particles &) = delete;
  /// Cannot be copied
  Particles &operator=(const Particles &) = delete;

  /**
   * Use the returned object for iterating over all ParticleData objects.
   *
   * You can use this object for use with range-based for:
   * \code
   * for (ParticleData &data : particles->data()) {
   *   ...
   * }
   * \endcode
   *
   * \returns Opaque object that provides begin/end functions for iteration of
   * all ParticleData objects.
   */
  MapIterationAdapter<ParticleDataMap> data() { return &data_; }
  /// const overload of the above
  MapIterationAdapter<const ParticleDataMap> data() const { return &data_; }

  /**
   * Return the specific data of a particle according to its id
   *
   * \throws std::out_of_range If there is no particle with the given \p id.
   */
  inline ParticleData &data(int id) { return data_.at(id); }
  /**
   * Return the specific data of a particle according to its id
   *
   * \throws std::out_of_range If there is no particle with the given \p id.
   */
  inline const ParticleData &data(int id) const;

  /// return the highest used id
  inline int id_max(void) const;
  /// inserts a new particle and returns its id
  inline int add_data(const ParticleData &particle_data);
  /// add a range of particles
  inline void create(size_t number, PdgCode pdg);
  /* add one particle and return pointer to it */
  inline ParticleData& create(const PdgCode pdg);
  /// remove a specific particle
  inline void remove(int id);
  /// size() of the ParticleData map
  inline size_t size(void) const;
  /// empty() check of the ParticleData map
  inline bool empty(void) const;
  /// check the existence of an element in the ParticleData map
  inline bool has_data(int id) const;
  /// return time of the computational frame
  inline double time(void) const;

  /// \ingroup exception
  struct LoadFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };
  /// \ingroup exception
  struct ParseError : public LoadFailure {
    using LoadFailure::LoadFailure;
  };

  /** Reset member data to the state the object had when the constructor
   * returned.
   */
  void reset();

 private:

  /// Highest id of a given particle
  int id_max_ = -1;
  /**
   * dynamic data of the particles a map between its id and data
   *
   * A map structure is used as particles decay and hence this
   * a swiss cheese over the runtime of SMASH. Also we want direct
   * lookup of the corresponding particle id with its data.
   */
  ParticleDataMap data_;
};

/* return the data of a specific particle */
inline const ParticleData &Particles::data(int particle_id) const {
  return data_.at(particle_id);
}

/* add a new particle data and return the id of the new particle */
inline int Particles::add_data(ParticleData const &particle_data) {
  id_max_++;
  data_.insert(std::make_pair(id_max_, particle_data));
  data_.at(id_max_).set_id(id_max_);
  return id_max_;
}

/* create a bunch of particles */
inline void Particles::create(size_t number, PdgCode pdgcode) {
  /* fixed pdgcode and no collision yet */
  ParticleData particle(ParticleType::find(pdgcode));
  particle.set_collision(0);
  for (size_t i = 0; i < number; i++) {
    id_max_++;
    particle.set_id(id_max_);
    data_.insert(std::make_pair(id_max_, particle));
  }
}

/* create a bunch of particles */
inline ParticleData& Particles::create(PdgCode pdgcode) {
  /* fixed pdgcode and no collision yet */
  ParticleData particle(ParticleType::find(pdgcode));
  particle.set_collision(0);
  id_max_++;
  particle.set_id(id_max_);
  data_.insert(std::make_pair(id_max_, particle));
  return data_.at(id_max_);
}

/* return the highest used id */
inline int Particles::id_max() const {
  return id_max_;
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

/* check the existence of an element in the ParticleData map */
inline bool Particles::has_data(int id) const {
  return data_.find(id) != data_.end();
}

/* return computation time which is reduced by the start up time */
inline double Particles::time() const {
  return data_.begin()->second.position().x0();
}

}  // namespace Smash

#endif  // SRC_INCLUDE_PARTICLES_H_
