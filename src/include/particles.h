/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARTICLES_H_
#define SRC_INCLUDE_PARTICLES_H_

#include <vector>

#include "macros.h"
#include "particledata.h"
#include "particletype.h"
#include "pdgcode.h"

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
   * Inserts the particle \p into the list of particles.
   * The argument \p will afterwards not be a valid copy of a particle of the
   * internal list. I.e.
   * \code
   * ParticleData pd(type);
   * particles.insert(pd);
   * particles.is_valid(pd); // returns false
   * \endcode
   */
  void insert(const ParticleData &p);

  /// Add \p number particles of the same type (\p pdg).
  void create(size_t number, PdgCode pdg);

  /// Add one particle of the given \p pdg code and return a reference to it
  ParticleData &create(const PdgCode pdg);

  /// Returns the current number of particles.
  size_t size() const { return data_.size() - dirty_.size() - 1; }

  /// empty() check of the ParticleData map
  bool is_empty() const { return data_.size() == 1; }

  /** return time of the computational frame
   *
   * \return computation time which is reduced by the start up time
   *
   * \fpPrecision Why \c double?
   */
  double time() const {
    assert(!is_empty());
    return data_.front().position().x0();
  }

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

  /**
   * Return whether the ParticleData copy is still a valid copy of the one
   * stored in the Particles object. If not, then the particle has interacted
   * between the copy and the call to is_valid.
   */
  bool is_valid(const ParticleData &copy) const {
    if (static_cast<int>(data_.size()) <= copy.index_) {
      return false;
    }
    return data_[copy.index_].id() ==
               copy.id()  // Check if the particles still exists. If it decayed
                          // or scattered inelastically it is gone.
           &&
           data_[copy.index_].id_process() ==
               copy.id_process();  // If the particle has scattered elastically,
                                   // its id_process has changed and we consider
                                   // it invalid.
  }

  /**
   * Returns a reference to the original particle in the list for the \p copy
   * that originated from the Particles list.
   *
   * \param copy This object must be a valid copy from the Particles list. On
   *             debug builds this is ensured via an assertion, on non-debug
   *             builds an invalid copy will lead to undefined behavior.
   */
  ParticleData &original(const ParticleData &copy) {
    assert(is_valid(copy));
    return data_[copy.index_];
  }
  /// const overload of the above
  const ParticleData &original(const ParticleData &copy) const {
    assert(is_valid(copy));
    return data_[copy.index_];
  }

  /**
   * Remove the given particle \p p from the list. The argument \p p must be a
   * valid copy obtained from Particles, i.e. a call to \ref is_valid must
   * return \c true.
   */
  void remove(const ParticleData &p) {
    assert(is_valid(p));
    const std::size_t index = p.index_;
    data_[index].set_id(-1);
    data_[index].index_ = ParticleData::invalid_index;
    dirty_.push_back(index);
  }

  /**
   * Replace the particles in \p to_remove with the particles in \p to_add in
   * the list of current particles. The particles in \p to_remove must be valid
   * copies obtained from Particles. The particles in \p to_add will not be
   * modified by this function call and therefore not be valid copies of the new
   * particles in the Particles list.
   */
  void replace(const ParticleList &to_remove, const ParticleList &to_add);

  /**
   * \internal
   * Iterator type that skips over the holes in data_.
   */
  template <typename Base>
  class iterator : public Base {
    // TODO(mkretz): operator[] is wrong! And it's not a RandomAccessIterator
    // either.
   public:
    typedef ParticleData value_type;
    typedef ParticleData *pointer;
    typedef ParticleData &reference;
    typedef const ParticleData *const_pointer;
    typedef const ParticleData &const_reference;

    iterator(Base it) : Base(it) {}  // NOLINT(runtime/explicit)

    iterator &operator++() {
      do {
        Base::operator++();
      } while ((*this)->index_ == ParticleData::invalid_index);
      return *this;
    }
    iterator operator++(int) {
      iterator old = *this;
      operator++();
      return old;
    }

    iterator &operator--() {
      do {
        Base::operator--();
      } while ((*this)->index_ == ParticleData::invalid_index);
      return *this;
    }
    iterator operator--(int) {
      iterator old = *this;
      operator--();
      return old;
    }
  };

  /// Returns a reference to the first particle in the list.
  ParticleData &front() { return *begin(); }
  /// const overload of the above
  const ParticleData &front() const { return *begin(); }

  /// Returns a reference to the last particle in the list.
  ParticleData &back() { return *(--end()); }
  /// const overload of the above
  const ParticleData &back() const { return *(--end()); }

  /**
   * Returns an iterator pointing to the first particle in the list. Use it to
   * iterate over all particles in the list.
   */
  iterator<ParticleList::iterator> begin() {
    auto it = data_.begin();
    if (size() != 0) {
      while (it->index_ == ParticleData::invalid_index) {
        ++it;
      }
    }
    return it;
  }
  /// const overload of the above
  iterator<ParticleList::const_iterator> begin() const {
    auto it = data_.begin();
    if (size() != 0) {
      while (it->index_ == ParticleData::invalid_index) {
        ++it;
      }
    }
    return it;
  }

  /**
   * Returns an iterator pointing behind the last particle in the list. Use it
   * to iterate over all particles in the list.
   */
  iterator<ParticleList::iterator> end() {
    return --data_.end();
  }
  /// const overload of the above
  iterator<ParticleList::const_iterator> end() const {
    return --data_.end();
  }

  /// Returns a const begin iterator.
  iterator<ParticleList::const_iterator> cbegin() const { return begin(); }
  /// Returns a const end iterator.
  iterator<ParticleList::const_iterator> cend() const { return end(); }

  /**
   * \ingroup logging
   * Print effective mass and type name for all particles to the stream.
   */
  friend std::ostream &operator<<(std::ostream &out, const Particles &p);

  /////////////////////// deprecated functions ///////////////////////////

  SMASH_DEPRECATED("use the begin() and end() iterators of Particles directly")
  std::vector<ParticleData> &data() { return data_; }
  SMASH_DEPRECATED("use the begin() and end() iterators of Particles directly")
  const std::vector<ParticleData> &data() const { return data_; }

  SMASH_DEPRECATED("don't reference particles by id") ParticleData
      &data(int id) {
    return data_.at(id);
  }

  SMASH_DEPRECATED("don't reference particles by id") const ParticleData
      &data(int id) const {
    for (auto &&x : data_) {
      if (x.id() == id) {
        return x;
      }
    }
    throw std::out_of_range("missing particle id");
  }

  SMASH_DEPRECATED("don't reference particles by id") int id_max(void) const {
    return id_max_;
  }

  SMASH_DEPRECATED("don't reference particles by id") void remove(int id) {
    for (auto it = data_.begin(); it != data_.end(); ++it) {
      if (it->id() == id) {
        data_.erase(it);
        return;
      }
    }
  }

  SMASH_DEPRECATED("don't reference particles by id") bool has_data(
      int id) const {
    for (auto &&x : data_) {
      if (x.id() == id) {
        return true;
      }
    }
    return false;
  }

  SMASH_DEPRECATED("use insert instead") int add_data(
      const ParticleData &particle_data) {
    insert(particle_data);
    return id_max_;
  }

  SMASH_DEPRECATED("use is_empty instead") bool empty() const {
    return is_empty();
  }

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
  ParticleList data_;

  /**
   * Stores the indexes in data_ that do not hold valid particle data and should
   * be reused when new particles are added.
   */
  std::vector<int> dirty_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_PARTICLES_H_
