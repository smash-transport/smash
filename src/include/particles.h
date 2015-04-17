/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_PARTICLES_H_
#define SRC_INCLUDE_PARTICLES_H_

#include <vector>
#include <type_traits>

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
   * Return a copy of all particles as a std::vector<ParticleData>.
   */
  ParticleList copy_to_vector() const {
    if (dirty_.empty()) {
      return {&data_[0], &data_[data_size_]};
    }
    return {begin(), end()};
  }

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
  size_t size() const { return data_size_ - dirty_.size(); }

  /// empty() check of the ParticleData map
  bool is_empty() const { return data_size_ == 0; }

  /** return time of the computational frame
   *
   * \return computation time which is reduced by the start up time
   *
   * \fpPrecision Why \c double?
   */
  double time() const {
    assert(!is_empty());
    return front().position().x0();
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
    if (data_size_ <= copy.index_) {
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
  void remove(const ParticleData &p);

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
  template <typename T>
  class iterator : public std::iterator<std::bidirectional_iterator_tag, T> {
    friend class Particles;

   public:
    using value_type = typename std::remove_const<T>::type;
    using pointer = typename std::add_pointer<T>::type;
    using reference = typename std::add_lvalue_reference<T>::type;
    using const_pointer = typename std::add_const<pointer>::type;
    using const_reference = typename std::add_const<reference>::type;

   private:
    iterator(pointer p) : ptr_(p) {}  // NOLINT(runtime/explicit)
    pointer ptr_;

   public:
    iterator &operator++() {
      do {
        ++ptr_;
      } while (ptr_->hole_);
      return *this;
    }
    iterator operator++(int) {
      iterator old = *this;
      operator++();
      return old;
    }

    iterator &operator--() {
      do {
        --ptr_;
      } while (ptr_->hole_);
      return *this;
    }
    iterator operator--(int) {
      iterator old = *this;
      operator--();
      return old;
    }

    reference operator*() { return *ptr_; }
    const_reference operator*() const { return *ptr_; }

    pointer operator->() { return ptr_; }
    const_pointer operator->() const { return ptr_; }

    bool operator==(const iterator &rhs) const { return ptr_ == rhs.ptr_; }
    bool operator!=(const iterator &rhs) const { return ptr_ != rhs.ptr_; }
    bool operator< (const iterator &rhs) const { return ptr_ <  rhs.ptr_; }
    bool operator> (const iterator &rhs) const { return ptr_ >  rhs.ptr_; }
    bool operator<=(const iterator &rhs) const { return ptr_ <= rhs.ptr_; }
    bool operator>=(const iterator &rhs) const { return ptr_ >= rhs.ptr_; }
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
  iterator<ParticleData> begin() {
    ParticleData *first = &data_[0];
    while (first->hole_) {
      ++first;
    }
    return first;
  }
  /// const overload of the above
  iterator<const ParticleData> begin() const {
    ParticleData *first = &data_[0];
    while (first->hole_) {
      ++first;
    }
    return first;
  }

  /**
   * Returns an iterator pointing behind the last particle in the list. Use it
   * to iterate over all particles in the list.
   */
  iterator<ParticleData> end() { return &data_[data_size_]; }
  /// const overload of the above
  iterator<const ParticleData> end() const { return &data_[data_size_]; }

  /// Returns a const begin iterator.
  iterator<const ParticleData> cbegin() const { return begin(); }
  /// Returns a const end iterator.
  iterator<const ParticleData> cend() const { return end(); }

  /**
   * \ingroup logging
   * Print effective mass and type name for all particles to the stream.
   */
  friend std::ostream &operator<<(std::ostream &out, const Particles &p);

  /////////////////////// deprecated functions ///////////////////////////

  SMASH_DEPRECATED("use the begin() and end() iterators of Particles directly")
  Particles &data() { return *this; }
  SMASH_DEPRECATED("use the begin() and end() iterators of Particles directly")
  const Particles &data() const { return *this; }

  SMASH_DEPRECATED("don't reference particles by id") ParticleData
      &data(int id) {
    for (ParticleData &x : *this) {
      if (x.id() == id) {
        return x;
      }
    }
    throw std::out_of_range("missing particle id");
  }

  SMASH_DEPRECATED("don't reference particles by id") const ParticleData
      &data(int id) const {
    for (const ParticleData &x : *this) {
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
    for (auto it = begin(); it != end(); ++it) {
      if (it->id() == id) {
        remove(*it);
        return;
      }
    }
  }

  SMASH_DEPRECATED("don't reference particles by id") bool has_data(
      int id) const {
    for (auto &&x : *this) {
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

  void increase_capacity(unsigned new_capacity);
  inline void ensure_capacity(unsigned to_add);
  inline void copy_in(ParticleData &to, const ParticleData &from);

  unsigned data_size_ = 0u;
  unsigned data_capacity_ = 100u;
  std::unique_ptr<ParticleData[]> data_;

  /**
   * Stores the indexes in data_ that do not hold valid particle data and should
   * be reused when new particles are added.
   */
  std::vector<unsigned> dirty_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_PARTICLES_H_
