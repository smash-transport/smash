/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/particles.h"

#include <iomanip>
#include <iostream>

namespace Smash {

Particles::Particles() : data_(new ParticleData[data_capacity_]) {
  for (unsigned i = 0; i < data_capacity_; ++i) {
    data_[i].index_ = i;
  }
}

inline void Particles::ensure_capacity(unsigned to_add) {
  if (data_size_ + to_add >= data_capacity_) {
    increase_capacity((data_capacity_ + to_add) * 2u);
    assert(data_size_ + to_add > data_capacity_);
  }
}

void Particles::increase_capacity(unsigned new_capacity) {
  assert(new_capacity > data_capacity_);
  data_capacity_ = new_capacity;
  std::unique_ptr<ParticleData[]> new_memory(new ParticleData[data_capacity_]);
  unsigned i = 0;
  for (; i < data_size_; ++i) {
    new_memory[i] = data_[i];
  }
  for (; i < data_capacity_; ++i) {
    new_memory[i].index_ = i;
  }
  std::swap(data_, new_memory);
}

inline void Particles::copy_in(ParticleData &to, const ParticleData &from) {
  to.id_ = ++id_max_;
  to.type_ = from.type_;
  to.momentum_ = from.momentum_;
  to.position_ = from.position_;
}

void Particles::insert(const ParticleData &p) {
  if (likely(dirty_.empty())) {
    ensure_capacity(1);
    copy_in(data_[data_size_], p);
    ++data_size_;
  } else {
    const auto offset = dirty_.back();
    dirty_.pop_back();
    copy_in(data_[offset], p);
    data_[offset].hole_ = false;
  }
}

void Particles::create(size_t number, PdgCode pdg) {
  const auto &type = ParticleType::find(pdg);
  while (number && !dirty_.empty()) {
    const auto offset = dirty_.back();
    dirty_.pop_back();
    data_[offset].id_ = ++id_max_;
    data_[offset].type_ = &type;
    data_[offset].hole_ = false;
    --number;
  }
  if (number) {
    ensure_capacity(number);
    const auto end_ptr = &data_[data_size_ + number];
    for (auto ptr = &data_[data_size_]; ptr < end_ptr; ++ptr) {
      ptr->id_ = ++id_max_;
      ptr->type_ = &type;
    }
    data_size_ += number;
  }
}

ParticleData &Particles::create(const PdgCode pdg) {
  const auto &type = ParticleType::find(pdg);
  ParticleData *ptr;
  if (likely(dirty_.empty())) {
    ensure_capacity(1);
    ptr = &data_[data_size_];
    ++data_size_;
  } else {
    const auto offset = dirty_.back();
    dirty_.pop_back();
    ptr = &data_[offset];
    ptr->hole_ = false;
  }
  ptr->id_ = ++id_max_;
  ptr->type_ = &type;
  return *ptr;
}

void Particles::remove(const ParticleData &p) {
  assert(is_valid(p));
  const unsigned index = p.index_;
  if (index == data_size_ - 1) {
    --data_size_;
  } else {
    data_[index].set_id(-1);
    data_[index].hole_ = true;
    dirty_.push_back(index);
  }
}

void Particles::replace(const ParticleList &to_remove,
                        const ParticleList &to_add) {
  std::size_t i = 0;
  for (; i < std::min(to_remove.size(), to_add.size()); ++i) {
    assert(is_valid(to_remove[i]));
    const auto index = to_remove[i].index_;
    copy_in(data_[index], to_add[i]);
  }
  for (; i < to_remove.size(); ++i) {
    remove(to_remove[i]);
  }
  for (; i < to_add.size(); ++i) {
    insert(to_add[i]);
  }
}

void Particles::reset() {
  id_max_ = -1;
  data_size_ = 0;
  for (auto index : dirty_) {
    data_[index].hole_ = false;
  }
  dirty_.clear();
}

std::ostream &operator<<(std::ostream &out, const Particles &particles) {
  using namespace std;
  out << particles.size() << " Particles:\n";
  for (unsigned i = 0; i < particles.data_size_; ++i) {
    const auto &p = particles.data_[i];
    if (p.id() < 0) {
      out << "------  ";
    } else {
      out << setw(5) << setprecision(3) << p.momentum().abs3()
          << p.type().name();
    }
    if ((i & 15) == 0) {
      out << '\n';
    }
  }
  return out;
}

}  // namespace Smash
