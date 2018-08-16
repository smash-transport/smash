/*
 *
 *    Copyright (c) 2013-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/particles.h"

#include <iomanip>
#include <iostream>

namespace smash {

Particles::Particles() : data_(new ParticleData[data_capacity_]) {
  for (unsigned i = 0; i < data_capacity_; ++i) {
    data_[i].index_ = i;
  }
}

inline void Particles::ensure_capacity(unsigned to_add) {
  if (data_size_ + to_add >= data_capacity_) {
    increase_capacity((data_capacity_ + to_add) * 2u);
    assert(data_size_ + to_add < data_capacity_);
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
  from.copy_to(to);
}

const ParticleData &Particles::insert(const ParticleData &p) {
  if (likely(dirty_.empty())) {
    ensure_capacity(1);
    ParticleData &in_vector = data_[data_size_];
    copy_in(in_vector, p);
    ++data_size_;
    return in_vector;
  } else {
    const auto offset = dirty_.back();
    dirty_.pop_back();
    copy_in(data_[offset], p);
    data_[offset].hole_ = false;
    return data_[offset];
  }
}

void Particles::create(size_t number, PdgCode pdg) {
  const ParticleData pd(ParticleType::find(pdg));
  while (number && !dirty_.empty()) {
    const auto offset = dirty_.back();
    dirty_.pop_back();
    pd.copy_to(data_[offset]);
    data_[offset].id_ = ++id_max_;
    data_[offset].type_ = pd.type_;
    data_[offset].hole_ = false;
    --number;
  }
  if (number) {
    ensure_capacity(number);
    const auto end_ptr = &data_[data_size_ + number];
    for (auto ptr = &data_[data_size_]; ptr < end_ptr; ++ptr) {
      pd.copy_to(*ptr);
      ptr->id_ = ++id_max_;
      ptr->type_ = pd.type_;
    }
    data_size_ += number;
  }
}

ParticleData &Particles::create(const PdgCode pdg) {
  const ParticleData pd(ParticleType::find(pdg));
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
  pd.copy_to(*ptr);
  ptr->id_ = ++id_max_;
  ptr->type_ = pd.type_;
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

void Particles::replace(const ParticleList &to_remove, ParticleList &to_add) {
  std::size_t i = 0;
  for (; i < std::min(to_remove.size(), to_add.size()); ++i) {
    assert(is_valid(to_remove[i]));
    const auto index = to_remove[i].index_;
    copy_in(data_[index], to_add[i]);
    to_add[i].id_ = data_[index].id_;
    to_add[i].index_ = index;
  }
  for (; i < to_remove.size(); ++i) {
    remove(to_remove[i]);
  }
  for (; i < to_add.size(); ++i) {
    const ParticleData &p = insert(to_add[i]);
    to_add[i].id_ = p.id_;
    to_add[i].index_ = p.index_;
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
  out << particles.size() << " Particles:\n";
  for (unsigned i = 0; i < particles.data_size_; ++i) {
    const auto &p = particles.data_[i];
    if (p.id() < 0) {
      out << "------  ";
    } else {
      out << std::setw(5) << std::setprecision(3) << p.momentum().abs3()
          << p.type().name();
    }
    if ((i & 15) == 0) {
      out << '\n';
    }
  }
  return out;
}

}  // namespace smash
