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

Particles::Particles() { data_.reserve(50);
  data_.emplace_back(ParticleType::list_all().front(), -1, -1);
}

void Particles::insert(const ParticleData &p) {
  if (dirty_.empty()) {
    data_.back() = p;
    data_.back().set_id(++id_max_);
    data_.back().index_ = data_.size() - 1;
    data_.emplace_back(ParticleType::list_all().front(), -1, -1);
  } else {
    const auto offset = dirty_.back();
    dirty_.pop_back();
    data_[offset] = p;
    data_[offset].set_id(++id_max_);
    data_[offset].index_ = offset;
  }
}

void Particles::create(size_t number, PdgCode pdg) {
  const auto &type = ParticleType::find(pdg);
  while (number && !dirty_.empty()) {
    const auto offset = dirty_.back();
    dirty_.pop_back();
    data_[offset] = ParticleData(type, ++id_max_, offset);
    --number;
  }
  if (number) {
    data_.reserve(data_.size() + number);
    data_.back() = ParticleData(type, ++id_max_, data_.size() - 1);
    for (auto i = number - 1; i; --i) {
      data_.emplace_back(type, ++id_max_, data_.size());
    }
    data_.emplace_back(ParticleType::list_all().front(), -1, -1);
  }
}

ParticleData &Particles::create(const PdgCode pdg) {
  const auto &type = ParticleType::find(pdg);
  if (dirty_.empty()) {
    data_.back() = ParticleData(type, ++id_max_, data_.size() - 1);
    data_.emplace_back(ParticleType::list_all().front(), -1, -1);
    return data_[data_.size() - 2];
  } else {
    const auto offset = dirty_.back();
    dirty_.pop_back();
    data_[offset] = ParticleData(type, ++id_max_, offset);
    return data_[offset];
  }
}

void Particles::replace(const ParticleList &to_remove,
                        const ParticleList &to_add) {
  std::size_t i = 0;
  for (; i < std::min(to_remove.size(), to_add.size()); ++i) {
    assert(is_valid(to_remove[i]));
    const auto index = to_remove[i].index_;
    data_[index] = to_add[i];
    data_[index].set_id(++id_max_);
    data_[index].index_ = index;
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
  data_.clear();
  data_.emplace_back(ParticleType::list_all().front(), -1, -1);
}

std::ostream &operator<<(std::ostream &out, const Particles &particles) {
  using namespace std;
  out << particles.size() << " Particles:\n";
  int n = 0;
  for (const auto &p : particles.data_) {
    if (p.id() < 0) {
      out << "------  ";
    } else {
      out << setw(5) << setprecision(3) << p.momentum().abs3()
          << p.type().name();
    }
    if ((++n & 15) == 0) {
      out << '\n';
    }
  }
  return out;
}

}  // namespace Smash
