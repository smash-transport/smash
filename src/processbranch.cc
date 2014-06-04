/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/processbranch.h"
#include "include/particledata.h"

namespace Smash {

ParticleList ProcessBranch::particle_list() const {
  ParticleList l;
  l.reserve(pdg_list_.size());
  for (auto pdgcode : pdg_list_) {
    l.push_back(ParticleData{ParticleType::find(pdgcode)});
  }
  return std::move(l);
}

}  // namespace Smash
