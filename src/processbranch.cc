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
#include "include/resonances.h"

namespace Smash {

ParticleList ProcessBranch::particle_list() const {
  ParticleList l;
  l.reserve(pdg_list_.size());
  for (auto pdgcode : pdg_list_) {
    l.push_back(ParticleData{ParticleType::find(pdgcode)});
  }
  return std::move(l);
}

float ProcessBranch::threshold() const {
  float thr = 0.;
  for (auto pdgcode : pdg_list_) {
    const ParticleType &t = ParticleType::find(pdgcode);
    if (t.is_stable()) {
      thr += t.mass();
    }
    else {
      thr += calculate_minimum_mass(pdgcode);
    }
  }
  return thr;
}

}  // namespace Smash
