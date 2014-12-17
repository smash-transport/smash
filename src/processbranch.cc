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
  l.reserve(particle_types_.size());
  for (const auto &type : particle_types_) {
    l.push_back(ParticleData{*type});
  }
  return std::move(l);
}

float ProcessBranch::threshold() const {
  float thr = 0.;
  /* Sum up the (minimum) masses of all final-state particles. */
  for (const auto &type : particle_types_) {
    thr += type->minimum_mass();
  }
  return thr;
}

}  // namespace Smash
