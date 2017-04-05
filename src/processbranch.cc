/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/processbranch.h"

#include <limits>

#include "include/particledata.h"

namespace Smash {

ParticleList ProcessBranch::particle_list() const {
  ParticleList l;
  l.reserve(particle_number());
  for (const auto &type : particle_types()) {
    l.push_back(ParticleData{*type});
  }
  return l;
}

float ProcessBranch::threshold() const {
  if (threshold_ < 0.f) {
    // Sum up the (minimum) masses of all final-state particles
    double thr = 0.;  // this requires double-precision to ensure that the sum is
                      // never smaller than the real sum would be without rounding
    for (const auto &type : particle_types()) {
      thr += type->minimum_mass();
    }
    const float rounded = thr;  // This may round up or down. Up is good. If down
                                // we must add one ULP via 'nextafter'.
    threshold_ =  rounded < thr
                  ? std::nextafter(rounded, std::numeric_limits<float>::max())
                  : rounded;
  }
  return threshold_;
}

std::ostream& operator<< (std::ostream& os, ProcessType process_type) {
  switch (process_type) {
    case ProcessType::None     : os << "None";     break;
    case ProcessType::Elastic  : os << "Elastic";  break;
    case ProcessType::TwoToOne : os << "TwoToOne"; break;
    case ProcessType::TwoToTwo : os << "TwoToTwo"; break;
    case ProcessType::String   : os << "String";   break;
    case ProcessType::Decay    : os << "Decay";    break;
    case ProcessType::Wall     : os << "Wall";     break;
    default: os.setstate(std::ios_base::failbit);
  }
  return os;
}

}  // namespace Smash
