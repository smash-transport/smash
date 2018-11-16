/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/processbranch.h"

#include <limits>

#include "smash/particledata.h"

namespace smash {

bool is_string_soft_process(ProcessType p) {
  return p == ProcessType::StringSoftSingleDiffractiveAX ||
         p == ProcessType::StringSoftSingleDiffractiveXB ||
         p == ProcessType::StringSoftDoubleDiffractive ||
         p == ProcessType::StringSoftAnnihilation ||
         p == ProcessType::StringSoftNonDiffractive;
}

ParticleList ProcessBranch::particle_list() const {
  ParticleList l;
  l.reserve(particle_number());
  for (const auto& type : particle_types()) {
    l.push_back(ParticleData{*type});
  }
  return l;
}

double ProcessBranch::threshold() const {
  if (threshold_ < 0.) {
    /* Sum up the (minimum) masses of all final-state particles
     * this requires double-precision to ensure that the sum is never
     * smaller than the real sum would be without rounding. */
    double thr = 0.;
    for (const auto& type : particle_types()) {
      thr += type->min_mass_kinematic();
    }
    /* This may round up or down. Up is good. If down
     * we must add one ULP via 'nextafter'. */
    const double rounded = thr;
    threshold_ =
        rounded < thr
            ? std::nextafter(rounded, std::numeric_limits<double>::max())
            : rounded;
  }
  return threshold_;
}

std::ostream& operator<<(std::ostream& os, const CollisionBranch& cbranch) {
  ProcessType ptype = cbranch.get_type();
  if (ptype == ProcessType::StringSoftSingleDiffractiveAX ||
      ptype == ProcessType::StringSoftSingleDiffractiveXB) {
    os << "1-diff";
  } else if (ptype == ProcessType::StringSoftDoubleDiffractive) {
    os << "2-diff";
  } else if (ptype == ProcessType::StringSoftAnnihilation) {
    os << "BBbar";
  } else if (ptype == ProcessType::StringSoftNonDiffractive) {
    os << "non-diff";
  } else if (ptype == ProcessType::StringHard) {
    os << "hard";
  } else if (ptype == ProcessType::TwoToOne || ptype == ProcessType::TwoToTwo ||
             ptype == ProcessType::Elastic || ptype == ProcessType::Decay) {
    ParticleTypePtrList ptype_list = cbranch.particle_types();
    /* Sorting ensures unique name for every channel
     * It avoids duplicates, such as Δ⁰Δ⁺⁺ and Δ⁺⁺Δ⁰,
     * which actually occur in SMASH, because of the way channels are added:
     * for example one channel can be added twice with halved cross-section. */
    std::sort(ptype_list.begin(), ptype_list.end());
    for (const auto& type : ptype_list) {
      os << type->name();
    }
  } else {
    os << ptype;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, ProcessType process_type) {
  switch (process_type) {
    case ProcessType::None:
      os << "None";
      break;
    case ProcessType::Elastic:
      os << "Elastic";
      break;
    case ProcessType::TwoToOne:
      os << "TwoToOne";
      break;
    case ProcessType::TwoToTwo:
      os << "TwoToTwo";
      break;
    case ProcessType::StringSoftSingleDiffractiveAX:
    case ProcessType::StringSoftSingleDiffractiveXB:
    case ProcessType::StringSoftDoubleDiffractive:
    case ProcessType::StringSoftAnnihilation:
    case ProcessType::StringSoftNonDiffractive:
      os << "Soft String Excitation";
      break;
    case ProcessType::StringHard:
      os << "Hard String via Pythia";
      break;
    case ProcessType::Decay:
      os << "Decay";
      break;
    case ProcessType::Wall:
      os << "Wall";
      break;
    default:
      os.setstate(std::ios_base::failbit);
  }
  return os;
}

}  // namespace smash
