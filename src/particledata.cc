/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/particledata.h"

#include <iomanip>
#include <iostream>
#include <vector>

#include "smash/constants.h"
#include "smash/iomanipulators.h"

namespace smash {

double ParticleData::effective_mass() const {
  const double m_pole = pole_mass();
  if (m_pole < really_small) {
    // prevent numerical problems with massless or very light particles
    return m_pole;
  } else {
    return momentum().abs();
  }
}

void ParticleData::set_history(int ncoll, uint32_t pid, ProcessType pt,
                               double time_last_coll,
                               const ParticleList &plist) {
  if (pt != ProcessType::Wall) {
    history_.collisions_per_particle = ncoll;
    history_.time_last_collision = time_last_coll;
  }
  history_.id_process = pid;
  history_.process_type = pt;
  switch (pt) {
    case ProcessType::Decay:
    case ProcessType::Wall:
      // only store one parent
      history_.p1 = plist[0].pdgcode();
      history_.p2 = 0x0;
      break;
    case ProcessType::Elastic:
    case ProcessType::TwoToOne:
    case ProcessType::TwoToTwo:
    case ProcessType::StringSoftSingleDiffractiveAX:
    case ProcessType::StringSoftSingleDiffractiveXB:
    case ProcessType::StringSoftDoubleDiffractive:
    case ProcessType::StringSoftAnnihilation:
    case ProcessType::StringSoftNonDiffractive:
    case ProcessType::StringHard:
      // store two parent particles
      history_.p1 = plist[0].pdgcode();
      history_.p2 = plist[1].pdgcode();
      break;
    case ProcessType::Thermalization:
    case ProcessType::None:
      // nullify parents
      history_.p1 = 0x0;
      history_.p2 = 0x0;
      break;
  }
}

double ParticleData::current_xsec_scaling_factor(double time_until_collision,
                                                 double power) const {
  double total_time = position_.x0() + time_until_collision;
  if (power <= 0.) {
    // use a step function to form particles
    if (total_time < formation_time_) {
      return cross_section_scaling_factor_;
    }
    return 1.;
  }
  if (formation_time_ <= total_time) {
    return 1.;
  }
  if (begin_formation_time_ >= total_time) {
    return cross_section_scaling_factor_;
  }
  return cross_section_scaling_factor_ +
         (1. - cross_section_scaling_factor_) *
             std::pow((total_time - begin_formation_time_) /
                          (formation_time_ - begin_formation_time_),
                      power);
}

std::ostream &operator<<(std::ostream &out, const ParticleData &p) {
  out.fill(' ');
  return out
#ifdef NDEBUG
         << std::setw(5) << p.type().pdgcode()
#else
         << p.type().name()
#endif
         << std::right << "{id:" << field<6> << p.id()
         << ", process:" << field<4> << p.id_process()
         << ", pos [fm]:" << p.position() << ", mom [GeV]:" << p.momentum()
         << ", formation time [fm]:" << p.formation_time()
         << ", cross section scaling factor:"
         << p.cross_section_scaling_factor() << "}";
}

std::ostream &operator<<(std::ostream &out, const ParticleList &particle_list) {
  auto column = out.tellp();
  out << '[';
  for (const auto &p : particle_list) {
    if (out.tellp() - column >= 201) {
      out << '\n';
      column = out.tellp();
      out << ' ';
    }
    out << std::setw(5) << std::setprecision(3) << p.momentum().abs3()
        << p.type().name();
  }
  return out << ']';
}

std::ostream &operator<<(std::ostream &out,
                         const PrintParticleListDetailed &particle_list) {
  bool first = true;
  out << '[';
  for (const auto &p : particle_list.list) {
    if (first) {
      first = false;
    } else {
      out << "\n ";
    }
    out << p;
  }
  return out << ']';
}

}  // namespace smash
