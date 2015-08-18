/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/isoparticletype.h"

#include <algorithm>

#include "include/forwarddeclarations.h"
#include "include/logging.h"

namespace Smash {

static IsoParticleTypeList iso_type_list;

IsoParticleType::IsoParticleType(const std::string &n, float m, float w, int i)
                                : name_(n), mass_(m), width_(w), isospin_(i) {}

const IsoParticleType& IsoParticleType::find(
                                                      const std::string &name) {
  const auto found = std::lower_bound(
      iso_type_list.begin(), iso_type_list.end(), name,
      [](const IsoParticleType &l, const std::string &r) {
        return l.name() < r;
      });
  if (found == iso_type_list.end() || found->name() != name) {
    throw ParticleNotFoundFailure("Isospin multiplet " + name + " not found!");
  }
  return *found;
}

IsoParticleType& IsoParticleType::find_private(const std::string &name) {
  auto found = std::lower_bound(
      iso_type_list.begin(), iso_type_list.end(), name,
      [](const IsoParticleType &l, const std::string &r) {
        return l.name() < r;
      });
  if (found == iso_type_list.end() || found->name() != name) {
    throw std::runtime_error("Isospin multiplet " + name +
                             " not found (privately)!");
  }
  return *found;
}

bool IsoParticleType::exists(const std::string &name) {
  const auto found = std::lower_bound(
      iso_type_list.begin(), iso_type_list.end(), name,
      [](const IsoParticleType &l, const std::string &r) {
        return l.name() < r;
      });
  if (found != iso_type_list.end()) {
    return found->name() == name;
  }
  return false;
}

/* Construct the name-string for an isospin multiplet from the given
 * name-string for the particle. */
static std::string multiplet_name(std::string name) {
  if (name.find("⁺⁺") != std::string::npos) {
    return name.substr(0, name.length() - sizeof("⁺⁺") + 1);
  } else if (name.find("⁺") != std::string::npos) {
    return name.substr(0, name.length() - sizeof("⁺") + 1);
  } else if (name.find("⁻⁻") != std::string::npos) {
    return name.substr(0, name.length() - sizeof("⁻⁻") + 1);
  } else if (name.find("⁻") != std::string::npos) {
    return name.substr(0, name.length() - sizeof("⁻") + 1);
  } else if (name.find("⁰") != std::string::npos) {
    return name.substr(0, name.length() - sizeof("⁰") + 1);
  } else {
    return name;
  }
}

bool IsoParticleType::has_anti_multiplet() const {
  if (states_[0]->has_antiparticle()) {
    ParticleTypePtr anti = states_[0]->get_antiparticle();
    return multiplet_name(states_[0]->name()) != multiplet_name(anti->name());
  } else {
    return false;
  }
}

const ParticleTypePtr IsoParticleType::find_state(
                                                        const std::string &n) {
  const IsoParticleType &multiplet = IsoParticleType::find(multiplet_name(n));
  auto found = std::find_if(
    multiplet.states_.begin(), multiplet.states_.end(),
    [&n](ParticleTypePtr p) {
      return p->name() == n;
    });
  if (found == multiplet.states_.end() || (*found)->name() != n) {
    throw std::runtime_error("Isospin state " + n + " not found!");
  }
  return *found;
}

void IsoParticleType::create_multiplet(const ParticleType &type) {
  const auto &log = logger<LogArea::ParticleType>();

  // create multiplet if it does not exist yet
  std::string multiname = multiplet_name(type.name());
  if (!exists(multiname)) {
    iso_type_list.emplace_back(multiname, type.mass(), type.width_at_pole(),
                               type.isospin());
    log.debug() << "Creating isospin multiplet " << multiname
                << " [ I = " << type.isospin() << "/2, m = " << type.mass()
                << ", Γ = " << type.width_at_pole() << " ]";
  }

  // sort the iso-type list by name
  std::sort(iso_type_list.begin(), iso_type_list.end(),
            [](const IsoParticleType &l, const IsoParticleType &r) {
              return l.name() < r.name();
            });

  // add the specific type to the multiplet
  IsoParticleType &multiplet = find_private(multiname);
  multiplet.add_state(type);
}

}  // namespace Smash
