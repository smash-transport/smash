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

IsoParticleType::IsoParticleType(std::string n, float m, float w, int i)
                                : name_(n), mass_(m), width_(w), isospin_(i) {}

SMASH_CONST const IsoParticleType& IsoParticleType::find(std::string name) {
  const auto found = std::lower_bound(
      iso_type_list.begin(), iso_type_list.end(), name,
      [](const IsoParticleType &l, const std::string &r) {
        return l.name() < r;
      });
  if (found == iso_type_list.end() || found->name() != name) {
    throw std::runtime_error("Isospin multiplet " + name + " not found!");
  }
  return *found;
}

IsoParticleType& IsoParticleType::find_private(std::string name) {
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

SMASH_CONST bool IsoParticleType::exists(std::string name) {
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
  char charges[] = "⁺⁻⁰";
  for (unsigned int i = 0; i < strlen(charges); ++i) {
    name.erase (std::remove(name.begin(), name.end(), charges[i]), name.end());
  }
  return name;
}

void IsoParticleType::create_multiplet(const ParticleType &type) {
  const auto &log = logger<LogArea::ParticleType>();

  // create multiplet if it does not exist yet
  std::string multiname = multiplet_name(type.name());
  if (!exists(multiname)) {
    iso_type_list.emplace_back(multiname, type.mass(), type.width_at_pole(),
                               type.isospin());
    log.info() << "Creating isospin multiplet " << multiname
               << " [ I = " << type.isospin() << "/2 ]";
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
