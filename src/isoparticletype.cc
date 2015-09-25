/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/isoparticletype.h"

#include <algorithm>

#include "include/forwarddeclarations.h"
#include "include/integrate.h"
#include "include/logging.h"
#include "include/resonances.h"

namespace Smash {

static IsoParticleTypeList iso_type_list;

IsoParticleType::IsoParticleType(const std::string &n, float m, float w, int i)
                                : name_(n), mass_(m), width_(w), isospin_(i) {}

const IsoParticleType& IsoParticleType::find(const std::string &name) {
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

IsoParticleType* IsoParticleType::find(const ParticleType &type) {
  std::string multiname = multiplet_name(type.name());
  IsoParticleType &multiplet = find_private(multiname);
  return &multiplet;
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


double IsoParticleType::get_integral_NR(double srts) {
  if (XS_NR_tabulation == nullptr) {
    // initialize tabulation
    /* TODO(weil): Move this lazy init to a global initialization function,
      * in order to avoid race conditions in multi-threading. */
    Integrator integrate;
    ParticleTypePtr type_res = states_[0];
    ParticleTypePtr nuc = IsoParticleType::find("N").get_states()[0];
    XS_NR_tabulation = make_unique<Tabulation>(
          type_res->minimum_mass() + nuc->mass(), 2.f, 100,
          [&](float sqrts) {
            return integrate(type_res->minimum_mass(), sqrts - nuc->mass(),
                             [&](float m) {
                               return spec_func_integrand_1res(m, sqrts,
                                                        nuc->mass(), *type_res);
                             });
          });
  }
  return XS_NR_tabulation->get_value_linear(srts);
}


double IsoParticleType::get_integral_DR(double srts) {
  if (XS_DR_tabulation == nullptr) {
    // initialize tabulation
    /* TODO(weil): Move this lazy init to a global initialization function,
      * in order to avoid race conditions in multi-threading. */
    Integrator2d integrate(1E4);
    ParticleTypePtr type_res = states_[0];
    ParticleTypePtr Delta = IsoParticleType::find("Δ").get_states()[0];
    XS_DR_tabulation = make_unique<Tabulation>(
          type_res->minimum_mass() + Delta->minimum_mass(), 2.5f, 100,
          [&](float sqrts) {
            return integrate(type_res->minimum_mass(),
                             sqrts - Delta->minimum_mass(),
                             Delta->minimum_mass(),
                             sqrts - type_res->minimum_mass(),
                             [&](float m1, float m2) {
                               return spec_func_integrand_2res(sqrts, m1, m2,
                                                            *type_res, *Delta);
                             });
          });
  }
  return XS_DR_tabulation->get_value_linear(srts);
}


}  // namespace Smash
