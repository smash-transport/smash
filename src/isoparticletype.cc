/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/isoparticletype.h"

#include "include/integrate.h"
#include "include/kinematics.h"
#include "include/logging.h"

namespace Smash {

static IsoParticleTypeList iso_type_list;

IsoParticleType::IsoParticleType(const std::string &n, float m, float w,
                                 unsigned int s)
                                : name_(n), mass_(m), width_(w), spin_(s) {}

const IsoParticleTypeList& IsoParticleType::list_all() { return iso_type_list; }

/// Helper function for IsoParticleType::try_find and friends.
static IsoParticleType* try_find_private(const std::string &name) {
  auto found = std::lower_bound(
      iso_type_list.begin(), iso_type_list.end(), name,
      [](const IsoParticleType &l, const std::string &r) {
        return l.name() < r;
      });
  if (found == iso_type_list.end() || found->name() != name) {
    return {};  // The default constructor creates an invalid pointer.
  }
  return &*found;
}

const IsoParticleType* IsoParticleType::try_find(const std::string &name) {
  return try_find_private(name);
}

const IsoParticleType& IsoParticleType::find(const std::string &name) {
  const auto found = try_find_private(name);
  if (!found) {
    throw ParticleNotFoundFailure("Isospin multiplet " + name + " not found!");
  }
  return *found;
}

IsoParticleType& IsoParticleType::find_private(const std::string &name) {
  auto found = try_find_private(name);
  if (!found) {
    throw ParticleNotFoundFailure("Isospin multiplet " + name
                                + " not found (privately)!");
  }
  return *found;
}

bool IsoParticleType::exists(const std::string &name) {
  const auto found = try_find_private(name);
  return found;
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

const ParticleTypePtr IsoParticleType::find_state(const std::string &n) {
  const IsoParticleType &multiplet = IsoParticleType::find(multiplet_name(n));
  auto found = std::find_if(
    multiplet.states_.begin(), multiplet.states_.end(),
    [&n](ParticleTypePtr p) {
      return p->name() == n;
    });
  if (found == multiplet.states_.end()) {
    throw std::runtime_error("Isospin state " + n + " not found!");
  }
  return *found;
}

IsoParticleType* IsoParticleType::find(const ParticleType &type) {
  std::string multiname = multiplet_name(type.name());
  IsoParticleType &multiplet = find_private(multiname);
  return &multiplet;
}

void IsoParticleType::add_state(const ParticleType &type) {
  states_.push_back(&type);

  // check if isospin symmetry is fulfilled
  const auto &log = logger<LogArea::ParticleType>();
  if (std::abs(mass() - type.mass()) > really_small) {
    log.warn() << "Isospin symmetry is broken by mass of " << type.name()
               << ": " << type.mass() << " vs. " << mass();
  }
  if (std::abs(width() - type.width_at_pole()) > really_small) {
    log.warn() << "Isospin symmetry is broken by width of " << type.name()
               << ": " << type.width_at_pole() << " vs. " << width();
  }
  if (spin() != type.spin()) {
    log.error() << "Isospin symmetry is broken by spin of " << type.name()
                << ": " << type.spin() << " vs. " << spin();
  }
}

void IsoParticleType::create_multiplet(const ParticleType &type) {
  const auto &log = logger<LogArea::ParticleType>();

  // create multiplet if it does not exist yet
  std::string multiname = multiplet_name(type.name());
  if (!exists(multiname)) {
    iso_type_list.emplace_back(multiname, type.mass(), type.width_at_pole(),
                               type.spin());
    log.debug() << "Creating isospin multiplet " << multiname
                << " [ m = " << type.mass()
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



static thread_local Integrator integrate;

double IsoParticleType::get_integral_NR(double sqrts) {
  if (XS_NR_tabulation_ == nullptr) {
    // initialize tabulation
    /* TODO(weil): Move this lazy init to a global initialization function,
     * in order to avoid race conditions in multi-threading. */
    ParticleTypePtr type_res = states_[0];
    ParticleTypePtr nuc = IsoParticleType::find("N").get_states()[0];
    XS_NR_tabulation_ = spectral_integral_semistable(integrate,
                                                     *type_res, *nuc, 2.0);
  }
  return XS_NR_tabulation_->get_value_linear(sqrts);
}

double IsoParticleType::get_integral_RK(double sqrts) {
  if (XS_RK_tabulation_ == nullptr) {
    // initialize tabulation
    /* TODO(weil): Move this lazy init to a global initialization function,
     * in order to avoid race conditions in multi-threading. */
    ParticleTypePtr type_res = states_[0];
    ParticleTypePtr kaon = IsoParticleType::find("K").get_states()[0];
    XS_RK_tabulation_ = spectral_integral_semistable(integrate,
                                                     *type_res, *kaon, 2.0);
  }
  return XS_RK_tabulation_->get_value_linear(sqrts);
}

static thread_local Integrator2d integrate2d(1E4);

double IsoParticleType::get_integral_RR(const ParticleType &type_res_2,
                                        double sqrts) {
  auto search = XS_RR_tabulations.find(find(type_res_2));
  if (search != XS_RR_tabulations.end()) {
    return search->second->get_value_linear(sqrts);
  }
  IsoParticleType* key = find(type_res_2);
  XS_RR_tabulations.emplace(key,
                            integrate_RR(find(type_res_2)->get_states()[0]));
  return XS_RR_tabulations.at(key)->get_value_linear(sqrts);
}

TabulationPtr IsoParticleType::integrate_RR(ParticleTypePtr &type_res_2) {
  ParticleTypePtr type_res_1 = states_[0];
  return make_unique<Tabulation>(
         type_res_1->min_mass_kinematic() +
         type_res_2->min_mass_kinematic(),
         3.f, 125,
         [&](float srts) {
            const auto result = integrate2d(type_res_1->minimum_mass_kinematic(),
                               srts - type_res_2->minimum_mass_kinematic(),
                               type_res_2->minimum_mass_kinematic(),
                               srts - type_res_1->minimum_mass_kinematic(),
                               [&](float m1, float m2) {
                                  return spec_func_integrand_2res(srts, m1, m2,
                                                      *type_res_1, *type_res_2);
                               });
            if (result.first == 0.) {
              return 0.;
            }
            const auto relerror = std::abs(result.second / result.first);
            const auto tol = 9e-2;
            if (relerror > tol) {
              std::stringstream error_msg;
              error_msg << "Integration error larger than " << tol*100 << "%: " << result.first << " +- " << result.second;
              throw std::runtime_error(error_msg.str());
            }
            return result.first;
         });
}

}  // namespace Smash
