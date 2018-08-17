/*
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/isoparticletype.h"

#include "smash/integrate.h"
#include "smash/kinematics.h"
#include "smash/logging.h"

namespace smash {

static IsoParticleTypeList iso_type_list;

IsoParticleType::IsoParticleType(const std::string &n, double m, double w,
                                 unsigned int s)
    : name_(n), mass_(m), width_(w), spin_(s) {}

const IsoParticleTypeList &IsoParticleType::list_all() { return iso_type_list; }

/// Helper function for IsoParticleType::try_find and friends.
static IsoParticleType *try_find_private(const std::string &name) {
  auto found =
      std::lower_bound(iso_type_list.begin(), iso_type_list.end(), name,
                       [](const IsoParticleType &l, const std::string &r) {
                         return l.name() < r;
                       });
  if (found == iso_type_list.end() || found->name() != name) {
    return {};  // The default constructor creates an invalid pointer.
  }
  return &*found;
}

const IsoParticleType *IsoParticleType::try_find(const std::string &name) {
  return try_find_private(name);
}

const IsoParticleType &IsoParticleType::find(const std::string &name) {
  const auto found = try_find_private(name);
  if (!found) {
    throw ParticleNotFoundFailure("Isospin multiplet " + name + " not found!");
  }
  return *found;
}

IsoParticleType &IsoParticleType::find_private(const std::string &name) {
  auto found = try_find_private(name);
  if (!found) {
    throw ParticleNotFoundFailure("Isospin multiplet " + name +
                                  " not found (privately)!");
  }
  return *found;
}

bool IsoParticleType::exists(const std::string &name) {
  const auto found = try_find_private(name);
  return found;
}

/**
 * Construct the name-string for an isospin multiplet from the given
 * name-string for the particle.
 *
 * \param[in] name name-string of the particle
 * \return the name-string for an isospin multiplet
 */
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
  auto found = std::find_if(multiplet.states_.begin(), multiplet.states_.end(),
                            [&n](ParticleTypePtr p) { return p->name() == n; });
  if (found == multiplet.states_.end()) {
    throw std::runtime_error("Isospin state " + n + " not found!");
  }
  return *found;
}

IsoParticleType *IsoParticleType::find(const ParticleType &type) {
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
                << " [ m = " << type.mass() << ", Γ = " << type.width_at_pole()
                << " ]";
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

static /*thread_local (see #3075)*/ Integrator integrate;

double IsoParticleType::get_integral_NR(double sqrts) {
  if (XS_NR_tabulation_ == nullptr) {
    // initialize tabulation
    /* TODO(weil): Move this lazy init to a global initialization function,
     * in order to avoid race conditions in multi-threading. */
    ParticleTypePtr type_res = states_[0];
    ParticleTypePtr nuc = IsoParticleType::find("N").get_states()[0];
    XS_NR_tabulation_ =
        spectral_integral_semistable(integrate, *type_res, *nuc, 2.0);
  }
  return XS_NR_tabulation_->get_value_linear(sqrts);
}

double IsoParticleType::get_integral_piR(double sqrts) {
  if (XS_piR_tabulation_ == nullptr) {
    // initialize tabulation
    /* TODO(weil): Move this lazy init to a global initialization function,
     * in order to avoid race conditions in multi-threading. */
    ParticleTypePtr type_res = states_[0];
    ParticleTypePtr pion = IsoParticleType::find("π").get_states()[0];
    XS_piR_tabulation_ =
        spectral_integral_semistable(integrate, *type_res, *pion, 2.0);
  }
  return XS_piR_tabulation_->get_value_linear(sqrts);
}

double IsoParticleType::get_integral_RK(double sqrts) {
  if (XS_RK_tabulation_ == nullptr) {
    // initialize tabulation
    /* TODO(weil): Move this lazy init to a global initialization function,
     * in order to avoid race conditions in multi-threading. */
    ParticleTypePtr type_res = states_[0];
    ParticleTypePtr kaon = IsoParticleType::find("K").get_states()[0];
    XS_RK_tabulation_ =
        spectral_integral_semistable(integrate, *type_res, *kaon, 2.0);
  }
  return XS_RK_tabulation_->get_value_linear(sqrts);
}

static /*thread_local (see #3075)*/ Integrator2dCuhre integrate2d;

double IsoParticleType::get_integral_RR(const ParticleType &type_res_2,
                                        double sqrts) {
  auto search = XS_RR_tabulations.find(find(type_res_2));
  if (search != XS_RR_tabulations.end()) {
    return search->second->get_value_linear(sqrts);
  }
  IsoParticleType *key = find(type_res_2);
  XS_RR_tabulations.emplace(key,
                            integrate_RR(find(type_res_2)->get_states()[0]));
  return XS_RR_tabulations.at(key)->get_value_linear(sqrts);
}

TabulationPtr IsoParticleType::integrate_RR(ParticleTypePtr &res2) {
  ParticleTypePtr res1 = states_[0];
  const double m1_min = res1->min_mass_kinematic();
  const double m2_min = res2->min_mass_kinematic();
  return make_unique<Tabulation>(m1_min + m2_min, 3., 125, [&](double srts) {
    const double m1_max = srts - m2_min;
    const double m2_max = srts - m1_min;
    const auto result =
        integrate2d(m1_min, m1_max, m2_min, m2_max, [&](double m1, double m2) {
          return spec_func_integrand_2res(srts, m1, m2, *res1, *res2);
        });
    return result.value();
  });
}

}  // namespace smash
