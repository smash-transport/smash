/*
 *
 *    Copyright (c) 2021 Christian Holm Christensen
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/hepmcinterface.h"

#include "HepMC3/Print.h"
#include "HepMC3/Setup.h"
#include "smash/logging.h"

namespace smash {

HepMcInterface::HepMcInterface(const std::string& name, const bool full_event,
                               const int total_N, const int proj_N)
    : OutputInterface(name),
      event_(HepMC3::Units::GEV, HepMC3::Units::MM),
      ion_(),
      xs_(),
      ip_(),
      total_N_(total_N),
      proj_N_(proj_N),
      full_event_(full_event) {
  logg[LOutput].debug() << "Name of output: " << name << " "
                        << (full_event_ ? "full event" : "final state only")
                        << " output" << std::endl;
  ion_ = std::make_shared<HepMC3::GenHeavyIon>();
  xs_ = std::make_shared<HepMC3::GenCrossSection>();

  HepMC3::Setup::set_debug_level(logg[LOutput].isEnabled<einhard::DEBUG>() ? 5
                                                                           : 0);
}

void HepMcInterface::at_eventstart(const Particles& particles,
                                   const int event_number,
                                   const EventInfo& event) {
  // Clear event and mapping and set event number
  clear();

  // Set header stuff on event
  ion_->impact_parameter = event.impact_parameter;
  xs_->set_cross_section(1, 1);  // Dummy values
  event_.set_event_number(event_number);
  event_.set_heavy_ion(ion_);
  coll_.resize(total_N_);

  // Create IP only if final state
  ip_ = std::make_shared<HepMC3::GenVertex>();
  event_.add_vertex(ip_);

  // Count up projectile and target
  smash::FourVector p_proj;
  smash::FourVector p_targ;
  AZ az_proj;
  AZ az_targ;
  bool is_coll = proj_N_ > 0 && total_N_ > 0;

  for (auto& data : particles) {
    if (is_coll) {
      if (!data.is_neutron() && !data.is_proton()) {
        throw std::runtime_error(
            "Particle of PID=" + std::to_string(data.pdgcode().get_decimal()) +
            " is not a valid HepMC beam particle!");
      }

      bool isproj = data.id() < proj_N_;
      smash::FourVector& p = (isproj ? p_proj : p_targ);
      AZ& az = (isproj ? az_proj : az_targ);
      p += data.momentum();
      az.second += data.is_proton();
      az.first++;
    }

    if (!full_event_ && is_coll) {
      continue;
    }

    // If we are not only adding final state particles, but all
    // interactions, or not in collider modus, then we add all the
    // incoming particles as outgoing (incoming for non-collider
    // modus) of the vertex.  In that way, we will keep track of
    // participants and spectators as well as make the event
    // structure consistent.
    auto op = make_register(data, Status::fnal);
    if (is_coll) {
      ip_->add_particle_out(op);
    } else {
      ip_->add_particle_in(op);
    }
  }
  // Make beam particles
  if (is_coll) {
    auto proj = make_gen(ion_pdg(az_proj), Status::beam, p_proj);
    auto targ = make_gen(ion_pdg(az_targ), Status::beam, p_targ);

    // Add to interaction point if we need to
    ip_->add_particle_in(proj);
    ip_->add_particle_in(targ);
  }
}

void HepMcInterface::at_interaction(const Action& action,
                                    const double /* density */) {
  HepMC3::GenVertexPtr vp;
  auto type = action.get_type();
  int status = get_status(type);

  if (full_event_) {
    FourVector v = action.get_interaction_point();
    vp = std::make_shared<HepMC3::GenVertex>(
        HepMC3::FourVector(v.x1(), v.x2(), v.x3(), v.x0()));
    vp->add_attribute("weight", std::make_shared<HepMC3::FloatAttribute>(
                                    action.get_total_weight()));
    vp->add_attribute(
        "partial_weight",
        std::make_shared<HepMC3::FloatAttribute>(action.get_partial_weight()));
    event_.add_vertex(vp);
  }

  bool has_proj = false;
  bool has_targ = false;
  // Check if we have initial particles on _both_ sides
  for (auto& i : action.incoming_particles()) {
    if (type != ProcessType::Elastic) {
      if (i.id() < proj_N_) {
        has_proj = true;
      } else if (i.id() < total_N_) {
        has_targ = true;
      }
    }
  }
  // Check if this a collision and count it
  bool is_coll = has_proj && has_targ;
  if (is_coll) {
    ncoll_++;
    if (type == ProcessType::StringHard) {
      ncoll_hard_++;
    }
  }

  // Now mark participants
  for (auto& i : action.incoming_particles()) {
    if (is_coll) {
      coll_[i.id()]++;
    }

    if (!full_event_)
      continue;

    // Create tree
    HepMC3::GenParticlePtr ip = find_or_make(i, status);
    ip->set_status(status);
    vp->add_particle_in(ip);
  }
  if (!full_event_)
    return;

  // Add outgoing particles
  for (auto& o : action.outgoing_particles()) {
    vp->add_particle_out(make_register(o, Status::fnal));
  }
}

void HepMcInterface::at_eventend(const Particles& particles,
                                 const int32_t /*event_number*/,
                                 const EventInfo& event) {
  // In case this was an empty event
  if (event.empty_event) {
    clear();
    return;
  }
  /* since the number of collisions is not an experimental obervable, we set it
   * to -1.
   */
  ion_->Ncoll_hard = -1;
  ion_->Ncoll = -1;
  /* This should be the number of participants in the projectile nucleus.
   * However, to avoid confusion with the Glauber model, we prefer to set it -1.
   */
  ion_->Npart_proj = -1;
  /* his should be the number of participants in the target nucleus.
   * However, to avoid confusion with the Glauber model, we prefer to set it -1.
   */
  ion_->Npart_targ = -1;
  // If we only do final state events, then take particle
  // Note, we should already have the particles if not only final
  // state
  // Take all passed particles and add as outgoing particles to event
  for (auto& p : particles) {
    if (!full_event_) {
      auto h = make_register(p, Status::fnal);
      ip_->add_particle_out(h);
    } else if (map_.find(p.id()) == map_.end()) {
      throw std::runtime_error("Dangling particle " + std::to_string(p.id()));
    }
  }
}

void HepMcInterface::clear() {
  event_.clear();
  map_.clear();
  ip_ = 0;
  coll_ = 0;
  ncoll_ = 0;
  ncoll_hard_ = 0;
}

int HepMcInterface::get_status(const ProcessType& t) const {
  if (t == ProcessType::Decay)
    return Status::dcy;

  // Make all other codes SMASH specific
  return static_cast<int>(t) + static_cast<int>(Status::off);
}

HepMC3::GenParticlePtr HepMcInterface::make_gen(int pid, int status,
                                                const smash::FourVector& mom,
                                                double mass) {
  auto p = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(mom.x1(), mom.x2(), mom.x3(), mom.x0()), pid, status);
  if (mass > 0)
    p->set_generated_mass(mass);
  return p;
}

HepMC3::GenParticlePtr HepMcInterface::make_register(const ParticleData& p,
                                                     int status) {
  int id = p.id();
  auto h = make_gen(p.pdgcode().get_decimal(), status, p.momentum(),
                    p.type().mass());
  map_[id] = h;

  return h;
}

HepMC3::GenParticlePtr HepMcInterface::find_or_make(const ParticleData& p,
                                                    int status,
                                                    bool force_new) {
  int id = p.id();
  if (!force_new) {
    auto it = map_.find(id);
    if (it != map_.end()) {
      return it->second;
    }
  }

  return make_register(p, status);
}

int HepMcInterface::ion_pdg(const AZ& az) const {
  if (az.second == 1 && az.first == 1)
    return 2212;  // Proton
  if (az.second == 1 && az.first == 0)
    return 2112;  // Neutron

  return 1'000'000'000 + az.second * 10'000 + az.first * 10;
}
}  // namespace smash
