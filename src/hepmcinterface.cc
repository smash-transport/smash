/*
 *
 *    Copyright (c) 2021 Christian Holm Christensen
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/hepmcinterface.h"
#include "smash/config.h"

#include "HepMC3/GenRunInfo.h"
#include "HepMC3/Print.h"
#include "HepMC3/Setup.h"
#include "smash/logging.h"

namespace smash {

HepMcInterface::HepMcInterface(const std::string& name, const bool full_event)
    : OutputInterface(name),
      event_(HepMC3::Units::GEV, HepMC3::Units::MM),
      ion_(),
      xs_(),
      ip_(),
      full_event_(full_event) {
  logg[LOutput].debug() << "Name of output: " << name << " "
                        << (full_event_ ? "full event" : "final state only")
                        << " output" << std::endl;
  ion_ = std::make_shared<HepMC3::GenHeavyIon>();
  ion_->set(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1.0, -1.0, -1.0, -1.0, -1.0,
            -1.0);

  xs_ = std::make_shared<HepMC3::GenCrossSection>();
  // GenRunInfo: HepMC3 class to store run-related information
  std::shared_ptr<HepMC3::GenRunInfo> run_info =
      std::make_shared<HepMC3::GenRunInfo>();
  std::vector<std::string> weightnames;
  weightnames.push_back("Default");
  run_info->set_weight_names(weightnames);
  HepMC3::GenRunInfo::ToolInfo tool;
  tool.name = "SMASH";
  tool.version = SMASH_VERSION_VERBOSE;
  tool.version = tool.version + GIT_BRANCH;
  tool.description = "";
  run_info->tools().push_back(tool);
  event_.set_run_info(run_info);
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

  // Create IP only if final state
  ip_ = std::make_shared<HepMC3::GenVertex>();
  event_.add_vertex(ip_);

  // Count up projectile and target
  smash::FourVector p_proj;
  smash::FourVector p_targ;
  AZ az_proj{0, 0};
  AZ az_targ{0, 0};
  bool is_coll = (event.impact_parameter >= 0.0);

  for (auto& data : particles) {
    if (is_coll) {
      if (!data.is_neutron() && !data.is_proton()) {
        throw std::runtime_error(
            "Particle of PID=" + std::to_string(data.pdgcode().get_decimal()) +
            " is not a valid HepMC beam particle!");
      }

      if (data.belongs_to() == BelongsTo::Projectile) {
        p_proj += data.momentum();
        az_proj.first++;
        az_proj.second += data.type().charge();
      } else if (data.belongs_to() == BelongsTo::Target) {
        p_targ += data.momentum();
        az_targ.first++;
        az_targ.second += data.type().charge();
      }
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

  coll_.resize(az_proj.first + az_targ.first);
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
    event_.add_vertex(vp);
    vp->add_attribute("weight", std::make_shared<HepMC3::FloatAttribute>(
                                    action.get_total_weight()));
    vp->add_attribute(
        "partial_weight",
        std::make_shared<HepMC3::FloatAttribute>(action.get_partial_weight()));
  }

  // Now mark participants
  if (full_event_) {
    for (auto& i : action.incoming_particles()) {
      // Create tree
      HepMC3::GenParticlePtr ip = find_or_make(i, status);
      ip->set_status(status);
      vp->add_particle_in(ip);
    }
  } else {
    return;
  }

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
  // Set the weights
  event_.weights() = std::vector<double>(1, 1);
  /* since the number of collisions is not an experimental obervable, we set it
   * to -1.
   */
  ion_->Ncoll_hard = -1;
  ion_->Ncoll = -1;
  /* This should be the number of participants in the projectile nucleus.
   * However, to avoid confusion with the Glauber model, we prefer to set it -1.
   */
  ion_->Npart_proj = -1;
  /* This should be the number of participants in the target nucleus.
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

  return 1000 * 1000 * 1000 + az.second * 10 * 1000 + az.first * 10;
}
}  // namespace smash
