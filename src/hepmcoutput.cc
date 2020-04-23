/*
 *
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/hepmcoutput.h"

#include "HepMC3/GenParticle.h"
#include "HepMC3/Print.h"

namespace smash {

HepMcOutput::HepMcOutput(const bf::path &path, std::string name,
            const OutputParameters &out_par) : OutputInterface(name),
            output_file_("data/"+name) {
  // TODO create output in correct path
}


void HepMcOutput::at_eventstart(const Particles &particles,
                                const int event_number) {
  // TODO What other meta-data is needed for event or vertex?
  current_event_ = HepMC3::GenEvent(HepMC3::Units::GEV,HepMC3::Units::MM);

  // Set cross section attribute
  std::shared_ptr<HepMC3::GenCrossSection> cross_section = std::make_shared<HepMC3::GenCrossSection>();
  current_event_.add_attribute("GenCrossSection",cross_section);
  // TODO Replace dummy XS
  cross_section->set_cross_section(1.0, 1.0);

  // Set heavy ion attribute
  std::shared_ptr<HepMC3::GenHeavyIon> heavy_ion = std::make_shared<HepMC3::GenHeavyIon>();
  current_event_.add_attribute("GenHeavyIon",heavy_ion);
  heavy_ion->set(-1, -1, -1, -1, -1, -1, -1, -1, -1, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0);

  current_event_.set_event_number(event_number);
  vertex_ =  HepMC3::make_shared<HepMC3::GenVertex>();
  current_event_.add_vertex(vertex_);
  for (const ParticleData &data : particles) {
    const FourVector mom = data.momentum();
    // TODO Does make_shared really makes sense here?
    HepMC3::GenParticlePtr p =  HepMC3::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(mom.x1(), mom.x2(), mom.x3(), mom.x0()), data.pdgcode().get_decimal(), status_code_for_beam_particles);
    vertex_->add_particle_in(p);
  }
}


void HepMcOutput::at_eventend(const Particles &particles, const int32_t /*event_number*/,
                 double /*impact_parameter*/, bool /*empty_event*/) {
  for (const ParticleData &data : particles) {
    const FourVector mom = data.momentum();
    HepMC3::GenParticlePtr p =  HepMC3::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(mom.x1(), mom.x2(), mom.x3(), mom.x0()), data.pdgcode().get_decimal(), status_code_for_final_particles);
    vertex_->add_particle_out(p);
  }

 output_file_.write_event(current_event_);
}


// Just have one vertex for whole event ATM.
void HepMcOutput::at_interaction(const Action &action, const double density) {
//  // COPIED EXAMPLE FROM HEPMC DOC:
//
//  // Add vertex with incoming and outgoing particles to event
//
//   GenParticlePtr p1 = make_shared<GenParticle>( FourVector( 0.0,    0.0,   7000.0,  7000.0  ),2212,  3 );
//   GenParticlePtr p2 = make_shared<GenParticle>( FourVector( 0.750, -1.569,   32.191,  32.238),   1,  3 );
//   GenParticlePtr p3 = make_shared<GenParticle>( FourVector( 0.0,    0.0,  -7000.0,  7000.0  ),2212,  3 );
//   GenParticlePtr p4 = make_shared<GenParticle>( FourVector(-3.047,-19.0,    -54.629,  57.920),  -2,  3 );
//
//   GenVertexPtr ver = make_shared<GenVertex>();
//   // Add for loop over incoming and outgoing particles
//   ver->add_particle_in (p1);
//   ver->add_particle_in (p2);
//   ver->add_particle_out(p3);
//   ver->add_particle_out(p4);
//   evt.add_vertex(ver);
}

}  // namespace smash
