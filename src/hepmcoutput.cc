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

  // TODO Is MM the correct unit?
  current_event_ = HepMC3::GenEvent(HepMC3::Units::GEV,HepMC3::Units::MM);

  /* Rivet needs a value for the cross section, but SMASH only knows about
   * nucleon cross sections, not about cross sections between nuclei,
   * so a dummy is used. */
  std::shared_ptr<HepMC3::GenCrossSection> cross_section = std::make_shared<HepMC3::GenCrossSection>();
  current_event_.add_attribute("GenCrossSection",cross_section);
  const double dummy_xs = 1.0;
  cross_section->set_cross_section(1.0, 1.0);

  current_event_.set_event_number(event_number);
  vertex_ =  HepMC3::make_shared<HepMC3::GenVertex>();
  current_event_.add_vertex(vertex_);


  // Add projectile and target

  // TODO Properly create two nuclei "particles" (do not hardcode information about nuclei!)

  const int na = 197;
  const int nz = 97;
  const int nl = 0;

  // Construct pdg nuclear code
  const int pdg_nuclear_code_prefix = 10 * 1E8;
  const int pdg_nuclear_code_lambda = nl * 1E7;
  const int pdg_nuclear_code_charge = nz * 1E4;
  const int pdg_nuclear_code_baryon = na * 1E1;
  const int pdg_nuclear_code_isomer =  0 * 1E0;  // we do not do isomers in SMASH

  const int pdg_nuclear_code = pdg_nuclear_code_prefix
                             + pdg_nuclear_code_lambda
                             + pdg_nuclear_code_charge
                             + pdg_nuclear_code_baryon
                             + pdg_nuclear_code_excite;

  // TODO Put nuclear code construction in proper function and into class (maybe pdgcode class, static function?)

  FourVector total_mom_proj;
  FourVector total_mom_target;
  int counter = 0;
  for (const ParticleData &data : particles) {
    if (data.id() < na) {
      total_mom_proj += data.momentum();
    } else {
      total_mom_target += data.momentum();
    }
    counter++;
  }

  if (counter != 2*na) {
    throw std::invalid_argument("HepMC currently only supported for AuAu nuclei ATM");
  }

  // TODO Does make_shared really makes sense here?
  HepMC3::GenParticlePtr projectile_p =  HepMC3::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(total_mom_proj.x1(), total_mom_proj.x2(), total_mom_proj.x3(), total_mom_proj.x0()), pdg_nuclear_code_au, status_code_for_beam_particles);
  vertex_->add_particle_in(projectile_p);
  HepMC3::GenParticlePtr target_p =  HepMC3::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(total_mom_target.x1(), total_mom_target.x2(), total_mom_target.x3(), total_mom_target.x0()), pdg_nuclear_code_au, status_code_for_beam_particles);
  vertex_->add_particle_in(target_p);

}


void HepMcOutput::at_eventend(const Particles &particles, const int32_t /*event_number*/,
                 double impact_parameter, bool /*empty_event*/) {

   // Set heavy ion attribute, only the impact parameter is known
   // TODO Test if it works to add the heavy ion attribute here
   std::shared_ptr<HepMC3::GenHeavyIon> heavy_ion = std::make_shared<HepMC3::GenHeavyIon>();
   current_event_.add_attribute("GenHeavyIon", heavy_ion);
   heavy_ion->set(-1, -1, -1, -1, -1, -1, -1, -1, -1, impact_parameter, -1.0, -1.0, -1.0, -1.0, -1.0);

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
