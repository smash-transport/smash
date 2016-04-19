#include "include/ScatterActionPhoton.h"

#include "include/constants.h"
#include "include/kinematics.h"
#include "include/cxx14compat.h"
#include "particletype.h"

namespace Smash {
    
    void ScatterActionPhoton::generate_final_state() override {
        /* Decide for a particular final state. */
        // modify so only photon producing actions can take place (still virtually of course) and weigh them correctly
        const CollisionBranch* proc = choose_channel <CollisionBranch>(collision_channels_photons_, cross_section_photons_);
        process_type_ = proc->get_type();
        outgoing_particles_ = proc->particle_list();

        /* The production point of the new particles.  */
        FourVector middle_point = get_interaction_point();

        /* 2->2 inelastic scattering */
        /* Sample the particle momenta in CM system. */
        sample_2body_phasespace();

        /* Set positions & boost to computational frame. */
        for (ParticleData &new_particle : outgoing_particles_) {
            new_particle.set_4position(middle_point);
            new_particle.boost_momentum(-beta_cm());
        }
    }

    
    virtual CollisionBranchList two_to_two_cross_sections() override {
      CollisionBranchList process_list;
      
      ParticleData &part_a = incoming_particles_[0];
      ParticleData &part_b = incoming_particles_[1];
    
      bool no_pion = false;
    
      if !(part_a.type().pdgcode().is_pion()) {
	if (part_b.type().pdgcode().is_pion()) {
	    ParticleData dummy = part_a;
	    part_a = part_b;
	    part_b = dummy;
	} else no_pion = true;
      } 
      
      if (no_pion) {
	return process_list;
      } else { 

	
	enum ReactionType {pi_pi, pi0_pi, piplus_rho0, pi_rho, pi0_rho, piplus_eta, no_reaction};
	ReactionType reac = no_reaction;
	if (part_a.type().charge()==0){
	  if (part_b.type().charge()!=0) {
	    if (part_b.type().pdgcode().is_pion()){
	      reac = pi0_pi;
	    }
	    if (part_b.type().pdgcode().is_rho()){
	      reac = pi0_rho;
	    }
	  }
	} else { //so part_a.type().charge()!=0
	  if (part_b.type().charge()==0){
	    if (part_b.type().pdgcode().is_pion()){
	      reac = pi0_pi;
	    }
	    if (part_b.type().pdgcode().is_rho()){
	      reac = piplus_rho0;
	    }
	    if (part_b.type().pdgcode()==0x221) { //corresponds to eta meson
	      reac = piplus_eta;
	    }
	  } else if (part_b.type().charge()==-part_a.type().charge()){
	    if (part_b.type().pdgcode().is_pion()){
	      reac = pi_pi;
	    }
	    if (part_b.type().pdgcode().is_rho()){
	      reac = pi_rho;
	    }
	  }
	}	
	
	const double s = mandelstam_s();
	const float &m1 = part_a.effective_mass();
	const float &m2 = part_b.effective_mass();
	float m3;
	
	// switch case .....
	
	std::array<double, 2> mandelstam_t = get_t_range(std::sqrt(s),m1,m2,m3,m4);
	double t1 = mandelstam_t[0];
	double t2 = mandelstam_t[1];
	
	ParticleTypePtr part_out, photon_out;
	// do a check according to incoming_particles_ and calculate the cross sections (xsection) for all possible reactions
	process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
      
      
      
	// create extra CollisionBranch only for photon producing reactions!
	add_processes<CollisionBranch>(std::move(process_list), collision_channels_photons_,cross_section_photons_);
      
	return process_list;
      }
    }
    
}  // namespace Smash