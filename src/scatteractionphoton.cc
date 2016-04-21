#include "include/ScatterActionPhoton.h"

#include "include/constants.h"
#include "include/kinematics.h"
#include "include/cxx14compat.h"
#include "particletype.h"

using std::sqrt;
using std::log;
using std::pow;

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

	const double s = mandelstam_s();
	const double sqrt_s = std::sqrt(s);
	const float &m1 = part_a.pole_mass();
	const float &m2 = part_b.pole_mass();
	float m3;
	const float p_cm_2 = p_cm_sqr(sqrt_s,m1,m2);
	ParticleTypePtr part_out, photon_out;
	
	enum ReactionType {pi_pi, pi0_pi, piplus_rho0, pi_rho, pi0_rho, piplus_eta, no_reaction};
	ReactionType reac = no_reaction;
	if (part_a.type().charge()==0){
	  if (part_b.type().charge()!=0) {
	    if (part_b.type().pdgcode().is_pion()){
	      reac = pi0_pi;
	      m3 = m_rho;
	    }
	    if (part_b.type().pdgcode().is_rho()){
	      reac = pi0_rho;
	      m3 = m_pi;
	    }
	  }
	} else { //so part_a.type().charge()!=0
	  if (part_b.type().charge()==0){
	    if (part_b.type().pdgcode().is_pion()){
	      reac = pi0_pi;
	      m3 = m_rho;
	    }
	    if (part_b.type().pdgcode().is_rho()){
	      reac = piplus_rho0;
	      m3 = m_pi;
	    }
	    if (part_b.type().pdgcode()==0x221) { //corresponds to eta meson
	      reac = piplus_eta;
	      m3 = m_pi;
	    }
	  } else if (part_b.type().charge()==-part_a.type().charge()){
	    if (part_b.type().pdgcode().is_pion()){
	      reac = pi_pi; // actually three reactions (afh) start with a
	      m3 = m_rho; // other cases: m3 = m_eta, m3 = 0
	    }
	    if (part_b.type().pdgcode().is_rho()){
	      reac = pi_rho;
	      m3 = m_rho;
	    }
	  }
	}	
	
	 std::array<double, 2> mandelstam_t = get_t_range(sqrt_s,m1,m2,m3,0);
	 double t1 = mandelstam_t[0];
	 double t2 = mandelstam_t[1];
	 double u1 = pow(m1,2)+pow(m2,2)+pow(m3,2)-s-t1;
	 double u2 = pow(m1,2)+pow(m2,2)+pow(m3,2)-s-t2;
	 float e,I0,I1,xsection;
	
	switch (reac){
	  case pi_pi:
	    xsection = alpha*g_rho_2/(4*s*p_cm_2);
            t1 += -m_pi_2;
            t2 += -m_pi_2;
            u1 += -m_pi_2;
            u2 += -m_pi_2;
            xsection=xsection*(2*(t2-t1)-DM*((s-2*m_pi_2)/(s-m_rho_2)*log(t2/t1)+m_pi_2*(t2-t1)/(t2*t1)+(s-2*m_pi_2)/(s-m_rho_2)*log(u1/u2)+m_pi_2*(u1-u2)/(u1*u2)));
	    process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
	    
	    //now the second possible reaction 
	    m3 = m_eta;
	    mandelstam_t = get_t_range(sqrt_s,m1,m2,m3,0);
	    t1 = mandelstam_t[0];
	    t2 = mandelstam_t[1];
	    u1 = pow(m1,2)+pow(m2,2)+pow(m3,2)-s-t1;
	    u2 = pow(m1,2)+pow(m2,2)+pow(m3,2)-s-t2;
	    xsection = 2*acos(0)*alpha*4.7*pow(m_rho,4)/(pow(s-m_rho_2,2)+pow(gamma_rho_tot,2)*m_rho_2)/(16*m_eta_2*pow(m_rho,4)*s*p_cm_2);
            xsection = xsection*((2*m_pi_2+m_eta_2-s)*s/2*(pow(t2,2)-pow(t1,2))-s/3*(pow(t2,3)-pow(t1,3))-(t2-t1)*m_pi_2*(pow(m_eta,4)+s*(m_pi_2-m_eta_2)));
	    process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
	    
	    // and the third possible reaction
	    m3 = 0;
	    mandelstam_t = get_t_range(sqrt_s,m1,m2,m3,0);
	    t1 = mandelstam_t[0];
	    t2 = mandelstam_t[1];
	    u1 = pow(m1,2)+pow(m2,2)+pow(m3,2)-s-t1;
	    u2 = pow(m1,2)+pow(m2,2)+pow(m3,2)-s-t2;
	    
	  break;
	  case pi0_pi:
	      
	  break;
	  case piplus_rho0:
	    
	  break;
	  case pi_rho:
	    
	  break;
	  case pi0_rho:
	    
	  break;
	  case piplus_eta:
	    
	  break;
	  default:
	    
	  break;  
	}

	// do a check according to incoming_particles_ and calculate the cross sections (xsection) for all possible reactions
	
      
       
	// create extra CollisionBranch only for photon producing reactions!
	add_processes<CollisionBranch>(std::move(process_list), collision_channels_photons_,cross_section_photons_);
      
	return process_list;
      }
    }
    
}  // namespace Smash