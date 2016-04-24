#include "include/scatteractionphoton.h"

#include "include/kinematics.h"
#include "include/cxx14compat.h"
#include "include/particletype.h"
#include "include/pdgcode.h"

using std::sqrt;
using std::pow;
using std::atan;

namespace Smash {
    
    void ScatterActionPhoton::generate_final_state() {
        /* Decide for a particular final state. */
        // modify so only photon producing actions can take place (still virtually of course) and weigh them correctly
        const CollisionBranch* proc = choose_channel <CollisionBranch>(collision_channels_photons_, cross_section_photons_);
        process_type_ = proc->get_type();
        outgoing_particles_ = proc->particle_list();
	weight_ = proc->weight();

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

    
    CollisionBranchList ScatterActionPhoton::two_to_two_cross_sections() {
      
      CollisionBranchList process_list;
      const float m_rho = ParticleType::find(0x113).mass();
      const float m_rho_2 = pow(m_rho,2);
      const float m_pi = ParticleType::find(0x111).mass();
      const float m_pi_2 = pow(m_pi,2);
      const float m_eta =ParticleType::find(0x221).mass(); 
      const float m_eta_2 = pow(m_eta,2);
      const float gamma_rho_tot = ParticleType::find(0x113).width_at_pole();
      const float g_rho_2 = 12*twopi*gamma_rho_tot*pow(m_rho,2)/pow(pow(m_rho,2)-4*pow(m_pi,2),3/2);
      const float DM = pow(m_rho,2)-4*pow(m_pi,2);

      ParticleData &part_a = incoming_particles_[0];
      ParticleData &part_b = incoming_particles_[1];
    
      bool no_pion = false;
    
      if (!part_a.type().pdgcode().is_pion()) {
	if (part_b.type().pdgcode().is_pion()) {
	    ParticleData dummy = part_a;
	    part_a = part_b;
	    part_b = dummy;
	} else no_pion = true;
      } 
      
      if (no_pion) {
	
      } else { // do a check according to incoming_particles_ and calculate the cross sections (xsection) for all possible reactions

	const double s = mandelstam_s();
	const double &m1 = part_a.pole_mass();
	const double &m2 = part_b.pole_mass();
	double m3 = 0.0; // will be fixed according to reaction outcome
	const double p_cm_2 = cm_momentum_squared();
	ParticleTypePtr part_out = &ParticleType::find(0x022);  
	ParticleTypePtr photon_out = &ParticleType::find(0x022);
	
	enum ReactionType {pi_pi, pi0_pi, piplus_rho0, pi_rho, pi0_rho, piplus_eta, no_reaction};
	ReactionType reac = no_reaction;
	if (part_a.type().charge()==0){
	  if (part_b.type().charge()!=0) {
	    if (part_b.type().pdgcode().is_pion()){
	      reac = pi0_pi;
	      part_out = &ParticleType::find(0x213); 
	    }
	    if (part_b.type().pdgcode().is_rho()){
	      reac = pi0_rho;
	      part_out = &ParticleType::find(0x211); 
	    }
	  }
	} else { //so part_a.type().charge()!=0
	  if (part_b.type().charge()==0){
	    if (part_b.type().pdgcode().is_pion()){
	      reac = pi0_pi;
	      part_out = &ParticleType::find(0x213);
	    }
	    if (part_b.type().pdgcode().is_rho()){
	      reac = piplus_rho0;
	      part_out = &ParticleType::find(0x211);
	    }
	    if (part_b.type().pdgcode()==0x221) { //corresponds to eta meson
	      reac = piplus_eta;
	      part_out = &ParticleType::find(0x211);
	    }
	  } else if (part_b.type().charge()==-part_a.type().charge()){
	    if (part_b.type().pdgcode().is_pion()){
	      reac = pi_pi; // actually three reactions (afh), start with h
	    }
	    if (part_b.type().pdgcode().is_rho()){
	      reac = pi_rho;
	      part_out = &ParticleType::find(0x111);
	    }
	  }
	}	
	
	m3 = part_out->mass();
	if ((sqrt_s()<=m3)||(sqrt_s()<=m1+m2)) reac = no_reaction;
	if (reac!=no_reaction){
	 std::array<double, 2> mandelstam_t = get_t_range(sqrt_s(),m1,m2,m3,0.0);
	 double t1 = mandelstam_t[1];
	 double t2 = mandelstam_t[0];

	//if (t1*t2<0) reac=no_reaction; //why does this case occur????
	
	 double u1 = pow(m1,2)+pow(m2,2)+pow(m3,2)-s-t1;
	 double u2 = pow(m1,2)+pow(m2,2)+pow(m3,2)-s-t2;

	 double e,I0,I1;
	 float xsection=0.0;
	
	 switch (reac){
	  case pi_pi: 
	    xsection = twopi*pow(alpha,2)/(s*p_cm_2);
            t1 += -m_pi_2; //is t+
            t2 += -m_pi_2;
            u1 = - s - t1;
            u2 = - s - t2;
            e = t2-t1+2*m_pi_2*((1-2*m_pi_2/s)*std::log(t2/t1)+m_pi_2*(t2-t1)/(t1*t2));
            e+= 2*m_pi_2*((1-2*m_pi_2/s)*std::log(u1/u2)+m_pi_2*(u1-u2)/(u1*u2));
            xsection = xsection*e;
	    process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
	    
	     //now the second possible reaction (produces eta)
	    part_out = &ParticleType::find(0x221);
	    m3 = part_out->mass();
	    if (sqrt_s()>m3){
	      mandelstam_t = get_t_range(sqrt_s(),m1,m2,m3,0.0);
	      t1 = mandelstam_t[1];
	      t2 = mandelstam_t[0];
	      u1 = pow(m1,2)+pow(m2,2)+pow(m3,2)-s-t1;
	      u2 = pow(m1,2)+pow(m2,2)+pow(m3,2)-s-t2;
	      xsection = M_PI*alpha*4.7*pow(m_rho,4)/(pow(s-m_rho_2,2)+pow(gamma_rho_tot,2)*m_rho_2)/(16*m_eta_2*pow(m_rho,4)*s*p_cm_2);
	      xsection = xsection*((2*m_pi_2+m_eta_2-s)*s/2*(pow(t2,2)-pow(t1,2))-s/3*(pow(t2,3)-pow(t1,3))-(t2-t1)*m_pi_2*(pow(m_eta,4)+s*(m_pi_2-m_eta_2)));
	      process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
	    
	       // and the third possible reaction (produces rho0)
	      part_out = &ParticleType::find(0x113);
	      m3 = part_out->mass();
	      if (sqrt_s()>m3){
		mandelstam_t = get_t_range(sqrt_s(),m1,m2,m3,0.0);
		t1 = mandelstam_t[1];
		t2 = mandelstam_t[0]; 
		xsection = alpha*g_rho_2/(4*s*p_cm_2);
		t1 += -m_pi_2;
		t2 += -m_pi_2;
		u1 += -m_pi_2;
		u2 += -m_pi_2;
		xsection=xsection*(2*(t2-t1)-DM*((s-2*m_pi_2)/(s-m_rho_2)*std::log(t2/t1)+m_pi_2*(t2-t1)/(t2*t1)+(s-2*m_pi_2)/(s-m_rho_2)*std::log(u1/u2)+m_pi_2*(u1-u2)/(u1*u2)));
		process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
	      }
	    }
	  break;
	  case pi0_pi:
	    xsection = -alpha*g_rho_2/(16*s*p_cm_2);
            e = 1.0/3.0 * (s-2*m_rho_2)/m_rho_2/pow(s-m_rho_2,2)*(pow(t2,3)-pow(t1,3))+1.0/3.0 * (s-2*m_rho_2)/m_rho_2/pow(s-m_rho_2,2)*(pow(u1,3)-pow(u2,3));
            e += 0.5*(s-6*m_rho_2)/m_rho_2/(s-m_rho_2)*(pow(t2,2)-pow(t1,2))+0.5*(s-6*m_rho_2)/m_rho_2/(s-m_rho_2)*(pow(u1,2)-pow(u2,2));
            t1+= -m_pi_2;
            t2+= -m_pi_2;
            u1+= -m_pi_2;
            u2+= -m_pi_2;
            e += (4*s*DM/pow(s-m_rho_2,2)+m_pi_2/m_rho_2-4.5)*(t2-t1+u1-u2);
            e +=  4*s*DM/(s-m_rho_2)*std::log(t2/t1*u1/u2);
            e += 4*m_pi_2*DM*((t2-t1)/(t2*t1)+(u1-u2)/(u1*u2));
            xsection = xsection*e; 
	    process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
	  break;
	  case piplus_rho0:
	    xsection = alpha*g_rho_2/(12*s*p_cm_2);
            t1 += -m_pi_2;
            t2 += -m_pi_2;
            xsection=xsection*(2*(t2-t1)-s*(DM)/pow(s-m_pi_2,2)*(t2-t1)-DM*((s-m_rho_2+m_pi_2)/(s-m_pi_2)*std::log(t2/t1)+m_pi_2*(t2-t1)/(t1*t2)));
	    process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
	  break;
	  case pi_rho:
	    xsection = -alpha*g_rho_2/(48*s*p_cm_2);
            t1 += -m_pi_2;
            t2 += -m_pi_2;
            u1 += -m_rho_2;
            u2 += -m_rho_2;
            e = 4*DM*(m_rho_2*(t2-t1)/(u1*u2)+m_pi_2*(t2-t1)/(t2*t1)+std::log(u1/u2*t2/t1)-m_rho_2/(s-m_pi_2)*std::log(t2/t1*u1/u2));
            e+= (s-m_pi_2)*(3.0+(s-m_pi_2)/m_rho_2)*std::log(u1/u2);
            e+= (t2-t1)*(s/m_rho_2-0.5-pow(s-m_pi_2,2)/(u1*u2));
            xsection = xsection*e;
	    process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
	  break;
	  case pi0_rho:
	    xsection = alpha*g_rho_2/(48*s*p_cm_2);
            u1 += -m_rho_2; //is u+
            u2 += -m_rho_2;
            e = (t2-t1)*(4.5-s/m_rho_2-4*s*DM/pow(s-m_pi_2,2)+(pow(s-m_pi_2,2)-4*m_rho_2*DM)/(u1*u2));
            e+= std::log(u1/u2)*(5*(s-m_pi_2)-pow(s-m_pi_2,2)/m_rho_2-4*DM*(s-m_pi_2+m_rho_2)/(s-m_pi_2));
            xsection = xsection*e; 
	    process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
	  break;
	  case piplus_eta:
	    xsection = M_PI*alpha*4.7/(16*m_eta_2*s*p_cm_2); 
            I0 = 1/(m_rho*gamma_rho_tot)*(atan((u1-m_rho_2)/(m_rho*gamma_rho_tot))-atan((u2-m_rho_2)/(m_rho*gamma_rho_tot)));
            I1 = std::log((pow(u2-m_rho_2,2)+m_rho_2*pow(gamma_rho_tot,2))/(pow(u1-m_rho_2,2)+m_rho_2*pow(gamma_rho_tot,2)));
            e = -m_pi_2*((t2+u2)*(s-m_pi_2)+pow(2*m_pi_2-s,2))*I0;
            e+= ((s-m_pi_2)*(m_pi_2+t2+u2)-2*m_pi_2*(s-2*m_pi_2))*((t2+u2-m_rho_2)*I0+0.5*I1);
            e+= -s*(t2-t1+(t2+u2-m_rho_2)*I1+pow(t2+u2-m_rho_2,2)*I0-m_rho_2*pow(gamma_rho_tot,2)*I0);
            xsection = xsection*e; 
	    process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
	  break;
	  case no_reaction:
	    //never reached
	  break;  
	 }
	}
	// add to extra CollisionBranch only for photon producing reactions!
	add_processes<CollisionBranch>(std::move(process_list), collision_channels_photons_,cross_section_photons_);
      
      }
      return process_list;
    }
    
}  // namespace Smash