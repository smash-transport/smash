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
      const PdgCode &pdg_a = incoming_particles_[0].type().pdgcode();
      const PdgCode &pdg_b = incoming_particles_[1].type().pdgcode();
      const float &m1 = incoming_particles_[0].effectiveMass();
      const double s = mandelstam_s();
      
      
      std::array<double, 2> mandelstam_t = get_t_range(std::sqrt(s),m1,m2,m3,m4);
      
      ParticleTypePtr part_out, photon_out;
      // do a check according to incoming_particles_ and calculate the cross sections (xsection) for all possible reactions
      process_list.push_back(make_unique<CollisionBranch>(*part_out, *photon_out, xsection,ProcessType::TwoToTwo));
      
      
      
      // create extra CollisionBranch only for photon producing reactions!
      add_processes<CollisionBranch>(std::move(process_list), collision_channels_photons_,cross_section_photons_);
      
      return process_list;
    }
    
}  // namespace Smash