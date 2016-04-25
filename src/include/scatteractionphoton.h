#ifndef SRC_INCLUDE_SCATTERACTIONPHOTON_H_
#define SRC_INCLUDE_SCATTERACTIONPHOTON_H_

#include "scatteraction.h"
#include "constants.h"

namespace Smash {

class ScatterActionPhoton : public ScatterAction {
  public:
	using ScatterAction::ScatterAction;
   // ScatterActionPhoton(const ParticleData &in_part1, const ParticleData &in_part2, float time_of_execution):
     //  ScatterAction(in_part1, in_part2, time_of_execution){}
    
    void generate_final_state() override; 
    float raw_weight_value() const override { return weight_; } 
    float cross_section() const override {
	if (cross_section_photons_<really_small) {
		return cross_section_photons_;
	} else return total_cross_section_;
    } 
    
    CollisionBranchList two_to_two_cross_sections() override;
     
  private:
    float weight_=0.0;
    /** List of possible collisions producing photons */
    CollisionBranchList collision_channels_photons_;
    float cross_section_photons_=0.0;
    enum ReactionType {pi_pi, pi0_pi, piplus_rho0, pi_rho, pi0_rho, piplus_eta, no_reaction};
    ReactionType reac = no_reaction;
    float diff_cross_section() const;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONPHOTON_H_

