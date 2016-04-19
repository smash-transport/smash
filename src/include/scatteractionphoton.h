#ifndef SRC_INCLUDE_SCATTERACTIONPHOTON_H_
#define SRC_INCLUDE_SCATTERACTIONPHOTON_H_

#include "scatteraction.h"

namespace Smash {

class ScatterActionPhoton : public ScatterAction {
  public:
    ScatterActionPhoton(const ParticleData &in_part1, const ParticleData &in_part2, float time_of_execution):
       ScatterAction(in_part1, in_part2, time_of_execution){}
    
    void generate_final_state() override; 
    float raw_weight_value() const override { return weight_; } 

  protected:
     virtual CollisionBranchList two_to_two_cross_sections() override;
     
  private:
    float weight_;
    /** List of possible collisions producing photons */
    // All photon reactions are added to a separarte CollisionList, which is later used to chose the scatteraction taking place
    // The total cross section is still the sum of both lists
    CollisionBranchList collision_channels_photons_;
    float cross_section_photons_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONPHOTON_H_

