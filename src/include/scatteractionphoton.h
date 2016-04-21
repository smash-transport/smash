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

    // fill with values:
	const float m_rho = ParticleType::find(0x113).mass();
	const float m_rho_2 = std::sqr(m_rho);
	const float m_pi = ParticleType::find(0x111).mass();
	const float m_pi_2 = std::sqr(m_pi);
	const float m_eta =ParticleType::find(0x221).mass(); 
	const float m_eta_2 = std::sqr(m_eta);
    	const float gamma_rho_tot = 0.1491;
    	const float g_rho_2 = 48*acos(0)*gamma_rho_tot*std::pow(m_rho,2)/std::pow(std::pow(m_rho,2)-4*std::pow(m_pi,2),3/2);
   	const float DM = std::pow(m_rho,2)-4*std::pow(m_pi,2);
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONPHOTON_H_

