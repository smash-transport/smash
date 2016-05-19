/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONPHOTON_H_
#define SRC_INCLUDE_SCATTERACTIONPHOTON_H_

#include "constants.h"
#include "scatteraction.h"

namespace Smash {

class ScatterActionPhoton : public ScatterAction {
 public:
  ScatterActionPhoton(const ParticleData &in_part1,
                      const ParticleData &in_part2, float time, int nofp)
      : ScatterAction(in_part1, in_part2, time),
        number_of_fractional_photons(nofp) {}
  void generate_final_state() override;
  float raw_weight_value() const override { return weight_; }
  float cross_section() const override {
    if (cross_section_photons_ < really_small) {
      return cross_section_photons_;
    } else
      return total_cross_section_;
  }
  CollisionBranchList two_to_two_cross_sections() override;

 private:
  int const number_of_fractional_photons;
  float weight_ = 0.0;
  /** List of possible collisions producing photons */
  CollisionBranchList collision_channels_photons_;
  float cross_section_photons_ = 0.0;
  const static int num_tab_pts = 200;
  enum ReactionType {
    pi_pi,
    pi0_pi,
    piplus_rho0,
    pi_rho,
    pi0_rho,
    piplus_eta,
    no_reaction
  };
  ReactionType reac = no_reaction;
  float pi_pi_rho0(const float M, const float s) const;
  float pi_pi0_rho(const float M, const float s) const;
  float diff_cross_section(float t) const;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONPHOTON_H_
