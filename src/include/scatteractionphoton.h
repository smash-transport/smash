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

#include <algorithm>

#include "constants.h"
#include "scatteraction.h"

namespace Smash {

class ScatterActionPhoton : public ScatterAction {
 public:
  ScatterActionPhoton(const ParticleList &in, float time, int nofp)
      : ScatterAction(in[0], in[1], time),
        number_of_fractional_photons(nofp) {}

  void generate_final_state() override;

  float raw_weight_value() const override { return weight_; }

  float cross_section() const override {
    if (cross_section_photons_ < really_small) {
      return cross_section_photons_;
    } else {
      return total_cross_section_ + cross_section_photons_;
    }
  }

  /** Overridden to effectively return the reaction channel. */
  virtual ProcessType get_type() const override {
    return static_cast<ProcessType>(reac);
  }

  /** Adds one dummy channel with a given cross-section. The intended use is to
   * add the hadronic cross-section from already performed hadronic action
   * without recomputing it. The photon action is never performed, so
   * this channel itself will never play any role. Only its cross-section will.
   */
  void add_dummy_hadronic_channels(float reaction_cross_section);

  /** To add only one reaction for testing purposes */
  void add_single_channel() {
    add_processes<CollisionBranch>(photon_cross_sections(),
                                  collision_channels_photons_,
                                  cross_section_photons_);
  }

  /// Tells if the given incoming particles may produce photon
  static bool is_photon_reaction(const ParticleList &in) {
    if (in.size() != 2) {
      return false;
    }
    // Turn all pi- and rho- into pi+ and rho+ to avoid listing
    // too many variants further, swap so that pion should be first.
    PdgCode a = (in[0].type().charge() >= 0) ? in[0].pdgcode() :
                 in[0].pdgcode().get_antiparticle();
    PdgCode b = (in[1].type().charge() >= 0) ? in[1].pdgcode() :
                 in[1].pdgcode().get_antiparticle();
    if (!a.is_pion()) {
      std::swap(a, b);
    }

    switch (pack(a.code(), b.code())) {
      case(pack(pdg::pi_p, pdg::pi_z)):
      case(pack(pdg::pi_z, pdg::pi_p)):
        // ReactionType::pi0_pi
      case(pack(pdg::pi_p, pdg::rho_z)):
        // ReactionType::piplus_rho0
      case(pack(pdg::pi_p, pdg::rho_p)):
        // ReactionType::pi_rho
      case(pack(pdg::pi_z, pdg::rho_p)):
        // ReactionType::pi0_rho
      case(pack(pdg::pi_p, pdg::eta)):
        // ReactionType::piplus_eta
      case(pack(pdg::pi_p, pdg::pi_p)):
        // ReactionType::pi_pi
        return true;
        break;
      default:
        return false;
    }
  }

 private:
  CollisionBranchList photon_cross_sections();

  int const number_of_fractional_photons;

  float weight_ = 0.0;

  /** List of possible collisions producing photons */
  CollisionBranchList collision_channels_photons_;

  float cross_section_photons_ = 0.0;

  const int num_tab_pts = 200;

  enum class ReactionType {
    no_reaction = 0,
    pi0_pi = 1,
    piplus_rho0 = 2,
    pi_rho = 3,
    pi0_rho = 4,
    piplus_eta = 5,
    pi_pi = 6
  };

  ReactionType reac = ReactionType::no_reaction;

  float pi_pi_rho0(const float M, const float s) const;

  float pi_pi0_rho(const float M, const float s) const;

  float diff_cross_section(float t, float m3) const;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONPHOTON_H_
