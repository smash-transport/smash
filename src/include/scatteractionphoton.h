/*
 *
 *    Copyright (c) 2016-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONPHOTON_H_
#define SRC_INCLUDE_SCATTERACTIONPHOTON_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>

#include "constants.h"
#include "scatteraction.h"

#include "angles.h"
#include "constants.h"
#include "cxx14compat.h"
#include "integrate.h"
#include "kinematics.h"
#include "particletype.h"
#include "pdgcode.h"
#include "photoncrosssections.h"
#include "pow.h"
#include "random.h"
#include "tabulation.h"

namespace Smash {

class ScatterActionPhoton : public ScatterAction {
 public:
  ScatterActionPhoton(const ParticleList &in, double time, int nofp)
      : ScatterAction(in[0], in[1], time),
        number_of_fractional_photons_(nofp),
        hadron_out_t_(outgoing_hadron_type(in)) {
    hadron_out_mass_ = hadron_out_t_->mass();
    reac_ = photon_reaction_type(in);
  }

  void generate_final_state() override;

  double raw_weight_value() const override { return weight_; }

  // returns the cross section of the underlying hadronic process (note: this
  // includes also the total cross section of the photon channels.)
  double cross_section() const override { return total_cross_section_; }

  // we have to override sample_masses from Action. In Action sample_masses
  // relies on the outgoing particles member. At the point where we need this
  // function, this member outgoing_particles_ of Action is not yet accessible.
  // ToDo: See if one can call the Action constructor in ScatterActionPhoton
  // constructor and set direct the outgoing particles (known at time of
  // constructor calling). But since all the inheritance stuff will change in
  // the near future we might as well go for now with the custom solution and
  // treat the photons in the special manner as they are meant.

  // actually we do not need a pair. second particle will always be a photon.
  double sample_out_hadron_mass(const ParticleTypePtr out_type);
  /** Overridden to effectively return the reaction channel. */
  ProcessType get_type() const override {
    return static_cast<ProcessType>(reac_);
  }

  /** Adds one dummy channel with a given cross-section. The intended use is to
   * add the hadronic cross-section from already performed hadronic action
   * without recomputing it. The photon action is never performed, so
   * this channel itself will never play any role. Only its cross-section will.
   */
  void add_dummy_hadronic_channels(double reaction_cross_section);

  /** To add only one reaction for testing purposes */
  void add_single_channel(bool from_check_collision = false) {
    add_processes<CollisionBranch>(photon_cross_sections(from_check_collision),
                                   collision_channels_photons_,
                                   cross_section_photons_);
  }

  enum class ReactionType {
    no_reaction,
    pi_z_pi_p_rho_p,
    pi_z_pi_m_rho_m,
    pi_p_rho_z_pi_p,
    pi_m_rho_z_pi_m,
    pi_m_rho_p_pi_z,
    pi_p_rho_m_pi_z,
    pi_z_rho_p_pi_p,
    pi_z_rho_m_pi_m,
    pi_p_pi_m_rho_z,
    pi_z_rho_z_pi_z
  };

  // the processes pi0 + rho => pi + y and pi + rho => pi0 + gamma can happen
  // via exchange of either (rho, a1, pi) or (omega). If MediatorType::SUM is
  // set,the cross section for both channels are added. If MediatorType::RHO /
  // OMEGA is set, only the channel where (rho, a1, pi) / (omega) is exchanged
  // will be used.
  enum class MediatorType { SUM, RHO, OMEGA };

  ReactionType reac_ = ReactionType::no_reaction;

  /// Tells if the given incoming particles may produce photon
  // static ReactionType is_photon_reaction(const ParticleList &in);
  static ReactionType photon_reaction_type(const ParticleList &in);

  // static bool is_photon_reaction(const double s, const ParticleList &in);

  static ParticleTypePtr outgoing_hadron_type(const ParticleList &in);

  static bool is_kinematically_possible(const double s, const ParticleList &in);

 private:
  CollisionBranchList collision_channels_photons_;

  int const number_of_fractional_photons_;

  ParticleTypePtr hadron_out_t_;

  double hadron_out_mass_;

  static constexpr MediatorType default_mediator_ = MediatorType::RHO;

  double weight_ = 0.0;

  /** List of possible collisions producing photons */

  double cross_section_photons_ = 0.0;

  double diff_cross_section(const double t, const double t2, const double t1,
                            const double m_rho,
                            MediatorType mediator = default_mediator_) const;

  double mediator_mass(ReactionType r) const;

  double mediator_mass() const;
  CollisionBranchList photon_cross_sections(
      bool from_check_collision = false, MediatorType mediator = default_mediator_);

  std::pair<double, double> diff_cross_section_single(const double t,
                                                      const double t2,
                                                      const double t1,
                                                      const double m_rho);

  std::pair<double, double> form_factor_single(const double E_photon);

  double form_factor(double E_photon);

  double diff_cross_section_w_ff(const double t, const double t2,
                                 const double t1, const double m_rho,
                                 const double E_photon);

  // conversion factor to millibarn
  const double to_mb = 0.3894;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONPHOTON_H_
