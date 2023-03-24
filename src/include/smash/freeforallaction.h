/*
 *
 *    Copyright (c) 2022-2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_FREEFORALLACTION_H_
#define SRC_INCLUDE_SMASH_FREEFORALLACTION_H_

#include "action.h"

namespace smash {

/**
 * \ingroup action
 * \brief Action class to create any incoming/outgoing particle combination
 * freely. This class is in particular designed to add and remove particles from
 * the evolution. This introduces violations of the conservation laws, but it is
 * needed for a concurrent running of transport and hydrodynamics.
 */
class FreeforallAction : public Action {
 public:
  /**
   * The FreeforallAction is able to add particles by providing an empty
   * particle list for \c in_part and a list with particles which are supposed to
   * be added to the evolution as \c out_part .
   * If particles should be removed a filled particle list with the corresponding
   * particles should be provided as \c in_part and an empty list as \c out_part .
   *
   * \param[in] in_part List of incoming particles
   * \param[in] out_part List of outgoing particles
   * \param[in] absolute_labframe_time Absolute time at which the action is
   * supposed to take place
   *
   * \note This action takes advantage of how the parent class works in order to
   * fake 0 &rarr; 1 and 1 &rarr; 0 processes. If particles need to be both added
   * and removed one needs to create and perform this action twice. Said
   * differently, the provided lists cannot be non-empty at the same time.
   */
  FreeforallAction(const ParticleList &in_part, const ParticleList &out_part,
                   double absolute_labframe_time)
      : Action(in_part, out_part, absolute_labframe_time,
               ProcessType::Freeforall) {}

  double get_total_weight() const { return 0.0; }
  double get_partial_weight() const { return 0.0; }

  /**
   * Generate a final state for the FreeforallAction in the sense that the
   * time and position of the particles in the list is scrolled back to the
   * action time. The outgoing particles are set in the constructor.
   */
  void generate_final_state() {
    // Set time for arbitrary outgoing particles to time of action
    for (auto &particle : outgoing_particles_) {
      const double t = particle.position().x0();
      const FourVector u(1.0, particle.velocity());
      particle.set_formation_time(t);
      particle.set_4position(particle.position() +
                             u * (time_of_execution() - t));
    }
  }

  /**
   * Function for debug output of incoming and outgoing particles from
   * freeforall action
   * \param[in] out Location of the output stream
   */
  void format_debug_output(std::ostream &out) const {
    out << "Freeforall action of " << incoming_particles_.size() << " to "
        << outgoing_particles_.size() << " particles.";
  }
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_FREEFORALLACTION_H_
