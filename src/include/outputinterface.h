/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_OUTPUTINTERFACE_H_
#define SRC_INCLUDE_OUTPUTINTERFACE_H_

#include "forwarddeclarations.h"

namespace Smash {
class Particles;

/**
 * \brief Abstraction of generic output
 * Any output should inherit this class. It provides virtual methods that will be called at predefined moments:
 * 1) At event start and event end
 * 2) After every N'th timestep
 * 3) Before and after collision
 */
class OutputInterface {
 public:
  virtual ~OutputInterface() = default;

  /**
   * Output launched at event start after initialization, when particles are generated, but not yet propagated.
   */
  virtual void at_eventstart(const Particles &, const int) = 0;

  /**
   * Output launched at event end. Event end is determined by maximal timestep option.
   */
  virtual void at_eventend(const Particles &, const int ) = 0;

  /**
   *Output before any collision or decay should be added when Collision class is set up.
   */
  virtual void before_collision() = 0;

  /**
   *Output after any collision or decay should be added when Collision class is set up.
   */
  virtual void after_collision() = 0;

  /**
   * Called whenever an action modified one or more particles.
   *
   * \param incoming_particles The list of particles before the Action was
   *                          performed.
   * \param outgoing_particles   The list of particles after the Action was
   *                          performed.
   */
  virtual void write_interaction(const ParticleList &/*incoming_particles*/,
                                 const ParticleList &/*outgoing_particles*/) {}

  /**
   * Output launched after every N'th timestep. N is controlled by an option.
   */
  virtual void after_Nth_timestep(const Particles &, const int, const int) = 0;

};

}  // namespace Smash

#endif  // SRC_INCLUDE_OUTPUTINTERFACE_H_
