/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_OUTPUTINTERFACE_H_
#define SRC_INCLUDE_OUTPUTINTERFACE_H_

#include "forwarddeclarations.h"
#include "density.h"
#include "lattice.h"
#include "macros.h"

namespace Smash {

/**
 * \ingroup output
 *
 * \brief Abstraction of generic output
 * Any output should inherit this class. It provides virtual methods that will
 * be called at predefined moments:
 * 1) At event start and event end
 * 2) After every N'th timestep
 * 3) At each interaction
 */
class OutputInterface {
 public:
  virtual ~OutputInterface() = default;

  /**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   * \param particles List of particles.
   * \param event_number Number of the current event.
   */
  virtual void at_eventstart(const Particles &particles,
                             const int event_number) = 0;

  /**
   * Output launched at event end. Event end is determined by maximal timestep
   * option.
   * \param particles List of particles.
   * \param event_number Number of the current event.
   */
  virtual void at_eventend(const Particles &particles,
                           const int event_number) = 0;

  /**
   * Called whenever an action modified one or more particles.
   *
   * \param action The action object, containing the initial and final state etc.
   * \param density The density at the interaction point.
   *
   * \fpPrecision Why \c double?
   */
  virtual void at_interaction(const Action &action, const double density) {
    SMASH_UNUSED(action);
    SMASH_UNUSED(density);
  }

  /**
   * Output launched after every N'th timestep. N is controlled by an option.
   * \param particles List of particles.
   * \param clock System clock.
   * \param dens_param Parameters for density calculation.
   */
  virtual void at_intermediate_time(const Particles &particles,
                                    const Clock &clock,
                                    const DensityParameters &dens_param) {
    SMASH_UNUSED(particles);
    SMASH_UNUSED(clock);
    SMASH_UNUSED(dens_param);
  }

  /**
   * Output to write thermodynamics from the lattice.
   * \param varname Variable name, used for file name etc.
   * \param lattice Lattice of tabulated values.
   */
  virtual void thermodynamics_output(const std::string &varname,
                            RectangularLattice<DensityOnLattice> &lattice) {
    SMASH_UNUSED(varname);
    SMASH_UNUSED(lattice);
  }
};

}  // namespace Smash

#endif  // SRC_INCLUDE_OUTPUTINTERFACE_H_
