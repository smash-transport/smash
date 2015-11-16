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

#include <map>

#include "forwarddeclarations.h"
#include "density.h"
#include "lattice.h"
#include "macros.h"
#include "processbranch.h"

namespace Smash {

/**
 * \ingroup output
 *
 * \brief Abstraction of generic output
 * Any output should inherit this class. It provides virtual methods that will
 * be called at predefined moments:
 * 1) At event start and event end
 * 2) After every N'th timestep
 * 3) Before and after collision
 */
class OutputInterface {
 public:
  virtual ~OutputInterface() = default;

  /**
   * Output launched at event start after initialization, when particles are
   * generated, but not yet propagated.
   */
  virtual void at_eventstart(const Particles &, const int) = 0;

  /**
   * Output launched at event end. Event end is determined by maximal timestep
   * option.
   */
  virtual void at_eventend(const Particles &, const int) = 0;

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
   */
  virtual void at_intermediate_time(const Particles &, const Clock &) {
  }

  /**
   * Output intended for writing out thermodynamics.
   * It is launched after every N'th timestep. N is controlled by an option.
   */
  virtual void thermodynamics_output(const Particles &particles,
                                     const ExperimentParameters &param,
                                     const DensityParameters &dens_param) {
    SMASH_UNUSED(particles);
    SMASH_UNUSED(param);
    SMASH_UNUSED(dens_param);
  }

  /**
   * Output to write thermodynamics from the lattice.
   */
  virtual void thermodynamics_output(const std::string varname,
                            RectangularLattice<DensityOnLattice> &lattice) {
    SMASH_UNUSED(varname);
    SMASH_UNUSED(lattice);
  }
};

}  // namespace Smash

#endif  // SRC_INCLUDE_OUTPUTINTERFACE_H_
