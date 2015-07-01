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
   * \param incoming_particles The list of particles before the Action was
   *                          performed.
   * \param outgoing_particles   The list of particles after the Action was
   *                          performed.
   * \param density The density at the interaction point
   * \param total_cross_section The total cross section of this interaction
   *                           or total width in case of decays
   * \param process_type Identifier for the type of process, e.g.
   *  		elastic scattering, resonance formation,...
   *
   * \fpPrecision Why \c double?
   */
  virtual void at_interaction(const ParticleList &incoming_particles,
                              const ParticleList &outgoing_particles,
                              const double density,
                              const double total_cross_section,
                              const ProcessType process_type) {
    SMASH_UNUSED(incoming_particles);
    SMASH_UNUSED(outgoing_particles);
    SMASH_UNUSED(density);
    SMASH_UNUSED(total_cross_section);
    SMASH_UNUSED(process_type);
  }

  /**
   * Output launched after every N'th timestep. N is controlled by an option.
   */
  virtual void at_intermediate_time(const Particles &, const int,
                                  const Clock &) = 0;

  /**
   * Output intended for writing out thermodynamics.
   * It is launched after every N'th timestep. N is controlled by an option.
   */
  virtual void thermodynamics_output(const Particles &particles,
                                     const ExperimentParameters &param) {
    SMASH_UNUSED(particles);
    SMASH_UNUSED(param);
  }

  /**
   * Output to write thermodynamics from the lattice.
   */
  virtual void thermodynamics_output(const std::string varname,
                            RectangularLattice<DensityOnLattice> &lattice,
                            const int event_number) {
    SMASH_UNUSED(varname);
    SMASH_UNUSED(lattice);
    SMASH_UNUSED(event_number);
  }

  // probably not necassary and will be removed
  // at_interaction could be used for dileptons in finder
  virtual void dileptons(const ParticleList &incoming_particles,
                         const ParticleList &outgoing_particles,
                         float shining_weight) {
    SMASH_UNUSED(incoming_particles);
    SMASH_UNUSED(outgoing_particles);
    SMASH_UNUSED(shining_weight);
  }

};

}  // namespace Smash

#endif  // SRC_INCLUDE_OUTPUTINTERFACE_H_
