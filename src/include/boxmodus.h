/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_BOXMODUS_H_
#define SRC_INCLUDE_BOXMODUS_H_

#include <stdint.h>
#include <cmath>
#include <list>

#include "include/crosssections.h"
#include "include/modusdefault.h"
#include "include/particles.h"

namespace Smash {

class BoxModus;
class Configuration;
struct ExperimentParameters;

/** BoxModus: Provides a modus for infinite matter calculations
 *
 * Matter is confined in a cubical box. Depending on the initial
 * condition, particles are either reflected on the boundaries or
 * inserted on opposite positions.
 *
 * To use this modus, chose
 * \code
 * General:
 *      MODUS: Box
 * \endcode
 * in the configuration file.
 *
 * Options for BoxModus go in the "Modi"→"Box" section of the
 * configuration:
 *
 * \code
 * Modi:
 *      Box:
 *              # definitions here
 * \endcode
 *
 * The following directives are understood:
 *
 * ...
 * TODO(baeuchle) needs to be completed.
 * 
 * 
 *
 *
 */
class BoxModus : public ModusDefault {
 public:
  /// Gathers all configuration variables for the Box.
  explicit BoxModus(Configuration modus_config,
           const ExperimentParameters &parameters);

  /** Prints some information about the initialization of BoxModus
   *
   * \see ModusDefalt::print_startup()
   */
  void print_startup();  // TODO(mkretz): needs to be discoverable from an
                         // outside "printer"

  /** creates initial conditions from the particles.
   */
  void initial_conditions(Particles *particles,
                          const ExperimentParameters &parameters);

  /** Enforces that all particles are inside the box
   *
   * \see enforce_periodic_boundaries
   */
  int sanity_check(Particles *particles);

  /** Propagates all particles through the box.
   *
   * \copydetails ModusDefault::propagate()
   * \param[in] output_list output objects
   *
   * In addition to the usual propagation, particles can touch the wall
   * in BoxModus and either be reflected or inserted opposite. In that
   * case, the OutputsList will be used.
   */
  void propagate(Particles *particles, const ExperimentParameters &parameters,
                                       const OutputsList &output_list);

 private:
  /// initial condition: ??? TODO(baeuchle): check what this is for.
  const int initial_condition_;
  /// length of the cube's edge in fm/c
  const float length_;
  /// Temperature of the Box in GeV
  const float temperature_;
  /// initial number density of the Box TODO(baeuchle): What is this
  /// for? Why isn't the density fixed by the temperature?
  float number_density_initial_ = 0.f;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_BOXMODUS_H_
