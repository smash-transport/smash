/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#ifndef SRC_INCLUDE_MODUSDEFAULT_H_
#define SRC_INCLUDE_MODUSDEFAULT_H_

#include <stdexcept>
#include <stdint.h>
#include <list>

#include "include/particles.h"

namespace Smash {

/* forward declarations */
class Particles;
class CrossSections;
struct ExperimentParameters;

/**
 * Baseclass for Modus classes that provides default function implementations.
 *
 * This is only a base class for actual Modus classes. Meaning there will never
 * be objects, references, or pointers to ModusDefault. Therefore, it does not
 * have - and will never need any virtual functions.
 *
 * The rules for adding functions to this class are as follows:
 * - This class is empty per default.
 * - You can add a function if you have a function that is different in at least
 *   two subclasses and equal in at least two subclasses.
 * - Code that is common to all goes into ExperimentImplementation.
 */
class ModusDefault {
 public:
  // never needs a virtual destructor

  // Missing functions for concrete Modus implementations:
  // void initial_conditions(Particles *particles);

  /** Enforces sensible positions for the particles.
   *
   * Currently, this is only needed for BoxModus; the other Modi do
   * nothing.
   *
   * \see BoxModus::sanity_check
   */
  int sanity_check(Particles * /*p*/) { return 0; }

  /** Propagates the positions of all particles on a straight line
   * through the current time step.
   *
   * For each particle, the position is shifted:
   * \f[\vec x^\prime = \vec x + \vec v \cdot \Delta t\f]
   * where \f$\vec x\f$ is the current position, \f$\vec v\f$ its
   * velocity and \f$\Delta t\f$ the duration of this timestep.
   *
   * \param[in,out] particles The particle list in the event
   * \param[in] parameters parameters for the experiment
   */
  void propagate(Particles *particles, const ExperimentParameters &parameters,
                                       const OutputsList &);

  /** BadInput is an error to throw if the configuration options are invalid.
   *
   **/
  struct BadInput : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };
};

}  // namespace Smash

#endif  // SRC_INCLUDE_MODUSDEFAULT_H_
