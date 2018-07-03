/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_BINARYOUTPUTPARTICLES_H_
#define SRC_INCLUDE_BINARYOUTPUTPARTICLES_H_

#include <string>

#include "binaryoutputcollisions.h"
#include "forwarddeclarations.h"
#include "outputparameters.h"

namespace smash {

/**
 * \ingroup output
 *
 * \brief Writes the particle list at specific times to the binary file
 *
 * This class writes the current particle list at a specific time t
 * to the binary output file. This specific time can
 * be: event start, event end or every next time interval \f$\Delta t \f$.
 * Writing (or not writing) the output at these moments is controlled by
 * different options. The time interval \f$\Delta t \f$ is also regulated by an
 * option. The output file is binary and has a block structure.
 *
 * Details of the output format can be found
 * on the wiki in the User Guide section, look for binary output.
 */

class BinaryOutputParticles : public BinaryOutputBase {
 public:
  /**
   * Create binary particle output.
   *
   * \param[in] path Output path.
   * \param[in] name Name of the ouput.
   * \param[in] out_par A structure containing the parameters of the output.
   */
  BinaryOutputParticles(const bf::path &path, std::string name,
                        const OutputParameters &out_par);

  /**
   * Writes the initial particle information of an event to the binary output.
   * \param[in] particles Current list of all particles.
   * \param[in] event_number Unused, needed since inherited.
   */
  void at_eventstart(const Particles &particles,
                     const int event_number) override;

  /**
   * Writes the final particle information of an event to the binary output.
   * \param[in] particles Current list of particles.
   * \param[in] event_number Number of event.
   * \param[in] impact_parameter Impact parameter of this event.
   */
  void at_eventend(const Particles &particles, const int event_number,
                   double impact_parameter) override;

  /**
   * Writes particles at each time interval; fixed by option OUTPUT_INTERVAL.
   * \param[in] particles Current list of particles.
   * \param[in] clock Unused, needed since inherited.
   * \param[in] dens_param Unused, needed since inherited.
   */
  void at_intermediate_time(const Particles &particles, const Clock &clock,
                            const DensityParameters &dens_param) override;

 private:
  /// Write only final particles (True) or both, inital and final (False).
  bool only_final_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_BINARYOUTPUTPARTICLES_H_
