/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/icoutput.h"

#include <fstream>

#include <boost/filesystem.hpp>

#include "smash/action.h"

namespace smash {

/*!\Userguide
 * \page IC_output_user_guide_ ASCII IC Output
 * The initial conditions output (SMASH_IC.dat) is used to create a hypersurface
 * of constant proper time. This output can be applied as an initial state
 * for hydrodynamic models. Curently, it is only available in ASCII format.
 * The provided particle data is printed in the computational frame.\n
 * \n
 * The proper time, at which the hypersurface is constructed can either be set
 * explicitly in the configuration file or determined from the collision system.
 * \n
 * By default, the proper time corresponds to the moment when the two nuclei
 * have entirely passed through each other:
 * \f$ \tau_0 = (r_a \ + \ r_b) \ \left(\left(\frac{\sqrt{s_\mathrm{NN}}}
 * {2 \ m_N}\right)^2 - 1\right)^{-1/2} \f$ \n
 * \n
 * The format of the file is the following: \n
 *
 * \n
 * **Header**
 * \code
 * # **smash_version** initial conditions: hypersurface of constant proper time
 * # tau x y eta mt px py Rap pdg ID charge
 * # fm fm fm none GeV GeV GeV none none none e
 * \endcode
 * The header consists of 3 lines starting with a '#', containing the following
 * information:
 * -# SMASH-version and the information, that 'initial conditions' are provided
 * -# The header with all column names
 * -# The units of all column quantities
 *
 * **Output block header**
 *
 * The initial conditions output is, similar to the OSCAR output, based on a
 * a block structure, where each block consists of 1 event. The header for a
 * new event is structured as follows:
 * \code
 * # event ev_num start
 * \endcode
 * where
 * \li \key ev_num: The number of the current event
 *
 * Note that 'event' and 'start' are no variables, but words that are
 * printed in the header. \n
 * \n
 * **Particle line**
 *
 * The particle lines are formatted as follows:
 * \code
 * tau x y eta mt px py Rap pdg ID charge
 * \endcode
 * where
 * \li \key tau: Proper time of the particle
 * \li \key x, \key y: Cartesian x and y coordinates of the particle
 * \li \key eta: Space-time rapidity of the particle
 * \li \key mt: Transverse mass of the particle
 * \li \key px, \key py: x and y components of the particle's momentum
 * \li \key Rap: Momentum space rapidity of the particle
 * \li \key pdg: PDG code of the particle (see http://pdg.lbl.gov/). It contains
 * all quantum numbers and uniquely identifies its type.
 * \li \key ID: Particle identifier in terms of an integer. It is unique for
 * every particle in the event.
 * \li \key charge: electric charge of the particle
 *
 * **Event end line**
 *
 * The end of an event is indicated by the following line:
 * \code
 * # event ev_num end
 * \endcode
 * where
 * \li \key ev_num: The number of the current event
 *
 * Note that 'event' and 'end' are no variables, but words that are
 * printed in the header. \n
 *
 */

ICOutput::ICOutput(const bf::path &path, const std::string &name,
                   const OutputParameters &out_par)
    : OutputInterface(name),
      file_{path / "SMASH_IC.dat", "w"},
      out_par_(out_par) {
  std::fprintf(file_.get(), "# %s initial conditions: hypersurface of constant proper time\n", VERSION_MAJOR);
  std::fprintf(file_.get(), "# tau x y eta mt px py Rap pdg ID charge \n");
  std::fprintf(file_.get(), "# fm fm fm none GeV GeV GeV none none none e \n");
}

ICOutput::~ICOutput() {}

void ICOutput::at_eventstart(const Particles &particles,
                             const int event_number) {
  std::fprintf(file_.get(), "# event %i start\n", event_number + 1);
  SMASH_UNUSED(particles);
};

void ICOutput::at_eventend(const Particles &particles, const int event_number,
                           double impact_parameter, bool empty_event) {
  const auto &log = logger<LogArea::HyperSurfaceCrossing>();
  std::fprintf(file_.get(), "# event %i end\n", event_number + 1);

  // If the runtime is too short some particles might not yet have
  // reached the hypersurface. Warning is printed.
  bool runtime_too_short = false;
  if (particles.size() != 0) {
    for (auto &p : particles) {
      double tau = p.position().tau();
      double t = p.position().x0();
      double z = p.position().x3();
      if ((tau < IC_proper_time_) || (fabs(t) < fabs(z))) {
        // If t < z, tau = sqrt(t^2 - z^2) returns NAN. Those particles are also
        // below the hypersurface and need further propagation
        runtime_too_short = true;
      }
    }
  }

  if (runtime_too_short) {
    log.warn(
        "End time might be too small. Hypersurface has not yet been crossed "
        "by ",
        particles.size(), " particle(s).");
  }
  SMASH_UNUSED(particles);
  SMASH_UNUSED(impact_parameter);
  SMASH_UNUSED(empty_event);
};

void ICOutput::at_intermediate_time(const Particles &particles,
                                    const Clock &clock,
                                    const DensityParameters &dens_param) {
  // Dummy, but virtual function needs to be declared.
  // It is never actually used.
  SMASH_UNUSED(particles);
  SMASH_UNUSED(clock);
  SMASH_UNUSED(dens_param);
}

void ICOutput::at_interaction(const Action &action, const double density) {
  SMASH_UNUSED(density);
  assert(action.get_type() == ProcessType::HyperSurfaceCrossing);
  assert(action.incoming_particles().size() == 1);

  ParticleData particle = action.incoming_particles()[0];

  // transverse mass
  double m_trans = sqrt(particle.type().mass() * particle.type().mass() +
                        particle.momentum()[1] * particle.momentum()[1] +
                        particle.momentum()[2] * particle.momentum()[2]);
  // momentum space rapidity
  double rapidity =
      0.5 * log((particle.momentum()[0] + particle.momentum()[3]) /
                (particle.momentum()[0] - particle.momentum()[3]));

  // write particle data
  std::fprintf(file_.get(), "%g %g %g %g %g %g %g %g %s %i %i \n",
               particle.position().tau(), particle.position()[1],
               particle.position()[2], particle.position().eta(), m_trans,
               particle.momentum()[1], particle.momentum()[2], rapidity,
               particle.pdgcode().string().c_str(), particle.id(),
               particle.type().charge());

  if (IC_proper_time_ < 0.0) {
    // First particle that is removed, overwrite negative default
    IC_proper_time_ =  particle.position().tau();
  } else {
    // Verify that all other particles have the same proper time
    double next_proper_time = particle.position().tau();
    if (!((next_proper_time - IC_proper_time_) < really_small))
      throw std::runtime_error("Hypersurface proper time changed during evolution.");
  }
}

}  // namespace smash
