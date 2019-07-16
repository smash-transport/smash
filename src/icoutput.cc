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

ICOutput::ICOutput(const bf::path &path, const std::string &name,
                   const OutputParameters &out_par)
    : OutputInterface(name),
      file_{path / "SMASH_IC.dat", "w"},
      out_par_(out_par) {
  const double prop_time = out_par_.IC_proper_time;
  std::fprintf(file_.get(), "# %s initial conditions\n", VERSION_MAJOR);
  std::fprintf(file_.get(), "# @ proper time: %7.4f fm \n", prop_time);
  std::fprintf(file_.get(), "# tau x y eta mt px py Rap pdg ID charge \n");
  std::fprintf(file_.get(), "# fm fm fm none GeV GeV GeV none none none e \n");
}

ICOutput::~ICOutput() {}

void ICOutput::at_eventstart(const Particles &particles,
                             const int event_number) {
  // dummy, but virtual function needs to be declared, this function is never
  // actually used
  SMASH_UNUSED(particles);
  SMASH_UNUSED(event_number);
};

void ICOutput::at_eventend(const Particles &particles, const int event_number,
                           double impact_parameter, bool empty_event) {
  // dummy, but virtual function needs to be declared, this function is never
  // actually used
  SMASH_UNUSED(particles);
  SMASH_UNUSED(event_number);
  SMASH_UNUSED(impact_parameter);
  SMASH_UNUSED(empty_event);
};

void ICOutput::at_intermediate_time(const Particles &particles,
                                    const Clock &clock,
                                    const DensityParameters &dens_param) {
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
}

}  // namespace smash
