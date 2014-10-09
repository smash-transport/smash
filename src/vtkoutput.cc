/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <memory>
#include <include/config.h>
#include "include/clock.h"
#include "include/filedeleter.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"
#include "include/vtkoutput.h"

namespace Smash {

VtkOutput::VtkOutput(bf::path path, Configuration&& /*conf*/)
  : base_path_(std::move(path)), vtk_output_counter_(0) {}
/*!\Userguide
 * \page input_vtk VTK
 *
 * Writes snapshots of simulated particles at fixed moments of time
 * to separate .vtk files. These fixed moments of time are event start,
 * event end and every next time interval \f$\Delta t\f$, where
 * \f$\Delta t\f$ is controlled by an option.
 * Produced output can be opened by paraview
 * and used for an easy visualization of the simulation.
 */

VtkOutput::~VtkOutput() {
}

void VtkOutput::at_eventstart(const Particles &particles,
                              const int event_number) {
  vtk_output_counter_ = 0;
  write(particles, event_number);
  vtk_output_counter_++;
}

void VtkOutput::at_eventend(const Particles &/*particles*/,
                            const int /*event_number*/) {
}

void VtkOutput::at_intermediate_time(const Particles &particles,
                                   const int event_number,
                                   const Clock& /*clock*/) {
  write(particles, event_number);
  vtk_output_counter_++;
}

void VtkOutput::write(const Particles &particles, const int event_number) {
  char filename[32];
  snprintf(filename, sizeof(filename), "pos_ev%05i_tstep%05i.vtk", event_number,
           vtk_output_counter_);
  FilePtr file_{fopen((base_path_ / filename).native().c_str(), "w")};

  /* Legacy VTK file format */
  fprintf(file_.get(), "# vtk DataFile Version 2.0\n");
  fprintf(file_.get(), "Generated from molecular-offset data %s\n",
                                                        VERSION_MAJOR);
  fprintf(file_.get(), "ASCII\n");

  /* Unstructured data sets are composed of points, lines, polygons, .. */
  fprintf(file_.get(), "DATASET UNSTRUCTURED_GRID\n");
  fprintf(file_.get(), "POINTS %zu double\n", particles.size());
  for (const auto &p : particles.data()) {
    fprintf(file_.get(), "%g %g %g\n", p.position().x1(),
            p.position().x2(), p.position().x3());
  }
  fprintf(file_.get(), "CELLS %zu %zu\n",
                       particles.size(), particles.size() * 2);
  for (size_t point_index = 0; point_index < particles.size(); point_index++) {
    fprintf(file_.get(), "1 %zu\n", point_index);
  }
  fprintf(file_.get(), "CELL_TYPES %zu\n", particles.size());
  for (size_t point_index = 0; point_index < particles.size(); point_index++) {
    fprintf(file_.get(), "1\n");
  }
  fprintf(file_.get(), "POINT_DATA %zu\n", particles.size());
  fprintf(file_.get(), "SCALARS pdg_codes int 1\n");
  fprintf(file_.get(), "LOOKUP_TABLE default\n");
  for (const auto &p : particles.data()) {
    fprintf(file_.get(), "%s\n", p.pdgcode().string().c_str());
  }
  fprintf(file_.get(), "VECTORS momentum double\n");
  for (const auto &p : particles.data()) {
    fprintf(file_.get(), "%g %g %g\n", p.momentum().x1(),
            p.momentum().x2(), p.momentum().x3());
  }
}

}  // namespace Smash
