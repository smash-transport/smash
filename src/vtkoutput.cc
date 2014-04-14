/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/vtkoutput.h"
#include "include/particles.h"
#include "include/filedeleter.h"
#include <memory>

namespace Smash {

VtkOutput::VtkOutput(boost::filesystem::path path)
    : base_path_(std::move(path)) {
}

VtkOutput::~VtkOutput() {
}

void VtkOutput::write_state(const Particles &particles) {
  char filename[32];
  snprintf(
      filename, sizeof(filename), "pos_0.%05i.vtk",
      static_cast<int>((particles.data().begin()->position().x0() - 1.0) * 10));
  std::unique_ptr<std::FILE> file_{
      fopen((base_path_ / filename).native().c_str(), "w")};

  /* Legacy VTK file format */
  fprintf(file_.get(), "# vtk DataFile Version 2.0\n");
  fprintf(file_.get(), "Generated from molecular-offset data\n");
  fprintf(file_.get(), "ASCII\n");

  /* Unstructured data sets are composed of points, lines, polygons, .. */
  fprintf(file_.get(), "DATASET UNSTRUCTURED_GRID\n");
  fprintf(file_.get(), "POINTS %zu double\n", particles.size());
  for (const auto &p : particles.data()) {
    fprintf(file_.get(), "%g %g %g\n", p.position().x1(),
            p.position().x2(), p.position().x3());
  }
  fprintf(file_.get(), "CELLS %zu %zu\n", particles.size(), particles.size() * 2);
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
    fprintf(file_.get(), "%x\n", p.pdgcode().code());
  }
  fprintf(file_.get(), "VECTORS momentum double\n");
  for (const auto &p : particles.data()) {
    fprintf(file_.get(), "%g %g %g\n", p.momentum().x1(),
            p.momentum().x2(), p.momentum().x3());
  }
}

}  // namespace Smash
