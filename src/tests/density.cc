/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "../include/configuration.h"
#include "../include/experiment.h"
#include "../include/modusdefault.h"

#include <boost/filesystem.hpp>

using namespace Smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "proton 0.938 0.0 2212\n");
}

static ParticleData create_proton(int id = -1) {
  return ParticleData{ParticleType::find(0x2212), id};
}

static ParticleData create_antiproton(int id = -1) {
  return ParticleData{ParticleType::find(-0x2212), id};
}


// create a particle list with particles moving in different directions.
static void create_particle_list(Particles &P) {
  // proton that moves along the x axis:
  ParticleData proton_x = create_proton();
  // proton that moves along the y axis:
  ParticleData proton_y = create_proton();
  // proton that moves along the z axis:
  ParticleData proton_z = create_proton();
  // proton that moves in neither of xyz planes, v = v0:
  ParticleData proton = create_proton();
  // proton that moves in neither of xyz planes, v = - v0:
  ParticleData antiproton = create_antiproton();

  double mass = 0.938;
  // set momenta:
  proton_x.set_4momentum(mass, ThreeVector(0.5, 0.0, 0.0));
  proton_y.set_4momentum(mass, ThreeVector(0.0, -.5, 0.0));
  proton_z.set_4momentum(mass, ThreeVector(0.0, 0.0, 0.5));
  proton.set_4momentum(mass, ThreeVector(0.4, 0.5, 0.7));
  antiproton.set_4momentum(mass, ThreeVector(-0.4, -0.5, -0.7));

  // set positions:
  proton_x.set_4position(FourVector(0.0, -5.0, 0.0, 0.0));
  proton_y.set_4position(FourVector(0.0, 0.0, 5.0, 0.0));
  proton_z.set_4position(FourVector(0.0, 0.0, 0.0, 5.0));
  proton.set_4position(FourVector(0.0,-4.0,-5.0,-7.0));
  antiproton.set_4position(FourVector(0.0, 4.0, 5.0, 7.0));

  // add particles (and make sure the particles get the correct ID):
  COMPARE(P.add_data(proton_x), 0);
  COMPARE(P.add_data(proton_y), 1);
  COMPARE(P.add_data(proton_z), 2);
  COMPARE(P.add_data(proton), 3);
  COMPARE(P.add_data(antiproton), 4);

  return;
}

TEST(density_compframe) {
  ModusDefault m;
  Particles Pdef;
  create_particle_list(Pdef);
  OutputsList out;
  // clock, output interval, cross-section, testparticles
  ExperimentParameters param{{0.f, 0.3f}, 1.f, 0.0, 1};
  double dx = 0.3;
  double dy = 0.3;
  double dz = 0.3;
  int nx = 20;
  int ny = 20;
  int nz = 20;
  double sigma = 0.8;
  ThreeVector r;
  FourVector jmu;
  FILE * pFile;
  ParticleList plist = ParticleList(Pdef.data().begin(), Pdef.data().end());
  /* for (const auto &p : plist) {
    printf("Momentum: %10.2f  %10.2f  %10.2f %10.2f\n", p.momentum().x0(),
                   p.momentum().x1(), p.momentum().x2(), p.momentum().x3());
  }
  */
  printf("Ready to write vtk files.\n");

  for (auto it = 0; it < 5; it++) {
    pFile = fopen(("bdens.vtk." + std::to_string(it)).c_str(), "w");
    // printf("Writing a header to file %d, pointer %p.\n", it, pFile);
    fprintf(pFile, "# vtk DataFile Version 2.0\n");
    fprintf(pFile, "baryon density\n");
    fprintf(pFile, "ASCII\n");
    fprintf(pFile, "DATASET STRUCTURED_POINTS\n");
    // printf("Continue writing a header to file %d.\n", it);
    fprintf(pFile, "DIMENSIONS %d %d %d\n", 2*nx+1, 2*ny+1, 2*nz+1);
    fprintf(pFile, "SPACING 1 1 1\n");
    fprintf(pFile, "ORIGIN %d %d %d\n", -nx, -ny, -nz);
    fprintf(pFile, "POINT_DATA %d\n", (2*nx+1)*(2*ny+1)*(2*nz+1) );
    fprintf(pFile, "SCALARS baryon_density float 1\n");
    fprintf(pFile, "LOOKUP_TABLE default\n");

    // printf("it = %d.\n", it);
    for (auto iz = -nz; iz <= nz; iz++) {
      for (auto iy = -ny; iy <= ny; iy++) {
        for (auto ix = -nx; ix <= nx; ix++) {
          r = ThreeVector(ix*dx, iy*dy, iz*dz);
          jmu = m.baryon_jmu(r, plist, sigma);
          fprintf(pFile,"%f ", jmu.x0());
        }
        fprintf(pFile,"\n");
      }
    }
    fclose(pFile);
    m.propagate(&Pdef, param, out);
  }

}
