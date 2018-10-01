/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This has to be the first include

#include "setup.h"

#include <smash/config.h>
#include <array>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <iostream>
#include <vector>

#include "../include/smash/clock.h"
#include "../include/smash/configuration.h"
#include "../include/smash/outputinterface.h"
#include "../include/smash/particles.h"
#include "../include/smash/random.h"
#include "../include/smash/vtkoutput.h"

using namespace smash;

static const double accuracy = 1.0e-4;
static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);
static auto random_value = random::make_uniform_distribution(-15.0, +15.0);

TEST(directory_is_created) {
  bf::create_directories(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

static void compare_threevector(const std::array<std::string, 3> &stringarray,
                                const ThreeVector &threevector) {
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(0).c_str()), threevector.x1(),
                         accuracy);
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(1).c_str()), threevector.x2(),
                         accuracy);
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(2).c_str()), threevector.x3(),
                         accuracy);
}

TEST(vtkoutputfile) {
  /* Create particles */
  Test::create_smashon_particletypes();

  Particles particles;
  const int number_of_particles = 5;
  for (int i = 0; i < number_of_particles; i++) {
    particles.insert(Test::smashon_random());
  }

  /* Create output object */
  std::unique_ptr<VtkOutput> vtkop =
      make_unique<VtkOutput>(testoutputpath, "Particles");
  int event_id = 0;
  /* Initial output */
  vtkop->at_eventstart(particles, event_id);
  const bf::path outputfilename = "pos_ev00000_tstep00000.vtk";
  const bf::path outputfilepath = testoutputpath / outputfilename;
  VERIFY(bf::exists(outputfilepath));
  /* Time step output */
  Clock clock(0.0, 1.0);
  DensityParameters dens_par(Test::default_parameters());
  vtkop->at_intermediate_time(particles, clock, dens_par);
  const bf::path outputfile2path =
      testoutputpath / "pos_ev00000_tstep00001.vtk";
  VERIFY(bf::exists(outputfile2path));

  {
    bf::fstream outputfile;
    outputfile.open(outputfilepath, std::ios_base::in);
    if (outputfile.good()) {
      std::string line, item;
      /* Check header */
      std::string output_header = "";
      std::string header =
          "# vtk DataFile Version 2.0\n"
          "Generated from molecular-offset data " VERSION_MAJOR
          "\n"
          "ASCII\n"
          "DATASET UNSTRUCTURED_GRID\n";
      do {
        std::getline(outputfile, line);
        output_header += line + '\n';
      } while (line != "DATASET UNSTRUCTURED_GRID" && !outputfile.eof());
      COMPARE(output_header, header);
      /* Check position information */
      outputfile >> item;
      COMPARE(item, "POINTS");
      outputfile >> item;
      COMPARE(std::stoul(item), particles.size());
      outputfile >> item;
      COMPARE(item, "double");
      COMPARE(particles.size(), size_t(number_of_particles));
      for (const auto &pd : particles) {
        std::array<std::string, 3> position_string;
        for (int j = 0; j < 3; j++) {
          outputfile >> item;
          position_string[j] = item;
        }
        compare_threevector(position_string, pd.position().threevec());
      }
      /* Check cell information */
      outputfile >> item;
      COMPARE(item, "CELLS");
      outputfile >> item;
      COMPARE(std::stoul(item), particles.size());
      outputfile >> item;
      COMPARE(std::stoul(item), particles.size() * 2);
      for (int i = 0; i < number_of_particles; i++) {
        outputfile >> item;
        COMPARE(std::atoi(item.c_str()), 1);
        outputfile >> item;
        COMPARE(std::atoi(item.c_str()), i);
      }
      outputfile >> item;
      COMPARE(item, "CELL_TYPES");
      outputfile >> item;
      COMPARE(std::stoul(item), particles.size());
      for (int i = 0; i < number_of_particles; i++) {
        outputfile >> item;
        COMPARE(std::atoi(item.c_str()), 1);
      }
      /* Check point data */
      outputfile >> item;
      COMPARE(item, "POINT_DATA");
      outputfile >> item;
      COMPARE(std::stoul(item), particles.size());
      outputfile >> item;
      COMPARE(item, "SCALARS");
      outputfile >> item;
      COMPARE(item, "pdg_codes");
      outputfile >> item;
      COMPARE(item, "int");
      outputfile >> item;
      COMPARE(item, "1");
      outputfile >> item;
      COMPARE(item, "LOOKUP_TABLE");
      outputfile >> item;
      COMPARE(item, "default");
      for (int i = 0; i < number_of_particles; i++) {
        outputfile >> item;
        COMPARE(item, "661");
      }
      /* Check SCALARS is_formed */
      outputfile >> item;
      COMPARE(item, "SCALARS");
      outputfile >> item;
      COMPARE(item, "is_formed");
      outputfile >> item;
      COMPARE(item, "int");
      outputfile >> item;
      COMPARE(item, "1");
      outputfile >> item;
      COMPARE(item, "LOOKUP_TABLE");
      outputfile >> item;
      COMPARE(item, "default");
      for (int i = 0; i < number_of_particles; i++) {
        outputfile >> item;
        COMPARE(item, "1");
      }
      /* Check SCALARS cross_section_scaling_factor */
      outputfile >> item;
      COMPARE(item, "SCALARS");
      outputfile >> item;
      COMPARE(item, "cross_section_scaling_factor");
      outputfile >> item;
      COMPARE(item, "double");
      outputfile >> item;
      COMPARE(item, "1");
      outputfile >> item;
      COMPARE(item, "LOOKUP_TABLE");
      outputfile >> item;
      COMPARE(item, "default");
      for (int i = 0; i < number_of_particles; i++) {
        outputfile >> item;
        COMPARE(item, "1");
      }
      /* Check SCALARS mass */
      outputfile >> item;
      COMPARE(item, "SCALARS");
      outputfile >> item;
      COMPARE(item, "mass");
      outputfile >> item;
      COMPARE(item, "double");
      outputfile >> item;
      COMPARE(item, "1");
      outputfile >> item;
      COMPARE(item, "LOOKUP_TABLE");
      outputfile >> item;
      COMPARE(item, "default");
      for (int i = 0; i < number_of_particles; i++) {
        outputfile >> item;
        COMPARE(item, "0.123");
      }
      /* Check momentum vectors */
      outputfile >> item;
      COMPARE(item, "VECTORS");
      outputfile >> item;
      COMPARE(item, "momentum");
      outputfile >> item;
      COMPARE(item, "double");
      COMPARE(particles.size(), size_t(number_of_particles));
      for (const auto &pd : particles) {
        std::array<std::string, 3> momentum_string;
        for (int j = 0; j < 3; j++) {
          outputfile >> item;
          momentum_string[j] = item;
        }
        compare_threevector(momentum_string, pd.momentum().threevec());
      }
    }
  }
  VERIFY(bf::remove(outputfilepath));
  VERIFY(bf::remove(outputfile2path));
}
