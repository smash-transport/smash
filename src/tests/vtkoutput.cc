/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include <cstdio>
#include <cstring>
#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include "../include/clock.h"
#include "../include/config.h"
#include "../include/configuration.h"
#include "../include/outputinterface.h"
#include "../include/particles.h"
#include "../include/random.h"
#include "../include/vtkoutput.h"

using namespace Smash;

static const float accuracy = 1.0e-4;
static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);
static auto random_value = Random::make_uniform_distribution(-15.0, +15.0);

TEST(directory_is_created) {
  bf::create_directory(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

static const float mass_smashon = 0.123;
static std::string mass_str = std::to_string(mass_smashon);
static std::string width_str = "1.200";
static std::string pdg_str = "661";
static std::string smashon_str = "smashon " + mass_str + " "
    + width_str + " " + pdg_str + "\n";
static const int zero = 0;

static ParticleData create_smashon_particle() {
  ParticleData particle = ParticleData{ParticleType::find(0x661)};
  particle.set_momentum(mass_smashon, random_value(), random_value(),
                        random_value());
  particle.set_position(FourVector(random_value(), random_value(),
                                   random_value(), random_value()));
  return particle;
}

static void compare_threevector(const std::array<std::string,3> &stringarray,
                               const ThreeVector &threevector) {
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(0).c_str()), threevector.x1(),
                         accuracy);
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(1).c_str()), threevector.x2(),
                         accuracy);
  COMPARE_ABSOLUTE_ERROR(std::atof(stringarray.at(2).c_str()), threevector.x3(),
                         accuracy);
}

TEST(outputfile) {
  /* Create particles */
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n" + smashon_str);

  Particles particles;
  const int number_of_particles = 5;
  for (int i = 0; i < number_of_particles; i++) {
    ParticleData particle = create_smashon_particle();
    particles.add_data(particle);
  }

  /* Create config file and object */
  std::string configfilename = "vtkconfig.yaml";
  std::ofstream cfgfile;
  cfgfile.open((testoutputpath / configfilename).native().c_str(),
               std::ios::out);
  cfgfile << "Options: None" << std::endl;
  cfgfile.close();
  Configuration&& op{testoutputpath, configfilename};
  /* Create output object */
  VtkOutput *vtkop = new VtkOutput(testoutputpath, std::move(op));
  int event_id = 0;
  /* Initial output */
  vtkop->at_eventstart(particles, event_id);
  std::string outputfilename = "pos_ev00000_tstep00000.vtk";
  VERIFY(bf::exists(testoutputpath / outputfilename));
  /* Time step output */
  Clock clock(0.0, 1.0);
  vtkop->at_intermediate_time(particles, event_id, clock);
  VERIFY(bf::exists(testoutputpath / "pos_ev00000_tstep00001.vtk"));

  std::fstream outputfile;
  outputfile.open((testoutputpath / outputfilename)
                  .native().c_str(), std::ios_base::in);
  if (outputfile.good()) {
    std::string line, item;
    /* Check header */
    std::string output_header = "";
    std::string header = "# vtk DataFile Version 2.0\n"
                         "Generated from molecular-offset data "
                         VERSION_MAJOR "\n"
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
    COMPARE(std::atoi(item.c_str()), particles.size());
    outputfile >> item;
    COMPARE(item, "double");
    for (int i = 0; i < number_of_particles; i++) {
      std::array<std::string, 3> position_string;
      for (int j = 0; j < 3; j++) {
        outputfile >> item;
        position_string[j] = item;
      }
      compare_threevector(position_string,
                          particles.data(i).position().threevec());
    }
    /* Check cell information */
    outputfile >> item;
    COMPARE(item, "CELLS");
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), particles.size());
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), particles.size() * 2);
    for (int i = 0; i < number_of_particles; i++) {
      outputfile >> item;
      COMPARE(std::atoi(item.c_str()), 1);
      outputfile >> item;
      COMPARE(std::atoi(item.c_str()), i);
    }
    outputfile >> item;
    COMPARE(item, "CELL_TYPES");
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), particles.size());
    for (int i = 0; i < number_of_particles; i++) {
      outputfile >> item;
      COMPARE(std::atoi(item.c_str()), 1);
    }
    /* Check point data */
    outputfile >> item;
    COMPARE(item, "POINT_DATA");
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), particles.size());
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
    /* Check momentum vectors */
    outputfile >> item;
    COMPARE(item, "VECTORS");
    outputfile >> item;
    COMPARE(item, "momentum");
    outputfile >> item;
    COMPARE(item, "double");
    for (int i = 0; i < number_of_particles; i++) {
      std::array<std::string, 3> momentum_string;
      for (int j = 0; j < 3; j++) {
        outputfile >> item;
        momentum_string[j] = item;
      }
      compare_threevector(momentum_string,
                          particles.data(i).momentum().threevec());
    }
  }
}
