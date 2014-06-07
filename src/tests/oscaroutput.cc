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
#include <iostream>
#include <fstream>
#include <string>
#include <boost/filesystem.hpp>
#include "../include/outputinterface.h"
#include "../include/oscarfullhistoryoutput.h"
#include "../include/oscarparticlelistoutput.h"
#include "../include/particles.h"
#include "../include/random.h"

using namespace Smash;

static const float accuracy = 1.0e-4;
static const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);
static auto random_value = Random::make_uniform_distribution(-15.0, +15.0);

TEST(directory_is_created) {
  bf::create_directory(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

static ParticleData create_smashon_particle() {
  return ParticleData{ParticleType::find(-0x331)};
}

static const float mass_smashon = 0.123;
static std::string mass_str = std::to_string(mass_smashon);
static std::string width_str = "1.200";
static std::string pdg_str = "-331";
static std::string smashon_str = "smashon " + mass_str + " "
    + width_str + " " + pdg_str + "\n";
static const int zero = 0;

TEST(fullhistory_format) {
  OscarFullHistoryOutput *oscfull = new OscarFullHistoryOutput(testoutputpath);
  VERIFY(bf::exists(testoutputpath / "full_event_history.oscar"));

  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n" + smashon_str);

  Particles particles({});

  ParticleData particle = create_smashon_particle();
  particle.set_momentum(mass_smashon, random_value(), random_value(),
                        random_value());
  particle.set_position(FourVector(random_value(), random_value(),
                                   random_value(), random_value()));
  particles.add_data(particle);

  int event_id = 0;
  oscfull->at_eventstart(particles, event_id);

  delete oscfull;

  std::fstream outputfile;
  outputfile.open((testoutputpath / "full_event_history.oscar")
                  .native().c_str(), std::ios_base::in);
  if (outputfile.good()) {
    std::string line, item;
    std::getline(outputfile, line);
    /* Check header */
    COMPARE(line, "# OSC1999A");
    std::getline(outputfile, line);
    COMPARE(line, "# full_event_history");
    std::getline(outputfile, line);
    COMPARE(line, "# smash");
    std::getline(outputfile, line);
    COMPARE(line, "# Block format:");
    std::getline(outputfile, line);
    COMPARE(line, "# nin nout event_number");
    std::getline(outputfile, line);
    COMPARE(line, "# id pdg 0 px py pz p0 mass x y z t");
    std::getline(outputfile, line);
    COMPARE(line, "# End of event: 0 0 event_number");
    std::getline(outputfile, line);
    COMPARE(line, "#");
    /* Check interaction block description line item by item */
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), 0);
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), particles.size());
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), event_id + 1);
    /* Check particle data line item by item */
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), 0);
    outputfile >> item;
    COMPARE(item, pdg_str);
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), zero);
    outputfile >> item;
    COMPARE_ABSOLUTE_ERROR(std::atof(item.c_str()), particle.momentum().x1(),
                           accuracy);
    outputfile >> item;
    COMPARE_ABSOLUTE_ERROR(std::atof(item.c_str()), particle.momentum().x2(),
                           accuracy);
    outputfile >> item;
    COMPARE_ABSOLUTE_ERROR(std::atof(item.c_str()), particle.momentum().x3(),
                           accuracy);
    outputfile >> item;
    COMPARE_ABSOLUTE_ERROR(std::atof(item.c_str()), particle.momentum().x0(),
                           accuracy);
    outputfile >> item;
    COMPARE(float(std::atof(item.c_str())), mass_smashon);
    outputfile >> item;
    COMPARE_ABSOLUTE_ERROR(std::atof(item.c_str()), particle.position().x1(),
                           accuracy);
    outputfile >> item;
    COMPARE_ABSOLUTE_ERROR(std::atof(item.c_str()), particle.position().x2(),
                           accuracy);
    outputfile >> item;
    COMPARE_ABSOLUTE_ERROR(std::atof(item.c_str()), particle.position().x3(),
                           accuracy);
    outputfile >> item;
    COMPARE_ABSOLUTE_ERROR(std::atof(item.c_str()), particle.position().x0(),
                           accuracy);
  }
}
