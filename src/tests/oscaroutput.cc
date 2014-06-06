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

using namespace Smash;

const bf::path testoutputpath = bf::absolute(SMASH_TEST_OUTPUT_PATH);

TEST(directory_is_created) {
  bf::create_directory(testoutputpath);
  VERIFY(bf::exists(testoutputpath));
}

static ParticleData create_smashon_particle(int id = -1) {
  return ParticleData{ParticleType::find(-0x331), id};
}

TEST(output_format) {
  OscarFullHistoryOutput *oscfull = new OscarFullHistoryOutput(testoutputpath);
  VERIFY(bf::exists(testoutputpath / "full_event_history.oscar"));

  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PDG\n"
      "smashon 0.123 1.2 -331\n");

  ParticleData particle = create_smashon_particle(0);
  particle.set_momentum(0.123, 1.1, 2.2, 3.3);
  particle.set_position(FourVector(7.7, 4.4, 5.5, 6.6));
  Particles particles({});
  particles.add_data(particle);
  int id = 0;
  oscfull->at_eventstart(particles, id);

  delete oscfull;

  std::fstream outputfile;
  outputfile.open((testoutputpath / "full_event_history.oscar")
                  .native().c_str(),
         std::ios_base::in);
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
    /* Check interaction block description line  */
    std::getline(outputfile, line);
    COMPARE(line, "0 1 1");
    /* Check particle data line item by item */
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), id);
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), -331);
    outputfile >> item;
    COMPARE(std::atoi(item.c_str()), 0);
    outputfile >> item;
    COMPARE(std::atof(item.c_str()), particle.momentum().x1());
    outputfile >> item;
    COMPARE(std::atof(item.c_str()), particle.momentum().x2());
    outputfile >> item;
    COMPARE(std::atof(item.c_str()), particle.momentum().x3());
    outputfile >> item;
    COMPARE(std::atof(item.c_str()), particle.momentum().x0());
    outputfile >> item;
    COMPARE(std::atof(item.c_str()), 0.123);
    outputfile >> item;
    COMPARE(std::atof(item.c_str()), particle.position().x1());
    outputfile >> item;
    COMPARE(std::atof(item.c_str()), particle.momentum().x2());
    outputfile >> item;
    COMPARE(std::atof(item.c_str()), particle.momentum().x3());
    outputfile >> item;
    COMPARE(std::atof(item.c_str()), particle.momentum().x0());
  }
}
