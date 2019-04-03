/*
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "smash/collidermodus.h"
#include "smash/customnucleus.h"
#include "smash/particledata.h"
#include "smash/particletype.h"
#include "smash/pdgcode.h"
#include "smash/pdgcode_constants.h"

namespace smash {

/*!\Userguide
 * \page projectile_and_target Projectile and Target
 * \n
 * Example: Configuring custom nuclei from external file
 * --------------
 * The following example illustrates how to configure a center-of-mass
 * heavy-ion collision with nuclei generated from an external file.
 * The nucleon positions are not sampled by smash but read in from an
 * external file. The given path and name
 * of the external file are made up and should be defined by the user
 * according to the used file.
 *\verbatim
 Modi:
     Collider:
         Projectile:
             Particles:    {2212: 79, 2112: 118}
             Custom:
                 File_Directory: "/home/username/custom_lists"
                 File_Name: "Au197_correlated.txt"
         Target:
             Particles:    {2212: 79, 2112: 118}
             Custom:
                 File_Directory: "/home/username/custom_lists"
                 File_Name: "Au197_correlated.txt"
         Sqrtsnn: 7.7
 \endverbatim
 *
 * The following example shows how an input file should be formatted:
 * <div class="fragment">
 * <div class="line"><span class="preprocessor"> 0.20100624   0.11402423
 * -2.40964466   0   0</span></div>
 * <div class="line"><span class="preprocessor"> 1.69072087  -3.21471918
 *  1.06050693   0   1</span></div>
 * <div class="line"><span class="preprocessor">-1.95791109  -3.51483782
 *  2.47294656   1   1</span></div>
 * <div class="line"><span class="preprocessor"> 0.43554894   4.35250733
 *  0.13331011   1   0</span></div>
 * </div>
 * The input file contains 5 columns (x, y, z, s, c). The first three
 * columns specify the spatial cordinates in fm. The fourth column denotes
 * the spin projection. The fifth contains the charge with 1 and 0 for
 * protons and neutrons respectively. In the example given the first line
 * defines a neutron and the second one a proton. Please make sure that
 * your file contains as many particles as you specified in the configuration.
 * For the example configuration your file needs to contain 79 protons and 
 * 118 neutrons in the first 197 lines. And the same number in the following
 * 197 lines. For every event you need to create 2 nuclei. It is crucial that
 * your file contains enough nucleons for the number of events you want to
 * simulate.
 *
 * \note
 * SMASH is shipped with an example configuration file to set up an collision 
 * with externally generated nucleon positions. This requires a particle list to
 * be read in. Both, the configuration file and the particle list, are located
 * in /input/custom_nucleus. To run SMASH with the provided example 
 * configuration and particle list, execute \n
 * \n
 * \verbatim
    ./smash -i INPUT_DIR/custom_nucleus/config.yaml
 \endverbatim
 * \n
 * Where 'INPUT_DIR' needs to be replaced by the path to the input directory
 * ('../input', if the build directory is located in the smash
 * folder).
 */

std::unique_ptr<std::ifstream> CustomNucleus::filestream_ = nullptr;
bool CustomNucleus::checkfileopen_ = false;

CustomNucleus::CustomNucleus(Configuration& config, int testparticles) {
  // Read in file directory from config
  std::string fd = config.take({"Custom", "File_Directory"});
  particle_list_file_directory_ = fd;
  // Read in file name from config
  std::string fn = config.take({"Custom", "File_Name"});
  particle_list_file_name_ = fn;

  if (particles_.size() != 0) {
    throw std::runtime_error(
        "Your Particle List is already filled before reading in from the "
        "external file."
        "Something went wrong. Please check your config.");
  }

  /*
   * Counts number of nucleons in one nucleus as it is specialized
   * by the user in the config file.
   * It is needed to read in the proper number of nucleons for one
   * nucleus and to restart at the listreading for the following
   * nucleus as one does not want to read configurations twice. */
  std::map<PdgCode, int> particle_list = config.take({"Particles"});
  for (const auto& particle : particle_list)
    number_of_nucleons_ += particle.second * testparticles;
  if (!checkfileopen_) {
    filestream_ = make_unique<std::ifstream>(
        streamfile(particle_list_file_directory_, particle_list_file_name_));
    // tracked that file was opened once
    checkfileopen_ = true;
  }
  custom_nucleon_ = readfile(*filestream_, number_of_nucleons_);

  fill_from_list(custom_nucleon_);
  // Inherited from nucleus class (see nucleus.h)
  set_parameters_automatic();
}

void CustomNucleus::fill_from_list(const std::vector<Nucleoncustom>& vec) {
  particles_.clear();
  index = 0;
  // checking if particle is proton or neutron
  for (const auto& it : vec) {
    PdgCode pdgcode;
    if (it.isospin == 1) {
      pdgcode = pdg::p;
    } else if (it.isospin == 0) {
      pdgcode = pdg::n;
    } else {
      throw std::runtime_error(
          "Your particles charges are not 1 = proton or 0 = neutron."
          "Check whether your list is correct or there is an error.");
    }
    // setting parameters for the particles in the particlelist in smash
    const ParticleType& current_type = ParticleType::find(pdgcode);
    double current_mass = current_type.mass();
    particles_.emplace_back(current_type);
    particles_.back().set_4momentum(current_mass, 0.0, 0.0, 0.0);
  }
}

ThreeVector CustomNucleus::distribute_nucleon() {
  if (index >= custom_nucleon_.size()) {
    custom_nucleon_ = readfile(*filestream_, number_of_nucleons_);

    fill_from_list(custom_nucleon_);
  }
  const auto& pos = custom_nucleon_.at(index);
  index++;
  return ThreeVector(pos.x, pos.y, pos.z);
}

std::string CustomNucleus::streamfile(const std::string& file_directory,
                                          const std::string& file_name) {
  if (file_directory.back() == '/') {
    return file_directory + file_name;
  } else {
    return file_directory + '/' + file_name;
  }
}

std::vector<Nucleoncustom> CustomNucleus::readfile(
    std::ifstream& infile, int particle_number) const {
  int A = particle_number;
  std::string line;
  std::vector<Nucleoncustom> nucleon;
  int i = 0;
  while (std::getline(infile, line)) {
    Nucleoncustom nuc;
    std::istringstream iss(line);
    if (!(iss >> nuc.x >> nuc.y >> nuc.z >> nuc.spinprojection >>
          nuc.isospin)) {
      throw std::runtime_error(
          "SMASH could not read in a line from your initial nuclei input file."
          "Check if your file has follwing format: x y z spinprojection "
          "isospin");
      break;
    }
    nucleon.push_back(nuc);
    // ensuring that only A particles are read in for one nucleus
    if (++i == A) {
      break;

}  // namespace smash
