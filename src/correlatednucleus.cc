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
#include "smash/correlatednucleus.h"
#include "smash/particledata.h"
#include "smash/particletype.h"
#include "smash/pdgcode.h"
#include "smash/pdgcode_constants.h"

namespace smash {

/*!\Userguide
 * \page projectile_and_target Projectile and Target
 *
 * \li \key File_Directory:
 * The directory where the external list with the nucleon configurations
 * is located. Please make sure to use an absolute path.
 *
 * \li \key File_Name:
 * The file name of the external list with the nucleon configurations.
 *
 * \n
 * Example: Configuring a correlated nucleus
 * --------------
 * The following example illustrates how to configure a COM 
 * heavy-ion collision with correlated nuclei. The given path and name 
 * of the external file are made up and should be defined by the user
 * according to the used file.
 *\verbatim
 Modi:
     Collider:
         Projectile:
             Particles:    {2212: 79, 2112: 118}
             Correlated:
                 File_Directory: /home/username/external_lists
                 File_Name: Au197_correlated.txt
         Target:
             Particles:    {2212: 79, 2112: 118}
             Correlated:
                 File_Directory: /home/username/external_lists
                 File_Name: Au197_correlated.txt
         Sqrtsnn: 7.7
 \endverbatim
 */

std::unique_ptr<std::ifstream> CorrelatedNucleus::filestream_ = nullptr;
bool CorrelatedNucleus::checkfileopen_ = false;

CorrelatedNucleus::CorrelatedNucleus(Configuration& config, int testparticles) {
  // Read in file directory from config
  std::string fd = config.take({"Correlated", "File_Directory"});
  particle_list_file_directory_ = fd;
  // Read in file name from config
  std::string fn = config.take({"Correlated", "File_Name"});
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
  corr_nucleon_ = readfile(*filestream_, number_of_nucleons_);

  fill_from_list(corr_nucleon_);
  // Inherited from nucleus class (see nucleus.h)
  set_parameters_automatic();
}

void CorrelatedNucleus::fill_from_list(const std::vector<Nucleoncorr>& vec) {
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

ThreeVector CorrelatedNucleus::distribute_nucleon() {
  if (index >= corr_nucleon_.size()) {
    corr_nucleon_ = readfile(*filestream_, number_of_nucleons_);

    fill_from_list(corr_nucleon_);
  }
  const auto& pos = corr_nucleon_.at(index);
  index++;
  return ThreeVector(pos.x, pos.y, pos.z);
}

std::string CorrelatedNucleus::streamfile(const std::string& file_directory,
                                          const std::string& file_name) {
  if (file_directory.back() == '/') {
    return file_directory + file_name;
  } else {
    return file_directory + '/' + file_name;
  }
}

std::vector<Nucleoncorr> CorrelatedNucleus::readfile(
    std::ifstream& infile, int particle_number) const {
  int A = particle_number;
  std::string line;
  std::vector<Nucleoncorr> nucleon;
  int i = 0;
  while (std::getline(infile, line)) {
    Nucleoncorr nuc;
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
    }
  }
  return nucleon;
}

}  // namespace smash
