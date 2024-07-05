/*
 *    Copyright (c) 2019-2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#include "smash/customnucleus.h"

#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "smash/constants.h"
#include "smash/input_keys.h"
#include "smash/logging.h"
#include "smash/particletype.h"
#include "smash/pdgcode.h"

namespace smash {
static constexpr int LCollider = LogArea::Collider::id;

std::unique_ptr<std::ifstream> CustomNucleus::filestream_shared_ = nullptr;

CustomNucleus::CustomNucleus(Configuration& config, int testparticles,
                             bool same_file) {
  const bool is_projectile = is_configuration_about_projectile(config);
  const auto& [file_dir_key, filename_key, particles_key] = [&is_projectile]() {
    return is_projectile
               ? std::make_tuple(
                     InputKeys::modi_collider_projectile_custom_fileDirectory,
                     InputKeys::modi_collider_projectile_custom_fileName,
                     InputKeys::modi_collider_projectile_particles)
               : std::make_tuple(
                     InputKeys::modi_collider_target_custom_fileDirectory,
                     InputKeys::modi_collider_target_custom_fileName,
                     InputKeys::modi_collider_target_particles);
  }();
  // Read in file directory from config
  const std::string particle_list_file_directory = config.take(file_dir_key);
  // Read in file name from config
  const std::string particle_list_file_name = config.take(filename_key);

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
   * nucleus as one does not want to read configurations twice.
   */
  std::map<PdgCode, int> particle_list = config.take(particles_key);
  for (const auto& particle : particle_list) {
    if (particle.first == pdg::p) {
      number_of_protons_ = particle.second * testparticles;
    } else if (particle.first == pdg::n) {
      number_of_neutrons_ = particle.second * testparticles;
    } else {
      throw std::runtime_error(
          "Your nucleus can only contain protons and/or neutrons."
          "Please check what particles you have specified in the config");
    }
    number_of_nucleons_ = number_of_protons_ + number_of_neutrons_;
  }
  /*
   * "if" statement makes sure the streams to the file are initialized
   * properly.
   */
  const std::string path =
      file_path(particle_list_file_directory, particle_list_file_name);
  if (same_file && !filestream_shared_) {
    filestream_shared_ = std::make_unique<std::ifstream>(path);
    used_filestream_ = &filestream_shared_;
  } else if (!same_file) {
    filestream_ = std::make_unique<std::ifstream>(path);
    used_filestream_ = &filestream_;
  } else {
    used_filestream_ = &filestream_shared_;
  }

  custom_nucleus_ = readfile(**used_filestream_);
  fill_from_list(custom_nucleus_);
  // Inherited from nucleus class (see nucleus.h)
  set_parameters_automatic();
}

void CustomNucleus::fill_from_list(const std::vector<Nucleoncustom>& vec) {
  particles_.clear();
  index_ = 0;
  // checking if particle is proton or neutron
  for (const auto& it : vec) {
    PdgCode pdgcode;
    if (it.isospin == 1) {
      pdgcode = pdg::p;
    } else if (it.isospin == 0) {
      pdgcode = pdg::n;
    } else {
      throw std::runtime_error(
          "Your particles charges are not 1 = proton or 0 = neutron.\n"
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
  /*
   * As only arrange_nucleons is called at the beginning of every
   * event it is important to have readfile and fill from list
   * called again when a new event starts. The constructor is only
   * called twice to initialize the first target and projectile.
   * Therefore this if statement is implemented.
   */
  if (index_ >= custom_nucleus_.size()) {
    custom_nucleus_ = readfile(**used_filestream_);
    fill_from_list(custom_nucleus_);
  }
  const auto& pos = custom_nucleus_.at(index_);
  index_++;
  ThreeVector nucleon_position(pos.x, pos.y, pos.z);
  // rotate nucleon about euler angle
  nucleon_position.rotate(euler_phi_, euler_theta_, euler_psi_);

  return nucleon_position;
}

void CustomNucleus::arrange_nucleons() {
  /* Randomly generate Euler angles for rotation everytime a new
   * custom nucleus is initialized. Therefore this is done 2 times per
   * event.
   */
  Nucleus::random_euler_angles();

  for (auto i = begin(); i != end(); i++) {
    // Initialize momentum
    i->set_4momentum(i->pole_mass(), 0.0, 0.0, 0.0);
    /* Sampling the Woods-Saxon, get the radial
     * position and solid angle for the nucleon. */
    ThreeVector pos = distribute_nucleon();
    // Set the position of the nucleon.
    i->set_4position(FourVector(0.0, pos));
  }
  // Recenter
  align_center();
}

void CustomNucleus::generate_fermi_momenta() {
  Nucleus::generate_fermi_momenta();
  logg[LCollider].warn() << "Fermi motion activated with a custom nucleus.\n";
  logg[LCollider].warn() << "Be aware that generating the Fermi momenta\n"
                         << "assumes nucleons distributed according to a\n"
                         << "Woods-Saxon distribution.";
}

std::string CustomNucleus::file_path(const std::string& file_directory,
                                     const std::string& file_name) {
  if (file_directory.back() == '/') {
    return file_directory + file_name;
  } else {
    return file_directory + '/' + file_name;
  }
}

std::vector<Nucleoncustom> CustomNucleus::readfile(
    std::ifstream& infile) const {
  int proton_counter = 0;
  int neutron_counter = 0;
  std::string line;
  std::vector<Nucleoncustom> custom_nucleus;
  // read in only A particles for one nucleus
  for (int i = 0; i < number_of_nucleons_; ++i) {
    std::getline(infile, line);
    // make sure the stream goes back to the beginning when it hits end of file
    if (infile.eof()) {
      infile.clear();
      infile.seekg(0, infile.beg);
      std::getline(infile, line);
    }
    Nucleoncustom nucleon;
    std::istringstream iss(line);
    if (!(iss >> nucleon.x >> nucleon.y >> nucleon.z >>
          nucleon.spinprojection >> nucleon.isospin)) {
      throw std::runtime_error(
          "SMASH could not read in a line from your initial nuclei input file."
          "\nCheck if your file has the following format: x y z "
          "spinprojection isospin");
    }
    if (nucleon.isospin == 1) {
      proton_counter++;
    } else if (nucleon.isospin == 0) {
      neutron_counter++;
    }
    custom_nucleus.push_back(nucleon);
  }
  if (proton_counter != number_of_protons_ ||
      neutron_counter != number_of_neutrons_) {
    throw std::runtime_error(
        "Number of protons and/or neutrons in the nuclei input file does not "
        "correspond to the number specified in the config.\nCheck the config "
        "and your input file.");
  } else {
    return custom_nucleus;
  }
}

}  // namespace smash
