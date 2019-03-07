/*
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_CORRELATEDNUCLEUS_H_
#define SRC_INCLUDE_CORRELATEDNUCLEUS_H_

#include <fstream>
#include <map>
#include <vector>

#include "collidermodus.h"
#include "nucleus.h"
#include "pdgcode.h"
#include "threevector.h"
namespace smash {

/**
 * Contains data for one nucleon that is read in from the list
 */
struct Nucleoncorr {
  double x;  // x-coordinate
  double y;  // y-coordinate
  double z;  // z-coordinate
  bool spinprojection;
  bool isospin;
};

/**
 * Inheriting from Nucleus-Class using modified Nucleon configurations.
 * Configurations are read in from external lists.
 */
class CorrelatedNucleus : public Nucleus {
 public:
  /// Constructor calling Collidermodus::readfile
  CorrelatedNucleus(Configuration& config, int testparticles);
  /// fill Particlelist
  void fill_from_list(const std::vector<Nucleoncorr>& vec);
  /// returns position of a nucleus
  ThreeVector distribute_nucleon() override;
  /**
   * Contains Data for one nucleus given in the particlelist
   */
  std::vector<Nucleoncorr> readfile(std::ifstream& infile,
                                    int particle_number) const;
  /** Directory where the nucleon configurations are located.
   *  Name is read in from manual input in the config.yaml
   */
  std::string particle_list_file_directory_;
  /** File name of the nucleon configurations.
   *  Name is read in from manual input in the config.yaml
   */
  std::string particle_list_file_name_;
  /** Number of Nucleons per Nucleus
   *  Set initally to zero to be modified in the constructor.
   *  Number is obtained by adding the proton and neutron numbers
   *  specified in the config.yaml
   */
  int number_of_nucleons_ = 0;
  /// Vector contianing Data for one nucleus given in the particlelist
  std::vector<Nucleoncorr> corr_nucleon_;
  /// Index needed to read out vector in distribute nucleon
  size_t index = 0;

 private:
  /// Generate the name of the stream file.
  static std::string streamfile(const std::string& file_directory,
                                const std::string& file_name);
  /// Variable carrying the output of the streamfile function
  /*
   * The unique_ptr is only required to work around a bug in GCC 4.8, because it
   * seems to be trying to use the non-existing copy-constructor of
   * `std::ifstream`. Newer compilers don't require this unneccessary
   * indirection.
   */
  static std::unique_ptr<std::ifstream> filestream_;
  /** Bool variable to check if the file was already opened. It ensures
   *  to read in every nucleus configuration given only once. If the bool
   *  is true the constructor uses the stream that was given the last time
   *  the constructor was called.
   */
  static bool checkfileopen_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_CORRELATEDNUCLEUS_H_
