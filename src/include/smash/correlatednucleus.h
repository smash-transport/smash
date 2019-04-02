/*
 *    Copyright (c) 2014-2019
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
  /// x-coordinate
  double x;
  /// y-coordinate
  double y;
  /// z-coordinate
  double z;
  /// spinprojection of the nucleon
  bool spinprojection;
  /// to differentiate between protons isospin=1 and neutrons isospin=0
  bool isospin;
};

/**
 * Inheriting from Nucleus-Class using modified Nucleon configurations.
 * Configurations are read in from external lists.
 */
class CorrelatedNucleus : public Nucleus {
 public:
  /**
   * Constructor that needs configuration parameters from input file
   * and the number of testparticles
   *
   * \param[in] config contains the parameters from the inputfile on the
   * numbers of particles with a certain PDG code and also the path where
   * the external particle list is located
   * \param[in] number of testparticles
   */
  CorrelatedNucleus(Configuration& config, int testparticles);
  /**
   * Fills Particlelist from vector containing data for one nucleus.
   * The data contains everything that is written in struct Nucleoncorr.
   *
   * \param[in] vector containing data from external list for one nucleus
   */
  void fill_from_list(const std::vector<Nucleoncorr>& vec);
  /// Returns position of a nucleon as given in the external file
  ThreeVector distribute_nucleon() override;
  /**
   * The returned vector contains Data for one nucleus given in the
   * particlelist.
   *
   * \param[in] infile is needed to read in from the external file
   * \param[in] particle_number ensures that only as many lines are read in
   * as the nucleus contains nucleons.
   */
  std::vector<Nucleoncorr> readfile(std::ifstream& infile,
                                    int particle_number) const;
  /**
   * Directory where the nucleon configurations are located.
   * Name is read in from manual input in the config.yaml
   */
  std::string particle_list_file_directory_;
  /**
   * File name of the nucleon configurations.
   * Name is read in from manual input in the config.yaml
   */
  std::string particle_list_file_name_;
  /**
   * Number of Nucleons per Nucleus
   * Set initally to zero to be modified in the constructor.
   * Number is obtained by adding the proton and neutron numbers
   * specified in the config.yaml
   */
  int number_of_nucleons_ = 0;
  /// Vector contianing Data for one nucleus given in the particlelist
  std::vector<Nucleoncorr> corr_nucleon_;
  /// Index needed to read out vector in distribute nucleon
  size_t index = 0;

 private:
  /**
   * Generates the name of the stream file.
   * \param[in] file_directory is the path to the external file
   * \param[in] file_name is the name of the external file
   */
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
