/*
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_CUSTOMNUCLEUS_H_
#define SRC_INCLUDE_SMASH_CUSTOMNUCLEUS_H_

#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "nucleus.h"
#include "pdgcode.h"
#include "threevector.h"
namespace smash {

/**
 * Contains data for one nucleon that is read in from the list
 */
struct Nucleoncustom {
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
class CustomNucleus : public Nucleus {
 public:
  /**
   * Constructor that needs configuration parameters from input file
   * and the number of testparticles
   *
   * \param[in] config contains the parameters from the inputfile on the
   * numbers of particles with a certain PDG code and also the path where
   * the external particle list is located
   * \param[in] testparticles represents the number of testparticles
   * \param[in] same_file specifies if target and projectile nucleus are
   * read in from the same file, which is important for the ifstream
   */
  CustomNucleus(Configuration& config, int testparticles, bool same_file);
  /**
   * Fills Particlelist from vector containing data for one nucleus.
   * The data contains everything that is written in struct Nucleoncustom.
   *
   * \param[in] vec vector containing data from external list for one nucleus
   */
  void fill_from_list(const std::vector<Nucleoncustom>& vec);
  /// Returns position of a nucleon as given in the external file
  ThreeVector distribute_nucleon() override;
  /// Sets the positions of the nucleons inside a nucleus.
  void arrange_nucleons() override;
  /**
   * The returned vector contains Data for one nucleus given in the
   * particlelist.
   *
   * \param[in] infile is needed to read in from the external file
   * \param[in] particle_number ensures that only as many lines are read in
   * as the nucleus contains nucleons.
   */
  std::vector<Nucleoncustom> readfile(std::ifstream& infile,
                                      int particle_number) const;
  /**
   * Generates the name of the stream file.
   * \param[in] file_directory is the path to the external file
   * \param[in] file_name is the name of the external file
   */
  std::string file_path(const std::string& file_directory,
                        const std::string& file_name);
  /**
   * Generates Fermi momenta as it is done in the mother class but in addition
   * prints a warning that the Fermi momenta are generated accoriding to
   * Woods-Saxon distributed nucleons.
   */
  void generate_fermi_momenta() override;
  /**
   * Number of Nucleons per Nucleus
   * Set initally to zero to be modified in the constructor.
   * Is obtained by adding the proton and neutron numbers
   * specified in the config.yaml
   */
  int number_of_nucleons_ = 0;
  /// Vector contianing Data for one nucleus given in the particlelist
  std::vector<Nucleoncustom> custom_nucleus_;
  /// Index needed to read out vector in distribute nucleon
  size_t index = 0;

 private:
  /**
   * Filestream variable used if projectile and target are read in from the
   * same file and they use the same static stream.
   */
  /*
   * The unique_ptr is only required to work around a bug in GCC 4.8, because it
   * seems to be trying to use the non-existing copy-constructor of
   * `std::ifstream`. Newer compilers don't require this unneccessary
   * indirection.
   */
  static std::unique_ptr<std::ifstream> filestream_shared_;
  /**
   * Filestream variable used if projectile and target are read in from
   * different files and they therefore use different streams.
   */
  std::unique_ptr<std::ifstream> filestream_;
  /// Pointer to the used filestream pointer
  std::unique_ptr<std::ifstream>* used_filestream_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_CUSTOMNUCLEUS_H_
