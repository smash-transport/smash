/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ROOTOUTPUT_H_
#define SRC_INCLUDE_ROOTOUTPUT_H_

#include <boost/filesystem.hpp>
#include <vector>
#include "outputinterface.h"
#include "TFile.h"
#include "TTree.h"

namespace Smash {
class Particles;

/**
  * \brief SMASH output to ROOT file
  * ----------------------------------------------------------------------------
  * SMASH supports ROOT (root.cern.ch) output as an option. By default SMASH
  * does not need ROOT for compilation. If one wants to produce ROOT output, one
  * has to compile SMASH as follows:
  * \code
  * cd build
  * cmake -D USE_ROOT=ON ..
  * make
  * \endcode
  * This class produces file smash_run.root, which contains
  * ROOT TTrees. Each TTree contains information about particles
  * at some moment of simulation from all SMASH events.
  * In addition there are TTrees with names at_eventstart and at_eventend,
  * which store information at event start (after initialization) and at 
  * event end (after the final moment of time). 
  * Particle information is stored in TBranches.
  * For each particle characteristic there is a separate branch.
  * Currently these are x,y,z (coordinates), p0,px,py,pz (4-momentum), 
  * id - unique identifier for every SMASH particle during one event,
  * pdgid - PDG code of particle, that characterizes its sort and
  * ev - number of event particle encountered in.
  * 
  * Here is an example of ROOT macro to read the ROOT output of SMASH:
  * \code
  * int rootfile_analysis_example(void) {
  *   // open SMASH output file to be read in
  *   TFile *input_file = TFile::Open("../build/data/0/smash_run.root");
  *   if (input_file->IsOpen()) {
  *     printf("Succesfully opened file %s\n", input_file->GetName());
  *   } else {
  *     printf("Error at opening file %s\n", input_file->GetName());
  *   }
  * 
  *   // Get a tree from file
  *   TTree *ev_end_tree = static_cast<TTree*>(input_file->Get("at_eventend"));
  * 
  *   // Get number of entries in a tree
  *   Int_t nentries = ev_end_tree->GetEntries();
  *   printf("Number of entries in a tree is %d\n", nentries);
  *
  *   // This draws p_T distribution 
  *   // ev_end_tree->Draw("sqrt(px*px + py*py)");

  *   // This draws 3D momentum space distribution
  *   ev_end_tree->Draw("px:py:pz");
  * 
  *   return 0;
  * }
  * \endcode
  * For examples of extracting info from .root file see root.cern.ch
  **/
class RootOutput : public OutputInterface {
 public:
  explicit RootOutput(boost::filesystem::path path);
  ~RootOutput();

  void at_eventstart(const Particles &particles,
                     const int event_number) override;
  void at_eventend(const Particles &particles,
                   const int event_number) override;
  void after_Nth_timestep(const Particles &particles,
                          const int event_number,
                          const Clock &) override;
  void write_interaction(const ParticleList &incoming_particles,
                         const ParticleList &outgoing_particles) override;

 private:
  const boost::filesystem::path base_path_;
  std::unique_ptr<TFile> root_out_file_;
  // TFile takes ownership of all TTrees.
  // That's why TTrees are not unique pointers.
  std::vector<TTree*> tree_list_;
  void particles_to_tree(const char* treename,
                         const char* treedescr,
                         const Particles &particles,
                         const int event_number);
};
}  // namespace Smash

#endif  // SRC_INCLUDE_ROOTOUTPUT_H_
