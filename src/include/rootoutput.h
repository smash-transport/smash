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
  * This class produces file smash_run.root, which contains a
  * ROOT TTree. TTree contains information about particles
  * during simulation from all SMASH events.
  * Output is happening in blocks. All particles in a block
  * are at the same time and in the same event. However, it is possible that
  * different blocks are at the same time and from the same event.
  * Particle information is stored in TBranches.
  * For each particle characteristic there is a separate branch.
  * Currently these are t,x,y,z (coordinates), p0,px,py,pz (4-momentum), 
  * pdgid - PDG code of particle, that characterizes its sort,
  * ev - number of event particle encountered in and
  * tcounter - number of output block in a given event.
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
  *   TTree *tree = static_cast<TTree*>(input_file->Get("particles"));
  * 
  *   // Get number of entries in a tree
  *   Int_t nentries = tree->GetEntries();
  *   printf("Number of entries in a tree is %d\n", nentries);
  *
  *   // This draws p_T distribution at initialization
  *   // tree->Draw("sqrt(px*px + py*py)","tcounter==0");

  *   // This draws 3D momentum space distribution at initialization
  *   tree->Draw("px:py:pz","tcounter==0");
  * 
  *   return 0;
  * }
  * \endcode
  * For examples of extracting info from .root file see root.cern.ch
  * To view ROOT file use TBrowser:
  * \code
  * root -l
  * new TBrowser
  * \endcode
  * If option write_collisions is set True, then in addition to particles 
  * TTree a collision TTree is created. Information about each collision is
  * written as one leaf: nin, nout - number of incoming and outgoing particles,
  * ev - event number, (t,x,y,z), (p0,px,py,pz) - arrays of dimension nin+nout
  * that contain coordinates and momenta.
  **/
class RootOutput : public OutputInterface {
 public:
  explicit RootOutput(boost::filesystem::path path, Options op);
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
  // That's why TTree is not a unique pointer.
  TTree* particles_tree_;
  TTree* collisions_tree_;
  void particles_to_tree(const Particles &particles,
                         const int event_number);
  void collisions_to_tree(const ParticleList &incoming,
                          const ParticleList &outgoing);
  // Counts number of output in a given event
  int output_counter_;
  int current_event_;

  static const int max_buffer_size_ = 10000;
  // Variables that serve as buffer for filling TTree
  std::array<double, max_buffer_size_> p0, px, py, pz, t, x, y, z;
  std::array<int, max_buffer_size_>    pdgcode;
  int npart, tcounter, ev, nin, nout;

  // Option to write collisions
  bool write_collisions_;
};
}  // namespace Smash

#endif  // SRC_INCLUDE_ROOTOUTPUT_H_
