/*
 *
 *    Copyright (c) 2014-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_ROOTOUTPUT_H_
#define SRC_INCLUDE_ROOTOUTPUT_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include "TFile.h"
#include "TTree.h"

#include "configuration.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "outputparameters.h"

namespace smash {
class Particles;

/**
 * \ingroup output
 *
 * \brief SMASH output to ROOT file
 * ----------------------------------------------------------------------------
 * SMASH supports ROOT output as an option (see http://root.cern.ch).
 * The ROOT framework needs to be installed when building SMASH, otherwise
 * ROOT support will be disabled.
 *
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
 * charge - electric charge of the particle,
 * ev - number of event particle encountered in,
 * tcounter - number of output block in a given event,
 * npart - number of particles,
 * impact_b - impact parameter, and
 * empty_event - whether there was no interaction between the projectile and
 * the target.
 *
 * Here is an example of a basic ROOT macro to read the ROOT output of SMASH:
 * \code
 * // file name: basic_macro.C
 *
 * #include <TFile.h>
 * #include <TTree.h>
 *
 * int rootfile_basic_example() {
 *   // open SMASH output file to be read in
 *   TFile *input_file = TFile::Open("../build/data/0/Particles.root");
 *   if (input_file->IsOpen()) {
 *     printf("Successfully opened file %s\n", input_file->GetName());
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
 *
 *   // This draws 3D momentum space distribution at initialization
 *   tree->Draw("px:py:pz","tcounter==0");
 *
 *   return 0;
 * }
 * \endcode
 * To execute this macro:
 * \code
 * root -l
 * .L basic_macro.C+
 * rootfile_basic_example()
 * \endcode
 *
 * To generate a ROOT macro for a more involved analysis:
 * \code
 * root -l Particles.root
 * particles->MakeClass("analysis")
 * \endcode
 * This creates analysis.h and analysis.C, the latter of which provides an
 * example of a basic loop over entries. The user can modify the Loop() function
 * and build a more complicated analysis. The following is an example of a macro
 * using an array of histograms to obtain the radial momentum distribution for
 * each of the entries in the Particles.root:
 * \code
 * #define analysis_cxx
 * #include "analysis.h"
 * #include <TH2.h>
 * #include <TStyle.h>
 * #include <TCanvas.h>
 * #include <iostream>
 *
 * void analysis::Loop()
 * {
 *    if (fChain == 0) return;
 *    Long64_t n_entries = fChain->GetEntriesFast();
 *
 *    // an array of histograms
 *    TH1D *h_p_avg[n_entries];
 *    // Each histogram needs to be declared with a unique name and title
 *    for (int i = 0; i < n_entries; i++){
 *      char h_p_avg_name[256];
 *      char h_p_avg_title[256];
 *
 *      sprintf(h_p_avg_name, "h_p_avg_entry_%d", i);
 *      sprintf(h_p_avg_title, "momentum distribution at entry %d", i);
 *      h_p_avg[i] = new TH1D(h_p_avg_name, h_p_avg_title, 50, 0, 1.0);
 *      h_p_avg[i]->Sumw2();
 *      h_p_avg[i]->SetStats(1);
 *    }
 *
 *    Long64_t nb = 0;
 *    // A loop over all entries
 *    for (Long64_t j_entry = 0; j_entry < n_entries; j_entry++){
 *      //Load the TTree data for that entry
 *      Long64_t i_entry = LoadTree(j_entry);
 *      if (i_entry < 0){
 *        std::cout << "Failed to load the TTree at j_entry = "
 * 		    << j_entry << std::endl;
 *        break;
 *      }
 *      nb = fChain->GetEntry(j_entry);
 *
 *      // A loop over all particles in the entry
 *      for (int i = 0; i < npart; i++) {
 *        const double p = sqrt( px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]
 * );
 *        // filling the j_entry-th histogram
 *        h_p_avg[j_entry]->Fill(p, 1.0/(p*p));
 *      }
 *    }
 *
 *    // drawing the histogram corresponding to the j_entry = 0 entry
 *    h_p_avg[0]->Draw();
 * }
 * /endcode
 * To run this analysis:
 * \code
 * root -l
 * .L analysis.C+
 * analysis t
 * t.Loop()
 * \endcode
 *
 * To quickly view a ROOT file from the command line:
 * \code
 * root -l Particles.root  // attaches the .root file
 * .ls  // lists objects contained in the .root file; here: particles
 * particles->Draw("p0", "tcounter == 0")
 * \endcode
 *
 * Viewing a ROOT file can be also done in a TBrowser:
 * \code
 * root -l
 * new TBrowser
 * \endcode
 *
 * For more examples of extracting info from .root file see root.cern.ch
 *
 * If option write_collisions is set True, then in addition to particles
 * TTree a collision TTree is created. Information about each collision is
 * written as one leaf: nin, nout - number of incoming and outgoing particles,
 * ev - event number, weight - total weight of the collision (wgt),
 * partial_weight - partial weight of the collision (par_wgt), (t,x,y,z),
 * (p0,px,py,pz) - arrays of dimension nin+nout
 * that contain coordinates and momenta.
 */
class RootOutput : public OutputInterface {
 public:
  /**
   * Construct ROOT output.
   *
   * \param[in] path Output path.
   * \param[in] name Name of the ouput.
   * \param[in] out_par A structure containing parameters of the output.
   */
  RootOutput(const bf::path &path, const std::string &name,
             const OutputParameters &out_par);

  /// Destructor
  ~RootOutput();

  /**
   * update event number and writes intermediate particles to a tree.
   * \param[in] particles Particles to be written to output.
   * \param[in] event_number event number to be used in ROOT output.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventstart(const Particles &particles, const int event_number,
                     const EventInfo &event) override;
  /**
   * update event number and impact parameter,
   * and writes intermediate particles to a tree.
   * \param[in] particles Particles to be written to output.
   * \param[in] event_number event number to be used in ROOT output.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventend(const Particles &particles, const int event_number,
                   const EventInfo &event) override;
  /**
   * Writes intermediate particles to a tree defined by treename,
   * if it is allowed (i.e., particles_only_final_ is No).
   * \param[in] particles Particles to be written to output.
   * \param[in] clock Unused, needed since inherited.
   * \param[in] dens_param Unused, needed since inherited.
   * \param[in] event Event info, see \ref event_info
   */
  void at_intermediate_time(const Particles &particles,
                            const std::unique_ptr<Clock> &clock,
                            const DensityParameters &dens_param,
                            const EventInfo &event) override;
  /**
   * Writes collisions to a tree defined by treename.
   * \param[in] action an Action object containing incoming, outgoing particles
   *            and type of interactions.
   * \param[in] density Unused, needed since inherited.
   */
  void at_interaction(const Action &action, const double density) override;

 private:
  /// Filename of output
  const bf::path filename_;
  /// Filename of output as long as simulation is still running.
  bf::path filename_unfinished_;
  /// Pointer to root output file.
  std::unique_ptr<TFile> root_out_file_;
  /**
   * TTree for particles output.
   *
   * TFile takes ownership of all TTrees.
   * That's why TTree is not a unique pointer.
   */
  TTree *particles_tree_;
  /**
   * TTree for collision output.
   *
   * TFile takes ownership of all TTrees.
   * That's why TTree is not a unique pointer.
   */
  TTree *collisions_tree_;
  /**
   * Writes particles to a tree defined by treename.
   * \param[in] particles Particles or ParticleList to be written to output.
   */
  template <typename T>
  void particles_to_tree(T &particles);
  /**
   * Writes collisions to a tree defined by treename.
   * \param[in] incoming Incoming particles to be written to output.
   * \param[in] outgoing Outgoing particles to be written to output.
   * \param[in] weight Total weight of the collision.
   * \param[in] partial_weight Partial weight of the collision
   */
  void collisions_to_tree(const ParticleList &incoming,
                          const ParticleList &outgoing, const double weight,
                          const double partial_weight);
  /// Number of output in a given event.
  int output_counter_ = 0;
  /// Number of current event.
  int current_event_ = 0;

  /*
   * Maximal buffer size.
   * When the number of particles N exceeds the buffer size B, data is flushed
   * to the ROOT file every B particles. This creates ceil(N/B) entries in the
   * ROOT Tree at every output.
   */
  static const int max_buffer_size_ = 500000;

  /** @name Buffer for filling TTree
   * See class documentation for definitions.
   */
  //@{
  /// Property that is written to ROOT output.
  std::vector<double> p0_, px_, py_, pz_, t_, x_, y_, z_, formation_time_,
      xsec_factor_, time_last_coll_;
  std::vector<int> pdgcode_, charge_, coll_per_part_, proc_id_origin_,
      proc_type_origin_, pdg_mother1_, pdg_mother2_;
  int npart_, tcounter_, ev_, nin_, nout_, test_p_;
  double wgt_, par_wgt_, impact_b_, modus_l_, current_t_;
  double E_kinetic_tot_, E_fields_tot_, E_tot_;
  bool empty_event_;
  //@}

  /// Option to write collisions tree.
  bool write_collisions_;

  /// Option to write particles tree.
  bool write_particles_;

  /// Option to write particles tree for initial conditions
  bool write_initial_conditions_;

  /// Print only final particles in the event, no intermediate output.
  OutputOnlyFinal particles_only_final_;

  /**
   * Root file cannot be read if it was not properly closed and finalized.
   * It can happen that SMASH simulation crashed and root file was not closed.
   * To save results of simulation in such case, "AutoSave" is
   * applied every N events. The autosave_frequency_ sets
   * this N (default N = 1000). Note that "AutoSave" operation is very
   * time-consuming, so the Autosave_Frequency is
   * always a compromise between safety and speed.
   */
  int autosave_frequency_;

  /// Whether extended particle output is on
  const bool part_extended_;
  /// Whether extended collisions output is on
  const bool coll_extended_;
  /// Whether extended ic output is on
  const bool ic_extended_;

  /**
   * Used for initializing buffer vectors to the maximum allowed size.
   */
  template <typename T>
  void resize_vector(std::vector<T> &vec) {
    for (int i = 0; i < max_buffer_size_; i++) {
      vec.push_back(0.0);
    }
  }

  /**
   * Basic initialization routine, creating the TTree objects
   * for particles and collisions.
   */
  void init_trees();
};

}  // namespace smash

#endif  // SRC_INCLUDE_ROOTOUTPUT_H_
