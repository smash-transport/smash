/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/clock.h"
#include "include/forwarddeclarations.h"
#include "include/inputfunctions.h"
#include "include/particles.h"
#include "include/rootoutput.h"
#include "TFile.h"
#include "TTree.h"

namespace Smash {

RootOutput::RootOutput(const bf::path &path, const std::string &name)
    : base_path_(std::move(path)),
      root_out_file_(
          new TFile((base_path_ / (name + ".root")).native().c_str(), "NEW")),
      output_counter_(0),
      write_collisions_(true), write_particles_(false),
      autosave_frequency_(1000) {
  init_trees();
}

/**
 * RootOuput constructor. Creates file smash_run.root in the output directory.
 */
RootOutput::RootOutput(const bf::path &path, Configuration&& conf)
    : base_path_(std::move(path)),
      root_out_file_(
          new TFile((base_path_ / "smash_run.root").native().c_str(), "NEW")),
      output_counter_(0),
      write_collisions_(conf.take({"Write_Collisions"}, false)),
      write_particles_(conf.take({"Write_Particles"}, true)),
      autosave_frequency_(conf.take({"Autosave_Frequency"}, 1000)) {
  /*!\Userguide
   * \page input_root ROOT
   * Enables output in a ROOT format. The latter is a structured binary format
   * invented at CERN. For more details see root.cern.ch. The resulting
   * output file can optionally contain two TTree's: the one containing
   * information about the particle list at fixed moments of time and
   * the other one containing information about the collision history.
   *
   * \key Enable (bool, optional, default = false):\n
   * true - Root output enabled\n
   * false - no Root output
   *
   * \key Write_Collisions (bool, optional, default = false): \n
   * true - information about collisions, decays and box wall
   * crossings will be written out \n
   * false - collision history information is suppressed
   *
   * \key Write_Particles (bool, optional, default = true): \n
   * true - particle list at output interval is written out \n
   * false - particle list is not written out
   *
   * \key Autosave_Frequency (int, optional, default = 1000): \n
   * Root file cannot be read if it was not properly closed and finalized.
   * It can happen that SMASH simulation crashed and root file was not closed.
   * To save results of simulation in such case, "AutoSave" is
   * applied every N events. The autosave_frequency option sets
   * this N (default N = 1000). Note that "AutoSave" operation is very
   * time-consuming, so the Autosave_Frequency is
   * always a compromise between safety and speed.
   *
   * For details on ROOT output format see \ref format_root
   */

  /*!\Userguide
   * \page format_root ROOT format
   * SMASH ROOT output has the same functionality as OSCAR output, but
   * root-files are faster to read, write and they need less disk space
   * for the same amount of information. This is reached due to complex
   * internal structure of ROOT file. ROOT file is not human-readable, but it
   * can be viewed using ROOT TBrowser. One can also open and read it using
   * ROOT functions. Full memory structure of the ROOT-file can be found
   * here http://root.cern.ch/root/html/TFile.html. We only desribe the logical
   * structure of the SMASH root output. Knowing the logical structure is
   * enough to read, write root-file and understand its view in TBrowser.
   *
   * Producing ROOT output requires ROOT installed (see http://root.cern.ch).
   *
   * SMASH produces one root file per run: \c smash_run.root. This file
   * contains TTree called \c particles and a TTree \c collisions.
   * The latter can be switched on and off by an option (see \ref input_root).
   * Particles tree contains the same information as OSCAR particles output
   * and collisions tree contains the same information as OSCAR collision output.
   *
   * Every physical quantity is in separate TBranch.\n
   * One entry in the particles TTree is:
   * \code
   * ev tcounter npart pdgcode[npart] t[npart] x[npart] y[npart] z[npart] ->
   * -> p0[npart] px[npart] py[npart] pz[npart]
   * \endcode
   * One tree entry is analogous to OSCAR output block, but maximal size of
   * ROOT entry is limited to 10000. This is done to limit buffer size needed
   * for root output. If number of particles in one block exceeds 10000,
   * then they are written in separate blocks with the same tcounter and ev.
   *
   * \li ev is event number
   * \li tcounter is number of output block in a given event in terms of OSCAR
   * \li npart is number of particles in the block
   * \li pdgcode is PDG id array
   * \li t, x, y, z are position arrays
   * \li p0, px, py, pz are 4-momenta arrays
   *
   * In collisions tree entries are organized in the same way, but additional
   * fields \c nin and \c nout are added to characterize number of incoming
   * and outgoing particles in the reaction, nin + nout = npart. Currently
   * writing initial and final configuration to collisions tree is not supported.
   *
   * See also \ref collisions_output_in_box_modus_.
   **/

  init_trees();
}

void RootOutput::init_trees() {
  if (write_particles_) {
    particles_tree_ = new TTree("particles", "particles");

    particles_tree_->Branch("npart", &npart, "npart/I");
    particles_tree_->Branch("ev", &ev, "ev/I");
    particles_tree_->Branch("tcounter", &tcounter, "tcounter/I");

    particles_tree_->Branch("pdgcode", &pdgcode[0], "pdgcode[npart]/I");

    particles_tree_->Branch("p0", &p0[0], "p0[npart]/D");
    particles_tree_->Branch("px", &px[0], "px[npart]/D");
    particles_tree_->Branch("py", &py[0], "py[npart]/D");
    particles_tree_->Branch("pz", &pz[0], "pz[npart]/D");

    particles_tree_->Branch("t", &t[0], "t[npart]/D");
    particles_tree_->Branch("x", &x[0], "x[npart]/D");
    particles_tree_->Branch("y", &y[0], "y[npart]/D");
    particles_tree_->Branch("z", &z[0], "z[npart]/D");
  }

  if (write_collisions_) {
    collisions_tree_ = new TTree("collisions", "collisions");

    collisions_tree_->Branch("nin", &nin, "nin/I");
    collisions_tree_->Branch("nout", &nout, "nout/I");
    collisions_tree_->Branch("npart", &npart, "npart/I");
    collisions_tree_->Branch("ev", &ev, "ev/I");
    collisions_tree_->Branch("weight", &wgt, "weight/D");

    collisions_tree_->Branch("pdgcode", &pdgcode[0], "pdgcode[npart]/I");

    collisions_tree_->Branch("p0", &p0[0], "p0[npart]/D");
    collisions_tree_->Branch("px", &px[0], "px[npart]/D");
    collisions_tree_->Branch("py", &py[0], "py[npart]/D");
    collisions_tree_->Branch("pz", &pz[0], "pz[npart]/D");

    collisions_tree_->Branch("t", &t[0], "t[npart]/D");
    collisions_tree_->Branch("x", &x[0], "x[npart]/D");
    collisions_tree_->Branch("y", &y[0], "y[npart]/D");
    collisions_tree_->Branch("z", &z[0], "z[npart]/D");
  }
}

/**
 * RootOutput destructor. Writes root objects (here TTrees) to file and closes it.
 */
RootOutput::~RootOutput() {
  // kOverwrite option prevents from writing extra TKey objects into root file
  root_out_file_->Write("", TObject::kOverwrite);
  root_out_file_->Close();
}

/**
 * Writes to tree "at_eventstart".
 */
void RootOutput::at_eventstart(const Particles &particles,
                               const int event_number) {
  // save event number
  current_event_ = event_number;

  if (write_particles_) {
    output_counter_ = 0;
    particles_to_tree(particles);
    output_counter_++;
  }
}

/**
 * Writes to tree "at_tstep_N", where N is timestep number counting from 1.
 */
void RootOutput::at_intermediate_time(const Particles &particles,
                                    const int /*event_number*/,
                                    const Clock &/*clock*/) {
  if (write_particles_) {
    particles_to_tree(particles);
    output_counter_++;
  }
}

/**
 * Writes to tree "at_eventend".
 */
void RootOutput::at_eventend(const Particles &/*particles*/,
                             const int /*event_number*/) {
  // Forced regular dump from operational memory to disk. Very demanding!
  // If program crashes written data will NOT be lost
  if (current_event_ > 0  && current_event_ % autosave_frequency_ == 0) {
    if (write_particles_) {
      particles_tree_->AutoSave("SaveSelf");
    }
    if (write_collisions_) {
      collisions_tree_->AutoSave("SaveSelf");
    }
  }
}

/**
 * Writes interactions to ROOT-file
 */
void RootOutput::at_interaction(const ParticleList &incoming,
                                const ParticleList &outgoing,
                                const double /*density*/,
                                const double weight,
                                ProcessType /*process_type*/) {
  if (write_collisions_) {
    collisions_to_tree(incoming, outgoing, weight);
  }
}

/**
 * Writes particles to a tree defined by treename.
 */
void RootOutput::particles_to_tree(const Particles &particles) {
  int i = 0;

  tcounter = output_counter_;
  ev = current_event_;

  for (const auto &p : particles) {
    // Buffer full - flush to tree, else fill with particles
    if (i >= max_buffer_size_) {
      npart = max_buffer_size_;
      i = 0;
      particles_tree_->Fill();
    } else {
      t[i] = p.position().x0();
      x[i] = p.position().x1();
      y[i] = p.position().x2();
      z[i] = p.position().x3();

      p0[i] = p.momentum().x0();
      px[i] = p.momentum().x1();
      py[i] = p.momentum().x2();
      pz[i] = p.momentum().x3();

      pdgcode[i] = p.pdgcode().get_decimal();

      i++;
    }
  }
  // Flush rest to tree
  if (i > 0) {
    npart = i;
    particles_tree_->Fill();
  }
}

void RootOutput::collisions_to_tree(const ParticleList &incoming,
                                    const ParticleList &outgoing,
                                    const double weight) {
  ev = current_event_;
  nin = incoming.size();
  nout = outgoing.size();
  npart = nin + nout;
  wgt = weight;

  int i = 0;

  // It is assumed that nin + nout < max_buffer_size_
  // This is true for any possible reaction for current buffer size: 10000
  // But if one wants initial/final particles written to collisions
  // then implementation should be updated.

  for (const auto &p : incoming) {
    t[i] = p.position().x0();
    x[i] = p.position().x1();
    y[i] = p.position().x2();
    z[i] = p.position().x3();

    p0[i] = p.momentum().x0();
    px[i] = p.momentum().x1();
    py[i] = p.momentum().x2();
    pz[i] = p.momentum().x3();

    pdgcode[i] = p.pdgcode().get_decimal();

    i++;
  }

  for (const auto &p : outgoing) {
    t[i] = p.position().x0();
    x[i] = p.position().x1();
    y[i] = p.position().x2();
    z[i] = p.position().x3();

    p0[i] = p.momentum().x0();
    px[i] = p.momentum().x1();
    py[i] = p.momentum().x2();
    pz[i] = p.momentum().x3();

    pdgcode[i] = p.pdgcode().get_decimal();

    i++;
  }

  collisions_tree_->Fill();
}
}  // namespace Smash
