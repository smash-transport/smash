/*
 *
 *    Copyright (c) 2014
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

/**
 * RootOuput constructor. Creates file smash_run.root in run directory
 * and TTree with the name "particles" in it.
 */
RootOutput::RootOutput(bf::path path, Configuration&& conf)
    : base_path_(std::move(path)),
      root_out_file_(
          new TFile((base_path_ / "smash_run.root").native().c_str(), "NEW")),
      output_counter_(0),
      write_collisions_(conf.has_value({"write_collisions"})
                                 ? conf.take({"write_collisions"}) : false),
      write_particles_(conf.has_value({"write_particles"})
                                 ? conf.take({"write_particles"}) : true),
      autosave_frequency_(conf.has_value({"autosave_frequency"})
                               ? conf.take({"autosave_frequency"}) : 1000) {
  /*!\Userguide
   * \page input_root ROOT
   * Enables output in a ROOT format. The latter is a structured binary format
   * invented at CERN. For more details see root.cern.ch. The resulting
   * output file can optionally contain two TTree's: the one containing
   * information about particles snapshots at fixed moments of time and
   * the one containing information about collisions.
   *
   *
   * \key write_collisions: \n
   * Optional, chooses if information about collisions, decays and box wall
   * crossings should be written out (true) or not (false, default).
   *
   * \key write_particles: \n
   * Optional, chooses if particles snapshots at fixed moments of time
   * should be written out (true, default) or not (false).
   *
   * \key autosave_frequency: \n
   * Root file cannot be read if it was not properly closed and finalized.
   * It can happen that SMASH simulation crashed and root file was not closed.
   * To save results of simulation in such case, "AutoSave" is
   * applied every N events. The autosave_frequency option sets
   * this N (default N = 1000). Note that "AutoSave" operation is very
   * time-consuming, so the autosave_frequency is
   * always a compromise between safety and speed.
   */

  /*!\Userguide
   * \page format_root ROOT format
   * Here there is INFO about ROOT format
   **/

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
  if (write_particles_) {
    output_counter_ = 0;
    particles_to_tree(particles, event_number);
    output_counter_++;
  }
}

/**
 * Writes to tree "at_tstep_N", where N is timestep number counting from 1.
 */
void RootOutput::at_intermediate_time(const Particles &particles,
                                    const int event_number,
                                    const Clock &/*clock*/) {
  if (write_particles_) {
    // This is needed by collision output
    current_event_ = event_number;

    particles_to_tree(particles, event_number);
    output_counter_++;
  }
}

/**
 * Writes to tree "at_eventend".
 */
void RootOutput::at_eventend(const Particles &/*particles*/,
                             const int event_number) {
  // Forced regular dump from operational memory to disk. Very demanding!
  // If program crashes written data will NOT be lost
  if (event_number > 0  && event_number % autosave_frequency_ == 0) {
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
                                   const ParticleList &outgoing) {
  if (write_collisions_) {
    collisions_to_tree(incoming, outgoing);
  }
}

/**
 * Writes particles to a tree defined by treename.
 */
void RootOutput::particles_to_tree(const Particles &particles,
                                   const int event_number) {
  int i = 0;

  tcounter = output_counter_;
  ev = event_number;

  for (const auto &p : particles.data()) {
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
                                    const ParticleList &outgoing) {
  ev = current_event_;
  nin = incoming.size();
  nout = outgoing.size();
  npart = nin + nout;

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
