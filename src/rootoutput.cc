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
#include "include/particles.h"
#include "include/rootoutput.h"
#include "TFile.h"
#include "TTree.h"

namespace Smash {

/**
 * RootOuput constructor. Creates file smash_run.root in run directory
 * and TTree with the name "particles" in it.
 */
RootOutput::RootOutput(bf::path path, Options op)
    : base_path_(std::move(path)),
      root_out_file_(
          new TFile((base_path_ / "smash_run.root").native().c_str(), "NEW")),
      output_counter_(0) {
  tree_ = new TTree("particles", "particles");

  tree_->Branch("npart", &npart, "npart/I");
  tree_->Branch("ev", &ev, "ev/I");
  tree_->Branch("tcounter", &tcounter, "tcounter/I");

  tree_->Branch("pdgcode", &pdgcode[0], "pdgcode[npart]/I");

  tree_->Branch("p0", &p0[0], "p0[npart]/D");
  tree_->Branch("px", &px[0], "px[npart]/D");
  tree_->Branch("py", &py[0], "py[npart]/D");
  tree_->Branch("pz", &pz[0], "pz[npart]/D");

  tree_->Branch("t", &t[0], "t[npart]/D");
  tree_->Branch("x", &x[0], "x[npart]/D");
  tree_->Branch("y", &y[0], "y[npart]/D");
  tree_->Branch("z", &z[0], "z[npart]/D");
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
  output_counter_ = 0;
  particles_to_tree(particles, event_number);
  output_counter_++;
}

/**
 * Writes to tree "at_tstep_N", where N is timestep number counting from 1.
 */
void RootOutput::after_Nth_timestep(const Particles &particles,
                                    const int event_number,
                                    const Clock &/*clock*/) {
  particles_to_tree(particles, event_number);
  output_counter_++;
}

/**
 * Writes to tree "at_eventend".
 */
void RootOutput::at_eventend(const Particles &/*particles*/,
                             const int event_number) {
  // Forced dump from operational memory to disk every 10 events
  // If program crashes written data will NOT be lost
  if (event_number%10 == 1) {
    tree_->AutoSave("SaveSelf");
  }
}

/**
 * Writes interactions to ROOT-file
 */
void RootOutput::write_interaction(const ParticleList &/*incoming_particles*/,
                                   const ParticleList &/*outgoing_particles*/) {
//  for (const auto &p : incoming_particles) {
//  }
//  for (const auto &p : outgoing_particles) {
//  }
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
      tree_->Fill();
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
    tree_->Fill();
  }
}
}  // namespace Smash
