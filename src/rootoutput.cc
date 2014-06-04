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
 * RootOuput constructor. Creates file smash_run.root in run directory.
 */
RootOutput::RootOutput(boost::filesystem::path path)
    : base_path_(std::move(path)),
      root_out_file_(
          new TFile((base_path_ / "smash_run.root").native().c_str(), "NEW")) {}

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
  particles_to_tree("at_eventstart", "Initial particles",
                    particles, event_number);
}

/**
 * Writes to tree "at_tstep_N", where N is timestep number counting from 1.
 */
void RootOutput::after_Nth_timestep(const Particles &particles,
                                    const int event_number,
                                    const Clock &clock) {
  char treename[32], treedescr[64];
  snprintf(treename, sizeof(treename),
           "at_time_%8.3f", clock.current_time());
  snprintf(treedescr, sizeof(treedescr),
           "Particles after time %8.3f", clock.current_time());
  particles_to_tree(treename, treedescr, particles, event_number);
}

/**
 * Writes to tree "at_eventend".
 */
void RootOutput::at_eventend(const Particles &particles,
                             const int event_number) {
  particles_to_tree("at_eventend", "Final particles",
                    particles, event_number);
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
void RootOutput::particles_to_tree(const char* treename,
                                   const char* treedescr,
                                   const Particles &particles,
                                   const int event_number) {
  TTree *curr_tree = static_cast<TTree*>(root_out_file_->Get(treename));
  if (curr_tree == NULL) {
    tree_list_.emplace_back(curr_tree = new TTree(treename, treedescr));
  }

  // Forced dump from operational memory to disk every 10 events
  // If program crashes written data will NOT be lost
  if (event_number%10 == 1) {
    curr_tree->AutoSave("SaveSelf");
  }


  double p0, px, py, pz;
  double x, y, z;
  int    id, pdgcode, ev;

  curr_tree->Branch("p0", &p0, "p0/D");
  curr_tree->Branch("px", &px, "px/D");
  curr_tree->Branch("py", &py, "py/D");
  curr_tree->Branch("pz", &pz, "pz/D");

  curr_tree->Branch("x", &x, "x/D");
  curr_tree->Branch("y", &y, "y/D");
  curr_tree->Branch("z", &z, "z/D");

  curr_tree->Branch("pdgcode", &pdgcode, "pdgcode/I");
  curr_tree->Branch("id", &id, "id/I");
  curr_tree->Branch("ev", &ev, "ev/I");

  for (const auto &p : particles.data()) {
    x = p.position().x1();
    y = p.position().x2();
    z = p.position().x3();

    p0 = p.momentum().x0();
    px = p.momentum().x1();
    py = p.momentum().x2();
    pz = p.momentum().x3();

    pdgcode = p.pdgcode().get_decimal();
    id = p.id();

    ev = event_number;

    curr_tree->Fill();
  }
}
}  // namespace Smash
