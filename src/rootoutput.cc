/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/rootoutput.h"
#include "include/particles.h"
//#include "include/filedeleter.h"
//#include <memory>
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
  root_out_file_->Write();
  root_out_file_->Close();
}

/**
 * Writes to tree "at_eventstart".
 */
void RootOutput::at_eventstart(const Particles &particles, const int event_number) {
  particles_to_tree("at_eventstart","Initial particles", particles, event_number);
}

void RootOutput::before_collision() {
}

void RootOutput::after_collision() {
}

/**
 * Writes to tree "at_tstep_N", where N is timestep number counting from 1.
 */
void RootOutput::after_Nth_timestep(const Particles &particles, const int event_number, const int timestep) {
  char treename[32], treedescr[64];
  snprintf(treename, sizeof(treename), "at_tstep_%07i", timestep + 1);
  snprintf(treedescr, sizeof(treedescr), "Particles after timestep %i", timestep + 1);
  particles_to_tree(treename, treedescr, particles, event_number);
}

/**
 * Writes to tree "at_eventend".
 */
void RootOutput::at_eventend(const Particles &particles, const int event_number) {
  particles_to_tree("at_eventend","Final particles", particles, event_number);
}

/**
 * Writes particles to a tree defined by treename.
 */
void RootOutput::particles_to_tree(const char* treename, const char* treedescr, const Particles &particles, const int event_number) {

  TTree *curr_tree = (TTree*)root_out_file_->Get(treename);
  if (curr_tree == NULL) {
    tree_list_.emplace_back(curr_tree = new TTree(treename, treedescr));
  }

  double p0,px,py,pz;
  double x,y,z;
  int    id, pdgcode, ev;

  curr_tree->Branch("p0",&p0,"p0/D");
  curr_tree->Branch("px",&px,"px/D");
  curr_tree->Branch("py",&py,"py/D");
  curr_tree->Branch("pz",&pz,"pz/D");

  curr_tree->Branch("x",&x,"x/D");
  curr_tree->Branch("y",&y,"y/D");
  curr_tree->Branch("z",&z,"z/D");

  curr_tree->Branch("pdgcode",&pdgcode,"pdgcode/I");
  curr_tree->Branch("id",&id,"id/I");
  curr_tree->Branch("ev",&ev,"ev/I");

  for (const auto &p : particles) {
    x = p.second.position().x1();
    y = p.second.position().x2();
    z = p.second.position().x3();

    p0 = p.second.momentum().x0();
    px = p.second.momentum().x1();
    py = p.second.momentum().x2();
    pz = p.second.momentum().x3();

    pdgcode = p.second.pdgcode();
    id = p.second.id();

    ev = event_number;

    curr_tree->Fill();
  }

}

}  // namespace Smash
