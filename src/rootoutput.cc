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

RootOutput::RootOutput(boost::filesystem::path path)
    : base_path_(std::move(path)),
      root_out_file(
          new TFile((base_path_ / "smash_run.root").native().c_str(), "NEW")) {}

RootOutput::~RootOutput() {
 root_out_file->Write();
 root_out_file->Close();
}


void RootOutput::at_runstart(){
}


void RootOutput::at_eventstart(const Particles &particles, const int evt_num){
 particles_to_tree("at_eventstart","Initial particles", particles, evt_num);
}


//void RootOutput::at_collision(const Collisions &collisions){}


void RootOutput::at_outtime(const Particles &particles, const int evt_num, const int timestep){
 char treename[32], treedescr[64];
 snprintf(treename, sizeof(treename), "at_tstep_%07i",timestep);
 snprintf(treedescr, sizeof(treedescr), "Particles after timestep %i",timestep);
 particles_to_tree(treename,treedescr, particles, evt_num);
}


void RootOutput::at_eventend(const Particles &particles, const int evt_num){
 particles_to_tree("at_eventend","Final particles", particles, evt_num);
}

void RootOutput::at_runend(){
}

void RootOutput::at_crash(){
 root_out_file->Write();
 root_out_file->Close();
}


void RootOutput::particles_to_tree(const char* treename, const char* treedescr, const Particles &particles, const int evt_num){

  TTree *curr_tree = (TTree*)root_out_file->Get(treename);
  if (curr_tree == NULL) {
    tree_list_.emplace_back(curr_tree = new TTree(treename, treedescr));
  }

   double p0,px,py,pz;
   double t,x,y,z;
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

    ev = evt_num;

    curr_tree->Fill();
 }

}




}  // namespace Smash
