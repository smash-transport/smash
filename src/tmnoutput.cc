/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/Tmnoutput.h"
#include "include/particles.h"
//#include "include/filedeleter.h"
//#include <memory>
#include "TFile.h"
#include "TTree.h"

namespace Smash {

TmnOutput::TmnOutput(boost::filesystem::path path)
    : base_path_(std::move(path)) {
}

TmnOutput::~TmnOutput() {
}


void TmnOutput::at_runstart(){
 root_out_file = new TFile( (base_path_ / "smash_run_Tmn.root").native().c_str(),"NEW");
}


void TmnOutput::at_eventstart(const Particles &particles, const int evt_num){
 Tmn_to_tree("at_eventstart","Initial Tmn", particles, evt_num);
}


//void TmnOutput::at_collision(const Collisions &collisions){}


void TmnOutput::at_outtime(const Particles &particles, const int evt_num, const int timestep){
 char treename[32], treedescr[64];
 snprintf(treename, sizeof(treename), "at_tstep_%07i",timestep);
 snprintf(treedescr, sizeof(treedescr), "Tmn after timestep %i",timestep);
 Tmn_to_tree(treename,treedescr, particles, evt_num);
}


void TmnOutput::at_eventend(const Particles &particles, const int evt_num){
 Tmn_to_tree("at_eventend","Final Tmn", particles, evt_num);
}

void TmnOutput::at_runend(){
 root_out_file->Write();
 root_out_file->Close();
}

void TmnOutput::at_crash(){
 root_out_file->Write();
 root_out_file->Close();
}


void TmnOutput::Tmn_to_tree(const char* treename, const char* treedescr, const Particles &particles, const int evt_num){

 curr_tree = (TTree*)root_out_file->Get(treename);
 if (curr_tree == NULL){
   curr_tree = new TTree(treename, treedescr);
 }

 curr_tree->Branch("t00",&t00,"t00/D");
 curr_tree->Branch("t01",&t01,"t01/D");
 curr_tree->Branch("t02",&t02,"t02/D");
 curr_tree->Branch("t03",&t03,"t03/D");
 curr_tree->Branch("t11",&t11,"t11/D");
 curr_tree->Branch("t12",&t12,"t12/D");
 curr_tree->Branch("t13",&t13,"t13/D");
 curr_tree->Branch("t22",&t22,"t22/D");
 curr_tree->Branch("t23",&t23,"t23/D");
 curr_tree->Branch("t33",&t33,"t33/D");

 curr_tree->Branch("ev",&ev,"ev/I");

 t00 = 0.0;
 t01 = 0.0;
 t02 = 0.0;
 t03 = 0.0;
 t11 = 0.0;
 t12 = 0.0;
 t13 = 0.0;
 t22 = 0.0;
 t23 = 0.0;
 t33 = 0.0;


 for (const auto &p : particles) {

    p0 = p.second.momentum().x0();
    px = p.second.momentum().x1();
    py = p.second.momentum().x2();
    pz = p.second.momentum().x3();
 
    t00+= p0;
    t01+= px;
    t02+= py;
    t03+= pz;
    t11+= px*px/p0;
    t12+= px*py/p0;
    t13+= px*pz/p0;
    t22+= py*py/p0;
    t23+= py*pz/p0;
    t33+= pz*pz/p0;

 }

 ev = evt_num;

 curr_tree->Fill();

}




}  // namespace Smash
