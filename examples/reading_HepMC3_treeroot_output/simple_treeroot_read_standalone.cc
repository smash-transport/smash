// -*- C++ -*-
// Example code kindly provided by the HepMC3
// developer Andrii Verbytskyi <andrii.verbytskyi@mpp.mpg.de>,
// with minor modifications.

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

#include <iostream>
#include <vector>

const Int_t kMaxparticles = 50000;
const Int_t kMaxvertices = 50000;

class SomeAnalysis {
 public:
  /// pointer to the analyzed TTree or TChain
  TChain *fChain;
  /// Event number in generated MC, HepMC trees only
  Int_t event_number;
  /// HepMC units GEV/MEV, HepMC trees only
  Int_t momentum_unit;
  /// HepMC units cm/mm, HepMC trees only
  Int_t length_unit;
  /// Number of HepMC particles, HepMC trees only
  Int_t particles_;
  /// HepMC PDG id, size[particles_], HepMC trees only
  Int_t particles_pid[kMaxparticles];
  /// HepMC particles status, size[particles_], HepMC trees only
  Int_t particles_status[kMaxparticles];
  /// Some flag, see HepMC docs, size[particles_], HepMC trees only
  Bool_t particles_is_mass_set[kMaxparticles];
  /// HepMC particles mass, size[particles_], HepMC trees only
  Double_t particles_mass[kMaxparticles];
  /// HepMC particles \f$ p_x \f$, size[particles_], HepMC trees only
  Double_t particles_momentum_m_v1[kMaxparticles];
  /// HepMC particles \f$ p_y \f$, size[particles_], HepMC trees only
  Double_t particles_momentum_m_v2[kMaxparticles];
  /// HepMC particles \f$ p_z \f$, size[particles_], HepMC trees only
  Double_t particles_momentum_m_v3[kMaxparticles];
  /// HepMC particles \f$ e \f$, size [particles_], HepMC trees only
  Double_t particles_momentum_m_v4[kMaxparticles];
  /// Number of HepMC vertices, HepMC trees only
  Int_t vertices_;
  /// HepMC vertex \f$x\f$, size[vertices_], HepMC trees only
  Double_t vertices_position_m_v1[kMaxvertices];
  /// HepMC vertex \f$y\f$, size[vertices_], HepMC trees only
  Double_t vertices_position_m_v2[kMaxvertices];
  /// HepMC vertex \f$z\f$, size[vertices_], HepMC trees only
  Double_t vertices_position_m_v3[kMaxvertices];
  /// HepMC vertex \f$t\f$, size[vertices_], HepMC trees only
  Double_t vertices_position_m_v4[kMaxvertices];
  /**
   * Relations between particles and vertices,
   * see HepMC documentation, HepMC trees only
   */
  std::vector<int> links1, links2;
  /// HepMC event attribute id, see HepMC documentation, HepMC trees only
  std::vector<int> attribute_id;
  /// HepMC event attribute name , see HepMC documentation, HepMC trees only
  std::vector<std::string> attribute_name;
  /// HepMC event attribute string, see HepMC documentation, HepMC trees only
  std::vector<std::string> attribute_string;
  /// SMASH internal total and partial weights of the interaction
  std::vector<double> weights;
  /// Name of the program that generated the data
  std::vector<std::string> tool_name;
  /// Version of the program that generated the data
  std::vector<std::string> tool_version;

  TBranch *b_hepmc3_event_event_number;
  TBranch *b_hepmc3_event_momentum_unit;
  TBranch *b_hepmc3_event_length_unit;
  TBranch *b_hepmc3_event_particles_;
  TBranch *b_particles_pid;
  TBranch *b_particles_status;
  TBranch *b_particles_is_mass_set;
  TBranch *b_particles_mass;
  TBranch *b_particles_momentum_m_v1;
  TBranch *b_particles_momentum_m_v2;
  TBranch *b_particles_momentum_m_v3;
  TBranch *b_particles_momentum_m_v4;
  TBranch *b_hepmc3_event_vertices_;
  TBranch *b_vertices_position_m_v1;
  TBranch *b_vertices_position_m_v2;
  TBranch *b_vertices_position_m_v3;
  TBranch *b_vertices_position_m_v4;
  TBranch *b_hepmc3_event_links1;
  TBranch *b_hepmc3_event_links2;
  TBranch *b_hepmc3_event_attribute_id;
  TBranch *b_hepmc3_event_attribute_name;
  TBranch *b_hepmc3_event_attribute_string;
  TBranch *b_hepmc3_event_weights;
  TBranch *b_tool_name;
  TBranch *b_tool_version;

  void Init(TChain *treeH) {
    if (!treeH)
      return;
    fChain = treeH;
    fChain->SetMakeClass(1);
    fChain->SetBranchAddress("event_number", &event_number,
                             &b_hepmc3_event_event_number);
    fChain->SetBranchAddress("momentum_unit", &momentum_unit,
                             &b_hepmc3_event_momentum_unit);
    fChain->SetBranchAddress("length_unit", &length_unit,
                             &b_hepmc3_event_length_unit);
    fChain->SetBranchAddress("particles", &particles_,
                             &b_hepmc3_event_particles_);
    fChain->SetBranchAddress("particles.pid", particles_pid, &b_particles_pid);
    fChain->SetBranchAddress("particles.status", particles_status,
                             &b_particles_status);
    fChain->SetBranchAddress("particles.is_mass_set", particles_is_mass_set,
                             &b_particles_is_mass_set);
    fChain->SetBranchAddress("particles.mass", particles_mass,
                             &b_particles_mass);
    fChain->SetBranchAddress("particles.momentum.m_v1", particles_momentum_m_v1,
                             &b_particles_momentum_m_v1);
    fChain->SetBranchAddress("particles.momentum.m_v2", particles_momentum_m_v2,
                             &b_particles_momentum_m_v2);
    fChain->SetBranchAddress("particles.momentum.m_v3", particles_momentum_m_v3,
                             &b_particles_momentum_m_v3);
    fChain->SetBranchAddress("particles.momentum.m_v4", particles_momentum_m_v4,
                             &b_particles_momentum_m_v4);
    fChain->SetBranchAddress("vertices", &vertices_, &b_hepmc3_event_vertices_);
    fChain->SetBranchAddress("vertices.position.m_v1", vertices_position_m_v1,
                             &b_vertices_position_m_v1);
    fChain->SetBranchAddress("vertices.position.m_v2", vertices_position_m_v2,
                             &b_vertices_position_m_v2);
    fChain->SetBranchAddress("vertices.position.m_v3", vertices_position_m_v3,
                             &b_vertices_position_m_v3);
    fChain->SetBranchAddress("vertices.position.m_v4", vertices_position_m_v4,
                             &b_vertices_position_m_v4);
    fChain->SetBranchAddress("links1", &links1, &b_hepmc3_event_links1);
    fChain->SetBranchAddress("links2", &links2, &b_hepmc3_event_links2);
    fChain->SetBranchAddress("attribute_id", &attribute_id,
                             &b_hepmc3_event_attribute_id);
    fChain->SetBranchAddress("attribute_name", &attribute_name,
                             &b_hepmc3_event_attribute_name);
    fChain->SetBranchAddress("attribute_string", &attribute_string,
                             &b_hepmc3_event_attribute_string);
    fChain->SetBranchAddress("weights", &weights, &b_hepmc3_event_weights);
    fChain->SetBranchAddress("tool_name", &tool_name, &b_tool_name);
    fChain->SetBranchAddress("tool_version", &tool_version, &b_tool_version);
  }
  SomeAnalysis(char *file) {
    TChain *TempChain = new TChain("hepmc3_tree");
    TempChain->Add(file);
    Init(TempChain);
  }
};

int main(int argc, char **argv) {
  if (argc != 2) {
    printf("Usage: %s <path to SMASH_HepMC_*.treeroot>\n", argv[0]);
    exit(1);
  }
  SomeAnalysis *A = new SomeAnalysis(argv[1]);
  for (int entry = 0; entry < A->fChain->GetEntries(); entry++) {
    A->fChain->GetEntry(entry);
    std::cout << "\ntool name: " << A->tool_name[0] << std::endl;
    std::cout << "tool version: " << A->tool_version[0] << std::endl;
    printf("Number of particles: %i\n", A->particles_);
    std::cout << "Fields: particle id, particle PDG id, px, py, pz, E"
              << std::endl;
    for (int np = 0; np < A->particles_; np++) {
      std::cout << np << "   " << A->particles_pid[np] << "   ";
      std::cout << A->particles_momentum_m_v1[np] << "   ";
      std::cout << A->particles_momentum_m_v2[np] << "   ";
      std::cout << A->particles_momentum_m_v3[np] << "   ";
      std::cout << A->particles_momentum_m_v4[np] << std::endl;
    }
    std::cout << std::endl;
  }
  delete A;
  return 0;
}
