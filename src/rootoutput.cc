/*
 *
 *    Copyright (c) 2014-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/rootoutput.h"

#include "TFile.h"
#include "TTree.h"
#include "smash/action.h"
#include "smash/clock.h"
#include "smash/forwarddeclarations.h"
#include "smash/particles.h"

namespace smash {

const int RootOutput::max_buffer_size_ = 500000;

/*!\Userguide
 * \page doxypage_output_root
 * SMASH ROOT output is a fast and disk-space efficient, but not human-readable
 * output. It is a custom format making information about the SMASH calculation
 * accessible with ROOT, mostly mirroring the information of the
 * \ref doxypage_output_oscar "OSCAR output" formats. This output is distinct
 * from the standarized \ref doxypage_output_hepmc "HepMC output" that is also
 * available in a ROOT format and more widely adopted.
 *
 * SMASH ROOT output files can be viewed using ROOT's TBrowser. One can also
 * access them using ROOT functions. The full memory structure of the ROOT
 * files can be found here: http://root.cern.ch/root/html/TFile.html. We only
 * describe the logical structure of the SMASH ROOT output. Knowing the logical
 * structure is enough to read and write ROOT files, be able to view them in
 * TBrowser, or write a ROOT macro to analyze them.
 *
 * Producing ROOT output requires ROOT installed (see http://root.cern.ch).
 *
 * Depending on configuration (see \ref doxypage_output "Output"), SMASH can
 * produce two ROOT files per run: \c Particles.root and \c Collisions.root.
 * These files contain a TTree called \c particles and a TTree called
 * \c collisions, respectively. The \c particles tree contains information about
 * the parameters of the run (such as the number of testparticles and event
 * number), information relating to individual particles (such as their position
 * or charge), and information about bulk observables in the system (kinetic
 * energy, mean field energy, and total energy). The \c collisions tree contains
 * information about each collision, such as number of incoming and outgoing
 * particles. It also has the full information about the incoming and outgoing\
 * particles of each collision.
 *
 * In case that the ROOT format is used for dilepton output
 * (see \ref doxypage_output_dileptons "Dileption output"), the ROOT file is
 * called \c Dileptons.root and only contains a \c collisions tree with all the
 * dilepton decays.
 *
 * Output is divided into ROOT entries (output blocks). In a ROOT entry,
 * information is stored in TBranches, where every physical quantity corresponds
 * to a separate TBranch. One entry (output block) in the \c particles TTree is
 * comprised of:
 * \code
 * ev ens tcounter npart test_p modus_l current_t impact_b empty_event
 * id[npart] pdgcode[npart] charge[npart] formation_time[npart]
 * time_last_collision[npart] t[npart] x[npart] y[npart] z[npart] p0[npart]
 * px[npart] py[npart] pz[npart] E_kinetic_tot E_fields_tot E_tot
 * \endcode
 * All particles in a ROOT entry (output block) are at the same time and in the
 * same event. The maximal number of particles in one entry is limited to
 * 500000. This is done to limit the buffer size needed for ROOT output. If the
 * number of particles in one ROOT entry (output block) exceeds 500000, then
 * they are written in separate entries with the same \c tcounter and \c ev.
 *
 * Each ROOT entry (output block) contains TBranches with overall information
 * about the event and the output block:
 * \li \c ev is event number
 * \li \c ens is ensemble number
 * \li \c tcounter is number of the output block in a given event
 * \li \c npart is number of particles in the block
 * \li \c test_p is number of testparticles per particle
 * \li \c modus_l is modus length
 * \li \c current_t is time associated with the output block, in fm
 * \li \c impact_b is the impact parameter of the event
 * \li \c empty_event indicates whether the projectile did not interact with the
 * target
 *
 * Then, each ROOT entry (output block) contains information about all particles
 * in the block. Particle characteristics such as position or charge are stored
 * in arrays which belong to separate TBranches. There are TBranches containing
 * arrays with information on the following quantities:
 * \li \c id is unique integer identifier of the particle
 * \li \c pdgcode is PDG id
 * \li \c charge is the electric charge
 * \li \c formation_time is particle formation time
 * \li \c time_last_collision is time of the last collision
 * \li \c p0, \c px, \c py, \c pz are components of particle 4-momentum
 * \li \c t, \c x, \c y, \c z are components of particle 4-position
 *
 * Finally, each ROOT entry (output block) contains information about some of
 * the bulk properties of the system:
 * \li \c E_kinetic_tot is total kinetic energy in the system
 * \li \c E_fields_tot is total mean field energy * test_p
 * \li \c E_tot is the sum of E_kinetic_tot and E_fields_tot
 *
 * In case of extended output
 * (see \ref input_output_content_specific_ "extended output") more fields are
 * added. Their description is the same that in case of OSCAR format, see
 * \ref extended_output_format_ "extended OSCAR format".
 * \note In contrast to the OSCAR format, the ROOT output includes
 * \c formation_time and \c time_last_collision in the standard (not extended)
 * output.
 *
 * If "Collisions:" section is present in the config file, then in addition to a
 * file Particles.root with particles TTree, another file Collisions.root is
 * created. Collisions.root contains information about each collision, written
 * as one leaf. Similarly to Particles.root, this includes global information
 * about the event:
 * \code
 * ev ens
 * \endcode
 * and information about particles participating in a collision:
 * \code
 * id[npart] pdgcode[npart] charge[npart] formation_time[npart]
 * time_last_collision[npart] p0[npart] px[npart] py[npart] pz[npart] t[npart]
 * x[npart] y[npart] z[npart]
 * \endcode
 * where all arrays are now of dimension nin+nout = npart.
 * Beyond these, Collisions.root contains a few additional fields:
 * \li \c nin is the number of incoming particles in the reaction,
 * \li \c nout is the number of outgoing particles in the reaction,
 * \li \c wgt is an action weight, whose meaning depends on the type of action:
 * for collisions it is the total cross section, for decays it is the total
 * decay width, and for dilepton decays it is the shining weight,
 * \li \c par_wgt is partial weight of the collision.
 *
 * Currently writing initial and final configurations to collisions tree is not
 * supported.
 *
 * See also \ref doxypage_output_collisions_box_modus.
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
 */
RootOutput::RootOutput(const std::filesystem::path &path,
                       const std::string &name, const OutputParameters &out_par)
    : OutputInterface(name),
      filename_(path / (name + ".root")),
      write_collisions_(name == "Collisions" || name == "Dileptons" ||
                        name == "Photons"),
      write_particles_(name == "Particles"),
      write_initial_conditions_(name == "SMASH_IC"),
      particles_only_final_(out_par.part_only_final),
      part_extended_(out_par.part_extended),
      coll_extended_(out_par.coll_extended),
      ic_extended_(out_par.ic_extended) {
  filename_unfinished_ = filename_;
  filename_unfinished_ += ".unfinished";
  root_out_file_ =
      std::make_unique<TFile>(filename_unfinished_.native().c_str(), "NEW");
  init_trees();
}

void RootOutput::init_trees() {
  if (write_particles_ || write_initial_conditions_) {
    particles_tree_ = new TTree("particles", "particles");

    particles_tree_->Branch("ev", &ev_, "ev/I");
    particles_tree_->Branch("ens", &ens_, "ens/I");
    particles_tree_->Branch("tcounter", &tcounter_, "tcounter/I");
    particles_tree_->Branch("npart", &npart_, "npart/I");
    particles_tree_->Branch("test_p", &test_p_, "test_p/I");
    particles_tree_->Branch("modus_l", &modus_l_, "modus_l/D");
    particles_tree_->Branch("current_t", &current_t_, "current_t/D");
    particles_tree_->Branch("impact_b", &impact_b_, "impact_b/D");
    particles_tree_->Branch("empty_event", &empty_event_, "empty_event/O");

    particles_tree_->Branch("id", &id_[0], "id[npart]/I");
    particles_tree_->Branch("pdgcode", &pdgcode_[0], "pdgcode[npart]/I");
    particles_tree_->Branch("charge", &charge_[0], "charge[npart]/I");
    particles_tree_->Branch("formation_time", &formation_time_[0],
                            "formation_time[npart]/D");
    particles_tree_->Branch("time_last_collision", &time_last_collision_[0],
                            "time_last_collision[npart]/D");

    particles_tree_->Branch("p0", &p0_[0], "p0[npart]/D");
    particles_tree_->Branch("px", &px_[0], "px[npart]/D");
    particles_tree_->Branch("py", &py_[0], "py[npart]/D");
    particles_tree_->Branch("pz", &pz_[0], "pz[npart]/D");

    particles_tree_->Branch("t", &t_[0], "t[npart]/D");
    particles_tree_->Branch("x", &x_[0], "x[npart]/D");
    particles_tree_->Branch("y", &y_[0], "y[npart]/D");
    particles_tree_->Branch("z", &z_[0], "z[npart]/D");

    particles_tree_->Branch("E_kinetic_tot", &E_kinetic_tot_,
                            "E_kinetic_tot/D");
    particles_tree_->Branch("E_fields_tot", &E_fields_tot_, "E_fields_tot/D");
    particles_tree_->Branch("E_tot", &E_tot_, "E_tot/D");

    if (part_extended_ || ic_extended_) {
      particles_tree_->Branch("ncoll", &coll_per_part_[0], "ncoll[npart]/I");
      particles_tree_->Branch("xsecfac", &xsec_factor_[0], "xsecfac[npart]/D");
      particles_tree_->Branch("proc_id_origin", &proc_id_origin_[0],
                              "proc_id_origin[npart]/I");
      particles_tree_->Branch("proc_type_origin", &proc_type_origin_[0],
                              "proc_type_origin[npart]/I");
      particles_tree_->Branch("pdg_mother1", &pdg_mother1_[0],
                              "pdg_mother1[npart]/I");
      particles_tree_->Branch("pdg_mother2", &pdg_mother2_[0],
                              "pdg_mother2[npart]/I");
      particles_tree_->Branch("baryon_number", &baryon_number_[0],
                              "baryon_number[npart]/I");
      particles_tree_->Branch("strangeness", &strangeness_[0],
                              "strangeness[npart]/I");
    }
  }

  if (write_collisions_) {
    collisions_tree_ = new TTree("collisions", "collisions");

    collisions_tree_->Branch("ev", &ev_, "ev/I");
    collisions_tree_->Branch("ens", &ens_, "ens/I");

    collisions_tree_->Branch("nin", &nin_, "nin/I");
    collisions_tree_->Branch("nout", &nout_, "nout/I");
    collisions_tree_->Branch("npart", &npart_, "npart/I");
    collisions_tree_->Branch("weight", &wgt_, "weight/D");
    collisions_tree_->Branch("partial_weight", &par_wgt_, "partial_weight/D");

    collisions_tree_->Branch("id", &id_[0], "id[npart]/I");
    collisions_tree_->Branch("pdgcode", &pdgcode_[0], "pdgcode[npart]/I");
    collisions_tree_->Branch("charge", &charge_[0], "charge[npart]/I");
    collisions_tree_->Branch("formation_time", &formation_time_[0],
                             "formation_time[npart]/D");
    collisions_tree_->Branch("time_last_collision", &time_last_collision_[0],
                             "time_last_collision[npart]/D");

    collisions_tree_->Branch("p0", &p0_[0], "p0[npart]/D");
    collisions_tree_->Branch("px", &px_[0], "px[npart]/D");
    collisions_tree_->Branch("py", &py_[0], "py[npart]/D");
    collisions_tree_->Branch("pz", &pz_[0], "pz[npart]/D");

    collisions_tree_->Branch("t", &t_[0], "t[npart]/D");
    collisions_tree_->Branch("x", &x_[0], "x[npart]/D");
    collisions_tree_->Branch("y", &y_[0], "y[npart]/D");
    collisions_tree_->Branch("z", &z_[0], "z[npart]/D");

    if (coll_extended_) {
      collisions_tree_->Branch("ncoll", &coll_per_part_[0], "ncoll[npart]/I");
      collisions_tree_->Branch("xsecfac", &xsec_factor_[0], "xsecfac[npart]/D");
      collisions_tree_->Branch("proc_id_origin", &proc_id_origin_[0],
                               "proc_id_origin[npart]/I");
      collisions_tree_->Branch("proc_type_origin", &proc_type_origin_[0],
                               "proc_type_origin[npart]/I");
      collisions_tree_->Branch("pdg_mother1", &pdg_mother1_[0],
                               "pdg_mother1[npart]/I");
      collisions_tree_->Branch("pdg_mother2", &pdg_mother2_[0],
                               "pdg_mother2[npart]/I");
      collisions_tree_->Branch("baryon_number", &baryon_number_[0],
                               "baryon_number[npart]/I");
      collisions_tree_->Branch("strangeness", &strangeness_[0],
                               "strangeness[npart]/I");
    }
  }
}

/**
 * RootOutput destructor. Writes root objects (here TTrees) to file and closes
 * it.
 */
RootOutput::~RootOutput() {
  // kOverwrite option prevents from writing extra TKey objects into root file
  root_out_file_->Write("", TObject::kOverwrite);
  root_out_file_->Close();
  std::filesystem::rename(filename_unfinished_, filename_);
}

void RootOutput::at_eventstart(const Particles &particles,
                               const EventLabel &event_label,
                               const EventInfo &event) {
  // save event and ensemble numbers
  current_event_ = event_label.event_number;
  current_ensemble_ = event_label.ensemble_number;

  modus_l_ = event.modus_length;
  test_p_ = event.test_particles;
  current_t_ = event.current_time;
  E_kinetic_tot_ = event.total_kinetic_energy;
  E_fields_tot_ = event.total_mean_field_energy;
  E_tot_ = event.total_energy;

  if (write_particles_ && particles_only_final_ == OutputOnlyFinal::No) {
    output_counter_ = 0;
    // This is to have only one output of positive impact parameter per event
    impact_b_ = -1.0;
    empty_event_ = false;
    particles_to_tree(particles);
    output_counter_++;
  }
}

void RootOutput::at_intermediate_time(const Particles &particles,
                                      const std::unique_ptr<Clock> &,
                                      const DensityParameters &,
                                      const EventLabel &event_label,
                                      const EventInfo &event) {
  current_ensemble_ = event_label.ensemble_number;
  modus_l_ = event.modus_length;
  test_p_ = event.test_particles;
  current_t_ = event.current_time;
  E_kinetic_tot_ = event.total_kinetic_energy;
  E_fields_tot_ = event.total_mean_field_energy;
  E_tot_ = event.total_energy;

  if (write_particles_ && particles_only_final_ == OutputOnlyFinal::No) {
    particles_to_tree(particles);
    output_counter_++;
  }
}

void RootOutput::at_eventend(const Particles &particles,
                             const EventLabel &event_label,
                             const EventInfo &event) {
  current_ensemble_ = event_label.ensemble_number;
  test_p_ = event.test_particles;
  modus_l_ = event.modus_length;
  current_t_ = event.current_time;
  impact_b_ = event.impact_parameter;
  empty_event_ = event.empty_event;

  E_kinetic_tot_ = event.total_kinetic_energy;
  E_fields_tot_ = event.total_mean_field_energy;
  E_tot_ = event.total_energy;

  if (write_particles_ &&
      !(event.empty_event &&
        particles_only_final_ == OutputOnlyFinal::IfNotEmpty)) {
    particles_to_tree(particles);
  }
  /* Calculate the autosave frequency to be used. In case multiple ensembles are
   * used this has to be at least 1, so set it to 1 if the number of ensembles
   * is larger than the initial value of the autosave frequency. The evaluation
   * of the actual frequency to be used has to be done only once, i.e. when the
   * autosave_frequency_ has still its meaningless initial negative value. */
  if (autosave_frequency_ < 0) {
    autosave_frequency_ = 1000 / event.n_ensembles;
    if (autosave_frequency_ == 0) {
      autosave_frequency_ = 1;
    }
  }
  /* Forced regular dump from operational memory to disk. Very demanding!
   * If program crashes written data will NOT be lost. */
  if (current_event_ > 0 && current_event_ % autosave_frequency_ == 0) {
    if (write_particles_ || write_initial_conditions_) {
      particles_tree_->AutoSave("SaveSelf");
    }
    if (write_collisions_) {
      collisions_tree_->AutoSave("SaveSelf");
    }
  }
}

void RootOutput::at_interaction(const Action &action,
                                const double /*density*/) {
  if (write_collisions_) {
    collisions_to_tree(action.incoming_particles(), action.outgoing_particles(),
                       action.get_total_weight(), action.get_partial_weight());
  }

  if (write_initial_conditions_ &&
      (action.get_type() == ProcessType::Fluidization ||
       action.get_type() == ProcessType::FluidizationNoRemoval)) {
    particles_to_tree(action.incoming_particles());
  }
}

template <typename T>
void RootOutput::particles_to_tree(T &particles) {
  int i = 0;

  ev_ = current_event_;
  ens_ = current_ensemble_;
  tcounter_ = output_counter_;
  bool exceeded_buffer_message = true;

  for (const auto &p : particles) {
    // Buffer full - flush to tree, else fill with particles
    if (i >= max_buffer_size_) {
      if (exceeded_buffer_message) {
        logg[LOutput].warn()
            << "\nThe number of particles N = " << particles.size()
            << " exceeds the maximum buffer size B = " << max_buffer_size_
            << ".\nceil(N/B) = "
            << std::ceil(particles.size() /
                         static_cast<double>(max_buffer_size_))
            << " separate ROOT Tree entries will be created at this output."
            << "\nMaximum buffer size (max_buffer_size_) can be changed in "
            << "rootoutput.h\n\n";
        exceeded_buffer_message = false;
      }
      npart_ = max_buffer_size_;
      i = 0;
      particles_tree_->Fill();
    } else {
      id_[i] = p.id();
      pdgcode_[i] = p.pdgcode().get_decimal();
      charge_[i] = p.type().charge();
      formation_time_[i] = p.formation_time();
      time_last_collision_[i] = p.get_history().time_last_collision;

      p0_[i] = p.momentum().x0();
      px_[i] = p.momentum().x1();
      py_[i] = p.momentum().x2();
      pz_[i] = p.momentum().x3();

      t_[i] = p.position().x0();
      x_[i] = p.position().x1();
      y_[i] = p.position().x2();
      z_[i] = p.position().x3();

      if (part_extended_ || ic_extended_) {
        const auto h = p.get_history();
        coll_per_part_[i] = h.collisions_per_particle;
        xsec_factor_[i] = p.xsec_scaling_factor();
        proc_id_origin_[i] = h.id_process;
        proc_type_origin_[i] = static_cast<int>(h.process_type);
        pdg_mother1_[i] = h.p1.get_decimal();
        pdg_mother2_[i] = h.p2.get_decimal();
        baryon_number_[i] = p.type().baryon_number();
        strangeness_[i] = p.type().strangeness();
      }

      i++;
    }
  }
  // Flush rest to tree
  if (i > 0) {
    npart_ = i;
    particles_tree_->Fill();
  }
}

void RootOutput::collisions_to_tree(const ParticleList &incoming,
                                    const ParticleList &outgoing,
                                    const double weight,
                                    const double partial_weight) {
  ev_ = current_event_;
  ens_ = current_ensemble_;
  nin_ = incoming.size();
  nout_ = outgoing.size();
  npart_ = nin_ + nout_;
  wgt_ = weight;
  par_wgt_ = partial_weight;

  int i = 0;

  /* It is assumed that nin + nout < max_buffer_size_
   * This is true for any possible reaction for current buffer size
   * But if one wants initial/final particles written to collisions
   * then implementation should be updated. */

  for (const ParticleList &plist : {incoming, outgoing}) {
    for (const auto &p : plist) {
      id_[i] = p.id();
      pdgcode_[i] = p.pdgcode().get_decimal();
      charge_[i] = p.type().charge();
      formation_time_[i] = p.formation_time();
      time_last_collision_[i] = p.get_history().time_last_collision;

      p0_[i] = p.momentum().x0();
      px_[i] = p.momentum().x1();
      py_[i] = p.momentum().x2();
      pz_[i] = p.momentum().x3();

      t_[i] = p.position().x0();
      x_[i] = p.position().x1();
      y_[i] = p.position().x2();
      z_[i] = p.position().x3();

      if (coll_extended_) {
        const auto h = p.get_history();
        coll_per_part_[i] = h.collisions_per_particle;
        xsec_factor_[i] = p.xsec_scaling_factor();
        proc_id_origin_[i] = h.id_process;
        proc_type_origin_[i] = static_cast<int>(h.process_type);
        pdg_mother1_[i] = h.p1.get_decimal();
        pdg_mother2_[i] = h.p2.get_decimal();
        baryon_number_[i] = p.type().baryon_number();
        strangeness_[i] = p.type().strangeness();
      }

      i++;
    }
  }

  collisions_tree_->Fill();
}
}  // namespace smash
