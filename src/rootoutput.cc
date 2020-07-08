/*
 *
 *    Copyright (c) 2014-2020
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
static constexpr int LHyperSurfaceCrossing = LogArea::HyperSurfaceCrossing::id;
static constexpr int LOutput = LogArea::Output::id;

/*!\Userguide
 * \page format_root ROOT Format
 * SMASH ROOT output shares functionalities with the OSCAR output, but ROOT
 * files are faster to read and write and they need less disk space for the
 * same amount of information. This is achieved due to an optimized internal
 * structure of ROOT files (and compression). ROOT files are not
 * human-readable, but they can be viewed using ROOT's TBrowser. One can also
 * access them using ROOT functions. The full memory structure of the ROOT
 * files can be found here: http://root.cern.ch/root/html/TFile.html. We only
 * desribe the logical structure of the SMASH ROOT output. Knowing the logical
 * structure is enough to read and write ROOT files, be able to view them in
 * TBrowser, or write a ROOT macro to analyze them.
 *
 * Producing ROOT output requires ROOT installed (see http://root.cern.ch).
 *
 * SMASH produces one ROOT file per run: \c smash_run.root. This file contains
 * a TTree called \c particles and a TTree called \c collisions, depending
 * on the required content (see \ref output_general_). The \c particles
 * tree contains information about the parameters of the run (such as the number
 * of testparticles and event number), information relating to individual
 * particles (such as their position or charge), and information about bulk
 * observables in the system (kinetic energy, mean field energy, and total
 * energy). The \c collisions tree contains the same information as OSCAR
 * collision output.
 *
 * In case that the ROOT format is used for dilepton output
 * (see \ref output_dileptons), the ROOT file is called \c Dileptons.root and
 * only contains a \c collisions tree with all the dilepton decays.
 *
 * Every physical quantity corresponds to a separate TBranch.
 * One entry in the \c particles TTree is:
 * \code
 * ev tcounter npart test_p modus_l current_t impact_b empty_event
 * pdgcode[npart] charge[npart] t[npart] x[npart] y[npart] z[npart] p0[npart]
 * px[npart] py[npart] pz[npart] E_kinetic_tot E_fields_tot E_tot
 * \endcode
 * One tree entry is analogous to an OSCAR output block, but the maximal
 * number of particles in one entry is limited to 500000. This is done to limit
 * the buffer size needed for ROOT output. If the number of particles in one
 * block exceeds 500000, then they are written in separate blocks with the same
 * \c tcounter and \c ev. The fields have the following meaning:
 *
 * \li \c ev is event number
 * \li \c tcounter is number of output block in a given event in terms of
 * OSCAR
 * \li \c npart is number of particles in the block
 * \li \c test_p is number of testparticles per particle
 * \li \c modus_l is modus length
 * \li \c current_t is time associated with the output block, in fm/c
 * \li \c impact_b is the impact parameter of the event
 * \li \c empty_event indicates whether the projectile did not interact with the
 * target
 * \li \c pdgcode is PDG id array
 * \li \c charge is the electric charge array
 * \li \c p0, \c px, \c py, \c pz are 4-momenta arrays
 * \li \c t, \c x, \c y, \c z are position arrays
 * \li \c E_kinetic_tot is total kinetic energy in the system
 * \li \c E_fields_tot is total mean field energy * test_p
 * \li \c E_total is the sum of E_kinetic_tot and E_fields_tot
 *
 * In case of extended output (see \ref output_content_specific_options_) more
 * fields are added. Their description is the same that in case of OSCAR
 * format, see \ref extended_output_format_.
 *
 * The entries in the \c collisions tree are organized in the same way, but
 * a few additional fields are present:
 * \li \c nin and \c nout are added to characterize number of incoming and
 *     outgoing particles in the reaction, with nin + nout = npart.
 * \li \c weight is an action weight, whose meaning depends on the type of
 *     action: For collisions it is the total cross section, for decays it is
 *     the total decay width and for dilepton decays it is the shining weight.
 *
 * Currently writing initial and final configuration to collisions tree is
 * not supported.
 *
 * See also \ref collisions_output_in_box_modus_.
 */
RootOutput::RootOutput(const bf::path &path, const std::string &name,
                       const OutputParameters &out_par)
    : OutputInterface(name),
      filename_(path / (name + ".root")),
      write_collisions_(name == "Collisions" || name == "Dileptons" ||
                        name == "Photons"),
      write_particles_(name == "Particles"),
      write_initial_conditions_(name == "SMASH_IC"),
      particles_only_final_(out_par.part_only_final),
      autosave_frequency_(1000),
      part_extended_(out_par.part_extended),
      coll_extended_(out_par.coll_extended),
      ic_extended_(out_par.ic_extended) {
  filename_unfinished_ = filename_;
  filename_unfinished_ += ".unfinished";
  root_out_file_ =
      make_unique<TFile>(filename_unfinished_.native().c_str(), "NEW");
  resize_vector(pdgcode_);
  resize_vector(charge_);
  resize_vector(p0_);
  resize_vector(px_);
  resize_vector(py_);
  resize_vector(pz_);
  resize_vector(t_);
  resize_vector(x_);
  resize_vector(y_);
  resize_vector(z_);
  resize_vector(formation_time_);
  resize_vector(xsec_factor_);
  resize_vector(time_last_coll_);
  resize_vector(coll_per_part_);
  resize_vector(proc_id_origin_);
  resize_vector(proc_type_origin_);
  resize_vector(pdg_mother1_);
  resize_vector(pdg_mother2_);
  init_trees();
}

void RootOutput::init_trees() {
  if (write_particles_ || write_initial_conditions_) {
    particles_tree_ = new TTree("particles", "particles");

    particles_tree_->Branch("ev", &ev_, "ev/I");
    particles_tree_->Branch("tcounter", &tcounter_, "tcounter/I");
    particles_tree_->Branch("npart", &npart_, "npart/I");
    particles_tree_->Branch("test_p", &test_p_, "test_p/I");
    particles_tree_->Branch("modus_l", &modus_l_, "modus_l/D");
    particles_tree_->Branch("current_t", &current_t_, "current_t/D");
    particles_tree_->Branch("impact_b", &impact_b_, "impact_b/D");
    particles_tree_->Branch("empty_event", &empty_event_, "empty_event/O");

    particles_tree_->Branch("pdgcode", &pdgcode_[0], "pdgcode[npart]/I");
    particles_tree_->Branch("charge", &charge_[0], "charge[npart]/I");

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
      particles_tree_->Branch("form_time", &formation_time_[0],
                              "form_time[npart]/D");
      particles_tree_->Branch("xsecfac", &xsec_factor_[0], "xsecfac[npart]/D");
      particles_tree_->Branch("proc_id_origin", &proc_id_origin_[0],
                              "proc_id_origin[npart]/I");
      particles_tree_->Branch("proc_type_origin", &proc_type_origin_[0],
                              "proc_type_origin[npart]/I");
      particles_tree_->Branch("time_last_coll", &time_last_coll_[0],
                              "time_last_coll[npart]/D");
      particles_tree_->Branch("pdg_mother1", &pdg_mother1_[0],
                              "pdg_mother1[npart]/I");
      particles_tree_->Branch("pdg_mother2", &pdg_mother2_[0],
                              "pdg_mother2[npart]/I");
    }
  }

  if (write_collisions_) {
    collisions_tree_ = new TTree("collisions", "collisions");

    collisions_tree_->Branch("nin", &nin_, "nin/I");
    collisions_tree_->Branch("nout", &nout_, "nout/I");
    collisions_tree_->Branch("npart", &npart_, "npart/I");
    collisions_tree_->Branch("ev", &ev_, "ev/I");
    collisions_tree_->Branch("weight", &wgt_, "weight/D");
    collisions_tree_->Branch("partial_weight", &par_wgt_, "partial_weight/D");

    collisions_tree_->Branch("pdgcode", &pdgcode_[0], "pdgcode[npart]/I");
    collisions_tree_->Branch("charge", &charge_[0], "charge[npart]/I");

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
      collisions_tree_->Branch("form_time", &formation_time_[0],
                               "form_time[npart]/D");
      collisions_tree_->Branch("xsecfac", &xsec_factor_[0], "xsecfac[npart]/D");
      collisions_tree_->Branch("proc_id_origin", &proc_id_origin_[0],
                               "proc_id_origin[npart]/I");
      collisions_tree_->Branch("proc_type_origin", &proc_type_origin_[0],
                               "proc_type_origin[npart]/I");
      collisions_tree_->Branch("time_last_coll", &time_last_coll_[0],
                               "time_last_coll[npart]/D");
      collisions_tree_->Branch("pdg_mother1", &pdg_mother1_[0],
                               "pdg_mother1[npart]/I");
      collisions_tree_->Branch("pdg_mother2", &pdg_mother2_[0],
                               "pdg_mother2[npart]/I");
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
  bf::rename(filename_unfinished_, filename_);
}

void RootOutput::at_eventstart(const Particles &particles,
                               const int event_number, const EventInfo &event) {
  // save event number
  current_event_ = event_number;

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
                                      const EventInfo &event) {
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
                             const int /*event_number*/,
                             const EventInfo &event) {
  modus_l_ = event.modus_length;
  test_p_ = event.test_particles;
  current_t_ = event.current_time;
  E_kinetic_tot_ = event.total_kinetic_energy;
  E_fields_tot_ = event.total_mean_field_energy;
  E_tot_ = event.total_energy;

  impact_b_ = event.impact_parameter;
  empty_event_ = event.empty_event;
  if (write_particles_ &&
      !(event.empty_event &&
        particles_only_final_ == OutputOnlyFinal::IfNotEmpty)) {
    particles_to_tree(particles);
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

  if (write_initial_conditions_) {
    // If the runtime is too short some particles might not yet have
    // reached the hypersurface. Warning is printed.
    if (particles.size() != 0) {
      logg[LHyperSurfaceCrossing].warn(
          "End time might be too small for initial conditions output. "
          "Hypersurface has not yet been crossed by ",
          particles.size(), " particle(s).");
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
      action.get_type() == ProcessType::HyperSurfaceCrossing) {
    particles_to_tree(action.incoming_particles());
  }
}

template <typename T>
void RootOutput::particles_to_tree(T &particles) {
  int i = 0;

  ev_ = current_event_;
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
            << ceil(particles.size() / static_cast<double>(max_buffer_size_))
            << " separate ROOT Tree entries will be created at this output."
            << "\nMaximum buffer size (max_buffer_size_) can be changed in "
            << "rootoutput.h\n\n";
        exceeded_buffer_message = false;
      }
      npart_ = max_buffer_size_;
      i = 0;
      particles_tree_->Fill();
    } else {
      pdgcode_[i] = p.pdgcode().get_decimal();
      charge_[i] = p.type().charge();

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
        formation_time_[i] = p.formation_time();
        xsec_factor_[i] = p.xsec_scaling_factor();
        time_last_coll_[i] = h.time_last_collision;
        coll_per_part_[i] = h.collisions_per_particle;
        proc_id_origin_[i] = h.id_process;
        proc_type_origin_[i] = static_cast<int>(h.process_type);
        pdg_mother1_[i] = h.p1.get_decimal();
        pdg_mother2_[i] = h.p2.get_decimal();
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
  nin_ = incoming.size();
  nout_ = outgoing.size();
  npart_ = nin_ + nout_;
  wgt_ = weight;
  par_wgt_ = partial_weight;

  int i = 0;

  /* It is assumed that nin + nout < max_buffer_size_
   * This is true for any possible reaction for current buffer size: 10000
   * But if one wants initial/final particles written to collisions
   * then implementation should be updated. */

  for (const ParticleList &plist : {incoming, outgoing}) {
    for (const auto &p : plist) {
      pdgcode_[i] = p.pdgcode().get_decimal();
      charge_[i] = p.type().charge();

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
        formation_time_[i] = p.formation_time();
        xsec_factor_[i] = p.xsec_scaling_factor();
        time_last_coll_[i] = h.time_last_collision;
        coll_per_part_[i] = h.collisions_per_particle;
        proc_id_origin_[i] = h.id_process;
        proc_type_origin_[i] = static_cast<int>(h.process_type);
        pdg_mother1_[i] = h.p1.get_decimal();
        pdg_mother2_[i] = h.p2.get_decimal();
      }

      i++;
    }
  }

  collisions_tree_->Fill();
}
}  // namespace smash
