/*
 *
 *    Copyright (c) 2014-2018
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

/*!\Userguide
 * \page format_root ROOT format
 * SMASH ROOT output has the same functionality as OSCAR output, but ROOT
 * files are faster to read and write and they need less disk space for the
 * same amount of information. This is achieved due to an optimized internal
 * structure of ROOT files (and compression). ROOT files are not
 * human-readable, but they can be viewed using ROOT's TBrowser. One can also
 * access them using ROOT functions. The full memory structure of the ROOT
 * files can be found here: http://root.cern.ch/root/html/TFile.html. We only
 * desribe the logical structure of the SMASH ROOT output. Knowing the logical
 * structure is enough to read and write ROOT files and understand their view
 * in TBrowser.
 *
 * Producing ROOT output requires ROOT installed (see http://root.cern.ch).
 *
 * SMASH produces one ROOT file per run: \c smash_run.root. This file contains
 * a TTree called \c particles and a TTree called \c collisions, depending
 * on the required content (see \ref output_general_). The \c particles
 * tree contains the same information as OSCAR particles output and the
 * \c collisions tree contains the same information as OSCAR collision output.
 *
 * In case that the ROOT format is used for dilepton output
 * (see \ref input_dileptons), the file is called \c Dileptons.root and
 * only contains a \c collisions tree with all the dilepton decays.
 *
 * Every physical quantity is in a separate TBranch.
 * One entry in the \c particles TTree is:
 * \code
 * ev tcounter npart impact_b pdgcode[npart] t[npart] x[npart] y[npart]
 * z[npart] p0[npart] px[npart] py[npart] pz[npart]
 * \endcode
 * One tree entry is analogous to an OSCAR output block, but the maximal
 * number of particles in one entry is limited to 10000. This is done to limit
 * the buffer size needed for ROOT output. If the number of particles in one
 * block exceeds 10000, then they are written in separate blocks with the same
 * \c tcounter and \c ev. The fields have the following meaning:
 *
 * \li \c ev is event number
 * \li \c tcounter is number of output block in a given event in terms of
 * OSCAR
 * \li \c npart is number of particles in the block
 * \li \c impact_b is the impact parameter of the event
 * \li \c pdgcode is PDG id array
 * \li \c t, \c x, \c y, \c z are position arrays
 * \li \c p0, \c px, \c py, \c pz are 4-momenta arrays
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
      particles_only_final_(out_par.part_only_final),
      autosave_frequency_(1000) {
  filename_unfinished_ = filename_;
  filename_unfinished_ += ".unfinished";
  root_out_file_ =
      make_unique<TFile>(filename_unfinished_.native().c_str(), "NEW");
  init_trees();
}

void RootOutput::init_trees() {
  if (write_particles_) {
    particles_tree_ = new TTree("particles", "particles");

    particles_tree_->Branch("npart", &npart, "npart/I");
    particles_tree_->Branch("impact_b", &impact_b, "impact_b/D");
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
    collisions_tree_->Branch("partial_weight", &par_wgt, "partial_weight/D");

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
                               const int event_number) {
  // save event number
  current_event_ = event_number;

  if (write_particles_ && !particles_only_final_) {
    output_counter_ = 0;
    // This is to have only one output of positive impact parameter per event
    impact_b = -1.0;
    particles_to_tree(particles);
    output_counter_++;
  }
}

void RootOutput::at_intermediate_time(const Particles &particles, const Clock &,
                                      const DensityParameters &) {
  if (write_particles_ && !particles_only_final_) {
    particles_to_tree(particles);
    output_counter_++;
  }
}

void RootOutput::at_eventend(const Particles &particles,
                             const int /*event_number*/,
                             double impact_parameter) {
  impact_b = impact_parameter;
  if (write_particles_) {
    particles_to_tree(particles);
  }
  /* Forced regular dump from operational memory to disk. Very demanding!
   * If program crashes written data will NOT be lost. */
  if (current_event_ > 0 && current_event_ % autosave_frequency_ == 0) {
    if (write_particles_) {
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
}

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
                                    const double weight,
                                    const double partial_weight) {
  ev = current_event_;
  nin = incoming.size();
  nout = outgoing.size();
  npart = nin + nout;
  wgt = weight;
  par_wgt = partial_weight;

  int i = 0;

  /* It is assumed that nin + nout < max_buffer_size_
   * This is true for any possible reaction for current buffer size: 10000
   * But if one wants initial/final particles written to collisions
   * then implementation should be updated. */

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
}  // namespace smash
