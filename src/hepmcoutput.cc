/*
 *
 *    Copyright (c) 2020-2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/hepmcoutput.h"

#include "HepMC3/Print.h"
#include "HepMC3/WriterAscii.h"

#ifdef SMASH_USE_HEPMC_ROOTIO
#include "HepMC3/WriterRootTree.h"
#endif

namespace smash {

/*!\Userguide
 * \page output_hepmc_ HepMC Output
 *
 * SMASH HepMC output is an implementation of the HepMC3 Event Record Library.
 * The aim is to provide a format compatible with other frameworks like Rivet
 * (https://rivet.hepforge.org). For resources regarding HepMC, see
 * http://hepmc.web.cern.ch/hepmc/ and https://arxiv.org/abs/1912.08005.
 *
 * The SMASH HepMC output can be:
 * - _HepMC_asciiv3_ plain human readable ASCII format
 * - _HepMC_treeroot_ ROOT Tree binary format, readable by ROOT
 *
 * You can find a snippet of the configuration for this output in \ref
 * configuring_output_.
 *
 * Producing HepMC output in asciiv3 format requires HepMC3 to be installed.
 * Download the tarball from http://hepmc.web.cern.ch/hepmc/
 * and follow the instructions or use the pre-compiled packages for your
 * OS distribution.
 * If the user wants to produce HepMC output in ROOT Tree format,
 * ROOT must be installed (https://root.cern.ch/), as well, and HepMC should
 * be compiled with ROOT IO support (`-DHEPMC3_ENABLE_ROOTIO:BOOL=ON`) or,
 * when using a binary precompiled distribution, the appropriate rootIO
 * package must be installed.
 *
 * \note
 * - Since some HepMC readers (e.g. Rivet) need a value for the
 * nuclei-nuclei cross section, a dummy cross section of 1.0 is written to the
 * output.
 * - To avoid confusion with the definition of these quantities within the
 * Glauber model, in the header of an event we set the numbers of participants
 * and of collisions to -1.
 * - If you use Fermi motion and want to read in the HepMC
 * ouput into Rivet, you need to disable the check for the beam particle
 * energies with the \key --ignore-beams option. When using the Rivet output
 * this check is disabled by default.
 * - The HepMC output is best suited to \ref input_modi_collider_ modus, where
 *   the two ion "particles" are constructed and included as the intial incoming
 *   particles in the particle list, as well.
 * - In other SMASH modi (box, sphere, list, etc.) only the inital
 *   and final hadrons are written out as the incoming and outgoing particles
 *   of a single vertex.
 * - Even though in the HepMC library root and treeroot outputs are distinct, in
 *   SMASH the extension of the HepMC treeroot output is simply .root because
 *   the ROOT browser tool does not recognize the .treeroot extension.
 *
 * \section output_hepmc_asciiv3_ ASCII HepMC Format
 *
 * HepMC generally structures each event into particles and vertices
 * connecting them, basically storing a graph of the event.

 * Two versions of HepMC are possible by specifying HepMC_asciiv3 or
 * HepMC_treeroot under \key Particles. The output only provides a particle
 * list of the final state. For this only one central vertex is used. All
 * initial state particles are incoming particles and all final state
 * particles are outgoing particles of this vertex. Scatterings
 * happening during the SMASH event are not recorded. For the collider
 * modus, the intial state particles are combined into two single
 * colliding nucleus "particles" with a nuclear pdg code.
 *
 * In certain circumstances one may be interested in the full
 * event structure.  To that end, the output module provides the
 * HepMC format also for the \key Collisions content. With this format, the
 * full event tree is written.  Furthermore, as above, for collider modus,
 * we lump all incoming nucleons into nuclei, but split them out immediately
 * afterwards to allow tracking of the individual nucleons.
 *
 * In all cases, the module will track the number of binary, inelastic
 * collisions between incident nucleons as well as keep track of
 * participating incident nucleons.
 *
 *
 * \section output_hepmc_root_ ROOT HepMC Format
 *
 * In this case the information about each event is inserted into a ROOT
 * Tree structure and saved in a binary file that can be read by ROOT.
 *
 * Here is an example of a basic ROOT macro that displays the structure of
 * the tree, assuming that the output file is in the same directory as the
 * macro:
 * \code
 * // file name: read_treeroot.C
 * #include <TFile.h>
 * #include <TTree.h>
 *
 * int read_hempc3_treeroot() {
 * // Opens a SMASH HepMC3 treeroot output file to be read in
 * TFile *input_file = TFile::Open("SMASH_HepMC_particles.root","read");
 * if (input_file->IsOpen()) {
 *   printf("Successfully opened file %s\n", input_file->GetName());
 * } else {
 *   printf("Error at opening file %s\n", input_file->GetName());
 * }
 *
 * // Shows the top level contents of the file
 * input_file->ls();
 *
 * // Gets a tree from file
 * TTree *tree = static_cast<TTree*>(input_file->Get("hepmc3_tree"));
 *
 * // Gets the number of entries (i.e. events) stored in the tree
 * Int_t nentries = tree->GetEntries();
 * std::cout << "\nThe number of entries in the tree is " << nentries <<
 std::endl;
 *
 * // Prints the branches of the tree
 * std::cout << "\n\nA bit more info:" << std::endl;
 *
 * tree->Print();
 *
 * // Extracts a bit more information
 * std::cout << "\n\nInfo about the branch hepmc3_event:" << std::endl;
 *
 * TBranch *b = tree->GetBranch("hepmc3_event");
 *
 * tree->Show(0);
 *
 * input_file->Close();
 *
 * return 0;
 * }
 * \endcode
 *
 * To load and execute the macro with ROOT:
 * \code
 * root
 * .L read_treeroot.C
 * read_hempc3_treeroot()
 * \endcode
 *
 * In the subdirectory `examples/reading_HepMC3_treeroot_output` under the
 * SMASH source code main directory there is a basic example of a C++ code
 * that reads the output without requiring HepMC3, but only ROOT.
 *
 * Among the examples of the HepMC3 library source code there is a converter
 * between different formats. The converter allows to transform the SMASH
 * HepMC3 output of a subtype (asciiv3 or treeroot) into the other.
 *
 **/
HepMcOutput::HepMcOutput(const bf::path &path, std::string name,
                         const bool full_event, std::string HepMC3_output_type)
    : HepMcInterface(name, full_event),
      filename_(path / (name + "." + HepMC3_output_type)) {
  filename_unfinished_ = filename_;
  filename_unfinished_ += +".unfinished";
#ifdef SMASH_USE_HEPMC_ROOTIO
  if (HepMC3_output_type == "asciiv3") {
#endif
    output_file_ = std::make_unique<HepMC3::WriterAscii>(
        filename_unfinished_.string(), event_.run_info());
    output_type_ = asciiv3;
#ifdef SMASH_USE_HEPMC_ROOTIO
  } else {
    output_file_ = std::make_unique<HepMC3::WriterRootTree>(
        filename_unfinished_.string(), event_.run_info());
    output_type_ = treeroot;
  }
#endif
}
HepMcOutput::~HepMcOutput() {
  logg[LOutput].debug() << "Renaming file " << filename_unfinished_ << " to "
                        << filename_ << std::endl;
  output_file_->close();
  bf::rename(filename_unfinished_, filename_);
}

void HepMcOutput::at_eventend(const Particles &particles,
                              const int32_t event_number,
                              const EventInfo &event) {
  HepMcInterface::at_eventend(particles, event_number, event);
  logg[LOutput].debug() << "Writing event " << event_number << " with "
                        << event_.particles().size() << " particles and "
                        << event_.vertices().size() << " vertices to output "
                        << std::endl;
  output_file_->write_event(event_);
}

}  // namespace smash
