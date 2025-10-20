/*
 *
 *    Copyright (c) 2020-2024
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
// clang-format off
/*!\Userguide
 * \page doxypage_output_hepmc
 *
 * SMASH HepMC output is an implementation of the HepMC3 Event Record Library.
 * The aim is to provide a format compatible with other frameworks like Rivet
 * (https://rivet.hepforge.org). For resources regarding HepMC, see
 * http://hepmc.web.cern.ch/hepmc/ and \iref{Buckley:2019xhk}.
 *
 * The SMASH HepMC output can be:
 * - _HepMC_asciiv3_ plain human readable ASCII format
 * - _HepMC_treeroot_ ROOT Tree binary format, readable by ROOT
 *
 * You can find a snippet of the configuration for this output in \ref
 * config_output_examples.
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
 * \attention
 * At the moment **it is not possible to produce HepMC output with multiple
 * parallel ensembles** and trying to do so will result in an error.
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
 * - The HepMC output is best suited to \ref doxypage_input_conf_modi_collider modus, where
 *   the two ion "particles" are constructed and included as the initial incoming
 *   particles in the particle list, as well.
 * - In other SMASH modi (box, sphere, list, etc.) only the inital
 *   and final hadrons are written out as the incoming and outgoing particles
 *   of a single vertex.
 * - Even though in the HepMC library root and treeroot outputs are distinct, in
 *   SMASH the extension of the HepMC treeroot output is simply .root because
 *   the ROOT browser tool does not recognize the .treeroot extension.
 *
 * \section output_particles_collisions_ Particles and Collisions
 *
 * HepMC generally structures each event into particles and vertices
 * connecting them, basically storing a graph of the event.
 * However, in SMASH it is possible to filter the amount of information
 * available in the output, depending on whether the HepMC_asciiv3 or
 * HepMC_treeroot output options are specified under \key %Particles or
 * \key Collisions (see  \ref doxypage_input_conf_output for an example).
 * - \key %Particles: the output only provides a particle
 *   <b>list of the final state</b>. In HepMC only one central vertex is created.
 *   All initial state particles are incoming particles and all final state
 *   particles are outgoing particles of this vertex. Scatterings
 *   happening during the SMASH event are not recorded. For the collider
 *   modus, the initial state particles are combined into two single
 *   colliding nucleus "particles" with a nuclear pdg code.
 * - \key Collisions:  with this format, the <b>full event tree</b> is written.
 *   Furthermore, as in the previous \key %Particles case, in collider modus
 *   we lump all incoming nucleons into nuclei, but split them out immediately
 *   afterwards to allow tracking of the individual nucleons.
 *
 * \attention
 * - If the HepMC output is intended to be part of a Rivet analysis that requires
 * access to the parent particles of the final state particles, then the \key Collisions
 * output should be used.
 * - By default, SMASH exclusively handles strong decays at the end of the evolution.
 * If one intends to include all non-strong decays as well,
 * the \key Ignore_Minimum_Decay_Width_For_Decays_At_The_End option in the
 * \ref doxypage_input_conf_collision_term section of the configuration should be enabled.
 *
 * \section output_hepmc_asciiv3_ ASCII HepMC Format
 *
 * In this case the information about each event is inserted into a plain,
 * human readable ascii text file.
 *
 * In most cases the HepMC output is read by using tools like Rivet,
 * however it can be useful to have a basic knowledge about the most important
 * pieces of information contained there. We refer to the official HepMC
 * documentation for more details.
 *
 * Here there is an example of the first lines of the HepMC_asciiv3 output, in
 * \key Collider modus and \key Particles output type:

 \verbatim

 HepMC::Version 3.02.05
 HepMC::Asciiv3-START_EVENT_LISTING
 W Default
 T SMASH\|SMASH-2.2.1-81-g450bc31amaster\|
 E 0 1 772
 U GEV MM
 W 1.0000000000000000000000e+00
 A 0 GenHeavyIon v0 -1 -1 -1 -1 -1 -1 -1 -1 -1 1.6090284 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0
 P 1 0 1000791970 -1.0269562977782698e-15 1.3877787807814457e-16 3.0788744605039813e+02 3.6105628807654466e+02 1.8859205636552139e+02 4
 P 2 0 1000791970 -5.6898930012039273e-16 9.7144514654701197e-17 -3.0817189202123001e+02 3.6138985479939856e+02 1.8876628968114440e+02 4
 V -1 0 [1,2]
 P 3 -1 2112 -1.3955349705352910e-01 -3.2052252367927581e-01 -1.3070952483948553e-01 1.0095240693561323e+00 9.3799999999999994e-01 1
 P 4 -1 2112 2.2504221328938548e-02 -8.1112446752570233e-02 1.8331677117778440e+00 2.0609302580389826e+00 9.3799999999999994e-01 1
 P 5 -1 111 -1.7932505641335714e-01 -2.3234517783752556e-01 7.9030334352518175e-02 3.3381364751890524e-01 1.3800000000000001e-01 1
 P 6 -1 -211 5.3428615692750892e-01 9.5409650270284543e-02 -1.6651817800882954e-01 5.8424053475982385e-01 1.3800000000000001e-01 1

 \endverbatim
 *
 * The first 4 lines provide information about the version of HepMC3 and
 * of the generator (i.e. SMASH) which have been used. These lines are not
 * repeated.
 *
 * Row 5: "E 0 1 772" marks the beginning of a new event, with the numbers
 * referring to the event ID, the number of vertexes and the number of
 * particles, respectively, which have been recorded.
 *
 * The content of line 8 _"A 0 GenHeavyIon ..."_ is described in
 * Appendix A.3.3 of \iref{Buckley:2019xhk}, however SMASH
 * prints only the value of the impact parameter (1.6090284 in the example).
 * This line closes the header, from now on, until the record of the next
 * event, there are only lines starting with P (for particle data) or V
 * (vertex data).
 *
 * Lines 9 and 10: these are the data regarding the two colliding ions,
 * which are treated like two big particles, with PDG code based on their
 * atomic and mass number (79 and 197, respectively, in this example with
 *  gold).
 *
 * Row 11: V, i.e. vertex. The first number is the ID of the vertex, which
 * decreases by 1 at each interaction. The square bracket contains the IDs
 * of the incoming particles.
 *
 * Lines from 12: P, i.e. particle. The contents of the columns are:
 * - "P", identifying a particle data line
 * - the internal particle ID
 * - the ID of the originating vertex
 * - the particle PDG ID
 * - particle momentum along x
 * - particle momentum along y
 * - particle momentum along z
 * - particle energy
 * - particle rest mass
 * - interaction code ( 1: final state, 2: decay, 4: beam particle)
 *
 * The HepMC3_asciiv3 with output under \key Collisions presents has a similar
 * structure, but, in addition to containing all the interaction vertexes and
 * the non final particles, it has the following differences after the header:
 * - a list of lines starting with "A", like
 *   "A -6311 partial_weight 99.6726",
 *   where the first number refers to a vertex ID and the partial_weights are
 *   cross sections, followed by a similar list with lines like
 *   "A -6311 weight 100.876",
 *   where the weights are total cross sections
 * - a list of particles and vertexes, in which the interaction code now
 *   can be also one of those used by SMASH, shifted by 100 (e.g. 103 or 145)
 *
 * In both cases, the HepMC output file ends with the line:
 * _HepMC::Asciiv3-END_EVENT_LISTING_, followed by an empty line.
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
 * HepMC3 output of a subtype (asciiv3 or treeroot) into the other without
 * loss of information.
 *
 **/

// clang-format on
HepMcOutput::HepMcOutput(const std::filesystem::path &path, std::string name,
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
  std::filesystem::rename(filename_unfinished_, filename_);
}

void HepMcOutput::at_eventend(const Particles &particles,
                              const EventLabel &event_label,
                              const EventInfo &event) {
  HepMcInterface::at_eventend(particles, event_label, event);
  logg[LOutput].debug() << "Writing event " << event_label.event_number
                        << " with " << event_.particles().size()
                        << " particles and " << event_.vertices().size()
                        << " vertices to output " << std::endl;
  output_file_->write_event(event_);
}

}  // namespace smash
