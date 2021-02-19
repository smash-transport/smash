/*
 *
 *    Copyright (c) 2014-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/hepmcoutput.h"

#include "HepMC3/Print.h"

namespace smash {

/*!\Userguide
 * \page hepmc_output_user_guide_ HepMC Output
 * SMASH HepMC output is an implementation of the HepMC3 Event Record Library.
 * The aim is to provide a format compatible with other frameworks like Rivet
 * (https://rivet.hepforge.org). For resources regarding HepMC, see
 * http://hepmc.web.cern.ch/hepmc/ and https://arxiv.org/abs/1912.08005.
 *
 * Producing HepMC output requires HepMC3 to be installed. Download the tarball
 * from http://hepmc.web.cern.ch/hepmc/.
 *
 * \section hepmc_output_user_guide_format_ ASCII HepMC Format
 *
 * HepMC generaly structures each event into particles and vertices
 * connecting them, basically storing a graph of the event. Since the
 * purpose of this output is to only provide a particle list of the
 * final state produced the output format is adapted accordingly for
 * the SMASH implementation: Only one central vertex is used. All
 * intial state particles are incoming particles and all final state
 * particles are outgoing particles of this vertex. Scatterings
 * happening during the SMASH event are not recorded. For the collider
 * modus, the intial state particles are combined into two single
 * colliding nucleus "particles" with a nuclear pdg code.
 *
 * This reduced event format is selected by the \key Format value \key ASCII.
 *
 * However, in certain circumstances one may be interested in the full
 * event structure.  To that end, the output module provides the
 * format \key ASCII-full.  With this format, the full event tree is
 * written.  Furthermore, as above, for collider modus, we lump all
 * incoming nucleons into nuclei, but split them out immediately
 * afterwards to allow tracking of the individual nucleons.
 *
 * In all cases, the module will track the number of binary, inelastic
 * collisions between incident nucleons as well as keep track of
 * participating incident nucleons.
 *
 * The output is written in Asciiv3, the HepMC3 native plain text format. See
 * https://arxiv.org/abs/1912.08005 for documentation of the format.
 *
 * \note Since some HepMC readers (e.g. Rivet) need a value for the
 * nuclei-nuclei cross section, a dummy cross section of 1.0 is written to the
 * output. Furthermore, if you use Fermi motion and want to read in the HepMC
 * ouput into Rivet, you need to disable the check for the beam particle
 * energies with the \key --ignore-beams option.
 */
HepMcOutput::HepMcOutput(const bf::path &path, std::string name,
                         const OutputParameters &out_par, const int total_N,
                         const int proj_N)
    : Base(name, out_par, total_N, proj_N),
      filename_(path / (name + ".asciiv3")) {
  filename_unfinished_ = filename_;
  filename_unfinished_ += +".unfinished";
  output_file_ =
      make_unique<HepMC3::WriterAscii>(filename_unfinished_.string());
}
HepMcOutput::~HepMcOutput() {
  logg[LOutput].debug() << "Renaming file " << filename_unfinished_ << " to "
                        << filename_ << std::endl;
  bf::rename(filename_unfinished_, filename_);
}

void HepMcOutput::at_eventend(const Particles &particles,
                              const int32_t event_number,
                              const EventInfo &event) {
  Base::at_eventend(particles, event_number, event);
  logg[LOutput].debug() << "Writing event " << event_number << " with "
                        << event_.particles().size() << " particles and "
                        << event_.vertices().size() << " vertices to output "
                        << std::endl;
  output_file_->write_event(event_);
}

}  // namespace smash
