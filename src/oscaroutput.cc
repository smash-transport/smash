/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/oscaroutput.h"

#include <boost/filesystem.hpp>
#include <string>

#include "include/clock.h"
#include "include/config.h"
#include "include/configuration.h"
#include "include/cxx14compat.h"
#include "include/forwarddeclarations.h"
#include "include/particles.h"

namespace Smash {

template <OscarOutputFormat Format, int Contents>
OscarOutput<Format, Contents>::OscarOutput(bf::path path, std::string name)
    : file_{std::fopen((path / (name + ".oscar")).native().c_str(), "w")} {
  /*!\Userguide
   * \page input_oscar_particlelist Oscar_Particlelist
   * Enables OSCAR particles output.
   * OSCAR particles output provides the particle list at the output intervals. The text format
   * is either OSCAR1999 or OSCAR2013, this is controlled by an option.
   * Fixed moments of output can be: event start, event end, every next
   * output interval \f$\Delta t \f$.
   * Writing (or not writing) output at these moments is controlled by options.
   * Output time interval \f$\Delta t \f$ is also regulated by an option.
   *
   * \key Enable (bool, optional, default = false):\n
   * true - OSCAR particle list output enabled\n
   * false - no OSCAR particle list output
   *
   * \key 2013_Format (bool, optional, default = false): \n
   * true - output will be in OSCAR2013 format\n
   * false - output will be in OSCAR1999 format
   *
   * \key Only_Final (bool, optional, default = true): \n
   * true - print only final particle list \n
   * false - particle list at output interval including initial time
   *
   *  Detailed specification of OSCAR particle list format can be found here:
   * \ref format_oscar_particlelist
   *
   * \page input_oscar_collisions Oscar_Collisions
   * Enables OSCAR collisions output. The latter saves information about
   * every collision, decay and box wall crossing in OSCAR1999 or OSCAR2013 format.
   * Optionally initial and final particle configurations can be written out.
   *
   * \key Enable (bool, optional, default = false):\n
   * true - OSCAR collision output enabled
   * false - no OSCAR collision output
   *
   * \key 2013_Format (bool, optional, default = false): \n
   * true - output will be in OSCAR2013 format\n
   * false - output will be in OSCAR1999 format
   *
   * \key Print_Start_End (bool, optional, default = false): \n
   * true - initial and final particle list is written out \n
   * false - initial and final particle list is not written out
   *
   * Detailed specification of OSCAR collisions format can be found here:
   * \ref format_oscar_collisions
   */

  /*!\Userguide
   * \page oscar_general_ OSCAR block structure
   * OSCAR outputs are a family of ASCII and binary formats that follow
   * OSCAR format conventions. \n
   * **All OSCAR outputs have the same general structure: header and arbitrary
   * number of event blocks.** Each event block consists of arbitrary number of
   * output blocks and special event end line that marks the end of event. One
   * output block consists of output block header and N particle lines, N is
   * specified in the output block  header. \n
   * File structure can be visualized in the following way:
   * \code
   * Header
   * Event block 1
   *   output block 1
   *       output block header
   *       particle line 1
   *       particle line 2
   *       ...
   *       particle line N
   *   output block 2
   *   ...
   *   output block k
   *   event end line
   * Event block 2
   * ...
   * \endcode
   * To fully characterise any OSCAR output one has to specify the following
   * formatting:
   * \li header
   * \li output block header
   * \li particle line
   * \li event end line
   *
   * Every OSCAR output can produce two types of files: collisions output and
   * particles output. In collisions output file and in particles output file
   * the above structure is the same, but meaning of blocks is different.
   * **In collision file one output block typically corresponds to one collision
   * / decay / box wall crossing, while in particles output one block
   * corresponds to the current particle list at one moment of time.**
   * Particles output may contain the particle list at event start
   * immediately after initialization, at event end (which is reached when
   * time is larger or equal than \c End_Time in configuration file) and
   * periodically during evolution, period is defined by \c Output_Interval
   * option in configuration file, see \ref input_general_.
   * Collisions output contains all collisions / decays / box wall crossings
   * and optionally initial and final configuration.
   */

  if (Format == OscarFormat2013) {
    std::fprintf(file_.get(), "#!OSCAR2013 %s %s ", name.c_str(),
                 VERSION_MAJOR);
    std::fprintf(file_.get(), "t x y z mass p0 px py pz pdg ID\n");
    std::fprintf(file_.get(),
                 "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none\n");
  } else {
    if (name == "particle_lists") {
      name = "final_id_p_x";  // FIXME: why is this necessary? I.e. what does
                              // the string on the second line tell, and why
                              // does it have to be this specific string?
    }
    std::fprintf(file_.get(), "# OSC1999A\n# %s\n# %s\n", name.c_str(),
                 VERSION_MAJOR);
    std::fprintf(file_.get(), "# Block format:\n");
    std::fprintf(file_.get(), "# nin nout event_number\n");
    std::fprintf(file_.get(), "# id pdg 0 px py pz p0 mass x y z t\n");
    std::fprintf(file_.get(), "# End of event: 0 0 event_number\n");
    std::fprintf(file_.get(), "#\n");
  }
}

template <OscarOutputFormat Format, int Contents>
inline void OscarOutput<Format, Contents>::write(const Particles &particles) {
  for (const ParticleData &data : particles) {
    write_particledata(data);
  }
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_eventstart(const Particles &particles,
                                                  const int event_number) {
  if (Contents & OscarAtEventstart) {
    if (Format == OscarFormat2013) {
      std::fprintf(file_.get(), "# event %i in %zu\n", event_number + 1,
                   particles.size());
    } else {
      /* OSCAR line prefix : initial particles; final particles; event id
       * First block of an event: initial = 0, final = number of particles
       */
      const size_t zero = 0;
      std::fprintf(file_.get(), "%zu %zu %i\n", zero, particles.size(),
                   event_number + 1);
    }
    write(particles);
  }
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_eventend(const Particles &particles,
                                                const int event_number) {
  if (Format == OscarFormat2013) {
    if (Contents & OscarParticlesAtEventend) {
      std::fprintf(file_.get(), "# event %i out %zu\n", event_number + 1,
                   particles.size());
      write(particles);
    }
    // Comment end of an event
    std::fprintf(file_.get(), "# event %i end 0\n", event_number + 1);
  } else {
    // OSCAR line prefix : initial particles; final particles; event id
    // Last block of an event: initial = number of particles, final = 0
    // Block ends with null interaction
    const size_t zero = 0;
    if (Contents & OscarParticlesAtEventend) {
      std::fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero,
                   event_number + 1);
      write(particles);
    }
    // Null interaction marks the end of an event
    std::fprintf(file_.get(), "%zu %zu %i\n", zero, zero, event_number + 1);
  }
  // Flush to disk
  std::fflush(file_.get());
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_interaction(
    const ParticleList &incoming_particles,
    const ParticleList &outgoing_particles,
    const double density,
    const double total_cross_section,
    const ProcessType process_type) {
  if (Contents & OscarInteractions) {
    if (Format == OscarFormat2013) {
      std::fprintf(
          file_.get(),
          "# interaction in %zu out %zu rho %12.7f weight %12.7g type %5i \n",
          incoming_particles.size(), outgoing_particles.size(), density,
          total_cross_section, static_cast<int>(process_type));
    } else {
      /* OSCAR line prefix : initial final
       * particle creation: 0 1
       * particle 2<->2 collision: 2 2
       * resonance formation: 2 1
       * resonance decay: 1 2
       * etc.
       */
      std::fprintf(file_.get(), "%zu %zu %12.7f %12.7f %5i \n",
                   incoming_particles.size(), outgoing_particles.size(),
                   density, total_cross_section,
                   static_cast<int>(process_type));
    }
    for (const auto &p : incoming_particles) {
      write_particledata(p);
    }
    for (const auto &p : outgoing_particles) {
      write_particledata(p);
    }
  }
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_intermediate_time(
    const Particles &particles, const int event_number,
    const Clock & /*clock*/) {
  if (Contents & OscarTimesteps) {
    if (Format == OscarFormat2013) {
      std::fprintf(file_.get(), "# event %i out %zu\n", event_number + 1,
                   particles.size());
    } else {
      const size_t zero = 0;
      std::fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero,
                   event_number + 1);
    }
    write(particles);
  }
}

  /*!\Userguide
   * \page format_oscar_particlelist Oscar particles format
   * The format follows general block structure of OSCAR format:
   * \ref oscar_general_. There are two kinds of this format -
   * OSCAR2013 and OSCAR1999. Information about OSCAR standard can be found at
   * https://karman.physics.purdue.edu/OSCAR and
   * http://phy.duke.edu/~jeb65/oscar2013. SMASH OSCAR particles output
   * produces \c particle_lists.oscar file. Format is flexible, options that
   * regulate output can be found at \ref input_oscar_particlelist
   * and at \ref input_general_. **Particle output always gives
   * the current particle list at a specific time.**
   * Oscar1999
   * ---------
   * This is ASCII (text) human-readable output according to OSCAR 1999
   * standard. Format specifics are the following:\n
   * **Header**
   * \code
   * # OSC1999A
   * # final_id_p_x
   * # smash <version>
   * # Block format:
   * # nin nout event_number
   * # id pdg 0 px py pz p0 mass x y z t
   * # End of event: 0 0 event_number
   * #
   * \endcode
   *
   * **Output block header**
   * \code
   * nin nout /(not guaranteed) event_number/
   * \endcode
   *
   * For initial particles block (nin, nout) = (0, npart), for intermediate
   * and final - (nin, nout) = (npart, 0). Here npart - total number of
   * particles. Output block header is followed by npart particle lines.
   *
   * **Particle line**
   * \code
   * id pdg 0 px py pz p0 mass x y z t
   * \endcode
   *
   * \li \c id is an integer particle identifier.
   *     It is unique for every particle in event.
   * \li \c pdg is a PDG code of the particle (see http://pdg.lbl.gov/).
   * It contains all the quantum numbers of the particle and uniquely
   * identifies its type.
   * \li \c px \c py \c pz \c p0 - 3-momentum and energy
   * \li \c x \c y \c z \c t - coordinates and time
   *
   * **Event end line**
   * \code
   * 0 0 event_number
   * \endcode
   *
   * Oscar2013
   * ---------
   *
   * This is ASCII (text) human-readable output according to OSCAR 2013
   * standard. Format specifics are the following:\n
   * **Header**
   * \code
   * #!OSCAR2013 final_particle_list t x y z mass p0 px py pz pdg ID
   * # Units: fm fm fm fm GeV GeV GeV GeV GeV none none
   * \endcode
   *
   * **Output block header**\n
   * At start of event:
   * \code
   * # event ev_num in npart
   * \endcode
   * At end of event or intermediate particle list output:
   * \code
   * # event ev_num out npart
   * \endcode
   *
   * **Particle line**
   * \code
   * t x y z mass p0 px py pz pdg ID
   * \endcode
   *
   * **Event end line**
   * \code
   * # event ev_num end 0
   * \endcode
   *
   * \page format_oscar_collisions Oscar collisions format
   * The format follows general block structure of OSCAR format:
   * \ref oscar_general_. There are two kinds of this format -
   * OSCAR2013 and OSCAR1999. Information about OSCAR standard can be found at
   * https://karman.physics.purdue.edu/OSCAR and
   * http://phy.duke.edu/~jeb65/oscar2013. SMASH OSCAR collisions output
   * produces \c full_event_history.oscar file. Format is flexible, options
   * that regulate output can be found at \ref input_oscar_collisions
   * and at \ref input_general_. **Collision output always gives
   * a list of collisions/decays/box wall crossings plus optionally
   * initial and final configuration.**
   *
   * See also \ref collisions_output_in_box_modus_.
   *
   * Oscar1999
   * ---------
   * Format specifics are the following:\n
   * **Header**
   * \code
   * # OSC1999A
   * # full_event_history
   * # smash <version>
   * # Block format:
   * # nin nout event_number
   * # id pdg 0 px py pz p0 mass x y z t
   * # End of event: 0 0 event_number
   * #
   * \endcode
   *
   * **Output block header**
   * \code
   * nin nout /(not guaranteed) event_number/
   * \endcode
   * nin, nout are numbers of incoming and outgoing particles in a given
   * reaction in collision output file. If initial and final configurations
   * are written to collision file then (nin nout) = (0 npart) in the initial
   * configuration and (nin nout) = (npart 0) in the final.
   *
   * **Particle line**
   * \code
   * id pdg 0 px py pz p0 mass x y z t
   * \endcode
   * **Event end line**
   * \code
   * # event ev_num end 0
   * \endcode
   *
   * Oscar2013
   * ---------
   *  Format specifics are the following:\n
   * **Header**
   * \code
   * #!OSCAR2013 full_event_history t x y z mass p0 px py pz pdg ID
   * # Units: fm fm fm fm GeV GeV GeV GeV GeV none none
   * \endcode
   *
   * **Output block header**\n
   * At start of event:
   * \code
   * # event ev_num in npart
   * \endcode
   * At end of event:
   * \code
   * # event ev_num out npart
   * \endcode
   * At interaction:
   * \code
   * # interaction in nin out nout
   * \endcode
   *
   * **Particle line**
   * \code
   * t x y z mass p0 px py pz pdg ID
   * \endcode
   *
   * **Event end line**
   * \code
   * # event ev_num end 0
   * \endcode
   **/
template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::write_particledata(
    const ParticleData &data) {
  if (Format == OscarFormat2013) {
    std::fprintf(file_.get(), "%g %g %g %g %g %.9g %.9g %.9g %.9g %s %i\n",
        data.position().x0(),
        data.position().x1(), data.position().x2(), data.position().x3(),
        data.effective_mass(), data.momentum().x0(),
        data.momentum().x1(), data.momentum().x2(), data.momentum().x3(),
        data.pdgcode().string().c_str(), data.id());
  } else {
    std::fprintf(file_.get(), "%i %s %i %g %g %g %g %g %g %g %g %g\n",
        data.id(), data.pdgcode().string().c_str(), 0,
        data.momentum().x1(), data.momentum().x2(), data.momentum().x3(),
        data.momentum().x0(), data.effective_mass(),
        data.position().x1(), data.position().x2(), data.position().x3(),
        data.position().x0());
  }
}

namespace {
template <int Contents>
std::unique_ptr<OutputInterface> create_select_format(bf::path path,
                                                      Configuration config,
                                                      std::string name) {
  const bool modern_format =
      config.has_value({"2013_Format"}) ? config.take({"2013_Format"}) : false;
  if (modern_format) {
    return make_unique<OscarOutput<OscarFormat2013, Contents>>(std::move(path),
                                                               std::move(name));
  } else {
    return make_unique<OscarOutput<OscarFormat1999, Contents>>(std::move(path),
                                                               std::move(name));
  }
}
}  // unnamed namespace

std::unique_ptr<OutputInterface> create_oscar_output(bf::path path,
                                                     Configuration config) {
  if (config.has_value({"Oscar_Particlelist", "Enable"})) {
    auto subconfig = config["Oscar_Particlelist"];
    const bool enabled = subconfig.take({"Enable"});
    if (!enabled) {
      config.take({"Oscar_Particlelist"});
    } else {
      const bool only_final = subconfig.has_value({"Only_Final"})
                                  ? subconfig.take({"Only_Final"})
                                  : true;
      if (only_final) {
        return create_select_format<OscarParticlesAtEventend>(
            std::move(path), std::move(subconfig), "particle_lists");
      } else {
        return create_select_format<OscarTimesteps | OscarAtEventstart |
                                    OscarParticlesAtEventend>(
            std::move(path), std::move(subconfig), "particle_lists");
      }
    }
  }
  if (config.has_value({"Oscar_Collisions", "Enable"})) {
    auto subconfig = config["Oscar_Collisions"];
    const bool enabled = subconfig.take({"Enable"});
    if (!enabled) {
      config.take({"Oscar_Collisions"});
    } else {
      const bool print_start_end = subconfig.has_value({"Print_Start_End"})
                                  ? subconfig.take({"Print_Start_End"})
                                  : false;
      if (print_start_end) {
        return create_select_format<OscarInteractions | OscarAtEventstart |
                                    OscarParticlesAtEventend>(
            std::move(path), std::move(subconfig), "full_event_history");
      } else {
        return create_select_format<OscarInteractions>(
            std::move(path), std::move(subconfig), "full_event_history");
      }
    }
  }
  return {};  // return a nullptr to signify the end of OSCAR outputs in the
              // config file
}

  /*!\Userguide
   * \page input_dileptons Dileptons
   * Enables Dilepton Output together with DecayActionsFinderDilepton.
   * Dilepton Output saves information about decays, which include Dileptons,
   * at every timestep. The output is formatted in the
   * \ref format_oscar_collisions (OSCAR2013 format).
   *
   * The treatment of Dilepton Decays is special:
   *
   * \li Dileptons are treted via the time integration method, also called
   * shining method as described in \iref{Schmidt:2008hm}, chapter 2D.
   * This means that, because dilepton decays are so rare , possible decays are
   * written in the ouput every single timestep without ever performing them
   * and afterwards you weight them properly with a "shining weight" to
   * compensate for the over production.
   * \li The shining weight can be found in the weight element of the ouput.
   * \li The shining method is implemented in the DecayActionsFinderDilepton,
   * which is enabled together with the dilepton output.
   *
   * \note If you want dilepton decays, you also have to modify decaymodes.txt.
   * Dilepton decays are commented out by default.
   *
   * \key Enable (bool, optional, default = false):\n
   * true - Dilepton Output and DecayActionsFinderDilepton enabled\n
   * false - no Dilepton Output and no DecayActionsFinderDilepton
   **/

   /*!\Userguide
   * \page format_dilepton_output Dilepton Output
   * The format follows \ref format_oscar_collisions in the OSCAR2013
   * version. Dilepton Output produces the \c DileptonOutput.oscar file.
   * Shining weights are found in the weight element. For further
   * documentation and input options see: \ref input_dileptons.
   **/

std::unique_ptr<OutputInterface> create_dilepton_output(bf::path path) {
  /* for now the Oscar Output in the 2013 format is sufficient
   * for dilepton output
   */
  return make_unique<OscarOutput<OscarFormat2013, OscarInteractions>>(
                                            std::move(path), "DileptonOutput");
}

}  // namespace Smash
