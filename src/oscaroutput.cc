/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "smash/oscaroutput.h"

#include <string>

#include <boost/filesystem.hpp>

#include "smash/action.h"
#include "smash/clock.h"
#include "smash/config.h"
#include "smash/configuration.h"
#include "smash/cxx14compat.h"
#include "smash/forwarddeclarations.h"
#include "smash/particles.h"

namespace smash {

template <OscarOutputFormat Format, int Contents>
OscarOutput<Format, Contents>::OscarOutput(const bf::path &path,
                                           const std::string &name)
    : OutputInterface(name),
      file_{path /
                (name + ".oscar" + ((Format == OscarFormat1999) ? "1999" : "")),
            "w"} {
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
   * option in configuration file, see
   * \ref output_content_specific_options_ "content-specific output options".
   * Collisions output contains all collisions / decays / box wall crossings
   * and optionally initial and final configuration.
   */
  if (Format == OscarFormat2013) {
    std::fprintf(file_.get(),
                 "#!OSCAR2013 %s t x y z mass "
                 "p0 px py pz pdg ID charge\n",
                 name.c_str());
    std::fprintf(file_.get(),
                 "# Units: fm fm fm fm "
                 "GeV GeV GeV GeV GeV none none none\n");
    std::fprintf(file_.get(), "# %s\n", VERSION_MAJOR);
  } else if (Format == OscarFormat2013Extended) {
    std::fprintf(file_.get(),
                 "#!OSCAR2013Extended %s t x y z mass p0 px py pz"
                 " pdg ID charge ncoll form_time xsecfac proc_id_origin"
                 " proc_type_origin time_last_coll pdg_mother1 pdg_mother2\n",
                 name.c_str());
    std::fprintf(file_.get(),
                 "# Units: fm fm fm fm GeV GeV GeV GeV GeV"
                 " none none none none fm none none none fm none none\n");
    std::fprintf(file_.get(), "# %s\n", VERSION_MAJOR);
  } else {
    const std::string &oscar_name =
        name == "particle_lists" ? "final_id_p_x" : name;
    // This is necessary because OSCAR199A requires
    // this particular string for particle output.

    std::fprintf(file_.get(), "# OSC1999A\n# %s\n# %s\n", oscar_name.c_str(),
                 VERSION_MAJOR);
    std::fprintf(file_.get(), "# Block format:\n");
    std::fprintf(file_.get(), "# nin nout event_number\n");
    std::fprintf(file_.get(), "# id pdg 0 px py pz p0 mass x y z t\n");
    std::fprintf(file_.get(),
                 "# End of event: 0 0 event_number"
                 " impact_parameter\n");
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
  current_event_ = event_number;
  if (Contents & OscarAtEventstart) {
    if (Format == OscarFormat2013 || Format == OscarFormat2013Extended) {
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
                                                const int event_number,
                                                double impact_parameter) {
  if (Format == OscarFormat2013 || Format == OscarFormat2013Extended) {
    if (Contents & OscarParticlesAtEventend) {
      std::fprintf(file_.get(), "# event %i out %zu\n", event_number + 1,
                   particles.size());
      write(particles);
    }
    // Comment end of an event
    std::fprintf(file_.get(), "# event %i end 0 impact %7.3f\n",
                 event_number + 1, impact_parameter);
  } else {
    /* OSCAR line prefix : initial particles; final particles; event id
     * Last block of an event: initial = number of particles, final = 0
     * Block ends with null interaction. */
    const size_t zero = 0;
    if (Contents & OscarParticlesAtEventend) {
      std::fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero,
                   event_number + 1);
      write(particles);
    }
    // Null interaction marks the end of an event
    std::fprintf(file_.get(), "%zu %zu %i %7.3f\n", zero, zero,
                 event_number + 1, impact_parameter);
  }
  // Flush to disk
  std::fflush(file_.get());
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_interaction(const Action &action,
                                                   const double density) {
  if (Contents & OscarInteractions) {
    if (Format == OscarFormat2013 || Format == OscarFormat2013Extended) {
      std::fprintf(file_.get(),
                   "# interaction in %zu out %zu rho %12.7f weight %12.7g"
                   " partial %12.7f type %5i\n",
                   action.incoming_particles().size(),
                   action.outgoing_particles().size(), density,
                   action.get_total_weight(), action.get_partial_weight(),
                   static_cast<int>(action.get_type()));
    } else {
      /* OSCAR line prefix : initial final
       * particle creation: 0 1
       * particle 2<->2 collision: 2 2
       * resonance formation: 2 1
       * resonance decay: 1 2
       * etc.*/
      std::fprintf(file_.get(), "%zu %zu %12.7f %12.7f %12.7f %5i\n",
                   action.incoming_particles().size(),
                   action.outgoing_particles().size(), density,
                   action.get_total_weight(), action.get_partial_weight(),
                   static_cast<int>(action.get_type()));
    }
    for (const auto &p : action.incoming_particles()) {
      write_particledata(p);
    }
    for (const auto &p : action.outgoing_particles()) {
      write_particledata(p);
    }
  }
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_intermediate_time(
    const Particles &particles, const Clock &, const DensityParameters &) {
  if (Contents & OscarTimesteps) {
    if (Format == OscarFormat2013 || Format == OscarFormat2013Extended) {
      std::fprintf(file_.get(), "# event %i out %zu\n", current_event_ + 1,
                   particles.size());
    } else {
      const size_t zero = 0;
      std::fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero,
                   current_event_ + 1);
    }
    write(particles);
  }
}

/*!\Userguide
 * \page format_oscar_particlelist OSCAR particles format
 * The format follows general block structure of OSCAR format:
 * \ref oscar_general_. There are two kinds of this format -
 * OSCAR2013 and OSCAR1999. Information about OSCAR standard can be found at
 * https://karman.physics.purdue.edu/OSCAR and
 * http://phy.duke.edu/~jeb65/oscar2013. SMASH OSCAR particles output
 * produces \c particle_lists.oscar file. Format is flexible, options that
 * regulate output can be found at
 * \ref output_content_specific_options_ "content-specific output options".
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
 * # End of event: 0 0 event_number impact_parameter
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
 * 0 0 event_number impact_parameter
 * \endcode
 *
 * \anchor oscar2013_format
 * Oscar2013
 * ---------
 *
 * This is ASCII (text) human-readable output according to OSCAR 2013
 * standard. Format specifics are the following:\n
 * **Header**
 * \code
 * #!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge
 * # Units: fm fm fm fm GeV GeV GeV GeV GeV none none none
 * # SMASH_version
 * \endcode
 *
 * For the extended version of this output the header is modified to read:\n
 * **Header**
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">#!OSCAR2013 particle_lists
 *   t x y z mass p0 px py pz pdg
 *   ID charge ncoll form_time xsecfac proc_id_origin proc_type_origin
 *   t_last_coll pdg_mother1 pdg_mother2</span></div>
 * <div class="line"><span class="preprocessor">\# Units: fm fm fm fm
 *   GeV GeV GeV GeV GeV
 *   none none none fm none none none fm none none</span></div>
 * <div class="line"><span class="preprocessor">\# SMASH_version</span></div>
 * </div>
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
 * For the extended version the particle line contains
 *
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">t x y z
 *  mass p0 px py pz pdg ID charge Ncoll formation_time
 *  xsecfac process_ID_origin process_type_origin t_last_coll
 *  PDG_mother1 PDG_mother2</span></div>
 * </div>
 *
 * **Event end line**
 * \code
 * # event ev_num end 0 impact impact_parameter
 * \endcode
 *
 * \page format_oscar_collisions OSCAR collisions format
 * The format follows general block structure of OSCAR format:
 * \ref oscar_general_. There are two kinds of this format -
 * OSCAR2013 and OSCAR1999. Information about OSCAR standard can be found at
 * https://karman.physics.purdue.edu/OSCAR and
 * http://phy.duke.edu/~jeb65/oscar2013. SMASH OSCAR collisions output
 * produces \c full_event_history.oscar file. Format is flexible, options
 * that regulate output can be found at
 * \ref output_content_specific_options_ "content-specific output options".
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
 * #!OSCAR2013 full_event_history t x y z mass p0 px py pz pdg ID charge
 * # Units: fm fm fm fm GeV GeV GeV GeV GeV none none
 * # SMASH_version
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
 * t x y z mass p0 px py pz pdg ID charge
 * \endcode
 *
 * **Event end line**
 * \code
 * # event ev_num end 0 impact impact_parameter
 * \endcode
 **/

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::write_particledata(
    const ParticleData &data) {
  const FourVector pos = data.position();
  const FourVector mom = data.momentum();
  if (Format == OscarFormat2013) {
    std::fprintf(file_.get(), "%g %g %g %g %g %.9g %.9g %.9g %.9g %s %i %i\n",
                 pos.x0(), pos.x1(), pos.x2(), pos.x3(), data.effective_mass(),
                 mom.x0(), mom.x1(), mom.x2(), mom.x3(),
                 data.pdgcode().string().c_str(), data.id(),
                 data.type().charge());
  } else if (Format == OscarFormat2013Extended) {
    const auto h = data.get_history();
    std::fprintf(
        file_.get(),
        "%g %g %g %g %g %.9g %.9g %.9g"
        " %.9g %s %i %i %i %g %g %i %i %g %s %s\n",
        pos.x0(), pos.x1(), pos.x2(), pos.x3(), data.effective_mass(), mom.x0(),
        mom.x1(), mom.x2(), mom.x3(), data.pdgcode().string().c_str(),
        data.id(), data.type().charge(), h.collisions_per_particle,
        data.formation_time(), data.cross_section_scaling_factor(),
        h.id_process, static_cast<int>(h.process_type), h.time_last_collision,
        h.p1.string().c_str(), h.p2.string().c_str());
  } else {
    std::fprintf(file_.get(), "%i %s %i %g %g %g %g %g %g %g %g %g\n",
                 data.id(), data.pdgcode().string().c_str(), 0, mom.x1(),
                 mom.x2(), mom.x3(), mom.x0(), data.effective_mass(), pos.x1(),
                 pos.x2(), pos.x3(), pos.x0());
  }
}

namespace {
/**
 * Helper function that creates the oscar output with the format selected by
 * create_oscar_output (except for dileptons and photons).
 *
 * \tparam Contents Determines what infomration will be written to the output
 * \param[in] modern_format Use the 1999 or 2013 format
 * \param[in] path Path of output
 * \param[in] out_par Output parameters that hold the output configuration
 * \param[in] name (File)name of ouput
 * \return Unique pointer to oscar output
 */
template <int Contents>
std::unique_ptr<OutputInterface> create_select_format(
    bool modern_format, const bf::path &path, const OutputParameters &out_par,
    const std::string &name) {
  bool extended_format = (Contents & OscarInteractions) ? out_par.coll_extended
                                                        : out_par.part_extended;
  if (modern_format && extended_format) {
    return make_unique<OscarOutput<OscarFormat2013Extended, Contents>>(path,
                                                                       name);
  } else if (modern_format) {
    return make_unique<OscarOutput<OscarFormat2013, Contents>>(path, name);
  } else {
    return make_unique<OscarOutput<OscarFormat1999, Contents>>(path, name);
  }
}
}  // unnamed namespace

std::unique_ptr<OutputInterface> create_oscar_output(
    const std::string &format, const std::string &content, const bf::path &path,
    const OutputParameters &out_par) {
  if (format != "Oscar2013" && format != "Oscar1999") {
    throw std::invalid_argument("Creating Oscar output: unknown format");
  }
  const bool modern_format = (format == "Oscar2013");
  if (content == "Particles") {
    if (out_par.part_only_final) {
      return create_select_format<OscarParticlesAtEventend>(
          modern_format, path, out_par, "particle_lists");
    } else {
      return create_select_format<OscarTimesteps | OscarAtEventstart |
                                  OscarParticlesAtEventend>(
          modern_format, path, out_par, "particle_lists");
    }
  } else if (content == "Collisions") {
    if (out_par.coll_printstartend) {
      return create_select_format<OscarInteractions | OscarAtEventstart |
                                  OscarParticlesAtEventend>(
          modern_format, path, out_par, "full_event_history");
    } else {
      return create_select_format<OscarInteractions>(
          modern_format, path, out_par, "full_event_history");
    }
  } else if (content == "Dileptons") {
    if (out_par.dil_extended) {
      return make_unique<
          OscarOutput<OscarFormat2013Extended, OscarInteractions>>(path,
                                                                   "Dileptons");
    } else {
      return make_unique<OscarOutput<OscarFormat2013, OscarInteractions>>(
          path, "Dileptons");
    }
  } else if (content == "Photons") {
    if (modern_format) {
      return make_unique<OscarOutput<OscarFormat2013, OscarInteractions>>(
          path, "Photons");
    } else {
      return make_unique<OscarOutput<OscarFormat1999, OscarInteractions>>(
          path, "Photons");
    }
  }

  throw std::invalid_argument("Create_oscar_output got unknown content.");
}

}  // namespace smash
