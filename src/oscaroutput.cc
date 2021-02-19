/*
 *
 *    Copyright (c) 2014-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "smash/oscaroutput.h"

#include <boost/filesystem.hpp>
#include <string>

#include "smash/action.h"
#include "smash/clock.h"
#include "smash/config.h"
#include "smash/cxx14compat.h"
#include "smash/forwarddeclarations.h"

namespace smash {
static constexpr int LHyperSurfaceCrossing = LogArea::HyperSurfaceCrossing::id;

template <OscarOutputFormat Format, int Contents>
OscarOutput<Format, Contents>::OscarOutput(const bf::path &path,
                                           const std::string &name)
    : OutputInterface(name),
      file_{path /
                (name + ".oscar" + ((Format == OscarFormat1999) ? "1999" : "")),
            "w"} {
  /*!\Userguide
   * \page oscar_general_ OSCAR Block Structure
   * OSCAR outputs are a family of ASCII and binary formats that follow
   * the OSCAR format conventions. \n
   * **All OSCAR outputs have the same general structure: a header and an
   * arbitrary number of event blocks.** Each event block consists of an
   * arbitrary number of output blocks and special event end lines that mark the
   * end of an event. One output block consists of an output block header and N
   * particle lines, where N is specified in the output block header. \n
   * The file structure can be visualized in the following way:
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
   * Each OSCAR output can produce two types of files: collisions output (see
   * \ref format_oscar_collisions) and particles output (see \ref
   * format_oscar_particlelist). In both output types, the above structure is
   * the same, however the meaning of the blocks is different. **In the
   * collision file one output block typically corresponds to one collision
   * / decay / box wall crossing, while in the particles output one block
   * corresponds to the current particle list at one moment of time.**
   * The particles output may contain the particle list at event start
   * (immediately after initialization), at event end (which is reached when
   * time is larger or equal than \c End_Time in configuration file) and
   * periodically during the evolution, the output period is defined by
   * the \c Output_Interval option in the configuration file, see
   * \ref output_content_specific_options_ "content-specific output options".
   * The collisions output contains all collisions / decays / box wall crossings
   * and optionally the initial and final configuration.
   */
  if (Format == OscarFormat2013) {
    std::fprintf(file_.get(),
                 "#!OSCAR2013 %s t x y z mass "
                 "p0 px py pz pdg ID charge\n",
                 name.c_str());
    std::fprintf(file_.get(),
                 "# Units: fm fm fm fm "
                 "GeV GeV GeV GeV GeV none none e\n");
    std::fprintf(file_.get(), "# %s\n", VERSION_MAJOR);
  } else if (Format == OscarFormat2013Extended) {
    std::fprintf(file_.get(),
                 "#!OSCAR2013Extended %s t x y z mass p0 px py pz"
                 " pdg ID charge ncoll form_time xsecfac proc_id_origin"
                 " proc_type_origin time_last_coll pdg_mother1 pdg_mother2\n",
                 name.c_str());
    std::fprintf(file_.get(),
                 "# Units: fm fm fm fm GeV GeV GeV GeV GeV"
                 " none none e none fm none none none fm none none\n");
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
                                                  const int event_number,
                                                  const EventInfo &) {
  current_event_ = event_number;
  if (Contents & OscarAtEventstart) {
    if (Format == OscarFormat2013 || Format == OscarFormat2013Extended) {
      std::fprintf(file_.get(), "# event %i in %zu\n", event_number,
                   particles.size());
    } else {
      /* OSCAR line prefix : initial particles; final particles; event id
       * First block of an event: initial = 0, final = number of particles
       */
      const size_t zero = 0;
      std::fprintf(file_.get(), "%zu %zu %i\n", zero, particles.size(),
                   event_number);
    }
    if (!(Contents & OscarParticlesIC)) {
      // We do not want the inital particle list to be printed in case of IC
      // output
      write(particles);
    }
  }
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_eventend(const Particles &particles,
                                                const int event_number,
                                                const EventInfo &event) {
  if (Format == OscarFormat2013 || Format == OscarFormat2013Extended) {
    if (Contents & OscarParticlesAtEventend ||
        (Contents & OscarParticlesAtEventendIfNotEmpty && !event.empty_event)) {
      std::fprintf(file_.get(), "# event %i out %zu\n", event_number,
                   particles.size());
      write(particles);
    }
    // Comment end of an event
    const char *empty_event_str = event.empty_event ? "no" : "yes";
    std::fprintf(
        file_.get(),
        "# event %i end 0 impact %7.3f scattering_projectile_target %s\n",
        event_number, event.impact_parameter, empty_event_str);
  } else {
    /* OSCAR line prefix : initial particles; final particles; event id
     * Last block of an event: initial = number of particles, final = 0
     * Block ends with null interaction. */
    const size_t zero = 0;
    if (Contents & OscarParticlesAtEventend ||
        (Contents & OscarParticlesAtEventendIfNotEmpty && !event.empty_event)) {
      std::fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero,
                   event_number);
      write(particles);
    }
    // Null interaction marks the end of an event
    std::fprintf(file_.get(), "%zu %zu %i %7.3f\n", zero, zero, event_number,
                 event.impact_parameter);
  }
  // Flush to disk
  std::fflush(file_.get());

  if (Contents & OscarParticlesIC) {
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
  } else if (Contents & OscarParticlesIC) {
    for (const auto &p : action.incoming_particles()) {
      write_particledata(p);
    }
  }
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_intermediate_time(
    const Particles &particles, const std::unique_ptr<Clock> &,
    const DensityParameters &, const EventInfo &) {
  if (Contents & OscarTimesteps) {
    if (Format == OscarFormat2013 || Format == OscarFormat2013Extended) {
      std::fprintf(file_.get(), "# event %i out %zu\n", current_event_,
                   particles.size());
    } else {
      const size_t zero = 0;
      std::fprintf(file_.get(), "%zu %zu %i\n", particles.size(), zero,
                   current_event_);
    }
    write(particles);
  }
}

/*!\Userguide
 * \page format_oscar_particlelist OSCAR Particles Format
 * The OSCAR particles format follows the general block structure of the OSCAR
 * format: \ref oscar_general_. We distinguish between two versions -
 * OSCAR2013 and OSCAR1999. Information about OSCAR standard can be found at
 * https://karman.physics.purdue.edu/OSCAR and
 * http://phy.duke.edu/~jeb65/oscar2013. \n
 * Enabling the OSCAR output for particles in the config.yaml file
 * (see \ref input_output_options_), a so-called \c particle_lists.oscar file is
 * produced when executing SMASH. It allows for a certain degree of flexibility,
 * see \ref output_content_specific_options_ "Content-specific output options"
 * for further details. \n
 * **Unless IC output is enabled, the Particle output always provides the
 * current particle list at a specific time.** See \ref input_ic for details
 * about the particles IC output. \n
 *
 * \n
 * Oscar1999
 * ---------
 * Oscar1999 is an ASCII (text) human-readable output following the OSCAR 1999
 * standard. The format specifics are the following:
 *
 * \n
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
 * The header consists of 8 lines starting with '#', of which the last one
 * is basically empty.
 * They contain the following information:
 * -# The specific OSCAR1999 version the formatting follows - OSCAR1999A
 * -# The substructure of each particle line: (id - momentum - coordinates)
 * -# The SMASH version with which the oputput was generated
 * -# - 7. Info on the block structure
 *
 * \n
 * **Output block header** \n
 * Each output block starts with a line indicating the numbers of ingoing and
 * outgoing particles as well the number of the event.
 * \code
 * nin nout event_number
 * \endcode
 * With
 * \li \key nin: Number of ingoing particles
 * \li \key nout: Number of outgoing particles
 * \li \key event_number: Number of the event
 *
 * For initial timesteps, (nin, nout) = (0, Nparticles), while (nin, nout) =
 * (Nparticles, 0) for intermediate and final timesteps. Nparticles is the
 * total number of particles at the specific timestep. It may differ from one
 * timestep to another if the test case allows more interactions than only
 * elastic scatterings. The output block header is followed by Nparticles
 * particle lines.
 *
 * \n
 * **Particle line** \n
 * The particle lines are formatted as follows:
 * \code
 * id pdg 0 px py pz p0 mass x y z t
 * \endcode
 *
 * Where
 * \li \key id: Particle identifier in terms of an integer.
 *     It is unique for every particle in the event.
 * \li \key pdg: PDG code of the particle (see http://pdg.lbl.gov/).
 * It contains all quantum numbers and uniquely identifies its type.
 * \li \key px, \key py, \key pz, \key p0: 3-momentum and energy
 * \li \key mass: Particle's rest-mass
 * \li \key x, \key y, \key z, \key t: Space-time coordinates
 *
 * \n
 * **Event end line** \n
 * The end of an event is indicated by the following line:
 * \code
 * 0 0 event_number impact_parameter
 * \endcode
 *
 * With
 * \li \key event_number: Number of the event
 * \li \key impact_parameter: Impact parameter of the collisions. In case of
 * a box or sphere setup, this value is 0.0.
 *
 * \n
 * \anchor oscar2013_format
 * Oscar2013
 * ---------
 *
 * Oscar2013 is an ASCII (text) human-readable output following the OSCAR 2013
 * standard. The format specifics are the following:\n
 * \n
 * **Header**
 * \code
 * #!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge
 * # Units: fm fm fm fm GeV GeV GeV GeV GeV none none e
 * # SMASH_version
 * \endcode
 * The header consists of 3 lines starting with '#'. They contain the following
 * information:
 * -# Output version (OSCAR2013) and the type of output (particle_lists),
 * followed by the substructure of the particle lines.
 * -# Units of the quantities in the particle lines
 * -# SMASH version
 *
 * \n
 * **Extended Output: Header** \n
 * If desired, the OSCAR2013 output can be extended
 * by additional particle properties. This requires enabling the extended
 * output in the configuration file, see the \key Extended switch in
 * \ref output_content_specific_options_ "content-specific output options" for
 * further details. The header of the extended OSCAR output is structured
 * identically to the non-extended version, but simply contains more columns
 * because of the additional entries:
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">#!OSCAR2013 particle_lists
 *   t x y z mass p0 px py pz pdg
 *   ID charge ncoll form_time xsecfac proc_id_origin proc_type_origin
 *   t_last_coll pdg_mother1 pdg_mother2</span></div>
 * <div class="line"><span class="preprocessor">\# Units: fm fm fm fm GeV GeV
 * GeV GeV GeV none none e none fm none none none fm none none</span></div>
 * <div class="line"><span class="preprocessor">\# SMASH_version</span></div>
 * </div>
 *
 * \n
 * **Output block header**\n
 * Just as the OSCAR1999 format, the OSCAR2013 format is based on a block
 * structure. The beginning of a new block is marked by either the start of a
 * new event or a new intermediate output (at the next timestep). \n
 * \n
 * Output block header for a new event:
 * \code
 * # event ev_num in Nparticles
 * \endcode
 * Where
 * \li \key ev_num: Event number
 * \li \key Nparticles: Number of particles initialized at the beginning of
 * the event
 *
 * Note that 'event' and 'in' are no variables, but words that are printed in
 * the header. \n
 * \n
 * Output block header for an intermediate output:
 * \code
 * # event ev_num out Nparticles
 * \endcode
 * Where
 * \li \key ev_num: Event number
 * \li \key Nparticles: Number of particles at the end of the timestep
 *
 * Note that 'event' and 'out' are no variables, but words that are printed in
 * the header. \n
 * \n
 *
 * **Particle line**\n
 * The particle lines are formatted as follows:
 * \code
 * t x y z mass p0 px py pz pdg ID charge
 * \endcode
 * Apart from the order, the entries are identical to those of the OSCAR1999
 * output, the only additional one is:
 * \li \key charge: the electric charge of the particle in units of the
 * elementary charge e.
 *
 * \n
 * For the extended version the particle line contains
 *
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">t x y z mass p0 px py pz pdg
 * ID charge ncoll form_time xsecfac proc_id_origin proc_type_origin
 * time_last_coll pdg_mother1 pdg_mother2</span></div>
 * </div>
 *
 * \anchor extended_output_format_
 * The additional particle properties available in the extended output format
 * are:
 * \li \key ncoll: Number of collisions the particle has undergone
 * \li \key form_time: Formation time of the particle
 * \li \key xsecfac: Cross section scaling factor (if the particles are
 * not yet fully formed at the time of interaction, the cross section for the
 * underlying process is scaled down by the cross section scaling factor)
 * \li \key proc_id_origin: ID of the process of the particle's last interaction
 * \li \key proc_type_origin: Type of the last process the particle has
 * undergone. The possible process types are listed in \subpage process_type.
 * \li \key t_last_coll: time of the particle's last interaction (except wall
 * crossings)
 * \li \key pdg_mother1: PDG code of the 1st mother particle (0 in case the
 * particle is sampled in a thermal bubble. It is not updated by elastic
 * scatterings.)
 * \li \key pdg_mother2: PDG code of the 2nd mother particle (0 in case the
 * particle results from the decay of a resonance or the appearance of a
 * thermal bubble. In the former case, \key pdg_mother1 is the PDG code of this
 * resonance. It is not updated by elastic scatterings.)
 *
 * The mother particles are also set in case of an elastic scattering process.
 *
 * \page process_type Process Types
 * The available process types are:
 * \li \key 0: No previous process yet, particle was created at initialization
 * \li \key 1: Elastic scattering
 * \li \key 2: Resonance formation (2 -> 1)
 * \li \key 3: Inelastic binary scattering (2 -> 2)
 * \li \key 5: Resonance decay
 * \li \key 6: Box wall crossing (due to periodic boundary conditions)
 * \li \key 7: Forced thermalization
 * \li \key 8: Hypersurface crossing
 * \li \key 9: Bremsstrahlung
 * \li \key 41: Soft string excitation, single diffractive AB -> AX
 * \li \key 42: Soft string excitation, single diffractive AB -> XB
 * \li \key 43: Soft string excitation, double diffractive
 * \li \key 44: Soft string N-Nbar annihilation
 * \li \key 45: Soft sring excitation, non-diffractive
 * \li \key 46: Hard string excitation
 * \li \key 47: Failed string process
 *
 * \page format_oscar_particlelist
 * \n
 * **Event end line**\n
 * The end of an event is indicated by the following line:
 * \code
 * # event ev_num end 0 impact impact_parameter empty yes_or_no
 * \endcode
 * Where
 * \li \key ev_num: Event number
 * \li \key Nparticles: Number of particles at the end of the timestep
 * \li \key impact_parameter: Impact parameter of the collision in case of a
 * collider setup, 0.0 otherwise.
 * \li \key yes_or_no: "no" if there was an interaction between the projectile
 * and the target, "yes" otherwise. For non-collider setups, this is always
 * "no".
 *
 * Note that 'event', 'end', 'impact' and 'empty' are no variables, but words
 * that are printed in the header. \n
 * \n
 * \page format_oscar_collisions OSCAR Collisions Format
 * The OSCAR particles format follows the general block structure of the OSCAR
 * format: \ref oscar_general_. We distinguish between two versions -
 * OSCAR2013 and OSCAR1999. Information about OSCAR standard can be found at
 * https://karman.physics.purdue.edu/OSCAR and
 * http://phy.duke.edu/~jeb65/oscar2013. \n
 * Enabling the OSCAR output for collisions in the config.yaml file
 * (see \ref input_output_options_), a so-called \c full_event_history.oscar
 * file is produced when executing SMASH. It allows for a certain degree of
 * flexibility, see \ref output_content_specific_options_
 * "Content-specific output options" for further details. \n
 * **Collision output always gives
 * a list of collisions/decays/box wall crossings plus optionally
 * initial and final configuration.**
 *
 * See also \ref collisions_output_in_box_modus_. \n
 *
 * Oscar1999
 * ---------
 * Oscar1999 is an ASCII (text) human-readable output following the OSCAR 1999
 * standard. The format specifics are the following: \n
 * \n
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
 * The header consists of 8 lines starting with '#', of which the last one
 * is basically empty.
 * They contain the following information:
 * -# The specific OSCAR1999 version the formatting follows - OSC1999A
 * -# The filename
 * -# The SMASH version with which the oputput was generated
 * -# - 7. Info on the block structure
 *
 * \n
 * **Output block header**\n
 * Each output block starts with a line of the following format:
 * \code
 * nin nout density tot_weight part_weight proc_type
 * \endcode
 * With
 * \li \key nin: Number of ingoing particles (initial state particles)
 * \li \key nout: Number of outgoing particles (final state particles)
 * \li \key density: Density at the interaction point
 * \li \key tot_weight: Total weight of the interaction. This is the total cross
 * section in case of a scattering and the total decay width in case of a decay.
 * If there is no weight for the specific process, e.g. a wall crossing, it's
 * value is 0.0.
 * \li \key part_weight: The partial weight of the interaction. This is the
 * specific weight for the chosen final state.
 * \li \key proc_type: The type of the underlying process. See
 * \ref process_type for possible types.
 *
 * If the \key Print_Start_End option is set (see \ref
 * output_content_specific_options_ "content-specific output options" for
 * details), (nin, nout) = (0, Nparticles) in the
 * initial timestep and (nin, nout) = (Nparticles, 0) in the final timestep.
 *
 * \n
 * **Particle line**\n
 * The particle lines are formatted as follows:
 * \code
 * id pdg 0 px py pz p0 mass x y z t
 * \endcode
 *
 * Where
 * \li \key id: Particle identifier in terms of an integer.
 *     It is unique for every particle in the event.
 * \li \key pdg: PDG code of the particle (see http://pdg.lbl.gov/).
 * It contains all quantum numbers and uniquely identifies its type.
 * \li \key px, \key py, \key pz, \key p0: 3-momentum and energy
 * \li \key mass: Particle's rest-mass
 * \li \key x, \key y, \key z, \key t: Space-time coordinates
 *
 * \n
 * **Event end line** \n
 * The end of an event is indicated by the following line:
 * \code
 * 0 0 event_number impact_parameter
 * \endcode
 *
 * With
 * \li \key event_number: Number of the event
 * \li \key impact_parameter: Impact parameter of the collisions. In case of
 * a box or sphere setup, this value is 0.0.
 *
 * \n
 * Oscar2013
 * ---------
 *  Oscar2013 is an ASCII (text) human-readable output following the OSCAR 2013
 * standard. The format specifics are the following:\n
 * \n
 * **Header**
 * \code
 * #!OSCAR2013 full_event_history t x y z mass p0 px py pz pdg ID charge
 * # Units: fm fm fm fm GeV GeV GeV GeV GeV none none
 * # SMASH_version
 * \endcode
 * The header consists of 3 lines starting with '#'. They contain the following
 * information:
 * -# Output version (OSCAR2013) and the type of output (particle_lists),
 * followed by the substructure of the particle lines.
 * -# Units of the quantities in the particle lines
 * -# SMASH version
 *
 * \n
 *
 * **Extended Output: Header** \n
 * If desired, the OSCAR2013 output can be extended
 * by additional particle properties. This requires enabling the extended
 * output in the configuration file, see the \key Extended switch in
 * \ref output_content_specific_options_ "content-specific output options" for
 * further details. The header of the extended OSCAR output is structured
 * identically to the non-extended version, but simply contains more columns
 * because of the additional entries:
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">#!OSCAR2013 particle_lists
 *   t x y z mass p0 px py pz pdg
 *   ID charge ncoll form_time xsecfac proc_id_origin proc_type_origin
 *   t_last_coll pdg_mother1 pdg_mother2</span></div>
 * <div class="line"><span class="preprocessor">\# Units: fm fm fm fm GeV GeV
 * GeV GeV GeV none none e none fm none none none fm none none</span></div>
 * <div class="line"><span class="preprocessor">\# SMASH_version</span></div>
 * </div>
 *
 * \n
 * **Output block header**\n
 * Just as the OSCAR1999 format, the OSCAR2013 format is based on a block
 * structure, where each block corresponds to one interaction. Each block starts
 * with a line formatted as follows:
 * <div class="fragment">
 * <div class="line"> <span class="preprocessor">
 *  \# interaction in nin out nout rho density weight tot_weight partial
 *part_weight type proc_type </span></div>
 * </div>
 * where
 * \li \key nin: Number of ingoing
 *particles (initial state particles) \li \key nout: Number of outgoing
 *particles (final state particles) \li \key density: Density at the interaction
 *point \li \key tot_weight: Total weight of the interaction. This is the total
 *cross section in case of a scattering and the total decay width in case of a
 *decay. If there is no weight for the specific process, e.g. a wall crossing,
 *it's value is 0.0. \li \key part_weight: The partial weight of the
 *interaction. This is the specific weight for the chosen final state. \li \key
 *proc_type: The type of the underlying process. See \ref process_type for
 *possible types.
 *
 * Note, that "interaction", "in", "out", "rho", "weight", "partial" and "type"
 * are no variables, but words that are printed.\n
 *
 * \n
 * **Particle line** \n
 * The particle lines are formatted as follows:
 * \code
 * t x y z mass p0 px py pz pdg ID charge
 * \endcode
 * Apart from the order, the entries are identical to those of the OSCAR1999
 * output, the only additional one is:
 * \li \key charge: the electric charge of the particle in units of the
 * elementary charge e.
 *
 * \n
 * For the extended version the particle line contains
 *
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">t x y z
 *  mass p0 px py pz pdg ID charge Ncoll formation_time
 *  xsecfac process_ID_origin process_type_origin t_last_coll
 *  PDG_mother1 PDG_mother2</span></div>
 * </div>
 *
 * The additional particle properties available in the extended output format
 * are:
 * \li \key ncoll: Number of collisions the particle has undergone
 * \li \key form_time: Formation time of the particle
 * \li \key xsecfac: Cross section scaling factor (if the particles are
 * not yet fully formed at the time of interaction, the cross section for the
 * underlying process is scaled down by the cross section scaling factor)
 * \li \key proc_id_origin: ID of the process of the particle's last interaction
 * \li \key proc_type_origin: Type of the last process the particle has
 * undergone. The possible process types are listed in \ref process_type.
 * \li \key t_last_coll: time of the particle's last interaction (except wall
 * crossings)
 * \li \key pdg_mother1: PDG code of the 1st mother particle
 * \li \key pdg_mother2: PDG code of the 2nd mother particle (0 in case the
 * particle results from the decay of a resonance, then \key pdg_mother1 is
 * the PDG code of this resonance)
 *
 * The mother particles are also set in case of an elastic scattering process.
 * \n
 * \n
 *
 * **Event end line** \n
 * The end of an event is indicated by the following line:
 * \code
 * # event ev_num end 0 impact impact_parameter
 * \endcode
 * where
 * \li \key ev_num: The event's number
 * \li \key impact_parameter: impact parameter of the collision in case of a
 * collider setup, otherwise 0.0.
 *
 * Note, that "event", "end" and "impact" are no variables, but words
 * that are printed.
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
        data.formation_time(), data.xsec_scaling_factor(), h.id_process,
        static_cast<int>(h.process_type), h.time_last_collision,
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
  } else if (modern_format && !extended_format) {
    return make_unique<OscarOutput<OscarFormat2013, Contents>>(path, name);
  } else if (!modern_format && !extended_format) {
    return make_unique<OscarOutput<OscarFormat1999, Contents>>(path, name);
  } else {
    // Only remaining possibility: (!modern_format && extended_format)
    logg[LOutput].warn() << "Creating Oscar output: "
                         << "There is no extended Oscar1999 format.";
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
    if (out_par.part_only_final == OutputOnlyFinal::Yes) {
      return create_select_format<OscarParticlesAtEventend>(
          modern_format, path, out_par, "particle_lists");
    } else if (out_par.part_only_final == OutputOnlyFinal::IfNotEmpty) {
      return create_select_format<OscarParticlesAtEventendIfNotEmpty>(
          modern_format, path, out_par, "particle_lists");

    } else {  // out_par.part_only_final == OutputOnlyFinal::No
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
    if (modern_format && out_par.dil_extended) {
      return make_unique<
          OscarOutput<OscarFormat2013Extended, OscarInteractions>>(path,
                                                                   "Dileptons");
    } else if (modern_format && !out_par.dil_extended) {
      return make_unique<OscarOutput<OscarFormat2013, OscarInteractions>>(
          path, "Dileptons");
    } else if (!modern_format && !out_par.dil_extended) {
      return make_unique<OscarOutput<OscarFormat1999, OscarInteractions>>(
          path, "Dileptons");
    } else if (!modern_format && out_par.dil_extended) {
      logg[LOutput].warn()
          << "Creating Oscar output: "
          << "There is no extended Oscar1999 (dileptons) format.";
    }
  } else if (content == "Photons") {
    if (modern_format && !out_par.photons_extended) {
      return make_unique<OscarOutput<OscarFormat2013, OscarInteractions>>(
          path, "Photons");
    } else if (modern_format && out_par.photons_extended) {
      return make_unique<
          OscarOutput<OscarFormat2013Extended, OscarInteractions>>(path,
                                                                   "Photons");
    } else if (!modern_format && !out_par.photons_extended) {
      return make_unique<OscarOutput<OscarFormat1999, OscarInteractions>>(
          path, "Photons");
    } else if (!modern_format && out_par.photons_extended) {
      logg[LOutput].warn()
          << "Creating Oscar output: "
          << "There is no extended Oscar1999 (photons) format.";
    }
  } else if (content == "Initial_Conditions") {
    if (modern_format && !out_par.ic_extended) {
      return make_unique<
          OscarOutput<OscarFormat2013, OscarParticlesIC | OscarAtEventstart>>(
          path, "SMASH_IC");
    } else if (modern_format && out_par.ic_extended) {
      return make_unique<OscarOutput<OscarFormat2013Extended,
                                     OscarParticlesIC | OscarAtEventstart>>(
          path, "SMASH_IC");
    } else if (!modern_format && !out_par.ic_extended) {
      return make_unique<
          OscarOutput<OscarFormat1999, OscarParticlesIC | OscarAtEventstart>>(
          path, "SMASH_IC");
    } else if (!modern_format && out_par.ic_extended) {
      logg[LOutput].warn()
          << "Creating Oscar output: "
          << "There is no extended Oscar1999 (initial conditions) format.";
    }
  }

  throw std::invalid_argument("Create_oscar_output got unknown content.");
}

}  // namespace smash
