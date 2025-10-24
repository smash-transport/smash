/*
 *
 *    Copyright (c) 2014-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/oscaroutput.h"

#include <filesystem>
#include <string>

#include "smash/action.h"
#include "smash/clock.h"
#include "smash/config.h"
#include "smash/forwarddeclarations.h"

namespace smash {

template <OscarOutputFormat Format, int Contents>
OscarOutput<Format, Contents>::OscarOutput(
    const std::filesystem::path &path, const std::string &name,
    const std::vector<std::string> quantities)
    : OutputInterface(name),
      file_{path / (name + ((Format == ASCII) ? ".dat" : ".oscar") +
                    ((Format == OscarFormat1999) ? "1999" : "")),
            "w"},
      formatter_{Format == ASCII ? quantities
                 : (Format == OscarFormat2013)
                     ? OutputDefaultQuantities::oscar2013
                 : (Format == OscarFormat2013Extended)
                     ? OutputDefaultQuantities::oscar2013extended
                     : OutputDefaultQuantities::oscar1999} {
  /*!\Userguide
   * \page doxypage_output_oscar
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
   * \ref doxypage_output_oscar_collisions) and particles output (see \ref
   * doxypage_output_oscar_particles). In both output types, the above
   * structure is the same, however the meaning of the blocks is different. **In
   * the collision file one output block typically corresponds to one collision
   * / decay / box wall crossing, while in the particles output one block
   * corresponds to the current particle list at one moment of time.**
   * The particles output may contain the particle list at event start
   * (immediately after initialization), at event end (which is reached when
   * time is larger or equal than \c End_Time in configuration file) and
   * periodically during the evolution, the output period is defined by
   * the \c Output_Interval option in the configuration file, see
   * \ref input_output_content_specific_ "content-specific output options".
   * The collisions output contains all collisions / decays / box wall crossings
   * and optionally the initial and final configuration.
   */
  if (Format != ASCII && !quantities.empty()) {
    throw std::logic_error(
        "Non-empty Quantities given alongside format other than ASCII.");
  }
  std::string format_name;
  if (Format == ASCII) {
    format_name = "ASCII";
  } else if (Format == OscarFormat2013) {
    format_name = "OSCAR2013";
  } else if (Format == OscarFormat2013Extended) {
    format_name = "OSCAR2013Extended";
  } else {
    format_name = "OSC1999A";
  }
  if (Format == ASCII || Format == OscarFormat2013 ||
      Format == OscarFormat2013Extended) {
    std::fprintf(file_.get(), "#!%s %s %s\n", format_name.c_str(), name.c_str(),
                 formatter_.quantities_line().c_str());
    std::fprintf(file_.get(), "# Units: %s\n", formatter_.unit_line().c_str());
    std::fprintf(file_.get(), "# %s\n", SMASH_VERSION);
  } else {
    const std::string &oscar_name =
        name == "particle_lists" ? "final_id_p_x" : name;
    // This is necessary because OSCAR1999A requires
    // this particular string for particle output.

    std::fprintf(file_.get(), "# %s\n# %s\n# %s\n", format_name.c_str(),
                 oscar_name.c_str(), SMASH_VERSION);
    std::fprintf(file_.get(), "# Block format:\n");
    if (oscar_name == "full_event_history") {
      std::fprintf(file_.get(),
                   "# nin nout density tot_weight part_weight proc_type\n");
    } else {
      std::fprintf(file_.get(), "# nin nout event_number ensemble_number\n");
    }
    std::fprintf(file_.get(), "# %s\n", formatter_.quantities_line().c_str());
    std::fprintf(
        file_.get(),
        "# End of event: 0 0 event_number ensemble_number impact_parameter\n");
    std::fprintf(file_.get(), "#\n");
  }
}

template <OscarOutputFormat Format, int Contents>
inline void OscarOutput<Format, Contents>::write(const Particles &particles) {
  write_in_chunk<ToASCII>(particles);
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_eventstart(const Particles &particles,
                                                  const EventLabel &event_label,
                                                  const EventInfo &) {
  // We do not want the inital particle list or number to be printed in case of
  // IC output
  if (Contents & OscarAtEventstart && !(Contents & OscarParticlesIC)) {
    if (Format == ASCII || Format == OscarFormat2013 ||
        Format == OscarFormat2013Extended) {
      std::fprintf(file_.get(), "# event %i ensemble %i in %zu\n",
                   event_label.event_number, event_label.ensemble_number,
                   particles.size());
    } else {
      /* OSCAR line prefix : initial particles; final particles; event id
       * First block of an event: initial = 0, final = number of particles
       */
      const size_t zero = 0;
      std::fprintf(file_.get(), "%zu %zu %i %i\n", zero, particles.size(),
                   event_label.event_number, event_label.ensemble_number);
    }
    write(particles);
  } else if (Contents & OscarParticlesIC) {
    if (Format == ASCII || Format == OscarFormat2013 ||
        Format == OscarFormat2013Extended) {
      std::fprintf(file_.get(), "# event %i ensemble %i start\n",
                   event_label.event_number, event_label.ensemble_number);
    } else if (Format == OscarFormat1999) {
      const size_t zero = 0;
      std::fprintf(file_.get(), "%zu %zu %i %i\n", zero, zero,
                   event_label.event_number, event_label.ensemble_number);
    }
  }
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_eventend(const Particles &particles,
                                                const EventLabel &event_label,
                                                const EventInfo &event) {
  if (Format == ASCII || Format == OscarFormat2013 ||
      Format == OscarFormat2013Extended) {
    if (Contents & OscarParticlesAtEventend ||
        (Contents & OscarParticlesAtEventendIfNotEmpty && !event.empty_event)) {
      std::fprintf(file_.get(), "# event %i ensemble %i out %zu\n",
                   event_label.event_number, event_label.ensemble_number,
                   particles.size());
      write(particles);
    }
    // Comment end of an event
    if (!(Contents & OscarParticlesIC)) {
      const char *empty_event_str = event.empty_event ? "no" : "yes";
      std::fprintf(file_.get(),
                   "# event %i ensemble %i end 0 impact %7.3f "
                   "scattering_projectile_target %s\n",
                   event_label.event_number, event_label.ensemble_number,
                   event.impact_parameter, empty_event_str);
    } else {
      std::fprintf(file_.get(), "# event %i ensemble %i end\n",
                   event_label.event_number, event_label.ensemble_number);
    }
  } else {
    /* OSCAR line prefix : initial particles; final particles; event id
     * Last block of an event: initial = number of particles, final = 0
     * Block ends with null interaction. */
    const size_t zero = 0;
    if (Contents & OscarParticlesAtEventend ||
        (Contents & OscarParticlesAtEventendIfNotEmpty && !event.empty_event)) {
      std::fprintf(file_.get(), "%zu %zu %i %i\n", particles.size(), zero,
                   event_label.event_number, event_label.ensemble_number);
      write(particles);
    }
    // Null interaction marks the end of an event
    std::fprintf(file_.get(), "%zu %zu %i %i %7.3f\n", zero, zero,
                 event_label.event_number, event_label.ensemble_number,
                 event.impact_parameter);
  }
  // Flush to disk
  std::fflush(file_.get());
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::at_interaction(const Action &action,
                                                   const double density) {
  if (Contents & OscarInteractions) {
    if (Format == ASCII || Format == OscarFormat2013 ||
        Format == OscarFormat2013Extended) {
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
    const DensityParameters &, const EventLabel &event_label,
    const EventInfo &) {
  if (Contents & OscarTimesteps) {
    if (Format == ASCII || Format == OscarFormat2013 ||
        Format == OscarFormat2013Extended) {
      std::fprintf(file_.get(), "# event %i ensemble %i out %zu\n",
                   event_label.event_number, event_label.ensemble_number,
                   particles.size());
    } else {
      const size_t zero = 0;
      std::fprintf(file_.get(), "%zu %zu %i %i\n", particles.size(), zero,
                   event_label.event_number, event_label.ensemble_number);
    }
    write(particles);
  }
}

/*!\Userguide
 * \page doxypage_output_oscar_particles
 * The OSCAR particles format follows the general block structure of the
 * \ref doxypage_output_oscar. We distinguish between two versions, OSCAR2013
 * and OSCAR1999. Additional information about OSCAR standard can be found
 * <a href="http://phy.duke.edu/~jeb65/oscar2013">here</a>. \n Enabling
 * the OSCAR output for particles in the config.yaml file (see \ref
 * doxypage_input_conf_output), a so-called \c particle_lists.oscar file is
 * produced when executing SMASH. It allows for a certain degree of flexibility,
 * see \ref input_output_content_specific_ "Content-specific output options" for
 * further details.
 *
 * **Unless IC output is enabled, the Particle output always provides the
 * current particle list at a specific time.** See \ref
 * doxypage_output_initial_conditions for details about the particles IC output.
 * Even though they are compatible, we do not recommend using Oscar1999 for the
 * Initial_Conditions output.
 * \n
 *
 * \anchor oscar2013_format
 * <h2> Oscar2013 </h2>
 *
 * Oscar2013 is an ASCII (text) human-readable output following the OSCAR 2013
 * standard. The format specifics are the following:\n
 * \n
 * **File header**
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
 * **File header for extended output**
 *
 * If desired, the OSCAR2013 output can be extended
 * by additional particle properties. This requires enabling the extended
 * output in the configuration file, see the \key Extended switch in
 * \ref input_output_content_specific_ "content-specific output options" for
 * further details. The header of the extended OSCAR output is structured
 * identically to the non-extended version, but simply contains more columns
 * because of the additional entries:
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">#!OSCAR2013 particle_lists
 *   t x y z mass p0 px py pz pdg
 *   ID charge ncoll form_time xsecfac proc_id_origin proc_type_origin
 *   time_last_coll pdg_mother1 pdg_mother2 baryon_number strangeness</span>
 * </div>
 * <div class="line"><span class="preprocessor">\# Units: fm fm fm fm GeV GeV
 * GeV GeV GeV none none e none fm none none none fm none none none none</span>
 * </div> <div class="line"><span class="preprocessor">\# SMASH_version</span>
 * </div>
 * </div>
 *
 * **Block header**
 *
 * The OSCAR2013 format is based on a block structure. The beginning of a new
 * block is marked by either the start of a new event or a new intermediate
 * output (at the next timestep).
 *
 * Output block header for a new event:
 * \code
 * # event ev_num ensemble ens_num in Nparticles
 * \endcode
 * Where
 * \li \key ev_num: Event number
 * \li \key ens_num: Ensemble number
 * \li \key Nparticles: Number of particles initialized at the beginning of
 * the event
 *
 * Note that `event`, `ensemble` and `in` are no variables, but words that are
 * printed in the header.
 *
 * Output block header for an intermediate output:
 * \code
 * # event ev_num ensemble ens_num out Nparticles
 * \endcode
 * Where
 * \li \key ev_num: Event number
 * \li \key ens_num: Ensemble number
 * \li \key Nparticles: Number of particles at the end of the timestep
 *
 * Note that `event`, `ensemble` and `out` are no variables, but words that are
 * printed in the header.
 *
 * **Particle line**
 *
 * The particle lines are formatted as follows:
 * \code
 * t x y z mass p0 px py pz pdg ID charge
 * \endcode
 *
 * where
 * \li \key t, \key x, \key y, \key z: Space-time coordinates
 * \li \key mass: Rest-mass
 * \li \key p0, \key px, \key py, \key pz: Energy and 3-momentum
 * \li \key pdg: PDG code of the particle (see http://pdg.lbl.gov/).
 * It contains all quantum numbers and uniquely identifies its type
 * \li \key ID: Particle identifier in terms of an integer, it is
 * unique for each particle in the event
 * \li \key charge: the electric charge of the particle in units of the
 * elementary charge e
 *
 * For the **extended output** the particle line contains
 *
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">t x y z mass p0 px py pz pdg
 * ID charge ncoll form_time xsecfac proc_id_origin proc_type_origin
 * time_last_coll pdg_mother1 pdg_mother2 baryon_number strangeness</span></div>
 * </div>
 *
 * \anchor extended_output_format_
 * The additional particle properties available in the extended output format
 * are:
 * \li \key ncoll: Number of collisions the particle has undergone
 * \li \key form_time: Formation time of the particle
 * \li \key xsecfac: Cross section scaling factor (if the particles are
 *     not yet fully formed at the time of interaction, the cross section for
 *     the underlying process is scaled down by the cross section scaling
 *     factor)
 * \li \key proc_id_origin: ID of the process of the particle's last interaction
 * \li \key proc_type_origin: Type of the last process the particle has
 *     undergone. The possible process types are listed in
 *     \ref doxypage_output_process_types
 * \li \key time_last_coll: time of the particle's last interaction (except wall
 *     crossings)
 * \li \key pdg_mother1: PDG code of the 1st mother particle (0 in case the
 *     particle is sampled in a thermal bubble. It is not updated by elastic
 *     scatterings)
 * \li \key pdg_mother2: PDG code of the 2nd mother particle (0 in case the
 *     particle results from the decay of a resonance or the appearance of a
 *     thermal bubble. In the former case, \key pdg_mother1 is the PDG code of
 *     this resonance. It is not updated by elastic scatterings)
 * \li \key baryon_number: Baryon number of the particle. 1 for baryons, -1 for
 *     anti-baryons and 0 for mesons
 * \li \key strangeness: Strangeness of the particle
 *
 * The mother particles are also set in case of an elastic scattering process.
 *
 * **Event end line**\n
 * The end of an event is indicated by the following line:
 * \code
 * # event ev_num ensemble ens_num end 0 impact impact_parameter empty yes_or_no
 * \endcode
 * Where
 * \li \key ev_num: Event number
 * \li \key ens_num: Ensemble number
 * \li \key impact_parameter: Impact parameter of the collision in case of a
 *          collider setup, 0.0 otherwise
 * \li \key yes_or_no: "no" if there was an interaction between the projectile
 * and the target, "yes" otherwise. For non-collider setups, this is always
 * "no"
 *
 * Note that `event`, `end`, `impact` and `empty` are no variables, but words
 * that are printed in the header.
 *
 * <h2> Oscar1999 </h2>
 *
 * Oscar1999 is an ASCII (text) human-readable output following the OSCAR 1999
 * standard. The format specifics are the following:
 *
 * **File header**
 * \code
 * # OSC1999A
 * # final_id_p_x
 * # smash <version>
 * # Block format:
 * # nin nout event_number ensemble_number
 * # id pdg 0 px py pz p0 mass x y z t
 * # End of event: 0 0 event_number ensemble_number impact_parameter
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
 * **Block header**
 *
 * Each output block starts with a line indicating the numbers of ingoing and
 * outgoing particles as well the numbers of the event and ensemble.
 * \code
 * nin nout event_number ensemble_number
 * \endcode
 * With
 * \li \key nin: Number of ingoing particles
 * \li \key nout: Number of outgoing particles
 * \li \key event_number: Number of the event
 * \li \key ensemble_number: Number of the ensemble
 *
 * For initial timesteps, (nin, nout) = (0, Nparticles), while (nin, nout) =
 * (Nparticles, 0) for intermediate and final timesteps. Nparticles is the
 * total number of particles at the specific timestep. It may differ from one
 * timestep to another if the test case allows more interactions than only
 * elastic scatterings. The output block header is followed by Nparticles
 * particle lines. In the Initial_Conditions output, Nparticles is always 0.
 *
 * **Particle line**
 *
 * The particle lines are formatted as follows:
 * \code
 * id pdg 0 px py pz p0 mass x y z t
 * \endcode
 *
 * Apart from the order, the entries are identical to those of the OSCAR2013
 * output.
 *
 * **Event end line**
 *
 * The end of an event is indicated by the following line:
 * \code
 * 0 0 event_number ensemble_number impact_parameter
 * \endcode
 *
 * With
 * \li \key event_number: Number of the event
 * \li \key ensemble_number: Number of the ensemble
 * \li \key impact_parameter: Impact parameter of the collisions. In case of
 * a box or sphere setup, this value is 0.0
 *
 * \page doxypage_output_process_types
 * The available process types are summarized in the following table.
 *
 * <table>
 * <tr><th>Process number<th>Description
 * \process_type{0} No previous process yet, particle was created at
 *                  initialization
 * \process_type{1} Elastic scattering
 * \process_type{2} Resonance formation (2 &rarr; 1)
 * \process_type{3} Inelastic binary scattering (2 &rarr; 2)
 * \process_type{4} Inelastic multi-particle scattering (2 &rarr; 3)
 * \process_type{5} Resonance decay
 * \process_type{6} Box wall crossing (due to periodic boundary conditions)
 * \process_type{7} Forced thermalization, many particles are replaced by a
 *                  thermalized ensemble
 * \process_type{8} Fluidization, particles that obey the fluidization
 *                  condition given in the configuration file are removed
 *                  from the evolution and printed to a separate output, to
 *                  serve as initial conditions for hybrid models.
 * \process_type{21} Fluidization as above, but particles are not removed from
 *                   the evolution. They are instead tagged as core.
 * \process_type{9}  Bremsstrahlung process: a + b &rarr; a + b + photon
 * \process_type{10} Inelastic multi-particle meson scattering (3 &rarr; 1)
 * \process_type{11} Inelastic multi-particle scattering (3 &rarr; 2)
 * \process_type{12} Inelastic multi-particle scattering (5 &rarr; 2)
 * \process_type{13} Inelastic multi-particle scattering (2 &rarr; 5)
 * \process_type{14} Inelastic multi-particle scattering (4 &rarr; 2)
 * \process_type{15} Inelastic multi-particle scattering (2 &rarr; 4)
 * \process_type{41} Soft string excitation, single diffractive AB &rarr; AX.
 *                   Both quark and anti-/di-quark taken from B.
 * \process_type{42} Soft string excitation, single diffractive AB &rarr; XB.
 *                   Both quark and anti-/di-quark taken from A. It makes sense
 *                   to distinguish it from AB &rarr; AX, because A and B can
 *                   be particles of different types, for example, a pion and a
 *                   proton. It matters then whether the pion or the proton
 *                   creates the string.
 * \process_type{43} Soft string excitation, double diffractive. Two strings are
 *                   formed, one from A and one from B.
 * \process_type{44} Soft string N-Nbar annihilation, a special case of
 *                   baryon-antibaryon annihilation. One pair qqbar annihilates
 *                   immediately and then two strings are formed.
 * \process_type{45} Soft string excitation, non-diffractive. Two strings are
 *                   formed both have ends in A and B.
 * \process_type{46} Hard string excitation, hard string process involving 2
 *                   &rarr; 2 QCD process by PYTHIA. Here quarks do not simply
 *                   form a string. They actually scatter on parton level first.
 * \process_type{47} Failed string process, Soft String NNbar annihilation
 *                   process can fail by lack of energy. This is a tag we add to
 *                   avoid mislabeling the events.
 * \process_type{90} Add or remove particle(s) process, which ignores
 *                   conservation laws. It can be thought of as a 0 &rarr; 1 or
 *                   a 1 &rarr; 0 process.
 * </table>
 *
 * \page doxypage_output_oscar_collisions
 * The OSCAR particles format follows the general block structure of the
 * \ref doxypage_output_oscar. We distinguish between two versions, OSCAR2013
 * and OSCAR1999. Additional information about OSCAR standard can be found
 * <a href="http://phy.duke.edu/~jeb65/oscar2013">here</a>. \n
 * Enabling the OSCAR output for collisions in the config.yaml file
 * (see \ref doxypage_input_conf_output), a so-called \c
 * full_event_history.oscar file is produced when executing SMASH. It allows for
 * a certain degree of flexibility, see \ref input_output_content_specific_
 * "Content-specific output options" for further details. \n
 * **Collision output always gives
 * a list of collisions/decays/box wall crossings plus optionally
 * initial and final configuration.**
 *
 * See also \ref doxypage_output_collisions_box_modus. \n
 *
 * \note The particle and event end lines for both OSCAR 2013 and 1999 formats
 * are identical as in the \ref doxypage_output_oscar_particles.
 *
 * <h2> Oscar2013 </h2>
 *
 *  Oscar2013 is an ASCII (text) human-readable output following the OSCAR 2013
 * standard. The format specifics are the following:\n
 * \n
 * **File header**
 * \code
 * #!OSCAR2013 full_event_history t x y z mass p0 px py pz pdg ID charge
 * # Units: fm fm fm fm GeV GeV GeV GeV GeV none none
 * # SMASH_version
 * \endcode
 * The header consists of 3 lines starting with '#'. They contain the following
 * information:
 * -# Output version (OSCAR2013) and the type of output (particle_lists),
 * followed by the substructure of the particle lines
 * -# Units of the quantities in the particle lines
 * -# SMASH version
 *
 * **File header for extended output** \n
 * If desired, the OSCAR2013 output can be extended
 * by additional particle properties. This requires enabling the extended
 * output in the configuration file, see the \key Extended switch in
 * \ref input_output_content_specific_ "content-specific output options" for
 * further details. The header of the extended OSCAR output is structured
 * identically to the non-extended version, but simply contains more columns
 * because of the additional entries:
 * <div class="fragment">
 * <div class="line"><span class="preprocessor">#!OSCAR2013 particle_lists
 *   t x y z mass p0 px py pz pdg
 *   ID charge ncoll form_time xsecfac proc_id_origin proc_type_origin
 *   time_last_coll pdg_mother1 pdg_mother2 baryon_number strangeness</span>
 * </div>
 * <div class="line"><span class="preprocessor">\# Units: fm fm fm fm GeV GeV
 * GeV GeV GeV none none e none fm none none none fm none none none
 * none</span></div> <div class="line"><span class="preprocessor">\#
 * SMASH_version</span></div>
 * </div>
 *
 * **Event block header**\n
 * The OSCAR2013 format is based on a block structure, where each block
 * corresponds to one interaction. Each block starts with a line formatted
 * as follows:
 * <div class="fragment">
 * <div class="line"> <span class="preprocessor">
 *  \# interaction in nin out nout rho density weight tot_weight partial
 *part_weight type proc_type </span></div>
 * </div>
 * where
 * \li \key nin: Number of ingoing particles (initial state particles)
 * \li \key nout: Number of outgoing particles (final state particles)
 * \li \key density: Density at the interaction point
 * \li \key tot_weight: Total weight of the interaction. This is the total
 *     cross section in case of a scattering and the total decay width in case
 *     of a decay. If there is no weight for the specific process, e.g. a wall
 *     crossing, it's value is 0.0
 * \li \key part_weight: The partial weight of the interaction. This is the
 *     specific weight for the chosen final state
 * \li \key proc_type: The type of the underlying process. See \ref
 *     doxypage_output_process_types for possible types
 *
 * Note, that "interaction", "in", "out", "rho", "weight", "partial" and "type"
 * are no variables, but words that are printed.
 *
 * <h2> Oscar1999 </h2>
 *
 * Oscar1999 is an ASCII (text) human-readable output following the OSCAR 1999
 * standard. The format specifics are the following:
 *
 * **File header**
 * \code
 * # OSC1999A
 * # full_event_history
 * # smash <version>
 * # Block format:
 * # nin nout density tot_weight part_weight proc_type
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
 * **Block header**\n
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
 * value is 0.0
 * \li \key part_weight: The partial weight of the interaction. This is the
 * specific weight for the chosen final state
 * \li \key proc_type: The type of the underlying process. See
 * \ref doxypage_output_process_types for possible types
 *
 * If the \key Print_Start_End option is set (see \ref
 * input_output_content_specific_ "content-specific output options" for
 * details), (nin, nout) = (0, Nparticles) in the
 * initial timestep and (nin, nout) = (Nparticles, 0) in the final timestep.
 **/

/*!\Userguide
 * \page doxypage_output_ascii
 * The \c ASCII format follows the general block structure of the \ref
 * doxypage_output_oscar, but offers more flexibility with the particle line
 * quantities written in the file. It is available for the \c %Particles,
 * \c Collisions, \c Dileptons, and \c Photons output contents (see \ref
 * doxypage_output), creating files with the extension <em>.dat</em>.
 * This format is useful to decrease storage usage.
 * \n
 *
 * <table>
 * <tr><th>Key(s)<th>C++ type <th> Description<th></tr>
 * <tr>
 * <td>\key t, \key x, \key y, \key z
 * <td> \c double <td>Space-time coordinates
 * </tr>
 * <tr>
 * <td>\key mass
 * <td>\c double
 * <td>Particle's rest-mass
 * </tr>
 * <tr>
 * <td>\key p0, \key px, \key py, \key pz
 * <td> \c double <td>Energy and 3-momentum
 * </tr>
 * <tr>
 * <td>\key pdg
 * <td>\c int32_t
 * <td>PDG code of the particle (see http://pdg.lbl.gov/). It contains all
 *     quantum numbers and uniquely identifies its type
 * </tr>
 * <tr>
 * <td>\key ID, \key id
 * <td> \c int32_t
 * <td>Particle identifier in terms of an integer. It is unique for every
 *     particle in the event. \key ID is used in the OSCAR 2013 standard, while
 *     \key id is used in OSCAR 1999
 * </tr>
 * <tr>
 * <td>\key charge
 * <td>\c int32_t
 * <td>Electric charge of the particle in units of the elementary charge \f$e\f$
 * </tr>
 * <tr>
 * <td>\key ncoll
 * <td>\c int32_t
 * <td>Number of collisions the particle has undergone
 * </tr>
 * <tr>
 * <td>\key form_time
 * <td>\c double
 * <td>Formation time of the particle
 * </tr>
 * <tr>
 * <td>\key xsecfac
 * <td>\c double
 * <td>Cross section scaling factor (if the particles are not yet fully formed
 *     at the time of interaction, the cross section for the underlying process
 *     is scaled down by the cross section scaling factor)
 * </tr>
 * <tr>
 * <td>\key proc_id_origin
 * <td>\c int32_t
 * <td>ID of the process of the particle's last interaction
 * </tr>
 * <tr>
 * <td>\key proc_type_origin
 * <td>\c int32_t
 * <td>Type of the last process the particle has undergone. The possible process
 *     types are listed in \ref doxypage_output_process_types
 * </tr>
 * <tr>
 * <td>
 * \key time_last_coll
 * <td>\c double
 * <td>Time of the particle's last interaction (except wall crossings)
 * </tr>
 * <tr>
 * <td>\key pdg_mother1
 * <td>\c int32_t
 * <td>PDG code of the 1st mother particle (0 in case the particle is sampled in
 *     a thermal bubble. It is not updated by elastic scatterings)
 * </tr>
 * <tr>
 * <td>\key pdg_mother2
 * <td>\c int32_t
 * <td>PDG code of the 2nd mother particle (0 in case the particle results from
 *     the decay of a resonance or the appearance of a thermal bubble. In the
 *     former case, \key pdg_mother1 is the PDG code of this resonance. It is
 *     not updated by elastic scatterings)
 * </tr>
 * <tr>
 * <td>\key baryon_number
 * <td>\c int32_t
 * <td>Baryon number of the particle: 1 for baryons, -1 for anti-baryons and 0
 *     for mesons
 * </tr>
 * <tr>
 * <td>\key strangeness
 * <td>\c int32_t
 * <td>Net-strangeness of the particles
 * </tr>
 * <tr>
 * <td>\key 0
 * <td> -
 * <td>Prints a column of 0 (for compatibility with OSCAR 1999)
 * </tr>
 * <tr>
 * <td>\key tau
 * <td>\c double
 * <td>Hyperbolic time \f$\tau=\sqrt{t^2-z^2}\f$
 * </tr>
 * <tr>
 * <td>\key eta,\key eta_s
 * <td>\c double
 * <td>Spacetime rapidity \f$\eta_s=\frac{1}{2}\log\frac{t+z}{t-z}\f$
 * </tr>
 * <tr>
 * <td>\key mt
 * <td>\c double
 * <td>Transverse mass \f$m_\perp=\sqrt{p_0^2-p_z^2}=\sqrt{m^2+p_x^2+p_y^2}\f$
 * </tr>
 * <tr>
 * <td>\key Rap,\key y_rap
 * <td>\c double
 * <td>Momentum rapidity
 *     \f$y_\mathrm{rap}=\frac{1}{2}\log\frac{p_0+p_z}{p_0-p_z}\f$
 * </tr>
 * </table>
 *
 * \attention Not all combinations of quantities are allowed and, in particular:
 *  - `id` and `ID` cannot be given together;
 *  - the same quantity cannot be repeated.
 *
 * **Example**
 *
 * If one is interested, for example, in the rate of production/annihilation of
 * resonances, for events divided into centrality classes based on the charged
 * particle yield at midrapidity, the config file would contain the following:
 *\verbatim
 Output:
     Particles:
         Format:     ["ASCII"]
         Quantities: ["p0","pz","pdg","charge"]
         Only_Final: IfNotEmpty
     Collisions:
         Format:     ["ASCII"]
         Quantities: ["t","pdg","ID","pdg_mother1","pdg_mother2"]
 \endverbatim
 *
 * Then, the output files would have the following headers:
 *
 * *particle_lists.dat*
 * \code
 * #!ASCII particle_lists p0 pz pdg charge
 * # Units: GeV GeV none e
 * # SMASH_version
 * \endcode
 *
 * *full_event_history.dat*
 * \code
 * #!ASCII full_event_history t pdg ID pdg_mother1 pdg_mother2
 * # Units: fm none none none none
 * # SMASH_version
 * \endcode
 **/

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::write_particledata(
    const ParticleData &data) {
  std::fprintf(file_.get(), "%s", formatter_.particle_line(data).c_str());
}

template <OscarOutputFormat Format, int Contents>
void OscarOutput<Format, Contents>::write(const ToASCII::type &buffer) {
  std::fprintf(file_.get(), "%s", buffer.c_str());
}

namespace {
/**
 * Helper function that creates the oscar output with the format selected by
 * create_oscar_output.
 *
 * \tparam Contents Determines what information will be written to the output
 * \param[in] path Path of output
 * \param[in] name (File)name of ouput
 * \param[in] modern_format Use the 1999 or 2013 format
 * \param[in] extended_format Whether the format is extended
 * \param[in] custom_format Whether the output has user-defined quantities
 * \param[in] quantities The user-defined quantities
 *
 * \return Unique pointer to oscar output
 */
template <int Contents>
std::unique_ptr<OutputInterface> create_selected_format(
    const std::filesystem::path &path, const std::string &name,
    bool modern_format, bool extended_format, bool custom_format,
    const std::vector<std::string> &quantities) {
  if (custom_format) {
    return std::make_unique<OscarOutput<ASCII, Contents>>(path, name,
                                                          quantities);
  } else {
    if (modern_format && extended_format) {
      return std::make_unique<OscarOutput<OscarFormat2013Extended, Contents>>(
          path, name);
    } else if (modern_format && !extended_format) {
      return std::make_unique<OscarOutput<OscarFormat2013, Contents>>(path,
                                                                      name);
    } else if (!modern_format && !extended_format) {
      return std::make_unique<OscarOutput<OscarFormat1999, Contents>>(path,
                                                                      name);
    } else {
      // Only remaining possibility: (!modern_format && extended_format)
      logg[LOutput].warn() << "There is no extended Oscar1999 format, creating "
                              "a regular Oscar1999 output instead.";
      return std::make_unique<OscarOutput<OscarFormat1999, Contents>>(path,
                                                                      name);
    }
  }
}
}  // unnamed namespace

std::unique_ptr<OutputInterface> create_oscar_output(
    const std::string &format, const std::string &content,
    const std::filesystem::path &path, const OutputParameters &out_par) {
  if (format != "Oscar2013" && format != "Oscar1999" && format != "ASCII") {
    throw std::invalid_argument("Creating Oscar output: unknown format");
  }
  const bool modern_format = (format == "Oscar2013");
  const bool custom_format = (format == "ASCII");
  const auto &quantities = custom_format ? out_par.quantities.at(content)
                                         : std::vector<std::string>{};

  if (content == "Particles") {
    if (out_par.part_only_final == OutputOnlyFinal::Yes) {
      return create_selected_format<OscarParticlesAtEventend>(
          path, "particle_lists", modern_format, out_par.part_extended,
          custom_format, quantities);
    } else if (out_par.part_only_final == OutputOnlyFinal::IfNotEmpty) {
      return create_selected_format<OscarParticlesAtEventendIfNotEmpty>(
          path, "particle_lists", modern_format, out_par.part_extended,
          custom_format, quantities);
    } else {  // out_par.part_only_final == OutputOnlyFinal::No
      return create_selected_format<OscarTimesteps | OscarAtEventstart |
                                    OscarParticlesAtEventend>(
          path, "particle_lists", modern_format, out_par.part_extended,
          custom_format, quantities);
    }
  } else if (content == "Collisions") {
    if (out_par.coll_printstartend) {
      return create_selected_format<OscarInteractions | OscarAtEventstart |
                                    OscarParticlesAtEventend>(
          path, "full_event_history", modern_format, out_par.coll_extended,
          custom_format, quantities);
    } else {
      return create_selected_format<OscarInteractions>(
          path, "full_event_history", modern_format, out_par.coll_extended,
          custom_format, quantities);
    }
  } else if (content == "Dileptons") {
    return create_selected_format<OscarInteractions>(
        path, "Dileptons", modern_format, out_par.dil_extended, custom_format,
        quantities);
  } else if (content == "Photons") {
    return create_selected_format<OscarInteractions>(
        path, "Photons", modern_format, out_par.photons_extended, custom_format,
        quantities);
  } else if (content == "Initial_Conditions") {
    return create_selected_format<OscarParticlesIC | OscarAtEventstart>(
        path, "SMASH_IC", modern_format, out_par.ic_extended, custom_format,
        quantities);
  }

  throw std::invalid_argument("Create_oscar_output got unknown content.");
}

}  // namespace smash
