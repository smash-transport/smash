/*
 *
 *    Copyright (c) 2012-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/experiment.h"

#include <cstdint>

#include "smash/actions.h"
#include "smash/boxmodus.h"
#include "smash/collidermodus.h"
#include "smash/cxx14compat.h"
#include "smash/fourvector.h"
#include "smash/listmodus.h"
#include "smash/spheremodus.h"

namespace smash {

/* ExperimentBase carries everything that is needed for the evolution */
ExperimentPtr ExperimentBase::create(Configuration config,
                                     const bf::path &output_path) {
  const auto &log = logger<LogArea::Experiment>();
  log.trace() << source_location;
  /*!\Userguide
   * \page input_general_ General
   * \key Modus (string, required): \n
   * Selects a modus for the calculation, e.g.\ infinite matter
   * calculation, collision of two particles or collision of nuclei. The modus
   * will be configured in \ref input_modi_. Recognized values are:
   * \li \key Collider - For collisions of nuclei or compound objects. See \ref
   *     \ColliderModus
   * \li \key Sphere - For calculations of the expansion of a thermalized
   * sphere. See \ref \SphereModus \li \key Box - For infinite matter
   * calculation in a rectangular box. See \ref \BoxModus \li \key List - For
   * given external particle list. See \ref \ListModus
   */

  /*!\Userguide
   * \page input_modi_ Modi
   * \li \subpage input_modi_collider_
   * \li \subpage input_modi_sphere_
   * \li \subpage input_modi_box_
   * \li \subpage input_modi_list_
   */
  const std::string modus_chooser = config.read({"General", "Modus"});
  log.debug() << "Modus for this calculation: " << modus_chooser;

  if (modus_chooser == "Box") {
    return make_unique<Experiment<BoxModus>>(config, output_path);
  } else if (modus_chooser == "List") {
    return make_unique<Experiment<ListModus>>(config, output_path);
  } else if (modus_chooser == "Collider") {
    return make_unique<Experiment<ColliderModus>>(config, output_path);
  } else if (modus_chooser == "Sphere") {
    return make_unique<Experiment<SphereModus>>(config, output_path);
  } else {
    throw InvalidModusRequest("Invalid Modus (" + modus_chooser +
                              ") requested from ExperimentBase::create.");
  }
}

/*!\Userguide
 * \page input_general_ General
 * \key Delta_Time (double, optional, default: 1.0): \n
 * Fixed time step at which the collision-finding grid is recreated, and, if
 * potentials are on, momenta are updated according to the equations of motion.
 * The collision-finding grid finds all the collisions from time
 * t_{beginning_of_timestep} until time t_{beginning_of_timestep} + Delta_Time,
 * and puts them into a vector. The collisions are then sorted in order of
 * occurrence, and particles are propagated from collision to collision. After
 * each performed collision, additional collisions are found for outgoing
 * particles and merged into the sorted vector.
 *
 * If potentials are on, the Delta_Time should be small enough, typically
 * around 0.1 fm/c. However, if potentials are off, it can be arbitrarily
 * large. In this case it only influences the runtime, but not physics.
 * If Time_Step_Mode = None is chosen, then the user-provided value of
 * Delta_Time is ignored and Delta_Time is set to the End_Time.
 *
 * \key Testparticles (int, optional, default = 1): \n
 * How many test particles per real particle should be simulated.
 *
 * \key Gaussian_Sigma (double, optional, default = 1.0): \n
 * Width of gaussians that represent Wigner density of particles, in fm.
 *
 * \key Gauss_Cutoff_In_Sigma (double, optional, default = 4.0): \n
 * Distance in sigma at which gaussian is considered 0.
 *
 * \page input_output_options_ Output Configuration
 *
 * Description of options
 * ---------------------
 * To produce a certain output content it is necessary to explicitly configure
 * it in the Output section of the configuration file. This means, that the
 * Output section needs to contain a subsection for the desired output.
 * Aditionally, there are general output configuration parameters. \n
 * \n
 * ### General output configuration parameters:
 * \key Output_Interval (double, optional, default = End_Time): \n
 * Defines the period of intermediate output of the status of the simulated
 * system in Standard Output and other output formats which support this
 * functionality.
 *
 * \key Density_Type (string, optional, default = "none"): \n
 * Determines which kind of density is printed into the headers of the
 * collision files.
 * Possible values:\n
 * \li \key "hadron" - Total hadronic density
 * \li \key "baryon" - Net baryon density
 * \li \key "baryonic isospin" - Baryonic isospin density
 * \li \key "pion" - Pion density
 * \li \key "none" - Do not calculate density, print 0.0
 *
 * \n
 * ### Format configuration independently of the specific output content
 * Further options are defined for every single output content
 * (see \ref output_contents_ "output contents" for the list of
 * possible contents). Independently of the content, it is always necessary
 * to provide the format in which the output should be generated.
 *
 * \key Format (list of formats, optional, default = [ ]):\n
 * List of formats for writing particular content.
 * Possible formats for every content are listed and described in
 * \ref output_contents_ "output contents". List of available formats is
 * \ref list_of_output_formats "here".
 * \n
 * Besides the universal \key Format option, there are also content-specific
 output
 * options that are listed below.
 *
 * ### Content-specific output options
 * \anchor output_content_specific_options_
 *
 * - \b Particles \n
 *   \key Extended (bool, optional, default = false, incompatible with
 *                  Oscar1999, VTK and Root formats): \n
 *   \li \key true - Print extended information for each particle \n
 *   \li \key false - Regular output for each particle
 *
 *   \key Only_Final (bool, optional, default = true, incompatible with
                      VTK format): \n
 *   \li \key true - Print only final particle list \n
 *   \li \key false - Particle list at output interval including initial time \n
 * \n
 * - \b Collisions (VTK not available) \n
 *   \key Extended (bool, optional, default = false, incompatible with
 *                  Oscar1999 and Root formats): \n
 *   \li \key true - Print extended information for each particle \n
 *   \li \key false - Regular output for each particle
 *
 *   \key Print_Start_End (bool, optional, default = false, incompatible with
 *                  Root format): \n
 *   \li \key true - Initial and final particle list is printed out \n
 *   \li \key false - Initial and final particle list is not printed out \n
 * \n
 * - \b Dileptons (Only Oscar1999, Oscar2013 and binary formats) \n
 *   \key Extended (bool, optional, default = false, incompatible with
 *                  Oscar1999 format): \n
 *   \li \key true - Print extended information for each particle \n
 *   \li \key false - Regular output for each particle \n
 * \n
 * - \b Photons (Only Oscar1999, Oscar2013 and binary formats) \n
 *   \key Fractions (int, required):
 *   Number of fractional photons sampled per single perturbatively produced
 *   photon. See \ref input_photons for further information. \n
 *   \key Extended (bool, optional, default = false, incompatible with
 *                  Oscar1999 format): \n
 *   \li \key true - Print extended information for each particle \n
 *   \li \key false - Regular output for each particle \n
 * \n
 * - \b Initial_Conditions (Oscar1999, Oscar2013, binary, ROOT and special ASCII
 * IC (\ref IC_output_user_guide_) formats)\n
 *   \key Proper_Time (double, optional, default = time after which nuclei
 *                     have passed through each other): Proper time at which
 *                     hypersurface is created \n
 *   \key Extended (bool, optional, default = false, incompatible with
 *                  Oscar1999, ROOT and ASCII format): \n
 *   \li \key true - Print extended information for each particle \n
 *   \li \key false - Regular output for each particle \n
 * \n
 * \anchor Thermodynamics
 * - \b Thermodynamics \n
 *   The user can print thermodynamical quantities:
 *   \li On the spatial lattice to vtk output. Note, that this output requires
 *       a lattice. This lattice needs to be enabled in the conguration file
 *       and is regulated by the options of \ref input_lattice_. See
 *       \ref output_vtk_lattice_ for further information.
 *   \li At a given point to ASCII output. See
 *       \ref thermodyn_output_user_guide_ for further information.
 *   \li Averaged over all particles to ASCII output. See
 *       \ref thermodyn_output_user_guide_ for further information.
 *
 *  \key Type (string, optional, default = \key "baryon"): \n
 *  Particle type taken into consideration, "baryon" corresponds to "net
 baryon".
 *   \li \key "hadron"
 *   \li \key "baryon"
 *   \li \key "baryonic isospin"
 *   \li \key "pion"
 *   \li \key "none"
 *   \li \key "total isospin"
 *
 *   \key Quantities (list of thermodynamic quantities, optional, default = [
 ]):\n
 *   List of thermodynamic quantities that are printed to the output. Possible
 *   quantities are:
 *   \li \key "rho_eckart" - Eckart rest frame density
 *   \li \key "tmn" - Energy-momentum tensor \f$T^{\mu\nu}(t,x,y,z) \f$
 *   \li \key "tmn_landau" - Energy-momentum tensor in the Landau rest frame.
 *      This tensor is computed by boosting \f$T^{\mu\nu}(t,x,y,z) \f$
 *      to the local rest frame, where \f$T^{0i} \f$ = 0.
 *   \li \key "landau_velocity" - Velocity of the Landau rest frame.
 *      The velocity is obtained from the energy-momentum tensor
 *      \f$T^{\mu\nu}(t,x,y,z) \f$ by solving the generalized eigenvalue
 *      equation \f$(T^{\mu\nu} - \lambda g^{\mu\nu})u_{\mu}=0 \f$.
 *
 *   \key Position (list of 3 doubles, optional, default = [0.0, 0.0, 0.0]): \n
 *   Point, at which thermodynamic quantities are computed.
 *
 *   \key Smearing (bool, optional, default = true): \n
 *   Using Gaussian smearing for computing thermodynamic quantities or not.
 *   This triggers whether thermodynamic quantities are evaluated at a fixed
 *   point (\key true) or summed over all particles (\key false).
 *   \li \key true - smearing applied
 *   \li \key false - smearing not applied
 *
 *   The contribution to the energy-momentum tensor and current (be it electric,
 *   baryonic or strange) from a single particle in its rest frame is:
 *   \f[\begin{eqnarray} j^{\mu} = B \frac{p_0^{\mu}}{p_0^0} W \\
 *   T^{\mu \nu} = \frac{p_0^{\mu}p_0^{\nu}}{p_0^0} W \end{eqnarray}\f]
 *   with B being the charge of interest and W being the weight given to this
 *   particle. Normally, if one computes thermodynamic quantities at a point,
 *   smearing should be applied, and then W takes on the following shape:
 *   \f[W = (2 \pi \sigma^2)^{-3/2} exp \left(- \frac{(\mathbf{r}
 *   - \mathbf{r_0(t)})^2}{2\sigma^2} \right)\f]
 *   It can however be useful to compute the thermo-
 *   dynamic quantities of all particles in a box with W = 1, which
 *   would correspond to \key "Smearing: false". Note that using this option
 *   changes the units of the thermodynamic quantities, as they are no longer
 *   spatially normalized. One should divide this quantity by
 *   by the volume of the box to restore units to the correct ones.
 *
 * \n
 * \page configuring_output_ Output Configuration
 * Example: Configuring the SMASH Output
 * --------------
 * The following example configures the output to be printed in an interval of
 * 1 fm and with the net baryon density being printed to the header.
 * The particles output is generated in "Oscar1999", VTK and "Root" format,
 * generating output for each time step. The collisions output is formatted
 * according to an extended "Oscar2013" format and the initial and final
 * particle lists are printed as well.
 *\verbatim
 Output:
     Output_Interval: 1.0
     Density_Type: "baryon"
     Particles:
         Format:    ["Oscar1999", "VTK", "Root"]
         Extended: False
         Only_Final: False
     Collisions:
         Format:    ["Oscar2013"]
         Extended: True
         Print_Start_End: True
 \endverbatim
 *
 * To further activate photons and dileptons in the SMASH simulation and to also
 * generate the output, the corresponding subsections need to be present in the
 * configuration file. In the following example, the dilepton output is
 * generated in extended "Oscar2013" and "Binary" format. The photon output
 * is printed in "Oscar2013" format while the calculation is performed with
 * 100 fractional photons.
 *\verbatim
     Dileptons:
         Format:    ["Oscar2013", "Binary"]
         Extended: True
     Photons:
         Format:    ["Oscar2013"]
         Fractions: 100
 \endverbatim
 *
 * Additionally, the thermodynsamics output can be activated. In this example,
 * thermodynamic output is activated for hadrons. The quanities that are printed
 * are the density in the Eckart rest frame and the energy momentum tensor in
 * the Landau rest frame. These quantities are printed at each time step for the
 * position (0,0,0). Gaussian smearing is not applied. The output is provided
 * in "ASCII" and "VTK" format.
 *\verbatim
     Thermodynamics:
         Format:    ["ASCII", "VTK"]
         Type: "hadron"
         Quantities:    ["rho_eckart", "tmn_landau"]
         Position:    [0.0, 0.0, 0.0]
         Smearing: False
 \endverbatim
 * SMASH can further be applied to extract initial conditions for hydrodynamic
 * simulations. The corresponding output, the "IC output" provides a
 * hypersurface of constant proper time. If desired, the proper time can be set
 * manually from the configuration file (in the following example at \f$\tau =
 * 0.5 \f$ fm). If not provided, the default proper time corresponds to the
 * moment when both nuclei have entirely passed through each other.\n
 * The initial conditions output can be enabled as follows:
 *\verbatim
     Initial_Conditions:
         Format:    ["ASCII", "Oscar1999", "Oscar2013", "Binary", "ROOT"]
         Extended: False
         Proper_Time: 0.5
 \endverbatim
 */

ExperimentParameters create_experiment_parameters(Configuration config) {
  const auto &log = logger<LogArea::Experiment>();
  log.trace() << source_location;

  const int ntest = config.take({"General", "Testparticles"}, 1);
  if (ntest <= 0) {
    throw std::invalid_argument("Testparticle number should be positive!");
  }

  const std::string modus_chooser = config.take({"General", "Modus"});
  // remove config maps of unused Modi
  config["Modi"].remove_all_but(modus_chooser);

  /* If this Delta_Time option is absent (this can be for timestepless mode)
   * just assign 1.0 fm/c, reasonable value will be set at event initialization
   */
  const double dt = config.take({"General", "Delta_Time"}, 1.);
  const double t_end = config.read({"General", "End_Time"});
  const double output_dt = config.take({"Output", "Output_Interval"}, t_end);
  auto config_coll = config["Collision_Term"];
  /* Elastic collisions between the nucleons with the square root s
   * below low_snn_cut are excluded. */
  const double low_snn_cut =
      config_coll.take({"Elastic_NN_Cutoff_Sqrts"}, 1.98);
  const auto proton = ParticleType::try_find(pdg::p);
  const auto pion = ParticleType::try_find(pdg::pi_z);
  if (proton && pion &&
      low_snn_cut > proton->mass() + proton->mass() + pion->mass()) {
    log.warn("The cut-off should be below the threshold energy",
             " of the process: NN to NNpi");
  }
  const bool potential_affect_threshold =
      config.take({"Lattice", "Potentials_Affect_Thresholds"}, false);
  return {{0., dt},
          {0.0, output_dt},
          ntest,
          config.take({"General", "Gaussian_Sigma"}, 1.),
          config.take({"General", "Gauss_Cutoff_In_Sigma"}, 4.),
          config_coll.take({"Two_to_One"}, true),
          config_coll.take({"Included_2to2"}, ReactionsBitSet().set()),
          config_coll.take({"Strings"}, modus_chooser != "Box"),
          config_coll.take({"Use_AQM"}, true),
          config_coll.take({"Strings_with_Probability"}, true),
          config_coll.take({"NNbar_Treatment"}, NNbarTreatment::Strings),
          config.has_value({"Output", "Photons"}),
          low_snn_cut,
          potential_affect_threshold};
}

std::string format_measurements(const Particles &particles,
                                uint64_t scatterings_this_interval,
                                const QuantumNumbers &conserved_initial,
                                SystemTimePoint time_start, double time) {
  const SystemTimeSpan elapsed_seconds = SystemClock::now() - time_start;

  const QuantumNumbers current_values(particles);
  const QuantumNumbers difference = conserved_initial - current_values;

  std::ostringstream ss;
  // clang-format off
  ss << field<8, 3> << time << field<13, 3> << difference.momentum().x0()
     << field<16, 3> << scatterings_this_interval
     << field<11, 3> << particles.size() << field<10, 3> << elapsed_seconds;
  // clang-format on
  return ss.str();
}

}  // namespace smash
