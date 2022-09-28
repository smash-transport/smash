/*
 *
 *    Copyright (c) 2013-2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/experiment.h"

#include <cstdint>

#include "smash/boxmodus.h"
#include "smash/collidermodus.h"
#include "smash/listmodus.h"
#include "smash/spheremodus.h"

namespace smash {

/* ExperimentBase carries everything that is needed for the evolution */
ExperimentPtr ExperimentBase::create(Configuration &config,
                                     const std::filesystem::path &output_path) {
  if (!std::filesystem::exists(output_path)) {
    throw NonExistingOutputPathRequest("The requested output path (" +
                                       output_path.string() +
                                       ") does not exist.");
  }
  logg[LExperiment].trace() << SMASH_SOURCE_LOCATION;
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
   * \li \key ListBox - For given external particle list in the Box.
   */

  /*!\Userguide
   * \page input_modi_ Modi
   * \li \subpage input_modi_collider_
   * \li \subpage input_modi_sphere_
   * \li \subpage input_modi_box_
   * \li \subpage input_modi_list_
   * \li \subpage input_modi_listbox_
   */
  const std::string modus_chooser = config.read({"General", "Modus"});
  logg[LExperiment].debug() << "Modus for this calculation: " << modus_chooser;

  if (modus_chooser == "Box") {
    return std::make_unique<Experiment<BoxModus>>(config, output_path);
  } else if (modus_chooser == "List") {
    return std::make_unique<Experiment<ListModus>>(config, output_path);
  } else if (modus_chooser == "ListBox") {
    return std::make_unique<Experiment<ListBoxModus>>(config, output_path);
  } else if (modus_chooser == "Collider") {
    return std::make_unique<Experiment<ColliderModus>>(config, output_path);
  } else if (modus_chooser == "Sphere") {
    return std::make_unique<Experiment<SphereModus>>(config, output_path);
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
 * \key Ensembles (int, optional, default = 1): \n
 * Number of parallel ensembles in the simulation.
 *
 * An ensemble is an instance of the system, and without mean-field potentials
 * it is practically equivalent to a completely separate and uncorrelated event.
 * Each ensemble is an independent simulation: initialization, collisions,
 * decays, box wall crossings, and propagation of particles is performed
 * independently within each ensemble.
 *
 * However, the densities and mean-field potentials are computed as averages
 * over all ensembles (within a given event). This process can be also viewed as
 * calculating densities and mean-fields by summing over particles in all
 * ensembles combined, where each particle carries a fraction 1/n_ensembles of
 * its "real" charge. Such technique is called *parallel ensemble* technique. It
 * increases statistics necessary for precise density calculation without
 * increasing the number of collisions, which is not the case in the *full
 * ensemble* method (see below). Because of this, the parallel ensembles
 * technique is computationally faster than the full ensemble technique.
 *
 * \key Testparticles (int, optional, default = 1): \n
 * Number of test-particles per real particle in the simulation.
 *
 * Amount of initial sampled particles is increased by this factor,
 * while all cross sections are decreased by this factor. In this
 * way mean free path does not change. Larger number of testparticles
 * helps to reduce spurious effects of geometric collision criterion
 * (see \iref{Cheng:2001dz}). It also reduces correlations related
 * to collisions and decays (but not the ones related to mean fields),
 * therefore the larger the number of testparticles, the closer the results
 * of the simulations should be to the solution of Boltzmann equation.
 * These advantages come at a cost of computational time.
 *
 * Testparticles are a way to increase statistics necessary for
 * precise density calculation, which is why they are needed for mean field
 * potentials. The technique of using testparticles for mean field
 * is called *full ensemble* technique. The number of collisions (and
 * consequently the simulation time) scales as square of the number of
 * testparticles, and that is why full ensemble is slower than parallel
 * ensemble.
 *
 * \key Derivatives_Mode (string, optional, default = "Covariant Gaussian"): \n
 * The mode of calculating the gradients, for example gradients of baryon
 * current. Currently SMASH supports two derivatives modes: "Covariant Gaussian"
 * and "Finite difference". Covariant Gaussian derivatives can be used when
 * Covariant Gaussian smearing is used; they are Lorentz covariant, but they do
 * not calculate the time derivative of the current properly. The "Finite
 * difference" mode requires using the lattice, and the derivatives are
 * calculated based on finite differences of a given quantity at adjacent
 * lattice nodes; this mode is more numerically efficient.

 * \key Rest_Frame_Density_Derivatives_Mode (string, optional, default = "Off"):
 * \n
 * The mode of calculating the gradients of currents, decides whether the rest
 * frame density derivatives are copmuted (these derivatives are needed for the
 * VDF potentials, but not for the Skyrme potentials).
 *
 * \key Smearing_Mode (string, optional, default = "Covariant Gaussian"): \n
 * The mode of smearing for density calculation.
 *
 * Smearing is necessary to ensure a smooth gradient calculation, and it can be
 * thought of as smoothing out charge density fluctuations due to the finite
 * number of test-particles used. In general, this is done by distributing the
 * contribution to charge density from a given particle according to some
 * prescription. For example, in Gaussian smearing the charge density of a
 * particle is given by a Gaussian with some chosen width, centered at the
 * position of the particle; the Gaussian is normalized such that integrating
 * over the entire space yields the charge of the particle. In result, the
 * particle's charge is "smeared" over the space around it. Note that the case
 * with no smearing is recovered when the charge contribution from each particle
 * is taken to be a Dirac delta function centered at the position of the
 * particle.
 *
 * Currently, SMASH supports three smearing modes:
 * 1) "Covariant Gaussian": This smearing represents the charge density of a
 * particle as a Gaussian centered at the position of a particle; the user can
 * specify the width and the cutoff of the Gaussian (the employed Gaussians, in
 * principle non-zero over the entire available space, are "cut off" at some
 * distance r_cut from the particle to improve calculation time). This smearing
 * is Lorentz covariant which results in correct density profiles of
 * relativistic systems. The downside of the smearing is its long computation
 * time, as well as the fact that when the density is added to lattice nodes, it
 * is done so by Euler approximation (using the density value at the lattice
 * node), which does not conserve the number of particles on the lattice.
 * 2) "Triangular": This smearing requires lattice; it represents the charge
 * density of a particle in a given space direction as a "triangle" peaking at
 * the particle's position and linearly decreasing over a specified range. The
 * user specifies the range of the smearing in units of lattice spacings. This
 * smearing is relatively fast, and it does conserve the number of particles on
 * the lattice (due to the fact that the Euler integration is exact for a linear
 * function).
 * 3) "Discrete": This smearing requires lattice; the easiest of all smearing
 * modes, it adds a specified portion of the particle's charge density to a node
 * closest to the particle's position, and distributes the remainder evenly
 * among the 6 nearest neighbor nodes. The user specifies the weight given to
 * the center node; for example, if this weight is 1/3, then each of the six
 * nearest neighbor nodes gets 1/9 of the particle's charge. This smearing is
 * extremely fast, but is also rather coarse and requires using a large number
 * of test-particles to produce smooth gradients.
 *
 *
 * \key Gaussian_Sigma (double, optional, default = 1.0): \n
 * Parameter for Covariant Gaussian smearing: Width of gaussians that represent
 * Wigner density of particles, in fm.
 *
 * \key Gauss_Cutoff_In_Sigma (double, optional, default = 4.0): \n
 * Parameter for Covariant Gaussian smearing: Distance in sigma at which
 * gaussian is considered 0.
 *
 * \key Triangular_Range (double, optional, default = 2.0): \n
 * Parameter for Triangular smearing: Half of the base of a symmetric triangle
 * that represents particle density, in units of lattice spacings.
 *
 * \key Discrete_Weight (double, optional, default = 0.333333): \n
 * Parameter for Discrete smearing: Weight given to particle density at the
 * the center node; cannot be smaller than 1./7. (the boundary case of 1./7.
 * results in an even distribution of particle's density over the center node
 * and 6 neighboring nodes).
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
 * \key Output_Times (doubles, optional, no default): \n
 * Explicitly defines the the times where output is generated in the form of
 * a list. Cannot be used in combination with Output_Interval. Output times
 * outside the simulation time are ignored. The following example will produce
 * output at event start, event end and at the specified times as long as they
 * are within the simulation time.
 *\verbatim
 Output:
     Output_Times: [-0.1, 0.0, 1.0, 2.0, 10.0]
 \endverbatim
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
 *                  Oscar1999, VTK, HepMC_asciiv3 and HepMC_treeroot formats):
 \n
 *   \li \key true - Print extended information for each particle \n
 *   \li \key false - Regular output for each particle
 *
 *   \key Only_Final (string, optional, default = Yes, incompatible with
                      VTK, HepMC_asciiv3 and HepMC_treeroot formats): \n
 *   \li \key Yes - Print only final particle list \n
 *   \li \key IfNotEmpty - Print only final particle list, but only if event
 *                         is not empty (i.e. any collisions happened between
 *                         projectile and target). Useful to save disk space. \n
 *   \li \key No - Particle list at output interval including initial time \n
 * \n
 * - \b Collisions (VTK not available) \n
 *   \key Extended (bool, optional, default = false, incompatible with
 *          Oscar1999, HepMC_asciiv3 and HepMC_treeroot formats): \n
 *   \li \key true - Print extended information for each particle \n
 *   \li \key false - Regular output for each particle
 *
 *   \key Print_Start_End (bool, optional, default = false, incompatible with
 *                  Root, HepMC_asciiv3 and HepMC_treeroot formats): \n
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
 *   \key Extended (bool, optional, default = false, incompatible with
 *                  Oscar1999 format): \n
 *   \li \key true - Print extended information for each particle \n
 *   \li \key false - Regular output for each particle \n
 * \n
 * - \b Initial_Conditions (Oscar1999, Oscar2013, binary, ROOT and special ASCII
 * IC (\ref IC_output_user_guide_) formats)\n
 *   \key Proper_Time (double, optional, default = nuclei passing time, if
 *   nuclei passing time > \key Lower_Bound, else \key Lower_Bound):
 *   Proper time at which hypersurface is created \n
 *   \key Lower_Bound (double, optional, default = 0.5 fm): Lower bound for the
 *    IC proper time if \key Proper_Time is not provided.\n
 *   \key Rapidity_Cut (double, optional, default = no cut): If set, employ a
 *                 rapidity cut for particles contributing to the initial
 *                 conditions for hydrodynamics. A positive value is expected
 *                 and the cut is employed symmetrically around 0. Only
 *                 particles characterized by
 *                 - \key Rapidity_Cut < y < \key Rapidity_Cut are printed to
 *                 the output file.
 *   \key pT_Cut (double, optional, default = no cut): If set, employ a
 *                 transverse momentum cut for particles contributing to the
 *                 initial conditions for hydrodynamics. A positive value is
 *                 expected. Only particles characterized by
 *                 0 < pT < \key pT_Cut are printed to the output file.
 *   \key Extended (bool, optional, default = false, incompatible with
 *                  Oscar1999, ROOT and ASCII format):\n
 *   \li \key true - Print extended information for each particle
 *   \li \key false - Regular output for each particle \n
 * \n
 * - \b Rivet (Only YODA format)\n
 *   See \ref rivet_output_user_guide_ for more information
 * \n
 * - \b Coulomb (Only VTK format)\n
 *   No content-specific output options \n
 * \n
 * \anchor Thermodynamics
 * - \b Thermodynamics \n
 *   The user can print thermodynamical quantities:
 *   \li On the spatial lattice to vtk output. Note, that this output requires
 *       a lattice. This lattice needs to be enabled in the conguration file
 *       and is regulated by the options of \ref input_lattice_. See
 *       \ref output_vtk_lattice_ for further information.
 *   \li On the spatial lattice to ASCII output. Note, that this output requires
 *       a lattice. This lattice needs to be enabled in the conguration file
 *       and is regulated by the options of \ref input_lattice_. See
 *       \ref thermodyn_lattice_output_ for further information.
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
 *   \li \key "j_QBS" - Electric (Q), baryonic (B) and strange (S) currents
 *      \f$j^{\mu}_{QBS}(t,x,y,z) \f$; note that all currents are given in
 *      units of "number of charges"; multiply the electric current by the
 *      elementary charge \f$\sqrt{4 \pi \alpha_{EM}} \f$ for charge units.
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
 *   by the volume of the box to restore units to the correct ones. \n
 *   \n
 *
 *   \anchor onlypart
 *   \key Only_Participants (bool, optional, default = false): \n
 *   If set to true, only participants are included in the computation of the
 *   energy momentum tensor and of the Eckart currents. In this context,
 *   a hadron is considered as a participant if it had at least one collision.
 *   When using Potentials this option must be either left unset or set to
 *   false. The reason behing this limitation is that in this case hadrons
 *   can influence the evolution of the system even without collisions.
 *
 * \n
 * \page configuring_output_ Output Configuration
 * **Example: Configuring the SMASH Output**\n
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
         Only_Final: No
     Collisions:
         Format:    ["Oscar2013"]
         Extended: True
         Print_Start_End: True
 \endverbatim
 *
 * In addition, the photon and dilepton output can be enabled as follows, where
 * the dilepton output is generated in extended "Oscar2013" and "Binary" format
 * and the photon output is printed in "Oscar2013" format.
 *\verbatim
     Dileptons:
         Format:    ["Oscar2013", "Binary"]
         Extended: True
     Photons:
         Format:    ["Oscar2013"]
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
 * simulations. The corresponding output provides the particle list on a
 * hypersurface of constant proper time. If desired, the proper time can be set
 * manually from the configuration file (in the following example at \f$\tau =
 * 1.5 \f$ fm). If not provided, the default proper time corresponds to the
 * moment when both nuclei have entirely passed through each other, while this
 * proper time is greater than 0.5 fm. Else it is set to \f$\tau = 0.5 \f$ fm.\n
 * The initial conditions output can be enabled as follows:
 *\verbatim
     Initial_Conditions:
         Format:    ["ASCII", "Oscar1999", "Oscar2013", "Binary", "ROOT"]
         Extended: False
         Proper_Time: 1.5
 \endverbatim
 * The HepMC_asciiv3 and/or HepMC_treeroot ouputs are enabled by specifying
 * these output options under Particles or Collisions depdening on the content
 * wanted.
 *\verbatim
 Output:
     Particles:
         Format:          ["HepMC_asciiv3","HepMC_treeroot"]
     Collisions:
         Format:          ["HepMC_asciiv3","HepMC_treeroot"]
 \endverbatim
 * If a lattice is configured and coulomb potentials are enabled, a VTK output
 * for the electric and magnetic fields is available. It can be obtained by
 * adding the following to the output section of the configuration:
 *\verbatim
     Coulomb:
         Format:   ["VTK"]
 \endverbatim
 */

ExperimentParameters create_experiment_parameters(Configuration &config) {
  logg[LExperiment].trace() << SMASH_SOURCE_LOCATION;

  const int ntest = config.take({"General", "Testparticles"}, 1);
  if (ntest <= 0) {
    throw std::invalid_argument("Testparticle number should be positive!");
  }

  // sets whether to consider only participants in thermodynamic outputs or not
  const bool only_participants =
      config.take({"Output", "Thermodynamics", "Only_Participants"}, false);

  if (only_participants && config.has_value({"Potentials"})) {
    throw std::invalid_argument(
        "Only_Participants option cannot be "
        "set to True when using Potentials.");
  }

  const std::string modus_chooser = config.take({"General", "Modus"});
  // remove config maps of unused Modi
  config.remove_all_entries_in_section_but_one(modus_chooser, {"Modi"});

  double box_length = -1.0;
  if (config.has_value({"Modi", "Box", "Length"})) {
    box_length = config.read({"Modi", "Box", "Length"});
  }

  if (config.has_value({"Modi", "ListBox", "Length"})) {
    box_length = config.read({"Modi", "ListBox", "Length"});
  }

  /* If this Delta_Time option is absent (this can be for timestepless mode)
   * just assign 1.0 fm/c, reasonable value will be set at event initialization
   */
  const double dt = config.take({"General", "Delta_Time"}, 1.);
  const double t_end = config.read({"General", "End_Time"});

  // Enforce a small time step, if the box modus is used
  if (box_length > 0.0 && dt > box_length / 10.0) {
    throw std::invalid_argument(
        "Please decrease the timestep size. "
        "A value of (dt < l_box / 10) is recommended in the boxmodus.");
  }

  // define output clock
  std::unique_ptr<Clock> output_clock = nullptr;
  if (config.has_value({"Output", "Output_Times"})) {
    if (config.has_value({"Output", "Output_Interval"})) {
      throw std::invalid_argument(
          "Please specify either Output_Interval or Output_Times");
    }
    std::vector<double> output_times = config.take({"Output", "Output_Times"});
    // Add an output time larger than the end time so that the next time is
    // always defined during the time evolution
    output_times.push_back(t_end + 1.);
    output_clock = std::make_unique<CustomClock>(output_times);
  } else {
    const double output_dt = config.take({"Output", "Output_Interval"}, t_end);
    output_clock = std::make_unique<UniformClock>(0.0, output_dt);
  }

  // Add proper error messages if photons are not configured properly.
  // 1) Missing Photon config section.
  if (config.has_value({"Output", "Photons"}) &&
      (!config.has_value({"Collision_Term", "Photons"}))) {
    throw std::invalid_argument(
        "Photon output is enabled although photon production is disabled. "
        "Photon production can be configured in the \"Photon\" subsection "
        "of the \"Collision_Term\".");
  }

  // 2) Missing Photon output section.
  bool missing_output_2to2 = false;
  bool missing_output_brems = false;
  if (!(config.has_value({"Output", "Photons"}))) {
    if (config.has_value({"Collision_Term", "Photons", "2to2_Scatterings"})) {
      missing_output_2to2 =
          config.read({"Collision_Term", "Photons", "2to2_Scatterings"});
    }
    if (config.has_value({"Collision_Term", "Photons", "Bremsstrahlung"})) {
      missing_output_brems =
          config.read({"Collision_Term", "Photons", "Bremsstrahlung"});
    }

    if (missing_output_2to2 || missing_output_brems) {
      throw std::invalid_argument(
          "Photon output is disabled although photon production is enabled. "
          "Please enable the photon output.");
    }
  }

  // Add proper error messages if dileptons are not configured properly.
  // 1) Missing Dilepton config section.
  if (config.has_value({"Output", "Dileptons"}) &&
      (!config.has_value({"Collision_Term", "Dileptons"}))) {
    throw std::invalid_argument(
        "Dilepton output is enabled although dilepton production is disabled. "
        "Dilepton production can be configured in the \"Dileptons\" subsection "
        "of the \"Collision_Term\".");
  }

  // 2) Missing Dilepton output section.
  bool missing_output_decays = false;
  if (!(config.has_value({"Output", "Dileptons"}))) {
    if (config.has_value({"Collision_Term", "Dileptons", "Decays"})) {
      missing_output_decays =
          config.read({"Collision_Term", "Dileptons", "Decays"});
    }

    if (missing_output_decays) {
      throw std::invalid_argument(
          "Dilepton output is disabled although dilepton production is "
          "enabled. "
          "Please enable the dilepton output.");
    }
  }

  auto config_coll = config.extract_sub_configuration(
      {"Collision_Term"}, Configuration::GetEmpty::Yes);
  /* Elastic collisions between the nucleons with the square root s
   * below low_snn_cut are excluded. */
  const double low_snn_cut =
      config_coll.take({"Elastic_NN_Cutoff_Sqrts"}, 1.98);
  const auto proton = ParticleType::try_find(pdg::p);
  const auto pion = ParticleType::try_find(pdg::pi_z);
  if (proton && pion &&
      low_snn_cut > proton->mass() + proton->mass() + pion->mass()) {
    logg[LExperiment].warn("The cut-off should be below the threshold energy",
                           " of the process: NN to NNpi");
  }
  const bool potential_affect_threshold =
      config.take({"Lattice", "Potentials_Affect_Thresholds"}, false);
  const double scale_xs = config_coll.take({"Cross_Section_Scaling"}, 1.0);

  const auto criterion =
      config_coll.take({"Collision_Criterion"}, CollisionCriterion::Covariant);

  if (config_coll.has_value({"Fixed_Min_Cell_Length"}) &&
      criterion != CollisionCriterion::Stochastic) {
    throw std::invalid_argument(
        "Only use a fixed minimal cell length with the stochastic collision "
        "criterion.");
  }
  if (config_coll.has_value({"Maximum_Cross_Section"}) &&
      criterion == CollisionCriterion::Stochastic) {
    throw std::invalid_argument(
        "Only use maximum cross section with the "
        "geometric collision criterion. Use Fixed_Min_Cell_Length to change "
        "the grid "
        "size for the stochastic criterion.");
  }

  /**
   * The maximum around 200 mb occurs in the Delta peak of the pi+p
   * cross section. Many SMASH cross sections diverge at the threshold,
   * these divergent parts are effectively cut off. If deuteron production
   * via d' is considered, then the default should be increased to 2000 mb
   * to function correctly (see \iref{Oliinychenko:2018ugs}). If the cross
   * sections are globally scaled, the maximum cross section is also scaled.
   */
  const double maximum_cross_section_default =
      ParticleType::exists("d'") ? 2000.0 : 200.0;

  bool cll_in_nucleus =
      config.take({"Modi", "Collider", "Collisions_Within_Nucleus"}, false);
  double maximum_cross_section = config_coll.take(
      {"Maximum_Cross_Section"}, maximum_cross_section_default);
  maximum_cross_section *= scale_xs;
  return {std::make_unique<UniformClock>(0.0, dt),
          std::move(output_clock),
          config.take({"General", "Ensembles"}, 1),
          ntest,
          config.take({"General", "Derivatives_Mode"},
                      DerivativesMode::CovariantGaussian),
          config.has_value({"Potentials", "VDF"})
              ? RestFrameDensityDerivativesMode::On
              : RestFrameDensityDerivativesMode::Off,
          config.take({"General", "Field_Derivatives_Mode"},
                      FieldDerivativesMode::ChainRule),
          config.take({"General", "Smearing_Mode"},
                      SmearingMode::CovariantGaussian),
          config.take({"General", "Gaussian_Sigma"}, 1.),
          config.take({"General", "Gauss_Cutoff_In_Sigma"}, 4.),
          config.take({"General", "Discrete_Weight"}, 1. / 3.0),
          config.take({"General", "Triangular_Range"}, 2.0),
          criterion,
          config_coll.take({"Two_to_One"}, true),
          config_coll.take({"Included_2to2"}, ReactionsBitSet().set()),
          config_coll.take({"Multi_Particle_Reactions"},
                           MultiParticleReactionsBitSet().reset()),
          config_coll.take({"Strings"}, modus_chooser != "Box"),
          config_coll.take({"Use_AQM"}, true),
          config_coll.take({"Resonance_Lifetime_Modifier"}, 1.),
          config_coll.take({"Strings_with_Probability"}, true),
          config_coll.take({"NNbar_Treatment"}, NNbarTreatment::Strings),
          low_snn_cut,
          potential_affect_threshold,
          box_length,
          maximum_cross_section,
          config_coll.take({"Fixed_Min_Cell_Length"}, 2.5),
          cll_in_nucleus,
          scale_xs,
          config_coll.take({"Additional_Elastic_Cross_Section"}, 0.0),
          only_participants,
          config_coll.take({"Include_Weak_And_EM_Decays_At_The_End"}, false)};
}

std::string format_measurements(const std::vector<Particles> &ensembles,
                                uint64_t scatterings_this_interval,
                                const QuantumNumbers &conserved_initial,
                                SystemTimePoint time_start, double time,
                                double E_mean_field,
                                double E_mean_field_initial) {
  const SystemTimeSpan elapsed_seconds = SystemClock::now() - time_start;

  const QuantumNumbers current_values(ensembles);
  const QuantumNumbers difference = current_values - conserved_initial;
  int total_particles = 0;
  for (const Particles &particles : ensembles) {
    total_particles += particles.size();
  }

  // Make sure there are no FPEs in case of IC output, were there will
  // eventually be no more particles in the system
  const double current_energy = current_values.momentum().x0();
  const double energy_per_part =
      (total_particles > 0) ? (current_energy + E_mean_field) / total_particles
                            : 0.0;

  std::ostringstream ss;
  // clang-format off
  ss << field<7, 3> << time
    // total kinetic energy in the system
     << field<11, 3> << current_energy
    // total mean field energy in the system
     << field<11, 3> << E_mean_field
    // total energy in the system
     << field<12, 3> << current_energy + E_mean_field
    // total energy per particle in the system
     << field<12, 6> << energy_per_part;
    // change in total energy per particle (unless IC output is enabled)
    if (total_particles == 0) {
     ss << field<13, 6> << "N/A";
    } else {
     ss << field<13, 6> << (difference.momentum().x0()
                            + E_mean_field - E_mean_field_initial)
                            / total_particles;
    }
    ss << field<14, 3> << scatterings_this_interval
     << field<10, 3> << total_particles
     << field<9, 3> << elapsed_seconds;
  // clang-format on
  return ss.str();
}

double calculate_mean_field_energy(
    const Potentials &potentials,
    RectangularLattice<smash::DensityOnLattice> &jmuB_lat,
    RectangularLattice<std::pair<ThreeVector, ThreeVector>> *em_lattice,
    const ExperimentParameters &parameters) {
  // basic parameters and variables
  const double V_cell = (jmuB_lat.cell_sizes())[0] *
                        (jmuB_lat.cell_sizes())[1] * (jmuB_lat.cell_sizes())[2];

  double E_mean_field = 0.0;
  double density_mean = 0.0;
  double density_variance = 0.0;

  /*
   * We anticipate having other options, like the vector DFT potentials, in the
   * future, hence we include checking which potentials are used.
   */
  if (potentials.use_skyrme()) {
    /*
     * Calculating the symmetry energy contribution to the total mean field
     * energy in the system is not implemented at this time.
     */
    if (potentials.use_symmetry() &&
        parameters.outputclock->current_time() == 0.0) {
      logg[LExperiment].warn()
          << "Note:"
          << "\nSymmetry energy is not included in the mean field calculation."
          << "\n\n";
    }

    /*
     * Skyrme potential parameters:
     * C1GeV are the Skyrme coefficients converted to GeV,
     * b1 are the powers of the baryon number density entering the expression
     * for the energy density of the system. Note that these exponents are
     * larger by 1 than those for the energy of a particle (which are used in
     * Potentials class). The formula for a total mean field energy due to a
     * Skyrme potential is E_MF = \sum_i (C_i/b_i) ( n_B^b_i )/( n_0^(b_i - 1) )
     * where nB is the local rest frame baryon number density and n_0 is the
     * saturation density. Then the single particle potential follows from
     * V = d E_MF / d n_B .
     */
    double C1GeV = (potentials.skyrme_a()) / 1000.0;
    double C2GeV = (potentials.skyrme_b()) / 1000.0;
    double b1 = 2.0;
    double b2 = (potentials.skyrme_tau()) + 1.0;

    /*
     * Note: calculating the mean field only works if lattice is used.
     * We iterate over the nodes of the baryon density lattice to sum their
     * contributions to the total mean field.
     */
    int number_of_nodes = 0;
    double lattice_mean_field_total = 0.0;

    for (auto &node : jmuB_lat) {
      number_of_nodes++;
      // the rest frame density
      double rhoB = node.rho();
      // the computational frame density
      const double j0B = node.jmu_net().x0();

      const double abs_rhoB = std::abs(rhoB);
      if (abs_rhoB < very_small_double) {
        continue;
      }
      density_mean += j0B;
      density_variance += j0B * j0B;

      /*
       * The mean-field energy for the Skyrme potential. Note: this expression
       * is only exact in the rest frame, and is expected to significantly
       * deviate from the correct value for systems that are considerably
       * relativistic. Note: symmetry energy is not taken into the account.
       *
       * TODO: Add symmetry energy.
       */
      double mean_field_contribution_1 = (C1GeV / b1) * std::pow(abs_rhoB, b1) /
                                         std::pow(nuclear_density, b1 - 1);
      double mean_field_contribution_2 = (C2GeV / b2) * std::pow(abs_rhoB, b2) /
                                         std::pow(nuclear_density, b2 - 1);

      lattice_mean_field_total +=
          V_cell * (mean_field_contribution_1 + mean_field_contribution_2);
    }

    // logging statistical properties of the density calculation
    density_mean = density_mean / number_of_nodes;
    density_variance = density_variance / number_of_nodes;
    double density_scaled_variance =
        std::sqrt(density_variance - density_mean * density_mean) /
        density_mean;
    logg[LExperiment].debug() << "\t\t\t\t\t";
    logg[LExperiment].debug()
        << "\n\t\t\t\t\t            density mean = " << density_mean;
    logg[LExperiment].debug()
        << "\n\t\t\t\t\t density scaled variance = " << density_scaled_variance;
    logg[LExperiment].debug()
        << "\n\t\t\t\t\t        total mean_field = "
        << lattice_mean_field_total * parameters.testparticles *
               parameters.n_ensembles
        << "\n";

    E_mean_field = lattice_mean_field_total;
  }  // if (potentials.use_skyrme())

  if (potentials.use_vdf()) {
    /*
     * Safety check:
     * Calculating the symmetry energy contribution to the total mean field
     * energy in the system is not implemented at this time.
     */
    if (potentials.use_symmetry() &&
        parameters.outputclock->current_time() == 0.0) {
      logg[LExperiment].error()
          << "\nSymmetry energy is not included in the VDF mean-field "
             "calculation"
          << "\nas VDF potentials haven't been fitted with symmetry energy."
          << "\n\n";
    }

    /*
     * The total mean-field energy density due to a VDF potential is
     * E_MF = \sum_i C_i rho^(b_i - 2) *
     *                 * [j_0^2 -  rho^2 * (b_i - 1)/b_i] / rho_0^(b_i - 1)
     * where j_0 is the local computational frame baryon density, rho is the
     * local rest frame baryon density, and rho_0 is the saturation density.
     */

    // saturation density of nuclear matter specified in the VDF parameters
    double rhoB_0 = potentials.saturation_density();

    /*
     * Note: calculating the mean field only works if lattice is used.
     * We iterate over the nodes of the baryon density lattice to sum their
     * contributions to the total mean field.
     */
    int number_of_nodes = 0;
    double lattice_mean_field_total = 0.0;

    for (auto &node : jmuB_lat) {
      number_of_nodes++;
      // the rest frame density
      double rhoB = node.rho();
      // the computational frame density
      const double j0B = node.jmu_net().x0();
      double abs_rhoB = std::abs(rhoB);
      density_mean += j0B;
      density_variance += j0B * j0B;

      /*
       * The mean-field energy for the VDF potential. This expression is correct
       * in any frame, and in the rest frame conforms to the Skyrme mean-field
       * energy (if same coefficients and powers are used).
       */
      // in order to prevent dividing by zero in case any b_i < 2.0
      if (abs_rhoB < very_small_double) {
        abs_rhoB = very_small_double;
      }
      double mean_field_contribution = 0.0;
      for (int i = 0; i < potentials.number_of_terms(); i++) {
        mean_field_contribution +=
            potentials.coeffs()[i] *
            std::pow(abs_rhoB, potentials.powers()[i] - 2.0) *
            (j0B * j0B -
             ((potentials.powers()[i] - 1.0) / potentials.powers()[i]) *
                 abs_rhoB * abs_rhoB) /
            std::pow(rhoB_0, potentials.powers()[i] - 1.0);
      }
      lattice_mean_field_total += V_cell * mean_field_contribution;
    }

    // logging statistical properties of the density calculation
    density_mean = density_mean / number_of_nodes;
    density_variance = density_variance / number_of_nodes;
    double density_scaled_variance =
        std::sqrt(density_variance - density_mean * density_mean) /
        density_mean;
    logg[LExperiment].debug() << "\t\t\t\t\t";
    logg[LExperiment].debug()
        << "\n\t\t\t\t\t            density mean = " << density_mean;
    logg[LExperiment].debug()
        << "\n\t\t\t\t\t density scaled variance = " << density_scaled_variance;
    logg[LExperiment].debug()
        << "\n\t\t\t\t\t        total mean_field = "
        << lattice_mean_field_total * parameters.testparticles *
               parameters.n_ensembles
        << "\n";

    E_mean_field = lattice_mean_field_total;
  }

  double electromagnetic_potential = 0.0;
  if (potentials.use_coulomb() && em_lattice) {
    // Use cell volume of electromagnetic fields lattice even though it should
    // be the same as for net-baryon density
    double V_cell_em = em_lattice->cell_sizes()[0] *
                       em_lattice->cell_sizes()[1] *
                       em_lattice->cell_sizes()[2];
    for (auto &fields : *em_lattice) {
      // Energy is 0.5 * int E^2 + B^2 dV
      electromagnetic_potential +=
          hbarc * 0.5 * V_cell_em * (fields.first.sqr() + fields.second.sqr());
    }
  }
  logg[LExperiment].debug() << "Total energy in electromagnetic field  = "
                            << electromagnetic_potential;
  E_mean_field += electromagnetic_potential;
  /*
   * E_mean_field is multiplied by the number of testparticles per particle and
   * the number of parallel ensembles because the total kinetic energy tracked
   * is that of all particles in the simulation, including test-particles and/or
   * ensembles, and so this way is more consistent.
   */
  E_mean_field =
      E_mean_field * parameters.testparticles * parameters.n_ensembles;

  return E_mean_field;
}

EventInfo fill_event_info(const std::vector<Particles> &ensembles,
                          double E_mean_field, double modus_impact_parameter,
                          const ExperimentParameters &parameters,
                          bool projectile_target_interact,
                          bool kinematic_cut_for_SMASH_IC) {
  const QuantumNumbers current_values(ensembles);
  const double E_kinetic_total = current_values.momentum().x0();
  const double E_total = E_kinetic_total + E_mean_field;

  EventInfo event_info{modus_impact_parameter,
                       parameters.box_length,
                       parameters.outputclock->current_time(),
                       E_kinetic_total,
                       E_mean_field,
                       E_total,
                       parameters.testparticles,
                       parameters.n_ensembles,
                       !projectile_target_interact,
                       kinematic_cut_for_SMASH_IC};
  return event_info;
}

}  // namespace smash
