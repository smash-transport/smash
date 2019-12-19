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
  logg[LExperiment].trace() << source_location;
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
  logg[LExperiment].debug() << "Modus for this calculation: " << modus_chooser;

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
 * \key Output_Times (doubles, optinal, no default): \n
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
 *   \key Proper_Time (double, optional, default = nuclei passing time): Proper
 *       time at which hypersurface is created \n
 *   \key Extended (bool, optional, default = false, incompatible with
 *                  Oscar1999, ROOT and ASCII format):\n
 *   \li \key true - Print extended information for each particle
 *   \li \key false - Regular output for each particle
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
 * simulations. The corresponding output provides the particle list on a
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
  logg[LExperiment].trace() << source_location;

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
    output_clock = make_unique<CustomClock>(output_times);
  } else {
    const double output_dt = config.take({"Output", "Output_Interval"}, t_end);
    output_clock = make_unique<UniformClock>(0.0, output_dt);
  }

  auto config_coll = config["Collision_Term"];
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
  return {make_unique<UniformClock>(0.0, dt),
          std::move(output_clock),
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
                                SystemTimePoint time_start, double time,
                                double E_mean_field,
                                double E_mean_field_initial) {
  const SystemTimeSpan elapsed_seconds = SystemClock::now() - time_start;

  const QuantumNumbers current_values(particles);
  const QuantumNumbers difference = current_values - conserved_initial;

  std::ostringstream ss;
  // clang-format off
  ss << field<7, 3> << time
    // total kinetic energy in the system
     << field<11, 3> << current_values.momentum().x0()
    // total mean field energy in the system
     << field<11, 3> << E_mean_field
    // total energy in the system
     << field<12, 3> << current_values.momentum().x0() + E_mean_field
    // total energy per particle in the system
     << field<12, 6>
     << (current_values.momentum().x0() + E_mean_field)/particles.size()
    // change in energy per particle
     << field<13, 6> << (difference.momentum().x0()
                         + E_mean_field - E_mean_field_initial)
                        /particles.size()
     << field<14, 3> << scatterings_this_interval
     << field<10, 3> << particles.size()
     << field<9, 3> << elapsed_seconds;
  // clang-format on
  return ss.str();
}

double calculate_mean_field_energy(
    const Potentials &potentials,
    RectangularLattice<smash::DensityOnLattice> &jmu_B_lat,
    const ExperimentParameters &parameters) {
  // basic parameters and variables
  const double V_cell = (jmu_B_lat.cell_sizes())[0] *
                        (jmu_B_lat.cell_sizes())[1] *
                        (jmu_B_lat.cell_sizes())[2];

  double E_mean_field = 0.0;
  double density_mean = 0.0;
  double density_variance = 0.0;

  /*
   * We anticipate having other options, like the vector DFT potentials, in the
   * future, hence we include checking which potentials are used.
   */
  if (potentials.use_skyrme()) {
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

    for (auto &node : jmu_B_lat) {
      number_of_nodes++;
      // the rest frame density
      double nB = node.density();
      // the computational frame density
      const double j0 = node.jmu_net().x0();

      density_mean += j0;
      density_variance += j0 * j0;

      /*
       * The mean-field energy for the Skyrme potential. Note: this expression
       * is only exact in the rest frame, and is expected to significantly 
       * deviate from the correct value for systems that are considerably 
       * relativistic. Note: symmetry energy is not taken into the account.
       * 
       * TODO: Add symmetry energy.
       */
      double mean_field_contribution_1 = (C1GeV/b1) * std::pow(nB, b1) /
                                          std::pow(nuclear_density, b1 - 1);
      double mean_field_contribution_2 = (C2GeV/b2) * std::pow(nB, b2) /
                                          std::pow(nuclear_density, b2 - 1);

      lattice_mean_field_total +=
          V_cell * (mean_field_contribution_1 + mean_field_contribution_2);
    }

    // logging statistical properties of the density calculation
    density_mean = density_mean / number_of_nodes;
    density_variance = density_variance / number_of_nodes;
    double density_scaled_variance =
        sqrt(density_variance - density_mean * density_mean) / density_mean;
    logg[LExperiment].debug() << "\t\t\t\t\t";
    logg[LExperiment].debug()
        << "\n\t\t\t\t\t            density mean = " << density_mean;
    logg[LExperiment].debug()
        << "\n\t\t\t\t\t density scaled variance = " << density_scaled_variance;
    logg[LExperiment].debug()
        << "\n\t\t\t\t\t        total mean_field = "
        << lattice_mean_field_total * parameters.testparticles << "\n";

    /*
     * E_mean_field is multiplied by the number of testparticles because the
     * total kinetic energy tracked is that of all particles, including
     * testparticles, and so this is more consistent with the general paradigm.
     */
    E_mean_field = lattice_mean_field_total * parameters.testparticles;
  }

  return E_mean_field;
}

}  // namespace smash
