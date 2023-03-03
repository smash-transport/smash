/*
 *
 *    Copyright (c) 2013-2023
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
 * \n
 * \page doxypage_output_conf_examples
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
   * just assign 1.0 fm, reasonable value will be set at event initialization
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
    output_clock = std::make_unique<UniformClock>(0.0, output_dt, t_end);
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
  /* Elastic collisions between the nucleons with the square root s
   * below low_snn_cut are excluded. */
  const double low_snn_cut =
      config.take({"Collision_Term", "Elastic_NN_Cutoff_Sqrts"}, 1.98);
  const auto proton = ParticleType::try_find(pdg::p);
  const auto pion = ParticleType::try_find(pdg::pi_z);
  if (proton && pion &&
      low_snn_cut > proton->mass() + proton->mass() + pion->mass()) {
    logg[LExperiment].warn("The cut-off should be below the threshold energy",
                           " of the process: NN to NNpi");
  }
  const bool potential_affect_threshold =
      config.take({"Lattice", "Potentials_Affect_Thresholds"}, false);
  const double scale_xs =
      config.take({"Collision_Term", "Cross_Section_Scaling"}, 1.0);

  const auto criterion = config.take({"Collision_Term", "Collision_Criterion"},
                                     CollisionCriterion::Covariant);

  if (config.has_value({"Collision_Term", "Fixed_Min_Cell_Length"}) &&
      criterion != CollisionCriterion::Stochastic) {
    throw std::invalid_argument(
        "Only use a fixed minimal cell length with the stochastic collision "
        "criterion.");
  }
  if (config.has_value({"Collision_Term", "Maximum_Cross_Section"}) &&
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

  double maximum_cross_section =
      config.take({"Collision_Term", "Maximum_Cross_Section"},
                  maximum_cross_section_default);
  maximum_cross_section *= scale_xs;
  return {
      std::make_unique<UniformClock>(0.0, dt, t_end),
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
      config.take({"Collision_Term", "Two_to_One"}, true),
      config.take({"Collision_Term", "Included_2to2"}, ReactionsBitSet().set()),
      config.take({"Collision_Term", "Multi_Particle_Reactions"},
                  MultiParticleReactionsBitSet().reset()),
      config.take({"Collision_Term", "Strings"}, modus_chooser != "Box"),
      config.take({"Collision_Term", "Resonance_Lifetime_Modifier"}, 1.),
      config.take({"Collision_Term", "NNbar_Treatment"},
                  NNbarTreatment::Strings),
      low_snn_cut,
      potential_affect_threshold,
      box_length,
      maximum_cross_section,
      config.take({"Collision_Term", "Fixed_Min_Cell_Length"}, 2.5),
      scale_xs,
      only_participants,
      config.take({"Collision_Term", "Include_Weak_And_EM_Decays_At_The_End"},
                  false)};
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

ParticleList check_particle_list(ParticleList &particle_list,
                                 int &n_warns_precision,
                                 int &n_warns_mass_consistency) {
  // Check if the particles in the list are valid, perfom the same checks
  // as in ListModus::try_create_particle()
  // (see function documentation in listmodus.h)
  constexpr int max_warns_precision = 10, max_warn_mass_consistency = 10;
  ParticleList plist_checked;
  int pdgcode_a = 0;
  for (auto &particle : particle_list) {
    try {
      pdgcode_a = particle.pdgcode().get_decimal();

      // Convert Kaon-L or Kaon-S into K0 or Anti-K0 used in SMASH
      if (pdgcode_a == 310 || pdgcode_a == 130) {
        pdgcode_a = (random::uniform_int(0, 1) == 0) ? 311 : -311;
      }

      ParticleData new_p{ParticleType::find(PdgCode::from_decimal(pdgcode_a))};
      const FourVector p = particle.momentum();
      const double mass = p.abs();

      // Check mass of particles in list
      if (particle.type().is_stable() &&
          std::abs(mass - particle.pole_mass()) > really_small) {
        if (n_warns_precision < max_warns_precision) {
          logg[LExperiment].warn()
              << "Provided mass of " << particle.type().name() << " = " << mass
              << " [GeV] is inconsistent with SMASH value = "
              << particle.pole_mass() << ". Forcing E = sqrt(p^2 + m^2)"
              << ", where m is SMASH mass.";
          n_warns_precision++;
        } else if (n_warns_precision == max_warns_precision) {
          logg[LExperiment].warn(
              "Further warnings about SMASH mass versus input mass"
              " inconsistencies will be suppressed.");
          n_warns_precision++;
        }
        new_p.set_4momentum(mass, ThreeVector(p.x1(), p.x2(), p.x3()));
      } else {
        new_p.set_4momentum(FourVector(p.x0(), p.x1(), p.x2(), p.x3()));
      }

      // On-shell condition consistency check
      if (std::abs(particle.momentum().sqr() - mass * mass) > really_small) {
        if (n_warns_mass_consistency < max_warn_mass_consistency) {
          logg[LExperiment].warn()
              << "Provided 4-momentum " << particle.momentum() << " and "
              << " mass " << mass << " do not satisfy E^2 - p^2 = m^2."
              << " This may originate from the lack of numerical"
              << " precision in the input. Setting E to sqrt(p^2 + m^2).";
          n_warns_mass_consistency++;
        } else if (n_warns_mass_consistency == max_warn_mass_consistency) {
          logg[LExperiment].warn(
              "Further warnings about E != sqrt(p^2 + m^2) will"
              " be suppressed.");
          n_warns_mass_consistency++;
        }
        new_p.set_4momentum(mass, ThreeVector(p.x1(), p.x2(), p.x3()));
      }
      // Set spatial coordinates, they will later be backpropagated if needed
      const FourVector r = particle.position();
      new_p.set_4position(FourVector(r.x0(), r.x1(), r.x2(), r.x3()));
      new_p.set_cross_section_scaling_factor(1.0);
      plist_checked.push_back(new_p);
    } catch (ParticleType::PdgNotFoundFailure &) {
      logg[LExperiment].warn()
          << "SMASH does not recognize pdg code " << pdgcode_a
          << " obtained from hadron list. This particle will be ignored.\n";
    }
  }
  return plist_checked;
}

}  // namespace smash
