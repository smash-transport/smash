/*
 *
 *    Copyright (c) 2013-2024
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

  const std::string modus_chooser = config.read(InputKeys::gen_modus);
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
 * \page doxypage_input_conf_output
 *
 * ---
 *
 * \section config_output_examples Examples for configuring the SMASH output
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
         Format: ["Oscar1999", "VTK", "Root"]
         Extended: False
         Only_Final: No
     Collisions:
         Format: ["Oscar2013"]
         Extended: True
         Print_Start_End: True
 \endverbatim
 *
 * In addition, the photon and dilepton output can be enabled as follows, where
 * the dilepton output is generated in extended "Oscar2013" and the
 * corresponding binary version of it, while the photon output is printed in
 * "Oscar2013" format only.
 *\verbatim
 Output:
     Dileptons:
         Format: ["Oscar2013", "Oscar2013_bin"]
         Extended: True
     Photons:
         Format: ["Oscar2013"]
 \endverbatim
 *
 * Additionally, the thermodynamics output can be activated. In this example,
 * thermodynamic output is activated for hadrons. The quanities that are printed
 * are the density in the Eckart rest frame and the energy momentum tensor in
 * the Landau rest frame. These quantities are printed at each time step for the
 * position (0,0,0). Gaussian smearing is not applied. The output is provided
 * in "ASCII" and "VTK" format.
 *\verbatim
 Output:
     Thermodynamics:
         Format: ["ASCII", "VTK"]
         Type: "hadron"
         Quantities: ["rho_eckart", "tmn_landau"]
         Position: [0.0, 0.0, 0.0]
         Smearing: False
 \endverbatim
 * SMASH can further be applied to extract initial conditions for hydrodynamic
 * simulations, either in a hypersurface of constant hyperbolic time or based on
 * the local energy density, as controlled by the key <tt>\ref key_MC_IC_type_
 * "Modi: Collider: Initial_Conditions: Type"</tt>. The "For_vHLLE" output
 * format is only available in conjunction with the iso-tau hypersurface
 * option.\n The initial conditions output can be enabled as follows (for
 * further possible formats, see the \ref list_of_output_formats
 * "corresponding documentation"):
 *\verbatim
 Output:
     Initial_Conditions:
         Format: ["For_vHLLE", "Oscar2013", "Oscar2013_bin", "Root"]
         Extended: False
 \endverbatim
 * The HepMC_asciiv3 and/or HepMC_treeroot ouputs are enabled by specifying
 * these output options under Particles or Collisions depdening on the content
 * wanted.
 *\verbatim
 Output:
     Particles:
         Format: ["HepMC_asciiv3","HepMC_treeroot"]
     Collisions:
         Format: ["HepMC_asciiv3","HepMC_treeroot"]
 \endverbatim
 * If a lattice is configured and coulomb potentials are enabled, a VTK output
 * for the electric and magnetic fields is available. It can be obtained by
 * adding the following to the output section of the configuration:
 *\verbatim
     Coulomb:
         Format: ["VTK"]
 \endverbatim
 */

ExperimentParameters create_experiment_parameters(Configuration &config) {
  logg[LExperiment].trace() << SMASH_SOURCE_LOCATION;

  const int ntest = config.take(InputKeys::gen_testparticles);
  if (ntest <= 0) {
    throw std::invalid_argument("Testparticle number should be positive!");
  }

  // sets whether to consider only participants in thermodynamic outputs or not
  const bool only_participants =
      config.take(InputKeys::output_thermodynamics_onlyParticipants);

  if (only_participants && config.has_section(InputSections::potentials)) {
    throw std::invalid_argument(
        "Only_Participants option cannot be "
        "set to True when using Potentials.");
  }

  const std::string modus_chooser = config.take(InputKeys::gen_modus);
  // remove config maps of unused Modi
  config.remove_all_entries_in_section_but_one(modus_chooser, {"Modi"});

  double box_length = -1.0;
  if (config.has_value(InputKeys::modi_box_length)) {
    box_length = config.read(InputKeys::modi_box_length);
  }
  if (config.has_value(InputKeys::modi_listBox_length)) {
    box_length = config.read(InputKeys::modi_listBox_length);
  }

  /* If this Delta_Time option is absent (this can be for timestepless mode)
   * just assign 1.0 fm, reasonable value will be set at event initialization
   */
  const double dt = config.take(InputKeys::gen_deltaTime);
  if (dt <= 0.) {
    throw std::invalid_argument("Delta_Time cannot be zero or negative.");
  }

  const double t_end = config.read(InputKeys::gen_endTime);
  if (t_end <= 0.) {
    throw std::invalid_argument("End_Time cannot be zero or negative.");
  }

  // Enforce a small time step, if the box modus is used
  if (box_length > 0.0 && dt > box_length / 10.0) {
    throw std::invalid_argument(
        "Please decrease the timestep size. "
        "A value of (dt <= l_box / 10) is necessary in the box modus.");
  }

  // define output clock
  std::unique_ptr<Clock> output_clock = nullptr;
  if (config.has_value(InputKeys::output_outputTimes)) {
    if (config.has_value(InputKeys::output_outputInterval)) {
      throw std::invalid_argument(
          "Please specify either Output_Interval or Output_Times");
    }
    std::vector<double> output_times =
        config.take(InputKeys::output_outputTimes);
    // Add an output time larger than the end time so that the next time is
    // always defined during the time evolution
    output_times.push_back(t_end + 1.);
    output_clock = std::make_unique<CustomClock>(output_times);
  } else {
    const double output_dt =
        config.take(InputKeys::output_outputInterval, t_end);
    if (output_dt <= 0.) {
      throw std::invalid_argument(
          "Output_Interval cannot be zero or negative.");
    }
    output_clock = std::make_unique<UniformClock>(0.0, output_dt, t_end);
  }

  // Add proper error messages if photons are not configured properly.
  // 1) Missing Photon config section.
  if (config.has_section(InputSections::o_photons) &&
      (!config.has_section(InputSections::c_photons))) {
    throw std::invalid_argument(
        "Photon output is enabled although photon production is disabled. "
        "Photon production can be configured in the \"Photon\" subsection "
        "of the \"Collision_Term\".");
  }

  // 2) Missing Photon output section.
  if (!(config.has_section(InputSections::o_photons))) {
    const bool missing_output_2to2 =
                   config.read(InputKeys::collTerm_photons_twoToTwoScatterings),
               missing_output_brems =
                   config.read(InputKeys::collTerm_photons_bremsstrahlung);
    if (missing_output_2to2 || missing_output_brems) {
      throw std::invalid_argument(
          "Photon output is disabled although photon production is enabled. "
          "Please enable the photon output.");
    }
  }

  // Add proper error messages if dileptons are not configured properly.
  // 1) Missing Dilepton config section.
  if (config.has_section(InputSections::o_dileptons) &&
      (!config.has_section(InputSections::c_dileptons))) {
    throw std::invalid_argument(
        "Dilepton output is enabled although dilepton production is disabled. "
        "Dilepton production can be configured in the \"Dileptons\" subsection "
        "of the \"Collision_Term\".");
  }

  // 2) Missing Dilepton output section.
  if (!(config.has_section(InputSections::o_dileptons))) {
    const bool missing_output_decays =
        config.read(InputKeys::collTerm_dileptons_decays);
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
      config.take(InputKeys::collTerm_elasticNNCutoffSqrts);
  const auto proton = ParticleType::try_find(pdg::p);
  const auto pion = ParticleType::try_find(pdg::pi_z);
  if (proton && pion &&
      low_snn_cut > proton->mass() + proton->mass() + pion->mass()) {
    logg[LExperiment].warn("The cut-off should be below the threshold energy",
                           " of the process: NN to NNpi");
  }
  const bool potential_affect_threshold =
      config.take(InputKeys::lattice_potentialsAffectThreshold);
  const double scale_xs = config.take(InputKeys::collTerm_crossSectionScaling);

  const auto criterion = config.take(InputKeys::collTerm_collisionCriterion);

  if (config.has_value(InputKeys::collTerm_fixedMinCellLength) &&
      criterion != CollisionCriterion::Stochastic) {
    throw std::invalid_argument(
        "Only use a fixed minimal cell length with the stochastic collision "
        "criterion.");
  }
  if (config.has_value(InputKeys::collTerm_maximumCrossSection) &&
      criterion == CollisionCriterion::Stochastic) {
    throw std::invalid_argument(
        "Only use maximum cross section with the geometric collision "
        "criterion. Use Fixed_Min_Cell_Length to change the grid size for the "
        "stochastic criterion.");
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

  double maximum_cross_section = config.take(
      InputKeys::collTerm_maximumCrossSection, maximum_cross_section_default);
  maximum_cross_section *= scale_xs;
  return {std::make_unique<UniformClock>(0.0, dt, t_end),
          std::move(output_clock),
          config.take(InputKeys::gen_ensembles),
          ntest,
          config.take(InputKeys::gen_derivativesMode),
          config.has_section(InputSections::p_vdf)
              ? RestFrameDensityDerivativesMode::On
              : RestFrameDensityDerivativesMode::Off,
          config.take(InputKeys::gen_fieldDerivativesMode),
          config.take(InputKeys::gen_smearingMode),
          config.take(InputKeys::gen_smearingGaussianSigma),
          config.take(InputKeys::gen_smearingGaussCutoffInSigma),
          config.take(InputKeys::gen_smearingDiscreteWeight),
          config.take(InputKeys::gen_smearingTriangularRange),
          criterion,
          config.take(InputKeys::collTerm_twoToOne),
          config.take(InputKeys::collTerm_includedTwoToTwo),
          config.take(InputKeys::collTerm_multiParticleReactions),
          config.take(InputKeys::collTerm_strings, modus_chooser != "Box"),
          config.take(InputKeys::collTerm_resonanceLifetimeModifier),
          config.take(InputKeys::collTerm_nnbarTreatment),
          low_snn_cut,
          potential_affect_threshold,
          box_length,
          maximum_cross_section,
          config.take(InputKeys::collTerm_fixedMinCellLength),
          scale_xs,
          only_participants,
          config.take(InputKeys::collTerm_ignoreDecayWidthAtTheEnd),
          config.take(InputKeys::collTerm_decayInitial),
          std::nullopt};
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

void validate_and_adjust_particle_list(ParticleList &particle_list) {
  static bool warn_mass_discrepancy = true;
  static bool warn_off_shell_particle = true;
  for (auto it = particle_list.begin(); it != particle_list.end();) {
    auto &particle = *it;
    auto pdgcode = particle.pdgcode();
    try {
      // Convert Kaon-L or Kaon-S into K0 or Anti-K0 used in SMASH
      if (pdgcode == 0x310 || pdgcode == 0x130) {
        pdgcode = (random::uniform_int(0, 1) == 0) ? pdg::K_z : pdg::Kbar_z;
      }
      /* ATTENTION: It would be wrong to directly assign here the return value
       * to 'particle', because this would potentially also change its id and
       * process number, which in turn, might lead to actions to be discarded.
       * Here, only the particle momentum has to be adjusted and this is done
       * creating a new particle and using its momentum to set 'particle' one.
       * The position and momentum of the particle are checked for nan values.
       */
      auto valid_smash_particle =
          create_valid_smash_particle_matching_provided_quantities(
              pdgcode, particle.effective_mass(), particle.position(),
              particle.momentum(), LExperiment, warn_mass_discrepancy,
              warn_off_shell_particle);
      particle.set_4position(valid_smash_particle.position());
      particle.set_4momentum(valid_smash_particle.momentum());
      particle.set_cross_section_scaling_factor(
          valid_smash_particle.xsec_scaling_factor());
      it++;
    } catch (ParticleType::PdgNotFoundFailure &) {
      logg[LExperiment].warn()
          << "SMASH does not recognize pdg code " << pdgcode
          << " obtained from hadron list. This particle will be ignored.\n";
      it = particle_list.erase(it);
    }
  }
}

}  // namespace smash
