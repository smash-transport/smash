/*
 *    Copyright (c) 2012-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/collidermodus.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "smash/configuration.h"
#include "smash/customnucleus.h"
#include "smash/experimentparameters.h"
#include "smash/fluidizationaction.h"
#include "smash/fourvector.h"
#include "smash/icparameters.h"
#include "smash/input_keys.h"
#include "smash/logging.h"
#include "smash/nucleus.h"
#include "smash/random.h"

namespace smash {
static constexpr int LCollider = LogArea::Collider::id;
static constexpr int LInitialConditions = LogArea::InitialConditions::id;

ColliderModus::ColliderModus(Configuration modus_config,
                             const ExperimentParameters &params) {
  Configuration modus_cfg = modus_config.extract_complete_sub_configuration(
      InputSections::m_collider);
  // Get the reference frame for the collision calculation.
  frame_ = modus_cfg.take(InputKeys::modi_collider_calculationFrame);

  Configuration proj_cfg = modus_cfg.extract_complete_sub_configuration(
      InputSections::m_c_projectile);
  Configuration targ_cfg =
      modus_cfg.extract_complete_sub_configuration(InputSections::m_c_target);
  /* Needed to check if projectile and target in customnucleus are read from
   * the same input file.*/
  bool same_file = false;
  // Set up the projectile nucleus
  if (proj_cfg.has_section(InputSections::m_c_p_deformed)) {
    projectile_ =
        create_deformed_nucleus(proj_cfg, params.testparticles, "projectile");
  } else if (proj_cfg.has_section(InputSections::m_c_p_custom)) {
    same_file = same_inputfile(proj_cfg, targ_cfg);
    projectile_ = std::make_unique<CustomNucleus>(
        proj_cfg, params.testparticles, same_file);
  } else if (proj_cfg.has_section(InputSections::m_c_p_alphaClustered)) {
    logg[LCollider].info() << "Projectile is alpha-clustered with woods-saxon "
                              "parameters for the He-clusters listed below.";
    projectile_ = create_alphaclustered_nucleus(proj_cfg, params.testparticles,
                                                "projectile");
  } else {
    projectile_ = std::make_unique<Nucleus>(proj_cfg, params.testparticles);
  }
  if (projectile_->size() < 1) {
    throw ColliderEmpty("Input Error: Projectile nucleus is empty.");
  }
  projectile_->set_label(BelongsTo::Projectile);

  // Set up the target nucleus
  if (targ_cfg.has_section(InputSections::m_c_t_deformed)) {
    target_ = create_deformed_nucleus(targ_cfg, params.testparticles, "target");
  } else if (targ_cfg.has_section(InputSections::m_c_t_custom)) {
    target_ = std::make_unique<CustomNucleus>(targ_cfg, params.testparticles,
                                              same_file);
  } else if (targ_cfg.has_section(InputSections::m_c_t_alphaClustered)) {
    logg[LCollider].info() << "Target is alpha-clustered with woods-saxon "
                              "parameters for the He-clusters listed below.";
    target_ =
        create_alphaclustered_nucleus(targ_cfg, params.testparticles, "target");
  } else {
    target_ = std::make_unique<Nucleus>(targ_cfg, params.testparticles);
  }
  if (target_->size() < 1) {
    throw ColliderEmpty("Input Error: Target nucleus is empty.");
  }
  target_->set_label(BelongsTo::Target);

  // Get the Fermi-Motion input (off, on, frozen)
  fermi_motion_ = modus_cfg.take(InputKeys::modi_collider_fermiMotion);
  if (fermi_motion_ == FermiMotion::On) {
    logg[LCollider].info() << "Fermi motion is ON.";
  } else if (fermi_motion_ == FermiMotion::Frozen) {
    logg[LCollider].info() << "FROZEN Fermi motion is on.";
  } else if (fermi_motion_ == FermiMotion::Off) {
    logg[LCollider].info() << "Fermi motion is OFF.";
  }

  // Get the total nucleus-nucleus collision energy. Since there is
  // no meaningful choice for a default energy, we require the user to
  // give one (and only one) energy input from the available options.
  int energy_input = 0;
  const double mass_projec = projectile_->mass();
  const double mass_target = target_->mass();
  // average mass of a particle in that nucleus
  const double mass_a =
      projectile_->mass() / projectile_->number_of_particles();
  const double mass_b = target_->mass() / target_->number_of_particles();
  // Option 1: Center of mass energy.
  if (modus_cfg.has_value(InputKeys::modi_collider_sqrtSNN)) {
    sqrt_s_NN_ = modus_cfg.take(InputKeys::modi_collider_sqrtSNN);
    // Check that input satisfies the lower bound (everything at rest).
    if (sqrt_s_NN_ <= mass_a + mass_b) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: sqrt(s_NN) is not larger than masses:\n" +
          std::to_string(sqrt_s_NN_) + " GeV <= " + std::to_string(mass_a) +
          " GeV + " + std::to_string(mass_b) + " GeV.");
    }
    // Set the total nucleus-nucleus collision energy.
    total_s_ = (sqrt_s_NN_ * sqrt_s_NN_ - mass_a * mass_a - mass_b * mass_b) *
                   mass_projec * mass_target / (mass_a * mass_b) +
               mass_projec * mass_projec + mass_target * mass_target;
    energy_input++;
  }
  /* Option 2: Total energy per nucleon of the projectile nucleus
   * (target at rest).  */
  if (modus_cfg.has_value(InputKeys::modi_collider_eTot)) {
    const double e_tot = modus_cfg.take(InputKeys::modi_collider_eTot);
    if (e_tot < 0) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: "
          "E_Tot must be nonnegative.");
    }
    // Set the total nucleus-nucleus collision energy.
    total_s_ = s_from_Etot(e_tot * projectile_->number_of_particles(),
                           mass_projec, mass_target);
    sqrt_s_NN_ = std::sqrt(s_from_Etot(e_tot, mass_a, mass_b));
    energy_input++;
  }
  /* Option 3: Kinetic energy per nucleon of the projectile nucleus
   * (target at rest).  */
  if (modus_cfg.has_value(InputKeys::modi_collider_eKin)) {
    const double e_kin = modus_cfg.take(InputKeys::modi_collider_eKin);
    if (e_kin < 0) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: "
          "E_Kin must be nonnegative.");
    }
    // Set the total nucleus-nucleus collision energy.
    total_s_ = s_from_Ekin(e_kin * projectile_->number_of_particles(),
                           mass_projec, mass_target);
    sqrt_s_NN_ = std::sqrt(s_from_Ekin(e_kin, mass_a, mass_b));
    energy_input++;
  }
  // Option 4: Momentum of the projectile nucleus (target at rest).
  if (modus_cfg.has_value(InputKeys::modi_collider_pLab)) {
    const double p_lab = modus_cfg.take(InputKeys::modi_collider_pLab);
    if (p_lab < 0) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: "
          "P_Lab must be nonnegative.");
    }
    // Set the total nucleus-nucleus collision energy.
    total_s_ = s_from_plab(p_lab * projectile_->number_of_particles(),
                           mass_projec, mass_target);
    sqrt_s_NN_ = std::sqrt(s_from_plab(p_lab, mass_a, mass_b));
    energy_input++;
  }
  // Option 5: Total energy per nucleon of _each_ beam
  if (proj_cfg.has_value(InputKeys::modi_collider_projectile_eTot) &&
      targ_cfg.has_value(InputKeys::modi_collider_target_eTot)) {
    const double e_tot_p =
        proj_cfg.take(InputKeys::modi_collider_projectile_eTot);
    const double e_tot_t = targ_cfg.take(InputKeys::modi_collider_target_eTot);
    if (e_tot_p < 0 || e_tot_t < 0) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: "
          "E_Tot must be nonnegative.");
    }
    total_s_ = s_from_Etot(e_tot_p * projectile_->number_of_particles(),
                           e_tot_t * target_->number_of_particles(),
                           mass_projec, mass_target);
    sqrt_s_NN_ = std::sqrt(s_from_Ekin(e_tot_p, e_tot_t, mass_a, mass_b));
    energy_input++;
  }
  // Option 6: Kinetic energy per nucleon of _each_ beam
  if (proj_cfg.has_value(InputKeys::modi_collider_projectile_eKin) &&
      targ_cfg.has_value(InputKeys::modi_collider_target_eKin)) {
    const double e_kin_p =
        proj_cfg.take(InputKeys::modi_collider_projectile_eKin);
    const double e_kin_t = targ_cfg.take(InputKeys::modi_collider_target_eKin);
    if (e_kin_p < 0 || e_kin_t < 0) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: "
          "E_Kin must be nonnegative.");
    }
    total_s_ = s_from_Ekin(e_kin_p * projectile_->number_of_particles(),
                           e_kin_t * target_->number_of_particles(),
                           mass_projec, mass_target);
    sqrt_s_NN_ = std::sqrt(s_from_Ekin(e_kin_p, e_kin_t, mass_a, mass_b));
    energy_input++;
  }
  // Option 7: Momentum per nucleon of _each_ beam
  if (proj_cfg.has_value(InputKeys::modi_collider_projectile_pLab) &&
      targ_cfg.has_value(InputKeys::modi_collider_target_pLab)) {
    const double p_lab_p =
        proj_cfg.take(InputKeys::modi_collider_projectile_pLab);
    const double p_lab_t = targ_cfg.take(InputKeys::modi_collider_target_pLab);
    if (p_lab_p < 0 || p_lab_t < 0) {
      throw ModusDefault::InvalidEnergy(
          "Input Error: "
          "P_Lab must be nonnegative.");
    }
    total_s_ = s_from_plab(p_lab_p * projectile_->number_of_particles(),
                           p_lab_t * target_->number_of_particles(),
                           mass_projec, mass_target);
    sqrt_s_NN_ = std::sqrt(s_from_plab(p_lab_p, p_lab_t, mass_a, mass_b));
    energy_input++;
  }
  if (energy_input == 0) {
    throw std::domain_error(
        "Input Error: Non-existent collision energy. "
        "Please provide one of Sqrtsnn/E_Kin/P_Lab.");
  }
  if (energy_input > 1) {
    throw std::invalid_argument(
        "Input Error: Redundant collision energy. "
        "Please provide only one of Sqrtsnn/E_Kin/P_Lab.");
  }

  /* Impact parameter setting: Either "Value", "Range", "Max" or "Sample".
   * Unspecified means 0 impact parameter.*/
  if (modus_cfg.has_value(InputKeys::modi_collider_impact_value)) {
    impact_ = modus_cfg.take(InputKeys::modi_collider_impact_value);
    imp_min_ = impact_;
    imp_max_ = impact_;
  } else {
    // If impact is not supplied by value, inspect sampling parameters:
    if (modus_cfg.has_value(InputKeys::modi_collider_impact_sample)) {
      sampling_ = modus_cfg.take(InputKeys::modi_collider_impact_sample);
      if (sampling_ == Sampling::Custom) {
        if (!(modus_cfg.has_value(InputKeys::modi_collider_impact_values) ||
              modus_cfg.has_value(InputKeys::modi_collider_impact_yields))) {
          throw std::invalid_argument(
              "Input Error: Need impact parameter spectrum for custom sampling."
              " Please provide Values and Yields.");
        }
        const std::vector<double> impacts =
            modus_cfg.take(InputKeys::modi_collider_impact_values);
        const std::vector<double> yields =
            modus_cfg.take(InputKeys::modi_collider_impact_yields);
        if (impacts.size() != yields.size()) {
          throw std::invalid_argument(
              "Input Error: Need as many impact parameter values as yields. "
              "Please make sure that Values and Yields have the same length.");
        }
        impact_interpolation_ = std::make_unique<InterpolateDataLinear<double>>(
            InterpolateDataLinear<double>(impacts, yields));

        const auto imp_minmax =
            std::minmax_element(impacts.begin(), impacts.end());
        imp_min_ = *imp_minmax.first;
        imp_max_ = *imp_minmax.second;
        yield_max_ = *std::max_element(yields.begin(), yields.end());
      }
    }
    if (modus_cfg.has_value(InputKeys::modi_collider_impact_range)) {
      const std::array<double, 2> range =
          modus_cfg.take(InputKeys::modi_collider_impact_range);
      imp_min_ = range[0];
      imp_max_ = range[1];
    }
    if (modus_cfg.has_value(InputKeys::modi_collider_impact_max)) {
      imp_min_ = 0.0;
      imp_max_ = modus_cfg.take(InputKeys::modi_collider_impact_max);
    }
  }
  /// \todo include a check that only one method of specifying impact is used
  // whether the direction of separation should be randomly sampled
  random_reaction_plane_ =
      modus_cfg.take(InputKeys::modi_collider_impact_randomReactionPlane);
  // Look for user-defined initial separation between nuclei.
  // The displacement is half the distance (both nuclei are shifted
  // initial_z_displacement_ away from origin)
  initial_z_displacement_ =
      modus_cfg.take(InputKeys::modi_collider_initialDistance) / 2.0;
  if (modus_cfg.has_section(InputSections::m_c_initialConditions)) {
    IC_for_hybrid_ = true;
    IC_parameters_ = std::make_unique<InitialConditionParameters>();
    IC_parameters_->type =
        modus_cfg.take(InputKeys::modi_collider_initialConditions_type);

    if (IC_parameters_->type == FluidizationType::ConstantTau) {
      FluidizationAction::remove_particle_ = true;
      if (modus_cfg.has_value(
              InputKeys::modi_collider_initialConditions_properTime)) {
        IC_parameters_->proper_time = modus_cfg.take(
            InputKeys::modi_collider_initialConditions_properTime);
      } else {
        IC_parameters_->lower_bound = modus_cfg.take(
            InputKeys::modi_collider_initialConditions_lowerBound);
        IC_parameters_->proper_time_scaling =
            modus_cfg.take(InputKeys::modi_collider_initialConditions_scaling);
      }
      IC_parameters_->rapidity_cut = modus_cfg.take(
          InputKeys::modi_collider_initialConditions_rapidityCut);
      IC_parameters_->pT_cut =
          modus_cfg.take(InputKeys::modi_collider_initialConditions_pTCut);
      validate_IC_kinematic_range();
    } else if (IC_parameters_->type == FluidizationType::Dynamic) {
      FluidizationAction::remove_particle_ = false;
      double threshold = modus_cfg.take(
          InputKeys::modi_collider_initialConditions_eDenThreshold);
      double min_time =
          modus_cfg.take(InputKeys::modi_collider_initialConditions_minTime);
      double max_time =
          modus_cfg.take(InputKeys::modi_collider_initialConditions_maxTime);
      int cells =
          modus_cfg.take(InputKeys::modi_collider_initialConditions_fluidCells);
      double form_time_fraction = modus_cfg.take(
          InputKeys::modi_collider_initialConditions_formTimeFraction);
      if (threshold <= 0 || max_time < min_time || min_time < 0 || cells < 2 ||
          form_time_fraction < 0) {
        logg[LInitialConditions].fatal()
            << "Bad parameters chosen for dynamic initial conditions. At least "
               "one of the following inequalities is violated:\n"
            << "  Energy_Density_Threshold = " << threshold << " > 0\n"
            << "  Maximum_Time = " << max_time << " > " << min_time
            << " = Minimum_Time > 0\n"
               "Fluidization_Cells = "
            << cells << " > 2\n"
            << " Formation_Time_Fraction < 0";
        throw std::invalid_argument("Please adjust the configuration file.");
      }

      IC_parameters_->fluidizable_processes = modus_cfg.take(
          InputKeys::modi_collider_initialConditions_fluidProcesses);

      double min_size = std::max(min_time, 40.);
      std::array<double, 3> length{2 * min_size, 2 * min_size, 2 * min_size};
      std::array<double, 3> origin{-min_size, -min_size, -min_size};
      std::array<int, 3> cell_array{cells, cells, cells};

      fluid_lattice_ =
          std::make_unique<RectangularLattice<EnergyMomentumTensor>>(
              length, cell_array, origin, false, LatticeUpdate::EveryTimestep);
      fluid_background_ = std::make_unique<std::map<int32_t, double>>();

      IC_parameters_->energy_density_threshold = threshold;
      IC_parameters_->min_time = min_time;
      IC_parameters_->max_time = max_time;
      IC_parameters_->num_fluid_cells = cells;
      logg[LInitialConditions].info()
          << "Preparing dynamic Initial Conditions with threshold " << threshold
          << " GeV/fm³ in energy density, between " << min_time << " and "
          << max_time << " fm.";
      IC_parameters_->formation_time_fraction = form_time_fraction;
      IC_parameters_->smearing_kernel_at_0 =
          std::pow(2 * M_PI * params.gaussian_sigma, -1.5);
      IC_parameters_->delay_initial_elastic = modus_cfg.take(
          InputKeys::modi_collider_initialConditions_delayInitialElastic);
    }
  }
}

void ColliderModus::validate_IC_kinematic_range() {
  bool bad_cuts = false;
  assert(IC_parameters_->rapidity_cut.has_value());
  assert(IC_parameters_->pT_cut.has_value());
  const double rapidity = IC_parameters_->rapidity_cut.value();
  const double pT = IC_parameters_->pT_cut.value();
  if (rapidity < 0.0) {
    logg[LInitialConditions].fatal()
        << "Rapidity cut for initial conditions configured as |y| < "
        << rapidity
        << " is unreasonable. \nPlease choose a positive, non-zero value or "
           "employ SMASH without rapidity cut.";
    bad_cuts = true;
  }
  if (pT < 0.0) {
    logg[LInitialConditions].fatal()
        << "Transverse momentum cut for initial conditions configured as pT < "
        << pT
        << " is unreasonable. \nPlease choose a positive, non-zero value or "
           "employ SMASH without pT cut.";
    bad_cuts = true;
  }
  if (bad_cuts) {
    throw std::runtime_error(
        "Kinematic cut for initial conditions malconfigured.");
  }

  std::ostringstream message{"Extracting iso-tau initial conditions ",
                             std::ios_base::ate};
  std::vector<std::string> cuts{};
  if (rapidity > 0.0) {
    cuts.emplace_back("|y| <= ");
    cuts.back() += std::to_string(rapidity);
  }
  if (pT > 0.0) {
    cuts.emplace_back("pT <= ");
    cuts.back() += std::to_string(pT) + " GeV.";
  }
  if (cuts.size() > 0) {
    message << "in kinematic range: " << join(cuts, "; ") << ".";
  } else {
    message << "without kinematic cuts.";
  }
  logg[LInitialConditions].info() << message.str();
}

std::ostream &operator<<(std::ostream &out, const ColliderModus &m) {
  return out << "-- Collider Modus:\n"
             << "sqrt(S) (nucleus-nucleus) = "
             << format(std::sqrt(m.total_s_), "GeV\n")
             << "sqrt(S) (nucleon-nucleon) = " << format(m.sqrt_s_NN_, "GeV\n")
             << "Projectile:\n"
             << *m.projectile_ << "\nTarget:\n"
             << *m.target_;
}

std::unique_ptr<DeformedNucleus> ColliderModus::create_deformed_nucleus(
    Configuration &nucleus_cfg, int ntest, const std::string &nucleus_type) {
  assert(has_projectile_or_target(nucleus_cfg));
  const bool is_projectile = is_about_projectile(nucleus_cfg);
  const auto &[automatic_key, beta2_key, beta3_key, beta4_key,
               gamma_key] = [&is_projectile]() {
    return is_projectile
               ? std::make_tuple(
                     InputKeys::modi_collider_projectile_deformed_automatic,
                     InputKeys::modi_collider_projectile_deformed_beta2,
                     InputKeys::modi_collider_projectile_deformed_beta3,
                     InputKeys::modi_collider_projectile_deformed_beta4,
                     InputKeys::modi_collider_projectile_deformed_gamma)
               : std::make_tuple(
                     InputKeys::modi_collider_target_deformed_automatic,
                     InputKeys::modi_collider_target_deformed_beta2,
                     InputKeys::modi_collider_target_deformed_beta3,
                     InputKeys::modi_collider_target_deformed_beta4,
                     InputKeys::modi_collider_target_deformed_gamma);
  }();

  bool automatic_deformation = nucleus_cfg.take(automatic_key);
  bool was_any_beta_given = nucleus_cfg.has_value(beta2_key) ||
                            nucleus_cfg.has_value(beta3_key) ||
                            nucleus_cfg.has_value(beta4_key);
  bool was_any_deformation_parameter_given =
      was_any_beta_given || nucleus_cfg.has_value(gamma_key);
  bool was_gamma_given_without_beta_2 =
      nucleus_cfg.has_value(gamma_key) && !nucleus_cfg.has_value(beta2_key);

  if (automatic_deformation && was_any_deformation_parameter_given) {
    throw std::invalid_argument(
        "Automatic deformation of " + nucleus_type +
        " nucleus requested, but deformation parameter(s) were provided as"
        " well. Please, check the 'Deformed' section in your input file.");
  } else if (!automatic_deformation && !was_any_beta_given) {
    throw std::invalid_argument(
        "Manual deformation of " + nucleus_type +
        " nucleus requested, but no deformation beta parameter was provided."
        " Please, check the 'Deformed' section in your input file.");
  } else if (!automatic_deformation && was_gamma_given_without_beta_2) {
    throw std::invalid_argument(
        "Manual deformation of " + nucleus_type +
        " nucleus requested, but 'Gamma' parameter was provided without "
        "providing a value of 'Beta_2' having hence no deformation effect. "
        "Please, check the 'Deformed' section in your input file.");
  } else {
    return std::make_unique<DeformedNucleus>(nucleus_cfg, ntest,
                                             automatic_deformation);
  }
}

std::unique_ptr<AlphaClusteredNucleus>
ColliderModus::create_alphaclustered_nucleus(Configuration &nucleus_cfg,
                                             int ntest,
                                             const std::string &nucleus_type) {
  const bool is_projectile = is_about_projectile(nucleus_cfg);
  const auto &[automatic_key, side_length_key] = [&is_projectile]() {
    return is_projectile
               ? std::make_tuple(
                     InputKeys::
                         modi_collider_projectile_alphaClustered_automatic,
                     InputKeys::
                         modi_collider_projectile_alphaClustered_sideLength)
               : std::make_tuple(
                     InputKeys::modi_collider_target_alphaClustered_automatic,
                     InputKeys::modi_collider_target_alphaClustered_sideLength);
  }();

  bool automatic_alphaclustering = nucleus_cfg.take(automatic_key);
  bool was_sidelength_given = nucleus_cfg.has_value(side_length_key);

  if (automatic_alphaclustering && was_sidelength_given) {
    throw std::invalid_argument(
        "Automatic alpha-clustering of " + nucleus_type +
        " nucleus requested, but a sidelength was provided as"
        " well. Please, check the 'Alpha_Clustered' section in your input "
        "file.");
  } else if (!automatic_alphaclustering && !was_sidelength_given) {
    throw std::invalid_argument(
        "Manual alpha-clustering of " + nucleus_type +
        " nucleus requested, but no sidelength was provided."
        " Please, check the 'Alpha_Clustered' section in your input file.");
  } else {
    return std::make_unique<AlphaClusteredNucleus>(nucleus_cfg, ntest,
                                                   automatic_alphaclustering);
  }
}

double ColliderModus::initial_conditions(Particles *particles,
                                         const ExperimentParameters &) {
  // Populate the nuclei with appropriately distributed nucleons.
  // If deformed, this includes rotating the nucleus.
  projectile_->arrange_nucleons();
  target_->arrange_nucleons();

  // Use the total mandelstam variable to get the frame-dependent velocity for
  // each nucleus. Position a is projectile, position b is target.
  double v_a, v_b;
  std::tie(v_a, v_b) =
      get_velocities(total_s_, projectile_->mass(), target_->mass());

  // If velocities are larger or equal to 1, throw an exception.
  if (v_a >= 1.0 || v_b >= 1.0) {
    throw std::domain_error(
        "Found velocity equal to or larger than 1 in "
        "ColliderModus::initial_conditions.\nConsider using "
        "the center of velocity reference frame.");
  }

  // Calculate the beam velocity of the projectile and the target, which will
  // be used to calculate the beam momenta in experiment.cc
  if (fermi_motion_ == FermiMotion::Frozen) {
    velocity_projectile_ = v_a;
    velocity_target_ = v_b;
  }

  // Generate Fermi momenta if necessary
  if (fermi_motion_ == FermiMotion::On ||
      fermi_motion_ == FermiMotion::Frozen) {
    // Frozen: Fermi momenta will be ignored during the propagation to
    // avoid that the nuclei will fly apart.
    projectile_->generate_fermi_momenta();
    target_->generate_fermi_momenta();
  } else if (fermi_motion_ == FermiMotion::Off) {
  } else {
    throw std::invalid_argument("Invalid Fermi_Motion input.");
  }

  // Boost the nuclei to the appropriate velocity.
  projectile_->boost(v_a);
  target_->boost(v_b);

  // Shift the nuclei into starting positions. Contracted spheres with
  // nuclear radii should touch exactly at t=0. Modus starts at negative
  // time corresponding to additional initial displacement.
  const double d_a = std::max(0., projectile_->get_diffusiveness());
  const double d_b = std::max(0., target_->get_diffusiveness());
  const double r_a = projectile_->get_nuclear_radius();
  const double r_b = target_->get_nuclear_radius();
  const double dz = initial_z_displacement_;

  const double simulation_time = -dz / std::abs(v_a);
  const double proj_z = -dz - std::sqrt(1.0 - v_a * v_a) * (r_a + d_a);
  const double targ_z =
      +dz * std::abs(v_b / v_a) + std::sqrt(1.0 - v_b * v_b) * (r_b + d_b);
  // rotation angle in the transverse plane
  const double phi =
      random_reaction_plane_ ? random::uniform(0.0, 2.0 * M_PI) : 0.0;

  projectile_->shift(proj_z, +impact_ / 2.0, simulation_time);
  target_->shift(targ_z, -impact_ / 2.0, simulation_time);

  // Put the particles in the nuclei into code particles.
  projectile_->copy_particles(particles);
  target_->copy_particles(particles);
  rotate_reaction_plane(phi, particles);
  return simulation_time;
}

void ColliderModus::rotate_reaction_plane(double phi, Particles *particles) {
  for (ParticleData &p : *particles) {
    ThreeVector pos = p.position().threevec();
    ThreeVector mom = p.momentum().threevec();
    pos.rotate_around_z(phi);
    mom.rotate_around_z(phi);
    p.set_3position(pos);
    p.set_3momentum(mom);
  }
}

void ColliderModus::sample_impact() {
  switch (sampling_) {
    case Sampling::Quadratic: {
      // quadratic sampling: Note that for bmin > bmax, this still yields
      // the correct distribution (however canonical() = 0 is then the
      // upper end, not the lower).
      impact_ = std::sqrt(imp_min_ * imp_min_ +
                          random::canonical() *
                              (imp_max_ * imp_max_ - imp_min_ * imp_min_));
    } break;
    case Sampling::Custom: {
      // rejection sampling based on given distribution
      assert(impact_interpolation_ != nullptr);
      double probability_random = 1;
      double probability = 0;
      double b;
      while (probability_random > probability) {
        b = random::uniform(imp_min_, imp_max_);
        probability = (*impact_interpolation_)(b) / yield_max_;
        assert(probability < 1.);
        probability_random = random::uniform(0., 1.);
      }
      impact_ = b;
    } break;
    case Sampling::Uniform: {
      // linear sampling. Still, min > max works fine.
      impact_ = random::uniform(imp_min_, imp_max_);
    }
  }
}

std::pair<double, double> ColliderModus::get_velocities(double s, double m_a,
                                                        double m_b) {
  double v_a = 0.0;
  double v_b = 0.0;
  // Frame dependent calculations of velocities. Assume v_a >= 0, v_b <= 0.
  switch (frame_) {
    case CalculationFrame::CenterOfVelocity:
      v_a = center_of_velocity_v(s, m_a, m_b);
      v_b = -v_a;
      break;
    case CalculationFrame::CenterOfMass: {
      // Compute center of mass momentum.
      double pCM = pCM_from_s(s, m_a, m_b);
      v_a = pCM / std::sqrt(m_a * m_a + pCM * pCM);
      v_b = -pCM / std::sqrt(m_b * m_b + pCM * pCM);
    } break;
    case CalculationFrame::FixedTarget:
      v_a = fixed_target_projectile_v(s, m_a, m_b);
      break;
    default:
      throw std::invalid_argument(
          "Invalid reference frame in "
          "ColliderModus::get_velocities.");
  }
  return std::make_pair(v_a, v_b);
}

std::string ColliderModus::custom_file_path(const std::string &file_directory,
                                            const std::string &file_name) {
  // make sure that path is correct even if the / at the end is missing
  if (file_directory.back() == '/') {
    return file_directory + file_name;
  } else {
    return file_directory + '/' + file_name;
  }
}

bool ColliderModus::same_inputfile(Configuration &proj_config,
                                   Configuration &targ_config) {
  /* Check if both nuclei are custom
   * Only check target as function is called after if statement for
   * projectile.
   */
  if (!targ_config.has_section(InputSections::m_c_t_custom)) {
    return false;
  }
  std::string projectile_file_directory = proj_config.read(
      InputKeys::modi_collider_projectile_custom_fileDirectory);
  std::string target_file_directory =
      targ_config.read(InputKeys::modi_collider_target_custom_fileDirectory);
  std::string projectile_file_name =
      proj_config.read(InputKeys::modi_collider_projectile_custom_fileName);
  std::string target_file_name =
      targ_config.read(InputKeys::modi_collider_target_custom_fileName);
  // Check if files are the same for projectile and target
  std::string proj_path =
      custom_file_path(projectile_file_directory, projectile_file_name);
  std::string targ_path =
      custom_file_path(target_file_directory, target_file_name);
  if (proj_path == targ_path) {
    return true;
  } else {
    return false;
  }
}

void ColliderModus::build_fluidization_lattice(
    const double t, const std::vector<Particles> &ensembles,
    const DensityParameters &dens_par) {
  if (fluid_lattice_ == nullptr) {
    throw std::logic_error(
        "Trying to build fluidization lattice with unset pointer in "
        "ColliderModus.");
  }
  if (t < IC_parameters_->min_time.value() ||
      t > IC_parameters_->max_time.value()) {
    return;
  }
  const double resizing_rate = 5;
  double side = fluid_lattice_->lattice_sizes()[0] / 2.;
  if (t > side) {
    side += resizing_rate;
    std::array<double, 3> new_length{2 * side, 2 * side, 2 * side};
    std::array<double, 3> new_origin{-side, -side, -side};
    fluid_lattice_->reset_and_resize(new_length, new_origin, std::nullopt);
    logg[LCollider].debug() << "Fluidization lattice resizing at " << t
                            << " fm to " << 2 * side << " fm";
  }

  update_lattice_accumulating_ensembles(
      fluid_lattice_.get(), LatticeUpdate::EveryTimestep, DensityType::Hadron,
      dens_par, ensembles, false);
}

}  // namespace smash
