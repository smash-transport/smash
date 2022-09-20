/*
 *    Copyright (c) 2017-2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_SMASH_OUTPUTPARAMETERS_H_
#define SRC_INCLUDE_SMASH_OUTPUTPARAMETERS_H_

#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "configuration.h"
#include "density.h"
#include "forwarddeclarations.h"
#include "logging.h"

namespace smash {
static constexpr int LExperiment = LogArea::Experiment::id;

/**
 * Helper structure for OutputParameters in order to store and hand over Rivet
 * parameters. OutputParameters has one member of this type.
 */
struct RivetOutputParameters {
  /// Logging in Rivet
  std::optional<std::map<std::string, std::string>> logs{std::nullopt};
  /// Paths to analyses libraries and data
  std::optional<std::vector<std::string>> paths{std::nullopt};
  /// Data files to pre-load e.g., for centrality configurations
  std::optional<std::vector<std::string>> preloads{std::nullopt};
  /// Analyses (including options) to add to run
  std::optional<std::vector<std::string>> analyses{std::nullopt};
  /// Weights to be enabled for processing
  std::optional<std::vector<std::string>> to_be_enabled_weights{std::nullopt};
  /// Weights to be disabled for processing
  std::optional<std::vector<std::string>> to_be_disabled_weights{std::nullopt};
  /// Cross sections
  std::optional<std::array<double, 2>> cross_sections{std::nullopt};
  /// Nominal weight name
  std::optional<std::string> nominal_weight_name{std::nullopt};
  /// Cap (maximum) on weights
  std::optional<double> cap_on_weights{std::nullopt};
  /// How to smear for NLO calculations
  std::optional<double> nlo_smearing{std::nullopt};
  /// Whether Rivet should not care about multi weights
  std::optional<bool> no_multi_weight{std::nullopt};
  /// Whether Rivet should ignore beams
  bool ignore_beams{true};
  /// Whether any weight parameter was specified
  bool any_weight_parameter_was_given{false};
};

/**
 * Helper structure for Experiment to hold output options and parameters.
 * Experiment has one member of this struct.
 */
struct OutputParameters {
  /// Default constructor, useful for tests
  OutputParameters()
      : td_position(ThreeVector()),
        td_dens_type(DensityType::None),
        td_rho_eckart(false),
        td_tmn(false),
        td_tmn_landau(false),
        td_v_landau(false),
        td_jQBS(false),
        td_smearing(true),
        td_only_participants(false),
        part_extended(false),
        part_only_final(OutputOnlyFinal::Yes),
        coll_extended(false),
        coll_printstartend(false),
        dil_extended(false),
        photons_extended(false),
        ic_extended(false),
        rivet_parameters{} {}

  /// Constructor from configuration
  explicit OutputParameters(Configuration conf) : OutputParameters() {
    logg[LExperiment].trace(SMASH_SOURCE_LOCATION);

    if (conf.has_value({"Thermodynamics"})) {
      auto subcon = conf["Thermodynamics"];
      if (subcon.has_value({"Position"})) {
        const std::array<double, 3> a = subcon.take({"Position"});
        td_position = ThreeVector(a[0], a[1], a[2]);
      }
      std::set<ThermodynamicQuantity> quan = subcon.take({"Quantities"});
      td_rho_eckart = (quan.count(ThermodynamicQuantity::EckartDensity) > 0);
      td_tmn = (quan.count(ThermodynamicQuantity::Tmn) > 0);
      td_tmn_landau = (quan.count(ThermodynamicQuantity::TmnLandau) > 0);
      td_v_landau = (quan.count(ThermodynamicQuantity::LandauVelocity) > 0);
      td_jQBS = (quan.count(ThermodynamicQuantity::j_QBS) > 0);
      td_dens_type = subcon.take({"Type"}, DensityType::Baryon);
      if (td_dens_type == DensityType::None &&
          (td_rho_eckart || td_tmn || td_tmn_landau || td_v_landau)) {
        logg[LExperiment].warn(
            "Requested Thermodynamics output with Density type None. ",
            "Change the density type to avoid output being dropped.");
      }
      td_smearing = subcon.take({"Smearing"}, true);
      td_only_participants = subcon.take({"Only_Participants"}, false);
    }

    if (conf.has_value({"Particles"})) {
      part_extended = conf.take({"Particles", "Extended"}, false);
      part_only_final =
          conf.take({"Particles", "Only_Final"}, OutputOnlyFinal::Yes);
    }

    if (conf.has_value({"Collisions"})) {
      coll_extended = conf.take({"Collisions", "Extended"}, false);
      coll_printstartend = conf.take({"Collisions", "Print_Start_End"}, false);
    }

    if (conf.has_value({"Dileptons"})) {
      dil_extended = conf.take({"Dileptons", "Extended"}, false);
    }

    if (conf.has_value({"Photons"})) {
      photons_extended = conf.take({"Photons", "Extended"}, false);
    }

    if (conf.has_value({"Initial_Conditions"})) {
      ic_extended = conf.take({"Initial_Conditions", "Extended"}, false);
    }

    if (conf.has_value({"Rivet"})) {
      logg[LOutput].debug() << "Reading Rivet section from configuration:\n"
                            << conf["Rivet"].to_string() << "\n";
      /*
       * std::optional<T> can be assigned from a value using the
       *   template<class U = T> optional& operator=( U&& value );
       * which is a perfect-forwarded assignment. However, in more complex cases
       * like here where we use a Configuration::Value object returned by
       * Configuration::take as value, it might be sometimes needed to
       * explicitly specify the type U. It is then advantageous to use
       * std::make_optional. When T is a built-in type, an assignment
       * would work (although not that mentioned above), but for simplicity
       * we always use make_optional.
       */
      if (conf.has_value({"Rivet", "Logging"})) {
        rivet_parameters.logs =
            std::make_optional<std::map<std::string, std::string>>(
                conf.take({"Rivet", "Logging"}));
      }
      if (conf.has_value({"Rivet", "Paths"})) {
        rivet_parameters.paths = std::make_optional<std::vector<std::string>>(
            conf.take({"Rivet", "Paths"}));
      }
      if (conf.has_value({"Rivet", "Preloads"})) {
        rivet_parameters.preloads =
            std::make_optional<std::vector<std::string>>(
                conf.take({"Rivet", "Preloads"}));
      }
      if (conf.has_value({"Rivet", "Analyses"})) {
        rivet_parameters.analyses =
            std::make_optional<std::vector<std::string>>(
                conf.take({"Rivet", "Analyses"}));
      }
      if (conf.has_value({"Rivet", "Cross_Section"})) {
        rivet_parameters.cross_sections =
            std::make_optional<std::array<double, 2>>(
                conf.take({"Rivet", "Cross_Section"}));
      }
      rivet_parameters.ignore_beams =
          conf.take({"Rivet", "Ignore_Beams"}, true);
      if (conf.has_value({"Weights"})) {
        rivet_parameters.any_weight_parameter_was_given = true;
        if (conf.has_value({"Rivet", "Weights", "Select"})) {
          rivet_parameters.to_be_enabled_weights =
              std::make_optional<std::vector<std::string>>(
                  conf.take({"Rivet", "Weights", "Select"}));
        }
        if (conf.has_value({"Rivet", "Weights", "Deselect"})) {
          rivet_parameters.to_be_disabled_weights =
              std::make_optional<std::vector<std::string>>(
                  conf.take({"Rivet", "Weights", "Deselect"}));
        }
        if (conf.has_value({"Rivet", "Weights", "Nominal"})) {
          rivet_parameters.nominal_weight_name =
              std::make_optional<std::string>(
                  conf.take({"Rivet", "Weights", "Nominal"}));
        }
        if (conf.has_value({"Rivet", "Weights", "Cap"})) {
          rivet_parameters.cap_on_weights = std::make_optional<double>(
              conf.take({"Rivet", "Weights", "Cap"}));
        }
        if (conf.has_value({"Rivet", "Weights", "NLO_Smearing"})) {
          rivet_parameters.nlo_smearing = std::make_optional<double>(
              conf.take({"Rivet", "Weights", "NLO_Smearing"}));
        }
        if (conf.has_value({"Rivet", "Weights", "No_Multi"})) {
          rivet_parameters.no_multi_weight = std::make_optional<bool>(
              conf.take({"Rivet", "Weights", "No_Multi"}));
        }
      }
      logg[LOutput].debug() << "After processing configuration:\n"
                            << conf["Rivet"].to_string() << "\n";
    }
  }

  /**
   * Pass correct extended flag to binary collision output constructor
   * \param[in] name (File)name of the output.
   * \return Extended flag for binary output.
   */
  bool get_coll_extended(std::string name) const {
    if (name == "Collisions") {
      return coll_extended;
    } else if (name == "Dileptons") {
      return dil_extended;
    } else if (name == "Photons") {
      return photons_extended;
    } else {
      return false;  // error
    }
  }

  /// Point, where thermodynamic quantities are calculated
  ThreeVector td_position;

  /// Type (e.g., baryon/pion/hadron) of thermodynamic quantity
  DensityType td_dens_type;

  /// Print out Eckart rest frame density of type td_dens_type or not?
  bool td_rho_eckart;

  /// Print out energy-momentum tensor of type td_dens_type or not?
  bool td_tmn;

  /**
   * Print out energy-momentum tensor in Landau rest frame
   * (of type td_dens_type) or not?
   */
  bool td_tmn_landau;

  /// Print out Landau velocity of type td_dens_type or not?
  bool td_v_landau;

  /// Print out QBS 4-currents or not?
  bool td_jQBS;

  /**
   * Whether smearing is on or off; WARNING : if smearing is off,
   * then final result is in GeV instead of GeV/fm3
   */
  bool td_smearing;

  /**
   * Flag reporting whether only participants are considered (true) or also
   * spectators (false)
   */
  bool td_only_participants;

  /// Extended format for particles output
  bool part_extended;

  /// Print only final particles in event
  OutputOnlyFinal part_only_final;

  /// Extended format for collisions output
  bool coll_extended;

  /// Print initial and final particles in event into collision output
  bool coll_printstartend;

  /// Extended format for dilepton output
  bool dil_extended;

  /// Extended format for photon output
  bool photons_extended;

  /// Extended initial conditions output
  bool ic_extended;

  /// Rivet specfic parameters
  RivetOutputParameters rivet_parameters;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_OUTPUTPARAMETERS_H_
