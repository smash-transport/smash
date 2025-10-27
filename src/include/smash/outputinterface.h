/*
 *
 *    Copyright (c) 2014-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_OUTPUTINTERFACE_H_
#define SRC_INCLUDE_SMASH_OUTPUTINTERFACE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "density.h"
#include "energymomentumtensor.h"
#include "forwarddeclarations.h"
#include "grandcan_thermalizer.h"
#include "lattice.h"
#include "macros.h"

namespace smash {
static constexpr int LOutput = LogArea::Output::id;

/**
 * \ingroup output
 *
 * \brief Structure to contain custom data for output
 *
 * \anchor event_info
 * This structure is intended to hold and conveniently pass information about
 * event such as impact parameter, total potential energy, and similar
 * auxiliary info.
 */
struct EventInfo {
  /// Impact parameter for collider modus, otherwise dummy
  double impact_parameter;
  /// Box length in case of box simulation, otherwise dummy
  double modus_length;
  /// Time in fm
  double current_time;
  /// Sum of kinetic energies of all particles
  double total_kinetic_energy;
  /// Total energy in the mean field
  double total_mean_field_energy;
  /// Kinetic + mean field energy
  double total_energy;
  /// Testparticle number, see Testparticles in \ref doxypage_input_conf_general
  int test_particles;
  /// Number of ensembles
  int n_ensembles;
  /// True if no collisions happened
  bool empty_event;
  /// Whether or not kinematic cuts are employed for SMASH IC
  bool impose_kinematic_cut_for_SMASH_IC;
};

/**
 * \ingroup output
 *
 * \brief Structure to contain information about the event and ensemble numbers
 *
 * \anchor event_label
 * This structure is intended to hold and conveniently pass information to
 * create a unique label for each block of output about a given ensemble of a
 * given event.
 */
struct EventLabel {
  /// The number of the event
  int32_t event_number;
  /// The number of the ensemble
  int32_t ensemble_number;
};

/**
 * \ingroup output
 *
 * \brief Abstraction of generic output
 *
 * Any output should inherit this class. It provides virtual methods that will
 * be called at predefined moments:
 * 1) At event start and event end: at_eventstart, at_eventend
 * 2) After every fixed time period: at_intermediate_time, thermodynamics_output
 * 3) At each interaction: at_interaction
 *
 * \attention This class provides more virtual methods than those needed in
 * different children classes. Although this is against the inheritance "is-a"
 * relationship, it somehow simplifies here the hierarchy, because we avoid
 * having many more interfaces. Furthermore all base virtual methods are
 * implemented as no-operations, i.e. empty. This implies that they can and will
 * be called if the interface is used polymorphically but the child does not
 * implement the called method. This happens e.g. in the Experiment class where
 * an array of pointers to the base class is initialized with different children
 * and then different methods are called on all array entries (some will do what
 * has to be done, but most will just do nothing).
 *
 * \note The parameters of most methods in this base class are not documented,
 * as irrelevant for the empty implementation. However, every child class which
 * overrides some methods documents them in detail. Refer to them for further
 * information.
 *
 * \attention Some methods take an \c EventLabel parameter and some just an
 * integer. The former contains both the number of event and the number of
 * ensemble while the latter is meant to be the number of the event, only.
 * Passing the number of event only is meant to emphasise that the method will
 * combine input from different ensembles, when multiple are used.
 */
class OutputInterface {
 public:
  /**
   * Construct output interface.
   * \param[in] name (File)name of output.
   */
  explicit OutputInterface(std::string name)
      : is_dilepton_output_(name == "Dileptons"),
        is_photon_output_(name == "Photons"),
        is_IC_output_(name.substr(0, 8) == "SMASH_IC") {}
  /**
   * Pure virtual destructor to make class abstract and prevent its
   * instantiation. It needs a definition which is done outside the class.
   */
  virtual ~OutputInterface() = 0;

  /**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   */
  virtual void at_eventstart(const Particles &, const EventLabel &,
                             const EventInfo &) {}

  /**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   */
  virtual void at_eventstart(const std::vector<Particles> &, int) {}

  /**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   */
  virtual void at_eventstart(const int, const ThermodynamicQuantity,
                             const DensityType,
                             RectangularLattice<DensityOnLattice>) {}

  /**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   */
  virtual void at_eventstart(const int, const ThermodynamicQuantity,
                             const DensityType,
                             RectangularLattice<EnergyMomentumTensor>) {}

  /**
   * Output launched at event end. Event end is determined by maximal time-step
   * option.
   */
  virtual void at_eventend(const ThermodynamicQuantity) {}

  /**
   * Output launched at event end. Event end is determined by maximal time-step
   * option.
   */
  virtual void at_eventend(const Particles &, const EventLabel &,
                           const EventInfo &) {}
  /**
   * Output launched at event end. Event end is determined by maximal time-step
   * option.
   */
  virtual void at_eventend(const std::vector<Particles> &, const int) {}

  /**
   * Called whenever an action modified one or more particles.
   */
  virtual void at_interaction(const Action &, const double) {}

  /**
   * Output launched after every N'th time-step. N is controlled by an option.
   */
  virtual void at_intermediate_time(const Particles &,
                                    const std::unique_ptr<Clock> &,
                                    const DensityParameters &,
                                    const EventLabel &, const EventInfo &) {}
  /**
   * Output launched after every N'th timestep. N is controlled by an option.
   */
  virtual void at_intermediate_time(const std::vector<Particles> &,
                                    const std::unique_ptr<Clock> &,
                                    const DensityParameters &) {}

  /**
   * Output to write thermodynamics from the lattice.
   * Used for vtk output.
   */
  virtual void thermodynamics_output(const ThermodynamicQuantity,
                                     const DensityType,
                                     RectangularLattice<DensityOnLattice> &) {}

  /**
   * Output to write energy-momentum tensor and related quantities from the
   * lattice. Used for vtk output.
   */
  virtual void thermodynamics_output(
      const ThermodynamicQuantity, const DensityType,
      RectangularLattice<EnergyMomentumTensor> &) {}

  /**
   * Output to write thermodynamics from the lattice.
   * Used for thermodynamic lattice output.
   */
  virtual void thermodynamics_lattice_output(
      RectangularLattice<DensityOnLattice> &, const double) {}

  /**
   * Output to write thermodynamics from the lattice.
   * Used for thermodynamic lattice output.
   */
  virtual void thermodynamics_lattice_output(
      RectangularLattice<DensityOnLattice> &, const double,
      const std::vector<Particles> &, const DensityParameters &) {}

  /**
   * Output to write energy-momentum tensor and related quantities from the
   * lattice. Used for thermodynamic lattice output.
   */
  virtual void thermodynamics_lattice_output(
      const ThermodynamicQuantity, RectangularLattice<EnergyMomentumTensor> &,
      const double) {}

  /**
   * Output to write energy-momentum tensor and related quantities from the
   * thermalizer class.
   * Only used for vtk output. Not connected to ThermodynamicOutput.
   */
  virtual void thermodynamics_output(const GrandCanThermalizer &) {}

  /**
   * Write fields in vtk output
   * Fields are a pair of threevectors for example electric and magnetic field
   */
  virtual void fields_output(
      const std::string, const std::string,
      RectangularLattice<std::pair<ThreeVector, ThreeVector>> &) {}

  /// Get, whether this is the dilepton output?
  bool is_dilepton_output() const { return is_dilepton_output_; }

  /// Get, whether this is the photon output?
  bool is_photon_output() const { return is_photon_output_; }

  /// Get, whether this is the IC output?
  bool is_IC_output() const { return is_IC_output_; }

  /**
   * Convert thermodynamic quantities to strings.
   * \param[in] tq Enum value of the thermodynamic quantity.
   * \return String description of the enumerator.
   */
  const char *to_string(const ThermodynamicQuantity tq) {
    switch (tq) {
      case ThermodynamicQuantity::EckartDensity:
        return "rho_eckart";
      case ThermodynamicQuantity::Tmn:
        return "tmn";
      case ThermodynamicQuantity::TmnLandau:
        return "tmn_landau";
      case ThermodynamicQuantity::LandauVelocity:
        return "v_landau";
      case ThermodynamicQuantity::j_QBS:
        return "j_QBS";
    }
    throw std::invalid_argument("Unknown thermodynamic quantity.");
  }

  /**
   * Convert density types to strings.
   * \param[in] dens_type enum value of the density type
   * \return String description of the enumerator.
   */
  const char *to_string(const DensityType dens_type) {
    switch (dens_type) {
      case DensityType::Hadron:
        return "hadron";
      case DensityType::Baryon:
        return "net_baryon";
      case DensityType::BaryonicIsospin:
        return "net_baryonI3";
      case DensityType::Pion:
        return "pion";
      case DensityType::Isospin3_tot:
        return "tot_isospin3";
      case DensityType::Charge:
        return "charge";
      case DensityType::Strangeness:
        return "strangeness";
      case DensityType::None:
        return "none";
    }
    throw std::invalid_argument("Unknown density type.");
  }

 protected:
  /// Is this the dilepton output?
  const bool is_dilepton_output_;

  /// Is this the photon output?
  const bool is_photon_output_;

  /// Is this the IC output?
  const bool is_IC_output_;
};

inline OutputInterface::~OutputInterface() = default;

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_OUTPUTINTERFACE_H_
