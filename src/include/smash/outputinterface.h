/*
 *
 *    Copyright (c) 2014-2022
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
 * \brief Abstraction of generic output
 *
 * Any output should inherit this class. It provides virtual methods that will
 * be called at predefined moments:
 * 1) At event start and event end: at_eventstart, at_eventend
 * 2) After every fixed time period: at_intermediate_time, thermodynamics_output
 * 3) At each interaction: at_interaction
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
        is_IC_output_(name == "SMASH_IC") {}
  virtual ~OutputInterface() = default;

  /**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   * \param particles List of particles.
   * \param event_number Number of the current event.
   * \param[in] info Event info, see \ref event_info
   */
  virtual void at_eventstart(const Particles &particles, const int event_number,
                             const EventInfo &info) {
    SMASH_UNUSED(particles);
    SMASH_UNUSED(event_number);
    SMASH_UNUSED(info);
  }
  /**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   * \param ensembles List of particles.
   * \param[in] event_number Number of the current event.
   */
  virtual void at_eventstart(const std::vector<Particles> &ensembles,
                             int event_number) {
    SMASH_UNUSED(ensembles);
    SMASH_UNUSED(event_number);
  }

  /**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   * \param[in] event_number Number of the current event.
   * \param[in] tq Thermodynamic quantity to deal with.
   * \param[in] dens_type Density type for the reference frame.
   * \param[in] lattice Lattice of tabulated values.
   */
  virtual void at_eventstart(const int event_number,
                             const ThermodynamicQuantity tq,
                             const DensityType dens_type,
                             RectangularLattice<DensityOnLattice> lattice) {
    SMASH_UNUSED(event_number);
    SMASH_UNUSED(tq);
    SMASH_UNUSED(dens_type);
    SMASH_UNUSED(lattice);
  }

  /**
   * Output launched atevent start after initialization, when particles are
   * generated but not yet propagated.
   * \param[in] event_number Number of the current event.
   * \param[in] tq Thermodynamic quantity to deal with.
   * \param[in] dens_type Density type for the reference frame.
   * \param[in] lattice Lattice of tabulated values.
   */
  virtual void at_eventstart(const int event_number,
                             const ThermodynamicQuantity tq,
                             const DensityType dens_type,
                             RectangularLattice<EnergyMomentumTensor> lattice) {
    SMASH_UNUSED(event_number);
    SMASH_UNUSED(tq);
    SMASH_UNUSED(dens_type);
    SMASH_UNUSED(lattice);
  }

  /**
   * Output launched at event end. Event end is determined by maximal timestep
   * option.
   * \param[in] event_number Number of the current event.
   * \param[in] tq Thermodynamic quantity to deal with
   * \param[in] dens_type Density type for the evaluation of thermodynamic
   *                      quantities
   */
  virtual void at_eventend(const int event_number,
                           const ThermodynamicQuantity tq,
                           const DensityType dens_type) {
    SMASH_UNUSED(event_number);
    SMASH_UNUSED(tq);
    SMASH_UNUSED(dens_type);
  }

  /**
   * Output launched at event end. Event end is determined by maximal timestep
   * option.
   * \param[in] tq Thermodynamic quantity to deal with.
   */
  virtual void at_eventend(const ThermodynamicQuantity tq) { SMASH_UNUSED(tq); }

  /**
   * Output launched at event end. Event end is determined by maximal timestep
   * option.
   * \param particles List of particles.
   * \param event_number Number of the current event.
   * \param[in] info Event info, see \ref event_info
   */
  virtual void at_eventend(const Particles &particles, const int event_number,
                           const EventInfo &info) {
    SMASH_UNUSED(particles);
    SMASH_UNUSED(event_number);
    SMASH_UNUSED(info);
  }
  /**
   * Output launched at event end. Event end is determined by maximal timestep
   * option.
   * \param ensembles List of particles.
   * \param event_number Number of the current event.
   */
  virtual void at_eventend(const std::vector<Particles> &ensembles,
                           const int event_number) {
    SMASH_UNUSED(ensembles);
    SMASH_UNUSED(event_number);
  }

  /**
   * Called whenever an action modified one or more particles.
   *
   * \param action The action object, containing the initial and final state
   * etc.
   * \param density The density at the interaction point.
   */
  virtual void at_interaction(const Action &action, const double density) {
    SMASH_UNUSED(action);
    SMASH_UNUSED(density);
  }

  /**
   * Output launched after every N'th timestep. N is controlled by an option.
   * \param particles List of particles.
   * \param clock System clock.
   * \param dens_param Parameters for density calculation.
   * \param[in] info Event info, see \ref event_info
   */
  virtual void at_intermediate_time(const Particles &particles,
                                    const std::unique_ptr<Clock> &clock,
                                    const DensityParameters &dens_param,
                                    const EventInfo &info) {
    SMASH_UNUSED(particles);
    SMASH_UNUSED(clock);
    SMASH_UNUSED(dens_param);
    SMASH_UNUSED(info);
  }
  /**
   * Output launched after every N'th timestep. N is controlled by an option.
   * \param ensembles List of particles.
   * \param clock System clock.
   * \param dens_param Parameters for density calculation.
   */
  virtual void at_intermediate_time(const std::vector<Particles> &ensembles,
                                    const std::unique_ptr<Clock> &clock,
                                    const DensityParameters &dens_param) {
    SMASH_UNUSED(ensembles);
    SMASH_UNUSED(clock);
    SMASH_UNUSED(dens_param);
  }

  /**
   * Output to write thermodynamics from the lattice.
   * \param tq Thermodynamic quantity to be written, used for file name etc.
   * \param dt Type of density, i.e. which particles to take into account.
   * \param lattice Lattice of tabulated values.
   *
   * Used for vtk output.
   */
  virtual void thermodynamics_output(
      const ThermodynamicQuantity tq, const DensityType dt,
      RectangularLattice<DensityOnLattice> &lattice) {
    SMASH_UNUSED(tq);
    SMASH_UNUSED(dt);
    SMASH_UNUSED(lattice);
  }

  /**
   * Output to write energy-momentum tensor and related quantities from the
   * lattice.
   * \param tq Thermodynamic quantity to be written: Tmn, Tmn_Landau, v_Landau.
   * \param dt Type of density, i.e. which particles to take into account.
   * \param lattice Lattice of tabulated values.
   *
   * Used for vtk output.
   */
  virtual void thermodynamics_output(
      const ThermodynamicQuantity tq, const DensityType dt,
      RectangularLattice<EnergyMomentumTensor> &lattice) {
    SMASH_UNUSED(tq);
    SMASH_UNUSED(dt);
    SMASH_UNUSED(lattice);
  }

  /**
   * Output to write thermodynamics from the lattice.
   * \param[in] lattice Lattice type DensityOnLattice of tabulated values.
   * \param[in] current_time Time of the simulation in the computational frame.
   *
   * Used for thermodynamic lattice output.
   */
  virtual void thermodynamics_lattice_output(
      RectangularLattice<DensityOnLattice> &lattice,
      const double current_time) {
    SMASH_UNUSED(lattice);
    SMASH_UNUSED(current_time);
  }

  /**
   * Output to write thermodynamics from the lattice.
   * \param[in] lattice Lattice type FourVector of tabulated values.
   * \param[in] current_time Time of the simulation in the computational frame.
   * \param[in] ensembles Particles, from which the 4-currents j_{Q,B,S} are
   *            computed
   * * \param[in] dens_param set of parameters, defining smearing.
   *            For more info about
   *            smearing see \ref doxypage_output_thermodyn.
   *
   * Used for thermodynamic lattice output.
   */
  virtual void thermodynamics_lattice_output(
      RectangularLattice<DensityOnLattice> &lattice, const double current_time,
      const std::vector<Particles> &ensembles,
      const DensityParameters &dens_param) {
    SMASH_UNUSED(lattice);
    SMASH_UNUSED(current_time);
    SMASH_UNUSED(ensembles);
    SMASH_UNUSED(dens_param);
  }

  /**
   * Output to write energy-momentum tensor and related quantities from the
   * lattice.
   * \param[in] tq Thermodynamic quantity to be written, used for file name etc.
   * \param[in] lattice Lattice type EnergyMomentumTensor of tabulated values.
   * \param[in] current_time Time of the simulation in the computational frame.
   *
   * Used for thermodynamic lattice output.
   */
  virtual void thermodynamics_lattice_output(
      const ThermodynamicQuantity tq,
      RectangularLattice<EnergyMomentumTensor> &lattice,
      const double current_time) {
    SMASH_UNUSED(tq);
    SMASH_UNUSED(lattice);
    SMASH_UNUSED(current_time);
  }

  /**
   * Output to write energy-momentum tensor and related quantities from the
   * thermalizer class.
   * \param gct Pointer to thermalizer
   *
   * Only used for vtk output. Not connected to ThermodynamicOutput.
   */
  virtual void thermodynamics_output(const GrandCanThermalizer &gct) {
    SMASH_UNUSED(gct);
  }

  /**
   * Write fields in vtk output
   *
   * Fields are a pair of threevectors for example electric and magnetic field
   *
   * \param[in] name1 Name of the first field
   * \param[in] name2 Name of the second field
   * \param[in] lat Lattice storing both fields
   */
  virtual void fields_output(
      const std::string name1, const std::string name2,
      RectangularLattice<std::pair<ThreeVector, ThreeVector>> &lat) {
    SMASH_UNUSED(name1);
    SMASH_UNUSED(name2);
    SMASH_UNUSED(lat);
  }

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

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_OUTPUTINTERFACE_H_
