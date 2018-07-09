/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_OUTPUTINTERFACE_H_
#define SRC_INCLUDE_OUTPUTINTERFACE_H_

#include <string>

#include "density.h"
#include "energymomentumtensor.h"
#include "forwarddeclarations.h"
#include "grandcan_thermalizer.h"
#include "lattice.h"
#include "macros.h"

namespace smash {

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
        is_photon_output_(name == "Photons") {}
  virtual ~OutputInterface() = default;

  /**
   * Output launched at event start after initialization, when particles are
   * generated but not yet propagated.
   * \param particles List of particles.
   * \param event_number Number of the current event.
   */
  virtual void at_eventstart(const Particles &particles,
                             const int event_number) = 0;

  /**
   * Output launched at event end. Event end is determined by maximal timestep
   * option.
   * \param particles List of particles.
   * \param event_number Number of the current event.
   * \param impact_parameter Distance between centers of nuclei in this event.
   *          Only makes sense for collider modus.
   */
  virtual void at_eventend(const Particles &particles, const int event_number,
                           double impact_parameter) = 0;

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
   */
  virtual void at_intermediate_time(const Particles &particles,
                                    const Clock &clock,
                                    const DensityParameters &dens_param) {
    SMASH_UNUSED(particles);
    SMASH_UNUSED(clock);
    SMASH_UNUSED(dens_param);
  }

  /**
   * Output to write thermodynamics from the lattice.
   * \param tq Thermodynamic quantity to be written, used for file name etc.
   * \param dt Type of density, i.e. which particles to take into account.
   * \param lattice Lattice of tabulated values.
   *
   * Only used for vtk output. Not connected to ThermodynamicOutput.
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
   * \param tq Thermodynamic quantity to be written: Tmn, Tmn_Landau, v_Landau
   * \param dt Type of density, i.e. which particles to take into account.
   * \param lattice Lattice of tabulated values.
   *
   * Only used for vtk output. Not connected to ThermodynamicOutput.
   */
  virtual void thermodynamics_output(
      const ThermodynamicQuantity tq, const DensityType dt,
      RectangularLattice<EnergyMomentumTensor> &lattice) {
    SMASH_UNUSED(tq);
    SMASH_UNUSED(dt);
    SMASH_UNUSED(lattice);
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

  /// Get, whether this is the dilepton output?
  bool is_dilepton_output() const { return is_dilepton_output_; }

  /// Get, whether this is the photon output?
  bool is_photon_output() const { return is_photon_output_; }

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
};

}  // namespace smash

#endif  // SRC_INCLUDE_OUTPUTINTERFACE_H_
