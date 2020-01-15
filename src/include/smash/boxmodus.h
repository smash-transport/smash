/*
 *    Copyright (c) 2012-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_BOXMODUS_H_
#define SRC_INCLUDE_BOXMODUS_H_

#include <map>
#include <memory>

#include "forwarddeclarations.h"
#include "modusdefault.h"

namespace smash {

/**
 * \ingroup modus
 * BoxModus: Provides a modus for infinite matter calculations
 *
 * Matter is confined in a cubical box. Depending on the initial
 * condition, particles are either reflected on the boundaries
 * (not implemented now) or inserted on opposite positions.
 *
 * To use this modus, choose
    Modus:         Box
 * \code
 * General:
 *      Modus: Box
 * \endcode
 * in the configuration file.
 *
 * Options for BoxModus go in the "Modi"→"Box" section of the
 * configuration:
 *
 * \code
 * Modi:
 *      Box:
 *              # definitions here
 * \endcode
 *
 * The following configuration options are understood:
 * \ref input_modi_box_
 */
class BoxModus : public ModusDefault {
 public:
  /**
   * Constructor
   *
   * Gathers all configuration variables for the Box.
   *
   * \param[in] modus_config The configuration object that sets all
   *                         initial conditions of the experiment.
   * \param[in] parameters Unused, but necessary because of templated
   *                       initialization
   */
  explicit BoxModus(Configuration modus_config,
                    const ExperimentParameters &parameters);

  /**
   * Generates initial state of the particles in the system according to
   * specified parameters: number of particles of each species, momentum
   * and coordinate space distributions. Subsequently makes the total
   * 3-momentum 0.
   *
   * \param[out] particles An empty list that gets filled up by this function
   * \param[in] parameters The initialization parameters of the box
   * \return The starting time of the simulation
   */
  double initial_conditions(Particles *particles,
                            const ExperimentParameters &parameters);

  /**
   * Enforces that all particles are inside the box at the beginning of
   * an event. It checks if the particles were placed correctly inside the box
   * at initialization and places them inside if they are not.
   *
   * \param[in] particles particles to check their position and possibly
   *            move it
   * \param[in] output_list output objects
   * \return The number of particles that were put back into the box
   *
   * In BoxModus if a particle crosses the wall of the box, it is
   * inserted from the opposite side. However these wall crossings are not
   * performed by this function but in the Experiment constructor when the
   * WallCrossActionsFinder are created. Wall crossings are written to
   * collision output: this is where OutputsList is used.
   */
  int impose_boundary_conditions(Particles *particles,
                                 const OutputsList &output_list = {});

  /// \copydoc smash::ModusDefault::create_grid
  Grid<GridOptions::PeriodicBoundaries> create_grid(
      const Particles &particles, double min_cell_length,
      double timestep_duration,
      CellSizeStrategy strategy = CellSizeStrategy::Optimal) const {
    return {{{0, 0, 0}, {length_, length_, length_}},
            particles,
            min_cell_length,
            timestep_duration,
            strategy};
  }

  /**
   * Creates GrandCanThermalizer. (Special Box implementation.)
   *
   * \param[in] conf configuration object
   * \return unique pointer to created thermalizer class
   */
  std::unique_ptr<GrandCanThermalizer> create_grandcan_thermalizer(
      Configuration &conf) const {
    const std::array<double, 3> lat_size = {length_, length_, length_};
    const std::array<double, 3> origin = {0., 0., 0.};
    const bool periodicity = true;
    return make_unique<GrandCanThermalizer>(conf, lat_size, origin,
                                            periodicity);
  }

  /// \copydoc smash::ModusDefault::max_timestep()
  double max_timestep(double max_transverse_distance_sqr) const {
    return 0.5 * std::sqrt(length_ * length_ - max_transverse_distance_sqr);
  }

  /// \return Length of the box
  double length() const { return length_; }
  /// \return equilibration time of the box
  double equilibration_time() const { return equilibration_time_; }
  /// \return whether the modus is box (also, trivially true)
  bool is_box() const { return true; }

 private:
  /// Initial momenta distribution: thermal or peaked momenta
  const BoxInitialCondition initial_condition_;
  /// Length of the cube's edge in fm/c
  const double length_;
  /// time after which output is written
  const double equilibration_time_;
  /// Temperature of the Box in GeV
  const double temperature_;
  /// Initial time of the box
  const double start_time_ = 0.;
  /**
   *  Whether to use a thermal initialization for all particles
   *  instead of specific numbers
   */
  const bool use_thermal_ = false;
  /**
   *  Baryon chemical potential for thermal initialization;
   *  only used if use_thermal_ is true
   */
  const double mub_;
  /**
   * Strange chemical potential for thermal initialization;
   * only used if use_thermal_ is true
   */
  const double mus_;
  /**
   * In case of thermal initialization: true -- account for resonance
   * spectral functions, while computing multiplicities and sampling masses,
   * false -- simply use pole masses.
   */
  const bool account_for_resonance_widths_;
  /**
   * Particle multiplicities at initialization;
   * required if use_thermal_ is false
   */
  const std::map<PdgCode, int> init_multipl_;
  /**
   * Average multiplicities in case of thermal initialization.
   * Saved to avoid recalculating at every event
   */
  std::map<PdgCode, double> average_multipl_;

  /**
   * Whether to insert a single high energy particle at the center of the
   * box (0,0,0). This particle will initially be moving along the x-axis.
   */
  const bool insert_jet_ = false;
  /**
   * Pdg of the particle to use as a jet; necessary if insert_jet_ is true,
   * unused otherwise
   */
  const PdgCode jet_pdg_;
  /**
   * Initial momentum of the jet particle; only used if insert_jet_ is true
   */
  const double jet_mom_;

  /**
   * \ingroup logging
   * Console output on startup of box specific parameters;
   * writes the initial state for the box to the output stream.
   *
   * \param[in] out The ostream into which to output
   * \param[in] m The BoxModus object to write into out
   */
  friend std::ostream &operator<<(std::ostream &out, const BoxModus &m);
};

}  // namespace smash

#endif  // SRC_INCLUDE_BOXMODUS_H_
