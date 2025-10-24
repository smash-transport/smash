/*
 *
 *    Copyright (c) 2014-2020,2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_ROOTOUTPUT_H_
#define SRC_INCLUDE_SMASH_ROOTOUTPUT_H_

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "configuration.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "outputparameters.h"

namespace smash {
class Particles;

/**
 * \ingroup output
 *
 * \brief <h2> SMASH output to ROOT file </h2>
 *
 * SMASH supports ROOT output as an option (see http://root.cern.ch).
 * The ROOT framework needs to be installed before building SMASH, otherwise
 * ROOT support will be disabled.
 *
 * This class produces file Particles.root, which contains a ROOT TTree. The
 * TTree contains information about particles from all SMASH events comprising a
 * simulation.
 * Output is performed in blocks. All particles in a block are at the same time
 * and in the same event. Due to the division of the output for large buffers,
 * it is possible that different blocks are at the same time and from the same
 * event.
 * Information is stored in TBranches. Each output block contains TBranches with
 * overall information about the event and the output block:
 * ev - event number in a given block,
 * ens - ensemble number in a given block,
 * tcounter - number of the output block in a given event,
 * npart - number of particles in the block,
 * test_p - number of testpartciles per particle,
 * modus_l - modus length,
 * current_t - time associated with the output block, in fm,
 * impact_b - impact parameter of the event,
 * empty_event - whether the projectile and the target did not collide/interact.
 * Then, each output block contains information about all particles in the
 * block. For each particle characteristic there is a separate TBranch. These
 * TBranches contain arrays of the following quantities:
 * id - a unique integer identifier of the particle,
 * pdgcode - PDG code of the particle, characterizing its type,
 * charge - electric charge of the particle,
 * formation_time - time when the particle was created,
 * time_last_collision - time of the last collision before the output time,
 * p0, px, py, pz - components of the 4-momentum,
 * t, x, y, z - coordinates.
 * In each array, the i-th entry corresponds to the i-th particle.
 * Finally, each output contains information about some of the bulk properties
 * of the system:
 * E_kinetic_tot - total kinetic energy in the system,
 * E_fields_tot - total mean field energy * test_p,
 * E_total - sum of E_kinetic_tot and E_fields_tot.
 *
 * It is possible to request an extended ROOT output, which in addition to these
 * characteristics also contains
 * ncoll - number of collisions the particle has undergone,
 * xsecfac - cross section scaling factor if the particles are not yet fully
 * formed at the time of interaction, the cross section for the underlying
 * process is scaled down by the cross section scaling factor),
 * proc_id_origin - ID of the process of the particle's last interaction,
 * proc_type_origin - type of the last process the particle has undergone (the
 * possible process types are listed in \ref doxypage_output_process_types),
 * pdg_mother1 - PDG code of the 1st mother particle (0 in case the particle is
 * sampled in a thermal bubble; it is not updated by elastic scatterings);
 * pdg_mother2 - PDG code of the 2nd mother particle (0 in case the particle
 * originates from the decay of a resonance or the appearance of a thermal
 * bubble; in the former case, \key pdg_mother1 is the PDG code of the mother
 * resonance; it is not updated by elastic scatterings),
 * baryon_number - baryon number of the particle,
 * strangeness - strangeness number of the particle.
 *
 * If "Collisions:" section is present in the config file, then in addition to a
 * file Particles.root with particles TTree, another file Collisions.root is
 * created. Collisions.root contains information about each collision, written
 * as one leaf. Similarly to Particles.root, this includes global information
 * about the event:
 * ev - event number,
 * ens - ensemble number,
 * and information about particles participating in a collision:
 * id - a unique integer identifier of the particle,
 * pdgcode - PDG code of the particle, characterizing its type,
 * charge - electric charge of the particle,
 * formation_time - time when the particle was created,
 * time_last_collision - time of the last collision before the output time,
 * p0, px, py, pz - components of the 4-momentum,
 * t, x, y, z - coordinates,
 * where all arrays (id, pdgcode, charge, formation_time, time_last_collision,
 * p0, px, py, pz, t, x, y, z) are now of dimension nin+nout = npart.
 * Beyond these, Collisions.root contains a few additional fields:
 * nin and nout - the number of incoming and outgoing particles in the reaction,
 * with nin + nout = npart,
 * wgt - an action weight, whose meaning depends on the type of action: for
 * collisions it is the total cross section, for decays it is the total decay
 * width, and for dilepton decays it is the shining weight,
 * par_wgt - partial weight of the collision.
 */
class RootOutput : public OutputInterface {
 public:
  /**
   * Construct ROOT output.
   *
   * \param[in] path Output path.
   * \param[in] name Name of the ouput.
   * \param[in] out_par A structure containing parameters of the output.
   */
  RootOutput(const std::filesystem::path &path, const std::string &name,
             const OutputParameters &out_par);

  /// Destructor
  ~RootOutput();

  /**
   * update event number and writes intermediate particles to a tree.
   * \param[in] particles Particles to be written to output.
   * \param[in] event_label Numbers of event and ensemble.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventstart(const Particles &particles, const EventLabel &event_label,
                     const EventInfo &event) override;
  /**
   * update event number and impact parameter,
   * and writes intermediate particles to a tree.
   * \param[in] particles Particles to be written to output.
   * \param[in] event_label Numbers of event and ensemble.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventend(const Particles &particles, const EventLabel &event_label,
                   const EventInfo &event) override;
  /**
   * Writes intermediate particles to a tree defined by treename,
   * if it is allowed (i.e., particles_only_final_ is No).
   * \param[in] particles Particles to be written to output.
   * \param[in] clock Unused, needed since inherited.
   * \param[in] dens_param Unused, needed since inherited.
   * \param[in] event_label Numbers of event and ensemble.
   * \param[in] event Event info, see \ref event_info
   */
  void at_intermediate_time(const Particles &particles,
                            const std::unique_ptr<Clock> &clock,
                            const DensityParameters &dens_param,
                            const EventLabel &event_label,
                            const EventInfo &event) override;
  /**
   * Writes collisions to a tree defined by treename.
   * \param[in] action an Action object containing incoming, outgoing particles
   *            and type of interactions.
   * \param[in] density Unused, needed since inherited.
   */
  void at_interaction(const Action &action, const double density) override;

 private:
  /// Filename of output
  const std::filesystem::path filename_;
  /// Filename of output as long as simulation is still running.
  std::filesystem::path filename_unfinished_;
  /// Pointer to root output file.
  std::unique_ptr<TFile> root_out_file_;
  /**
   * TTree for particles output.
   *
   * TFile takes ownership of all TTrees.
   * That's why TTree is not a unique pointer.
   */
  TTree *particles_tree_;
  /**
   * TTree for collision output.
   *
   * TFile takes ownership of all TTrees.
   * That's why TTree is not a unique pointer.
   */
  TTree *collisions_tree_;
  /**
   * Writes particles to a tree defined by treename.
   * \param[in] particles Particles or ParticleList to be written to output.
   */
  template <typename T>
  void particles_to_tree(T &particles);
  /**
   * Writes collisions to a tree defined by treename.
   * \param[in] incoming Incoming particles to be written to output.
   * \param[in] outgoing Outgoing particles to be written to output.
   * \param[in] weight Total weight of the collision.
   * \param[in] partial_weight Partial weight of the collision
   */
  void collisions_to_tree(const ParticleList &incoming,
                          const ParticleList &outgoing, const double weight,
                          const double partial_weight);
  /// Number of output in a given event.
  int output_counter_ = 0;
  /// Number of current event.
  int current_event_ = 0;
  /// Number of current ensemble.
  int current_ensemble_ = 0;

  /**
   * Maximal buffer size.
   * When the number of particles N exceeds the buffer size B, data is flushed
   * to the ROOT file every B particles. This creates ceil(N/B) entries in the
   * ROOT Tree at every output.
   */
  static const int max_buffer_size_;

  /** @name Buffer for filling TTree
   * See class documentation for definitions.
   */
  //@{
  /// Property that is written to ROOT output.
  int ev_{};
  int ens_{};
  int tcounter_{};
  int npart_{};
  int test_p_{};
  double modus_l_{};
  double current_t_{};
  double impact_b_{};
  bool empty_event_{};
  std::vector<int> id_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> pdgcode_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> charge_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<double> formation_time_ =
    std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> time_last_collision_ =
    std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> p0_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> px_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> py_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> pz_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> t_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> x_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> y_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> z_ = std::vector<double>(max_buffer_size_, 0.0);
  double E_kinetic_tot_{};
  double E_fields_tot_{};
  double E_tot_{};
  std::vector<int> coll_per_part_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<double> xsec_factor_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<int> proc_id_origin_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> proc_type_origin_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> pdg_mother1_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> pdg_mother2_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> baryon_number_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> strangeness_ = std::vector<int>(max_buffer_size_, 0);
  int nin_{};
  int nout_{};
  double wgt_{};
  double par_wgt_{};
  //@}

  /// Option to write collisions tree.
  bool write_collisions_;

  /// Option to write particles tree.
  bool write_particles_;

  /// Option to write particles tree for initial conditions
  bool write_initial_conditions_;

  /// Print only final particles in the event (no intermediate output).
  OutputOnlyFinal particles_only_final_;

  /**
   * ROOT file cannot be read if it was not properly closed and finalized.
   * It can happen that SMASH simulation crashed and ROOT file was not closed.
   * To save results of simulation in such case, "AutoSave" is applied every N
   * events. If multiple ensembles are used, N is divided by the number of
   * ensembles per event. This makes sense especially in case of a large number
   * of ensembles and, at the same time, it ensures that all ensembles belonging
   * to the same event are saved. The autosave_frequency_ sets this N (at the
   * moment N=1000 is hard-coded). Note that "AutoSave" operation is very
   * time-consuming, so the autosave frequency is always a compromise between
   * safety and speed.
   */
  int autosave_frequency_ = -1;

  /// Whether extended particle output is on
  const bool part_extended_;
  /// Whether extended collisions output is on
  const bool coll_extended_;
  /// Whether extended ic output is on
  const bool ic_extended_;

  /**
   * Basic initialization routine, creating the TTree objects
   * for particles and collisions.
   */
  void init_trees();
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_ROOTOUTPUT_H_
