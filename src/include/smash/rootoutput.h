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
 * This class produces file Particles.root, which contains a
 * ROOT TTree. TTree contains information about particles from all SMASH events
 * comprising a simulation.
 * Output is happening in blocks. All particles in a block
 * are at the same time and in the same event. However, it is possible that
 * different blocks are at the same time and from the same event.
 * Particle information is stored in TBranches.
 * For each particle characteristic there is a separate branch.
 * Currently these are t,x,y,z (coordinates), p0,px,py,pz (4-momentum),
 * pdgcode - PDG code of the particle, characterizing its type,
 * charge - electric charge of the particle,
 * ev - event number in a given block,
 * tcounter - number of the output block in a given event,
 * npart - number of particles in the block,
 * test_p - number of testpartciles per particle,
 * modus_l - modus length,
 * current_t - time associated with the output block, in fm,
 * impact_b - impact parameter of the event,
 * empty_event - whether there was no interaction between the projectile and
 * the target,
 * E_kinetic_tot - total kinetic energy in the system,
 * E_fields_tot - total mean field energy * test_p,
 * E_total - sum of E_kinetic_tot and E_fields_tot.
 *
 * This class also produces file Collisions.root, organized in the same way,
 * with a few additional fields:
 * nin and nout - characterize number of incoming and outgoing particles in the
 * reaction, with nin + nout = npart,
 * weight - an action weight, whose meaning depends on the type of action: For
 * collisions it is the total cross section, for decays it is the total decay
 * width and for dilepton decays it is the shining weight.
 *
 * If "Collisions:" section is present, then in addition to a file
 * Particles.root with particles TTree, another file Collisions.root is
 * created. It contains information about each collision, written as one leaf:
 * nin, nout - number of incoming and outgoing particles, ev - event number,
 * weight - total weight of the collision (wgt), partial_weight - partial
 * weight of the collision (par_wgt), (t,x,y,z),
 * (p0,px,py,pz) - arrays of dimension nin+nout
 * that contain coordinates and momenta.
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
  std::vector<double> p0_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> px_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> py_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> pz_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> t_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> x_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> y_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> z_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> formation_time_ =
      std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> xsec_factor_ = std::vector<double>(max_buffer_size_, 0.0);
  std::vector<double> time_last_coll_ =
      std::vector<double>(max_buffer_size_, 0.0);
  std::vector<int> pdgcode_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> charge_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> coll_per_part_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> proc_id_origin_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> proc_type_origin_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> pdg_mother1_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> pdg_mother2_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> baryon_number_ = std::vector<int>(max_buffer_size_, 0);
  std::vector<int> strangeness_ = std::vector<int>(max_buffer_size_, 0);
  int npart_, tcounter_, ev_, ens_, nin_, nout_, test_p_;
  double wgt_, par_wgt_, impact_b_, modus_l_, current_t_;
  double E_kinetic_tot_, E_fields_tot_, E_tot_;
  bool empty_event_;
  //@}

  /// Option to write collisions tree.
  bool write_collisions_;

  /// Option to write particles tree.
  bool write_particles_;

  /// Option to write particles tree for initial conditions
  bool write_initial_conditions_;

  /// Print only final particles in the event, no intermediate output.
  OutputOnlyFinal particles_only_final_;

  /**
   * ROOT file cannot be read if it was not properly closed and finalized.
   * It can happen that SMASH simulation crashed and ROOT file was not closed.
   * To save results of simulation in such case, "AutoSave" is applied every N
   * events. If multiple ensembles are used, N is divided by the number of
   * ensembles per event. This makes sense especially in case of a large number
   * of ensembles and, at the same time, it ensures that all ensembles belonging
   * to the same event are saved. The autosave_frequency_ sets this N (at the
   * moment N=1000 hard-coded). Note that "AutoSave" operation is very
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
