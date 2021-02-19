
/*
 *
 *    Copyright (c) 2021 Christian Holm Christensen
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_HEPMCINTERFACE_H_
#define SRC_INCLUDE_SMASH_HEPMCINTERFACE_H_

#include <memory>
#include <string>
#include <valarray>

#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenHeavyIon.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "action.h"
#include "forwarddeclarations.h"
#include "outputinterface.h"
#include "outputparameters.h"
namespace smash {

/**
 * \ingroup output
 *
 * \brief Base class for output handlers that need the HepMC3 structure
 *
 * This class can write the full event info or just the initial state
 * (i.e., beam particles) and final state (i.e., final state
 * particle).
 *
 * The class serves as a base class for output routines that utilizes
 * the HepMC event format (currently HepMcOutput and RivetOutput).
 *
 * A techincal point: We need to generate HepMC::GenParticle objects,
 * and we need to keep track of which HepMC::GenParticle corresponds
 * to which smash::ParticleData.  We therefor set up a map from the
 * smash::ParticleData identifier (integer) to HepMC::GenParticlePtr.
 * We use that map to keep track of the particles and interaction
 * points.

 * In case full event history:
 *
 * At each interaction we create a new vertex, and add the incoming
 * particles as "in" particles to that vertex.  We also set the
 * appropriate "state" of the incoming particles.  That is, if the
 * interaction corresponded to a decay, then the HepMC state is set to
 * two (2).  All other interactions are set to ten (10) plus the SMASH
 * interaction code (since there is no standard for these codes other
 * than 1: final state, 2: particle has decayed, 4: particle is beam
 * particle).
 *
 * For outgoing particles this is a bit different.  In case of elastic
 * scatterings, SMASH will keep the incoming particle around as an
 * outgoing particle. This is not how the HepMC event record is
 * invisioned.  Indeed, the particle has changed momentum and that
 * should be recorded in the event record.  In that case, we therefore
 * generate a new particle which we add as outgoing particle.
 *
 * All outgoing particles of a vertex have their initial status set to
 * one (1 - final state), but it can be changed later due to other
 * interactions.
 *
 * If the outgoing particle corresponds to an incoming particle, and
 * in particular if the incoming particle is a beam particle (or later
 * fragment thereof), we need to fragment the ion so that the outgoing
 * particle is dissociated from the incoming particle.  To that end,
 * we check if the outgoing particle was part of the beam particles.
 * If so, we remove the outgoing particle from the register of
 * identifiers that make up the beam particle, and create a new beam
 * particle.  This ensures that the HepMC event record is sound (in
 * most cases).
 */
class HepMcInterface : public OutputInterface {
 public:
  using AZ = std::pair<int, int>;
  /**
   * Create HepMC particle event in memory.
   *
   * \param[in] name    Name of output
   * \param[in] out_par Unused, needed since inhertied.
   * \param[in] total_N Total number of particles in both nuclei.
   * \param[in] proj_N  Number of particles in projectile.
   */
  HepMcInterface(const std::string& name, const OutputParameters& out_par,
                 const int total_N, const int proj_N);
  /**
   * Add the initial particles information of an event to the
   * central vertex.  Construct projectile and target particles with
   * nuclear pdg code if collider.
   *
   * \param[in] particles    Current list of all particles.
   * \param[in] event_number Current event number
   * \param[in] event        Event information
   * \throw std::runtime_error if nuclei with non-nucleon particle
   *                   (like hypernuclei) are tried to be constructed
   */
  void at_eventstart(const Particles& particles, const int event_number,
                     const EventInfo& event) override;
  /**
   * Writes collisions to event.
   *
   * \param[in] action an Action object containing incoming,
   *                   outgoing particles and type of interactions.
   * \param[in] density Unused, needed since inherited.
   */
  void at_interaction(const Action& action, const double density) override;
  /**
   * Add the final particles information of an event to the central vertex.
   * Store impact paramter and write event.
   *
   * \param[in] particles Current list of particles.
   * \param[in] event_number Number of event.
   * \param[in] event Event info, see \ref event_info
   */
  void at_eventend(const Particles& particles, const int32_t event_number,
                   const EventInfo& event) override;

 protected:
  /** HepMC status codes */
  enum Status {
    beam = 4,  // Beam particle
    fnal = 1,  // final state `final` is a reserved word
    dcy = 2,   // Decay
    off = 10
  };
  /** Type of mapping from SMASH ID to HepMC ID */
  using IdMap = std::map<int, HepMC3::GenParticlePtr>;
  /** Counter of collitions per incoming particle */
  using CollCounter = std::valarray<int>;
  /** Clear before an event */
  void clear();
  /** Convert SMASH process type to HepMC status */
  int get_status(const ProcessType& t) const;
  /**
   * Make an HepMC particle
   *
   * \param[in] pid     Particle type identifier
   * \param[in] status  Status code of particle
   * \param[in] mom     Four momentum of particle
   * \param[in] mass    Generator mass of particle
   *
   * \return A shared pointer to a HepMC::GenParticle object
   */
  HepMC3::GenParticlePtr make_gen(int pid, int status,
                                  const smash::FourVector& mom,
                                  double mass = -1);
  /**
   * Find particle in mapping or generate it.
   *
   * \param[in] p         ParticleData object
   * \param[in] status    HepMC status code
   *
   * \return r The existing or generated particle
   */
  HepMC3::GenParticlePtr make_register(const ParticleData& p,
                                       int status = Status::fnal);
  /**
   * Find particle in mapping or generate it.
   *
   * \param[in] p         ParticleData object
   * \param[in] status    HepMC status code
   * \param[in] force_new Create new particle even if could found
   *
   * \return r The existing or generated particle
   */
  HepMC3::GenParticlePtr find_or_make(const ParticleData& p,
                                      int status = Status::fnal,
                                      bool force_new = false);
  /**
   * Encode ion PDG
   *
   * \param[in] az Pair of Atomic weight and number
   *
   * \return PDG code of ion
   */
  int ion_pdg(const AZ& az) const;
  /** The event */
  HepMC3::GenEvent event_;
  /** The heavy-ion structure */
  HepMC3::GenHeavyIonPtr ion_;
  /** Dummy cross-section */
  HepMC3::GenCrossSectionPtr xs_;
  /** The interaction point */
  HepMC3::GenVertexPtr ip_;
  /** Mapping from ID to particle */
  IdMap map_;
  /**
   * Collision counter.  For each incoming particle we keep track of
   * how many colisions that particle partook in \f$ c_i\f$.  In the
   * end,
   *
   * \f{eqnarray}{
   N_{\mathrm{part}} &=$ \sum_i \begin{cases} 1 & c_i>0\\ 0 &c_i=0\end{case}\\
   N_{\mathrm{coll}} &=& \sum_i^{N_{\mathrm{pro}}} c_i\\
   \f}
   *
   */
  CollCounter coll_;
  /** counter of binary collisions (e.g., where both incoming
      particles are from the beams. */
  int ncoll_;
  /** counter of hard binary collisions (e.g., where both incoming
      particles are from the beams. */
  int ncoll_hard_;
  /**
   * Total number of nucleons in projectile and target,
   * needed for converting nuclei to single particles.
   */
  const int total_N_;
  /**
   * Total number of nucleons in projectile,
   * needed for converting nuclei to single particles.
   */
  const int proj_N_;
  /** Whether to only write final-state particles */
  bool only_final_;
  /** Extended particle information - not used? */
  bool part_extended_;
  /** Collsition (interaction) extended information - not used? */
  bool coll_extended_;
};
}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_RIVETOUTPUT_H_
