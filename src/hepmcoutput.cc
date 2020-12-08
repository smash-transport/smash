/*
 *
 *    Copyright (c) 2014-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/hepmcoutput.h"

#include "HepMC3/Print.h"

namespace smash {

/*!\Userguide
 * \page hepmc_output_user_guide_ HepMC Output
 * SMASH HepMC output is an implementation of the HepMC3 Event Record Library.
 * The aim is to provide a format compatible with other frameworks like Rivet
 * (https://rivet.hepforge.org). For resources regarding HepMC, see
 * http://hepmc.web.cern.ch/hepmc/ and https://arxiv.org/abs/1912.08005.
 *
 * Producing HepMC output requires HepMC3 to be installed. Download the tarball
 * from http://hepmc.web.cern.ch/hepmc/.
 *
 * \section hepmc_output_user_guide_format_ ASCII HepMC Format
 *
 * HepMC generaly structures each event into particles and vertices connecting
 * them, basically storing a graph of the event. Since the purpose of this
 * output is to only provide a particle list of the final state produced the
 * output format is adapted accordingly for the SMASH implementation: Only one
 * central vertex is used. All intial state particles are incoming particles and
 * all final state particles are outgoing particles of this vertex. Scatterings
 * happening during the SMASH event are not recorded. For the collider modus,
 * the intial state particles are combined into two single colliding nucleus
 * "particles" with a nuclear pdg code.
 *
 * The output is written in Asciiv3, the HepMC3 native plain text format. See
 * https://arxiv.org/abs/1912.08005 for documentation of the format.
 *
 * \note Since some HepMC readers (e.g. Rivet) need a value for the
 * nuclei-nuclei cross section, a dummy cross section of 1.0 is written to the
 * output. Furthermore, if you use Fermi motion and want to read in the HepMC
 * ouput into Rivet, you need to disable the check for the beam particle
 * energies with the \key --ignore-beams option.
 */

const int HepMcOutput::status_code_for_beam_particles = 4;
const int HepMcOutput::status_code_for_final_particles = 1;

HepMcOutput::HepMcOutput(const bf::path &path, std::string name,
                         const OutputParameters & /*out_par*/,
                         const int total_N, const int proj_N)
    : OutputInterface(name),
      filename_(path / (name + ".asciiv3")),
      total_N_(total_N),
      proj_N_(proj_N) {
  filename_unfinished_ = filename_;
  filename_unfinished_ += +".unfinished";
  output_file_ =
      make_unique<HepMC3::WriterAscii>(filename_unfinished_.string());
}

HepMcOutput::~HepMcOutput() { bf::rename(filename_unfinished_, filename_); }

int HepMcOutput::construct_nuclear_pdg_code(int na, int nz) const {
  const int pdg_nuclear_code_prefix = 10 * 1E8;
  // Hypernuclei not supported here
  const int pdg_nuclear_code_lambda = 0 * 1E7;
  const int pdg_nuclear_code_charge = nz * 1E4;
  const int pdg_nuclear_code_baryon = na * 1E1;
  // SMASH does not do isomers
  const int pdg_nuclear_code_isomer = 0 * 1E0;

  return pdg_nuclear_code_prefix + pdg_nuclear_code_lambda +
         pdg_nuclear_code_charge + pdg_nuclear_code_baryon +
         pdg_nuclear_code_isomer;
}

void HepMcOutput::at_eventstart(const Particles &particles,
                                const int event_number, const EventInfo &) {
  current_event_ =
      make_unique<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);

  /* Rivet needs a value for the cross section, but SMASH only knows about
   * nucleon cross sections, not about cross sections between nuclei,
   * so a dummy is used. */
  std::shared_ptr<HepMC3::GenCrossSection> cross_section =
      std::make_shared<HepMC3::GenCrossSection>();
  current_event_->add_attribute("GenCrossSection", cross_section);
  const double dummy_xs = 1.0;
  cross_section->set_cross_section(dummy_xs, dummy_xs);

  current_event_->set_event_number(event_number);
  vertex_ = std::make_shared<HepMC3::GenVertex>();
  current_event_->add_vertex(vertex_);

  if (proj_N_ > 0) {
    // Collider modus: Construct and write projectile and target as two intial
    // particles
    int targ_N = total_N_ - proj_N_;
    int proj_Z = 0;
    int targ_Z = 0;
    FourVector total_mom_proj;
    FourVector total_mom_targ;
    for (const ParticleData &data : particles) {
      if (data.id() < proj_N_) {
        total_mom_proj += data.momentum();
        if (data.is_proton()) {
          proj_Z += 1;
        } else if (!data.is_neutron()) {
          throw std::invalid_argument(
              "HepMC output in SMASH only supports colliding nuclei consisting "
              "of p and n.");
        }
      } else {
        total_mom_targ += data.momentum();
        if (data.is_proton()) {
          targ_Z += 1;
        } else if (!data.is_neutron()) {
          throw std::invalid_argument(
              "HepMC output in SMASH only supports colliding nuclei consisting "
              "of p and n.");
        }
      }
    }

    const int proj_nuclear_pdg_code =
        construct_nuclear_pdg_code(proj_N_, proj_Z);
    const int targ_nuclear_pdg_code =
        construct_nuclear_pdg_code(targ_N, targ_Z);

    HepMC3::GenParticlePtr projectile_p = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector(total_mom_proj.x1(), total_mom_proj.x2(),
                           total_mom_proj.x3(), total_mom_proj.x0()),
        proj_nuclear_pdg_code, status_code_for_beam_particles);
    vertex_->add_particle_in(projectile_p);
    HepMC3::GenParticlePtr target_p = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector(total_mom_targ.x1(), total_mom_targ.x2(),
                           total_mom_targ.x3(), total_mom_targ.x0()),
        targ_nuclear_pdg_code, status_code_for_beam_particles);
    vertex_->add_particle_in(target_p);

  } else {
    // Other modi (not collider): Write all inital particles into output
    for (const ParticleData &data : particles) {
      const FourVector mom = data.momentum();
      HepMC3::GenParticlePtr p = std::make_shared<HepMC3::GenParticle>(
          HepMC3::FourVector(mom.x1(), mom.x2(), mom.x3(), mom.x0()),
          data.pdgcode().get_decimal(), status_code_for_beam_particles);
      vertex_->add_particle_in(p);
    }
  }
}

void HepMcOutput::at_eventend(const Particles &particles,
                              const int32_t /*event_number*/,
                              const EventInfo &event) {
  // Set heavy ion attribute, only the impact parameter is known
  std::shared_ptr<HepMC3::GenHeavyIon> heavy_ion =
      std::make_shared<HepMC3::GenHeavyIon>();
  current_event_->add_attribute("GenHeavyIon", heavy_ion);
  // Impact paramter has to converted to mm (x1E-12), since fm not a supported
  // unit in HepMC, -1(.0) are placeholders
  heavy_ion->set(-1, -1, -1, -1, -1, -1, -1, -1, -1,
                 event.impact_parameter * 1E-12, -1.0, -1.0, -1.0, -1.0, -1.0);

  for (const ParticleData &data : particles) {
    const FourVector mom = data.momentum();
    HepMC3::GenParticlePtr p = std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector(mom.x1(), mom.x2(), mom.x3(), mom.x0()),
        data.pdgcode().get_decimal(), status_code_for_final_particles);
    vertex_->add_particle_out(p);
  }

  output_file_->write_event(*current_event_);
}

}  // namespace smash
