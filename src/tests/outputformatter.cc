/*
 *
 *    Copyright (c) 2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/outputformatter.h"

#include "setup.h"

using namespace smash;

TEST_CATCH(invalid_quanitity, std::invalid_argument) {
  std::vector<std::string> invalid_quantities = {"gibberish"};
  OutputFormatter<ASCII> outputformatter(invalid_quantities);
}

TEST(valid_line_maker) {
  ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
      "π⁺ 0.1380 0      -  211\n");

  PdgCode pdg = 0x211;
  ParticleData p{ParticleType::find(pdg)};

  std::vector<std::string> valid_quantities = {"t",
                                               "x",
                                               "y",
                                               "z",
                                               "mass",
                                               "p0",
                                               "px",
                                               "py",
                                               "pz",
                                               "pdgcode",
                                               "ID",
                                               "charge",
                                               "ncoll",
                                               "form_time",
                                               "xsecfac",
                                               "proc_id_origin",
                                               "proc_type_origin",
                                               "t_last_coll",
                                               "pdg_mother1",
                                               "pdg_mother2",
                                               "baryon_number",
                                               "strangeness",
                                               "spin"};

  OutputFormatter<ASCII> formatter(valid_quantities);
  std::string correct_line;

  for (int i = 0; i < 4; ++i) {
    correct_line += std::to_string(p.position()[i]) + ",";
  }
  correct_line += std::to_string(p.effective_mass()) + ",";
  for (int i = 0; i < 4; ++i) {
    correct_line += std::to_string(p.momentum()[i]) + ",";
  }
  correct_line += p.pdgcode().string() + ",";
  correct_line += std::to_string(p.id()) + ",";
  correct_line += std::to_string(p.type().charge()) + ",";
  correct_line += std::to_string(p.get_history().collisions_per_particle) + ",";
  correct_line += std::to_string(p.formation_time()) + ",";
  correct_line += std::to_string(p.get_history().id_process) + ",";
  correct_line +=
      std::to_string(static_cast<int>(p.get_history().process_type)) + ",";
  correct_line += std::to_string(p.get_history().time_last_collision) + ",";
  correct_line += p.get_history().p1.string() + ",";
  correct_line += p.get_history().p2.string() + ",";
  correct_line += std::to_string(p.pdgcode().baryon_number()) + ",";
  correct_line += std::to_string(p.pdgcode().strangeness()) + ",";
  correct_line += std::to_string(p.spin_projection()) + ",";
  correct_line.pop_back();

  std::string tested_line = formatter.data_line(p);

  VERIFY(correct_line == tested_line);
}
